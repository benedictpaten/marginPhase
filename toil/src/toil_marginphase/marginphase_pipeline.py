#!/usr/bin/env python2.7
from __future__ import print_function

import argparse
import os
import multiprocessing
import subprocess
import sys
import textwrap
import tarfile
from urlparse import urlparse
import math
import shutil
import glob
import time
import datetime
import logging

import yaml
from bd2k.util.files import mkdir_p
from bd2k.util.processes import which
from toil.job import Job
from toil.lib.docker import dockerCall, dockerCheckOutput, _fixPermissions
from toil_lib import require, UserError
from toil_lib.files import tarball_files, copy_files
from toil_lib.jobs import map_job
from toil_lib.urls import download_url, s3am_upload

# formats
SCHEMES = ('http', 'file', 's3', 'ftp')

# filenames
DEFAULT_CONFIG_NAME = 'config-toil-marginphase.yaml'
DEFAULT_MANIFEST_NAME = 'manifest-toil-marginphase.tsv'

# docker images
DOCKER_SAMTOOLS = "quay.io/ucsc_cgl/samtools:1.3--256539928ea162949d8a65ca5c79a72ef557ce7c"


def parse_samples(path_to_manifest):
    """
    Parses samples, specified in either a manifest or listed with --samples

    :param str path_to_manifest: Path to configuration file
    :return: Samples and their attributes as defined in the manifest
    :rtype: list[list]
    """

    samples = []
    with open(path_to_manifest, 'r') as f:
        for line in f.readlines():
            if line.isspace() or line.startswith('#'):
                continue
            sample = line.strip().split('\t')

            # validate structure
            if len(sample) != 2:
                raise UserError('Bad manifest format! Expected 2 tab-separated columns, got: {}'.format(sample))

            # extract sample parts
            uuid, url = sample

            # validation?

            sample = [uuid, url]
            samples.append(sample)
    return samples


def prepare_input(job, sample, config):

    # job prep
    config = argparse.Namespace(**vars(config))
    uuid, url = sample
    config.uuid = uuid
    work_dir = job.fileStore.getLocalTempDir()
    start = time.time()
    job.fileStore.logToMaster("START:{}:{}".format(config.uuid, datetime.datetime.now()))

    # todo global resource estimation
    #config.cores = min(config.maxCores, multiprocessing.cpu_count())
    #config.memory
    #config.disk


    download_url(job, url=url, work_dir=work_dir)
    filename = os.path.basename(url)
    data_bam_location = os.path.join("/data", filename)

    # index the bam
    docker_params = ["index", data_bam_location]
    job.fileStore.logToMaster("Running {} with parameters: {}".format(DOCKER_SAMTOOLS, docker_params))
    dockerCall(job, tool=DOCKER_SAMTOOLS, workDir=work_dir, parameters=docker_params)
    _fixPermissions(DOCKER_SAMTOOLS, work_dir)  # todo tmpfix

    # sanity check
    workdir_bai_location = os.path.join(work_dir, filename + ".bai")
    if not os.path.isfile(workdir_bai_location):
        raise UserError("BAM index file not created for {}: {}".format(filename, workdir_bai_location))

    # get start and end location
    get_idx_cmd = [
        ["samtools", "view", data_bam_location],
        ["head", "-1"],
        ["awk", "'{print $4}'"]
    ]
    start_idx_str = dockerCheckOutput(job, tool=DOCKER_SAMTOOLS, workDir=work_dir, parameters=get_idx_cmd)
    get_idx_cmd[1][0] = "tail"
    end_idx_str = dockerCheckOutput(job, tool=DOCKER_SAMTOOLS, workDir=work_dir, parameters=get_idx_cmd)
    job.fileStore.logToMaster("{}: start_pos:{}, end_pos:{}".format(config.uuid, start_idx_str, end_idx_str))
    start_idx = int(start_idx_str) - 1
    end_idx = int(end_idx_str) + 1

    # get reads from positions
    positions = list()
    idx = start_idx
    while idx < end_idx:
        chunk_start = idx - (config.partition_margin if idx != start_idx else 0)
        idx += config.partition_size
        chunk_end = idx + config.partition_margin
        positions.append([chunk_start, chunk_end])

    # enqueue jobs
    idx = 0
    return_values = list()
    for position in positions:
        margin_phase_job = job.addChildJobFn(run_margin_phase, config, "file_id", idx)
        return_values.append(margin_phase_job.rv())

    # enqueue consolidation
    job.addFollowOnJobFn(consolidate_output, config, return_values)

    # log
    _log_time(job, "prepare_input", start, config.uuid)


def run_margin_phase(job, config, chunk_file_id, chunk_idx):
    job.fileStore.logToMaster("run_margin_phase : {} : chunk_{}".format(config.uuid, chunk_idx))
    chunk_identifier = "{}.{}".format(config.uuid, chunk_idx)
    return chunk_identifier

def consolidate_output(job, config, output_list):
    job.fileStore.logToMaster("consolidate_output : {}".format(config.uuid))
    for elem in output_list:
        job.fileStore.logToMaster(str(elem))



def _log_time(job, function_name, start_time, sample_identifier=''):
    job.fileStore.logToMaster("TIME:{}:{}:{}".format(sample_identifier, function_name, int(time.time() - start_time)))


def _get_default_docker_params(work_dir):
    return ['--rm','--log-driver','none', '-v', '{}:/data'.format(work_dir)]

def generate_config():
    return textwrap.dedent("""
        # MarginPhase Pipeline configuration file
        # This configuration file is formatted in YAML. Simply write the value (at least one space) after the colon.
        # Edit the values in this configuration file and then rerun the pipeline: "toil-marginphase run"
        #
        # URLs can take the form: http://, ftp://, file://, s3://
        # Local inputs follow the URL convention: file:///full/path/to/input
        # S3 URLs follow the convention: s3://bucket/directory/file.txt
        #
        # Comments (beginning with #) do not need to be removed. Optional parameters left blank are treated as false.
        ##############################################################################################################
        # Required: URL {scheme} to reference contig
        # TODO: this should be altered to accept multiple references
        reference: file:///your/reference/here.fa

        # Required: Output location of sample. Can be full path to a directory or an s3:// URL
        # Warning: S3 buckets must exist prior to upload or it will fail.
        output-dir:

        # Required: Size of each bam partition
        partition-size: 2000000

        # Required: Margin to apply on each partition
        partition-margin: 5000

    """.format(scheme=[x + '://' for x in SCHEMES])[1:])


def generate_manifest():
    return textwrap.dedent("""
        #   Edit this manifest to include information pertaining to each sample to be run.
        #
        #   There are 4 tab-separated columns: filetype, paired/single, UUID, URL(s) to sample
        #
        #   UUID        This should be a unique identifier for the sample to be processed
        #   filetype    Filetype of the sample. Options: "tar", "fq", "bam"
        #   paired      Indicates whether the data is paired or single-ended. Options:  "paired" or "single"
        #   URL         A URL {scheme} pointing to the sample or a full path to a directory
        #
        #   If samples are in .tar or .tar.gz format, the files must be in the root of the tarball (not in a folder)
        #   If samples are submitted as multiple fastq files, provide comma separated URLs.
        #   If samples are submitted as multiple bam files, provide comma separated URLs.
        #   Samples must have the same extension - do not mix and match gzip and non-gzipped sample pairs.
        #
        #   Samples in "tar" format must have one of these extensions: .tar .tar.gz
        #   Files in the tarballs must have one of these extensions: .fastq .fq
        #   Paired files in tarballs must follow name conventions ending in an R1/R2 or _1/_2
        #
        #   Samples in "fq" format must have one of these extensions: fastq.gz, fastq, fq.gz, fq
        #   Paired samples in "fq" format must follow name conventions ending in an R1/R2 or _1/_2.  Ordering in manifest does not matter.
        #
        #   Samples in "bam" format will be converted to fastq with picard's SamToFastq utility.
        #   The paired/single parameter will determine whether the pipleine requests single or paired output from picard
        #
        #   When multiple samples are submitted, they will be concatinated into a single (or pair of) fastq file(s) before alignment
        #
        #   If a full path to a directory is provided for a sample, every file inside needs to be a fastq(.gz).
        #   Do not have other superfluous files / directories inside or the pipeline will complain.
        #
        #   Examples of several combinations are provided below. Lines beginning with # are ignored.
        #
        #   UUID_1  file:///path/to/file.bam
        #   UUID_2  s3://path/to/file.bam
        #
        #   Place your samples below, one per line.
        """.format(scheme=[x + '://' for x in SCHEMES])[1:])


def generate_file(file_path, generate_func):
    """
    Checks file existance, generates file, and provides message

    :param str file_path: File location to generate file
    :param function generate_func: Function used to generate file
    """
    require(not os.path.exists(file_path), file_path + ' already exists!')
    with open(file_path, 'w') as f:
        f.write(generate_func())
    print('\t{} has been generated in the current working directory.'.format(os.path.basename(file_path)))


def main():
    """
    Computational Genomics Lab, Genomics Institute, UC Santa Cruz
    MarginPhase pipeline

    =======================================
    Dependencies
    Curl:       apt-get install curl
    Docker:     wget -qO- https://get.docker.com/ | sh
    Toil:       pip install toil
    Boto:       pip install boto (OPTIONAL)
    """

    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(dest='command')

    # Generate subparsers
    subparsers.add_parser('generate-config', help='Generates an editable config in the current working directory.')
    subparsers.add_parser('generate-manifest', help='Generates an editable manifest in the current working directory.')
    subparsers.add_parser('generate', help='Generates a config and manifest in the current working directory.')

    # Run subparser
    parser_run = subparsers.add_parser('run', help='Runs the MarginPhase pipeline')
    group = parser_run.add_mutually_exclusive_group()
    parser_run.add_argument('--config', default=DEFAULT_CONFIG_NAME, type=str,
                            help='Path to the (filled in) config file, generated with "generate-config". '
                                 '\nDefault value: "%(default)s"')
    group.add_argument('--manifest', default=DEFAULT_MANIFEST_NAME, type=str,
                       help='Path to the (filled in) manifest file, generated with "generate-manifest". '
                            '\nDefault value: "%(default)s"')

    # If no arguments provided, print full help menu
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # Add Toil options
    Job.Runner.addToilOptions(parser_run)
    args = parser.parse_args()

    # Parse subparsers related to generation of config and manifest
    cwd = os.getcwd()
    if args.command == 'generate-config' or args.command == 'generate':
        generate_file(os.path.join(cwd, DEFAULT_CONFIG_NAME), generate_config)
    if args.command == 'generate-manifest' or args.command == 'generate':
        generate_file(os.path.join(cwd, DEFAULT_MANIFEST_NAME), generate_manifest)

    # Pipeline execution
    elif args.command == 'run':
        # sanity check
        require(os.path.exists(args.config), '{} not found. Please run '
                                             '"toil-marginphase generate-config"'.format(args.config))
        require(os.path.exists(args.manifest), '{} not found and no samples provided. Please '
                                               'run "toil-marginphase generate-manifest"'.format(args.manifest))

        # get samples
        samples = parse_samples(path_to_manifest=args.manifest)

        # Parse config
        parsed_config = {x.replace('-', '_'): y for x, y in yaml.load(open(args.config).read()).iteritems()}
        config = argparse.Namespace(**parsed_config)
        config.maxCores = int(args.maxCores) if args.maxCores else sys.maxint
        config.maxDisk = int(args.maxDisk) if args.maxDisk else sys.maxint
        config.maxMemory = args.maxMemory if args.maxMemory else str(sys.maxint)
        config.memory = args.maxMemory if args.maxMemory else str(sys.maxint)

        # Config sanity checks
        require(config.output_dir, 'No output location specified')
        if urlparse(config.output_dir).scheme != "s3":
            mkdir_p(config.output_dir)
        if not config.output_dir.endswith('/'):
            config.output_dir += '/'
        require(config.reference, 'No reference contig specified')
        require(config.partition_size, "Configuration parameter partition_size is required")
        require(config.partition_margin, "Configuration parameter partition_margin is required")

        # Program checks
        for program in ['docker']:
            require(next(which(program), None), program + ' must be installed on every node.'.format(program))

        # Start the workflow
        Job.Runner.startToil(Job.wrapJobFn(map_job, prepare_input, samples, config), args)


if __name__ == '__main__':
    try:
        main()
    except UserError as e:
        print(e.message, file=sys.stderr)
        sys.exit(1)
