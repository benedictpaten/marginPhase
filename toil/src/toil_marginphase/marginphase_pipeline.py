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
from contextlib import closing

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
DOCKER_MARGIN_PHASE = "quay.io/ucsc_cgl/margin_phase:latest"

# resource
MP_CPU = 2
MP_MEM_BAM_FACTOR = 1024 #todo account for learning iterations
MP_MEM_REF_FACTOR = 2
MP_DSK_BAM_FACTOR = 2.5
MP_DSK_REF_FACTOR = 1.1


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

    # download references
    #ref genome
    download_url(job, url=config.reference_contig, work_dir=work_dir)
    ref_genome_filename = os.path.basename(config.reference_contig)
    ref_genome_fileid = job.fileStore.writeGlobalFile(os.path.join(work_dir, ref_genome_filename))
    config.reference_genome_fileid = ref_genome_fileid
    ref_genome_size = os.stat(os.path.join(work_dir, ref_genome_filename)).st_size
    #ref vcf
    download_url(job, url=config.reference_vcf, work_dir=work_dir)
    ref_vcf_filename = os.path.basename(config.reference_vcf)
    ref_vcf_fileid = job.fileStore.writeGlobalFile(os.path.join(work_dir, ref_vcf_filename))
    config.reference_vcf_fileid = ref_vcf_fileid
    #params
    download_url(job, url=config.params, work_dir=work_dir)
    params_filename = os.path.basename(config.params)
    params_fileid = job.fileStore.writeGlobalFile(os.path.join(work_dir, params_filename))
    config.params_fileid = params_fileid


    # download bam
    download_url(job, url=url, work_dir=work_dir)
    bam_filename = os.path.basename(url)
    data_bam_location = os.path.join("/data", bam_filename)

    # index the bam
    docker_params = ["index", data_bam_location]
    job.fileStore.logToMaster("Running {} with parameters: {}".format(DOCKER_SAMTOOLS, docker_params))
    dockerCall(job, tool=DOCKER_SAMTOOLS, workDir=work_dir, parameters=docker_params)
    _fixPermissions(DOCKER_SAMTOOLS, work_dir)  # todo tmpfix

    # sanity check
    workdir_bai_location = os.path.join(work_dir, bam_filename + ".bai")
    if not os.path.isfile(workdir_bai_location):
        raise UserError("BAM index file not created for {}: {}".format(bam_filename, workdir_bai_location))

    # get start and end location
    get_idx_cmd = [
        ["samtools", "view", data_bam_location],
        ["head", "-n", "1"],
        [os.path.join("/data", write_select_column_script(work_dir))]
    ]
    start_idx_str = dockerCheckOutputExcept141(job, tool=DOCKER_SAMTOOLS, work_dir=work_dir, parameters=get_idx_cmd).strip()
    get_idx_cmd[1][0] = "tail"
    end_idx_str = dockerCheckOutputExcept141(job, tool=DOCKER_SAMTOOLS, work_dir=work_dir, parameters=get_idx_cmd).strip()
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
    job.fileStore.logToMaster("{}: Enqueueing {} jobs".format(config.uuid, len(positions)))
    idx = 0
    enqueued_jobs = 0
    return_values = list()
    total_mem = 0
    for position in positions:
        #prep
        chunk_position_description = "{}:{}-{}".format(config.contig_name, position[0], position[1])
        bam_split_command = ["view", "-b", data_bam_location, chunk_position_description]
        chunk_name = "{}.{}.bam".format(config.uuid, idx)
        #write chunk
        chunk_location = os.path.join(work_dir, chunk_name)
        with open(chunk_location, 'w') as out:
            job.fileStore.logToMaster("Running {} with parameters: {}".format(DOCKER_SAMTOOLS, bam_split_command))
            dockerCall(job, tool=DOCKER_SAMTOOLS, workDir=work_dir, parameters=bam_split_command, outfile=out)
        #document read count
        chunk_size = os.stat(chunk_location).st_size
        read_count= get_bam_read_count(job, work_dir, chunk_name)
        job.fileStore.logToMaster("{}: chunk from {} for idx {} is {}b ({}mb) and has {} reads"
                                  .format(config.uuid, chunk_position_description, idx, chunk_size,
                                          int(chunk_size / 1024 / 1024), read_count))
        # enqueue marginPhase job
        if read_count > 0:
            chunk_fileid = job.fileStore.writeGlobalFile(chunk_location)
            mp_cores = int(min(MP_CPU, config.maxCores))
            mp_mem = int(min(chunk_size * MP_MEM_BAM_FACTOR + ref_genome_size * MP_MEM_REF_FACTOR, config.maxMemory))
            mp_disk = int(min(chunk_size * MP_DSK_BAM_FACTOR + ref_genome_size * MP_DSK_REF_FACTOR, config.maxDisk))
            job.fileStore.logToMaster("{}:{} requesting {} cores, {}b ({}mb} disk, {}b ({}gb) mem"
                                      .format(config.uuid, idx, mp_cores,mp_disk, int(mp_disk / 1024 / 1024 ),
                                              mp_mem, int(mp_mem / 1024 / 1024 / 1024)))
            total_mem += mp_mem
            mp_mem = str(int(mp_mem / 1024)) + "K"
            mp_disk = str(int(mp_disk) / 1024) + "K"
            margin_phase_job = job.addChildJobFn(run_margin_phase, config, chunk_fileid, idx,
                                                 memory=mp_mem, cores=mp_cores, disk=mp_disk)
            return_values.append(margin_phase_job.rv())
            enqueued_jobs += 1
        idx += 1

    job.fileStore.logToMaster("{}: Enqueued {} jobs, requested total of {}gb ({}b) mem"
                              .format(config.uuid, int(total_mem/1024/1024/1024), total_mem))

    # enqueue consolidation job
    job.addFollowOnJobFn(consolidate_output, config, return_values)

    # log
    _log_time(job, "prepare_input", start, config.uuid)


def write_select_column_script(work_dir, column=4):
    # I feel bad for doing this, but I can't send single quotes into a toil command without them becoming escaped
    # so I'm just creating a script which does that
    filename = "select_column_{}.sh".format(column)
    file_location = os.path.join(work_dir, filename)
    with open(file_location, 'w') as out:
        print("#!/usr/bin/awk -f", file=out)
        print("{print $%d}" % column, file=out)
    os.chmod(file_location, 1023) #rwxrwxrwx
    return filename


def dockerCheckOutputExcept141(job, tool, work_dir, parameters):
    # there's something strange with the return code for commands which stop reading from stdin (like "head")
    # and so we need to ignore the returncode
    try:
        return dockerCheckOutput(job, tool, parameters=parameters, workDir=work_dir)
    except subprocess.CalledProcessError, e:
        if e.returncode == 141:
            return e.output
        else:
            raise e

def get_bam_read_count(job, work_dir, bam_name):
    params = [
        ["samtools", "view", os.path.join("/data", bam_name)],
        ["wc", "-l"]
    ]
    line_count_str = dockerCheckOutput(job, DOCKER_SAMTOOLS, params, work_dir)
    return int(line_count_str)


def run_margin_phase(job, config, chunk_file_id, chunk_idx):
    # prep
    start = time.time()
    work_dir = job.fileStore.getLocalTempDir()
    chunk_identifier = "{}.{}".format(config.uuid, chunk_idx)
    chunk_name = "{}.in.bam".format(chunk_identifier)
    chunk_location = os.path.join(work_dir, chunk_name)
    job.fileStore.logToMaster("run_margin_phase:{}:{}:{}".format(config.uuid, chunk_idx, datetime.datetime.now()))

    # download bam chunk
    job.fileStore.readGlobalFile(chunk_file_id, chunk_location)
    if not os.path.isfile(chunk_location):
        raise UserError("Failed to download chunk {} from {}".format(chunk_name, chunk_file_id))

    # download references
    #ref genome
    genome_reference_name = "reference.fa"
    genome_reference_location = os.path.join(work_dir, genome_reference_name)
    job.fileStore.readGlobalFile(config.reference_genome_fileid, genome_reference_location)
    if not os.path.isfile(genome_reference_location):
        raise UserError("Failed to download genome reference {} from {}"
                        .format(os.path.basename(config.reference_genome), config.reference_genome_fileid))
    #ref vcf
    vcf_reference_name = "reference.vcf"
    vcf_reference_location = os.path.join(work_dir, vcf_reference_name)
    job.fileStore.readGlobalFile(config.reference_vcf_fileid, vcf_reference_location)
    if not os.path.isfile(genome_reference_location):
        raise UserError("Failed to download vcf reference {} from {}"
                        .format(os.path.basename(config.reference_vcf), config.reference_vcf_fileid))
    # params
    params_name = "params.json"
    params_location = os.path.join(work_dir, params_name)
    job.fileStore.readGlobalFile(config.params_fileid, params_location)
    if not os.path.isfile(params_location):
        raise UserError("Failed to download params {} from {}"
                        .format(os.path.basename(config.params), config.params_fileid))

    # run marginPhase
    params = [os.path.join("/data", chunk_name), os.path.join("/data", genome_reference_name),
              "-p", os.path.join("/data", params_name), "-o", os.path.join("/data","{}.out".format(chunk_identifier)),
              "-r",  os.path.join("/data", vcf_reference_name)]
    job.fileStore.logToMaster("Running {} with parameters: {}".format(DOCKER_MARGIN_PHASE, params))
    dockerCall(job, tool=DOCKER_MARGIN_PHASE, workDir=work_dir, parameters=params)
    os.rename(os.path.join(work_dir, "marginPhase.log"), os.path.join(work_dir, "{}.log".format(chunk_identifier)))

    # document output
    job.fileStore.logToMaster("Output files for {}:".format(chunk_identifier))
    output_file_locations = glob.glob(os.path.join(work_dir, "{}*".format(chunk_identifier)))
    for f in output_file_locations:
        job.fileStore.logToMaster("\t{}".format(os.path.basename(f)))

    # tarball the output and save
    tarball_name = "{}.tar.gz".format(chunk_identifier)
    tarball_files(tar_name=tarball_name, file_paths=output_file_locations, output_dir=work_dir)
    output_fileid = job.fileStore.writeGlobalFile(os.path.join(work_dir, tarball_name))

    # log
    _log_time(job, "run_margin_phase", start, chunk_identifier)
    return output_fileid

def consolidate_output(job, config, all_output):
    #prep
    start = time.time()
    work_dir = job.fileStore.getLocalTempDir()
    out_tar = os.path.join(work_dir, '{}.tar.gz'.format(config.uuid))
    job.fileStore.logToMaster("consolidate_output:{}:{}".format(config.uuid, datetime.datetime.now()))
    job.fileStore.logToMaster("consolidating {} files".format(len(all_output)))

    # build tarball
    out_tars = [out_tar]
    with tarfile.open(out_tar, 'w:gz') as f_out:
        idx = 0
        for tar in all_output:
            tar_file = os.path.join(work_dir, "{}.tar.gz".format(idx))
            job.fileStore.readGlobalFile(tar, tar_file)
            out_tars.append(tar_file)
            with tarfile.open(tar_file, 'r') as f_in:
                for tarinfo in f_in:
                    with closing(f_in.extractfile(tarinfo)) as f_in_file:
                        f_out.addfile(tarinfo, fileobj=f_in_file)
            idx += 1

    # Move to output location
    if urlparse(config.output_dir).scheme == 's3':
        job.fileStore.logToMaster('Uploading {} to S3: {}'.format(out_tar, config.output_dir))
        s3am_upload(fpath=out_tar, s3_dir=config.output_dir, num_cores=config.cores)
    else:
        job.fileStore.logToMaster('Moving {} to output dir: {}'.format(out_tar, config.output_dir))
        mkdir_p(config.output_dir)
        copy_files(file_paths=[out_tar], output_dir=config.output_dir)

    # log
    _log_time(job, "consolidate_output", start, config.uuid)
    job.fileStore.logToMaster("END:{}:{}".format(config.uuid, datetime.datetime.now()))


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
        reference-contig: file:///your/reference/here.fa

        # Required: URL {scheme} to reference VCF
        reference-vcf: s3://your/reference/here.vcf

        # Required: URL {scheme} to marginPhase params JSON file
        params: ftp://your/params/here.vcf

        # Required: Name of contig
        # TODO: this should be done better
        contig-name: chr3

        # Required: Output location of sample. Can be full path to a directory or an s3:// URL
        # Warning: S3 buckets must exist prior to upload or it will fail.
        # Warning: Do not use "file://" syntax if output directory is local location
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
        #   URL         A URL {scheme} pointing to the sample or a full path to a directory
        #
        #   If a full path to a directory is provided for a sample, every file inside needs to be a fastq(.gz).
        #   Do not have other superfluous files / directories inside or the pipeline will complain.
        #
        #   Examples of several combinations are provided below. Lines beginning with # are ignored.
        #
        #   UUID_1  file:///path/to/file.bam    chr3
        #   UUID_2  s3://path/to/file.bam   chrX
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
        config.maxCores = int(args.maxCores) if args.maxCores else sys.maxsize
        config.maxDisk = int(args.maxDisk) if args.maxDisk else sys.maxint
        config.maxMemory = args.maxMemory if args.maxMemory else str(sys.maxint)

        # Config sanity checks
        require(config.output_dir, 'No output location specified')
        if urlparse(config.output_dir).scheme != "s3":
            mkdir_p(config.output_dir)
        if not config.output_dir.endswith('/'):
            config.output_dir += '/'
        require(config.reference_contig, 'No reference contig specified')
        require(config.contig_name, 'No contig name specified')
        require(config.reference_vcf, 'No reference vcf specified')
        require(config.params, 'No params specified')
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
