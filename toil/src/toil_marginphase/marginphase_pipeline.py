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
from docker.errors import ContainerError
import pysam

import yaml
from bd2k.util.files import mkdir_p
from bd2k.util.processes import which
from toil import physicalMemory
from toil.job import Job
from toil.lib.docker import apiDockerCall, dockerCall
from toil_lib import require, UserError
from toil_lib.files import tarball_files, copy_files
from toil_lib.jobs import map_job
from toil_lib.urls import download_url, s3am_upload

# input / output schemes
SCHEMES = ('http', 'file', 's3', 'ftp')

# filenames
DEFAULT_CONFIG_NAME = 'config-toil-marginphase.yaml'
DEFAULT_MANIFEST_NAME = 'manifest-toil-marginphase.tsv'

# docker images
DOCKER_SAMTOOLS_IMG = "quay.io/ucsc_cgl/samtools"
DOCKER_SAMTOOLS_TAG = "1.8--cba1ddbca3e1ab94813b58e40e56ab87a59f1997"
DOCKER_MARGIN_PHASE_IMG_DEFAULT = "tpesout/margin_phase"
DOCKER_MARGIN_PHASE_TAG_DEFAULT = "latest"
DOCKER_CPECAN_IMG_DEFAULT = "tpesout/cpecan"
DOCKER_CPECAN_TAG_DEFAULT = "latest"

# resource
MP_CPU = 16
MP_MEM_BAM_FACTOR = 128 #todo account for learning iterations
MP_MEM_REF_FACTOR = 2
MP_DSK_BAM_FACTOR = 8 #input bam chunk, output (in sam fmt), vcf etc
MP_DSK_CPECAN_FACTOR = 4
MP_DSK_REF_FACTOR = 2

# for debugging
DEBUG = True
DOCKER_LOGGING = True

#chunk_info_keys
CI_CHUNK_BOUNDARY_START = "chunk_bounary_start_pos" #position where chunk split occured
CI_CHUNK_BOUNDARY_END = "chunk_bounary_end_pos" #position where chunk split occured
CI_CHUNK_START = "chunk_start_pos" #chunk boundary modified by the margin
CI_CHUNK_END = "chunk_end_pos" #chunk boundary modified by the margin
CI_READ_COUNT = "read_count" #how many reads were in the chunk
CI_CHUNK_SIZE = "chunk_size" #chunk size in bytes
CI_REF_FA_SIZE = "ref_fa_size" #chunk size in bytes
CI_CHUNK_INDEX = "chunk_index" #index of the chunk
CI_OUTPUT_FILE_ID = "output_file_id"

# read classification
CR_BEFORE_START = "before_start"
CR_SPAN_START = "span_start"
CR_WHOLLY_WITHIN = "wholly_within"
CR_SPAN_END = "span_end"
CR_AFTER_END = "after_end"
CR_START_HAP1 = "start_hap1"
CR_START_HAP2 = "start_hap2"
CR_END_HAP1 = "end_hap1"
CR_END_HAP2 = "end_hap2"
CR_START_PHASE_BLOCK = "start_phase_block"
CR_END_PHASE_BLOCK = "end_phase_block"

# merge_chunks phase block classification
PB_HAP1_READS = "pb_hap1"
PB_HAP2_READS = "pb_hap2"
PB_LENGTH = "pb_len"
PB_BLOCK_ID = "pb_block_id"
PB_START_POS = PB_BLOCK_ID
PB_END_POS = "pb_end_pos"

# merge_chunks read classification
RD_ID = "read_id"
RD_ALN_START = "rd_aln_start"
RD_ALN_END = "rd_aln_end"
RD_PHASE_BLOCKS = "rd_chunk_pbs"
RD_HAPLOTYPE_TAG = "rd_haplotype_tag"

# merge_chunks read phase block classification
RPB_BLOCK_ID = "rpb_bid"
RPB_READ_START = "rpb_rs"
RPB_READ_LENGTH = "rpb_rl"
RPB_IS_HAP1 = "rpb_h1"

# merge actions
ACTION_INVERT = "invert"

# helper functions
percent=lambda s, b: int(100.0 * s / b) if b != 0 else 0

# output naming conventions
SAM_UNIFIED_SUFFIX = "out.sam"
VCF_SUFFIX = "out.vcf"

# tags
TAG_GENOTYPE = "GT"
TAG_PHASE_SET = "PS"
TAG_MARGIN_PHASE_IDENTIFIER = "MPI"
TAG_HAPLOTYPE = "ht"
TAG_CHUNK_ID = "mp"

# cpecan locations - todo this is kind of a hack to not have to specify in config/manifest
CPECAN_NANOPORE_HMM = "/opt/cPecan/hmm/nanopore.hmm"
CPECAN_PACBIO_HMM = "/opt/cPecan/hmm/pacbio.s1-gc5.hmm"

# for retries
MAX_RETRIES = 1
CONTINUE_AFTER_FAILURE = False

def parse_samples(config, path_to_manifest):
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
            if len(sample) < 2:
                raise UserError('Bad manifest format! Required at least 2 tab-separated columns, got: {}'.format(sample))
            if len(sample) > 6:
                raise UserError('Bad manifest format! Required at most 6 tab-separated columns, got: {}'.format(sample))

            # extract sample parts
            uuid = sample[0]
            url = sample[1]
            contig_name, reference_url, params_url = "", "", ""
            if len(sample) > 2: contig_name = sample[2]
            if len(sample) > 3: reference_url = sample[3]
            if len(sample) > 4: params_url = sample[4]

            # fill defaults
            if len(contig_name) == 0: contig_name = config.default_contig
            if len(reference_url) == 0: reference_url = config.default_reference
            if len(params_url) == 0: params_url = config.default_params

            sample = [uuid, url, contig_name, reference_url, params_url]
            samples.append(sample)
    return samples


def prepare_input(job, sample, config):

    # job prep
    config = argparse.Namespace(**vars(config))
    uuid, url, contig_name, reference_url, params_url = sample
    config.uuid = uuid
    config.contig_name = contig_name
    config.reference_url = reference_url
    config.params_url = params_url
    work_dir = job.fileStore.getLocalTempDir()
    start = time.time()
    job.fileStore.logToMaster("{}:START:{}".format(config.uuid, datetime.datetime.now()))
    job.fileStore.logToMaster("{}:prepare_input: Preparing input with URL:{}, contig:{}, reference_url:{}, params_url:{}"
                              .format(uuid, url, contig_name, reference_url, params_url))

    # todo global resource estimation
    config.maxCores = min(config.maxCores, multiprocessing.cpu_count())
    config.defaultCores = min(MP_CPU, config.maxCores)
    config.maxMemory = min(config.maxMemory, int(physicalMemory() * .95))
    #config.disk

    # download references
    #ref fasta
    download_url(reference_url, work_dir=work_dir)
    ref_genome_filename = os.path.basename(reference_url)
    ref_genome_fileid = job.fileStore.writeGlobalFile(os.path.join(work_dir, ref_genome_filename))
    config.reference_genome_fileid = ref_genome_fileid
    ref_genome_size = os.stat(os.path.join(work_dir, ref_genome_filename)).st_size
    #params
    download_url(params_url, work_dir=work_dir)
    params_filename = os.path.basename(params_url)
    params_fileid = job.fileStore.writeGlobalFile(os.path.join(work_dir, params_filename))
    config.params_fileid = params_fileid

    # download bam
    download_url(url, work_dir=work_dir)
    bam_filename = os.path.basename(url)
    data_bam_location = os.path.join("/data", bam_filename)

    # index the bam
    _index_bam(job, config, work_dir, bam_filename)

    # sanity check
    workdir_bai_location = os.path.join(work_dir, bam_filename + ".bai")
    if not os.path.isfile(workdir_bai_location):
        raise UserError("BAM index file not created for {}: {}".format(bam_filename, workdir_bai_location))

    # get start and end location
    start_idx = sys.maxint
    end_idx = 0
    with closing(pysam.AlignmentFile(bam_filename, 'rb' if bam_filename.endswith("bam") else 'r')) as aln:
        for read in aln.fetch():
            align_start = read.reference_start
            align_end = read.reference_end
            start_idx = min([start_idx, align_start])
            end_idx = max([end_idx, align_end])
    job.fileStore.logToMaster("{}:prepare_input: start_pos:{}, end_pos:{}".format(config.uuid, start_idx, end_idx))

    # get reads from positions
    chunk_infos = list()
    idx = start_idx
    while idx < end_idx:
        ci = dict()
        ci[CI_CHUNK_BOUNDARY_START] = idx
        chunk_start = idx - config.partition_margin
        ci[CI_CHUNK_START] = chunk_start
        idx += config.partition_size
        ci[CI_CHUNK_BOUNDARY_END] = idx
        chunk_end = idx + config.partition_margin
        ci[CI_CHUNK_END] = chunk_end
        chunk_infos.append(ci)

    # enqueue jobs
    job.fileStore.logToMaster("{}:prepare_input: Enqueueing {} jobs".format(config.uuid, len(chunk_infos)))
    idx = 0
    enqueued_jobs = 0
    returned_tarballs = list()
    for ci in chunk_infos:
        #prep
        ci[CI_CHUNK_INDEX] = idx
        chunk_start = ci[CI_CHUNK_START]
        chunk_end = ci[CI_CHUNK_END]
        chunk_position_description = "{}:{}-{}".format(config.contig_name, chunk_start, chunk_end)
        bam_split_command = ["view", "-b", data_bam_location, chunk_position_description]
        chunk_name = "{}.{}.bam".format(config.uuid, idx)
        #write chunk
        chunk_location = os.path.join(work_dir, chunk_name)
        with open(chunk_location, 'w') as out:
            if DOCKER_LOGGING:
                job.fileStore.logToMaster("{}: Running {} with parameters: {}".format(config.uuid, "{}:{}".format(DOCKER_SAMTOOLS_IMG, DOCKER_SAMTOOLS_TAG), bam_split_command))
            dockerCall(job, tool="{}:{}".format(DOCKER_SAMTOOLS_IMG, DOCKER_SAMTOOLS_TAG), workDir=work_dir,
                       parameters=bam_split_command, outfile=out)
        #document read count
        chunk_size = os.stat(chunk_location).st_size
        ci[CI_CHUNK_SIZE] = chunk_size
        ci[CI_REF_FA_SIZE] = ref_genome_size
        read_count= prepare_input__get_bam_read_count(job, work_dir, chunk_name)
        ci[CI_READ_COUNT] = read_count
        job.fileStore.logToMaster("{}:prepare_input: chunk from {} for idx {} is {}b ({}mb) and has {} reads"
                                  .format(config.uuid, chunk_position_description, idx, chunk_size,
                                          int(chunk_size / 1024 / 1024), read_count))
        if config.intermediate_file_location is not None:
            copy_files(file_paths=[chunk_location], output_dir=config.intermediate_file_location)


        # enqueue marginPhase job
        if read_count > 0:
            chunk_fileid = job.fileStore.writeGlobalFile(chunk_location)
            mp_cores = config.defaultCores
            mp_mem = int(min(int(chunk_size * MP_MEM_BAM_FACTOR + ref_genome_size * MP_MEM_REF_FACTOR),
                             config.maxMemory))
            mp_disk = int(min(int(chunk_size * MP_DSK_BAM_FACTOR + ref_genome_size * MP_DSK_REF_FACTOR +
                                  (0 if config.cpecan_probabilities else MP_DSK_CPECAN_FACTOR) * chunk_size),
                              config.maxDisk))
            job.fileStore.logToMaster("{}:{}:prepare_input: requesting {} cores, {}b ({}mb) disk, {}b ({}gb) mem"
                                      .format(config.uuid, idx, mp_cores, mp_disk, int(mp_disk / 1024 / 1024 ),
                                              mp_mem, int(mp_mem / 1024 / 1024 / 1024)))
            mp_mem = str(int(mp_mem / 1024)) + "K"
            mp_disk = str(int(mp_disk) / 1024) + "K"
            margin_phase_job = job.addChildJobFn(run_margin_phase, config, chunk_fileid, ci,
                                                 memory=mp_mem, cores=mp_cores, disk=mp_disk)
            returned_tarballs.append(margin_phase_job.rv())
            enqueued_jobs += 1
        idx += 1

    job.fileStore.logToMaster("{}:prepare_input: Enqueued {} jobs".format(config.uuid, enqueued_jobs))

    # enqueue merging and consolidation job
    merge_job = job.addFollowOnJobFn(merge_chunks, config, returned_tarballs)
    merge_job.addFollowOnJobFn(consolidate_output, config, merge_job.rv())

    # log
    _log_time(job, "prepare_input", start, config.uuid)


def prepare_input__get_bam_read_count(job, work_dir, bam_name):
    params = [
        ["samtools", "view", os.path.join("/data", bam_name)],
        ["wc", "-l"]
    ]
    line_count_str = apiDockerCall(job, "{}:{}".format(DOCKER_SAMTOOLS_IMG, DOCKER_SAMTOOLS_TAG), parameters=params,
                                   working_dir=work_dir, detach=False, stdout=True)
    return int(line_count_str)


def run_margin_phase(job, config, chunk_file_id, chunk_info):
    # prep
    start = time.time()
    work_dir = job.fileStore.getLocalTempDir()
    chunk_idx = chunk_info[CI_CHUNK_INDEX]
    chunk_identifier = "{}.{}".format(config.uuid, chunk_idx)
    chunk_name = "{}.in.bam".format(chunk_identifier)
    chunk_location = os.path.join(work_dir, chunk_name)
    job.fileStore.logToMaster("{}:run_margin_phase:{}:{}".format(config.uuid, chunk_idx, datetime.datetime.now()))

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
    # params
    params_name = "params.json"
    params_location = os.path.join(work_dir, params_name)
    job.fileStore.readGlobalFile(config.params_fileid, params_location)
    if not os.path.isfile(params_location):
        raise UserError("Failed to download params {} from {}"
                        .format(os.path.basename(config.params), config.params_fileid))

    # do we want to run cPecan?
    cpecan_prob_location = None
    if config.cpecan_probabilities:
        cpecan_prob_location = run_margin_phase__run_cpecan_alignment(job, config, chunk_identifier, work_dir,
                                                                      chunk_name, genome_reference_name)

    # run marginPhase
    params = [os.path.join("/data", chunk_name), os.path.join("/data", genome_reference_name),
              "-p", os.path.join("/data", params_name), "-o", os.path.join("/data","{}.out".format(chunk_identifier)),
              '--tag', "{},{}-{}".format(chunk_idx, chunk_info[CI_CHUNK_BOUNDARY_START], chunk_info[CI_CHUNK_BOUNDARY_END])]
    if cpecan_prob_location is not None:
        params.extend(['--singleNuclProbDir', os.path.join("/data", cpecan_prob_location)])
    docker_margin_phase = "{}:{}".format(config.margin_phase_image, config.margin_phase_tag)
    if DOCKER_LOGGING:
        job.fileStore.logToMaster("{}:run_margin_phase: Running {} with parameters: {}".format(chunk_identifier, docker_margin_phase, params))
    apiDockerCall(job, docker_margin_phase, working_dir=work_dir, parameters=params, user="root")
    os.rename(os.path.join(work_dir, "marginPhase.log"), os.path.join(work_dir, "{}.log".format(chunk_identifier)))

    # document output
    job.fileStore.logToMaster("{}:run_margin_phase: Output files after marginPhase:".format(chunk_identifier))
    output_file_locations = glob.glob(os.path.join(work_dir, "{}*".format(chunk_identifier)))
    found_vcf, found_sam = False, False
    for f in output_file_locations:
        job.fileStore.logToMaster("{}:\t\t{}".format(chunk_identifier, os.path.basename(f)))
        if f.endswith(VCF_SUFFIX): found_vcf = True
        if f.endswith(SAM_UNIFIED_SUFFIX): found_sam = True
    if cpecan_prob_location is not None:
        cpecan_tarball = glob.glob(os.path.join(work_dir, cpecan_prob_location, "*.tar.gz"))
        if len(cpecan_tarball) == 0:
            # todo why has tarball_files failed in this location?
            job.fileStore.logToMaster("{}:run_margin_phase: Found no cpecan output tarball! Trying alt location"
                                      .format(chunk_identifier))
            cpecan_tarball = glob.glob(os.path.join(work_dir, "*.tar.gz"))

        if len(cpecan_tarball) == 0:
            job.fileStore.logToMaster("{}:run_margin_phase: Found no cpecan output tarball!".format(chunk_identifier))
        elif len(cpecan_tarball) > 1:
            job.fileStore.logToMaster("{}:run_margin_phase: Found {} cpecan output tarballs: {}!"
                                      .format(chunk_identifier, len(cpecan_tarball), cpecan_tarball))
        else:
            job.fileStore.logToMaster("{}:run_margin_phase: Saving cpecan output tarball".format(chunk_identifier))
            output_file_locations.append(cpecan_tarball[0])

    # tarball the output and save
    tarball_name = "{}.tar.gz".format(chunk_identifier)
    tarball_files(tar_name=tarball_name, file_paths=output_file_locations, output_dir=work_dir)

    # validate output, retry if not
    if not (found_sam and found_vcf):
        if "retry_attempts" not in config:
            config.retry_attempts = 1
        else:
            config.retry_attempts += 1
            if config.retry_attempts > MAX_RETRIES:
                error = "{}:run_margin_phase: Failed to generate appropriate output files {} times".format(chunk_identifier, MAX_RETRIES)
                job.fileStore.logToMaster(error)
                # this enables us to "recover" in the face of failure during a run
                if CONTINUE_AFTER_FAILURE:
                    output_file_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, tarball_name))
                    chunk_info[CI_OUTPUT_FILE_ID] = output_file_id
                    return chunk_info
                raise UserError(error)

        job.fileStore.logToMaster("{}:run_margin_phase: missing output files.  Attepmting retry {}"
                                  .format(chunk_identifier, config.retry_attempts))
        job.fileStore.logToMaster("{}:run_margin_phase: failed job log file:".format(chunk_identifier))
        with open(os.path.join(work_dir, "{}.log".format(chunk_identifier)), 'r') as input:
            for line in input:
                job.fileStore.logToMaster("{}:run_margin_phase:\t\t{}".format(chunk_identifier, line.strip()))
        # new job
        mp_cores = int(min(MP_CPU, config.maxCores))
        mp_mem = int(config.maxMemory) # on retry, try full mem
        mp_disk = int(min(int(chunk_info[CI_CHUNK_SIZE] * MP_DSK_BAM_FACTOR +
                              chunk_info[CI_REF_FA_SIZE] * MP_DSK_REF_FACTOR), config.maxDisk)) * 2 # on retry, double
        mp_mem = str(int(mp_mem / 1024)) + "K"
        mp_disk = str(int(mp_disk) / 1024) + "K"
        retry_job = job.addChildJobFn(run_margin_phase, config, chunk_file_id, chunk_info,
                                      memory=mp_mem, cores=mp_cores, disk=mp_disk)
        # save failed output
        if config.intermediate_file_location is not None:
            tarball_fail_name = "{}.FAILURE.{}.tar.gz".format(chunk_identifier, config.retry_attempts)
            os.rename(os.path.join(work_dir, tarball_name), os.path.join(work_dir, tarball_fail_name))
            copy_files(file_paths=[os.path.join(work_dir, tarball_fail_name)], output_dir=config.intermediate_file_location)

        return retry_job.rv()

    # if successfull, save output
    if config.intermediate_file_location is not None:
        copy_files(file_paths=[os.path.join(work_dir, tarball_name)], output_dir=config.intermediate_file_location)
    output_file_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, tarball_name))
    chunk_info[CI_OUTPUT_FILE_ID] = output_file_id

    # log
    _log_time(job, "run_margin_phase", start, chunk_identifier)
    return chunk_info


def run_margin_phase__run_cpecan_alignment(job, config, chunk_identifier, work_dir, alignment_filename, reference_filename):
    # prep
    start = time.time()
    job.fileStore.logToMaster("{}:run_margin_phase:run_cpecan_alignment:{}".format(chunk_identifier, datetime.datetime.now()))
    job.fileStore.logToMaster("{}:run_margin_phase:run_cpecan_alignment: Running cPecan positional probabilities on {}".format(chunk_identifier, alignment_filename))

    # index bam
    _index_bam(job, config, work_dir, alignment_filename)

    # build cPecan args
    out_dir_name = "cPecan_out"
    params = [['python', '/opt/cPecan/marginPhaseIntegration.py',
               '--ref', os.path.join("/data", reference_filename),
               '--alignment_file', os.path.join("/data", alignment_filename),
               '--workdir_directory', '/data/tmp',
               '--output_directory', os.path.join("/data", out_dir_name),
               '--realign_exe', '/opt/sonLib/bin/cPecanRealign',
               '--validate',
               '--threads', str(config.defaultCores) #is there a better way to read current allotted toil cores?
              ]]
    hmm_location = run_margin_phase__infer_cpecan_hmm_location(chunk_identifier)
    if hmm_location is not None: params[0].extend(['--realign_hmm', hmm_location])
    docker_cpecan = "{}:{}".format(config.cpecan_image, config.cpecan_tag)
    if DOCKER_LOGGING:
        job.fileStore.logToMaster("{}:run_margin_phase:run_cpecan_alignment: Running {} with parameters: {}".format(chunk_identifier, docker_cpecan, params))
    cpecan_output = apiDockerCall(job, docker_cpecan, working_dir=work_dir, parameters=params, user="root")
    if DEBUG:
        for line in (cpecan_output if type(cpecan_output) == list else cpecan_output.split('\n')):
            job.fileStore.logToMaster("{}:run_margin_phase:run_cpecan_alignment: \t{}".format(chunk_identifier, line))

    # document output
    output_files = glob.glob(os.path.join(work_dir, out_dir_name, "*.tsv".format(chunk_identifier)))
    job.fileStore.logToMaster("{}:run_margin_phase:run_cpecan_alignment: cPecan generated {} output files".format(chunk_identifier, len(output_files)))

    # tarball the output and save
    tarball_name = "{}.nuc_pos_prob.tar.gz".format(chunk_identifier)
    try:
        tarball_files(tar_name=tarball_name, file_paths=output_files, output_dir=os.path.join(work_dir, out_dir_name))
    except Exception, e:
        job.fileStore.logToMaster(
            "{}:run_margin_phase:run_cpecan_alignment: {} error making cPecan tarball: {}"
                .format(chunk_identifier, type(e), e))
        tarball_files(tar_name=tarball_name, file_paths=output_files, output_dir=work_dir)
        job.fileStore.logToMaster(
            "{}:run_margin_phase:run_cpecan_alignment: created tarball in work_dir: {}"
                .format(chunk_identifier, os.path.join(work_dir)))

    # cleanup
    _log_time(job, "run_cpecan_alignment", start, config.uuid)
    return out_dir_name


def run_margin_phase__infer_cpecan_hmm_location(identifier):
    if "np" in identifier and "pb" not in identifier:
        return CPECAN_NANOPORE_HMM
    if "np" not in identifier and "pb" in identifier:
        return CPECAN_PACBIO_HMM
    return None


def merge_chunks(job, config, chunk_infos):
    # prep
    start = time.time()
    work_dir = job.fileStore.getLocalTempDir()
    job.fileStore.logToMaster("{}:merge_chunks:{}".format(config.uuid, datetime.datetime.now()))
    job.fileStore.logToMaster("{}:merge_chunks: Merging {} chunks".format(config.uuid, len(chunk_infos)))
    if config.minimal_output:
        job.fileStore.logToMaster("{}:merge_chunks: Minimal output is configured, will only save full chromosome vcf and merged BAMs"
                                  .format(config.uuid))

    # work directory for tar management
    # output files
    merged_chunks_directory = os.path.join(work_dir, "merged")
    os.mkdir(merged_chunks_directory)
    full_merged_vcf_file = os.path.join(merged_chunks_directory, "{}.merged.vcf".format(config.uuid))
    full_merged_sam_file = os.path.join(merged_chunks_directory, "{}.merged.sam".format(config.uuid))

    # sort by chunk index and validate
    chunk_infos.sort(key=(lambda x: x[CI_CHUNK_INDEX]))
    idx = 0
    missing_indices = []
    for ci in chunk_infos:
        while ci[CI_CHUNK_INDEX] > idx:
            missing_indices.append(idx)
            idx += 1
        idx += 1
    if len(missing_indices) > 0:
        job.fileStore.logToMaster("{}:merge_chunks: Found {} missing indices: {}"
                                  .format(config.uuid, len(missing_indices), missing_indices))

    # prep for iteration
    prev_chunk_workdir = ""
    prev_chunk_sam_file = None
    prev_chunk_vcf_file = None
    prev_chunk = {CI_CHUNK_INDEX: "start"}
    prev_written_reads = set()
    prev_vcf_split_pos = None
    prev_vcf_phase_action = None

    # iterate over all chunks
    for chunk in chunk_infos:

        # get current chunk info/files
        chunk_idx = chunk[CI_CHUNK_INDEX]
        chunk_boundary = chunk[CI_CHUNK_BOUNDARY_START]
        merging_step_identifier = "{}:{}-{}".format(config.uuid, prev_chunk[CI_CHUNK_INDEX], chunk[CI_CHUNK_INDEX])
        curr_chunk_workdir = os.path.join(work_dir, "tmp-{}".format(chunk_idx))
        curr_chunk_sam_file, curr_chunk_vcf_file = merge_chunks__extract_chunk_tarball(
            job, config, curr_chunk_workdir, chunk)

        # error out if missing files
        if curr_chunk_sam_file is None or curr_chunk_vcf_file is None:
            error = "{}:merge_chunks:{}: Missing expected output file, sam:{}, vcf:{}, chunk_info:{}".format(
                config.uuid, chunk_idx, curr_chunk_sam_file, curr_chunk_vcf_file, chunk)
            job.fileStore.logToMaster(error)
            if CONTINUE_AFTER_FAILURE:
                # prev chunk info is maintained, and will be written during next chunk
                continue
            raise UserError(error)

        # skip writing the first chunk
        if prev_chunk_sam_file is None:
            curr_written_reads = set()
            curr_vcf_split_pos = 0
            curr_vcf_phase_action = dict()

        # write the rest of the chunks
        else:
            # get chunk splitting
            prev_reads, curr_reads, curr_vcf_split_pos, curr_vcf_phase_action =\
                merge_chunks__determine_chunk_splitting(job, merging_step_identifier, prev_chunk_sam_file,
                                                        curr_chunk_sam_file, chunk_boundary)

            # write sam
            curr_written_reads = merge_chunks__append_sam_reads(job, merging_step_identifier, prev_chunk_sam_file,
                                                                full_merged_sam_file, prev_reads, prev_written_reads)
            if len(curr_reads) > 0:
                curr_written_right_reads = merge_chunks__append_sam_reads(job, merging_step_identifier, curr_chunk_sam_file,
                                                                full_merged_sam_file, curr_reads, prev_written_reads)
                curr_written_reads = curr_written_reads.union(curr_written_right_reads)

            # write vcf
            merge_chunks__append_vcf_calls(job, merging_step_identifier, prev_chunk_vcf_file, full_merged_vcf_file,
                                           prev_vcf_split_pos, curr_vcf_split_pos, prev_vcf_phase_action,
                                           mp_identifier=prev_chunk[CI_CHUNK_INDEX])


        # cleanup
        if os.path.isdir(prev_chunk_workdir):
            shutil.rmtree(prev_chunk_workdir)

        # iterate
        prev_chunk = chunk
        prev_chunk_workdir = curr_chunk_workdir
        prev_chunk_sam_file = curr_chunk_sam_file
        prev_chunk_vcf_file = curr_chunk_vcf_file
        prev_written_reads = curr_written_reads
        prev_vcf_split_pos = curr_vcf_split_pos
        prev_vcf_phase_action = curr_vcf_phase_action

    # write the final reads and calls
    merging_step_identifier = "{}:{}-{}".format(config.uuid, prev_chunk[CI_CHUNK_INDEX], "end")
    merge_chunks__append_sam_reads(job, merging_step_identifier, prev_chunk_sam_file,
                                   full_merged_sam_file, {None: None}, prev_written_reads)
    merge_chunks__append_vcf_calls(job, merging_step_identifier, prev_chunk_vcf_file, full_merged_vcf_file,
                                   prev_vcf_split_pos, sys.maxint, prev_vcf_phase_action,
                                   mp_identifier=prev_chunk[CI_CHUNK_INDEX])

    # tarball the output and save
    job.fileStore.logToMaster("{}:merge_chunks: Output files for merge:".format(config.uuid))
    output_file_locations = glob.glob(os.path.join(merged_chunks_directory, "*.*"))
    output_file_locations.sort()
    tmp = output_file_locations
    output_file_locations = list()
    for f in tmp:
        if os.path.isdir(f):
            job.fileStore.logToMaster("{}:merge_chunks:\t\t{} (skipped, directory)".format(config.uuid, os.path.basename(f)))
        else:
            job.fileStore.logToMaster("{}:merge_chunks:\t\t{}".format(config.uuid, os.path.basename(f)))
            output_file_locations.append(f)
    tarball_name = "{}.merged.tar.gz".format(config.uuid)
    tarball_files(tar_name=tarball_name, file_paths=output_file_locations, output_dir=work_dir)
    output_file_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, tarball_name))
    # we need to return the input list of chunk infos for consolidation
    chunk_infos.append({CI_OUTPUT_FILE_ID: output_file_id, CI_CHUNK_INDEX: "merged"})

    _log_time(job, "merge_chunks", start, config.uuid)
    return chunk_infos
















def merge_chunks__determine_chunk_splitting(job, chunk_identifier, left_chunk_location, right_chunk_location, chunk_boundary):

    # get read and phase block info
    l_reads, l_phase_blocks = merge_chunks__read_chunk(job, chunk_identifier, left_chunk_location)
    r_reads, r_phase_blocks = merge_chunks__read_chunk(job, chunk_identifier, right_chunk_location)

    # organize chunk comparison
    all_phase_blocks, perfect_matches, inverted_matches, shared_phase_blocks = merge_chunks__organize_reads_and_blocks(
        job, chunk_identifier, l_reads, l_phase_blocks, r_reads, r_phase_blocks)

    # recommend strategy
    split_position, phase_block, invert_right, decision_summary = merge_chunks__recommend_merge_strategy(
        job, chunk_identifier, chunk_boundary, perfect_matches, inverted_matches, shared_phase_blocks)

    # implement strategy
    left_reads_writing, right_reads_writing, vcf_split_pos, vcf_right_phase_action = merge_chunks__specify_split_action(
        job, chunk_identifier, split_position, phase_block, invert_right,
        l_reads, l_phase_blocks, r_reads, r_phase_blocks)

    # log summarization
    job.fileStore.logToMaster(("{}:merge_chunks:determine_chunk_splitting: "
                              "read merge action: write {} from left, {} from right")
                              .format(chunk_identifier, len(left_reads_writing), len(right_reads_writing)))
    job.fileStore.logToMaster(("{}:merge_chunks:determine_chunk_splitting: call merge action: "
                              "split at {}, right action {}")
                              .format(chunk_identifier, vcf_split_pos, vcf_right_phase_action))

    # return
    return left_reads_writing, right_reads_writing, vcf_split_pos, vcf_right_phase_action


def merge_chunks__organize_reads_and_blocks(job, chunk_identifier, l_reads, l_phase_blocks, r_reads, r_phase_blocks):

    # get all phase blocks with same start pos
    all_phase_blocks = list(set(l_phase_blocks.keys()).union(set(r_phase_blocks.keys())))
    all_phase_block_count = len(all_phase_blocks)

    # get all phase blocks with same start pos
    shared_phase_blocks = list(set(l_phase_blocks.keys()).intersection(set(r_phase_blocks.keys())))
    shared_phase_blocks.sort()

    # get phase blocks by start and end pos
    l_pb_uniq_ids = list(map(lambda x: "{}-{}".format(x[PB_START_POS], x[PB_END_POS]), l_phase_blocks.values()))
    r_pb_uniq_ids = set(map(lambda x: "{}-{}".format(x[PB_START_POS], x[PB_END_POS]), r_phase_blocks.values()))

    # find all matches
    shared_phase_blocks = list()
    perfect_matches = list()
    inverted_matches = list()
    for l_pb_uniq in l_pb_uniq_ids:
        if l_pb_uniq in r_pb_uniq_ids:
            # we know phase block positions align
            shared_phase_blocks.append(l_pb_uniq)
            shared_phase_block = int(l_pb_uniq.split("-")[0])
            phase_block_median_pos = int((int(l_pb_uniq.split("-")[0]) + int(l_pb_uniq.split("-")[1]))/2.0)

            # get all reads in phase block (get from both blocks, in case some were skipped)
            reads_in_phase_block = set()
            for reads in [l_phase_blocks[shared_phase_block][PB_HAP1_READS],l_phase_blocks[shared_phase_block][PB_HAP2_READS],
                          r_phase_blocks[shared_phase_block][PB_HAP1_READS],r_phase_blocks[shared_phase_block][PB_HAP2_READS]]:
                for read in reads:
                    reads_in_phase_block.add(read)

            # get reads spannign median position
            reads_spanning_median = list(filter(
                lambda x: l_reads[x][RD_ALN_START] <= phase_block_median_pos and
                          l_reads[x][RD_ALN_END] >= phase_block_median_pos,
                reads_in_phase_block))

            # perfect match?
            perfect_matched_reads = list(filter(
                lambda x: l_reads[x][RD_HAPLOTYPE_TAG] == r_reads[x][RD_HAPLOTYPE_TAG], reads_spanning_median
            ))
            if len(perfect_matched_reads) == len(reads_spanning_median):
                perfect_matches.append(l_pb_uniq)
                continue

            # inverted match?
            inverted_matched_reads = list(filter(
                lambda x: l_reads[x][RD_HAPLOTYPE_TAG].replace("h1","<TMP>").replace("h2","h1").replace("<TMP>","h2")
                          == r_reads[x][RD_HAPLOTYPE_TAG], reads_spanning_median
            ))
            if len(inverted_matched_reads) == len(reads_spanning_median):
                inverted_matches.append(l_pb_uniq)
                continue

    # loggit
    job.fileStore.logToMaster("{}:merge_chunks:organize_reads_and_blocks: Found {} distinct phase blocks"
        .format(chunk_identifier, all_phase_block_count))
    job.fileStore.logToMaster("{}:merge_chunks:organize_reads_and_blocks: Found {} ({}%) perfect matches"
        .format(chunk_identifier, len(perfect_matches), percent(len(perfect_matches), all_phase_block_count)))
    job.fileStore.logToMaster("{}:merge_chunks:organize_reads_and_blocks: Found {} ({}%) inverted matches"
        .format(chunk_identifier, len(inverted_matches), percent(len(inverted_matches), all_phase_block_count)))
    job.fileStore.logToMaster("{}:merge_chunks:organize_reads_and_blocks: Found {} ({}%) matched phase starts"
        .format(chunk_identifier, len(shared_phase_blocks), percent(len(shared_phase_blocks), all_phase_block_count)))

    # return what we found
    return all_phase_blocks, perfect_matches, inverted_matches, shared_phase_blocks


def merge_chunks__recommend_merge_strategy(job, chunk_identifier, chunk_boundary, perfect_matches,
                                           inverted_matches, shared_phase_blocks):

    # helper function
    def pick_closest_elem_to_chunk_boundary(elems):
        elems.sort(key=lambda x: abs(chunk_boundary - sum(map(int, str(x).split("-"))) / len(str(x).split("-"))))
        return elems[0]

    # description of merge strategy
    split_position, phase_block, invert_right, decision_summary = None, None, None, None

    # case1: perfect match:
    #   READS: left of and spanning split_pos from left chunk, starting after split_pos from right chunk
    #   CALLS: left of split pos from left VCF, right of split pos from right VCF
    if len(perfect_matches) > 0:
        parts = map(int, pick_closest_elem_to_chunk_boundary(perfect_matches).split("-"))
        split_position = int(sum(parts) / len(parts))
        phase_block = parts[0]
        invert_right = False
        decision_summary = "PERFECT_MATCH"
        job.fileStore.logToMaster("{}:merge_chunks:recommend_merge_strategy: Found perfect match at "
                                  "pos {} in phase block {}".format(chunk_identifier, split_position, phase_block))


    # case2: perfect match but inverted haploptyes
    #   READS: left of and spanning split_pos from left chunk, starting after split_pos from right chunk
    #       reverse haplotype of included reads in phase_block from right chunk
    #   CALLS: left of split pos from left VCF, right of split pos from right VCF
    #       reverse phasing of calls in phase block from right chunk
    elif len(inverted_matches) > 0:
        parts = map(int, pick_closest_elem_to_chunk_boundary(inverted_matches).split("-"))
        split_position = int(sum(parts) / len(parts))
        phase_block = parts[0]
        invert_right = True
        decision_summary = "INVERT_MATCH"
        job.fileStore.logToMaster(
            "{}:merge_chunks:recommend_merge_strategy: Found inverted match at "
            "pos {} in phase block {}".format(chunk_identifier, split_position, phase_block))

    # case3: found a phase block starting at the same posistion in each chunk
    #   READS: finishing before split_pos from left chunk, starting after split pos from right chunk
    #       reads spanning split_pos get hap info from left before split_pos, and hap info from right after and including
    #   CALLS: left of split pos from left VCF, right of split pos from right VCF
    elif len(shared_phase_blocks) > 0:
        phase_block = pick_closest_elem_to_chunk_boundary(shared_phase_blocks)
        split_position = phase_block
        invert_right = False
        decision_summary = "PHASE_START_MATCH"
        job.fileStore.logToMaster("{}:merge_chunks:recommend_merge_strategy: Found phase block start match "
                                  "at {}".format(chunk_identifier, phase_block))

    # case4: no matching phase blocks
    #   READS: finishing before split_pos from left chunk, reads spanning split_pos get phasing finishing left of
    #       split_pos from left chunk, phasing in phase blocks spanning split_pos get split in two, phasing in phase
    #       blocks starting after split_pos from right chunk
    #   CALLS: starting left of split_pos from left VCF, starting after split_pos from right VCF, calls from right
    #       phase block spanning split_pos get new phase_block
    else:
        phase_block = None
        split_position = chunk_boundary
        invert_right = False
        decision_summary = "NO_MATCH"
        job.fileStore.logToMaster("{}:merge_chunks:recommend_merge_strategy: Found no match, creating "
                                  "new phase block at {}".format(chunk_identifier, split_position))

    # return data
    return split_position, phase_block, invert_right, decision_summary


def merge_chunks__specify_split_action(job, chunk_identifier, split_position, phase_block, invert_right,
                                       l_reads, l_phase_blocks, r_reads, r_phase_blocks):

    # describes read inclusion and modifications to haplotype string
    left_reads_writing = dict()
    right_reads_writing = dict()
    right_phase_block_conversion = dict()

    for read in l_reads.values():
        # this read belongs wholly to the right chunk
        if read[RD_ALN_START] > split_position:
            continue

        # this read belongs wholly to the left chunk
        elif read[RD_ALN_END] <= split_position:
            left_reads_writing[read[RD_ID]] = None

        # this read spans the split_position (and needs analysis)
        elif read[RD_ALN_START] <= split_position and read[RD_ALN_END] > split_position:

            # case4: new phase block created at split pos
            if phase_block is None:
                l_read = read
                r_read = r_reads[read[RD_ID]]
                new_hap_str, old_right_haplotype = merge_chunks__create_new_phase_block_at_position(split_position,
                                                                                                    l_read, r_read)
                left_reads_writing[read[RD_ID]] = [read[RD_HAPLOTYPE_TAG], new_hap_str]
                right_phase_block_conversion[old_right_haplotype] = split_position

                # santity check
                if len(right_phase_block_conversion) > 1:
                    raise UserError("SANITY_CHECK_FAILURE:{}: got inconsistent phase blocks ({}) spanning {} for read {}"
                                    .format(chunk_identifier, right_phase_block_conversion.keys(), split_position, read[RD_ID]))

            # case3: take hap info before split_pos from left, after right.  phase block exists at split_pos
            elif phase_block == split_position:
                l_read = l_reads[read[RD_ID]]
                r_read = r_reads[read[RD_ID]]
                haps = list(filter(lambda x: x[RPB_BLOCK_ID] < split_position, l_read[RD_PHASE_BLOCKS]))
                haps.extend(list(filter(lambda x: x[RPB_BLOCK_ID] >= split_position, r_read[RD_PHASE_BLOCKS])))
                haps.sort(key=lambda x: x[RPB_BLOCK_ID])
                new_hap_str = ";".join(map(merge_chunks__encode_phase_info, haps))
                left_reads_writing[read[RD_ID]] = [read[RD_HAPLOTYPE_TAG], new_hap_str]

            # case2, case1:
            else:
                left_reads_writing[read[RD_ID]] = None

    # get right reads we care about (reads in the spanning-the-split_pos phase chunk)
    # (everthing before split_pos comes from left chunk, everything after is unchanged)
    analysis_read_ids = set()
    if phase_block is None:
        if len(right_phase_block_conversion) == 0:
            job.fileStore.logToMaster("{}:merge_chunks:specify_split_action: No reads spanning {} were found!"
                .format(chunk_identifier, split_position))
        else:
            analysis_phase_block_id = list(right_phase_block_conversion.keys())[0]
            analysis_phase_block = r_phase_blocks[analysis_phase_block_id]
            analysis_read_ids = analysis_phase_block[PB_HAP1_READS].union(analysis_phase_block[PB_HAP2_READS])
    else:
        analysis_phase_block = r_phase_blocks[phase_block]
        analysis_read_ids = analysis_phase_block[PB_HAP1_READS].union(analysis_phase_block[PB_HAP2_READS])

    for read_id in analysis_read_ids:
        read = r_reads[read_id]

        # this read belongs wholly to the left chunk
        if read[RD_ALN_END] <= split_position:
            continue

        # this was analyzed with the left reads
        elif read_id in left_reads_writing:
            continue

        # now we need to analyize - we know these reads start after split_pos
        # case4
        if phase_block is None:
            if len(right_phase_block_conversion) == 0:
                raise UserError("SANITY_CHECK_FAILURE:{}: new phase block determined, but no conversion for read {}"
                                .format(chunk_identifier, read_id))
            pb_from = list(right_phase_block_conversion.keys())[0]
            pb_to = right_phase_block_conversion[pb_from]
            new_hap_str = read[RD_HAPLOTYPE_TAG].replace("p{},".format(pb_from), "p{},".format(pb_to))
            right_reads_writing[read_id] = [read[RD_HAPLOTYPE_TAG], new_hap_str]

        # case2
        elif invert_right:
            h1_str = "h1,p{}".format(phase_block)
            h1_tmp = "h1,p<TMP>"
            h2_str = "h2,p{}".format(phase_block)
            new_hap_str = read[RD_HAPLOTYPE_TAG].replace(h1_str, h1_tmp).replace(h2_str, h1_str).replace(h1_tmp, h2_str)
            right_reads_writing[read_id] = [read[RD_HAPLOTYPE_TAG], new_hap_str]

        # case1, case3
        else:
            pass

    # summarize vcf
    vcf_split_position = split_position
    vcf_right_phase_conversion = dict()

    if phase_block is None: # case 4
        if len(right_phase_block_conversion) != 0:
            # this also requires a new phase block when started
            vcf_right_phase_conversion = right_phase_block_conversion
        else: # no reads span this, so no action
            pass
    elif invert_right: # case 2
        vcf_right_phase_conversion = {phase_block: ACTION_INVERT}
    else:
        # case1 or case 3: no action
        pass

    # finish
    return left_reads_writing, right_reads_writing, vcf_split_position, vcf_right_phase_conversion


def merge_chunks__parse_phase_info(phase_info):
    parts = phase_info.split(",")
    if len(parts) != 4 or not (parts[0].startswith('h') and parts[1].startswith('p') and parts[2].startswith('r')
                               and parts[3].startswith('l')):
        return None, None, None, None
    parts = map(lambda x: int(x[1:]), parts)
    haplotype, phase_block, read_start, read_length = parts[0], parts[1], parts[2], parts[3]
    return haplotype, phase_block, read_start, read_length


def merge_chunks__encode_phase_info(read_data):
    haplotype = 0 if read_data[RPB_IS_HAP1] is None else (1 if read_data[RPB_IS_HAP1] else 2)
    phase_block = read_data[RPB_BLOCK_ID]
    read_start = read_data[RPB_READ_START]
    read_length = read_data[RPB_READ_LENGTH]
    return "h{},p{},r{},l{}".format(haplotype, phase_block, read_start, read_length)


def merge_chunks__save_read_info(all_reads, read):
    # get info
    read_id = read.query_name
    align_start = read.reference_start
    align_end = read.reference_end

    # save storage data
    read_data = dict()
    all_reads[read_id] = read_data
    read_data[RD_ID] = read_id
    read_data[RD_ALN_START] = align_start
    read_data[RD_ALN_END] = align_end
    read_data[RD_HAPLOTYPE_TAG] = read.get_tag(TAG_HAPLOTYPE)
    read_data[RD_PHASE_BLOCKS] = list()

    # return data
    return read_data


def merge_chunks__save_phase_block_info(job, chunk_identifier, phase_blocks, phase_info, read_id):
    # get info
    haplotype, phase_block, read_start, read_length = merge_chunks__parse_phase_info(phase_info)
    if haplotype is None:
        raise UserError("SANITY_CHECK_FAILURE:{}: malformed phase info for read {}: {}"
                        .format(chunk_identifier, read_id, phase_info))
    if haplotype == 0:
        return None

    # get storage data
    if phase_block in phase_blocks:
        phase_block_data = phase_blocks[phase_block]
    else:
        phase_block_data = dict()
        phase_blocks[phase_block] = phase_block_data
        phase_block_data[PB_BLOCK_ID] = phase_block
        phase_block_data[PB_HAP1_READS] = set()
        phase_block_data[PB_HAP2_READS] = set()
        phase_block_data[PB_END_POS] = None
        phase_blocks[phase_block] = phase_block_data

    # read pb info
    read_phase_block_info = dict()
    read_phase_block_info[RPB_BLOCK_ID] = phase_block
    read_phase_block_info[RPB_READ_START] = read_start
    read_phase_block_info[RPB_READ_LENGTH] = read_length
    read_phase_block_info[RPB_IS_HAP1] = None

    # save read
    if haplotype == 1:
        phase_block_data[PB_HAP1_READS].add(read_id)
        read_phase_block_info[RPB_IS_HAP1] = True
    elif haplotype == 2:
        phase_block_data[PB_HAP2_READS].add(read_id)
        read_phase_block_info[RPB_IS_HAP1] = False
    else:
        raise UserError("SANITY_CHECK_FAILURE:{}: unexpected haplotype in phase_info for read {}: {}"
                        .format(chunk_identifier, read_id, phase_info))

    # return read phase data
    return read_phase_block_info


def merge_chunks__read_chunk(job, chunk_identifier, chunk_location):
    # log
    job.fileStore.logToMaster("{}:merge_chunks:read_chunk: reading chunk {}"
                              .format(chunk_identifier, os.path.basename(chunk_location)))

    # data storage
    phase_blocks = dict()
    reads = dict()
    read_count = 0
    failed_reads = 0

    with closing(pysam.AlignmentFile(chunk_location, 'rb' if chunk_location.endswith("bam") else 'r')) as aln:

        start = time.time()
        for read in aln.fetch():
            # get read data
            read_id = read.query_name
            read_count += 1

            # find haplotype tag
            for tag in [TAG_HAPLOTYPE, TAG_CHUNK_ID]:
                if not read.has_tag(tag):
                    job.fileStore.logToMaster("{}:merge_chunks:read_chunk: read {} had no {} tag"
                                              .format(chunk_identifier, read_id, tag))
                    failed_reads += 1
                    continue

            # save read data
            read_data = merge_chunks__save_read_info(reads, read)

            # save haplotpye data
            haplotype_tags = read.get_tag(TAG_HAPLOTYPE).split(";")
            for pb_tag in haplotype_tags:
                rpb_info = merge_chunks__save_phase_block_info(job, chunk_identifier, phase_blocks, pb_tag, read_id)
                if rpb_info is not None: read_data[RD_PHASE_BLOCKS].append(rpb_info)

        job.fileStore.logToMaster("{}:merge_chunks:read_chunk: read {} reads ({}s)"
                                  .format(chunk_identifier, read_count, int(time.time() - start)))

    # finish phase block analysis
    phase_block_ids = list(phase_blocks.keys())
    phase_block_ids.sort()
    prev_pb = None
    for pb_id in phase_block_ids:
        curr_pb = phase_blocks[pb_id]
        if prev_pb is not None: prev_pb[PB_END_POS] = curr_pb[PB_START_POS]
        prev_pb = curr_pb
    # we aren't going to use this last one anyway
    prev_pb[PB_END_POS] = prev_pb[PB_START_POS]

    # return chunk data
    return reads, phase_blocks


def merge_chunks__create_new_phase_block_at_position(split_position, l_read, r_read):
    # get all documented haplotpyes
    l_haps = list(filter(lambda x: x[RPB_BLOCK_ID] < split_position, l_read[RD_PHASE_BLOCKS]))
    r_haps = list(filter(lambda x: x[RPB_BLOCK_ID] >= split_position, r_read[RD_PHASE_BLOCKS]))

    # data we want at the end
    haps = list()
    old_right_haplotype = None

    # get desired (and modified) haplotypes from l_read
    for hap in l_haps:
        if hap[RPB_BLOCK_ID] >= split_position:
            continue  # belongs to r_read
        elif hap[RPB_BLOCK_ID] + hap[RPB_READ_LENGTH] < split_position:
            haps.append(hap)  # before split
        else:  # this read needs to be split
            new_hap = {
                RPB_IS_HAP1: hap[RPB_IS_HAP1],
                RPB_BLOCK_ID: hap[RPB_BLOCK_ID],
                RPB_READ_START: hap[RPB_READ_START],
                RPB_READ_LENGTH: split_position - hap[RPB_BLOCK_ID]
            }
            haps.append(new_hap)

    # get desired (and modified) haplotypes from r_read
    for hap in r_haps:
        if hap[RPB_BLOCK_ID] >= split_position:
            haps.append(hap)  # after split
        elif hap[RPB_BLOCK_ID] + hap[RPB_READ_LENGTH] < split_position:
            continue  # belongs to l_read
        else:  # this read needs to be split
            split_diff = split_position - hap[RPB_BLOCK_ID]
            new_hap = {
                RPB_IS_HAP1: hap[RPB_IS_HAP1],
                RPB_BLOCK_ID: split_position,
                RPB_READ_START: hap[RPB_READ_START] + split_diff,
                RPB_READ_LENGTH: hap[RPB_READ_LENGTH] - split_diff
            }
            haps.append(new_hap)
            # sanity check and save old haplotype
            if old_right_haplotype is not None:
                raise UserError("SANITY_CHECK_FAILURE: " +
                                "found multiple phase_blocks ({}, {}) spanning split_position {} for read {}:".format(
                                    old_right_haplotype, hap[RPB_BLOCK_ID], split_position, l_read[RD_ID]))
            old_right_haplotype = hap[RPB_BLOCK_ID]

    # sanity check
    if old_right_haplotype is None:
        raise UserError("SANITY_CHECK_FAILURE: found no phase_blocks spanning split_position {} for read {}:"
            .format(split_position, l_read[RD_ID]))

    # save haploptyes
    haps.sort(key=lambda x: x[RPB_BLOCK_ID])
    new_hap_str = ";".join(map(merge_chunks__encode_phase_info, haps))
    return new_hap_str, old_right_haplotype


def merge_chunks__append_sam_reads(job, chunk_identifier, input_sam_file, output_sam_file, included_reads,
                                   excluded_reads):
    # data to track
    header_lines = 0
    written_read_ids = set()
    reads_not_written = 0
    explicitly_excluded = 0
    modified_writes = 0

    # include header?
    include_header = not os.path.isfile(output_sam_file)
    if include_header:
        job.fileStore.logToMaster("{}: output sam file does not exist, writing header".format(chunk_identifier))

    # read and write
    with open(output_sam_file, 'a') as output, open(input_sam_file, 'r') as input:
        for line in input:
            # header
            if line.startswith("@"):
                if include_header:
                    output.write(line)
                    header_lines += 1
                continue
            read_id = line.split("\t")[0]

            # already written in a previous chunk
            if read_id in excluded_reads:
                reads_not_written += 1
                explicitly_excluded += 1
                continue

            # should be written in this chunk
            elif read_id in included_reads:
                if included_reads[read_id] is not None:
                    for substitution in included_reads[read_id]:
                        line = line.replace(substitution[0], substitution[1])
                    modified_writes += 1
                output.write(line)
                written_read_ids.add(read_id)

            # edge case for inlusion of all reads in last chunk
            elif None in included_reads:
                output.write(line)
                written_read_ids.add(read_id)

            # if not explicitly told to write, read is not written
            else:
                reads_not_written += 1

    if include_header:
        job.fileStore.logToMaster("{}: wrote {} header lines from {} to {}".format(chunk_identifier, header_lines,
                                          os.path.basename(input_sam_file), os.path.basename(output_sam_file)))
    job.fileStore.logToMaster("{}: wrote {} reads ({} modified) with {} excluded ({} explicitly so) from {} to {}"
                              .format(chunk_identifier, len(written_read_ids), modified_writes,
                                      reads_not_written, explicitly_excluded,
                                      os.path.basename(input_sam_file), os.path.basename(output_sam_file)))
    return written_read_ids


def merge_chunks__append_vcf_calls(job, chunk_identifier, input_vcf_file, output_vcf_file, start_pos, end_pos,
                                   vcf_conversion, mp_identifier=None):

    # include header?
    include_header = not os.path.isfile(output_vcf_file)
    if include_header:
        job.fileStore.logToMaster("{}:merge_chunks:append_vcf_calls: output vcf file does not exist, writing header"
                                  .format(chunk_identifier))

    # data we track
    written_header_lines = 0
    written_lines = 0
    lines_outside_boundaries = 0

    # read and write
    with open(output_vcf_file, 'a') as output, open(input_vcf_file, 'r') as input:
        first_analyzed_line = True # may need to manage the phase set (only for the first phase of a chunk)
        for line in input:
            if line.startswith("#"):  # header
                if include_header:
                    output.write(line)
                    written_header_lines += 1
                continue

            # break line into parts
            line = line.rstrip().split("\t")
            position = int(line[1])
            if position < start_pos or position >= end_pos: #only include positions in given range
                lines_outside_boundaries += 1
                continue

            # get info and tags (and positions)
            info = line[-1].split(":")
            tags = line[-2].split(":")
            genotype_tag_idx = None
            phase_block_tag_idx = None
            idx = 0
            for tag in tags:
                if tag == TAG_GENOTYPE: genotype_tag_idx = idx
                if tag == TAG_PHASE_SET: phase_block_tag_idx = idx
                idx += 1
            if genotype_tag_idx is None:
                raise UserError("{}:SANITY_CHECK_FAILURE: malformed vcf {} phasing line (no {} tag): {}"
                                .format(chunk_identifier, os.path.basename(input_vcf_file), TAG_GENOTYPE, line))
            if phase_block_tag_idx is None:
                raise UserError("{}:SANITY_CHECK_FAILURE: malformed vcf {} phasing line (no {} tag): {}"
                                .format(chunk_identifier, os.path.basename(input_vcf_file), TAG_PHASE_SET, line))

            # phase block
            initial_phase_block = int(info[phase_block_tag_idx])
            reverse_phasing = False
            updated_phase_block = str(initial_phase_block)
            if initial_phase_block in vcf_conversion:
                if vcf_conversion[initial_phase_block] == ACTION_INVERT:
                    reverse_phasing = True
                else:
                    updated_phase_block = str(vcf_conversion[initial_phase_block])
            # set phase block (may be same as old phase block)
            info[phase_block_tag_idx] = updated_phase_block

            # get genotype
            genotype = info[genotype_tag_idx]
            continued_phase_block = "|" in genotype
            new_phase_block = "/" in genotype
            if not (continued_phase_block ^ new_phase_block):
                raise UserError("{}:SANITY_CHECK_FAILURE: Malformed vcf {} phasing line (unexpected genotype): {}"
                                .format(chunk_identifier, os.path.basename(input_vcf_file), line))
            genotype = genotype.split("|") if continued_phase_block else genotype.split("/")
            # make updates to genotype (potentially)
            if reverse_phasing:
                genotype.reverse()
            if first_analyzed_line and initial_phase_block != updated_phase_block:
                new_phase_block = True
            # save genotype
            genotype = ("/" if new_phase_block else "|").join(genotype)
            info[genotype_tag_idx] = genotype


            # add identifier (if appropriate)
            if mp_identifier is not None:
                info.append(mp_identifier)
                tags.append(TAG_MARGIN_PHASE_IDENTIFIER)
                line[-2] = ":".join(tags)

            # cleanup
            line[-1] = ":".join(info)
            line = "\t".join(line) + "\n"
            output.write(line)
            written_lines += 1
            first_analyzed_line = False

    # loggit
    job.fileStore.logToMaster(
        "{}:merge_chunks:append_vcf_calls: wrote {} lines ({} skipped from being outside boundaries {}-{}) from {} to {}"
            .format(chunk_identifier, written_lines, lines_outside_boundaries, start_pos, end_pos,
                    os.path.basename(input_vcf_file), os.path.basename(output_vcf_file)))


def merge_chunks__extract_chunk_tarball(job, config, tar_work_dir, chunk):
    # prep
    os.mkdir(tar_work_dir)
    tar_file = os.path.join(tar_work_dir, "chunk.tar.gz")

    # get file
    job.fileStore.readGlobalFile(chunk[CI_OUTPUT_FILE_ID], tar_file, mutable=True)
    with tarfile.open(tar_file, 'r') as tar:
        tar.extractall(tar_work_dir)

    # find desired files
    sam_unified, vcf = None, None
    for name in os.listdir(tar_work_dir):
        if name.endswith(VCF_SUFFIX): vcf = name
        elif name.endswith(SAM_UNIFIED_SUFFIX): sam_unified = name
    sam_unified_file = None if sam_unified is None else os.path.join(tar_work_dir, sam_unified)
    vcf_file = None if vcf is None else os.path.join(tar_work_dir, vcf)

    # return file locations
    return sam_unified_file, vcf_file


def consolidate_output(job, config, chunk_infos):
    #prep
    start = time.time()
    work_dir = job.fileStore.getLocalTempDir()
    out_tar = os.path.join(work_dir, '{}.tar.gz'.format(config.uuid))
    job.fileStore.logToMaster("{}:consolidate_output:{}".format(config.uuid, datetime.datetime.now()))
    job.fileStore.logToMaster("{}: consolidating {} files".format(config.uuid, len(chunk_infos)))

    # build tarball
    out_tars = [out_tar]
    output_file_count = 0
    with tarfile.open(out_tar, 'w:gz') as f_out:
        for ci in chunk_infos:
            file_id = ci[CI_OUTPUT_FILE_ID]
            tar_file = os.path.join(work_dir, "{}.tar.gz".format(ci[CI_CHUNK_INDEX]))
            job.fileStore.readGlobalFile(file_id, tar_file)
            out_tars.append(tar_file)
            with tarfile.open(tar_file, 'r') as f_in:
                for tarinfo in f_in:
                    if config.minimal_output and (
                            (tarinfo.name.endswith("bam") or
                                tarinfo.name.endswith("sam") or
                                tarinfo.name.endswith("bai"))
                            and "merged" not in tarinfo.name):
                        job.fileStore.logToMaster("{}: (Minimal Output) Skipping output file: {}".format(
                            config.uuid, tarinfo.name))
                        continue
                    if config.minimal_cpecan_output and tarinfo.name.endswith("gz"):
                        job.fileStore.logToMaster("{}: (Minimal cPecan Output) Skipping output file: {}".format(
                            config.uuid, tarinfo.name))
                        continue
                    job.fileStore.logToMaster("{}: file {}".format(config.uuid, tarinfo.name))
                    with closing(f_in.extractfile(tarinfo)) as f_in_file:
                        f_out.addfile(tarinfo, fileobj=f_in_file)
                        output_file_count += 1
    job.fileStore.logToMaster("{}: Consolidated {} files in {} tarballs".format(
        config.uuid, output_file_count, len(out_tars)))

    # Move to output location
    if urlparse(config.output_dir).scheme == 's3':
        job.fileStore.logToMaster('{}: Uploading {} to S3: {}'.format(config.uuid, out_tar, config.output_dir))
        s3am_upload(fpath=out_tar, s3_dir=config.output_dir, num_cores=config.maxCores)
    else:
        job.fileStore.logToMaster('{}: Moving {} to output dir: {}'.format(config.uuid, out_tar, config.output_dir))
        mkdir_p(config.output_dir)
        copy_files(file_paths=[out_tar], output_dir=config.output_dir)

    # log
    _log_time(job, "consolidate_output", start, config.uuid)
    job.fileStore.logToMaster("{}:END:{}".format(config.uuid, datetime.datetime.now()))


def _index_bam(job, config, work_dir, bam_filename):
    data_bam_location = os.path.join("/data", bam_filename)
    docker_params = ["index", data_bam_location]
    if DOCKER_LOGGING:
        job.fileStore.logToMaster("{}: Running {} with parameters: {}".format(config.uuid, "{}:{}".format(DOCKER_SAMTOOLS_IMG, DOCKER_SAMTOOLS_TAG), docker_params))
    apiDockerCall(job, "{}:{}".format(DOCKER_SAMTOOLS_IMG, DOCKER_SAMTOOLS_TAG), working_dir=work_dir,
                  parameters=docker_params, user="root")


def _log_time(job, function_name, start_time, sample_identifier=''):
    job.fileStore.logToMaster("{}:TIME:{}:{}".format(sample_identifier, function_name, int(time.time() - start_time)))


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

        # Required: Output location of sample. Can be full path to a directory or an s3:// URL
        # Warning: S3 buckets must exist prior to upload or it will fail.
        # Warning: Do not use "file://" syntax if output directory is local location
        output-dir: /tmp

        # Required: Size of each bam partition
        partition-size: 2000000

        # Required: Margin to apply on each partition
        partition-margin: 5000

        # Required: Minimum ratio of reads appearing in cross-chunk boundary to trigger a merge
        min-merge-ratio: .8

        # Optional: Identifier for marginPhase Docker image
        margin-phase-image: tpesout/margin_phase

        # Optional: Tag for marginPhase Docker image
        margin-phase-tag: latest

        # Optional: Perform cPecan single nucleotide probabilities calculation
        cpecan-probabilities: False

        # Optional: Identifier for cPecan Docker image
        cpecan-image: tpesout/cpecan

        # Optional: Tag for cpecan Docker image
        cpecan-tag: latest

        # Optional: URL {scheme} for default FASTA reference
        default-reference: file://path/to/reference.fa

        # Optional: Default contig name (must match sample URL's contig and reference fasta)
        default-contig: chr1

        # Optional: URL {scheme} for default parameters file
        default-params: file://path/to/reference.fa

        # Optional: Don't include BAM or SAM in output
        minimal-output: False

        # Optional: Don't include cPecan probabilities in output
        minimal-cpecan-output: False

        # Optional: for debugging, this will save intermediate files to the output directory (only works for file:// scheme)
        save-intermediate-files: False

    """.format(scheme=[x + '://' for x in SCHEMES])[1:])


def generate_manifest():
    return textwrap.dedent("""
        #   Edit this manifest to include information pertaining to each sample to be run.
        #
        #   There are 6 tab-separated columns: UUID, URL, contig name, reference fasta URL, parameters URL, reference vcf URL
        #
        #   UUID            Required    A unique identifier for the sample to be processed.
        #   URL             Required    A URL ['http://', 'file://', 's3://', 'ftp://'] pointing to the sample bam.  It must belong to a single contig
        #   CONTIG_NAME     Optional    Contig name (must match the contig in the URL and the reference)
        #   REFERENCE_URL   Optional    A URL ['http://', 'file://', 's3://', 'ftp://'] pointing to reference fasta file
        #   PARAMS_URL      Optional    A URL ['http://', 'file://', 's3://', 'ftp://'] pointing to parameters file for the run
        #
        #   For the four optional values, there must be a value specified either in this manifest, or in the configuration
        #   file.  Any value specified in the manifest overrides whatever is specified in the config file
        #
        #   Examples of several combinations are provided below. Lines beginning with # are ignored.
        #
        #   UUID_1\tfile:///path/to/file.bam
        #   UUID_2\ts3://path/to/file.bam\tchrX\ts3://path/to/chrX.reference.fa
        #   UUID_3\ts3://path/to/file.bam\tchr4\tfile:///path/to/chr4.reference.fa\tfile:///path/to/params.json
        #   UUID_4\tfile:///path/to/file.bam\t\t\tfile:///path/to/params.json
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

        # Parse config
        parsed_config = {x.replace('-', '_'): y for x, y in yaml.load(open(args.config).read()).iteritems()}
        config = argparse.Namespace(**parsed_config)
        config.maxCores = int(args.maxCores) if args.maxCores else sys.maxsize
        config.defaultCores = int(min(MP_CPU, config.maxCores))
        config.maxDisk = int(args.maxDisk) if args.maxDisk else sys.maxint
        config.maxMemory = sys.maxint
        # fix parsing of GB to int
        if args.maxMemory:
            args.maxMemory = args.maxMemory.upper()
            if args.maxMemory.endswith('B'):
                args.maxMemory = args.maxMemory.rstrip('B')
            # actual parsing
            if args.maxMemory.endswith('G'):
                config.maxMemory = int(args.maxMemory.rstrip('G')) * 1024 * 1024 * 1024
            elif args.maxMemory.endswith('M'):
                config.maxMemory = int(args.maxMemory.rstrip('M')) * 1024 * 1024
            elif args.maxMemory.endswith('K'):
                config.maxMemory = int(args.maxMemory.rstrip('K')) * 1024
            else:
                config.maxMemory = int(args.maxMemory)

        # Config sanity checks
        require(config.output_dir, 'No output location specified')
        if urlparse(config.output_dir).scheme != "s3":
            mkdir_p(config.output_dir)
        if not config.output_dir.endswith('/'):
            config.output_dir += '/'
        require(config.partition_size, "Configuration parameter partition-size is required")
        require(config.partition_margin, "Configuration parameter partition-margin is required")
        require(config.min_merge_ratio, "Configuration parameter min-merge-ratio is required")
        require(config.min_merge_ratio > .5 and config.min_merge_ratio <= 1,
                "Configuration parameter min-merge-ratio must be in range (.5,1]")
        if config.save_intermediate_files != True or urlparse(config.output_dir).scheme == "s3":
            config.intermediate_file_location = None
        else:
            intermediate_location = os.path.join(config.output_dir, "intermediate")
            mkdir_p(intermediate_location)
            config.intermediate_file_location = intermediate_location
        if "margin_phase_image" not in config or len(config.margin_phase_image) == 0:
            config.margin_phase_image = DOCKER_MARGIN_PHASE_IMG_DEFAULT
        if "margin_phase_tag" not in config or len(config.margin_phase_tag) == 0:
            config.margin_phase_tag = DOCKER_MARGIN_PHASE_TAG_DEFAULT
        if "cpecan_image" not in config or len(config.cpecan_image) == 0:
            config.cpecan_image = DOCKER_CPECAN_IMG_DEFAULT
        if "cpecan_tag" not in config or len(config.cpecan_tag) == 0:
            config.cpecan_tag = DOCKER_CPECAN_TAG_DEFAULT
        if "unittest" not in config:
            config.unittest = False
        if "minimal_output" not in config:
            config.minimal_output = False
        if "minimal_cpecan_output" not in config:
            config.minimal_cpecan_output = False
        if "cpecan_probabilities" not in config:
            config.cpecan_probabilities = False

        # get samples
        samples = parse_samples(config, args.manifest)

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
