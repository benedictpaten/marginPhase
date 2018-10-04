from __future__ import print_function

import os
import time

from toil.lib.docker import apiDockerCall, dockerCall
from toil_lib import UserError


###########################################################
#                   WORKFLOW DEFINITIONS                  #
###########################################################

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


###########################################################
#                 TOIL UTILITY FUNCTIONS                  #
###########################################################


def docker_call(job, config, work_dir, params, image, tag, outfile=None):
    tagged_image = "{}:{}".format(image, tag)
    if DOCKER_LOGGING:
        log(job, "Running '{}' with parameters: {}".format(tagged_image, params), config.uuid, 'docker')
    if outfile is None:
        output = apiDockerCall(job, tagged_image, working_dir=work_dir, parameters=params, user="root")
    else:
        output = dockerCall(job, tool=tagged_image, workDir=work_dir, parameters=params, outfile=outfile)
    return output


def require_docker_file_output(job, config, work_dir, output_filenames, function_id, log_filename=None,
                               max_directory_contents=None, max_log_lines=None):
    missing_filenames = list(filter(lambda x: not os.path.isfile(os.path.join(work_dir, x)), output_filenames))
    if len(missing_filenames) > 0:
        # document missing
        log(job, "Missing files after docker call: ", config.uuid, function_id)
        for missing in missing_filenames:
            log(job, "\t{}".format(missing), config.uuid, function_id)

        # document contents
        directory_contents = os.listdir(work_dir)
        if max_directory_contents is not None and len(directory_contents) > max_directory_contents:
            directory_contents = directory_contents[0:max_directory_contents]
            directory_contents.append("[{} items total]".format(len(directory_contents)))
        log(job, "Current files in work_dir: {}".format(work_dir), config.uuid, function_id)
        for missing in directory_contents:
            log(job, "\t{}".format(missing), config.uuid, function_id)

        # document log
        if log_filename is not None:
            log_location = os.path.join(work_dir, log_filename)
            if os.path.isfile(log_location):
                log(job, "Log file contents: {}".format(log_filename), config.uuid, function_id)
                log_lines = 0
                with open(log_location) as log_stream:
                    for ll in log_stream:
                        if max_log_lines is None or log_lines < max_log_lines:
                            log(job, "\t{}".format(ll.rstrip()), config.uuid, function_id)
                        log_lines += 1
                if max_log_lines is not None and log_lines <= max_log_lines:
                    log(job, "\t[{} lines total]".format(log_lines), config.uuid, function_id)
        else:
            log(job, "Log file {} was not found".format(log_filename), config.uuid, function_id)

        # die
        raise UserError("Missing files after running {} on {}: {}".format(function_id, config.uuid, missing_filenames))


def log_debug_from_docker(job, log_file_location, identifier, function, input_file_locations=None):
    # sanity check
    if not os.path.isfile(log_file_location):
        log(job, "Logfile missing!", identifier, function)
        return

    # file size logging
    if input_file_locations is not None:
        if type(input_file_locations) == str: input_file_locations = [input_file_locations]
        for input_file_location in input_file_locations:
            log(job, "DEBUG_INPUT_FILESIZE:{}:{}".format(
                os.stat(input_file_location).st_size, os.path.basename(input_file_location)),
                identifier, function)

    # any logging from logfile
    with open(log_file_location) as log_file:
        for line in log_file:
            line = line.strip()
            if line.startswith("DEBUG"):
                log(job, line, identifier, function)


def log(job, message, identifier=None, function=None):
    if function is not None:
        message = "{}:{}".format(function, message)
    if identifier is not None:
        message = "{}:{}".format(identifier, message)
    job.log(message)


def log_time(job, function_name, start_time, sample_identifier=''):
    log(job, "TIME:{}".format(int(time.time() - start_time)), sample_identifier, function_name)

