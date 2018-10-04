#!/usr/bin/env python2.7
from __future__ import print_function

import argparse
import os
import multiprocessing
import sys
import textwrap
import tarfile
from urlparse import urlparse
import shutil
import glob
import time
import datetime

import yaml
from bd2k.util.files import mkdir_p
from bd2k.util.processes import which
from toil import physicalMemory
from toil.job import Job
from toil_lib import require, UserError
from toil_lib.files import tarball_files, copy_files
from toil_lib.jobs import map_job
from toil_lib.urls import download_url, s3am_upload

from marginphase_chunk_merging import *
from marginphase_core import *



###########################################################
#                    WORKFLOW STEPS                       #
###########################################################


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
    log(job, "{}".format(datetime.datetime.now()), config.uuid, 'START')
    log(job, "Preparing input with URL:{}, contig:{}, reference_url:{}, params_url:{}"
                              .format(url, contig_name, reference_url, params_url), uuid, 'prepare_input')

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
    workdir_bam_location = os.path.join(work_dir, bam_filename)

    # index the bam
    _index_bam(job, config, work_dir, bam_filename)

    # sanity check
    workdir_bai_location = os.path.join(work_dir, bam_filename + ".bai")
    if not os.path.isfile(workdir_bai_location):
        raise UserError("BAM index file not created for {}: {}".format(bam_filename, workdir_bai_location))

    # get start and end location
    start_idx = sys.maxint
    end_idx = 0
    with closing(pysam.AlignmentFile(workdir_bam_location, 'rb' if bam_filename.endswith("bam") else 'r')) as aln:
        for read in aln.fetch():
            align_start = read.reference_start
            align_end = read.reference_end
            start_idx = min([start_idx, align_start])
            end_idx = max([end_idx, align_end])
    log(job, "start_pos:{}, end_pos:{}".format(config.uuid, start_idx, end_idx), uuid, 'prepare_input')

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
    log(job, "Enqueueing {} jobs".format(len(chunk_infos)), uuid, 'prepare_input')
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
            docker_call(job, config, work_dir, bam_split_command, DOCKER_SAMTOOLS_IMG, DOCKER_SAMTOOLS_TAG, outfile=out)

        #document read count
        chunk_size = os.stat(chunk_location).st_size
        ci[CI_CHUNK_SIZE] = chunk_size
        ci[CI_REF_FA_SIZE] = ref_genome_size
        read_count= prepare_input__get_bam_read_count(job, work_dir, chunk_name)
        ci[CI_READ_COUNT] = read_count
        log(job, "chunk from {} for idx {} is {}b ({}mb) and has {} reads".format(
            chunk_position_description, idx, chunk_size, int(chunk_size / 1024 / 1024), read_count),
            uuid, 'prepare_input')
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
            log(job, "requesting {} cores, {}b ({}mb) disk, {}b ({}gb) mem".format(
                mp_cores, mp_disk, int(mp_disk / 1024 / 1024 ), mp_mem, int(mp_mem / 1024 / 1024 / 1024)),
                "{}.{}".format(uuid, idx), 'prepare_input')
            mp_mem = str(int(mp_mem / 1024)) + "K"
            mp_disk = str(int(mp_disk) / 1024) + "K"
            margin_phase_job = job.addChildJobFn(run_margin_phase, config, chunk_fileid, ci,
                                                 memory=mp_mem, cores=mp_cores, disk=mp_disk)
            returned_tarballs.append(margin_phase_job.rv())
            enqueued_jobs += 1
        idx += 1

    log(job, "Enqueued {} jobs".format(enqueued_jobs), uuid, 'prepare_input')

    # enqueue merging and consolidation job
    merge_job = job.addFollowOnJobFn(merge_chunks, config, returned_tarballs)
    merge_job.addFollowOnJobFn(consolidate_output, config, merge_job.rv())

    # log
    log_time(job, "prepare_input", start, config.uuid)


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
    log(job, str(datetime.datetime.now()), chunk_identifier, 'run_margin_phase')

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
              os.path.join("/data", params_name), "-o", os.path.join("/data","{}.out".format(chunk_identifier)),
              '--tag', "{},{}-{}".format(chunk_idx, chunk_info[CI_CHUNK_BOUNDARY_START], chunk_info[CI_CHUNK_BOUNDARY_END])]
    if cpecan_prob_location is not None:
        params.extend(['--singleNuclProbDir', os.path.join("/data", cpecan_prob_location)])
    docker_call(job, config, work_dir, params, config.margin_phase_image, config.margin_phase_tag)
    os.rename(os.path.join(work_dir, "marginPhase.log"), os.path.join(work_dir, "{}.log".format(chunk_identifier)))

    # document output
    log(job, "Output files after marginPhase:", chunk_identifier, 'run_margin_phase')
    output_file_locations = glob.glob(os.path.join(work_dir, "{}*".format(chunk_identifier)))
    found_vcf, found_sam = False, False
    for f in output_file_locations:
        log(job, "\t\t{}".format(os.path.basename(f)), chunk_identifier, 'run_margin_phase')
        if f.endswith(VCF_SUFFIX): found_vcf = True
        if f.endswith(SAM_UNIFIED_SUFFIX): found_sam = True
    if cpecan_prob_location is not None:
        cpecan_tarball = glob.glob(os.path.join(work_dir, cpecan_prob_location, "*.tar.gz"))
        if len(cpecan_tarball) == 0:
            # todo why has tarball_files failed in this location?
            log(job, "Found no cpecan output tarball! Trying alt location.", chunk_identifier, 'run_margin_phase')
            cpecan_tarball = glob.glob(os.path.join(work_dir, "*.tar.gz"))

        if len(cpecan_tarball) == 0:
            log(job, "Found no cpecan output tarball!", chunk_identifier, 'run_margin_phase')
        elif len(cpecan_tarball) > 1:
            log(job, "Found {} cpecan output tarballs: {}!".format(len(cpecan_tarball), cpecan_tarball),
                chunk_identifier, 'run_margin_phase')
        else:
            log(job, "Saving cpecan output tarball", chunk_identifier, 'run_margin_phase')
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
                log(job, "", chunk_identifier, 'run_margin_phase')
                error = "Failed to generate appropriate output files {} times".format(MAX_RETRIES)
                log(job, error, chunk_identifier, 'run_margin_phase')
                # this enables us to "recover" in the face of failure during a run
                if CONTINUE_AFTER_FAILURE:
                    output_file_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, tarball_name))
                    chunk_info[CI_OUTPUT_FILE_ID] = output_file_id
                    return chunk_info
                raise UserError("{}:{}".format(chunk_identifier, error))

        log(job, "Missing output files.  Attepmting retry {}".format(config.retry_attempts),
            chunk_identifier, 'run_margin_phase')
        log(job, "Failed job log file:", chunk_identifier, 'run_margin_phase')
        log(job, "", chunk_identifier, 'run_margin_phase')
        with open(os.path.join(work_dir, "{}.log".format(chunk_identifier)), 'r') as input:
            for line in input:
                log(job, "\t\t{}".format(line.strip()), chunk_identifier, 'run_margin_phase')

        # new job
        retry_job = job.addChildJobFn(run_margin_phase, config, chunk_file_id, chunk_info,
                                      memory=str(int(config.maxMemory / 1024)) + "K", cores=job.cores, disk=job.disk)
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
    log_time(job, "run_margin_phase", start, chunk_identifier)
    return chunk_info


def run_margin_phase__run_cpecan_alignment(job, config, chunk_identifier, work_dir, alignment_filename, reference_filename):
    # prep
    start = time.time()
    fcn_identifier = "run_margin_phase:run_cpecan_alignment"
    log(job, "{}".format(datetime.datetime.now()), chunk_identifier, fcn_identifier)
    log(job, "Running cPecan positional probabilities on {}".format(alignment_filename),
        chunk_identifier, fcn_identifier)

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
               '--threads', str(job.cores) #is there a better way to read current allotted toil cores?
              ]]
    hmm_location = run_margin_phase__infer_cpecan_hmm_location(chunk_identifier)
    if hmm_location is not None: params[0].extend(['--realign_hmm', hmm_location])
    cpecan_output = docker_call(job, config, work_dir, params, config.cpecan_image, config.cpecan_tag)
    if DEBUG:
        log(job, "CPecan output:", chunk_identifier, fcn_identifier)
        for line in (cpecan_output if type(cpecan_output) == list else cpecan_output.split('\n')):
            log(job, "\t{}".format(line.strip()), chunk_identifier, fcn_identifier)

    # document output
    output_files = glob.glob(os.path.join(work_dir, out_dir_name, "*.tsv".format(chunk_identifier)))
    log(job, "cPecan generated {} output files".format(len(output_files)), chunk_identifier, fcn_identifier)

    # tarball the output and save
    tarball_name = "{}.nuc_pos_prob.tar.gz".format(chunk_identifier)
    try:
        tarball_files(tar_name=tarball_name, file_paths=output_files, output_dir=os.path.join(work_dir, out_dir_name))
    except Exception, e:
        log(job, "{} error making cPecan tarball: {}".format(type(e), e), chunk_identifier, fcn_identifier)
        tarball_files(tar_name=tarball_name, file_paths=output_files, output_dir=work_dir)
        log(job, "created tarball in work_dir: {}".format(os.path.join(work_dir)), chunk_identifier, fcn_identifier)

    # cleanup
    log_time(job, "run_cpecan_alignment", start, config.uuid)
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
    uuid = config.uuid
    work_dir = job.fileStore.getLocalTempDir()
    log(job, "{}".format(datetime.datetime.now()), uuid, 'merge_chunks')
    log(job, "Merging {} chunks".format(len(chunk_infos)), uuid, 'merge_chunks')
    if config.minimal_output:
        log(job, "Minimal output is configured, will only save full chromosome vcf and merged BAMs",
            uuid, 'merge_chunks')

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
        log(job, "Found {} missing indices: {}".format(len(missing_indices), missing_indices), uuid, 'merge_chunks')

    # prep for iteration
    merge_decisions = dict()
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
        log(job, "merging {} and {} across boundary {}".format(prev_chunk[CI_CHUNK_INDEX], chunk_idx, chunk_boundary),
            uuid, 'merge_chunks')

        # error out if missing files
        if curr_chunk_sam_file is None or curr_chunk_vcf_file is None:
            error = "{}: Missing expected output file, sam:{}, vcf:{}, chunk_info:{}".format(
                chunk_idx, curr_chunk_sam_file, curr_chunk_vcf_file, chunk)
            log(job, error, uuid, 'merge_chunks')
            job.fileStore.logToMaster(error)
            if CONTINUE_AFTER_FAILURE:
                # prev chunk info is maintained, and will be written during next chunk
                continue
            raise UserError("{}:{}".format(uuid, error))

        # skip writing the first chunk
        if prev_chunk_sam_file is None:
            curr_written_reads = set()
            curr_vcf_split_pos = 0
            curr_vcf_phase_action = dict()

        # write the rest of the chunks
        else:
            # get chunk splitting
            prev_reads, curr_reads, curr_vcf_split_pos, curr_vcf_phase_action, decision_summary =\
                merge_chunks__determine_chunk_splitting(job, merging_step_identifier, prev_chunk_sam_file,
                                                        curr_chunk_sam_file, chunk_boundary)
            merge_decisions[decision_summary] =\
                merge_decisions[decision_summary] + 1 if decision_summary in merge_decisions else 1

            # write sam
            curr_written_reads = merge_chunks__append_sam_reads(job, merging_step_identifier, prev_chunk_sam_file,
                                                                full_merged_sam_file, prev_reads, prev_written_reads)
            if len(curr_reads) > 0:
                curr_written_right_reads = merge_chunks__append_sam_reads(job, merging_step_identifier, curr_chunk_sam_file,
                                                                full_merged_sam_file, curr_reads, curr_written_reads)
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

    # loggit
    log(job, "Finished merge with following matches:", uuid, 'merge_chunks')
    job.fileStore.logToMaster("{}:merge_chunks: ".format(config.uuid))
    for decision, count in merge_decisions.items():
        log(job, "\t\t{}: \t{}".format(decision, count), uuid, 'merge_chunks')

    # tarball the output and save
    log(job, "Output files for merge:".format(), uuid, 'merge_chunks')
    output_file_locations = glob.glob(os.path.join(merged_chunks_directory, "*.*"))
    output_file_locations.sort()
    tmp = output_file_locations
    output_file_locations = list()
    for f in tmp:
        if os.path.isdir(f):
            log(job, "\t\t{} (skipped, directory)".format(os.path.basename(f)), uuid, 'merge_chunks')
        else:
            log(job, "\t\t{}".format(os.path.basename(f)), uuid, 'merge_chunks')
            output_file_locations.append(f)
    tarball_name = "{}.merged.tar.gz".format(config.uuid)
    tarball_files(tar_name=tarball_name, file_paths=output_file_locations, output_dir=work_dir)
    output_file_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, tarball_name))
    # we need to return the input list of chunk infos for consolidation
    chunk_infos.append({CI_OUTPUT_FILE_ID: output_file_id, CI_CHUNK_INDEX: "merged"})

    log_time(job, "merge_chunks", start, config.uuid)
    return chunk_infos


def consolidate_output(job, config, chunk_infos):
    #prep
    start = time.time()
    uuid = config.uuid
    work_dir = job.fileStore.getLocalTempDir()
    out_tar = os.path.join(work_dir, '{}.tar.gz'.format(config.uuid))

    log(job, "{}".format(datetime.datetime.now()), uuid, 'consolidate_output')
    log(job, "consolidating {} files".format(len(chunk_infos)), uuid, 'consolidate_output')

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
                        log(job, "(Minimal Output) Skipping output file: {}".format(tarinfo.name),
                            uuid, 'consolidate_output')
                        continue
                    if config.minimal_cpecan_output and tarinfo.name.endswith("gz"):
                        log(job, "(Minimal cPecan Output) Skipping output file: {}".format(tarinfo.name),
                            uuid, 'consolidate_output')
                        continue
                    log(job, "file {}".format(tarinfo.name), uuid, 'consolidate_output')
                    with closing(f_in.extractfile(tarinfo)) as f_in_file:
                        f_out.addfile(tarinfo, fileobj=f_in_file)
                        output_file_count += 1
    log(job, "Consolidated {} files in {} tarballs".format(output_file_count, len(out_tars)), uuid, 'consolidate_output')

    # Move to output location
    if urlparse(config.output_dir).scheme == 's3':
        log(job, "Uploading {} to S3: {}".format(out_tar, config.output_dir), uuid, 'consolidate_output')
        s3am_upload(fpath=out_tar, s3_dir=config.output_dir, num_cores=config.maxCores)
    else:
        log(job, "Moving {} to output dir: {}".format(out_tar, config.output_dir), uuid, 'consolidate_output')
        mkdir_p(config.output_dir)
        copy_files(file_paths=[out_tar], output_dir=config.output_dir)

    # log
    log_time(job, "consolidate_output", start, config.uuid)
    log(job, "{}".format(datetime.datetime.now()), uuid, 'END')


def _index_bam(job, config, work_dir, bam_filename):
    docker_call(job, config, work_dir, ["index", os.path.join("/data", bam_filename)],
                DOCKER_SAMTOOLS_IMG, DOCKER_SAMTOOLS_TAG)


###########################################################
#             WORKFLOW MANAGEMENT FUNCTIONS               #
###########################################################

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
