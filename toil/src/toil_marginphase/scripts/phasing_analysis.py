#!/usr/bin/env python
from __future__ import print_function
import argparse
import glob
import sys
import numpy as np
import pysam
import time
from contextlib import closing


TAG_HAPLOTYPE = "ht"
TAG_CHUNK_ID = "mp"

PB_HAP1 = "pb_hap1"
PB_HAP2 = "pb_hap2"
PB_LENGTH = "pb_len"
PB_BLOCK_ID = "pb_block_id"

CI_CHUNK_IDX = "ci_chunk_idx"
CI_CHUNK_START = "ci_chunk_start"
CI_CHUNK_END = "ci_chunk_end"
CI_PHASE_BLOCKS = "ci_phase_blocks"
CI_READS = "ci_reads"
CI_FIRST_PB = "ci_first_pb"
CI_LAST_PB = "ci_last_pb"

RD_ID = "read_id"
RD_ALN_START = "rd_aln_start"
RD_ALN_END = "rd_aln_end"
RD_CHUNK_PHASE_BLOCKS = "rd_chunk_pbs"

RPB_BLOCK_ID = "rpb_bid"
RPB_READ_START = "rpb_rs"
RPB_READ_LENGTH = "rpb_rl"
RPB_IS_HAP1 = "rpb_h1"

def parse_args(args = None):
    parser = argparse.ArgumentParser("Provides statistics on a BAM/SAM file")
    parser.add_argument('--input_glob', '-i', dest='input_glob', required=True, type=str,
                        help='Glob matching SAM or BAM file(s) for analysis')
    # parser.add_argument('--generic_stats', '-g', dest='generic_stats', action='store_true', default=False,
    #                     help='Print generic stats for all files')
    # parser.add_argument('--read_length', '-l', dest='read_length', action='store_true', default=False,
    #                     help='Print statistics on read length for all files')
    # parser.add_argument('--read_depth', '-d', dest='read_depth', action='store_true', default=False,
    #                     help='Print statistics on read depth for all files')
    # parser.add_argument('--verbose', '-v', dest='verbose', action='store_true', default=False,
    #                     help='Print histograms for length and depth')
    # parser.add_argument('--silent', '-V', dest='silent', action='store_true', default=False,
    #                     help='Print nothing')
    # parser.add_argument('--depth_spacing', '-s', dest='depth_spacing', action='store', default=1000, type=int,
    #                     help='How far to sample read data')
    # parser.add_argument('--depth_range', '-r', dest='depth_range', action='store', default=None,
    #                     help='Whether to only calculate depth within a range, ie: \'100000-200000\'')
    # parser.add_argument('--filter_secondary', '-f', dest='filter_secondary', action='store_true', default=False,
    #                     help='Filter secondary alignments out (only include primary)')
    # parser.add_argument('--filter_primary', '-F', dest='filter_primary', action='store_true', default=False,
    #                     help='Filter primary alignments out (only include secondary)')
    # parser.add_argument('--min_alignment_threshold', '-a', dest='min_alignment_threshold', action='store', default=None,
    #                     type=int, help='Minimum alignment threshold, below which reads are not included')

    return parser.parse_args() if args is None else parser.parse_args(args)

def log(msg, depth=0):
    print("\t" * depth + str(msg))


def parse_chunk_id(chunk_id):
    parts = chunk_id.split(",")
    if len(parts) != 2 or len(parts[1].split("-")) != 2:
        log("malformed chunk_id: {}".format(chunk_id))
        return None, None, None
    positions = parts[1].split("-")
    chunk_idx, chunk_start, chunk_end = int(parts[0]), int(positions[0]), int(positions[1])
    return chunk_idx, chunk_start, chunk_end


def parse_phase_info(phase_info):
    parts = phase_info.split(",")
    if len(parts) != 4 or not (parts[0].startswith('h') and parts[1].startswith('p') and parts[2].startswith('r')
                               and parts[3].startswith('l')):
        log("malformed phase_info: {}".format(phase_info))
        return None, None, None, None
    parts = map(lambda x: int(x[1:]), parts)
    haplotype, phase_block, read_start, read_length = parts[0], parts[1], parts[2], parts[3]
    return haplotype, phase_block, read_start, read_length


def save_read_info(all_reads, read):
    # get info
    read_id = read.query_name
    align_start = read.reference_start
    align_end = read.reference_end

    # get storage data
    if read_id in all_reads:
        read_data = all_reads[read_id]
    else:
        read_data = dict()
        all_reads[read_id] = read_data
        read_data[RD_ID] = read_id
        read_data[RD_ALN_START] = align_start
        read_data[RD_ALN_END] = align_end
        read_data[RD_CHUNK_PHASE_BLOCKS] = dict()

    # return data
    return read_data


def save_phase_block_info(phase_blocks, phase_info, read_id):
    # get info
    haplotype, phase_block, read_start, read_length = parse_phase_info(phase_info)
    if haplotype == 0:
        return None

    # get storage data
    if phase_block in phase_blocks:
        phase_block_data = phase_blocks[phase_block]
    else:
        phase_block_data = dict()
        phase_blocks[phase_block] = phase_block_data
        phase_block_data[PB_BLOCK_ID] = phase_block
        phase_block_data[PB_HAP1] = set()
        phase_block_data[PB_HAP2] = set()
        phase_blocks[phase_block] = phase_block_data

    # read pb info
    read_phase_block_info = dict()
    read_phase_block_info[RPB_BLOCK_ID] = phase_block
    read_phase_block_info[RPB_READ_START] = read_start
    read_phase_block_info[RPB_READ_LENGTH] = read_length
    read_phase_block_info[RPB_IS_HAP1] = None

    # save read
    if haplotype == 1:
        phase_block_data[PB_HAP1].add(read_id)
        read_phase_block_info[RPB_IS_HAP1] = True
    elif haplotype == 2:
        phase_block_data[PB_HAP2].add(read_id)
        read_phase_block_info[RPB_IS_HAP1] = False
    else:
        log("unknown haplotype in phase_info for read {}: {}".format(read_id, phase_info))

    # return read phase data
    return read_phase_block_info


def save_chunk_info(chunks_by_idx, chunk_info):

    # parse data
    chunk_idx, chunk_start, chunk_end = parse_chunk_id(chunk_info)

    # get chunk_info_data
    if chunk_idx in chunks_by_idx:
        chunk_data = chunks_by_idx[chunk_idx]
        if chunk_start != chunk_data[CI_CHUNK_START]:
            log("\t\tmismatched chunks: {} / {}".format(chunk_start, chunk_data[CI_CHUNK_START]))
        if chunk_end != chunk_data[CI_CHUNK_END]:
            log("\t\tmismatched chunks: {} / {}".format(chunk_end, chunk_data[CI_CHUNK_END]))
    else:
        chunk_data = dict()
        chunks_by_idx[chunk_idx] = chunk_data
        chunk_data[CI_CHUNK_IDX] = chunk_idx
        chunk_data[CI_CHUNK_START] = chunk_start
        chunk_data[CI_CHUNK_END] = chunk_end
        chunk_data[CI_PHASE_BLOCKS] = dict()
        chunk_data[CI_READS] = set()
        chunk_data[CI_FIRST_PB] = None
        chunk_data[CI_LAST_PB] = None
        log("analyzing chunk {}: {} - {}".format(chunk_idx, chunk_start, chunk_end), depth=2)

    # return data
    return chunk_data


def main(args = None):
    # get our arguments
    args = parse_args() if args is None else parse_args(args)

    # get filenames, sanity check
    in_alignments = glob.glob(args.input_glob)
    if len(in_alignments) == 0:
        log("No files matching {}".format(args.input_glob))
        return
    else:
        log("Analyzing {} files".format(len(in_alignments)))

    # for storing all information
    all_analysis = dict()

    # iterate over all alignments
    for aln_filename in in_alignments:
        # sanity check
        if not (aln_filename.endswith("sam") or aln_filename.endswith("bam")):
            print("Matched file {} has unexpected filetype".format(aln_filename))
            continue

        # log
        log("{}:".format(aln_filename))

        # data storage
        aln_analysis = dict()
        all_analysis[aln_filename] = aln_analysis
        chunks_by_idx = dict()
        reads_by_id = dict()
        read_count = 0

        # examine all reads
        with closing(pysam.AlignmentFile(aln_filename, 'rb' if aln_filename.endswith("bam") else 'r')) as aln:
            # log
            log_depth = 1
            log("reading", depth=1)
            start = time.time()
            for read in aln.fetch():
                # get read data
                read_id = read.query_name
                read_count += 1

                # find haplotype tag
                for tag in [TAG_HAPLOTYPE, TAG_CHUNK_ID]:
                    if not read.has_tag(tag):
                        log("read {} had no {} tag".format(read_id, tag), depth=2)
                        continue

                # save read data
                read_data = save_read_info(reads_by_id, read)

                # get chunk data
                chunk_data = save_chunk_info(chunks_by_idx, read.get_tag(TAG_CHUNK_ID))
                chunk_idx = chunk_data[CI_CHUNK_IDX]
                if chunk_idx not in read_data[RD_CHUNK_PHASE_BLOCKS]:
                    read_data[RD_CHUNK_PHASE_BLOCKS][chunk_idx] = list()

                # save haplotpye data
                haplotype_tags = read.get_tag(TAG_HAPLOTYPE).split(";")
                for pb in haplotype_tags:
                    rpb_info = save_phase_block_info(chunk_data[CI_PHASE_BLOCKS], pb, read_id)
                    if rpb_info is not None: read_data[RD_CHUNK_PHASE_BLOCKS][chunk_idx].append(rpb_info)

            log("read {} reads ({}s)".format(read_count, int(time.time() - start)), depth=2)

        ### chunk and blockanalysis
        log("chunk/block analysis", depth=1)

        #prep
        chunk_idxs = list(chunks_by_idx.keys())
        chunk_idxs.sort()
        prev_chunk = None

        # analysis per chunk
        for chunk_idx in chunk_idxs:
            # log
            log("chunk {}:".format(chunk_idx), depth=2)

            # chunk_data
            chunk = chunks_by_idx[chunk_idx]
            phase_blocks = chunk[CI_PHASE_BLOCKS]
            phase_block_ids = list(phase_blocks.keys())
            phase_block_ids.sort()
            chunk[CI_FIRST_PB] = phase_block_ids[0]
            chunk[CI_LAST_PB] = phase_block_ids[-1]

            # phase block analysis
            prev_pb_id = None
            def print_pb(pb):
                log("pb:%9d  length:%9d  hap1:%6d  hap2:%6d" %
                    (pb[PB_BLOCK_ID], pb[PB_LENGTH], len(pb[PB_HAP1]), len(pb[PB_HAP2])), depth=3)
            for phase_block_id in phase_block_ids:
                pb = phase_blocks[phase_block_id]
                if pb[PB_BLOCK_ID] > chunk[CI_CHUNK_END]: break
                # print last one and get length
                if prev_pb_id is not None:
                    prev_pb = phase_blocks[prev_pb_id]
                    prev_pb[PB_LENGTH] = pb[PB_BLOCK_ID] - prev_pb[PB_BLOCK_ID]
                    print_pb(prev_pb)
                prev_pb_id = phase_block_id

            if prev_pb_id is not None:
                prev_pb = phase_blocks[prev_pb_id]
                prev_pb[PB_LENGTH] = chunk[CI_CHUNK_END] - prev_pb[PB_BLOCK_ID]
                print_pb(prev_pb)

        ### read analysis
        log("border analysis", depth=1)
        double_reads = list(filter(lambda x: len(x[RD_CHUNK_PHASE_BLOCKS]) > 1, reads_by_id.values()))
        double_chunk_reads = dict()
        for read in double_reads:
            chunks = list(read[RD_CHUNK_PHASE_BLOCKS].keys())
            chunks.sort()
            chunk_join_id = "-".join(map(str, chunks))
            if chunk_join_id not in double_chunk_reads: double_chunk_reads[chunk_join_id] = list()
            double_chunk_reads[chunk_join_id].append(read)

        double_chunk_ids = list(double_chunk_reads.keys())
        double_chunk_ids.sort(key=lambda x: int(x.split("-")[0]))
        for double_chunk_id in double_chunk_ids:
            chunk_ids = double_chunk_id.split("-")
            chunk_ids = list(map(int, chunk_ids))
            all_phase_blocks = set()
            prev_phase_blocks = dict()
            next_phase_blocks = dict()
            for read in double_chunk_reads[double_chunk_id]:
                for rpb in read[RD_CHUNK_PHASE_BLOCKS][chunk_ids[0]]:
                    rpb_id = rpb[RPB_BLOCK_ID]
                    all_phase_blocks.add(rpb_id)
                    prev_phase_blocks[rpb_id] = 1 if rpb_id not in prev_phase_blocks else (prev_phase_blocks[rpb_id] + 1)
                for rpb in read[RD_CHUNK_PHASE_BLOCKS][chunk_ids[1]]:
                    rpb_id = rpb[RPB_BLOCK_ID]
                    all_phase_blocks.add(rpb_id)
                    next_phase_blocks[rpb_id] = 1 if rpb_id not in next_phase_blocks else (next_phase_blocks[rpb_id] + 1)

            # analyze
            all_phase_blocks = list(all_phase_blocks)
            all_phase_blocks.sort()
            prev_str = "%3d:\t" % chunk_ids[0]
            next_str = "%3d:\t" % chunk_ids[1]
            for pb in all_phase_blocks:
                prev_str += ("%9d: %3d" % (pb, prev_phase_blocks[pb]) if pb in prev_phase_blocks else " " * 14) + "\t"
                next_str += ("%9d: %3d" % (pb, next_phase_blocks[pb]) if pb in next_phase_blocks else " " * 14) + "\t"

            log("chunk_boundary: {} ({})".format(double_chunk_id, chunks_by_idx[chunk_ids[0]][CI_CHUNK_END]), depth=2)
            log(prev_str, depth=3)
            log(next_str, depth=3)




    pass




if __name__ == "__main__":
    main()
