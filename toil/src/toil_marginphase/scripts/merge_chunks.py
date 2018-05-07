#!/usr/bin/env python
from __future__ import print_function
import argparse
import glob
import sys
import numpy as np
import pysam
import time
import os
from contextlib import closing


TAG_HAPLOTYPE = "ht"
TAG_CHUNK_ID = "mp"

PB_HAP1_READS = "pb_hap1"
PB_HAP2_READS = "pb_hap2"
PB_LENGTH = "pb_len"
PB_BLOCK_ID = "pb_block_id"
PB_START_POS = PB_BLOCK_ID
PB_END_POS = "pb_end_pos"

RD_ID = "read_id"
RD_ALN_START = "rd_aln_start"
RD_ALN_END = "rd_aln_end"
RD_PHASE_BLOCKS = "rd_chunk_pbs"
RD_HAPLOTYPE_TAG = "rd_haplotype_tag"

RPB_BLOCK_ID = "rpb_bid"
RPB_READ_START = "rpb_rs"
RPB_READ_LENGTH = "rpb_rl"
RPB_IS_HAP1 = "rpb_h1"

percent=lambda s, b: int(100.0 * s / b)

def parse_args(args = None):
    parser = argparse.ArgumentParser("Merges phasing for two overlapping BAM/SAM files")
    parser.add_argument('--chunkLeft', '-l', dest='chunkLeft', required=True, type=str,
                        help='Location of left chunk BAM/SAM')
    parser.add_argument('--chunkRight', '-r', dest='chunkRight', required=True, type=str,
                        help='Location of right chunk BAM/SAM')

    return parser.parse_args() if args is None else parser.parse_args(args)


def log(msg, depth=0):
    print("\t" * depth + str(msg))


def merge_chunks__parse_chunk_id(chunk_id):
    parts = chunk_id.split(",")
    if len(parts) != 2 or len(parts[1].split("-")) != 2:
        log("malformed chunk_id: {}".format(chunk_id))
        return None, None, None
    positions = parts[1].split("-")
    chunk_idx, chunk_start, chunk_end = int(parts[0]), int(positions[0]), int(positions[1])
    return chunk_idx, chunk_start, chunk_end


def merge_chunks__parse_phase_info(phase_info):
    parts = phase_info.split(",")
    if len(parts) != 4 or not (parts[0].startswith('h') and parts[1].startswith('p') and parts[2].startswith('r')
                               and parts[3].startswith('l')):
        log("malformed phase_info: {}".format(phase_info))
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


def merge_chunks__save_phase_block_info(phase_blocks, phase_info, read_id):
    # get info
    haplotype, phase_block, read_start, read_length = merge_chunks__parse_phase_info(phase_info)
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
        log("unknown haplotype in phase_info for read {}: {}".format(read_id, phase_info))

    # return read phase data
    return read_phase_block_info


def merge_chunks__read_chunk(chunk_location):
    # log
    log("{}:".format(chunk_location))

    # data storage
    phase_blocks = dict()
    reads = dict()
    read_count = 0
    failed_reads = 0

    with closing(pysam.AlignmentFile(chunk_location, 'rb' if chunk_location.endswith("bam") else 'r')) as aln:

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
                    failed_reads += 1
                    continue

            # save read data
            read_data = merge_chunks__save_read_info(reads, read)

            # save haplotpye data
            haplotype_tags = read.get_tag(TAG_HAPLOTYPE).split(";")
            for pb_tag in haplotype_tags:
                rpb_info = merge_chunks__save_phase_block_info(phase_blocks, pb_tag, read_id)
                if rpb_info is not None: read_data[RD_PHASE_BLOCKS].append(rpb_info)

        log("read {} reads ({}s)".format(read_count, int(time.time() - start)), depth=2)

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
            if old_right_haplotype is not None: raise Exception("SANITY_CHECK_FAIL: " +
                                "found multiple phase_blocks ({}, {}) spanning split_position {} for read {}:".format(
                                    old_right_haplotype, hap[RPB_BLOCK_ID], split_position, l_read[RD_ID]))
            old_right_haplotype = hap[RPB_BLOCK_ID]

    # sanity check
    if old_right_haplotype is None: raise Exception(
        "SANITY_CHECK_FAIL: found no phase_blocks spanning split_position {} for read {}:".format(
            split_position, l_read[RD_ID]))

    # save haploptyes
    haps.sort(key=lambda x: x[RPB_BLOCK_ID])
    new_hap_str = ";".join(map(merge_chunks__encode_phase_info, haps))
    return new_hap_str, old_right_haplotype


def merge_chunks__organize_reads_and_blocks(l_reads, l_phase_blocks, r_reads, r_phase_blocks):

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
    log("Found {} distinct phase blocks".format(all_phase_block_count))
    log("Found {} ({}%) perfect matches".format(
        len(perfect_matches), percent(len(perfect_matches), all_phase_block_count)))
    log("Found {} ({}%) inverted matches".format(
        len(inverted_matches), percent(len(inverted_matches), all_phase_block_count)))
    log("Found {} ({}%) matched phase starts".format(
        len(shared_phase_blocks), percent(len(shared_phase_blocks), all_phase_block_count)))

    # return what we found
    return all_phase_blocks, perfect_matches, inverted_matches, shared_phase_blocks


def merge_chunks__recommend_merge_strategy(chunk_boundary, perfect_matches, inverted_matches, shared_phase_blocks):

    # helper function
    def pick_closest_elem(elems,center):
        elems.sort(key=lambda x: abs(center - sum(map(int, str(x).split("-"))) / len(str(x).split("-"))))
        return elems[0]

    split_position, phase_block, invert_right, decision_summary = None, None, None, None

    # case1: perfect match:
    #   READS: left of and spanning split_pos from left chunk, starting after split_pos from right chunk
    #   CALLS: left of split pos from left VCF, right of split pos from right VCF
    if len(perfect_matches) > 0:
        parts = map(int, pick_closest_elem(perfect_matches, chunk_boundary).split("-"))
        split_position = int(np.mean(parts))
        phase_block = parts[0]
        invert_right = False
        decision_summary = "PERFECT_MATCH"
        log("Found perfect match at pos {} in phase block {}".format(split_position, phase_block))


    # case2: perfect match but inverted haploptyes
    #   READS: left of and spanning split_pos from left chunk, starting after split_pos from right chunk
    #       reverse haplotype of included reads in phase_block from right chunk
    #   CALLS: left of split pos from left VCF, right of split pos from right VCF
    #       reverse phasing of calls in phase block from right chunk
    elif len(inverted_matches) > 0:
        parts = map(int, pick_closest_elem(inverted_matches, chunk_boundary).split("-"))
        split_position = int(np.mean(parts))
        phase_block = parts[0]
        invert_right = True
        decision_summary = "INVERT_MATCH"
        log("Found inverted match at pos {} in phase block {}".format(split_position, phase_block))

    # case3: found a phase block starting at the same posistion in each chunk
    #   READS: finishing before split_pos from left chunk, starting after split pos from right chunk
    #       reads spanning split_pos get hap info from left before split_pos, and hap info from right after and including
    #   CALLS: left of split pos from left VCF, right of split pos from right VCF
    elif len(shared_phase_blocks) > 0:
        phase_block = pick_closest_elem(shared_phase_blocks)
        split_position = phase_block
        invert_right = False
        decision_summary = "PHASE_START_MATCH"
        log("Found phase block start match at {}".format(phase_block))

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
        log("Found no match, creating new phase block at {}".format(split_position))

    # return data
    return split_position, phase_block, invert_right, decision_summary


def merge_chunks__specify_split_action(split_position, phase_block, invert_right,
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
                left_reads_writing[read[RD_ID]] = new_hap_str
                right_phase_block_conversion[old_right_haplotype] = split_position

                # santity check
                if len(right_phase_block_conversion) > 1:
                    raise Exception("SANITY_CHECK_FAIL: got inconsistent phase blocks ({}) spanning {} for read {}"
                                    .format(right_phase_block_conversion.keys(), split_position, read[RD_ID]))

            # case3: take hap info before split_pos from left, after right.  phase block exists at split_pos
            elif phase_block == split_position:
                l_read = l_reads[read[RD_ID]]
                r_read = r_reads[read[RD_ID]]
                haps = list(filter(lambda x: x[RPB_BLOCK_ID] < split_position, l_read[RD_PHASE_BLOCKS]))
                haps.extend(list(filter(lambda x: x[RPB_BLOCK_ID] >= split_position, r_read[RD_PHASE_BLOCKS])))
                haps.sort(key=lambda x: x[RPB_BLOCK_ID])
                new_hap_str = ";".join(map(merge_chunks__encode_phase_info, haps))
                left_reads_writing[read[RD_ID]] = new_hap_str

            # case2, case1:
            else:
                left_reads_writing[read[RD_ID]] = True

    # get right reads we care about (reads in the spanning-the-split_pos phase chunk)
    # (everthing before split_pos comes from left chunk, everything after is unchanged)
    analysis_read_ids = set()
    if phase_block is None:
        if len(right_phase_block_conversion) == 0:
            log("No reads spanning {} were found!".format(split_position))
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
                raise Exception("SANITY_CHECK_FAIL: new phase block determined, but no conversion for read {}"
                                .format(read_id))
            pb_from = list(right_phase_block_conversion.keys())[0]
            pb_to = right_phase_block_conversion[pb_from]
            new_hap_str = read[RD_HAPLOTYPE_TAG].replace("p{},".format(pb_from), "p{},".format(pb_to))
            right_reads_writing[read_id] = new_hap_str

        # case2
        elif invert_right:
            h1_str = "h1,p{}".format(phase_block)
            h1_tmp = "h1,p<TMP>"
            h2_str = "h2,p{}".format(phase_block)
            new_hap_str = read[RD_HAPLOTYPE_TAG].replace(h1_str, h1_tmp).replace(h2_str, h1_str).replace(h1_tmp, h2_str)
            right_reads_writing[read_id] = new_hap_str

        # case1, case3
        else:
            pass

    # summarize vcf
    vcf_split_position = split_position
    vcf_right_phase_action = None
    if invert_right:
        right_phase_action = "INVERT"
    elif phase_block is None:
        if len(right_phase_block_conversion) != 0:
            right_phase_action = right_phase_block_conversion
        else:
            # no reads span this, so no action
            pass
    else:
        # case1 or case 3: no action
        pass

    # finish
    return left_reads_writing, right_reads_writing, vcf_split_position, vcf_right_phase_action


def merge_chunks__determine_chunk_splitting(args = None):


    # get read and phase block info
    l_reads, l_phase_blocks = merge_chunks__read_chunk(args.chunkLeft)
    r_reads, r_phase_blocks = merge_chunks__read_chunk(args.chunkRight)

    # organize chunk comparison
    all_phase_blocks, perfect_matches, inverted_matches, shared_phase_blocks = merge_chunks__organize_reads_and_blocks(
        l_reads, l_phase_blocks, r_reads, r_phase_blocks)

    # todo
    chunk_boundary = np.median(all_phase_blocks)


    # recommend strategy
    split_position, phase_block, invert_right, decision_summary = merge_chunks__recommend_merge_strategy(
        chunk_boundary, perfect_matches, inverted_matches, shared_phase_blocks)


    # implement strategy
    left_reads_writing, right_reads_writing, vcf_split_pos, vcf_right_phase_action = merge_chunks__specify_split_action(
        split_position, phase_block, invert_right, l_reads, l_phase_blocks, r_reads, r_phase_blocks)

    # log summarization
    log("read merge action: write {} from left, {} from right".format(len(left_reads_writing), len(right_reads_writing)))
    log("call merge action: split at {}, right action {}".format(vcf_split_pos, vcf_right_phase_action))


    """
    CASES
    1) perfect split  
            split_position
            phase_block != split_position
            invert_right = False
        take reads left of and spanning split_position from prev chunk
        take reads starting after split_position from curr chunk
        take calls based on split_position
        no modification
        ITERATION:
            reads written from prev
    2) inverted split
            split_position
            phase_block != split_position
            invert_right = True
        take reads left of and spanning split_position from prev chunk
        take reads starting after split_position from curr chunk
        take calls based on split_position
        reverse haplotype string of phase_block from read in curr chunk
        reverse phasing on calls in vcf
        ITERATION:
            reads written from prev
            reads written with reversed phase block
    3) shared phase start
            split_position
            phase_block = split_position
            invert_right = False
        take reads finishing left of split_position from prev chunk
        take reads spanning start of split_position from prev and curr chunk, 
                keep haplotyping from prev chunk until split_pos
                modify haplotyping to include phasing after split pos
        take reads starting after phase from curr chunk
        take calls based on split_position
        ITERATION:
            reads written from prev
            reads written from curr
    4) new phase block
            split_position
            phase_block = None
            invert_right = False
        take reads finishing left of split_position
        take reads spanning split_position phase from prev and curr chunk, 
                keep haplotyping from prev chunk until split_pos
                modify haplotyping to include phasing after split pos
        ITERATION:
            reads written from prev
            reads written from curr
    """



if __name__ == "__main__":

    # get our arguments
    args = parse_args()
    assert False not in [map(os.path.isfile, [args.chunkLeft, args.chunkRight])]
    merge_chunks__determine_chunk_splitting(args.chunkLeft, args.chunkRight)
