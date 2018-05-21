#!/usr/bin/env python
from __future__ import print_function
import argparse
import glob
import sys
import numpy as np
import pysam
import time
from contextlib import closing
import bam_stats




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

AN_READS = "reads"
AN_PHASE_BLOCKS = "pb"
AN_UNIQ_PHASE_BLOCKS = "uniq_pb"
AN_DEPTH_ANALYSIS = "depth"
AN_CHUNK_BOUNDARIES = "chunk_boundaries"
AN_CHUNK_SIZE = "chunk_size"

DEPTH_SPACING = 100


percent=lambda s, b: int(100.0 * s / b)


def parse_args(args = None):
    parser = argparse.ArgumentParser("Provides statistics on a BAM/SAM file")
    parser.add_argument('--input_glob', '-i', dest='input_glob', required=True, type=str,
                        help='Glob matching SAM or BAM file(s) for analysis')

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


def read_alignment_file(alignment_location):
    # log
    log("{}:".format(alignment_location))

    # data storage
    phase_blocks = dict()
    reads = dict()
    chunks = set()
    read_count = 0
    failed_reads = 0
    duplicated_reads = 0
    max_read_end = 0
    chrom = None

    with closing(pysam.AlignmentFile(alignment_location, 'rb' if alignment_location.endswith("bam") else 'r')) as aln:

        start = time.time()
        for read in aln.fetch():
            # get read data
            read_id = read.query_name
            max_read_end = max(max_read_end, read.reference_end)
            if chrom is None: chrom = read.reference_name
            if chrom != read.reference_name:
                raise Exception("{} has multiple contigs: {}, {}".format(alignment_location, chrom, read.reference_name))
            read_count += 1

            # find haplotype tag
            for tag in [TAG_HAPLOTYPE, TAG_CHUNK_ID]:
                if not read.has_tag(tag):
                    log("read {} had no {} tag".format(read_id, tag), depth=2)
                    failed_reads += 1
                    continue

            # see if duplicate
            if read_id in reads:
                duplicated_reads += 1
                continue

            # save read data
            read_data = save_read_info(reads, read)

            # save haplotpye data
            haplotype_tags = read.get_tag(TAG_HAPLOTYPE).split(";")
            for pb_tag in haplotype_tags:
                rpb_info = save_phase_block_info(phase_blocks, pb_tag, read_id)
                if rpb_info is not None: read_data[RD_PHASE_BLOCKS].append(rpb_info)

            # save chunk data
            chunks.add(read.get_tag(TAG_CHUNK_ID))

        if duplicated_reads > 0:
            log("found {} ({}%) duplicated reads".format(duplicated_reads, percent(duplicated_reads, read_count)),
                depth=1)
        log("read {} reads ({}s)".format(read_count, int(time.time() - start)), depth=1)

    # finish phase block analysis
    phase_block_ids = list(phase_blocks.keys())
    phase_block_ids.sort()
    prev_pb = None
    for pb_id in phase_block_ids:
        curr_pb = phase_blocks[pb_id]
        if prev_pb is not None: prev_pb[PB_END_POS] = curr_pb[PB_START_POS]
        prev_pb = curr_pb
    prev_pb[PB_END_POS] = max_read_end

    # finish chunk data
    chunk_boundaries = set()
    for chunk in chunks:
        _, chunk_start, chunk_end = parse_chunk_id(chunk)
        chunk_boundaries.add(chunk_start)
        chunk_boundaries.add(chunk_end)
    chunk_boundaries = list(chunk_boundaries)
    chunk_boundaries.sort()

    # return chunk data
    return reads, phase_blocks, chunk_boundaries, chrom


def get_depth_analysis(aln_filename, chrom):
    # get depth from bam_stats
    args = ['-i', aln_filename, '--silent', '--depth_spacing', str(DEPTH_SPACING)]
    bam_summaries, length_summaries, depth_summaries = bam_stats.main(args)
    return depth_summaries[aln_filename][chrom][bam_stats.D_ALL_DEPTH_MAP]


def get_average_phase_depth(depth_info, phase_block):
    phase_block = map(int, phase_block.split("-"))
    depths = [depth_info[int(x/DEPTH_SPACING)] for x in range(phase_block[0], phase_block[1])]
    return np.mean(depths)


def main(args = None):
    # get our arguments
    args = parse_args() if args is None else parse_args(args)

    # get filenames, sanity check
    in_alignments = glob.glob(args.input_glob)
    if len(in_alignments) == 0:
        log("No files matching {}".format(args.input_glob))
        return
    else:
        log("Analyzing {} files\n".format(len(in_alignments)))

    # for storing all information
    all_analysis = dict()
    depth_analysis = None

    # iterate over all alignments
    for aln_filename in in_alignments:
        # sanity check
        if not (aln_filename.endswith("sam") or aln_filename.endswith("bam")):
            print("Matched file {} has unexpected filetype".format(aln_filename))
            continue

        # get data
        aln_analysis = dict()
        all_analysis[aln_filename] = aln_analysis

        reads, phase_blocks, chunk_boundaries, chrom = read_alignment_file(aln_filename)
        aln_analysis[AN_READS] = reads
        aln_analysis[AN_CHUNK_BOUNDARIES] = chunk_boundaries
        aln_analysis[AN_PHASE_BLOCKS] = phase_blocks
        uniq_phase_blocks = set(map(lambda x: "{}-{}".format(x[PB_START_POS], x[PB_END_POS]), phase_blocks.values()))
        aln_analysis[AN_UNIQ_PHASE_BLOCKS]= uniq_phase_blocks
        if depth_analysis is None: depth_analysis = get_depth_analysis(aln_filename, chrom)

    log("\nAnalysis:")
    # how many files agree on phase block boundaries
    all_uniq_phase_blocks = dict()
    for aln_data in all_analysis.values():
        for uniq_pb in aln_data[AN_UNIQ_PHASE_BLOCKS]:
            if uniq_pb not in all_uniq_phase_blocks: all_uniq_phase_blocks[uniq_pb] = 0
            all_uniq_phase_blocks[uniq_pb] += 1

    # stats grouped by "agreement levels"
    total_uniq_pb = len(all_uniq_phase_blocks)
    # overall agreement count
    pb_representation = {x:0 for x in range(1,len(all_analysis)+1)}
    # length based on agreement count
    pb_rep_length = {x:list() for x in range(1,len(all_analysis)+1)}
    # depth based on agreement count
    pb_rep_depth = {x:list() for x in range(1,len(all_analysis)+1)}
    # distance to chunk boundary based on agreement count
    pb_start_chunk_dist = {x:list() for x in range(1,len(all_analysis)+1)}
    pb_center_chunk_dist = {x:list() for x in range(1,len(all_analysis)+1)}
    pb_end_chunk_dist = {x:list() for x in range(1,len(all_analysis)+1)}
    def min_dist_ratio_to_boundary(boundaries, location):
        half_chunk_size = (boundaries[1] - boundaries[0]) / 2.0
        closest_below, closest_above = None, None
        for boundary in boundaries:
            if boundary <= location: closest_below = abs(boundary - location)
            if boundary >= location: closest_above = abs(boundary - location)
            if closest_above is not None: break
        return min(closest_below, closest_above) / half_chunk_size

    # calculate these stats
    for phase_block, count in all_uniq_phase_blocks.items():
        # how many pb's are shared among 'count' files
        pb_representation[count] += 1
        # length of this phase block
        pb_rep_length[count].append(int(phase_block.split("-")[1]) - int(phase_block.split("-")[0]))
        # average depth of this phase block
        pb_rep_depth[count].append(get_average_phase_depth(depth_analysis, phase_block))
        # calculating how close to chunk splits these are
        for analysis in all_analysis.values():
            if phase_block not in analysis[AN_UNIQ_PHASE_BLOCKS]: continue
            pb_start = int(phase_block.split("-")[0])
            pb_end = int(phase_block.split("-")[1])
            pb_center = (pb_start + pb_end) / 2.0
            pb_start_chunk_dist[count].append(min_dist_ratio_to_boundary(analysis[AN_CHUNK_BOUNDARIES], pb_start))
            pb_center_chunk_dist[count].append(min_dist_ratio_to_boundary(analysis[AN_CHUNK_BOUNDARIES], pb_center))
            pb_end_chunk_dist[count].append(min_dist_ratio_to_boundary(analysis[AN_CHUNK_BOUNDARIES], pb_end))

    # buckets of distance
    DIST_BUCKETS = 6
    pb_start_dist_buckets = {y:{x:0 for x in range(DIST_BUCKETS)} for y in range(1,len(all_analysis)+1)}
    pb_center_dist_buckets = {y:{x:0 for x in range(DIST_BUCKETS)} for y in range(1,len(all_analysis)+1)}
    pb_end_dist_buckets = {y:{x:0 for x in range(DIST_BUCKETS)} for y in range(1,len(all_analysis)+1)}
    for count, values in pb_start_chunk_dist.items():
        for value in values:
            pb_start_dist_buckets[count][int(DIST_BUCKETS * value)] += 1
    for count, values in pb_center_chunk_dist.items():
        for value in values:
            pb_center_dist_buckets[count][int(DIST_BUCKETS * value)] += 1
    for count, values in pb_end_chunk_dist.items():
        for value in values:
            pb_end_dist_buckets[count][int(DIST_BUCKETS * value)] += 1
    def print_hist(buckets, depth):
        keys = buckets.keys()
        keys.sort()
        max_val = max(buckets.values())
        for key in keys:
            hash_cnt = int(12.0 * buckets[key] / max_val)
            log("{}: {}{} ({})"
                .format("boundary" if key == 0 else ("  center" if key + 1 == DIST_BUCKETS else "        "),
                        "#" * hash_cnt, " " * (12 - hash_cnt), buckets[key]), depth=depth)

    # report information organized based on consensus
    log("Unique PBs: {}".format(total_uniq_pb), depth=1)
    for x in reversed(range(1,len(all_analysis)+1)):
        log("")
        log("Extant in {} file: {} ({}%)".format(x, pb_representation[x], percent(pb_representation[x], total_uniq_pb)),
            depth=1)
        log("Length AVG: %8d" % np.mean(pb_rep_length[x]), depth=2)
        log("Length STD: %8d" % np.std(pb_rep_length[x]), depth=2)
        log("AvgDepth AVG:  %.2f" % np.mean(pb_rep_depth[x]), depth=2)
        log("AvgDepth STD:  %.2f" % np.std(pb_rep_depth[x]), depth=2)
        log("Avg Dist Ratio - phase block start to closest chunk boundary:  %.3f (std %.3f)" %
            (np.mean(pb_start_chunk_dist[x]), np.std(pb_start_chunk_dist[x])), depth=2)
        log("Avg Dist Ratio - phase block center to closest chunk boundary: %.3f (std %.3f)" %
            (np.mean(pb_center_chunk_dist[x]), np.std(pb_center_chunk_dist[x])), depth=2)
        log("Avg Dist Ratio - phase block end to closest chunk boundary:    %.3f (std %.3f)" %
            (np.mean(pb_end_chunk_dist[x]), np.std(pb_end_chunk_dist[x])), depth=2)
        log("PB start to chunk boundary:", depth=2)
        print_hist(pb_start_dist_buckets[x], 3)
        log("PB center to chunk boundary:", depth=2)
        print_hist(pb_center_dist_buckets[x], 3)
        log("PB end to chunk boundary:", depth=2)
        print_hist(pb_end_dist_buckets[x], 3)

    # do analysis of unique phase blocks
    singleton_phase_blocks = {y:dict() for y in
                              filter(lambda x: all_uniq_phase_blocks[x] == 1, all_uniq_phase_blocks.keys())}
    SPB_SPANNED_PHASE_BLOCKS = "spanned_pbs"
    SPB_WHOLLY_CONTAINS = "wholly_contain"
    SPB_CONTAINED_WHOLLY_WITHIN = "wholly_within"

    for spb in singleton_phase_blocks.keys():
        spanned_pb_lists = list()
        spb_contained_wholly = list()
        spb_wholly_contains = list()
        spb_start = int(spb.split("-")[0])
        spb_end = int(spb.split("-")[1])
        for aln_analysis in all_analysis.values():
            if spb in aln_analysis[AN_UNIQ_PHASE_BLOCKS]: continue
            spanned_pbs = list()
            spb_wholly_contains_value = 0
            spb_contained_wholly_value = 0
            for upb in aln_analysis[AN_UNIQ_PHASE_BLOCKS]:
                upb_start = int(upb.split("-")[0])
                upb_end = int(upb.split("-")[1])
                if upb_end > spb_start and spb_end > upb_start:
                    spanned_pbs.append(upb)
                    if upb_end <= spb_end and upb_start >= spb_start:
                        spb_wholly_contains.append(upb_end - upb_start)
                    elif spb_end <= upb_end and spb_start >= upb_start:
                        spb_contained_wholly.append(upb_end - upb_start)
            spanned_pb_lists.append(spanned_pbs)
            spb_wholly_contains.append(spb_wholly_contains_value)
            spb_contained_wholly.append(spb_contained_wholly_value)
        singleton_phase_blocks[spb][SPB_SPANNED_PHASE_BLOCKS] = spanned_pb_lists
        singleton_phase_blocks[spb][SPB_WHOLLY_CONTAINS] = list(filter(lambda x: x != 0, spb_wholly_contains))
        singleton_phase_blocks[spb][SPB_CONTAINED_WHOLLY_WITHIN] = list(filter(lambda x: x != 0, spb_contained_wholly))

    for spb in singleton_phase_blocks.keys():
        pass

    print(str(singleton_phase_blocks))



if __name__ == "__main__":
    main()
