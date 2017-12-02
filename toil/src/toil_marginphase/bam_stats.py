#!/usr/bin/env python
from __future__ import print_function
import argparse
import glob
import json
import sys
import numpy as np
import os
import pysam
import math

# read keys
R_ID = "id"
R_START_POS = "start_pos"
R_END_POS = "end_pos"
R_LENGTH = "length"
R_SECONDARY = "secondary_alignment"
R_MAPPING_QUALITY = "mapping_quality"

#length summary
L_READ_COUNT = "read_count"
L_FILTERED_READ_COUNT = "filtered_read_count"
L_LOG_LENGTH_BUCKETS = "log_length_buckets"
L_MIN = "min_length"
L_MAX = "max_length"
L_AVG = "avg_length"
L_STD = "std_lenght"
L_LOG_BASE = "log_base"
L_LOG_MAX = "log_max"
L_ALL_LENGTHS = "all_lengths"

# depth summary
D_MAX = "max_depth"
D_MIN = "min_depth"
D_AVG = "avg_depth"
D_STD = "std_depth"
D_ALL_DEPTHS = "all_depths"
D_ALL_DEPTH_BINS = "all_depth_bins"
D_SPACING = "depth_spacing"
D_START_IDX = "depth_start_idx"
D_RANGE = "depth_range"


def parse_args(args = None):
    parser = argparse.ArgumentParser("Provides statistics on a BAM/SAM file")
    parser.add_argument('--input_glob', '-i', dest='input_glob', default="*.bam", type=str,
                        help='Glob matching SAM or BAM file(s)')
    parser.add_argument('--read_length', '-l', dest='read_length', action='store_true', default=False,
                        help='Print statistics on read length for all files')
    parser.add_argument('--read_depth', '-d', dest='read_depth', action='store_true', default=False,
                        help='Print statistics on read depth for all files')
    parser.add_argument('--verbose', '-v', dest='verbose', action='store_true', default=False,
                        help='Print histograms for length and depth')
    parser.add_argument('--silent', '-V', dest='silent', action='store_true', default=False,
                        help='Print nothing')
    parser.add_argument('--depth_spacing', '-s', dest='depth_spacing', action='store', default=1000, type=int,
                        help='How far to sample read data')
    parser.add_argument('--depth_range', '-r', dest='depth_range', action='store', default=None,
                        help='Whether to only calculate depth within a range, ie: \'100000-200000\'')
    parser.add_argument('--filter_secondary', '-f', dest='filter_secondary', action='store_true', default=False,
                        help='Filter secondary alignments')
    parser.add_argument('--min_alignment_threshold', '-a', dest='min_alignment_threshold', action='store', default=None,
                        type=int, help='Minimum alignment threshold, below which reads are not included')

    return parser.parse_args() if args is None else parser.parse_args(args)


def get_read_summary(read):
    summary = {
        R_START_POS: read.reference_start,
        R_END_POS: read.reference_end,
        R_LENGTH: read.reference_length,
        R_ID: read.query_name,
        R_SECONDARY: read.is_secondary,
        R_MAPPING_QUALITY: read.mapping_quality
    }
    return summary



def get_read_length_summary(read_summaries, length_log_base=2, length_log_max=32):
    # what we look for
    log_length_bins = [0 for _ in range(length_log_max)]
    all_lengths = []

    for read_summary in read_summaries:
        log_length_bins[int(math.log(read_summary[R_LENGTH], length_log_base))] += 1
        all_lengths.append(read_summary[R_LENGTH])

    summary = {
        L_LOG_LENGTH_BUCKETS: log_length_bins,
        L_MAX: max(all_lengths),
        L_MIN: min(all_lengths),
        L_AVG: np.mean(all_lengths),
        L_STD: np.std(all_lengths),
        L_LOG_BASE: length_log_base,
        L_LOG_MAX: length_log_max,
        L_ALL_LENGTHS: all_lengths
    }

    return summary


def print_read_length_summary(summary, verbose=False):
    print("\tREAD LENGTHS")
    print("\t\tmin: {}".format(summary[L_MIN]))
    print("\t\tmax: {}".format(summary[L_MAX]))
    print("\t\tavg: {}".format(int(summary[L_AVG])))
    print("\t\tstd: {}".format(int(summary[L_STD])))

    if verbose:
        length_log_max = summary[L_LOG_MAX]
        log_length_bins = summary[L_LOG_LENGTH_BUCKETS]
        print("\t\tread length log_{}:".format(summary[L_LOG_BASE]))
        max_bucket = max(list(filter(lambda x: log_length_bins[x] != 0, [x for x in range(length_log_max)])))
        min_bucket = min(list(filter(lambda x: log_length_bins[x] != 0, [x for x in range(length_log_max)])))
        max_bucket_size = max(log_length_bins)
        for bucket in range(min_bucket, max_bucket + 1):
            id = "%3d:" % bucket
            count = log_length_bins[bucket]
            pound_count = int(32.0 * count / max_bucket_size)
            print("\t\t\t{} {} {}".format(id, "#"*pound_count, count))


def get_read_depth_summary(read_summaries, spacing=1000, included_range=None):
    S, E = 's', 'e'

    # get reads which start or end on spacing interval
    positions = []
    for summary in read_summaries:
        positions.append((S, int(summary[R_START_POS]/spacing)))
        positions.append((E, int(summary[R_END_POS]/spacing)))

    # sort them: we iterate from the high end by popping off
    positions.sort(key=lambda x: x[1])
    start_idx = positions[0][1]
    end_idx = positions[-1][1]

    # data we care about
    depths = [0 for _ in range(end_idx - start_idx + 1)]

    # iterate over all read starts and ends
    depth = 0
    idx = end_idx
    while idx >= start_idx:
        curr = positions.pop()
        while curr[1] == idx:
            if curr[0] == E: depth += 1
            if curr[0] == S: depth -= 1
            # iterate
            if len(positions) == 0: break
            else: curr = positions.pop()
        positions.append(curr)
        # save and iterate
        depths[idx - start_idx] = depth
        idx -= 1

    assert depth == 0
    assert len(positions) == 1

    # check range before outputting summary
    if included_range is not None:
        # get range
        included_range = list(map(int, included_range.split("-")))
        if len(included_range) != 2:
            raise Exception("Malformed depth range: '{}'".format("-".join(included_range)))
        range_start = int(included_range[0]/spacing)
        range_end = int(included_range[1]/spacing)
        # sanity check
        if range_start > end_idx or range_end < start_idx or range_start >= range_end:
            raise Exception("Range {} outside of bounds of chunks: {}".format("-".join(included_range),
                                                                              "-".join([start_idx, end_idx])))
        # get appropriate depths
        new_depths = list()
        for depth_idx in range(end_idx - start_idx + 1):
            if start_idx + depth_idx < range_start: continue
            if start_idx + depth_idx > range_end: break
            new_depths.append(depths[depth_idx])
        # update values
        depths = new_depths
        start_idx = max(start_idx, range_start)
        assert len(depths) > 0

    # get read depth log value
    log_depth_bins = [0 for _ in range(16)]
    for depth in depths:
        if depth == 0:
            log_depth_bins[0] += 1
        else:
            log_depth_bins[int(math.log(depth, 2))] += 1

    # get depth summary
    summary = {
        D_MAX: max(depths),
        D_MIN: min(depths),
        D_AVG: np.mean(depths),
        D_STD: np.std(depths),
        D_ALL_DEPTHS: depths,
        D_ALL_DEPTH_BINS: log_depth_bins,
        D_SPACING: spacing,
        D_START_IDX: start_idx,
        D_RANGE: included_range
    }

    return summary


def print_read_depth_summary(summary, verbose=False):
    print("\tREAD DEPTHS")
    print("\t\tmax: {}".format(summary[D_MAX]))
    print("\t\tmin: {}".format(summary[D_MIN]))
    print("\t\tavg: {}".format(summary[D_AVG]))
    print("\t\tstd: {}".format(summary[D_STD]))
    log_depth_bins = summary[D_ALL_DEPTH_BINS]
    total_depths = sum(log_depth_bins)
    log_depth_pairs = [(i, log_depth_bins[i]) for i in range(len(log_depth_bins))]
    log_depth_pairs.sort(key=lambda x: x[1], reverse=True)
    print("\t\tmost frequent read depths [floor(log2(depth))]:")
    for i in range(0,min(len(list(filter(lambda x: x[1] != 0, log_depth_pairs))), 3)):
        print("\t\t\t#{}: depth:{} count:{} ({}%)".format(i + 1, log_depth_pairs[i][0], log_depth_pairs[i][1],
                                                          int(100.0 * log_depth_pairs[i][1] / total_depths)))

    if verbose:
        print("\t\tdepths with spacing {}{}:".format(summary[D_SPACING],
                                           "" if summary[D_RANGE] is None else ", and range {}".format(summary[D_RANGE])))
        idx = summary[D_START_IDX]
        for depth in summary[D_ALL_DEPTHS]:
            id = "%4d:" % idx
            pound_count = int(32.0 * depth / summary[D_MAX])
            print("\t\t\t{} {} {}".format(id, '#' * pound_count, depth))
            idx += 1

        print("\t\tread depth log_2 at above intervals:")
        max_bucket = max(list(filter(lambda x: log_depth_bins[x] != 0, [x for x in range(16)])))
        min_bucket = min(list(filter(lambda x: log_depth_bins[x] != 0, [x for x in range(16)])))
        max_bucket_size = max(log_depth_bins)
        for bucket in range(min_bucket, max_bucket + 1):
            id = "%3d:" % bucket
            count = log_depth_bins[bucket]
            pound_count = int(32.0 * count / max_bucket_size)
            print("\t\t\t{} {} {}".format(id, "#" * pound_count, count))


def main(args = None):
    # get our arguments
    args = parse_args() if args is None else parse_args(args)

    # get filenames, sanity check
    in_alignments = glob.glob(args.input_glob)
    if len(in_alignments) == 0:
        print("No files matching {}".format(args.input_glob))
        return 1
    else:
        if not args.silent: print("Analyzing {} files".format(len(in_alignments)))

    # data we care about
    length_summaries = dict()
    depth_summaries = dict()

    # iterate over all alignments
    for alignment_filename in in_alignments:
        # sanity check
        if not (alignment_filename.endswith("sam") or alignment_filename.endswith("bam")):
            print("Matched file {} has unexpected filetype".format(alignment_filename))
            continue

        # data we're gathering
        read_summaries = list()

        # get read data we care about
        samfile = None
        read_count = -1
        try:
            if not args.silent: print("Read {}:".format(alignment_filename))
            samfile = pysam.AlignmentFile(alignment_filename, 'rb' if alignment_filename.endswith("bam") else 'r')
            for read in samfile.fetch():
                read_count += 1
                read_summaries.append(get_read_summary(read))
            if not args.silent: print("read_count: {}".format(read_count))
        finally:
            if samfile is not None: samfile.close()

        # filter if appropriate
        if args.filter_secondary:
            read_summaries = list(filter(lambda x: not x[R_SECONDARY], read_summaries))
        if args.min_alignment_threshold is not None:
            read_summaries = list(filter(lambda x: x[R_MAPPING_QUALITY] >= args.min_alignment_threshold, read_summaries))
        if args.filter_secondary or args.min_alignment_threshold is not None:
            if not args.silent: print("filtered read_count: {} ".format(len(read_summaries)))

        # summarize
        length_summaries[alignment_filename] = get_read_length_summary(read_summaries)
        length_summaries[alignment_filename][L_READ_COUNT] = read_count
        length_summaries[alignment_filename][L_FILTERED_READ_COUNT] = len(read_summaries)
        depth_summaries[alignment_filename] = get_read_depth_summary(read_summaries,
                                                                     spacing=args.depth_spacing,
                                                                     included_range=args.depth_range)
        # print
        if args.read_length:
            if not args.silent: print_read_length_summary(length_summaries[alignment_filename], verbose=args.verbose)
        if args.read_depth:
            if not args.silent: print_read_depth_summary(depth_summaries[alignment_filename], verbose=args.verbose)

    return length_summaries, depth_summaries




if __name__ == "__main__":
    main()
