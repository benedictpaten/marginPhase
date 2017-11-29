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

#length summary
L_LOG_LENGTH_BUCKETS = "log_length_buckets"
L_MIN = "min_length"
L_MAX = "max_length"
L_AVG = "avg_length"
L_STD = "std_lenght"

# depth summary
D_MAX = "max_depth"
D_AVG = "avg_depth"
D_STD = "std_depth"
D_ALL_DEPTHS = "all_depths"
D_SPACING = "depth_spacin"
D_START_IDX = "depth_start_idx"



def parse_args():
    parser = argparse.ArgumentParser("Provides statistics on a BAM/SAM file")
    parser.add_argument('--input_glob', '-i', dest='input_glob', default="*.bam", type=str,
                       help='Glob matching SAM or BAM file(s)')
    parser.add_argument('--verbose_read_length', '-l', dest='read_length', action='store_true', default=False,
                       help='Print statistics on read length for all files')
    parser.add_argument('--verbose_read_depth', '-d', dest='read_depth', action='store_true', default=False,
                       help='Print statistics on read depth for all files')

    return parser.parse_args()


def get_read_summary(read):
    summary = {
        R_START_POS: read.reference_start,
        R_END_POS: read.reference_end,
        R_LENGTH: read.reference_length,
        R_ID: read.query_name
    }
    return summary


def get_read_length_summary(read_summaries, print_it=True, length_log_base=2, length_log_max=32):
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
        L_STD: np.std(all_lengths)
    }

    if print_it: print_read_length_summary(summary, length_log_base, length_log_max)

    return summary

def print_read_length_summary(summary, length_log_base=2, length_log_max=32):
    print("\tREAD LENGTHS")
    print("\tmin: {}".format(summary[L_MIN]))
    print("\tmax: {}".format(summary[L_MAX]))
    print("\tavg: {}".format(int(summary[L_AVG])))
    print("\tstd: {}".format(int(summary[L_STD])))

    log_length_bins = summary[L_LOG_LENGTH_BUCKETS]
    print("\tread length log_{}:".format(length_log_base))
    max_bucket = max(list(filter(lambda x: log_length_bins[x] != 0, [x for x in range(length_log_max)])))
    min_bucket = min(list(filter(lambda x: log_length_bins[x] != 0, [x for x in range(length_log_max)])))
    max_bucket_size = max(log_length_bins)
    for bucket in range(min_bucket, max_bucket + 1):
        id = "%3d:" % bucket
        count = log_length_bins[bucket]
        pound_count = int(32.0 * count / max_bucket_size)
        print("\t\t{} {} {}".format(id, "#"*pound_count, count))


def get_read_depth_summary(read_summaries, print_it=True, spacing = 1000):
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

    # get depth summary
    summary = {
        D_MAX: max(depths),
        D_AVG: np.mean(depths),
        D_STD: np.std(depths),
        D_ALL_DEPTHS: depths,
        D_SPACING: spacing,
        D_START_IDX: start_idx
    }

    #print it
    if print_it: print_depth_summary(summary)

    return summary


def print_depth_summary(summary):
    print("\tREAD DEPTHS")
    print("\tmax: {}".format(summary[D_MAX]))
    print("\tavg: {}".format(summary[D_AVG]))
    print("\tstd: {}".format(summary[D_STD]))
    print("\tdepths with spacing {}:".format(summary[D_SPACING]))
    idx = summary[D_START_IDX]
    for depth in summary[D_ALL_DEPTHS]:
        id = "%4d:" % idx
        pound_count = int(32.0 * depth / summary[D_MAX])
        print("\t\t{} {} {}".format(id, '#' * pound_count, depth))
        idx += 1







def main():
    # get our arguments
    args = parse_args()

    # get filenames, sanity check
    in_alignments = glob.glob(args.input_glob)
    if len(in_alignments) == 0:
        print("No files matching {}".format(args.input_glob))
        return 1
    else:
        print("Analyzing {} files".format(len(in_alignments)))

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
        try:
            print("Handling read {}:\n\t".format(alignment_filename), end="")
            samfile = pysam.AlignmentFile(alignment_filename, 'rb' if alignment_filename.endswith("bam") else 'r')
            read_count = 0
            for read in samfile.fetch():
                read_summaries.append(get_read_summary(read))
                read_count += 1
                if read_count % 1024 == 0:
                    print(". ", end="")
            print(" {}".format(read_count))
        finally:
            if samfile is not None: samfile.close()

        # summarize
        length_summaries[alignment_filename] = get_read_length_summary(read_summaries, print_it=args.read_length)
        depth_summaries[alignment_filename] = get_read_depth_summary(read_summaries, print_it=args.read_depth)




if __name__ == "__main__":
    main()
