#!/usr/bin/env python
from __future__ import print_function
import argparse
import glob
import gzip
import sys
import numpy as np
from scipy import stats
import os
import matplotlib.pyplot as plt
import math
import pickle
import bam_stats

chrom_sort = lambda x: int(x.replace("chr", ""))
xor = lambda a, b: (a and not b) or (not a and b)


def parse_args():
    parser = argparse.ArgumentParser("Produces BED file of bam where read depth is at or above a threshold")
    parser.add_argument('--input_glob', '-i', dest='input_glob', default="*.bam", type=str,
                        help='Glob matching SAM or BAM file(s) to produce bed file for')
    parser.add_argument('--verbose', '-v', dest='verbose', action='store_true', default=False,
                       help='Print extra information')
    parser.add_argument('--output', '-o', dest='output', default="-", type=str,
                       help="If set, output will be written to this file, otherwise stdout. "
                            "Can use python formatting with parameters:\n"
                            "file:    input filename (less extension), "
                            "chr:     chromosome, "
                            "depth:   'read_depth' parameter, "
                            "spacing: 'sampling_spacing' parameter "
                        )
    parser.add_argument('--sampling_spacing', '-s', dest='sampling_spacing', default=500, type=int,
                       help='Sampling spacing for read depth')
    parser.add_argument('--read_depth', '-d', dest='read_depth', required=True, type=int,
                       help='Minimum depth for file to be included in bed file')
    parser.add_argument('--bam_file', '-b', dest='bam_file', default=None, type=str,
                       help='Filename of input bam')
    parser.add_argument('--invert_bed', '-n', dest='invert_bed', default=False, action='store_true',
                       help='Produce BED file where read depth is BELOW \'read_depth\' parameter')

    return parser.parse_args()


def log(msg):
    print(msg, file=sys.stderr)


def get_genome_read_depths(bam_file, spacing, verbose=False):

    if not os.path.isfile(bam_file):
        raise Exception("Genome bam {} does not exist!".format(bam_file))
    log(bam_file)

    args = [
        "--input_glob", bam_file,
        "--depth_spacing", str(spacing),
        '--filter_secondary',
        '--silent'
    ]
    bam_summaries, length_summaries, depth_summaries = bam_stats.main(args)

    positions = dict()
    chromosomes = list(depth_summaries[bam_file].keys())
    chromosomes.sort(key=chrom_sort)

    for chromosome in chromosomes:
        if verbose:
            log("{}: \t{}: \tread_count: {}".format(bam_file, chromosome, bam_summaries[bam_file][chromosome][bam_stats.B_FILTERED_READ_COUNT]))
            log("{}: \t{}: \tmax_depth:  {}".format(bam_file, chromosome, depth_summaries[bam_file][chromosome][bam_stats.D_MAX]))
            log("{}: \t{}: \tmin_depth:  {}".format(bam_file, chromosome, depth_summaries[bam_file][chromosome][bam_stats.D_MIN]))
            log("{}: \t{}: \tavg_depth:  {}".format(bam_file, chromosome, depth_summaries[bam_file][chromosome][bam_stats.D_AVG]))
            log("{}: \t{}: \tstd_depth:  {}".format(bam_file, chromosome, depth_summaries[bam_file][chromosome][bam_stats.D_STD]))
        read_depths = depth_summaries[bam_file][chromosome][bam_stats.D_ALL_DEPTHS]
        read_depth_positions = depth_summaries[bam_file][chromosome][bam_stats.D_ALL_DEPTH_POSITIONS]
        assert(len(read_depths) == len(read_depth_positions))
        assert(len(read_depths) != 0)
        positions[chromosome] = {pos:depth for pos, depth in zip(read_depth_positions, read_depths)}

    return positions


def write_to_bedfile(depth_map, chromosome, output, args):
    # args
    verbose = args.verbose
    spacing = args.sampling_spacing
    min_depth = args.read_depth
    write_blocks_above_threshold = not args.invert_bed

    # prep
    max_idx = max(depth_map.keys())
    min_idx = min(depth_map.keys())
    block_start = 0
    block_total = 0
    block_below_thresh = True
    lines_written = 0

    def write_bedline(block_start, block_end, block_sum):
        start_idx = block_start * spacing
        end_idx = (block_end + 1) * spacing - 1
        avg_depth = int(1.0 * block_sum / (block_end - block_start + 1))
        output.write("{}\t{}\t{}\tavg_depth:{}\n".format(chromosome, start_idx, end_idx, avg_depth))

    # iterate over blocks
    for idx in range(min_idx, max_idx + 1):
        depth = depth_map[idx]
        below_thresh = depth < min_depth

        # are these different? (ie, have we gone from a string of belows to an above, or string of aboves to a below)
        if xor(block_below_thresh, below_thresh):
            # should we write? (ie, the block we're writing is below and we want to write the aboves)
            if xor(block_below_thresh, write_blocks_above_threshold):
                write_bedline(block_start, idx - 1, block_total)
                lines_written += 1

            # update block
            block_start = idx
            block_total = 0
            block_below_thresh = below_thresh

        block_total += depth

    # maybe write the last block
    if xor(block_below_thresh, write_blocks_above_threshold):
        write_bedline(block_start, idx, block_total)
        lines_written += 1

    # finish
    return lines_written


def main():
    # prep
    args = parse_args()
    stdout_output = args.output == '-'
    def get_output_filename(file=None, chr=None):
        assert(args.output is not None)
        return args.output.format(file=file, chr=chr, depth=args.read_depth, spacing=args.sampling_spacing)
    all_output_files = list()

    # get files
    in_alignments = glob.glob(args.input_glob)
    if len(in_alignments) == 0:
        log("No files matching {}".format(args.input_glob))
        return 1
    log("Analyzing {} files".format(len(in_alignments)))
    in_alignments.sort()

    # get depths
    all_file_depths = dict()
    for alignment_file in in_alignments:
        chrom_depths = get_genome_read_depths(alignment_file, args.sampling_spacing, args.verbose)
        all_file_depths[alignment_file] = chrom_depths

    # output beds
    output_file = sys.stdout if stdout_output else None
    output_filename = "STDOUT" if stdout_output else None
    try:
        for alignment_file in in_alignments:
            file_identifier = os.path.splitext(os.path.basename(alignment_file))[0]
            chromosomes = list(all_file_depths[alignment_file].keys())
            chromosomes.sort(key=chrom_sort)
            for chromosome in chromosomes:
                # get outstream
                if not stdout_output:
                    output_filename = get_output_filename(file_identifier, chromosome)
                    output_file = open(output_filename, 'w' if output_filename not in all_output_files else "a")
                    all_output_files.append(output_filename)

                # write beds
                depth_map = all_file_depths[alignment_file][chromosome]
                lines_written = write_to_bedfile(depth_map, chromosome, output_file, args)
                if args.verbose:
                    log("{}: wrote {} lines".format(output_filename, lines_written))

                # close outstream
                if not stdout_output:
                    output_file.close()
    finally:
        if output_file is not None and not stdout_output: output_file.close()


    log("\nFin.")










if __name__ == "__main__":
    main()