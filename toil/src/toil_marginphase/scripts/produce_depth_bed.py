#!/usr/bin/env python
from __future__ import print_function

import argparse
import glob
import sys

import bam_stats
import os

xor = lambda a, b: (a and not b) or (not a and b)

AVG_DEPTH_PARAM = "avg_depth"
MIN_DEPTH_PARAM = "min_depth"

def parse_args():
    parser = argparse.ArgumentParser("Produces BED file of bam where read depth is at or above a threshold")
    parser.add_argument('--input_glob', '-i', dest='input_glob', default="*.bam", type=str,
                        help='Glob matching SAM or BAM file(s) to produce bed file for')
    parser.add_argument('--verbose', '-v', dest='verbose', action='store_true', default=False,
                       help='Print extra information')
    parser.add_argument('--output', '-o', dest='output', default="-", type=str,
                       help="If set, output will be written to this file, otherwise stdout. "
                            "Can use python formatting with parameters:\n"
                            "file: input filename (less extension), "
                            "chr: chromosome, "
                            "depth_threshold: 'depth_threshold' parameter, "
                            "depth_increment: 'depth_increment' parameter, "
                            "spacing: 'sampling_spacing' parameter "
                        )

    parser.add_argument('--filter_secondary', '-f', dest='filter_secondary', action='store_true', default=False,
                        help='Filter secondary alignments out (only include primary)')
    parser.add_argument('--filter_primary', '-F', dest='filter_primary', action='store_true', default=False,
                        help='Filter primary alignments out (only include secondary)')
    parser.add_argument('--sampling_spacing', '-s', dest='sampling_spacing', default=500, type=int,
                       help='Sampling spacing for read depth')
    parser.add_argument('--depth_increment', '-d', dest='depth_increment', default=None, type=int,
                       help='If set, will report min depth per \'sampling_spacing\' bases rounded down to '
                            'a multiple of \'depth_increment\' ')
    parser.add_argument('--depth_threshold', '-t', dest='depth_threshold', type=int, default=None,
                       help='Depth threshold for region ffffffto be included in bed file')
    parser.add_argument('--invert_bed', '-n', dest='invert_bed', default=False, action='store_true',
                       help='Produce BED file where read depth is BELOW \'depth_threshold\' parameter')

    args = parser.parse_args()

    if args.depth_increment is None and args.depth_threshold is None:
        raise Exception("Either 'depth_increment' or 'depth_threshold' must be specified.")
    if args.depth_increment is not None and args.depth_threshold is not None:
        raise Exception("Only one of 'depth_increment' or 'depth_threshold' may be specified.")

    return args


def log(msg):
    print(msg, file=sys.stderr)


def chrom_sort(chr):
    new_chr = chr.replace("chr", "").replace("X", "23").replace("Y", "24")
    return int(new_chr) if new_chr.isdigit() else 100


def get_genome_read_depths(bam_file, args):
    spacing = args.sampling_spacing
    verbose = args.verbose

    if not os.path.isfile(bam_file):
        raise Exception("Genome bam {} does not exist!".format(bam_file))
    log(bam_file)

    bs_args = [
        "--input_glob", bam_file,
        "--depth_spacing", str(spacing),
        '--silent'
    ]
    if args.filter_secondary: bs_args.append( '--filter_secondary')
    if args.filter_primary: bs_args.append( '--filter_primary')

    bam_summaries, _, depth_summaries = bam_stats.main(bs_args)
    if bam_stats.GENOME_KEY in depth_summaries[bam_file]: del depth_summaries[bam_file][bam_stats.GENOME_KEY]

    positions = dict()
    chromosomes = list(depth_summaries[bam_file].keys())
    chromosomes.sort(key=chrom_sort)

    for chromosome in chromosomes:
        if verbose:
            log("{}: \t{}: \tread_count: {}".format(bam_file, chromosome, bam_summaries[bam_file][chromosome][
                bam_stats.B_FILTERED_READ_COUNT]))
            log("{}: \t{}: \tmax_depth:  {}".format(bam_file, chromosome, depth_summaries[bam_file][chromosome][
                bam_stats.D_MAX]))
            log("{}: \t{}: \tmin_depth:  {}".format(bam_file, chromosome, depth_summaries[bam_file][chromosome][
                bam_stats.D_MIN]))
            log("{}: \t{}: \tavg_depth:  {}".format(bam_file, chromosome, depth_summaries[bam_file][chromosome][
                bam_stats.D_AVG]))
            log("{}: \t{}: \tstd_depth:  {}".format(bam_file, chromosome, depth_summaries[bam_file][chromosome][
                bam_stats.D_STD]))

        positions[chromosome] = depth_summaries[bam_file][chromosome][bam_stats.D_ALL_DEPTH_MAP]
        pass

    return positions


def write_to_depth_threshold_bedfile(depth_map, chromosome, output, args):
    # args
    spacing = args.sampling_spacing
    depth_threshold = args.depth_threshold
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
        end_idx = (block_end + 1) * spacing
        avg_depth = int(1.0 * block_sum / (block_end - block_start + 1))
        output.write("{}\t{}\t{}\t{}:{}\n".format(chromosome, start_idx, end_idx, AVG_DEPTH_PARAM, avg_depth))

    # iterate over blocks
    for idx in range(min_idx, max_idx + 1):
        depth = depth_map[idx]
        below_thresh = depth < depth_threshold

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


def write_to_depth_summary_bedfile(depth_map, chromosome, output, args):
    # args
    spacing = args.sampling_spacing
    depth_increment = args.depth_increment

    # prep
    max_idx = max(depth_map.keys())
    min_idx = min(depth_map.keys())
    block_start = 0
    block_depth = 0
    lines_written = 0

    def get_depth_by_increment(depth):
        return depth_increment * int(1.0 * depth / depth_increment)

    def write_bedline(block_start, block_end, block_min):
        start_idx = block_start * spacing
        end_idx = (block_end + 1) * spacing
        output.write("{}\t{}\t{}\t{}:{}\n".format(chromosome, start_idx, end_idx, MIN_DEPTH_PARAM,block_min))

    # iterate over blocks
    for idx in range(min_idx, max_idx + 1):
        current_depth = get_depth_by_increment(depth_map[idx])

        # is curret pos depth different than current block depth
        if block_depth != current_depth:
            write_bedline(block_start, idx - 1, block_depth)
            lines_written += 1

            # update block
            block_start = idx
            block_depth = current_depth


    # write the last block
    write_bedline(block_start, idx, block_depth)
    lines_written += 1

    # finish
    return lines_written


def main():
    # prep
    args = parse_args()
    stdout_output = args.output == '-'
    def get_output_filename(file=None, chr=None):
        assert(args.output is not None)
        return args.output.format(file=file, chr=chr, depth_threshold=args.depth_threshold,
                                  spacing=args.sampling_spacing, depth_increment=args.depth_increment)
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
        chrom_depths = get_genome_read_depths(alignment_file, args)
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
                if args.depth_threshold is not None:
                    lines_written = write_to_depth_threshold_bedfile(depth_map, chromosome, output_file, args)
                elif args.depth_increment is not None:
                    lines_written = write_to_depth_summary_bedfile(depth_map, chromosome, output_file, args)
                else:
                    raise Exception("PROGRAMMER ERROR: argument sanity checks failed: {}".format(args))
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