#!/usr/bin/env python
from __future__ import print_function
import argparse
import glob
import sys
import numpy as np
import pysam

MISSING_QUAL_CHAR = '"' # got from calling 'samtools bam2fq' on bam w/o qualities

def parse_args():
    parser = argparse.ArgumentParser("Produce FASTQ file from soft-clipped reads")
    parser.add_argument('--input', '-i', dest='input', default=None, required=True, type=str,
                       help='Alignment file')
    parser.add_argument('--output', '-o', dest='output', default=None, required=False, type=str,
                       help='Output FQ (defaults to stdout)')

    return parser.parse_args()



def main():
    args = parse_args()
    alignment_filename = args.input

    # get read data we care about
    samfile = None
    outfile = None
    read_count = -1
    try:
        print("Reading {}:".format(alignment_filename))
        samfile = pysam.AlignmentFile(alignment_filename, 'rb' if alignment_filename.endswith("bam") else 'r')
        outfile = sys.stdout if args.output is None else open(args.output, 'w')
        for read in samfile.fetch():
            read_count += 1
            read_name = read.query_name
            sequence = read.query_alignment_sequence
            qualities = read.query_alignment_qualities
            if qualities is None: qualities = MISSING_QUAL_CHAR * len(sequence)
            outfile.write("@{}\n{}\n+\n{}\n".format(read_name, sequence, qualities))
            assert(len(sequence) == len(qualities), "Read {} (#{}) has sequence len {} and qual len {}"
                   .format(read_name, read_count, len(sequence), len(qualities)))
    finally:
        if samfile is not None: samfile.close()
        if outfile is not None and args.output is not None: outfile.close()


if __name__ == "__main__":
    main()