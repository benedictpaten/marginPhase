#!/usr/bin/env python
from __future__ import print_function
import argparse
import gzip
import sys
import os


def parse_args():
    parser = argparse.ArgumentParser("Produce a BED based on VCF calls")
    parser.add_argument('--input_vcf', '-i', dest='input_vcf', required=True, type=str,
                       help='VCF to generate calls from (supports .gz and .vcf)')
    parser.add_argument('--output_file', '-o', dest='output_file', required=False, default=None, type=str,
                       help='Write output to file (otherwise stdout)')
    parser.add_argument('--bp_width', '-b', dest='bp_width', default=50, type=int,
                       help='Produce a bed with this many bases around each variant')
    parser.add_argument('--include_vcf', dest='include_vcf', action='store_true', default=False,
                       help='Include VCF contents in bed file')

    return parser.parse_args()


def main():
    args = parse_args()

    input_vcf = args.input_vcf
    if not os.path.isfile(input_vcf):
        raise Exception("Input VCF {} does not exist".format(input_vcf))
    width = args.bp_width

    input_open_fcn = open if not input_vcf.endswith("gz") else gzip.open
    with input_open_fcn(input_vcf, 'r') as input:
        with (open(args.output_file, 'w') if args.output_file is not None else sys.stdout) as output:
            for line in input:
                if type(line) == bytes: line = line.decode("utf-8")
                if line.startswith("#"): continue
                parts = line.split("\t")
                contig = parts[0]
                pos = int(parts[1])

                output.write("{}\t{}\t{}".format(contig, max(0, pos - width), pos + width))
                if args.include_vcf:
                    output.write("\t{}".format(line[2:]))
                else:
                    output.write("\n")


if __name__ == "__main__":
    main()