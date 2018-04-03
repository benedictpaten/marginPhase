#!/usr/bin/env python
from __future__ import print_function
import argparse
import glob
import gzip
import math
import numpy as np
import os
import random


def parse_args():
    parser = argparse.ArgumentParser("Makes plots and reports details of vcfeval output")
    parser.add_argument('--input_qual_file', '-i', dest='input_qual_file', required=True, type=str,
                       help='Glob matching bed files for analysis')
    parser.add_argument('--sample_ratio', '-s', dest='sample_ratio', type=float, default=None,
                        help="Sample only this ratio of lines (for reduced memory usage on a real big BAM)")

    return parser.parse_args()


def main():
    args = parse_args()

    quals = list()
    qual_buckets = dict()
    total = 0
    with open(args.input_qual_file) as input:
        for line in input:
            total += 1
            qual = int(line)
            if qual not in qual_buckets: qual_buckets[qual] = 0
            qual_buckets[qual] += 1
            if args.sample_ratio is None or random.uniform(0, 1) < args.sample_ratio:
                quals.append(qual)


    print("Count:  {}".format(len(quals)))
    if args.sample_ratio is not None:
        print("Total:  {}".format(total))
    print("Median: {}".format(np.median(quals)))
    print("Mean:   {}".format(np.mean(quals)))
    print("StdDev: {}".format(np.std(quals)))

    print("\nHistogram:")
    qual_keys = list(qual_buckets.keys())
    qual_keys.sort()
    max_qual_value = max(qual_buckets.values())
    qual_hash_size = max_qual_value / 32.0
    for qual in qual_keys:
        value = qual_buckets[qual]
        print("%3d: %s (%9d)" % (qual, ("#" * int(value / qual_hash_size)), value))

    print("\nLog Histogram:")
    max_qual_value = math.log10(max_qual_value)
    qual_hash_size = max_qual_value / 32.0
    for qual in qual_keys:
        value = math.log10(qual_buckets[qual])
        print("%3d: %s (%2.3f)" % (qual, ("#" * int(value / qual_hash_size)), value))



if __name__ == "__main__":
    main()