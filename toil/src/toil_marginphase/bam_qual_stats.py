#!/usr/bin/env python
from __future__ import print_function
import argparse
import glob
import gzip
import math
import numpy as np
import os
import matplotlib.pyplot as plt


from matplotlib import rcParams


bed_details = lambda start, end: {START_POS:int(start), END_POS:int(end)}

chrom_sort = lambda x: int(x.replace("chr", "").replace("X", "23").replace("Y", "24"))
percent = lambda part, whole: int(100.0 * part / whole) if whole != 0 else "--"
log = print

BUCKET_COUNT_DEFAULT = 100


DEBUG = False
DEBUG_VERBOSE = DEBUG and False

IMG_DIR = "img"
if DEBUG: IMG_DIR = "test_" + IMG_DIR

def parse_args():
    parser = argparse.ArgumentParser("Makes plots and reports details of vcfeval output")
    parser.add_argument('--input_qual_file', '-i', dest='input_qual_file', required=True, type=str,
                       help='Glob matching bed files for analysis')

    return parser.parse_args()


def main():
    args = parse_args()

    quals = list()
    with open(args.input_qual_file) as input:
        for line in input:
            quals.append(int(line))

    print("Count:  {}".format(len(quals)))
    print("Median: {}".format(np.median(quals)))
    print("Mean:   {}".format(np.mean(quals)))
    print("StdDev: {}".format(np.std(quals)))



if __name__ == "__main__":
    main()