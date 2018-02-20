#!/usr/bin/env python
from __future__ import print_function
import argparse
import glob
import gzip
import math
import numpy as np
import sys




def parse_args():
    parser = argparse.ArgumentParser("Makes plots and reports details of vcfeval output")
    parser.add_argument('--input_bed', '-i', dest='input_bed', required=True, type=str,
                       help='Input bed file')
    parser.add_argument('--chroms_to_remove', '-x', dest='chroms_to_remove', default="", type=str,
                       help='Comma-separated list of chromosomes to remove')
    return parser.parse_args()


def main():
    args = parse_args()
    args.chroms_to_remove = args.chroms_to_remove.split(",")

    depths = dict()
    with open(args.input_bed) as bed_in:
        for line in bed_in:
            line = line.split()
            if line[0] in args.chroms_to_remove: continue
            start = int(line[1])
            end = int(line[2])
            depth = int(line[3].lstrip("min_depth:"))
            if depth not in depths:
                depths[depth] = 0
            depths[depth] += end - start

    max_distance = max(depths.values())
    distance_increment = max_distance / 32
    depth_keys = list(depths.keys())
    depth_keys.sort()
    print("\nHistogram:", file=sys.stderr)
    for depth in depth_keys:
        distance = depths[depth]
        print("%6d: (%1.2f)\t%s (%d)" % (depth, distance, "#"*int(1.0*distance/distance_increment),
                                         distance), file=sys.stderr)

    max_log_distance = math.log10(max_distance)
    log_distance_increment = max_log_distance / 32
    depth_keys = list(depths.keys())
    depth_keys.sort()
    print("\nLog Histogram:", file=sys.stderr)
    for depth in depth_keys:
        distance = depths[depth]
        log_distance = math.log10(distance)
        print("%6d: (%1.2f)\t%s (%d)" % (depth, log_distance, "#"*int(1.0*log_distance/log_distance_increment),
                                         distance), file=sys.stderr)

    median_depth = None
    upper_depth = depth_keys.pop()
    upper_distance = depths[upper_depth]
    lower_depth = depth_keys[0]
    depth_keys.remove(lower_depth)
    lower_distance = depths[lower_depth]
    while True:
        if upper_distance > lower_distance:
            if len(depth_keys) == 0:
                median_depth = upper_depth
                break
            upper_distance -= lower_distance
            lower_depth = depth_keys[0]
            depth_keys.remove(lower_depth)
            lower_distance = depths[lower_depth]
        elif upper_distance < lower_distance:
            if len(depth_keys) == 0:
                median_depth = lower_depth
                break
            lower_distance -= upper_distance
            upper_depth = depth_keys.pop()
            upper_distance = depths[upper_depth]
        else:
            if len(depth_keys) >= 2:
                upper_depth = depth_keys.pop()
                upper_distance = depths[upper_depth]
                lower_depth = depth_keys[0]
                depth_keys.remove(lower_depth)
                lower_distance = depths[lower_depth]
            elif len(depth_keys) == 1:
                median_depth = depth_keys.pop()
                break
            else:
                median_depth = (upper_depth + lower_depth) / 2.0
                break

    print("\nMedian Depth: {}\n".format(median_depth), file=sys.stderr)
    print(median_depth)







if __name__ == "__main__":
    main()