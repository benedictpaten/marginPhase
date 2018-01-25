#!/usr/bin/env python
from __future__ import print_function
import argparse
import glob
import gzip
import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import pickle
import time
import math


# bed detail indices
START_POS = "start"
END_POS = "end"
CENTER_POS = "center"
LENGTH = "length"
ABOVE_CENT = "above_cent"
DIST_TO_CENT = "cent_dist"
DIST_TO_CENT_RATIO = "centromere_dist_ratio"
DIST_TO_CHROM_END = "chrom_end_dist"
DIST_TO_CHROM_END_RATIO = "chrom_end_dist_ratio"

bed_details = lambda start, end: {START_POS:int(start), END_POS:int(end)}


chrom_sort = lambda x: int(x.replace("chr", "").replace("X", "23").replace("Y", "24"))
percent = lambda part, whole: int(100.0 * part / whole) if whole != 0 else "--"
log = print


def parse_args():
    parser = argparse.ArgumentParser("Makes plots and reports details of vcfeval output")
    parser.add_argument('--input_bed_glob', '-i', dest='input_bed_glob', default="*.bed", type=str,
                       help='Glob matching bed files for analysis')
    parser.add_argument('--chrom_bed', '-r', dest='chrom_bed', default="hg38.chrom.bed", type=str,
                       help='Chromosome bed file (describing size of chromosomes)')
    parser.add_argument('--centromere_bed', '-e', dest='centromere_bed', default="hg38.cen.bed", type=str,
                       help='Centromere bed file (describing location of centromere for each chromosome)')
    parser.add_argument('--chrom_to_remove', '-x', dest='chrom_to_remove', default="", type=str,
                       help='Comma-separated list of chromosomes to remove')
    parser.add_argument('--plot', '-p', dest='plot', action='store_true', default=False,
                       help='Make plots of data')
    parser.add_argument('--verbose', '-v', dest='verbose', action='store_true', default=False,
                       help='Print extra information')

    return parser.parse_args()


def read_bed_file(bed_filename, args):
    #prep
    bed_contents = dict()
    open_fcn = open if not bed_filename.endswith("gz") else gzip.open

    # read in bed
    with open_fcn(bed_filename, 'r') as input:
        for line in input:
            line = line.strip().split()
            if len(line) < 3: raise Exception("BED file malformed: {}, line: {}".format(bed_filename, "\t".join(line)))
            chrom = line[0]
            if chrom in args.chrom_to_remove: continue
            if chrom in bed_contents:
                chrom_bed = bed_contents[chrom]
            else:
                chrom_bed = list()
                bed_contents[chrom] = chrom_bed
            chrom_bed.append(bed_details(int(line[1]), int(line[2])))

    # sort and sanity
    for chrom_bed in bed_contents.values():
        chrom_bed.sort(key=lambda x: x[START_POS])
        last_bedframe = None
        for bedframe in chrom_bed:
            assert last_bedframe is None or last_bedframe[END_POS] <= bedframe[START_POS]

    return bed_contents


def convert_bed_to_centromere_map(bed_contents):
    centromere_map = dict()
    for chrom in bed_contents.keys():
        chrom_cent_data = bed_contents[chrom]
        centromere_map[chrom] = {
            START_POS: min(map(lambda x: x[START_POS], chrom_cent_data)),
            END_POS: min(map(lambda x: x[END_POS], chrom_cent_data)),
        }
    return centromere_map


def convert_bed_to_chrom_size_map(bed_contents):
    chrom_size_map = dict()
    for chrom in bed_contents.keys():
        chrom_size_data = bed_contents[chrom]
        if len(chrom_size_data) != 1:
            raise Exception("Chromosome size data has {} entries for {}".format(len(chrom_size_data), chrom))
        chrom_size_map[chrom] = {
            START_POS: min(map(lambda x: x[START_POS], chrom_size_data)),
            END_POS: min(map(lambda x: x[END_POS], chrom_size_data)),
        }
    return chrom_size_map


def add_bed_information(bed_contents, centromeres, chrom_sizes):
    for chrom in bed_contents.keys():
        chrom_cent = centromeres[chrom]
        chrom_size = chrom_sizes[chrom]
        chrom_calls = bed_contents[chrom]
        for call in chrom_calls:
            call[LENGTH] = call[END_POS] - call[START_POS]
            call[CENTER_POS] = call[START_POS] + int(call[LENGTH] / 2)
            if call[CENTER_POS] < chrom_cent[START_POS]:
                call[ABOVE_CENT] = False
                call[DIST_TO_CENT] = chrom_cent[START_POS] - call[CENTER_POS]
                call[DIST_TO_CHROM_END] = call[CENTER_POS] - chrom_size[START_POS]
                arm_length = chrom_cent[START_POS] - chrom_size[START_POS]
                call[DIST_TO_CENT_RATIO] = 1.0 * call[DIST_TO_CENT] / arm_length
                call[DIST_TO_CHROM_END_RATIO] = 1.0 * call[DIST_TO_CHROM_END] / arm_length
            elif call[CENTER_POS] > chrom_cent[END_POS]:
                call[ABOVE_CENT] = True
                call[DIST_TO_CENT] = call[CENTER_POS] - chrom_cent[END_POS]
                call[DIST_TO_CHROM_END] = chrom_size[END_POS] - call[CENTER_POS]
                arm_length = chrom_size[END_POS] - chrom_cent[END_POS]
                call[DIST_TO_CENT_RATIO] = 1.0 * call[DIST_TO_CENT] / arm_length
                call[DIST_TO_CHROM_END_RATIO] = 1.0 * call[DIST_TO_CHROM_END] / arm_length
            else:
                call[ABOVE_CENT] = None
                call[DIST_TO_CENT] = 0
                call[DIST_TO_CHROM_END] = min(chrom_cent[START_POS] - chrom_size[START_POS],
                                              chrom_size[END_POS] - chrom_cent[END_POS])
                call[DIST_TO_CENT_RATIO] = 0.0
                call[DIST_TO_CHROM_END_RATIO] = 1.0


def print_reports(reports):
    max_key_len = max(list(map(lambda x: len(x[0]), reports)))
    for report in reports:
        log(("%"+str(max_key_len)+"s \t") % report[0], " \t".join(report[1:]))


def plot_dist_ratios(all_data, key, bucket_count=20, title=None, save_name=None):

    raw_buckets = [0 for _ in range(bucket_count)]
    for data in all_data:
        datum = data[key]
        assert datum >= 0 and datum <= 1
        raw_buckets[int(datum * bucket_count)] += 1

    total_count = sum(raw_buckets)
    expected_count = int(total_count / bucket_count)
    count_buckets = list(map(lambda x: x - expected_count, raw_buckets))

    expected_ratio = 1.0 / bucket_count
    # ratio_buckets = list(map(lambda x: (1.0 * x / total_count) - expected_ratio, raw_buckets))
    ratio_buckets = list(map(lambda x: (1.0 * x / total_count) , raw_buckets))


    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].bar([i for i in range(bucket_count)], count_buckets, color='blue')
    if title is not None: axarr[0].set_title(title)
    axarr[0].set_xticks([i for i in range(bucket_count)])
    axarr[0].set_ylabel('Count (Expected {})'.format(expected_count))
    axarr[1].bar([i for i in range(bucket_count)], ratio_buckets, color='red')
    axarr[1].set_xticks([i for i in range(bucket_count)])
    axarr[1].set_xticklabels(list(map(str, [1.0*i/bucket_count for i in range(bucket_count)])))
    axarr[1].set_ylabel('Ratio')
    axarr[1].set_xlabel('Buckets')

    if save_name is not None: plt.savefig(save_name)
    plt.show()


def plot_dist_ratios_with_sum_key(all_data, key, sum_key, bucket_count=20, title=None, save_name=None):

    raw_buckets = [0 for _ in range(bucket_count)]
    for data in all_data:
        datum = data[key]
        assert datum >= 0 and datum <= 1
        raw_buckets[int(datum * bucket_count)] += data[sum_key]

    total_count = sum(raw_buckets)
    expected_count = int(total_count / bucket_count)
    count_buckets = list(map(lambda x: x - expected_count, raw_buckets))

    ratio_buckets = list(map(lambda x: (1.0 * x / total_count) if total_count != 0 else 0 , raw_buckets))


    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].bar([i for i in range(bucket_count)], count_buckets, color='green')
    if title is not None: axarr[0].set_title(title)
    axarr[0].set_xticks([i for i in range(bucket_count)])
    axarr[0].set_ylabel('Position Count x Length (Expected {})'.format(expected_count))
    axarr[1].bar([i for i in range(bucket_count)], ratio_buckets, color='orange')
    axarr[1].set_xticks([i for i in range(bucket_count)])
    axarr[1].set_xticklabels(list(map(str, [1.0*i/bucket_count for i in range(bucket_count)])))
    axarr[1].set_ylabel('Ratio')
    axarr[1].set_xlabel('Buckets')

    if save_name is not None: plt.savefig(save_name)
    # plt.show()
    plt.close()



def summarize_call_data(identifier, call_list, verbose):

        report = [
            "%-9s  " % identifier
        ]
        if len(call_list) == 0: return report

        def ratio_within_threshold(key, threshold, data):
            whole = len(data)
            part = len(list(filter(lambda x: x <= threshold, list(map(lambda x: x[key], data)))))
            return 1.0 * part / whole

        report.extend([
            "count: %6d" % len(call_list),
            "len_avg: %5d" % int(np.mean(map(lambda x: x[LENGTH], call_list))),
            "len_std: %4.3f" % np.std(map(lambda x: x[LENGTH], call_list)),
            "cent_dist_ratio: %.4f" % np.mean(map(lambda x: x[DIST_TO_CENT_RATIO], call_list)),
            "cent_dist_lt_10: %.3f" % ratio_within_threshold(DIST_TO_CENT_RATIO, .1, call_list),
            "chrom_dist_ratio: %.4f" % np.mean(map(lambda x: x[DIST_TO_CHROM_END_RATIO], call_list)),
            "chrom_dist_lt_10: %.3f" % ratio_within_threshold(DIST_TO_CHROM_END_RATIO, .1, call_list),
        ])

        if verbose:
            report.extend([
                "cent_dist_avg: %9d" % int(np.mean(map(lambda x: x[DIST_TO_CENT], call_list))),
                "cent_dist_std: %9.f" % np.std(map(lambda x: x[DIST_TO_CENT], call_list)),
                "chrom_dist_avg: %9d" % int(np.mean(map(lambda x: x[DIST_TO_CHROM_END], call_list))),
                "chrom_dist_std: %9.f" % np.std(map(lambda x: x[DIST_TO_CHROM_END], call_list)),
            ])
        return report


def main():
    args = parse_args()
    args.chrom_to_remove = args.chrom_to_remove.split(",")

    files = glob.glob(args.input_bed_glob)
    if args.chrom_bed in files:
        files.remove(args.chrom_bed)
    if args.centromere_bed in files:
        files.remove(args.centromere_bed)
    if len(files) == 0:
        raise Exception("No files matching {}".format(files))
    log("Found {} files".format(len(files)))

    chrom_sizes = convert_bed_to_chrom_size_map(read_bed_file(args.chrom_bed, args))
    centromeres = convert_bed_to_centromere_map(read_bed_file(args.centromere_bed, args))

    reports = list()
    file_areas = dict()

    # make image dir
    if args.plot:
        if not os.path.isdir("img"):
            os.mkdir("img")


    for file in files:
        file_name = ".".join(os.path.basename(file).split(".")[:-1])
        bed_contents = read_bed_file(file, args)
        add_bed_information(bed_contents, centromeres, chrom_sizes)
        reports.append([file])
        all_areas = list()
        all_chroms = bed_contents.keys()
        all_chroms.sort(key=chrom_sort)
        for chrom in all_chroms:
            reports.append(summarize_call_data(chrom, bed_contents[chrom], args.verbose))
            all_areas.extend(bed_contents[chrom])
            if args.plot:
                if not os.path.isdir("img/{}".format(file_name)):
                    os.mkdir("img/{}".format(file_name))
                plot_dist_ratios_with_sum_key(list(filter(lambda x: x[ABOVE_CENT], bed_contents[chrom])), DIST_TO_CENT_RATIO, LENGTH,
                                              title="{}: {} Upper Arm\n(0.0 - {}:{}, 1.0 - {}:{})".format(file_name, chrom, chrom, centromeres[chrom][END_POS], chrom, chrom_sizes[chrom][END_POS]),
                                              save_name="img/{}/CENT_ARM_LOC.BY_LENGTH.{}.{}.UPPER.png".format(file_name, file_name, chrom))
                plot_dist_ratios_with_sum_key(list(filter(lambda x: not x[ABOVE_CENT], bed_contents[chrom])), DIST_TO_CENT_RATIO, LENGTH,
                                 title="{}: {} Lower Arm\n(0.0 - {}:{}, 1.0 - {}:{})" .format(file_name, chrom, chrom, centromeres[chrom][START_POS], chrom, chrom_sizes[chrom][START_POS]),
                                              save_name="img/{}/CENT_ARM_LOC.BY_LENGTH.{}.{}.LOWER.png".format(file_name, file_name, chrom))

        reports.append(summarize_call_data("GENOME", all_areas, True))
        file_areas[file] = all_areas


    log("\nReporting:\n")
    print_reports(reports)


    if args.plot:
        log("\nPlotting\n")
        for file in files:
            areas = file_areas[file]
            file_name = ".".join(os.path.basename(file).split(".")[:-1])
            # plot_dist_ratios(areas, DIST_TO_CENT_RATIO,
            #                  title="{}\nPosition along chrom arm (0 is centromere, 1 is telomere)"
            #                  .format(file_name), save_name="CENT_ARM_LOC.{}.png".format(file_name))
            plot_dist_ratios_with_sum_key(areas, DIST_TO_CENT_RATIO, LENGTH,
                             title="{}\n(0.0 - centromere, 1.0 - telomere)"
                             .format(file_name), save_name="img/CENT_ARM_LOC.BY_LENGTH.{}.png".format(file_name))

    log("\nFin.")










if __name__ == "__main__":
    main()