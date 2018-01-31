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
rcParams.update({'figure.autolayout': True})

# bed detail indices
START_POS = "start"
END_POS = "end"
CENTER_POS = "center"
LENGTH = "length"
ABOVE_CENT = "above_cent"
DIST_TO_CENT = "cent_dist"
DIST_TO_CHROM_END = "chrom_end_dist"

RATIO_CENTER_TO_CENTROMERE = "ratio_cent_to_centro"
RATIO_CLOSEST_TO_CENTROMERE = "ratio_start_to_centro"
RATIO_FARTHEST_FROM_CENTROMERE = "ratio_end_to_centro"

ARM_LENGTH = "arm_length"
NAME = "name"



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
    parser.add_argument('--input_bed_glob', '-i', dest='input_bed_glob', default="*.bed", type=str,
                       help='Glob matching bed files for analysis')
    parser.add_argument('--compare_beds_to', '-c', dest='compare_beds_to', default=None, type=str,
                       help='Compare matched BEDs to this file for coverage')
    parser.add_argument('--chrom_bed', '-r', dest='chrom_bed', default="hg38.chrom.bed", type=str,
                       help='Chromosome bed file (describing size of chromosomes)')
    parser.add_argument('--centromere_bed', '-e', dest='centromere_bed', default="hg38.cen.bed", type=str,
                       help='Centromere bed file (describing location of centromere for each chromosome)')
    parser.add_argument('--chrom_to_remove', '-x', dest='chrom_to_remove', default="", type=str,
                       help='Comma-separated list of chromosomes to remove')
    parser.add_argument('--plot', '-p', dest='plot', action='store_true', default=False,
                       help='Make plots of data')
    parser.add_argument('--reduced_memory', '-m', dest='reduced_memory', action='store_true', default=False,
                       help='Reduce memory usage')
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
        centromere_map[chrom][CENTER_POS] = int((centromere_map[chrom][END_POS] + centromere_map[chrom][START_POS]) / 2)
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


def position_as_chrom_arm_ratio(chrom_size, chrom_cent, pos):
    if pos > chrom_cent[CENTER_POS]:
        arm_length = chrom_size[END_POS] - chrom_cent[END_POS]
        return 1.0 * (pos - chrom_cent[END_POS]) / arm_length
    else:
        arm_length = chrom_cent[START_POS] - chrom_size[START_POS]
        return 1.0 * (chrom_cent[START_POS] - pos) / arm_length

def add_bed_information(bed_contents, centromeres, chrom_sizes, args):
    for chrom in list(bed_contents.keys()):
        chrom_cent = centromeres[chrom]
        chrom_size = chrom_sizes[chrom]
        chrom_areas = bed_contents[chrom]
        for area in chrom_areas:
            area[LENGTH] = area[END_POS] - area[START_POS]
            area[CENTER_POS] = area[START_POS] + int((area[END_POS] - area[START_POS]) / 2)
            if area[CENTER_POS] < chrom_cent[CENTER_POS]:
                area[ABOVE_CENT] = False
                area[DIST_TO_CENT] = chrom_cent[START_POS] - area[CENTER_POS]
                area[DIST_TO_CHROM_END] = area[CENTER_POS] - chrom_size[START_POS]
                area[ARM_LENGTH] = chrom_cent[START_POS] - chrom_size[START_POS]
                area[RATIO_CENTER_TO_CENTROMERE] = position_as_chrom_arm_ratio(chrom_size, chrom_cent, area[CENTER_POS])
                area[RATIO_CLOSEST_TO_CENTROMERE] = position_as_chrom_arm_ratio(chrom_size, chrom_cent, area[END_POS])
                area[RATIO_FARTHEST_FROM_CENTROMERE] = position_as_chrom_arm_ratio(chrom_size, chrom_cent, area[START_POS])
            elif area[CENTER_POS] > chrom_cent[CENTER_POS]:
                area[ABOVE_CENT] = True
                area[DIST_TO_CENT] = area[CENTER_POS] - chrom_cent[END_POS]
                area[DIST_TO_CHROM_END] = chrom_size[END_POS] - area[CENTER_POS]
                area[ARM_LENGTH] = chrom_size[END_POS] - chrom_cent[END_POS]
                area[RATIO_CENTER_TO_CENTROMERE] = position_as_chrom_arm_ratio(chrom_size, chrom_cent, area[CENTER_POS])
                area[RATIO_CLOSEST_TO_CENTROMERE] = position_as_chrom_arm_ratio(chrom_size, chrom_cent, area[START_POS])
                area[RATIO_FARTHEST_FROM_CENTROMERE] = position_as_chrom_arm_ratio(chrom_size, chrom_cent, area[END_POS])
            else:
                raise Exception("Got call centered at centromere: {} ({})".format(area, chrom_cent))
            assert area[RATIO_CLOSEST_TO_CENTROMERE] < area[RATIO_FARTHEST_FROM_CENTROMERE] or \
                   (area[RATIO_CLOSEST_TO_CENTROMERE] > area[RATIO_FARTHEST_FROM_CENTROMERE]
                    and area[RATIO_CLOSEST_TO_CENTROMERE] < 0
                    and area[RATIO_FARTHEST_FROM_CENTROMERE] < 0)



def print_reports(reports):
    max_key_len = max(list(map(lambda x: len(x[0]), reports)))
    for report in reports:
        log(("%"+str(max_key_len)+"s \t") % report[0], " \t".join(report[1:]))

# def create_bins_from_start_and_end_pos(start_pos, end_pos, bin_count):
#     forward = start_pos < end_pos
#     coverage_length = abs(end_pos - start_pos)
#     bin_size = coverage_length / bin_count
#
#     bins = list()
#     prev_bin_end = None
#     bin_idx = 0
#     while bin_idx < bin_count:
#         bin_start = start_pos if prev_bin_end is None else prev_bin_end
#         bin_end = (bin_start + bin_size) if forward else (bin_start - bin_size)
#         bins.append({
#             START_POS: bin_start,
#             END_POS: bin_end,
#             FORWARD: forward,
#             NAME: bin_idx
#         })
#         prev_bin_end = bin_end
#         bin_idx += 1
#
#     if forward:
#         assert (end_pos - bin_count) <= prev_bin_end and (end_pos + bin_count) >= prev_bin_end
#     else:
#         assert (end_pos - bin_count) >= prev_bin_end and (end_pos + bin_count) <= prev_bin_end
#
#     return bins

def create_bins_from_start_and_end_pos(start_pos, end_pos, bin_count, prepended_bins = 0, appended_bins=0):
    coverage_length = end_pos - start_pos
    bin_size = coverage_length / bin_count

    bins = list()
    prev_bin_end = None
    bin_idx = 0
    while bin_idx < bin_count:
        bin_start = start_pos if prev_bin_end is None else prev_bin_end
        bin_end = bin_start + bin_size
        bins.append({
            START_POS: bin_start,
            END_POS: bin_end,
            NAME: str(bin_idx)
        })
        prev_bin_end = bin_end
        bin_idx += 1

    assert (end_pos - bin_count) < prev_bin_end and (end_pos + bin_count) > prev_bin_end
    assert (end_pos - bin_count) < prev_bin_end < (end_pos + bin_count)
    assert (end_pos + bin_count) > prev_bin_end > (end_pos - bin_count)

    if prepended_bins > 0:
        new_bins = []
        while prepended_bins > 0:
            bins.append({
                START_POS:-1 * bin_size * (prepended_bins + 1) + start_pos,
                END_POS: -1 * bin_size * (prepended_bins) + start_pos,
                NAME: str(-1 * prepended_bins)
            })
            prepended_bins -= 1
        bins = new_bins.extend(bins)

    if appended_bins > 0:
        new_bins = []
        bin_idx = bin_count
        while bin_idx < bin_count + appended_bins:
            bins.append({
                START_POS: bin_size * (bin_idx) + start_pos,
                END_POS: bin_size * (bin_idx + 1) + start_pos,
                NAME: str(bin_idx)
            })
            bin_idx += 1
        bins.extend(new_bins)

    return bins


def plot_coverage_in_bins(all_data, buckets, start_fcn, end_fcn, ticks=10, title=None, save_name=None,
                          non_canonical_bucket_fcn = lambda x: x[START_POS] >= 0.0 and x[END_POS <= 1.0]):
    bucket_count = len(buckets)

    #which bins does this fit in
    def find_buckets(data):
        return list(filter(lambda x: start_fcn(data) < x[END_POS] and end_fcn(data) > x[START_POS], buckets))

    all_output_buckets = {x:0 for x in map(lambda x: x[NAME], buckets)}
    for data in all_data:
        buckets_for_data = find_buckets(data)
        for bucket in buckets_for_data:
            bucket_addition = min(bucket[END_POS], end_fcn(data)) - max(bucket[START_POS], start_fcn(data))
            all_output_buckets[bucket[NAME]] += bucket_addition

    bucket_value_by_size = list()
    for bucket in buckets:
        name = bucket[NAME]
        value = all_output_buckets[name] / (bucket[END_POS] - bucket[START_POS])
        bucket_value_by_size.append(value)
    oversize_buckets = len(list(filter(lambda x: x > 1.0, bucket_value_by_size)))
    if oversize_buckets > 0:
        print("\t\tWARN: got {}/{} ({}%) oversize buckets in {}".format(oversize_buckets, len(buckets),
                                                                  percent(oversize_buckets, len(buckets)),
                                                                        save_name))

    total_count = sum(bucket_value_by_size)
    average_count = total_count / len(buckets)
    count_buckets = list(map(lambda x: x - average_count, bucket_value_by_size))

    ratio_buckets = list(map(lambda x: (1.0 * x / total_count) if total_count != 0 else 0 , bucket_value_by_size))

    f, axarr = plt.subplots(2, sharex=True)
    if title is not None: axarr[0].set_title(title)
    axarr[0].bar([i for i in range(bucket_count)], count_buckets, color='green')
    axarr[0].set_ylabel('Bucket Coverage\n(Average %.3f)' % average_count)
    axarr[0].set_xticks([i*(bucket_count/ticks) for i in range(ticks)])
    axarr[1].bar([i for i in range(bucket_count)], ratio_buckets, color='orange')
    axarr[1].set_xticks([i*(bucket_count/ticks) for i in range(ticks)])
    axarr[1].set_xticklabels(list(map(str, [1.0*i/ticks for i in range(ticks)])))
    axarr[1].set_ylabel('Ratio\n(of total coverage)')
    axarr[1].set_xlabel('Buckets')

    if save_name is not None: plt.savefig(save_name)
    # plt.show()
    plt.close()


def plot_compare_coverage_information(top_bucket, bottom_bucket, buckets, start_fcn, end_fcn, ticks=10, title=None, save_name=None,
                                      non_canonical_bucket_fcn = lambda x: x[START_POS] >= 0.0 and x[END_POS <= 1.0]):
    bucket_count = len(buckets)

    #which bins does this fit in
    def find_buckets(data):
        return list(filter(lambda x: start_fcn(data) < x[END_POS] and end_fcn(data) > x[START_POS], buckets))

    top_output_buckets = {x:0 for x in map(lambda x: x[NAME], buckets)}
    for data in top_bucket:
        buckets_for_data = find_buckets(data)
        for bucket in buckets_for_data:
            bucket_addition = min(bucket[END_POS], end_fcn(data)) - max(bucket[START_POS], start_fcn(data))
            top_output_buckets[bucket[NAME]] += bucket_addition
    bottom_output_buckets = {x:0 for x in map(lambda x: x[NAME], buckets)}
    for data in bottom_bucket:
        buckets_for_data = find_buckets(data)
        for bucket in buckets_for_data:
            bucket_addition = min(bucket[END_POS], end_fcn(data)) - max(bucket[START_POS], start_fcn(data))
            bottom_output_buckets[bucket[NAME]] += bucket_addition

    top_bucket_value_by_size = list()
    bottom_bucket_value_by_size = list()
    for bucket in buckets:
        name = bucket[NAME]
        top = top_output_buckets[name] / (bucket[END_POS] - bucket[START_POS])
        bottom = bottom_output_buckets[name] / (bucket[END_POS] - bucket[START_POS])
        top_bucket_value_by_size.append(top)
        bottom_bucket_value_by_size.append(bottom)
    oversize_top_buckets = len(list(filter(lambda x: x > 1.0, top_bucket_value_by_size)))
    if oversize_top_buckets > 0:
        print("\t\tWARN: got {}/{} ({}%) oversize TOP buckets in {}".format(oversize_top_buckets, len(buckets),
                                                                        percent(oversize_top_buckets, len(buckets)),
                                                                        save_name))
    oversize_bottom_buckets = len(list(filter(lambda x: x > 1.0, bottom_bucket_value_by_size)))
    if oversize_bottom_buckets > 0:
        print("\t\tWARN: got {}/{} ({}%) oversize BOTTOM buckets in {}".format(oversize_bottom_buckets, len(buckets),
                                                                        percent(oversize_bottom_buckets, len(buckets)),
                                                                        save_name))

    top_average = 1.0 * sum(top_bucket_value_by_size) / len(buckets)
    if top_average == 0: top_average = 1.0
    bottom_average = 1.0 * sum(bottom_bucket_value_by_size) / len(buckets)
    if bottom_average == 0: bottom_average = 1.0
    count_buckets = map(lambda x: x[0] - x[1], zip(top_bucket_value_by_size, bottom_bucket_value_by_size))
    ratio_buckets = map(lambda x: (x[0]/top_average) - (x[1]/bottom_average),
                        zip(top_bucket_value_by_size, bottom_bucket_value_by_size))


    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].bar([i for i in range(bucket_count)], map(lambda x: 0 if x < 0 else x, count_buckets), color='blue')
    axarr[0].bar([i for i in range(bucket_count)], map(lambda x: 0 if x > 0 else x, count_buckets), color='red')
    if title is not None: axarr[0].set_title(title)
    axarr[0].set_xticks([i*(bucket_count/ticks) for i in range(ticks)])
    axarr[0].set_ylabel('Coverage Difference\n(ratio of bucket coverage)')
    axarr[1].bar([i for i in range(bucket_count)], map(lambda x: 0 if x < 0 else x, ratio_buckets), color='blue')
    axarr[1].bar([i for i in range(bucket_count)], map(lambda x: 0 if x > 0 else x, ratio_buckets), color='red')
    axarr[1].set_xticks([i*(bucket_count/ticks) for i in range(ticks)])
    axarr[1].set_xticklabels(list(map(str, [1.0*i/ticks for i in range(ticks)])))
    axarr[1].set_ylabel('Normalized Difference\n(by mean coverage)')
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
            "cent_dist_ratio: %.4f" % np.mean(map(lambda x: x[RATIO_CENTER_TO_CENTROMERE], call_list)),
            "cent_dist_lt_10: %.3f" % ratio_within_threshold(RATIO_CENTER_TO_CENTROMERE, .1, call_list),
        ])

        if verbose:
            report.extend([
                "cent_dist_avg: %9d" % int(np.mean(map(lambda x: x[DIST_TO_CENT], call_list))),
                "cent_dist_std: %9.f" % np.std(map(lambda x: x[DIST_TO_CENT], call_list)),
                "chrom_dist_avg: %9d" % int(np.mean(map(lambda x: x[DIST_TO_CHROM_END], call_list))),
                "chrom_dist_std: %9.f" % np.std(map(lambda x: x[DIST_TO_CHROM_END], call_list)),
            ])
        return report


def get_file_identifier(filename):
    file_name = ".".join(os.path.basename(filename).split(".")[:-1])
    file_identifier = file_name
    if len(file_name) > 36:
        file_identifier = "{}..{}".format(file_name[:16],file_name[-16:])
    return file_identifier



def main():
    args = parse_args()
    args.chrom_to_remove = args.chrom_to_remove.split(",")

    if args.compare_beds_to is not None and not os.path.isfile(args.compare_beds_to):
        raise Exception("Compared BED not found: {}".format(args.compare_beds_to))
    files = glob.glob(args.input_bed_glob)
    if args.chrom_bed in files:
        files.remove(args.chrom_bed)
    if args.centromere_bed in files:
        files.remove(args.centromere_bed)
    if len(files) == 0:
        raise Exception("No files matching {}".format(args.input_bed_glob))
    log("Found {} files".format(len(files)))

    chrom_sizes = convert_bed_to_chrom_size_map(read_bed_file(args.chrom_bed, args))
    centromeres = convert_bed_to_centromere_map(read_bed_file(args.centromere_bed, args))

    # make image dir
    if args.plot:
        if not os.path.isdir(IMG_DIR):
            os.mkdir(IMG_DIR)

    if args.compare_beds_to is not None:
        tmp = [args.compare_beds_to]
        tmp.extend(files)
        files = tmp

    compare_bed_identifier = None
    compare_bed_chroms = None
    compare_bed_all_areas = None
    for file in files:
        print("\nAnalyzing {}".format(file))
        if args.compare_beds_to is not None and compare_bed_identifier is None:
            print("\tCOMPARE BED")
        file_identifier = get_file_identifier(file)
        bed_contents = read_bed_file(file, args)
        print("\tFound {} records".format(sum([len(x) for x in bed_contents.values()])))
        add_bed_information(bed_contents, centromeres, chrom_sizes, args)
        reports = list()
        reports.append([file])
        all_areas = list()
        all_chroms = bed_contents.keys()
        all_chroms.sort(key=chrom_sort)

        #plot prep
        bucket_ranges = create_bins_from_start_and_end_pos(0.0, 1.0, BUCKET_COUNT_DEFAULT)
        above_cent_filter = lambda x: list(filter(lambda y: y[ABOVE_CENT], x))
        below_cent_filter = lambda x: list(filter(lambda y: not y[ABOVE_CENT], x))
        ratio_start_fcn = lambda x: x[RATIO_CLOSEST_TO_CENTROMERE]
        ratio_end_fcn = lambda x: x[RATIO_FARTHEST_FROM_CENTROMERE]

        print("\tAnalyzing chroms")
        for chrom in all_chroms:
            reports.append(summarize_call_data(chrom, bed_contents[chrom], args.verbose))
            all_areas.extend(bed_contents[chrom])
            # plot chromosomes
            if args.plot:
                if not os.path.isdir("{}/{}".format(IMG_DIR, file_identifier)):
                    os.mkdir("{}/{}".format(IMG_DIR, file_identifier))

                # coverage by chrom arm
                plot_coverage_in_bins(
                    above_cent_filter(bed_contents[chrom]), bucket_ranges, ratio_start_fcn, ratio_end_fcn,
                    title="{}: {} Upper Arm\n(0.0 - {}:{}, 1.0 - {}:{})"
                        .format(file_identifier, chrom, chrom, centromeres[chrom][END_POS], chrom, chrom_sizes[chrom][END_POS]),
                    save_name="{}/{}/CENT_ARM_LOC_BY_LENGTH.b{}.{}.{}.UPPER.png"
                        .format(IMG_DIR, file_identifier, BUCKET_COUNT_DEFAULT, file_identifier, chrom))
                plot_coverage_in_bins(
                    below_cent_filter(bed_contents[chrom]), bucket_ranges, ratio_start_fcn, ratio_end_fcn,
                    title="{}: {} Lower Arm\n(0.0 - {}:{}, 1.0 - {}:{})"
                        .format(file_identifier, chrom, chrom, centromeres[chrom][START_POS], chrom, chrom_sizes[chrom][START_POS]),
                    save_name="{}/{}/CENT_ARM_LOC_BY_LENGTH.b{}.{}.{}.LOWER.png"
                        .format(IMG_DIR, file_identifier, BUCKET_COUNT_DEFAULT, file_identifier, chrom))

                # comparative coverage by chrom arm
                if compare_bed_identifier is not None:
                    plot_compare_coverage_information(
                        above_cent_filter(bed_contents[chrom]), above_cent_filter(compare_bed_chroms[chrom]),
                        bucket_ranges, ratio_start_fcn, ratio_end_fcn,
                        title="{} (Top/Blue)\n{} (Bottom/Red)\n{} Upper Arm\n(0.0 - {}:{}, 1.0 - {}:{})"
                            .format(file_identifier, compare_bed_identifier, chrom, chrom, centromeres[chrom][END_POS], chrom, chrom_sizes[chrom][END_POS]),
                        save_name="{}/{}/COMPARE_CENT_ARM.b{}.T-{}.B-{}.{}.UPPER.png"
                            .format(IMG_DIR, file_identifier, BUCKET_COUNT_DEFAULT, file_identifier, compare_bed_identifier, chrom))
                    plot_compare_coverage_information(
                        below_cent_filter(bed_contents[chrom]), below_cent_filter(compare_bed_chroms[chrom]),
                        bucket_ranges, ratio_start_fcn, ratio_end_fcn,
                        title="{} (Top/Blue)\n{} (Bottom/Red)\n{} Lower Arm\n(0.0 - {}:{}, 1.0 - {}:{})"
                            .format(file_identifier, compare_bed_identifier, chrom, chrom, centromeres[chrom][END_POS], chrom, chrom_sizes[chrom][END_POS]),
                        save_name="{}/{}/COMPARE_CENT_ARM.b{}.T-{}.B-{}.{}.LOWER.png"
                            .format(IMG_DIR, file_identifier, BUCKET_COUNT_DEFAULT, file_identifier, compare_bed_identifier, chrom))


        print("\tAnalyzing genome")
        reports.append(summarize_call_data("GENOME", all_areas, True))
        log("\tReporting")
        print_reports(reports)

        if args.plot:
            plot_coverage_in_bins(
                    all_areas, bucket_ranges, ratio_start_fcn, ratio_end_fcn,
                title="{}\n(0.0 - centromere, 1.0 - telomere)".format(file_identifier),
                save_name="{}/CENT_ARM_LOC_BY_LENGTH.b{}.{}.png".format(IMG_DIR, BUCKET_COUNT_DEFAULT, file_identifier))
            if compare_bed_identifier is not None:
                plot_compare_coverage_information(
                    all_areas, compare_bed_all_areas,
                    bucket_ranges, ratio_start_fcn, ratio_end_fcn,
                    title="{} (Top/Blue)\n{} (Bottom/Red)\n(0.0 - centromere, 1.0 - telomere)"
                        .format(file_identifier, compare_bed_identifier),
                    save_name="{}/COMPARE_CENT_ARM.b{}.T-{}.B-{}.genome.png"
                        .format(IMG_DIR, BUCKET_COUNT_DEFAULT, file_identifier, compare_bed_identifier))

        # set the compare bed options if appropriate
        if args.compare_beds_to is not None and compare_bed_identifier is None:
            compare_bed_identifier = file_identifier
            compare_bed_chroms = bed_contents
            compare_bed_all_areas = all_areas



    log("\nFin.")






if __name__ == "__main__":
    main()