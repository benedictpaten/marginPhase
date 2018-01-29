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


# def plot_dist_ratios_by_length(all_data, key, bucket_count=BUCKET_COUNT_DEFAULT, ticks=10,
#                                title=None, save_name=None):
#     assert ticks <= BUCKET_COUNT_DEFAULT
#
#     raw_buckets = [[] for _ in range(bucket_count)]
#     raw_bucket_data = [[] for _ in range(bucket_count)]
#     bucket_sizes = list()
#     for data in all_data:
#         datum = data[key]
#         assert datum >= 0 and datum <= 1
#         bucket = int(datum * bucket_count)
#
#         # simple and naive
#         # raw_buckets[bucket] += data[sum_key]
#
#         # complicated and accurate
#         bucket_size = int(1.0 * data[ARM_LENGTH] / bucket_count)
#         bucket_sizes.append(bucket_size)
#         sum_data = data[END_POS] - data[START_POS]
#         center_distance = 0
#         max_bucket_dist = 0 if sum_data <= bucket_size else math.ceil(1.0 * (sum_data - bucket_size) / (bucket_size*2))
#         while sum_data > 0:
#             if center_distance == 0:
#                 add_to_bucket = min(sum_data, bucket_size)
#                 raw_buckets[bucket].append(add_to_bucket)
#                 raw_bucket_data[bucket].append([0, max_bucket_dist, add_to_bucket, data])
#                 sum_data -= bucket_size
#             else:
#                 add_to_bucket = bucket_size if sum_data > bucket_size * 2 else sum_data / 2
#                 raw_buckets[max(bucket - center_distance, 0)].append(add_to_bucket)
#                 raw_buckets[min(bucket + center_distance, bucket_count - 1)].append(add_to_bucket)
#                 raw_bucket_data[max(bucket - center_distance, 0)]\
#                     .append([-1 * center_distance, max_bucket_dist, add_to_bucket, data])
#                 raw_bucket_data[min(bucket + center_distance, bucket_count - 1)]\
#                     .append([center_distance, max_bucket_dist, add_to_bucket, data])
#                 sum_data -= bucket_size * 2
#             center_distance += 1
#         assert center_distance - 1 == max_bucket_dist
#
#     bucket_sums = list(map(lambda x: sum(x), raw_buckets))
#     mean_bucket_size = np.mean(bucket_sizes)
#     oversize_buckets = list(filter(lambda x: bucket_sums[x] > mean_bucket_size, [i for i in range(bucket_count)]))
#     if DEBUG:
#         print(title)
#         print("\tbucket size: mean {} (std {})".format(mean_bucket_size, np.std(bucket_sizes)))
#         print("\tbuckets over their size: {}/{}".format(len(oversize_buckets), bucket_count))
#         for idx in oversize_buckets:
#             print("\t\t{}: {}/{} (+{}, {})".format(idx, bucket_sums[idx], int(mean_bucket_size),
#                                                    int(bucket_sums[idx] - mean_bucket_size),
#                                                    1.0 * bucket_sums[idx] / mean_bucket_size ))
#             if DEBUG_VERBOSE:
#                 for b in raw_bucket_data[idx]:
#                     dpos = 3
#                     print("\t\t\t%d (%d)\t%9d (%1.4f)\t%s" % (b[0], b[1], b[2], 1.0 * b[2] / mean_bucket_size,
#                                                               {START_POS:b[dpos][START_POS], END_POS:b[dpos][END_POS],
#                                                                RATIO_CENTER_TO_CENTROMERE:("%.4f" % b[dpos][RATIO_CENTER_TO_CENTROMERE]),
#                                                                CENTER_POS:b[dpos][CENTER_POS], LENGTH:b[dpos][LENGTH]}))
#
#
#     total_count = sum(bucket_sums)
#     average_count = int(total_count / bucket_count)
#     count_buckets = list(map(lambda x: x - average_count, bucket_sums))
#
#     ratio_buckets = list(map(lambda x: (1.0 * x / total_count) if total_count != 0 else 0 , bucket_sums))
#
#
#     f, axarr = plt.subplots(2, sharex=True)
#     axarr[0].bar([i for i in range(bucket_count)], count_buckets, color='green')
#     if title is not None: axarr[0].set_title(title)
#     axarr[0].set_xticks([i for i in range(bucket_count)])
#     axarr[0].set_ylabel('Count x Length (Average {})'.format(average_count))
#     axarr[1].bar([i for i in range(bucket_count)], ratio_buckets, color='orange')
#     axarr[1].set_xticks([i*(bucket_count/ticks) for i in range(ticks)])
#     axarr[1].set_xticklabels(list(map(str, [1.0*i/ticks for i in range(ticks)])))
#     axarr[1].set_ylabel('Ratio')
#     axarr[1].set_xlabel('Buckets')
#
#     if save_name is not None: plt.savefig(save_name)
#     # plt.show()
#     plt.close()


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

    total_count = sum(bucket_value_by_size)
    average_count = total_count / len(buckets)
    count_buckets = list(map(lambda x: x - average_count, bucket_value_by_size))

    ratio_buckets = list(map(lambda x: (1.0 * x / total_count) if total_count != 0 else 0 , bucket_value_by_size))

    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].bar([i for i in range(bucket_count)], count_buckets, color='green')
    if title is not None: axarr[0].set_title(title)
    axarr[0].set_xticks([i*(bucket_count/ticks) for i in range(ticks)])
    axarr[0].set_ylabel('Count x Length (Average %.3f)' % average_count)
    axarr[1].bar([i for i in range(bucket_count)], ratio_buckets, color='orange')
    axarr[1].set_xticks([i*(bucket_count/ticks) for i in range(ticks)])
    axarr[1].set_xticklabels(list(map(str, [1.0*i/ticks for i in range(ticks)])))
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

    for file in files:
        print("\nAnalyzing {}".format(file))
        file_identifier = get_file_identifier(file)
        bed_contents = read_bed_file(file, args)
        print("\tFound {} records".format(sum([len(x) for x in bed_contents.values()])))
        add_bed_information(bed_contents, centromeres, chrom_sizes, args)
        reports = list()
        reports.append([file])
        all_areas = list()
        all_chroms = bed_contents.keys()
        all_chroms.sort(key=chrom_sort)
        print("\tAnalyzing chroms")
        for chrom in all_chroms:
            reports.append(summarize_call_data(chrom, bed_contents[chrom], args.verbose))
            all_areas.extend(bed_contents[chrom])
            # plot chromosomes
            if args.plot:
                if not os.path.isdir("{}/{}".format(IMG_DIR, file_identifier)):
                    os.mkdir("{}/{}".format(IMG_DIR, file_identifier))
                plot_coverage_in_bins(
                    list(filter(lambda x: x[ABOVE_CENT], bed_contents[chrom])),
                    create_bins_from_start_and_end_pos(0.0, 1.0, BUCKET_COUNT_DEFAULT),
                    lambda x: x[RATIO_CLOSEST_TO_CENTROMERE], lambda x: x[RATIO_FARTHEST_FROM_CENTROMERE],
                    title="{}: {} Upper Arm\n(0.0 - {}:{}, 1.0 - {}:{})"
                        .format(file_identifier, chrom, chrom, centromeres[chrom][END_POS], chrom, chrom_sizes[chrom][END_POS]),
                    save_name="{}/{}/CENT_ARM_LOC_BY_LENGTH.b{}.{}.{}.UPPER.png"
                        .format(IMG_DIR, file_identifier, BUCKET_COUNT_DEFAULT, file_identifier, chrom))
                plot_coverage_in_bins(
                    list(filter(lambda x: not x[ABOVE_CENT], bed_contents[chrom])),
                    create_bins_from_start_and_end_pos(0.0, 1.0, BUCKET_COUNT_DEFAULT),
                    lambda x: x[RATIO_CLOSEST_TO_CENTROMERE], lambda x: x[RATIO_FARTHEST_FROM_CENTROMERE],
                    title="{}: {} Lower Arm\n(0.0 - {}:{}, 1.0 - {}:{})"
                        .format(file_identifier, chrom, chrom, centromeres[chrom][START_POS], chrom, chrom_sizes[chrom][START_POS]),
                    save_name="{}/{}/CENT_ARM_LOC_BY_LENGTH.b{}.{}.{}.LOWER.png"
                        .format(IMG_DIR, file_identifier, BUCKET_COUNT_DEFAULT, file_identifier, chrom))
                # plot_dist_ratios_by_length(
                #     list(filter(lambda x: x[ABOVE_CENT], bed_contents[chrom])), DIST_TO_CENT_RATIO,
                #     title="{}: {} Upper Arm\n(0.0 - {}:{}, 1.0 - {}:{})"
                #         .format(file_identifier, chrom, chrom, centromeres[chrom][END_POS], chrom, chrom_sizes[chrom][END_POS]),
                #     save_name="{}/{}/CENT_ARM_LOC_BY_LENGTH.b{}.{}.{}.UPPER.png"
                #         .format(IMG_DIR, file_identifier, BUCKET_COUNT_DEFAULT, file_identifier, chrom))
                # plot_dist_ratios_by_length(
                #     list(filter(lambda x: not x[ABOVE_CENT], bed_contents[chrom])), DIST_TO_CENT_RATIO,
                #     title="{}: {} Lower Arm\n(0.0 - {}:{}, 1.0 - {}:{})"
                #         .format(file_identifier, chrom, chrom, centromeres[chrom][START_POS], chrom, chrom_sizes[chrom][START_POS]),
                #     save_name="{}/{}/CENT_ARM_LOC_BY_LENGTH.b{}.{}.{}.LOWER.png"
                #         .format(IMG_DIR, file_identifier, BUCKET_COUNT_DEFAULT, file_identifier, chrom))
        print("\tAnalyzing genome")
        reports.append(summarize_call_data("GENOME", all_areas, True))
        if args.plot:
            plot_coverage_in_bins(
                    all_areas, create_bins_from_start_and_end_pos(0.0, 1.0, BUCKET_COUNT_DEFAULT),
                    lambda x: x[RATIO_CLOSEST_TO_CENTROMERE], lambda x: x[RATIO_FARTHEST_FROM_CENTROMERE],
                title="{}\n(0.0 - centromere, 1.0 - telomere)".format(file_identifier),
                save_name="{}/CENT_ARM_LOC_BY_LENGTH.b{}.{}.png".format(IMG_DIR, BUCKET_COUNT_DEFAULT, file_identifier))
        log("\n\tReporting:\n")
        print_reports(reports)



    log("\nFin.")






if __name__ == "__main__":
    main()