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

# file names
FP_NAME = "fp.vcf"
FN_NAME = "fn.vcf"
TP_NAME = "tp.vcf"

# misc
MPI_TAG = "MPI"
DEPTH_TAG = "DP"

# chunk details
START_POS_KEY = 'start_pos'
END_POS_KEY = 'end_pos'
FP_COUNT_KEY = 'fp_count'
FN_COUNT_KEY = 'fn_count'
TP_COUNT_KEY = 'tp_count'
MPI_KEY = 'mpi'
CALL_DETAILS_LIST_KEY = 'call_details_list'

# chunk summaries
SENSITIVITY_KEY = "sensitivity"
PRECISION_KEY = "precision"
FMEASURE_KEY = "fmeasure"
AVG_QUAL_KEY = "avg_qual"
AVG_TP_QUAL_KEY = "avg_tp_qual"
AVG_FP_QUAL_KEY = "avg_fp_qual"
AVG_DEPTH_KEY = "avg_depth"
AVG_TP_DEPTH_KEY = "avg_tp_depth"
AVG_FP_DEPTH_KEY = "avg_fp_depth"

# read depth analysis
ORIGINAL_READ_COUNT_KEY = "original_read_count"
FILTERED_READ_COUNT_KEY = "filtered_read_count"
AVG_READ_DEPTH_KEY = "avg_read_depth"
STD_READ_DEPTH_KEY = "std_read_depth"
MAX_READ_DEPTH_KEY = "max_read_depths"
MIN_READ_DEPTH_KEY = "min_read_depths"
ALL_READ_DEPTH_KEY = "all_read_depths"
ALL_READ_DEPTH_POSITIONS_KEY = "all_read_depth_positions"
READ_DEPTH_MAP_KEY = "read_depth_map"

# AVG_READ_LENGTH_KEY = "avg_read_length"
# STD_READ_LENGTH_KEY = "std_read_length"
# ALL_READ_LENGTH_KEY = "all_read_length"
# MAX_READ_LENGTH_KEY = "max_read_length"
# MIN_READ_LENGTH_KEY = "min_read_length"

# call details indices
C_POS_IDX = 0
C_QUAL_IDX = 1
C_VCF_DEPTH_IDX = 2
C_TYPE_IDX = 3
C_BAM_DEPTH_IDX = 4
call_details = lambda pos, qual, depth, type: [pos, qual, depth, type, -1]

# mpi lookup position
M_START_POS_IDX = 0
M_END_POS_IDX = 1
M_MPI_IDX = 2
mpi_details = lambda chunk: [chunk[START_POS_KEY], chunk[END_POS_KEY], chunk[MPI_KEY]]

# plot info indices
P_COLOR_IDX = 0
P_LABEL_IDX = 1
P_DATA_IDX = 2
P_KEY_IDX = 3
plot_details = lambda data, key, label, color: [color, label, data, key]

# short functions
percent = lambda part, whole: int(100.0 * part / whole)
calculate_fmeasure = lambda precision, sensitivity: 2 * (precision * sensitivity) / (precision + sensitivity)

# config
DEPTH_SAMPLE_SPACING = 500
QUICK_REPORT = False

def parse_args():
    parser = argparse.ArgumentParser("Makes plots and details of vcf output")
    parser.add_argument('--vcf_directory', '-d', dest='vcf_directory', default=".", type=str,
                       help='Directory where VCF files will live')
    parser.add_argument('--unzipped', '-z', dest='unzipped', action='store_true', default=False,
                       help='Whether files are unzipped (*.vcf, not *.vcf.gz)')
    parser.add_argument('--verbose', '-v', dest='verbose', action='store_true', default=False,
                       help='Print extra information')
    parser.add_argument('--chunks', '-c', dest='chunks', default=None, type=str,
                       help='Specify comma-separated MPIs of specific chunks to analyze')
    parser.add_argument('--plot', '-p', dest='plot', action='store_true', default=False,
                       help='Make plots of data')
    parser.add_argument('--qual_min', '-q', dest='qual_min', default=None, type=int,
                       help='FILTER: Filter VCF qualities below threshold: TP->FN, FP->/dev/null')
    parser.add_argument('--vcf_read_depth_min', '-r', dest='vcf_read_depth_min', default=None, type=int,
                       help='FILTER: Filter VCF read depth below threshold: TP->FN, FP->/dev/null')
    parser.add_argument('--bam_read_depth_min', '-R', dest='bam_read_depth_min', default=None, type=int,
                       help='FILTER: Filter calls where read depth below threshold in bam. '
                            'Used in conjuction with --chunk_bam_format option. TP,FN,FP -> /dev/null')
    parser.add_argument('--remove_chunks', '-C', dest='remove_chunks', default=None, type=str,
                       help='FILTER: Remove comma-separated MPIs (or ranges) of chunks from analysis. ex: \'0-4,22\'')
    # parser.add_argument('--chunk_bam_format', '-b', dest='chunk_bam_format', default=None, type=str,
    #                    help='Python string.format() string for input bam, used in '
    #                         'conjunction with --chunks option. ex: "./SMPL1.{mpi}.in.bam"')
    parser.add_argument('--genome_bam_file', '-b', dest='genome_bam_file', default=None, type=str,
                       help='Filename of input bam (will only output depth statistics)')
    parser.add_argument('--pickle_filename', '-s', dest='pickle_filename', default=None, type=str,
                        help="Filename for pickling or unpickling read depth data.  Depth info will be unpickled from "
                             "here if the file exists, otherwise data in 'genome_bam_file' parameter will be saved.")
    # parser.add_argument('--unpickle', '-S', dest='unpickle', default=False, action='store_true',
    #                     help="If set, will unpickle read depth data from 'pickle_filename' argument.")

    return parser.parse_args()


def read_file(chunk_map, input, increment_key, args, has_mpi_tag=True):
    # get range mapping if necessary
    mpi_lookup_list = None
    current_mpi = None
    if not has_mpi_tag:
        # get positional summary
        mpi_lookup_list = list()
        for chunk in chunk_map.values():
            mpi_lookup_list.append(mpi_details(chunk))
        mpi_lookup_list.sort(key=lambda x: x[M_START_POS_IDX])
        # create mpis for what we're missing
        missing_mpis = list()
        idx = 0
        for mpi in mpi_lookup_list:
            missing_mpis.append([idx, mpi[M_START_POS_IDX] - 1, None])
            idx = mpi[M_END_POS_IDX] + 1
        missing_mpis.append([idx, sys.maxint, None])
        # add missing mpis to list of mpis
        for mpi in missing_mpis:
            mpi_lookup_list.append(mpi)
    def is_in_mpi(curr, pos):
        return pos >= curr[M_START_POS_IDX] and pos <= curr[M_END_POS_IDX]
    def lookup_mpi(pos):
        for mpi in mpi_lookup_list:
            if is_in_mpi(mpi, pos): return mpi
        assert False

    mpi_idx = 3
    depth_idx = 1
    linenr = -1
    for line in input:
        # prep
        linenr += 1
        if line.startswith("#"): continue
        line = line.split("\t")
        if len(line) < 10:
            raise Exception("Line {} is malformed in {}:\n\t[{}]".format(linenr, increment_key, '\t'.join(line)))

        # data we care about
        chr = line[0]
        pos = int(line[1])
        qual = None
        depth = None
        tags = line[8].split(":")
        smpl = line[9].split(":")
        if has_mpi_tag:
            # qual only  meaningful on our calls
            qual = int(line[5])
            # get mpi for call
            if len(tags) < mpi_idx or tags[mpi_idx] != MPI_TAG:
                mpi_idx = None
                for i in range(0, len(tags)):
                    if tags[i] == MPI_TAG: mpi_idx = i
                if mpi_idx is None: raise Exception("Line {} is missing MPI tag in {}:\n\t[{}]"
                                                    .format(linenr, increment_key, '\t'.join(line)))
            mpi = int(smpl[mpi_idx])
            # get depth
            if len(tags) < depth_idx or tags[depth_idx] != DEPTH_TAG:
                depth_idx = None
                for i in range(0, len(tags)):
                    if tags[i] == DEPTH_TAG: mpi_idx = i
                if depth_idx is None: raise Exception("Line {} is missing DP tag in {}:\n\t[{}]"
                                                    .format(linenr, increment_key, '\t'.join(line)))
            depth = int(smpl[depth_idx])
        else:
            if current_mpi is None or not is_in_mpi(current_mpi, pos):
                current_mpi = lookup_mpi(pos)
            mpi = current_mpi[M_MPI_IDX]
            if mpi is None:
                continue

        # init chunk (if necessary)
        if mpi not in chunk_map:
            chunk_map[mpi] = {
                START_POS_KEY: pos,
                END_POS_KEY: pos,
                FP_COUNT_KEY: 0,
                FN_COUNT_KEY: 0,
                TP_COUNT_KEY: 0,
                MPI_KEY: mpi,
                CALL_DETAILS_LIST_KEY: list()
            }
        chunk = chunk_map[mpi]

        # save information
        chunk[START_POS_KEY] = min(pos, chunk[START_POS_KEY])
        chunk[END_POS_KEY] = max(pos, chunk[END_POS_KEY])
        chunk[CALL_DETAILS_LIST_KEY].append(call_details(pos, qual, depth, increment_key))
        chunk[increment_key] += 1


def summarize_chunk(chunk):
    chunk[SENSITIVITY_KEY] = 1.0 * chunk[TP_COUNT_KEY] / max(1, chunk[TP_COUNT_KEY] + chunk[FN_COUNT_KEY])
    chunk[PRECISION_KEY] = 1.0 * chunk[TP_COUNT_KEY] / max(1, chunk[TP_COUNT_KEY] + chunk[FP_COUNT_KEY])
    chunk[FMEASURE_KEY] = calculate_fmeasure(max(.001, chunk[SENSITIVITY_KEY]), max(.001, chunk[PRECISION_KEY]))
    tps = list(filter(lambda x: x[C_TYPE_IDX] == TP_COUNT_KEY, chunk[CALL_DETAILS_LIST_KEY]))
    chunk[AVG_TP_QUAL_KEY] = -1 if len(tps) == 0 else int(np.mean(map(lambda x: x[C_QUAL_IDX], tps)))
    chunk[AVG_TP_DEPTH_KEY] = -1 if len(tps) == 0 else int(np.mean(map(lambda x: x[C_VCF_DEPTH_IDX], tps)))
    fps = list(filter(lambda x: x[C_TYPE_IDX] == FP_COUNT_KEY, chunk[CALL_DETAILS_LIST_KEY]))
    chunk[AVG_FP_QUAL_KEY] = -1 if len(fps) == 0 else int(np.mean(map(lambda x: x[C_QUAL_IDX], fps)))
    chunk[AVG_FP_DEPTH_KEY] = -1 if len(fps) == 0 else int(np.mean(map(lambda x: x[C_VCF_DEPTH_IDX], fps)))
    our_calls = []
    our_calls.extend(tps)
    our_calls.extend(fps)
    chunk[AVG_QUAL_KEY] = -1 if len(our_calls) == 0 else int(np.mean(map(lambda x: x[C_QUAL_IDX], our_calls)))
    chunk[AVG_DEPTH_KEY] = -1 if len(our_calls) == 0 else int(np.mean(map(lambda x: x[C_VCF_DEPTH_IDX], our_calls)))


def get_genome_read_depths(args):

    if not os.path.isfile(args.genome_bam_file):
        raise Exception("Genome bam {} does not exist!".format(args.genome_bam_file))

    import bam_stats
    length_summaries, depth_summaries = bam_stats.main([
        "--input_glob", args.genome_bam_file,
        "--depth_spacing", str(DEPTH_SAMPLE_SPACING),
        '--filter_secondary',
        '--silent'
    ])
    read_depths = list(depth_summaries.values())[0][bam_stats.D_ALL_DEPTHS]
    read_depth_positions = list(depth_summaries.values())[0][bam_stats.D_ALL_DEPTH_POSITIONS]
    assert(len(read_depths) == len(read_depth_positions))
    assert(len(read_depths) != 0)

    return {pos:depth for pos, depth in zip(read_depth_positions, read_depths)}


def fill_read_depth_summaries_for_chunks(chunk_map, genome_read_depths):

    for chunk in list(chunk_map.values()):
        chunk_depth_map = {}
        start_pos = int(1.0 * chunk[START_POS_KEY] / DEPTH_SAMPLE_SPACING)
        end_pos = int(1.0 * chunk[END_POS_KEY] / DEPTH_SAMPLE_SPACING)
        read_depths = list()
        read_depth_positions = list()
        for i in range(start_pos, end_pos + 1):
            depth = genome_read_depths[i] if i in genome_read_depths else None
            chunk_depth_map[i] = depth
            read_depth_positions.append(i)
            read_depths.append(depth)

        for call in chunk[CALL_DETAILS_LIST_KEY]:
            call_pos = int(1.0 * call[C_POS_IDX] / DEPTH_SAMPLE_SPACING)
            call[C_BAM_DEPTH_IDX] = chunk_depth_map[call_pos]


        chunk[AVG_READ_DEPTH_KEY] = np.mean(read_depths)
        chunk[STD_READ_DEPTH_KEY] = np.std(read_depths)
        chunk[MAX_READ_DEPTH_KEY] = max(read_depths)
        chunk[MIN_READ_DEPTH_KEY] = min(read_depths)
        chunk[ALL_READ_DEPTH_KEY] = read_depths
        chunk[ALL_READ_DEPTH_POSITIONS_KEY] = read_depth_positions
        chunk[READ_DEPTH_MAP_KEY] = chunk_depth_map


def filter_calls_by_bam_read_depth(chunk, read_analysis, min_depth):
    #prep
    call_data = chunk[CALL_DETAILS_LIST_KEY]
    depth_data = read_analysis[ALL_READ_DEPTH_KEY]
    depth_positions = read_analysis[ALL_READ_DEPTH_POSITIONS_KEY]
    depth_map = {pos:depth for pos,depth in zip(depth_positions, depth_data)}

    # filter calls
    new_calls = list()
    for call in call_data:
        call_pos = int(1.0 * call[C_POS_IDX] / DEPTH_SAMPLE_SPACING)
        if depth_map[call_pos] >= min_depth:
            new_calls.append(call)

    # create new chunk
    filtered_chunk = {
        START_POS_KEY: chunk[START_POS_KEY],
        END_POS_KEY: chunk[END_POS_KEY],
        FP_COUNT_KEY: len(list(filter(lambda x: x[C_TYPE_IDX] == FP_COUNT_KEY, new_calls))),
        FN_COUNT_KEY: len(list(filter(lambda x: x[C_TYPE_IDX] == FN_COUNT_KEY, new_calls))),
        TP_COUNT_KEY: len(list(filter(lambda x: x[C_TYPE_IDX] == TP_COUNT_KEY, new_calls))),
        MPI_KEY: chunk[MPI_KEY],
        CALL_DETAILS_LIST_KEY: new_calls
    }
    summarize_chunk(filtered_chunk)

    # return
    return filtered_chunk


def print_reports(reports):
    max_key_len = max(list(map(lambda x: len(x[0]), reports)))
    for report in reports:
        print(("%-"+str(max_key_len)+"s \t") % report[0], " \t".join(report[1:]))


def plot_sensitivity_x_precision(chunk_map):

    sensitivity, precision = [], []
    for chunk in chunk_map.values():
        if chunk[SENSITIVITY_KEY] == 0 or chunk[PRECISION_KEY] == 0: continue
        sensitivity.append(chunk[SENSITIVITY_KEY])
        precision.append(chunk[PRECISION_KEY])

    plt.scatter(sensitivity, precision, c="g", alpha=0.5)
    plt.xlabel("Sensitivity")
    plt.ylabel("Precision")
    plt.show()
    return plt


def plot_filtered_sensitivity_x_precision(chunk_map, threshold_values={10:'r', 25:'g', 100:'b'},
                                                  call_filter_key=C_QUAL_IDX, call_filter_name="Quality",
                                                  save_name=None):
    mpis = chunk_map.keys()
    base_tp, base_fp, base_fn = {},{},{}

    # get bases
    for mpi in mpis:
        chunk = chunk_map[mpi]
        base_tp[mpi] = len(list(filter(lambda x: x[C_TYPE_IDX] == TP_COUNT_KEY, chunk[CALL_DETAILS_LIST_KEY])))
        base_fp[mpi] = len(list(filter(lambda x: x[C_TYPE_IDX] == FP_COUNT_KEY, chunk[CALL_DETAILS_LIST_KEY])))
        base_fn[mpi] = len(list(filter(lambda x: x[C_TYPE_IDX] == FN_COUNT_KEY, chunk[CALL_DETAILS_LIST_KEY])))

    # per threshold computation
    threshold_keys = list(threshold_values.keys())
    threshold_keys.sort()
    all_tp,all_fp,all_fn = {t:{} for t in threshold_keys},{t:{} for t in threshold_keys},{t:{} for t in threshold_keys}
    for threshold in threshold_keys:
        for mpi in mpis:
            chunk = chunk_map[mpi]
            orig_tp = base_tp[mpi]
            orig_fp = base_fp[mpi]
            orig_fn = base_fn[mpi]
            tp_below_thresh = len(list(filter(lambda x: x[C_TYPE_IDX] == TP_COUNT_KEY and x[call_filter_key] < threshold,
                                              chunk[CALL_DETAILS_LIST_KEY])))
            fp_below_thresh = len(list(filter(lambda x: x[C_TYPE_IDX] == FP_COUNT_KEY and x[call_filter_key] < threshold,
                                              chunk[CALL_DETAILS_LIST_KEY])))
            new_tp = orig_tp - tp_below_thresh
            new_fp = orig_fp - fp_below_thresh
            new_fn = orig_fn + tp_below_thresh

            all_tp[threshold][mpi] = new_tp
            all_fp[threshold][mpi] = new_fp
            all_fn[threshold][mpi] = new_fn

    # per threshold plotting
    for threshold in threshold_keys:
        x, y = [], []
        for mpi in mpis:
            new_tp = all_tp[threshold][mpi]
            new_fp = all_fp[threshold][mpi]
            new_fn = all_fn[threshold][mpi]

            if (new_tp + new_fn) > 0 and (new_tp + new_fp) > 0:
                x.append(1.0 * new_tp / (new_tp + new_fn))
                y.append(1.0 * new_tp / (new_tp + new_fp))

        # calculate totals (to report on legend)
        total_thresh_tp = sum(map(lambda mpi: all_tp[threshold][mpi], mpis))
        total_thresh_fp = sum(map(lambda mpi: all_fp[threshold][mpi], mpis))
        total_thresh_fn = sum(map(lambda mpi: all_fn[threshold][mpi], mpis))
        total_sensitivity = 1.0 * total_thresh_tp / (total_thresh_tp + total_thresh_fn)
        total_precision = 1.0 * total_thresh_tp / (total_thresh_tp + total_thresh_fp)

        label = "%3d (overall: S %s, P %s)" % (threshold, ("%.2f" % total_sensitivity).lstrip('0'),
                                              ("%.2f" % total_precision).lstrip('0'))
        plt.scatter(x, y, alpha=.3, color=threshold_values[threshold], label=label)

    plt.title("Per-Chunk Sensitivity and Precision for {} Thresholds".format(call_filter_name))
    plt.xlabel("Sensitivity")
    plt.ylabel("Precision")
    plt.legend(loc=2)

    if save_name is not None: plt.savefig(save_name+"_"+("-".join(map(str, threshold_keys))))
    plt.show()


def multi_plot_chunks_by_parameter(x_plot_details, chunk_map_y, y_param, title=None, x_label=None, y_label=None,
                             save_name=None, show=True, legend_loc=1):

    for plot_detail in x_plot_details:
        data, key, color, label = plot_detail[P_DATA_IDX], plot_detail[P_KEY_IDX], \
                                  plot_detail[P_COLOR_IDX], plot_detail[P_LABEL_IDX]
        x, y = [], []
        for mpi in chunk_map_y.keys():
            if mpi not in chunk_map_y: continue
            x.append(data[mpi][key] if type(key) == str else key(data[mpi]) )
            y.append(chunk_map_y[mpi][y_param] if type(y_param) == str else y_param(chunk_map_y[mpi]) )
        plt.scatter(x, y, c=color, alpha=0.5, label=label)

    plt.legend(loc=legend_loc)
    plt.ylabel(y_param if y_label is None else y_label)
    if x_label is not None: plt.xlabel(x_label)
    if title is not None: plt.title(title)
    if save_name is not None: plt.savefig(save_name)
    if show: plt.show()


def plot_chunks_by_parameter(chunk_map_x, x_param, chunk_map_y, y_param, title=None, x_label=None, y_label=None,
                             color='r', save_name=None, show=True):

    x, y = [], []
    for mpi in chunk_map_x.keys():
        if mpi not in chunk_map_y: continue
        x.append(chunk_map_x[mpi][x_param] if type(x_param) == str else x_param(chunk_map_x[mpi]) )
        y.append(chunk_map_y[mpi][y_param] if type(y_param) == str else y_param(chunk_map_y[mpi]) )

    plt.scatter(x, y, c=color, alpha=0.5)
    plt.xlabel(x_param if x_label is None else x_label)
    plt.ylabel(y_param if y_label is None else y_label)
    plt.title("{} by {}".format(x_param, y_param) if title is None else title)
    if save_name is not None: plt.savefig(save_name)
    if show: plt.show()


def plot_all_call_qual(all_calls, bucket_size=10, logscale=True, title=None):

    tp = list(filter(lambda x: x[C_TYPE_IDX] == TP_COUNT_KEY, all_calls))
    fp = list(filter(lambda x: x[C_TYPE_IDX] == FP_COUNT_KEY, all_calls))

    bucket_count = int(100 / bucket_size) + 1
    bucket_names = []
    for i in range(bucket_count):
        start = i*bucket_size
        end = min((i+1)*bucket_size - 1, 100)
        bucket_names.append(str(end) if start == end else "{}-{}".format(start, end))
    tp_buckets = [0 for _ in range(bucket_count)]
    fp_buckets = [0 for _ in range(bucket_count)]

    def fill_buckets(calls, buckets):
        for call in calls:
            buckets[int(call[C_QUAL_IDX] / bucket_size)] += 1

    fill_buckets(tp, tp_buckets)
    fill_buckets(fp, fp_buckets)

    max_y = max(max(tp_buckets), max(fp_buckets))

    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].bar([i for i in range(bucket_count)], tp_buckets, color='blue')
    axarr[0].set_title('{}True Positives'.format("" if title is None else "{}: ".format(title)))
    axarr[0].set_xticks([i for i in range(bucket_count)])
    axarr[1].bar([i for i in range(bucket_count)], fp_buckets, color='red')
    axarr[1].set_title('{}False Positives'.format("" if title is None else "{}: ".format(title)))
    axarr[1].set_xticks([i for i in range(bucket_count)])
    axarr[1].set_xticklabels(bucket_names)
    if logscale:
        axarr[0].set_yscale("log", nonposy='clip')
        axarr[1].set_yscale("log", nonposy='clip')
        max_y *= 1.1
    else:
        max_y *= 1.05
    axarr[0].set_ylim(ymax=max_y)
    axarr[1].set_ylim(ymax=max_y)

    plt.show()
    return plt


def main():
    args = parse_args()
    do_filtering = args.qual_min is not None or \
                   args.vcf_read_depth_min is not None or \
                   args.bam_read_depth_min is not None or \
                   args.remove_chunks is not None
    do_depth_analysis = (args.genome_bam_file is not None) or (args.pickle_filename is not None)

    # get files
    fp_file = os.path.join(args.vcf_directory, FP_NAME) + ("" if args.unzipped else ".gz")
    fn_file = os.path.join(args.vcf_directory, FN_NAME) + ("" if args.unzipped else ".gz")
    tp_file = os.path.join(args.vcf_directory, TP_NAME) + ("" if args.unzipped else ".gz")
    # sanity check
    missing_files = [f for f in list(filter(lambda x: not os.path.isfile(x), [fp_file, fn_file, tp_file]))]
    if len(missing_files): raise Exception("Could not find requisite files: \n\t{}".format("\n\t".join(missing_files)))

    # data we care about
    chunk_map = dict()

    # document true positives
    print("Reading vcfeval files")
    open_fcn = open if args.unzipped else gzip.open
    with open_fcn(fp_file, 'r') as fp_in:
        read_file(chunk_map, fp_in, FP_COUNT_KEY, args)
    with open_fcn(tp_file, 'r') as tp_in:
        read_file(chunk_map, tp_in, TP_COUNT_KEY, args)
    with open_fcn(fn_file, 'r') as fn_in:
        read_file(chunk_map, fn_in, FN_COUNT_KEY, args, has_mpi_tag=False)

    # holistic analysis
    print("Analyzing all chunks")
    chunk_map_keys = list(chunk_map.keys())
    chunk_map_keys.sort()
    chunks_len = len(chunk_map_keys)
    missing_chunks = list(filter(lambda x: x not in chunk_map, [i for i in range(min(chunk_map_keys), max(chunk_map_keys))]))
    missing_chunks_len = len(missing_chunks)
    total_fp = sum(map(lambda x: x[FP_COUNT_KEY], chunk_map.values()))
    total_tp = sum(map(lambda x: x[TP_COUNT_KEY], chunk_map.values()))
    total_fn = sum(map(lambda x: x[FN_COUNT_KEY], chunk_map.values()))
    mean_calls_per_chunk = int(np.mean(map(lambda x: len(x[CALL_DETAILS_LIST_KEY]), chunk_map.values())))
    specific_chunks = [] if args.chunks is None else \
        list(map(int, args.chunks.split(","))) if args.chunks != "*" else chunk_map_keys

    # per-chunk analysis
    print("Analyzing each chunk")
    for chunk in chunk_map.values():
        summarize_chunk(chunk)

    # chunk ranking
    chunk_ranking = list(chunk_map.values())
    chunk_ranking.sort(key=lambda x: x[FMEASURE_KEY])
    chunk_ranking = list(map(lambda x: x[MPI_KEY], chunk_ranking))
    chunk_ranking_third = int(len(chunk_ranking) / 3.0)
    best_chunks = chunk_ranking[-1*min(3,len(chunk_ranking)):]
    worst_chunks = chunk_ranking[0:min(3,len(chunk_ranking))]

    # read depth analysis
    if do_depth_analysis:
        print("Read depth analysis")

        # get depths and sanity check
        genome_read_depths = None
        if args.pickle_filename is not None and os.path.isfile(args.pickle_filename):
            print("\tUnpickling read depth data from file {}".format(args.pickle_filename))
            with open(args.pickle_filename, "rb") as pickle_input:
                genome_read_depths = pickle.load(pickle_input)
            if len(genome_read_depths) == 0: print("\tUnpickled empty data!")
            else: print("\tUnpickled successfully")
        elif args.genome_bam_file is not None:
            print("Getting read depths from {}".format(args.genome_bam_file))
            genome_read_depths = get_genome_read_depths(args)
            if args.pickle_filename is not None:
                print("\tPickling read depth data into file {}".format(args.pickle_filename))
                with open(args.pickle_filename, "wb") as pickle_output:
                    pickle.dump(genome_read_depths, pickle_output)
        if genome_read_depths is None:
            raise Exception("PROGRAMMER ERROR: read depth data improperly configured")

        # populate chunks with depth information
        fill_read_depth_summaries_for_chunks(chunk_map, genome_read_depths)

    # get modified filter based on read depth (if appropriate)
    filtered_chunk_map = None
    if args.bam_read_depth_min is not None:
        filtered_chunk_map = dict()
        for chunk_id in chunk_map_keys:
            filtered_chunk_map[chunk_id] = filter_calls_by_bam_read_depth(
                chunk_map[chunk_id], chunk_map[chunk_id], args.bam_read_depth_min)
    #TODO coalesce the chunk_map stuff in this function

    # read analysis
    # read_analysis_chunks = {}
    # read_analysis_chunk_ids = []

    # overall filtering analysis
    # these functions are used to get TPs and FPs from a set of calls (from the arguments)
    true_positive_filter = lambda x: \
        (x[C_TYPE_IDX] == TP_COUNT_KEY) and \
        (x[C_QUAL_IDX] < args.qual_min if args.qual_min is not None else True) and \
        (x[C_VCF_DEPTH_IDX] < args.vcf_read_depth_min if args.vcf_read_depth_min is not None else True)
    false_positive_filter = lambda x: \
        (x[C_TYPE_IDX] == FP_COUNT_KEY) and \
        (x[C_QUAL_IDX] < args.qual_min if args.qual_min is not None else True) and \
        (x[C_VCF_DEPTH_IDX] < args.vcf_read_depth_min if args.vcf_read_depth_min is not None else True)
    if do_filtering:
        print("Filtering calls")
        # prep (if necessary)
        if filtered_chunk_map is None: filtered_chunk_map = dict()
        # get chunks to not include (if appropriate)
        chunks_to_remove = list()
        if args.remove_chunks is not None:
            remove_chunks = args.remove_chunks.split(",")
            for remove_chunk in remove_chunks:
                if '-' in remove_chunk:
                    remove_chunk = remove_chunk.split("-")
                    for i in range(int(remove_chunk[0]), int(remove_chunk[1])+1):
                        if i in chunk_map_keys: chunks_to_remove.append(i)

                else:
                    if i in chunk_map_keys: chunks_to_remove.append(int(remove_chunk))
            chunks_to_remove.sort()

        # get all calls
        all_calls = list()
        for mpi in chunk_map.keys():
            # some chunks are just skipped
            if mpi in chunks_to_remove: continue
            # get appropriate chunk (filtered by read depth in bam or not)
            chunk = filtered_chunk_map[mpi] if mpi in filtered_chunk_map else chunk_map[mpi]
            # get calls
            for call in chunk[CALL_DETAILS_LIST_KEY]:
                all_calls.append(call)
        # get baselines
        filtered_fp = len(list(filter(lambda x: x[C_TYPE_IDX] == FP_COUNT_KEY, all_calls)))
        filtered_tp = len(list(filter(lambda x: x[C_TYPE_IDX] == TP_COUNT_KEY, all_calls)))
        filtered_fn = len(list(filter(lambda x: x[C_TYPE_IDX] == FN_COUNT_KEY, all_calls)))
        if args.qual_min is not None or args.vcf_read_depth_min is not None:
            tp_below_thresh = len(list(filter(true_positive_filter, all_calls)))
            fp_below_thresh = len(list(filter(false_positive_filter, all_calls)))
            filtered_tp -= tp_below_thresh
            filtered_fn += tp_below_thresh
            filtered_fp -= fp_below_thresh

    # report
    print("\nReporting\n")
    reports = list()

    #analysis on whole dataset
    reports.append(["--------------------"])
    reports.append(["-      Summary     -"])
    reports.append(["--------------------"])
    reports.append(["Missing chunks:", "{} ({}%)".format(missing_chunks_len, percent(missing_chunks_len, missing_chunks_len + chunks_len))])
    if args.verbose:
        reports.append(["", str(missing_chunks)])
    reports.append(['Average Chunk Calls:', str(mean_calls_per_chunk)])
    reports.append(["Total TPs:", "{} ({}%)".format(total_tp, percent(total_tp, total_tp + total_fp + total_fn))])
    reports.append(["Total FPs:", "{} ({}%)".format(total_fp, percent(total_fp, total_tp + total_fp + total_fn))])
    reports.append(["Total FNs:", "{} ({}%)".format(total_fp, percent(total_fn, total_tp + total_fp + total_fn))])
    sensitivity = 1.0 * total_tp / (total_tp + total_fn)
    precision = 1.0 * total_tp / (total_tp + total_fp)
    fmeasure = calculate_fmeasure(sensitivity, precision)
    reports.append(["Sensitivity:", "%.5f" % sensitivity])
    reports.append(["Precision:", "%.5f" % precision])
    reports.append(["FMeasure:", "%.5f" % fmeasure])

    # analysis with quality filtering
    if do_filtering:
        reports.append([""])
        reports.append(["--------------------"])
        reports.append(["- Filtered Summary -"])
        reports.append(["--------------------"])
        reports.append(["Minimum Quality (VCF):", str(args.qual_min)])
        reports.append(["Minimum Depth (VCF):", str(args.vcf_read_depth_min)])
        reports.append(["Minimum Depth (BAM):", str(args.bam_read_depth_min)])
        reports.append(["Removed Chunks:", str(chunks_to_remove)])
        total_calls = sum([total_fp,total_tp,total_fn])
        filtered_calls = sum([filtered_tp,filtered_fp,filtered_fn])
        reports.append(["Remaining Calls:", "{} ({}%)".format(filtered_calls, percent(filtered_calls, total_calls))])
        reports.append(["Filtered TPs:", "{} ({}%)".format(filtered_tp, percent(filtered_tp, filtered_calls))])
        reports.append(["Filtered FPs:", "{} ({}%)".format(filtered_fp, percent(filtered_fp, filtered_calls))])
        reports.append(["Filtered FNs:", "{} ({}%)".format(filtered_fn, percent(filtered_fn, filtered_calls))])
        filt_sensitivity = 1.0 * filtered_tp / (filtered_tp + filtered_fn)
        filt_precision = 1.0 * filtered_tp / (filtered_tp + filtered_fp)
        filt_fmeasure = calculate_fmeasure(filt_sensitivity, filt_precision)
        reports.append(["Sensitivity:", "%.5f " % filt_sensitivity, "diff: %.5f" % (filt_sensitivity - sensitivity)])
        reports.append(["Precision:", "%.5f " % filt_precision, "diff: %.5f" % (filt_precision - precision)])
        reports.append(["FMeasure:", "%.5f " % filt_fmeasure, "diff: %.5f" % (filt_fmeasure - fmeasure)])

    if QUICK_REPORT:
        print_reports(reports)
        sys.exit(0)

    # analysis on a per-chunk basis
    reports.append([""])
    reports.append(["--------------------"])
    reports.append(["-  Chunk Summary   -"])
    reports.append(["--------------------"])
    def get_chunk_report(chunk, extended=True):
        sensitivity = chunk[SENSITIVITY_KEY]
        precision = chunk[PRECISION_KEY]
        fmeasure = chunk[FMEASURE_KEY]
        report = ["%8d %s%s" % (mpi,
                                  "+" if mpi in best_chunks else " ",
                                  "-" if mpi in worst_chunks else " "),
                  "fmeasure: %.3f %s " % (fmeasure,  "+  " if fmeasure > .85 else ( " = " if fmeasure > .7 else "  -")),
                  "sensitivity: %.3f %s " % (sensitivity,
                                             "+  " if sensitivity > .8 else (
                                                 " = " if sensitivity > .5 else "  -")),
                  "precision: %.3f %s " % (precision,
                                           "+  " if precision > .85 else (" = " if precision > .7 else "  -"))
            ]
        if extended:
            report.extend([
                      "total_calls: %5d  " % len(chunk[CALL_DETAILS_LIST_KEY]),
                      "avg_qual: %3d" % chunk[AVG_QUAL_KEY],
                      "avg_tp_qual: %3d" % chunk[AVG_TP_QUAL_KEY],
                      "avg_fp_qual: %3d" % chunk[AVG_FP_QUAL_KEY],
                      "avg_depth: %3d" % chunk[AVG_DEPTH_KEY],
                      "avg_tp_depth: %3d" % chunk[AVG_TP_DEPTH_KEY],
                      "avg_fp_depth: %3d" % chunk[AVG_FP_DEPTH_KEY],
                      "range: %9d-%-9d (%7d)" % (chunk[START_POS_KEY], chunk[END_POS_KEY],
                                                 chunk[END_POS_KEY] - chunk[START_POS_KEY])
                  ])
        return report
    # get reports
    for mpi in chunk_map_keys:
        chunk = chunk_map[mpi]
        if args.verbose or mpi in best_chunks or mpi in worst_chunks or mpi in specific_chunks:
            reports.append(get_chunk_report(chunk))

    # specific chunks
    def get_read_analysis_chunk_report(chunk):
        report = [
            "%8d  %s%s%s " % (chunk[MPI_KEY],
                              "+" if chunk[MPI_KEY] in chunk_ranking[chunk_ranking_third*2:] else " ",
                              "=" if chunk[MPI_KEY] in chunk_ranking[chunk_ranking_third:chunk_ranking_third*2] else " ",
                              "-" if chunk[MPI_KEY] in chunk_ranking[0:chunk_ranking_third] else " ")
        ]
        if len(chunk) > 1:
            report.extend( [
                # "read count: %6d" % chunk[ORIGINAL_READ_COUNT_KEY],
                # "filtered RC: %6d (%d%%)" % (chunk[FILTERED_READ_COUNT_KEY],
                #                              percent(chunk[FILTERED_READ_COUNT_KEY],
                #                                      chunk[ORIGINAL_READ_COUNT_KEY])),
                "depth avg: %3d" % int(chunk[AVG_READ_DEPTH_KEY]),
                "depth std: %3d" % int(chunk[STD_READ_DEPTH_KEY]),
                "depth max: %4d" % int(chunk[MAX_READ_DEPTH_KEY]),
                "depth min: %3d" % int(chunk[MIN_READ_DEPTH_KEY]),

                # "length avg: %5d" % int(chunk[AVG_READ_LENGTH_KEY]),
                # "length std: %5d" % int(chunk[STD_READ_LENGTH_KEY]),
                # "length max: %7d" % int(chunk[MAX_READ_LENGTH_KEY]),
                # "length min: %3d" % int(chunk[MIN_READ_LENGTH_KEY]),
                ])
        return report
    if do_depth_analysis:
        reports.append([""])
        reports.append(["--------------------"])
        reports.append(["-  Read  Analysis  -"])
        reports.append(["--------------------"])
        for mpi in chunk_map_keys:
            if args.verbose or mpi in best_chunks or mpi in worst_chunks or mpi in specific_chunks:
                chunk = chunk_map[mpi]
                reports.append(get_read_analysis_chunk_report(chunk))

    if do_filtering:
        reports.append([""])
        reports.append(["--------------------"])
        reports.append(["- Filtered  Chunks -"])
        reports.append(["--------------------"])
        for mpi in chunk_map_keys:
            if mpi in chunks_to_remove: continue
            if args.verbose or mpi in best_chunks or mpi in worst_chunks or mpi in specific_chunks:
                chunk = chunk_map[mpi]
                if filtered_chunk_map is not None and mpi in filtered_chunk_map:
                    chunk = filtered_chunk_map[mpi]
                tp_below_thresh = len(list(filter(true_positive_filter, chunk[CALL_DETAILS_LIST_KEY])))
                fp_below_thresh = len(list(filter(false_positive_filter, chunk[CALL_DETAILS_LIST_KEY])))
                chunk[TP_COUNT_KEY] -= tp_below_thresh
                chunk[FN_COUNT_KEY] += tp_below_thresh
                chunk[FP_COUNT_KEY] -= fp_below_thresh
                reports.append(get_chunk_report(chunk, False))


    # print them
    print_reports(reports)

    # make plots
    if args.plot:
        print("\nPlotting")

        # functions
        depth_count_below_threshold = lambda x, y: len(list(filter(lambda z: z < y, x[ALL_READ_DEPTH_KEY])))
        total_depth_below_threshold = lambda x, y: sum(list(map(lambda z: min(z, y), x[ALL_READ_DEPTH_KEY])))

        ### generic
        # plot_all_call_qual(all_calls, title="All Calls", logscale=True)
        # plot_sensitivity_x_precision(chunk_map)
        plot_filtered_sensitivity_x_precision(chunk_map, threshold_values={10:'r', 25:'g', 100:'b'},
                                              call_filter_key=C_QUAL_IDX, call_filter_name="Quality",
                                              save_name="Chunk_QualFiltered_SensPrec")
        plot_filtered_sensitivity_x_precision(chunk_map, threshold_values={0:'w', 5:'r', 20:'g', 25:'b'},
                                              call_filter_key=C_VCF_DEPTH_IDX, call_filter_name="VCF Read Depth",
                                              save_name="Chunk_VcfDepthFiltered_SensPrec")
        plot_filtered_sensitivity_x_precision(chunk_map, threshold_values={0:'w', 10:'r', 20:'g', 30:'b'},
                                              call_filter_key=C_BAM_DEPTH_IDX, call_filter_name="BAM Read Depth",
                                              save_name="Chunk_BamDepthFiltered_SensPrec")

        ### multi plot
        # plot_information = [
        #     plot_details(color='r', data=chunk_map, key=PRECISION_KEY, label="Precision"),
        #     plot_details(color='b', data=chunk_map, key=SENSITIVITY_KEY, label="Sensitivity"),
        # #    plot_details(color='g', data=chunk_map, key=HMEAN_KEY, label="Mean")
        # ]
        # multi_plot_chunks_by_parameter(plot_information, specific_chunks, STD_READ_DEPTH_KEY, legend_loc=1,
        #                          title="Depth StdDev", y_label= "Depth StdDev", save_name="All3_vs_DepthStd")
        # multi_plot_chunks_by_parameter(plot_information, specific_chunks, MIN_READ_DEPTH_KEY, legend_loc=2,
        #                          title="Depth Minimum", y_label="Depth Min", save_name="All3_vs_DepthMin")
        # multi_plot_chunks_by_parameter(plot_information, specific_chunks, lambda x: depth_count_below_threshold(x, 5),
        #                          title="Depths below 5", legend_loc=3,
        #                          y_label="Depths below 5 ", save_name="All3_vs_DepthCntLt5")
        # multi_plot_chunks_by_parameter(plot_information, specific_chunks, lambda x: total_depth_below_threshold(x, 5)/33.0,
        #                          title="Avg Depths Below 5", legend_loc=4,
        #                          y_label="Avg depth below 5", save_name="All3_vs_DepthAvgLt5")
        # multi_plot_chunks_by_parameter(plot_information, specific_chunks, lambda x: total_depth_below_threshold(x, 10)/33.0,
        #                          title="Avg Depths Below 10", legend_loc=4,
        #                          y_label="Avg depth below 10", save_name="All3_vs_DepthAvgLt10")
        # multi_plot_chunks_by_parameter(plot_information, specific_chunks, lambda x: total_depth_below_threshold(x, 20)/33.0,
        #                          title="Avg Depths Below 20", legend_loc=4,
        #                          y_label="Avg depth below 20", save_name="All3_vs_DepthAvgLt20")
        # multi_plot_chunks_by_parameter(plot_information, specific_chunks, lambda x: sum(x[ALL_READ_DEPTH_KEY])/33.0,
        #                          title="Avg Depths", legend_loc=1,
        #                          y_label="Avg depth", save_name="All3_vs_DepthAvg")


        ### precision
        # plot_chunks_by_parameter(chunk_map, PRECISION_KEY, specific_chunks, STD_READ_DEPTH_KEY,
        #                          title="Precision vs Depth StdDev", save_name="Prec_vs_DepthStd")
        # plot_chunks_by_parameter(chunk_map, PRECISION_KEY, specific_chunks, MIN_READ_DEPTH_KEY,
        #                          title="Precision vs Min Depth", save_name="Prec_vs_DepthMin", color='r')
        # plot_chunks_by_parameter(chunk_map, PRECISION_KEY, specific_chunks, lambda x: depth_count_below_threshold(x, 5),
        #                          title="Precision vs Depths below 5", x_label="Precision",
        #                          y_label="Depths below 5 ", save_name="Prec_vs_DepthCntLt5", color='r')
        # plot_chunks_by_parameter(chunk_map, PRECISION_KEY, specific_chunks, lambda x: total_depth_below_threshold(x, 5)/33.0,
        #                          title="Precision vs Avg Depths Below 5", x_label="Precision",
        #                          y_label="Avg depth below 5", save_name="Prec_vs_DepthAvgLt5", color='r')
        # plot_chunks_by_parameter(chunk_map, PRECISION_KEY, specific_chunks, lambda x: total_depth_below_threshold(x, 10)/33.0,
        #                          title="Precision vs Avg Depths Below 10", x_label="Precision",
        #                          y_label="Avg depth below 10", save_name="Prec_vs_DepthAvgLt10", color='r')
        # plot_chunks_by_parameter(chunk_map, PRECISION_KEY, specific_chunks, lambda x: total_depth_below_threshold(x, 20)/33.0,
        #                          title="Precision vs Avg Depths Below 20", x_label="Precision",
        #                          y_label="Avg depth below 20", save_name="Prec_vs_DepthAvgLt20", color='r')
        # plot_chunks_by_parameter(chunk_map, PRECISION_KEY, specific_chunks, lambda x: sum(x[ALL_READ_DEPTH_KEY])/33.0,
        #                          title="Precision vs Avg Depths", x_label="Precision",
        #                          y_label="Avg depth", save_name="Prec_vs_DepthAvg", color='r')
        #
        ### sensitivity
        # plot_chunks_by_parameter(chunk_map, SENSITIVITY_KEY, specific_chunks, STD_READ_DEPTH_KEY,
        #                          title="Sensitivity vs Depth StdDev", save_name="Sens_vs_DepthStd", color='b')
        # plot_chunks_by_parameter(chunk_map, SENSITIVITY_KEY, specific_chunks, MIN_READ_DEPTH_KEY,
        #                          title="Sensitivity vs Min Depth", save_name="Sens_vs_DepthMin", color='b')
        # plot_chunks_by_parameter(chunk_map, SENSITIVITY_KEY, specific_chunks, lambda x: depth_count_below_threshold(x, 5),
        #                          title="Sensitivity vs Depths below 5", x_label="Sensitivity",
        #                          y_label="Depths below 5 ", save_name="Sens_vs_DepthCntLt5", color='b')
        # plot_chunks_by_parameter(chunk_map, SENSITIVITY_KEY, specific_chunks, lambda x: total_depth_below_threshold(x, 5)/33.0,
        #                          title="Sensitivity vs Avg Depths Below 5", x_label="Sensitivity",
        #                          y_label="Avg depth below 5", save_name="Sens_vs_DepthAvgLt5", color='b')
        # plot_chunks_by_parameter(chunk_map, SENSITIVITY_KEY, specific_chunks, lambda x: total_depth_below_threshold(x, 10)/33.0,
        #                          title="Sensitivity vs Avg Depths Below 10", x_label="Sensitivity",
        #                          y_label="Avg depth below 10", save_name="Sens_vs_DepthAvgLt10", color='b')
        # plot_chunks_by_parameter(chunk_map, SENSITIVITY_KEY, specific_chunks, lambda x: total_depth_below_threshold(x, 20)/33.0,
        #                          title="Sensitivity vs Avg Depths Below 20", x_label="Sensitivity",
        #                          y_label="Avg depth below 20", save_name="Sens_vs_DepthAvgLt20", color='b')
        # plot_chunks_by_parameter(chunk_map, SENSITIVITY_KEY, specific_chunks, lambda x: sum(x[ALL_READ_DEPTH_KEY])/33.0,
        #                          title="Sensitivity vs Avg Depths", x_label="Sensitivity",
        #                          y_label="Avg depth", save_name="Sens_vs_DepthAvg", color='b')
        #
        ### sensitivity and precision
        # plot_chunks_by_parameter(chunk_map, HMEAN_KEY, specific_chunks, STD_READ_DEPTH_KEY,
        #                          title="Sensitivity&Precision vs Depth StdDev", save_name="SnP_vs_DepthStd", color='g')
        # plot_chunks_by_parameter(chunk_map, HMEAN_KEY, specific_chunks, MIN_READ_DEPTH_KEY,
        #                          title="Sensitivity&Precision vs Min Depth", save_name="SnP_vs_DepthMin", color='g')
        # plot_chunks_by_parameter(chunk_map, HMEAN_KEY, specific_chunks, lambda x: depth_count_below_threshold(x, 5),
        #                          title="Sensitivity&Precision vs Depths below 5", x_label="Harmonic Mean",
        #                          y_label="Depths below 5", save_name="SnP_vs_DepthCntLt5", color='g')
        # plot_chunks_by_parameter(chunk_map, HMEAN_KEY, specific_chunks, lambda x: total_depth_below_threshold(x, 5)/33.0,
        #                          title="Sensitivity&Precision vs Avg Depths Below 5", x_label="Harmonic Mean",
        #                          y_label="Avg depth below 5", save_name="SnP_vs_DepthAvgLt5", color='g')
        # plot_chunks_by_parameter(chunk_map, HMEAN_KEY, specific_chunks, lambda x: total_depth_below_threshold(x, 10)/33.0,
        #                          title="Sensitivity&Precision vs Avg Depths Below 10", x_label="Harmonic Mean",
        #                          y_label="Avg depth below 10", save_name="SnP_vs_DepthAvgLt10", color='g')
        # plot_chunks_by_parameter(chunk_map, HMEAN_KEY, specific_chunks, lambda x: total_depth_below_threshold(x, 20)/33.0,
        #                          title="Sensitivity&Precision vs Avg Depths Below 20", x_label="Harmonic Mean",
        #                          y_label="Avg depth below 20", save_name="SnP_vs_DepthAvgLt20", color='g')
        # plot_chunks_by_parameter(chunk_map, HMEAN_KEY, specific_chunks, lambda x: sum(x[ALL_READ_DEPTH_KEY])/33.0,
        #                          title="Sensitivity&Precision vs Avg Depths", x_label="Harmonic Mean",
        #                          y_label="Avg depth", save_name="SnP_vs_DepthAvg", color='g')
        #
        # if len(read_analysis_chunk_ids) > 0:
        #     for mpi in read_analysis_chunk_ids:
        #         plot_all_call_qual(chunk_map[mpi][CALL_DETAILS_LIST_KEY], title="Chunk {}".format(mpi), logscale=False)



    print("\nFin.")










if __name__ == "__main__":
    main()