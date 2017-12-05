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
HMEAN_KEY = "hmean"
AVG_QUAL_KEY = "avg_qual"
AVG_TP_QUAL_KEY = "avg_tp_qual"
AVG_FP_QUAL_KEY = "avg_fp_qual"


# chunk analysis
ORIGINAL_READ_COUNT_KEY = "original_read_count"
FILTERED_READ_COUNT_KEY = "filtered_read_count"
AVG_READ_DEPTH_KEY = "avg_read_depth"
STD_READ_DEPTH_KEY = "std_read_depth"
ALL_READ_DEPTH_KEY = "all_read_depths"
MAX_READ_DEPTH_KEY = "max_read_depths"
MIN_READ_DEPTH_KEY = "min_read_depths"

AVG_READ_LENGTH_KEY = "avg_read_length"
STD_READ_LENGTH_KEY = "std_read_length"
ALL_READ_LENGTH_KEY = "all_read_length"
MAX_READ_LENGTH_KEY = "max_read_length"
MIN_READ_LENGTH_KEY = "min_read_length"

# call details indices
C_POS_IDX = 0
C_QUAL_IDX = 1
C_TYPE_IDX = 2
call_details = lambda pos, qual, type: [pos, qual, type]

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


def parse_args():
    parser = argparse.ArgumentParser("Makes plots and details of vcf output")
    parser.add_argument('--vcf_directory', '-d', dest='vcf_directory', default=".", type=str,
                       help='Directory where VCF files will live')
    parser.add_argument('--unzipped', '-z', dest='unzipped', action='store_true', default=False,
                       help='Whether files are unzipped (*.vcf, not *.vcf.gz)')
    parser.add_argument('--verbose', '-v', dest='verbose', action='store_true', default=False,
                       help='Print extra information')
    parser.add_argument('--plot', '-p', dest='plot', action='store_true', default=False,
                       help='Make plots of data')
    parser.add_argument('--qual_min', '-q', dest='qual_min', default=None, type=int,
                       help='Filter qualities below threshold: TP->FN, FP->/dev/null')
    parser.add_argument('--chunks', '-c', dest='chunks', default=None, type=str,
                       help='Specify comma-separated MPIs of specific chunks to analyze')
    parser.add_argument('--in_bam_format', '-b', dest='in_bam_format', default=None, type=str,
                       help='Python string.format() string for input bam, used (optionally) in '
                            'conjunction with --chunks option. ex: "./SMPL1.{mpi}.in.bam"')
    parser.add_argument('--pickle_filename', '-s', dest='pickle_filename', default=None, type=str,
                        help="Filename for pickling specific_chunk read data.  No pickling will happen if unset.")
    parser.add_argument('--unpickle', '-S', dest='unpickle', default=False, action='store_true',
                        help="If set, will unpickle specific_chunk data from 'pickle_filename' argument.")

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
        qual = int(line[5])
        tags = line[8].split(":")
        smpl = line[9].split(":")
        if has_mpi_tag:
            if len(tags) < mpi_idx or tags[mpi_idx] != MPI_TAG:
                mpi_idx = None
                for i in range(0, len(tags)):
                    if tags[i] == MPI_TAG: mpi_idx = i
                if mpi_idx is None: raise Exception("Line {} is missing MPI tag in {}:\n\t[{}]"
                                                    .format(linenr, increment_key, '\t'.join(line)))
            mpi = int(smpl[mpi_idx])
        else:
            if current_mpi is None or not is_in_mpi(current_mpi, pos):
                current_mpi = lookup_mpi(pos)
            mpi = current_mpi[M_MPI_IDX]
            if mpi is None:
                continue;


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
        chunk[END_POS_KEY] = pos
        chunk[CALL_DETAILS_LIST_KEY].append(call_details(pos, qual, increment_key))
        chunk[increment_key] += 1


def analyze_specific_chunk(chunk, args):
    analysis = {MPI_KEY: chunk[MPI_KEY]}
    if args.in_bam_format is None:
        return analysis
    bam_location = args.in_bam_format.format(mpi=chunk[MPI_KEY])
    if not os.path.isfile(bam_location):
        print("\tCould not find bam for chunk {}: {}".format(chunk[MPI_KEY], bam_location))
        return analysis

    print("\tAnalyzing {}".format(bam_location))

    import bam_stats
    length_summaries, depth_summaries = bam_stats.main([
        "--input_glob", bam_location,
        "--depth_range", "{}-{}".format(chunk[START_POS_KEY], chunk[END_POS_KEY]),
        "--depth_spacing", str(int((chunk[END_POS_KEY] - chunk[START_POS_KEY]) / 32)),
        '--filter_secondary',
        '--min_alignment_threshold', str(20),
        '--silent'
    ])
    read_lengths = list(length_summaries.values())[0][bam_stats.L_ALL_LENGTHS]
    read_depths = list(depth_summaries.values())[0][bam_stats.D_ALL_DEPTHS]

    analysis[AVG_READ_DEPTH_KEY] = np.mean(read_depths)
    analysis[STD_READ_DEPTH_KEY] = np.std(read_depths)
    analysis[ALL_READ_DEPTH_KEY] = read_depths
    analysis[MAX_READ_DEPTH_KEY] = max(read_depths)
    analysis[MIN_READ_DEPTH_KEY] = min(read_depths)

    analysis[AVG_READ_LENGTH_KEY] = np.mean(read_lengths)
    analysis[STD_READ_LENGTH_KEY] = np.std(read_lengths)
    analysis[ALL_READ_LENGTH_KEY] = read_lengths
    analysis[MAX_READ_LENGTH_KEY] = max(read_lengths)
    analysis[MIN_READ_LENGTH_KEY] = min(read_lengths)

    analysis[ORIGINAL_READ_COUNT_KEY] = list(length_summaries.values())[0][bam_stats.L_READ_COUNT]
    analysis[FILTERED_READ_COUNT_KEY] = list(length_summaries.values())[0][bam_stats.L_FILTERED_READ_COUNT]

    return analysis


def print_reports(reports):
    max_key_len = max(list(map(lambda x: len(x[0]), reports)))
    for report in reports:
        print(("%-"+str(max_key_len)+"s\t") % report[0], "\t".join(report[1:]))


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

# {10:'r', 20:'g', 25:'c', 30:'b', 40:'m', 100:'k'}
# {10:'r', 20:'g', 30:'b', 100:'k'}
# {10:'r', 25:'g', 40:'b', 100:'k'}
# {10:'r', 25:'g', 100:'b'}
def plot_quality_filtered_sensitivity_x_precision(chunk_map, threshold_values={10:'r', 25:'g', 100:'b'}, save_name=None):
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
            tp_below_thresh = len(list(filter(lambda x: x[C_TYPE_IDX] == TP_COUNT_KEY and x[C_QUAL_IDX] < threshold,
                                              chunk[CALL_DETAILS_LIST_KEY])))
            fp_below_thresh = len(list(filter(lambda x: x[C_TYPE_IDX] == FP_COUNT_KEY and x[C_QUAL_IDX] < threshold,
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

    plt.title("Per-Chunk Sensitivity and Precision for Quality Thresholds")
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

    # get all calls
    all_calls = list()
    for chunk in chunk_map.values():
        for call in chunk[CALL_DETAILS_LIST_KEY]:
            all_calls.append(call)

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

    # all calls analysis
    print("Analyzing all calls")
    filtered_fp = total_fp
    filtered_tp = total_tp
    filtered_fn = total_fn
    if args.qual_min is not None:
        tp_below_thresh = len(list(filter(
            lambda x: x[C_QUAL_IDX] < args.qual_min and x[C_TYPE_IDX] == TP_COUNT_KEY, all_calls)))
        fp_below_thresh = len(list(filter(
            lambda x: x[C_QUAL_IDX] < args.qual_min and x[C_TYPE_IDX] == FP_COUNT_KEY, all_calls)))
        filtered_tp -= tp_below_thresh
        filtered_fn += tp_below_thresh
        filtered_fp -= fp_below_thresh

    # per-chunk analysis
    print("Analyzing each chunk")
    for chunk in chunk_map.values():
        chunk[SENSITIVITY_KEY] = 1.0 * chunk[TP_COUNT_KEY] / max(1, chunk[TP_COUNT_KEY] + chunk[FP_COUNT_KEY])
        chunk[PRECISION_KEY] = 1.0 * chunk[TP_COUNT_KEY] / max(1, chunk[TP_COUNT_KEY] + chunk[FN_COUNT_KEY])
        chunk[HMEAN_KEY] = stats.hmean([max(.001, chunk[SENSITIVITY_KEY]), max(.001, chunk[PRECISION_KEY])])
        chunk[AVG_QUAL_KEY] = int(np.mean(map(lambda x: x[C_QUAL_IDX], chunk[CALL_DETAILS_LIST_KEY])))
        tps = list(filter(lambda x: x[C_TYPE_IDX] == TP_COUNT_KEY, chunk[CALL_DETAILS_LIST_KEY]))
        chunk[AVG_TP_QUAL_KEY] = -1 if len(tps) == 0 else int(np.mean(map(lambda x: x[C_QUAL_IDX], tps)))
        fps = list(filter(lambda x: x[C_TYPE_IDX] == FP_COUNT_KEY, chunk[CALL_DETAILS_LIST_KEY]))
        chunk[AVG_FP_QUAL_KEY] = -1 if len(fps) == 0 else int(np.mean(map(lambda x: x[C_QUAL_IDX], fps)))

    # all chunks analysis
    chunk_ranking = list(chunk_map.values())
    chunk_ranking.sort(key=lambda x: x[HMEAN_KEY])
    chunk_ranking = list(map(lambda x: x[MPI_KEY], chunk_ranking))
    chunk_ranking_third = int(len(chunk_ranking) / 3.0)
    best_chunks = chunk_ranking[-1*min(3,len(chunk_ranking)):]
    worst_chunks = chunk_ranking[0:min(3,len(chunk_ranking))]

    #specific chunk analysis
    specific_chunks = {}
    specific_chunk_ids = []
    if args.chunks is not None:
        print("Analyzing specific chunks")
        specific_chunk_ids = map(int, args.chunks.split(",")) if args.chunks != '*' else chunk_map_keys
        # unpickle
        if args.unpickle and os.path.isfile(args.pickle_filename):
            print("\tUnpickling chunk data from file {}".format(args.pickle_filename))
            specific_chunks = pickle.load(open(args.pickle_filename, "rb"))
            print("\tUnpickled {} chunks".format(len(specific_chunks)))
            unpickled_chunk_ids = list(map(lambda x: x[MPI_KEY], specific_chunks.values()))
            for specific_chunk_id in specific_chunk_ids:
                if specific_chunk_id not in unpickled_chunk_ids:
                    print("\t\tUnpickled data does not include chunk {}".format(specific_chunk_id))
                    specific_chunks[specific_chunk_id] = {MPI_KEY: specific_chunk_id}
        else:
            if args.unpickle: #pickle_filname didn't exist
                print("\tPickle filename '{}' doesn't exist to unpickle from.  Importing data manually."
                      .format(args.pickle_filename))
            for mpi in specific_chunk_ids:
                specific_chunks[mpi] = analyze_specific_chunk(chunk_map[mpi], args)

        # pickle data if appropriate (pickle filename is set and we didn't just unpickle it)
        if args.pickle_filename is not None and not (args.unpickle and os.path.isfile(args.pickle_filename)):
            print("\tPickling chunk data into file {}".format(args.pickle_filename))
            pickle.dump(specific_chunks, open(args.pickle_filename, "wb"))

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
    reports.append(["Sensitivity:", "{}".format(1.0 * total_tp / (total_tp + total_fp))])
    reports.append(["Precision:", "{}".format(1.0 * total_tp / (total_tp + total_fn))])

    # analysis with quality filtering
    if args.qual_min is not None:
        reports.append(["--------------------"])
        reports.append(["- Filtered Summary -"])
        reports.append(["--------------------"])
        reports.append(["Minimum Quality:", str(args.qual_min)])
        reports.append(["Total TPs:", "{} ({}%)".format(filtered_tp, percent(filtered_tp, filtered_tp + filtered_fp + filtered_fn))])
        reports.append(["Total FPs:", "{} ({}%)".format(total_fp, percent(total_fp, filtered_tp + filtered_fp + filtered_fn))])
        reports.append(["Total FNs:", "{} ({}%)".format(total_fp, percent(filtered_fn, filtered_tp + filtered_fp + filtered_fn))])
        reports.append(["Sensitivity:", "{}".format(1.0 * filtered_tp / (filtered_tp + filtered_fp))])
        reports.append(["Precision:", "{}".format(1.0 * filtered_tp / (filtered_tp + filtered_fn))])

    # analysis on a per-chunk basis
    reports.append(["--------------------"])
    reports.append(["-  Chunk Summary   -"])
    reports.append(["--------------------"])
    def get_chunk_report(chunk):
        sensitivity = chunk[SENSITIVITY_KEY]
        precision = chunk[PRECISION_KEY]
        hmean = chunk[HMEAN_KEY]
        report = ["%8d %s%s%s" % (mpi,
                                  "+" if mpi in best_chunks else " ",
                                  "-" if mpi in worst_chunks else " ",
                                  "*" if mpi in specific_chunks.keys() else " ",),
                  "hmean: %.3f %s " % (hmean,  "+  " if hmean > .85 else ( " = " if hmean > .7 else "  -")),
                  "sensitivity: %.3f %s " % (sensitivity,
                                             "+  " if sensitivity > .8 else (
                                                 " = " if sensitivity > .5 else "  -")),
                  "precision: %.3f %s " % (precision,
                                           "+  " if precision > .85 else (" = " if precision > .7 else "  -")),
                  "total_calls: %5d  " % len(chunk[CALL_DETAILS_LIST_KEY]),
                  "avg_qual: %3d" % chunk[AVG_QUAL_KEY],
                  "avg_tp_qual: %3d" % chunk[AVG_TP_QUAL_KEY],
                  "avg_fp_qual: %3d" % chunk[AVG_FP_QUAL_KEY],
                  "range: %9d-%-9d (%7d)" % (chunk[START_POS_KEY], chunk[END_POS_KEY],
                                             chunk[END_POS_KEY] - chunk[START_POS_KEY])
                  ]
        return report
    # get reports
    for mpi in chunk_map_keys:
        chunk = chunk_map[mpi]
        if args.verbose or mpi in best_chunks or mpi in worst_chunks or mpi in specific_chunks.keys():
            reports.append(get_chunk_report(chunk))

    # specific chunks
    def get_specific_chunk_report(chunk):
        report = [
            "%8d  %s%s%s %s" % (chunk[MPI_KEY],
                              "+" if chunk[MPI_KEY] in chunk_ranking[chunk_ranking_third*2:] else " ",
                              "=" if chunk[MPI_KEY] in chunk_ranking[chunk_ranking_third:chunk_ranking_third*2] else " ",
                              "-" if chunk[MPI_KEY] in chunk_ranking[0:chunk_ranking_third] else " ",
                              "*" if chunk[MPI_KEY] in specific_chunks.keys() and args.chunks != '*' else " ")
        ]
        if len(chunk) > 1:
            report.extend( [
                "read count: %6d" % chunk[ORIGINAL_READ_COUNT_KEY],
                "filtered RC: %6d (%d%%)" % (chunk[FILTERED_READ_COUNT_KEY],
                                             percent(chunk[FILTERED_READ_COUNT_KEY],
                                                     chunk[ORIGINAL_READ_COUNT_KEY])),
                "depth avg: {}".format(int(chunk[AVG_READ_DEPTH_KEY])),
                "depth std: {}".format(int(chunk[STD_READ_DEPTH_KEY])),
                "depth max: {}".format(int(chunk[MAX_READ_DEPTH_KEY])),
                "depth min: {}".format(int(chunk[MIN_READ_DEPTH_KEY])),

                "length avg: {}".format(int(chunk[AVG_READ_LENGTH_KEY])),
                "length std: {}".format(int(chunk[STD_READ_LENGTH_KEY])),
                "length max: {}".format(int(chunk[MAX_READ_LENGTH_KEY])),
                "length min: {}".format(int(chunk[MIN_READ_LENGTH_KEY])),
                ])
        return report
    if len(specific_chunk_ids) > 0:
        reports.append(["--------------------"])
        reports.append(["- Specific  Chunk  -"])
        reports.append(["--------------------"])
        for mpi in specific_chunk_ids:
            chunk = specific_chunks[mpi]
            reports.append(get_specific_chunk_report(chunk))

    # print them
    print_reports(reports)

    # make plots
    if args.plot:
        print("Plotting")

        depth_count_below_threshold = lambda x, y: len(list(filter(lambda z: z < y, x[ALL_READ_DEPTH_KEY])))
        total_depth_below_threshold = lambda x, y: sum(list(map(lambda z: min(z, y), x[ALL_READ_DEPTH_KEY])))

        # generic
        # plot_all_call_qual(chunk_map, title="All Calls", logscale=True)
        # plot_sensitivity_x_precision(chunk_map)
        plot_quality_filtered_sensitivity_x_precision(chunk_map, save_name="Chunk_Filtered_SensPrec")

        # multi plot
        # plot_information = [
        #     plot_details(color='r', data=chunk_map, key=PRECISION_KEY, label="Precision"),
        #     plot_details(color='b', data=chunk_map, key=SENSITIVITY_KEY, label="Sensitivity"),
        #     plot_details(color='g', data=chunk_map, key=HMEAN_KEY, label="Mean")
        # ]
        # multi_plot_chunks_by_parameter(plot_information, specific_chunks, STD_READ_DEPTH_KEY, legend_loc=1,
        #                          title="Depth StdDev", y_label= "Depth StdDev", save_name="All3_vs_DepthStd")
        # multi_plot_chunks_by_parameter(plot_information, specific_chunks, MIN_READ_DEPTH_KEY, legend_loc=2,
        #                          title="Depth Minimum", y_label="Depth Min", save_name="All3_vs_DepthMin")
        # multi_plot_chunks_by_parameter(plot_information, specific_chunks, lambda x: depth_count_below_threshold(x, 5),
        #                          title="Depths below 5", legend_loc=3,
        #                          y_label="Depths below 5 (sampling 32 times)", save_name="All3_vs_DepthCntLt5")
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


        # # precision
        # plot_chunks_by_parameter(chunk_map, PRECISION_KEY, specific_chunks, STD_READ_DEPTH_KEY,
        #                          title="Precision vs Depth StdDev", save_name="Prec_vs_DepthStd")
        # plot_chunks_by_parameter(chunk_map, PRECISION_KEY, specific_chunks, MIN_READ_DEPTH_KEY,
        #                          title="Precision vs Min Depth", save_name="Prec_vs_DepthMin", color='r')
        # plot_chunks_by_parameter(chunk_map, PRECISION_KEY, specific_chunks, lambda x: depth_count_below_threshold(x, 5),
        #                          title="Precision vs Depths below 5", x_label="Precision",
        #                          y_label="Depths below 5 (sampling 32 times)", save_name="Prec_vs_DepthCntLt5", color='r')
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
        # # sensitivity
        # plot_chunks_by_parameter(chunk_map, SENSITIVITY_KEY, specific_chunks, STD_READ_DEPTH_KEY,
        #                          title="Sensitivity vs Depth StdDev", save_name="Sens_vs_DepthStd", color='b')
        # plot_chunks_by_parameter(chunk_map, SENSITIVITY_KEY, specific_chunks, MIN_READ_DEPTH_KEY,
        #                          title="Sensitivity vs Min Depth", save_name="Sens_vs_DepthMin", color='b')
        # plot_chunks_by_parameter(chunk_map, SENSITIVITY_KEY, specific_chunks, lambda x: depth_count_below_threshold(x, 5),
        #                          title="Sensitivity vs Depths below 5", x_label="Sensitivity",
        #                          y_label="Depths below 5 (sampling 32 times)", save_name="Sens_vs_DepthCntLt5", color='b')
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
        # # sensitivity and precision
        # plot_chunks_by_parameter(chunk_map, HMEAN_KEY, specific_chunks, STD_READ_DEPTH_KEY,
        #                          title="Sensitivity&Precision vs Depth StdDev", save_name="SnP_vs_DepthStd", color='g')
        # plot_chunks_by_parameter(chunk_map, HMEAN_KEY, specific_chunks, MIN_READ_DEPTH_KEY,
        #                          title="Sensitivity&Precision vs Min Depth", save_name="SnP_vs_DepthMin", color='g')
        # plot_chunks_by_parameter(chunk_map, HMEAN_KEY, specific_chunks, lambda x: depth_count_below_threshold(x, 5),
        #                          title="Sensitivity&Precision vs Depths below 5", x_label="Harmonic Mean",
        #                          y_label="Depths below 5 (sampling 32 times)", save_name="SnP_vs_DepthCntLt5", color='g')
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

        if len(specific_chunk_ids) > 0:
            for mpi in specific_chunk_ids:
                # plot_all_call_qual(chunk_map[mpi][CALL_DETAILS_LIST_KEY], title="Chunk {}".format(mpi), logscale=False)
                pass



    print("\nFin.")










if __name__ == "__main__":
    main()