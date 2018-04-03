#!/usr/bin/env python
from __future__ import print_function
import argparse
import glob
import json
import sys
import numpy as np
import os

LEARNED_PARAM_KEY = "readErrorModel"
MIN_DIST_OFF = .00000001

def parse_args():
    parser = argparse.ArgumentParser("Consolidates output learned parameters for chunked marginPhase run")
    parser.add_argument('--params', '-p', dest='param_glob', default="*.json", type=str,
                       help='Matches for the params we\'re getting values for')
    parser.add_argument('--output', '-o', dest='output_file', default=None, type=str,
                       help='Output file. Default:stdout')
    parser.add_argument('--alphabet_size', '-s', dest='alphabet_size', default=5, type=int,
                       help='Size of alphabet for parameters')

    return parser.parse_args()


def log(msg, include_stdout_also = False):
    print(msg, file=sys.stderr)
    if include_stdout_also:
        print(msg, file=sys.stdout)


def print_params(param_list, alphabet_size, file=sys.stderr):
    file.write("    [\n")
    for i in range(alphabet_size):
        file.write(" " * 8)
        for j in range(alphabet_size):
            idx = i * alphabet_size + j
            file.write("    %.8f" % param_list[idx])
            if idx != alphabet_size * alphabet_size - 1: file.write(",")
        file.write("\n")
    file.write("    ]\n")


def main():
    args = parse_args()

    in_params = glob.glob(args.param_glob)
    if len(in_params) == 0:
        log("No files matching {}".format(args.param_glob))
        return 1
    else:
        log("Analyzing {} files".format(len(in_params)))

    param_size = args.alphabet_size ** 2
    all_param_values = [[] for _ in range(param_size)]

    missing_key = 0
    wrong_length = 0
    for in_param in in_params:
        with open(in_param) as input:
            param_str = ""
            for line in input:
                param_str += " ".join(line.strip().split())
            if param_str.endswith(",}"):
                param_str = param_str.replace(",}", " }")
            params = json.loads(param_str)
            if LEARNED_PARAM_KEY not in params:
                if missing_key == 0: log("Key '{}' not in {}".format(LEARNED_PARAM_KEY, in_param))
                missing_key += 1
                continue
            learned_param = params[LEARNED_PARAM_KEY]
            if len(learned_param) != param_size:
                if wrong_length == 0: log("Expected size {}x{}, got {} in {}".format(
                    args.alphabet_size,args.alphabet_size,len(learned_param), in_param))
                wrong_length += 1
                continue
            for i in range(param_size):
                all_param_values[i].append(learned_param[i])

    if missing_key != 0: log("{} files were missing keys".format(missing_key))
    if wrong_length != 0: log("{} files had the wrong length".format(wrong_length))
    if len(all_param_values[0]) == 0:
        log("No valid files matching {}".format(args.param_glob))
        return 1
    else:
        log("Consolidating params from {} files".format(len(all_param_values[0])))

    stddev_params = [np.std(x) for x in all_param_values]
    output_params = [np.mean(x) for x in all_param_values]
    for i in range(args.alphabet_size):
        row_total = sum(output_params[(i*args.alphabet_size):((i+1)*args.alphabet_size)])
        for j in range(args.alphabet_size):
            idx = i*args.alphabet_size + j
            output_params[idx] = output_params[idx] / row_total
        new_row_total = sum(output_params[(i*args.alphabet_size):((i+1)*args.alphabet_size)])
        assert abs(new_row_total-1) < MIN_DIST_OFF
    assert abs(sum(output_params)- args.alphabet_size) < MIN_DIST_OFF

    log("\nSTD DEV (beware if these are large):")
    print_params(stddev_params, args.alphabet_size)

    log("\nOutput params:")
    print_params(output_params, args.alphabet_size, file=(sys.stdout if args.output_file is None else sys.stderr))

    if args.output_file is not None:
        with open(args.output_file, 'w') as output:
            print_params(output_params, args.alphabet_size, file=output)

    log("\nFin.")









if __name__ == "__main__":
    main()