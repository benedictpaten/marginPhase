from __future__ import print_function
import argparse
import glob
import numpy as np
import os

def parse_args():
    parser = argparse.ArgumentParser("Summarize VCFEval summaries")
    parser.add_argument('--input', '-i', dest='input', default=None, required=True, type=str,
                       help='Glob matching summaries')
    parser.add_argument('--verbose', '-v', dest='verbose', default=False, action='store_true',
                       help='Print values')

    return parser.parse_args()


def main():
    args = parse_args()
    files = glob.glob(args.input)

    if len(files) == 0:
        print("Found {} files matching {}".format(len(files), args.input))
        return

    sensitivities = []
    precisions = []
    s_map = {}
    p_map = {}
    for file in files:
        with open(file, 'r') as input:
            for line in input:
                if line.startswith("-"): continue
                if line.startswith("Threshold"): continue
                line = line.split()
                precisions.append(float(line[5]))
                sensitivities.append(float(line[6]))
                p_map[os.path.basename(file)] = float(line[5])
                s_map[os.path.basename(file)] = float(line[6])

    print("{} ({} files):".format( args.input, len(files)))
    print("\tPrecision\t{} ({})".format(np.mean(precisions), np.std(precisions)))
    print("\tSensitivity\t{} ({})".format(np.mean(sensitivities), np.std(sensitivities)))
    if args.verbose:
        keys = p_map.keys()
        keys.sort()
        max_len = max(map(len, keys))
        for k in keys:
            print(("\t\t%" + str(max_len) + "s: \tP: %.3f \tS: %.3f") % (k, p_map[k], s_map[k]))


if __name__ == "__main__":
    main()