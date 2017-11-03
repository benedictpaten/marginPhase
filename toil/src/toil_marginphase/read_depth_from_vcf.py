from __future__ import print_function
import argparse
import glob
import os

def parse_args():
    parser = argparse.ArgumentParser("Computes read depth for output vcf glob from marginPhase run")
    parser.add_argument('--vcf', '-v', dest='vcf_glob', default="*.vcf", type=str,
                       help='The VCF(s) we\'re getting read depth from')

    return parser.parse_args()

def handle_vcf(vcf_file):
    depths = list()
    with open(vcf_file) as vcf:
        for line in vcf:
            if line.startswith("#"): continue
            line = line.split()
            tag_names = line[8].split(":")
            values = line[9].split(":")
            for i in range(len(tag_names)):
                if tag_names[i] == 'DP':
                    depths.append(int(values[i]))
                    break

    length = len(depths)
    total = 1.0*sum(depths)
    avg = total/length
    print("{}:".format(os.path.basename(vcf_file)))
    print("\tCALLS:     {}".format(length))
    print("\tAVG depth: {}".format(avg))

    return total, length

def main():
    args = parse_args()

    vcf_files = glob.glob(args.vcf_glob)
    vcf_files.sort()

    totals = list()
    lengths = list()
    for vcf_file in vcf_files:
        t, l = handle_vcf(vcf_file)
        totals.append(t)
        lengths.append(l)

    total_lengths = sum(lengths)
    print("\nOVERALL: ({})".format(args.vcf_glob))
    print("\tCALLS:     {}".format(total_lengths))
    print("\tAVG depth: {}".format(sum(totals)/total_lengths))
    print()


if __name__ == "__main__":
    main()