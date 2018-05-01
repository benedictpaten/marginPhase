#!/usr/bin/env python
from __future__ import print_function
import argparse
import glob
import os
import subprocess
import sys



CHR = "c"
START = "s"
END = "e"
DESC = "d"

def parse_args():
    parser = argparse.ArgumentParser("Split BAM by region")
    parser.add_argument('--input_bam_glob', '-i', dest='input_bam_glob', required=True, type=str,
                       help='Glob matching input BAMs (will perform for all bams)')
    parser.add_argument('--coordinate_tsv', '-c', dest='coordinate_tsv', required=True, type=str,
                       help='Coordinates for splitting ($CHROM\\t$START\\t$END)')
    parser.add_argument('--output_location', '-o', dest='output_location', default=".", type=str,
                       help='Location where output files are put')
    parser.add_argument('--description_column', '-d', dest='description_column', default=None, type=int,
                       help='0-based index of description field in TSV (not required)')

    return parser.parse_args()


def get_output_filename(input_file_location, output_directory, coordinates):

    input_file_name = os.path.basename(input_file_location)
    input_file_parts = input_file_name.split(".")
    output_file_name = "{}.{}_{}-{}".format(".".join(input_file_parts[0:-1]), coordinates[CHR], coordinates[START],
                                            coordinates[END])
    if coordinates[DESC] is not None:
        output_file_name += "." + coordinates[DESC]
    output_file_name += "." + input_file_parts[-1]
    return os.path.join(output_directory, output_file_name)


def main():
    args = parse_args()
    assert False not in [len(args.input_bam_glob) > 0, os.path.isfile(args.coordinate_tsv), os.path.isdir(args.output_location)]

    coords = list()
    with open(args.coordinate_tsv) as tsv_in:
        header=True
        for line in tsv_in:
            if header:
                header = False
                continue
            line = line.split("\t")
            coords.append({
                CHR:line[0],
                START:int(line[1]),
                END:int(line[2]),
                DESC: None if args.description_column is None else "_".join(line[args.description_column].split())
            })

    for file in glob.glob(args.input_bam_glob):
        for coord in coords:
            outfile = get_output_filename(file, args.output_location, coord)
            print("{}:\n\tloc:  {}:{}-{}\n\tdesc: {}\n\tout:  {}".format(file, coord[CHR], coord[START], coord[END],
                                                                         coord[DESC], outfile), file=sys.stderr)
            samtools_args = ['samtools', 'view', '-hb', file, "{}:{}-{}".format(coord[CHR], coord[START], coord[END])]
            with open(outfile, 'w') as output:
                subprocess.check_call(samtools_args, stdout=output)

    print("Fin.", file=sys.stderr)






if __name__ == "__main__":
    main()