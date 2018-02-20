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

VCF_1 = 1
VCF_2 = 2

VCF_LINE = "l"
VCF_POS = "p"
VCF_CHR = "c"

CHR_IDX = 0
POS_IDX = 1
ID_IDX = 2
REF_IDX = 3
ALT_IDX = 4
QUAL_IDX = 5
FILTER_IDX = 6
INFO_IDX = 7
FORMAT_IDX = 8
SAMPLE_IDX = 9

INFO_SEP = ";"
INFO_ELEM_SEP = "="
FORMAT_SEP = ":"

chrom_sort = lambda x: int(x.replace("chr", "").replace("X", "23").replace("Y", "24"))
percent = lambda part, whole: int(100.0 * part / whole) if whole != 0 else "--"
EQUIVALENT_FILTERS = [
    ["PASS", "."]
]

def parse_args():
    parser = argparse.ArgumentParser("Merges two VCFs")
    parser.add_argument('--vcf_1', '-1', dest='vcf_1', type=str, required=True,
                       help='First VCF file for merging')
    parser.add_argument('--vcf_2', '-2', dest='vcf_2', type=str, required=True,
                       help='Second VCF file for merging')
    parser.add_argument('--output', '-o', dest='output', default="-", type=str,
                       help="If set, output will be written to this file, otherwise stdout. ")

    #tags
    parser.add_argument('--tag_prefix_1', '-t1', dest='tag_prefix_1', default="", type=str,
                       help='Prepend this parameter to the tags in \'vcf_1\'')
    parser.add_argument('--tag_prefix_2', '-t2', dest='tag_prefix_2', default="", type=str,
                       help='Prepend this parameter to the tags in \'vcf_2\'')
    parser.add_argument('--tag_suffix_1', '-s1', dest='tag_suffix_1', default="1", type=str,
                       help='Append this parameter to the tags in \'vcf_1\'')
    parser.add_argument('--tag_suffix_2', '-s2', dest='tag_suffix_2', default="2", type=str,
                       help='Append this parameter to the tags in \'vcf_2\'')

    # misc
    parser.add_argument('--sample_name', '-s', dest='sample_name', default=None, type=str,
                       help='New name to use for sample column')
    parser.add_argument('--verbose', '-v', dest='verbose', action='store_true', default=False,
                       help='Print extra information')


    # parser.add_argument('--quiet', '-V', dest='quiet', action='store_true', default=False,
    #                    help='Print less information')
    # parser.add_argument('--chunks', '-c', dest='chunks', default=None, type=str,
    #                    help='Specify comma-separated MPIs of specific chunks to analyze')
    # parser.add_argument('--plot', '-p', dest='plot', action='store_true', default=False,
    #                    help='Make plots of data')
    # parser.add_argument('--chromosomes', '-o', dest='chromosomes', default=None, type=str,
    #                    help='Limit calls in VCF and BED files to specified comma-separated chromosomes')
    #
    # parser.add_argument('--qual_min', '-q', dest='qual_min', default=None, type=int,
    #                    help='FILTER: Filter VCF qualities below threshold: TP->FN, FP->/dev/null')
    # parser.add_argument('--vcf_read_depth_min', '-d', dest='vcf_read_depth_min', default=None, type=int,
    #                    help='FILTER: Filter VCF read depth below threshold: TP->FN, FP->/dev/null')
    # parser.add_argument('--bam_read_depth_min', '-D', dest='bam_read_depth_min', default=None, type=int,
    #                    help='FILTER: Filter calls where read depth below threshold in bam. '
    #                         'Used in conjuction with --chunk_bam_format option. TP,FN,FP -> /dev/null')
    # parser.add_argument('--remove_chunks', '-C', dest='remove_chunks', default=None, type=str,
    #                    help='FILTER: Remove comma-separated MPIs (or ranges) of chunks from analysis. ex: \'0-4,22\'')
    # parser.add_argument('--inclusion_beds', '-e', dest='inclusion_beds', default=None, type=str,
    #                    help='FILTER: only include calls from within comma-separated bed files')
    # parser.add_argument('--exclusion_beds', '-E', dest='exclusion_beds', default=None, type=str,
    #                    help='FILTER: exclude calls from within comma-separated bed files')
    #
    # parser.add_argument('--genome_bam_glob', '-b', dest='genome_bam_glob', default=None, type=str,
    #                    help='Filename glob for input BAMs (will output depth statistics)')
    # parser.add_argument('--pickle_filename', '-s', dest='pickle_filename', default=None, type=str,
    #                     help="Filename for pickling or unpickling read depth data.  Depth info will be unpickled from "
    #                          "here if the file exists, otherwise data in 'genome_bam_glob' parameter will be saved.")

    return parser.parse_args()


def read_file(header, calls, vcf_id, input):
    # prep
    header[vcf_id] = list()
    header = header[vcf_id]

    linenr = -1
    for line in input:
        # prep
        linenr += 1
        if line.startswith("#"):
            header.append(line.rstrip())
            continue

        # get line info
        line_parts = line.strip().split("\t")
        if len(line_parts) < 10:
            raise Exception("Line {} is malformed:\n\t[{}]".format(linenr, '\t'.join(line)))
        chrom = line_parts[0]
        pos = int(line_parts[1])

        # get chrom map
        if chrom not in calls:
            chrom_map = dict()
            calls[chrom] = chrom_map
        else:
            chrom_map = calls[chrom]

        # get pos map
        if pos not in chrom_map:
            pos_map = dict()
            chrom_map[pos] = pos_map
        else:
            pos_map = chrom_map[pos]

        #sanity and save
        if vcf_id in pos_map: raise Exception("VCF {} already found at {} {}!\n\tNew Line: {}\n\tOld Line: {}"
                                              .format(vcf_id, chrom, pos, line, pos_map[vcf_id]))
        pos_map[vcf_id] = line

    return linenr


def merge_header(args, headers):
    misc_lines =  list()
    command_lines = list()
    command_lines.append("CL=" + " ".join(sys.argv))
    contig_lines = list()
    filter_lines = list()
    info_lines = list()
    format_lines = list()
    sample_name = args.sample_name

    # how to place header lines into approprate places
    def organize_header_line(args, orig_line, is_vcf1):
        #clean lines
        line = orig_line.lstrip("#")
        #command lines
        if line.upper().startswith("CL"):
            command_lines.append(line)
        # contig
        elif line.upper().startswith("CONTIG"):
            if line not in contig_lines: contig_lines.append(line)
        # filter
        elif line.upper().startswith("FILTER"):
            if line not in filter_lines: filter_lines.append(line)
        # info needs to be modified with prefix and suffix
        elif line.upper().startswith("INFO"):
            line = line.lstrip("INFO").lstrip("info").lstrip("=").lstrip("<").rstrip(">")
            line_parts = line.split(",")
            for i in range(len(line_parts)):
                if line_parts[i].startswith("ID"):
                    id_parts = line_parts[i].split("=")
                    id_parts[1] = "{}{}{}".format(args.tag_prefix_1 if is_vcf1 else args.tag_prefix_2, id_parts[1],
                                                  args.tag_suffix_1 if is_vcf1 else args.tag_suffix_2)
                    line_parts[i] = "=".join(id_parts)
            line = ",".join(line_parts)
            info_lines.append("INFO=<{}>".format(line))
        # format needs to be modified with prefix and suffix
        elif line.upper().startswith("FORMAT"):
            line = line.lstrip("FORMAT").lstrip("format").lstrip("=").lstrip("<").rstrip(">")
            line_parts = line.split(",")
            for i in range(len(line_parts)):
                if line_parts[i].startswith("ID"):
                    id_parts = line_parts[i].split("=")
                    id_parts[1] = "{}{}{}".format(args.tag_prefix_1 if is_vcf1 else args.tag_prefix_2, id_parts[1],
                                                  args.tag_suffix_1 if is_vcf1 else args.tag_suffix_2)
                    line_parts[i] = "=".join(id_parts)
            line = ",".join(line_parts)
            format_lines.append("FORMAT=<{}>".format(line))
        # everything else is misc
        else:
            if line not in misc_lines:
                misc_lines.append(line)

    # for each of the vcf files, get all the lines (vcf_1 first) and handle them
    for key in [VCF_1, VCF_2]:
        for line in headers[key]:
            # header lines
            if line.startswith("##"):
                organize_header_line(args, line, VCF_1 == key)
            # describes the columns
            elif line.startswith("#"):
                line = line.split("\t")
                if len(line) <= SAMPLE_IDX:
                    raise Exception("Malformed header in VCF {}: {}".format(1 if VCF_1 == key else 2, "\t".join(line)))
                if sample_name is None:
                    sample_name = line[SAMPLE_IDX]
            # should never happen
            else:
                raise Exception("PROGRAMMER ERROR: header managment improperly configured")

    # build output header with lines appropriately ordered
    header_out = list()
    header_out.extend(misc_lines)
    header_out.extend(command_lines)
    header_out.extend(contig_lines)
    header_out.extend(filter_lines)
    header_out.extend(info_lines)
    header_out.extend(format_lines)
    header_out = ["##{}".format(x) for x in header_out]
    # add the last column line
    header_out.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}".format(sample_name))

    # return
    return header_out


def merge_lines(args, line1, line2):
    line1 = line1.strip().split('\t')
    line2 = line2.strip().split('\t')

    assert line1[CHR_IDX] == line2[CHR_IDX]
    assert line1[POS_IDX] == line2[POS_IDX]
    differences = list()
    reverse_vcf2_phase = False
    if line1[REF_IDX] != line2[REF_IDX]:
        differences.append("REF")
    if line1[ALT_IDX] != line2[ALT_IDX]:
        line1_alts = set(line1[ALT_IDX].split(","))
        line2_alts = set(line2[ALT_IDX].split(","))
        if len(line1_alts) != len(line2_alts) or len(list(filter(lambda x: x not in line2_alts, line1_alts))) > 0:
            differences.append("ALT")
        else:
            reverse_vcf2_phase = True
    if line1[FILTER_IDX] != line2[FILTER_IDX]:
        equiv = False
        for equiv_filter in EQUIVALENT_FILTERS:
            if line1[FILTER_IDX] in equiv_filter and line2[FILTER_IDX] in equiv_filter:
                equiv = True
                break
        if not equiv:
            differences.append("FILTER")
    if len(differences) != 0:
        return False, "Different {}\n\tVCF1:\t{}\n\tVCF2:\t{}".format(",".join(differences), line1, line2)

    # get easy elements
    chr = line1[CHR_IDX]
    pos = line1[POS_IDX]
    id = line1[ID_IDX] if line1[ID_IDX] == line2[ID_IDX] else '.'
    ref = line1[REF_IDX]
    alt = line1[ALT_IDX]
    qual = str(min([int(line1[QUAL_IDX]),int(line2[QUAL_IDX])]))
    filtr = line1[FILTER_IDX]

    # get more complicated elements
    info = list()
    for info_elem in line1[INFO_IDX].split(INFO_SEP):
        tag, val = info_elem.split(INFO_ELEM_SEP)
        new_tag = "{}{}{}".format(args.tag_prefix_1, tag, args.tag_suffix_1)
        info.append(INFO_ELEM_SEP.join([new_tag,val]))
    for info_elem in line2[INFO_IDX].split(INFO_SEP):
        tag, val = info_elem.split(INFO_ELEM_SEP)
        new_tag = "{}{}{}".format(args.tag_prefix_2, tag, args.tag_suffix_2)
        info.append(INFO_ELEM_SEP.join([new_tag,val]))
    info = INFO_SEP.join(info)

    format = list()
    vcf_phase_idx = None
    for format_tag in line1[FORMAT_IDX].split(FORMAT_SEP):
        format.append("{}{}{}".format(args.tag_prefix_1, format_tag, args.tag_suffix_1))
    for format_tag in line2[FORMAT_IDX].split(FORMAT_SEP):
        if reverse_vcf2_phase and format_tag == "GT":
            vcf_phase_idx = len(format)
        format.append("{}{}{}".format(args.tag_prefix_2, format_tag, args.tag_suffix_2))
    format_len = len(format)
    format = FORMAT_SEP.join(format)

    sample = list()
    for sample_val in line1[SAMPLE_IDX].split(FORMAT_SEP):
        sample.append(sample_val)
    for sample_val in line2[SAMPLE_IDX].split(FORMAT_SEP):
        if vcf_phase_idx is not None and len(sample) == vcf_phase_idx:
            if len(sample_val) != 3:
                return False, "Malformed GT value:\n\tVCF1:\t{}\n\tVCF2:\t{}".format(line1, line2)
            sample_val = "{}{}{}".format(sample_val[2],sample_val[1],sample_val[0])
        sample.append(sample_val)
    sample_len = len(sample)
    sample = FORMAT_SEP.join(sample)

    # sanity check
    if format_len != sample_len:
        return False, "FORMAT and SAMPLE have different lengths:\n\tVCF1:\t{}\n\tVCF2:\t{}".format(line1, line2)

    # get new line
    new_line = "\t".join([chr, pos, id, ref, alt, qual, filtr, info, format, sample])
    return True, new_line

def main():
    # prep
    overall_start = time.clock()
    args = parse_args()
    header = dict()
    calls = dict()

    # get files and sanity check
    print("Reading VCF files", file=sys.stderr)
    start = time.clock()
    vcf1_file = args.vcf_1
    vcf2_file = args.vcf_2
    if not os.path.isfile(vcf1_file):
        raise Exception("Could not find file {}".format(vcf1_file))
    if not os.path.isfile(vcf2_file):
        raise Exception("Could not find file {}".format(vcf2_file))

    # read files
    with (gzip.open if vcf1_file.endswith("gz") else open)(vcf1_file, 'r') as fp_in:
        line_cnt = read_file(header, calls, VCF_1, fp_in)
        print("\t{}: {} lines".format(vcf1_file, line_cnt), file=sys.stderr)
    with (gzip.open if vcf2_file.endswith("gz") else open)(vcf2_file, 'r') as fp_in:
        line_cnt = read_file(header, calls, VCF_2, fp_in)
        print("\t{}: {} lines".format(vcf2_file, line_cnt), file=sys.stderr)
    print("\t{}s".format(int(time.clock() - start)), file=sys.stderr)

    # header
    print("Handling header", file=sys.stderr)
    start = time.clock()
    header_out = merge_header(args, header)
    print("\tmerged {} header lines into {}".format(sum([len(header[x]) for x in header.keys()]), len(header_out)),
          file=sys.stderr)
    print("\t{}s".format(int(time.clock() - start)), file=sys.stderr)

    # merge calls
    print("Merging calls", file=sys.stderr)
    start = time.clock()
    calls_in_vcf1 = list()
    calls_in_vcf2 = list()
    call_count_in_both = 0
    errors = list()
    calls_out = list()
    chroms = calls.keys()
    chroms.sort(key=chrom_sort)
    # get all calls
    for chrom in chroms:
        chrom_calls = calls[chrom]
        positions = list(chrom_calls.keys())
        positions.sort(key=lambda x: int(x))
        # get all positions
        for pos in positions:
            pos_calls = chrom_calls[pos]
            if len(pos_calls) == 1:
                # missing one call
                if VCF_1 in pos_calls:
                    calls_in_vcf1.append(pos_calls[VCF_1])
                elif VCF_2 in pos_calls:
                    calls_in_vcf2.append(pos_calls[VCF_2])
                else: assert False
            elif len(pos_calls) == 2:
                # have both calls
                call_count_in_both += 1
                success, value = merge_lines(args, pos_calls[VCF_1], pos_calls[VCF_2])
                if success:
                    calls_out.append(value)
                else:
                    errors.append(value)
    call_count_in_vcf1 = len(calls_in_vcf1)
    call_count_in_vcf2 = len(calls_in_vcf2)
    call_count_total = sum([call_count_in_vcf1, call_count_in_vcf2, call_count_in_both])
    print("\ttotal calls:        {}".format(call_count_total), file=sys.stderr)
    print("\tcalls only in VCF1: {} ({}%)".format(call_count_in_vcf1, percent(call_count_in_vcf1, call_count_total)), file=sys.stderr)
    if args.verbose:
        for call in calls_in_vcf1:
            print("\t\t{}".format(call.strip()), file=sys.stderr)
    print("\tcalls only in VCF2: {} ({}%)".format(call_count_in_vcf2, percent(call_count_in_vcf2, call_count_total)), file=sys.stderr)
    if args.verbose:
        for call in calls_in_vcf2:
            print("\t\t{}".format(call.strip()), file=sys.stderr)
    print("\tcalls in both:      {} ({}%)".format(call_count_in_both, percent(call_count_in_both, call_count_total)), file=sys.stderr)
    print("\t    successes:      {} ({}%)".format(len(calls_out), percent(len(calls_out), call_count_in_both)), file=sys.stderr)
    print("\t    failures:       {} ({}%)".format(len(errors), percent(len(errors), call_count_in_both)), file=sys.stderr)
    if args.verbose:
        for err in errors:
            print("\t\t{}".format(err.replace("\n", "\n\t\t")), file=sys.stderr)
    print("\t{}s".format(int(time.clock() - start)), file=sys.stderr)

    print("Writing output")
    start = time.clock()
    lines_out = 0
    with (open(args.output, 'w') if args.output != "-" else sys.stdout) as out:
        for h in header_out:
            print(h, file=out)
            lines_out += 1
        for c in calls_out:
            print(c, file=out)
            lines_out += 1
        print("\t{}: wrote {} lines".format(args.output if args.output != "-" else "STDOUT", lines_out),file=sys.stderr)
    print("\t{}s".format(int(time.clock() - start)))

    print("\nFin ({}s).".format(int(time.clock() - overall_start)),file=sys.stderr)


if __name__ == "__main__":
    main()