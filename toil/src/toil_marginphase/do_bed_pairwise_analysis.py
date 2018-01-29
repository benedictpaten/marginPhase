#!/usr/bin/env python

from __future__ import print_function
import subprocess
import os

BEDS = [
    # 'np_pb_intersect.bed',
    # 'np_pb_int_gte10.bed',
    'np_pb_int_gte5.bed',
    # 'np_pb_int_gte15.bed',
    # 'gatk_gvcf_callable.bed',
    'cloci_call_exces.bed',
    # 'giab_high_conf.bed'
]

GENOME = 'hg38.genome'


THAT_HACKY_INTERSECT_SCRIPT_NAME = "do_that_hacky_intersect.sh"
def write_that_hacky_intersect_script(filename=THAT_HACKY_INTERSECT_SCRIPT_NAME):
    with open(filename, 'w') as output:
        print("#!/bin/bash", file=output)
        print("bedtools intersect -a $1 -b $2 | bedtools sort | bedtools merge >$3", file=output)
    # os.chmod(filename, 0777)
    assert os.path.isfile(filename)


def get_intersect_name(bed1, bed2):
    bed1 = bed1.rstrip(".bed")
    bed2 = bed2.rstrip(".bed")
    name_parts = ("{}.{}".format(bed1, bed2)).split(".")
    name_parts.sort()
    return "{}.bed".format(".".join(name_parts))


def create_intersect_bed(bed1, bed2, output):
    if os.path.isfile(output): return
    print("Creating: {}".format(output))
    cmd = ['bash', THAT_HACKY_INTERSECT_SCRIPT_NAME, bed1, bed2, output]
    subprocess.check_call(cmd)
    assert os.path.isfile(output)


THAT_HACKY_COVERAGE_SCRIPT_NAME = "do_that_hacky_coverage.sh"
def write_that_hacky_coverage_script(filename=THAT_HACKY_COVERAGE_SCRIPT_NAME):
    with open(filename, 'w') as output:
        print("#!/bin/bash", file=output)
        print("bedtools genomecov -g {} -i $1 >$2".format(GENOME), file=output)
    # os.chmod(filename, 0777)
    assert os.path.isfile(filename)


def get_genome_coverage(bed, silent=False):
    coverage_file = "{}.coverage".format(bed.rstrip(".bed"))
    if not silent: print("{}:".format(coverage_file))
    if not os.path.isfile(coverage_file):
        cmd = ['bash', THAT_HACKY_COVERAGE_SCRIPT_NAME, bed, coverage_file]
        subprocess.check_call(cmd)
    assert os.path.isfile(coverage_file)

    included = None
    with open(coverage_file, 'r') as genomecov:
        for line in genomecov:
            if 'genome' not in line: continue
            line = line.split()
            assert len(line) == 5
            assert line[1] in ["0", "1"]
            if line[1] == '1':
                included = float(line[4])
                if not silent: print("\t%.5f" % included)
    assert included is not None
    return included


def main():
    write_that_hacky_intersect_script()
    write_that_hacky_coverage_script()

    # do first pass
    prev_run = []
    print("\nDEPTH: 0")
    for file in BEDS:
        get_genome_coverage(file)
        prev_run.append(file)

    # do next passes
    depth = 1
    while depth < len(BEDS):
        print("\nDEPTH: {}".format(depth))
        current_run = set()
        for file1 in BEDS:
            for file2 in prev_run:
                parts = (file2.rstrip(".bed")).split(".")
                if file1.rstrip(".bed") in parts: continue
                intersect = get_intersect_name(file1, file2)
                if intersect in current_run: continue
                create_intersect_bed(file1, file2, intersect)
                get_genome_coverage(intersect)
                current_run.add(intersect)
        prev_run = list(current_run)
        depth += 1

    bit_to_bed = dict()
    bed_to_bit = dict()
    all_bits = list()
    bit = 1
    for file in BEDS:
        bed_to_bit[file] = bit
        bit_to_bed[bit] = file
        all_bits.append(bit)
        bit *= 2
    all_bits.reverse()

    all_coverages = dict()
    bitstring_map = dict()
    bit_ids = dict()
    idx = 1
    while idx < 2 ** len(BEDS):
        current_files = list()
        bitmap = idx
        current_bitstring = ""
        for bit in all_bits:
            current_bitstring += "1" if bitmap & bit else "0"
            if bitmap & bit > 0:
                current_files.append(bit_to_bed[bit].rstrip('.bed'))
                bitmap -= bit
        if bitmap != 0:
            raise Exception("Bad logic for idx {}".format(idx))
        assert len(current_files) > 0

        if len(current_files) == 1:
            bit_ids[current_bitstring] = current_files[0]
        bitstring_map[current_bitstring] = idx
        current_files.sort()
        intersection_bed = "{}.bed".format(".".join(current_files))
        all_coverages[idx] = get_genome_coverage(intersection_bed, silent=True)

        idx += 1

    bitstrings = list(bit_ids.keys())
    bitstrings.sort()
    print("\nBED Keys:")
    for bs in bitstrings:
        print("{}\t{}".format(bs, bit_ids[bs]))

    bitstrings = list(bitstring_map.keys())
    bitstrings.sort()
    print("\nBED Intersections:")
    for bs in bitstrings:
        print("{}\t{}".format(bs, all_coverages[bitstring_map[bs]]))







if __name__ == "__main__":
    main()