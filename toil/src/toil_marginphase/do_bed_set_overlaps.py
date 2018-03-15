#!/usr/bin/env python

from __future__ import print_function
import subprocess
import os
import shutil
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})


COVERAGE_BEDS = [
    'giab_high_conf.bed',
    'gatk_gvcf_callable.bed',
    # 'int_gte10_lt2xm.bed',
    'int_gte5_lt2xm.bed',
    # 'pb_gte20.bed',
    # 'np_gte20.bed',
    # 'pb_gte15.bed',
    # 'np_gte15.bed',
    # 'pb_gte10.bed',
    # 'np_gte10.bed',
    'pb_gte5_lt92.bed',
    'np_gte5_lt74.bed',
    # 'pb_gte5.bed',
    # 'np_gte5.bed',
]

FEATURE_BEDS = [
    'hg38_alu.bed',
    'hg38_gencode24.bed',
    'hg38_gencode27.bed',
    'hg38_tandemRepeats.bed',
    'hg38_repeatMasker.bed',
    'hg38_segmentalDups.bed',
]

GENOME = 'genomic_features/hg38.genome'

BASELINE_COLOR = 'DarkSlateGray'
LONG_READ_COLOR = 'MediumPurple'
SHORT_READ_COLOR = 'IndianRed'

IDENTIFIER_KEY='id'
COLOR_KEY='color'
LABEL_KEY='label'
VALUE_KEY='value'
BED_LOC_KEY='bed_loc'

LABEL_IDX=0
COLOR_IDX=1

BED_INFO = {
    'giab_high_conf.bed': ['GIAB High Confidence', SHORT_READ_COLOR],
    'gatk_gvcf_callable.bed': ['GATK GVCF Callable', SHORT_READ_COLOR],
    'int_gte10_lt2xm.bed': ['NP/PB Intersection (Depth 10)', LONG_READ_COLOR],
    'int_gte5_lt2xm.bed': ['NP/PB Intersection', LONG_READ_COLOR],
    'pb_gte20.bed': ['PacBio (Depth 20)', LONG_READ_COLOR],
    'np_gte20.bed': ['Nanopore (Depth 20)', LONG_READ_COLOR],
    # 'pb_gte15.bed': ['', LONG_READ_COLOR],
    # 'np_gte15.bed': ['', LONG_READ_COLOR],
    # 'pb_gte10.bed': ['', LONG_READ_COLOR],
    # 'np_gte10.bed': ['', LONG_READ_COLOR],
    'pb_gte5_lt92.bed': ['PacBio', LONG_READ_COLOR],
    'np_gte5_lt74.bed': ['Nanopore', LONG_READ_COLOR],
    'pb_gte5.bed': ['PacBio (Depth 5)', LONG_READ_COLOR],
    'np_gte5.bed': ['Nanopore (Depth 5)', LONG_READ_COLOR],
    'hg38_alu.bed': ['ALUs', BASELINE_COLOR],
    'hg38_gencode24.bed': ['Gencode v24', BASELINE_COLOR],
    'hg38_repeatMasker.bed': ['Repeat Masker', BASELINE_COLOR],
    'hg38_segmentalDups.bed': ['Segmental Dups', BASELINE_COLOR],
    'hg38_gencode27.bed': ['Gencode v27', BASELINE_COLOR],
    'hg38_tandemRepeats.bed': ['Tandem Repeats', BASELINE_COLOR],
}


COVERAGE_BED_DIR = "genomic_features/coverage_beds"
FEATURE_BED_DIR = "genomic_features"
OUTPUT_BED_DIR = "genomic_features/coverage_output"
PLOT_BED_DIR = "genomic_features/plots"


def get_bed_intersect_location(coverage, feature):
    return os.path.join(OUTPUT_BED_DIR, "{}.{}.bed".format(coverage.rstrip(".bed"), feature.rstrip(".bed")))
def get_coverage_bed_location(coverage):
    return os.path.join(OUTPUT_BED_DIR, coverage)
def get_feature_bed_location(feature):
    return os.path.join(OUTPUT_BED_DIR, feature)


THAT_HACKY_INTERSECT_SCRIPT_NAME = "do_that_hacky_intersect.sh"
def write_that_hacky_intersect_script(filename=THAT_HACKY_INTERSECT_SCRIPT_NAME):
    with open(filename, 'w') as output:
        print("#!/bin/bash", file=output)
        print("bedtools intersect -a $1 -b $2 | bedtools sort | bedtools merge >$3", file=output)
    # os.chmod(filename, 0777)
    assert os.path.isfile(filename)


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


def get_genome_coverage(bed):
    coverage_file = "{}.coverage".format(bed.rstrip(".bed"))
    if not os.path.isfile(coverage_file):
        print("Creating: {}".format(coverage_file))
        cmd = ['bash', THAT_HACKY_COVERAGE_SCRIPT_NAME, bed, coverage_file]
        subprocess.check_call(cmd)
    assert os.path.isfile(coverage_file)

    included = None
    with open(coverage_file, 'r') as genomecov:
        for line in genomecov:
            if 'genome' not in line: continue
            line = line.split()
            assert len(line) == 5
            if int(line[1]) > 0:
                if included is None: included = 0.0
                included += float(line[4])
    assert included is not None
    return included


def format_plot_element(bed_identifier, label, color, bed_location):
    return {
        IDENTIFIER_KEY: bed_identifier,
        COLOR_KEY: color,
        LABEL_KEY: label,
        BED_LOC_KEY: bed_location,
        VALUE_KEY: get_genome_coverage(bed_location) if bed_location is not None else 1.0
    }


def plot_feature_coverage(coverage_beds, feature_bed, show=False, save=True):

    # get elements we plot
    elements = list()
    # add our baseline
    if feature_bed is None:
        elements.append(format_plot_element('genome', "Genome", BASELINE_COLOR, None))
    else:
        elements.append(format_plot_element(feature_bed, BED_INFO[feature_bed][LABEL_IDX], BASELINE_COLOR,
                                            get_feature_bed_location(feature_bed)))
    # add all our coverage beds
    for bed in coverage_beds:
        elements.append(format_plot_element(bed, BED_INFO[bed][LABEL_IDX], BED_INFO[bed][COLOR_IDX],
                            get_bed_intersect_location(bed, feature_bed) if feature_bed is not None
                            else get_coverage_bed_location(bed)))
    #sort
    elements.sort(key=lambda x: x[VALUE_KEY])

    xs = []
    ys = []
    for i,v in enumerate(elements):
        xs.append(i)
        ys.append(v[VALUE_KEY])

    fig, ax = plt.subplots(figsize=(12, 3))
    ax.barh(xs, ys, color=list(map(lambda x: x[COLOR_KEY], elements)))
    ax.set_yticks(xs)
    ax.set_yticklabels(list(map(lambda x: x[LABEL_KEY], elements)))
    for x, v in enumerate(ys):
        ax.text(max(ys) * 1.005, x-.1, ("%.3f" % v), color=list(map(lambda x: x[COLOR_KEY], elements))[x])

    plt.title("{} Coverage".format("Genome" if feature_bed is None else BED_INFO[feature_bed][LABEL_IDX]))
    if show:
        plt.show()
    if save:
        fig.savefig(os.path.join(PLOT_BED_DIR, "{}.png".format(feature_bed if feature_bed is not None else 'genome')), dpi=300)


def main():
    write_that_hacky_intersect_script()
    write_that_hacky_coverage_script()

    for coverage in COVERAGE_BEDS:
        if not os.path.isfile(get_coverage_bed_location(coverage)):
            shutil.copy(os.path.join(COVERAGE_BED_DIR, coverage), get_coverage_bed_location(coverage))
        get_genome_coverage(get_coverage_bed_location(coverage))
    plot_feature_coverage(COVERAGE_BEDS, None)

    spacing = str(max(max(map(len, FEATURE_BEDS)), 4 + max(map(len, COVERAGE_BEDS))))
    for feature in FEATURE_BEDS:
        print("\n#############################################")
        if not os.path.isfile(get_feature_bed_location(feature)):
            shutil.copy(os.path.join(FEATURE_BED_DIR, feature), get_feature_bed_location(feature))
        their_cov = get_genome_coverage(os.path.join(OUTPUT_BED_DIR, feature))
        print(('%-'+spacing+'s:  %.5f') % (feature, their_cov))
        for coverage in COVERAGE_BEDS:
            intersect = get_bed_intersect_location(coverage, feature)
            create_intersect_bed(get_coverage_bed_location(coverage), get_feature_bed_location(feature), intersect)
            overlap = get_genome_coverage(intersect)
            print(('%-'+spacing+'s:  %.5f    %2.1f%%') % ("   "+coverage, overlap, 100.0 * overlap / their_cov))

        plot_feature_coverage(COVERAGE_BEDS, feature)



if __name__ == "__main__":
    main()