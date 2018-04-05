#!/usr/bin/env python

from __future__ import print_function
import subprocess
import os
import shutil
import math
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rcParams
rcParams.update({'font.size': 16})
rcParams.update({'figure.autolayout': True})


COVERAGE_BEDS = [
    # 'int_gte5_lt2xm.bed',
    # 'pb_gte20.bed',
    # 'np_gte20.bed',
    # 'pb_gte15.bed',
    # 'np_gte15.bed',
    # 'pb_gte10.bed',
    # 'np_gte10.bed',
    # 'pb_gte5.bed',
    # 'np_gte5.bed',
    # 'pb_noreadfilt_gte5.bed',
    # 'np_noreadfilt_gte5.bed',
    'uni_gte1.bed',
    # 'uni_gte5.bed',
    'pb_gte1.bed',
    'np_gte1.bed',
    # 'pb_gte5_lt92.bed',
    # 'np_gte5_lt74.bed',
    'int_gte10_lt2xm.bed',
    'gatk_cloci_callable.bed',
    'gatk_gvcf_callable.bed',
    'giab_high_conf.bed',
]
COVERAGE_BEDS.reverse()

FEATURE_BEDS = [
    # 'hg38_alu.bed',
    # 'hg38_gencode24.bed',
    # 'hg38_gencode27.bed',
    'hg38_gencode27_basic.bed',
    'hg38_gencode27_exome.bed',
    'hg38_repeatMasker.bed',
    # 'Satellite.bed',
    # 'SimpleRepeat.bed',
    # 'LINE.bed',
    # 'SINE.bed',
    'XINE.bed',
    # 'Alu.bed',
    'hg38_segmentalDups.bed',
    'hg38_tandemRepeats.bed',
    'Telomere.bed',
    'Centromere.bed',
    # 'genomic_regions_definitions_modeledcentromere.bed',
    # 'INC_gencode.INC_repeatmask.bed',
]



GENOME = 'genomic_features/hg38.genome'

BASELINE_COLOR = 'DarkSlateGray'
BASELINE_COLOR = '#393939'
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
    'int_gte5_lt2xm.bed': ['NP/PB Intersection (Depth 5)', LONG_READ_COLOR],
    'pb_gte5.bed': ['PacBio (Depth 5)', LONG_READ_COLOR],
    'np_gte5.bed': ['Nanopore (Depth 5)', LONG_READ_COLOR],
    'pb_gte20.bed': ['PacBio (Depth 20)', LONG_READ_COLOR],
    'np_gte20.bed': ['Nanopore (Depth 20)', LONG_READ_COLOR],
    # 'pb_gte15.bed': ['', LONG_READ_COLOR],
    # 'np_gte15.bed': ['', LONG_READ_COLOR],
    # 'pb_gte10.bed': ['', LONG_READ_COLOR],
    # 'np_gte10.bed': ['', LONG_READ_COLOR],
    # 'pb_noreadfilt_gte5.bed': ['PacBio (All Reads)', '#C00000'],
    # 'np_noreadfilt_gte5.bed': ['Nanopore (All Reads)', '#C00060'],

    # 'uni_gte1.bed': ['NP/PB Union (Depth 1)', LONG_READ_COLOR],
    # 'uni_gte5.bed': ['NP/PB Union (Depth 5)', LONG_READ_COLOR],
    # 'pb_gte1.bed': ['PacBio (Depth 1)', LONG_READ_COLOR],
    # 'np_gte1.bed': ['Nanopore (Depth 1)', LONG_READ_COLOR],
    # 'pb_gte5_lt92.bed': ['PacBio', LONG_READ_COLOR],
    # 'np_gte5_lt74.bed': ['Nanopore', LONG_READ_COLOR],
    # 'int_gte10_lt2xm.bed': ['NP/PB Intersection', LONG_READ_COLOR],
    # 'giab_high_conf.bed': ['GIAB High Confidence', SHORT_READ_COLOR],
    # 'gatk_gvcf_callable.bed': ['GATK GVCF Callable', SHORT_READ_COLOR],
    # 'gatk_cloci_callable.bed': ['GATK CLoci Callable', SHORT_READ_COLOR],

    # 'uni_gte5.bed': ['NP/PB Union (Depth 5)', '#F00060'],
    # 'pb_gte5_lt92.bed': ['PacBio (Depth 5)', '#C00000'],
    # 'np_gte5_lt74.bed': ['Nanopore (Depth 5)', '#C00060'],

    'uni_gte1.bed': ['Long Read Mappable', '#FF0F3F'],
    'pb_gte1.bed': ['PacBio Mappable', '#CF2F0F'],
    'np_gte1.bed': ['Nanopore Mappable', '#CF0F3F'],
    'int_gte10_lt2xm.bed': ['Long Read Callable', '#9F0F3F'],
    'giab_high_conf.bed': ['GIAB High Confidence', '#0F0F9F'],
    'gatk_gvcf_callable.bed': ['GATK Callable', '#3F0FCF'],
    'gatk_cloci_callable.bed': ['Short Read Mappable', '#6F0FFF'],

    'hg38_alu.bed': ['Alu', BASELINE_COLOR],
    'hg38_gencode24.bed': ['Gencode v24', BASELINE_COLOR],
    'hg38_repeatMasker.bed': ['Repeat Masker', BASELINE_COLOR],
    'hg38_segmentalDups.bed': ['Segmental Duplications', BASELINE_COLOR],
    'hg38_gencode27.bed': ['Gencode v27', BASELINE_COLOR],
    'hg38_gencode27_exome.bed': ['Gencode v27 (Exome)', BASELINE_COLOR],
    'hg38_gencode27_basic.bed': ['Gencode v27', BASELINE_COLOR],
    'hg38_tandemRepeats.bed': ['Tandem Repeats', BASELINE_COLOR],
    'Alu.bed': ['Alu', BASELINE_COLOR],
    'Centromere.bed': ['Centromeres', BASELINE_COLOR],
    'genomic_regions_definitions_modeledcentromere.bed': ['Centromeres', BASELINE_COLOR],
    'LINE.bed': ['LINE', BASELINE_COLOR],
    'Satellite.bed': ['Satellite DNA', BASELINE_COLOR],
    'SimpleRepeat.bed': ['Satellite DNA', BASELINE_COLOR],
    'SINE.bed': ['SINE', BASELINE_COLOR],
    'XINE.bed': ['LINEs and SINEs', BASELINE_COLOR],
    'Telomere.bed': ['Telomeres', BASELINE_COLOR],

    'INC_gencode.INC_repeatmask.bed': ['GenExon RMSK', BASELINE_COLOR],

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
            if included is None: included = 0.0
            line = line.split()
            assert len(line) == 5
            if int(line[1]) > 0:
                included += float(line[4])
    assert included is not None, "malformed file: {}".format(coverage_file)
    return included


def format_plot_element(bed_identifier, label, color, bed_location):
    return {
        IDENTIFIER_KEY: bed_identifier,
        COLOR_KEY: color,
        LABEL_KEY: label,
        BED_LOC_KEY: bed_location,
        VALUE_KEY: get_genome_coverage(bed_location) if bed_location is not None else 1.0
    }


def plot_all_feature_coverage(coverage_beds, feature_beds, show=False, save=True, sort=True):

    # get genome elements
    genome_elements = []
    for bed in coverage_beds:
        genome_elements.append(format_plot_element(bed, BED_INFO[bed][LABEL_IDX], BED_INFO[bed][COLOR_IDX],
                                                   get_coverage_bed_location(bed)))
    genome_elements.append(format_plot_element('genome', "Genome", BASELINE_COLOR, None))

    # get feature elements
    # feature_elements = [format_plot_element('genome', "Genome", BASELINE_COLOR, None)]
    # feature_elements = []
    # for bed in feature_beds:
    #     feature_elements.append(format_plot_element(bed, BED_INFO[bed][LABEL_IDX], BED_INFO[bed][COLOR_IDX],
    #                                                 get_feature_bed_location(bed)))
    # feature_elements.reverse()

    # get coverage elements
    coverage_elements = {f:[] for f in feature_beds}
    for feature_bed in feature_beds:
        for coverage_bed in coverage_beds:
            coverage_elements[feature_bed].append(
                format_plot_element(coverage_bed, BED_INFO[coverage_bed][LABEL_IDX], BED_INFO[coverage_bed][COLOR_IDX],
                                    get_bed_intersect_location(coverage_bed, feature_bed)))
        coverage_elements[feature_bed].append(format_plot_element(feature_bed, BED_INFO[feature_bed][LABEL_IDX],
                                                                  BASELINE_COLOR, get_feature_bed_location(feature_bed)))


    # create figure
    fig = plt.figure(1)
    # set up subplot grid
    n_row_start = 1
    n_cols = 2
    n_rows = int(n_row_start + math.ceil(len(feature_beds) / 2.0))
    grid_size = (n_rows, n_cols)
    gridspec.GridSpec(n_rows, n_cols)

    # # genome feature coverage
    # plt.subplot2grid(grid_size, (0, 0), colspan=2, rowspan=n_row_start - 1)
    # plt.title('Genome Feature Coverage')
    #
    # # legend
    # from matplotlib.lines import Line2D
    # legend_elements = [Line2D([0], [0], color=BED_INFO[x][COLOR_IDX], lw=16, label=BED_INFO[x][LABEL_IDX]) for x in coverage_beds]
    # legend_elements.reverse()
    # plt.legend(handles=legend_elements, bbox_to_anchor=(0.6, 0.8))
    #
    # # data
    # xs, ys, elements = [], [] , feature_elements
    # for i, v in enumerate(elements):
    #     xs.append(i)
    #     ys.append(v[VALUE_KEY])
    # plt.barh(xs, ys, color=list(map(lambda x: x[COLOR_KEY], elements)))
    # plt.yticks(xs, list(map(lambda x: x[LABEL_KEY], elements)))
    # plt.xlim(0,1.07)
    # plt.xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    # for x, v in enumerate(ys):
    #     plt.text(1.02, x - .1, ("%.3f" % v), color=list(map(lambda x: x[COLOR_KEY], elements))[x])


    # genome read coverage
    plt.subplot2grid(grid_size, (n_row_start - 1, 0), colspan=2, rowspan=1)
    plt.title('Whole Genome')
    # data
    xs, ys, elements = [], [] , genome_elements
    if sort: elements.sort(key=lambda x: x[VALUE_KEY])
    for i, v in enumerate(elements):
        xs.append(i)
        ys.append(v[VALUE_KEY])

    plt.barh(xs, ys, color=list(map(lambda x: x[COLOR_KEY], elements)))
    plt.yticks(xs, list(map(lambda x: x[LABEL_KEY], elements)))
    plt.xlim(0,1.07)
    plt.xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    for x, v in enumerate(ys):
        plt.text(1.02, x - .1, ("%.3f" % v), color=list(map(lambda x: x[COLOR_KEY], elements))[x])


    # specific feature coverage
    for n, feature in enumerate(feature_beds):
        y = int(n_row_start + math.floor(n / 2.0))
        x = int(n % 2)

        # subplot managment
        plt.subplot2grid(grid_size, (y, x))
        plt.title(BED_INFO[feature][LABEL_IDX])

        # data
        elements = coverage_elements[feature]
        if sort: elements.sort(key=lambda x: x[VALUE_KEY])
        max_value = get_genome_coverage(get_feature_bed_location(feature))
        xs = []
        ys = []
        for i, v in enumerate(elements):
            xs.append(i)
            ys.append(1.0 * v[VALUE_KEY] / max_value)

        plt.barh(xs, ys, color=list(map(lambda x: x[COLOR_KEY], elements)))
        # plt.xlim(0,1.17)
        # plt.yticks(xs, ['' for _ in xs])
        plt.xlim(0,1.35)
        plt.yticks(xs, list(map(lambda x: x[LABEL_KEY], elements)))
        plt.xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        for i, v in enumerate(ys):
            element = elements[i]
            # plt.text(1.02, i - .1, ("%-2.1f%%   (%-2.1f%%)" % (100.0 * v, 100.0 * element[VALUE_KEY])), color=element[COLOR_KEY])
            plt.text(1.02, i - .1, ("%.3f   (%.3f)" % (v, element[VALUE_KEY])), color=element[COLOR_KEY])

    # sizing and saving
    # fig.tight_layout()
    fig.set_size_inches(h=3*n_rows, w=10*n_cols)
    if show:
        plt.show()
    if save:
        fig.savefig(os.path.join(PLOT_BED_DIR, 'all_features.png'), dpi=200)


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
    ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xlim(0, max(ys) * 1.25)
    for x, v in enumerate(ys):
        ax.text(max(ys) * 1.005, x-.1, ("%.3f" % v), color=list(map(lambda x: x[COLOR_KEY], elements))[x])


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
    # plot_feature_coverage(COVERAGE_BEDS, None)

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

    plot_all_feature_coverage(COVERAGE_BEDS, FEATURE_BEDS)


if __name__ == "__main__":
    main()