from __future__ import print_function
import argparse
import glob
import sys
import os

def parse_args():
    parser = argparse.ArgumentParser("Consolidates output logs for chunked marginPhase run")
    parser.add_argument('--log_dir', '-l', dest='log_dir', default='.', type=str,
                       help='Location where marginPhase .log files are located.')

    return parser.parse_args()

def log(msg, include_stdout_also = False):
    print(msg, file=sys.stderr)
    if include_stdout_also:
        print(msg, file=sys.stdout)

class LogDetails():
    def __init__(self, log_path):
        self.log_path = log_path
        self.log_name = os.path.basename(log_path)
        self.chunk_idx = int(self.log_name.split(".")[-2])
        
        # read in lines, save results
        self.log_contents = []
        lines_for_analysis = []
        with open(self.log_name, 'r') as contents:
            found_results = False
            for line in contents:
                if len(line.rstrip()) == 0: continue
                if found_results:
                    self.log_contents.append(line.rstrip())
                    lines_for_analysis.append(line.rstrip())
                elif line.find("RESULTS") != -1:
                    found_results = True
                    
        # stats
        self.genome_fragments = -1
        self.average_fragment_len = -1.0
        self.sensitivity = -1.0
        self.sensitivity_no_indels = -1.0
        self.sens_true_positive = -1
        self.sens_ref_true_positive = -1
        self.sens_true_positive_no_indels = -1
        self.sens_ref_true_positive_no_indels = -1
        self.homozygous_variants_in_sample = -1
        self.false_negatives = -1
        self.specificity = -1.0
        self.spec_true_negatives = -1
        self.spec_total_reference_positions = -1
        self.incorrect_positives = -1
        self.false_positives = -1
        self.false_positives_with_gaps = -1
        self.false_negatives_bad_partition = -1
        self.false_negatives_indel_missed = -1
        self.switch_error_numerator = -1
        self.switch_error_denominator = -1
        self.average_distance_between_switch_errors = -1
        self.uncertain_phasing_spots = -1
        self.load_details(lines_for_analysis)

    def load_details(self, lines_for_analysis):
        lines_for_analysis.reverse()
        line = lines_for_analysis.pop().strip().split()
        self.genome_fragments = int(line[5])
        self.average_fragment_len = float(line[11])

        line = lines_for_analysis.pop().strip().split()
        self.sensitivity = float(line[1].rstrip(","))
        self.sensitivity_no_indels = float(line[4])

        line = lines_for_analysis.pop().strip().split()
        self.sens_true_positive = int(line[8])
        self.sens_ref_true_positive = int(line[11])
        self.sens_true_positive_no_indels = int(line[13])
        self.sens_ref_true_positive_no_indels = int(line[16].rstrip(")"))

        line = lines_for_analysis.pop().strip().split()
        self.homozygous_variants_in_sample = int(line[4])

        line = lines_for_analysis.pop().strip().split()
        self.false_negatives = int(line[2].rstrip(","))

        line = lines_for_analysis.pop().strip().split()
        self.specificity = float(line[1])

        line = lines_for_analysis.pop().strip().split()
        self.spec_true_negatives = int(line[8])
        self.spec_total_reference_positions = int(line[11].rstrip(")"))

        line = lines_for_analysis.pop().strip().split()
        self.incorrect_positives = int(line[2])

        line = lines_for_analysis.pop().strip().split()
        self.false_positives = int(line[2].rstrip(","))
        self.false_positives_with_gaps = int(line[5])

        line = lines_for_analysis.pop().strip().split()
        line = lines_for_analysis.pop().strip().split()
        self.false_negatives_bad_partition = int(line[2])

        line = lines_for_analysis.pop().strip().split()
        self.false_negatives_indel_missed = int(line[2])

        line = lines_for_analysis.pop().strip().split()
        line = lines_for_analysis.pop().strip().split()
        self.switch_error_numerator = int(line[4].lstrip("("))
        self.switch_error_denominator = int(line[7].rstrip(","))

        line = lines_for_analysis.pop().strip().split()
        self.average_distance_between_switch_errors = None if line[5].find("nan") != -1 else int(line[5])

        line = lines_for_analysis.pop().strip().split()
        self.uncertain_phasing_spots = int(line[3])

"""
----- RESULTS -----
There were a total of 13 genome fragments. Average length = 2824.615479
Sensitivity: 0.447368, 	 without indels: 0.478873
	(= fraction of true positives compared to reference, 	34 out of 76 / 34 out of 71)
	Homozygous variants in sample: 27
	False negatives: 28
Specificity: 0.999554
	(= fraction of true negatives compared to reference, 	35869 out of  35885)
	Incorrect positives: 0
	False positives: 16,	with gaps: 6
False negatives:
	Partition bad: 23 		(0.821429)
	Indel missed: 5 		(0.178571)
Phasing:
	Switch error rate: 0.000000 	 (0 out of 30, fraction correct: 1.000000)
	Average distance between switch errors: -nan
	Uncertain phasing spots: 4
"""

def main():
    args = parse_args()

    logfiles = glob.glob(os.path.join(args.log_dir, "*.log"))
    log_count = len(logfiles)
    if log_count == 0:
        log("Found no log files matching %s" % os.path.join(args.log_dir, "*.log"))
        sys.exit(0)
    else:
        log("Fond %d log files in %s" % (log_count, args.log_dir))

    log_details = list()
    for l in logfiles:
        log_details.append(LogDetails(l))
    log_details.sort(key=lambda x: x.chunk_idx)

    # total_reference_positions = sum(list(map(lambda x: x., log_details)))
    spec_total_reference_positions = sum(list(map(lambda x: x.spec_total_reference_positions, log_details)))
    spec_true_negatives = sum(list(map(lambda x: x.spec_true_negatives, log_details)))

    sens_true_positive = sum(list(map(lambda x: x.sens_true_positive, log_details)))
    sens_ref_true_positive = sum(list(map(lambda x: x.sens_ref_true_positive, log_details)))

    sens_true_positive_no_indels = sum(list(map(lambda x: x.sens_true_positive_no_indels, log_details)))
    sens_ref_true_positive_no_indels = sum(list(map(lambda x: x.sens_ref_true_positive_no_indels, log_details)))

    switch_error_numerator = sum(list(map(lambda x: x.switch_error_numerator, log_details)))
    switch_error_denominator = sum(list(map(lambda x: x.switch_error_denominator, log_details)))

    uncertain_phasing_spots = sum(list(map(lambda x: x.uncertain_phasing_spots, log_details)))

    for l in log_details:
        print("\n%3d:\t%s" % (l.chunk_idx, l.log_name))
        for line in l.log_contents:
            print("\t%s" % line)
    log("\n\nSUMMARY:", True)
    log("\tAnalyzed %d chunks" % len(log_details), True)
    log("Specificity: %f" % (1.0 * spec_true_negatives / spec_total_reference_positions), True)
    log("\t%d / %d" % (spec_true_negatives, spec_total_reference_positions), True)
    log("Sensitivity: %f " % (1.0 * sens_true_positive / sens_ref_true_positive), True)
    log("\t%d / %d " % (sens_true_positive, sens_ref_true_positive), True)
    log("Sensitivity: %f  (no indels) " % (1.0 * sens_true_positive_no_indels / sens_ref_true_positive_no_indels), True)
    log("\t%d / %d " % (sens_true_positive_no_indels, sens_ref_true_positive_no_indels), True)
    log("Switch Error: %f" % (1.0 * switch_error_numerator / switch_error_denominator), True)
    log("\t%d / %d" % (switch_error_numerator, switch_error_denominator), True)
    log("Uncertain phasing spots: %d" % uncertain_phasing_spots, True)






# self.genome_fragments = -1
# self.average_fragment_len = -1.0
# self.sensitivity = -1.0
# self.sensitivity_no_indels = -1.0
# self.sens_true_positive = -1
# self.sens_ref_true_positive = -1
# self.sens_true_positive_no_indels = -1
# self.sens_ref_true_positive_no_indels = -1
# self.homozygous_variants_in_sample = -1
# self.false_negatives = -1
# self.specificity = -1.0
# self.spec_true_negatives = -1
# self.spec_total_reference_positions = -1
# self.incorrect_positives = -1
# self.false_positives = -1
# self.false_positives_with_gaps = -1
# self.false_negatives_bad_partition = -1
# self.false_negatives_indel_missed = -1
# self.switch_error_numerator = -1
# self.switch_error_denominator = -1
# self.average_distance_between_switch_errors = -1
# self.uncertain_phasing_spots = -1




if __name__ == "__main__":
    main()