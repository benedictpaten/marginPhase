/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com) & Arthur Rand (arand@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <getopt.h>
#include <stdio.h>
#include <ctype.h>
#include "stRPHmm.h"

void usage() {
    fprintf(stderr, "marginPhase [options] BAM_FILE\n");
    fprintf(stderr,
            "Phases the reads in an interval of a BAM file reporting a gVCF file "
            "giving genotypes and haplotypes for region.\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
}

int main(int argc, char *argv[]) {
    // Parameters / arguments
    char * logLevelString = NULL;
    char *bamFile = NULL;
    char *vcfFile = NULL;

    // Parse the options
    while (1) {
        static struct option long_options[] = {
                { "logLevel", required_argument, 0, 'a' },
                { "help", no_argument, 0, 'h' },
                { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:h", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
        case 'a':
            logLevelString = stString_copy(optarg);
            st_setLogLevelFromString(logLevelString);
            break;
        case 'h':
            usage();
            return 0;
        default:
            usage();
            return 1;
        }
    }

    st_setLogLevelFromString(logLevelString);

    // Parse reads for interval
    st_logInfo("Parsing input reads\n");

    stProfileSeq *stProfileSeq_constructEmptyProfile(char *referenceName,
            int64_t referenceStart, int64_t length);

    // Parse any model parameters
    st_logInfo("Parsing model parameters\n");

    stRPHmmParameters *stRPHmmParameters_construct(uint16_t *hetSubModel,
            double *hetSubModelSlow,
            uint16_t *readErrorSubModel,
            double *readErrorSubModelSlow,
            bool maxNotSumTransitions,
            int64_t maxPartitionsInAColumn,
            int64_t maxCoverageDepth);

    void setSubstitutionProb(uint16_t *logSubMatrix, double *logSubMatrixSlow,
            int64_t sourceCharacterIndex,
            int64_t derivedCharacterIndex, double prob);

    // Create HMMs
    st_logInfo("Creating read partitioning HMMs\n");

    // Here we will call

    stList *getRPHmms(stList *profileSeqs, stRPHmmParameters *params);

    // Create HMM traceback and genome fragments
    st_logInfo("Creating genome fragments\n");

    void stRPHmm_forwardBackward(stRPHmm *hmm);

    stList *stRPHmm_forwardTraceBack(stRPHmm *hmm);

    stSet *stRPHmm_partitionSequencesByStatePath(stRPHmm *hmm, stList *path, bool partition1);

    stGenomeFragment *stGenomeFragment_construct(stRPHmm *hmm, stList *path);

    // Write out VCF
    st_logInfo("Writing out VCF\n");

    // Optionally write out two BAMs, one for each read partition
    st_logInfo("Writing out BAM partitions\n");


    //while(1);

    return 0;
}

