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
    char *paramsFile = NULL;

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

    stList *profileSequences = NULL;
    /*
     * TODO: Use htslib to parse the reads within an input interval of a reference sequence of a bam file
     * and create a list of profile sequences using
     * stProfileSeq *stProfileSeq_constructEmptyProfile(char *referenceName,
            int64_t referenceStart, int64_t length);
       where for each position you turn the character into a profile probability, as shown in the tests

       In future we can use the mapq scores for reads to adjust these profiles, or for signal level alignments
       use the posterior probabilities.
     */

    // Parse any model parameters
    st_logInfo("Parsing model parameters\n");

    stRPHmmParameters *params = NULL;
    /*
     * TODO: Get model parameters. I suggest we make a simple json or yaml file to hold these parameters.
     * Minimally we need a heterozygozity rate (the fraction of reference positions that are different between the
     * haplotypes (excluding gaps)
     * We also need to figure out what to do with gap positions (which are just treated as an additional character)
     * Once these are read in we need to construct (as shown in the tests) the different matrices.
     * We will need:
     * stRPHmmParameters *stRPHmmParameters_construct(uint16_t *hetSubModel,
            double *hetSubModelSlow,
            uint16_t *readErrorSubModel,
            double *readErrorSubModelSlow,
            bool maxNotSumTransitions,
            int64_t maxPartitionsInAColumn,
            int64_t maxCoverageDepth);
       And:
       void setSubstitutionProb(uint16_t *logSubMatrix, double *logSubMatrixSlow,
            int64_t sourceCharacterIndex,
            int64_t derivedCharacterIndex, double prob);
     */

    // Create HMMs
    st_logInfo("Creating read partitioning HMMs\n");

    stList *hmms = getRPHmms(profileSequences, params);

    // Create HMM traceback and genome fragments

    // For each read partitioning HMM
    for(int64_t i=0; i<stList_length(hmms); i++) {
        stRPHmm *hmm = stList_get(hmms, i);

        st_logInfo("Creating genome fragment for reference sequence: %s, start: %" PRIi64 ", length: %" PRIi64 "\n",
                    hmm->referenceName, hmm->refStart, hmm->refLength);

        // Run the forward-backward algorithm
        stRPHmm_forwardBackward(hmm);

        // Now compute a high probability path through the hmm
        stList *path = stRPHmm_forwardTraceBack(hmm);

        /*
         * TODO: BENEDICT: Break up the path in places where we can not be sure about the phasing
         */

        // Compute the genome fragment
        stGenomeFragment *gF = stGenomeFragment_construct(hmm, path);

        // Write out VCF
        st_logInfo("Writing out VCF for fragment\n");
        /*
         * TODO: Convert the genome fragment into a portion of a VCF file (we'll need to write the header out earlier)
         * We can express the genotypes and (using phase sets) the phasing relationships.
         */

        // Optionally write out two BAMs, one for each read partition
        st_logInfo("Writing out BAM partitions for fragment\n");
        /*
         * TODO: Optionally, write out a new BAM file expressing the partition (which reads belong in which partition)
         * Not sure if we need to write out multiple files or if we can add a per read flag to express this information.
         */
    }

    // Cleanup
    stList_destruct(hmms);

    //while(1); // Use this for testing for memory leaks

    return 0;
}

