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
    fprintf(stderr, "-b --bamFile : Input BAM file\n");
    fprintf(stderr, "-v --vcfFile : Output VCF file\n");
    fprintf(stderr, "-p --params : Input params file\n");
    fprintf(stderr, "-n --refSeqName : Name of reference sequence\n");
    fprintf(stderr, "-s --intervalStart : Starting position of interval to read\n");
    fprintf(stderr, "-e --intervalEnd : Ending position of interval to read\n");
}


int main(int argc, char *argv[]) {
    // Parameters / arguments
    char * logLevelString = NULL;
    char *bamFile = NULL;
    char *vcfFile = NULL;
    char *paramsFile = "params.json";

    char *refSeqName = NULL;
    int32_t intervalStart = -1;
    int32_t intervalEnd = -1;

    // When done testing, set random seed
    // st_randomSeed();

    // Parse the options
    while (1) {
        static struct option long_options[] = {
                { "logLevel", required_argument, 0, 'a' },
                { "help", no_argument, 0, 'h' },
                { "bamFile", required_argument, 0, 'b'},
                { "vcfFile", required_argument, 0, 'v'},
                { "params", required_argument, 0, 'p'},
                { "refSeqName", required_argument, 0, 'n'},
                { "intervalStart", required_argument, 0, 's'},
                { "intervalEnd", required_argument, 0, 'e'},
                { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:v:p:n:s:e:h", long_options, &option_index);

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
        case 'b':
            bamFile = stString_copy(optarg);
            break;
        case 'v':
            vcfFile = stString_copy(optarg);
            break;
        case 'p':
            paramsFile = stString_copy(optarg);
            break;
        case 'n':
            refSeqName = stString_copy(optarg);
            break;
        case 's':
            intervalStart = atoi(optarg);
            break;
        case 'e':
            intervalEnd = atoi(optarg);
            break;
        default:
            usage();
            return 1;
        }
    }

    // Parse any model parameters
    st_logInfo("Parsing model parameters\n");

    char *alphabet[ALPHABET_SIZE];
    char *wildcard;
    stRPHmmParameters *params = parseParameters(paramsFile, alphabet, &wildcard);

    st_logDebug("Alphabet array: \n");
    for (int j = 0; j < ALPHABET_SIZE; j++ ) {
        st_logDebug("%d: %s \t%d\n", j, alphabet[j], strlen(alphabet[j]));
    }
    st_logDebug("Wildcard: %s\n", wildcard);


    // Parse reads for interval
    st_logInfo("Parsing input reads\n");

    stList *profileSequences = stList_construct();
    parseReads(profileSequences, bamFile, alphabet, wildcard, refSeqName, intervalStart, intervalEnd);



    // Create HMMs
    st_logInfo("Creating read partitioning HMMs\n");
    stList *hmms = getRPHmms(profileSequences, params);

    // Break up the hmms where the phasing is uncertain
    st_logInfo("Breaking apart HMMs where the phasing is uncertain\n");

    stList *l = stList_construct3(0, (void (*)(void *))stRPHmm_destruct2);
    while(stList_length(hmms) > 0) {
        stList_appendAll(l, stRPHMM_splitWherePhasingIsUncertain(stList_pop(hmms)));
    }
    hmms = l;

    // Create HMM traceback and genome fragments

    st_logDebug("Created %d hmms \n", stList_length(hmms));
    // For each read partitioning HMM
    for(int64_t i=0; i<stList_length(hmms); i++) {
        stRPHmm *hmm = stList_get(hmms, i);

        st_logInfo("Creating genome fragment for reference sequence: %s, start: %" PRIi64 ", length: %" PRIi64 "\n",
                    hmm->referenceName, hmm->refStart, hmm->refLength);

        // Run the forward-backward algorithm
        stRPHmm_forwardBackward(hmm);

        // Now compute a high probability path through the hmm
        stList *path = stRPHmm_forwardTraceBack(hmm);

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

