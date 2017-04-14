/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com) & Arthur Rand (arand@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <getopt.h>
#include <stdio.h>
#include <ctype.h>
#include "stRPHmm.h"
#include "sam.h"
#include "bgzf.h"

void usage() {
    fprintf(stderr, "marginPhase [options] BAM_FILE\n");
    fprintf(stderr,
            "Phases the reads in an interval of a BAM file reporting a gVCF file "
            "giving genotypes and haplotypes for region.\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
    fprintf(stderr, "-b --bamFile : Input BAM file\n");
    fprintf(stderr, "-v --vcfFile : Output VCF file\n");
    fprintf(stderr, "-p --paramsFile : Input params file\n");
    fprintf(stderr, "-n --refSeqName : Name of reference sequence\n");
    fprintf(stderr, "-s --intervalStart : Starting position of interval to read\n");
    fprintf(stderr, "-e --intervalEnd : Ending position of interval to read\n");
}

char baseToAlphabet(char b) {
    // As far as I can tell, it wasn't actually specified which of '0,1,2,3' correspond
    // to 'A,C,G,T', so I chose arbitrarily
    if (b == 'A') return FIRST_ALPHABET_CHAR;
    if (b == 'C') return FIRST_ALPHABET_CHAR + 1;
    if (b == 'G') return FIRST_ALPHABET_CHAR + 2;
    if (b == 'T') return FIRST_ALPHABET_CHAR + 3;
    return  FIRST_ALPHABET_CHAR - 1;
}

void parseReads(stList *profileSequences, char *bamFile, char *refSeqName, int32_t intervalStart, int32_t intervalEnd) {

    /*
     * TODO: Use htslib to parse the reads within an input interval of a reference sequence of a bam file
     * and create a list of profile sequences using
     * stProfileSeq *stProfileSeq_constructEmptyProfile(char *referenceName,
            int64_t referenceStart, int64_t length);
       where for each position you turn the character into a profile probability, as shown in the tests

       In future we can use the mapq scores for reads to adjust these profiles, or for signal level alignments
       use the posterior probabilities.
     */

    st_logDebug("Bam file: %s \n", bamFile);
    BGZF *in = bgzf_open(bamFile,"r"); //open bam file
    bam_hdr_t *bamHdr = bam_hdr_read(in); //read header
    bam1_t *aln = bam_init1(); //initialize an alignment

    bool readWholeFile = false;
    if (!refSeqName) {
        st_logInfo("No reference sequence name given, reading whole file.\n");
        readWholeFile = true;
    }
    if (intervalStart < 0) {
        intervalStart = 0;
        st_logInfo("Reading from start of chromosome.\n");
    }
    if (intervalEnd < 0) {
        st_logInfo("Reading until end of chromosome.\n", intervalEnd);
    }
    int32_t readCount = 0;

    while(bam_read1(in,aln) > 0){

        int32_t pos = aln->core.pos +1; //left most position of alignment in zero based coordianate (+1)
        char *chr = bamHdr->target_name[aln->core.tid] ; //contig name (chromosome)
        uint32_t len = aln->core.l_qseq; //length of the read.
        uint8_t *q = bam_get_seq(aln);  // DNA sequence

        // TODO: is this the right way to deal with the reference sequence / input interval?
        if(readWholeFile || strcmp(chr, refSeqName) == 0) {
            // Should reads that cross boundaries be counted?
            if (pos >= intervalStart && (intervalEnd < 0 || pos + len <= intervalEnd)) {
                readCount++;

                stProfileSeq *seq = stProfileSeq_constructEmptyProfile(chr, pos, len);

                for(int64_t i=0; i<len; i++) {
                    // For each position turn character into profile probability
                    // As is, this makes the probability 1 for the base read by the file, and 0 otherwise
                    // Should this be modified to take into account error rates?
                    char b = baseToAlphabet(seq_nt16_str[bam_seqi(q,i)]);
                    seq->profileProbs[i*ALPHABET_SIZE + b - FIRST_ALPHABET_CHAR] = ALPHABET_MAX_PROB;
                }
                stList_append(profileSequences, seq);
            }
        }


        bool withinInputInterval = true;
        if (withinInputInterval) {

            // stProfileSeq_print(seq, stderr, true);
        }
    }
    st_logDebug("Number of reads found with name: %s in interval %d - %d : %d\n", refSeqName, intervalStart, intervalEnd, readCount);

    bam_destroy1(aln);
    bgzf_close(in);
}

int main(int argc, char *argv[]) {
    // Parameters / arguments

    char * logLevelString = NULL;
    char *bamFile = NULL;
    char *vcfFile = NULL;
    char *paramsFile = NULL;

    char *refSeqName = NULL;
    int32_t intervalStart = -1;
    int32_t intervalEnd = -1;

    // Parse the options
    while (1) {
        static struct option long_options[] = {
                { "logLevel", required_argument, 0, 'a' },
                { "help", no_argument, 0, 'h' },
                { "bamFile", required_argument, 0, 'b'},
                { "vcfFile", required_argument, 0, 'v'},
                { "paramsFile", required_argument, 0, 'p'},
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

    // Parse reads for interval
    st_logInfo("Parsing input reads\n");

    stList *profileSequences = stList_construct();
    parseReads(profileSequences, bamFile, refSeqName, intervalStart, intervalEnd);


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

    // Break up the hmms where the phasing is uncertain
    st_logInfo("Breaking apart HMMs where the phasing is uncertain\n");

    stList *l = stList_construct3(0, (void (*)(void *))stRPHmm_destruct2);
    while(stList_length(hmms) > 0) {
        stList_appendAll(l, stRPHMM_splitWherePhasingIsUncertain(stList_pop(hmms)));
    }
    hmms = l;

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

