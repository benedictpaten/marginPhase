/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com) & Arthur Rand (arand@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <getopt.h>
#include <stdio.h>
#include <ctype.h>
#include <memory.h>

#include "stRPHmm.h"

/*
 * Functions used to print debug output.
 */

void printSequenceStats(FILE *fH, stList *profileSequences) {
    /*
     * Print stats about the set of profile sequences.
     */
    int64_t totalLength=0;
    for(int64_t i=0; i<stList_length(profileSequences); i++) {
        stProfileSeq *profileSeq = stList_get(profileSequences, i);
        totalLength += profileSeq->length;
    }
    fprintf(fH, "Got %" PRIi64 " profile sequences, with total length: %" PRIi64 ", average length: %f\n",
            stList_length(profileSequences), totalLength, ((float)totalLength)/stList_length(profileSequences));
}

double getExpectedNumberOfMatches(uint64_t *haplotypeString, int64_t start, int64_t length, stProfileSeq *profileSeq) {
    /*
     * Returns the expected number of positions in the profile sequence that are identical to the given haplotype string.
     */
    double totalExpectedMatches = 0.0;

    for(int64_t i=0; i<profileSeq->length; i++) {
        // Get base in the haplotype sequence
        int64_t j = i + profileSeq->refStart - start;
        if(j >= 0 && j < length) {
            uint64_t hapBase = haplotypeString[j];
            assert(hapBase < ALPHABET_SIZE);

            // Expectation of a match
            totalExpectedMatches += getProb(&(profileSeq->profileProbs[i * ALPHABET_SIZE]), hapBase);
        }
    }
    return totalExpectedMatches;
}

double getExpectedIdentity(uint64_t *haplotypeString, int64_t start, int64_t length, stSet *profileSeqs) {
    /*
     * Returns the expected fraction of positions in the profile sequences that match their corresponding position in the
     * given haplotype string.
     */
    double totalExpectedNumberOfMatches = 0.0;
    int64_t totalLength = 0;
    stSetIterator *it = stSet_getIterator(profileSeqs);
    stProfileSeq *pSeq;
    while((pSeq = stSet_getNext(it)) != NULL) {
        totalExpectedNumberOfMatches += getExpectedNumberOfMatches(haplotypeString, start, length, pSeq);
        totalLength += pSeq->length;
    }
    return totalExpectedNumberOfMatches/totalLength;
}

double getIdentityBetweenHaplotypes(uint64_t *hap1String, uint64_t *hap2String, int64_t length) {
    /*
     * Returns the fraction of positions in two haplotypes that are identical.
     */
    int64_t totalMatches = 0;
    for(int64_t i=0; i<length; i++) {
        if(hap1String[i] == hap2String[i]) {
            totalMatches++;
        }
    }
    return ((double)totalMatches)/length;
}

double getIdentityBetweenHaplotypesExcludingIndels(uint64_t *hap1String, uint64_t *hap2String, int64_t length) {
    /*
     * Returns the fraction of positions in two haplotypes that are identical.
     */
    int64_t totalMatches = 0;
    int64_t numGaps = 0;
    for(int64_t i=0; i<length; i++) {
        if(hap1String[i] == hap2String[i]) {
            totalMatches++;
        } else if (hap1String[i] == ALPHABET_SIZE-1 || hap2String[i] == ALPHABET_SIZE-1) {
            numGaps++;
        }
    }
    return ((double)totalMatches)/(length - numGaps);
}

void getExpectedMatcheBetweenProfileSeqs(stProfileSeq *pSeq1, stProfileSeq *pSeq2, int64_t *totalAlignedPositions, double *totalExpectedMatches) {
    /*
     * Calculates the number of base overlaps and expected base matches between two profile sequences.
     */

    for(int64_t i=0; i<pSeq1->length; i++) {
        // Establish if the coordinate is in both sequences
        int64_t j = i + pSeq1->refStart - pSeq2->refStart;
        if(j >= 0 && j < pSeq2->length) {
            (*totalAlignedPositions)++;

            // Calculate expectation of match
            for(int64_t k=0; k<ALPHABET_SIZE; k++) {
                double e1 = getProb(&(pSeq1->profileProbs[i * ALPHABET_SIZE]), k);
                double e2 = getProb(&(pSeq2->profileProbs[j * ALPHABET_SIZE]), k);
                assert(e1 * e2 <= 1.0);
                *totalExpectedMatches += e1 * e2;
            }
        }
    }
}

void printAvgIdentityBetweenProfileSequences(FILE *fH, stList *profileSequences, int64_t maxSeqs) {
    /*
     * Prints the average base identity between pairwise base overlaps between the given profile sequences
     */

    double totalExpectedMatches = 0.0;
    int64_t totalAlignedPositions = 0;

    for(int64_t i=0; i<stList_length(profileSequences) && i<maxSeqs; i++) {
        for(int64_t j=i+1; j<stList_length(profileSequences) && j<maxSeqs; j++) {
            getExpectedMatcheBetweenProfileSeqs(stList_get(profileSequences, i), stList_get(profileSequences, j),
                                                &totalAlignedPositions, &totalExpectedMatches);

        }
    }

    fprintf(fH, "Avg. pairwise identity between profile sequences: %f measured at %" PRIi64 " overlapping sites\n",
            totalExpectedMatches/totalAlignedPositions, totalAlignedPositions);
}

double *getHaplotypeBaseComposition(uint64_t *hapString, int64_t length) {
    /*
     * Get the count of each alphabet character in the haplotype sequence, returned
     * as an array.
     */
    double *baseCounts = st_calloc(ALPHABET_SIZE, sizeof(double));
    for(int64_t i=0; i<length; i++) {
        baseCounts[hapString[i]] += 1;
    }
    return baseCounts;
}

double *getExpectedProfileSequenceBaseComposition(stSet *profileSeqs) {
    /*
     * Get the expected count of each alphabet character in the profile sequences, returned
     * as an array.
     */
    double *baseCounts = st_calloc(ALPHABET_SIZE, sizeof(double));
    stSetIterator *it = stSet_getIterator(profileSeqs);
    stProfileSeq *pSeq;
    while((pSeq = stSet_getNext(it)) != NULL) {
        for(int64_t i=0; i<pSeq->length; i++) {
            for(int64_t j=0; j<ALPHABET_SIZE; j++) {
                baseCounts[j] += getProb(&(pSeq->profileProbs[i*ALPHABET_SIZE]), j);
            }
        }
    }
    return baseCounts;
}

void printBaseComposition(FILE *fH, double *baseCounts) {
    /*
     * Print the counts/fraction of each alphabet character.
     */
    double totalCount = 0;
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        totalCount += baseCounts[i];
    }
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        fprintf(fH, "Base %" PRIi64 " count: %f fraction: %f\n", i, baseCounts[i], baseCounts[i]/totalCount);
    }
}

void addProfileSeqIdsToSet(stSet *pSeqs, stSet *readIds) {
    stSetIterator *it = stSet_getIterator(pSeqs);
    stProfileSeq *pSeq;
    while((pSeq = stSet_getNext(it)) != NULL) {
        stSet_insert(readIds, pSeq->readId);
    }
    stSet_destructIterator(it);
}

/*
 * Main functions
 */


void usage() {
fprintf(stderr, "marginPhase BAM_FILE REFERENCE_FASTA [options]\n");
    fprintf(stderr,
            "Phases the reads in an interval of a BAM file (BAM_FILE) reporting a gVCF file "
            "giving genotypes and haplotypes for region.\n"
            "REFERENCE_FASTA is the reference sequence for the region in fasta format.\n");
    fprintf(stderr, "-a --logLevel   : Set the log level [default = info]\n");
    fprintf(stderr, "-h --help       : Print this help screen\n");
    fprintf(stderr, "-o --outSamBase : Output SAM Base (\"example\" -> \"example1.sam\", \"example2.sam\")\n");
    fprintf(stderr, "-v --vcfFile    : Output VCF file\n");
    fprintf(stderr, "-p --params     : Input params file\n");
    fprintf(stderr, "-l --iterationsOfParameterLearning : Iterations of parameter learning\n");
}

int main(int argc, char *argv[]) {
    // Parameters / arguments
    char *logLevelString = stString_copy("info");
    char *bamInFile = NULL;
    char *referenceFastaFile = NULL;

    char *samOutBase = NULL;
    char *vcfOutFile = "output.vcf";
    char *paramsFile = "params.json";
    int64_t iterationsOfParameterLearning = 0;

    // TODO: When done testing, optionally set random seed using st_randomSeed();

    if(argc < 3) {
        usage();
        return 0;
    }

    // Parse the options
    while (1) {
        static struct option long_options[] = {
                { "logLevel", required_argument, 0, 'a' },
                { "help", no_argument, 0, 'h' },
                { "outSamBase", required_argument, 0, 'o'},
                { "vcfFile", required_argument, 0, 'v'},
                { "params", required_argument, 0, 'p'},
                { "iterationsOfParameterLearning", required_argument, 0, 'l'},
                { "reference", required_argument, 0, 'r'},
                { 0, 0, 0, 0 } };

        int option_index = 0;
        int key = getopt_long(argc-2, &argv[2], "a:v:p:l:h", long_options, &option_index);

        if (key == -1) {
            break;
        }

        int i;

        switch (key) {
        case 'a':
            free(logLevelString);
            logLevelString = stString_copy(optarg);
            break;
        case 'h':
            usage();
            return 0;
        case 'o':
            samOutBase = stString_copy(optarg);
            break;
        case 'v':
            vcfOutFile = stString_copy(optarg);
            break;
        case 'p':
            paramsFile = stString_copy(optarg);
            break;
        case 'l':
            i = sscanf(optarg, "%" PRIi64 "", &iterationsOfParameterLearning);
            assert(i == 1);
            assert(iterationsOfParameterLearning >= 0);
            break;
        default:
            usage();
            return 0;
        }
    }
    st_setLogLevelFromString(logLevelString);
    free(logLevelString);

    bamInFile = stString_copy(argv[1]);
    referenceFastaFile = stString_copy(argv[2]);

    // Parse any model parameters
    st_logInfo("> Parsing model parameters from file: %s\n", paramsFile);
    stBaseMapper *baseMapper = stBaseMapper_construct();
    stRPHmmParameters *params = parseParameters(paramsFile, baseMapper);
    free(paramsFile);

    // Print a report of the parsed parameters
    if(st_getLogLevel() == debug) {
        stRPHmmParameters_printParameters(params, stderr);
    }

    // Parse reads for interval
    st_logInfo("> Parsing input reads from file: %s\n", bamInFile);
    stList *profileSequences = stList_construct3(0, (void (*)(void *))stProfileSeq_destruct);
    parseReads(profileSequences, bamInFile, baseMapper);

    // Print some stats about the input sequences
    if(st_getLogLevel() == debug) {
        printSequenceStats(stderr, profileSequences);
        printAvgIdentityBetweenProfileSequences(stderr, profileSequences, 100);
    }

    // Getting reference sequence priors
    st_logInfo("> Parsing prior probabilities on positions from reference sequences: %s\n", referenceFastaFile);
    stHash *referenceNamesToReferencePriors = createReferencePriorProbabilities(referenceFastaFile, profileSequences, baseMapper, params);

    // Learn the parameters for the input data
    st_logInfo("> Learning parameters for HMM model (%" PRIi64 ")\n", iterationsOfParameterLearning);
    stRPHmmParameters_learnParameters(params, profileSequences, referenceNamesToReferencePriors, iterationsOfParameterLearning);

    // Print a report of the parsed parameters
    if(st_getLogLevel() == debug && iterationsOfParameterLearning > 0) {
        st_logDebug("> Learned parameters\n");
        stRPHmmParameters_printParameters(params, stderr);
    }

    // Create HMMs
    st_logInfo("> Creating read partitioning HMMs\n");
    stList *hmms = getRPHmms(profileSequences, referenceNamesToReferencePriors, params);
    st_logDebug("Got %" PRIi64 " hmms before splitting\n", stList_length(hmms));

    // Break up the hmms where the phasing is uncertain
    st_logInfo("> Breaking apart HMMs where the phasing is uncertain\n");

    stList *l = stList_construct3(0, (void (*)(void *))stRPHmm_destruct2);
    while(stList_length(hmms) > 0) {
        stList_appendAll(l, stRPHMM_splitWherePhasingIsUncertain(stList_pop(hmms)));
    }
    hmms = l;
    st_logDebug("\tCreated %d hmms after splitting at uncertain regions of phasing\n", stList_length(hmms));

    // Start VCF generation
    vcfFile *vcfOutFP = vcf_open(vcfOutFile, "w");
    bcf_hdr_t *hdr = writeVcfHeader(vcfOutFP, l, referenceFastaFile);

    // Prep for BAM outputs
    stSet *read1Ids = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, NULL);
    stSet *read2Ids = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, NULL);

    // For each read partitioning HMM
    for(int64_t i=0; i<stList_length(hmms); i++) {
        stRPHmm *hmm = stList_get(hmms, i);

        st_logInfo("> Creating genome fragment for reference sequence: %s, start: %" PRIi64 ", length: %" PRIi64 "\n",
                    hmm->referenceName, hmm->refStart, hmm->refLength);

        // Run the forward-backward algorithm
        stRPHmm_forwardBackward(hmm);

        // Now compute a high probability path through the hmm
        stList *path = stRPHmm_forwardTraceBack(hmm);

        // Compute the genome fragment
        stGenomeFragment *gF = stGenomeFragment_construct(hmm, path);

        // Get the reads which mapped to each path
        stSet *reads1 = stRPHmm_partitionSequencesByStatePath(hmm, path, true);
        stSet *reads2 = stRPHmm_partitionSequencesByStatePath(hmm, path, false);

        addProfileSeqIdsToSet(reads1, read1Ids);
        addProfileSeqIdsToSet(reads2, read2Ids);

        if(st_getLogLevel() == debug) {
            st_logDebug("\nThere are %" PRIi64 " reads covered by the hmm, bipartitioned into sets of %" PRIi64 " and %" PRIi64 " reads\n",
                    stList_length(hmm->profileSeqs), stSet_size(reads1), stSet_size(reads2));

            // Print the similarity between the two imputed haplotypes sequences
            st_logDebug("The haplotypes have an %f identity\n", getIdentityBetweenHaplotypes(gF->haplotypeString1, gF->haplotypeString2, gF->length));
            st_logDebug("\tIdentity excluding indels: %f \n\n", getIdentityBetweenHaplotypesExcludingIndels(gF->haplotypeString1, gF->haplotypeString2, gF->length));

            // Print the base composition of the haplotype sequences
            double *hap1BaseCounts = getHaplotypeBaseComposition(gF->haplotypeString1, gF->length);
            st_logDebug("The base composition of haplotype 1:\n");
            printBaseComposition(stderr, hap1BaseCounts);
            free(hap1BaseCounts);

            double *hap2BaseCounts = getHaplotypeBaseComposition(gF->haplotypeString2, gF->length);
            st_logDebug("The base composition of haplotype 2:\n");
            printBaseComposition(stderr, hap2BaseCounts);
            free(hap2BaseCounts);

            // Print the base composition of the reads
            double *reads1BaseCounts =getExpectedProfileSequenceBaseComposition(reads1);
            st_logDebug("The base composition of reads1 set:\n");
            printBaseComposition(stderr, reads1BaseCounts);
            free(reads1BaseCounts);

            double *reads2BaseCounts =getExpectedProfileSequenceBaseComposition(reads2);
            st_logDebug("The base composition of reads2 set:\n");
            printBaseComposition(stderr, reads2BaseCounts);
            free(reads2BaseCounts);

            st_logDebug("Genome fragment info: refStart = %" PRIi64 ", length = %" PRIi64 "\n", gF->refStart, gF->length);

            // Print some summary stats about the differences between haplotype sequences and the bipartitioned reads
            st_logDebug("hap1 vs. reads1 identity: %f\n",
                    getExpectedIdentity(gF->haplotypeString1, gF->refStart, gF->length, reads1));
            st_logDebug("hap1 vs. reads2 identity: %f\n",
                    getExpectedIdentity(gF->haplotypeString1, gF->refStart, gF->length, reads2));

            st_logDebug("hap2 vs. reads2 identity: %f\n",
                    getExpectedIdentity(gF->haplotypeString2, gF->refStart, gF->length, reads2));
            st_logDebug("hap2 vs. reads1 identity: %f\n",
                    getExpectedIdentity(gF->haplotypeString2, gF->refStart, gF->length, reads1));
        }

        // Write out VCF
        st_logInfo("> Writing out VCF for fragment into file: %s\n", vcfOutFile);

        // Write two vcfs, one using the reference fasta file and one not
        writeVcfFragment(vcfOutFP, hdr, gF, referenceFastaFile, baseMapper, true);

        // Cleanup
        stGenomeFragment_destruct(gF);
        stSet_destruct(reads1);
        stSet_destruct(reads2);
        stList_destruct(path);
    }

    // Cleanup vcf
    vcf_close(vcfOutFP);
    bcf_hdr_destroy(hdr);

    // Write out two BAMs, one for each read partition
    if(samOutBase != NULL) {
        st_logInfo("> Writing out SAM files for each partition into files: %s.1.sam and %s.1.sam\n", samOutBase, samOutBase);
        writeSplitSams(bamInFile, samOutBase, read1Ids, read2Ids);
    }

    // Cleanup
    stHash_destruct(referenceNamesToReferencePriors);
    stList_destruct(hmms);
    stList_destruct(profileSequences);
    stSet_destruct(read1Ids);
    stSet_destruct(read2Ids);

    stBaseMapper_destruct(baseMapper);
    stRPHmmParameters_destruct(params);

    // TODO: only free these if they need to be
//    free(bamInFile);
//    free(samOutBase);
//    free(vcfOutFile);
//    free(referenceName);

    //while(1); // Use this for testing for memory leaks

    return 0;
}

