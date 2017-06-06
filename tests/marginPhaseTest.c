/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */


#include <htslib/vcf.h>
#include <htslib/sam.h>
#include "CuTest.h"
#include "sonLib.h"
#include "stRPHmm.h"
#include "../externalTools/sonLib/C/impl/sonLibListPrivate.h"

void printSequenceStats(stList *profileSequences) {
    /*
     * Print stats about the set of profile sequences.
     */
    int64_t totalLength=0;
    for(int64_t i=0; i<stList_length(profileSequences); i++) {
        stProfileSeq *profileSeq = stList_get(profileSequences, i);
        totalLength += profileSeq->length;
    }
    st_logDebug("Got %" PRIi64 " profile sequences, with total length: %" PRIi64 ", average length: %f\n",
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
    stSet_destructIterator(it);
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

void printAvgIdentityBetweenProfileSequences(stList *profileSequences, int64_t maxSeqs) {
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

    st_logDebug("Avg. pairwise identity between profile sequences: %f measured at %" PRIi64 " overlapping sites\n",
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
    stSet_destructIterator(it);
    return baseCounts;
}





void printBaseComposition(double *baseCounts) {
    /*
     * Print the counts/fraction of each alphabet character.
     */
    double totalCount = 0;
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        totalCount += baseCounts[i];
    }
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        st_logDebug("Base %" PRIi64 " count: %f fraction: %f\n", i, baseCounts[i], baseCounts[i]/totalCount);
    }
}


void printPartitionComposition(int64_t evalPos, stRPHmm *hmm, stSet *reads1, stSet *reads2) {
    double *read1BaseCounts = getProfileSequenceBaseCompositionAtPosition(reads1, evalPos);
    double *read2BaseCounts = getProfileSequenceBaseCompositionAtPosition(reads2, evalPos);
    printColumnAtPosition(hmm, evalPos);
    st_logDebug("\tPartition 1: \n");
    printBaseComposition2(read1BaseCounts);
    st_logDebug("\tPartition 2: \n");
    printBaseComposition2(read2BaseCounts);
    free(read1BaseCounts);
    free(read2BaseCounts);
}


int genotypingTest(char *paramsFile, char *bamFile, char *vcfOutFile,
        char *referenceFile, char *vcfReference, int64_t iterationsOfParameterLearning, bool verbose) {

    // Run margin phase
    char *logString = verbose ? "--logLevel DEBUG" : "--logLevel INFO";
    char *command = stString_print("./marginPhase --bamFile %s --referenceFasta %s %s --params %s --vcfFile %s "
            "--iterationsOfParameterLearning %" PRIi64 " --referenceVcf %s",
            bamFile, referenceFile, logString,
            paramsFile, vcfOutFile, iterationsOfParameterLearning, vcfReference);
    st_logInfo("> Running margin phase on %s\n", bamFile);
    return st_system(command);

    // TODO : Do VCF comparison using VCF eval
}

void test_5kbGenotyping(CuTest *testCase) {

    char *paramsFile = "../tests/params.json";
    char *referenceFile = "../tests/hg19.chr3.9mb.fa";
    char *vcfOutFile = "test_5kb.vcf";
    char *vcfOutFileDiff = "test_5kb_diff.vcf";
    char *samOutBase = "test_100kb";
    double falseNegativeThreshold = 0.20;
    int64_t iterationsOfParameterLearning = 3;
    bool verbose = true;

    char *bamFile = "../tests/NA12878.pb.chr3.5kb.bam";
    // TODO: create vcf specifically for this 5 kb region

    char *vcfReference = "../tests/NA12878.PG.chr3.100kb.0.vcf";


    st_logInfo("Testing haplotype inference on %s\n", bamFile);
    int i = genotypingTest(paramsFile, bamFile, vcfOutFile,
                            referenceFile, vcfReference, iterationsOfParameterLearning, verbose);
    CuAssertTrue(testCase, i == 0);
}

void test_100kbGenotyping(CuTest *testCase) {

    char *paramsFile = "../tests/params.json";
    char *referenceFile = "../tests/hg19.chr3.9mb.fa";
    char *vcfOutFile = "test_100kb.vcf";
    char *vcfOutFileDiff = "test_100kb_diff.vcf";
    char *samOutBase = "test_100kb";
    double falseNegativeThreshold = 0.15;
    int64_t iterationsOfParameterLearning = 1;
    bool verbose = true;

    char *bamFile = "../tests/NA12878.pb.chr3.100kb.1.bam";
//    char *bamFile = "../tests/NA12878.ihs.chr3.100kb.4.bam";
    char *vcfReference = "../tests/NA12878.PG.chr3.100kb.1.vcf";
//    char *vcfReference = "../tests/HG001.GRCh37.chr3.100kb.vcf";
//    char *bamFile = "../tests/NA12878.pb.chr3.2mb.bam";
//    char *vcfReference = "../tests/HG001.GRCh37.chr3.2mb.vcf";

    st_logInfo("Testing haplotype inference on %s\n", bamFile);

//    stGenotypeResults *results = genotypingTest(paramsFile, bamFile, vcfOutFile, vcfOutFileDiff,
//                                                     referenceFile, vcfReference, samOutBase,
//                                                     falseNegativeThreshold, iterationsOfParameterLearning);
//    printGenotypeResults(results);
//
//    free(results);


    int i = genotypingTest(paramsFile, bamFile, vcfOutFile,
                        referenceFile, vcfReference, iterationsOfParameterLearning, verbose);
    CuAssertTrue(testCase, i == 0);

}

void test_multiple100kbGenotyping(CuTest *testCase) {

    st_setLogLevelFromString("info");

    char *paramsFile = "../tests/params.json";
    char *referenceFile = "../tests/hg19.chr3.9mb.fa";
    char *vcfOutFile = "test_100kb.vcf";
    char *vcfOutFileDiff = "test_100kb_diff.vcf";
    char *samOutBase = "test_100kb_0";
    double falseNegativeThreshold = 0.15;
    int64_t iterationsOfParameterLearning = 1;
    bool verbose = false;

    char *bamFile = "../tests/NA12878.pb.chr3.100kb.0.bam";
    char *vcfReference = "../tests/NA12878.PG.chr3.100kb.0.vcf";
    st_logInfo("Testing haplotype inference on %s\n", bamFile);

    int i = genotypingTest(paramsFile, bamFile, vcfOutFile,
                           referenceFile, vcfReference, iterationsOfParameterLearning, verbose);
    CuAssertTrue(testCase, i == 0);


    bamFile = "../tests/NA12878.pb.chr3.100kb.1.bam";
    vcfReference = "../tests/NA12878.PG.chr3.100kb.1.vcf";
    samOutBase = "test_100kb_1";
    st_logInfo("Testing haplotype inference on %s\n", bamFile);

    i = genotypingTest(paramsFile, bamFile, vcfOutFile,
                       referenceFile, vcfReference, iterationsOfParameterLearning, verbose);
    CuAssertTrue(testCase, i == 0);

    bamFile = "../tests/NA12878.pb.chr3.100kb.2.bam";
    vcfReference = "../tests/NA12878.PG.chr3.100kb.2.vcf";
    samOutBase = "test_100kb_2";
    st_logInfo("Testing haplotype inference on %s\n", bamFile);

    i = genotypingTest(paramsFile, bamFile, vcfOutFile,
                       referenceFile, vcfReference, iterationsOfParameterLearning, verbose);
    CuAssertTrue(testCase, i == 0);


    bamFile = "../tests/NA12878.pb.chr3.100kb.3.bam";
    vcfReference = "../tests/NA12878.PG.chr3.100kb.3.vcf";
    samOutBase = "test_100kb_3";
    st_logInfo("Testing haplotype inference on %s\n", bamFile);

    i = genotypingTest(paramsFile, bamFile, vcfOutFile,
                       referenceFile, vcfReference, iterationsOfParameterLearning, verbose);
    CuAssertTrue(testCase, i == 0);

    bamFile = "../tests/NA12878.pb.chr3.100kb.4.bam";
    vcfReference = "../tests/NA12878.PG.chr3.100kb.4.vcf";
    samOutBase = "test_100kb_4";
    st_logInfo("Testing haplotype inference on %s\n", bamFile);

    i = genotypingTest(paramsFile, bamFile, vcfOutFile,
                       referenceFile, vcfReference, iterationsOfParameterLearning, verbose);
    CuAssertTrue(testCase, i == 0);

}




CuSuite *marginPhaseTestSuite(void) {

    CuSuite* suite = CuSuiteNew();

    st_setLogLevelFromString("debug");
//    st_setLogLevelFromString("info");

//    SUITE_ADD_TEST(suite, test_5kbGenotyping);
    SUITE_ADD_TEST(suite, test_100kbGenotyping);
//    SUITE_ADD_TEST(suite, test_multiple100kbGenotyping);

    return suite;
}
