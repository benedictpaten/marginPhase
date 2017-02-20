/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "stRPHmm.h"

#define RANDOM_TEST_NO 100

char *getRandomSequence(int64_t referenceLength) {
    /*
     *
     */
}

char *permuteSequence(char *referenceSeq) {
    /*
     *
     */

}

stProfileSeq *getRandomProfileSeq(char *hapSeq, int64_t readLength) {
    /*
     *
     */

}

static void test_systemTest(CuTest *testCase, int64_t minReferenceSeqNumber, int64_t maxReferenceSeqNumber,
        int64_t minReferenceLength, int64_t maxReferenceLength, int64_t minCoverage, int64_t maxCoverage,
        int64_t minReadLength, int64_t maxReadLength, double posteriorProbabilityThreshold,
        int64_t minColumnDepthToFilter) {
    /*
     * System level test
     *
     * Repeatedly:
     *  Creates reference sequence
     *  Generates two haplotypes
     *  Generates sequences from each haplotype
     *  Creates read HMMs
     *  Check that HMMs represent overlapping sequences as expected
     *  Check forward-backward algorithm is behaving
     *  Create traceBack
     *  Reports how close predicted partition is to actual read partition
     */

    // Print info about test parameters

    for(int64_t test=0; test<RANDOM_TEST_NO; test++) {
        // Creates reference sequences
        // Generates two haplotypes for each reference sequence
        // Generates sequences from each haplotype

        int64_t referenceSeqNumber = st_randomInt(minReferenceSeqNumber, maxReferenceSeqNumber);
        stList *referenceSeqs = stList_construct3(0, free);
        stList *hapSeqs1 = stList_construct3(0, free);
        stList *hapSeqs2 = stList_construct3(0, free);
        stList *profileSeqs1 = stList_construct3(0, (void (*)(void *))stProfileSeq_destruct);
        stList *profileSeqs2 = stList_construct3(0, (void (*)(void *))stProfileSeq_destruct);

        // For each reference sequence
        for(int64_t i=0; i<<referenceSeqNumber; i++) {
            // Generate random reference sequence
            int64_t referenceLength = st_randomInt(minReferenceLength, maxReferenceLength);
            char *referenceSeq = getRandomSequence(referenceLength);

            stList_append(referenceSeqs, referenceSeq);

            // Create haplotype sequences for reference
            char *haplotypeSeq1 = permuteSequence(referenceSeq);
            char *haplotypeSeq2 = permuteSequence(referenceSeq);

            stList_append(hapSeqs1, haplotypeSeq1);
            stList_append(hapSeqs2, haplotypeSeq2);

            // Create read sequences to given coverage
            int64_t coverage = st_randomInt(minCoverage, maxCoverage);
            int64_t totalBasesToSimulate = coverage * referenceLength;
            while(totalBasesToSimulate > 0) {
                // Randomly pick a haplotype sequence to template from
                char *hapSeq = haplotypeSeq1;
                if(st_random() > 0.5) {
                    hapSeq = haplotypeSeq2;
                }
                int64_t readLength = st_randomInt(minReadLength, maxReadLength);
                stList_append(getRandomProfileSeq(hapSeq, readLength));
                totalBasesToSimulate -= readLength;
            }
        }

        stList *profileSeqs = stList_construct();
        stList_appendAll(profileSeqs, profileSeqs1);
        stList_appendAll(profileSeqs, profileSeqs2);


        // Print info about simulated sequences


        // Creates read HMMs
        stList *hmms = getRPHmms(profileSeqs, posteriorProbabilityThreshold,
                                minColumnDepthToFilter, maxCoverage);

        // For each hmm
        for(int64_t i=0; i<stList_length(hmms); i++) {
            stRPHmm *hmm = stList_get(hmms, i);

            // Print HMM info

            // Check it doesn't overlap with another hmm
            for(int64_t j=i+1; j<stList_length(hmms); j++) {
                stRPHmm *hmm2 = stList_get(hmms, j);
                CuAssertTrue(testCase, !stRPHmm_overlapOnReference(hmm, hmm2));
            }

            // Check that each sequence contained in it is in the reference coordinate span
            for(int64_t j=0; stList_length(hmm->profileSeqs); j++) {
                stProfileSeq *profileSeq = stList_get(hmm->profileSeqs, j);

                // Check reference coordinate containment
                CuAssertTrue(testCase, stString_eq(profileSeq->referenceName, hmm->referenceName));
                CuAsseetTrue(testCase, hmm->refStart <= profileSeq->refStart);
                CuAsseetTrue(testCase, hmm->refStart + hmm->refLength >= profileSeq->refStart + profileSeq->length);
            }
        }

        // Check that HMMs represent the overlapping sequences as expected

        // For each sequence, check it is contained in an hmm
        for(int64_t i=0; i<stList_length(profileSeqs); i++) {
            stProfileSeq *profileSeq = stList_get(profileSeqs, i);
            bool containedInHmm = 0;
            for(int64_t j=0; j<stList_length(hmms); j++) {
                stRPHmm *hmm = stList_get(hmms, j);
                // If are on the same reference sequence
                if(stString_eq(profileSeq->referenceName, hmm->referenceName)) {
                    // If overlapping
                    if(hmm->refStart <= profileSeq->refStart && hmm->refStart + hmm->refLength > profileSeq->refStart) {

                        // Must be contained in the hmm
                        CuAssertTrue(testCase, stList_contains(hmm->profileSeqs, profileSeq));

                        // Must not be partially overlapping
                        CuAssertTrue(testCase, hmm->refStart + hmm->refLength >= profileSeq->refStart + profileSeq->length);

                        // Is not contained in another hmm.
                        CuAssertTrue(testCase, !containedInHmm);
                        containedInHmm = 1;
                    }
                }
            }
        }

        // Check columns of each hmm / parameters of hmm
        for(int64_t i=0; i<stList_length(hmms); i++) {
            stRPHmm *hmm = stList_get(hmms, i);
            int64_t columnNumber = 0;
            int64_t maxDepth = 0;
            int64_t refStart = hmm->refStart;

            // Fpr each column
            stRPColumn *column = hmm->firstColumn;
            // Check pointer
            CuAssertTrue(testCase, column->pColumn == NULL);

            while(column != NULL) {
                columnNumber++;

                // Check coordinates
                CuAssertIntEquals(testCase, refStart, column->refStart);
                CuAssertTrue(testCase, column->length > 0);
                refStart += column->length;

                // Check sequences in column
                CuAssertTrue(testCase, column->depth >= 0);
                if(column->depth > maxDepth) {
                    maxDepth = column->depth;
                }
                for(int64_t j=0; j<column->depth; j++) {
                    stProfileSeq *profileSeq = column->seqHeaders[j];
                    // Check is part of hmm
                    CuAssertTrue(testCase, stList_contains(hmm->profileSeqs, profileSeq));
                    // Check column is contained in reference intervsl of the profile sequence
                    CuAssertTrue(testCase, profileSeq->refStart <= column->refStart);
                    CuAssertTrue(testCase, profileSeq->refStart + profileSeq->length >= column->refStart + column->length);
                    // Check the corresponding profile probability array is correct
                    CuAssertPtrEquals(testCase, &(profileSeq->profileProbs[column->refStart-profileSeq->refStart]), column->seqs[j]);
                }

                // Check cells in column
                stRPCell *cell = column->head;
                while(cell != NULL) {
                    // Check that partition is properly specified
                    CuAssertIntEquals(testCase, cell->partition >> column->depth, 0);

                    cell = cell->nCell;
                }

                if(column->nColumn == NULL) {
                    CuAssertPtrEquals(testCase, hmm->lastColumn, column);
                    break;
                }

                // Check merge column
                stRPMergeColumn *mColumn = column->nColumn;

                // Check back link
                CuAssertPtrEquals(testCase, column, mColumn->pColumn);

                for(int64_t j=0; j<column->depth; j++) {
                    stProfileSeq *profileSeq = column->seqHeaders[j];
                    // If the sequence ends at the boundary of the column check it is masked out
                    if(profileSeq->refStart + profileSeq->length == column->refStart + column->length) {
                        CuAssertIntEquals(testCase, (mColumn->maskFrom >> j) & 1, 0);
                    }
                    // Otherwise check it is not masked out
                    else {
                        CuAssertIntEquals(testCase, (mColumn->maskFrom >> j) & 1, 1);
                    }
                }

                column = column->nColumn->nColumn;

                // Check back link pointer
                CuAssertPtrEquals(testCase, mColumn, column->pColumn);

                for(int64_t j=0; j<column->depth; j++) {
                    stProfileSeq *profileSeq = column->seqHeaders[j];
                    // If the sequence starts at the start of the column check it is masked out
                    if(profileSeq->refStart == column->refStart) {
                        CuAssertIntEquals(testCase, (mColumn->maskTo >> j) & 1, 0);
                    }
                    // Otherwise check it is not masked out
                    else {
                        CuAssertIntEquals(testCase, (mColumn->maskTo >> j) & 1, 1);
                    }
                }

                // Check merge cells are same in both the from and to sets
                stList *mCellsFrom = stHash_getValues(mColumn->mergeCellsFrom);
                stList *mCellsTo = stHash_getValues(mColumn->mergeCellsTo);
                stSet *mCellsFromSet = stList_getSet(mCellsFrom);
                stSet *mCellsToSet = stList_getSet(mCellsTo);
                CuAssertTrue(testCase, stSet_equals(mCellsFromSet, mCellsToSet));
                stSet_destruct(mCellsFromSet);
                stSet_destruct(mCellsToSet);
                stList_destruct(mCellsFrom);
                stList_destruct(mCellsTo);

                // check merge cells
                stRPMergeCell *mCell;
                stHashIterator *it = stHash_getIterator(mColumn->mergeCellsFrom);
                while((mCell = stHash_getNext(mCell)) != NULL) {
                    // Check partitions
                    CuAssertTrue(testCase, mCell->fromPartition & mColumn->maskFrom == mCell->fromPartition);
                    CuAssertTrue(testCase, mCell->toPartition & mColumn->maskTo == mCell->toPartition);
                }
                stHash_destructIterator(it);

            }

            // Check recorded column number and max depth match up with what we observed
            CuAssertIntEquals(testCase, columnNumber, hmm->columnNumber);
            CuAssertIntEquals(testCase, maxDepth, hmm->maxDepth);
            CuAssertTrue(testCase, hmm-maxDepth <= MAX_READ_PARTITIONING_DEPTH);
        }

        // Check forward-backward algorithm is behaving
        for(int64_t i=0; i<stList_length(hmms); i++) {
            stRPHmm *hmm = stList_get(hmms, i);

            // Run the forward and backward algorithms
            stRPHmm_forward(hmm);
            stRPHmm_backward(hmm);

            // Check the forward/backward/total probabilities
            double forward = hmm->forwardProbability;
            double backward = hmm->backwardProbability;
            double totalProb = hmm->forwardProbability + hmm->backwardProbability;

            // Check each column / cell
            stRPColumn *column = hmm->firstColumn;
            while(1) {
                // Check column
                CuAssertDblEquals(testCase, forward, column->forwardProb, 0.01);
                CuAssertDblEquals(testCase, backward, column->backwardProb, 0.01);

                // Check posterior probabilities
                stRPCell *cell = column->head;
                double totalProb = 0.0;
                while(cell != NULL) {
                    double posteriorProb = stRPCell_posteriorProb(cell, column);
                    CuAssertTrue(testCase, posteriorProb > 0.0);
                    CuAssertTrue(testCase, posteriorProb <= 1.0);
                    totalProb += posteriorProb;
                    cell = cell->nCell;
                }

                CuAssertDblEquals(testCase, 1.0, totalProb, 0.01);

                if(column->nColumn == NULL) {
                    break;
                }

                // Check merge column
                stRPMergeColumn *mColumn = column->nColumn;

                // Check posterior probabilities of merge cells
                stRPMergeCell *mCell;
                stHashIterator *it = stHash_getIterator(mColumn->mergeCellsFrom);
                totalProb = 0.0;
                double totalProb = 0.0;
                while((mCell = stHash_getNext(it)) != NULL) {
                    double posteriorProb = stRPMergeCell_posteriorProb(mCell);
                    CuAssertTrue(testCase, posteriorProb > 0.0);
                    CuAssertTrue(testCase, posteriorProb <= 1.0);
                    totalProb += posteriorProb;
                }
                CuAssertDblEquals(testCase, 1.0, totalProb, 0.01);

                column = mColumn->nColumn;
            }
        }

        // Create traceBacks
        for(int64_t i=0; i<stList_length(hmms); i++) {
            stRPHmm *hmm = stList_get(hmms, i);

            stList *traceBackPath = stRPHmm_forwardTraceBack(hmm);

            // Check each cell in traceback


            stSet *predictedHaplotype1Seqs = stRPHmm_partitionSequencesByStatePath(hmm, traceBackPath);

            // Reports how close predicted partition is to actual read partition
            stSet *actualHaplotype1Seqs = stList_getSet(profileSeqs1);

            stSet *hapSeq1Overlap = stSet_getIntersection(predictedHaplotype1Seqs, actualHaplotype1Seqs);
            double precision = ((float)stSet_size(predictedHaplotype1Seqs) - stSet_size(hapSeq1Overlap))/stSet_size(predictedHaplotype1Seqs);
            double recall = ((float)stSet_size(hapSeq1Overlap))/stSet_size(actualHaplotype1Seqs);

            fprintf(STD_ERR, " There were %" PRIi64 " hap1 seqs and %" PRIi64 "hap2 seqs, got precision: %f and recall: %f", stList_length(profileSeqs1),
                    stList_length(profileSeqs2), (float)precision, (float)recall);
        }

    }
}

static void test_systemSingleReferenceFullLengthReads(CuTest *testCase) {
    int64_t minReferenceSeqNumber = 1;
    int64_t maxReferenceSeqNumber = 1;
    int64_t minReferenceLength = 1000;
    int64_t maxReferenceLength = 1000;
    int64_t minCoverage = 20;
    int64_t maxCoverage = 20;
    int64_t minReadLength = 1000;
    int64_t maxReadLength = 1000;
    double posteriorProbabilityThreshold = 0.1;
    int64_t minColumnDepthToFilter = 10;

    test_systemTest(testCase, minReferenceSeqNumber, maxReferenceSeqNumber,
            minReferenceLength, maxReferenceLength, minCoverage, maxCoverage,
            minReadLength, maxReadLength, posteriorProbabilityThreshold,
            minColumnDepthToFilter);
}

static void test_systemSingleReferenceFixedLengthReads(CuTest *testCase) {
    int64_t minReferenceSeqNumber = 1;
    int64_t maxReferenceSeqNumber = 1;
    int64_t minReferenceLength = 1000;
    int64_t maxReferenceLength = 1000;
    int64_t minCoverage = 20;
    int64_t maxCoverage = 20;
    int64_t minReadLength = 100;
    int64_t maxReadLength = 100;
    double posteriorProbabilityThreshold = 0.1;
    int64_t minColumnDepthToFilter = 10;

    test_systemTest(testCase, minReferenceSeqNumber, maxReferenceSeqNumber,
            minReferenceLength, maxReferenceLength, minCoverage, maxCoverage,
            minReadLength, maxReadLength, posteriorProbabilityThreshold,
            minColumnDepthToFilter);
}

static void test_systemSingleReference(CuTest *testCase) {
    int64_t minReferenceSeqNumber = 1;
    int64_t maxReferenceSeqNumber = 1;
    int64_t minReferenceLength = 1000;
    int64_t maxReferenceLength = 1000;
    int64_t minCoverage = 20;
    int64_t maxCoverage = 20;
    int64_t minReadLength = 10;
    int64_t maxReadLength = 300;
    double posteriorProbabilityThreshold = 0.1;
    int64_t minColumnDepthToFilter = 10;

    test_systemTest(testCase, minReferenceSeqNumber, maxReferenceSeqNumber,
            minReferenceLength, maxReferenceLength, minCoverage, maxCoverage,
            minReadLength, maxReadLength, posteriorProbabilityThreshold,
            minColumnDepthToFilter);
}

static void test_systemMultipleReferences(CuTest *testCase) {
    int64_t minReferenceSeqNumber = 2;
    int64_t maxReferenceSeqNumber = 5;
    int64_t minReferenceLength = 1000;
    int64_t maxReferenceLength = 2000;
    int64_t minCoverage = 5;
    int64_t maxCoverage = 20;
    int64_t minReadLength = 10;
    int64_t maxReadLength = 300;
    double posteriorProbabilityThreshold = 0.1;
    int64_t minColumnDepthToFilter = 10;

    test_systemTest(testCase, minReferenceSeqNumber, maxReferenceSeqNumber,
            minReferenceLength, maxReferenceLength, minCoverage, maxCoverage,
            minReadLength, maxReadLength, posteriorProbabilityThreshold,
            minColumnDepthToFilter);
}


CuSuite *stRPHmmTestSuite(void) {
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, test_systemSingleReferenceFullLengthReads);
    SUITE_ADD_TEST(suite, test_systemSingleReferenceFixedLengthReads);
    SUITE_ADD_TEST(suite, test_systemSingleReference);
    SUITE_ADD_TEST(suite, test_systemMultipleReferences);

    return suite;
}
