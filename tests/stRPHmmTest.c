/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "stRPHmm.h"
#include <math.h>
#include <time.h>

#define RANDOM_TEST_NO 2

char getRandomBase() {
    /*
     * Returns an ascii character starting from ascii symbol 48: '0', '1', '2', ... alphabetsize
     */
    return st_randomInt(FIRST_ALPHABET_CHAR, FIRST_ALPHABET_CHAR+ALPHABET_SIZE);
}

char *getRandomSequence(int64_t referenceLength) {
    /*
     * Creates a random sequence of form [ACTG]*referenceLength
     */
    char *randomSequence = st_malloc(sizeof(char) * (referenceLength+1));
    for(int64_t i=0; i<referenceLength; i++) {
        randomSequence[i] = getRandomBase();
    }
    randomSequence[referenceLength] = '\0';

    return randomSequence;
}

char *permuteSequence(char *referenceSeq, double hetRate) {
    /*
     * Takes a random sequence and returns a copy of it, permuting randomly each position with rate hetRate.
     */
    referenceSeq = stString_copy(referenceSeq);
    int64_t strLength = strlen(referenceSeq);
    for(int64_t i=0; i<strLength; i++) {
        if(st_random() < hetRate) {
            referenceSeq[i] = getRandomBase();
        }
    }
    return referenceSeq;
}

stProfileSeq *getRandomProfileSeq(char *referenceName, char *hapSeq, int64_t hapLength,
        int64_t readLength, double readErrorRate) {
    /*
     * Creates a random "read" of the given length from a random sub-interval
     * of the input sequence with error rate errRate.
     */
    assert(hapLength-readLength >= 0);
    int64_t start = st_randomInt(0, hapLength-readLength+1);

    stProfileSeq *pSeq = stProfileSeq_constructEmptyProfile(referenceName, " ", start,
            readLength);

    for(int64_t i=0; i<readLength; i++) {
        // Haplotype base or error at random
        char b = st_random() < readErrorRate ? getRandomBase() : hapSeq[start+i];
        assert(b - FIRST_ALPHABET_CHAR >= 0);
        assert(b - FIRST_ALPHABET_CHAR < ALPHABET_SIZE);
        // Fill in the profile probabilities according to the chosen base
        pSeq->profileProbs[i*ALPHABET_SIZE + b - FIRST_ALPHABET_CHAR] = ALPHABET_MAX_PROB;
    }

    return pSeq;
}

static stRPHmmParameters *getHmmParams(int64_t maxPartitionsInAColumn,
        double hetRate, double readErrorRate,
        bool maxNotSumTransitions, int64_t minReadCoverageToSupportPhasingBetweenHeterozygousSites) {
    // Substitution models
    uint16_t *hetSubModel = st_calloc(ALPHABET_SIZE*ALPHABET_SIZE, sizeof(uint16_t));
    uint16_t *readErrorSubModel = st_calloc(ALPHABET_SIZE*ALPHABET_SIZE, sizeof(uint16_t));
    double *readErrorSubModelSlow = st_calloc(ALPHABET_SIZE*ALPHABET_SIZE, sizeof(double));
    double *hetSubModelSlow = st_calloc(ALPHABET_SIZE*ALPHABET_SIZE, sizeof(double));

    // Fill in substitition matrix with simple symmetric probs matching
    // the read error rate
    assert(readErrorRate <= 1.0);
    assert(readErrorRate >= 0.0);
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        for(int64_t j=0; j<ALPHABET_SIZE; j++) {
            setSubstitutionProb(readErrorSubModel, readErrorSubModelSlow, i, j, i==j ?
                    1.0-readErrorRate : readErrorRate/(ALPHABET_SIZE-1));
            setSubstitutionProb(hetSubModel, hetSubModelSlow, i, j, i==j ?
                    1.0-hetRate : hetRate/(ALPHABET_SIZE-1));
        }
    }

    return stRPHmmParameters_construct(hetSubModel, hetSubModelSlow,
            readErrorSubModel, readErrorSubModelSlow,
            maxNotSumTransitions, maxPartitionsInAColumn, MAX_READ_PARTITIONING_DEPTH,
            minReadCoverageToSupportPhasingBetweenHeterozygousSites);
}

static void simulateReads(stList *referenceSeqs, stList *hapSeqs1, stList *hapSeqs2,
        stList *profileSeqs1, stList *profileSeqs2,
        int64_t minReferenceSeqNumber, int64_t maxReferenceSeqNumber,
        int64_t minReferenceLength, int64_t maxReferenceLength,
        int64_t minCoverage, int64_t maxCoverage,
        int64_t minReadLength, int64_t maxReadLength,
        double hetRate, double readErrorRate) {
    /*
     * Simulate reference sequence, haplotypes and derived reads, represented as profile
     * sequences, placing the results in the argument lists.
     *
     *
     */
    int64_t referenceSeqNumber = st_randomInt(minReferenceSeqNumber, maxReferenceSeqNumber+1);

    // For each reference sequence
    for(int64_t i=0; i<referenceSeqNumber; i++) {
        // Generate random reference sequence
        int64_t referenceLength = st_randomInt(minReferenceLength, maxReferenceLength+1);
        char *referenceSeq = getRandomSequence(referenceLength);
        // Reference name
        char *referenceName = stString_print("Reference_%" PRIi64 "", i);

        stList_append(referenceSeqs, referenceSeq);

        // Create haplotype sequences for reference
        char *haplotypeSeq1 = permuteSequence(referenceSeq, hetRate);
        char *haplotypeSeq2 = permuteSequence(referenceSeq, hetRate);

        stList_append(hapSeqs1, haplotypeSeq1);
        stList_append(hapSeqs2, haplotypeSeq2);

        // Print info about simulated sequences
        /*fprintf(stderr, "Ref seq: %s\n", referenceSeq);
        fprintf(stderr, "  Hap 1: %s\n", haplotypeSeq1);
        fprintf(stderr, "  Hap 2: %s\n", haplotypeSeq2);*/

        int64_t diffs = 0;
        for(int64_t j=0; j<strlen(haplotypeSeq1); j++) {
            if(haplotypeSeq1[j] != haplotypeSeq2[j]) {
                diffs++;
            }
        }
        fprintf(stderr, " Differences between hap1 and hap2: %" PRIi64 " length %"
                PRIi64 " het rate per haplotype base: %f actual diff rate: %f\n", diffs,
                (int64_t)strlen(haplotypeSeq1), hetRate, (double)diffs/strlen(haplotypeSeq1));

        // Create read sequences to given coverage
        int64_t coverage = st_randomInt(minCoverage, maxCoverage+1);
        int64_t totalBasesToSimulate = coverage * referenceLength;
//        fprintf(stderr, "Total bases to simulate: %" PRIi64 " for coverage: %" PRIi64 "\n", totalBasesToSimulate, coverage);
        while(totalBasesToSimulate > 0) {
            // Randomly pick a haplotype sequence to template from
            char *hapSeq = haplotypeSeq1;
            stList *readSeqs = profileSeqs1;
            if(st_random() > 0.5) {
                hapSeq = haplotypeSeq2;
                readSeqs = profileSeqs2;
            }
            int64_t readLength = st_randomInt(minReadLength, maxReadLength+1);
            stProfileSeq *pSeq = getRandomProfileSeq(referenceName, hapSeq,
                    referenceLength, readLength, readErrorRate);
            stList_append(readSeqs, pSeq);
            totalBasesToSimulate -= readLength;
            //fprintf(stderr, "Simulating read from haplotype: %s\n", hapSeq);
            //stProfileSeq_print(pSeq, stderr, 1);
        }

        // Cleanup
        free(referenceName);
    }
}

static void test_systemTest(CuTest *testCase, int64_t minReferenceSeqNumber, int64_t maxReferenceSeqNumber,
        int64_t minReferenceLength, int64_t maxReferenceLength, int64_t minCoverage, int64_t maxCoverage,
        int64_t minReadLength, int64_t maxReadLength,
        int64_t maxPartitionsInAColumn, double hetRate, double readErrorRate,
        bool maxNotSumTransitions, bool splitHmmsWherePhasingUncertain,
        int64_t minReadCoverageToSupportPhasingBetweenHeterozygousSites,
        bool printHmm) {
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
    fprintf(stderr, " System test parameters:\n"
            "\tminReferenceSequenceNumber: %" PRIi64 "\n"
            "\tmaxReferenceSequenceNumber: %" PRIi64 "\n"
            "\tminReferenceLength: %" PRIi64 "\n"
            "\tmaxReferenceLength: %" PRIi64 "\n"
            "\tminCoverage: %" PRIi64 "\n"
            "\tmaxCoverage: %" PRIi64 "\n"
            "\tminReadLength: %" PRIi64 "\n"
            "\tmaxReadLength: %" PRIi64 "\n"
            "\tmaxPartitionsInAColumn: %" PRIi64 "\n"
            "\thetRate: %f\n"
            "\treadErrorRate: %f\n"
            "\tmaxNotSumTransitions: %i\n"
            "\tsplitHmmsWherePhasingUncertain: %i\n",
            minReferenceSeqNumber, maxReferenceSeqNumber,
            minReferenceLength, maxReferenceLength, minCoverage, maxCoverage,
            minReadLength, maxReadLength,
            maxPartitionsInAColumn, (float)hetRate, (float)readErrorRate,
            maxNotSumTransitions, splitHmmsWherePhasingUncertain);

    int64_t totalProfile1SeqsOverAllTests = 0;
    int64_t totalProfile2SeqsOverAllTests = 0;
    int64_t totalPartitionErrorsOverAllTests = 0;

    time_t startTime = time(NULL);

    for(int64_t test=0; test<RANDOM_TEST_NO; test++) {

        fprintf(stderr, "Starting test iteration: #%" PRIi64 "\n", test);

        stRPHmmParameters *params = getHmmParams(maxPartitionsInAColumn,
                hetRate, readErrorRate, maxNotSumTransitions,
                minReadCoverageToSupportPhasingBetweenHeterozygousSites);

        stList *referenceSeqs = stList_construct3(0, free);
        stList *hapSeqs1 = stList_construct3(0, free);
        stList *hapSeqs2 = stList_construct3(0, free);
        stList *profileSeqs1 = stList_construct3(0, (void (*)(void *))stProfileSeq_destruct);
        stList *profileSeqs2 = stList_construct3(0, (void (*)(void *))stProfileSeq_destruct);

        // Creates reference sequences
        // Generates two haplotypes for each reference sequence
        // Generates profile sequences from each haplotype
        simulateReads(referenceSeqs, hapSeqs1, hapSeqs2,
                        profileSeqs1, profileSeqs2,
                        minReferenceSeqNumber, maxReferenceSeqNumber,
                        minReferenceLength, maxReferenceLength,
                        minCoverage, maxCoverage,
                        minReadLength, maxReadLength,
                        hetRate, readErrorRate);

        stList *profileSeqs = stList_construct();
        stList_appendAll(profileSeqs, profileSeqs1);
        stList_appendAll(profileSeqs, profileSeqs2);
        stList_shuffle(profileSeqs); // Ensure we don't have the reads already partitioned!

        // Set representations of the profile sequences
        stSet *profileSeqs1Set = stList_getSet(profileSeqs1);
        stSet *profileSeqs2Set = stList_getSet(profileSeqs2);

        fprintf(stderr, "Running get hmms with %" PRIi64 " profile sequences \n", stList_length(profileSeqs));

        // Creates read HMMs
        stList *hmms = getRPHmms(profileSeqs, params);

        // Split hmms where phasing is uncertain
        if(splitHmmsWherePhasingUncertain) {
            stList *splitHmms = stList_construct3(0,  (void (*)(void *))stRPHmm_destruct2);
            while(stList_length(hmms) > 0) {
                stList *l = stRPHMM_splitWherePhasingIsUncertain(stList_pop(hmms));
                stList_appendAll(splitHmms, l);
                stList_setDestructor(l, NULL);
                stList_destruct(l);
            }
            stList_destruct(hmms);
            hmms = splitHmms;
        }

        fprintf(stderr, "Got %" PRIi64 " hmms\n", stList_length(hmms));

        // For each hmm
        for(int64_t i=0; i<stList_length(hmms); i++) {
            stRPHmm *hmm = stList_get(hmms, i);

            // Check it doesn't overlap with another hmm
            for(int64_t j=i+1; j<stList_length(hmms); j++) {
                stRPHmm *hmm2 = stList_get(hmms, j);
                CuAssertTrue(testCase, !stRPHmm_overlapOnReference(hmm, hmm2));
            }

            // Check that each sequence contained in it is in the reference coordinate span
            for(int64_t j=0; j<stList_length(hmm->profileSeqs); j++) {
                stProfileSeq *profileSeq = stList_get(hmm->profileSeqs, j);
                // Check reference coordinate containment
                CuAssertTrue(testCase, stString_eq(profileSeq->referenceName, hmm->referenceName));
                if(!splitHmmsWherePhasingUncertain) { // The profile sequence is only guaranteed to be wholly contained
                    // in the hmm if the hmm was not split at points where the phasing is uncertain
                    CuAssertTrue(testCase, hmm->refStart <= profileSeq->refStart);
                    CuAssertTrue(testCase, hmm->refStart + hmm->refLength >= profileSeq->refStart + profileSeq->length);
                }
                else {
                    // Must overlap
                    assert(hmm->refStart + hmm->refLength > profileSeq->refStart);
                    assert(profileSeq->refStart + profileSeq->length > hmm->refStart);
                }
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

                        if(!splitHmmsWherePhasingUncertain) { // The profile sequence is only guaranteed to be wholly contained
                                            // in the hmm if the hmm was not split at points where the phasing is uncertain
                            // Must not be partially overlapping
                            CuAssertTrue(testCase, hmm->refStart + hmm->refLength >= profileSeq->refStart + profileSeq->length);
                        }
                        else {
                            // Must overlap
                            assert(hmm->refStart + hmm->refLength > profileSeq->refStart);
                            assert(profileSeq->refStart + profileSeq->length > hmm->refStart);
                        }

                        // Is not contained in another hmm.
                        CuAssertTrue(testCase, !containedInHmm);
                        containedInHmm = 1;
                    }
                }
            }
            CuAssertTrue(testCase, containedInHmm);
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
                    CuAssertPtrEquals(testCase, &(profileSeq->profileProbs[(column->refStart-profileSeq->refStart) * ALPHABET_SIZE]), column->seqs[j]);
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
                // Clean up
                stSet_destruct(mCellsFromSet);
                stSet_destruct(mCellsToSet);
                stList_destruct(mCellsFrom);
                stList_destruct(mCellsTo);

                // Check merge cells
                stRPMergeCell *mCell;
                stHashIterator *it = stHash_getIterator(mColumn->mergeCellsFrom);
                while((mCell = stHash_getNext(it)) != NULL) {
                    // Check partitions
                    CuAssertTrue(testCase, (mCell->fromPartition & mColumn->maskFrom) == mCell->fromPartition);
                    CuAssertTrue(testCase, (mCell->toPartition & mColumn->maskTo) == mCell->toPartition);
                }
                stHash_destructIterator(it);

            }

            // Check recorded column number and max depth match up with what we observed
            CuAssertIntEquals(testCase, columnNumber, hmm->columnNumber);
            CuAssertIntEquals(testCase, maxDepth, hmm->maxDepth);
            CuAssertTrue(testCase, hmm->maxDepth <= MAX_READ_PARTITIONING_DEPTH);
        }

        // Check forward-backward algorithm is behaving
        for(int64_t i=0; i<stList_length(hmms); i++) {
            stRPHmm *hmm = stList_get(hmms, i);

            // Run the forward and backward algorithms
            stRPHmm_forwardBackward(hmm);

            // Print HMM info, which will contain forward and backward probs
            if(printHmm) {
                stRPHmm_print(hmm, stderr, 1, 1);
            }

            // Check the forward and backward probs are close
            CuAssertDblEquals(testCase, hmm->forwardLogProb, hmm->backwardLogProb, 0.1);

            // Check each column / cell
            stRPColumn *column = hmm->firstColumn;
            while(1) {
                // Check column total probability is same as total probability in model
                CuAssertDblEquals(testCase, hmm->forwardLogProb,
                        column->totalLogProb, 0.1);

                // Check posterior probabilities
                stRPCell *cell = column->head;
                double totalProb = 0.0;
                while(cell != NULL) {
                    double posteriorProb = stRPCell_posteriorProb(cell, column);
                    CuAssertTrue(testCase, posteriorProb >= 0.0);
                    CuAssertTrue(testCase, posteriorProb <= 1.0);
                    totalProb += posteriorProb;
                    cell = cell->nCell;
                }

                if(!maxNotSumTransitions) {
                    CuAssertDblEquals(testCase, 1.0, totalProb, 0.1);
                }
                if(column->nColumn == NULL) {
                    break;
                }

                // Check merge column
                stRPMergeColumn *mColumn = column->nColumn;

                // Check posterior probabilities of merge cells
                stRPMergeCell *mCell;
                stHashIterator *it = stHash_getIterator(mColumn->mergeCellsFrom);
                totalProb = 0.0;
                while((mCell = stHash_getNext(it)) != NULL) {
                    double posteriorProb = stRPMergeCell_posteriorProb(mCell, mColumn);
                    CuAssertTrue(testCase, posteriorProb >= 0.0);
                    CuAssertTrue(testCase, posteriorProb <= 1.0);
                    totalProb += posteriorProb;
                }
                stHash_destructIterator(it);
                if(!maxNotSumTransitions) {
                    CuAssertDblEquals(testCase, 1.0, totalProb, 0.1);
                }
                column = mColumn->nColumn;
            }
        }

        // Create traceBacks
        int64_t totalPartitionError = 0; // The number of switch operations required to make the
        // partitions predicted by the models perfect

        for(int64_t i=0; i<stList_length(hmms); i++) {
            stRPHmm *hmm = stList_get(hmms, i);

            stList *traceBackPath = stRPHmm_forwardTraceBack(hmm);

            // Check each cell in traceback
            stRPColumn *column = hmm->firstColumn;
            for(int64_t j=0; j<stList_length(traceBackPath); j++) {
                stRPCell *cell = stList_get(traceBackPath, j);

                // Must belong to the given column
                stRPCell *cell2 = column->head;
                bool inColumn = 0;
                while((cell2 != NULL)) {
                    if(cell == cell2) {
                        inColumn = 1;
                        break;
                    }
                    cell2 = cell2->nCell;
                }
                CuAssertTrue(testCase, inColumn);

                // Must be compatible with previous cell (i.e. point to the same merge cell)
                if(j > 0) {
                    stRPCell *pCell = stList_get(traceBackPath, j-1);
                    stRPMergeCell *mCell = stRPMergeColumn_getPreviousMergeCell(cell, column->pColumn);
                    stRPMergeCell *mCell2 = stRPMergeColumn_getNextMergeCell(pCell, column->pColumn);
                    CuAssertPtrEquals(testCase, mCell, mCell2);
                }

                // Get next column
                if(j+1 < stList_length(traceBackPath)) {

                    CuAssertTrue(testCase, column->nColumn != NULL);
                    column = column->nColumn->nColumn;
                }
                else { // Verify we've reached the expected end
                    CuAssertTrue(testCase, column == hmm->lastColumn);
                }
            }

            stSet *profileSeqsPartition1 = stRPHmm_partitionSequencesByStatePath(hmm, traceBackPath, 1);
            stSet *profileSeqsPartition2 = stRPHmm_partitionSequencesByStatePath(hmm, traceBackPath, 0);

            // Check all the sequences accounted for
            CuAssertIntEquals(testCase, stList_length(hmm->profileSeqs),
                    stSet_size(profileSeqsPartition1)+stSet_size(profileSeqsPartition2));

            /*
             * Comparing a given HMMs partition to the true read partition there are four set
             * cardinalities we care about:
             *
             * p11 - the cardinality of the intersection of the first subpartition of the hmm
             * and the profile sequences for the first haplotype sequence
             *
             * p12 - the cardinality of the intersection of the first subpartition of the hmm
             * and the profile sequences for the second haplotype sequence
             *
             * p21 - the cardinality of the intersection of the second subpartition of the hmm
             * and the profile sequences for the first haplotype sequence
             *
             * p22 - the cardinality of the intersection of the second subpartition of the hmm
             * and the profile sequences for the second haplotype sequence
             *
             * Given these quantities the partitionErrors is the minimum number of "movements" of elements
             * in the hmm partition to ensure that the HMM partition is a subset of the true partition.
             */

            stSet *x = stSet_getIntersection(profileSeqsPartition1, profileSeqs1Set);
            int64_t p11 = stSet_size(x);
            stSet_destruct(x);

            int64_t p12 = stSet_size(profileSeqsPartition1) - p11;
            assert(p12 >= 0);

            x = stSet_getIntersection(profileSeqsPartition2, profileSeqs2Set);
            int64_t p21 = stSet_size(x);
            stSet_destruct(x);

            int64_t p22 = stSet_size(profileSeqsPartition2) - p21;
            assert(p22 >= 0);

            int64_t partitionErrors = p12 + p22 < p11 + p21 ? p12 + p22 : p11 + p21;

            fprintf(stderr, "For HMM %" PRIi64 " there were %" PRIi64 " hap1 seqs and %" PRIi64 " hap2 seqs, "
                    "got partition errors: %" PRIi64 "\n", i,
                    stSet_size(profileSeqsPartition1),
                    stSet_size(profileSeqsPartition2), partitionErrors);

            /*for(int64_t k=0; k<stList_length(hapSeqs1); k++) {
                fprintf(stderr, "Hap1: %s\n", (char *)stList_get(hapSeqs1, k));
            }
            for(int64_t k=0; k<stList_length(hapSeqs2); k++) {
                fprintf(stderr, "Hap2: %s\n", (char *)stList_get(hapSeqs2, k));
            }
            printPartition(stderr, profileSeqsPartition1, profileSeqsPartition2);*/

            totalPartitionError += partitionErrors;

            totalProfile1SeqsOverAllTests += stList_length(profileSeqs1);
            totalProfile2SeqsOverAllTests += stList_length(profileSeqs2);
            totalPartitionErrorsOverAllTests += partitionErrors;

            /*
             * Test genome fragment code
             */
            stGenomeFragment *gF = stGenomeFragment_construct(hmm, traceBackPath);

            // Check coordinates
            CuAssertStrEquals(testCase, hmm->referenceName, gF->referenceName);
            CuAssertIntEquals(testCase, hmm->refStart, gF->refStart);
            CuAssertIntEquals(testCase, hmm->refLength, gF->length);

            // Get the haplotype sequences
            stList *tokens = stString_splitByString(hmm->referenceName, "_");
            assert(stList_length(tokens) == 2);
            assert(stString_eq(stList_get(tokens, 0), "Reference"));
            int64_t refSeqIndex = stSafeStrToInt64(stList_peek(tokens));
            stList_destruct(tokens);
            char *hap1Seq = stList_get(hapSeqs1, refSeqIndex);
            char *hap2Seq = stList_get(hapSeqs2, refSeqIndex);

            int64_t correctGenotypes = 0;
            int64_t totalHets = 0;
            int64_t correctHets = 0;
            double probsOfCorrectGenotypes = 0.0;
            double probsOfIncorrectGenotypes = 0.0;
            int64_t hap1ToPredictedHap1Diffs = 0;
            int64_t hap1ToPredictedHap2Diffs = 0;
            int64_t hap2ToPredictedHap1Diffs = 0;
            int64_t hap2ToPredictedHap2Diffs = 0;
            int64_t hap1ToPredictedHap1HetDiffs = 0;
            int64_t hap1ToPredictedHap2HetDiffs = 0;
            int64_t hap2ToPredictedHap1HetDiffs = 0;
            int64_t hap2ToPredictedHap2HetDiffs = 0;
            // For each character
            for(int64_t j=0; j<gF->length; j++) {
                // Check genotype
                CuAssertTrue(testCase, gF->genotypeString[j] <= ALPHABET_SIZE * ALPHABET_SIZE);
                char hap1Char = hap1Seq[j + gF->refStart] - FIRST_ALPHABET_CHAR;
                char hap2Char = hap2Seq[j + gF->refStart] - FIRST_ALPHABET_CHAR;
                uint64_t trueGenotype = hap1Char < hap2Char ? hap1Char * ALPHABET_SIZE + hap2Char : hap2Char * ALPHABET_SIZE + hap1Char;

                totalHets += hap1Char != hap2Char ? 1 : 0;

                if(gF->genotypeString[j] == trueGenotype) {
                    correctGenotypes++;
                    correctHets += hap1Char != hap2Char ? 1 : 0;
                    probsOfCorrectGenotypes += gF->genotypeProbs[j];
                }
                else {
                    probsOfIncorrectGenotypes += gF->genotypeProbs[j];
                }

                // Check genotype posterior probability
                CuAssertTrue(testCase, gF->genotypeProbs[j] >= 0.0);
                CuAssertTrue(testCase, gF->genotypeProbs[j] <= 1.0);

                // Check haplotypes
                CuAssertTrue(testCase, gF->haplotypeString1[j] <= ALPHABET_SIZE);
                CuAssertTrue(testCase, gF->haplotypeString2[j] <= ALPHABET_SIZE);
                if(gF->haplotypeString1[j] < gF->haplotypeString2[j]) {
                    CuAssertTrue(testCase, gF->haplotypeString1[j] * ALPHABET_SIZE + gF->haplotypeString2[j] == gF->genotypeString[j]);
                }
                else {
                    CuAssertTrue(testCase, gF->haplotypeString2[j] * ALPHABET_SIZE + gF->haplotypeString1[j] == gF->genotypeString[j]);
                }

                hap1ToPredictedHap1Diffs += gF->haplotypeString1[j] == hap1Char ? 0 : 1;
                hap1ToPredictedHap2Diffs += gF->haplotypeString2[j] == hap1Char ? 0 : 1;
                hap2ToPredictedHap1Diffs += gF->haplotypeString1[j] == hap2Char ? 0 : 1;
                hap2ToPredictedHap2Diffs += gF->haplotypeString2[j] == hap2Char ? 0 : 1;

                if(hap1Char != hap2Char) {
                    hap1ToPredictedHap1HetDiffs += gF->haplotypeString1[j] == hap1Char ? 0 : 1;
                    hap1ToPredictedHap2HetDiffs += gF->haplotypeString2[j] == hap1Char ? 0 : 1;
                    hap2ToPredictedHap1HetDiffs += gF->haplotypeString1[j] == hap2Char ? 0 : 1;
                    hap2ToPredictedHap2HetDiffs += gF->haplotypeString2[j] == hap2Char ? 0 : 1;
                }

                // Check haplotype posterior probabilities
                CuAssertTrue(testCase, gF->haplotypeProbs1[j] >= 0.0);
                CuAssertTrue(testCase, gF->haplotypeProbs1[j] <= 1.0);
                CuAssertTrue(testCase, gF->haplotypeProbs2[j] >= 0.0);
                CuAssertTrue(testCase, gF->haplotypeProbs2[j] <= 1.0);
            }

            // Pick the best pairing of the haplotypes to report
            if(hap1ToPredictedHap1Diffs + hap2ToPredictedHap2Diffs >
                    hap2ToPredictedHap1Diffs + hap1ToPredictedHap2Diffs) {
                hap1ToPredictedHap1Diffs = hap2ToPredictedHap1Diffs;
                hap2ToPredictedHap2Diffs = hap1ToPredictedHap2Diffs;

                hap1ToPredictedHap1HetDiffs = hap2ToPredictedHap1HetDiffs;
                hap2ToPredictedHap2HetDiffs = hap1ToPredictedHap2HetDiffs;
            }

            fprintf(stderr, "For an HMM of length: %" PRIi64 " got %" PRIi64 " genotypes correct, rate: %f\n"
                    " Got %" PRIi64 " hets, correct: %" PRIi64 ", rate: %f\n"
                    " Got %" PRIi64 " hap1 differences, rate: %f\n"
                    " Got %" PRIi64 " hap2 differences, rate: %f\n"
                    " Got %" PRIi64 " hap1 het differences, rate: %f\n"
                    " Got %" PRIi64 " hap2 het differences, rate: %f\n",
                    hmm->refLength,
                    correctGenotypes, (float)correctGenotypes/hmm->refLength,
                    totalHets, correctHets, (float)correctHets/totalHets,
                    hap1ToPredictedHap1Diffs, (float)hap1ToPredictedHap1Diffs/hmm->refLength,
                    hap2ToPredictedHap2Diffs, (float)hap2ToPredictedHap2Diffs/hmm->refLength,
                    hap1ToPredictedHap1HetDiffs, (float)hap1ToPredictedHap1HetDiffs/totalHets,
                    hap2ToPredictedHap2HetDiffs, (float)hap2ToPredictedHap2HetDiffs/totalHets);

            fprintf(stderr, "Avg posterior prob. of correct genotype call: %f, avg. posterior"
                    " prob. of incorrect genotype call: %f\n", probsOfCorrectGenotypes/correctGenotypes,
                    probsOfIncorrectGenotypes/(hmm->refLength - correctGenotypes));

            // Cleanup
            stGenomeFragment_destruct(gF);
            stList_destruct(traceBackPath);
            stSet_destruct(profileSeqsPartition1);
            stSet_destruct(profileSeqsPartition2);
        }

        fprintf(stderr, " For %" PRIi64 " hap 1 sequences and %" PRIi64 " hap 2 sequences there were %" PRIi64
                " hmms and %" PRIi64 " partition errors, with %f partition errors per hmm\n",
                stList_length(profileSeqs1), stList_length(profileSeqs2),
                stList_length(hmms), totalPartitionError, (float)totalPartitionError/stList_length(hmms));

        // Cleanup
        stList_destruct(hmms);
        stList_destruct(profileSeqs);
        stList_destruct(referenceSeqs);
        stList_destruct(hapSeqs1);
        stList_destruct(hapSeqs2);
        stSet_destruct(profileSeqs1Set);
        stSet_destruct(profileSeqs2Set);
        stList_destruct(profileSeqs1);
        stList_destruct(profileSeqs2);
        stRPHmmParameters_destruct(params);
    }

    int64_t totalTime = time(NULL) - startTime;

    fprintf(stderr, " For %i tests there were avg. %f hap 1 sequences and avg. %f hap 2 sequences there were "
            " avg. %f partition errors in %" PRIi64 " seconds \n",
                    RANDOM_TEST_NO,
                    (float)totalProfile1SeqsOverAllTests/RANDOM_TEST_NO,
                    (float)totalProfile2SeqsOverAllTests/RANDOM_TEST_NO,
                    (float)totalPartitionErrorsOverAllTests/RANDOM_TEST_NO,
                    totalTime);
}

void test_systemSingleReferenceFullLengthReads(CuTest *testCase) {
    int64_t minReferenceSeqNumber = 1;
    int64_t maxReferenceSeqNumber = 1;
    int64_t minReferenceLength = 1000;
    int64_t maxReferenceLength = 1000;
    int64_t minCoverage = 30;
    int64_t maxCoverage = 30;
    int64_t minReadLength = 1000;
    int64_t maxReadLength = 1000;
    int64_t maxPartitionsInAColumn = 50;
    double hetRate = 0.001;
    double readErrorRate = 0.1;
    bool maxNotSumTransitions = 0;
    bool splitHmmsWherePhasingUncertain = 0;
    int64_t minReadCoverageToSupportPhasingBetweenHeterozygousSites = 0;
    bool printHmm = 1;

    test_systemTest(testCase, minReferenceSeqNumber, maxReferenceSeqNumber,
            minReferenceLength, maxReferenceLength, minCoverage, maxCoverage,
            minReadLength, maxReadLength, maxPartitionsInAColumn, hetRate, readErrorRate,
            maxNotSumTransitions, splitHmmsWherePhasingUncertain,
            minReadCoverageToSupportPhasingBetweenHeterozygousSites, printHmm);
}

void test_systemSingleReferenceFixedLengthReads(CuTest *testCase) {
    int64_t minReferenceSeqNumber = 1;
    int64_t maxReferenceSeqNumber = 1;
    int64_t minReferenceLength = 40000;
    int64_t maxReferenceLength = 40000;
    int64_t minCoverage = 30;
    int64_t maxCoverage = 30;
    int64_t minReadLength = 3000;
    int64_t maxReadLength = 3000;
    int64_t maxPartitionsInAColumn = 50;
    double hetRate = 0.0007;
    double readErrorRate = 0.05;
    bool maxNotSumTransitions = 0;
    bool splitHmmsWherePhasingUncertain = 1;
    int64_t minReadCoverageToSupportPhasingBetweenHeterozygousSites = 3;
    bool printHmm = 0;

    test_systemTest(testCase, minReferenceSeqNumber, maxReferenceSeqNumber,
            minReferenceLength, maxReferenceLength, minCoverage, maxCoverage,
            minReadLength, maxReadLength, maxPartitionsInAColumn, hetRate, readErrorRate,
            maxNotSumTransitions, splitHmmsWherePhasingUncertain,
            minReadCoverageToSupportPhasingBetweenHeterozygousSites, printHmm);
}

void test_systemSingleReference(CuTest *testCase) {
    int64_t minReferenceSeqNumber = 1;
    int64_t maxReferenceSeqNumber = 1;
    int64_t minReferenceLength = 1000;
    int64_t maxReferenceLength = 1000;
    int64_t minCoverage = 30;
    int64_t maxCoverage = 30;
    int64_t minReadLength = 10;
    int64_t maxReadLength = 300;
    int64_t maxPartitionsInAColumn = 50;
    double hetRate = 0.01;
    double readErrorRate = 0.01;
    bool maxNotSumTransitions = 0;
    bool splitHmmsWherePhasingUncertain = 1;
    int64_t minReadCoverageToSupportPhasingBetweenHeterozygousSites = 15;
    bool printHmm = 1;

    test_systemTest(testCase, minReferenceSeqNumber, maxReferenceSeqNumber,
            minReferenceLength, maxReferenceLength, minCoverage, maxCoverage,
            minReadLength, maxReadLength, maxPartitionsInAColumn,
            hetRate, readErrorRate,
            maxNotSumTransitions, splitHmmsWherePhasingUncertain,
            minReadCoverageToSupportPhasingBetweenHeterozygousSites, printHmm);
}

void test_systemMultipleReferences(CuTest *testCase) {
    int64_t minReferenceSeqNumber = 2;
    int64_t maxReferenceSeqNumber = 5;
    int64_t minReferenceLength = 500;
    int64_t maxReferenceLength = 1000;
    int64_t minCoverage = 5;
    int64_t maxCoverage = 30;
    int64_t minReadLength = 10;
    int64_t maxReadLength = 300;
    int64_t maxPartitionsInAColumn = 50;
    double hetRate = 0.01;
    double readErrorRate = 0.01;
    bool maxNotSumTransitions = 0;
    bool splitHmmsWherePhasingUncertain = 1;
    int64_t minReadCoverageToSupportPhasingBetweenHeterozygousSites = 15;
    bool printHmm = 0;

    test_systemTest(testCase, minReferenceSeqNumber, maxReferenceSeqNumber,
            minReferenceLength, maxReferenceLength, minCoverage, maxCoverage,
            minReadLength, maxReadLength, maxPartitionsInAColumn, hetRate, readErrorRate,
            maxNotSumTransitions, splitHmmsWherePhasingUncertain,
            minReadCoverageToSupportPhasingBetweenHeterozygousSites, printHmm);
}

void test_popCount64(CuTest *testCase) {
    CuAssertIntEquals(testCase, popcount64(0), 0);
    CuAssertIntEquals(testCase, popcount64(1), 1);
    CuAssertIntEquals(testCase, popcount64(2), 1);
    CuAssertIntEquals(testCase, popcount64(3), 2);
    CuAssertIntEquals(testCase, popcount64(0xFF), 8);
    CuAssertIntEquals(testCase, popcount64(0xFFFFFFFF), 32);
    CuAssertIntEquals(testCase, popcount64(0xFFFFFFFFFFFFFFFF), 64);
    CuAssertIntEquals(testCase, popcount64(0x1111111111111111), 16);
}

static double getExpectedInstanceNumberSimple(uint8_t **seqs, uint64_t partition,
        int64_t depth, int64_t length,
        int64_t position, int64_t characterIndex) {
    uint64_t expectation = 0;
    for(int64_t i=0; i<depth; i++) {
        if((partition & ((uint64_t)1 << i)) != 0) {
            expectation += seqs[i][position*ALPHABET_SIZE + characterIndex];
        }
    }
    return (double)expectation/255;
}

static uint64_t getRandomPartition(int64_t depth) {
    //return 0xFFFFFFFFFFFFFFFF;
    return st_randomInt(0, 0xFFFFFFFFFFFFFFFF/2);
}

void test_bitCountVectors(CuTest *testCase) {
    for(int64_t depth=0; depth<64; depth++) {
        for(int64_t test=0; test<100; test++) {
            // Make column as set of uint8_t sequences
            int64_t length = st_randomInt(0, 10);
            uint8_t **seqs = st_malloc(sizeof(uint8_t *) * depth);
            for(int64_t i=0; i<depth; i++) {
                seqs[i] = st_calloc(length * ALPHABET_SIZE, sizeof(uint8_t));
                // Initialise probs randomly
                for(int64_t j=0; j<ALPHABET_SIZE*length; j++) {
                    seqs[i][j] = st_randomInt(0, 255);
                }
            }

//            fprintf(stderr, "Test: %" PRIi64 " depth: %" PRIi64 " alphabet_size: %i length: %" PRIi64 "\n",
//                    test, depth, ALPHABET_SIZE, length);

            // Calculate the bit vectors
            uint64_t *countBitVectors = calculateCountBitVectors(seqs, depth, length);

            // Partition
            uint64_t partition = getRandomPartition(depth);

            // Test we get the expected output
            for(int64_t i=0; i<length; i++) {
                for(int64_t j=0; j<ALPHABET_SIZE; j++) {
                    CuAssertDblEquals(testCase,
                            getExpectedInstanceNumberSimple(seqs, partition, depth,
                                    length, i, j),
                            (double)getExpectedInstanceNumber(countBitVectors, depth,
                                    partition, i, j)/ALPHABET_MAX_PROB, 0.0000001);

                }
            }

            // Cleanup
            for(int64_t i=0; i<depth; i++) {
                free(seqs[i]);
            }
            free(seqs);
            free(countBitVectors);
        }
    }
}

void buildComponent(stRPHmm *hmm1, stSortedSet *component, stSet *seen) {
    stSet_insert(seen, hmm1);
    stSortedSetIterator *it = stSortedSet_getIterator(component);
    stRPHmm *hmm2;
    while((hmm2 = stSortedSet_getNext(it)) != NULL) {
        if(stRPHmm_overlapOnReference(hmm1, hmm2) && stSet_search(seen, hmm2) == NULL) {
            buildComponent(hmm2, component, seen);
        }
    }
    stSortedSet_destructIterator(it);
}

void test_getOverlappingComponents(CuTest *testCase) {
    int64_t minReferenceSeqNumber = 1;
    int64_t maxReferenceSeqNumber = 10;
    int64_t minReferenceLength = 1000;
    int64_t maxReferenceLength = 1000;
    int64_t minCoverage = 1;
    int64_t maxCoverage = 30;
    int64_t minReadLength = 100;
    int64_t maxReadLength = 100;
    int64_t maxPartitionsInAColumn = 100;
    double hetRate = 0.02;
    double readErrorRate = 0.01;
    bool maxNotSumTransitions = 0;

    for(int64_t test=0; test<RANDOM_TEST_NO; test++) {
        fprintf(stderr, "Starting test iteration: #%" PRIi64 "\n", test);

        stRPHmmParameters *params = getHmmParams(maxPartitionsInAColumn,
                        hetRate, readErrorRate,
                        maxNotSumTransitions, 0);

        stList *referenceSeqs = stList_construct3(0, free);
        stList *hapSeqs1 = stList_construct3(0, free);
        stList *hapSeqs2 = stList_construct3(0, free);
        stList *profileSeqs1 = stList_construct3(0, (void (*)(void *))stProfileSeq_destruct);
        stList *profileSeqs2 = stList_construct3(0, (void (*)(void *))stProfileSeq_destruct);

        // Creates reference sequences
        // Generates two haplotypes for each reference sequence
        // Generates profile sequences from each haplotype
        simulateReads(referenceSeqs, hapSeqs1, hapSeqs2,
                        profileSeqs1, profileSeqs2,
                        minReferenceSeqNumber, maxReferenceSeqNumber,
                        minReferenceLength, maxReferenceLength,
                        minCoverage, maxCoverage,
                        minReadLength, maxReadLength,
                        hetRate, readErrorRate);

        // Make simple hmms
        stSortedSet *readHmms = stSortedSet_construct3(stRPHmm_cmpFn, NULL);
        for(int64_t i=0; i<stList_length(profileSeqs1); i++) {
            stRPHmm *hmm = stRPHmm_construct(stList_get(profileSeqs1, i), params);
            CuAssertTrue(testCase, stSortedSet_search(readHmms, hmm) == NULL);
            stSortedSet_insert(readHmms, hmm);
            CuAssertTrue(testCase, stSortedSet_search(readHmms, hmm) == hmm);
        }
        for(int64_t i=0; i<stList_length(profileSeqs2); i++) {
            stRPHmm *hmm = stRPHmm_construct(stList_get(profileSeqs2, i), params);
            CuAssertTrue(testCase, stSortedSet_search(readHmms, hmm) == NULL);
            stSortedSet_insert(readHmms, hmm);
            CuAssertTrue(testCase, stSortedSet_search(readHmms, hmm) == hmm);
        }
        CuAssertIntEquals(testCase, stSortedSet_size(readHmms), stList_length(profileSeqs1) + stList_length(profileSeqs2));

        // Organise HMMs into "tiling paths" consisting of sequences of hmms that do not overlap
        stList *readHmmsList = stSortedSet_getList(readHmms);
        stList *tilingPaths = getTilingPaths(readHmms);

        // Check tiling paths - that every hmm is in a tiling path and that hmms do not overlap
        // and that tiling paths are minimal
        stSet *seen = stSet_construct();
        for(int64_t i=0; i<stList_length(tilingPaths); i++) {
            stList *tilingPath = stList_get(tilingPaths, i);
            for(int64_t j=0; j<stList_length(tilingPath); j++) {
                stRPHmm *hmm1 = stList_get(tilingPath, j);

                // The hmm1 is in only one tiling path
                CuAssertTrue(testCase, stSet_search(seen, hmm1) == NULL);
                stSet_insert(seen, hmm1);

                // The hmm1 does not overlap with any other hmm1 in the tiling path
                // and precedes all subsequent hmms in tiling path
                for(int64_t k=j+1; k<stList_length(tilingPath); k++) {
                    stRPHmm *hmm2 = stList_get(tilingPath, k);
                    CuAssertTrue(testCase, !stRPHmm_overlapOnReference(hmm1, hmm2));
                    CuAssertTrue(testCase, stRPHmm_cmpFn(hmm1, hmm2) < 0);
                }

                // The next hmm1 in the tiling path is as close as possible while not overlapping
                for(int64_t k=i+1; k<stList_length(tilingPaths); k++) {
                    stList *tilingPath2 = stList_get(tilingPaths, k);
                    for(int64_t l=0; l<stList_length(tilingPath2); l++) {
                        stRPHmm *hmm2 = stList_get(tilingPath2, l);
                        // If not overlapping and hmm2 occurs after hmm1
                        if(!stRPHmm_overlapOnReference(hmm1, hmm2) && stRPHmm_cmpFn(hmm1, hmm2) < 0) {
                            // There must be an hmm following hmm1 in the tiling path
                            CuAssertTrue(testCase, j+1 < stList_length(tilingPath));

                            stRPHmm *hmm3 = stList_get(tilingPath, j+1);
                            // hmm3 must precede this next hmm in the sort
                            CuAssertTrue(testCase, stRPHmm_cmpFn(hmm3, hmm2) <= 0);
                        }
                    }
                }
            }
        }
        // Check we accounted for all the hmms
        CuAssertIntEquals(testCase, stSet_size(seen), stList_length(readHmmsList));
        stSet_destruct(seen);

        // Get components
        while(stList_length(tilingPaths) > 1) {
            stList *tilingPath2 = stList_pop(tilingPaths);
            stList *tilingPath1 = stList_pop(tilingPaths);

            stSet *components = getOverlappingComponents(tilingPath1, tilingPath2);
            // Check that all the hmms in the tiling paths are in one component
            // Check that within a component hmms overlap
            stHash *hmmToComponent = stHash_construct();
            stSetIterator *componentsIt = stSet_getIterator(components);
            stSortedSet *component;
            while((component = stSet_getNext(componentsIt)) != NULL) {
                //
                stSortedSetIterator *componentIt = stSortedSet_getIterator(component);
                stRPHmm *hmm;
                while((hmm = stSortedSet_getNext(componentIt)) != NULL) {

                    CuAssertTrue(testCase, stHash_search(hmmToComponent, hmm) == NULL);
                    stHash_insert(hmmToComponent, hmm, component);
                }
                stSortedSet_destructIterator(componentIt);

                // Check that hmms in component are single overlap component
                seen = stSet_construct();
                buildComponent(stSortedSet_getFirst(component), component, seen);
                CuAssertIntEquals(testCase, stSet_size(seen), stSortedSet_size(component));
                stSet_destruct(seen);
            }
            stSet_destructIterator(componentsIt);

            CuAssertIntEquals(testCase, stHash_size(hmmToComponent), stList_length(tilingPath1)+stList_length(tilingPath2));

            // Check that no hmms overlap between components
            stList *hmmsInComponents = stHash_getKeys(hmmToComponent);
            for(int64_t i=0; i<stList_length(hmmsInComponents); i++) {
                stRPHmm *hmm1 = stList_get(hmmsInComponents, i);
                for(int64_t j=i+1; j<stList_length(hmmsInComponents); j++) {
                    stRPHmm *hmm2 = stList_get(hmmsInComponents, j);

                    if(stRPHmm_overlapOnReference(hmm1, hmm2)) {
                        CuAssertTrue(testCase, stHash_search(hmmToComponent, hmm1) == stHash_search(hmmToComponent, hmm2));
                    }
                }
            }

            // Cleanup
            stList_destruct(tilingPath1);
            stList_destruct(tilingPath2);
            stHash_destruct(hmmToComponent);
            stList_destruct(hmmsInComponents);
            stSet_destruct(components);
        }

        // Cleanup
        for(int64_t i=0; i<stList_length(readHmmsList); i++) {
            stRPHmm_destruct(stList_get(readHmmsList, i), 1);
        }
        stList_destruct(readHmmsList);
        while(stList_length(tilingPaths) > 0) {
            stList_destruct(stList_pop(tilingPaths));
        }
        stList_destruct(tilingPaths);
        stList_destruct(referenceSeqs);
        stList_destruct(hapSeqs1);
        stList_destruct(hapSeqs2);
        stList_destruct(profileSeqs1);
        stList_destruct(profileSeqs2);
        stRPHmmParameters_destruct(params);
    }
}

/*
 * Functions to test emission function
 */

void test_emissionLogProbability(CuTest *testCase) {
    int64_t minReferenceSeqNumber = 1;
    int64_t maxReferenceSeqNumber = 10;
    int64_t minReferenceLength = 1000;
    int64_t maxReferenceLength = 1000;
    int64_t minCoverage = 10;
    int64_t maxCoverage = 10;
    int64_t minReadLength = 10;
    int64_t maxReadLength = 1000;
    int64_t maxPartitionsInAColumn = 100;
    double hetRate = 0.02;
    double readErrorRate = 0.01;
    bool maxNotSumTransitions = 0;

    for(int64_t test=0; test<RANDOM_TEST_NO; test++) {
        fprintf(stderr, "Starting test iteration: #%" PRIi64 "\n", test);

        stRPHmmParameters *params = getHmmParams(maxPartitionsInAColumn,
                        hetRate, readErrorRate,
                        maxNotSumTransitions, 0);

        stList *referenceSeqs = stList_construct3(0, free);
        stList *hapSeqs1 = stList_construct3(0, free);
        stList *hapSeqs2 = stList_construct3(0, free);
        stList *profileSeqs1 = stList_construct3(0, (void (*)(void *))stProfileSeq_destruct);
        stList *profileSeqs2 = stList_construct3(0, (void (*)(void *))stProfileSeq_destruct);

        // Creates reference sequences
        // Generates two haplotypes for each reference sequence
        // Generates profile sequences from each haplotype
        simulateReads(referenceSeqs, hapSeqs1, hapSeqs2,
                        profileSeqs1, profileSeqs2,
                        minReferenceSeqNumber, maxReferenceSeqNumber,
                        minReferenceLength, maxReferenceLength,
                        minCoverage, maxCoverage,
                        minReadLength, maxReadLength,
                        hetRate, readErrorRate);


        // Creates read HMMs
        stList *profileSeqs = stList_copy(profileSeqs1, NULL);
        stList_appendAll(profileSeqs, profileSeqs2);
        stList *hmms = getRPHmms(profileSeqs, params);

        // For each hmm
        while(stList_length(hmms) > 0) {
            stRPHmm *hmm = stList_pop(hmms);

            // For each column
            stRPColumn *column = hmm->firstColumn;
            while(1) {
                // Get bit count vectors
                uint64_t *bitCountVectors = calculateCountBitVectors(
                        column->seqs, column->depth, column->length);

                // For each cell
                stRPCell *cell = column->head;
                do {
                    // Check slow and fast way to calculate emission probabilities
                    // are equivalent
                    double e1 = emissionLogProbabilitySlow(column, cell,
                            bitCountVectors, params, 1);
                    double e2 = emissionLogProbability(column,
                                                        cell, bitCountVectors, params);
                    //st_uglyf("Boo %f %f\n", e1, e2);
                    CuAssertDblEquals(testCase, e1, e2, 0.1);
                } while((cell = cell->nCell) != NULL);

                // Clean up
                free(bitCountVectors);

                if(column->nColumn == NULL) {
                    break;
                }
                column = column->nColumn->nColumn;
            }

            // Clean up hmm
            stRPHmm_destruct(hmm, 1);
        }

        // Clean up
        stList_destruct(profileSeqs);
        stList_destruct(hmms);
        stList_destruct(referenceSeqs);
        stList_destruct(hapSeqs1);
        stList_destruct(hapSeqs2);
        stList_destruct(profileSeqs1);
        stList_destruct(profileSeqs2);
        stRPHmmParameters_destruct(params);
    }
}

CuSuite *stRPHmmTestSuite(void) {
    CuSuite* suite = CuSuiteNew();

    // System level tests
    SUITE_ADD_TEST(suite, test_systemSingleReferenceFullLengthReads);
    SUITE_ADD_TEST(suite, test_systemSingleReferenceFixedLengthReads);
    SUITE_ADD_TEST(suite, test_systemSingleReference);
    SUITE_ADD_TEST(suite, test_systemMultipleReferences);

    // Constituent function tests
    SUITE_ADD_TEST(suite, test_popCount64);
    SUITE_ADD_TEST(suite, test_bitCountVectors);
    SUITE_ADD_TEST(suite, test_getOverlappingComponents);
    SUITE_ADD_TEST(suite, test_emissionLogProbability);

    return suite;
}
