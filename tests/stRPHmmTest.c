/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "stRPHmm.h"
#include <math.h>
#include <time.h>

#define RANDOM_TEST_NO 1

char getRandomBase(int64_t alphabetSize) {
    /*
     * Returns an ascii character starting from ascii symbol 48: '0', '1', '2', ... alphabetsize
     */
    return st_randomInt(FIRST_ALPHABET_CHAR, FIRST_ALPHABET_CHAR+alphabetSize);
}

char *getRandomSequence(int64_t referenceLength, int64_t alphabetSize) {
    /*
     * Creates a random sequence of form [ACTG]*referenceLength
     */
    char *randomSequence = st_malloc(sizeof(char) * (referenceLength+1));
    for(int64_t i=0; i<referenceLength; i++) {
        randomSequence[i] = getRandomBase(alphabetSize);
    }
    randomSequence[referenceLength] = '\0';

    return randomSequence;
}

char *permuteSequence(char *referenceSeq, double hetRate, int64_t alphabetSize) {
    /*
     * Takes a random sequence and returns a copy of it, permuting randomly each position with rate hetRate.
     */
    referenceSeq = stString_copy(referenceSeq);
    int64_t strLength = strlen(referenceSeq);
    for(int64_t i=0; i<strLength; i++) {
        if(st_random() < hetRate) {
            referenceSeq[i] = getRandomBase(alphabetSize);
        }
    }
    return referenceSeq;
}

stProfileSeq *getRandomProfileSeq(char *referenceName, char *hapSeq, int64_t hapLength,
        int64_t readLength, double readErrorRate, int64_t alphabetSize) {
    /*
     * Creates a random "read" of the given length from a random sub-interval
     * of the input sequence with error rate errRate.
     */
    assert(hapLength-readLength >= 0);
    int64_t start = st_randomInt(0, hapLength-readLength+1);

    stProfileSeq *pSeq = stProfileSeq_constructEmptyProfile(referenceName, start,
            readLength, alphabetSize);

    for(int64_t i=0; i<readLength; i++) {
        // Haplotype base or error at random
        char b = st_random() < readErrorRate ? getRandomBase(alphabetSize) : hapSeq[start+i];
        assert(b - FIRST_ALPHABET_CHAR >= 0);
        assert(b - FIRST_ALPHABET_CHAR < alphabetSize);
        // Fill in the profile probabilities according to the chosen base
        pSeq->profileProbs[i*alphabetSize + b - FIRST_ALPHABET_CHAR] = ALPHABET_MAX_PROB;
    }

    return pSeq;
}

static stRPHmmParameters *getHmmParams(int64_t maxPartitionsInAColumn,
        double hetRate, double readErrorRate, int64_t alphabetSize,
        bool maxNotSumEmissions, bool maxNotSumTransitions) {
    // Substitution models
    stSubModel *hetSubModel = stSubModel_constructEmptyModel(alphabetSize);
    stSubModel *readErrorSubModel = stSubModel_constructEmptyModel(alphabetSize);

    // Fill in substitition matrix with simple symmetric probs matching
    // the read error rate
    assert(readErrorRate <= 1.0);
    assert(readErrorRate >= 0.0);
    assert(alphabetSize >= 2);
    for(int64_t i=0; i<alphabetSize; i++) {
        for(int64_t j=0; j<alphabetSize; j++) {
            stSubModel_setSubstitutionProb(readErrorSubModel, i, j, i==j ?
                    1.0-readErrorRate : readErrorRate/(alphabetSize-1));
            stSubModel_setSubstitutionProb(hetSubModel, i, j, i==j ?
                                    1.0-hetRate : hetRate/(alphabetSize-1));
        }
    }

    return stRPHmmParameters_construct(hetSubModel, readErrorSubModel, maxNotSumEmissions,
                    maxNotSumTransitions, maxPartitionsInAColumn, MAX_READ_PARTITIONING_DEPTH);
}

static void simulateReads(stList *referenceSeqs, stList *hapSeqs1, stList *hapSeqs2,
        stList *profileSeqs1, stList *profileSeqs2,
        int64_t minReferenceSeqNumber, int64_t maxReferenceSeqNumber,
        int64_t minReferenceLength, int64_t maxReferenceLength,
        int64_t minCoverage, int64_t maxCoverage,
        int64_t minReadLength, int64_t maxReadLength,
        double hetRate, double readErrorRate, int64_t alphabetSize) {
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
        char *referenceSeq = getRandomSequence(referenceLength, alphabetSize);
        // Reference name
        char *referenceName = stString_print("Reference_%" PRIi64 "", i);

        stList_append(referenceSeqs, referenceSeq);

        // Create haplotype sequences for reference
        char *haplotypeSeq1 = permuteSequence(referenceSeq, hetRate, alphabetSize);
        char *haplotypeSeq2 = permuteSequence(referenceSeq, hetRate, alphabetSize);

        stList_append(hapSeqs1, haplotypeSeq1);
        stList_append(hapSeqs2, haplotypeSeq2);

        // Print info about simulated sequences
        fprintf(stderr, "Ref seq: %s\n", referenceSeq);
        fprintf(stderr, "  Hap 1: %s\n", haplotypeSeq1);
        fprintf(stderr, "  Hap 2: %s\n", haplotypeSeq2);

        // Create read sequences to given coverage
        int64_t coverage = st_randomInt(minCoverage, maxCoverage+1);
        int64_t totalBasesToSimulate = coverage * referenceLength;
        fprintf(stderr, "Total bases to simulate: %" PRIi64 " for coverage: %" PRIi64 "\n", totalBasesToSimulate, coverage);
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
                    referenceLength, readLength, readErrorRate, alphabetSize);
            stList_append(readSeqs, pSeq);
            totalBasesToSimulate -= readLength;
            fprintf(stderr, "Simulating read from haplotype: %s\n", hapSeq);
            //stProfileSeq_print(pSeq, stderr, 1);
        }

        // Cleanup
        free(referenceName);
    }
}

static void test_systemTest(CuTest *testCase, int64_t minReferenceSeqNumber, int64_t maxReferenceSeqNumber,
        int64_t minReferenceLength, int64_t maxReferenceLength, int64_t minCoverage, int64_t maxCoverage,
        int64_t minReadLength, int64_t maxReadLength,
        int64_t maxPartitionsInAColumn, double hetRate, double readErrorRate, int64_t alphabetSize,
        bool maxNotSumEmissions, bool maxNotSumTransitions) {
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
            "\talphabetSize: %" PRIi64 "\n"
            "\tmaxNotSumEmissions: %i\n"
            "\tmaxNotSumTransitions: %i\n",
            minReferenceSeqNumber, maxReferenceSeqNumber,
            minReferenceLength, maxReferenceLength, minCoverage, maxCoverage,
            minReadLength, maxReadLength,
            maxPartitionsInAColumn, (float)hetRate, (float)readErrorRate, alphabetSize,
            maxNotSumEmissions, maxNotSumTransitions);

    int64_t totalProfile1SeqsOverAllTests = 0;
    int64_t totalProfile2SeqsOverAllTests = 0;
    int64_t totalPartitionErrorsOverAllTests = 0;

    time_t startTime = time(NULL);

    for(int64_t test=0; test<RANDOM_TEST_NO; test++) {

        fprintf(stderr, "Starting test iteration: #%" PRIi64 "\n", test);

        stRPHmmParameters *params = getHmmParams(maxPartitionsInAColumn,
                hetRate, readErrorRate, alphabetSize,
                maxNotSumEmissions, maxNotSumTransitions);

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
                        hetRate, readErrorRate, alphabetSize);

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
                CuAssertTrue(testCase, hmm->refStart <= profileSeq->refStart);
                CuAssertTrue(testCase, hmm->refStart + hmm->refLength >= profileSeq->refStart + profileSeq->length);
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
                    CuAssertPtrEquals(testCase, &(profileSeq->profileProbs[(column->refStart-profileSeq->refStart) * alphabetSize]), column->seqs[j]);
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
            stRPHmm_print(hmm, stderr, 1, 1);

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
                    CuAssertDblEquals(testCase, 1.0, totalProb, 0.03);
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
                    stList_length(profileSeqs1),
                    stList_length(profileSeqs2), partitionErrors);

            for(int64_t k=0; k<stList_length(hapSeqs1); k++) {
                fprintf(stderr, "Hap1: %s\n", stList_get(hapSeqs1, k));
            }
            for(int64_t k=0; k<stList_length(hapSeqs2); k++) {
                fprintf(stderr, "Hap2: %s\n", stList_get(hapSeqs2, k));
            }
            printPartition(stderr, profileSeqsPartition1, profileSeqsPartition2);

            totalPartitionError += partitionErrors;

            totalProfile1SeqsOverAllTests += stList_length(profileSeqs1);
            totalProfile2SeqsOverAllTests += stList_length(profileSeqs2);
            totalPartitionErrorsOverAllTests += partitionErrors;

            // Cleanup
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

static void test_systemSingleReferenceFullLengthReads(CuTest *testCase) {
    int64_t minReferenceSeqNumber = 1;
    int64_t maxReferenceSeqNumber = 1;
    int64_t minReferenceLength = 1000;
    int64_t maxReferenceLength = 1000;
    int64_t minCoverage = 10;
    int64_t maxCoverage = 10;
    int64_t minReadLength = 1000;
    int64_t maxReadLength = 1000;
    int64_t maxPartitionsInAColumn = 100;
    double hetRate = 0.005;
    double readErrorRate = 0.1;
    int64_t alphabetSize = 2;
    bool maxNotSumEmissions = 1;
    bool maxNotSumTransitions = 0;

    test_systemTest(testCase, minReferenceSeqNumber, maxReferenceSeqNumber,
            minReferenceLength, maxReferenceLength, minCoverage, maxCoverage,
            minReadLength, maxReadLength, maxPartitionsInAColumn, hetRate, readErrorRate, alphabetSize,
            maxNotSumEmissions, maxNotSumTransitions);
}

static void test_systemSingleReferenceFixedLengthReads(CuTest *testCase) {
    return;
    int64_t minReferenceSeqNumber = 1;
    int64_t maxReferenceSeqNumber = 1;
    int64_t minReferenceLength = 1000;
    int64_t maxReferenceLength = 1000;
    int64_t minCoverage = 10;
    int64_t maxCoverage = 10;
    int64_t minReadLength = 200;
    int64_t maxReadLength = 200;
    int64_t maxPartitionsInAColumn = 100;
    double hetRate = 0.01;
    double readErrorRate = 0.1;
    int64_t alphabetSize = 4;
    bool maxNotSumEmissions = 1;
    bool maxNotSumTransitions = 0;

    test_systemTest(testCase, minReferenceSeqNumber, maxReferenceSeqNumber,
            minReferenceLength, maxReferenceLength, minCoverage, maxCoverage,
            minReadLength, maxReadLength, maxPartitionsInAColumn, hetRate, readErrorRate, alphabetSize,
            maxNotSumEmissions, maxNotSumTransitions);
}

static void test_systemSingleReference(CuTest *testCase) {
    return;
    int64_t minReferenceSeqNumber = 1;
    int64_t maxReferenceSeqNumber = 1;
    int64_t minReferenceLength = 1000;
    int64_t maxReferenceLength = 1000;
    int64_t minCoverage = 30;
    int64_t maxCoverage = 30;
    int64_t minReadLength = 10;
    int64_t maxReadLength = 300;
    int64_t maxPartitionsInAColumn = 100;
    double hetRate = 0.02;
    double readErrorRate = 0.01;
    int64_t alphabetSize = 3;
    bool maxNotSumEmissions = 1;
    bool maxNotSumTransitions = 0;

    test_systemTest(testCase, minReferenceSeqNumber, maxReferenceSeqNumber,
            minReferenceLength, maxReferenceLength, minCoverage, maxCoverage,
            minReadLength, maxReadLength, maxPartitionsInAColumn,
            hetRate, readErrorRate, alphabetSize,
            maxNotSumEmissions, maxNotSumTransitions);
}

static void test_systemMultipleReferences(CuTest *testCase) {
    return;
    int64_t minReferenceSeqNumber = 2;
    int64_t maxReferenceSeqNumber = 5;
    int64_t minReferenceLength = 500;
    int64_t maxReferenceLength = 1000;
    int64_t minCoverage = 5;
    int64_t maxCoverage = 20;
    int64_t minReadLength = 10;
    int64_t maxReadLength = 300;
    int64_t maxPartitionsInAColumn = 100;
    double hetRate = 0.02;
    double readErrorRate = 0.01;
    int64_t alphabetSize = 5;
    bool maxNotSumEmissions = 1;
    bool maxNotSumTransitions = 0;

    test_systemTest(testCase, minReferenceSeqNumber, maxReferenceSeqNumber,
            minReferenceLength, maxReferenceLength, minCoverage, maxCoverage,
            minReadLength, maxReadLength, maxPartitionsInAColumn, hetRate, readErrorRate, alphabetSize,
            maxNotSumEmissions, maxNotSumTransitions);
}

static void test_popCount64(CuTest *testCase) {
    CuAssertIntEquals(testCase, popcount64(0), 0);
    CuAssertIntEquals(testCase, popcount64(1), 1);
    CuAssertIntEquals(testCase, popcount64(2), 1);
    CuAssertIntEquals(testCase, popcount64(3), 2);
    CuAssertIntEquals(testCase, popcount64(0xFF), 8);
    CuAssertIntEquals(testCase, popcount64(0xFFFFFFFF), 32);
    CuAssertIntEquals(testCase, popcount64(0xFFFFFFFFFFFFFFFF), 64);
    CuAssertIntEquals(testCase, popcount64(0x1111111111111111), 16);
}

static double getExpectedInstanceNumberSimple(uint8_t **seqs, uint64_t partition, int64_t depth, int64_t length,
        int64_t alphabetSize, int64_t position, int64_t characterIndex) {
    uint64_t expectation = 0;
    for(int64_t i=0; i<depth; i++) {
        if((partition & ((uint64_t)1 << i)) != 0) {
            expectation += seqs[i][position*alphabetSize + characterIndex];
        }
    }
    return (double)expectation/255;
}

static uint64_t getRandomPartition(int64_t depth) {
    //return 0xFFFFFFFFFFFFFFFF;
    return st_randomInt(0, 0xFFFFFFFFFFFFFFFF/2);
}

static void test_bitCountVectors(CuTest *testCase) {
    return;
    for(int64_t depth=0; depth<64; depth++) {
        for(int64_t test=0; test<100; test++) {
            // Make column as set of uint8_t sequences
            int64_t alphabetSize = st_randomInt(1, 10);
            int64_t length = st_randomInt(0, 10);
            uint8_t **seqs = st_malloc(sizeof(uint8_t *) * depth);
            for(int64_t i=0; i<depth; i++) {
                seqs[i] = st_calloc(length * alphabetSize, sizeof(uint8_t));
                // Initialise probs randomly
                for(int64_t j=0; j<alphabetSize*length; j++) {
                    seqs[i][j] = st_randomInt(0, 255);
                }
            }

            fprintf(stderr, "Test: %" PRIi64 " depth: %" PRIi64 " alphabetSize: %" PRIi64 " length: %" PRIi64 "\n",
                    test, depth, alphabetSize, length);

            // Calculate the bit vectors
            uint64_t *countBitVectors = calculateCountBitVectors(seqs, depth, length, alphabetSize);

            // Partition
            uint64_t partition = getRandomPartition(depth);

            // Test we get the expected output
            for(int64_t i=0; i<length; i++) {
                for(int64_t j=0; j<alphabetSize; j++) {
                    CuAssertDblEquals(testCase,
                            getExpectedInstanceNumberSimple(seqs, partition, depth,
                                    length, alphabetSize, i, j),
                            (double)getExpectedInstanceNumber(countBitVectors, depth,
                                    partition, i, j, alphabetSize)/ALPHABET_MAX_PROB, 0.0000001);

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

static void test_getOverlappingComponents(CuTest *testCase) {
    return;
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
    int64_t alphabetSize = 2;
    bool maxNotSumEmissions = 1;
    bool maxNotSumTransitions = 0;

    for(int64_t test=0; test<RANDOM_TEST_NO; test++) {
        fprintf(stderr, "Starting test iteration: #%" PRIi64 "\n", test);

        stRPHmmParameters *params = getHmmParams(maxPartitionsInAColumn,
                        hetRate, readErrorRate, alphabetSize,
                        maxNotSumEmissions, maxNotSumTransitions);

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
                        hetRate, readErrorRate, alphabetSize);

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

double getLogProbOfReadCharactersSlow(stSubModel *alphabet, uint64_t *expectedInstanceNumbers,
        int64_t sourceCharacterIndex) {
    /*
     * Get the log probability of a given source character given the expected number of instances of
     * each character in the reads.
     */
    double logCharacterProb = stSubModel_getSubstitutionProbSlow(alphabet, sourceCharacterIndex, 0) *
            ((double)expectedInstanceNumbers[0])/ALPHABET_MAX_PROB;

    for(int64_t i=1; i<alphabet->alphabetSize; i++) {
        logCharacterProb += stSubModel_getSubstitutionProbSlow(alphabet, sourceCharacterIndex, i) *
                ((double)expectedInstanceNumbers[i])/ALPHABET_MAX_PROB;
    }

    return logCharacterProb;
}

void columnIndexLogHapProbabilitySlow(stRPColumn *column, uint64_t index,
        uint64_t partition, uint64_t *bitCountVectors, stRPHmmParameters *params, double *rootCharacterProbs) {
    /*
     * Get the probabilities of the "root" characters for a given read sub-partition and a haplotype.
     */
    // For each possible read character calculate the expected number of instances in the
    // partition and store counts in an array
    uint64_t expectedInstanceNumbers[params->alphabetSize];
    for(int64_t i=0; i<params->alphabetSize; i++) {
        expectedInstanceNumbers[i] = getExpectedInstanceNumber(bitCountVectors,
                               column->depth, partition, index, i, params->alphabetSize);
    }

    // Calculate the probability of the read characters for each possible haplotype character
    double characterProbsHap[params->alphabetSize];
    for(int64_t i=0; i<params->alphabetSize; i++) {
        characterProbsHap[i] = getLogProbOfReadCharactersSlow(params->readErrorSubModel, expectedInstanceNumbers, i);
    }

    // Calculate the probability of haplotype characters and read characters for each root character
    for(int64_t i=0; i<params->alphabetSize; i++) {
        rootCharacterProbs[i] = characterProbsHap[0] +
                stSubModel_getSubstitutionProbSlow(params->hetSubModel, i, 0);
        for(int64_t j=1; j<params->alphabetSize; j++) {
            rootCharacterProbs[i] =
                    logAddP(rootCharacterProbs[i],
                            characterProbsHap[j] + stSubModel_getSubstitutionProbSlow(params->hetSubModel, i, j), params->maxNotSumEmissions);
        }
    }
}

double columnIndexLogProbabilitySlow(stRPColumn *column, uint64_t index,
        uint64_t partition, uint64_t *bitCountVectors,
        stRPHmmParameters *params) {
    /*
     * Get the probability of a the characters in a given position within a column for a given partition.
     */
    // Get the sum of log probabilities of the derived characters over the possible source characters
    double rootCharacterProbsHap1[params->alphabetSize];
    columnIndexLogHapProbabilitySlow(column, index,
            partition, bitCountVectors, params, rootCharacterProbsHap1);
    double rootCharacterProbsHap2[params->alphabetSize];
    columnIndexLogHapProbabilitySlow(column, index,
            ~partition, bitCountVectors, params, rootCharacterProbsHap2);

    // Combine the probabilities to calculate the overall probability of a given position in a column
    double logColumnProb = rootCharacterProbsHap1[0] + rootCharacterProbsHap2[0];
    for(int64_t i=1; i<params->alphabetSize; i++) {
        logColumnProb = logAddP(logColumnProb, rootCharacterProbsHap1[i] + rootCharacterProbsHap2[i], params->maxNotSumEmissions);
    }

    return logColumnProb; // + log(1.0/params->alphabetSize);
}

double emissionLogProbabilitySlow(stRPColumn *column,
        stRPCell *cell, uint64_t *bitCountVectors,
        stRPHmmParameters *params) {
    /*
     * Get the log probability of a set of reads for a given column.
     */
    assert(column->length > 0);
    double logPartitionProb = columnIndexLogProbabilitySlow(column, 0,
            cell->partition, bitCountVectors, params);
    for(int64_t i=1; i<column->length; i++) {
        logPartitionProb += columnIndexLogProbabilitySlow(column, i,
                cell->partition, bitCountVectors, params);
    }
    return logPartitionProb;
}

void test_emissionLogProbability(CuTest *testCase) {

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
