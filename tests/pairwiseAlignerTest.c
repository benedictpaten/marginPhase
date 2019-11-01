/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "pairwiseAligner.h"
//#include "multipleAligner.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "randomSequences.h"

void test_diagonal(CuTest *testCase) {
    //Construct an example diagonal.
    int64_t xL = 10, yL = 20, xU = 30, yU = 0; //Coordinates of the upper and lower
    //pairs in x,y coordinates
    Diagonal d = diagonal_construct(xL + yL, xL - yL, xU - yU);
    CuAssertIntEquals(testCase, diagonal_getXay(d), xL + yL);
    CuAssertIntEquals(testCase, diagonal_getMinXmy(d), xL - yL);
    CuAssertIntEquals(testCase, diagonal_getMaxXmy(d), xU - yU);
    CuAssertIntEquals(testCase, diagonal_getWidth(d), (xU - yU - (xL - yL)) / 2 + 1);
    CuAssertIntEquals(testCase, diagonal_getXCoordinate(xL + yL, xL - yL), xL);
    CuAssertIntEquals(testCase, diagonal_getYCoordinate(xL + yL, xL - yL), yL);
    CuAssertTrue(testCase, diagonal_equals(d, d));
    CuAssertTrue(testCase, !diagonal_equals(d, diagonal_construct(0, 0, 0)));
    //A bogus diagonal is one such that |xay + xmy| % 2 != 0 or such that xmyR < xmyL.
    //Try constructing bogus diagonals, should throw exceptions.
    stTry
        {
            diagonal_construct(10, 5, 5);
            CuAssertTrue(testCase, 0);
        }
        stCatch(PAIRWISE_ALIGNMENT_EXCEPTION_ID)
            {
                st_logInfo(stExcept_getMsg(PAIRWISE_ALIGNMENT_EXCEPTION_ID));
            }stTryEnd
    stTry
        {
            diagonal_construct(10, 5, 5);
            CuAssertTrue(testCase, 0);
        }
        stCatch(PAIRWISE_ALIGNMENT_EXCEPTION_ID)
            {
                st_logInfo(stExcept_getMsg(PAIRWISE_ALIGNMENT_EXCEPTION_ID));
            }stTryEnd
    stTry
        {
            diagonal_construct(10, 6, 4);
            CuAssertTrue(testCase, 0);
        }
        stCatch(PAIRWISE_ALIGNMENT_EXCEPTION_ID)
            {
                st_logInfo(stExcept_getMsg(PAIRWISE_ALIGNMENT_EXCEPTION_ID));
            }stTryEnd
}

static bool testDiagonalsEqual(Diagonal d1, Diagonal d2) {
    bool b = diagonal_equals(d1, d2);
    if (!b) {
        st_logCritical("Diagonals not equal: d1: %s, d2: %s \n", diagonal_getString(d1), diagonal_getString(d2));
    }
    return b;
}

void test_bands(CuTest *testCase) {
    stList *anchorPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    ///stList_append(anchorPairs, stIntTuple_construct2( 0, 0));
    stList_append(anchorPairs, stIntTuple_construct2(1, 0));
    stList_append(anchorPairs, stIntTuple_construct2(2, 1));
    stList_append(anchorPairs, stIntTuple_construct2(3, 3));
    /////stList_append(anchorPairs, stIntTuple_construct2( 5, 4));
    //Start the traversal
    int64_t lX = 6, lY = 5;
    Band *band = band_construct(anchorPairs, lX, lY, 2);
    BandIterator *bandIt = bandIterator_construct(band);

    //Forward pass
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(0, 0, 0)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(1, -1, 1)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(2, -2, 2)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(3, -1, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(4, -2, 4)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(5, -1, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(6, -2, 4)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(7, -3, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(8, -2, 2)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(9, -1, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(10, 0, 2)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(11, 1, 1)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(11, 1, 1)));

    //Go backward
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(11, 1, 1)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(10, 0, 2)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(9, -1, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(8, -2, 2)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(7, -3, 3)));
    //Now walk forward a bit
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(7, -3, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(8, -2, 2)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(9, -1, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(10, 0, 2)));
    //Now carry on back again
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(10, 0, 2)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(9, -1, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(8, -2, 2)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(7, -3, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(6, -2, 4)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(5, -1, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(4, -2, 4)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(3, -1, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(2, -2, 2)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(1, -1, 1)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(0, 0, 0)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(0, 0, 0)));

    //Now forward again
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(0, 0, 0)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(1, -1, 1)));
    //Now back again
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(1, -1, 1)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(0, 0, 0)));

    //Cleanup
    bandIterator_destruct(bandIt);
    band_destruct(band);
    stList_destruct(anchorPairs);
}

void test_logAdd(CuTest *testCase) {
    for (int64_t test = 0; test < 100000; test++) {
        double i = st_random();
        double j = st_random();
        double k = i + j;
        double l = exp(logAdd(log(i), log(j)));
        //st_logInfo("I got %f %f\n", k, l);
        CuAssertTrue(testCase, l < k + 0.001);
        CuAssertTrue(testCase, l > k - 0.001);
    }
}

void test_symbol(CuTest *testCase) {
	Alphabet *a = alphabet_constructNucleotide();
    Symbol cA[9] = { 0, 1, 2, 3, 4, 3, 4, 1, 2 };
    SymbolString cA2 = symbolString_construct("AcGTntNCG", 9, a);
    for (int64_t i = 0; i < 9; i++) {
        CuAssertTrue(testCase, cA[i] == cA2.sequence[i]);
    }
    symbolString_destruct(cA2);
    alphabet_destruct(a);
}

void test_cell(CuTest *testCase) {
    StateMachine *sM = stateMachine3_constructNucleotide(threeState);
    double lowerF[sM->stateNumber], middleF[sM->stateNumber], upperF[sM->stateNumber], currentF[sM->stateNumber];
    double lowerB[sM->stateNumber], middleB[sM->stateNumber], upperB[sM->stateNumber], currentB[sM->stateNumber];
    for (int64_t i = 0; i < sM->stateNumber; i++) {
        middleF[i] = sM->startStateProb(sM, i);
        middleB[i] = LOG_ZERO;
        lowerF[i] = LOG_ZERO;
        lowerB[i] = LOG_ZERO;
        upperF[i] = LOG_ZERO;
        upperB[i] = LOG_ZERO;
        currentF[i] = LOG_ZERO;
        currentB[i] = sM->endStateProb(sM, i);
    }
    Symbol cX = 0, cY = 3;  // A and T
    //Do forward
    cell_calculateForward(sM, lowerF, NULL, NULL, middleF, cX, cY, NULL);
    cell_calculateForward(sM, upperF, middleF, NULL, NULL, cX, cY, NULL);
    cell_calculateForward(sM, currentF, lowerF, middleF, upperF, cX, cY, NULL);
    //Do backward
    cell_calculateBackward(sM, currentB, lowerB, middleB, upperB, cX, cY, NULL);
    cell_calculateBackward(sM, upperB, middleB, NULL, NULL, cX, cY, NULL);
    cell_calculateBackward(sM, lowerB, NULL, NULL, middleB, cX, cY, NULL);
    double totalProbForward = cell_dotProduct2(currentF, sM, sM->endStateProb);
    double totalProbBackward = cell_dotProduct2(middleB, sM, sM->startStateProb);
    st_logInfo("Total probability for cell test, forward %f and backward %f\n", totalProbForward, totalProbBackward);
    CuAssertDblEquals(testCase, totalProbForward, totalProbBackward, 0.00001); //Check the forward and back probabilities are about equal
    stateMachine_destruct(sM);
}

void test_dpDiagonal(CuTest *testCase) {
    StateMachine *sM = stateMachine3_constructNucleotide(threeState);
    Diagonal diagonal = diagonal_construct(3, -1, 1);
    DpDiagonal *dpDiagonal = dpDiagonal_construct(diagonal, sM->stateNumber);

    //Get cell
    double *c1 = dpDiagonal_getCell(dpDiagonal, -1);
    CuAssertTrue(testCase, c1 != NULL);
    double *c2 = dpDiagonal_getCell(dpDiagonal, 1);
    CuAssertTrue(testCase, c2 != NULL);
    CuAssertTrue(testCase, dpDiagonal_getCell(dpDiagonal, 3) == NULL);
    CuAssertTrue(testCase, dpDiagonal_getCell(dpDiagonal, -3) == NULL);

    dpDiagonal_initialiseValues(dpDiagonal, sM, sM->endStateProb); //Test initialise values
    double totalProb = LOG_ZERO;
    for (int64_t i = 0; i < sM->stateNumber; i++) {
        CuAssertDblEquals(testCase, c1[i], sM->endStateProb(sM, i), 0.0);
        CuAssertDblEquals(testCase, c2[i], sM->endStateProb(sM, i), 0.0);
        totalProb = logAdd(totalProb, 2 * c1[i]);
        totalProb = logAdd(totalProb, 2 * c2[i]);
    }

    DpDiagonal *dpDiagonal2 = dpDiagonal_clone(dpDiagonal);
    CuAssertTrue(testCase, dpDiagonal_equals(dpDiagonal, dpDiagonal2));

    CuAssertDblEquals(testCase, totalProb, dpDiagonal_dotProduct(dpDiagonal, dpDiagonal2), 0.001); //Check it runs

    dpDiagonal_destruct(dpDiagonal);
    dpDiagonal_destruct(dpDiagonal2);
}

void test_dpMatrix(CuTest *testCase) {
    int64_t lX = 3, lY = 2;
    DpMatrix *dpMatrix = dpMatrix_construct(lX + lY, 5);

    CuAssertIntEquals(testCase, dpMatrix_getActiveDiagonalNumber(dpMatrix), 0);

    for (int64_t i = -1; i <= lX + lY + 10; i++) {
        CuAssertTrue(testCase, dpMatrix_getDiagonal(dpMatrix, i) == NULL);
    }

    for (int64_t i = 0; i <= lX + lY; i++) {
        DpDiagonal *dpDiagonal = dpMatrix_createDiagonal(dpMatrix, diagonal_construct(i, -i, i));
        CuAssertTrue(testCase, dpDiagonal == dpMatrix_getDiagonal(dpMatrix, i));
        CuAssertIntEquals(testCase, dpMatrix_getActiveDiagonalNumber(dpMatrix), i + 1);
    }

    for (int64_t i = lX + lY; i >= 0; i--) {
        dpMatrix_deleteDiagonal(dpMatrix, i);
        CuAssertTrue(testCase, dpMatrix_getDiagonal(dpMatrix, i) == NULL);
        CuAssertIntEquals(testCase, dpMatrix_getActiveDiagonalNumber(dpMatrix), i);
    }

    CuAssertIntEquals(testCase, dpMatrix_getActiveDiagonalNumber(dpMatrix), 0);

    dpMatrix_destruct(dpMatrix);
}

void test_diagonalDPCalculations(CuTest *testCase) {
    //Sets up a complete matrix for the following example and checks the total marginal
    //probability and the posterior probabilities of the matches

	Alphabet *a = alphabet_constructNucleotide();
    const char *sX = "AGCG";
    const char *sY = "AGTTCG";
    int64_t lX = strlen(sX);
    int64_t lY = strlen(sY);
    SymbolString sX2 = symbolString_construct(sX, lX, a);
    SymbolString sY2 = symbolString_construct(sY, lY, a);
    StateMachine *sM = stateMachine3_constructNucleotide(threeState);
    DpMatrix *dpMatrixForward = dpMatrix_construct(lX + lY, sM->stateNumber);
    DpMatrix *dpMatrixBackward = dpMatrix_construct(lX + lY, sM->stateNumber);
    stList *anchorPairs = stList_construct();
    Band *band = band_construct(anchorPairs, lX, lY, 2);
    BandIterator *bandIt = bandIterator_construct(band);

    //Initialise matrices
    for (int64_t i = 0; i <= lX + lY; i++) {
        Diagonal d = bandIterator_getNext(bandIt);
        //initialisation
        dpDiagonal_zeroValues(dpMatrix_createDiagonal(dpMatrixBackward, d));
        dpDiagonal_zeroValues(dpMatrix_createDiagonal(dpMatrixForward, d));
    }
    dpDiagonal_initialiseValues(dpMatrix_getDiagonal(dpMatrixForward, 0), sM, sM->startStateProb);
    dpDiagonal_initialiseValues(dpMatrix_getDiagonal(dpMatrixBackward, lX + lY), sM, sM->endStateProb);

    //Forward algorithm
    for (int64_t i = 1; i <= lX + lY; i++) {
        //Do the forward calculation
        diagonalCalculationForward(sM, i, dpMatrixForward, sX2, sY2);
    }

    //Backward algorithm
    for (int64_t i = lX + lY; i > 0; i--) {
        //Do the backward calculation
        diagonalCalculationBackward(sM, i, dpMatrixBackward, sX2, sY2);
    }

    //Calculate total probabilities
    double totalProbForward = cell_dotProduct2(
            dpDiagonal_getCell(dpMatrix_getDiagonal(dpMatrixForward, lX + lY), lX - lY), sM, sM->endStateProb);
    double totalProbBackward = cell_dotProduct2(dpDiagonal_getCell(dpMatrix_getDiagonal(dpMatrixBackward, 0), 0), sM,
            sM->startStateProb);
    st_logInfo("Total forward and backward prob %f %f\n", (float) totalProbForward, (float) totalProbBackward);

    CuAssertDblEquals(testCase, totalProbForward, totalProbBackward, 0.001); //Check the forward and back probabilities are about equal

    //Test calculating the posterior probabilities along the diagonals of the matrix.
    for (int64_t i = 0; i <= lX + lY; i++) {
        //Calculate the total probs
        double totalDiagonalProb = diagonalCalculationTotalProbability(sM, i, dpMatrixForward, dpMatrixBackward, sX2,
                sY2);
        CuAssertDblEquals(testCase, totalProbForward, totalDiagonalProb, 0.01); //Check the forward and back probabilities are about equal
    }

    //Now do the posterior probabilities
    stList *alignedPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    void *extraArgs[1] = { alignedPairs };
    for (int64_t i = 1; i <= lX + lY; i++) {
        PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
        p->threshold = 0.2;
        diagonalCalculationPosteriorMatchProbs(sM, i, dpMatrixForward, dpMatrixBackward, sX2, sY2, totalProbForward, p,
                extraArgs);
        pairwiseAlignmentBandingParameters_destruct(p);
    }

    stSortedSet *alignedPairsSet = stSortedSet_construct3((int (*)(const void *, const void *)) stIntTuple_cmpFn,
            (void (*)(void *)) stIntTuple_destruct);
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(0, 0));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(1, 1));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(2, 4));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(3, 5));

    for (int64_t i = 0; i < stList_length(alignedPairs); i++) {
        stIntTuple *pair = stList_get(alignedPairs, i);
        int64_t x = stIntTuple_get(pair, 1), y = stIntTuple_get(pair, 2);
        st_logInfo("Pair %f %" PRIi64 " %" PRIi64 "\n", (float) stIntTuple_get(pair, 0) / PAIR_ALIGNMENT_PROB_1, x, y);
        CuAssertTrue(testCase, stSortedSet_search(alignedPairsSet, stIntTuple_construct2( x, y)) != NULL);
    }
    CuAssertIntEquals(testCase, stList_length(alignedPairs), 4);

}

stList *getRandomAnchorPairs(int64_t lX, int64_t lY) {
    stList *anchorPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    int64_t x = -1;
    int64_t y = -1;
    while (1) {
        x += st_randomInt(1, 20);
        y += st_randomInt(1, 20);
        int64_t expansion = 2 * st_randomInt(0, 5);
        if (x >= lX || y >= lY) {
            break;
        }
        assert(x >= 0 && x < lX);
        assert(y >= 0 && y < lY);
        stList_append(anchorPairs, stIntTuple_construct3(x, y, expansion));
    }
    return anchorPairs;
}

static void checkAlignedPairs(CuTest *testCase, stList *blastPairs, int64_t lX, int64_t lY, bool xGaps, bool yGaps) {
    st_logInfo("I got %" PRIi64 " pairs to check\n", stList_length(blastPairs));
    stSortedSet *pairs = stSortedSet_construct3((int (*)(const void *, const void *)) stIntTuple_cmpFn,
            (void (*)(void *)) stIntTuple_destruct);
    for (int64_t i = 0; i < stList_length(blastPairs); i++) {
        stIntTuple *j = stList_get(blastPairs, i);
        CuAssertTrue(testCase, stIntTuple_length(j) == 3);

        int64_t x = stIntTuple_get(j, 1);
        int64_t y = stIntTuple_get(j, 2);
        int64_t score = stIntTuple_get(j, 0);
        CuAssertTrue(testCase, score > 0);
        CuAssertTrue(testCase, score <= PAIR_ALIGNMENT_PROB_1);

        if(yGaps) {
        	CuAssertTrue(testCase, x >= -1);
        }
        else {
        	CuAssertTrue(testCase, x >= 0);
        }

        if(xGaps) {
        	CuAssertTrue(testCase, y >= -1);
        }
        else {
        	CuAssertTrue(testCase, y >= 0);
        }

        CuAssertTrue(testCase, x < lX);
        CuAssertTrue(testCase, y < lY);

        //Check is unique
        stIntTuple *pair = stIntTuple_construct2(x, y);
        CuAssertTrue(testCase, stSortedSet_search(pairs, pair) == NULL);
        stSortedSet_insert(pairs, pair);
    }
    stSortedSet_destruct(pairs);
}

static void checkAlignedPairsAreOrdered(CuTest *testCase, stList *alignedPairs) {
	/*
	 * Checks that the list of aligned pairs are ordered in sequence space so that no aligned pair i
	 * proceeds another j in the sequence such that i's coordinates in the two sequences both precede j's.
	 */
	for(int64_t i=0; i<stList_length(alignedPairs); i++) {
		stIntTuple *aPair = stList_get(alignedPairs, i);
		int64_t x = stIntTuple_get(aPair, 1);
		int64_t y = stIntTuple_get(aPair, 2);

		for(int64_t j=0; j<i; j++) {
			stIntTuple *pPair = stList_get(alignedPairs, j);
			int64_t x2 = stIntTuple_get(pPair, 1);
			int64_t y2 = stIntTuple_get(pPair, 2);

			CuAssertTrue(testCase, x2 < x || y2 < y);
		}
	}
}

void test_getAlignedPairsWithBanding(CuTest *testCase) {
	Alphabet *a = alphabet_constructNucleotide();
    for (int64_t test = 0; test < 100; test++) {
        //Make a pair of sequences
        char *sX = getRandomSequence(st_randomInt(0, 100));
        char *sY = evolveSequence(sX); //stString_copy(seqX);
        int64_t lX = strlen(sX);
        int64_t lY = strlen(sY);
        st_logInfo("Sequence X to align: %s END\n", sX);
        st_logInfo("Sequence Y to align: %s END\n", sY);
        SymbolString sX2 = symbolString_construct(sX, lX, a);
        SymbolString sY2 = symbolString_construct(sY, lY, a);
        //Now do alignment
        PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
        p->traceBackDiagonals = st_randomInt(1, 10);
        p->minDiagsBetweenTraceBack = p->traceBackDiagonals + st_randomInt(2, 10);
        p->diagonalExpansion = st_randomInt(0, 10) * 2;
        p->dynamicAnchorExpansion = st_random() > 0.5;
        StateMachine *sM = stateMachine3_constructNucleotide(threeState);
        stList *anchorPairs = getRandomAnchorPairs(lX, lY);

        stList *alignedPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
        void *extraArgs[1] = { alignedPairs };
        getPosteriorProbsWithBanding(sM, anchorPairs, sX2, sY2, p, 0, 0, diagonalCalculationPosteriorMatchProbs,
                extraArgs);
        //Check the aligned pairs.
        checkAlignedPairs(testCase, alignedPairs, lX, lY, 0, 0);

        //Cleanup
        stateMachine_destruct(sM);
        free(sX);
        free(sY);
        free(sX2.sequence);
        free(sY2.sequence);
        stList_destruct(alignedPairs);
    }
    alphabet_destruct(a);
}

static void checkBlastPairs(CuTest *testCase, stList *blastPairs, int64_t lX, int64_t lY, int64_t diagonalExpansion, bool checkNonOverlapping) {
    st_logInfo("I got %" PRIi64 " pairs to check\n", stList_length(blastPairs));
    int64_t pX = -1;
    int64_t pY = -1;
    for (int64_t i = 0; i < stList_length(blastPairs); i++) {
        stIntTuple *j = stList_get(blastPairs, i);
        CuAssertTrue(testCase, stIntTuple_length(j) == 3);

        int64_t x = stIntTuple_get(j, 0);
        int64_t y = stIntTuple_get(j, 1);
        int64_t expansion = stIntTuple_get(j, 2);

        CuAssertTrue(testCase, x >= 0);
        CuAssertTrue(testCase, y >= 0);
        CuAssertTrue(testCase, x < lX);
        CuAssertTrue(testCase, y < lY);
        if (checkNonOverlapping) {
            CuAssertTrue(testCase, x > pX);
            CuAssertTrue(testCase, y > pY);
        }
        pX = x;
        pY = y;

        CuAssertTrue(testCase, expansion == diagonalExpansion);
        CuAssertTrue(testCase, expansion >= 0);
        CuAssertTrue(testCase, expansion % 2 == 0);

    }
}

void test_getSplitPoints(CuTest *testCase) {
    int64_t matrixSize = 2000 * 2000;

    stList *anchorPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);

    //Test a small region, which produces no splits
    int64_t lX = 3000;
    int64_t lY = 1000;
    stList *splitPoints = getSplitPoints(anchorPairs, lX, lY, matrixSize, 0, 0);
    CuAssertIntEquals(testCase, 1, stList_length(splitPoints));
    CuAssertTrue(testCase, stIntTuple_equalsFn(stList_get(splitPoints, 0), stIntTuple_construct4(0, 0, lX, lY)));
    stList_destruct(splitPoints);

    //Test with one really big matrix with no anchors
    lX = 20000;
    lY = 25000;
    splitPoints = getSplitPoints(anchorPairs, lX, lY, matrixSize, 1, 1);
    CuAssertIntEquals(testCase, 0, stList_length(splitPoints));
    stList_destruct(splitPoints);

    splitPoints = getSplitPoints(anchorPairs, lX, lY, matrixSize, 1, 0);
    CuAssertIntEquals(testCase, 1, stList_length(splitPoints));
    CuAssertTrue(testCase,
            stIntTuple_equalsFn(stList_get(splitPoints, 0), stIntTuple_construct4(18000, 23000, lX, lY)));
    stList_destruct(splitPoints);

    splitPoints = getSplitPoints(anchorPairs, lX, lY, matrixSize, 0, 1);
    CuAssertIntEquals(testCase, 1, stList_length(splitPoints));
    CuAssertTrue(testCase, stIntTuple_equalsFn(stList_get(splitPoints, 0), stIntTuple_construct4(0, 0, 2000, 2000)));
    stList_destruct(splitPoints);

    splitPoints = getSplitPoints(anchorPairs, lX, lY, matrixSize, 0, 0);
    CuAssertIntEquals(testCase, 2, stList_length(splitPoints));
    CuAssertTrue(testCase, stIntTuple_equalsFn(stList_get(splitPoints, 0), stIntTuple_construct4(0, 0, 2000, 2000)));
    CuAssertTrue(testCase,
            stIntTuple_equalsFn(stList_get(splitPoints, 1), stIntTuple_construct4(18000, 23000, lX, lY)));
    stList_destruct(splitPoints);

    //Now test with some more points
    stList_append(anchorPairs, stIntTuple_construct2(2000, 2000)); //This should not create a split
    stList_append(anchorPairs, stIntTuple_construct2(4002, 4001)); //This should cause a split
    stList_append(anchorPairs, stIntTuple_construct2(5000, 5000)); //This should not cause a split
    stList_append(anchorPairs, stIntTuple_construct2(8000, 6000)); //Neither should this (it is maximum sized)
    stList_append(anchorPairs, stIntTuple_construct2(9000, 9000)); //Or this (it is maximum sized)
    stList_append(anchorPairs, stIntTuple_construct2(10000, 14000)); //This should create a split
    stList_append(anchorPairs, stIntTuple_construct2(15000, 15000)); //This should also create a split
    stList_append(anchorPairs, stIntTuple_construct2(16000, 16000)); //This should not, but there will be a split with the end.

    splitPoints = getSplitPoints(anchorPairs, lX, lY, matrixSize, 0, 0);

    for (int64_t i = 0; i < stList_length(splitPoints); i++) {
        stIntTuple *j = stList_get(splitPoints, i);
        st_logInfo("I got split point: x1: %" PRIi64 " y1: %" PRIi64 " x2: %" PRIi64 " y2: %" PRIi64 "\n",
                stIntTuple_get(j, 0), stIntTuple_get(j, 1), stIntTuple_get(j, 2), stIntTuple_get(j, 3));
    }

    CuAssertIntEquals(testCase, 5, stList_length(splitPoints));
    CuAssertTrue(testCase, stIntTuple_equalsFn(stList_get(splitPoints, 0), stIntTuple_construct4(0, 0, 3001, 3001)));
    CuAssertTrue(testCase,
            stIntTuple_equalsFn(stList_get(splitPoints, 1), stIntTuple_construct4(3002, 3001, 9500, 11001)));
    CuAssertTrue(testCase,
            stIntTuple_equalsFn(stList_get(splitPoints, 2), stIntTuple_construct4(9501, 12000, 12001, 14500)));
    CuAssertTrue(testCase,
            stIntTuple_equalsFn(stList_get(splitPoints, 3), stIntTuple_construct4(13000, 14501, 18000, 18001)));
    CuAssertTrue(testCase,
            stIntTuple_equalsFn(stList_get(splitPoints, 4), stIntTuple_construct4(18001, 23000, 20000, 25000)));

    stList_destruct(splitPoints);
    stList_destruct(anchorPairs);
}

void test_getAlignedPairs(CuTest *testCase) {
    for (int64_t test = 0; test < 100; test++) {
        //Make a pair of sequences
        char *sX = getRandomSequence(st_randomInt(0, 100));
        char *sY = evolveSequence(sX); //stString_copy(seqX);
        int64_t lX = strlen(sX);
        int64_t lY = strlen(sY);
        st_logInfo("Sequence X to align: %s END\n", sX);
        st_logInfo("Sequence Y to align: %s END\n", sY);

        //Now do alignment
        PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
        StateMachine *sM = stateMachine3_constructNucleotide(threeState);

        stList *alignedPairs = getAlignedPairs(sM, sX, sY, p, 0, 0);

        //Check the aligned pairs.
        checkAlignedPairs(testCase, alignedPairs, lX, lY, 0, 0);

        //Cleanup
        stateMachine_destruct(sM);
        free(sX);
        free(sY);
        stList_destruct(alignedPairs);
    }
}

void test_getAlignedPairsWithRaggedEnds(CuTest *testCase) {
    for (int64_t test = 0; test < 1000; test++) {
        //Make a pair of sequences
        int64_t coreLength = 100, randomPortionLength = 100;
        char *sX = getRandomSequence(coreLength);
        char *randomPrefix = getRandomSequence(randomPortionLength);
        char *randomSuffix = getRandomSequence(randomPortionLength);
        char *sY = stString_print("%s%s%s", randomPrefix, sX, randomSuffix); //x with an extra bit at the end.

        st_logInfo("Sequence X to align: %s END\n", sX);
        st_logInfo("Sequence Y to align: %s END\n", sY);

        //Now do alignment
        PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
        StateMachine *sM = stateMachine3_constructNucleotide(threeState);
        stList *alignedPairs = getAlignedPairs(sM, sX, sY, p, 1, 1);
        alignedPairs = filterPairwiseAlignmentToMakePairsOrdered(alignedPairs, sX, sY, p);

        //Check the aligned pairs.
        checkAlignedPairs(testCase, alignedPairs, strlen(sX), strlen(sY), 0, 0);
        CuAssertIntEquals(testCase, stList_length(alignedPairs), coreLength);
        for (int64_t i = 0; i < stList_length(alignedPairs); i++) {
            stIntTuple *j = stList_get(alignedPairs, i);
            CuAssertTrue(testCase, stIntTuple_length(j) == 3);

            int64_t x = stIntTuple_get(j, 1);
            int64_t y = stIntTuple_get(j, 2);
            CuAssertIntEquals(testCase, x + randomPortionLength, y);
        }

        //Cleanup
        stateMachine_destruct(sM);
        free(sX);
        free(sY);
        free(randomPrefix);
        free(randomSuffix);
        stList_destruct(alignedPairs);
        pairwiseAlignmentBandingParameters_destruct(p);
    }
}

/*
 * Test indel posterior prob calculating methods
 */

int64_t getGapScore(int64_t *gapProbsX, int64_t *gapProbsY, int64_t xS, int64_t yS, int64_t xE, int64_t yE,
		PairwiseAlignmentParameters *p) {
	double s = 0;
	for(int64_t i=xS+1; i<xE; i++) {
		s += gapProbsX[i];
	}
	for(int64_t i=yS+1; i<yE; i++) {
		s += gapProbsY[i];
	}
	return s * p->gapGamma;
}

static stList *getMEAlignment(stList *alignedPairs, int64_t seqLengthX, int64_t seqLengthY,
		stList *gapXPairs, stList *gapYPairs, PairwiseAlignmentParameters *p, double *alignmentScore) {
	double *scores = st_calloc(stList_length(alignedPairs), sizeof(double));
	int64_t *backPointers = st_calloc(stList_length(alignedPairs), sizeof(int64_t));

	int64_t *gapProbsX = st_calloc(seqLengthX, sizeof(int64_t));
	for(int64_t i=0; i<stList_length(gapXPairs); i++) {
		stIntTuple *gPair = stList_get(gapXPairs, i);
		gapProbsX[stIntTuple_get(gPair, 1)] += stIntTuple_get(gPair, 0);
	}

	int64_t *gapProbsY = st_calloc(seqLengthY, sizeof(int64_t));
	for(int64_t i=0; i<stList_length(gapYPairs); i++) {
		stIntTuple *gPair = stList_get(gapYPairs, i);
		gapProbsY[stIntTuple_get(gPair, 2)] += stIntTuple_get(gPair, 0);
	}

	for(int64_t i=0; i<stList_length(alignedPairs); i++) {
		stIntTuple *aPair = stList_get(alignedPairs, i);
		int64_t matchProb = stIntTuple_get(aPair, 0);
		int64_t x = stIntTuple_get(aPair, 1);
		int64_t y = stIntTuple_get(aPair, 2);
		scores[i] = matchProb + getGapScore(gapProbsX, gapProbsY, -1, -1, x, y, p);
		backPointers[i] = -1;

		for(int64_t j=0; j<stList_length(alignedPairs); j++) {
			stIntTuple *pPair = stList_get(alignedPairs, j);
			int64_t x2 = stIntTuple_get(pPair, 1);
			int64_t y2 = stIntTuple_get(pPair, 2);

			if(x2 < x && y2 < y) {
				double s = scores[j] + matchProb + getGapScore(gapProbsX, gapProbsY, x2, y2, x, y, p);

				if(s > scores[i]) {
					scores[i] = s;
					backPointers[i] = j;
				}
			}
		}
	}

	double maxScore = getGapScore(gapProbsX, gapProbsY, -1, -1, seqLengthX, seqLengthY, p);
	int64_t maxScoreIndex = -1;
	for(int64_t i=0; i<stList_length(alignedPairs); i++) {
		stIntTuple *aPair = stList_get(alignedPairs, i);
		int64_t x = stIntTuple_get(aPair, 1);
		int64_t y = stIntTuple_get(aPair, 2);

		double s = scores[i] + getGapScore(gapProbsX, gapProbsY, x, y, seqLengthX, seqLengthY, p);
		if(s > maxScore) {
			maxScore = s;
			maxScoreIndex = i;
		}
	}

	stList *filteredAlignment = stList_construct();
	while(maxScoreIndex != -1) {
		stList_append(filteredAlignment, stList_get(alignedPairs, maxScoreIndex));
		maxScoreIndex = backPointers[maxScoreIndex];
	}
	stList_reverse(filteredAlignment);

	free(scores);
	free(backPointers);
	free(gapProbsX);
	free(gapProbsY);

	*alignmentScore = maxScore;
	return filteredAlignment;
}

void checkAlignedPairsAreTotallyOrdered(CuTest *testCase, stList *alignedPairs, int64_t lX, int64_t lY) {
	int64_t x = -1, y = -1;
	for(int64_t i=0; i<stList_length(alignedPairs); i++) {
		stIntTuple *aPair = stList_get(alignedPairs, i);
		int64_t x2 = stIntTuple_get(aPair, 1);
		int64_t y2 = stIntTuple_get(aPair, 2);
		int64_t score = stIntTuple_get(aPair, 0);

		CuAssertTrue(testCase, x < x2);
		CuAssertTrue(testCase, y < y2);

		x = x2;
		y = y2;

		CuAssertTrue(testCase, score >= 0);
		CuAssertTrue(testCase, score <= PAIR_ALIGNMENT_PROB_1);
	}
	CuAssertTrue(testCase, x < lX);
	CuAssertTrue(testCase, y < lY);
}

/*void printAlignment(stList *alignedPairs, char *sX, char *sY) {
	stList *seqX = stList_construct();
	stList *seqY = stList_construct();
	int64_t x = -1, y = -1;
	for(int64_t i=0; i<stList_length(alignedPairs); i++) {
		stIntTuple *aPair = stList_get(alignedPairs, i);
		int64_t x2 = stIntTuple_get(aPair, 1);
		int64_t y2 = stIntTuple_get(aPair, 2);
		while(x2 - x > 1) {
			stList_append(seqX, stString_print("%c", sX[++x]));
			stList_append(seqY, stString_print("-"));
		}
		while(y2 - y > 1) {
			stList_append(seqY, stString_print("%c", sY[++y]));
			stList_append(seqX, stString_print("-"));
		}
		stList_append(seqX, stString_print("%c", sX[x2]));
		stList_append(seqY, stString_print("%c", sY[y2]));
		x = x2;
		y = y2;
	}
	fprintf(stderr, "SeaX: %s\n", stString_join2("", seqX));
	fprintf(stderr, "SeaY: %s\n", stString_join2("", seqY));
}*/

void checkAlignmentNoFurtherLeftShifts(CuTest *testCase, stList *alignedPairs, char *sX, char *sY) {
	int64_t x = -1, y = -1;
	for(int64_t i=0; i<stList_length(alignedPairs); i++) {
		stIntTuple *aPair = stList_get(alignedPairs, 0);
		int64_t x2 = stIntTuple_get(aPair, 1);
		int64_t y2 = stIntTuple_get(aPair, 2);
		if(x2 - x > 1 && y2 > 0) { // Possible shift
			CuAssertTrue(testCase, toupper(sY[y2-1]) != toupper(sX[x2-1])); // Check no shift is possible
		}
		if(y2 - y > 1 && x2 > 0) { // Possible shift
			CuAssertTrue(testCase, toupper(sY[y2-1]) != toupper(sX[x2-1])); // Check no shift is possible
		}
		x = x2;
		y = y2;
	}
}

void test_getAlignedPairsWithIndels(CuTest *testCase) {
    for (int64_t test = 0; test < 1000; test++) {
        //Make a pair of sequences
        char *sX = getRandomSequence(st_randomInt(0, 100));
        char *sY = evolveSequence(sX); //stString_copy(seqX);

        int64_t lX = strlen(sX);
        int64_t lY = strlen(sY);
        st_logInfo("Sequence X to align: %s END\n", sX);
        st_logInfo("Sequence Y to align: %s END\n", sY);

        //Now do alignment
        PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
        StateMachine *sM = stateMachine3_constructNucleotide(threeState);

        stList *alignedPairs = NULL, *gapXPairs = NULL, *gapYPairs = NULL;

        getAlignedPairsWithIndels(sM, sX, sY, p, &alignedPairs, &gapXPairs, &gapYPairs, st_random() > 0.5, st_random() > 0.5);

        //Check the aligned pairs.
        checkAlignedPairs(testCase, alignedPairs, lX, lY, 0, 0);

        // Check the inserts in Y
        checkAlignedPairs(testCase, gapXPairs, lX, lY, 1, 0);

        // Check the inserts in X
        checkAlignedPairs(testCase, gapYPairs, lX, lY, 0, 1);

        checkAlignedPairsAreOrdered(testCase, alignedPairs);
        checkAlignedPairsAreOrdered(testCase, gapXPairs);
        checkAlignedPairsAreOrdered(testCase, gapYPairs);

        // Get the MEA alignment, setting gap gamma to 0 so that our match based MEA algorithm is equivalent
        double alignmentScore;
        stList *filteredAlignment = getMaximalExpectedAccuracyPairwiseAlignment(alignedPairs, gapXPairs, gapYPairs, lX, lY, &alignmentScore, p);

        // Get bad algorithm MEA
        double expectedAlignmentScore;
        stList *filteredAlignment2 = getMEAlignment(alignedPairs, lX, lY, gapXPairs, gapYPairs, p, &expectedAlignmentScore);

        st_logInfo("Started with : %i pairs, ended with %i aligned pairs, expected %i aligned pairs\n",
        		(int)stList_length(alignedPairs), (int)stList_length(filteredAlignment), (int)stList_length(filteredAlignment2));

        checkAlignedPairs(testCase, filteredAlignment, lX, lY, 0, 0);
        checkAlignedPairsAreTotallyOrdered(testCase, filteredAlignment, lX, lY);

        st_logInfo("Scores: %f %f\n", expectedAlignmentScore/PAIR_ALIGNMENT_PROB_1, alignmentScore/PAIR_ALIGNMENT_PROB_1);

        CuAssertDblEquals(testCase, expectedAlignmentScore/PAIR_ALIGNMENT_PROB_1, alignmentScore/PAIR_ALIGNMENT_PROB_1, 0.0001);
        CuAssertIntEquals(testCase, stList_length(filteredAlignment2), stList_length(filteredAlignment));

        // Now do left shift alignment
        stList *shiftedAlignedPairs = leftShiftAlignment(filteredAlignment, sX, sY);

        //printAlignment(filteredAlignment, sX, sY);
        //fprintf(stderr, "After shift\n");
        //printAlignment(shiftedAlignedPairs, sX, sY);

        CuAssertTrue(testCase, stList_length(shiftedAlignedPairs) >= stList_length(filteredAlignment));
        checkAlignedPairs(testCase, shiftedAlignedPairs, lX, lY, 0, 0);
        checkAlignedPairsAreTotallyOrdered(testCase, shiftedAlignedPairs, lX, lY);
        checkAlignmentNoFurtherLeftShifts(testCase, shiftedAlignedPairs, sX, sY);

        //Cleanup
        stateMachine_destruct(sM);
        free(sX);
        free(sY);
        stList_destruct(alignedPairs);
        stList_destruct(gapXPairs);
        stList_destruct(gapYPairs);
        stList_destruct(filteredAlignment);
        stList_destruct(filteredAlignment2);
        stList_destruct(shiftedAlignedPairs);
        pairwiseAlignmentBandingParameters_destruct(p);
    }
}

void test_leftShiftAlignment(CuTest *testCase) {
	stList *alignedPairs = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
	char *seqX = "GATTTACATC";
	char *seqY = "GATTACAATCTG";

	// Make following alignment
	// 01234567-8--9
	// GATTTACA-T--C
	// GATT-ACAATCTG
	// 0123-45678901

	// expect the following output
	// 01234567---89
	// GATTTACA---TC
	// GA-TTACAATCTG
	// 01-2345678901

	int64_t alignedPairsX[] = { 0, 1, 2, 3, 5, 6, 7, 8, 9 };
	int64_t alignedPairsY[] = { 0, 1, 2, 3, 4, 5, 6, 8, 11 };

	for(int64_t i=0; i<9; i++) {
		stList_append(alignedPairs, stIntTuple_construct3(1, alignedPairsX[i], alignedPairsY[i]));
	}

	// Run left shift
	stList *leftShiftedAlignment = leftShiftAlignment(alignedPairs, seqX, seqY);

	//for(int64_t i=0; i<stList_length(leftShiftedAlignment); i++) {
	//	stIntTuple *aPair = stList_get(leftShiftedAlignment, i);
	//	fprintf(stderr, "FUC %i %i\n", (int)stIntTuple_get(aPair, 1), (int)stIntTuple_get(aPair, 2));
	//}

	// Test we get what we expect
	CuAssertIntEquals(testCase, 9, stList_length(leftShiftedAlignment));

	int64_t shiftedAlignedPairsX[] = { 0, 1, 3, 4, 5, 6, 7, 8, 9 };
	int64_t shiftedAlignedPairsY[] = { 0, 1, 2, 3, 4, 5, 6, 10, 11 };

	for(int64_t i=0; i<9; i++) {
		stIntTuple *aPair = stList_get(leftShiftedAlignment, i);
		CuAssertIntEquals(testCase, shiftedAlignedPairsX[i], stIntTuple_get(aPair, 1));
		CuAssertIntEquals(testCase, shiftedAlignedPairsY[i], stIntTuple_get(aPair, 2));
	}

	// Cleanup
	stList_destruct(alignedPairs);
	stList_destruct(leftShiftedAlignment);
}

/*
 * EM training tests.
 */

void test_hmm(CuTest *testCase, StateMachineType stateMachineType, EmissionType emissionType, int64_t alphabet_size) {
    //Expectation object
    Hmm *hmm = hmm_constructEmpty(0.0, stateMachineType, emissionType);

    //Add some transition expectations
    for (int64_t from = 0; from < hmm->stateNumber; from++) {
        for (int64_t to = 0; to < hmm->stateNumber; to++) {
            hmm_addToTransitionExpectation(hmm, from, to, from * hmm->stateNumber + to);
        }
    }


    //Add some emission expectations
    for (int64_t state = 0; state < hmm->stateNumber; state++) {
        for (int64_t x = 0; x < hmm->emissionNoPerState; x++) {
        	hmm_addToEmissionsExpectation(hmm, state, x, state * hmm->emissionNoPerState + x);
        }
    }

    /*
    //Write to a file
    char *tempFile = stString_print("./temp%" PRIi64 ".hmm", st_randomInt(0, INT64_MAX));
    CuAssertTrue(testCase, !stFile_exists(tempFile)); //Quick check that we don't write over anything.
    FILE *fH = fopen(tempFile, "w");
    hmm_write(hmm, fH);
    fclose(fH);
    hmm_destruct(hmm);

    //Load from a file
    hmm = hmm_loadFromFile(tempFile);
    stFile_rmrf(tempFile);*/

    //Check the transition expectations
    for (int64_t from = 0; from < hmm->stateNumber; from++) {
        for (int64_t to = 0; to < hmm->stateNumber; to++) {
            CuAssertTrue(testCase, hmm_getTransition(hmm, from, to) == from * hmm->stateNumber + to);
        }
    }

    //Check the emission expectations
    for (int64_t state = 0; state < hmm->stateNumber; state++) {
    	for (int64_t x = 0; x < hmm->emissionNoPerState; x++) {
    		CuAssertTrue(testCase, hmm_getEmissionsExpectation(hmm, state, x) == state * hmm->emissionNoPerState + x);
        }
    }

    //Normalise
    hmm_normalise(hmm);

    //Recheck transitions
    for (int64_t from = 0; from < hmm->stateNumber; from++) {
        for (int64_t to = 0; to < hmm->stateNumber; to++) {
            double z = from * hmm->stateNumber * hmm->stateNumber + (hmm->stateNumber * (hmm->stateNumber - 1)) / 2;
            CuAssertDblEquals(testCase, (from * hmm->stateNumber + to) / z, hmm_getTransition(hmm, from, to), 0.0);
        }
    }

    //Recheck the emissions
    for (int64_t state = 0; state < hmm->stateNumber; state++) {
        for (int64_t x = 0; x < hmm->emissionNoPerState; x++) {
			double z = hmm->emissionNoPerState * hmm->emissionNoPerState * state
					+ (hmm->emissionNoPerState * (hmm->emissionNoPerState - 1))
							/ 2;
			CuAssertTrue(testCase,
					hmm_getEmissionsExpectation(hmm, state, x) == (state * hmm->emissionNoPerState + x)/z);
        }
    }

    //Clean up
    hmm_destruct(hmm);
}

void test_hmm_3State(CuTest *testCase) {
    test_hmm(testCase, threeState, nucleotideEmissions, 4);
}

void test_hmm_3StateAsymmetric(CuTest *testCase) {
    test_hmm(testCase, threeStateAsymmetric, nucleotideEmissions, 4);
}

void test_em(CuTest *testCase, StateMachineType stateMachineType, EmissionType emissionType) {
    for (int64_t test = 0; test < 100; test++) {
        //Make a pair of sequences
        char *sX = getRandomSequence(st_randomInt(10, 100));
        char *sY = evolveSequence(sX); //stString_copy(seqX);
        st_logInfo("Sequence X to align: %s END\n", sX);
        st_logInfo("Sequence Y to align: %s END\n", sY);

        //Now do alignment
        PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();

        //Currently starts from random model and iterates.
        double pLikelihood = -INFINITY;
        Hmm *hmm = hmm_constructEmpty(0.0, stateMachineType, emissionType);
        hmm_randomise(hmm);
        StateMachine *sM = hmm_getStateMachine(hmm);
        hmm_destruct(hmm);

        for (int64_t iteration = 0; iteration < 10; iteration++) {
            hmm = hmm_constructEmpty(0.000000000001, stateMachineType, emissionType); //The tiny pseudo count prevents overflow
            getExpectations(sM, hmm, sX, sY, p, 0, 0);
            hmm_normalise(hmm);
            //Log stuff
            for (int64_t from = 0; from < sM->stateNumber; from++) {
                for (int64_t to = 0; to < sM->stateNumber; to++) {
                    st_logInfo("Transition from %" PRIi64 " to %" PRIi64 " has expectation %f\n", from, to,
                            hmm_getTransition(hmm, from, to));
                }
            }
            for (int64_t x = 0; x < hmm->emissionNoPerState; x++) {
            	st_logInfo("Emission x %" PRIi64 " has expectation %f\n", x,
                            hmm_getEmissionsExpectation(hmm, sM->matchState, x));
            }

            st_logInfo("->->-> Got expected likelihood %f for trial %" PRIi64 " and  iteration %" PRIi64 "\n",
                    hmm->likelihood, test, iteration);
            //assert(pLikelihood <= hmm->likelihood * 0.95);
            //CuAssertTrue(testCase, pLikelihood <= hmm->likelihood * 0.95);
            pLikelihood = hmm->likelihood;
            stateMachine_destruct(sM);
            sM = hmm_getStateMachine(hmm);
            hmm_destruct(hmm);
        }

        //Cleanup
        pairwiseAlignmentBandingParameters_destruct(p);
        free(sX);
        free(sY);
        stateMachine_destruct(sM);
    }
}

void test_em_3StateAsymmetric(CuTest *testCase) {
    test_em(testCase, threeStateAsymmetric, nucleotideEmissions);
}

void test_em_3State(CuTest *testCase) {
    test_em(testCase, threeState, nucleotideEmissions);
}

void test_computeForwardProbability(CuTest *testCase) {
	for (int64_t test = 0; test < 1000; test++) {
		// Make a pair of sequences
		char *sX = getRandomSequence(st_randomInt(10, 100));
		char *sY = evolveSequence(sX); //stString_copy(seqX);
		st_logInfo("Sequence X to align: %s END\n", sX);
		st_logInfo("Sequence Y to align: %s END\n", sY);

		// Now do alignment
		PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
		StateMachine *sM = stateMachine3_constructNucleotide(threeState);

		// Forward probability
		stList *anchorPairs = stList_construct();
		bool raggedLeftEnd = st_random() > 0.5;
		bool raggedRightEnd = st_random() > 0.5;
		double logForwardProb = computeForwardProbability(sX, sY, anchorPairs, p, sM, raggedLeftEnd, raggedRightEnd);
		double logForwardProbIdentity = computeForwardProbability(sX, sX, anchorPairs, p, sM, raggedLeftEnd, raggedRightEnd);

		//st_uglyf("Boo:\n\t%s\n\t%s\t%f\t%f\n\n", sX, sY, logForwardProb, logForwardProbIdentity);

		CuAssertTrue(testCase, logForwardProb <= LOG_ONE);
		CuAssertTrue(testCase, logForwardProb > LOG_ZERO);
		CuAssertTrue(testCase, logForwardProb <= logForwardProbIdentity);

		// Cleanup
		stateMachine_destruct(sM);
		free(sX);
		free(sY);
		pairwiseAlignmentBandingParameters_destruct(p);
	}
}

CuSuite* pairwiseAlignmentTestSuite(void) {
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, test_diagonal);
    SUITE_ADD_TEST(suite, test_bands);
    SUITE_ADD_TEST(suite, test_logAdd);
    SUITE_ADD_TEST(suite, test_symbol);
    SUITE_ADD_TEST(suite, test_cell);
    SUITE_ADD_TEST(suite, test_dpDiagonal);
    SUITE_ADD_TEST(suite, test_dpMatrix);
    SUITE_ADD_TEST(suite, test_diagonalDPCalculations);
    SUITE_ADD_TEST(suite, test_getAlignedPairsWithBanding);
    SUITE_ADD_TEST(suite, test_getSplitPoints);
    SUITE_ADD_TEST(suite, test_getAlignedPairs);
    SUITE_ADD_TEST(suite, test_getAlignedPairsWithRaggedEnds);
    SUITE_ADD_TEST(suite, test_getAlignedPairsWithIndels);
    SUITE_ADD_TEST(suite, test_hmm_3State);
    SUITE_ADD_TEST(suite, test_hmm_3StateAsymmetric);
    SUITE_ADD_TEST(suite, test_em_3State);
    SUITE_ADD_TEST(suite, test_em_3StateAsymmetric);
    SUITE_ADD_TEST(suite, test_leftShiftAlignment);
    SUITE_ADD_TEST(suite, test_computeForwardProbability);

    return suite;
}
