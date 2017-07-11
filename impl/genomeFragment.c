/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "stRPHmm.h"

stGenomeFragment *stGenomeFragment_construct(stRPHmm *hmm, stList *path) {
    /*
     * Returns an genome fragment inferred from the hmm and given path through it.
     */
    stGenomeFragment *gF = st_calloc(1, sizeof(stGenomeFragment));

    // Set coordinates
    gF->referenceName = stString_copy(hmm->referenceName);
    gF->refStart = hmm->refStart;
    gF->length = hmm->refLength;
    gF->refCoords = st_calloc(hmm->refLength, sizeof(int64_t));
    gF->insertionsBeforePosition = st_calloc(hmm->refLength, sizeof(int64_t));

    // Allocate genotype arrays
    gF->genotypeString = st_calloc(gF->length, sizeof(uint64_t));
    gF->genotypeProbs = st_calloc(gF->length, sizeof(float));

    // Allocate haplotype arrays
    gF->haplotypeString1 = st_calloc(gF->length, sizeof(uint64_t));
    gF->haplotypeProbs1 = st_calloc(gF->length, sizeof(float));
    gF->haplotypeString2 = st_calloc(gF->length, sizeof(uint64_t));
    gF->haplotypeProbs2 = st_calloc(gF->length, sizeof(float));

    // For each cell in the hmm
    stRPColumn *column = hmm->firstColumn;
    for(int64_t i=0; i<stList_length(path)-1; i++) {
        stRPCell *cell = stList_get(path, i);
        assert(cell != NULL);

        // Calculate the predicted genotype/haplotypes for the given cell
        fillInPredictedGenome(gF, cell, column, hmm->referencePriorProbs, (stRPHmmParameters *)hmm->parameters);

        column = column->nColumn->nColumn;
    }

    // Get predictions for the last column
    assert(column != NULL);
    assert(column->nColumn == NULL);
    fillInPredictedGenome(gF, stList_peek(path), column, hmm->referencePriorProbs, (stRPHmmParameters *)hmm->parameters);

    return gF;
}

void stGenomeFragment_destruct(stGenomeFragment *genomeFragment) {
    // Coordinates
    free(genomeFragment->referenceName);
    free(genomeFragment->refCoords);
    free(genomeFragment->insertionsBeforePosition);

    // Genotypes
    free(genomeFragment->genotypeString);
    free(genomeFragment->genotypeProbs);

    // Haplotypes
    free(genomeFragment->haplotypeString1);
    free(genomeFragment->haplotypeString2);
    free(genomeFragment->haplotypeProbs1);
    free(genomeFragment->haplotypeProbs2);

    free(genomeFragment);
}

void stGenomeFragment_setInsertionCounts(stGenomeFragment *gF) {
    // Add in insertion counts
    gF->insertionsBeforePosition[0] = 0;
    for (int64_t i = 1; i < gF->length; i++) {
        if (gF->refCoords[i] == gF->refCoords[i-1]) {
            gF->insertionsBeforePosition[i] = gF->insertionsBeforePosition[i-1]+1;
        } else {
            gF->insertionsBeforePosition[i] = gF->insertionsBeforePosition[i-1];
        }
    }
    if (gF->insertionsBeforePosition[gF->length - 1] != 0) {
        st_logInfo("* Total insertions in genome fragment: %d \n",
                   gF->insertionsBeforePosition[gF->length -1]);
    }
//    int64_t positionOfInterest = 8098619;
//    if (positionOfInterest >= gF->refCoords[0] && positionOfInterest <= gF->refCoords[gF->length -1]) {
//        st_logInfo("Insertions before position %d: %d\n", positionOfInterest,
//                   gF->insertionsBeforePosition[positionOfInterest - gF->refCoords[0]]);
//        for (int64_t i = 1; i < gF->length; i++) {
//            if (gF->insertionsBeforePosition[i] != gF->insertionsBeforePosition[i-1] && (gF->insertionsBeforePosition[i] < 50 || gF->insertionsBeforePosition[i] > 8480)) {
//                st_logInfo("Insertions at position %d: %d\n", i+gF->refCoords[0],
//                           gF->insertionsBeforePosition[i]);
//            }
//        }
//    }
}
