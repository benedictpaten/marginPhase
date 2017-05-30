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

