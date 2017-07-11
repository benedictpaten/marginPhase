/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <htslib/vcf.h>
#include "stRPHmm.h"

/*
 * Functions for reference prior probabilities
 */
static int stHash_intPtrEqualKey(const void *key1, const void *key2) {
    return (*(int64_t *) key1) == (*(int64_t  *)key2);
}

stReferencePriorProbs *stReferencePriorProbs_constructEmptyProfile(char *referenceName,
        int64_t referenceStart, int64_t length) {
    /*
     * Creates an empty object, with all the profile probabilities set to 0.
     */
    stReferencePriorProbs *referencePriorProbs = st_calloc(1, sizeof(stReferencePriorProbs));
    referencePriorProbs->referenceName = stString_copy(referenceName);
    referencePriorProbs->refStart = referenceStart;
    referencePriorProbs->length = length;
    referencePriorProbs->profileProbs = st_calloc(length*ALPHABET_SIZE, sizeof(uint16_t));
    referencePriorProbs->profileSequence = st_calloc(length, sizeof(uint8_t));
    referencePriorProbs->baseCounts = st_calloc(length*ALPHABET_SIZE, sizeof(double));
    referencePriorProbs->insertionCounts = st_calloc(length, sizeof(int64_t));
    referencePriorProbs->gapSizes = st_calloc(length, sizeof(int64_t));
    referencePriorProbs->refCoords = st_calloc(length, sizeof(int64_t));
    referencePriorProbs->refCoordMap = stHash_construct3(stHash_stringKey, stHash_intPtrEqualKey, NULL, NULL);
    return referencePriorProbs;
}

void stReferencePriorProbs_destruct(stReferencePriorProbs *referencePriorProbs) {
    /*
     * Destructor for reference prior
     */
    free(referencePriorProbs->profileProbs);
    free(referencePriorProbs->referenceName);
    free(referencePriorProbs->profileSequence);
    free(referencePriorProbs->insertionCounts);
    free(referencePriorProbs->gapSizes);
    free(referencePriorProbs->refCoords);
    stHash_destruct(referencePriorProbs->refCoordMap);
    free(referencePriorProbs);
}

static stReferencePriorProbs *getNext(stList *profileSequences) {
    /*
     * Construct an empty stReferencePriorProbs for the next interval of a reference sequence.
     */
    assert(stList_length(profileSequences) > 0);
    stProfileSeq *pSeq = stList_pop(profileSequences);
    char *refName = pSeq->referenceName;
    int64_t refStart = pSeq->refStart;
    int64_t refEnd = pSeq->refStart + pSeq->length;
    while(stList_length(profileSequences) > 0) {
        stProfileSeq *pSeq = stList_peek(profileSequences);

        if(!stString_eq(pSeq->referenceName, refName)) {
            break;
        }
        stList_pop(profileSequences);

        assert(pSeq->refStart <= refStart);
        refStart = pSeq->refStart;
        refEnd = pSeq->refStart + pSeq->length > refEnd ? pSeq->refStart + pSeq->length : refEnd;
    }

    return stReferencePriorProbs_constructEmptyProfile(refName, refStart, refEnd-refStart);
}

void setRefCoordinates(stReferencePriorProbs *rProbs, stList *profileSequences, bool includeInsertions) {
    /*
     * Sets the reference coordinates in the set of profile sequences at each position.
     */
    // Walk through rProbs, setting coordinates as you go
    int64_t rProbsIndex = 0;
    for(int64_t i=0; i<stList_length(profileSequences); i++) {
        stProfileSeq *pSeq = stList_get(profileSequences, i);
        for (int64_t j = 0; j < pSeq->length; j++) {
            if (rProbsIndex == j + pSeq->refStart - rProbs->refStart) {
                rProbs->refCoords[rProbsIndex] = pSeq->refCoords[j];
                rProbsIndex++;
            }
//            int64_t k = includeInsertions?
//                        findIndexForRefCoord(rProbs->refCoords, pSeq->refCoords, j, rProbs->length)
//                                : j + pSeq->refStart - rProbs->refStart;
//            // TODO: Will update the same index many times...
//            rProbs->refCoords[k] = pSeq->refCoords[j];
        }
    }
    for (int64_t i = 0; i < rProbs->length; i++) {
        if (rProbs->refCoords[i] == 0) {
            st_logInfo("rProbs->refCoords[ %d ] = 0 \n", rProbs->refStart + i);
        }
    }
}

void setReadBaseCounts(stReferencePriorProbs *rProbs, stList *profileSequences, bool includeInsertions) {
    /*
     * Set the base counts observed in the set of profile sequences at each position.
     */
    for(int64_t i=0; i<stList_length(profileSequences); i++) {
        stProfileSeq *pSeq = stList_get(profileSequences, i);
        for (int64_t j = 0; j < pSeq->length; j++) {
            int64_t k = includeInsertions ?
                        findIndexForRefCoord(pSeq->refCoords, rProbs->refCoords, j, pSeq->length)
                                          : j + pSeq->refStart - rProbs->refStart;
            if(k >= 0 && k < rProbs->length) {
                for (int64_t l = 0; l < ALPHABET_SIZE; l++) {
                    rProbs->baseCounts[k*ALPHABET_SIZE + l] +=
                            getProb(&(pSeq->profileProbs[j*ALPHABET_SIZE]), l);
                }
            }
//            rProbs->refCoords[k] = pSeq->refCoords[j];
        }
    }
}

void setInsertions(stReferencePriorProbs *rProbs, stList *profileSequences) {
    /*
     * Set the counts of insertions observed at each position in the set of profile sequences.
     */
    int64_t numDiffInsertionSizes = 0;
    int64_t numSameInsertionSizes = 0;
    for (int64_t i = 0; i < stList_length(profileSequences); i++) {
        stProfileSeq *pSeq = stList_get(profileSequences, i);
        for (int64_t j = 0; j < pSeq->length; j++) {
            int64_t k = j + pSeq->refStart - rProbs->refStart;
            if(k >= 0 && k < rProbs->length) {
                rProbs->insertionCounts[k] += (pSeq->insertions[j] > 1 ? 1 : pSeq->insertions[j]);

                if ((rProbs->gapSizes[k] > 0 && pSeq-> insertions[j] > 0) &&
                        (rProbs->gapSizes[k] != pSeq->insertions[j]))
                {
                    numDiffInsertionSizes++;
                } else {
                    numSameInsertionSizes++;
                }
                // Update gap size
                // If there's already a gap, choose the smaller gap size
                // TODO: choose most common size instead?
                if (pSeq->insertions[j] > 0) {
                    if (rProbs->gapSizes[k] == 0 ||
                            (rProbs->gapSizes[k] > 0 && pSeq->insertions[j] < rProbs->gapSizes[k])) {
                        rProbs->gapSizes[k] = pSeq->insertions[j];
                    }
                }
            }
        }
    }
//    st_logInfo("Number of times gap size differed: %d \n", numDiffInsertionSizes);
//    st_logInfo("Number of times gap size same: %d \n", numSameInsertionSizes);
}


stHash *createEmptyReferencePriorProbabilities(stList *profileSequences, bool includeInsertions) {
    /*
     * Create a set of stReferencePriorProbs that cover the reference intervals included in the profile sequences.
     * Each has a flat probability over reference positions.
     * The return value is encoded as a map from the reference sequence name (as a string)
     * to the stReferencePriorProbs.
     */

    // Make map from reference sequence names to reference priors
    stHash *referenceNamesToReferencePriors = stHash_construct3(stHash_stringKey, stHash_stringEqualKey,
            NULL, (void (*)(void *))stReferencePriorProbs_destruct);

    // Copy and sort profile sequences in order
    // TODO need destructors?
    profileSequences = stList_copy(profileSequences, (void (*)(void *))stProfileSeq_destruct);
    stList *profileSequencesCopy = stList_copy(profileSequences, (void (*)(void *))stProfileSeq_destruct);
    stList_sort(profileSequences, stRPProfileSeq_cmpFn);

    uint16_t flatValue = scaleToLogIntegerSubMatrix(log(1.0/ALPHABET_SIZE));
    while(stList_length(profileSequences) > 0) {
        // Get next reference prior prob interval
        stReferencePriorProbs *rProbs = getNext(profileSequences);

        // Fill in probability profile
        for(int64_t i=0; i<rProbs->length; i++) {
            for(int64_t j=0; j<ALPHABET_SIZE; j++) {
                rProbs->profileProbs[i*ALPHABET_SIZE + j] = flatValue;
            }
        }
        setRefCoordinates(rProbs, profileSequencesCopy, includeInsertions);
        for (int64_t i = 0; i < rProbs->length; i++) {
            if (rProbs->refCoords[i] == 0 && i < 14400) {
                st_logInfo("1 refCoords[ %d ] = 0\n", i+rProbs->refStart);
            }
        }
        setReadBaseCounts(rProbs, profileSequencesCopy, includeInsertions);
        for (int64_t i = 0; i < rProbs->length; i++) {
            if (rProbs->refCoords[i] == 0 && i < 14400) {
                st_logInfo(" 2 refCoords[ %d ] = 0\n", i+rProbs->refStart);
            }
        }
        setInsertions(rProbs, profileSequencesCopy);

        assert(stHash_search(referenceNamesToReferencePriors, rProbs->referenceName) == NULL);
        stHash_insert(referenceNamesToReferencePriors, rProbs->referenceName, rProbs);

        for (int64_t i = 0; i < rProbs->length; i++) {
            if (rProbs->refCoords[i] == 0 && i < 14400) {
                st_logInfo("  3 refCoords[ %d ] = 0\n", i+rProbs->refStart);
            }
        }
    }



    // Cleanup
    stList_destruct(profileSequences);

    return referenceNamesToReferencePriors;
}

stHash *createReferencePriorProbabilities(char *referenceFastaFile, stList *profileSequences,
        stBaseMapper *baseMapper, stRPHmmParameters *params, bool includeInsertions) {
    /*
     * Create a set of stReferencePriorProbs that cover the reference intervals included in the profile sequences.
     * The return value is encoded as a map from the reference sequence name (as a string)
     * to the stReferencePriorProbs.
     */

    // Make map from reference sequence names to reference priors
    stHash *referenceNamesToReferencePriors = createEmptyReferencePriorProbabilities(profileSequences, includeInsertions);

    // Load reference fasta index
    faidx_t *fai = fai_load(referenceFastaFile);
    if ( !fai ) {
        st_errAbort("Could not load fai index of %s.  Maybe you should run 'samtools faidx %s'\n",
                       referenceFastaFile, referenceFastaFile);
    }

    stHashIterator *hashIt = stHash_getIterator(referenceNamesToReferencePriors);
    char *referenceName;
    while((referenceName = stHash_getNext(hashIt)) != NULL) {
        stReferencePriorProbs *rProbs = stHash_search(referenceNamesToReferencePriors, referenceName);

        // Now get the corresponding reference sequence
        int seqLen;
        char *referenceSeq = fai_fetch(fai, rProbs->referenceName, &seqLen);
        if ( seqLen < 0 ) {
            st_errAbort("Failed to fetch reference sequence %s\n", rProbs->referenceName);
        }

        // Build probability profile
        assert(seqLen >= rProbs->length + rProbs->refStart);
        for(int64_t i=0; i<rProbs->length; i++) {
//            uint8_t refChar = stBaseMapper_getValueForChar(baseMapper, referenceSeq[i+rProbs->refStart-1]);
            uint8_t refChar = stBaseMapper_getValueForChar(baseMapper,
                                                           referenceSeq[rProbs->refCoords[i] - 1]);
            assert(refChar >= 0 && refChar < ALPHABET_SIZE);
            rProbs->profileSequence[i] = refChar;
            for(int64_t j=0; j<ALPHABET_SIZE; j++) {
                rProbs->profileProbs[i*ALPHABET_SIZE + j] = *getSubstitutionProb(params->hetSubModel, refChar, j);
            }
        }
    }

    // Cleanup
    fai_destroy(fai);
    stHash_destructIterator(hashIt);

    return referenceNamesToReferencePriors;
}

int64_t findIndexForRefCoord(int64_t *coords, int64_t *refCoords, int64_t index, int64_t coordsLength) {
    /*
     * Given the index in a profile sequence refSeq, find the index that corresponds
     * to the same location in another progile sequence, seq.
     */
    int64_t refCoord = refCoords[index];

    // If in an insertion, find which index in the gap it is
    int64_t tempIndex = index;
    int64_t offsetL = 0;
    while (refCoords[tempIndex] == refCoords[tempIndex - 1]) {
        offsetL++;
        tempIndex--;
    }
    tempIndex = index;
    int64_t offsetR = 0;
    while (refCoords[tempIndex] == refCoords[tempIndex + 1]) {
        offsetR++;
        tempIndex++;
    }

    int64_t coordIndex = refCoords[index] - coords[0];
    if (coordIndex < 0 || coordIndex >= coordsLength) return coordIndex;

    int64_t currentCoord =  coords[coordIndex];

    while (currentCoord < refCoord) {
        coordIndex++;
        if (coordIndex < 0 || coordIndex >= coordsLength) return coordIndex;
        currentCoord = coords[coordIndex];
        if (currentCoord == refCoord) coordIndex += offsetL;
    }

    while (currentCoord > refCoord) {
        coordIndex--;
        if (coordIndex < 0 || coordIndex >= coordsLength) return coordIndex;
        currentCoord = coords[coordIndex];
        if (currentCoord == refCoord) coordIndex -= offsetR;
    }


//    if (refCoord != coords[coordIndex]) {
// TODO delete
//            st_logInfo("Coords don't match! Ref[%d] = %d, Seq[%d] = %d\n", index+refSeq->refStart, refCoord, coordIndex+seq->refStart, seq->refCoords[coordIndex]);
//    }

    return coordIndex;
}