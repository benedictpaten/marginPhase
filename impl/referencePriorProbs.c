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

int stHash_intPtrEqualKey(const void *key1, const void *key2) {
    /*
     * Function to test if two int hash keys are equal.
     */
    return (*(int64_t *) key1) == (*(int64_t  *)key2);
}

int64_t getRefCoordIndex(int64_t refCoord, stReferencePriorProbs *rProbs, bool last) {
    /*
     * Given a reference coordinate, find which index in rProbs it corresponds to.
     * Return -1 if not in rProbs
     */
    int64_t *idxPtr = stHash_search(rProbs->refCoordMap, &refCoord);
    int64_t index;
    if (idxPtr != NULL) {
        index = *idxPtr;
        if (last && index < rProbs->length - 1) {
            while(rProbs->refIndexes[index]->refCoord == rProbs->refIndexes[index+1]->refCoord) {
                index++;
            }
            assert(rProbs->refIndexes[index]->refCoord != rProbs->refIndexes[index+1]->refCoord);
        } else {
            if (index > 0) assert(rProbs->refIndexes[index]->refCoord != rProbs->refIndexes[index-1]->refCoord);
        }
        return index;
    } else {
        return -1;
    }
}

stReferencePriorProbs *stReferencePriorProbs_constructEmptyProfile(char *referenceName,
        int64_t referenceStart, int64_t refLength, int64_t numInsertions) {
    /*
     * Creates an empty object, with all the profile probabilities set to 0.
     */
    stReferencePriorProbs *referencePriorProbs = st_calloc(1, sizeof(stReferencePriorProbs));
    referencePriorProbs->referenceName = stString_copy(referenceName);
    referencePriorProbs->refStart = referenceStart;
    int64_t length = refLength + numInsertions;
    referencePriorProbs->length = length;
    referencePriorProbs->profileProbs = st_calloc(length*ALPHABET_SIZE, sizeof(uint16_t));
    referencePriorProbs->profileSequence = st_calloc(length, sizeof(uint8_t));
    referencePriorProbs->baseCounts = st_calloc(length*ALPHABET_SIZE, sizeof(double));
    referencePriorProbs->insertionCounts = st_calloc(length, sizeof(int64_t));
    referencePriorProbs->gapSizes = st_calloc(length, sizeof(int64_t));
    referencePriorProbs->insertionsBeforePosition = st_calloc(length, sizeof(int64_t));
    referencePriorProbs->refCoordMap = stHash_construct3(stHash_stringKey, stHash_intPtrEqualKey, NULL, NULL);
    referencePriorProbs->refIndexes = st_calloc(length, sizeof(stRefIndex));
    referencePriorProbs->referencePositionsIncluded = st_calloc(length, sizeof(bool));
    for(int64_t i=0; i<length; i++) {
        referencePriorProbs->referencePositionsIncluded[i] = true;
    }

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
    free(referencePriorProbs->insertionsBeforePosition);
    stHash_destruct(referencePriorProbs->refCoordMap);
    free(referencePriorProbs->refIndexes);
    free(referencePriorProbs->baseCounts);
    free(referencePriorProbs->referencePositionsIncluded);
    free(referencePriorProbs);
}

stRefIndex *stRefIndex_construct(int64_t refCoord, int64_t index) {
    /*
     * Constructor for reference index
     */
    stRefIndex *refIndex = st_calloc(1, sizeof(stRefIndex));
    refIndex->refCoord = refCoord;
    refIndex->index = index;
    return refIndex;
}


static stReferencePriorProbs *getNext(stList *profileSequences, stList *profileSequencesOnReference, int64_t numInsertions) {
    /*
     * Construct an empty stReferencePriorProbs for the next interval of a reference sequence.
     */
    assert(stList_length(profileSequences) > 0);
    stProfileSeq *pSeq = stList_pop(profileSequences);
    char *refName = pSeq->referenceName;
    int64_t refStart = pSeq->refStart;
    int64_t refEnd = pSeq->refEnd;

    // Add to list of profile sequences contained within the reference interval
    // being considered
    stList_append(profileSequencesOnReference, pSeq);

    while(stList_length(profileSequences) > 0) {
        stProfileSeq *pSeq = stList_peek(profileSequences);

        if(!stString_eq(pSeq->referenceName, refName)) {
            break;
        }
        stList_pop(profileSequences);

        stList_append(profileSequencesOnReference, pSeq);

        assert(pSeq->refStart <= refStart);
        refStart = pSeq->refStart;
        refEnd = pSeq->refEnd > refEnd ? pSeq->refEnd : refEnd;
    }

    stReferencePriorProbs *rProbs =  stReferencePriorProbs_constructEmptyProfile(refName, refStart, refEnd-refStart+1, numInsertions);
    return rProbs;
}

void setRefCoordinates(stReferencePriorProbs *rProbs, stList *profileSequences) {
    /*
     * Sets the reference coordinates in the set of profile sequences at each position.
     */
    // First position
    rProbs->refIndexes[0] = stRefIndex_construct(rProbs->refStart, 0);
    if (stHash_search(rProbs->refCoordMap, &rProbs->refIndexes[0]->refCoord) == NULL) {
        stHash_insert(rProbs->refCoordMap, &rProbs->refIndexes[0]->refCoord, &rProbs->refIndexes[0]->index);
    }
    // Walk through rProbs, setting coordinates as you go
    int64_t rProbsIndex = 1;
    int64_t nextSeqStartingIndex = 1;
    stProfileSeq *pSeq = stList_get(profileSequences, 0);
    for(int64_t i=1; i<stList_length(profileSequences); i++) {
        if (nextSeqStartingIndex < 0 || nextSeqStartingIndex >= pSeq->length) {
            stProfileSeq *pSeq2 = stList_get(profileSequences, i);
            nextSeqStartingIndex = findCorrespondingRefCoordIndex(pSeq->length-1, pSeq->refIndexes, pSeq2->refCoordMap);
            if (nextSeqStartingIndex >= 0 && nextSeqStartingIndex < pSeq2->length) {
                pSeq = pSeq2;
                nextSeqStartingIndex++;
            }
        } else {
            for (int64_t j = nextSeqStartingIndex; j < pSeq->length; j++) {
                rProbs->refIndexes[rProbsIndex] =
                        stRefIndex_construct(pSeq->refIndexes[j]->refCoord, rProbsIndex);

                // Only insert if not already in the table (to make sure index is first on an insertion)
                if (stHash_search(rProbs->refCoordMap, &rProbs->refIndexes[rProbsIndex]->refCoord) == NULL) {
                    stHash_insert(rProbs->refCoordMap, &rProbs->refIndexes[rProbsIndex]->refCoord, &rProbs->refIndexes[rProbsIndex]->index);
                }
                rProbsIndex++;
            }
            // Find the next index to look at in the next sequence, so we don't repeat coordinates
            stProfileSeq *pSeq2 = stList_get(profileSequences, i);
            nextSeqStartingIndex = findCorrespondingRefCoordIndex(pSeq->length-1, pSeq->refIndexes, pSeq2->refCoordMap);

            // If the next index is actually in the next read, use it
            if (nextSeqStartingIndex >= 0 && nextSeqStartingIndex < pSeq2->length) {
                pSeq = pSeq2;
                nextSeqStartingIndex++;
            }
        }
    }


    rProbs->refEnd = rProbs->refIndexes[rProbs->length-1]->refCoord;
}

void setReadBaseCounts(stReferencePriorProbs *rProbs, stList *profileSequences, int64_t numInsertions) {
    /*
     * Set the base counts observed in the set of profile sequences at each reference position.
     */
    for(int64_t i=0; i<stList_length(profileSequences); i++) {
        stProfileSeq *pSeq = stList_get(profileSequences, i);
        assert(stString_eq(rProbs->referenceName, pSeq->referenceName));
        for (int64_t j = 0; j < pSeq->length; j++) {
            int64_t k = numInsertions > 0 ?
                        findCorrespondingRefCoordIndex(j, pSeq->refIndexes, rProbs->refCoordMap)
                                          : j + pSeq->refStart - rProbs->refStart;
            if(k >= 0 && k < rProbs->length) {
                for (int64_t l = 0; l < ALPHABET_SIZE; l++) {
                    rProbs->baseCounts[k*ALPHABET_SIZE + l] +=
                            getProb(&(pSeq->profileProbs[j*ALPHABET_SIZE]), l);
                }
            }
        }
    }
}

void setInsertions(stReferencePriorProbs *rProbs, stList *profileSequences) {
    /*
     * Set the counts of insertions observed at each position in the set of profile sequences.
     */

    for (int64_t i = 0; i < stList_length(profileSequences); i++) {
        stProfileSeq *pSeq = stList_get(profileSequences, i);
        for (int64_t j = 0; j < pSeq->length; j++) {
            int64_t k = findCorrespondingRefCoordIndex(j, pSeq->refIndexes, rProbs->refCoordMap);
            if(k >= 0 && k < rProbs->length) {
                rProbs->insertionCounts[k] += (pSeq->insertions[j] > 1 ? 1 : pSeq->insertions[j]);

                // Update gap size
                // If there's already a gap, choose the smaller gap size
                // TODO: choose most common size instead? or floor of average?
                if (pSeq->insertions[j] > 0) {
                    if (rProbs->gapSizes[k] == 0 ||
                            (rProbs->gapSizes[k] > 0 && pSeq->insertions[j] < rProbs->gapSizes[k])) {
                        rProbs->gapSizes[k] = pSeq->insertions[j];
                    }
                }
            }
        }
    }
    // Walk back through and set field insertionsBeforePosition in rProbs
    for (int64_t i = 1; i < rProbs->length; i++) {
        if (rProbs->insertionCounts[i] > 1) {
            rProbs->insertionsBeforePosition[i] = rProbs->insertionsBeforePosition[i-1] + 1;
        } else {
            rProbs->insertionsBeforePosition[i] = rProbs->insertionsBeforePosition[i-1];
        }
    }
}


stHash *createEmptyReferencePriorProbabilities(stList *profileSequences, int64_t numInsertions) {
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
    stList_sort(profileSequences, stRPProfileSeq_cmpFn);

    uint16_t flatValue = scaleToLogIntegerSubMatrix(log(1.0/ALPHABET_SIZE));
    while(stList_length(profileSequences) > 0) {
        // Get next reference prior prob interval
        stList *l = stList_construct();
        stReferencePriorProbs *rProbs = getNext(profileSequences, l, numInsertions);
        stList_sort(l, stRPProfileSeq_cmpFn);

        // Fill in probability profile
        for(int64_t i=0; i<rProbs->length; i++) {
            for(int64_t j=0; j<ALPHABET_SIZE; j++) {
                rProbs->profileProbs[i*ALPHABET_SIZE + j] = flatValue;
            }
        }
        // Set the coordinates and counts of each base observed at each reference position
        setRefCoordinates(rProbs, l);
        setReadBaseCounts(rProbs, l, numInsertions);
        setInsertions(rProbs, l);

        stList_destruct(l);

        assert(stHash_search(referenceNamesToReferencePriors, rProbs->referenceName) == NULL);
        stHash_insert(referenceNamesToReferencePriors, rProbs->referenceName, rProbs);
    }
    // Cleanup
    stList_destruct(profileSequences);

    return referenceNamesToReferencePriors;
}

stHash *createReferencePriorProbabilities(char *referenceFastaFile, stList *profileSequences,
        stBaseMapper *baseMapper, stRPHmmParameters *params, int64_t numInsertions) {
    /*
     * Create a set of stReferencePriorProbs that cover the reference intervals included in the profile sequences.
     * The return value is encoded as a map from the reference sequence name (as a string)
     * to the stReferencePriorProbs.
     */

    // Make map from reference sequence names to reference priors
    stHash *referenceNamesToReferencePriors =
            createEmptyReferencePriorProbabilities(profileSequences, numInsertions);

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
        assert(seqLen >= rProbs->refEnd + 1);
        for(int64_t i=0; i<rProbs->length; i++) {
            int64_t refSeqIndex = rProbs->refIndexes[i]->refCoord;
            uint8_t refChar = stBaseMapper_getValueForChar(baseMapper, referenceSeq[refSeqIndex- 1]);
            assert(refChar >= 0 && refChar < ALPHABET_SIZE);
            rProbs->profileSequence[i] = refChar;
            for(int64_t j=0; j<ALPHABET_SIZE; j++) {
                rProbs->profileProbs[i*ALPHABET_SIZE + j] =
                        *getSubstitutionProb(params->hetSubModel, refChar, j);
            }
        }
    }

    // Cleanup
    fai_destroy(fai);
    stHash_destructIterator(hashIt);

    return referenceNamesToReferencePriors;
}

int64_t numTotalInsertionColumns(stReferencePriorProbs *rProbs, int64_t threshold) {
    /*
     * Returns the total number of columns that will be inserted into the model.
     */
    int64_t numColumns = 0;
    for (int64_t i = 0; i < rProbs->length; i++) {
        if (rProbs->insertionCounts[i] >= threshold) {
            numColumns += rProbs->gapSizes[i];
        }
    }
    return numColumns;
}

int64_t stReferencePriorProbs_setReferencePositionFilter(stReferencePriorProbs *rProbs, stRPHmmParameters *params) {
    /*
     * Sets the rProbs->referencePositionsIncluded array to exclude positions likely to be homozygous.
     *
     * In a column, let x be the number of occurrences of the most frequent non-reference base.
     * If x is less than params.minNonReferenceBaseFilter then the column is filtered out.
     */
    int64_t filteredPositions = 0;
    for(int64_t i=0; i<rProbs->length; i++) {
        double *baseCounts = &rProbs->baseCounts[i*ALPHABET_SIZE];

        // Get total expected coverage at base and most frequent base in column
        double totalCount = baseCounts[0];
        int64_t mostFrequentBase = 0;
        double mostFrequentBaseCount = totalCount;
        for(int64_t j=1; j<ALPHABET_SIZE; j++) {
            totalCount += baseCounts[j];
            if(baseCounts[j] > mostFrequentBaseCount) {
                mostFrequentBase = j;
                mostFrequentBaseCount = baseCounts[j];
            }
        }

        // Maximum number of instances of a non-reference base
        double secondMostFrequentBaseCount = 0.0;
        for(int64_t j=0; j<ALPHABET_SIZE; j++) {
            if(j != mostFrequentBase && baseCounts[j] > secondMostFrequentBaseCount) {
                secondMostFrequentBaseCount = baseCounts[j];
            }
        }

        if(secondMostFrequentBaseCount < params->minSecondMostFrequentBaseFilter) {
            filteredPositions++;
            // Mask the position
            assert(rProbs->referencePositionsIncluded[i] == true);
            rProbs->referencePositionsIncluded[i] = false;
        }
        else if(st_getLogLevel() == debug) {
            stList *l = stList_construct3(0, free);
            for(int64_t k=0; k<ALPHABET_SIZE; k++) {
                stList_append(l, stString_print("%f", baseCounts[k]));
            }
            char *baseCountString = stString_join2(" ", l);
//            st_logDebug("Including position: %s %" PRIi64 " with coverage depth: %f and base counts: %s\n",
//                    rProbs->referenceName, rProbs->refIndexes[i]->refCoord, totalCount, baseCountString);
            stList_destruct(l);
            free(baseCountString);
        }
    }

    return filteredPositions;
}

int64_t filterHomozygousReferencePositions(stHash *referenceNamesToReferencePriors, stRPHmmParameters *params, int64_t *totalPositions) {
    /*
     * Sets a column filter to ignore columns likely to be homozygous reference.
     *
     * Returns total number of positions filtered out, sets totalPositions to the total number of
     * reference positions.
     */
    int64_t filteredPositions = 0;
    *totalPositions = 0;
    stHashIterator *hashIt = stHash_getIterator(referenceNamesToReferencePriors);
    char *referenceName;
    while((referenceName = stHash_getNext(hashIt)) != NULL) {
        stReferencePriorProbs *rProbs = stHash_search(referenceNamesToReferencePriors, referenceName);
        filteredPositions += stReferencePriorProbs_setReferencePositionFilter(rProbs, params);
        *totalPositions += rProbs->length;
    }

    // Cleanup
    stHash_destructIterator(hashIt);

    return filteredPositions;
}

int64_t insertionsInRange(stReferencePriorProbs *rProbs, int64_t leftIndex, int64_t rightIndex, int64_t threshold) {
    /*
     * Calculates the number of insertions in a given range on indexes on the reference,
     * up to but not including the end index.
     */
    int64_t numInsertions = 0;
    for (int64_t i = leftIndex; i < rightIndex; i++) {
        if (rProbs->insertionCounts[i] >= threshold) {
            numInsertions += rProbs->gapSizes[i];
        }
    }
    return numInsertions;
}
