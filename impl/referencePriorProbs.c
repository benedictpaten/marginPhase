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
    referencePriorProbs->refCoords = st_calloc(length, sizeof(int64_t));
    referencePriorProbs->refCoordMap = stHash_construct3(stHash_stringKey, stHash_intPtrEqualKey, NULL, NULL);

//    referencePriorProbs->refIndexes = st_calloc(length, sizeof(stRefIndex));
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
    free(referencePriorProbs->refCoords);
    stHash_destruct(referencePriorProbs->refCoordMap);

//    free(referencePriorProbs->refIndexes);
    free(referencePriorProbs);
}

static stReferencePriorProbs *getNext(stList *profileSequences, int64_t numInsertions) {
    /*
     * Construct an empty stReferencePriorProbs for the next interval of a reference sequence.
     */
    assert(stList_length(profileSequences) > 0);
    stProfileSeq *pSeq = stList_pop(profileSequences);
    char *refName = pSeq->referenceName;
    int64_t refStart = pSeq->refStart;
    int64_t refEnd = pSeq->refEnd;

    while(stList_length(profileSequences) > 0) {
        stProfileSeq *pSeq = stList_peek(profileSequences);

        if(!stString_eq(pSeq->referenceName, refName)) {
            break;
        }
        stList_pop(profileSequences);

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
    // Walk through rProbs, setting coordinates as you go
    // TODO: figure out how to do this without calloc...
    int64_t *indexes = st_calloc(rProbs->length, sizeof(int64_t));

    // First position
    rProbs->refCoords[0] = rProbs->refStart;
    indexes[0] = 0;
    if (stHash_search(rProbs->refCoordMap, &rProbs->refStart) == NULL) {
        stHash_insert(rProbs->refCoordMap, &rProbs->refStart, &indexes[0]);
    }
    int64_t rProbsIndex = 1;
    int64_t nextSeqStartingIndex = 1;
    stProfileSeq *pSeq = stList_get(profileSequences, 0);
    for(int64_t i=1; i<stList_length(profileSequences); i++) {
        if (nextSeqStartingIndex <= 0 || nextSeqStartingIndex >= pSeq->length) {
            stProfileSeq *pSeq2 = stList_get(profileSequences, i);
            nextSeqStartingIndex = findCorrespondingRefCoordIndex(pSeq->length-1, pSeq->refCoords, pSeq2->refCoordMap) + 1;
            if (nextSeqStartingIndex > 0 && nextSeqStartingIndex < pSeq2->length) {
                pSeq = pSeq2;
            }
        } else {
            // Only look at a profile sequence until you reach the start of the next
            for (int64_t j = nextSeqStartingIndex; j < pSeq->length; j++) {
                rProbs->refCoords[rProbsIndex] = pSeq->refCoords[j];
                indexes[rProbsIndex] = rProbsIndex;
                // Only insert if not already in the table
                if (stHash_search(rProbs->refCoordMap, &pSeq->refCoords[j]) == NULL) {
                    stHash_insert(rProbs->refCoordMap, &pSeq->refCoords[j], &indexes[rProbsIndex]);
                }
                rProbsIndex++;
            }
            stProfileSeq *pSeq2 = stList_get(profileSequences, i);
            nextSeqStartingIndex = findCorrespondingRefCoordIndex(pSeq->length-1, pSeq->refCoords, pSeq2->refCoordMap) + 1;

            if (nextSeqStartingIndex > 0 && nextSeqStartingIndex < pSeq2->length) {
                pSeq = pSeq2;
            }
        }
    }
    rProbs->refEnd = rProbs->refCoords[rProbs->length-1];
}

void setReadBaseCounts(stReferencePriorProbs *rProbs, stList *profileSequences, int64_t numInsertions) {
    /*
     * Set the base counts observed in the set of profile sequences at each position.
     */
    for(int64_t i=0; i<stList_length(profileSequences); i++) {
        stProfileSeq *pSeq = stList_get(profileSequences, i);
        for (int64_t j = 0; j < pSeq->length; j++) {
            int64_t k = numInsertions > 0 ?
                        findCorrespondingRefCoordIndex(j, pSeq->refCoords, rProbs->refCoordMap)
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

void setInsertions(stReferencePriorProbs *rProbs, stList *profileSequences, int64_t insertionCountThreshold) {
    /*
     * Set the counts of insertions observed at each position in the set of profile sequences.
     */

    int64_t numDiffInsertionSizes = 0;
    int64_t numSameInsertionSizes = 0;
    for (int64_t i = 0; i < stList_length(profileSequences); i++) {
        stProfileSeq *pSeq = stList_get(profileSequences, i);
        for (int64_t j = 0; j < pSeq->length; j++) {
            int64_t k = findCorrespondingRefCoordIndex(j, pSeq->refCoords, rProbs->refCoordMap);
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
//    st_logInfo("Number of times gap size differed: %d \n", numDiffInsertionSizes);
//    st_logInfo("Number of times gap size same: %d \n", numSameInsertionSizes);
}


stHash *createEmptyReferencePriorProbabilities(stList *profileSequences, int64_t numInsertions, stRPHmmParameters *params) {
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
        stReferencePriorProbs *rProbs = getNext(profileSequences, numInsertions);

        // Fill in probability profile
        for(int64_t i=0; i<rProbs->length; i++) {
            for(int64_t j=0; j<ALPHABET_SIZE; j++) {
                rProbs->profileProbs[i*ALPHABET_SIZE + j] = flatValue;
            }
        }
        setRefCoordinates(rProbs, profileSequencesCopy);
        setReadBaseCounts(rProbs, profileSequencesCopy, numInsertions);
        setInsertions(rProbs, profileSequencesCopy, params->insertionCountThreshold);

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
            createEmptyReferencePriorProbabilities(profileSequences, numInsertions, params);

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
            uint8_t refChar =
                    stBaseMapper_getValueForChar(baseMapper, referenceSeq[rProbs->refCoords[i] - 1]);
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
