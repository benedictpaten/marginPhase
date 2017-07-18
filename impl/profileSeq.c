/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "stRPHmm.h"

/*
 * Functions for profile sequence
 */

stProfileSeq *stProfileSeq_constructEmptyProfile(char *referenceName, char *readId,
                                                 int64_t referenceStart, int64_t length) {
    /*
     * Creates an empty profile sequence, with all the profile probabilities set to 0.
     */
    stProfileSeq *seq = st_malloc(sizeof(stProfileSeq));
    seq->referenceName = stString_copy(referenceName);
    seq->readId = stString_copy(readId);
    seq->refStart = referenceStart;
    seq->refEnd = referenceStart + length - 1;
    seq->length = length;
    seq->profileProbs = st_calloc(length*ALPHABET_SIZE, sizeof(uint8_t));
    seq->refCoords = st_calloc(length, sizeof(int64_t));
    seq->insertions = st_calloc(length, sizeof(int64_t));
    seq->numInsertions = 0;
    seq->insertionSeqs = st_calloc(length, sizeof(int64_t *));
    seq->refCoordMap = stHash_construct3(stHash_stringKey, stHash_intPtrEqualKey, NULL, NULL);
    return seq;
}

void stProfileSeq_destruct(stProfileSeq *seq) {
    /*
     * Cleans up memory for profile sequence.
     */
    for (int64_t i = 0; i < seq->length; i++) {
        free(seq->insertionSeqs[i]);
    }
    free(seq->insertionSeqs);
//    stHash_destruct(seq->refCoordMap);
//    free(seq->refCoords);
    free(seq->insertions);
    free(seq->profileProbs);
    free(seq->readId);
    free(seq->referenceName);
    free(seq);
}

float getProb(uint8_t *p, int64_t characterIndex) {
    /*
     * Gets probability of a given character as a float.
     */
    return ((float)p[characterIndex])/255;
}

void stProfileSeq_print(stProfileSeq *seq, FILE *fileHandle, bool includeProbs) {
    /*
     * Prints a debug representation of a profile sequence.
     */
    char profileString[seq->length+1];
    profileString[seq->length] = '\0';
    for(int64_t i=0; i<seq->length; i++) {
        uint8_t *p = &seq->profileProbs[i * ALPHABET_SIZE];
        float maxProb = getProb(p, 0);
        int64_t maxChar = 0;
        for(int64_t j=1; j<ALPHABET_SIZE; j++) {
            float prob = getProb(p, j);
            if(prob > maxProb) {
                maxProb = prob;
                maxChar = j;
            }
        }
        profileString[i] = maxChar + FIRST_ALPHABET_CHAR;
    }
    fprintf(fileHandle, "\tSEQUENCE REF_NAME: %s REF_START %"
            PRIi64 " REF_LENGTH: %" PRIi64 " ML_STRING: %s\n",
            seq->referenceName, seq->refStart, seq->length, profileString);
    if(includeProbs) {
        // FIXME
        for(int64_t i=0; i<seq->length; i++) {
//        for(int64_t i=0; i<26; i++) {
            int64_t o = seq->refCoords[i];
            fprintf(fileHandle, "\t%"PRIi64 "", o);
            uint8_t *p = &seq->profileProbs[i * ALPHABET_SIZE];
            // Print individual character probs
            fprintf(fileHandle, "\t\tPOS: %" PRIi64 "", i);
            for(int64_t j=0; j<ALPHABET_SIZE; j++) {
                fprintf(fileHandle, "\t%f\t", getProb(p, j));
            }
            fprintf(fileHandle, "\n");
        }
    }

}

void printSeqs(FILE *fileHandle, stSet *profileSeqs) {
    /*
     * Prints a set of profile seqs.
     */
    stSetIterator *seqIt = stSet_getIterator(profileSeqs);
    stProfileSeq *pSeq;
    while((pSeq = stSet_getNext(seqIt)) != NULL) {
        stProfileSeq_print(pSeq, fileHandle, 0);
    }
    stSet_destructIterator(seqIt);
}

void printPartition(FILE *fileHandle, stSet *profileSeqs1, stSet *profileSeqs2) {
    /*
     * Print a partition.
     */
    fprintf(fileHandle, "First partition\n");
    printSeqs(fileHandle, profileSeqs1);
    fprintf(fileHandle, "Second partition\n");
    printSeqs(fileHandle, profileSeqs2);
}

static int cmpint64(int64_t i, int64_t j) {
    return i > j ? 1 : i < j ? -1 : 0;
}

inline int stRPProfileSeq_cmpFn(const void *a, const void *b) {
    /*
     * Compares two profile sequences by coordinate on the reference.
     * If two profile sequences have the same reference coordinates
     * will return equal only if they are the same profile sequence, with the same memory
     * address, otherwise compares pointers for equal profile sequences.
     */
    stProfileSeq *pSeq1 = (stProfileSeq *)a, *pSeq2 = (stProfileSeq *)b;
    int i = strcmp(pSeq1->referenceName, pSeq2->referenceName);
    if(i == 0) {
        i = cmpint64(pSeq1->refStart,  pSeq2->refStart);
        if(i == 0) {
            i = cmpint64(pSeq1->length,  pSeq2->length);
            if(i == 0) {
                i = pSeq1 > pSeq2 ? 1 : (pSeq1 < pSeq2 ? -1 : 0);
            }
        }
    }
    return i;
}

void addProfileSeqIdsToSet(stSet *pSeqs, stSet *readIds) {
    /*
     * Add the set of readIDs of a set of profile sequences to a separate set.
     */
    stSetIterator *it = stSet_getIterator(pSeqs);
    stProfileSeq *pSeq;
    while((pSeq = stSet_getNext(it)) != NULL) {
        stSet_insert(readIds, pSeq->readId);
    }
    stSet_destructIterator(it);
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

int64_t numInsertionColumnsInSeq(stProfileSeq *pSeq, stReferencePriorProbs *rProbs, int64_t threshold) {
    /*
     * Counts the number of columns that will be added to the profile sequence for insertions.
     * The number of insertions seen at any location must be at least the threshold to be added.
     */
    int64_t numColumns = 0;
    for (int64_t i = 0; i < pSeq->length; i++) {
        int64_t j = findCorrespondingRefCoordIndex(i, pSeq->refCoords, rProbs->refCoordMap);
//        int64_t j = i + pSeq->refStart - rProbs->refStart;
        if (j >= 0 && j < rProbs->length) {
            if (rProbs->insertionCounts[j] >= threshold) {
                numColumns += rProbs->gapSizes[j];
            }
        }
    }
    return numColumns;
}

void addInsertedBases(stReferencePriorProbs *rProbs, stProfileSeq *seq1, stProfileSeq *seq2, int64_t seq1Index, int64_t seq2Index, int64_t rProbsIndex, int64_t *indexes) {
    /*
     * Adds bases to be inserted (stored in seq1 insertionSeqs) into seq2.
     */
    for (int64_t i = 0; i < rProbs->gapSizes[rProbsIndex]; i++) {
        if (i < seq1->insertions[seq1Index]) {
            int64_t b = seq1->insertionSeqs[seq1Index][i];
            seq2->profileProbs[(seq2Index+1) * ALPHABET_SIZE + b] = ALPHABET_MAX_PROB;
        } else {
            // Chosen gap size was bigger, add extra gaps
            seq2->profileProbs[(seq2Index+1) * ALPHABET_SIZE + (ALPHABET_SIZE - 1)]
                    = ALPHABET_MAX_PROB;
        }
        // Reference coordinates for the sequence all refer back to the start of the gap
        seq2->refCoords[seq2Index+1] = seq1Index + seq1->refStart;
        seq2->insertions[seq2Index+1] = 1;
        seq2Index++;
    }
}

stProfileSeq *stProfileSeq_constructProfileWithInsertions(stProfileSeq *pSeq, stReferencePriorProbs *rProbs, int64_t threshold) {
    /*
     * Construct a new profile sequence that is a copy of pSeq, but with insertion columns (above
     * the correct threshold) added.
     */
    int64_t columnsToAdd = numInsertionColumnsInSeq(pSeq, rProbs, threshold);
    int64_t newLength = pSeq->length + columnsToAdd;
    stProfileSeq *insertionSeq = stProfileSeq_constructEmptyProfile(pSeq->referenceName, pSeq->readId, pSeq->refStart, newLength);
    int64_t insertionSeqIndex = 0;
    int64_t *indexes = st_calloc(newLength, sizeof(int64_t));
    indexes[0] = 0;

    for (int64_t i = 0; i < pSeq->length; i++) {

        int64_t rProbsIndex = findCorrespondingRefCoordIndex(i, pSeq->refCoords, rProbs->refCoordMap);
        if(rProbsIndex >= 0 && rProbsIndex < rProbs->length) {
            // Add existing character in column
            for (int64_t k = 0; k < ALPHABET_SIZE; k++) {
                insertionSeq->profileProbs[insertionSeqIndex * ALPHABET_SIZE + k] =
                        pSeq->profileProbs[i * ALPHABET_SIZE + k];
            }

            insertionSeq->refCoords[insertionSeqIndex] = pSeq->refCoords[i];
            indexes[insertionSeqIndex] = insertionSeqIndex;
            stHash_insert(insertionSeq->refCoordMap,
                          &insertionSeq->refCoords[insertionSeqIndex], &indexes[insertionSeqIndex]);
            insertionSeq->insertions[insertionSeqIndex] = 0;
            if(pSeq->insertions[i] > 0 && rProbs->insertionCounts[rProbsIndex] >= threshold) {
                // Insertion in both this sequence and others
                addInsertedBases(rProbs, pSeq, insertionSeq,
                                 i, insertionSeqIndex, rProbsIndex, indexes);
                insertionSeqIndex += rProbs->gapSizes[rProbsIndex] + 1;
            } else if (pSeq->insertions[i] > 0) {
                // Number of insertions didn't exceed threshold, don't add this spot
                insertionSeqIndex++;
            } else if (rProbs->insertionCounts[rProbsIndex] >= threshold) {
                // Insertions in other sequences, but not in this one - make it a gap
                // TODO: this assumes gap is always last character in alphabet
                for (int64_t j = 0; j < rProbs->gapSizes[rProbsIndex]; j++) {
                    insertionSeq->profileProbs[(insertionSeqIndex+1) * ALPHABET_SIZE + (ALPHABET_SIZE - 1)]
                            = ALPHABET_MAX_PROB;
                    insertionSeq->refCoords[insertionSeqIndex+1] = pSeq->refCoords[i];
                    insertionSeq->insertions[insertionSeqIndex+1] = 1;
                    insertionSeqIndex++;
                }
                insertionSeqIndex++;
            } else {
                // No insertion in this profile sequence nor in others
                insertionSeqIndex++;
            }
        }
    }
    insertionSeq->refEnd = insertionSeq->refCoords[insertionSeq->length - 1];
    return insertionSeq;
}

stList *addInsertionColumnsToSeqs(stList *profileSequences, stHash *referenceNamesToReferencePriors, int64_t threshold, int64_t *numInsertions) {
    /*
     * Modify the set of profile sequences to insert new columns wherever there is likely an
     * insertion relative to the reference.
     */
    stList *insertionPSeqs = stList_construct3(0, (void (*)(void *))stProfileSeq_destruct);
    stHashIterator *hashIt = stHash_getIterator(referenceNamesToReferencePriors);
    char *referenceName;
    while ((referenceName = stHash_getNext(hashIt)) != NULL) {
        stReferencePriorProbs *rProbs = stHash_search(referenceNamesToReferencePriors, referenceName);
        for (int64_t i = 0; i < stList_length(profileSequences); i++) {
            stProfileSeq *pSeq = stList_get(profileSequences, i);
            stProfileSeq *insertionSeq = stProfileSeq_constructProfileWithInsertions(pSeq, rProbs, threshold);
            stList_append(insertionPSeqs, insertionSeq);
        }
        (*numInsertions) += numTotalInsertionColumns(rProbs, threshold);
    }
    // Destroy old profile sequences
    stList_destruct(profileSequences);
    return insertionPSeqs;
}

int64_t gapSizeAtIndex(int64_t *refCoords, int64_t index) {
    /*
     * Returns the offset from the index in the profile sequence to the first spot with the
     * same reference coordinate (the size of the current gap up until the index)
     */
    int64_t gapSize = 0;
    while(index > 0 && refCoords[index] == refCoords[index - 1]) {
        gapSize++;
        index--;
    }
    return gapSize;
}


int64_t findCorrespondingRefCoordIndex(int64_t index1, int64_t *refCoords1, stHash *refCoordMap2) {
    /*
     * Given an index for an object that has reference coordinates,
     * return the index of the corresponding location found in another object with reference coordinates.
     */
    int64_t *idxPtr = stHash_search(refCoordMap2, &refCoords1[index1]);
    int64_t index2;
    if (idxPtr != NULL) {
        index2 = *idxPtr;
        int64_t gapSize = gapSizeAtIndex(refCoords1, index1);
        index2 += gapSize;
        return index2;
    } else {
        return -1;
    }
}

