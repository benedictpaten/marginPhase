/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "margin.h"

/*
 * Functions for profile sequence
 */

stProfileSeq *stProfileSeq_constructEmptyProfile(stReference *ref, char *readId,
		int64_t referenceStart, int64_t length) {
    /*
     * Creates an empty profile sequence, with all the profile probabilities set to 0.
     */

    stProfileSeq *seq = st_calloc(1, sizeof(stProfileSeq));
    seq->ref = ref;
    seq->readId = stString_copy(readId);
    seq->refStart = referenceStart;
    seq->length = length;
    seq->alleleOffset = ref->sites[referenceStart].alleleOffset;
    uint64_t lastAllele = referenceStart+length < ref->length ? ref->sites[referenceStart+length].alleleOffset : ref->totalAlleles;
    seq->profileProbs = st_calloc(lastAllele-seq->alleleOffset, sizeof(uint8_t));
    return seq;
}

void stProfileSeq_destruct(stProfileSeq *seq) {
    /*
     * Cleans up memory for profile sequence.
     */

    free(seq->profileProbs);
    free(seq->readId);
    free(seq);
}

uint8_t *stProfileSeq_getProb(stProfileSeq *seq, uint64_t site, uint64_t allele) {
    /*
     * Gets probability of a given character as a float.
     */
    return &(seq->profileProbs[seq->ref->sites[site].alleleOffset - seq->alleleOffset + allele]); //((float)p[characterIndex])/255;
}

void stProfileSeq_print(stProfileSeq *seq, FILE *fileHandle) {
    /*
     * Prints a debug representation of a profile sequence.
     */

    char profileString[seq->length+1];
    profileString[seq->length] = '\0';
    for(int64_t i=0; i<seq->length; i++) {
    	stSite *site = &(seq->ref->sites[i+seq->refStart]);
        uint8_t maxProb = seq->profileProbs[site->alleleOffset - seq->alleleOffset];
        int64_t maxAllele = 0;
        for(int64_t j=1; j<site->alleleNumber; j++) {
        	uint8_t prob = seq->profileProbs[site->alleleOffset - seq->alleleOffset + j];
            if(prob < maxProb) {
                maxProb = prob;
                maxAllele = j;
            }
        }
        //TODO this prints weird boxes 'cause this isn't a char
        profileString[i] = maxAllele;
    }

    fprintf(fileHandle, "\tSEQUENCE REF_NAME: %s REF_START %"
            PRIi64 " REF_LENGTH: %" PRIi64 " ML_STRING: %s\n",
            seq->ref->referenceName, seq->refStart, seq->length, profileString);
}

void printSeqs(FILE *fileHandle, stSet *profileSeqs) {
    /*
     * Prints a set of profile seqs.
     */

    stSetIterator *seqIt = stSet_getIterator(profileSeqs);
    stProfileSeq *pSeq;
    while((pSeq = stSet_getNext(seqIt)) != NULL) {
        stProfileSeq_print(pSeq, fileHandle);
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
    int i = strcmp(pSeq1->ref->referenceName, pSeq2->ref->referenceName);
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
    stSetIterator *it = stSet_getIterator(pSeqs);
    stProfileSeq *pSeq;
    while((pSeq = stSet_getNext(it)) != NULL) {
        stSet_insert(readIds, pSeq->readId);
    }
    stSet_destructIterator(it);
}

