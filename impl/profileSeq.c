/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "margin.h"

/*
 * Functions for profile sequence
 */

stProfileSeq *stProfileSeq_constructEmptyProfile(char *referenceName, char *readId, int64_t referenceStart,
        int64_t length) {
    /*
     * Creates an empty profile sequence, with all the profile probabilities set to 0.
     */

    stProfileSeq *seq = st_malloc(sizeof(stProfileSeq));
    seq->referenceName = stString_copy(referenceName);
    seq->readId = stString_copy(readId);
    seq->refStart = referenceStart;
    seq->length = length;
    seq->profileProbs = st_calloc(length*ALPHABET_SIZE, sizeof(uint8_t));
    return seq;
}

void stProfileSeq_destruct(stProfileSeq *seq) {
    /*
     * Cleans up memory for profile sequence.
     */

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
        for(int64_t i=0; i<seq->length; i++) {
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
    stSetIterator *it = stSet_getIterator(pSeqs);
    stProfileSeq *pSeq;
    while((pSeq = stSet_getNext(it)) != NULL) {
        stSet_insert(readIds, pSeq->readId);
    }
    stSet_destructIterator(it);
}

/*
 * Create profile sequence by summing over the alignment of a read to a given reference sequence,
 * giving the probability of each base for each reference position.
 */
stProfileSeq *stProfileSeq_constructFromPosteriorProbs(char *refName, char *refSeq, int64_t refLength,
													   char *readId, char *readSeq, stList *anchorAlignment,
													   Params *params) {

	// Generate the posterior probabilities
	stList *matches = NULL, *inserts = NULL, *deletes = NULL;
	getAlignedPairsWithIndelsCroppingReference(refSeq, refLength, readSeq, anchorAlignment, &matches, &inserts, &deletes, params->polishParams);

	// Get min and max reference coordinates
	int64_t refStart, refEnd;
	if(stList_length(matches) > 0) {
		refStart=refLength, refEnd=0;
		for(int64_t i=0; i<stList_length(matches); i++) {
			stIntTuple *aPair = stList_get(matches, i);
			int64_t j = stIntTuple_get(aPair, 1);
			refStart = j < refStart ? j : refStart;
			refEnd = j > refEnd ? j : refEnd;
		}
	}
	else {
		refStart = 0;
		refEnd = refLength-1;
	}
	assert(refEnd+1 >= refStart);

	int64_t *probs = st_calloc((refEnd+1-refStart)*ALPHABET_SIZE, sizeof(int64_t));

	// Add posterior match probabilities
	for(int64_t i=0; i<stList_length(matches); i++) {
		stIntTuple *aPair = stList_get(matches, i);
		char base = readSeq[stIntTuple_get(aPair, 2)];
		assert(stIntTuple_get(aPair, 1) >= refStart);
		assert(stIntTuple_get(aPair, 1) <= refEnd);
		probs[(stIntTuple_get(aPair, 1)-refStart) * ALPHABET_SIZE + stBaseMapper_getValueForChar(params->baseMapper, base)] += stIntTuple_get(aPair, 0);
	}

	// Add posterior delete probs
	if(params->phaseParams->gapCharactersForDeletions) {
		for(int64_t i=0; i<stList_length(deletes); i++) {
			stIntTuple *dPair = stList_get(deletes, i);
			int64_t j = stIntTuple_get(dPair, 1);
			if(j >= refStart && j <= refEnd) {
				probs[(stIntTuple_get(dPair, 1)-refStart) * ALPHABET_SIZE + stBaseMapper_getValueForChar(params->baseMapper, '-')] += stIntTuple_get(dPair, 0);
			}
		}
	}

	// Make empty profile sequence
	stProfileSeq *pSeq = stProfileSeq_constructEmptyProfile(refName, readId, refStart,
		        											refEnd-refStart+1);
	for(int64_t i=0; i<refEnd+1-refStart; i++) { // Add the probs
		int64_t totalProb = 0;
		for(int64_t j=0; j<ALPHABET_SIZE; j++) {
			totalProb += probs[i*ALPHABET_SIZE + j];
		}
		for(int64_t j=0; j<ALPHABET_SIZE; j++) {
			pSeq->profileProbs[i*ALPHABET_SIZE + j] = totalProb > 0 ?
					round(ALPHABET_MAX_PROB * ((double)probs[i*ALPHABET_SIZE + j]/totalProb)) : ALPHABET_MIN_PROB; //round((double)ALPHABET_MAX_PROB / ALPHABET_SIZE);
		}
	}

	// Cleanup
	free(probs);
	stList_destruct(matches);
	stList_destruct(inserts);
	stList_destruct(deletes);

	return pSeq;
}
