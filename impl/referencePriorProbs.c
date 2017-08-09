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
    referencePriorProbs->referenceSequence = st_calloc(length, sizeof(uint8_t));
    referencePriorProbs->baseCounts = st_calloc(length*ALPHABET_SIZE, sizeof(double));
    referencePriorProbs->referencePositionsIncluded = st_calloc(length, sizeof(bool));
    for(int64_t i=0; i<length; i++) {
        referencePriorProbs->referencePositionsIncluded[i] = true;
    }
    return referencePriorProbs;
}

void stReferencePriorProbs_destruct(stReferencePriorProbs *referencePriorProbs) {
    free(referencePriorProbs->profileProbs);
    free(referencePriorProbs->referenceName);
    free(referencePriorProbs->referenceSequence);
    free(referencePriorProbs->baseCounts);
    free(referencePriorProbs->referencePositionsIncluded);
    free(referencePriorProbs);
}

static stReferencePriorProbs *getNext(stList *profileSequences, stList *profileSequencesOnReference) {
    /*
     * Construct an empty stReferencePriorProbs for the next interval of a reference sequence.
     */
    assert(stList_length(profileSequences) > 0);
    stProfileSeq *pSeq = stList_pop(profileSequences);
    char *refName = pSeq->referenceName;
    int64_t refStart = pSeq->refStart;
    int64_t refEnd = pSeq->refStart + pSeq->length;

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
        refEnd = pSeq->refStart + pSeq->length > refEnd ? pSeq->refStart + pSeq->length : refEnd;
    }

    return stReferencePriorProbs_constructEmptyProfile(refName, refStart, refEnd-refStart);
}

void setReadBaseCounts(stReferencePriorProbs *rProbs, stList *profileSequences) {
    /*
     * Set the base counts observed in the set of profile sequences at each reference position.
     */
    for(int64_t i=0; i<stList_length(profileSequences); i++) {
        stProfileSeq *pSeq = stList_get(profileSequences, i);
        assert(stString_eq(rProbs->referenceName, pSeq->referenceName));
        for (int64_t j = 0; j < pSeq->length; j++) {
            int64_t k = j + pSeq->refStart - rProbs->refStart;
            assert(k >= 0 && k < rProbs->length);
            for (int64_t l = 0; l < ALPHABET_SIZE; l++) {
                rProbs->baseCounts[k*ALPHABET_SIZE + l] += getProb(&(pSeq->profileProbs[j*ALPHABET_SIZE]), l);
            }
        }
    }
}

stHash *createEmptyReferencePriorProbabilities(stList *profileSequences) {
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
    profileSequences = stList_copy(profileSequences, NULL);
    stList_sort(profileSequences, stRPProfileSeq_cmpFn);

    uint16_t flatValue = scaleToLogIntegerSubMatrix(log(1.0/ALPHABET_SIZE));
    while(stList_length(profileSequences) > 0) {
        // Get next reference prior prob interval
        stList *l = stList_construct();
        stReferencePriorProbs *rProbs = getNext(profileSequences, l);

        // Fill in probability profile
        for(int64_t i=0; i<rProbs->length; i++) {
            for(int64_t j=0; j<ALPHABET_SIZE; j++) {
                rProbs->profileProbs[i*ALPHABET_SIZE + j] = flatValue;
            }
        }

        // Set the counts of each base observed at each reference position
        setReadBaseCounts(rProbs, l);
        stList_destruct(l);

        assert(stHash_search(referenceNamesToReferencePriors, rProbs->referenceName) == NULL);
        stHash_insert(referenceNamesToReferencePriors, rProbs->referenceName, rProbs);
    }

    // Cleanup
    stList_destruct(profileSequences);

    return referenceNamesToReferencePriors;
}

stHash *createReferencePriorProbabilities(char *referenceFastaFile, stList *profileSequences,
        stBaseMapper *baseMapper, stRPHmmParameters *params) {
    /*
     * Create a set of stReferencePriorProbs that cover the reference intervals included in the profile sequences.
     * The return value is encoded as a map from the reference sequence name (as a string)
     * to the stReferencePriorProbs.
     */

    // Make map from reference sequence names to reference priors
    stHash *referenceNamesToReferencePriors = createEmptyReferencePriorProbabilities(profileSequences);

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
            uint8_t refChar = stBaseMapper_getValueForChar(baseMapper, referenceSeq[i+rProbs->refStart-1]);
            assert(refChar >= 0 && refChar < ALPHABET_SIZE);
            rProbs->referenceSequence[i] = refChar;
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

double *stReferencePriorProbs_estimateReadErrorProbs(stHash *referenceNamesToReferencePriors, stRPHmmParameters *params) {
    /*
     * Estimate the read error substitution matrix from the alignment of the reads to the reference.
     */
    
    // Make read error substitution matrix
    double *readErrorSubModel = getEmptyReadErrorSubstitutionMatrix(params);

    stHashIterator *rProbsIt = stHash_getIterator(referenceNamesToReferencePriors);
    char *refSeq;
    while((refSeq = stHash_getNext(rProbsIt)) != NULL) {
        stReferencePriorProbs *rProbs = stHash_search(referenceNamesToReferencePriors, refSeq);
        for(int64_t i=0; i<rProbs->length; i++) {
            // Reference character
            uint8_t refChar = rProbs->referenceSequence[i];

            // Now get read bases
            double *baseCounts = &rProbs->baseCounts[i*ALPHABET_SIZE];

            //
            for(int64_t j=0; j<ALPHABET_SIZE; j++) {
                readErrorSubModel[ALPHABET_SIZE*refChar + j] += baseCounts[j];
            }
        }
    }
    stHash_destructIterator(rProbsIt);
    
    normaliseSubstitutionMatrix(readErrorSubModel);

    return readErrorSubModel;
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
        int64_t secondMostFrequentBase = 0;
        for(int64_t j=0; j<ALPHABET_SIZE; j++) {
            if(j != mostFrequentBase && baseCounts[j] > secondMostFrequentBaseCount) {
                secondMostFrequentBaseCount = baseCounts[j];
                secondMostFrequentBase = j;
            }
        }

        // Filter columns with either of two options:
        // (1) If the second most frequent character in the column occurs less than
        // params->minSecondMostFrequentBaseFilter times; this is
        // intentionally a simple and somewhat robust filter
        // (2) If the -log probability of observing the second most frequent character as a result
        // of independent substitutions from the
        // most frequent character is less than params->minSecondMostFrequentBaseLogProbFilter.
        // This second options takes
        // into account read error substitution probabilities. It can be used to ignore columns with
        // larger numbers of gaps, for example, if
        // gaps occur frequently, while including columns with rarer apparent substitution patterns
        if(secondMostFrequentBaseCount < params->minSecondMostFrequentBaseFilter ||
           -*getSubstitutionProbSlow(params->readErrorSubModelSlow, mostFrequentBase, secondMostFrequentBase) * secondMostFrequentBaseCount <
           params->minSecondMostFrequentBaseLogProbFilter) {
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
//                        rProbs->referenceName, rProbs->refStart + i, totalCount, baseCountString);
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
