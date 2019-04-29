//
// Created by tpesout on 3/29/19.
//

#include "margin.h"
#include "helenFeatures.h"

#ifdef _HDF5
#include <hdf5.h>
#endif


PoaFeatureSimpleWeight *PoaFeature_SimpleWeight_construct(int64_t refPos, int64_t insPos) {
    PoaFeatureSimpleWeight *feature = st_calloc(1, sizeof(PoaFeatureSimpleWeight));
    feature->refPosition = refPos;
    feature->insertPosition = insPos;
    feature->label = '\0';
    feature->nextInsert = NULL;
    return feature;
}

void PoaFeature_SimpleWeight_destruct(PoaFeatureSimpleWeight *feature) {
    if (feature->nextInsert != NULL) {
        PoaFeature_SimpleWeight_destruct(feature->nextInsert);
    }
    free(feature);
}


PoaFeatureRleWeight *PoaFeature_RleWeight_construct(int64_t refPos, int64_t insPos) {
    PoaFeatureRleWeight *feature = st_calloc(1, sizeof(PoaFeatureRleWeight));
    feature->refPosition = refPos;
    feature->insertPosition = insPos;
    feature->labelChar = '\0';
    feature->labelRunLength = 0;
    feature->predictedRunLength = 0;
    feature->nextInsert = NULL;
    return feature;
}

void PoaFeature_RleWeight_destruct(PoaFeatureRleWeight *feature) {
    if (feature->nextInsert != NULL) {
        PoaFeature_RleWeight_destruct(feature->nextInsert);
    }
    free(feature);
}


PoaFeatureSplitRleWeight *PoaFeature_SplitRleWeight_construct(int64_t refPos, int64_t insPos, int64_t rlPos,
        int64_t maxRunLength) {
    PoaFeatureSplitRleWeight *feature = st_calloc(1, sizeof(PoaFeatureSplitRleWeight));
    feature->refPosition = refPos;
    feature->insertPosition = insPos;
    feature->runLengthPosition = rlPos;
    feature->labelChar = '\0';
    feature->labelRunLength = 0;
    feature->nextRunLength = NULL;
    feature->nextInsert = NULL;
    feature->maxRunLength = maxRunLength;
    feature->weights = st_calloc(((SYMBOL_NUMBER - 1) * (1 + maxRunLength) + 1) * 2, sizeof(double));
    return feature;
}
void PoaFeature_SplitRleWeight_destruct(PoaFeatureSplitRleWeight *feature) {
    if (feature->nextRunLength != NULL) {
        PoaFeature_SplitRleWeight_destruct(feature->nextRunLength);
    }
    if (feature->nextInsert != NULL) {
        PoaFeature_SplitRleWeight_destruct(feature->nextInsert);
    }
    free(feature->weights);
    free(feature);
}

int PoaFeature_SimpleWeight_charIndex(Symbol character, bool forward) {
    int pos = character * 2 + (forward ? POS_STRAND_IDX : NEG_STRAND_IDX);
    assert(pos < POAFEATURE_SIMPLE_WEIGHT_TOTAL_SIZE);
    return pos;
}
int PoaFeature_SimpleWeight_gapIndex(bool forward) {
    int pos = POAFEATURE_SYMBOL_GAP_POS * 2 + (forward ? POS_STRAND_IDX : NEG_STRAND_IDX);
    assert(pos < POAFEATURE_SIMPLE_WEIGHT_TOTAL_SIZE);
    return pos;
}
int PoaFeature_RleWeight_charIndex(Symbol character, int64_t runLength, bool forward) {
    assert(runLength > 0);
    assert(runLength <= POAFEATURE_MAX_RUN_LENGTH);
    runLength -= 1;
    int pos = (character * POAFEATURE_MAX_RUN_LENGTH + runLength) * 2 + (forward ? POS_STRAND_IDX : NEG_STRAND_IDX);
    assert(pos < POAFEATURE_RLE_WEIGHT_TOTAL_SIZE);
    return pos;
}
int PoaFeature_RleWeight_gapIndex(bool forward) {
    int pos = (POAFEATURE_SYMBOL_GAP_POS * POAFEATURE_MAX_RUN_LENGTH) * 2 + (forward ? POS_STRAND_IDX : NEG_STRAND_IDX);
    assert(pos < POAFEATURE_RLE_WEIGHT_TOTAL_SIZE);
    return pos;
}
int PoaFeature_SplitRleWeight_charIndex(int64_t maxRunLength, Symbol character, int64_t runLength, bool forward) {
    assert(runLength >= 0);
    assert(runLength <= maxRunLength);
    int pos = (character * ((int)maxRunLength + 1) + runLength) * 2 + (forward ? POS_STRAND_IDX : NEG_STRAND_IDX);
    return pos;
}
int PoaFeature_SplitRleWeight_gapIndex(int64_t maxRunLength, bool forward) {
    int pos = ((SYMBOL_NUMBER - 1) * ((int)maxRunLength + 1)) * 2 + (forward ? POS_STRAND_IDX : NEG_STRAND_IDX);
    return pos;

}

stList *poa_getSimpleWeightFeatures(Poa *poa, stList *bamChunkReads) {

    // initialize feature list
    stList *featureList = stList_construct3(0, (void (*)(void *)) PoaFeature_SimpleWeight_destruct);
    for(int64_t i=1; i<stList_length(poa->nodes); i++) {
        stList_append(featureList, PoaFeature_SimpleWeight_construct(i - 1, 0));
    }

    // for logging (of errors)
    char *logIdentifier = getLogIdentifier();

    // iterate over all positions
    for(int64_t i=0; i<stList_length(featureList); i++) {

        // get feature and node
        PoaFeatureSimpleWeight* feature = stList_get(featureList, i);
        PoaNode *node = stList_get(poa->nodes, i + 1); //skip the first poa node, as it's always an 'N', so featureIdx and poaIdx are off by one

        // get weights for all bases, strands
        double totalWeight, totalPositiveWeight, totalNegativeWeight;
        double *baseWeights = poaNode_getStrandSpecificBaseWeights(node, bamChunkReads,
                                                                   &totalWeight, &totalPositiveWeight, &totalNegativeWeight);
        for(int64_t j=0; j<SYMBOL_NUMBER; j++) {
            double positiveStrandBaseWeight = baseWeights[j * 2 + POS_STRAND_IDX];
            double negativeStrandBaseWeight = baseWeights[j * 2 + NEG_STRAND_IDX];
            feature->weights[PoaFeature_SimpleWeight_charIndex(j, TRUE)] += positiveStrandBaseWeight;
            feature->weights[PoaFeature_SimpleWeight_charIndex(j, FALSE)] += negativeStrandBaseWeight;
        }
        free(baseWeights);

        // Deletes
        if (stList_length(node->deletes) > 0) {

            // iterate over all deletes
            for (int64_t j = 0; j < stList_length(node->deletes); j++) {
                PoaDelete *delete = stList_get(node->deletes, j);

                // Deletes start AFTER the current position, need to add counts/weights to nodes after the current node
                for (int64_t k = 1; k < delete->length; k++) {
                    if (i + k >= stList_length(poa->nodes)) {
                        st_logInfo(" %s Encountered DELETE that occurs after end of POA!\n", logIdentifier);
                        break;
                    }
                    PoaFeatureSimpleWeight *delFeature = stList_get(featureList, i + k);

                    delFeature->weights[PoaFeature_SimpleWeight_gapIndex(TRUE)] += delete->weightForwardStrand;
                    delFeature->weights[PoaFeature_SimpleWeight_gapIndex(FALSE)] += delete->weightReverseStrand;
                }
            }
        }

        // Inserts
        if (stList_length(node->inserts) > 0) {

            // iterate over all inserts
            for (int64_t j = 0; j < stList_length(node->inserts); j++) {
                PoaInsert *insert = stList_get(node->inserts, j);

                // get feature iterator
                PoaFeatureSimpleWeight *prevFeature = feature;
                for (int64_t k = 0; k < strlen(insert->insert); k++) {
                    // get current feature (or create if necessary)
                    PoaFeatureSimpleWeight *currFeature = prevFeature->nextInsert;
                    if (currFeature == NULL) {
                        currFeature = PoaFeature_SimpleWeight_construct(i, k + 1);
                        prevFeature->nextInsert = currFeature;
                    }

                    Symbol c = symbol_convertCharToSymbol(insert->insert[k]);

                    // add weights
                    currFeature->weights[PoaFeature_SimpleWeight_charIndex(c, TRUE)] += insert->weightForwardStrand;
                    currFeature->weights[PoaFeature_SimpleWeight_charIndex(c, FALSE)] += insert->weightReverseStrand;

                    // iterate
                    prevFeature = currFeature;
                }
            }
        }
    }

    free(logIdentifier);
    return featureList;
}

stList *poa_getRleWeightFeatures(Poa *poa, stList *bamChunkReads, stList *rleStrings, RleString *consensusRleString) {

    // initialize feature list
    stList *featureList = stList_construct3(0, (void (*)(void *)) PoaFeature_RleWeight_destruct);
    for(int64_t i=1; i<stList_length(poa->nodes); i++) {
        stList_append(featureList, PoaFeature_RleWeight_construct(i - 1, 0));
    }

    // for logging (of errors)
    char *logIdentifier = getLogIdentifier();

    // iterate over all positions
    for(int64_t i=0; i<stList_length(featureList); i++) {

        // get feature and node
        PoaFeatureRleWeight* feature = stList_get(featureList, i);
        PoaNode *node = stList_get(poa->nodes, i + 1); //skip the first poa node, as it's always an 'N', so featureIdx and poaIdx are off by one

        // handle each observations
        for (int64_t o = 0; o < stList_length(node->observations); o++) {
            // get objects we need
            PoaBaseObservation *observation = stList_get(node->observations, o);
            BamChunkRead *bamChunkRead = stList_get(bamChunkReads, observation->readNo);
            RleString *rleString = stList_get(rleStrings, observation->readNo);

            // save weight based on character and runLength
            Symbol character = symbol_convertCharToSymbol(rleString->rleString[observation->offset]);
            int64_t runLength = rleString->repeatCounts[observation->offset];
            if (runLength == 0) continue;
            if (runLength > POAFEATURE_MAX_RUN_LENGTH) runLength = POAFEATURE_MAX_RUN_LENGTH;
            feature->weights[PoaFeature_RleWeight_charIndex(character, runLength, bamChunkRead->forwardStrand)] += observation->weight;
        }
        feature->predictedRunLength = consensusRleString->repeatCounts[i];


        // Deletes
        if (stList_length(node->deletes) > 0) {

            // iterate over all deletes
            for (int64_t d = 0; d < stList_length(node->deletes); d++) {
                PoaDelete *delete = stList_get(node->deletes, d);

                // Deletes start AFTER the current position, need to add counts/weights to nodes after the current node
                for (int64_t k = 1; k < delete->length; k++) {
                    if (i + k >= stList_length(poa->nodes)) {
                        st_logInfo(" %s Encountered DELETE that occurs after end of POA!\n", logIdentifier);
                        break;
                    }
                    PoaFeatureRleWeight *delFeature = stList_get(featureList, i + k);

                    delFeature->weights[PoaFeature_RleWeight_gapIndex(TRUE)] += delete->weightForwardStrand;
                    delFeature->weights[PoaFeature_RleWeight_gapIndex(FALSE)] += delete->weightReverseStrand;
                }
            }
        }

        // Inserts
        if (stList_length(node->inserts) > 0) {

            // iterate over all inserts
            for (int64_t n = 0; n < stList_length(node->inserts); n++) {
                PoaInsert *insert = stList_get(node->inserts, n);

                // handle each observation
                for (int64_t o = 0; o < stList_length(insert->observations); o++) {
                    // get objects we need
                    PoaBaseObservation *observation = stList_get(insert->observations, o);
                    BamChunkRead *bamChunkRead = stList_get(bamChunkReads, observation->readNo);
                    RleString *rleString = stList_get(rleStrings, observation->readNo);

                    // get feature iterator
                    PoaFeatureRleWeight *prevFeature = feature;

                    // iterate over all positions in insert
                    for (int64_t k = 0; k < strlen(insert->insert); k++) {
                        // get current feature (or create if necessary)
                        PoaFeatureRleWeight *currFeature = prevFeature->nextInsert;
                        if (currFeature == NULL) {
                            currFeature = PoaFeature_RleWeight_construct(i, k + 1);
                            prevFeature->nextInsert = currFeature;
                        }

                        // get character and runLength
                        int64_t stringPos = observation->offset + k;
                        assert(stringPos < rleString->length);
                        Symbol character = symbol_convertCharToSymbol(rleString->rleString[stringPos]);
                        int64_t runLength = rleString->repeatCounts[stringPos];
                        if (runLength == 0) continue;
                        if (runLength > POAFEATURE_MAX_RUN_LENGTH) runLength = POAFEATURE_MAX_RUN_LENGTH;

                        // save weight for position
                        currFeature->weights[PoaFeature_RleWeight_charIndex(character, runLength,
                                                                            bamChunkRead->forwardStrand)] += observation->weight;
                    }
                }
            }
        }
    }

    free(logIdentifier);
    return featureList;
}


void poa_addRunLengthFeaturesForObservations(PoaFeatureSplitRleWeight *baseFeature, stList *observations,
        stList *bamChunkReads, stList *rleStrings, const int64_t maxRunLength, int64_t observationOffset) {


    PoaFeatureSplitRleWeight *currFeature = baseFeature;
    int64_t currentRunLengthIndex = 0;
    bool beforeMaxObservedRunLength = TRUE;

    while (beforeMaxObservedRunLength) {
        beforeMaxObservedRunLength = FALSE;

        // examine each observation
        for (int64_t i = 0; i < stList_length(observations); i++) {

            // get data
            PoaBaseObservation *observation = stList_get(observations, i);
            RleString *rleString = stList_get(rleStrings, observation->readNo);
            BamChunkRead *bamChunkRead = stList_get(bamChunkReads, observation->readNo);
            Symbol symbol = symbol_convertCharToSymbol(rleString->rleString[observation->offset + observationOffset]);
            int64_t runLength = rleString->repeatCounts[observation->offset + observationOffset];
            bool forward = bamChunkRead->forwardStrand;

            // get correct run length
            runLength -= currentRunLengthIndex * maxRunLength;
            if (runLength < 0) {
                runLength = 0;
            } else if (runLength > maxRunLength) {
                runLength = maxRunLength;
                beforeMaxObservedRunLength = TRUE;
            }

            int64_t pos = PoaFeature_SplitRleWeight_charIndex(maxRunLength, symbol, runLength, forward);
            currFeature->weights[pos] += observation->weight;

        }

        // update currFeature if we're going ot run again
        if (beforeMaxObservedRunLength) {
            currentRunLengthIndex++;
            PoaFeatureSplitRleWeight *prevFeature = currFeature;
            currFeature = PoaFeature_SplitRleWeight_construct(baseFeature->refPosition, baseFeature->insertPosition,
                                                              currentRunLengthIndex, maxRunLength);
            prevFeature->nextRunLength = currFeature;
            currFeature->weights[PoaFeature_SplitRleWeight_gapIndex(maxRunLength, TRUE)] =
                    baseFeature->weights[PoaFeature_SplitRleWeight_gapIndex(maxRunLength, TRUE)];
            currFeature->weights[PoaFeature_SplitRleWeight_gapIndex(maxRunLength, FALSE)] =
                    baseFeature->weights[PoaFeature_SplitRleWeight_gapIndex(maxRunLength, FALSE)];
        }
    }
}


stList *poa_getSplitRleWeightFeatures(Poa *poa, stList *bamChunkReads, stList *rleStrings, const int64_t maxRunLength) {
    // initialize feature list
    stList *featureList = stList_construct3(0, (void (*)(void *)) PoaFeature_SplitRleWeight_destruct);
    for(int64_t i=1; i<stList_length(poa->nodes); i++) {
        stList_append(featureList, PoaFeature_SplitRleWeight_construct(i - 1, 0, 0, maxRunLength));
    }

    // for logging (of errors)
    char *logIdentifier = getLogIdentifier();

    // iterate over all positions
    for(int64_t i=0; i<stList_length(featureList); i++) {

        // get feature and node
        PoaFeatureSplitRleWeight* feature = stList_get(featureList, i);
        PoaNode *node = stList_get(poa->nodes, i + 1); //skip the first poa node, as it's always an 'N', so featureIdx and poaIdx are off by one

        // save run length nodes
        poa_addRunLengthFeaturesForObservations(feature, node->observations, bamChunkReads, rleStrings, maxRunLength, 0);

        // Deletes
        if (stList_length(node->deletes) > 0) {

            // iterate over all deletes
            for (int64_t d = 0; d < stList_length(node->deletes); d++) {
                PoaDelete *delete = stList_get(node->deletes, d);

                // Deletes start AFTER the current position, need to add counts/weights to nodes after the current node
                for (int64_t k = 1; k < delete->length; k++) {
                    if (i + k >= stList_length(poa->nodes)) {
                        st_logInfo(" %s Encountered DELETE that occurs after end of POA!\n", logIdentifier);
                        break;
                    }
                    PoaFeatureSplitRleWeight *delFeature = stList_get(featureList, i + k);

                    delFeature->weights[PoaFeature_SplitRleWeight_gapIndex(maxRunLength, TRUE)] += delete->weightForwardStrand;
                    delFeature->weights[PoaFeature_SplitRleWeight_gapIndex(maxRunLength, FALSE)] += delete->weightReverseStrand;
                }
            }
        }

        // Inserts
        if (stList_length(node->inserts) > 0) {

            // iterate over all inserts
            for (int64_t n = 0; n < stList_length(node->inserts); n++) {
                PoaInsert *insert = stList_get(node->inserts, n);

                // handle each insert base
                PoaFeatureSplitRleWeight *prevFeature = feature;
                for (int64_t o = 0; o < strlen(insert->insert); o++) {

                    // get feature iterator
                    PoaFeatureSplitRleWeight *currFeature = prevFeature->nextInsert;
                    if (currFeature == NULL) {
                        currFeature = PoaFeature_SplitRleWeight_construct(i, o + 1, 0, maxRunLength);
                        prevFeature->nextInsert = currFeature;
                    }

                    // save insert run lengths
                    poa_addRunLengthFeaturesForObservations(currFeature, insert->observations, bamChunkReads, rleStrings,
                            maxRunLength, o);
                }
            }
        }
    }

    free(logIdentifier);
    return featureList;
}


void printMEAAlignment(char *X, char *Y, int64_t lX, int64_t lY, stList *alignedPairs, int64_t *Xrl, int64_t *Yrl) {
    // should we do run lengths
    bool handleRunLength = Xrl != NULL && Yrl != NULL;

    // build strings to print
    int64_t bufferLen = lX + lY;
    char *alnXStr = st_calloc(bufferLen, sizeof(char));
    char *alnYStr = st_calloc(bufferLen, sizeof(char));
    char *alnDesc = st_calloc(bufferLen, sizeof(char));
    char *rlXStr = NULL;
    char *rlYStr = NULL;
    if (handleRunLength) {
        rlXStr = st_calloc(bufferLen, sizeof(char));
        rlYStr = st_calloc(bufferLen, sizeof(char));
    }

    // stats to track
    int64_t nuclMatches = 0;
    int64_t nuclMismatches = 0;
    int64_t nuclXInserts = 0;
    int64_t nuclYInserts = 0;
    int64_t rlMatches = 0;
    int64_t rlMismatches = 0;


    // iterate over alignment
    stListIterator *alignmentItor = stList_getIterator(alignedPairs);
    stIntTuple *currAlign = stList_getNext(alignmentItor);
    int64_t posX = stIntTuple_get(currAlign, 1);
    int64_t posY = stIntTuple_get(currAlign, 2);
    int64_t outStrPos = 0;

    while (TRUE) {
        if (currAlign == NULL) break;
        int64_t currAlignPosX = stIntTuple_get(currAlign, 1);
        int64_t currAlignPosY = stIntTuple_get(currAlign, 2);
        // Y gap / X insert
        if (posX < currAlignPosX) {
            alnXStr[outStrPos] = X[posX];
            alnYStr[outStrPos] = '_';
            alnDesc[outStrPos] = ' ';
            if (handleRunLength) {
                rlXStr[outStrPos] = (char) ('0' + Xrl[posX]);
                rlYStr[outStrPos] = ' ';
            }
            posX++;
            nuclXInserts++;
        }

        // X gap / Y insert
        else if (posY < currAlignPosY) {
            alnXStr[outStrPos] = '_';
            alnYStr[outStrPos] = Y[posY];
            alnDesc[outStrPos] = ' ';
            if (handleRunLength) {
                rlXStr[outStrPos] = ' ';
                rlYStr[outStrPos] = (char) ('0' + Yrl[posY]);
            }
            posY++;
            nuclYInserts++;
        }

        // match
        else if (posX == currAlignPosX && posY == currAlignPosY) {
            alnXStr[outStrPos] = X[posX];
            alnYStr[outStrPos] = Y[posY];
            if (handleRunLength) {
                rlXStr[outStrPos] = (char) ('0' + Xrl[posX]);
                rlYStr[outStrPos] = (char) ('0' + Yrl[posY]);
            }
            if (X[posX] == Y[posY]) {
                nuclMatches++;
                alnDesc[outStrPos] = '|';
                if (handleRunLength) {
                    if (Xrl[posX] == Yrl[posY]) {
                        rlMatches++;
                    } else {
                        rlMismatches++;
                        alnDesc[outStrPos] = ':';
                    }
                }
            } else {
                alnDesc[outStrPos] = ' ';
                nuclMismatches++;
            }
            posX++;
            posY++;
            currAlign = stList_getNext(alignmentItor);
        }

        // should never happen
        else {
            assert(FALSE);
        }

        outStrPos++;
    }

    // print
    fprintf(stderr, "\n");
    if (handleRunLength) fprintf(stderr, "%s\n", rlXStr);
    fprintf(stderr, "%s\n", alnXStr);
    fprintf(stderr, "%s\n", alnDesc);
    fprintf(stderr, "%s\n", alnYStr);
    if (handleRunLength) fprintf(stderr, "%s\n", rlYStr);
    fprintf(stderr, "Matches:    %"PRId64"\n", nuclMatches);
    if (handleRunLength) {
        fprintf(stderr, "  RL Match: %"PRId64"\n", rlMatches);
        fprintf(stderr, "  RL Miss:  %"PRId64"\n", rlMismatches);
    }
    fprintf(stderr, "Mismatches: %"PRId64"\n", nuclMismatches);
    fprintf(stderr, "X Inserts:  %"PRId64"\n", nuclXInserts);
    fprintf(stderr, "Y Inserts:  %"PRId64"\n", nuclYInserts);
    fprintf(stderr, "\n");

    // cleanup
    free(alnXStr);
    free(alnYStr);
    free(alnDesc);
    if (handleRunLength) {
        free(rlXStr);
        free(rlYStr);
    }
    stList_destructIterator(alignmentItor);
}


void poa_annotateHelenFeaturesWithTruth(stList *features, HelenFeatureType featureType, stList *trueRefAlignment,
                                        RleString *trueRefRleString, int64_t *firstMatchedFeaure,
                                        int64_t *lastMatchedFeature) {
    /*
     Each index in features represents a character in the final consensus string (which in turn was a single node in the
     final POA which was used to generate the consensus).  This consensus was aligned against the true reference (which
     is reflected in the trueRefRleString).  Items in the true refAlignment are stIntTuples with index 0 being the
     alignment weight (discarded), index 1 being the position in the consensus string (and features), and index 2 being
     the position in the trueRefRleString.  So we can iterate over features and the true alignment to assign truth
     labels to each feature.
     */
    static int FEATURE_POS = 1;
    static int REFERENCE_POS = 2;
    *firstMatchedFeaure = -1;
    *lastMatchedFeature = -1;

    // iterate over true ref alignment
    stListIterator *trueRefAlignItor = stList_getIterator(trueRefAlignment);
    stIntTuple *currRefAlign = stList_getNext(trueRefAlignItor);

    // iterate over features
    int64_t trueRefPos = stIntTuple_get(currRefAlign, REFERENCE_POS);
    for (int64_t featureRefPos = 0; featureRefPos < stList_length(features); featureRefPos++) {
        void *feature = stList_get(features, featureRefPos);
        void *prevFeature = NULL;
        int64_t trueRunLength = -1;
        PoaFeatureSplitRleWeight* rlFeature = NULL;

        int64_t featureInsPos = 0;
        while (feature != NULL) {

            // no more ref bases, everything is gaps
            if (currRefAlign == NULL) {
                switch (featureType) {
                    case HFEAT_SIMPLE_WEIGHT:
                        ((PoaFeatureSimpleWeight*)feature)->label = '_';
                        feature = ((PoaFeatureSimpleWeight*)feature)->nextInsert;
                        break;
                    case HFEAT_RLE_WEIGHT:
                    case HFEAT_NUCL_AND_RL_WEIGHT:
                        ((PoaFeatureRleWeight*)feature)->labelChar = '_';
                        ((PoaFeatureRleWeight*)feature)->labelRunLength = 0;
                        feature = ((PoaFeatureRleWeight*)feature)->nextInsert;
                        break;
                    case HFEAT_SPLIT_RLE_WEIGHT:
                        rlFeature = ((PoaFeatureSplitRleWeight*)feature);
                        while (rlFeature != NULL) {
                            rlFeature->labelChar = '_';
                            rlFeature->labelRunLength = 0;
                            rlFeature = rlFeature->nextRunLength;
                        }
                        feature = ((PoaFeatureSplitRleWeight*)feature)->nextInsert;
                        break;
                    default:
                        st_errAbort("Unhandled FeatureType!\n");
                }
                continue;
            }

            // sanity checks
            assert(stIntTuple_get(currRefAlign, FEATURE_POS) >= featureRefPos && stIntTuple_get(currRefAlign, REFERENCE_POS) >= trueRefPos);

            // match
            if (stIntTuple_get(currRefAlign, FEATURE_POS) == featureRefPos && stIntTuple_get(currRefAlign, REFERENCE_POS) == trueRefPos) {
                // save label (based on feature type)
                switch (featureType) {
                    case HFEAT_SIMPLE_WEIGHT:
                        ((PoaFeatureSimpleWeight *) feature)->label = trueRefRleString->rleString[trueRefPos];
                        break;
                    case HFEAT_RLE_WEIGHT:
                    case HFEAT_NUCL_AND_RL_WEIGHT:
                        ((PoaFeatureRleWeight*)feature)->labelChar = trueRefRleString->rleString[trueRefPos];
                        trueRunLength = trueRefRleString->repeatCounts[trueRefPos];
                        if (trueRunLength > POAFEATURE_MAX_RUN_LENGTH) trueRunLength = POAFEATURE_MAX_RUN_LENGTH;
                        ((PoaFeatureRleWeight*)feature)->labelRunLength = trueRunLength;
                        break;
                    case HFEAT_SPLIT_RLE_WEIGHT:
                        rlFeature = ((PoaFeatureSplitRleWeight*)feature);
                        trueRunLength = trueRefRleString->repeatCounts[trueRefPos];
                        while (rlFeature != NULL) {
                            rlFeature->labelChar = trueRefRleString->rleString[trueRefPos];
                            if (trueRunLength <= 0) {
                                rlFeature->labelRunLength = 0;
                            } else if (trueRunLength > rlFeature->maxRunLength) {
                                rlFeature->labelRunLength = rlFeature->maxRunLength;
                            } else {
                                rlFeature->labelRunLength = trueRunLength;
                            }
                            trueRunLength -= rlFeature->maxRunLength;
                            rlFeature = rlFeature->nextRunLength;
                        }
                        break;
                    default:
                        st_errAbort("Unhandled FeatureType!\n");
                }

                // iterate
                trueRefPos++;
                currRefAlign = stList_getNext(trueRefAlignItor);
                // handle first and last match
                if (featureInsPos == 0) {
                    if (*firstMatchedFeaure == -1) {
                        *firstMatchedFeaure = featureRefPos;
                    }
                    *lastMatchedFeature = featureRefPos;
                }
            }

            // insert
            else if (trueRefPos < stIntTuple_get(currRefAlign, REFERENCE_POS)) {
                // apply label
                switch (featureType) {
                    case HFEAT_SIMPLE_WEIGHT:
                        ((PoaFeatureSimpleWeight*)feature)->label = trueRefRleString->rleString[trueRefPos];
                        break;
                    case HFEAT_RLE_WEIGHT:
                    case HFEAT_NUCL_AND_RL_WEIGHT:
                        ((PoaFeatureRleWeight*)feature)->labelChar = trueRefRleString->rleString[trueRefPos];
                        trueRunLength = trueRefRleString->repeatCounts[trueRefPos];
                        if (trueRunLength > POAFEATURE_MAX_RUN_LENGTH) trueRunLength = POAFEATURE_MAX_RUN_LENGTH;
                        ((PoaFeatureRleWeight*)feature)->labelRunLength = trueRunLength;
                        break;
                    case HFEAT_SPLIT_RLE_WEIGHT:
                        rlFeature = ((PoaFeatureSplitRleWeight*)feature);
                        trueRunLength = trueRefRleString->repeatCounts[trueRefPos];
                        while (rlFeature != NULL) {
                            rlFeature->labelChar = trueRefRleString->rleString[trueRefPos];
                            if (trueRunLength <= 0) {
                                rlFeature->labelRunLength = 0;
                            } else if (trueRunLength > rlFeature->maxRunLength) {
                                rlFeature->labelRunLength = rlFeature->maxRunLength;
                            } else {
                                rlFeature->labelRunLength = trueRunLength;
                            }
                            trueRunLength -= rlFeature->maxRunLength;
                            rlFeature = rlFeature->nextRunLength;
                        }
                        break;
                    default:
                        st_errAbort("Unhandled FeatureType!\n");
                }
                trueRefPos++;
            }

            // delete
            else if (featureRefPos < stIntTuple_get(currRefAlign, FEATURE_POS)) {
                // apply label
                switch (featureType) {
                    case HFEAT_SIMPLE_WEIGHT:
                        ((PoaFeatureSimpleWeight*)feature)->label = '_';
                        break;
                    case HFEAT_RLE_WEIGHT:
                    case HFEAT_NUCL_AND_RL_WEIGHT:
                        ((PoaFeatureRleWeight*)feature)->labelChar = '_';
                        ((PoaFeatureRleWeight*)feature)->labelRunLength = 0;
                        break;
                    case HFEAT_SPLIT_RLE_WEIGHT:
                        rlFeature = ((PoaFeatureSplitRleWeight*)feature);
                        while (rlFeature != NULL) {
                            rlFeature->labelChar = '_';
                            rlFeature->labelRunLength = 0;
                            rlFeature = rlFeature->nextRunLength;
                        }
                        break;
                    default:
                        st_errAbort("Unhandled FeatureType!\n");
                }
            }

            // programmer error
            else {
                st_errAbort("Unhandled case annotating features with true reference characters!\n");
            }

            // always iterate over insert features
            prevFeature = feature;
            switch (featureType) {
                case HFEAT_SIMPLE_WEIGHT:
                    feature = ((PoaFeatureSimpleWeight*)feature)->nextInsert;
                    break;
                case HFEAT_RLE_WEIGHT:
                case HFEAT_NUCL_AND_RL_WEIGHT:
                    feature = ((PoaFeatureRleWeight*)feature)->nextInsert;
                    break;
                case HFEAT_SPLIT_RLE_WEIGHT:
                    feature = ((PoaFeatureSplitRleWeight*)feature)->nextInsert;
                    break;
                default:
                    st_errAbort("Unhandled FeatureType!\n");
            }
            featureInsPos++;
        }

        // this catches any true inserts which are not present in the poa / feature list
        while (currRefAlign != NULL &&
                featureRefPos < stIntTuple_get(currRefAlign, FEATURE_POS) &&
                trueRefPos < stIntTuple_get(currRefAlign, REFERENCE_POS)) {
            trueRefPos++;
        }
    }

    stList_destructIterator(trueRefAlignItor);
}

void poa_writeHelenFeatures(HelenFeatureType type, Poa *poa, stList *bamChunkReads, stList *rleStrings,
        char *outputFileBase, BamChunk *bamChunk, stList *trueRefAlignment, RleString *consensusRleString,
        RleString *trueRefRleString, bool fullFeatureOutput, int64_t maxRunLength) {
    // prep
    int64_t firstMatchedFeature = -1;
    int64_t lastMatchedFeature = -1;
    stList *features = NULL;
    bool outputLabels = trueRefAlignment != NULL && trueRefRleString != NULL;

    // handle differently based on type
    switch (type) {
        case HFEAT_SIMPLE_WEIGHT :
            // get features
            features = poa_getSimpleWeightFeatures(poa, bamChunkReads);
            firstMatchedFeature = 0;
            lastMatchedFeature = stList_length(features) - 1;

            // get truth (if we have it)
            if (outputLabels) {
                poa_annotateHelenFeaturesWithTruth(features, type, trueRefAlignment, trueRefRleString,
                                                   &firstMatchedFeature, &lastMatchedFeature);
            }

            // write it out
            if (fullFeatureOutput) {
                writeSimpleWeightHelenFeaturesTSV(outputFileBase, bamChunk, outputLabels, features,
                                                  firstMatchedFeature, lastMatchedFeature);
            }

            #ifdef _HDF5
            writeSimpleWeightHelenFeaturesHDF5(outputFileBase, bamChunk, outputLabels, features,
                                                            firstMatchedFeature, lastMatchedFeature);
            #endif

            break;

        case HFEAT_RLE_WEIGHT:
        case HFEAT_NUCL_AND_RL_WEIGHT:
            // get features
            features = poa_getRleWeightFeatures(poa, bamChunkReads, rleStrings, consensusRleString);
            firstMatchedFeature = 0;
            lastMatchedFeature = stList_length(features) - 1;

            // get truth (if we have it)
            if (outputLabels) {
                poa_annotateHelenFeaturesWithTruth(features, type, trueRefAlignment, trueRefRleString,
                                                   &firstMatchedFeature, &lastMatchedFeature);
            }

            // write it out
            if (fullFeatureOutput) {
                writeRleWeightHelenFeaturesTSV(outputFileBase, bamChunk, outputLabels, features,
                                               firstMatchedFeature, lastMatchedFeature);
            }

            #ifdef _HDF5
            if (type == HFEAT_RLE_WEIGHT) {
                writeRleWeightHelenFeaturesHDF5(outputFileBase, bamChunk, outputLabels, features,
                                                firstMatchedFeature, lastMatchedFeature);
            } else {
                writeNucleotideAndRleWeightHelenFeaturesHDF5(outputFileBase, bamChunk, outputLabels, features,
                                                             firstMatchedFeature, lastMatchedFeature);
            }
            #endif
            break;

        case HFEAT_SPLIT_RLE_WEIGHT:
            // get features
            features = poa_getSplitRleWeightFeatures(poa, bamChunkReads, rleStrings, maxRunLength);
            firstMatchedFeature = 0;
            lastMatchedFeature = stList_length(features) - 1;

            // get truth (if we have it)
            if (outputLabels) {
                poa_annotateHelenFeaturesWithTruth(features, type, trueRefAlignment, trueRefRleString,
                                                   &firstMatchedFeature, &lastMatchedFeature);
            }

            #ifdef _HDF5
            writeSplitRleWeightHelenFeaturesHDF5(outputFileBase, bamChunk, outputLabels, features,
                                                 firstMatchedFeature, lastMatchedFeature, maxRunLength);
            #endif
            break;
        default:
            st_errAbort("Unhandled HELEN feature type!\n");
    }

    //cleanup
    stList_destruct(features);
}


void writeSimpleWeightHelenFeaturesTSV(char *outputFileBase, BamChunk *bamChunk, bool outputLabels, stList *features,
                                       int64_t featureStartIdx, int64_t featureEndIdxInclusive) {

    char *outputFile = stString_print("%s.tsv", outputFileBase);
    FILE *fH = fopen(outputFile, "w");

    // print header
    fprintf(fH, "##contig:%s\n", bamChunk->refSeqName);
    fprintf(fH, "##contigStartPos:%"PRId64"\n", bamChunk->chunkBoundaryStart);
    fprintf(fH, "##contigEndPos:%"PRId64"\n", bamChunk->chunkBoundaryEnd);
    fprintf(fH, "#refPos\tinsPos");
    if (outputLabels) fprintf(fH, "\tlabel");
    for (int64_t i = 0; i < SYMBOL_NUMBER_NO_N; i++) {
        fprintf(fH, "\t%c_fwd\t%c_rev", symbol_convertSymbolToChar((Symbol)i), symbol_convertSymbolToChar((Symbol)i));
    }
    fprintf(fH, "\tgap_fwd\tgap_rev\n");

    // iterate over features
    for (int64_t i = featureStartIdx; i <= featureEndIdxInclusive; i++) {
        PoaFeatureSimpleWeight *feature = stList_get(features, i);

        // iterate over all inserts for each assembly position
        while (feature != NULL) {

            // position and label
            fprintf(fH, "%"PRId64, feature->refPosition);
            fprintf(fH, "\t%"PRId64, feature->insertPosition);
            if (outputLabels) {
                fprintf(fH, "\t%c", feature->label);
            }

            // print weights
            for (int64_t j = 0; j < SYMBOL_NUMBER_NO_N; j++) {
                fprintf(fH, "\t%7.4f", feature->weights[PoaFeature_SimpleWeight_charIndex((Symbol) j, TRUE)] / PAIR_ALIGNMENT_PROB_1);
                fprintf(fH, "\t%7.4f", feature->weights[PoaFeature_SimpleWeight_charIndex((Symbol) j, FALSE)] / PAIR_ALIGNMENT_PROB_1);
            }
            fprintf(fH, "\t%7.4f", feature->weights[PoaFeature_SimpleWeight_gapIndex(TRUE)] / PAIR_ALIGNMENT_PROB_1);
            fprintf(fH, "\t%7.4f\n", feature->weights[PoaFeature_SimpleWeight_gapIndex(FALSE)] / PAIR_ALIGNMENT_PROB_1);

            // iterate
            feature = feature->nextInsert;
        }
    }

    fclose(fH);
    free(outputFile);
}


void writeRleWeightHelenFeaturesTSV(char *outputFileBase, BamChunk *bamChunk, bool outputLabels, stList *features,
                                       int64_t featureStartIdx, int64_t featureEndIdxInclusive) {

    char *outputFile = stString_print("%s.tsv", outputFileBase);
    FILE *fH = fopen(outputFile, "w");

    // print header
    fprintf(fH, "##contig:%s\n", bamChunk->refSeqName);
    fprintf(fH, "##contigStartPos:%"PRId64"\n", bamChunk->chunkBoundaryStart);
    fprintf(fH, "##contigEndPos:%"PRId64"\n", bamChunk->chunkBoundaryEnd);
    fprintf(fH, "#refPos\tinsPos");
    if (outputLabels) {
        fprintf(fH, "\tlabelChar");
        fprintf(fH, "\tlabelRunLength");
    }
    for (int64_t symbol = 0; symbol < SYMBOL_NUMBER_NO_N; symbol++) {
        for (int64_t runLength = 1; runLength <= POAFEATURE_MAX_RUN_LENGTH; runLength++) {
            fprintf(fH, "\t%c_%"PRId64"_fwd\t%c_%"PRId64"_rev", symbol_convertSymbolToChar((Symbol)symbol), runLength,
                    symbol_convertSymbolToChar((Symbol)symbol), runLength);
        }
    }
    fprintf(fH, "\tgap_fwd\tgap_rev\n");

    // iterate over features
    for (int64_t i = featureStartIdx; i <= featureEndIdxInclusive; i++) {
        PoaFeatureRleWeight *feature = stList_get(features, i);

        // iterate over all inserts for each assembly position
        while (feature != NULL) {

            // position and label
            fprintf(fH, "%"PRId64, feature->refPosition);
            fprintf(fH, "\t%"PRId64, feature->insertPosition);
            if (outputLabels) {
                fprintf(fH, "\t%c", feature->labelChar);
                fprintf(fH, "\t%"PRId64, feature->labelRunLength);
            }

            // print weights
            for (int64_t symbol = 0; symbol < SYMBOL_NUMBER_NO_N; symbol++) {
                for (int64_t runLength = 1; runLength <= POAFEATURE_MAX_RUN_LENGTH; runLength++) {
                    fprintf(fH, "\t%7.4f", feature->weights[PoaFeature_RleWeight_charIndex((Symbol) symbol, runLength, TRUE)] / PAIR_ALIGNMENT_PROB_1);
                    fprintf(fH, "\t%7.4f", feature->weights[PoaFeature_RleWeight_charIndex((Symbol) symbol, runLength, FALSE)] / PAIR_ALIGNMENT_PROB_1);
                }
            }

            for (int64_t j = 0; j < SYMBOL_NUMBER_NO_N; j++) {
            }
            fprintf(fH, "\t%7.4f", feature->weights[PoaFeature_RleWeight_gapIndex(TRUE)] / PAIR_ALIGNMENT_PROB_1);
            fprintf(fH, "\t%7.4f\n", feature->weights[PoaFeature_RleWeight_gapIndex(FALSE)] / PAIR_ALIGNMENT_PROB_1);

            // iterate
            feature = feature->nextInsert;
        }
    }

    fclose(fH);
    free(outputFile);
}

#ifdef _HDF5

#define HDF5_FEATURE_SIZE 1000

double **getTwoDArrayDouble(int64_t rowCount, int64_t columnCount, bool zeroValues) {
    double **array  = st_calloc(rowCount, sizeof(double*));
    array[0] = (double*) st_calloc(columnCount * rowCount, sizeof(double) );
    if (zeroValues) {
        for (int64_t i = 0; i < columnCount * rowCount; i++) {
            (*array)[i] = 0.0;
        }
    }
    for (int64_t i=1; i < rowCount; i++) {
        array[i] = array[0] + i*columnCount;
    }
    return array;
}

float **getTwoDArrayFloat(int64_t rowCount, int64_t columnCount) {
    float **array  = st_calloc(rowCount, sizeof(float*));
    array[0] = (float*) st_calloc(columnCount * rowCount, sizeof(float) );
    for (int64_t i=1; i < rowCount; i++) {
        array[i] = array[0] + i*columnCount;
    }
    return array;
}

uint32_t **getTwoDArrayUInt32(int64_t rowCount, int64_t columnCount) {
    uint32_t **array  = st_calloc(rowCount, sizeof(uint32_t*));
    array[0] = (uint32_t*) st_calloc(columnCount * rowCount, sizeof(uint32_t) );
    for (int64_t i=1; i < rowCount; i++) {
        array[i] = array[0] + i*columnCount;
    }
    return array;
}

uint8_t **getTwoDArrayUInt8(int64_t rowCount, int64_t columnCount) {
    uint8_t **array  = st_calloc(rowCount, sizeof(uint8_t*));
    array[0] = (uint8_t*) st_calloc(columnCount * rowCount, sizeof(uint8_t) );
    for (int64_t i=1; i < rowCount; i++) {
        array[i] = array[0] + i*columnCount;
    }
    return array;
}

char **getTwoDArrayChar(int64_t rowCount, int64_t columnCount) {
    char **array  = st_calloc(rowCount, sizeof(char*));
    array[0] = (char*) st_calloc(columnCount * rowCount, sizeof(char) );
    for (int64_t i=1; i < rowCount; i++) {
        array[i] = array[0] + i*columnCount;
    }
    return array;
}

void writeSimpleWeightHelenFeaturesHDF5(char *outputFileBase, BamChunk *bamChunk, bool outputLabels, stList *features,
                                        int64_t featureStartIdx, int64_t featureEndIdxInclusive) {

    herr_t      status = 0;

    /*
     * Get feature data set up
     */

    // count features, create feature array
    uint64_t featureCount = 0;
    for (int64_t i = featureStartIdx; i <= featureEndIdxInclusive; i++) {
        PoaFeatureSimpleWeight *feature = stList_get(features, i);
        while (feature != NULL) {
            featureCount++;
            feature = feature->nextInsert;
        }
    }
    if (featureCount < HDF5_FEATURE_SIZE && outputLabels) {
        char *logIdentifier = getLogIdentifier();
        st_logInfo(" %s Feature count %"PRId64" less than minimum of %d\n", logIdentifier, featureCount, HDF5_FEATURE_SIZE);
        free(logIdentifier);
        return;
    }


    // get all feature data into an array
    uint32_t **positionData = getTwoDArrayUInt32(featureCount, 2);
    int64_t columnCount = SYMBOL_NUMBER * 2; //{A, C, T, G, Gap} x {fwd, rev}
    float **rleWeightData = getTwoDArrayFloat(featureCount, columnCount);
    char **labelCharacterData = NULL;
    if (outputLabels) {
        labelCharacterData = getTwoDArrayChar(featureCount, 1);
    }

    // add all data to features
    featureCount = 0;
    for (int64_t i = featureStartIdx; i <= featureEndIdxInclusive; i++) {
        PoaFeatureSimpleWeight *feature = stList_get(features, i);
        while (feature != NULL) {
            positionData[featureCount][0] = (uint32_t) feature->refPosition;
            positionData[featureCount][1] = (uint32_t) feature->insertPosition;

            for (int64_t symbol = 0; symbol < SYMBOL_NUMBER_NO_N; symbol++) {
                int64_t pos = PoaFeature_SimpleWeight_charIndex((Symbol) symbol, TRUE);
                rleWeightData[featureCount][pos] = (float) (feature->weights[pos] / PAIR_ALIGNMENT_PROB_1);
                pos = PoaFeature_SimpleWeight_charIndex((Symbol) symbol, FALSE);
                rleWeightData[featureCount][pos] = (float) (feature->weights[pos] / PAIR_ALIGNMENT_PROB_1);
            }

            // weights include 'N' index which is not included in features
            int64_t pos = PoaFeature_SimpleWeight_gapIndex(TRUE);
            rleWeightData[featureCount][pos - 2] = (float) (feature->weights[pos] / PAIR_ALIGNMENT_PROB_1);
            pos = PoaFeature_SimpleWeight_gapIndex(FALSE);
            rleWeightData[featureCount][pos - 2] = (float) (feature->weights[pos] / PAIR_ALIGNMENT_PROB_1);

            if (outputLabels) {
                labelCharacterData[featureCount][0] = feature->label;
            }

            featureCount++;
            feature = feature->nextInsert;
        }
    }

    /*
     * Get hdf5 data set up
     */

    hid_t int64Type = H5Tcopy(H5T_NATIVE_UINT32);
    status |= H5Tset_order(int64Type, H5T_ORDER_LE);
    hid_t uint32Type = H5Tcopy(H5T_NATIVE_UINT32);
    status |= H5Tset_order(uint32Type, H5T_ORDER_LE);
    hid_t floatType = H5Tcopy(H5T_NATIVE_FLOAT);
    status |= H5Tset_order(floatType, H5T_ORDER_LE);
    hid_t uint8Type = H5Tcopy(H5T_NATIVE_UINT8);
    status |= H5Tset_order(uint8Type, H5T_ORDER_LE);
    hid_t stringType = H5Tcopy (H5T_C_S1);
    status |= H5Tset_size (stringType, strlen(bamChunk->refSeqName) + 1);

    // so that we can produce chunks smaller than HDF5_FEATURE_SIZE (not used during training)
    hsize_t featureSize = (hsize_t) (featureCount < HDF5_FEATURE_SIZE ? featureCount : HDF5_FEATURE_SIZE);

    hsize_t metadataDimension[1] = {1};
    hsize_t postionDimension[2] = {featureSize, 2};
    hsize_t labelCharacterDimension[2] = {featureSize, 1};
    hsize_t rleWeightDimension[2] = {featureSize, (hsize_t) columnCount};

    hid_t metadataSpace = H5Screate_simple(1, metadataDimension, NULL);
    hid_t positionSpace = H5Screate_simple(2, postionDimension, NULL);
    hid_t labelCharacterSpace = H5Screate_simple(2, labelCharacterDimension, NULL);
    hid_t rleWeightSpace = H5Screate_simple(2, rleWeightDimension, NULL);

    /*
     * Write features to files
     */

    // each file must have exactly 1000 features
    int64_t totalFeatureFiles = (int64_t) (featureCount / HDF5_FEATURE_SIZE) + (featureCount % HDF5_FEATURE_SIZE == 0 ? 0 : 1);
    int64_t featureOffset = 0;
    if (featureCount >= HDF5_FEATURE_SIZE) {
        featureOffset = (int64_t) ((HDF5_FEATURE_SIZE * totalFeatureFiles - featureCount) / (int64_t) (featureCount / HDF5_FEATURE_SIZE));
    }

    for (int64_t featureIndex = 0; featureIndex < totalFeatureFiles; featureIndex++) {
        // get start pos
        int64_t chunkFeatureStartIdx = (HDF5_FEATURE_SIZE * featureIndex) - (featureOffset * featureIndex);
        if (featureIndex + 1 == totalFeatureFiles && featureCount >= HDF5_FEATURE_SIZE) {
            chunkFeatureStartIdx = featureCount - HDF5_FEATURE_SIZE;
        }

        // create file
        char *outputFile = stString_print("%s.%"PRId64".h5", outputFileBase, featureIndex);
        hid_t file = H5Fcreate (outputFile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        // write metadata
        hid_t contigDataset = H5Dcreate (file, "contig", stringType, metadataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (contigDataset, stringType, H5S_ALL, H5S_ALL, H5P_DEFAULT, bamChunk->refSeqName);
        hid_t contigStartDataset = H5Dcreate (file, "contig_start", int64Type, metadataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (contigStartDataset, int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bamChunk->chunkBoundaryStart);
        hid_t contigEndDataset = H5Dcreate (file, "contig_end", int64Type, metadataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (contigEndDataset, int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bamChunk->chunkBoundaryEnd);
        hid_t chunkIndexDataset = H5Dcreate (file, "feature_chunk_idx", int64Type, metadataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (chunkIndexDataset, int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &featureIndex);

        // write position info
        hid_t positionDataset = H5Dcreate (file, "position", uint32Type, positionSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (positionDataset, uint32Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, positionData[chunkFeatureStartIdx]);

        // write rle data
        hid_t rleWeightDataset = H5Dcreate (file, "image", floatType, rleWeightSpace, H5P_DEFAULT, H5P_DEFAULT,
                                            H5P_DEFAULT);
        status |= H5Dwrite (rleWeightDataset, floatType, H5S_ALL, H5S_ALL, H5P_DEFAULT, rleWeightData[chunkFeatureStartIdx]);

        // if labels, add all these too
        if (outputLabels) {
            hid_t labelCharacterDataset = H5Dcreate (file, "label_base", uint8Type, labelCharacterSpace, H5P_DEFAULT,
                                                     H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite (labelCharacterDataset, uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                                labelCharacterData[chunkFeatureStartIdx]);
            status |= H5Dclose (labelCharacterDataset);
        }

        // cleanup
        status |= H5Dclose (contigDataset);
        status |= H5Dclose (contigStartDataset);
        status |= H5Dclose (contigEndDataset);
        status |= H5Dclose (chunkIndexDataset);
        status |= H5Dclose (positionDataset);
        status |= H5Dclose (rleWeightDataset);
        status |= H5Fclose (file);
        free(outputFile);
    }

    // cleanup
    free(rleWeightData[0]);
    free(rleWeightData);
    free(positionData[0]);
    free(positionData);
    status |= H5Tclose (int64Type);
    status |= H5Tclose (uint32Type);
    status |= H5Tclose (floatType);
    status |= H5Tclose (uint8Type);
    status |= H5Tclose (stringType);
    status |= H5Sclose (metadataSpace);
    status |= H5Sclose (rleWeightSpace);
    status |= H5Sclose (positionSpace);
    status |= H5Sclose (labelCharacterSpace);
    if (outputLabels) {
        free(labelCharacterData[0]);
        free(labelCharacterData);
    }

    if (status) {
        char *logIdentifier = getLogIdentifier();
        st_logInfo(" %s Error writing HELEN features to HDF5 files: %s\n", logIdentifier, outputFileBase);
        free(logIdentifier);
    }
}

#define MAX_TOTAL_WEIGHT 64.0
uint8_t convertTotalWeightToUInt8(double totalWeight) {
    // convert to "depth space"
    totalWeight /= PAIR_ALIGNMENT_PROB_1;
    // cap at 64
    if (totalWeight > MAX_TOTAL_WEIGHT) totalWeight = MAX_TOTAL_WEIGHT;
    // convert to uint8
    return (uint8_t) (totalWeight / MAX_TOTAL_WEIGHT * (UINT8_MAX - 1));
}
uint8_t normalizeWeightToUInt8(double totalWeight, double weight) {
    return (uint8_t) (weight / totalWeight * (UINT8_MAX - 1));
}

void writeRleWeightHelenFeaturesHDF5(char *outputFileBase, BamChunk *bamChunk, bool outputLabels, stList *features,
                                     int64_t featureStartIdx, int64_t featureEndIdxInclusive) {

    herr_t      status = 0;

    /*
     * Get feature data set up
     */

    // count features, create feature array
    uint64_t featureCount = 0;
    for (int64_t i = featureStartIdx; i <= featureEndIdxInclusive; i++) {
        PoaFeatureRleWeight *feature = stList_get(features, i);
        while (feature != NULL) {
            featureCount++;
            feature = feature->nextInsert;
        }
    }
    if (featureCount < HDF5_FEATURE_SIZE && outputLabels) {
        char *logIdentifier = getLogIdentifier();
        st_logInfo(" %s Feature count %"PRId64" less than minimum of %d\n", logIdentifier, featureCount, HDF5_FEATURE_SIZE);
        free(logIdentifier);
        return;
    }

    // get all feature data into an array
    uint32_t **positionData = getTwoDArrayUInt32(featureCount, 2);
    uint8_t **predictedRunLengthData = getTwoDArrayUInt8(featureCount, 1);
    int64_t rleNucleotideColumnCount = POAFEATURE_RLE_WEIGHT_TOTAL_SIZE - POAFEATURE_MAX_RUN_LENGTH * 2; // don't output 'N' chars
    uint8_t **normalizationData = getTwoDArrayUInt8(featureCount, 1);
    uint8_t **normalizedRleNucleotideWeights = getTwoDArrayUInt8(featureCount, rleNucleotideColumnCount);
    char **labelCharacterData = NULL;
    uint8_t **labelRunLengthData = NULL;
    if (outputLabels) {
        labelCharacterData = getTwoDArrayChar(featureCount, 1);
        labelRunLengthData = getTwoDArrayUInt8(featureCount, 1);
    }

    // add all data to features
    featureCount = 0;
    for (int64_t i = featureStartIdx; i <= featureEndIdxInclusive; i++) {
        PoaFeatureRleWeight *feature = stList_get(features, i);
        while (feature != NULL) {
            positionData[featureCount][0] = (uint32_t) feature->refPosition;
            positionData[featureCount][1] = (uint32_t) feature->insertPosition;
            predictedRunLengthData[featureCount][0] = (uint8_t) feature->predictedRunLength;

            // get total weight
            double totalWeight = 0;
            int64_t pos;
            for (int64_t symbol = 0; symbol < SYMBOL_NUMBER_NO_N; symbol++) {
                for (int64_t runLength = 1; runLength <= POAFEATURE_MAX_RUN_LENGTH; runLength++) {
                    pos = PoaFeature_RleWeight_charIndex((Symbol) symbol, runLength, TRUE);
                    totalWeight += feature->weights[pos];
                    pos = PoaFeature_RleWeight_charIndex((Symbol) symbol, runLength, FALSE);
                    totalWeight += feature->weights[pos];
                }
            }
            pos = PoaFeature_RleWeight_gapIndex(TRUE);
            totalWeight += feature->weights[pos];
            pos = PoaFeature_RleWeight_gapIndex(FALSE);
            totalWeight += feature->weights[pos];
            normalizationData[featureCount][0] = convertTotalWeightToUInt8(totalWeight);

            for (int64_t symbol = 0; symbol < SYMBOL_NUMBER_NO_N; symbol++) {
                for (int64_t runLength = 1; runLength <= POAFEATURE_MAX_RUN_LENGTH; runLength++) {
                    pos = PoaFeature_RleWeight_charIndex((Symbol) symbol, runLength, TRUE);
                    normalizedRleNucleotideWeights[featureCount][pos] = normalizeWeightToUInt8(totalWeight, feature->weights[pos]);
                    pos = PoaFeature_RleWeight_charIndex((Symbol) symbol, runLength, FALSE);
                    normalizedRleNucleotideWeights[featureCount][pos] = normalizeWeightToUInt8(totalWeight, feature->weights[pos]);
                }
            }
            pos = PoaFeature_RleWeight_gapIndex(TRUE);
            normalizedRleNucleotideWeights[featureCount][pos - POAFEATURE_MAX_RUN_LENGTH * 2] =
                    normalizeWeightToUInt8(totalWeight, feature->weights[pos]);
            pos = PoaFeature_RleWeight_gapIndex(FALSE);
            normalizedRleNucleotideWeights[featureCount][pos - POAFEATURE_MAX_RUN_LENGTH * 2] =
                    normalizeWeightToUInt8(totalWeight, feature->weights[pos]);

            if (outputLabels) {
                labelCharacterData[featureCount][0] = feature->labelChar;
                labelRunLengthData[featureCount][0] = (uint8_t) feature->labelRunLength;
            }

            featureCount++;
            feature = feature->nextInsert;
        }
    }

    /*
     * Get hdf5 data set up
     */


    hid_t int64Type = H5Tcopy(H5T_NATIVE_UINT32);
    status |= H5Tset_order(int64Type, H5T_ORDER_LE);
    hid_t uint32Type = H5Tcopy(H5T_NATIVE_UINT32);
    status |= H5Tset_order(uint32Type, H5T_ORDER_LE);
    hid_t uint8Type = H5Tcopy(H5T_NATIVE_UINT8);
    status |= H5Tset_order(uint8Type, H5T_ORDER_LE);
    hid_t stringType = H5Tcopy (H5T_C_S1);
    status |= H5Tset_size (stringType, strlen(bamChunk->refSeqName) + 1);

    // so that we can produce chunks smaller than HDF5_FEATURE_SIZE (not used during training)
    hsize_t featureSize = (hsize_t) (featureCount < HDF5_FEATURE_SIZE ? featureCount : HDF5_FEATURE_SIZE);

    hsize_t metadataDimension[1] = {1};
    hsize_t postionDimension[2] = {featureSize, 2};
    hsize_t predictedRunLengthDimension[2] = {featureSize, 1};
    hsize_t labelCharacterDimension[2] = {featureSize, 1};
    hsize_t labelRunLengthDimension[2] = {featureSize, 1};
    hsize_t normalizationDimension[2] = {featureSize, 1};
    hsize_t rleNucleotideDimension[2] = {featureSize, (hsize_t) rleNucleotideColumnCount};

    hid_t metadataSpace = H5Screate_simple(1, metadataDimension, NULL);
    hid_t positionSpace = H5Screate_simple(2, postionDimension, NULL);
    hid_t predictedRunLengthSpace = H5Screate_simple(1, predictedRunLengthDimension, NULL);
    hid_t labelCharacterSpace = H5Screate_simple(2, labelCharacterDimension, NULL);
    hid_t labelRunLengthSpace = H5Screate_simple(2, labelRunLengthDimension, NULL);
    hid_t normalizationSpace = H5Screate_simple(2, normalizationDimension, NULL);
    hid_t rleNucleotideSpace = H5Screate_simple(2, rleNucleotideDimension, NULL);


    /*
     * Write features to files
     */

    // each file must have exactly 1000 features
    int64_t totalFeatureFiles = (int64_t) (featureCount / HDF5_FEATURE_SIZE) + (featureCount % HDF5_FEATURE_SIZE == 0 ? 0 : 1);
    int64_t featureOffset = 0;
    if (featureCount >= HDF5_FEATURE_SIZE) {
        featureOffset = (int64_t) ((HDF5_FEATURE_SIZE * totalFeatureFiles - featureCount) / (int64_t) (featureCount / HDF5_FEATURE_SIZE));
    }
    for (int64_t featureIndex = 0; featureIndex < totalFeatureFiles; featureIndex++) {
        // get start pos
        int64_t chunkFeatureStartIdx = (HDF5_FEATURE_SIZE * featureIndex) - (featureOffset * featureIndex);
        if (featureIndex + 1 == totalFeatureFiles && featureCount >= HDF5_FEATURE_SIZE) {
            chunkFeatureStartIdx = featureCount - HDF5_FEATURE_SIZE;
        }

        // create file
        char *outputFile = stString_print("%s.%"PRId64".h5", outputFileBase, featureIndex);
        hid_t file = H5Fcreate (outputFile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        // write metadata
        hid_t contigDataset = H5Dcreate (file, "contig", stringType, metadataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (contigDataset, stringType, H5S_ALL, H5S_ALL, H5P_DEFAULT, bamChunk->refSeqName);
        hid_t contigStartDataset = H5Dcreate (file, "contig_start", int64Type, metadataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (contigStartDataset, int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bamChunk->chunkBoundaryStart);
        hid_t contigEndDataset = H5Dcreate (file, "contig_end", int64Type, metadataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (contigEndDataset, int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bamChunk->chunkBoundaryEnd);
        hid_t chunkIndexDataset = H5Dcreate (file, "feature_chunk_idx", int64Type, metadataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (chunkIndexDataset, int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &featureIndex);

        // write position info
        hid_t positionDataset = H5Dcreate (file, "position", uint32Type, positionSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (positionDataset, uint32Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, positionData[chunkFeatureStartIdx]);

        // write rle data
        hid_t predictedRunLengthDataset = H5Dcreate (file, "bayesian_run_length_prediction", uint8Type,
                                                     predictedRunLengthSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (predictedRunLengthDataset, uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                            predictedRunLengthData[chunkFeatureStartIdx]);
        hid_t rleWeightDataset = H5Dcreate (file, "image", uint8Type, rleNucleotideSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (rleWeightDataset, uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                            normalizedRleNucleotideWeights[chunkFeatureStartIdx]);
        hid_t rleNormalizationDataset = H5Dcreate (file, "normalization", uint8Type, normalizationSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (rleNormalizationDataset, uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                            normalizationData[chunkFeatureStartIdx]);

        // if labels, add all these too
        if (outputLabels) {
            hid_t labelCharacterDataset = H5Dcreate (file, "label_base", uint8Type, labelCharacterSpace, H5P_DEFAULT,
                                                     H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite (labelCharacterDataset, uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                                labelCharacterData[chunkFeatureStartIdx]);
            hid_t labelRunLengthDataset = H5Dcreate (file, "label_run_length", uint8Type, labelRunLengthSpace,
                                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite (labelRunLengthDataset, uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                                labelRunLengthData[chunkFeatureStartIdx]);

            status |= H5Dclose (labelCharacterDataset);
            status |= H5Dclose (labelRunLengthDataset);
        }

        // cleanup
        status |= H5Dclose (contigDataset);
        status |= H5Dclose (contigStartDataset);
        status |= H5Dclose (contigEndDataset);
        status |= H5Dclose (chunkIndexDataset);
        status |= H5Dclose (positionDataset);
        status |= H5Dclose (predictedRunLengthDataset);
        status |= H5Dclose (rleWeightDataset);
        status |= H5Dclose (rleNormalizationDataset);
        status |= H5Fclose (file);
        free(outputFile);
    }

    // cleanup
    free(predictedRunLengthData[0]);
    free(predictedRunLengthData);
    free(normalizedRleNucleotideWeights[0]);
    free(normalizedRleNucleotideWeights);
    free(normalizationData[0]);
    free(normalizationData);
    free(positionData[0]);
    free(positionData);
    status |= H5Tclose (int64Type);
    status |= H5Tclose (uint32Type);
    status |= H5Tclose (uint8Type);
    status |= H5Tclose (stringType);
    status |= H5Sclose (metadataSpace);
    status |= H5Sclose (positionSpace);
    status |= H5Sclose (predictedRunLengthSpace);
    status |= H5Sclose (rleNucleotideSpace);
    status |= H5Sclose (normalizationSpace);
    status |= H5Sclose (labelRunLengthSpace);
    status |= H5Sclose (labelCharacterSpace);
    if (outputLabels) {
        free(labelCharacterData[0]);
        free(labelCharacterData);
        free(labelRunLengthData[0]);
        free(labelRunLengthData);
    }

    if (status) {
        char *logIdentifier = getLogIdentifier();
        st_logInfo(" %s Error writing HELEN features to HDF5 files: %s\n", logIdentifier, outputFileBase);
        free(logIdentifier);
    }
}


void writeNucleotideAndRleWeightHelenFeaturesHDF5(char *outputFileBase, BamChunk *bamChunk, bool outputLabels,
        stList *features, int64_t featureStartIdx, int64_t featureEndIdxInclusive) {

    herr_t      status = 0;

    /*
     * Get feature data set up
     */

    // count features, create feature array
    uint64_t featureCount = 0;
    for (int64_t i = featureStartIdx; i <= featureEndIdxInclusive; i++) {
        PoaFeatureRleWeight *feature = stList_get(features, i);
        while (feature != NULL) {
            featureCount++;
            feature = feature->nextInsert;
        }
    }
    if (featureCount < HDF5_FEATURE_SIZE && outputLabels) {
        char *logIdentifier = getLogIdentifier();
        st_logInfo(" %s Feature count %"PRId64" less than minimum of %d\n", logIdentifier, featureCount, HDF5_FEATURE_SIZE);
        free(logIdentifier);
        return;
    }

    // get all feature data into an array
    uint32_t **positionData = getTwoDArrayUInt32(featureCount, 2);
    uint8_t **predictedRunLengthData = getTwoDArrayUInt8(featureCount, 1);
    int64_t runLengthColumnCount = (POAFEATURE_MAX_RUN_LENGTH + 1) * 2; // gaps are rl 0, fwd/rev
    int64_t nucleotideColumnCount = (SYMBOL_NUMBER) * 2; // actg gap, fwd/rev
    uint8_t **normalizationData = getTwoDArrayUInt8(featureCount, 1);
    double **normalizedNucleotidePrep = getTwoDArrayDouble(featureCount, nucleotideColumnCount, TRUE);
    double **normalizedRunLengthPrep = getTwoDArrayDouble(featureCount, runLengthColumnCount, TRUE);
    uint8_t **normalizedNucleotideAndRunLengthData = getTwoDArrayUInt8(featureCount, nucleotideColumnCount + runLengthColumnCount);
    char **labelCharacterData = NULL;
    uint8_t **labelRunLengthData = NULL;
    if (outputLabels) {
        labelCharacterData = getTwoDArrayChar(featureCount, 1);
        labelRunLengthData = getTwoDArrayUInt8(featureCount, 1);
    }

    // add all data to features
    featureCount = 0;
    for (int64_t i = featureStartIdx; i <= featureEndIdxInclusive; i++) {
        PoaFeatureRleWeight *feature = stList_get(features, i);
        while (feature != NULL) {
            positionData[featureCount][0] = (uint32_t) feature->refPosition;
            positionData[featureCount][1] = (uint32_t) feature->insertPosition;
            predictedRunLengthData[featureCount][0] = (uint8_t) feature->predictedRunLength;

            // get total weight
            double totalWeight = 0;
            int64_t pos;

            for (int64_t symbol = 0; symbol < SYMBOL_NUMBER_NO_N; symbol++) {
                for (int64_t runLength = 1; runLength <= POAFEATURE_MAX_RUN_LENGTH; runLength++) {
                    //fwd
                    pos = PoaFeature_RleWeight_charIndex((Symbol) symbol, runLength, TRUE);
                    totalWeight += feature->weights[pos];
                    normalizedNucleotidePrep[featureCount][symbol * 2 + POS_STRAND_IDX] += feature->weights[pos];
                    normalizedRunLengthPrep[featureCount][runLength * 2 + POS_STRAND_IDX] += feature->weights[pos];
                    // rev
                    pos = PoaFeature_RleWeight_charIndex((Symbol) symbol, runLength, FALSE);
                    totalWeight += feature->weights[pos];
                    normalizedNucleotidePrep[featureCount][symbol * 2 + NEG_STRAND_IDX] += feature->weights[pos];
                    normalizedRunLengthPrep[featureCount][runLength * 2 + NEG_STRAND_IDX] += feature->weights[pos];
                }
            }
            // fwd gap
            pos = PoaFeature_RleWeight_gapIndex(TRUE);
            totalWeight += feature->weights[pos];
            normalizedNucleotidePrep[featureCount][SYMBOL_NUMBER_NO_N * 2 + POS_STRAND_IDX] += feature->weights[pos];
            normalizedRunLengthPrep[featureCount][POS_STRAND_IDX] += feature->weights[pos];
            // rev gap
            pos = PoaFeature_RleWeight_gapIndex(FALSE);
            totalWeight += feature->weights[pos];
            normalizedNucleotidePrep[featureCount][SYMBOL_NUMBER_NO_N * 2 + NEG_STRAND_IDX] += feature->weights[pos];
            normalizedRunLengthPrep[featureCount][NEG_STRAND_IDX] += feature->weights[pos];

            // save total weight
            normalizationData[featureCount][0] = convertTotalWeightToUInt8(totalWeight);

            // covert to total weight and save
            for (int64_t normNuclPos = 0; normNuclPos < nucleotideColumnCount; normNuclPos++) {
                normalizedNucleotideAndRunLengthData[featureCount][normNuclPos] =
                        normalizeWeightToUInt8(totalWeight, normalizedNucleotidePrep[featureCount][normNuclPos]);
            }
            for (int64_t normRunLenPos = 0; normRunLenPos < runLengthColumnCount; normRunLenPos++) {
                normalizedNucleotideAndRunLengthData[featureCount][nucleotideColumnCount + normRunLenPos] =
                        normalizeWeightToUInt8(totalWeight, normalizedRunLengthPrep[featureCount][normRunLenPos]);
            }

            // save labels if appropriate
            if (outputLabels) {
                labelCharacterData[featureCount][0] = feature->labelChar;
                labelRunLengthData[featureCount][0] = (uint8_t) feature->labelRunLength;
            }

            featureCount++;
            feature = feature->nextInsert;
        }
    }

    /*
     * Get hdf5 data set up
     */


    hid_t int64Type = H5Tcopy(H5T_NATIVE_UINT32);
    status |= H5Tset_order(int64Type, H5T_ORDER_LE);
    hid_t uint32Type = H5Tcopy(H5T_NATIVE_UINT32);
    status |= H5Tset_order(uint32Type, H5T_ORDER_LE);
    hid_t uint8Type = H5Tcopy(H5T_NATIVE_UINT8);
    status |= H5Tset_order(uint8Type, H5T_ORDER_LE);
    hid_t stringType = H5Tcopy (H5T_C_S1);
    status |= H5Tset_size (stringType, strlen(bamChunk->refSeqName) + 1);

    // so that we can produce chunks smaller than HDF5_FEATURE_SIZE (not used during training)
    hsize_t featureSize = (hsize_t) (featureCount < HDF5_FEATURE_SIZE ? featureCount : HDF5_FEATURE_SIZE);

    hsize_t metadataDimension[1] = {1};
    hsize_t postionDimension[2] = {featureSize, 2};
    hsize_t predictedRunLengthDimension[2] = {featureSize, 1};
    hsize_t labelCharacterDimension[2] = {featureSize, 1};
    hsize_t labelRunLengthDimension[2] = {featureSize, 1};
    hsize_t normalizationDimension[2] = {featureSize, 1};
    hsize_t nucleotideAndRunLengthDimension[2] = {featureSize, (hsize_t) (nucleotideColumnCount + runLengthColumnCount)};

    hid_t metadataSpace = H5Screate_simple(1, metadataDimension, NULL);
    hid_t positionSpace = H5Screate_simple(2, postionDimension, NULL);
    hid_t predictedRunLengthSpace = H5Screate_simple(1, predictedRunLengthDimension, NULL);
    hid_t labelCharacterSpace = H5Screate_simple(2, labelCharacterDimension, NULL);
    hid_t labelRunLengthSpace = H5Screate_simple(2, labelRunLengthDimension, NULL);
    hid_t normalizationSpace = H5Screate_simple(2, normalizationDimension, NULL);
    hid_t nucleotideAndRunLengthSpace = H5Screate_simple(2, nucleotideAndRunLengthDimension, NULL);


    /*
     * Write features to files
     */

    // each file must have exactly 1000 features
    int64_t totalFeatureFiles = (int64_t) (featureCount / HDF5_FEATURE_SIZE) + (featureCount % HDF5_FEATURE_SIZE == 0 ? 0 : 1);
    int64_t featureOffset = 0;
    if (featureCount >= HDF5_FEATURE_SIZE) {
        featureOffset = (int64_t) ((HDF5_FEATURE_SIZE * totalFeatureFiles - featureCount) / (int64_t) (featureCount / HDF5_FEATURE_SIZE));
    }
    for (int64_t featureIndex = 0; featureIndex < totalFeatureFiles; featureIndex++) {
        // get start pos
        int64_t chunkFeatureStartIdx = (HDF5_FEATURE_SIZE * featureIndex) - (featureOffset * featureIndex);
        if (featureIndex + 1 == totalFeatureFiles && featureCount >= HDF5_FEATURE_SIZE) {
            chunkFeatureStartIdx = featureCount - HDF5_FEATURE_SIZE;
        }

        // create file
        char *outputFile = stString_print("%s.%"PRId64".h5", outputFileBase, featureIndex);
        hid_t file = H5Fcreate (outputFile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        // write metadata
        hid_t contigDataset = H5Dcreate (file, "contig", stringType, metadataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (contigDataset, stringType, H5S_ALL, H5S_ALL, H5P_DEFAULT, bamChunk->refSeqName);
        hid_t contigStartDataset = H5Dcreate (file, "contig_start", int64Type, metadataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (contigStartDataset, int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bamChunk->chunkBoundaryStart);
        hid_t contigEndDataset = H5Dcreate (file, "contig_end", int64Type, metadataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (contigEndDataset, int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bamChunk->chunkBoundaryEnd);
        hid_t chunkIndexDataset = H5Dcreate (file, "feature_chunk_idx", int64Type, metadataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (chunkIndexDataset, int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &featureIndex);

        // write position info
        hid_t positionDataset = H5Dcreate (file, "position", uint32Type, positionSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (positionDataset, uint32Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, positionData[chunkFeatureStartIdx]);

        // write rle data
        hid_t predictedRunLengthDataset = H5Dcreate (file, "bayesian_run_length_prediction", uint8Type,
                                                     predictedRunLengthSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (predictedRunLengthDataset, uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                            predictedRunLengthData[chunkFeatureStartIdx]);
        hid_t nucleotideAndRunLengthDataset = H5Dcreate (file, "image", uint8Type, nucleotideAndRunLengthSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (nucleotideAndRunLengthDataset, uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                            normalizedNucleotideAndRunLengthData[chunkFeatureStartIdx]);
        hid_t normalizationDataset = H5Dcreate (file, "normalization", uint8Type, normalizationSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (normalizationDataset, uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                            normalizationData[chunkFeatureStartIdx]);

        // if labels, add all these too
        if (outputLabels) {
            hid_t labelCharacterDataset = H5Dcreate (file, "label_base", uint8Type, labelCharacterSpace, H5P_DEFAULT,
                                                     H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite (labelCharacterDataset, uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                                labelCharacterData[chunkFeatureStartIdx]);
            hid_t labelRunLengthDataset = H5Dcreate (file, "label_run_length", uint8Type, labelRunLengthSpace,
                                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite (labelRunLengthDataset, uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                                labelRunLengthData[chunkFeatureStartIdx]);

            status |= H5Dclose (labelCharacterDataset);
            status |= H5Dclose (labelRunLengthDataset);
        }

        // cleanup
        status |= H5Dclose (contigDataset);
        status |= H5Dclose (contigStartDataset);
        status |= H5Dclose (contigEndDataset);
        status |= H5Dclose (chunkIndexDataset);
        status |= H5Dclose (positionDataset);
        status |= H5Dclose (predictedRunLengthDataset);
        status |= H5Dclose (nucleotideAndRunLengthDataset);
        status |= H5Dclose (normalizationDataset);
        status |= H5Fclose (file);
        free(outputFile);
    }

    // cleanup
    free(predictedRunLengthData[0]);
    free(predictedRunLengthData);
    free(normalizedNucleotidePrep[0]);
    free(normalizedNucleotidePrep);
    free(normalizedNucleotideAndRunLengthData[0]);
    free(normalizedRunLengthPrep[0]);
    free(normalizedRunLengthPrep);
    free(normalizationData[0]);
    free(normalizationData);
    free(positionData[0]);
    free(positionData);
    status |= H5Tclose (int64Type);
    status |= H5Tclose (uint32Type);
    status |= H5Tclose (uint8Type);
    status |= H5Tclose (stringType);
    status |= H5Sclose (metadataSpace);
    status |= H5Sclose (positionSpace);
    status |= H5Sclose (predictedRunLengthSpace);
    status |= H5Sclose (nucleotideAndRunLengthSpace);
    status |= H5Sclose (normalizationSpace);
    status |= H5Sclose (labelRunLengthSpace);
    status |= H5Sclose (labelCharacterSpace);
    if (outputLabels) {
        free(labelCharacterData[0]);
        free(labelCharacterData);
        free(labelRunLengthData[0]);
        free(labelRunLengthData);
    }

    if (status) {
        char *logIdentifier = getLogIdentifier();
        st_logInfo(" %s Error writing HELEN features to HDF5 files: %s\n", logIdentifier, outputFileBase);
        free(logIdentifier);
    }
}

void writeSplitRleWeightHelenFeaturesHDF5(char *outputFileBase, BamChunk *bamChunk, bool outputLabels, stList *features,
                                          int64_t featureStartIdx, int64_t featureEndIdxInclusive, const int64_t maxRunLength) {

    herr_t      status = 0;

    /*
     * Get feature data set up
     */

    // count features, create feature array
    uint64_t featureCount = 0;
    for (int64_t i = featureStartIdx; i <= featureEndIdxInclusive; i++) {
        PoaFeatureSplitRleWeight *feature = stList_get(features, i);
        while (feature != NULL) {
            PoaFeatureSplitRleWeight *rlFeature = feature;
            while (rlFeature != NULL) {
                featureCount++;
                rlFeature = rlFeature->nextRunLength;
            }
            feature = feature->nextInsert;
        }
    }
    if (featureCount < HDF5_FEATURE_SIZE && outputLabels) {
        char *logIdentifier = getLogIdentifier();
        st_logInfo(" %s Feature count %"PRId64" less than minimum of %d\n", logIdentifier, featureCount, HDF5_FEATURE_SIZE);
        free(logIdentifier);
        return;
    }

    // get all feature data into an array
    uint32_t **positionData = getTwoDArrayUInt32(featureCount, 3);
    int64_t rleNucleotideColumnCount = ((SYMBOL_NUMBER - 1) * (maxRunLength + 1) + 1) * 2;
    uint8_t **normalizationData = getTwoDArrayUInt8(featureCount, 1);
    uint8_t **imageData = getTwoDArrayUInt8(featureCount, rleNucleotideColumnCount);
    uint8_t **labelCharacterData = NULL;
    uint8_t **labelRunLengthData = NULL;
    if (outputLabels) {
        labelCharacterData = getTwoDArrayUInt8(featureCount, 1);
        labelRunLengthData = getTwoDArrayUInt8(featureCount, 1);
    }

    // add all data to features
    featureCount = 0;
    for (int64_t i = featureStartIdx; i <= featureEndIdxInclusive; i++) {
        PoaFeatureSplitRleWeight *refFeature = stList_get(features, i);

        // total weight is calculated for the very first refPos feature, used for all inserts and run lengths
        double totalWeight = 0;
        for (int64_t j = 0; j < rleNucleotideColumnCount; j++) {
            totalWeight += refFeature->weights[j];
        }

        // iterate over all insert features
        PoaFeatureSplitRleWeight *insFeature = refFeature;
        while (insFeature != NULL) {

            // iterate over all run length features
            PoaFeatureSplitRleWeight *rlFeature = insFeature;
            while (rlFeature != NULL) {
                // position
                positionData[featureCount][0] = (uint32_t) rlFeature->refPosition;
                positionData[featureCount][1] = (uint32_t) rlFeature->insertPosition;
                positionData[featureCount][2] = (uint32_t) rlFeature->runLengthPosition;

                // normalization
                normalizationData[featureCount][0] = convertTotalWeightToUInt8(totalWeight);

                // copy weights over (into normalized uint8 space)
                for (int64_t j = 0; j < rleNucleotideColumnCount; j++) {
                    imageData[featureCount][j] =
                            normalizeWeightToUInt8(totalWeight, rlFeature->weights[j]);
                }

                // labels
                if (outputLabels) {
                    Symbol label = symbol_convertCharToSymbol(rlFeature->labelChar);
                    labelCharacterData[featureCount][0] = (uint8_t) (label == n ? 0 : label + 1);
                    labelRunLengthData[featureCount][0] = (uint8_t) (label == n ? 0 : rlFeature->labelRunLength);
                    if (labelRunLengthData[featureCount][0] > maxRunLength) {
                        st_errAbort("Encountered run length of %d (max %"PRId64") in chunk %s:%"PRId64"-%"PRId64,
                                    labelRunLengthData[featureCount][0], maxRunLength, bamChunk->refSeqName,
                                    bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
                    }
                }

                // iterate
                featureCount++;
                rlFeature = rlFeature->nextRunLength;
            }
            insFeature = insFeature->nextInsert;
        }
    }

    /*
     * Get hdf5 data set up
     */


    hid_t int64Type = H5Tcopy(H5T_NATIVE_UINT32);
    status |= H5Tset_order(int64Type, H5T_ORDER_LE);
    hid_t uint32Type = H5Tcopy(H5T_NATIVE_UINT32);
    status |= H5Tset_order(uint32Type, H5T_ORDER_LE);
    hid_t uint8Type = H5Tcopy(H5T_NATIVE_UINT8);
    status |= H5Tset_order(uint8Type, H5T_ORDER_LE);
    hid_t stringType = H5Tcopy (H5T_C_S1);
    status |= H5Tset_size (stringType, strlen(bamChunk->refSeqName) + 1);

    // so that we can produce chunks smaller than HDF5_FEATURE_SIZE (not used during training)
    hsize_t featureSize = (hsize_t) (featureCount < HDF5_FEATURE_SIZE ? featureCount : HDF5_FEATURE_SIZE);

    hsize_t metadataDimension[1] = {1};
    hsize_t postionDimension[2] = {featureSize, 3};
    hsize_t labelCharacterDimension[2] = {featureSize, 1};
    hsize_t labelRunLengthDimension[2] = {featureSize, 1};
    hsize_t normalizationDimension[2] = {featureSize, 1};
    hsize_t imageDimension[2] = {featureSize, (hsize_t) rleNucleotideColumnCount};

    hid_t metadataSpace = H5Screate_simple(1, metadataDimension, NULL);
    hid_t positionSpace = H5Screate_simple(2, postionDimension, NULL);
    hid_t labelCharacterSpace = H5Screate_simple(2, labelCharacterDimension, NULL);
    hid_t labelRunLengthSpace = H5Screate_simple(2, labelRunLengthDimension, NULL);
    hid_t normalizationSpace = H5Screate_simple(2, normalizationDimension, NULL);
    hid_t imageSpace = H5Screate_simple(2, imageDimension, NULL);


    /*
     * Write features to files
     */

    // each file must have exactly 1000 features
    int64_t totalFeatureFiles = (int64_t) (featureCount / HDF5_FEATURE_SIZE) + (featureCount % HDF5_FEATURE_SIZE == 0 ? 0 : 1);
    int64_t featureOffset = 0;
    if (featureCount >= HDF5_FEATURE_SIZE) {
        featureOffset = (int64_t) ((HDF5_FEATURE_SIZE * totalFeatureFiles - featureCount) / (int64_t) (featureCount / HDF5_FEATURE_SIZE));
    }
    for (int64_t featureIndex = 0; featureIndex < totalFeatureFiles; featureIndex++) {
        // get start pos
        int64_t chunkFeatureStartIdx = (HDF5_FEATURE_SIZE * featureIndex) - (featureOffset * featureIndex);
        if (featureIndex + 1 == totalFeatureFiles && featureCount >= HDF5_FEATURE_SIZE) {
            chunkFeatureStartIdx = featureCount - HDF5_FEATURE_SIZE;
        }

        // create file
        char *outputFile = stString_print("%s.%"PRId64".h5", outputFileBase, featureIndex);
        hid_t file = H5Fcreate (outputFile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        // write metadata
        hid_t contigDataset = H5Dcreate (file, "contig", stringType, metadataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (contigDataset, stringType, H5S_ALL, H5S_ALL, H5P_DEFAULT, bamChunk->refSeqName);
        hid_t contigStartDataset = H5Dcreate (file, "contig_start", int64Type, metadataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (contigStartDataset, int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bamChunk->chunkBoundaryStart);
        hid_t contigEndDataset = H5Dcreate (file, "contig_end", int64Type, metadataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (contigEndDataset, int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bamChunk->chunkBoundaryEnd);
        hid_t chunkIndexDataset = H5Dcreate (file, "feature_chunk_idx", int64Type, metadataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (chunkIndexDataset, int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &featureIndex);

        // write position info
        hid_t positionDataset = H5Dcreate (file, "position", uint32Type, positionSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (positionDataset, uint32Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, positionData[chunkFeatureStartIdx]);

        // write rle data
        hid_t imageDataset = H5Dcreate (file, "image", uint8Type, imageSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (imageDataset, uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                            imageData[chunkFeatureStartIdx]);
        hid_t normalizationDataset = H5Dcreate (file, "normalization", uint8Type, normalizationSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (normalizationDataset, uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                            normalizationData[chunkFeatureStartIdx]);

        // if labels, add all these too
        if (outputLabels) {
            hid_t labelCharacterDataset = H5Dcreate (file, "label_base", uint8Type, labelCharacterSpace, H5P_DEFAULT,
                                                     H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite (labelCharacterDataset, uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                                labelCharacterData[chunkFeatureStartIdx]);
            hid_t labelRunLengthDataset = H5Dcreate (file, "label_run_length", uint8Type, labelRunLengthSpace,
                                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite (labelRunLengthDataset, uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                                labelRunLengthData[chunkFeatureStartIdx]);

            status |= H5Dclose (labelCharacterDataset);
            status |= H5Dclose (labelRunLengthDataset);
        }

        // cleanup
        status |= H5Dclose (contigDataset);
        status |= H5Dclose (contigStartDataset);
        status |= H5Dclose (contigEndDataset);
        status |= H5Dclose (chunkIndexDataset);
        status |= H5Dclose (positionDataset);
        status |= H5Dclose (imageDataset);
        status |= H5Dclose (normalizationDataset);
        status |= H5Fclose (file);
        free(outputFile);
    }

    // cleanup
    free(imageData[0]);
    free(imageData);
    free(normalizationData[0]);
    free(normalizationData);
    free(positionData[0]);
    free(positionData);
    status |= H5Tclose (int64Type);
    status |= H5Tclose (uint32Type);
    status |= H5Tclose (uint8Type);
    status |= H5Tclose (stringType);
    status |= H5Sclose (metadataSpace);
    status |= H5Sclose (positionSpace);
    status |= H5Sclose (imageSpace);
    status |= H5Sclose (normalizationSpace);
    status |= H5Sclose (labelRunLengthSpace);
    status |= H5Sclose (labelCharacterSpace);
    if (outputLabels) {
        free(labelCharacterData[0]);
        free(labelCharacterData);
        free(labelRunLengthData[0]);
        free(labelRunLengthData);
    }

    if (status) {
        char *logIdentifier = getLogIdentifier();
        st_logInfo(" %s Error writing HELEN features to HDF5 files: %s\n", logIdentifier, outputFileBase);
        free(logIdentifier);
    }
}

#endif