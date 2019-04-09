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
    feature->nextInsert = NULL;
    return feature;
}

void PoaFeature_RleWeight_destruct(PoaFeatureRleWeight *feature) {
    if (feature->nextInsert != NULL) {
        PoaFeature_RleWeight_destruct(feature->nextInsert);
    }
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

                    //TODO += delete->weightXXXStrand / delete->length?
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
                    //TODO += insert->weightXXXStrand / strlen(insert->insert)?
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

stList *poa_getWeightedRleFeatures(Poa *poa, stList *bamChunkReads, stList *rleStrings) {

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

                    //TODO += delete->weightXXXStrand / delete->length?
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
                        //TODO += insert->weightXXXStrand / strlen(insert->insert)?
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
                        ((PoaFeatureRleWeight*)feature)->labelChar = '_';
                        ((PoaFeatureRleWeight*)feature)->labelRunLength = 1; //TODO 0?
                        feature = ((PoaFeatureRleWeight*)feature)->nextInsert;
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
                        ((PoaFeatureRleWeight*)feature)->labelChar = trueRefRleString->rleString[trueRefPos];
                        int64_t trueRunLength = trueRefRleString->repeatCounts[trueRefPos];
                        if (trueRunLength > POAFEATURE_MAX_RUN_LENGTH) trueRunLength = POAFEATURE_MAX_RUN_LENGTH;
                        ((PoaFeatureRleWeight*)feature)->labelRunLength = trueRunLength;
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
                        ((PoaFeatureRleWeight*)feature)->labelChar = trueRefRleString->rleString[trueRefPos];
                        ((PoaFeatureRleWeight*)feature)->labelRunLength = trueRefRleString->repeatCounts[trueRefPos];
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
                        ((PoaFeatureRleWeight*)feature)->labelChar = '_';
                        ((PoaFeatureRleWeight*)feature)->labelRunLength = 1; //TODO 0?
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
                    feature = ((PoaFeatureRleWeight*)feature)->nextInsert;
                    break;
                default:
                    st_errAbort("Unhandled FeatureType!\n");
            }
            featureInsPos++;
        }

        // this catches any true inserts which are not present in the poa / feature list
        while (currRefAlign != NULL && featureRefPos < stIntTuple_get(currRefAlign, FEATURE_POS) && trueRefPos < stIntTuple_get(currRefAlign, REFERENCE_POS)) {
            // DO NOT make new empty feature and save truth
            // TODO remove this once you're sure
            //PoaFeatureSimpleWeight *newFeature = PoaFeature_SimpleWeight_construct(featureRefPos, featureInsPos);
            //newFeature->label = trueRefRleString->rleString[trueRefPos];

            // DO NOT save and DO iterate
            //prevFeature->nextInsert = newFeature;
            //prevFeature = newFeature;
            //featureInsPos++;
            trueRefPos++;
        }
    }

    stList_destructIterator(trueRefAlignItor);
}

void poa_writeHelenFeatures(HelenFeatureType type, Poa *poa, stList *bamChunkReads, stList *rleStrings,
        char *outputFileBase, BamChunk *bamChunk, stList *trueRefAlignment, RleString *trueRefRleString,
        bool fullFeatureOutput) {
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
            // get features
            features = poa_getWeightedRleFeatures(poa, bamChunkReads, rleStrings);
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
            writeRleWeightHelenFeaturesHDF5(outputFileBase, bamChunk, outputLabels, features,
                                                            firstMatchedFeature, lastMatchedFeature);
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

double **getTwoDArrayDouble(int64_t rowCount, int64_t columnCount) {
    double **array  = malloc(rowCount*sizeof(double*));
    array[0] = (double*) malloc(columnCount * rowCount * sizeof(double) );
    for (int64_t i=1; i < rowCount; i++) {
        array[i] = array[0] + i*columnCount;
    }
    return array;
}

int64_t **getTwoDArrayInt64(int64_t rowCount, int64_t columnCount) {
    int64_t **array  = malloc(rowCount*sizeof(int64_t*));
    array[0] = (int64_t*) malloc(columnCount * rowCount * sizeof(int64_t) );
    for (int64_t i=1; i < rowCount; i++) {
        array[i] = array[0] + i*columnCount;
    }
    return array;
}

char **getTwoDArrayChar(int64_t rowCount, int64_t columnCount) {
    char **array  = malloc(rowCount*sizeof(char*));
    array[0] = (char*) malloc(columnCount * rowCount * sizeof(char) );
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
    if (featureCount < HDF5_FEATURE_SIZE) {
        char *logIdentifier = getLogIdentifier();
        st_logInfo(" %s Feature count %"PRId64" less than minimum of %d\n", logIdentifier, featureCount, HDF5_FEATURE_SIZE);
        free(logIdentifier);
        return;
    }


    // get all feature data into an array
    int64_t **positionData = getTwoDArrayInt64(featureCount, 2);
    int64_t rleColumnCount = SYMBOL_NUMBER * 2; //{A, C, T, G, Gap} x {fwd, rev}
    double **rleWeightData = getTwoDArrayDouble(featureCount, rleColumnCount);
    char **labelCharacterData = NULL;
    if (outputLabels) {
        labelCharacterData = getTwoDArrayChar(featureCount, 1);
    }

    // add all data to features
    featureCount = 0;
    for (int64_t i = featureStartIdx; i <= featureEndIdxInclusive; i++) {
        PoaFeatureSimpleWeight *feature = stList_get(features, i);
        while (feature != NULL) {
            positionData[featureCount][0] = feature->refPosition;
            positionData[featureCount][1] = feature->insertPosition;

            for (int64_t symbol = 0; symbol < SYMBOL_NUMBER_NO_N; symbol++) {
                int64_t pos = PoaFeature_SimpleWeight_charIndex((Symbol) symbol, TRUE);
                rleWeightData[featureCount][pos] = feature->weights[pos] / PAIR_ALIGNMENT_PROB_1;
                pos = PoaFeature_SimpleWeight_charIndex((Symbol) symbol, FALSE);
                rleWeightData[featureCount][pos] = feature->weights[pos] / PAIR_ALIGNMENT_PROB_1;
            }

            // weights include 'N' index which is not included in features
            int64_t pos = PoaFeature_SimpleWeight_gapIndex(TRUE);
            rleWeightData[featureCount][pos - 2] = feature->weights[pos] / PAIR_ALIGNMENT_PROB_1;
            pos = PoaFeature_SimpleWeight_gapIndex(FALSE);
            rleWeightData[featureCount][pos - 2] = feature->weights[pos] / PAIR_ALIGNMENT_PROB_1;

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

    hid_t int64Type = H5Tcopy(H5T_NATIVE_INT64);
    status |= H5Tset_order(int64Type, H5T_ORDER_LE);
    hid_t doubleType = H5Tcopy(H5T_NATIVE_DOUBLE);
    status |= H5Tset_order(doubleType, H5T_ORDER_LE);
    hid_t charType = H5Tcopy(H5T_NATIVE_CHAR);
    status |= H5Tset_order(charType, H5T_ORDER_LE);
    hid_t stringType = H5Tcopy (H5T_C_S1);
    status |= H5Tset_size (stringType, strlen(bamChunk->refSeqName) + 1);

    hsize_t metadataDimension[1] = {1};
    hsize_t postionDimension[2] = {(hsize_t)HDF5_FEATURE_SIZE, 2};
    hsize_t labelCharacterDimension[2] = {(hsize_t)HDF5_FEATURE_SIZE, 1};
    hsize_t rleWeightDimension[2] = {(hsize_t)HDF5_FEATURE_SIZE, (hsize_t)rleColumnCount};

    hid_t metadataSpace = H5Screate_simple(1, metadataDimension, NULL);
    hid_t positionSpace = H5Screate_simple(2, postionDimension, NULL);
    hid_t labelCharacterSpace = H5Screate_simple(2, labelCharacterDimension, NULL);
    hid_t rleWeightSpace = H5Screate_simple(2, rleWeightDimension, NULL);

    /*
     * Write features to files
     */

    // each file must have exactly 1000 features
    int64_t totalFeatureFiles = (int64_t) (featureCount / HDF5_FEATURE_SIZE) + (featureCount % HDF5_FEATURE_SIZE == 0 ? 0 : 1);
    int64_t featureOffset = (int64_t) ((HDF5_FEATURE_SIZE * totalFeatureFiles - featureCount) / (int64_t) (featureCount / HDF5_FEATURE_SIZE));

    for (int64_t featureIndex = 0; featureIndex < totalFeatureFiles; featureIndex++) {
        // get start pos
        int64_t chunkFeatureStartIdx = (HDF5_FEATURE_SIZE * featureIndex) - (featureOffset * featureIndex);
        if (featureIndex + 1 == totalFeatureFiles) chunkFeatureStartIdx = featureCount - HDF5_FEATURE_SIZE;

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
        hid_t positionDataset = H5Dcreate (file, "position", int64Type, positionSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (positionDataset, int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, positionData[chunkFeatureStartIdx]);

        // write rle data
        hid_t rleWeightDataset = H5Dcreate (file, "image", doubleType, rleWeightSpace, H5P_DEFAULT, H5P_DEFAULT,
                                            H5P_DEFAULT);
        status |= H5Dwrite (rleWeightDataset, doubleType, H5S_ALL, H5S_ALL, H5P_DEFAULT, rleWeightData[chunkFeatureStartIdx]);

        // if labels, add all these too
        if (outputLabels) {
            hid_t labelCharacterDataset = H5Dcreate (file, "label_base", charType, labelCharacterSpace, H5P_DEFAULT,
                                                     H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite (labelCharacterDataset, charType, H5S_ALL, H5S_ALL, H5P_DEFAULT,
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
    status |= H5Tclose (doubleType);
    status |= H5Tclose (charType);
    status |= H5Tclose (stringType);
    status |= H5Sclose (metadataSpace);
    status |= H5Sclose (rleWeightSpace);
    status |= H5Sclose (positionSpace);
    status |= H5Sclose (labelCharacterSpace);
    if (outputLabels) {
        free(labelCharacterData[0]);
        free(labelCharacterData);
    }

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
    if (featureCount < HDF5_FEATURE_SIZE) {
        char *logIdentifier = getLogIdentifier();
        st_logInfo(" %s Feature count %"PRId64" less than minimum of %d\n", logIdentifier, featureCount, HDF5_FEATURE_SIZE);
        free(logIdentifier);
        return;
    }

    // get all feature data into an array
    int64_t **positionData = getTwoDArrayInt64(featureCount, 2);
    int64_t rleColumnCount = POAFEATURE_RLE_WEIGHT_TOTAL_SIZE - POAFEATURE_MAX_RUN_LENGTH * 2; // don't output 'N' chars
    double **rleWeightData = getTwoDArrayDouble(featureCount, rleColumnCount);
    char **labelCharacterData = NULL;
    int64_t **labelRunLengthData = NULL;
    if (outputLabels) {
        labelCharacterData = getTwoDArrayChar(featureCount, 1);
        labelRunLengthData = getTwoDArrayInt64(featureCount, 1);
    }

    // add all data to features
    featureCount = 0;
    for (int64_t i = featureStartIdx; i <= featureEndIdxInclusive; i++) {
        PoaFeatureRleWeight *feature = stList_get(features, i);
        while (feature != NULL) {
            positionData[featureCount][0] = feature->refPosition;
            positionData[featureCount][1] = feature->insertPosition;

            for (int64_t symbol = 0; symbol < SYMBOL_NUMBER_NO_N; symbol++) {
                for (int64_t runLength = 1; runLength <= POAFEATURE_MAX_RUN_LENGTH; runLength++) {
                    int64_t pos = PoaFeature_RleWeight_charIndex((Symbol) symbol, runLength, TRUE);
                    rleWeightData[featureCount][pos] = feature->weights[pos] / PAIR_ALIGNMENT_PROB_1;
                    pos = PoaFeature_RleWeight_charIndex((Symbol) symbol, runLength, FALSE);
                    rleWeightData[featureCount][pos] = feature->weights[pos] / PAIR_ALIGNMENT_PROB_1;
                }
            }
            int64_t pos = PoaFeature_RleWeight_gapIndex(TRUE);
            rleWeightData[featureCount][pos - POAFEATURE_MAX_RUN_LENGTH * 2] = feature->weights[pos] / PAIR_ALIGNMENT_PROB_1;
            pos = PoaFeature_RleWeight_gapIndex(FALSE);
            rleWeightData[featureCount][pos - POAFEATURE_MAX_RUN_LENGTH * 2] = feature->weights[pos] / PAIR_ALIGNMENT_PROB_1;

            if (outputLabels) {
                labelCharacterData[featureCount][0] = feature->labelChar;
                labelRunLengthData[featureCount][0] = feature->labelRunLength;
            }

            featureCount++;
            feature = feature->nextInsert;
        }
    }

    /*
     * Get hdf5 data set up
     */


    hid_t int64Type = H5Tcopy(H5T_NATIVE_INT64);
    status |= H5Tset_order(int64Type, H5T_ORDER_LE);
    hid_t doubleType = H5Tcopy(H5T_NATIVE_DOUBLE);
    status |= H5Tset_order(doubleType, H5T_ORDER_LE);
    hid_t charType = H5Tcopy(H5T_NATIVE_CHAR);
    status |= H5Tset_order(charType, H5T_ORDER_LE);
    hid_t stringType = H5Tcopy (H5T_C_S1);
    status |= H5Tset_size (stringType, strlen(bamChunk->refSeqName) + 1);

    hsize_t metadataDimension[1] = {1};
    hsize_t postionDimension[2] = {(hsize_t)HDF5_FEATURE_SIZE, 2};
    hsize_t labelCharacterDimension[2] = {(hsize_t)HDF5_FEATURE_SIZE, 1};
    hsize_t labelRunLengthDimension[2] = {(hsize_t)HDF5_FEATURE_SIZE, 1};
    hsize_t rleDimension[2] = {(hsize_t)HDF5_FEATURE_SIZE, (hsize_t)rleColumnCount};

    hid_t metadataSpace = H5Screate_simple(1, metadataDimension, NULL);
    hid_t positionSpace = H5Screate_simple(2, postionDimension, NULL);
    hid_t labelCharacterSpace = H5Screate_simple(2, labelCharacterDimension, NULL);
    hid_t labelRunLengthSpace = H5Screate_simple(2, labelRunLengthDimension, NULL);
    hid_t rleWeightSpace = H5Screate_simple(2, rleDimension, NULL);


    /*
     * Write features to files
     */

    // each file must have exactly 1000 features
    int64_t totalFeatureFiles = (int64_t) (featureCount / HDF5_FEATURE_SIZE) + (featureCount % HDF5_FEATURE_SIZE == 0 ? 0 : 1);
    int64_t featureOffset = (int64_t) ((HDF5_FEATURE_SIZE * totalFeatureFiles - featureCount) / (int64_t) (featureCount / HDF5_FEATURE_SIZE));

    for (int64_t featureIndex = 0; featureIndex < totalFeatureFiles; featureIndex++) {
        // get start pos
        int64_t chunkFeatureStartIdx = (HDF5_FEATURE_SIZE * featureIndex) - (featureOffset * featureIndex);
        if (featureIndex + 1 == totalFeatureFiles) chunkFeatureStartIdx = featureCount - HDF5_FEATURE_SIZE;

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
        hid_t positionDataset = H5Dcreate (file, "position", int64Type, positionSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (positionDataset, int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, positionData[chunkFeatureStartIdx]);

        // write rle data
        hid_t rleWeightDataset = H5Dcreate (file, "image", doubleType, rleWeightSpace, H5P_DEFAULT, H5P_DEFAULT,
                                               H5P_DEFAULT);
        status |= H5Dwrite (rleWeightDataset, doubleType, H5S_ALL, H5S_ALL, H5P_DEFAULT, rleWeightData[chunkFeatureStartIdx]);

        // if labels, add all these too
        if (outputLabels) {
            hid_t labelCharacterDataset = H5Dcreate (file, "label_base", charType, labelCharacterSpace, H5P_DEFAULT,
                    H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite (labelCharacterDataset, charType, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    labelCharacterData[chunkFeatureStartIdx]);
            hid_t labelRunLengthDataset = H5Dcreate (file, "label_run_length", int64Type, labelCharacterSpace,
                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite (labelRunLengthDataset, int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
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
    status |= H5Tclose (doubleType);
    status |= H5Tclose (charType);
    status |= H5Tclose (stringType);
    status |= H5Sclose (metadataSpace);
    status |= H5Sclose (positionSpace);
    status |= H5Sclose (rleWeightSpace);
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