//
// Created by tpesout on 3/29/19.
//

#include "margin.h"
#include "helenFeatures.h"

#ifdef _HDF5
#include <hdf5.h>
#endif


PoaFeatureSimpleCharacterCount *PoaFeature_SimpleCharacterCount_construct(int64_t refPos, int64_t insPos) {
    PoaFeatureSimpleCharacterCount *scc = st_calloc(1, sizeof(PoaFeatureSimpleCharacterCount));
    scc->refPosition = refPos;
    scc->insertPosition = insPos;
    scc->label = '\0';
    scc->nextInsert = NULL;
    return scc;
}
void PoaFeature_SimpleCharacterCount_destruct(PoaFeatureSimpleCharacterCount *scc) {
    if (scc->nextInsert != NULL) {
        PoaFeature_SimpleCharacterCount_destruct(scc->nextInsert);
    }
    free(scc);
}

int64_t PoaFeature_SimpleCharacterCount_getTotalCount(PoaFeatureSimpleCharacterCount *scc) {
    int64_t total = 0;
    for (int64_t i = 0; i < POAFEATURE_TOTAL_SIZE; i++) {
        total += scc->counts[i];
    }
    return total == 0 ? 1 : total;
}

double PoaFeature_SimpleCharacterCount_getTotalWeight(PoaFeatureSimpleCharacterCount *scc) {
    double total = 0;
    for (int64_t i = 0; i < POAFEATURE_TOTAL_SIZE; i++) {
        total += scc->weights[i];
    }
    return total == 0 ? 1.0 : total;
}

stList *poa_getSimpleCharacterCountFeatures(Poa *poa, stList *bamChunkReads) {

    // initialize feature list
    stList *featureList = stList_construct3(0, (void (*)(void *)) PoaFeature_SimpleCharacterCount_destruct);
    for(int64_t i=1; i<stList_length(poa->nodes); i++) {
        stList_append(featureList, PoaFeature_SimpleCharacterCount_construct(i-1, 0));
    }

    // for logging (of errors)
    char *logIdentifier = getLogIdentifier();

    // iterate over all positions
    for(int64_t i=0; i<stList_length(featureList); i++) {

        // get feature and node
        PoaFeatureSimpleCharacterCount* feature = stList_get(featureList, i);
        PoaNode *node = stList_get(poa->nodes, i + 1); //skip the first poa node, as it's always an 'N', so featureIdx and poaIdx are off by one

        // only include nucleotide counts once for each read
        stHash* readToMaxMatchWeight = stHash_construct();
        for (int64_t j = 0; j < stList_length(node->observations); j++) {
            PoaBaseObservation *obs = stList_get(node->observations, j);

            stHash_insert(readToMaxMatchWeight, (void *) obs->readNo, (void *)
                    ((int64_t) obs->weight > (int64_t) stHash_search(readToMaxMatchWeight, (void *) obs->readNo) ?
                     (int64_t) obs->weight : (int64_t) stHash_search(readToMaxMatchWeight, (void *) obs->readNo)));
        }

        // iterate over all observations in node (to get character counts)
        for(int64_t j=0; j<stList_length(node->observations); j++) {
            PoaBaseObservation *baseObs = stList_get(node->observations, j);

            // only get max weight observations for each readId
            if ((int64_t) baseObs->weight != (int64_t) stHash_search(readToMaxMatchWeight, (void *) baseObs->readNo)) {
                continue;
            } else {
                // for rare case of equal weight
                stHash_insert(readToMaxMatchWeight, (void*) baseObs->readNo, 0);
            }

            // count strand and nucleotide
            BamChunkRead *read = stList_get(bamChunkReads, baseObs->readNo);
            char base = read->nucleotides[baseObs->offset];
            feature->counts[symbol_convertCharToSymbol(base) * 2 + (read->forwardStrand ? POS_STRAND_IDX : NEG_STRAND_IDX)] += 1;
        }

        // get weights for all bases, strands
        double totalWeight, totalPositiveWeight, totalNegativeWeight;
        double *baseWeights = poaNode_getStrandSpecificBaseWeights(node, bamChunkReads,
                                                                   &totalWeight, &totalPositiveWeight, &totalNegativeWeight);
        for(int64_t j=0; j<SYMBOL_NUMBER; j++) {
            double positiveStrandBaseWeight = baseWeights[j * 2 + POS_STRAND_IDX];
            double negativeStrandBaseWeight = baseWeights[j * 2 + NEG_STRAND_IDX];
            feature->weights[j * 2 + POS_STRAND_IDX] += positiveStrandBaseWeight;
            feature->weights[j * 2 + NEG_STRAND_IDX] += negativeStrandBaseWeight;
        }
        free(baseWeights);

        // Deletes
        if (stList_length(node->deletes) > 0) {

            // only count each read once for deletes
            stHash* readToMaxDelWeight = stHash_construct();
            for (int64_t j = 0; j < stList_length(node->deletes); j++) {
                PoaInsert *delete = stList_get(node->deletes, j);
                for (int64_t k = 0; k < stList_length(delete->observations); k++) {
                    PoaBaseObservation *obs = stList_get(delete->observations, k);

                    stHash_insert(readToMaxDelWeight, (void *) obs->readNo, (void *)
                            ((int64_t) obs->weight > (int64_t) stHash_search(readToMaxDelWeight, (void *) obs->readNo) ?
                             (int64_t) obs->weight : (int64_t) stHash_search(readToMaxDelWeight, (void *) obs->readNo)));
                }
            }

            // iterate over all deletes
            for (int64_t j = 0; j < stList_length(node->deletes); j++) {
                PoaDelete *delete = stList_get(node->deletes, j);

                // get alignment counts
                int64_t posObvsCount = 0;
                int64_t negObvsCount = 0;
                for (int64_t k = 0; k < stList_length(delete->observations); k++) {
                    PoaBaseObservation *obs = stList_get(delete->observations, k);
                    // only get max weight observations for each readId
                    if ((int64_t) obs->weight != (int64_t) stHash_search(readToMaxDelWeight, (void *) obs->readNo)) {
                        continue;
                    } else {
                        // for rare case of equal weight
                        stHash_insert(readToMaxDelWeight, (void*) obs->readNo, 0);
                    }

                    if (((BamChunkRead *) stList_get(bamChunkReads, obs->readNo))->forwardStrand) {
                        posObvsCount++;
                    } else {
                        negObvsCount++;
                    }
                }

                // Deletes start AFTER the current position, need to add counts/weights to nodes after the current node
                for (int64_t k = 1; k < delete->length; k++) {
                    if (i + k >= stList_length(poa->nodes)) {
                        st_logInfo(" %s Encountered DELETE that occurs after end of POA!\n", logIdentifier);
                        break;
                    }
                    PoaFeatureSimpleCharacterCount *delFeature = stList_get(featureList, i + k);
                    delFeature->weights[POAFEATURE_SYMBOL_GAP_POS * 2 + POS_STRAND_IDX] += delete->weightForwardStrand;
                    delFeature->weights[POAFEATURE_SYMBOL_GAP_POS * 2 + NEG_STRAND_IDX] += delete->weightReverseStrand;
                    delFeature->counts[POAFEATURE_SYMBOL_GAP_POS * 2 + POS_STRAND_IDX] += posObvsCount;
                    delFeature->counts[POAFEATURE_SYMBOL_GAP_POS * 2 + NEG_STRAND_IDX] += negObvsCount;
                }
            }

            // cleanup
            stHash_destruct(readToMaxDelWeight);
        }

        // Inserts
        if (stList_length(node->inserts) > 0) {

            // only include nucleotide counts once for each read
            stHash* readToMaxInsWeight = stHash_construct();
            for (int64_t j = 0; j < stList_length(node->inserts); j++) {
                PoaInsert *insert = stList_get(node->inserts, j);
                for (int64_t k = 0; k < stList_length(insert->observations); k++) {
                    PoaBaseObservation *obs = stList_get(insert->observations, k);

                    stHash_insert(readToMaxInsWeight, (void *) obs->readNo, (void *)
                            ((int64_t) obs->weight > (int64_t) stHash_search(readToMaxInsWeight, (void *) obs->readNo) ?
                             (int64_t) obs->weight : (int64_t) stHash_search(readToMaxInsWeight, (void *) obs->readNo)));
                }
            }

            // iterate over all inserts
            for (int64_t j = 0; j < stList_length(node->inserts); j++) {
                PoaInsert *insert = stList_get(node->inserts, j);

                // get counts for observations
                int64_t fwdObservations = 0;
                int64_t revObservations = 0;
                for (int64_t k = 0; k < stList_length(insert->observations); k++) {
                    PoaBaseObservation *obs = stList_get(insert->observations, k);
                    // only get max weight observations for each readId
                    if ((int64_t) obs->weight != (int64_t) stHash_search(readToMaxInsWeight, (void *) obs->readNo)) {
                        continue;
                    } else {
                        // for rare case of equal weight
                        stHash_insert(readToMaxInsWeight, (void*) obs->readNo, 0);
                    }

                    BamChunkRead *bcr = stList_get(bamChunkReads, obs->readNo);
                    if (bcr->forwardStrand) {
                        fwdObservations++;
                    } else {
                        revObservations++;
                    }

                    // sanity check
                    //char *obsStr = stString_getSubString(bcr->nucleotides, obs->offset, strlen(insert->insert));
                }

                // get feature iterator
                PoaFeatureSimpleCharacterCount *prevFeature = feature;
                for (int64_t k = 0; k < strlen(insert->insert); k++) {
                    // get current feature (or create if necessary)
                    PoaFeatureSimpleCharacterCount *currFeature = prevFeature->nextInsert;
                    if (currFeature == NULL) {
                        currFeature = PoaFeature_SimpleCharacterCount_construct(i, k + 1);
                        prevFeature->nextInsert = currFeature;
                    }

                    Symbol c = symbol_convertCharToSymbol(insert->insert[k]);

                    // add weights
                    currFeature->weights[c * 2 + POS_STRAND_IDX] += insert->weightForwardStrand;
                    currFeature->weights[c * 2 + NEG_STRAND_IDX] += insert->weightReverseStrand;

                    // add counts
                    currFeature->counts[c * 2 + POS_STRAND_IDX] += fwdObservations;
                    currFeature->counts[c * 2 + NEG_STRAND_IDX] += revObservations;

                    // iterate
                    prevFeature = currFeature;
                }
            }

            // cleanup
            stHash_destruct(readToMaxInsWeight);
        }
        stHash_destruct(readToMaxMatchWeight);
    }

    free(logIdentifier);
    return featureList;
}

void poa_annotateSimpleCharacterCountFeaturesWithTruth(stList *features, stList *trueRefAlignment, RleString *trueRefRleString) {
    /*
     Each index in features represents a character in the final consensus string (which in turn was a single node in the
     final POA which was used to generate the consensus).  This consensus was aligned against the true reference (which
     is reflected in the trueRefRleString).  Items in the true refAlignment are stIntTuples with index 0 being the
     alignment weight (discarded), index 1 being the position in the consensus string (and features), and index 2 being
     the position in the trueRefRleString.  So we can iterate over features and the true alignment to assign truth
     labels to each feature.
     */
    int FPOS = 1;
    int RPOS = 2;
    bool TP_DEBUG = FALSE;
    if (TP_DEBUG) fprintf(stderr, "\nAnnotating features with truth:\n");

    // iterate over true ref alignment
    stListIterator *trueRefAlignItor = stList_getIterator(trueRefAlignment);
    stIntTuple *currRefAlign = stList_getNext(trueRefAlignItor);

    // iterate over features
    int64_t trueRefPos = stIntTuple_get(currRefAlign, RPOS);
    for (int64_t featureRefPos = 0; featureRefPos < stList_length(features); featureRefPos++) {
        PoaFeatureSimpleCharacterCount *feature = stList_get(features, featureRefPos);
        PoaFeatureSimpleCharacterCount *prevFeature = NULL;
        if (TP_DEBUG) fprintf(stderr, "\n");

        int64_t featureInsPos = 0;
        while (feature != NULL) {
            // no more ref bases, everything is gaps
            if (currRefAlign == NULL) {
                if (TP_DEBUG) fprintf(stderr, "(\t    )\t%"PRId64",%"PRId64"\t%"PRId64"\t-> DEL\n", featureRefPos,
                        featureInsPos, trueRefPos);
                feature->label = '_';
                feature = feature->nextInsert;
                continue;
            }

            // debug and sanity checks
            if (TP_DEBUG) fprintf(stderr, "(%4"PRId64"\t%4"PRId64")\t%4"PRId64",%"PRId64"\t%4"PRId64,
                    (int64_t ) stIntTuple_get(currRefAlign, FPOS), (int64_t ) stIntTuple_get(currRefAlign, RPOS),
                    featureRefPos, featureInsPos, trueRefPos);
            assert(stIntTuple_get(currRefAlign, FPOS) >= featureRefPos && stIntTuple_get(currRefAlign, RPOS) >= trueRefPos);

            // match
            if (stIntTuple_get(currRefAlign, FPOS) == featureRefPos && stIntTuple_get(currRefAlign, RPOS) == trueRefPos) {
                if (TP_DEBUG) fprintf(stderr, "\t -> MATCH\n");
                feature->label = trueRefRleString->rleString[trueRefPos];
                trueRefPos++;
                currRefAlign = stList_getNext(trueRefAlignItor);
            }
                // insert
            else if (trueRefPos < stIntTuple_get(currRefAlign, RPOS)) {
                if (TP_DEBUG) fprintf(stderr, "\t -> INS\n");
                feature->label = trueRefRleString->rleString[trueRefPos];
                trueRefPos++;
            }
                // delete
            else if (featureRefPos < stIntTuple_get(currRefAlign, FPOS)) {
                if (TP_DEBUG) fprintf(stderr, "\t -> DEL\n");
                feature->label = '_';
            }
                // programmer error
            else {
                st_errAbort("Unhandled case annotating features with true reference characters!\n");
            }

            // always iterate over insert features
            prevFeature = feature;
            feature = feature->nextInsert;
            featureInsPos++;
        }

        // this catches any true inserts which are not present in the poa / feature list
        while (currRefAlign != NULL && featureRefPos < stIntTuple_get(currRefAlign, FPOS) && trueRefPos < stIntTuple_get(currRefAlign, RPOS)) {
            if (TP_DEBUG) fprintf(stderr, "(%4"PRId64"\t%4"PRId64")\t%4"PRId64",%"PRId64"\t%4"PRId64,
                    (int64_t ) stIntTuple_get(currRefAlign, FPOS), (int64_t ) stIntTuple_get(currRefAlign, RPOS),
                    featureRefPos, featureInsPos, trueRefPos);
            if (TP_DEBUG) fprintf(stderr, "\t -> INS (new)\n");
            // make new empty feature and save truth
            PoaFeatureSimpleCharacterCount *newFeature = PoaFeature_SimpleCharacterCount_construct(featureRefPos, featureInsPos);
            newFeature->label = trueRefRleString->rleString[trueRefPos];

            // save and iterate
            prevFeature->nextInsert = newFeature;
            prevFeature = newFeature;
            trueRefPos++;
            featureInsPos++;
        }
    }

    if (TP_DEBUG) fprintf(stderr, "Finished annotating with truth!\n\n");
    stList_destructIterator(trueRefAlignItor);
}

void poa_writeHelenFeatures(HelenFeatureType type, Poa *poa, stList *bamChunkReads, char *outputFileBase,
                            BamChunk *bamChunk, stList *trueRefAlignment, RleString *trueRefRleString) {
    // prep
    stList *features = NULL;
    bool outputLabels = trueRefAlignment != NULL && trueRefRleString != NULL;

    // handle differently based on type
    switch (type) {
        case HFEAT_SIMPLE_COUNT :
        case HFEAT_SIMPLE_WEIGHT :
            // get features
            features = poa_getSimpleCharacterCountFeatures(poa, bamChunkReads);
            if (outputLabels) {
                poa_annotateSimpleCharacterCountFeaturesWithTruth(features, trueRefAlignment, trueRefRleString);
            }

            char *outputFile = stString_print("%s.tsv", outputFileBase);
            writeSimpleHelenFeaturesTSV(outputFile, bamChunk, outputLabels, features, type);
            free(outputFile);

            #ifdef _HDF5
            outputFile = stString_print("%s.h5", outputFileBase);
            int status = writeSimpleHelenFeaturesHDF5(outputFile, bamChunk, outputLabels, features, type);
            if (status) {
                st_logInfo(" Error writing HELEN features to %s\n", outputFile);
            }
            free(outputFile);
            #endif

            break;
        default:
            st_errAbort("Unhandled HELEN feature type!\n");
    }

    //cleanup
    stList_destruct(features);
}


void writeSimpleHelenFeaturesTSV(char *outputFile, BamChunk *bamChunk, bool outputLabels, stList *features,
                                 HelenFeatureType type) {

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
    for (int64_t i = 0; i < stList_length(features); i++) {
        PoaFeatureSimpleCharacterCount *feature = stList_get(features, i);

        // iterate over all inserts for each assembly position
        while (feature != NULL) {

            // position and label
            fprintf(fH, "%"PRId64, feature->refPosition);
            fprintf(fH, "\t%"PRId64, feature->insertPosition);
            if (outputLabels) {
                fprintf(fH, "\t%c", feature->label);
            }

            // print counts or weights
            if (type == HFEAT_SIMPLE_COUNT) {
                for (int64_t j = 0; j < SYMBOL_NUMBER_NO_N; j++) {
                    fprintf(fH, "\t%"PRId64, feature->counts[j * 2 + POS_STRAND_IDX]);
                    fprintf(fH, "\t%"PRId64, feature->counts[j * 2 + NEG_STRAND_IDX]);
                }
                fprintf(fH, "\t%"PRId64, feature->counts[POAFEATURE_SYMBOL_GAP_POS * 2 + POS_STRAND_IDX]);
                fprintf(fH, "\t%"PRId64"\n", feature->counts[POAFEATURE_SYMBOL_GAP_POS * 2 + NEG_STRAND_IDX]);
            } else {
                for (int64_t j = 0; j < SYMBOL_NUMBER_NO_N; j++) {
                    fprintf(fH, "\t%7.4f", feature->weights[j * 2 + POS_STRAND_IDX] / PAIR_ALIGNMENT_PROB_1);
                    fprintf(fH, "\t%7.4f", feature->weights[j * 2 + NEG_STRAND_IDX] / PAIR_ALIGNMENT_PROB_1);
                }
                fprintf(fH, "\t%7.4f", feature->weights[POAFEATURE_SYMBOL_GAP_POS * 2 + POS_STRAND_IDX] / PAIR_ALIGNMENT_PROB_1);
                fprintf(fH, "\t%7.4f\n", feature->weights[POAFEATURE_SYMBOL_GAP_POS * 2 + NEG_STRAND_IDX] / PAIR_ALIGNMENT_PROB_1);
            }

            // iterate
            feature = feature->nextInsert;
        }
    }

    fclose(fH);
}

#ifdef _HDF5
typedef struct {
    int64_t     refPos;
    int64_t     insPos;
} helen_features_position_hdf5_record_t;

typedef struct {
    char        label;
} helen_features_label_hdf5_record_t;

typedef struct {
    double      wAFwd;
    double      wARev;
    double      wCFwd;
    double      wCRev;
    double      wGFwd;
    double      wGRev;
    double      wTFwd;
    double      wTRev;
    double      wGapFwd;
    double      wGapRev;
} helen_features_simple_weight_hdf5_record_t;      /* Compound type */

int writeSimpleHelenFeaturesHDF5(char *outputFile, BamChunk *bamChunk, bool outputLabels, stList *features, HelenFeatureType type) {

    herr_t      status = 0;

    // count features, create feature array
    uint64_t featureCount = 0;
    for (int64_t i = 0; i < stList_length(features); i++) {
        PoaFeatureSimpleCharacterCount *feature = stList_get(features, i);
        while (feature != NULL) {
            featureCount++;
            feature = feature->nextInsert;
        }
    }
    helen_features_position_hdf5_record_t *positionData =
            st_calloc(featureCount, sizeof(helen_features_simple_weight_hdf5_record_t));
    helen_features_simple_weight_hdf5_record_t *simpleWeightData =
            st_calloc(featureCount, sizeof(helen_features_simple_weight_hdf5_record_t));
    helen_features_label_hdf5_record_t *labelData = NULL;
    if (outputLabels) {
        labelData = st_calloc(featureCount, sizeof(helen_features_simple_weight_hdf5_record_t));
    }

    // add all data to features
    featureCount = 0;
    for (int64_t i = 0; i < stList_length(features); i++) {
        PoaFeatureSimpleCharacterCount *feature = stList_get(features, i);
        while (feature != NULL) {
            positionData[featureCount].refPos = feature->refPosition;
            positionData[featureCount].insPos = feature->insertPosition;
            simpleWeightData[featureCount].wAFwd = feature->weights[symbol_convertCharToSymbol('A') * 2 + POS_STRAND_IDX] / PAIR_ALIGNMENT_PROB_1;
            simpleWeightData[featureCount].wARev = feature->weights[symbol_convertCharToSymbol('A') * 2 + NEG_STRAND_IDX] / PAIR_ALIGNMENT_PROB_1;
            simpleWeightData[featureCount].wCFwd = feature->weights[symbol_convertCharToSymbol('C') * 2 + POS_STRAND_IDX] / PAIR_ALIGNMENT_PROB_1;
            simpleWeightData[featureCount].wCRev = feature->weights[symbol_convertCharToSymbol('C') * 2 + NEG_STRAND_IDX] / PAIR_ALIGNMENT_PROB_1;
            simpleWeightData[featureCount].wGFwd = feature->weights[symbol_convertCharToSymbol('G') * 2 + POS_STRAND_IDX] / PAIR_ALIGNMENT_PROB_1;
            simpleWeightData[featureCount].wGRev = feature->weights[symbol_convertCharToSymbol('G') * 2 + NEG_STRAND_IDX] / PAIR_ALIGNMENT_PROB_1;
            simpleWeightData[featureCount].wTFwd = feature->weights[symbol_convertCharToSymbol('T') * 2 + POS_STRAND_IDX] / PAIR_ALIGNMENT_PROB_1;
            simpleWeightData[featureCount].wTRev = feature->weights[symbol_convertCharToSymbol('T') * 2 + NEG_STRAND_IDX] / PAIR_ALIGNMENT_PROB_1;
            simpleWeightData[featureCount].wGapFwd = feature->weights[POAFEATURE_SYMBOL_GAP_POS * 2 + POS_STRAND_IDX] / PAIR_ALIGNMENT_PROB_1;
            simpleWeightData[featureCount].wGapRev = feature->weights[POAFEATURE_SYMBOL_GAP_POS * 2 + NEG_STRAND_IDX] / PAIR_ALIGNMENT_PROB_1;
            if (outputLabels) {
                labelData[featureCount].label = feature->label;
            }

            featureCount++;
            feature = feature->nextInsert;
        }
    }

    // create file
    hid_t file = H5Fcreate (outputFile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // create datatypes
    hid_t simpleWeightType = H5Tcreate (H5T_COMPOUND, sizeof (helen_features_simple_weight_hdf5_record_t));
    status |= H5Tinsert (simpleWeightType, "a_fwd", HOFFSET (helen_features_simple_weight_hdf5_record_t, wAFwd), H5T_NATIVE_DOUBLE);
    status |= H5Tinsert (simpleWeightType, "a_rev", HOFFSET (helen_features_simple_weight_hdf5_record_t, wARev), H5T_NATIVE_DOUBLE);
    status |= H5Tinsert (simpleWeightType, "c_fwd", HOFFSET (helen_features_simple_weight_hdf5_record_t, wCFwd), H5T_NATIVE_DOUBLE);
    status |= H5Tinsert (simpleWeightType, "c_rev", HOFFSET (helen_features_simple_weight_hdf5_record_t, wCRev), H5T_NATIVE_DOUBLE);
    status |= H5Tinsert (simpleWeightType, "g_fwd", HOFFSET (helen_features_simple_weight_hdf5_record_t, wGFwd), H5T_NATIVE_DOUBLE);
    status |= H5Tinsert (simpleWeightType, "g_rev", HOFFSET (helen_features_simple_weight_hdf5_record_t, wGRev), H5T_NATIVE_DOUBLE);
    status |= H5Tinsert (simpleWeightType, "t_fwd", HOFFSET (helen_features_simple_weight_hdf5_record_t, wTFwd), H5T_NATIVE_DOUBLE);
    status |= H5Tinsert (simpleWeightType, "t_rev", HOFFSET (helen_features_simple_weight_hdf5_record_t, wTRev), H5T_NATIVE_DOUBLE);
    status |= H5Tinsert (simpleWeightType, "gap_fwd", HOFFSET (helen_features_simple_weight_hdf5_record_t, wGapFwd), H5T_NATIVE_DOUBLE);
    status |= H5Tinsert (simpleWeightType, "gap_rev", HOFFSET (helen_features_simple_weight_hdf5_record_t, wGapRev), H5T_NATIVE_DOUBLE);

    hid_t positionType = H5Tcreate (H5T_COMPOUND, sizeof (helen_features_position_hdf5_record_t));
    status |= H5Tinsert (positionType, "refPos", HOFFSET (helen_features_position_hdf5_record_t, refPos), H5T_NATIVE_INT64);
    status |= H5Tinsert (positionType, "insPos", HOFFSET (helen_features_position_hdf5_record_t, insPos), H5T_NATIVE_INT64);

    hid_t labelType = H5Tcreate (H5T_COMPOUND, sizeof (helen_features_label_hdf5_record_t));
    status |= H5Tinsert (labelType, "label", HOFFSET (helen_features_label_hdf5_record_t, label), H5T_NATIVE_CHAR);

    // size of dataset
    hsize_t dimension = featureCount;
    hid_t space = H5Screate_simple (1, &dimension, NULL);

    // create and write datasets
    hid_t positionDataset = H5Dcreate (file, "position", positionType, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t simpleWeightDataset = H5Dcreate (file, "simpleWeight", simpleWeightType, space, H5P_DEFAULT, H5P_DEFAULT,
                                           H5P_DEFAULT);
    status |= H5Dwrite (simpleWeightDataset, simpleWeightType, H5S_ALL, H5S_ALL, H5P_DEFAULT, simpleWeightData);
    status |= H5Dwrite (positionDataset, positionType, H5S_ALL, H5S_ALL, H5P_DEFAULT, positionData);

    // if labels, add all these too
    if (outputLabels) {
        hid_t labelDataset = H5Dcreate (file, "label", labelType, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (labelDataset, labelType, H5S_ALL, H5S_ALL, H5P_DEFAULT, labelData);
        status |= H5Dclose (labelDataset);
    }

    // cleanup
    status |= H5Dclose (positionDataset);
    status |= H5Dclose (simpleWeightDataset);
    status |= H5Sclose (space);
    status |= H5Tclose (positionType);
    status |= H5Tclose (simpleWeightType);
    status |= H5Tclose (labelType);
    status |= H5Fclose (file);

    return (int) status;
}
#endif