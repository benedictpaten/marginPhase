//
// Created by tpesout on 3/29/19.
//

#ifdef _HDF5

#include "margin.h"
#include "htsIntegration.h"
#include "helenFeatures.h"
#include "ssw.h"
#include <hdf5.h>
#include <omp.h>

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

//todo deprecated
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


//todo deprecated
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



PoaFeatureChannelRleWeight *PoaFeature_ChannelRleWeight_construct(int64_t refPos, int64_t insPos, int64_t rlPos,
                                                                  int64_t maxRunLength) {
    PoaFeatureChannelRleWeight *feature = st_calloc(1, sizeof(PoaFeatureChannelRleWeight));
    feature->refPosition = refPos;
    feature->insertPosition = insPos;
    feature->runLengthPosition = rlPos;
    feature->labelChar = '\0';
    feature->labelRunLength = 0;
    feature->nextRunLength = NULL;
    feature->nextInsert = NULL;
    feature->maxRunLength = maxRunLength;
    feature->nucleotideWeights = st_calloc((SYMBOL_NUMBER) * 2, sizeof(double));
    feature->runLengthWeights = st_calloc((SYMBOL_NUMBER - 1) * (1 + maxRunLength) * 2, sizeof(double));
    return feature;
}
void PoaFeature_ChannelRleWeight_destruct(PoaFeatureChannelRleWeight *feature) {
    if (feature->nextRunLength != NULL) {
        PoaFeature_ChannelRleWeight_destruct(feature->nextRunLength);
    }
    if (feature->nextInsert != NULL) {
        PoaFeature_ChannelRleWeight_destruct(feature->nextInsert);
    }
    free(feature->nucleotideWeights);
    free(feature->runLengthWeights);
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

int PoaFeature_ChannelRleWeight_charNuclIndex(Symbol character, bool forward) {
    int pos = character * 2 + (forward ? POS_STRAND_IDX : NEG_STRAND_IDX);
    return pos;
}
int PoaFeature_ChannelRleWeight_gapNuclIndex(bool forward) {
    int pos = (SYMBOL_NUMBER - 1) * 2 + (forward ? POS_STRAND_IDX : NEG_STRAND_IDX);
    return pos;
}
int PoaFeature_ChannelRleWeight_charRLIndex(int64_t maxRunLength, Symbol character, int64_t runLength, bool forward) {
    assert(runLength >= 0);
    assert(runLength <= maxRunLength);
    int pos = (character * ((int)maxRunLength + 1) + runLength) * 2 + (forward ? POS_STRAND_IDX : NEG_STRAND_IDX);
    return pos;
}

void handleHelenFeatures(
        // global params
        HelenFeatureType helenFeatureType, BamChunker *trueReferenceBamChunker,
        int64_t splitWeightMaxRunLength, void **helenHDF5Files, bool fullFeatureOutput,
        char *trueReferenceBam, Params *params,

        // chunk params
        char *logIdentifier, int64_t chunkIdx, BamChunk *bamChunk, Poa *poa, stList *rleReads, stList *rleNucleotides,
        char *polishedConsensusString, RleString *polishedRleConsensus) {

    st_logInfo(">%s Performing feature generation for chunk.\n", logIdentifier);

    // get filename
    char *helenFeatureOutfileBase = NULL;
    switch (helenFeatureType) {
        case HFEAT_SIMPLE_WEIGHT:
            helenFeatureOutfileBase = stString_print("simpleWeight.C%05"PRId64".%s-%"PRId64"-%"PRId64,
                                                     chunkIdx, bamChunk->refSeqName,
                                                     bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
            break;
        case HFEAT_SPLIT_RLE_WEIGHT:
            // name of folder, not of file
            helenFeatureOutfileBase = stString_print("splitRleWeight.C%05"PRId64".%s-%"PRId64"-%"PRId64,
                                                     chunkIdx, bamChunk->refSeqName,
                                                     bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
            break;
        case HFEAT_CHANNEL_RLE_WEIGHT:
            // name of folder, not of file
            helenFeatureOutfileBase = stString_print("channelRleWeight.C%05"PRId64".%s-%"PRId64"-%"PRId64,
                                                     chunkIdx, bamChunk->refSeqName,
                                                     bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
            break;
        default:
            st_errAbort("Unhandled HELEN feature type!\n");
    }

    // necessary to annotate poa with truth (if true reference BAM has been specified)
    stList *trueRefAlignment = NULL;
    RleString *trueRefRleString = NULL;
    bool validReferenceAlignment = FALSE;

    // get reference chunk
    if (trueReferenceBam != NULL) {
        // get alignment of true ref to assembly
        stList *trueRefReads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
        stList *unused = stList_construct3(0, (void (*)(void *)) stList_destruct);
        // construct new chunk
        BamChunk *trueRefBamChunk = bamChunk_copyConstruct(bamChunk);
        trueRefBamChunk->parent = trueReferenceBamChunker;
        // get true ref as "read"
        uint32_t trueAlignmentCount = convertToReadsAndAlignments(trueRefBamChunk, trueRefReads, unused);

        // poor man's "do we have a unique alignment"
        if (trueAlignmentCount == 1) {
            BamChunkRead *trueRefRead = stList_get(trueRefReads, 0);

            stList *trueRefAlignmentRawSpace = alignConsensusAndTruth(polishedConsensusString, trueRefRead->nucleotides);
            if (st_getLogLevel() == debug) {
                printMEAAlignment(polishedConsensusString, trueRefRead->nucleotides,
                                  strlen(polishedConsensusString), strlen(trueRefRead->nucleotides),
                                  trueRefAlignmentRawSpace, NULL, NULL);
            }


            // convert to rleSpace if appropriate
            if (params->polishParams->useRunLengthEncoding) {
                trueRefRleString = rleString_construct(trueRefRead->nucleotides);
                trueRefAlignment = runLengthEncodeAlignment2(trueRefAlignmentRawSpace, polishedRleConsensus,
                                                             trueRefRleString, 1, 2, 0);
                if (st_getLogLevel() == debug) {
                    printMEAAlignment(polishedRleConsensus->rleString, trueRefRleString->rleString,
                                      strlen(polishedRleConsensus->rleString),
                                      strlen(trueRefRleString->rleString),
                                      trueRefAlignment, polishedRleConsensus->repeatCounts,
                                      trueRefRleString->repeatCounts);
                }
                stList_destruct(trueRefAlignmentRawSpace);
            } else {
                trueRefRleString = rleString_constructNoRLE(trueRefRead->nucleotides);
                trueRefAlignment = trueRefAlignmentRawSpace;
            }


            // we found a single alignment of reference
            double refLengthRatio = 1.0 * trueRefRleString->length / polishedRleConsensus->length;
            double alnLengthRatio = 1.0 * stList_length(trueRefAlignment) / polishedRleConsensus->length;
            int refLengthRatioHundredthsOffOne = abs((int) (100 * (1.0 - refLengthRatio)));
            int alnLengthRatioHundredthsOffOne = abs((int) (100 * (1.0 - alnLengthRatio)));
            if (stList_length(trueRefAlignment) > 0 && refLengthRatioHundredthsOffOne < 10 &&
                alnLengthRatioHundredthsOffOne < 10) {
                validReferenceAlignment = TRUE;
            } else {
                st_logInfo(" %s True reference alignment QC failed:  polished length %"PRId64", true ref length"
                           " ratio (true/polished) %f, aligned pairs length ratio (true/polished): %f\n",
                           logIdentifier, polishedRleConsensus->length, refLengthRatio, alnLengthRatio);
            }
        }

        stList_destruct(trueRefReads);
        stList_destruct(unused);
        bamChunk_destruct(trueRefBamChunk);
    }

    // either write it, or note that we failed to find a valid reference alignment
    if (trueReferenceBam != NULL && !validReferenceAlignment) {
        st_logInfo(" %s No valid reference alignment was found, skipping HELEN feature output.\n", logIdentifier);
    } else {
        st_logInfo(" %s Writing HELEN features with filename base: %s\n", logIdentifier, helenFeatureOutfileBase);

        // write the actual features (type dependent)
        poa_writeHelenFeatures(helenFeatureType, poa, rleReads, rleNucleotides, helenFeatureOutfileBase,
                               bamChunk, trueRefAlignment, polishedRleConsensus, trueRefRleString, fullFeatureOutput,
                               splitWeightMaxRunLength, (HelenFeatureHDF5FileInfo**) helenHDF5Files);

        // write the polished chunk in fasta format
        if (fullFeatureOutput) {
            char *chunkPolishedRefFilename = stString_print("%s.fa", helenFeatureOutfileBase);
            char *chunkPolishedRefContigName = stString_print("%s\t%"PRId64"\t%"PRId64"\t%s",
                                                              bamChunk->refSeqName,
                                                              bamChunk->chunkBoundaryStart,
                                                              bamChunk->chunkBoundaryEnd,
                                                              helenFeatureOutfileBase);
            FILE *chunkPolishedRefOutFh = fopen(chunkPolishedRefFilename, "w");
            fastaWrite(polishedConsensusString, chunkPolishedRefContigName, chunkPolishedRefOutFh);
            fclose(chunkPolishedRefOutFh);
            free(chunkPolishedRefFilename);
            free(chunkPolishedRefContigName);
        }
    }

    // cleanup
    free(helenFeatureOutfileBase);
    if (trueRefAlignment != NULL) stList_destruct(trueRefAlignment);
    if (trueRefRleString != NULL) rleString_destruct(trueRefRleString);
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


void poa_addSplitRunLengthFeaturesForObservations(PoaFeatureSplitRleWeight *baseFeature, stList *observations,
                                                  stList *bamChunkReads, stList *rleStrings, const int64_t maxRunLength,
                                                  int64_t observationOffset) {


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
            if (currFeature->nextRunLength != NULL) {
                currFeature = currFeature->nextRunLength;
            } else {
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
        poa_addSplitRunLengthFeaturesForObservations(feature, node->observations, bamChunkReads, rleStrings,
                                                     maxRunLength, 0);

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
                    poa_addSplitRunLengthFeaturesForObservations(currFeature, insert->observations, bamChunkReads,
                                                                 rleStrings,
                                                                 maxRunLength, o);
                }
            }
        }
    }

    free(logIdentifier);
    return featureList;
}



void poa_addChannelRunLengthFeaturesForObservations(PoaFeatureChannelRleWeight *baseFeature, stList *observations,
                                                    stList *bamChunkReads, stList *rleStrings, const int64_t maxRunLength,
                                                    int64_t observationOffset) {

    PoaFeatureChannelRleWeight *currFeature = baseFeature;
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

            // save weight for both nucleotide totals and runLenght totals
            int64_t nuclPos = PoaFeature_ChannelRleWeight_charNuclIndex(symbol, forward);
            int64_t runLengthPos = PoaFeature_ChannelRleWeight_charRLIndex(maxRunLength, symbol, runLength, forward);
            currFeature->nucleotideWeights[nuclPos] += observation->weight;
            currFeature->runLengthWeights[runLengthPos] += observation->weight;

        }

        // update currFeature if we're going ot run again
        if (beforeMaxObservedRunLength) {
            currentRunLengthIndex++;
            if (currFeature->nextRunLength != NULL) {
                currFeature = currFeature->nextRunLength;
            } else {
                PoaFeatureChannelRleWeight *prevFeature = currFeature;
                currFeature = PoaFeature_ChannelRleWeight_construct(baseFeature->refPosition, baseFeature->insertPosition,
                                                                    currentRunLengthIndex, maxRunLength);
                prevFeature->nextRunLength = currFeature;
                currFeature->nucleotideWeights[PoaFeature_ChannelRleWeight_gapNuclIndex(TRUE)] =
                        baseFeature->nucleotideWeights[PoaFeature_ChannelRleWeight_gapNuclIndex(TRUE)];
                currFeature->nucleotideWeights[PoaFeature_ChannelRleWeight_gapNuclIndex(FALSE)] =
                        baseFeature->nucleotideWeights[PoaFeature_ChannelRleWeight_gapNuclIndex(FALSE)];
            }
        }
    }
}

stList *poa_getChannelRleWeightFeatures(Poa *poa, stList *bamChunkReads, stList *rleStrings, int64_t maxRunLength) {
    // initialize feature list
    stList *featureList = stList_construct3(0, (void (*)(void *)) PoaFeature_ChannelRleWeight_destruct);
    for(int64_t i=1; i<stList_length(poa->nodes); i++) {
        stList_append(featureList, PoaFeature_ChannelRleWeight_construct(i - 1, 0, 0, maxRunLength));
    }

    // for logging (of errors)
    char *logIdentifier = getLogIdentifier();

    // iterate over all positions
    for(int64_t i=0; i<stList_length(featureList); i++) {

        // get feature and node
        PoaFeatureChannelRleWeight* feature = stList_get(featureList, i);
        PoaNode *node = stList_get(poa->nodes, i + 1); //skip the first poa node, as it's always an 'N', so featureIdx and poaIdx are off by one

        // save run length nodes
        poa_addChannelRunLengthFeaturesForObservations(feature, node->observations, bamChunkReads, rleStrings,
                                                       maxRunLength, 0);

        // Deletes
        if (stList_length(node->deletes) > 0) {

            // iterate over all deletes
            for (int64_t d = 0; d < stList_length(node->deletes); d++) {
                PoaDelete *delete = stList_get(node->deletes, d);

                // Deletes start AFTER the current position, need to add counts/weights to nodes after the current node
                for (int64_t k = 1; k < delete->length; k++) {
                    if (i + k >= stList_length(poa->nodes)) {
                        st_logCritical(" %s Encountered DELETE that occurs after end of POA!\n", logIdentifier);
                        break;
                    }
                    PoaFeatureChannelRleWeight *delFeature = stList_get(featureList, i + k);

                    delFeature->nucleotideWeights[PoaFeature_ChannelRleWeight_gapNuclIndex(TRUE)] +=
                            delete->weightForwardStrand;
                    delFeature->nucleotideWeights[PoaFeature_ChannelRleWeight_gapNuclIndex(FALSE)] +=
                            delete->weightReverseStrand;
                }
            }
        }

        // Inserts
        if (stList_length(node->inserts) > 0) {

            // iterate over all inserts
            for (int64_t n = 0; n < stList_length(node->inserts); n++) {
                PoaInsert *insert = stList_get(node->inserts, n);

                // handle each insert base
                PoaFeatureChannelRleWeight *prevFeature = feature;
                for (int64_t o = 0; o < strlen(insert->insert); o++) {

                    // get feature iterator
                    PoaFeatureChannelRleWeight *currFeature = prevFeature->nextInsert;
                    if (currFeature == NULL) {
                        currFeature = PoaFeature_ChannelRleWeight_construct(i, o + 1, 0, maxRunLength);
                        prevFeature->nextInsert = currFeature;
                    }

                    // save insert run lengths
                    poa_addChannelRunLengthFeaturesForObservations(currFeature, insert->observations, bamChunkReads,
                                                                   rleStrings, maxRunLength, o);
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
    char *logIdentifier = getLogIdentifier();

    // iterate over true ref alignment
    stListIterator *trueRefAlignItor = stList_getIterator(trueRefAlignment);
    stIntTuple *currTrueRefAlign = stList_getNext(trueRefAlignItor);

    // iterate over features
    int64_t trueRefPos = stIntTuple_get(currTrueRefAlign, REFERENCE_POS);
    for (int64_t featureRefPos = 0; featureRefPos < stList_length(features); featureRefPos++) {
        void *feature = stList_get(features, featureRefPos);
        int64_t trueRunLength = -1;
        PoaFeatureSplitRleWeight* srlFeature = NULL;
        PoaFeatureChannelRleWeight* crlFeature = NULL;

        int64_t featureInsPos = 0;
        while (feature != NULL) {

            // no more ref bases, everything is gaps
            if (currTrueRefAlign == NULL) {
                switch (featureType) {
                    case HFEAT_SIMPLE_WEIGHT:
                        ((PoaFeatureSimpleWeight*)feature)->label = '_';
                        feature = ((PoaFeatureSimpleWeight*)feature)->nextInsert;
                        break;
                    case HFEAT_SPLIT_RLE_WEIGHT:
                        srlFeature = ((PoaFeatureSplitRleWeight*)feature);
                        while (srlFeature != NULL) {
                            srlFeature->labelChar = '_';
                            srlFeature->labelRunLength = 0;
                            srlFeature = srlFeature->nextRunLength;
                        }
                        feature = ((PoaFeatureSplitRleWeight*)feature)->nextInsert;
                        break;
                    case HFEAT_CHANNEL_RLE_WEIGHT:
                        crlFeature = ((PoaFeatureChannelRleWeight*)feature);
                        while (crlFeature != NULL) {
                            crlFeature->labelChar = '_';
                            crlFeature->labelRunLength = 0;
                            crlFeature = crlFeature->nextRunLength;
                        }
                        feature = ((PoaFeatureChannelRleWeight*)feature)->nextInsert;
                        break;
                    default:
                        st_errAbort("Unhandled FeatureType!\n");
                }
                continue;
            }

            // sanity checks
            assert(stIntTuple_get(currTrueRefAlign, FEATURE_POS) >= featureRefPos && stIntTuple_get(currTrueRefAlign, REFERENCE_POS) >= trueRefPos);

            // match
            if (stIntTuple_get(currTrueRefAlign, FEATURE_POS) == featureRefPos && stIntTuple_get(currTrueRefAlign, REFERENCE_POS) == trueRefPos) {
                st_logDebug(" %s LABEL MATCH  %c trueRefPos:%"PRId64" featureRefPos:%"PRId64" featureInsPos:%"PRId64"\n",
                           logIdentifier, featureInsPos == 0 ? ' ' : 'I', trueRefPos, featureRefPos, featureInsPos);
                // save label (based on feature type)
                switch (featureType) {
                    case HFEAT_SIMPLE_WEIGHT:
                        ((PoaFeatureSimpleWeight *) feature)->label = trueRefRleString->rleString[trueRefPos];
                        break;
                    case HFEAT_SPLIT_RLE_WEIGHT:
                        srlFeature = ((PoaFeatureSplitRleWeight*)feature);
                        trueRunLength = trueRefRleString->repeatCounts[trueRefPos];
                        while (srlFeature != NULL) {
                            srlFeature->labelChar = trueRefRleString->rleString[trueRefPos];
                            if (trueRunLength <= 0) {
                                srlFeature->labelRunLength = 0;
                            } else if (trueRunLength > srlFeature->maxRunLength) {
                                srlFeature->labelRunLength = srlFeature->maxRunLength;
                            } else {
                                srlFeature->labelRunLength = trueRunLength;
                            }
                            trueRunLength -= srlFeature->maxRunLength;
                            srlFeature = srlFeature->nextRunLength;
                        }
                        break;
                    case HFEAT_CHANNEL_RLE_WEIGHT:
                        crlFeature = ((PoaFeatureChannelRleWeight*)feature);
                        trueRunLength = trueRefRleString->repeatCounts[trueRefPos];
                        while (crlFeature != NULL) {
                            crlFeature->labelChar = trueRefRleString->rleString[trueRefPos];
                            if (trueRunLength <= 0) {
                                crlFeature->labelRunLength = 0;
                            } else if (trueRunLength > crlFeature->maxRunLength) {
                                crlFeature->labelRunLength = crlFeature->maxRunLength;
                            } else {
                                crlFeature->labelRunLength = trueRunLength;
                            }
                            trueRunLength -= crlFeature->maxRunLength;
                            crlFeature = crlFeature->nextRunLength;
                        }
                        break;
                    default:
                        st_errAbort("Unhandled FeatureType!\n");
                }

                // iterate
                trueRefPos++;
                currTrueRefAlign = stList_getNext(trueRefAlignItor);
                // handle first and last match
                if (featureInsPos == 0) {
                    if (*firstMatchedFeaure == -1) {
                        *firstMatchedFeaure = featureRefPos;
                    }
                    *lastMatchedFeature = featureRefPos;
                }
            }

            // insert
            else if (trueRefPos < stIntTuple_get(currTrueRefAlign, REFERENCE_POS)) {
                st_logDebug(" %s LABEL INSERT %c trueRefPos:%"PRId64" featureRefPos:%"PRId64" featureInsPos:%"PRId64"\n",
                           logIdentifier, featureInsPos == 0 ? ' ' : 'I', trueRefPos, featureRefPos, featureInsPos);
                // apply label
                switch (featureType) {
                    case HFEAT_SIMPLE_WEIGHT:
                        ((PoaFeatureSimpleWeight*)feature)->label = trueRefRleString->rleString[trueRefPos];
                        break;
                    case HFEAT_SPLIT_RLE_WEIGHT:
                        srlFeature = ((PoaFeatureSplitRleWeight*)feature);
                        trueRunLength = trueRefRleString->repeatCounts[trueRefPos];
                        while (srlFeature != NULL) {
                            srlFeature->labelChar = trueRefRleString->rleString[trueRefPos];
                            if (trueRunLength <= 0) {
                                srlFeature->labelRunLength = 0;
                            } else if (trueRunLength > srlFeature->maxRunLength) {
                                srlFeature->labelRunLength = srlFeature->maxRunLength;
                            } else {
                                srlFeature->labelRunLength = trueRunLength;
                            }
                            trueRunLength -= srlFeature->maxRunLength;
                            srlFeature = srlFeature->nextRunLength;
                        }
                        break;
                    case HFEAT_CHANNEL_RLE_WEIGHT:
                        crlFeature = ((PoaFeatureChannelRleWeight*)feature);
                        trueRunLength = trueRefRleString->repeatCounts[trueRefPos];
                        while (crlFeature != NULL) {
                            crlFeature->labelChar = trueRefRleString->rleString[trueRefPos];
                            if (trueRunLength <= 0) {
                                crlFeature->labelRunLength = 0;
                            } else if (trueRunLength > crlFeature->maxRunLength) {
                                crlFeature->labelRunLength = crlFeature->maxRunLength;
                            } else {
                                crlFeature->labelRunLength = trueRunLength;
                            }
                            trueRunLength -= crlFeature->maxRunLength;
                            crlFeature = crlFeature->nextRunLength;
                        }
                        break;
                    default:
                        st_errAbort("Unhandled FeatureType!\n");
                }
                trueRefPos++;
            }

            // delete
            else if (featureRefPos < stIntTuple_get(currTrueRefAlign, FEATURE_POS)) {
                st_logDebug(" %s LABEL DELETE %c trueRefPos:%"PRId64" featureRefPos:%"PRId64" featureInsPos:%"PRId64"\n",
                           logIdentifier, featureInsPos == 0 ? ' ' : 'I', trueRefPos, featureRefPos, featureInsPos);
                // apply label
                switch (featureType) {
                    case HFEAT_SIMPLE_WEIGHT:
                        ((PoaFeatureSimpleWeight*)feature)->label = '_';
                        break;
                    case HFEAT_SPLIT_RLE_WEIGHT:
                        srlFeature = ((PoaFeatureSplitRleWeight*)feature);
                        while (srlFeature != NULL) {
                            srlFeature->labelChar = '_';
                            srlFeature->labelRunLength = 0;
                            srlFeature = srlFeature->nextRunLength;
                        }
                        break;
                    case HFEAT_CHANNEL_RLE_WEIGHT:
                        crlFeature = ((PoaFeatureChannelRleWeight*)feature);
                        while (crlFeature != NULL) {
                            crlFeature->labelChar = '_';
                            crlFeature->labelRunLength = 0;
                            crlFeature = crlFeature->nextRunLength;
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
            switch (featureType) {
                case HFEAT_SIMPLE_WEIGHT:
                    feature = ((PoaFeatureSimpleWeight*)feature)->nextInsert;
                    break;
                case HFEAT_SPLIT_RLE_WEIGHT:
                    feature = ((PoaFeatureSplitRleWeight*)feature)->nextInsert;
                    break;
                case HFEAT_CHANNEL_RLE_WEIGHT:
                    feature = ((PoaFeatureChannelRleWeight*)feature)->nextInsert;
                    break;
                default:
                    st_errAbort("Unhandled FeatureType!\n");
            }
            featureInsPos++;
        }

        // this catches any true inserts which are not present in the poa / feature list
        while (currTrueRefAlign != NULL &&
                featureRefPos < stIntTuple_get(currTrueRefAlign, FEATURE_POS) &&
                trueRefPos < stIntTuple_get(currTrueRefAlign, REFERENCE_POS)) {
            trueRefPos++;
        }
    }

    stList_destructIterator(trueRefAlignItor);
    free(logIdentifier);
}

void poa_writeHelenFeatures(HelenFeatureType type, Poa *poa, stList *bamChunkReads, stList *rleStrings,
        char *outputFileBase, BamChunk *bamChunk, stList *trueRefAlignment, RleString *consensusRleString,
        RleString *trueRefRleString, bool fullFeatureOutput, int64_t maxRunLength,
        HelenFeatureHDF5FileInfo** helenHDF5Files) {
    // prep
    int64_t firstMatchedFeature = -1;
    int64_t lastMatchedFeature = -1;
    stList *features = NULL;
    bool outputLabels = trueRefAlignment != NULL && trueRefRleString != NULL;
    int64_t threadIdx = omp_get_thread_num();

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

            writeSimpleWeightHelenFeaturesHDF5(helenHDF5Files[threadIdx], outputFileBase, bamChunk, outputLabels,
                    features, firstMatchedFeature, lastMatchedFeature);

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

            writeSplitRleWeightHelenFeaturesHDF5(helenHDF5Files[threadIdx],
                    outputFileBase, bamChunk, outputLabels, features, firstMatchedFeature, lastMatchedFeature,
                    maxRunLength);
            break;

        case HFEAT_CHANNEL_RLE_WEIGHT:
            // get features
            features = poa_getChannelRleWeightFeatures(poa, bamChunkReads, rleStrings, maxRunLength);
            firstMatchedFeature = 0;
            lastMatchedFeature = stList_length(features) - 1;

            // get truth (if we have it)
            if (outputLabels) {
                poa_annotateHelenFeaturesWithTruth(features, type, trueRefAlignment, trueRefRleString,
                                                   &firstMatchedFeature, &lastMatchedFeature);
            }

            writeChannelRleWeightHelenFeaturesHDF5(helenHDF5Files[threadIdx],
                                                 outputFileBase, bamChunk, outputLabels, features, firstMatchedFeature,
                                                 lastMatchedFeature, maxRunLength);
            break;

        default:
            st_errAbort("Unhandled HELEN feature type!\n");
    }

    //cleanup
    stList_destruct(features);
}

// this function taken from https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/blob/master/src/example.c
stList *alignConsensusAndTruth(char *consensusStr, char *truthStr) {

    int64_t l, m, k;
    uint8_t match = 2, mismatch = 2, gap_open = 3, gap_extension = 1;	// default parameters for genome sequence alignment

    int64_t consensusLen = strlen(consensusStr);
    int64_t truthLen = strlen(truthStr);

    s_profile* profile;

    int8_t* num = st_calloc(consensusLen, sizeof(int8_t));	// the read sequence represented in numbers
    int8_t* ref_num = st_calloc(truthLen, sizeof(int8_t));  // the read sequence represented in numbers
    s_align* result;

    /* This table is used to transform nucleotide letters into numbers. */
    static const int8_t nt_table[128] = {
            4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
            4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
            4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
            4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
            4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
            4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
            4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
            4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
    };

    // initialize scoring matrix for genome sequences
    //  A  C  G  T	N (or other ambiguous code)
    //  2 -2 -2 -2 	0	A
    // -2  2 -2 -2 	0	C
    // -2 -2  2 -2 	0	G
    // -2 -2 -2  2 	0	T
    //	0  0  0  0  0	N (or other ambiguous code)
    int8_t* mat = (int8_t*)calloc(25, sizeof(int8_t));
    for (l = k = 0; l < 4; ++l) {
        for (m = 0; m < 4; ++m) mat[k++] = (int8_t ) (l == m ? match : - mismatch);	/* weight_match : -weight_mismatch */
        mat[k++] = 0; // ambiguous base: no penalty
    }
    for (m = 0; m < 5; ++m) mat[k++] = 0;

    for (m = 0; m < consensusLen; ++m) num[m] = nt_table[(int) consensusStr[m]];
    profile = ssw_init(num, (int32_t) consensusLen, mat, 5, 2);
    for (m = 0; m < truthLen; ++m) ref_num[m] = nt_table[(int) truthStr[m]];

    // Only the 8 bit of the flag is setted. ssw_align will always return the best alignment beginning position and cigar.
    result = ssw_align (profile, ref_num, (int32_t) truthLen, gap_open, gap_extension, 1, 0, 0, 15);

    // Convert from cigar to aligned pairs
    int32_t consensusPos = result->read_begin1;
    int32_t truthPos = result->ref_begin1;
    stList *alignedPairs = stList_construct3(0, (void(*)(void*)) stIntTuple_destruct);
    if (result->cigar) {
        for (int32_t cigIdx = 0; cigIdx < result->cigarLen; cigIdx++) {
            char letter = cigar_int_to_op(result->cigar[cigIdx]);
            uint32_t length = cigar_int_to_len(result->cigar[cigIdx]);

            if (letter == 'M') {
                for (int32_t matchIdx = 0; matchIdx < length; matchIdx++) {
                    stList_append(alignedPairs, stIntTuple_construct3(0, consensusPos, truthPos));
                    consensusPos++;
                    truthPos++;
                }
            } else if (letter == 'I') {
                consensusPos += length;
            } else if (letter == 'D') {
                truthPos += length;
            }
        }
    }

    // cleanup
    align_destroy(result);
    init_destroy(profile);
    free(mat);
    free(ref_num);
    free(num);

    return alignedPairs;
}

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

uint8_t ***getThreeDArrayUInt8(int64_t depthCount, int64_t rowCount, int64_t columnCount) {
    uint8_t ***array  = st_calloc(depthCount, sizeof(uint8_t**));
    array[0] = (uint8_t**) st_calloc(depthCount * rowCount, sizeof(uint8_t*) );
    array[0][0] = (uint8_t*) st_calloc(depthCount * rowCount * columnCount, sizeof(uint8_t) );
    for (int64_t i=0; i < depthCount; i++) {
        array[i] = array[0] + i*rowCount;
        for (int64_t j=0; j < rowCount; j++) {
            array[i][j] = (array[0][0] + i*rowCount*columnCount + j*columnCount);
        }
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


void writeSimpleWeightHelenFeaturesHDF5(HelenFeatureHDF5FileInfo* hdf5FileInfo, char *outputFileBase,
        BamChunk *bamChunk, bool outputLabels, stList *features, int64_t featureStartIdx, int64_t featureEndIdxInclusive) {

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

    hid_t status;
    hid_t stringType = H5Tcopy (H5T_C_S1);
    status = H5Tset_size (stringType, strlen(bamChunk->refSeqName) + 1);

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

        // create group
        char *outputGroup = stString_print("images/%s.%"PRId64, outputFileBase, featureIndex);
        hid_t group = H5Gcreate (hdf5FileInfo->file, outputGroup, hdf5FileInfo->groupPropertyList, H5P_DEFAULT, H5P_DEFAULT);

        // write metadata
        hid_t contigDataset = H5Dcreate (group, "contig", stringType, metadataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (contigDataset, stringType, H5S_ALL, H5S_ALL, H5P_DEFAULT, bamChunk->refSeqName);
        hid_t contigStartDataset = H5Dcreate (group, "contig_start", hdf5FileInfo->int64Type, metadataSpace,
                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (contigStartDataset, hdf5FileInfo->int64Type,
                H5S_ALL, H5S_ALL, H5P_DEFAULT, &bamChunk->chunkBoundaryStart);
        hid_t contigEndDataset = H5Dcreate (group, "contig_end", hdf5FileInfo->int64Type, metadataSpace,
                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (contigEndDataset, hdf5FileInfo->int64Type,
                H5S_ALL, H5S_ALL, H5P_DEFAULT, &bamChunk->chunkBoundaryEnd);
        hid_t chunkIndexDataset = H5Dcreate (group, "feature_chunk_idx", hdf5FileInfo->int64Type, metadataSpace,
                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (chunkIndexDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &featureIndex);

        // write position info
        hid_t positionDataset = H5Dcreate (group, "position", hdf5FileInfo->uint32Type, positionSpace,
                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (positionDataset, hdf5FileInfo->uint32Type,
                H5S_ALL, H5S_ALL, H5P_DEFAULT, positionData[chunkFeatureStartIdx]);

        // write rle data
        hid_t rleWeightDataset = H5Dcreate (group, "image", hdf5FileInfo->floatType, rleWeightSpace,
                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (rleWeightDataset, hdf5FileInfo->floatType, H5S_ALL, H5S_ALL, H5P_DEFAULT, rleWeightData[chunkFeatureStartIdx]);

        // if labels, add all these too
        if (outputLabels) {
            hid_t labelCharacterDataset = H5Dcreate (group, "label_base", hdf5FileInfo->uint8Type, labelCharacterSpace,
                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite (labelCharacterDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
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
        status |= H5Gclose (group);
        free(outputGroup);
    }


    // cleanup
    free(rleWeightData[0]);
    free(rleWeightData);
    free(positionData[0]);
    free(positionData);
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

//
//
//void writeSimpleWeightHelenFeaturesHDF5(char *outputFileBase, BamChunk *bamChunk, bool outputLabels, stList *features,
//                                        int64_t featureStartIdx, int64_t featureEndIdxInclusive) {
//
//    herr_t      status = 0;
//
//    /*
//     * Get feature data set up
//     */
//
//    // count features, create feature array
//    uint64_t featureCount = 0;
//    for (int64_t i = featureStartIdx; i <= featureEndIdxInclusive; i++) {
//        PoaFeatureSimpleWeight *feature = stList_get(features, i);
//        while (feature != NULL) {
//            featureCount++;
//            feature = feature->nextInsert;
//        }
//    }
//    if (featureCount < HDF5_FEATURE_SIZE && outputLabels) {
//        char *logIdentifier = getLogIdentifier();
//        st_logInfo(" %s Feature count %"PRId64" less than minimum of %d\n", logIdentifier, featureCount, HDF5_FEATURE_SIZE);
//        free(logIdentifier);
//        return;
//    }
//
//
//    // get all feature data into an array
//    uint32_t **positionData = getTwoDArrayUInt32(featureCount, 2);
//    int64_t columnCount = SYMBOL_NUMBER * 2; //{A, C, T, G, Gap} x {fwd, rev}
//    float **rleWeightData = getTwoDArrayFloat(featureCount, columnCount);
//    char **labelCharacterData = NULL;
//    if (outputLabels) {
//        labelCharacterData = getTwoDArrayChar(featureCount, 1);
//    }
//
//    // add all data to features
//    featureCount = 0;
//    for (int64_t i = featureStartIdx; i <= featureEndIdxInclusive; i++) {
//        PoaFeatureSimpleWeight *feature = stList_get(features, i);
//        while (feature != NULL) {
//            positionData[featureCount][0] = (uint32_t) feature->refPosition;
//            positionData[featureCount][1] = (uint32_t) feature->insertPosition;
//
//            for (int64_t symbol = 0; symbol < SYMBOL_NUMBER_NO_N; symbol++) {
//                int64_t pos = PoaFeature_SimpleWeight_charIndex((Symbol) symbol, TRUE);
//                rleWeightData[featureCount][pos] = (float) (feature->weights[pos] / PAIR_ALIGNMENT_PROB_1);
//                pos = PoaFeature_SimpleWeight_charIndex((Symbol) symbol, FALSE);
//                rleWeightData[featureCount][pos] = (float) (feature->weights[pos] / PAIR_ALIGNMENT_PROB_1);
//            }
//
//            // weights include 'N' index which is not included in features
//            int64_t pos = PoaFeature_SimpleWeight_gapIndex(TRUE);
//            rleWeightData[featureCount][pos - 2] = (float) (feature->weights[pos] / PAIR_ALIGNMENT_PROB_1);
//            pos = PoaFeature_SimpleWeight_gapIndex(FALSE);
//            rleWeightData[featureCount][pos - 2] = (float) (feature->weights[pos] / PAIR_ALIGNMENT_PROB_1);
//
//            if (outputLabels) {
//                labelCharacterData[featureCount][0] = feature->label;
//            }
//
//            featureCount++;
//            feature = feature->nextInsert;
//        }
//    }
//
//    /*
//     * Get hdf5 data set up
//     */
//
//    hid_t int64Type = H5Tcopy(H5T_NATIVE_UINT32);
//    status |= H5Tset_order(int64Type, H5T_ORDER_LE);
//    hid_t uint32Type = H5Tcopy(H5T_NATIVE_UINT32);
//    status |= H5Tset_order(uint32Type, H5T_ORDER_LE);
//    hid_t floatType = H5Tcopy(H5T_NATIVE_FLOAT);
//    status |= H5Tset_order(floatType, H5T_ORDER_LE);
//    hid_t uint8Type = H5Tcopy(H5T_NATIVE_UINT8);
//    status |= H5Tset_order(uint8Type, H5T_ORDER_LE);
//    hid_t stringType = H5Tcopy (H5T_C_S1);
//    status |= H5Tset_size (stringType, strlen(bamChunk->refSeqName) + 1);
//
//    // so that we can produce chunks smaller than HDF5_FEATURE_SIZE (not used during training)
//    hsize_t featureSize = (hsize_t) (featureCount < HDF5_FEATURE_SIZE ? featureCount : HDF5_FEATURE_SIZE);
//
//    hsize_t metadataDimension[1] = {1};
//    hsize_t postionDimension[2] = {featureSize, 2};
//    hsize_t labelCharacterDimension[2] = {featureSize, 1};
//    hsize_t rleWeightDimension[2] = {featureSize, (hsize_t) columnCount};
//
//    hid_t metadataSpace = H5Screate_simple(1, metadataDimension, NULL);
//    hid_t positionSpace = H5Screate_simple(2, postionDimension, NULL);
//    hid_t labelCharacterSpace = H5Screate_simple(2, labelCharacterDimension, NULL);
//    hid_t rleWeightSpace = H5Screate_simple(2, rleWeightDimension, NULL);
//
//    /*
//     * Write features to files
//     */
//
//    // each file must have exactly 1000 features
//    int64_t totalFeatureFiles = (int64_t) (featureCount / HDF5_FEATURE_SIZE) + (featureCount % HDF5_FEATURE_SIZE == 0 ? 0 : 1);
//    int64_t featureOffset = 0;
//    if (featureCount >= HDF5_FEATURE_SIZE) {
//        featureOffset = (int64_t) ((HDF5_FEATURE_SIZE * totalFeatureFiles - featureCount) / (int64_t) (featureCount / HDF5_FEATURE_SIZE));
//    }
//
//    for (int64_t featureIndex = 0; featureIndex < totalFeatureFiles; featureIndex++) {
//        // get start pos
//        int64_t chunkFeatureStartIdx = (HDF5_FEATURE_SIZE * featureIndex) - (featureOffset * featureIndex);
//        if (featureIndex + 1 == totalFeatureFiles && featureCount >= HDF5_FEATURE_SIZE) {
//            chunkFeatureStartIdx = featureCount - HDF5_FEATURE_SIZE;
//        }
//
//        // create file
//        char *outputFile = stString_print("%s.%"PRId64".h5", outputFileBase, featureIndex);
//        hid_t file = H5Fcreate (outputFile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
//
//        // write metadata
//        hid_t contigDataset = H5Dcreate (file, "contig", stringType, metadataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//        status |= H5Dwrite (contigDataset, stringType, H5S_ALL, H5S_ALL, H5P_DEFAULT, bamChunk->refSeqName);
//        hid_t contigStartDataset = H5Dcreate (file, "contig_start", int64Type, metadataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//        status |= H5Dwrite (contigStartDataset, int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bamChunk->chunkBoundaryStart);
//        hid_t contigEndDataset = H5Dcreate (file, "contig_end", int64Type, metadataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//        status |= H5Dwrite (contigEndDataset, int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bamChunk->chunkBoundaryEnd);
//        hid_t chunkIndexDataset = H5Dcreate (file, "feature_chunk_idx", int64Type, metadataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//        status |= H5Dwrite (chunkIndexDataset, int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &featureIndex);
//
//        // write position info
//        hid_t positionDataset = H5Dcreate (file, "position", uint32Type, positionSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//        status |= H5Dwrite (positionDataset, uint32Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, positionData[chunkFeatureStartIdx]);
//
//        // write rle data
//        hid_t rleWeightDataset = H5Dcreate (file, "image", floatType, rleWeightSpace, H5P_DEFAULT, H5P_DEFAULT,
//                                            H5P_DEFAULT);
//        status |= H5Dwrite (rleWeightDataset, floatType, H5S_ALL, H5S_ALL, H5P_DEFAULT, rleWeightData[chunkFeatureStartIdx]);
//
//        // if labels, add all these too
//        if (outputLabels) {
//            hid_t labelCharacterDataset = H5Dcreate (file, "label_base", uint8Type, labelCharacterSpace, H5P_DEFAULT,
//                                                     H5P_DEFAULT, H5P_DEFAULT);
//            status |= H5Dwrite (labelCharacterDataset, uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
//                                labelCharacterData[chunkFeatureStartIdx]);
//            status |= H5Dclose (labelCharacterDataset);
//        }
//
//        // cleanup
//        status |= H5Dclose (contigDataset);
//        status |= H5Dclose (contigStartDataset);
//        status |= H5Dclose (contigEndDataset);
//        status |= H5Dclose (chunkIndexDataset);
//        status |= H5Dclose (positionDataset);
//        status |= H5Dclose (rleWeightDataset);
//        status |= H5Fclose (file);
//        free(outputFile);
//    }
//
//    // cleanup
//    free(rleWeightData[0]);
//    free(rleWeightData);
//    free(positionData[0]);
//    free(positionData);
//    status |= H5Tclose (int64Type);
//    status |= H5Tclose (uint32Type);
//    status |= H5Tclose (floatType);
//    status |= H5Tclose (uint8Type);
//    status |= H5Tclose (stringType);
//    status |= H5Sclose (metadataSpace);
//    status |= H5Sclose (rleWeightSpace);
//    status |= H5Sclose (positionSpace);
//    status |= H5Sclose (labelCharacterSpace);
//    if (outputLabels) {
//        free(labelCharacterData[0]);
//        free(labelCharacterData);
//    }
//
//    if (status) {
//        char *logIdentifier = getLogIdentifier();
//        st_logInfo(" %s Error writing HELEN features to HDF5 files: %s\n", logIdentifier, outputFileBase);
//        free(logIdentifier);
//    }
//}

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


void writeSplitRleWeightHelenFeaturesHDF5(HelenFeatureHDF5FileInfo* hdf5FileInfo, char *outputFileBase, BamChunk *bamChunk,
                                          bool outputLabels, stList *features, int64_t featureStartIdx, int64_t featureEndIdxInclusive,
                                          const int64_t maxRunLength) {

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

    hid_t stringType = H5Tcopy (H5T_C_S1);
    H5Tset_size (stringType, strlen(bamChunk->refSeqName) + 1);


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

        // create group
        char *outputGroup = stString_print("images/%s.%"PRId64, outputFileBase, featureIndex);
        hid_t group = H5Gcreate (hdf5FileInfo->file, outputGroup, hdf5FileInfo->groupPropertyList, H5P_DEFAULT, H5P_DEFAULT);

        // write metadata
        hid_t contigDataset = H5Dcreate (group, "contig", stringType, metadataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (contigDataset, stringType, H5S_ALL, H5S_ALL, H5P_DEFAULT, bamChunk->refSeqName);
        hid_t contigStartDataset = H5Dcreate (group, "contig_start", hdf5FileInfo->int64Type, metadataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (contigStartDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bamChunk->chunkBoundaryStart);
        hid_t contigEndDataset = H5Dcreate (group, "contig_end", hdf5FileInfo->int64Type, metadataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (contigEndDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bamChunk->chunkBoundaryEnd);
        hid_t chunkIndexDataset = H5Dcreate (group, "feature_chunk_idx", hdf5FileInfo->int64Type, metadataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (chunkIndexDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &featureIndex);

        // write position info
        hid_t positionDataset = H5Dcreate (group, "position", hdf5FileInfo->uint32Type, positionSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (positionDataset, hdf5FileInfo->uint32Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, positionData[chunkFeatureStartIdx]);

        // write rle data
        hid_t imageDataset = H5Dcreate (group, "image", hdf5FileInfo->uint8Type, imageSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (imageDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                            imageData[chunkFeatureStartIdx]);
        hid_t normalizationDataset = H5Dcreate (group, "normalization", hdf5FileInfo->uint8Type, normalizationSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (normalizationDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                            normalizationData[chunkFeatureStartIdx]);

        // if labels, add all these too
        if (outputLabels) {
            hid_t labelCharacterDataset = H5Dcreate (group, "label_base", hdf5FileInfo->uint8Type, labelCharacterSpace, H5P_DEFAULT,
                                                     H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite (labelCharacterDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                                labelCharacterData[chunkFeatureStartIdx]);
            hid_t labelRunLengthDataset = H5Dcreate (group, "label_run_length", hdf5FileInfo->uint8Type, labelRunLengthSpace,
                                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite (labelRunLengthDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
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
        status |= H5Gclose (group);
        free(outputGroup);
    }

    // cleanup
    free(imageData[0]);
    free(imageData);
    free(normalizationData[0]);
    free(normalizationData);
    free(positionData[0]);
    free(positionData);
    status |= H5Sclose (metadataSpace);
    status |= H5Sclose (positionSpace);
    status |= H5Sclose (imageSpace);
    status |= H5Sclose (normalizationSpace);
    status |= H5Sclose (labelRunLengthSpace);
    status |= H5Sclose (labelCharacterSpace);
    status |= H5Tclose (stringType);
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



void writeChannelRleWeightHelenFeaturesHDF5(HelenFeatureHDF5FileInfo* hdf5FileInfo, char *outputFileBase,
        BamChunk *bamChunk, bool outputLabels, stList *features, int64_t featureStartIdx,
        int64_t featureEndIdxInclusive, const int64_t maxRunLength) {

    herr_t      status = 0;

    /*
     * Get feature data set up
     */

    // count features, create feature array
    uint64_t featureCount = 0;
    for (int64_t i = featureStartIdx; i <= featureEndIdxInclusive; i++) {
        PoaFeatureChannelRleWeight *feature = stList_get(features, i);
        while (feature != NULL) {
            PoaFeatureChannelRleWeight *rlFeature = feature;
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

    // sizes
    int64_t nucleotideColumnCount = SYMBOL_NUMBER * 2; //ACGTGap x {fwd,bwd}
    int64_t runLengthColumnCount = (maxRunLength + 1) * 2; //(runLenght + 0) x {fwd,bwd}

    // get all feature data into an array
    uint32_t **positionData = getTwoDArrayUInt32(featureCount, 3);
    uint8_t **normalizationData = getTwoDArrayUInt8(featureCount, 1);
    uint8_t **nucleotideData = getTwoDArrayUInt8(featureCount, nucleotideColumnCount);
    uint8_t ***runLengthData = getThreeDArrayUInt8(featureCount, runLengthColumnCount, (SYMBOL_NUMBER - 1));
    //uint8_t ***runLengthData = getThreeDArrayUInt8(featureCount, (SYMBOL_NUMBER - 1), runLengthColumnCount);
    uint8_t **labelCharacterData = NULL;
    uint8_t **labelRunLengthData = NULL;
    if (outputLabels) {
        labelCharacterData = getTwoDArrayUInt8(featureCount, 1);
        labelRunLengthData = getTwoDArrayUInt8(featureCount, 1);
    }

    // add all data to features
    featureCount = 0;
    for (int64_t i = featureStartIdx; i <= featureEndIdxInclusive; i++) {
        PoaFeatureChannelRleWeight *refFeature = stList_get(features, i);

        // total weight is calculated for the very first refPos feature, used for all inserts and run lengths
        double totalWeight = 0;
        for (int64_t j = 0; j < nucleotideColumnCount; j++) {
            totalWeight += refFeature->nucleotideWeights[j];
        }

        // iterate over all insert features
        PoaFeatureChannelRleWeight *insFeature = refFeature;
        while (insFeature != NULL) {

            // iterate over all run length features
            PoaFeatureChannelRleWeight *rlFeature = insFeature;
            while (rlFeature != NULL) {
                // position
                positionData[featureCount][0] = (uint32_t) rlFeature->refPosition;
                positionData[featureCount][1] = (uint32_t) rlFeature->insertPosition;
                positionData[featureCount][2] = (uint32_t) rlFeature->runLengthPosition;

                // normalization
                normalizationData[featureCount][0] = convertTotalWeightToUInt8(totalWeight);

                // copy weights over (into normalized uint8 space)
                for (int64_t c = 0; c < SYMBOL_NUMBER - 1; c++) {
                    // overall nucl count
                    nucleotideData[featureCount][c * 2 + POS_STRAND_IDX] = normalizeWeightToUInt8(totalWeight,
                            rlFeature->nucleotideWeights[PoaFeature_ChannelRleWeight_charNuclIndex(c, TRUE)]);
                    nucleotideData[featureCount][c * 2 + NEG_STRAND_IDX] = normalizeWeightToUInt8(totalWeight,
                            rlFeature->nucleotideWeights[PoaFeature_ChannelRleWeight_charNuclIndex(c, FALSE)]);

                    // run length counts
                    for (int64_t r = 0; r <= maxRunLength; r++) {
                        runLengthData[featureCount][r * 2 + POS_STRAND_IDX][c] =
                                normalizeWeightToUInt8(totalWeight, rlFeature->runLengthWeights[
                                        PoaFeature_ChannelRleWeight_charRLIndex(maxRunLength, c, r, TRUE)]);
                        runLengthData[featureCount][r * 2 + NEG_STRAND_IDX][c] =
                                normalizeWeightToUInt8(totalWeight, rlFeature->runLengthWeights[
                                        PoaFeature_ChannelRleWeight_charRLIndex(maxRunLength, c, r, FALSE)]);

                        //runLengthData[featureCount][c][r * 2 + POS_STRAND_IDX] =
                        //        normalizeWeightToUInt8(totalWeight, rlFeature->runLengthWeights[
                        //                PoaFeature_ChannelRleWeight_charRLIndex(maxRunLength, c, r, TRUE)]);
                        //runLengthData[featureCount][c][r * 2 + NEG_STRAND_IDX] =
                        //        normalizeWeightToUInt8(totalWeight, rlFeature->runLengthWeights[
                        //                PoaFeature_ChannelRleWeight_charRLIndex(maxRunLength, c, r, FALSE)]);
                    }
                }
                // gap counts
                nucleotideData[featureCount][SYMBOL_NUMBER_NO_N * 2 + 0 + POS_STRAND_IDX] = normalizeWeightToUInt8(
                        totalWeight, rlFeature->nucleotideWeights[PoaFeature_ChannelRleWeight_gapNuclIndex(TRUE)]);
                nucleotideData[featureCount][SYMBOL_NUMBER_NO_N * 2 + NEG_STRAND_IDX] = normalizeWeightToUInt8(
                        totalWeight, rlFeature->nucleotideWeights[PoaFeature_ChannelRleWeight_gapNuclIndex(FALSE)]);

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


    // so that we can produce chunks smaller than HDF5_FEATURE_SIZE (not used during training)
    hsize_t featureSize = (hsize_t) (featureCount < HDF5_FEATURE_SIZE ? featureCount : HDF5_FEATURE_SIZE);

    hsize_t metadataDimension[1] = {1};
    hsize_t postionDimension[2] = {featureSize, 3};
    hsize_t labelCharacterDimension[2] = {featureSize, 1};
    hsize_t labelRunLengthDimension[2] = {featureSize, 1};
    hsize_t normalizationDimension[2] = {featureSize, 1};
    hsize_t nucleotideDimension[2] = {featureSize, (hsize_t) nucleotideColumnCount};
    hsize_t runLengthDimension[3] = {featureSize, (hsize_t) runLengthColumnCount, (hsize_t) SYMBOL_NUMBER - 1};
//    hsize_t runLengthDimension[3] = {featureSize, (hsize_t) SYMBOL_NUMBER - 1, (hsize_t) runLengthColumnCount};

    hid_t metadataSpace = H5Screate_simple(1, metadataDimension, NULL);
    hid_t positionSpace = H5Screate_simple(2, postionDimension, NULL);
    hid_t labelCharacterSpace = H5Screate_simple(2, labelCharacterDimension, NULL);
    hid_t labelRunLengthSpace = H5Screate_simple(2, labelRunLengthDimension, NULL);
    hid_t normalizationSpace = H5Screate_simple(2, normalizationDimension, NULL);
    hid_t nucleotideSpace = H5Screate_simple(2, nucleotideDimension, NULL);
    hid_t runLengthSpace = H5Screate_simple(3, runLengthDimension, NULL);

    hid_t stringType = H5Tcopy (H5T_C_S1);
    H5Tset_size (stringType, strlen(bamChunk->refSeqName) + 1);


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

        // create group
        char *outputGroup = stString_print("images/%s.%"PRId64, outputFileBase, featureIndex);
        hid_t group = H5Gcreate (hdf5FileInfo->file, outputGroup, hdf5FileInfo->groupPropertyList, H5P_DEFAULT, H5P_DEFAULT);

        // write metadata
        hid_t contigDataset = H5Dcreate (group, "contig", stringType, metadataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (contigDataset, stringType, H5S_ALL, H5S_ALL, H5P_DEFAULT, bamChunk->refSeqName);
        hid_t contigStartDataset = H5Dcreate (group, "contig_start", hdf5FileInfo->int64Type, metadataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (contigStartDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bamChunk->chunkBoundaryStart);
        hid_t contigEndDataset = H5Dcreate (group, "contig_end", hdf5FileInfo->int64Type, metadataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (contigEndDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bamChunk->chunkBoundaryEnd);
        hid_t chunkIndexDataset = H5Dcreate (group, "feature_chunk_idx", hdf5FileInfo->int64Type, metadataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (chunkIndexDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &featureIndex);

        // write position info
        hid_t positionDataset = H5Dcreate (group, "position", hdf5FileInfo->uint32Type, positionSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (positionDataset, hdf5FileInfo->uint32Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, positionData[chunkFeatureStartIdx]);

        // write nucl, rl, and norm data
        hid_t nucleotideDataset = H5Dcreate (group, "nucleotide", hdf5FileInfo->uint8Type, nucleotideSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (nucleotideDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                            nucleotideData[chunkFeatureStartIdx]);
        hid_t runLengthDataset = H5Dcreate (group, "runLengths", hdf5FileInfo->uint8Type, runLengthSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (runLengthDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                            runLengthData[chunkFeatureStartIdx][0]);
        hid_t normalizationDataset = H5Dcreate (group, "normalization", hdf5FileInfo->uint8Type, normalizationSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite (normalizationDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                            normalizationData[chunkFeatureStartIdx]);

        // if labels, add all these too
        if (outputLabels) {
            hid_t labelCharacterDataset = H5Dcreate (group, "label_base", hdf5FileInfo->uint8Type, labelCharacterSpace, H5P_DEFAULT,
                                                     H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite (labelCharacterDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                                labelCharacterData[chunkFeatureStartIdx]);
            hid_t labelRunLengthDataset = H5Dcreate (group, "label_run_length", hdf5FileInfo->uint8Type, labelRunLengthSpace,
                                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite (labelRunLengthDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
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
        status |= H5Dclose (nucleotideDataset);
        status |= H5Dclose (runLengthDataset);
        status |= H5Dclose (normalizationDataset);
        status |= H5Gclose (group);
        free(outputGroup);
    }

    // cleanup
    free(nucleotideData[0]);
    free(nucleotideData);
    free(runLengthData[0][0]);
    free(runLengthData[0]);
    free(runLengthData);
    free(normalizationData[0]);
    free(normalizationData);
    free(positionData[0]);
    free(positionData);
    status |= H5Sclose (metadataSpace);
    status |= H5Sclose (positionSpace);
    status |= H5Sclose (nucleotideSpace);
    status |= H5Sclose (runLengthSpace);
    status |= H5Sclose (normalizationSpace);
    status |= H5Sclose (labelRunLengthSpace);
    status |= H5Sclose (labelCharacterSpace);
    status |= H5Tclose (stringType);
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



HelenFeatureHDF5FileInfo* HelenFeatureHDF5FileInfo_construct(char *filename) {
    HelenFeatureHDF5FileInfo *fileInfo = st_calloc(1, sizeof(HelenFeatureHDF5FileInfo));
    fileInfo->filename = stString_copy(filename);
    fileInfo->file = H5Fcreate (filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    fileInfo->int64Type = H5Tcopy(H5T_NATIVE_UINT32);
    H5Tset_order(fileInfo->int64Type, H5T_ORDER_LE);
    fileInfo->uint32Type = H5Tcopy(H5T_NATIVE_UINT32);
    H5Tset_order(fileInfo->uint32Type, H5T_ORDER_LE);
    fileInfo->uint8Type = H5Tcopy(H5T_NATIVE_UINT8);
    H5Tset_order(fileInfo->uint8Type, H5T_ORDER_LE);
    fileInfo->floatType = H5Tcopy(H5T_NATIVE_FLOAT);
    H5Tset_order(fileInfo->floatType, H5T_ORDER_LE);
    fileInfo->groupPropertyList = H5Pcreate (H5P_LINK_CREATE);
    H5Pset_create_intermediate_group (fileInfo->groupPropertyList, 1);
    return fileInfo;
}

void HelenFeatureHDF5FileInfo_destruct(HelenFeatureHDF5FileInfo *fileInfo) {
    free(fileInfo->filename);
    H5Tclose(fileInfo->int64Type);
    H5Tclose(fileInfo->uint32Type);
    H5Tclose(fileInfo->uint8Type);
    H5Tclose(fileInfo->floatType);
    H5Pclose(fileInfo->groupPropertyList);
    H5Fclose(fileInfo->file);
    free(fileInfo);
}

HelenFeatureHDF5FileInfo** openHelenFeatureHDF5FilesByThreadCount(char *filenameBase, int64_t threadCount) {
    HelenFeatureHDF5FileInfo** infoArray = st_calloc(threadCount, sizeof(HelenFeatureHDF5FileInfo*));
    for (int64_t i = 0; i < threadCount; i++) {
        char *filename = stString_print("%s.T%02"PRId64".h5", filenameBase, i);
        infoArray[i] = HelenFeatureHDF5FileInfo_construct(filename);
        free(filename);
    }
    return infoArray;
}

#endif