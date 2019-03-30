//
// Created by tpesout on 3/29/19.
//

#ifndef MARGINPHASE_HELENFEATURES_H
#define MARGINPHASE_HELENFEATURES_H

#include "margin.h"

/*
 * Simple feature (no run lengths)
 */
typedef enum {
    HFEAT_NONE=0,
    HFEAT_SIMPLE_COUNT=1,
    HFEAT_SIMPLE_WEIGHT=2
} HelenFeatureType;
#define POAFEATURE_SYMBOL_GAP_POS SYMBOL_NUMBER
#define POAFEATURE_TOTAL_SIZE ((SYMBOL_NUMBER + 1) * 2) // {A, C, G, T, N, gap} x {fwd, bkwd}
typedef struct _poaFeatureSimpleCharacterCount PoaFeatureSimpleCharacterCount;
struct _poaFeatureSimpleCharacterCount {
    int64_t refPosition;
    int64_t insertPosition;
    char label;
    int64_t counts[POAFEATURE_TOTAL_SIZE];
    double weights[POAFEATURE_TOTAL_SIZE];
    PoaFeatureSimpleCharacterCount* nextInsert; //so we can model all inserts after a position
};

PoaFeatureSimpleCharacterCount *PoaFeature_SimpleCharacterCount_construct(int64_t refPos, int64_t insPos);
void PoaFeature_SimpleCharacterCount_destruct(PoaFeatureSimpleCharacterCount *scc);
int64_t PoaFeature_SimpleCharacterCount_getTotalCount(PoaFeatureSimpleCharacterCount *scc);
double PoaFeature_SimpleCharacterCount_getTotalWeight(PoaFeatureSimpleCharacterCount *scc);
stList *poa_getSimpleCharacterCountFeatures(Poa *poa, stList *bamChunkReads);
void poa_writeHelenFeatures(HelenFeatureType type, Poa *poa, stList *bamChunkReads, char *outputFileBase,
                            BamChunk *bamChunk, stList *trueRefAlignment, RleString *trueRefRleString);

void writeSimpleHelenFeaturesTSV(char *outputFile, BamChunk *bamChunk, bool outputLabels, stList *features,
                                 HelenFeatureType type);
int writeSimpleHelenFeaturesHDF5(char *outputFileBase, BamChunk *bamChunk, bool outputLabels, stList *features,
                                 HelenFeatureType type);

#endif //MARGINPHASE_HELENFEATURES_H
