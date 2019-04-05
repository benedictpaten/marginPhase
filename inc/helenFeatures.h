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
    HFEAT_SIMPLE_WEIGHT=1,
} HelenFeatureType;
#define POAFEATURE_SYMBOL_GAP_POS SYMBOL_NUMBER
#define POAFEATURE_SIMPLE_WEIGHT_TOTAL_SIZE ((SYMBOL_NUMBER + 1) * 2) // {A, C, G, T, N, gap} x {fwd, bkwd}
typedef struct _poaFeatureSimpleWeight PoaFeatureSimpleWeight;
struct _poaFeatureSimpleWeight {
    int64_t refPosition;
    int64_t insertPosition;
    char label;
    double weights[POAFEATURE_SIMPLE_WEIGHT_TOTAL_SIZE];
    PoaFeatureSimpleWeight* nextInsert; //so we can model all inserts after a position
};

PoaFeatureSimpleWeight *PoaFeature_SimpleWeight_construct(int64_t refPos, int64_t insPos);
void PoaFeature_SimpleCharacterCount_destruct(PoaFeatureSimpleWeight *scc);
stList *poa_getSimpleWeightFeatures(Poa *poa, stList *bamChunkReads);
void poa_writeHelenFeatures(HelenFeatureType type, Poa *poa, stList *bamChunkReads, char *outputFileBase,
                            BamChunk *bamChunk, stList *trueRefAlignment, RleString *trueRefRleString);

void writeSimpleWeightHelenFeaturesTSV(char *outputFileBase, BamChunk *bamChunk, bool outputLabels, stList *features,
                                       HelenFeatureType type, int64_t featureStartIdx, int64_t featureEndIdxInclusive);
int writeSimpleWeightHelenFeaturesHDF5(char *outputFileBase, BamChunk *bamChunk, bool outputLabels, stList *features,
                                       HelenFeatureType type, int64_t featureStartIdx, int64_t featureEndIdxInclusive);

#endif //MARGINPHASE_HELENFEATURES_H
