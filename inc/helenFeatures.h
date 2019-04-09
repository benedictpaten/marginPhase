//
// Created by tpesout on 3/29/19.
//

#ifndef MARGINPHASE_HELENFEATURES_H
#define MARGINPHASE_HELENFEATURES_H

#include "margin.h"

typedef enum {
    HFEAT_NONE=0,
    HFEAT_SIMPLE_WEIGHT=1,
    HFEAT_RLE_WEIGHT=2,
} HelenFeatureType;

#define POAFEATURE_SYMBOL_GAP_POS SYMBOL_NUMBER
#define POAFEATURE_SIMPLE_WEIGHT_TOTAL_SIZE ((SYMBOL_NUMBER + 1) * 2) // {A,C,G,T,N,gap} x {fwd,bkwd}
#define POAFEATURE_MAX_RUN_LENGTH 20
#define POAFEATURE_RLE_WEIGHT_TOTAL_SIZE ((SYMBOL_NUMBER * POAFEATURE_MAX_RUN_LENGTH + 1) * 2 ) // ({A,C,G,T,N} x {rlesize} + {gap}) x {fwd,bkwd}

typedef struct _poaFeatureSimpleWeight PoaFeatureSimpleWeight;
struct _poaFeatureSimpleWeight {
    int64_t refPosition;
    int64_t insertPosition;
    char label;
    double weights[POAFEATURE_SIMPLE_WEIGHT_TOTAL_SIZE];
    PoaFeatureSimpleWeight* nextInsert; //so we can model all inserts after a position
};

typedef struct _poaFeatureRleWeight PoaFeatureRleWeight;
struct _poaFeatureRleWeight {
    int64_t refPosition;
    int64_t insertPosition;
    char labelChar;
    int64_t labelRunLength;
    double weights[POAFEATURE_RLE_WEIGHT_TOTAL_SIZE];
    PoaFeatureRleWeight* nextInsert; //so we can model all inserts after a position
};

int PoaFeature_SimpleWeight_charIndex(Symbol character, bool forward);
int PoaFeature_SimpleWeight_gapIndex(bool forward);
int PoaFeature_RleWeight_charIndex(Symbol character, int64_t runLength, bool forward);
int PoaFeature_RleWeight_gapIndex(bool forward);

PoaFeatureSimpleWeight *PoaFeature_SimpleWeight_construct(int64_t refPos, int64_t insPos);
void PoaFeature_SimpleWeight_destruct(PoaFeatureSimpleWeight *feature);

PoaFeatureRleWeight *PoaFeature_RleWeight_construct(int64_t refPos, int64_t insPos);
void PoaFeature_RleWeight_destruct(PoaFeatureRleWeight *feature);

stList *poa_getSimpleWeightFeatures(Poa *poa, stList *bamChunkReads);
stList *poa_getWeightedRleFeatures(Poa *poa, stList *bamChunkReads, stList *rleStrings);

void poa_writeHelenFeatures(HelenFeatureType type, Poa *poa, stList *bamChunkReads, stList *rleStrings,
        char *outputFileBase, BamChunk *bamChunk, stList *trueRefAlignment, RleString *trueRefRleString,
        bool fullFeatureOutput);

void writeSimpleWeightHelenFeaturesTSV(char *outputFileBase, BamChunk *bamChunk, bool outputLabels, stList *features,
                                       int64_t featureStartIdx, int64_t featureEndIdxInclusive);

void writeSimpleWeightHelenFeaturesHDF5(char *outputFileBase, BamChunk *bamChunk, bool outputLabels, stList *features,
                                        int64_t featureStartIdx, int64_t featureEndIdxInclusive);

void writeRleWeightHelenFeaturesTSV(char *outputFileBase, BamChunk *bamChunk, bool outputLabels, stList *features,
                                    int64_t featureStartIdx, int64_t featureEndIdxInclusive);

void writeRleWeightHelenFeaturesHDF5(char *outputFileBase, BamChunk *bamChunk, bool outputLabels, stList *features,
                                     int64_t featureStartIdx, int64_t featureEndIdxInclusive);

#endif //MARGINPHASE_HELENFEATURES_H
