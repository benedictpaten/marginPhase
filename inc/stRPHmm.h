/*
 * model.h
 *
 *  Created on: Feb 4, 2017
 *      Author: benedictpaten
 */

#ifndef ST_RP_HMM_H_
#define ST_RP_HMM_H_

#include <stdint.h>
#include <stdlib.h>
#include <stddef.h>
#include "sonLib.h"

/*
 * Function documentation is in the .c file
 */

typedef struct _stProfileSeq stProfileSeq;
typedef struct _stRPHmm stRPHmm;
typedef struct _stRPHmmParameters stRPHmmParameters;
typedef struct _stRPColumn stRPColumn;
typedef struct _stRPCell stRPCell;
typedef struct _stRPMergeColumn stRPMergeColumn;
typedef struct _stRPMergeCell stRPMergeCell;
typedef struct _stSubModel stSubModel;

// The maximum read depth the model can support
#define MAX_READ_PARTITIONING_DEPTH 64

/*
 * Alphabet
 */

#define ALPHABET_MAX_PROB 255
#define ALPHABET_MIN_PROB 0
#define ALPHABET_CHARACTER_BITS 8

struct _stSubModel {
    int64_t alphabetSize;
    double *logSubMatrix;
};

stSubModel *stSubModel_constructEmptyModel(int64_t alphabetSize);

void stSubModel_destruct(stSubModel *alphabet);

double stSubModel_getSubstitutionProb(stSubModel *alphabet, int64_t sourceCharacterIndex,
        int64_t derivedCharacterIndex);

void stSubModel_setSubstitutionProb(stSubModel *alphabet, int64_t sourceCharacterIndex,
        int64_t derivedCharacterIndex, double prob);

char * intToBinaryString(uint64_t i);

/*
 * Profile sequence
 */

#define FIRST_ALPHABET_CHAR 48 // Ascii symbol '0'

struct _stProfileSeq {
    char *referenceName;
    int64_t refStart;
    int64_t length;
    int64_t alphabetSize;
    // The probability of alphabet characters, as specified by stSubModel
    // Each is expressed as an 8 bit unsigned int, with 0x00 representing 0 prob and
    // 0xFF representing 1.0 and each step between representing a linear step in probability of
    // 1.0/255
    uint8_t *profileProbs;
};

stProfileSeq *stProfileSeq_constructEmptyProfile(char *referenceName,
        int64_t referenceStart, int64_t length, int64_t alphabetSize);

void stProfileSeq_destruct(stProfileSeq *seq);

void stProfileSeq_print(stProfileSeq *seq, FILE *fileHandle, bool includeProbs);

float getProb(uint8_t *p, int64_t characterIndex);

void printSeqs(FILE *fileHandle, stSet *profileSeqs);

void printPartition(FILE *fileHandle, stSet *profileSeqs1, stSet *profileSeqs2);

/*
 * Emission probabilities
 */

double emissionLogProbability(stRPColumn *column, stRPCell *cell, uint64_t *bitCountVectors,
                                stRPHmmParameters *params);

// Constituent functions tested and used to do bit twiddling

int popcount64(uint64_t x);

double getExpectedInstanceNumber(uint64_t *bitCountVectors, uint64_t depth, uint64_t partition,
        int64_t position, int64_t characterIndex, int64_t alphabetSize);

uint64_t *calculateCountBitVectors(uint8_t **seqs, int64_t depth, int64_t length, int64_t alphabetSize);

/*
 * Read partitioning hmm
 */

struct _stRPHmmParameters {
    /*
     * Parameters used for the HMM computation
     */
    int64_t alphabetSize;
    stSubModel *hetSubModel;
    stSubModel *readErrorSubModel;
    bool maxNotSumEmissions;
    bool maxNotSumTransitions;
    int64_t maxPartitionsInAColumn;
    int64_t maxCoverageDepth;
};

stRPHmmParameters *stRPHmmParameters_construct(stSubModel *hetSubModel,
        stSubModel *readErrorSubModel,
        bool maxNotSumEmissions,
        bool maxNotSumTransitions,
        int64_t maxPartitionsInAColumn,
        int64_t maxCoverageDepth);

void stRPHmmParameters_destruct(stRPHmmParameters *params);

struct _stRPHmm {
    char *referenceName;
    int64_t refStart;
    int64_t refLength;
    stList *profileSeqs; // List of stProfileSeq
    int64_t columnNumber; // Number of columns, excluding merge columns
    int64_t maxDepth;
    stRPColumn *firstColumn;
    stRPColumn *lastColumn;
    const stRPHmmParameters *parameters;
    //Forward/backward probability calculation things
    double forwardLogProb;
    double backwardLogProb;
};

stList *filterProfileSeqsToMaxCoverageDepth(stList *profileSeqs, int64_t maxDepth);

stList *getRPHmms(stList *profileSeqs, stRPHmmParameters *params);

stList *getTilingPaths(stSortedSet *hmms);

stSet *getOverlappingComponents(stList *tilingPath1, stList *tilingPath2);

stRPHmm *stRPHmm_construct(stProfileSeq *profileSeq, stRPHmmParameters *params);

void stRPHmm_destruct(stRPHmm *hmm, bool destructColumns);

void stRPHmm_destruct2(stRPHmm *hmm);

bool stRPHmm_overlapOnReference(stRPHmm *hmm1, stRPHmm *hmm2);

stRPHmm *stRPHmm_createCrossProductOfTwoAlignedHmm(stRPHmm *hmm1, stRPHmm *hmm2);

void stRPHmm_alignColumns(stRPHmm *hmm1, stRPHmm *hmm2);

stRPHmm *stRPHmm_fuse(stRPHmm *leftHmm, stRPHmm *rightHmm);

void stRPHmm_forward(stRPHmm *hmm);

void stRPHmm_forwardBackward(stRPHmm *hmm);

void stRPHmm_prune(stRPHmm *hmm);

void stRPHmm_print(stRPHmm *hmm, FILE *fileHandle, bool includeColumns, bool includeCells);

stList *stRPHmm_forwardTraceBack(stRPHmm *hmm);

stSet *stRPHmm_partitionSequencesByStatePath(stRPHmm *hmm, stList *path, bool partition1);

int stRPHmm_cmpFn(const void *a, const void *b);

/*
 * Column of read partitioning hmm
 */

struct _stRPColumn {
    int64_t refStart;
    int64_t length;
    int64_t depth;
    stProfileSeq **seqHeaders;
    uint8_t **seqs;
    stRPCell *head;
    stRPMergeColumn *nColumn, *pColumn;
    double totalLogProb;
};

stRPColumn *stRPColumn_construct(int64_t refStart, int64_t length, int64_t depth,
        stProfileSeq **seqHeaders, uint8_t **seqs);

void stRPColumn_destruct(stRPColumn *column);

void stRPColumn_print(stRPColumn *column, FILE *fileHandle, bool includeCells);

void stRPColumn_split(stRPColumn *column, int64_t firstHalfLength, stRPHmm *hmm);

/*
 * State of read partitioning hmm
 */

struct _stRPCell {
    uint64_t partition;
    double forwardLogProb, backwardLogProb;
    stRPCell *nCell;
};

stRPCell *stRPCell_construct(int64_t partition);

void stRPCell_destruct(stRPCell *cell);

void stRPCell_print(stRPCell *cell, FILE *fileHandle);

double stRPCell_posteriorProb(stRPCell *cell, stRPColumn *column);

/*
 * Merge column of read partitioning hmm
 */

struct _stRPMergeColumn {
    uint64_t maskFrom;
    uint64_t maskTo;
    stHash *mergeCellsFrom;
    stHash *mergeCellsTo;
    stRPColumn *nColumn, *pColumn;
};

stRPMergeColumn *stRPMergeColumn_construct(uint64_t maskFrom, uint64_t maskTo);

void stRPMergeColumn_destruct(stRPMergeColumn *mColumn);

void stRPMergeColumn_print(stRPMergeColumn *mColumn, FILE *fileHandle, bool includeCells);

stRPMergeCell *stRPMergeColumn_getNextMergeCell(stRPCell *cell, stRPMergeColumn *mergeColumn);

stRPMergeCell *stRPMergeColumn_getPreviousMergeCell(stRPCell *cell, stRPMergeColumn *mergeColumn);

int64_t stRPMergeColumn_numberOfPartitions(stRPMergeColumn *mColumn);

/*
 * Merge cell of read partitioning hmm
 */

struct _stRPMergeCell {
    uint64_t fromPartition;
    uint64_t toPartition;
    double forwardLogProb, backwardLogProb;
};

stRPMergeCell *stRPMergeCell_construct(uint64_t fromPartition,
        uint64_t toPartition, stRPMergeColumn *mColumn);

void stRPMergeCell_destruct(stRPMergeCell *mCell);

void stRPMergeCell_print(stRPMergeCell *mCell, FILE *fileHandle);

double stRPMergeCell_posteriorProb(stRPMergeCell *mCell, stRPMergeColumn *mColumn);


#endif /* ST_RP_HMM_H_ */
