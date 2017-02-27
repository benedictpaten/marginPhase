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
typedef struct _stProfileProb stProfileProb;
typedef struct _stRPHmm stRPHmm;
typedef struct _stRPColumn stRPColumn;
typedef struct _stRPCell stRPCell;
typedef struct _stRPMergeColumn stRPMergeColumn;
typedef struct _stRPMergeCell stRPMergeCell;

// The maximum read depth the model can support
#define MAX_READ_PARTITIONING_DEPTH 64

/*
 * Profile sequence
 */

struct _stProfileSeq {
    char *referenceName;
    int64_t refStart;
    int64_t length;
    stProfileProb *profileProbs;
};

stProfileSeq *stProfileSeq_constructEmptyProfile(char *referenceName, int64_t referenceStart, int64_t length);

void stProfileSeq_destruct(stProfileSeq *seq);

void stProfileSeq_print(stProfileSeq *seq, FILE *fileHandle, bool includeSequence);


/*
 * Element of a profile sequence
 */

#define NUCLEOTIDE_ALPHABET_SIZE 8
#define NUCLEOTIDE_BITS sizeof(uint8_t)
#define NUCLEOTIDE_GAP 0
#define NUCLEOTIDE_A 1
#define NUCLEOTIDE_C 2
#define NUCLEOTIDE_G 3
#define NUCLEOTIDE_T 4
#define NUCLEOTIDE_METHYL_C 5
#define NUCLEOTIDE_HYDROXYMETHYL_C 6
#define NUCLEOTIDE_METHYL_A 7
#define NUCLEOTIDE_MAX_PROB 255
#define NUCLEOTIDE_MIN_PROB 0

struct _stProfileProb {
    // The probability of, in order, -, A, C, T, G, methyl-C, hydroxymethyl-C and methyl-A
    // Each is expressed as an 8 bit unsigned int, with 0x00 representing 0 prob and
    // 0xFF representing 1.0 and each step between representing a linear step in probability of
    // 1.0/255
    uint8_t probs[NUCLEOTIDE_ALPHABET_SIZE];
};

float stProfileProb_prob(stProfileProb *p, int64_t characterIndex);

/*
 * Emission probabilities
 */

double emissionLogProbability(stRPColumn *column, stRPCell *cell, uint64_t *bitCountVectors, double *logSubMatrix);

double *getLogSubstitutionMatrix();

/*
 * Read partitioning hmm
 */

struct _stRPHmm {
    char *referenceName;
    int64_t refStart;
    int64_t refLength;
    stList *profileSeqs; // List of stProfileSeq
    int64_t columnNumber; // Number of columns, excluding merge columns
    int64_t maxDepth;
    stRPColumn *firstColumn;
    stRPColumn *lastColumn;

    //Forward/backward probability calculation things
    double forwardLogProb;
    double backwardLogProb;
    double *logSubMatrix;
};

stList *filterProfileSeqsToMaxCoverageDepth(stList *profileSeqs, int64_t maxDepth);

stList *getRPHmms(stList *profileSeqs, double posteriorProbabilityThreshold,
        int64_t minColumnDepthToFilter, int64_t maxCoverageDepth, double *logSubMatrix);

stRPHmm *stRPHmm_construct(stProfileSeq *profileSeq,
        double *logSubMatrix);

void stRPHmm_destruct(stRPHmm *hmm);

bool stRPHmm_overlapOnReference(stRPHmm *hmm1, stRPHmm *hmm2);

stRPHmm *stRPHmm_createCrossProductOfTwoAlignedHmm(stRPHmm *hmm1, stRPHmm *hmm2);

void stRPHmm_alignColumns(stRPHmm *hmm1, stRPHmm *hmm2);

stRPHmm *stRPHmm_fuse(stRPHmm *leftHmm, stRPHmm *rightHmm);

void stRPHmm_forward(stRPHmm *hmm);

void stRPHmm_backward(stRPHmm *hmm);

void stRPHmm_prune(stRPHmm *hmm, double posteriorProbabilityThreshold, int64_t minColumnDepthToFilter);

void stRPHmm_print(stRPHmm *hmm, FILE *fileHandle, bool includeColumns, bool includeCells);

stList *stRPHmm_forwardTraceBack(stRPHmm *hmm);

stSet *stRPHmm_partitionSequencesByStatePath(stRPHmm *hmm, stList *path);

/*
 * Column of read partitioning hmm
 */

struct _stRPColumn {
    int64_t refStart;
    int64_t length;
    int64_t depth;
    stProfileSeq **seqHeaders;
    stProfileProb **seqs;
    stRPCell *head;
    stRPMergeColumn *nColumn, *pColumn;
    double forwardLogProb, backwardLogProb;
};

stRPColumn *stRPColumn_construct(int64_t refStart, int64_t length, int64_t depth,
        stProfileSeq **seqHeaders, stProfileProb **seqs);

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

int64_t stRPMergeColumn_depth(stRPMergeColumn *mColumn);

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
