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

/*
 * Profile sequence
 */

struct _stProfileSeq {
    char *referenceName;
    int64_t refStart;
    int64_t length;
    stProfileProb *profileProbs;
};

stProfileSeq *stProfileSeq_constructEmptyProfile(int64_t length);

void stProfileSeq_destruct(stProfileSeq *seq);

/*
 * Element of a profile sequence
 */

struct _stProfileProb {
    double probA;
    double probC;
    double probG;
    double probT;
};

/*
 * Read partitioning hmm
 */

struct _stRPHmm {
    char *referenceName;
    int64_t refStart;
    int64_t refLength;
    stList *profileSeqs; // List of stProfileSeq
    int64_t columnNumber;
    int64_t maxDepth;
    stRPColumn *firstColumn;
    stRPColumn *lastColumn;

    //Forward/backward probability calculation things
    double forwardProbability;
    double backwardProbability;
    double (*emissionProbability)(stRPColumn *, stRPCell *);
};

stList *getRPHmms(stList *profileSeqs, double posteriorProbabilityThreshold);

stRPHmm *stRPHmm_construct(stProfileSeq *profileSeq);

void stRPHmm_destruct(stRPHmm *hmm);

bool stRPHmm_overlapOnReference(stRPHmm *hmm1, stRPHmm *hmm2);

stRPHmm *stRPHmm_createCrossProductHmm(stRPHmm *hmm1, stRPHmm *hmm2);

void stRPHmm_alignColumns(stRPHmm *hmm1, stRPHmm *hmm2);

stRPHmm *stRPHmm_fuse(stRPHmm *leftHmm, stRPHmm *rightHmm);

void stRPHmm_forward(stRPHmm *hmm);

void stRPHmm_backward(stRPHmm *hmm);

void stRPHmm_prune(stRPHmm *hmm, double posteriorProbabilityThreshold);

/*
 * Column of read partitioning hmm
 */

struct _stRPColumn {
    int64_t refStart;
    int64_t length;
    int64_t depth;
    stProfileProb **seqs;
    stRPCell *head;
    stRPMergeColumn *nColumn, *pColumn;
    double forwardProb, backwardProb;
};

stRPColumn *stRPColumn_construct(int64_t refStart, int64_t length, int64_t depth, stProfileProb **seqs);

void stRPColumn_destruct(stRPColumn *column);

void stRPColumn_split(stRPColumn *column, int64_t firstHalfLength, stRPHmm *hmm);

/*
 * State of read partitioning hmm
 */

struct _stRPCell {
    int32_t partition;
    double forwardProb, backwardProb;
    stRPCell *nCell;
};

stRPCell *stRPCell_construct(int64_t partition);

void stRPCell_destruct(stRPCell *cell, bool);

double stRPCell_posteriorProb(stRPCell *cell, stRPColumn *column);

/*
 * Merge column of read partitioning hmm
 */

struct _stRPMergeColumn {
    int32_t maskFrom;
    int32_t maskTo;
    stHash *mergeCellsFrom;
    stHash *mergeCellsTo;
    stRPColumn *nColumn, *pColumn;
};

stRPMergeCell *stRPMergeColumn_getNextMergeCell(stRPCell *cell, stRPMergeColumn *mergeColumn);

stRPMergeCell *stRPMergeColumn_getPreviousMergeCell(stRPCell *cell, stRPMergeColumn *mergeColumn);

stRPMergeColumn *stRPMergeColumn_construct(int32_t maskFrom, int32_t maskTo);

void stRPMergeColumn_destruct(stRPMergeColumn *mColumn);

void stRPMergeColumn_addCell(stRPMergeColumn *mColumn, int32_t fromPartition, int32_t toPartition);

double stRPMergeCell_posteriorProb(stRPMergeCell *mCell, stRPMergeColumn *mColumn);

/*
 * Merge cell of read partitioning hmm
 */

struct _stRPMergeCell {
    int32_t fromPartition;
    int32_t toPartition;
    double forwardProb, backwardProb;
};

void stRPMergeCell_destruct(stRPMergeCell *mCell);

#endif /* ST_RP_HMM_H_ */
