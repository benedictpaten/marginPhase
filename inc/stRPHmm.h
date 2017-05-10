/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef ST_RP_HMM_H_
#define ST_RP_HMM_H_

#include <stdint.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <float.h>

#include <htslib/vcf.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>

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
typedef struct _stGenomeFragment stGenomeFragment;

/*
 * Overall coordination functions
 */

stList *filterProfileSeqsToMaxCoverageDepth(stList *profileSeqs, int64_t maxDepth);

stList *getRPHmms(stList *profileSeqs, stRPHmmParameters *params);

stList *getTilingPaths(stSortedSet *hmms);

stSet *getOverlappingComponents(stList *tilingPath1, stList *tilingPath2);

/*
 * Math
 */
#define ST_MATH_LOG_ZERO -INFINITY
#define ST_MATH_LOG_ONE 0.0

double logAddP(double a, double b, bool maxNotSum);

/*
 * Alphabet
 */

#define ALPHABET_SIZE 5
#define ALPHABET_MAX_PROB 255
#define ALPHABET_MIN_PROB 0
#define ALPHABET_CHARACTER_BITS 8
#define ALPHABET_MIN_SUBSTITUTION_PROB 65535 // 2^16 -1

// Each value is expressed as an unsigned integer scaled linearly from 0 to 2^16-1, with 0 = log(1) and 2^16-1 = -7 = log(0.0000001)
uint16_t scaleToLogIntegerSubMatrix(double logProb);
double invertScaleToLogIntegerSubMatrix(int64_t i);

void setSubstitutionProb(uint16_t *logSubMatrix, double *logSubMatrixSlow,
        int64_t sourceCharacterIndex,
        int64_t derivedCharacterIndex, double prob);

/*
 * Binary partition stuff
 */

// The maximum read depth the model can support
#define MAX_READ_PARTITIONING_DEPTH 64

char * intToBinaryString(uint64_t i);

uint64_t makeAcceptMask(int64_t depth);

uint64_t mergePartitionsOrMasks(uint64_t partition1, uint64_t partition2,
        uint64_t depthOfPartition1, uint64_t depthOfPartition2);

uint64_t maskPartition(uint64_t partition, uint64_t mask);

bool seqInHap1(uint64_t partition, int64_t seqIndex);

/*
 * Profile sequence
 */

#define FIRST_ALPHABET_CHAR 48 // Ascii symbol '0'

struct _stProfileSeq {
    char *referenceName;
    char *readId;
    int64_t refStart;
    int64_t length;
    // The probability of alphabet characters, as specified by uint16_t
    // Each is expressed as an 8 bit unsigned int, with 0x00 representing 0 prob and
    // 0xFF representing 1.0 and each step between representing a linear step in probability of
    // 1.0/255
    uint8_t *profileProbs;
};

stProfileSeq *stProfileSeq_constructEmptyProfile(char *referenceName, char *readId,
                                                 int64_t referenceStart, int64_t length);

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

double emissionLogProbabilitySlow(stRPColumn *column,
        stRPCell *cell, uint64_t *bitCountVectors,
        stRPHmmParameters *params, bool maxNotSum);

void fillInPredictedGenome(stGenomeFragment *gF, stRPCell *cell,
        stRPColumn *column, stRPHmmParameters *params);

// Constituent functions tested and used to do bit twiddling

int popcount64(uint64_t x);

uint64_t getExpectedInstanceNumber(uint64_t *bitCountVectors, uint64_t depth, uint64_t partition,
        int64_t position, int64_t characterIndex);

uint64_t *calculateCountBitVectors(uint8_t **seqs, int64_t depth, int64_t length);

/*
 * Read partitioning hmm
 */

struct _stRPHmmParameters {
    /*
     * Parameters used for the HMM computation
     */
    uint16_t *hetSubModel;
    uint16_t *readErrorSubModel;
    double *hetSubModelSlow;
    double *readErrorSubModelSlow;
    bool maxNotSumTransitions;
    int64_t maxPartitionsInAColumn;
    int64_t maxCoverageDepth;
    int64_t minReadCoverageToSupportPhasingBetweenHeterozygousSites;
};

stRPHmmParameters *stRPHmmParameters_construct(uint16_t *hetSubModel,
        double *hetSubModelSlow,
        uint16_t *readErrorSubModel,
        double *readErrorSubModelSlow,
        bool maxNotSumTransitions,
        int64_t maxPartitionsInAColumn,
        int64_t maxCoverageDepth,
        int64_t minReadCoverageToSupportPhasingBetweenHeterozygousSites);

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

stRPHmm *stRPHmm_construct(stProfileSeq *profileSeq, stRPHmmParameters *params);

void stRPHmm_destruct(stRPHmm *hmm, bool destructColumns);

void stRPHmm_destruct2(stRPHmm *hmm);

bool stRPHmm_overlapOnReference(stRPHmm *hmm1, stRPHmm *hmm2);

stRPHmm *stRPHmm_createCrossProductOfTwoAlignedHmm(stRPHmm *hmm1, stRPHmm *hmm2);

void stRPHmm_alignColumns(stRPHmm *hmm1, stRPHmm *hmm2);

stRPHmm *stRPHmm_fuse(stRPHmm *leftHmm, stRPHmm *rightHmm);

void stRPHmm_forwardBackward(stRPHmm *hmm);

void stRPHmm_prune(stRPHmm *hmm);

void stRPHmm_print(stRPHmm *hmm, FILE *fileHandle, bool includeColumns, bool includeCells);

stList *stRPHmm_forwardTraceBack(stRPHmm *hmm);

stSet *stRPHmm_partitionSequencesByStatePath(stRPHmm *hmm, stList *path, bool partition1);

int stRPHmm_cmpFn(const void *a, const void *b);

stRPHmm *stRPHmm_split(stRPHmm *hmm, int64_t splitPoint);

void stRPHmm_resetColumnNumberAndDepth(stRPHmm *hmm);

stList *stRPHMM_splitWherePhasingIsUncertain(stRPHmm *hmm);

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

stSet *stRPColumn_getSequencesInCommon(stRPColumn *column1, stRPColumn *column2);

stSet *stRPColumn_getColumnSequencesAsSet(stRPColumn *column);

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

/*
 * String to represent genotype and haplotype inference from an HMM
 */

struct _stGenomeFragment {
    // A string where each element represents the predicted genotype at the corresponding
    // position.
    // A genotype is represented by an integer in the range [0, ALPHABET_SIZE**2)
    // A genotype expresses two characters. For two characters x, y represented by two integers
    // in [0, ALPHABET_SIZE) then the genotype is expressed as x * ALPHABET_SIZE + y if x <= y
    // else y * ALPHABET_SIZE + x
    uint64_t *genotypeString;

    // An array of genotype posterior probabilities,
    // each between 0 and 1, for the corresponding genotypes
    // in the genotype string
    float *genotypeProbs;

    // Strings representing the predicted haplotypes, where each element is an alphabet character
    // index in [0, ALPHABET_SIZE)
    uint64_t *haplotypeString1;
    uint64_t *haplotypeString2;

    // An array of haplotype posterior probabilities,
    // each between 0 and 1, for the corresponding haplotypes
    // in the haplotype strings
    float *haplotypeProbs1;
    float *haplotypeProbs2;

    // The reference coordinates of the genotypes
    char *referenceName;
    int64_t refStart;
    int64_t length;
};

stGenomeFragment *stGenomeFragment_construct(stRPHmm *hmm, stList *path);
void stGenomeFragment_destruct(stGenomeFragment *genomeFragment);

// Struct for alphabet and mapping bases to numbers
struct _stBaseMapper {
    char *baseToNum;
    int *numToBase;
    char *wildcard;
    int size;
};
typedef struct _stBaseMapper stBaseMapper;
stBaseMapper* stBaseMapper_construct();
void stBaseMapper_addBases(stBaseMapper *bm, char *bases);
void stBaseMapper_setWildcard(stBaseMapper* bm, char *wildcard);
int stBaseMapper_getBaseForValue(stBaseMapper *bm, int value);
int stBaseMapper_getValueForBase(stBaseMapper *bm, char base);

// Parsing stuff
stRPHmmParameters *parseParameters(char *paramsFile, stBaseMapper *baseMapper);
void parseReads(stList *profileSequences, char *bamFile, stBaseMapper *baseMapper);

// File writing
void writeVcfFragment(vcfFile *out, bcf_hdr_t *bcf_hdr, stGenomeFragment *gF, char *referenceName, stBaseMapper *baseMapper, bool includeReference);
bcf_hdr_t* writeVcfHeader(vcfFile *out, stList *genomeFragments);
bcf_hdr_t* writeSplitSams(char *bamInFile, char *bamOutBase, stSet *haplotype1Ids, stSet *haplotype2Ids);


#endif /* ST_RP_HMM_H_ */
