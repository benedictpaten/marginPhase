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
typedef struct _stReferencePriorProbs stReferencePriorProbs;
typedef struct _stBaseMapper stBaseMapper;
typedef struct _stGenotypeResults stGenotypeResults;
typedef struct _stReferencePositionFilter stReferencePositionFilter;

/*
 * Overall coordination functions
 */


stList *filterReadsByCoverageDepth(stList *profileSeqs, stRPHmmParameters *params, stList *filteredProfileSeqs,
        stList *discardedProfileSeqs, stHash *referenceNamesToReferencePriors);

stList *getRPHmms(stList *profileSeqs, stHash *referenceNamesToReferencePriors, stRPHmmParameters *params);

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

uint16_t *getSubstitutionProb(uint16_t *matrix, int64_t from, int64_t to);

double *getSubstitutionProbSlow(double *matrix, int64_t from, int64_t to);

/*
 * Binary partition stuff
 */

// The maximum read depth the model can support
#define MAX_READ_PARTITIONING_DEPTH 64

char * intToBinaryString(uint64_t i);

uint64_t makeAcceptMask(uint64_t depth);

uint64_t mergePartitionsOrMasks(uint64_t partition1, uint64_t partition2,
        uint64_t depthOfPartition1, uint64_t depthOfPartition2);

uint64_t maskPartition(uint64_t partition, uint64_t mask);

bool seqInHap1(uint64_t partition, int64_t seqIndex);

uint64_t invertPartition(uint64_t partition, uint64_t depth);

/*
 * Profile sequence
 */

#define FIRST_ALPHABET_CHAR 48 // Ascii symbol '0'

struct _stProfileSeq {
    char *referenceName;
    char *readId;
    int64_t refStart;
    int64_t length;
    // The probability of alphabet characters, as specified by uint8_t
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

int stRPProfileSeq_cmpFn(const void *a, const void *b);

/*
 * Prior over reference positions
 */

struct _stReferencePriorProbs {
    char *referenceName;
    int64_t refStart;
    int64_t length;
    // The log probability of alphabet characters, as specified by uint16_t
    // see scaleToLogIntegerSubMatrix()
    // and invertScaleToLogIntegerSubMatrix() to see how probabilities are stored
    uint16_t *profileProbs;
    // The reference sequence
    uint8_t *referenceSequence;
    // Read counts for the bases seen in reads
    double *baseCounts;
    // Filter array of positions in the reference, used
    // to ignore some columns in the alignment
    bool *referencePositionsIncluded;
};

stReferencePriorProbs *stReferencePriorProbs_constructEmptyProfile(char *referenceName, int64_t referenceStart, int64_t length);

void stReferencePriorProbs_destruct(stReferencePriorProbs *seq);

stHash *createEmptyReferencePriorProbabilities(stList *profileSequences);

stHash *createReferencePriorProbabilities(char *referenceFastaFile, stList *profileSequences,
        stBaseMapper *baseMapper, stRPHmmParameters *params);

int64_t filterHomozygousReferencePositions(stHash *referenceNamesToReferencePriors, stRPHmmParameters *params, int64_t *totalPositions);

double *stReferencePriorProbs_estimateReadErrorProbs(stHash *referenceNamesToReferencePriors, stRPHmmParameters *params);

/*
 * Emission probabilities
 */

double emissionLogProbability(stRPColumn *column, stRPCell *cell, uint64_t *bitCountVectors,
                                stReferencePriorProbs *referencePriorProbs,
                                stRPHmmParameters *params);

double emissionLogProbabilitySlow(stRPColumn *column,
        stRPCell *cell, uint64_t *bitCountVectors, stReferencePriorProbs *referencePriorProbs,
        stRPHmmParameters *params, bool maxNotSum);

void fillInPredictedGenome(stGenomeFragment *gF, stRPCell *cell,
        stRPColumn *column, stReferencePriorProbs *referencePriorProbs, stRPHmmParameters *params);

// Constituent functions tested and used to do bit twiddling

int popcount64(uint64_t x);

uint64_t getExpectedInstanceNumber(uint64_t *bitCountVectors, uint64_t depth, uint64_t partition,
        int64_t position, int64_t characterIndex);

uint64_t *calculateCountBitVectors(uint8_t **seqs, int64_t depth, int64_t *activePositions, int64_t totalActivePositions);

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

    // Filters on the number of states in a column
    // Used to prune the hmm
    int64_t minPartitionsInAColumn;
    int64_t maxPartitionsInAColumn;
    double minPosteriorProbabilityForPartition;

    // MaxCoverageDepth is the maximum depth of profileSeqs to allow at any base. If the coverage depth is higher
    // than this then some profile seqs are randomly discarded.
    int64_t maxCoverageDepth;
    int64_t minReadCoverageToSupportPhasingBetweenHeterozygousSites;
    // Training

    // Pseudo counts used to make training of substitution matrices a bit more robust
    double offDiagonalReadErrorPseudoCount;
    double onDiagonalReadErrorPseudoCount;
    // Before doing any training estimate the read error substitution parameters empirically
    bool estimateReadErrorProbsEmpirically;
    // Number of iterations of training
    int64_t trainingIterations;
    // Whether or not to filter out poorly matching reads after one round and try again
    bool filterBadReads;
    double filterMatchThreshold;

    // Use a prior for the reference sequence
    bool useReferencePrior;

    // Verbosity options for printing
    bool verboseTruePositives;
    bool verboseFalsePositives;

    // Ensure symmetry in the HMM such that the inverted partition of each partition is included in the HMM
    bool includeInvertedPartitions;

    // Options to filter which positions in the reference sequence are included in the computation
    bool filterLikelyHomozygousSites;
    double minSecondMostFrequentBaseFilter; // See stReferencePriorProbs_setReferencePositionFilter
    double minSecondMostFrequentBaseLogProbFilter; //  See stReferencePriorProbs_setReferencePositionFilter

    // Whether or not to make deletions gap characters (otherwise, profile probs will be flat)
    bool gapCharactersForDeletions;

    // Any read that has one of the following sam flags is ignored when parsing the reads from the SAM/BAM file.
    // This allows the ability to optionally ignore, for example, secondary alignments.
    uint16_t filterAReadWithAnyOneOfTheseSamFlagsSet;

    // Whether or not to do the vcf comparison within marginPhase
    bool compareVCFs;
};

void stRPHmmParameters_destruct(stRPHmmParameters *params);

void stRPHmmParameters_learnParameters(stRPHmmParameters *params, stList *profileSequences,
        stHash *referenceNamesToReferencePriors);

void stRPHmmParameters_printParameters(stRPHmmParameters *params, FILE *fH);

void stRPHmmParameters_setReadErrorSubstitutionParameters(stRPHmmParameters *params, double *readErrorSubModel);

void normaliseSubstitutionMatrix(double *subMatrix);

double *getEmptyReadErrorSubstitutionMatrix(stRPHmmParameters *params);

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
    // Prior over reference bases
    stReferencePriorProbs *referencePriorProbs;
    // Filter used to mask column positions from consideration
    stReferencePositionFilter *referencePositionFilter;
};

stRPHmm *stRPHmm_construct(stProfileSeq *profileSeq, stReferencePriorProbs *referencePriorProbs, stRPHmmParameters *params);

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

void printBaseComposition2(double *baseCounts);
double *getColumnBaseComposition(stRPColumn *column, int64_t pos);
void printColumnAtPosition(stRPHmm *hmm, int64_t pos);
double *getProfileSequenceBaseCompositionAtPosition(stSet *profileSeqs, int64_t pos);

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
    // Record of which positions in the column are not filtered out
    int64_t *activePositions; // List of positions that are not filtered out, relative to the start of the column in reference coordinates
    int64_t totalActivePositions; // The length of activePositions
};

stRPColumn *stRPColumn_construct(int64_t refStart, int64_t length, int64_t depth,
        stProfileSeq **seqHeaders, uint8_t **seqs, stReferencePriorProbs *rProbs);

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
    uint8_t *referenceSequence;
};

stGenomeFragment *stGenomeFragment_construct(stRPHmm *hmm, stList *path);
void stGenomeFragment_destruct(stGenomeFragment *genomeFragment);

// Struct for alphabet and mapping bases to numbers
struct _stBaseMapper {
    uint8_t *charToNum;
    char *numToChar;
    char *wildcard;
    uint8_t size;
};

stBaseMapper* stBaseMapper_construct();
void stBaseMapper_destruct(stBaseMapper *bm);
void stBaseMapper_addBases(stBaseMapper *bm, char *bases);
void stBaseMapper_setWildcard(stBaseMapper* bm, char *wildcard);
char stBaseMapper_getCharForValue(stBaseMapper *bm, int value);
uint8_t stBaseMapper_getValueForChar(stBaseMapper *bm, char base);

// Parsing stuff
stRPHmmParameters *parseParameters(char *paramsFile, stBaseMapper *baseMapper);
int64_t parseReads(stList *profileSequences, char *bamFile, stBaseMapper *baseMapper, stRPHmmParameters *params);
void countIndels(uint32_t *cigar, uint32_t ncigar, int64_t *numInsertions, int64_t *numDeletions);
// Verbosity for what's printed.  To add more verbose options, you need to update:
//  usage, setVerbosity, struct _stRPHmmParameters, stRPHmmParameters_printParameters, writeParamFile
#define LOG_TRUE_POSITIVES 1
#define LOG_FALSE_POSITIVES 2
void setVerbosity(stRPHmmParameters *params, int64_t bitstring);

// File writing
void writeVcfFragment(vcfFile *out, bcf_hdr_t *bcf_hdr, stGenomeFragment *gF, char *referenceName, stBaseMapper *baseMapper, bool differencesOnly);
bcf_hdr_t* writeVcfHeader(vcfFile *out, stList *genomeFragments, char *referenceName);
void writeParamFile(char *outputFilename, stRPHmmParameters *params);

/*
 * Stores information about relevant test results.
 */
struct _stGenotypeResults {
    // Variants in reference
    int64_t negatives;
    int64_t positives;
    int64_t indelsInRef;
    int64_t homozygousVariantsInRef;

    // Variants in evaluated vcf
    int64_t truePositives;
    int64_t falsePositives;
    int64_t trueNegatives;
    int64_t falseNegatives;

    // Stats for specific types of variants
    int64_t truePositiveGaps;
    int64_t falsePositiveGaps;
    int64_t falseNegativeGaps;
    int64_t truePositiveHomozygous;

    // Types of errors
    int64_t error_badPartition;
    int64_t error_homozygousInRef;
    int64_t error_incorrectVariant;

    // Phasing
    int64_t switchErrors;
    float switchErrorDistance;
    int64_t uncertainPhasing;
};

void compareVCFs(FILE *fh, stList *hmms, char *vcf_toEval, char *vcf_ref,
                 stBaseMapper *baseMapper, stGenotypeResults *results, stRPHmmParameters *params);
void printGenotypeResults(stGenotypeResults *results);


void writeSplitBams(char *bamInFile, char *bamOutBase, stSet *haplotype1Ids, stSet *haplotype2Ids);
void addProfileSeqIdsToSet(stSet *pSeqs, stSet *readIds);


#endif /* ST_RP_HMM_H_ */
