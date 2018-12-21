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
#include <inttypes.h>
#include <ctype.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h>
#include <util.h>
#include <stdio.h>
#include <string.h>

#include <htslib/vcf.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/bgzf.h>
#include <htslib/hts.h>

#include "sonLib.h"
#include "hashTableC.h"
#include "pairwiseAligner.h"
#include "randomSequences.h"
#include "multipleAligner.h"

/*
 * More function documentation is in the .c files
 */

/*
 * Phasing structs
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
 * Polisher structs
 */
typedef struct _stReadHaplotypeSequence stReadHaplotypeSequence;
typedef struct hashtable stReadHaplotypePartitionTable;
typedef struct _repeatSubMatrix RepeatSubMatrix;
typedef struct _polishParams PolishParams;
typedef struct _Poa Poa;
typedef struct _poaNode PoaNode;
typedef struct _poaInsert PoaInsert;
typedef struct _poaDelete PoaDelete;
typedef struct _poaBaseObservation PoaBaseObservation;
typedef struct _rleString RleString;
typedef struct _refMsaView MsaView;
/*
 * Combined params object
 */
typedef struct _params Params;


/*
 * Combined parameter object for phase, polish, view, etc.
 */

struct _params {
	PolishParams *polishParams;
	stRPHmmParameters *phaseParams;
	stBaseMapper *baseMapper;
};

Params *params_readParams(FILE *fp);

void params_destruct(Params *params);

void params_printParameters(Params *params, FILE *fh);

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

// Each value is expressed as an unsigned integer scaled linearly from 0 to 2^16-1,
// with 0 = log(1) and 2^16-1 = -7 = log(0.0000001)
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

uint64_t flipAReadsPartition(uint64_t partition, uint64_t readIndex);

/*
 * _stProfileSeq
 * Struct for profile sequence
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

stProfileSeq *stProfileSeq_constructFromPosteriorProbs(char *refName, char *refSeq, int64_t refLength,
													   char *readId, char *readSeq, stList *anchorAlignment,
													   Params *params);

void stProfileSeq_destruct(stProfileSeq *seq);

void stProfileSeq_print(stProfileSeq *seq, FILE *fileHandle, bool includeProbs);

float getProb(uint8_t *p, int64_t characterIndex);

void printSeqs(FILE *fileHandle, stSet *profileSeqs);

void printPartition(FILE *fileHandle, stSet *profileSeqs1, stSet *profileSeqs2);

int stRPProfileSeq_cmpFn(const void *a, const void *b);



/*
 * _stReferencePriorProbs
 * Struct for prior over reference positions
 */
struct _stReferencePriorProbs {
    char *referenceName;
    int64_t refStart;
    int64_t length;

    // The log probability of alphabet characters, as specified by uint16_t
    // see scaleToLogIntegerSubMatrix()
    // and invertScaleToLogIntegerSubMatrix() to see how probabilities are stored
    uint16_t *profileProbs;
    uint8_t *referenceSequence; // The reference sequence
    // Read counts for the bases seen in reads
    double *baseCounts;
    // Filter array of positions in the reference, used to ignore some columns in the alignment
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
 * Emission probabilities methods
 */
double emissionLogProbability(stRPColumn *column, stRPCell *cell, uint64_t *bitCountVectors,
                                stReferencePriorProbs *referencePriorProbs,
                                stRPHmmParameters *params);

double emissionLogProbabilitySlow(stRPColumn *column,
        stRPCell *cell, uint64_t *bitCountVectors, stReferencePriorProbs *referencePriorProbs,
        stRPHmmParameters *params, bool maxNotSum);

void fillInPredictedGenome(stGenomeFragment *gF, uint64_t partition,
        stRPColumn *column, stReferencePriorProbs *referencePriorProbs, stRPHmmParameters *params);

/*
 * Constituent functions tested and used to do bit twiddling
*/
int popcount64(uint64_t x);

uint64_t getExpectedInstanceNumber(uint64_t *bitCountVectors, uint64_t depth, uint64_t partition,
        int64_t position, int64_t characterIndex);

uint64_t *calculateCountBitVectors(uint8_t **seqs, int64_t depth, int64_t *activePositions, int64_t totalActivePositions);

/*
 * _stRPHmmParameters
 * Struct for hmm parameters
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

    // MaxCoverageDepth is the maximum depth of profileSeqs to allow at any base.
    // If the coverage depth is higher than this then some profile seqs are randomly discarded.
    int64_t maxCoverageDepth;
    int64_t minReadCoverageToSupportPhasingBetweenHeterozygousSites;

    // Training

    // Number of iterations of training
    int64_t trainingIterations;
    // Pseudo counts used to make training of substitution matrices a bit more robust
    double offDiagonalReadErrorPseudoCount;
    double onDiagonalReadErrorPseudoCount;
    // Before doing any training estimate the read error substitution parameters empirically
    bool estimateReadErrorProbsEmpirically;

    // Whether or not to filter out poorly matching reads after one round and try again
    bool filterBadReads;
    double filterMatchThreshold;

    // Use a prior for the reference sequence
    bool useReferencePrior;

    // Verbosity options for printing
    bool verboseTruePositives;
    bool verboseFalsePositives;
    bool verboseFalseNegatives;

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

    // Filter out any reads with a MAPQ score less than or equal to this.
    int64_t mapqFilter;

    // Number of rounds of iterative refinement to attempt to improve the partition.
    int64_t roundsOfIterativeRefinement;

    // Whether or not to write a gvcf as output
    bool writeGVCF;

    // What types of file formats of split reads to output
    bool writeSplitSams;
    bool writeUnifiedSam;
//    bool writeSplitBams;
};

void stRPHmmParameters_destruct(stRPHmmParameters *params);

void stRPHmmParameters_learnParameters(stRPHmmParameters *params, stList *profileSequences,
        stHash *referenceNamesToReferencePriors);

void stRPHmmParameters_printParameters(stRPHmmParameters *params, FILE *fH);

void stRPHmmParameters_setReadErrorSubstitutionParameters(stRPHmmParameters *params, double *readErrorSubModel);

void normaliseSubstitutionMatrix(double *subMatrix);

double *getEmptyReadErrorSubstitutionMatrix(stRPHmmParameters *params);

/*
 * _stRPHmm
 * Struct for read partitioning hmm
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

void logHmm(stRPHmm *hmm, stSet *reads1, stSet *reads2, stGenomeFragment *gF);

/*
 * _stRPColumn
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
 * _stRPCell
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
 * _stRPMergeColumn
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
 * _stRPMergeCell
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
 * _stGenomeFragment
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

    // The reference coordinates of the genotypes & other read info
    char *referenceName;
    int64_t refStart;
    int64_t length;
    uint8_t *referenceSequence;

    // Depth and allele counts
    uint8_t *hap1Depth;
    uint8_t *hap2Depth;
    uint8_t *alleleCountsHap1;
    uint8_t *alleleCountsHap2;
    uint8_t *allele2CountsHap1;
    uint8_t *allele2CountsHap2;
};

stGenomeFragment *stGenomeFragment_construct(stRPHmm *hmm, stList *path);

void stGenomeFragment_destruct(stGenomeFragment *genomeFragment);

void stGenomeFragment_refineGenomeFragment(stGenomeFragment *gF, stSet *reads1, stSet *reads2,
        stRPHmm *hmm, stList *path, int64_t maxIterations);

double getLogProbOfReadGivenHaplotype(uint64_t *haplotypeString, int64_t start, int64_t length,
                                      stProfileSeq *profileSeq, stRPHmmParameters *params);

/*
 * _stBaseMapper
 * Struct for alphabet and mapping bases to numbers
 */
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

/*
 * Parsing methods
 */

stRPHmmParameters *parseParameters(char *paramsFile, stBaseMapper *baseMapper);

int64_t parseReads(stList *profileSequences, char *bamFile, stBaseMapper *baseMapper, stRPHmmParameters *params);

int64_t parseReadsWithSingleNucleotideProbs(stList *profileSequences, char *bamFile, stBaseMapper *baseMapper,
                                            stRPHmmParameters *params, char *signalAlignDirectory, bool onlySignalAlign);

int64_t getAlignedReadLength(bam1_t *aln);
int64_t getAlignedReadLength2(bam1_t *aln, int64_t *start_softclip, int64_t *end_softclip);
int64_t getAlignedReadLength3(bam1_t *aln, int64_t *start_softclip, int64_t *end_softclip, bool boundaryAtMatch);

void countIndels(uint32_t *cigar, uint32_t ncigar, int64_t *numInsertions, int64_t *numDeletions);

// Verbosity for what's printed.  To add more verbose options, you need to update:
//  usage, setVerbosity, struct _stRPHmmParameters, stRPHmmParameters_printParameters, writeParamFile
#define LOG_TRUE_POSITIVES 1
#define LOG_FALSE_POSITIVES 2
#define LOG_FALSE_NEGATIVES 4
void setVerbosity(stRPHmmParameters *params, int64_t bitstring);

/*
 * File writing methods
 */

void writeVcfFragment(vcfFile *out, bcf_hdr_t *bcf_hdr, stGenomeFragment *gF, char *referenceName,
                      stBaseMapper *baseMapper, bool gvcf);

bcf_hdr_t* writeVcfHeader(vcfFile *out, stList *genomeFragments, char *referenceName);

void writeParamFile(char *outputFilename, stRPHmmParameters *params);

/*
 * _stGenotypeResults
 * Struct which stores information about relevant test results.
 */
struct _stGenotypeResults {

    // Variants in reference
    int64_t negatives;
    int64_t positives;
    int64_t homozygousVariantsInRef;
    int64_t homozygousVariantsInRef_Insertions;
    int64_t homozygousVariantsInRef_Deletions;
    int64_t hetsInRef;
    int64_t hetsInRef_Insertions;
    int64_t hetsInRef_Deletions;

    // Variants in evaluated vcf
    int64_t truePositives;
    int64_t falsePositives;
    int64_t trueNegatives;
    int64_t falseNegatives;

    // Stats for specific types of variants
    int64_t truePositiveIndels;
    int64_t falsePositiveIndels;
    int64_t truePositiveHomozygous;
    int64_t truePositiveHet;
    int64_t truePositiveHomozygousIndels;
    int64_t truePositiveHetIndels;

    // Types of errors
    int64_t error_missedHet;
    int64_t error_missedHet_Insertions;
    int64_t error_missedHet_Deletions;
    int64_t error_homozygousInRef;
    int64_t error_homozygous_Insertions;
    int64_t error_homozygous_Deletions;

    // Phasing
    int64_t switchErrors;
    float switchErrorDistance;
    int64_t uncertainPhasing;
};
void printGenotypeResults(stGenotypeResults *results);

/*
 * VCF comparison methods
 */

void compareVCFs(FILE *fh, stList *hmms, char *vcf_toEval, char *vcf_ref,
                 stBaseMapper *baseMapper, stGenotypeResults *results, stRPHmmParameters *params);

void compareVCFsBasic(FILE *fh, char *vcf_toEval, char *vcf_ref, stGenotypeResults *results);

void compareVCFs_debugWithBams(char *vcf_toEval, char *vcf_ref, char *bamFile1, char *bamFile2, char *referenceFasta,
                               stBaseMapper *baseMapper, stGenotypeResults *results, stRPHmmParameters *params);

// Tag definitions (for haplotype output)
#define HAPLOTYPE_TAG "ht"
#define MARGIN_PHASE_TAG "mp"

/*
 * _stReadHaplotypeSequence
 * Struct for tracking haplotypes for read
 */

struct _stReadHaplotypeSequence {
    int64_t readStart;
    int64_t phaseBlock;
    int64_t length;
    int8_t haplotype;
    void *next;
};
stReadHaplotypeSequence *stReadHaplotypeSequence_construct(int64_t readStart, int64_t phaseBlock, int64_t length,
                                                           int8_t haplotype);

char *stReadHaplotypeSequence_toString(stReadHaplotypeSequence *rhs);

char *stReadHaplotypeSequence_toStringEmpty();

void stReadHaplotypeSequence_destruct(stReadHaplotypeSequence * rhs);

/*
 * stReadHaplotypePartitionTable
 * Tracking haplotypes for all reads
 */

stReadHaplotypePartitionTable *stReadHaplotypePartitionTable_construct(int64_t initialSize);

void stReadHaplotypePartitionTable_add(stReadHaplotypePartitionTable *hpt, char *readName, int64_t readStart,
                                       int64_t phaseBlock, int64_t length, int8_t haplotype);

void stReadHaplotypePartitionTable_destruct(stReadHaplotypePartitionTable *hpt);

void populateReadHaplotypePartitionTable(stReadHaplotypePartitionTable *hpt, stGenomeFragment *gF, stRPHmm *hmm,
                                         stList *path);

// Output file writing methods
void writeHaplotypedSam(char *bamInFile, char *bamOutBase, stReadHaplotypePartitionTable *readHaplotypePartitions,
                        char *marginPhaseTag);

void writeSplitSams(char *bamInFile, char *bamOutBase, stReadHaplotypePartitionTable *readHaplotypePartitions,
                    char *marginPhaseTag);

void addProfileSeqIdsToSet(stSet *pSeqs, stSet *readIds);

/*
 * Polish functions
 */

/*
 * Parameter object for polish algorithm
 */

struct _polishParams {
	bool useRunLengthEncoding;
	double referenceBasePenalty; // used by poa_getConsensus to weight against picking the reference base
	double minPosteriorProbForAlignmentAnchor; // used by by poa_getAnchorAlignments to determine which alignment pairs
	// to use for alignment anchors during poa_realignIterative
	Hmm *hmm; // Pair hmm used for aligning reads to the reference.
	StateMachine *sM; // Statemachine derived from the hmm
	PairwiseAlignmentParameters *p; // Parameters object used for aligning
	RepeatSubMatrix *repeatSubMatrix; // Repeat submatrix
	// chunking configuration
	bool includeSoftClipping;
	uint64_t chunkSize;
	uint64_t chunkBoundary;
	double candidateVariantWeight; // The fraction (from 0 to 1) of the average position coverage needed to define a candidate variant
	uint64_t columnAnchorTrim; // The min distance between a column anchor and a candidate variant
	uint64_t maxConsensusStrings; // The maximum number of different consensus strings to consider for a substring.
	uint64_t maxPoaConsensusIterations; // Maximum number of poa_consensus / realignment iterations
	uint64_t minPoaConsensusIterations; // Minimum number of poa_consensus / realignment iterations
	uint64_t maxRealignmentPolishIterations; // Maximum number of poa_polish iterations
	uint64_t minRealignmentPolishIterations; // Minimum number of poa_polish iterations
};

PolishParams *polishParams_readParams(FILE *fileHandle);

void polishParams_printParameters(PolishParams *polishParams, FILE *fh);

void polishParams_destruct(PolishParams *polishParams);

/*
 * Basic data structures for representing a POA alignment.
 */

struct _Poa {
	char *refString; // The reference string
	stList *nodes;
};

struct _poaNode {
	stList *inserts; // Inserts that happen immediately after this position
	stList *deletes; // Deletes that happen immediately after this position
	char base; // Char representing base, e.g. 'A', 'C', etc.
	double *baseWeights; // Array of length SYMBOL_NUMBER, encoding the weight given go each base, using the Symbol enum
	stList *observations; // Individual events representing event, a list of PoaObservations
};

struct _poaInsert {
	char *insert; // String representing characters of insert e.g. "GAT", etc.
	double weightForwardStrand;
	double weightReverseStrand;
};

struct _poaDelete {
	int64_t length; // Length of delete
	double weightForwardStrand;
	double weightReverseStrand;
};

struct _poaBaseObservation {
	int64_t readNo;
	int64_t offset;
	double weight;
};

/*
 * Poa functions.
 */

double poaInsert_getWeight(PoaInsert *insert);

double poaDelete_getWeight(PoaDelete *delete);

/*
 * Creates a POA representing the given reference sequence, with one node for each reference base and a
 * prefix 'N' base to represent place to add inserts/deletes that precede the first position of the reference.
 */
Poa *poa_getReferenceGraph(char *reference);

/*
 * Adds to given POA the matches, inserts and deletes from the alignment of the given read to the reference.
 * Adds the inserts and deletes so that they are left aligned.
 */
void poa_augment(Poa *poa, char *read, bool readStrand, int64_t readNo, stList *matches, stList *inserts, stList *deletes);

/*
 * Creates a POA representing the reference and the expected inserts / deletes and substitutions from the
 * alignment of the given set of reads aligned to the reference. Anchor alignments is a set of pairwise
 * alignments between the reads and the reference sequence. There is one alignment for each read. See
 * poa_getAnchorAlignments. The anchorAlignments can be null, in which case no anchors are used.
 */
Poa *poa_realign(stList *reads, stList *alignments, char *reference, PolishParams *polishParams);

/*
 * Generates a set of anchor alignments for the reads aligned to a consensus sequence derived from the poa.
 * These anchors can be used to restrict subsequent alignments to the consensus to generate a new poa.
 * PoaToConsensusMap is a map from the positions in the poa reference sequence to the derived consensus
 * sequence. See poa_getConsensus for description of poaToConsensusMap. If poaToConsensusMap is NULL then
 * the alignment is just the reference sequence of the poa.
 */
stList *poa_getAnchorAlignments(Poa *poa, int64_t *poaToConsensusMap, int64_t noOfReads,
							    PolishParams *polishParams);

/*
 * Generates a set of maximal expected alignments for the reads aligned to the the POA reference sequence.
 * Unlike the draft anchor alignments, these are designed to be complete, high quality alignments.
 */
stList *poa_getReadAlignmentsToConsensus(Poa *poa, stList *reads, PolishParams *polishParams);

/*
 * Prints representation of the POA.
 */
void poa_print(Poa *poa, FILE *fH, float indelSignificanceThreshold, float strandBalanceRatio);

/*
 * Prints some summary stats on the POA.
 */
void poa_printSummaryStats(Poa *poa, FILE *fH);

/*
 * Creates a consensus reference sequence from the POA. poaToConsensusMap is a pointer to an
 * array of integers of length str(poa->refString), giving the index of the reference positions
 * alignment to the consensus sequence, or -1 if not aligned. It is initialised as a
 * return value of the function.
 */
char *poa_getConsensus(Poa *poa, int64_t **poaToConsensusMap, PolishParams *polishParams);

Poa *poa_polish(Poa *poa, stList *bamChunkReads, PolishParams *params);

char *poa_polish2(Poa *poa, stList *reads, PolishParams *params,
				  int64_t **poaToConsensusMap);

/*
 * Iteratively used poa_realign and poa_getConsensus to refine the median reference sequence
 * for the given reads and the starting reference.
 */
Poa *poa_realignIterative(stList *bamChunkReads, stList *alignments, char *reference, PolishParams *polishParams);

/*
 * Ad poa_realignIterative, but allows the specification of the min and max number of realignment cycles,
 * also, can switch between the "poa_polish" and the "poa_consensus" algorithm using hmmNotRealign (poa_consensus
 * if non-zero).
 */
Poa *poa_realignIterative2(stList *bamChunkReads,
						   stList *anchorAlignments, char *reference,
						   PolishParams *polishParams, bool hmmNotRealign,
						   int64_t minIterations, int64_t maxIterations);

/*
 * As poa_realignIterative, but takes a starting poa. Input poa is destroyed by function.
 */
Poa *poa_realignIterative3(Poa *poa, stList *bamChunkReads,
						   PolishParams *polishParams, bool hmmMNotRealign,
						   int64_t minIterations, int64_t maxIterations);

/*
 * Convenience function that iteratively polishes sequence using poa_consensus and then poa_polish for
 * a specified number of iterations.
 */
Poa *poa_realignAll(stList *bamChunkReads, stList *anchorAlignments, char *reference,
						  PolishParams *polishParams);

/*
 * Greedily evaluate the top scoring indels.
 */
Poa *poa_checkMajorIndelEditsGreedily(Poa *poa, stList *reads, PolishParams *polishParams);

void poa_destruct(Poa *poa);

/*
 * Finds shift, expressed as a reference coordinate, that the given substring str can
 * be shifted left in the refString, starting from a match at refStart.
 */
int64_t getShift(char *refString, int64_t refStart, char *str, int64_t length);

/*
 * Get sum of weights for reference bases in poa - proxy to agreement of reads
 * with reference.
 */
double poa_getReferenceNodeTotalMatchWeight(Poa *poa);

/*
 * Get sum of weights for delete in poa - proxy to delete disagreement of reads
 * with reference.
 */
double poa_getDeleteTotalWeight(Poa *poa);

/*
 * Get sum of weights for inserts in poa - proxy to insert disagreement of reads
 * with reference.
 */
double poa_getInsertTotalWeight(Poa *poa);

/*
 * Get sum of weights for non-reference bases in poa - proxy to disagreement of read positions
 * aligned with reference.
 */
double poa_getReferenceNodeTotalDisagreementWeight(Poa *poa);

/*
 * Functions for run-length encoding/decoding with POAs
 */

// Data structure for representing RLE strings
struct _rleString {
	char *rleString; //Run-length-encoded (RLE) string
	int64_t *repeatCounts; // Count of repeat for each position in rleString
	int64_t *rleToNonRleCoordinateMap; // For each position in the RLE string the corresponding, left-most position
	// in the expanded non-RLE string
	int64_t *nonRleToRleCoordinateMap; // For each position in the expanded non-RLE string the corresponding position
	// in the RLE string
	int64_t length; // Length of the rleString
	int64_t nonRleLength; // Length of the expanded non-rle string
};

RleString *rleString_construct(char *string);

void rleString_destruct(RleString *rlString);

/*
 * Generates the expanded non-rle string.
 */
char *rleString_expand(RleString *rleString);

// Data structure for storing log-probabilities of observing
// one repeat count given another
struct _repeatSubMatrix {
	double *logProbabilities;
	int64_t maximumRepeatLength;
	int64_t maxEntry;
};

/*
 * Reads the repeat count matrix from a given input file.
 */

RepeatSubMatrix *repeatSubMatrix_constructEmpty();

void repeatSubMatrix_destruct(RepeatSubMatrix *repeatSubMatrix);

/*
 * Gets the log probability of observing a given repeat conditioned on an underlying repeat count and base.
 */
double repeatSubMatrix_getLogProb(RepeatSubMatrix *repeatSubMatrix, Symbol base, bool strand,
								  int64_t observedRepeatCount, int64_t underlyingRepeatCount);

/*
 * As gets, but returns the address.
 */
double *repeatSubMatrix_setLogProb(RepeatSubMatrix *repeatSubMatrix, Symbol base, bool strand, int64_t observedRepeatCount, int64_t underlyingRepeatCount);

/*
 * Gets the log probability of observing a given set of repeat observations conditioned on an underlying repeat count and base.
 */
double repeatSubMatrix_getLogProbForGivenRepeatCount(RepeatSubMatrix *repeatSubMatrix, Symbol base,
        stList *observations, stList *rleReads, stList *bamChunkReads, int64_t underlyingRepeatCount);

/*
 * Gets the maximum likelihood underlying repeat count for a given set of observed read repeat counts.
 * Puts the ml log probility in *logProbabilty.
 */
int64_t repeatSubMatrix_getMLRepeatCount(RepeatSubMatrix *repeatSubMatrix, Symbol base, stList *observations,
        stList *rleReads, stList *bamChunkReads, double *logProbability);

/*
 * Takes a POA done in run-length space and returns the expanded consensus string in
 * non-run-length space as an RleString.
 */
RleString *expandRLEConsensus(Poa *poa, stList *rlReads, stList *bamChunkReads, RepeatSubMatrix *repeatSubMatrix);

/*
 * Translate a sequence of aligned pairs (as stIntTuples) whose coordinates are monotonically increasing
 * in both underlying sequences (seqX and seqY) into an equivalent run-length encoded space alignment.
 */
stList *runLengthEncodeAlignment(stList *alignment, RleString *seqX, RleString *seqY);

/*
 * Make edited string with given insert. Edit start is the index of the position to insert the string.
 */
char *addInsert(char *string, char *insertString, int64_t editStart);

/*
 * Make edited string with given insert. Edit start is the index of the first position to delete from the string.
 */
char *removeDelete(char *string, int64_t deleteLength, int64_t editStart);

/*
 * Generates aligned pairs and indel probs, but first crops reference to only include sequence from first
 * to last anchor position.
 */
void getAlignedPairsWithIndelsCroppingReference(char *reference, int64_t refLength,
		char *read, stList *anchorPairs,
		stList **matches, stList **inserts, stList **deletes, PolishParams *polishParams);

/*
 * Functions for processing BAMs
 */

// TODO: MOVE BAMCHUNKER TO PARSER .c

typedef struct _bamChunker {
    // file locations
	char *bamFile;
    // configuration
    uint64_t chunkSize;
	uint64_t chunkBoundary;
	bool includeSoftClip;
	PolishParams *params;
	// internal data
    stList *chunks;
    uint64_t chunkCount;
    int64_t itorIdx;
} BamChunker;

typedef struct _bamChunk {
	char *refSeqName;          // name of contig
	int64_t chunkBoundaryStart;  // the first 'position' where we have an aligned read
	int64_t chunkStart;        // the actual boundary of the chunk, calculations from chunkMarginStart to chunkStart
	//  should be used to initialize the probabilities at chunkStart
	int64_t chunkEnd;          // same for chunk end
	int64_t chunkBoundaryEnd;    // no reads should start after this position
	BamChunker *parent;        // reference to parent (may not be needed)
} BamChunk;

typedef struct _bamChunkRead {
	char *readName;          	// read name
	char *nucleotides;			// nucleotide string
	int64_t readLength;
	bool forwardStrand;			// whether the alignment is matched to the forward strand
	BamChunk *parent;        	// reference to parent chunk
} BamChunkRead;

BamChunker *bamChunker_construct(char *bamFile, PolishParams *params);
BamChunker *bamChunker_construct2(char *bamFile, char *region, PolishParams *params);
void bamChunker_destruct(BamChunker *bamChunker);
BamChunk *bamChunker_getNext(BamChunker *bamChunker);

BamChunk *bamChunk_construct();
BamChunk *bamChunk_construct2(char *refSeqName, int64_t chunkBoundaryStart, int64_t chunkStart, int64_t chunkEnd,
                              int64_t chunkBoundaryEnd, BamChunker *parent);
void bamChunk_destruct(BamChunk *bamChunk);

BamChunkRead *bamChunkRead_construct();
BamChunkRead *bamChunkRead_construct2(char *readName, char *nucleotides, bool forwardStrand, BamChunk *parent);
BamChunkRead *bamChunkRead_constructRLECopy(BamChunkRead  *read, RleString *rle);
void bamChunkRead_destruct(BamChunkRead *bamChunkRead);

/*
 * Converts chunk of aligned reads into list of reads and alignments.
 */
uint32_t convertToReadsAndAlignments(BamChunk *bamChunk, stList *reads, stList *alignments);

/*
 * Remove overlap between two overlapping strings. Returns max weight of split point.
 */
int64_t removeOverlap(char *prefixString, char *suffixString, int64_t approxOverlap, PolishParams *polishParams,
				      int64_t *prefixStringCropEnd, int64_t *suffixStringCropStart);

/*
 * View functions
 */

struct _refMsaView {
	int64_t refLength; // The length of the reference sequence
	char *refSeq; // The reference sequence - this is not copied by the constructor
	char *refSeqName; // The reference sequence name - this is not copied by the constructor, and can be NULL
	int64_t seqNo; // The number of non-ref sequences aligned to the reference
	stList *seqs; // The non-ref sequences - - this is not copied by the constructor
	stList *seqNames; // The non-ref sequence names - this is not copied by the constructor, and can be NULL
	int64_t *seqCoordinates; // A matrix giving the coordinates of the non-reference sequence
	// as aligned to the reference sequence
	int64_t *maxPrecedingInsertLengths; // The maximum length of an insert in
	// any of the sequences preceding the reference positions
	int64_t **precedingInsertCoverages; // The number of sequences with each given indel position
};

/*
 * Get the coordinate in the given sequence aligned to the given reference position. Returns -1 if
 * no sequence position is aligned to the reference position.
 */
int64_t msaView_getSeqCoordinate(MsaView *view, int64_t refCoordinate, int64_t seqIndex);

/*
 * Gets the length of any insert in the given sequence preceding the given reference position. If
 * the sequence is not aligned at the given reference position returns 0.
 */
int64_t msaView_getPrecedingInsertLength(MsaView *view, int64_t rightRefCoordinate, int64_t seqIndex);

/*
 * Gets the first position in the sequence of an insert preceding the given reference position. If there
 * is no such insert returns -1.
 */
int64_t msaView_getPrecedingInsertStart(MsaView *view, int64_t rightRefCoordinate, int64_t seqIndex);

/*
 * Gets the maximum length of an indel preceding the given reference position
 */
int64_t msaView_getMaxPrecedingInsertLength(MsaView *view, int64_t rightRefCoordinate);

/*
 * Get the number of sequences with an insertion at a given position. IndelOffset if the position, from 0, of the
 indel from left-to-right.
 */
int64_t msaView_getPrecedingCoverageDepth(MsaView *view, int64_t rightRefCoordinate, int64_t indelOffset);

/*
 * Get the maximum length of an insertion at a given position with a minimum of reads supporting it.
 */
int64_t msaView_getMaxPrecedingInsertLengthWithGivenCoverage(MsaView *view, int64_t rightRefCoordinate, int64_t minCoverage);

/*
 * Builds an MSA view for the given reference and aligned sequences.
 * Does not copy the strings or string names, just holds references.
 */
MsaView *msaView_construct(char *refSeq, char *refName,
						   stList *refToSeqAlignments, stList *seqs, stList *seqNames);

void msaView_destruct(MsaView * view);

/*
 * Prints a quick view of the MSA for debugging/browsing.
 */
void msaView_print(MsaView *view, int64_t minInsertCoverage, FILE *fh);

/*
 * Prints the repeat counts of the MSA.
 */
void msaView_printRepeatCounts(MsaView *view, int64_t minInsertCoverage,
							   RleString *refString, stList *rleStrings, FILE *fh);

/*
 * Phase to polish functions
 */

void phaseReads(char *reference, int64_t referenceLength, stList *reads, stList *anchorAlignments,
				stList **readsPartition1, stList **readsPartition2, Params *params);

#endif /* ST_RP_HMM_H_ */
