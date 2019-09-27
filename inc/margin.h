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
typedef struct _stSite stSite;
typedef struct _stReference stReference;
typedef struct _stProfileSeq stProfileSeq;
typedef struct _stRPHmm stRPHmm;
typedef struct _stRPHmmParameters stRPHmmParameters;
typedef struct _stRPColumn stRPColumn;
typedef struct _stRPCell stRPCell;
typedef struct _stRPMergeColumn stRPMergeColumn;
typedef struct _stRPMergeCell stRPMergeCell;
typedef struct _stGenomeFragment stGenomeFragment;

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
};

Params *params_readParams(FILE *fp);

void params_destruct(Params *params);

void params_printParameters(Params *params, FILE *fh);

int64_t getAlignedReadLength3(bam1_t *aln, int64_t *start_softclip, int64_t *end_softclip, bool boundaryAtMatch);

/*
 * Overall coordination functions
 */

stList *filterReadsByCoverageDepth(stList *profileSeqs, stRPHmmParameters *params, stList *filteredProfileSeqs,
        stList *discardedProfileSeqs);

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
#define ALLELE_LOG_PROB_BITS 8


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
 * Reference / site definition
 */

struct _stSite {
	uint64_t alleleNumber; // Number of alleles at the site in the reference
	uint64_t alleleOffset; // The index of the first allele in this site
	// in a sequence of all alleles in the reference, ordered first by site then
	// by order in the site.
	uint16_t *substitutionLogProbs; // Log probabilities of substitutions between the alleles
	uint16_t *allelePriorLogProbs; // Prior log-prob on alleles, allows upweighting of reference allele
};

uint16_t *stSite_getSubstitutionProb(stSite *s, int64_t from, int64_t to);

struct _stReference {
	char *referenceName;
	uint64_t length; // Number of sites
	uint64_t totalAlleles; // Total number of alleles across all sites
	stSite *sites;
};

void stReference_destruct(stReference *ref);

/*
 * _stProfileSeq
 * Struct for profile sequence
 */

struct _stProfileSeq {
    stReference *ref;
    char *readId;
    uint64_t refStart; // The first site in the reference
    uint64_t length; // Number of reference sites

    uint64_t alleleOffset; // The index of the first allele in this sequence
    // in a sequence of all alleles in the reference, ordered first by site then
    // by order in the site.

    // The log-probability of alleles, as specified by uint8_t
    // Each is expressed as an 8 bit unsigned int, with the value from 0 to -255
    uint8_t *profileProbs;
};

stProfileSeq *stProfileSeq_constructEmptyProfile(stReference *ref, char *readId,
                                                 int64_t referenceStart, int64_t length);

void stProfileSeq_destruct(stProfileSeq *seq);

void stProfileSeq_print(stProfileSeq *seq, FILE *fileHandle);

uint8_t *stProfileSeq_getProb(stProfileSeq *seq, uint64_t site, uint64_t allele);

int stRPProfileSeq_cmpFn(const void *a, const void *b);

/*
 * Emission probabilities methods
 */
double emissionLogProbability(stRPColumn *column, stRPCell *cell, uint64_t *bitCountVectors,
                                stReference *reference,
                                stRPHmmParameters *params);

void fillInPredictedGenome(stGenomeFragment *gF, uint64_t partition,
        stRPColumn *column, stRPHmmParameters *params);

/*
 * Constituent functions tested and used to do bit twiddling
*/
int popcount64(uint64_t x);

uint64_t getLogProbOfAllele(uint64_t *bitCountVectors, uint64_t depth, uint64_t partition,
							uint64_t siteOffset, uint64_t allele);

uint64_t *calculateCountBitVectors(uint8_t **seqs, stReference *ref,
								   uint64_t firstSite, uint64_t length, uint64_t depth);

/*
 * _stRPHmmParameters
 * Struct for hmm parameters
 */
struct _stRPHmmParameters {
    /*
     * Parameters used for the HMM computation
     */
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

    // Ensure symmetry in the HMM such that the inverted partition of each partition is included in the HMM
    bool includeInvertedPartitions;

    // Number of rounds of iterative refinement to attempt to improve the partition.
    int64_t roundsOfIterativeRefinement;
};

void stRPHmmParameters_destruct(stRPHmmParameters *params);

void stRPHmmParameters_printParameters(stRPHmmParameters *params, FILE *fH);

/*
 * _stRPHmm
 * Struct for read partitioning hmm
 */
struct _stRPHmm {
	stReference *ref;
    int64_t refStart; // First site in reference
    int64_t refLength; // Number of sites in the reference
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

void logHmm(stRPHmm *hmm, stGenomeFragment *gF);

/*
 * _stRPColumn
 * Column of read partitioning hmm
 */
struct _stRPColumn {
    int64_t refStart; // First site in the reference
    int64_t length; // Number of sites
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
	// The reference coordinates of the genotypes & other read info
	stReference *reference; // The reference this fragment refers to
	uint64_t refStart; // First site in the reference
	uint64_t length; // The number of sites

	// Reads
	stSet *reads1; // The reads in the first partition
	stSet *reads2; // The reads in the second partition

    // A string where each element represents the predicted genotype at the corresponding
    // position.
    // A genotype is represented by an integer in the range [0, allele_number**2), where
	// where allele_number is the number alleles at the given site
    // A genotype expresses two characters. For two characters x, y represented by two integers
    // in [0, allele_number) then the genotype is expressed as x * allele_number + y if x <= y
    // else y * allele_number + x
    uint64_t *genotypeString;

    // Strings representing the predicted haplotypes, where each element is a reference to an allele
    uint64_t *haplotypeString1;
    uint64_t *haplotypeString2;

    // An array of genotype posterior probabilities,
    // each between 0 and 1, for the corresponding genotypes
    // in the genotype string
    float *genotypeProbs;
};

stGenomeFragment *stGenomeFragment_construct(stRPHmm *hmm, stList *path);

void stGenomeFragment_destruct(stGenomeFragment *genomeFragment);

void stGenomeFragment_refineGenomeFragment(stGenomeFragment *gF,
        stRPHmm *hmm, stList *path, int64_t maxIterations);

double getLogProbOfReadGivenHaplotype(uint64_t *haplotypeString, int64_t start, int64_t length,
                                      stProfileSeq *profileSeq, stReference *ref);

/*
 * Parsing methods
 */

stRPHmmParameters *parseParameters(char *paramsFile);


/*
 * Polish functions
 */

/*
 * Parameter object for polish algorithm
 */

struct _polishParams {
	bool useRunLengthEncoding;
	double referenceBasePenalty; // used by poa_getConsensus to weight against picking the reference base
	double *minPosteriorProbForAlignmentAnchors; // used by by poa_getAnchorAlignments to determine which alignment pairs
	// to use for alignment anchors during poa_realignIterative, of the form of even-length array of form
	// [ min_posterio_anchor_prob_1, diagonal_expansion_1,  min_posterio_anchor_prob_2, diagonal_expansion_2, ... ]
	int64_t minPosteriorProbForAlignmentAnchorsLength;  // Length of array minPosteriorProbForAlignmentAnchors
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

	uint64_t minReadsToCallConsensus; // Min reads to choose between consensus sequences for a region
	uint64_t filterReadsWhileHaveAtLeastThisCoverage; // Only filter read substrings if we have at least this coverage
	// at a locus
	double minAvgBaseQuality; // Minimum average base quality to include a substring for consensus finding
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
Poa *poa_realign(stList *bamChunkReads, stList *alignments, char *reference, PolishParams *polishParams);

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
stList *poa_getReadAlignmentsToConsensus(Poa *poa, stList *bamChunkReads, PolishParams *polishParams);

/*
 * Prints representation of the POA.
 */
void poa_print(Poa *poa, FILE *fH,
			  stList *bamChunkReads,
			  float indelSignificanceThreshold, float strandBalanceRatio);

/*
 * Prints a tab separated version of the POA graph.
 */
void poa_printTSV(Poa *poa, FILE *fH,
		stList *bamChunkReads,
		float indelSignificanceThreshold, float strandBalanceRatio);

/*
 * Print repeat count observations.
 */
void poa_printRepeatCounts(Poa *poa, FILE *fH, stList *bamChunkReads);

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

char *poa_polish2(Poa *poa, stList *bamChunkReads, PolishParams *params,
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
Poa *poa_checkMajorIndelEditsGreedily(Poa *poa, stList *bamChunkReads, PolishParams *polishParams);

void poa_destruct(Poa *poa);

double *poaNode_getStrandSpecificBaseWeights(PoaNode *node, stList *bamChunkReads,
											 double *totalWeight, double *totalPositiveWeight, double *totalNegativeWeight);

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

// Data structure for representing RLE strings
struct _rleString {
	char *rleString; //Run-length-encoded (RLE) string
	uint64_t *repeatCounts; // Count of repeat for each position in rleString
	uint64_t length; // Length of the rleString
	uint64_t nonRleLength; // Length of the expanded, non-rle string
};

RleString *rleString_construct(char *string);

RleString *rleString_construct_no_rle(char *string);

void rleString_destruct(RleString *rlString);

/*
 * Generates the expanded non-rle string.
 */
char *rleString_expand(RleString *rleString);

/*
 * Gets an array giving the position in the rleString of a corresponding position in the expanded string.
 */
uint64_t *rleString_getNonRleToRleCoordinateMap(RleString *rleString);

uint8_t *rleString_rleQualities(RleString *rleString, uint8_t *qualities);

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
        stList *observations, stList *bamChunkReads, int64_t underlyingRepeatCount);

/*
 * Gets the maximum likelihood underlying repeat count for a given set of observed read repeat counts.
 * Puts the ml log probility in *logProbabilty.
 */
int64_t repeatSubMatrix_getMLRepeatCount(RepeatSubMatrix *repeatSubMatrix, Symbol base, stList *observations,
        stList *bamChunkReads, double *logProbability);

/*
 * Takes a POA done in run-length space and returns the expanded consensus string in
 * non-run-length space as an RleString.
 */
RleString *expandRLEConsensus(Poa *poa, stList *bamChunkReads, RepeatSubMatrix *repeatSubMatrix);

/*
 * Translate a sequence of aligned pairs (as stIntTuples) whose coordinates are monotonically increasing
 * in both underlying sequences (seqX and seqY) into an equivalent run-length encoded space alignment.
 */
stList *runLengthEncodeAlignment(stList *alignment,
		uint64_t *seqXNonRleToRleCoordinateMap, uint64_t *seqYNonRleToRleCoordinateMap);

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
	RleString *rleRead; 		// rle read
	uint8_t *qualities;			// quality scores. will be NULL if not given, else will be of length rleRead->length
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
BamChunkRead *bamChunkRead_construct2(char *readName, char *nucleotides, uint8_t *qualities, bool forwardStrand, BamChunk *parent);
void bamChunkRead_destruct(BamChunkRead *bamChunkRead);

/*
 * Generates the expanded non-rle version of bam chunk read nucleotide sequence.
 */
char *bamChunkRead_rleExpand(BamChunkRead *read);

typedef struct _bamChunkReadSubstring {
	BamChunkRead *read; // The parent read from which the substring arises from
	uint64_t start; // The 0 based offset of the start position in the parent read (inclusive)
	uint64_t length; // The length of the substring
	char *readSubstring; // TODO: Test not generating this but getting via method
	double qualValue;
} BamChunkReadSubstring;

/*
 * Converts chunk of aligned reads into list of reads and alignments.
 */
uint32_t convertToReadsAndAlignments(BamChunk *bamChunk, RleString *reference, stList *reads, stList *alignments);

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
 * Bubble graphs
 */

typedef struct _bubble {
	uint64_t refStart; //First inclusive position
	uint64_t length; // Length of the reference sub-sequence covered by the bubble
	char *refAllele; // The current reference allele
	uint64_t alleleNo; // Number of alleles
	char **alleles; // Array of allele strings
	uint64_t readNo; // Number of reads overlapping bubble
	BamChunkReadSubstring **reads; // Array of read substrings aligned to the bubble
	float *alleleReadSupports; // An array of log-likelihoods giving the support of
	// each allele for each read
	uint64_t alleleOffset; // The index of the first allele in this bubble
	// in a sequence of all alleles in the bubble graph, ordered first by bubble then
	// by order in the bubble.
} Bubble;

typedef struct _bubbleGraph {
	char *refString; // The reference string
	uint64_t refLength; // The length of the reference string
	uint64_t bubbleNo; // The number of bubbles
	Bubble *bubbles; // An array of bubbles
} BubbleGraph;

/*
 * Get a consensus path through bubble graph by picking the highest
 * likelihood allele at each bubble. Returned as a string of bg->refLength integers,
 * each denoting a chosen allele.
 */
uint64_t *bubbleGraph_getConsensusPath(BubbleGraph *bg, PolishParams *polishParams);

/*
 * Get a consensus sequences from the bubble graph by picking the highest
 * likelihood allele at each bubble.
 */
char *bubbleGraph_getConsensusString(BubbleGraph *bg, uint64_t *consensusPath, int64_t **poaToConsensusMap, PolishParams *polishParams);

/*
 * Create a bubble graph from a POA.
 */
BubbleGraph *bubbleGraph_constructFromPoa(Poa *poa, stList *bamChunkReads, PolishParams *params);

void bubbleGraph_destruct(BubbleGraph *bg);

/*
 * Prints a quick view of the bubble graph for debugging/browsing.
 */
void bubbleGraph_print(BubbleGraph *bg, FILE *fh);

/*
 * The the index in b->alleles of the allele with highest likelihood
 * given the reads
 */
int64_t bubble_getIndexOfHighestLikelihoodAllele(Bubble *b);

/*
 * Gets the likelihood of a given allele giving rise to the reads.
 */
double bubble_getLogLikelihoodOfAllele(Bubble *b, int64_t allele);

/*
 * Gets a set of profile sequences for the reads aligned to the bubble graph.
 * Allows them to then be phased.
 */
stList *bubbleGraph_getProfileSeqs(BubbleGraph *bg, stReference *ref);

/*
 * Gets an stReference that can be used for phasing.
 */
stReference *bubbleGraph_getReference(BubbleGraph *bg);

/*
 * Phase bubble graph.
 */
stGenomeFragment *bubbleGraph_phaseBubbleGraph(BubbleGraph *bg, Params *params);

/*
 * Get Poa from bubble graph.
 */
Poa *bubbleGraph_getNewPoa(BubbleGraph *bg, uint64_t *consensusPath, Poa *poa, stList *reads, Params *params);

#endif /* ST_RP_HMM_H_ */
