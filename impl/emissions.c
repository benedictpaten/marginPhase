/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "stRPHmm.h"

/*
 * Character alphabet and substitutions
 */

uint16_t *getSubstitutionProb(uint16_t *matrix, int64_t from, int64_t to) {
    /*
     * Gets the (log) substitution probability of getting the derived (to) character given the source (from/haplotype) character.
     */
    return &matrix[from * ALPHABET_SIZE + to];
}

double *getSubstitutionProbSlow(double *matrix, int64_t from, int64_t to) {
    /*
     * As getSubstitutionProb.
     */
    return &matrix[from * ALPHABET_SIZE + to];
}

uint16_t scaleToLogIntegerSubMatrix(double logProb) {
    /*
     * Convert log probability into scaled form for substitution matrix.
     */
    assert(logProb <= 0);
    if(logProb < -10) {
        st_errAbort("Attempting to set a substitution probability smaller than x=0.00001 (log(x) ~= -12)");
    }
    return round(ALPHABET_MIN_SUBSTITUTION_PROB * (-logProb/12.0));
}

double invertScaleToLogIntegerSubMatrix(int64_t i) {
    /*
     * Invert scaled form to log probability.
     */
    return (12.0 * ((double)-i))/ALPHABET_MIN_SUBSTITUTION_PROB;
}

void setSubstitutionProb(uint16_t *logSubMatrix, double *logSubMatrixSlow,
        int64_t sourceCharacterIndex,
        int64_t derivedCharacterIndex, double prob) {
    /*
     * Sets the substitution probability, scaling it appropriately by taking the log and then storing as integer (see definition)
     */
    if(prob <= 0 || prob > 1.0) {
        st_errAbort("Attempting to set substitution probability out of 0-1 range");
    }
    *getSubstitutionProb(logSubMatrix, sourceCharacterIndex, derivedCharacterIndex) = scaleToLogIntegerSubMatrix(log(prob));
    *getSubstitutionProbSlow(logSubMatrixSlow, sourceCharacterIndex, derivedCharacterIndex) = log(prob);
}

/*
 * Emission probabilities with optimization to make fast
 */

/*
 * Following implement Hamming weight for uint64_t ints, taken from
 * https://en.wikipedia.org/wiki/Hamming_weight
 */

//types and constants used in the functions below
//uint64_t is an unsigned 64-bit integer variable type (defined in C99 version of C language)
const uint64_t m1  = 0x5555555555555555; //binary: 0101...
const uint64_t m2  = 0x3333333333333333; //binary: 00110011..
const uint64_t m4  = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
const uint64_t h01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...

inline int popcount64NoBuiltIn(uint64_t x) {
    /*
     * Returns Hamming weight of input unsigned integer.
     */
    //This uses fewer arithmetic operations than any other known
    //implementation on machines with fast multiplication.
    //This algorithm uses 12 arithmetic operations, one of which is a multiply.
    x -= (x >> 1) & m1;             //put count of each 2 bits into those 2 bits
    x = (x & m2) + ((x >> 2) & m2); //put count of each 4 bits into those 4 bits
    x = (x + (x >> 4)) & m4;        //put count of each 8 bits into those 8 bits
    return (x * h01) >> 56;  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ...
}

inline int popcount64(uint64_t x) {
    /*
     * Returns Hamming weight of input unsigned integer.
     */
    //return popcount64NoBuiltIn(x); // Use this line if the builtin is unavailable
    return __builtin_popcountll(x);
}

static inline uint64_t *retrieveBitCountVector(uint64_t *bitCountVector,
        int64_t position, int64_t characterIndex, int64_t bit) {
    /*
     * Returns a pointer to a bit count vector for a given position (offset in the column),
     * character index and bit.
     */
    return &bitCountVector[position * ALPHABET_CHARACTER_BITS * ALPHABET_SIZE + characterIndex * ALPHABET_CHARACTER_BITS + bit];
}

uint64_t calculateBitCountVector(uint8_t **seqs, int64_t depth,
        int64_t position, int64_t characterIndex, int64_t bit) {
    /*
     * Calculates the bit count vector for a given position, character index and bit.
     */
    uint64_t bitCountVector = 0;
    for(int64_t i=0; i<depth; i++) {
        uint8_t *p = &(seqs[i][ALPHABET_SIZE * position]);
        bitCountVector |= ((((uint64_t)p[characterIndex] >> bit) & 1) << i);
    }

    return bitCountVector;
}

uint64_t *calculateCountBitVectors(uint8_t **seqs, int64_t depth, int64_t length) {
    /*
     * Calculates the bit count vector for every position, character and bit in the column.
     */

    // Array of bit vectors, for each position, for each character and for each bit in uint8_t
    uint64_t *bitCountVectors = st_malloc(length * ALPHABET_SIZE *
            ALPHABET_CHARACTER_BITS * sizeof(uint64_t));

    // For each position
    for(int64_t i=0; i<length; i++) {
        // For each character
        for(int64_t j=0; j<ALPHABET_SIZE; j++) {
            // For each bit
            for(int64_t k=0; k<ALPHABET_CHARACTER_BITS; k++) {
                *retrieveBitCountVector(bitCountVectors, i, j, k) =
                        calculateBitCountVector(seqs, depth, i, j, k);
            }
        }
    }

    return bitCountVectors;
}

uint64_t getExpectedInstanceNumber(uint64_t *bitCountVectors, uint64_t depth, uint64_t partition,
        int64_t position, int64_t characterIndex) {
    /*
     * Returns the number of instances of a character, given by characterIndex, at the given position within the column for
     * the given partition. Returns value scaled between 0 and ALPHABET_MAX_PROB, where the return value divided by ALPHABET_MAX_PROB
     * is the expected number of instances of the given character in the given subpartition of the column.
     */
    uint64_t *j = retrieveBitCountVector(bitCountVectors, position, characterIndex, 0);
    uint64_t expectedCount = popcount64(j[0] & partition);

    for(int64_t i=1; i<ALPHABET_CHARACTER_BITS; i++) {
        expectedCount += (popcount64(j[i] & partition) << i);
    }

    assert(expectedCount >= 0.0);
    assert((double)expectedCount / ALPHABET_MAX_PROB <= depth);
    return expectedCount;
}

static inline uint64_t getLogProbOfReadCharacters(uint16_t *logSubMatrix, uint64_t *expectedInstanceNumbers,
        int64_t sourceCharacterIndex) {
    /*
     * Get the log probability of a given source character given the expected number of instances of
     * each character in the reads.
     */
    uint16_t *j = getSubstitutionProb(logSubMatrix, sourceCharacterIndex, 0);
    uint64_t logCharacterProb = j[0] * expectedInstanceNumbers[0];

    for(int64_t i=1; i<ALPHABET_SIZE; i++) {
        logCharacterProb += j[i] * expectedInstanceNumbers[i];
    }

    return logCharacterProb;
}

static inline uint64_t minP(uint64_t a, uint64_t b) {
    return a < b ? a : b;
}

static inline void columnIndexLogHapProbability(stRPColumn *column, uint64_t index,
        uint64_t partition, uint64_t *bitCountVectors,
        stRPHmmParameters *params,
        uint64_t *rootCharacterProbs) {
    /*
     * Get the probabilities of the "root" characters for a given read sub-partition and a haplotype.
     */
    // For each possible read character calculate the expected number of instances in the
    // partition and store counts in an array
    uint64_t expectedInstanceNumbers[ALPHABET_SIZE];
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        expectedInstanceNumbers[i] = getExpectedInstanceNumber(bitCountVectors,
                               column->depth, partition, index, i);
    }

    // Calculate the probability of the read characters for each possible haplotype character
    uint64_t characterProbsHap[ALPHABET_SIZE];
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        characterProbsHap[i] = getLogProbOfReadCharacters(params->readErrorSubModel, expectedInstanceNumbers, i);
    }

    // Calculate the probability of haplotype characters and read characters for each root character
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        uint16_t *j = getSubstitutionProb(params->hetSubModel, i, 0);
        rootCharacterProbs[i] = characterProbsHap[0] + j[0] * ALPHABET_MAX_PROB;
        for(int64_t k=1; k<ALPHABET_SIZE; k++) {
            rootCharacterProbs[i] =
                    minP(rootCharacterProbs[i],
                         characterProbsHap[k] + j[k] * ALPHABET_MAX_PROB);
        }
    }
}

static inline uint64_t columnIndexLogProbability(stRPColumn *column, uint64_t index,
        uint64_t partition, uint64_t *bitCountVectors, uint16_t *referencePriorProbs,
        stRPHmmParameters *params) {
    /*
     * Get the probability of the characters in a given position within a column for a given partition.
     */
    // Get the sum of log probabilities of the derived characters over the possible source characters
    uint64_t rootCharacterProbsHap1[ALPHABET_SIZE];
    columnIndexLogHapProbability(column, index,
            partition, bitCountVectors, params, rootCharacterProbsHap1);
    uint64_t rootCharacterProbsHap2[ALPHABET_SIZE];
    columnIndexLogHapProbability(column, index,
            ~partition, bitCountVectors, params, rootCharacterProbsHap2);

    // Combine the probabilities to calculate the overall probability of a given position in a column
    uint64_t logColumnProb = rootCharacterProbsHap1[0] + rootCharacterProbsHap2[0] + referencePriorProbs[0] * ALPHABET_MAX_PROB;
    for(int64_t i=1; i<ALPHABET_SIZE; i++) {
        logColumnProb = minP(logColumnProb, rootCharacterProbsHap1[i] + rootCharacterProbsHap2[i] + referencePriorProbs[i] * ALPHABET_MAX_PROB); // + (i == ALPHABET_SIZE-1 ? scaleToLogIntegerSubMatrix(0.001) : 0));
    }

    return logColumnProb;
}

double emissionLogProbability(stRPColumn *column,
        stRPCell *cell, uint64_t *bitCountVectors, stReferencePriorProbs *referencePriorProbs,
        stRPHmmParameters *params) {
    /*
     * Get the log probability of a set of reads for a given column.
     */
    assert(column->length > 0);
    uint16_t *rProbs = &referencePriorProbs->profileProbs[(column->refStart - referencePriorProbs->refStart) * ALPHABET_SIZE];
    uint64_t logPartitionProb = columnIndexLogProbability(column, 0,
            cell->partition, bitCountVectors, rProbs, params);

    for(int64_t i=1; i<column->length; i++) {
        rProbs = &rProbs[ALPHABET_SIZE]; // Move to the next column of the reference prior
        logPartitionProb = logPartitionProb + columnIndexLogProbability(column, i,
                cell->partition, bitCountVectors, rProbs, params);
    }

    return invertScaleToLogIntegerSubMatrix(logPartitionProb)/ALPHABET_MAX_PROB;
}

/*
 * Functions for calculating genotypes/haplotypes
 */

double getLogProbOfReadCharactersSlow(double *logSubMatrix, uint64_t *expectedInstanceNumbers,
        int64_t sourceCharacterIndex) {
    /*
     * Get the log probability of a given source character given the expected number of instances of
     * each character in the reads.
     */
    double logCharacterProb = *getSubstitutionProbSlow(logSubMatrix, sourceCharacterIndex, 0) *
            ((double)expectedInstanceNumbers[0]);

    for(int64_t i=1; i<ALPHABET_SIZE; i++) {
        logCharacterProb += *getSubstitutionProbSlow(logSubMatrix, sourceCharacterIndex, i) *
                ((double)expectedInstanceNumbers[i]);
    }

    return logCharacterProb/ALPHABET_MAX_PROB;
}

void columnIndexLogHapProbabilitySlow(stRPColumn *column, uint64_t index,
        uint64_t partition, uint64_t *bitCountVectors, stRPHmmParameters *params, double *characterProbsHap) {
    /*
     * Get the probabilities of the haplotype characters for a given read sub-partition and a haplotype.
     */
    // For each possible read character calculate the expected number of instances in the
    // partition and store counts in an array
    uint64_t expectedInstanceNumbers[ALPHABET_SIZE];
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        expectedInstanceNumbers[i] = getExpectedInstanceNumber(bitCountVectors,
                               column->depth, partition, index, i);
    }

    // Calculate the probability of the read characters for each possible haplotype character
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        characterProbsHap[i] = getLogProbOfReadCharactersSlow(params->readErrorSubModelSlow, expectedInstanceNumbers, i);
    }
}

void calculateRootCharacterProbs(double *characterProbsHap, stRPHmmParameters *params,
        double *rootCharacterProbs, bool maxNotSum) {
    /*
     * Calculate the probability of haplotype characters and read characters for each root character
     * given the probability of each individual haplotype character
     */
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        rootCharacterProbs[i] = characterProbsHap[0] +
                *getSubstitutionProbSlow(params->hetSubModelSlow, i, 0);
        for(int64_t j=1; j<ALPHABET_SIZE; j++) {
            rootCharacterProbs[i] =
                    logAddP(rootCharacterProbs[i],
                            characterProbsHap[j] +
                            *getSubstitutionProbSlow(params->hetSubModelSlow, i, j),
                            maxNotSum);
        }
    }
}

void columnIndexLogRootHapProbabilitySlow(stRPColumn *column, uint64_t index,
        uint64_t partition, uint64_t *bitCountVectors,
        stRPHmmParameters *params, double *rootCharacterProbs, bool maxNotSum) {
    /*
     * Get the probabilities of the "root" characters for a given read sub-partition and a haplotype.
     */
    double characterProbsHap[ALPHABET_SIZE];

    columnIndexLogHapProbabilitySlow(column, index,
            partition, bitCountVectors, params, characterProbsHap);

    calculateRootCharacterProbs(characterProbsHap, params, rootCharacterProbs, maxNotSum);
}

double columnIndexLogProbabilitySlow(stRPColumn *column, uint64_t index,
        uint64_t partition, uint64_t *bitCountVectors, uint16_t *referencePriorProbs,
        stRPHmmParameters *params, bool maxNotSum) {
    /*
     * Get the probability of a the characters in a given position within a column for a given partition.
     */
    // Get the sum of log probabilities of the derived characters over the possible source characters
    double rootCharacterProbsHap1[ALPHABET_SIZE];
    columnIndexLogRootHapProbabilitySlow(column, index,
            partition, bitCountVectors, params, rootCharacterProbsHap1, maxNotSum);
    double rootCharacterProbsHap2[ALPHABET_SIZE];
    columnIndexLogRootHapProbabilitySlow(column, index,
            ~partition, bitCountVectors, params, rootCharacterProbsHap2, maxNotSum);

    // Combine the probabilities to calculate the overall probability of a given position in a column
    double logColumnProb = rootCharacterProbsHap1[0] + rootCharacterProbsHap2[0] +
            invertScaleToLogIntegerSubMatrix(referencePriorProbs[0]);
    for(int64_t i=1; i<ALPHABET_SIZE; i++) {
        logColumnProb = logAddP(logColumnProb, rootCharacterProbsHap1[i] + rootCharacterProbsHap2[i] +
                invertScaleToLogIntegerSubMatrix(referencePriorProbs[i]), maxNotSum);
    }

    return logColumnProb; // + log(1.0/params->alphabetSize);
}

double emissionLogProbabilitySlow(stRPColumn *column,
        stRPCell *cell, uint64_t *bitCountVectors, stReferencePriorProbs *referencePriorProbs,
        stRPHmmParameters *params, bool maxNotSum) {
    /*
     * Get the log probability of a set of reads for a given column.
     */
    assert(column->length > 0);
    uint16_t *rProbs = &referencePriorProbs->profileProbs[(column->refStart - referencePriorProbs->refStart) * ALPHABET_SIZE];
    double logPartitionProb = columnIndexLogProbabilitySlow(column, 0,
            cell->partition, bitCountVectors, rProbs, params, maxNotSum);
    for(int64_t i=1; i<column->length; i++) {
        rProbs = &rProbs[ALPHABET_SIZE]; // Move to the next column of the reference prior
        logPartitionProb += columnIndexLogProbabilitySlow(column, i,
            cell->partition, bitCountVectors, rProbs, params, maxNotSum);
    }
    return logPartitionProb;
}

uint64_t getMLHapChar(double *characterProbsHap,
        stRPHmmParameters *params,
        int64_t rootChar) {
    /*
     * Return the haplotype character with maximum probability.
     */
    int64_t maxProbHapChar = 0;

    double maxHapProb = characterProbsHap[0] +
                *getSubstitutionProbSlow(params->hetSubModelSlow,
                rootChar, 0);

    for(int64_t i=1; i<ALPHABET_SIZE; i++) {
        double hapProb = characterProbsHap[i] +
                *getSubstitutionProbSlow(params->hetSubModelSlow,
                rootChar, i);
        if(hapProb > maxHapProb) {
            maxHapProb = hapProb;
            maxProbHapChar = i;
        }
    }

    return maxProbHapChar;
}

double getHaplotypeProb(double characterReadProb,
            uint64_t hapChar, double *rootCharacterProbsOtherHap, stRPHmmParameters *params,
            uint16_t *referencePriorProbs) {
    /*
     * Return the probability of the tree given that the haplotype character was hapChar
     */

    double logHapProb = ST_MATH_LOG_ZERO;

    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        logHapProb = stMath_logAdd(logHapProb, characterReadProb +
                *getSubstitutionProbSlow(params->hetSubModelSlow,
                i, hapChar) + rootCharacterProbsOtherHap[i] + invertScaleToLogIntegerSubMatrix(referencePriorProbs[i]));
    }

    return logHapProb;
}

void fillInPredictedGenomePosition(stGenomeFragment *gF, stRPCell *cell,
        stRPColumn *column, stRPHmmParameters *params,
        stReferencePriorProbs *referencePriorProbs,
        uint64_t *bitCountVectors, uint64_t index) {
    /*
     * Computes the most probable haplotype characters / genotype and associated posterior
     * probabilities for a given position within a cell/column.
     */

    uint16_t *rProbs = &referencePriorProbs->profileProbs[(column->refStart-referencePriorProbs->refStart + index)*ALPHABET_SIZE];

    // Get the haplotype characters that are most probable given the root character
    double characterProbsHap1[ALPHABET_SIZE];
    double characterProbsHap2[ALPHABET_SIZE];

    columnIndexLogHapProbabilitySlow(column, index,
                    cell->partition, bitCountVectors,
                    params, characterProbsHap1);

    columnIndexLogHapProbabilitySlow(column, index,
                        ~cell->partition, bitCountVectors,
                        params, characterProbsHap2);

    // Get the root character with maximum posterior probability.

    // Get the sum of log probabilities of the derived characters over the possible source characters
    double rootCharacterProbsHap1[ALPHABET_SIZE];
    double rootCharacterProbsHap2[ALPHABET_SIZE];

    calculateRootCharacterProbs(characterProbsHap1, params,
            rootCharacterProbsHap1, 0);
    calculateRootCharacterProbs(characterProbsHap2, params,
                rootCharacterProbsHap2, 0);

    // Combine the probabilities to calculate the overall probability of a given partition of
    // read characters and the root character with maximum posterior prob
    double logColumnProbSum = rootCharacterProbsHap1[0] + rootCharacterProbsHap2[0] +
            invertScaleToLogIntegerSubMatrix(rProbs[0]);
    double logColumnProbMax = logColumnProbSum;
    int64_t maxProbRootChar = 0;
    for(int64_t i=1; i<ALPHABET_SIZE; i++) {
        double logColumnProb = rootCharacterProbsHap1[i] + rootCharacterProbsHap2[i] +
                invertScaleToLogIntegerSubMatrix(rProbs[i]);
        logColumnProbSum = stMath_logAdd(logColumnProbSum, logColumnProb);
        if(logColumnProb > logColumnProbMax) {
            logColumnProbMax = logColumnProb;
            maxProbRootChar = i;
        }
    }

    int64_t j = column->refStart + index - gF->refStart;

    // Get the haplotype characters with highest posterior probability.
    uint64_t hapChar1  = getMLHapChar(characterProbsHap1, params, maxProbRootChar);
    gF->haplotypeString1[j] = hapChar1;
    uint64_t hapChar2 = getMLHapChar(characterProbsHap2, params, maxProbRootChar);
    gF->haplotypeString2[j] = hapChar2;

    // Calculate haplotype probabilities
    gF->haplotypeProbs1[j] = exp(getHaplotypeProb(characterProbsHap1[hapChar1],
            hapChar1, rootCharacterProbsHap2, params, rProbs) - logColumnProbSum);
    gF->haplotypeProbs2[j] = exp(getHaplotypeProb(characterProbsHap2[hapChar2],
            hapChar2, rootCharacterProbsHap1, params, rProbs) - logColumnProbSum);

    // Get combined genotype
    gF->genotypeString[j] = hapChar1 < hapChar2 ? hapChar1 * ALPHABET_SIZE + hapChar2 :
            hapChar2 * ALPHABET_SIZE + hapChar1;

    // Calculate genotype posterior probability
    double genotypeProb = ST_MATH_LOG_ZERO;
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        genotypeProb = stMath_logAdd(genotypeProb,
                characterProbsHap1[hapChar1] + characterProbsHap2[hapChar2] +
                *getSubstitutionProbSlow(params->hetSubModelSlow, i, hapChar1) +
                *getSubstitutionProbSlow(params->hetSubModelSlow, i, hapChar2) +
                invertScaleToLogIntegerSubMatrix(rProbs[i]));
    }
    gF->genotypeProbs[j] = exp(genotypeProb - logColumnProbSum);
}

void fillInPredictedGenome(stGenomeFragment *gF, stRPCell *cell,
        stRPColumn *column, stReferencePriorProbs *referencePriorProbs, stRPHmmParameters *params) {
    /*
     * Computes the most probable haplotype characters / genotypes and associated posterior
     * probabilities for a given interval defined by a cell/column. Fills in these values in the
     * genome fragment argument.
     */

    // Calculate the bit vectors
    uint64_t *bitCountVectors = calculateCountBitVectors(column->seqs, column->depth,
            column->length);

    assert(column->length > 0);

    for(int64_t i=0; i<column->length; i++) {
        fillInPredictedGenomePosition(gF, cell, column, params,
          referencePriorProbs, bitCountVectors, i);
    }

    // Cleanup
    free(bitCountVectors);
}
