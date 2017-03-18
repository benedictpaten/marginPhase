/*
 * stRPHmm.c
 *
 *  Created on: Feb 4, 2017
 *      Author: benedictpaten
 */

#include "stRPHmm.h"
#include "sonLib.h"
#include <float.h>
#include <math.h>

#define ST_MATH_LOG_ZERO -INFINITY
#define ST_MATH_LOG_ONE 0.0

static double logAddP(double a, double b, bool maxNotSum) {
    /*
     * Local function for doing addition of logs or (if doing Viterbi style calculation), to take the max.
     */
    return maxNotSum ? (a > b ? a : b) : stMath_logAdd(a, b);
}

/*
 * Functions for manipulating read partitions described in binary
 */

uint64_t mergePartitionsOrMasks(uint64_t partition1, uint64_t partition2,
        uint64_t depthOfPartition1, uint64_t depthOfPartition2) {
    /*
     * Take two read partitions or masks and merge them together
     */
    assert(depthOfPartition1 + depthOfPartition2 <= MAX_READ_PARTITIONING_DEPTH);
    return (partition2 << depthOfPartition1) | partition1;
}

uint64_t maskPartition(uint64_t partition, uint64_t mask) {
    /*
     * Mask a read partition
     */
    return partition & mask;
}

bool seqInHap1(uint64_t partition, int64_t seqIndex) {
    /*
     * Returns non-zero if the sequence indexed by seqIndex is in the first haplotype,
     * rather than the second, according to the given partition.
     */
    assert(seqIndex < MAX_READ_PARTITIONING_DEPTH);
    return (partition >> seqIndex) & 1;
}

uint64_t makeAcceptMask(int64_t depth) {
    /*
     * Returns a mask to the given sequence depth that includes all the sequences
     */
    assert(depth <= MAX_READ_PARTITIONING_DEPTH);
    return ~(0xFFFFFFFFFFFFFFFF << depth);
}

char * intToBinaryString(uint64_t i) {
    /*
     * Converts the unsigned int to a binary string.
     */
    int64_t bits = sizeof(uint64_t)*8;
    char * str = st_malloc((bits + 1) * sizeof(char));
    str[bits] = '\0'; //terminate the string

    // Decode the bits in from low to high in order
    // so that 14 will end up as 1110 and 15 will end up as
    // 1111 (plus some prefix bits)
    for(int64_t bit=0; bit < bits; i >>= 1) {
        str[(sizeof(uint64_t)*8)-++bit] = i & 1 ? '1' : '0';
    }

    return str;
}

/*
 * Functions to create a set of read partitioning HMMs that include a given input set of reads.
 */

int cmpint(int64_t i, int64_t j) {
    return i > j ? 1 : i < j ? -1 : 0;
}

int stRPHmm_cmpFn(const void *a, const void *b) {
    /*
     * Compares two read partitioning HMMs by coordinate on the reference.
     * Will return equal only if they are the same HMM, with the same memory
     * address, otherwise compares pointers for equal HMMs.
     */
    stRPHmm *hmm1 = (stRPHmm *)a, *hmm2 = (stRPHmm *)b;
    int i = strcmp(hmm1->referenceName, hmm2->referenceName);
    if(i == 0) {
        i = cmpint(hmm1->refStart,  hmm2->refStart);
        if(i == 0) {
            i = cmpint(hmm1->refLength,  hmm2->refLength);
            if(i == 0) {
                i = hmm1 > hmm2 ? 1 : (hmm1 < hmm2 ? -1 : 0);
            }
        }
    }
    return i;
}

stRPHmm *getNextClosestNonoverlappingHmm(stRPHmm *hmm1, stSortedSet *readHmms) {
    /*
     * Returns the HMM from the set readHmms that does not overlap hmm1
     * but whose start coordinate is closest to
     * the end coordinate of hmm1. If does not exist returns NULL.
     */

    // Iterator in the set starting from hmm1
    assert(stSortedSet_search(readHmms, hmm1) == hmm1);
    stSortedSetIterator *it = stSortedSet_getIteratorFrom(readHmms, hmm1);
    stRPHmm *hmm2 = stSortedSet_getNext(it);
    assert(hmm2 == hmm1);

    // For each hmm in readHmms whose coordinate is >= than hmm1's
    while((hmm2 = stSortedSet_getNext(it)) != NULL) {
        // Compare the hmm coordinates just to check that hmm2 has a coordinate >= to hmm1s
        int i = stRPHmm_cmpFn(hmm1, hmm2);
        assert(i <= 0);

        // If hmm1 and hmm2 are on different references, then hmm2 is the closest non-overlapping
        // hmm to hmm1 in reference space
        i = strcmp(hmm1->referenceName, hmm2->referenceName);
        if(i != 0) {
            break;
        }

        // If hmm2 does not overlap hmm1 it must be the closest non-overlapping hmm to hmm1
        if(hmm1->refStart + hmm1->refLength <= hmm2->refStart) {
            break;
        }
    }

    // Cleanup
    stSortedSet_destructIterator(it);

    return hmm2;
}

stSortedSet *makeComponent(stRPHmm *hmm, stSet *components, stHash *componentsHash) {
    /*
     * Create a component containing hmm and add the component to components.
     */
    stSortedSet *component = stSortedSet_construct3(stRPHmm_cmpFn, NULL);
    stSortedSet_insert(component, hmm);
    stSet_insert(components, component);
    assert(stHash_search(componentsHash, hmm) == NULL);
    stHash_insert(componentsHash, hmm, component);
    return component;
}

stSet *getOverlappingComponents(stList *tilingPath1, stList *tilingPath2) {
    /*
     * Two hmms overlap if their reference coordinate intervals overlaps.
     * The transitive closure of the overlap relation
     * partitions a set of hmms into connected components.
     * This function returns this partition for the hmms in tilingPath1
     * and tilingPath2, each of which is a set of hmms sorted by reference
     * coordinate and which do not overlap in reference
     * coordinates. Each component is a stSortedSet.
     */

    // A map of hmms to components
    stHash *componentsHash = stHash_construct();

    // The set of components
    stSet *components = stSet_construct2((void (*)(void *))stSortedSet_destruct);

    // The "lagging" index of the hmm in tilingPath2 that could possibly overlap hmm1
    int64_t j = 0;

    // For each hmm in tilingPath1, in order
    for(int64_t i=0; i<stList_length(tilingPath1); i++) {
        stRPHmm *hmm1 = stList_get(tilingPath1, i);

        // Start with the component being undefined
        stSortedSet *component = NULL;

        // The "leading" index of the hmm in tilingPath2 that could possibly overlap hmm1
        int64_t k = 0;

        // While there exists an hmm in tilingPath2 that precedes or overlaps with hmm1
        while(j+k<stList_length(tilingPath2)) {
            stRPHmm *hmm2 = stList_get(tilingPath2, j+k); // Note the j+k

            // If hmm1 and hmm2 overlap
            if(stRPHmm_overlapOnReference(hmm1, hmm2)) {
                // The leading index is increased
                k++;

                // If component is still NULL
                if(component == NULL) {

                    // Look for a component for hmm2
                    component = stHash_search(componentsHash, hmm2);

                    // If hmm2 has no component make one
                    if(component == NULL) {
                        component = makeComponent(hmm2, components, componentsHash);
                    }

                    // Add hmm1 to the component
                    assert(stSortedSet_search(component, hmm1) == NULL);
                    assert(stHash_search(componentsHash, hmm1) == NULL);
                    stSortedSet_insert(component, hmm1);
                    stHash_insert(componentsHash, hmm1, component);
                }
                // Otherwise component is defined
                else {
                    // Add hmm2 to the component
                    assert(stSortedSet_search(component, hmm2) == NULL);
                    assert(stHash_search(componentsHash, hmm2) == NULL); // Impossible to be defined,
                    // as implies that two
                    // hmms in tilingPath2 each both overlap two hmms in tilingPath1.
                    stSortedSet_insert(component, hmm2);
                    stHash_insert(componentsHash, hmm2, component);
                }
            }
            // Else hmm1 and hmm2 do not overlap
            else {
                // If hmm1 occurs before hmm2 in the reference ordering
                if(stRPHmm_cmpFn(hmm1, hmm2) < 0) {

                    // If has no component, make a trivial component containing just hmm1
                    // (it doesn't overlap with any other hmm)
                    if(component == NULL) {
                        component = makeComponent(hmm1, components, componentsHash);
                    }

                    // Done with hmm1
                    break;
                }
                // else hmm2 occurs before hmm1 in the reference ordering
                else {

                    // Add hmm2 to a trivial component if it does not overlap an HMM in tiling path1
                    if(stHash_search(componentsHash, hmm2) == NULL) {
                        makeComponent(hmm2, components, componentsHash);
                    }

                    // Increase the lagging index as hmm1 and proceding hmms can not overlap with hmm2
                    j++;
                }
            }
        }

        if(component == NULL) {
            //
            assert(stHash_search(componentsHash, hmm1) == NULL);
            makeComponent(hmm1, components, componentsHash);
        }
    }

    // For any remaining hmms in tilingPath2 that have not been placed in a component
    // put them in a component
    while(j < stList_length(tilingPath2)) {
        stRPHmm *hmm2 = stList_get(tilingPath2, j++);
        if(stHash_search(componentsHash, hmm2) == NULL) {
            makeComponent(hmm2, components, componentsHash);
        }
    }

    // Cleanup
    stHash_destruct(componentsHash);

    return components;
}

stList *getTilingPaths(stSortedSet *hmms) {
    /*
     * Takes set of hmms ordered by reference coordinate (see stRPHmm_cmpFn) and returns
     * a list of tiling paths. Each tiling path consisting of maximal sequences of hmms
     * that do not overlap. Destroys sortedSet in the process.
     */
    stList *tilingPaths = stList_construct();
    while(stSortedSet_size(hmms) > 0) {

        // Make an empty tiling path and add to set of tiling paths built so far
        stList *tilingPath = stList_construct();
        stList_append(tilingPaths, tilingPath);

        // Get the hmm with lowest reference coordinate and add to the tiling path
        stRPHmm *hmm = stSortedSet_getFirst(hmms);
        assert(hmm != NULL);
        assert(stSortedSet_search(hmms, hmm) == hmm);
        stList_append(tilingPath, hmm);

        // While it exists, get the next closest non-overlapping hmm
        // and add to the tiling path progressively, removing the chain of hmms from the
        // set of hmms left to tile
        stRPHmm *hmm2;
        while((hmm2 = getNextClosestNonoverlappingHmm(hmm, hmms)) != NULL) {
            stSortedSet_remove(hmms, hmm);
            stList_append(tilingPath, hmm2);
            hmm = hmm2;
            assert(stSortedSet_search(hmms, hmm) == hmm);
        }
        stSortedSet_remove(hmms, hmm);
    }

    // Cleanup the input set
    stSortedSet_destruct(hmms);

    return tilingPaths;
}

stRPHmm *fuseTilingPath(stList *tilingPath) {
    /*
     * Fuse together the hmms in the tiling path into one hmm.
     * Destroys the tiling path and cleans it up.
     */
    stRPHmm *rightHmm = stList_pop(tilingPath);

    // While there remain other hmms in the list fuse them together
    while(stList_length(tilingPath) > 0) {
        stRPHmm *leftHmm = stList_pop(tilingPath);
        rightHmm = stRPHmm_fuse(leftHmm, rightHmm);
    }

    // Cleanup
    stList_destruct(tilingPath);

    return rightHmm;
}

stList *mergeTwoTilingPaths(stList *tilingPath1, stList *tilingPath2) {
    /*
     *  Takes two lists, tilingPath1 and tilingPath2, each of which is a set of hmms
     *  ordered by reference coordinates and
     *  non-overlapping in reference coordinates.
     *  Merges together the hmms and returns a single tiling path as a result in the
     *  same format as the input lists.
     *  Destroys the input tilingPaths in the process and cleans them up.
     */
    // Partition of the hmms into overlapping connected components
    stSet *components = getOverlappingComponents(tilingPath1, tilingPath2);
    // Cleanup the input tiling paths
    stList_destruct(tilingPath1);
    stList_destruct(tilingPath2);

    // The output tiling path, which starts out empty
    stList *newTilingPath = stList_construct();

    // Fuse the hmms

    // For each component of overlapping hmms
    stList *componentsList = stSet_getList(components);
    for(int64_t i=0; i<stList_length(componentsList); i++) {
        stSortedSet *component = stList_get(componentsList, i);
        stSet_remove(components, component);

        // Make two sub-tiling paths (there can only be two maximal paths, by definition)
        stList *tilingPaths = getTilingPaths(component);

        stRPHmm *hmm = NULL;

        if(stList_length(tilingPaths) == 2) {
            stList *subTilingPath1 = stList_get(tilingPaths, 0);
            stList *subTilingPath2 = stList_get(tilingPaths, 1);

            // Fuse the hmms in each sub tiling path
            stRPHmm *hmm1 = fuseTilingPath(subTilingPath1);
            stRPHmm *hmm2 = fuseTilingPath(subTilingPath2);

            // Align
            stRPHmm_alignColumns(hmm1, hmm2);

            // Merge
            hmm = stRPHmm_createCrossProductOfTwoAlignedHmm(hmm1, hmm2);
            stRPHmm_destruct(hmm1, 1);
            stRPHmm_destruct(hmm2, 1);

            // Prune
            stRPHmm_forwardBackward(hmm);
            stRPHmm_prune(hmm);
        }
        else { // Case that component is just one hmm that does not
            // overlap anything else
            assert(stList_length(tilingPaths) == 1);
            stList *subTilingPath1 = stList_get(tilingPaths, 0);
            assert(stList_length(subTilingPath1) == 1);

            hmm = stList_pop(subTilingPath1);
            stList_destruct(subTilingPath1);
        }

        // Add to output tiling path
        stList_append(newTilingPath, hmm);

        stList_destruct(tilingPaths);
    }

    //Cleanup

    stList_destruct(componentsList);
    stSet_destruct(components);

    // Sort new tiling path
    stList_sort(newTilingPath, stRPHmm_cmpFn);

    return newTilingPath;
}

stList *mergeTilingPaths(stList *tilingPaths) {
    /*
     * Like mergeTwoTilingPaths(), except instead of just two tiling paths it takes a list.
     * Destroys the tiling path as it goes.
     */

    // If no tiling paths in input warn and return an empty tiling path
    if(stList_length(tilingPaths) == 0) {
        st_logCritical("WARNING: Zero tiling paths to merge\n");
        stList_destruct(tilingPaths);
        return stList_construct();
    }

    // If only one tiling path in the input, the output is just the single input tiling path
    if(stList_length(tilingPaths) == 1) {
        stList *tilingPath = stList_get(tilingPaths, 0);
        stList_destruct(tilingPaths);
        return tilingPath;
    }

    stList *tilingPath1;
    stList *tilingPath2;

    // If there are more than two tiling paths
    // split the problem into two recursively until there are just two remaining
    // tiling paths
    if(stList_length(tilingPaths) > 2) {

        // Recursively turn the first half of the tiling paths into one tiling path
        stList *tilingPaths1 = stList_construct();
        for(int64_t i=0; i<stList_length(tilingPaths)/2; i++) {
            stList_append(tilingPaths1, stList_get(tilingPaths, i));
        }
        tilingPath1 = mergeTilingPaths(tilingPaths1);

        // Recursively turn the other half of the tiling paths into the other tiling path
        stList *tilingPaths2 = stList_construct();
        for(int64_t i=stList_length(tilingPaths)/2; i < stList_length(tilingPaths); i++) {
            stList_append(tilingPaths2, stList_get(tilingPaths, i));
        }
        tilingPath2 = mergeTilingPaths(tilingPaths2);
    }
    // Otherwise the number of tiling paths is two
    else {
        tilingPath1 = stList_get(tilingPaths, 0);
        tilingPath2 = stList_get(tilingPaths, 1);
    }

    // Merge together the two tiling paths and return result
    assert(tilingPath1 != NULL);
    assert(tilingPath2 != NULL);
    stList_destruct(tilingPaths);
    return mergeTwoTilingPaths(tilingPath1, tilingPath2);
}

stList *getRPHmms(stList *profileSeqs, stRPHmmParameters *params) {
    /*
     * Takes a set of profile sequences (stProfileSeq) and returns a list of read partitioning
     * hmms (stRPHmm) ordered and non-overlapping in reference coordinates.
     *
     * PosteriorProbabilityThreshold is the probability threshold used to keep cells during pruning.
     * MinColumnDepth is the size of a column to need before applying pruning.
     * MaxCoverageDepth is the maximum depth of profileSeqs to allow at any base. If the coverage depth is higher
     * than this then some profile seqs are randomly discarded.
     */
    // Create a read partitioning HMM for every sequence and put in ordered set, ordered by reference coordinate
    stSortedSet *readHmms = stSortedSet_construct3(stRPHmm_cmpFn, NULL);
    for(int64_t i=0; i<stList_length(profileSeqs); i++) {
        stRPHmm *hmm = stRPHmm_construct(stList_get(profileSeqs, i), params);
        stSortedSet_insert(readHmms, hmm);
    }
    assert(stSortedSet_size(readHmms) == stList_length(profileSeqs));

    // Organise HMMs into "tiling paths" consisting of sequences of hmms that do not overlap
    stList *tilingPaths = getTilingPaths(readHmms);

    // Eliminate HMMs that cause the maximum coverage depth to exceed a threshold
    while(stList_length(tilingPaths) > params->maxCoverageDepth) {
        stList *tilingPath = stList_pop(tilingPaths);
        for(int64_t i=0; i<stList_length(tilingPath); i++) {
            stRPHmm_destruct(stList_get(tilingPath, i), 1);
        }
        stList_destruct(tilingPath);
    }

    // Merge together the tiling paths into one merged tiling path, merging the individual hmms when
    // they overlap on the reference
    stList *finalTilingPath = mergeTilingPaths(tilingPaths);
    stList_setDestructor(finalTilingPath, (void (*)(void *))stRPHmm_destruct2);

    return finalTilingPath;
}

/*
 * Functions for profile sequence
 */

stProfileSeq *stProfileSeq_constructEmptyProfile(char *referenceName, int64_t referenceStart,
        int64_t length, int64_t alphabetSize) {
    /*
     * Creates an empty profile sequence, with all the profile probabilities set to 0.
     */
    stProfileSeq *seq = st_malloc(sizeof(stProfileSeq));
    seq->referenceName = stString_copy(referenceName);
    seq->refStart = referenceStart;
    seq->length = length;
    seq->alphabetSize = alphabetSize;
    seq->profileProbs = st_calloc(length*alphabetSize, sizeof(uint8_t));
    return seq;
}

void stProfileSeq_destruct(stProfileSeq *seq) {
    /*
     * Cleans up memory for profile sequence.
     */
    free(seq->profileProbs);
    free(seq->referenceName);
    free(seq);
}

float getProb(uint8_t *p, int64_t characterIndex) {
    /*
     * Gets probability of a given character as a float.
     */
    return ((float)p[characterIndex])/255;
}

void stProfileSeq_print(stProfileSeq *seq, FILE *fileHandle, bool includeProbs) {
    /*
     * Prints a debug representation of a profile sequence.
     */
    char profileString[seq->length+1];
    profileString[seq->length] = '\0';
    for(int64_t i=0; i<seq->length; i++) {
        uint8_t *p = &seq->profileProbs[i * seq->alphabetSize];
        float maxProb = getProb(p, 0);
        int64_t maxChar = 0;
        for(int64_t j=1; j<seq->alphabetSize; j++) {
            float prob = getProb(p, j);
            if(prob < maxProb) {
                maxProb = prob;
                maxChar = j;
            }
        }
        profileString[i] = maxChar + FIRST_ALPHABET_CHAR;
    }
    fprintf(fileHandle, "\tSEQUENCE REF_NAME: %s REF_START %"
            PRIi64 " REF_LENGTH: %" PRIi64 " ML_STRING: %s\n",
            seq->referenceName, seq->refStart, seq->length, profileString);
    if(includeProbs) {
        for(int64_t i=0; i<seq->length; i++) {
            uint8_t *p = &seq->profileProbs[i * seq->alphabetSize];
            // Print individual character probs
            fprintf(fileHandle, "\t\tPOS: %" PRIi64 "", i);
            for(int64_t j=0; j<seq->alphabetSize; j++) {
                fprintf(fileHandle, "\t%f\t", getProb(p, j));
            }
            fprintf(fileHandle, "\n");
        }
    }
}

void printSeqs(FILE *fileHandle, stSet *profileSeqs) {
    /*
     * Prints a set of profile seqs.
     */
    stSetIterator *seqIt = stSet_getIterator(profileSeqs);
    stProfileSeq *pSeq;
    while((pSeq = stSet_getNext(seqIt)) != NULL) {
        stProfileSeq_print(pSeq, fileHandle, 0);
    }
    stSet_destructIterator(seqIt);
}

void printPartition(FILE *fileHandle, stSet *profileSeqs1, stSet *profileSeqs2) {
    /*
     * Print a partition.
     */
    fprintf(fileHandle, "First partition\n");
    printSeqs(fileHandle, profileSeqs1);
    fprintf(fileHandle, "Second partition\n");
    printSeqs(fileHandle, profileSeqs2);
}

/*
 * Character alphabet and substitutions
 */

static inline uint16_t *subProb(uint16_t *matrix, int64_t from, int64_t to, int64_t alphabetSize) {
    return &matrix[from * alphabetSize + to];
}

uint16_t stSubModel_getSubstitutionProb(stSubModel *alphabet, int64_t sourceCharacterIndex,
        int64_t derivedCharacterIndex) {
    /*
     * Gets the (log) substitution probability of getting the derived character given the source (haplotype) character.
     */
    return *subProb(alphabet->logSubMatrix, sourceCharacterIndex, derivedCharacterIndex, alphabet->alphabetSize);
}

uint16_t scaleToLogIntegerSubMatrix(double logProb) {
    /*
     * Convert log probability into scaled form for substitution matrix.
     */
    assert(logProb <= 0);
    if(logProb < -7) {
        st_errAbort("Attempting to set a substitution probability smaller than x=0.0000001 (log(x) = -7)");
    }
    return ALPHABET_MIN_SUBSTITUTION_PROB * (-logProb/7.0);
}

double invertScaleToLogIntegerSubMatrix(int64_t i) {
    /*
     * Invert scaled form to log probability.
     */
    return (7 * ((double)-i))/ALPHABET_MIN_SUBSTITUTION_PROB;
}

void stSubModel_setSubstitutionProb(stSubModel *alphabet, int64_t sourceCharacterIndex,
        int64_t derivedCharacterIndex, double prob) {
    /*
     * Sets the substitution probability, scaling it appropriately by taking the log and then storing as integer (see definition)
     */
    if(prob <= 0 || prob > 1.0) {
        st_errAbort("Attempting to set substitution probability out of 0-1 range");
    }
    *subProb(alphabet->logSubMatrix, sourceCharacterIndex, derivedCharacterIndex, alphabet->alphabetSize) = scaleToLogIntegerSubMatrix(log(prob));
}

stSubModel *stSubModel_constructEmptyModel(int64_t alphabetSize) {
    /*
     * Creates an empty substitution matrix model
     */
    stSubModel *alphabet = st_malloc(sizeof(stSubModel));
    alphabet->alphabetSize = alphabetSize;
    alphabet->logSubMatrix = st_calloc(alphabetSize * alphabetSize, sizeof(uint64_t));
    alphabet->logSubMatrixSlow = st_calloc(alphabetSize * alphabetSize, sizeof(double));

    return alphabet;
}

void stSubModel_destruct(stSubModel *alphabet) {
    free(alphabet->logSubMatrix);
    free(alphabet);
}

/*
 * Emission probabilities
 */

/*
 * Following implement Hamming weight for uint64_t ints, taken from
 * https://en.wikipedia.org/wiki/Hamming_weight
 * TODO: Fiddle with built in popcount instruction
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

uint64_t *retrieveBitCountVector(uint64_t *bitCountVector,
        int64_t position, int64_t characterIndex, int64_t bit, int64_t alphabetSize) {
    /*
     * Returns a pointer to a bit count vector for a given position (offset in the column),
     * character index and bit.
     */
    return &bitCountVector[position * ALPHABET_CHARACTER_BITS * alphabetSize + characterIndex * ALPHABET_CHARACTER_BITS + bit];
}

uint64_t calculateBitCountVector(uint8_t **seqs, int64_t depth,
        int64_t position, int64_t characterIndex, int64_t bit, int64_t alphabetSize) {
    /*
     * Calculates the bit count vector for a given position, character index and bit.
     */
    uint64_t bitCountVector = 0;
    for(int64_t i=0; i<depth; i++) {
        uint8_t *p = &(seqs[i][alphabetSize * position]);
        bitCountVector |= ((((uint64_t)p[characterIndex] >> bit) & 1) << i);
    }

    return bitCountVector;
}

uint64_t *calculateCountBitVectors(uint8_t **seqs, int64_t depth, int64_t length, int64_t alphabetSize) {
    /*
     * Calculates the bit count vector for every position, character and bit in the column.
     */

    // Array of bit vectors, for each position, for each character and for each bit in uint8_t
    uint64_t *bitCountVectors = st_malloc(length * alphabetSize *
            ALPHABET_CHARACTER_BITS * sizeof(uint64_t));

    // For each position
    for(int64_t i=0; i<length; i++) {
        // For each character
        for(int64_t j=0; j<alphabetSize; j++) {
            // For each bit
            for(int64_t k=0; k<ALPHABET_CHARACTER_BITS; k++) {
                *retrieveBitCountVector(bitCountVectors, i, j, k, alphabetSize) =
                        calculateBitCountVector(seqs, depth, i, j, k, alphabetSize);
            }
        }
    }

    return bitCountVectors;
}

inline uint64_t getExpectedInstanceNumber(uint64_t *bitCountVectors, uint64_t depth, uint64_t partition,
        int64_t position, int64_t characterIndex, int64_t alphabetSize) {
    /*
     * Returns the number of instances of a character, given by characterIndex, at the given position within the column for
     * the given partition. Returns value scaled between 0 and ALPHABET_MAX_PROB, where the return value divided by ALPHABET_MAX_PROB
     * is the expected number of instances of the given character in the given subpartition of the column.
     */
    uint64_t expectedCount = 0;
    for(int64_t i=0; i<ALPHABET_CHARACTER_BITS; i++) {
        uint64_t j = *retrieveBitCountVector(bitCountVectors, position, characterIndex, i, alphabetSize);
        expectedCount += (popcount64(j & partition) << i);
    }

    assert(expectedCount >= 0.0);
    assert((double)expectedCount / ALPHABET_MAX_PROB <= depth);
    return expectedCount;
}

uint64_t getLogProbOfReadCharacters(stSubModel *alphabet, uint64_t *expectedInstanceNumbers,
        int64_t sourceCharacterIndex) {
    /*
     * Get the log probability of a given source character given the expected number of instances of
     * each character in the reads.
     */
    uint64_t logCharacterProb = stSubModel_getSubstitutionProb(alphabet, sourceCharacterIndex, 0) *
            expectedInstanceNumbers[0];

    for(int64_t i=1; i<alphabet->alphabetSize; i++) {
        logCharacterProb += stSubModel_getSubstitutionProb(alphabet, sourceCharacterIndex, i) *
                expectedInstanceNumbers[i];
    }

    return logCharacterProb;
}

static uint64_t minP(uint64_t a, uint64_t b) {
    return a < b ? a : b;
}

void columnIndexLogHapProbability(stRPColumn *column, uint64_t index,
        uint64_t partition, uint64_t *bitCountVectors, stRPHmmParameters *params,
        uint64_t *rootCharacterProbs) {
    /*
     * Get the probabilities of the "root" characters for a given read sub-partition and a haplotype.
     */
    // For each possible read character calculate the expected number of instances in the
    // partition and store counts in an array
    uint64_t expectedInstanceNumbers[params->alphabetSize];
    for(int64_t i=0; i<params->alphabetSize; i++) {
        expectedInstanceNumbers[i] = getExpectedInstanceNumber(bitCountVectors,
                               column->depth, partition, index, i, params->alphabetSize);
    }

    // Calculate the probability of the read characters for each possible haplotype character
    uint64_t characterProbsHap[params->alphabetSize];
    for(int64_t i=0; i<params->alphabetSize; i++) {
        characterProbsHap[i] = getLogProbOfReadCharacters(params->readErrorSubModel, expectedInstanceNumbers, i);
    }

    // Calculate the probability of haplotype characters and read characters for each root character
    for(int64_t i=0; i<params->alphabetSize; i++) {
        rootCharacterProbs[i] = characterProbsHap[0] +
                stSubModel_getSubstitutionProb(params->hetSubModel, i, 0) * ALPHABET_MAX_PROB;
        for(int64_t j=1; j<params->alphabetSize; j++) {
            rootCharacterProbs[i] =
                    minP(rootCharacterProbs[i],
                            characterProbsHap[j] + stSubModel_getSubstitutionProb(params->hetSubModel, i, j) * ALPHABET_MAX_PROB);
        }
    }
}

uint64_t columnIndexLogProbability(stRPColumn *column, uint64_t index,
        uint64_t partition, uint64_t *bitCountVectors,
        stRPHmmParameters *params) {
    /*
     * Get the probability of a the characters in a given position within a column for a given partition.
     */
    // Get the sum of log probabilities of the derived characters over the possible source characters
    uint64_t rootCharacterProbsHap1[params->alphabetSize];
    columnIndexLogHapProbability(column, index,
            partition, bitCountVectors, params, rootCharacterProbsHap1);
    uint64_t rootCharacterProbsHap2[params->alphabetSize];
    columnIndexLogHapProbability(column, index,
            ~partition, bitCountVectors, params, rootCharacterProbsHap2);

    // Combine the probabilities to calculate the overall probability of a given position in a column
    uint64_t logColumnProb = rootCharacterProbsHap1[0] + rootCharacterProbsHap2[0];
    for(int64_t i=1; i<params->alphabetSize; i++) {
        logColumnProb = minP(logColumnProb, rootCharacterProbsHap1[i] + rootCharacterProbsHap2[i]);
    }

    return logColumnProb;
}

double emissionLogProbability(stRPColumn *column,
        stRPCell *cell, uint64_t *bitCountVectors,
        stRPHmmParameters *params) {
    /*
     * Get the log probability of a set of reads for a given column.
     */
    assert(column->length > 0);
    uint64_t logPartitionProb = columnIndexLogProbability(column, 0,
            cell->partition, bitCountVectors, params);
    for(int64_t i=1; i<column->length; i++) {
        logPartitionProb += columnIndexLogProbability(column, i,
                cell->partition, bitCountVectors, params);
    }
    return invertScaleToLogIntegerSubMatrix(logPartitionProb)/ALPHABET_MAX_PROB;
}

/*
 * Functions for the read partitioning hmm object stRPHmm.
 */

stRPHmmParameters *stRPHmmParameters_construct(stSubModel *hetSubModel,
        stSubModel *readErrorSubModel,
        bool maxNotSumEmissions,
        bool maxNotSumTransitions,
        int64_t maxPartitionsInAColumn,
        int64_t maxCoverageDepth) {
    /*
     * Create an parameters object for an HMM.
     */
    stRPHmmParameters *params = st_malloc(sizeof(stRPHmmParameters));

    params->alphabetSize = hetSubModel->alphabetSize;
    params->hetSubModel = hetSubModel;
    params->readErrorSubModel = readErrorSubModel;
    params->maxNotSumEmissions = maxNotSumEmissions;
    params->maxNotSumTransitions = maxNotSumTransitions;
    params->maxPartitionsInAColumn = maxPartitionsInAColumn;
    params->maxCoverageDepth = maxCoverageDepth;

    // Checks
    if(params->maxCoverageDepth > MAX_READ_PARTITIONING_DEPTH) {
            st_errAbort("The maximum coverage depth %" PRIi64 " is greater than the maximum allowed by the model: %"
                    PRIi64 "\n", params->maxCoverageDepth, MAX_READ_PARTITIONING_DEPTH);
    }
    if(hetSubModel->alphabetSize != readErrorSubModel->alphabetSize) {
        st_errAbort("Read and haplotype substitution models have different alphabet sizes");
    }
    if(maxPartitionsInAColumn <= 0) {
        st_errAbort("Minimum number of partitions in a column is not a positive integer");
    }

    return params;
}

void stRPHmmParameters_destruct(stRPHmmParameters *params) {
    stSubModel_destruct(params->hetSubModel);
    stSubModel_destruct(params->readErrorSubModel);
    free(params);
}

stRPHmm *stRPHmm_construct(stProfileSeq *profileSeq, stRPHmmParameters *params) {
    /*
     * Create a read partitioning HMM representing the single sequence profile.
     */

    stRPHmm *hmm = st_calloc(1, sizeof(stRPHmm));

    //  Set reference coordinates
    hmm->referenceName = stString_copy(profileSeq->referenceName);
    hmm->refStart = profileSeq->refStart;
    hmm->refLength = profileSeq->length;

    // Add the single profile sequence to the list of the hmm's sequences
    hmm->profileSeqs = stList_construct();
    stList_append(hmm->profileSeqs, profileSeq);

    hmm->parameters = params; // Parameters for the model for computation, this is shared by different HMMs

    hmm->columnNumber = 1; // The number of columns in the model, initially just 1
    hmm->maxDepth = 1; // The maximum number of states in a column, initially just 1

    // Create the first column of the model
    stProfileSeq **seqHeaders = st_malloc(sizeof(stProfileSeq *));
    seqHeaders[0] = profileSeq;
    uint8_t **seqs = st_malloc(sizeof(uint8_t *));
    seqs[0] = profileSeq->profileProbs;
    stRPColumn *column = stRPColumn_construct(hmm->refStart, hmm->refLength, 1, seqHeaders, seqs);
    hmm->firstColumn = column;
    hmm->lastColumn = column;

    // Add two cells to the column to represent the two possible partitions of the single profile sequence
    stRPCell *cell = stRPCell_construct(1);
    column->head = cell;
    cell->nCell = stRPCell_construct(0);

    return hmm;
}

void stRPHmm_destruct(stRPHmm *hmm, bool destructColumns) {
    /*
     * Free memory owned by the hmm, including columns.
     */
    free(hmm->referenceName);
    stList_destruct(hmm->profileSeqs);

    if(destructColumns) {
        // Cleanup the columns of the hmm
        stRPColumn *column = hmm->firstColumn;
        while(1) {
            stRPMergeColumn *mColumn = column->nColumn;
            stRPColumn_destruct(column);
            if(mColumn == NULL) {
                break;
            }
            column = mColumn->nColumn;
            stRPMergeColumn_destruct(mColumn);
        }
    }

    free(hmm);
}

void stRPHmm_destruct2(stRPHmm *hmm) {
    /*
     * Cleans up hmm and columns
     */
    stRPHmm_destruct(hmm, 1);
}

stList *stRPHmm_forwardTraceBack(stRPHmm *hmm) {
    /*
     * Traces back through the forward matrix picking the most probable path.
     * (yes, this is non-symmetric)
     * Returns the result as a list of cells, one from each column.
     */
    stList *path = stList_construct();

    stRPColumn *column = hmm->lastColumn;

    // Pick cell in the last column with highest probability
    stRPCell *cell = column->head;
    double maxProb = cell->forwardLogProb;
    stRPCell *maxCell = cell;
    while((cell = cell->nCell) != NULL) {
        if(cell->forwardLogProb > maxProb) {
            maxProb = cell->forwardLogProb;
            maxCell = cell;
        }
    }

    stList_append(path, maxCell); // Add chosen cell to output

    // Walk back through previous columns
    while(column->pColumn != NULL) {
        // Get previous merge cell
        stRPMergeCell *mCell = stRPMergeColumn_getPreviousMergeCell(maxCell, column->pColumn);
        assert(mCell != NULL);

        // Switch to previous column
        column = column->pColumn->pColumn;

        // Walk through cells in the previous column to find the one with the
        // highest forward probability that transitions to maxCell
        cell = column->head;
        maxCell = NULL;
        maxProb = ST_MATH_LOG_ZERO;
        do {
            // If compatible and has greater probability
            if(stRPMergeColumn_getNextMergeCell(cell, column->nColumn) == mCell && cell->forwardLogProb > maxProb) {
                maxProb = cell->forwardLogProb;
                maxCell = cell;
            }
        } while((cell = cell->nCell) != NULL);

        assert(maxCell != NULL);
        stList_append(path, maxCell);
    }

    stList_reverse(path); // So cells go in order

    return path;
}

stSet *stRPHmm_partitionSequencesByStatePath(stRPHmm *hmm, stList *path, bool partition1) {
    /*
     * For an hmm and path through the hmm (e.g. computed with stRPHmm_forwardTraceBack) returns the
     * set of sequences in the hmm that are predicted to come from one given haplotype.
     */

    stSet *seqsInHap1 = stSet_construct();

    // For each cell/column pair
    stRPColumn *column = hmm->firstColumn;
    for(int64_t i=0; i<stList_length(path); i++) {
        stRPCell *cell = stList_get(path, i);

        // Get sequences in first or second partition
        for(int64_t j=0; j<column->depth; j++) {
            if((seqInHap1(cell->partition, j) && partition1) ||
                    (!seqInHap1(cell->partition, j) && !partition1)) {
                stSet_insert(seqsInHap1, column->seqHeaders[j]);
            }
        }

        if(column->nColumn != NULL) {
            column = column->nColumn->nColumn;
        }
    }

    return seqsInHap1;
}

void stRPHmm_print(stRPHmm *hmm, FILE *fileHandle, bool includeColumns, bool includeCells) {
    /*
     * Prints a debug friendly representation of the state of an hmm.
     */
    //Header line
    fprintf(fileHandle, "HMM REF_NAME: %s REF_START: %" PRIi64 " REF_LENGTH %" PRIi64
            " COLUMN_NUMBER %" PRIi64 " MAX_DEPTH: %" PRIi64 " FORWARD_PROB: %f BACKWARD_PROB: %f\n,",
            hmm->referenceName, hmm->refStart, hmm->refLength,
            hmm->columnNumber, hmm->maxDepth,
            (float)hmm->forwardLogProb, (float)hmm->backwardLogProb);

    if(includeColumns) {
        stRPColumn *column = hmm->firstColumn;
        int64_t i=0;
        while(1) {
            fprintf(fileHandle, "Column %" PRIi64 "\n", i++);

            // Print the column
            stRPColumn_print(column, fileHandle, includeCells);

            if(column->nColumn == NULL) {
                break;
            }

            // Print the merge column
            stRPMergeColumn_print(column->nColumn, fileHandle, includeCells);

            column = column->nColumn->nColumn;
        }
    }
}

stRPHmm *stRPHmm_fuse(stRPHmm *leftHmm, stRPHmm *rightHmm) {
    /*
     * Fuses together two hmms, such that leftHmm and rightHMM are on the same reference sequence and non-overlapping and
     * left hmm preceds right hmm on the reference sequence.
     * Returns fused hmm, destroys input hmms in the process.
     */

    // Checks
    if(!stString_eq(leftHmm->referenceName, rightHmm->referenceName)) {
        st_errAbort("Attempting to fuse two hmms not on the same reference sequence");
    }

    if(stRPHmm_overlapOnReference(leftHmm, rightHmm)) {
        st_errAbort("Attemping to fuse two hmms that overlap in reference coordinates");
    }

    if(leftHmm->refStart >= rightHmm->refStart) {
        st_errAbort("Left hmm does not precede right hmm in reference coordinates for merge");
    }

    // Create a new empty hmm
    stRPHmm *hmm = st_malloc(sizeof(stRPHmm));
    // Set the reference interval
    hmm->referenceName = stString_copy(leftHmm->referenceName);
    hmm->refStart = leftHmm->refStart;
    hmm->refLength = rightHmm->refStart + rightHmm->refLength - leftHmm->refStart;
    // Create the combined list of profile seqs
    hmm->profileSeqs = stList_copy(leftHmm->profileSeqs, NULL);
    stList_appendAll(hmm->profileSeqs, rightHmm->profileSeqs);
    // Set column number
    hmm->columnNumber = leftHmm->columnNumber + rightHmm->columnNumber;
    // Max depth
    hmm->maxDepth = leftHmm->maxDepth > rightHmm->maxDepth ? leftHmm->maxDepth : rightHmm->maxDepth;
    // Parameters
    if(leftHmm->parameters != rightHmm->parameters) {
        st_errAbort("HMM parameters differ in fuse function, panic.");
    }
    hmm->parameters = leftHmm->parameters;

    // Make columns to fuse left hmm and right hmm's columns
    stRPMergeColumn *mColumn = stRPMergeColumn_construct(0, 0);

    // Links
    leftHmm->lastColumn->nColumn = mColumn;
    mColumn->pColumn = leftHmm->lastColumn;

    // Add merge cell to connect the cells in the two columns
    stRPMergeCell_construct(0, 0, mColumn);

    int64_t gapLength = rightHmm->refStart - (leftHmm->refStart + leftHmm->refLength);
    assert(gapLength >= 0);
    if(gapLength > 0) {
        // Make column in the gap
        stRPColumn *column = stRPColumn_construct(leftHmm->refStart + leftHmm->refLength,
                gapLength, 0, NULL, NULL);

        // Links
        mColumn->nColumn = column;
        column->pColumn = mColumn;
        // Make cell for empty column
        column->head = stRPCell_construct(0);

        // Add right merge column
        mColumn = stRPMergeColumn_construct(0, 0);

        // Add merge cell to connect the cells in the two columns
        stRPMergeCell_construct(0, 0, mColumn);

        // Links
        column->nColumn = mColumn;
        mColumn->pColumn = column;
        // Increase the column number to account for the introduced gap column
        hmm->columnNumber += 1;
    }
    mColumn->nColumn = rightHmm->firstColumn;
    rightHmm->firstColumn->pColumn = mColumn;

    // Initialise first/last columns of fused hmm
    hmm->firstColumn = leftHmm->firstColumn;
    hmm->lastColumn = rightHmm->lastColumn;

    // Cleanup
    stRPHmm_destruct(leftHmm, 0);
    stRPHmm_destruct(rightHmm, 0);

    return hmm;
}

void stRPHmm_alignColumns(stRPHmm *hmm1, stRPHmm *hmm2) {
    /*
     * Align the input hmms, modifying them in place, so that they each
     *  (1) span the same reference interval,
     *  (2) have the same number of columns, and
     *  (3) so that for all i, column i in each model span the same interval.
     */
    assert(hmm1 != hmm2);

    // If the two hmms don't overlap in reference space then complain
    if(!stRPHmm_overlapOnReference(hmm1, hmm2)) {
        st_errAbort("Attempting to align two HMMs that do not overlap in reference coordinate space");
    }

    // If hmm1 starts after hmm2 then call the other way around
    if(hmm1->refStart > hmm2->refStart) {
        stRPHmm_alignColumns(hmm2, hmm1);
        return;
    }

    // If hmm1 starts before hmm2 add an empty prefix interval to hmm2
    // so they have the same start coordinate
    if(hmm1->refStart < hmm2->refStart) {
        // Create column
        stRPColumn *column = stRPColumn_construct(hmm1->refStart,
                hmm2->refStart - hmm1->refStart,
                0, NULL, NULL);
        // Add cell
        column->head = stRPCell_construct(0);
        // Create merge column
        stRPMergeColumn *mColumn = stRPMergeColumn_construct(0,0);
        // Add merge cell
        stRPMergeCell_construct(0, 0, mColumn);
        // Create links
        hmm2->firstColumn->pColumn = mColumn;
        mColumn->nColumn = hmm2->firstColumn;
        mColumn->pColumn = column;
        column->nColumn = mColumn;
        assert(column->pColumn == NULL);
        hmm2->firstColumn = column;
        //Adjust start and length of hmm2 interval
        hmm2->refLength += hmm2->refStart - hmm1->refStart;
        hmm2->refStart = hmm1->refStart;
        // Increase column number
        hmm2->columnNumber++;
    }

    // If hmm1 has a shorter reference interval length than hmm2 then call the function
    // with the hmms reversed.
    if(hmm1->refLength < hmm2->refLength) {
        stRPHmm_alignColumns(hmm2, hmm1);
        return;
    }

    // If hmm1 has a longer reference interval than hmm2 append an empty suffix
    // interval to hmm2 to make them the same length.
    if(hmm1->refLength > hmm2->refLength) {
        // Create column
        stRPColumn *column = stRPColumn_construct(hmm2->lastColumn->refStart + hmm2->lastColumn->length,
                hmm1->refLength - hmm2->refLength, 0, NULL, NULL);
        // Add cell
        column->head = stRPCell_construct(0);
        // Create merge column
        stRPMergeColumn *mColumn = stRPMergeColumn_construct(0, 0);
        // Add merge cell
        stRPMergeCell_construct(0, 0, mColumn);
        // Create links
        hmm2->lastColumn->nColumn = mColumn;
        mColumn->pColumn = hmm2->lastColumn;
        mColumn->nColumn = column;
        column->pColumn = mColumn;
        assert(column->nColumn == NULL);
        hmm2->lastColumn = column;
        //Adjust start and length of hmm2 interval
        hmm2->refLength = hmm1->refLength;
        // Increase column number
        hmm2->columnNumber++;
    }

    // Quick coordinate checks
    assert(hmm1->refStart == hmm2->refStart);
    assert(hmm1->refLength == hmm2->refLength);
    assert(hmm1->firstColumn->refStart == hmm1->refStart);
    assert(hmm2->firstColumn->refStart == hmm2->refStart);
    assert(hmm1->lastColumn->refStart + hmm1->lastColumn->length == hmm1->refStart + hmm1->refLength);
    assert(hmm2->lastColumn->refStart + hmm2->lastColumn->length == hmm2->refStart + hmm2->refLength);

    // At this point both hmms have the same reference interval

    // While one hmm has a shorter reference interval than the other split the other interval
    // otherwise move on to the next
    stRPColumn *column1 = hmm1->firstColumn;
    stRPColumn *column2 = hmm2->firstColumn;
    while(1) {
        assert(column1->refStart == column2->refStart);

        if(column1->length > column2->length) {
            stRPColumn_split(column1, column2->length, hmm1);
            assert(column1->nColumn->nColumn->refStart == column1->refStart + column2->length);
        }
        else if(column1->length < column2->length) {
            stRPColumn_split(column2, column1->length, hmm2);
        }

        assert(column1->refStart == column2->refStart);
        assert(column1->length == column2->length); // Now have equal length/start

        // There are no more columns, so break
        if(column1->nColumn == NULL) {
            assert(hmm1->lastColumn == column1);
            assert(column2->nColumn == NULL);
            assert(hmm2->lastColumn == column2);
            break;
        }

        column1 = column1->nColumn->nColumn;
        assert(column2->nColumn != NULL);
        column2 = column2->nColumn->nColumn;
        assert(column1 != NULL);
        assert(column2 != NULL);
    }

    assert(hmm1->columnNumber == hmm2->columnNumber);
}

stRPHmm *stRPHmm_createCrossProductOfTwoAlignedHmm(stRPHmm *hmm1, stRPHmm *hmm2) {
    /*
     *  For two aligned hmms (see stRPHmm_alignColumns) returns a new hmm that represents the
     *  cross product of all the states of the two input hmms.
     */

    // Do sanity checks that the two hmms have been aligned
    if(!stString_eq(hmm1->referenceName, hmm2->referenceName)) {
        st_errAbort("Trying to create cross product of two HMMs "
                "on different reference sequences");
    }
    if(hmm1->refStart != hmm2->refStart) {
        st_errAbort("Trying to create cross product of two HMMs "
                "with different reference interval starts");
    }
    if(hmm1->refLength != hmm2->refLength) {
        st_errAbort("Trying to create cross product of two HMMs "
                "with different reference interval length");
    }
    if(hmm1->columnNumber != hmm2->columnNumber) {
        st_errAbort("Trying to create cross product of two HMMs "
                "with different column numbers");
    }

    // Create a new empty hmm
    stRPHmm *hmm = st_calloc(1, sizeof(stRPHmm));
    // Set the reference interval
    hmm->referenceName = stString_copy(hmm1->referenceName);
    hmm->refStart = hmm1->refStart;
    hmm->refLength = hmm1->refLength;
    // Create the combined list of profile seqs
    hmm->profileSeqs = stList_copy(hmm1->profileSeqs, NULL);
    stList_appendAll(hmm->profileSeqs, hmm2->profileSeqs);
    // Set column number
    hmm->columnNumber = hmm1->columnNumber;
    // Set substitution matrices
    if(hmm1->parameters != hmm1->parameters) {
        st_errAbort("Hmm parameters differ in fuse function, panic.");
    }
    hmm->parameters = hmm1->parameters;

    // For each pair of corresponding columns
    stRPColumn *column1 = hmm1->firstColumn;
    stRPColumn *column2 = hmm2->firstColumn;
    assert(column1 != NULL);
    assert(column2 != NULL);
    stRPMergeColumn *mColumn = NULL;

    while(1) {
        // Check columns aligned
        assert(column1->refStart == column2->refStart);
        assert(column1->length == column2->length);

        // Create the new column

        // Depth
        int64_t newColumnDepth = column1->depth+column2->depth;
        if(newColumnDepth > hmm->maxDepth) {
            hmm->maxDepth = newColumnDepth;
        }

        // Seq headers
        stProfileSeq **seqHeaders = st_malloc(sizeof(stProfileSeq *) * newColumnDepth);
        memcpy(seqHeaders, column1->seqHeaders, sizeof(stProfileSeq *) * column1->depth);
        memcpy(&seqHeaders[column1->depth], column2->seqHeaders, sizeof(stProfileSeq *) * column2->depth);

        // Profiles
        uint8_t **seqs = st_malloc(sizeof(uint8_t *) * newColumnDepth);
        memcpy(seqs, column1->seqs, sizeof(uint8_t *) * column1->depth);
        memcpy(&seqs[column1->depth], column2->seqs, sizeof(uint8_t *) * column2->depth);

        stRPColumn *column = stRPColumn_construct(column1->refStart, column1->length,
                newColumnDepth, seqHeaders, seqs);

        // If the there is a previous column
        if(mColumn != NULL) {
            mColumn->nColumn = column;
            column->pColumn = mColumn;
        }
        else {
            hmm->firstColumn = column;
            assert(column->pColumn == NULL);
        }

        // Create cross product of columns
        stRPCell **pCell = &column->head;
        stRPCell *cell1 = column1->head;
        do {
            stRPCell *cell2 = column2->head;
            do {
                stRPCell *cell = stRPCell_construct(mergePartitionsOrMasks(cell1->partition, cell2->partition,
                        column1->depth, column2->depth));
                // Link cells
                *pCell = cell;
                pCell = &cell->nCell;
            } while((cell2 = cell2->nCell) != NULL);
        } while((cell1 = cell1->nCell) != NULL);

        // Get the next merged column
        stRPMergeColumn *mColumn1 = column1->nColumn;
        stRPMergeColumn *mColumn2 = column2->nColumn;

        // If column is NULL, we have reached the last column
        // and we can exit
        if(mColumn1 == NULL) {
            assert(mColumn2 == NULL);
            assert(hmm1->lastColumn == column1);
            assert(hmm2->lastColumn == column2);

            // Set the last column pointer
            hmm->lastColumn = column;
            break;
        }

        // Create new merged column
        uint64_t fromMask = mergePartitionsOrMasks(mColumn1->maskFrom, mColumn2->maskFrom,
                mColumn1->pColumn->depth, mColumn2->pColumn->depth);
        uint64_t toMask = mergePartitionsOrMasks(mColumn1->maskTo, mColumn2->maskTo,
                        mColumn1->nColumn->depth, mColumn2->nColumn->depth);
        mColumn = stRPMergeColumn_construct(fromMask, toMask);

        // Connect links
        mColumn->pColumn = column;
        column->nColumn = mColumn;

        // Create cross product of merged columns
        stHashIterator *cellIt1 = stHash_getIterator(mColumn1->mergeCellsFrom);
        stRPMergeCell *mCell1;
        while((mCell1 = stHash_getNext(cellIt1)) != NULL) {
            stHashIterator *cellIt2 = stHash_getIterator(mColumn2->mergeCellsFrom);
            stRPMergeCell *mCell2;
            while((mCell2 = stHash_getNext(cellIt2)) != NULL) {
                uint64_t fromPartition = mergePartitionsOrMasks(mCell1->fromPartition,
                        mCell2->fromPartition,
                        mColumn1->pColumn->depth, mColumn2->pColumn->depth);

                uint64_t toPartition = mergePartitionsOrMasks(mCell1->toPartition,
                        mCell2->toPartition,
                        mColumn1->nColumn->depth, mColumn2->nColumn->depth);

                stRPMergeCell_construct(fromPartition, toPartition, mColumn);
            }
            stHash_destructIterator(cellIt2);
        }
        stHash_destructIterator(cellIt1);

        // Get next column
        column1 = mColumn1->nColumn;
        column2 = mColumn2->nColumn;
        assert(column1 != NULL);
        assert(column2 != NULL);
    }

    return hmm;
}

void stRPHmm_initialiseForwardProbs(stRPHmm *hmm) {
    /*
     * Initialize the forward matrix.
     */
    // Initialise total forward probability
    hmm->forwardLogProb = ST_MATH_LOG_ZERO;

    // Iterate through columns from first to last
    stRPColumn *column = hmm->firstColumn;
    while(1) {
        // Initialise cells in the column
        stRPCell *cell = column->head;
        do {
            cell->forwardLogProb = ST_MATH_LOG_ZERO;
        } while((cell = cell->nCell) != NULL);

        if(column->nColumn == NULL) {
            break;
        }

        // Initialise cells in the next merge column
        stList *mergeCells = stHash_getValues(column->nColumn->mergeCellsFrom);
        for(int64_t i=0; i<stList_length(mergeCells); i++) {
            stRPMergeCell *mergeCell = stList_get(mergeCells, i);
            mergeCell->forwardLogProb = ST_MATH_LOG_ZERO;
        }
        stList_destruct(mergeCells);

        column = column->nColumn->nColumn;
    }
}

void stRPHmm_forward(stRPHmm *hmm) {
    /*
     * Forward algorithm for hmm.
     */
    // Initialise state values
    stRPHmm_initialiseForwardProbs(hmm);

    stRPColumn *column = hmm->firstColumn;

    // Iterate through columns from first to last
    while(1) {
        // Get the bit count vectors for the column
        uint64_t *bitCountVectors = calculateCountBitVectors(column->seqs, column->depth,
                column->length, hmm->parameters->alphabetSize);

        // Iterate through states in column
        stRPCell *cell = column->head;
        do {
            // If the previous merge column exists then propagate forward probability from merge state
            if(column->pColumn != NULL) {
                stRPMergeCell *mCell = stRPMergeColumn_getPreviousMergeCell(cell, column->pColumn);
                cell->forwardLogProb = mCell->forwardLogProb;
            }
            // Otherwise initialize probability with log(1.0)
            else {
                cell->forwardLogProb = ST_MATH_LOG_ONE;
            }

            // Emission prob
            cell->forwardLogProb += emissionLogProbability(column, cell, bitCountVectors,
                    (stRPHmmParameters *)hmm->parameters);

            // If the next merge column exists then propagate forward probability to the merge state
            if(column->nColumn != NULL) {
                // Add to the next merge cell
                stRPMergeCell *mCell = stRPMergeColumn_getNextMergeCell(cell, column->nColumn);
                mCell->forwardLogProb = logAddP(mCell->forwardLogProb, cell->forwardLogProb, hmm->parameters->maxNotSumTransitions);
            }
            else {
                // Else propagate probability to total forward probability of model
                hmm->forwardLogProb = logAddP(hmm->forwardLogProb, cell->forwardLogProb, hmm->parameters->maxNotSumTransitions);
            }
        }
        while((cell = cell->nCell) != NULL);

        // Cleanup the bit count vectors
        free(bitCountVectors);

        if(column->nColumn == NULL) {
            break;
        }
        column = column->nColumn->nColumn;
    }
}

void stRPHmm_initialiseBackwardProbs(stRPHmm *hmm) {
    /*
     * Initialize the backward matrix.
     */
    // Initialize total backward probability
    hmm->backwardLogProb = ST_MATH_LOG_ZERO;

    // Iterate through columns
    stRPColumn *column = hmm->firstColumn;
    while(1) {
        // Set total log prob
        column->totalLogProb = ST_MATH_LOG_ZERO;

        // Initialize cells in the column
        stRPCell *cell = column->head;
        do {
            cell->backwardLogProb = ST_MATH_LOG_ZERO;
        } while((cell = cell->nCell) != NULL);

        if(column->nColumn == NULL) {
            break;
        }

        // Initialize cells in the next merge column
        stList *mergeCells = stHash_getValues(column->nColumn->mergeCellsFrom);
        for(int64_t i=0; i<stList_length(mergeCells); i++) {
            stRPMergeCell *mergeCell = stList_get(mergeCells, i);
            mergeCell->backwardLogProb = ST_MATH_LOG_ZERO;
        }
        stList_destruct(mergeCells);

        column = column->nColumn->nColumn;
    }
}

void stRPHmm_backward(stRPHmm *hmm) {
    /*
     * Backward algorithm for hmm.
     */
    stRPColumn *column = hmm->lastColumn;

    // Initialize backward probabilities
    stRPHmm_initialiseBackwardProbs(hmm);

    // Iterate through columns from last to first
    while(1) {
        // Get the bit count vectors for the column
        uint64_t *bitCountVectors = calculateCountBitVectors(column->seqs, column->depth,
                column->length, hmm->parameters->alphabetSize);

        // Iterate through states in column
        stRPCell *cell = column->head;
        do {
            // If the next merge column exists then propagate backward probability from merge state
            if(column->nColumn != NULL) {
                stRPMergeCell *mCell = stRPMergeColumn_getNextMergeCell(cell, column->nColumn);
                cell->backwardLogProb = mCell->backwardLogProb;
            }
            // Otherwise initialize probability with log(1.0)
            else {
                cell->backwardLogProb = ST_MATH_LOG_ONE;
            }

            // Add to column total probability
            column->totalLogProb = logAddP(column->totalLogProb, cell->forwardLogProb + cell->backwardLogProb, hmm->parameters->maxNotSumTransitions);

            // Total backward prob to propagate
            double probabilityToPropagateLogProb = cell->backwardLogProb + emissionLogProbability(column, cell,
                    bitCountVectors,  (stRPHmmParameters *)hmm->parameters);

            // If the previous merge column exists then propagate backward probability to the merge state
            if(column->pColumn != NULL) {
                // Add to the previous merge cell
                stRPMergeCell *mCell = stRPMergeColumn_getPreviousMergeCell(cell, column->pColumn);
                mCell->backwardLogProb = logAddP(mCell->backwardLogProb, probabilityToPropagateLogProb,
                            hmm->parameters->maxNotSumTransitions);
            }
            else {
                hmm->backwardLogProb = logAddP(hmm->backwardLogProb, probabilityToPropagateLogProb,
                        hmm->parameters->maxNotSumTransitions);
            }
        }
        while((cell = cell->nCell) != NULL);

        // Cleanup the bit count vectors
        free(bitCountVectors);

        if(column->pColumn == NULL) {
            break;
        }
        column = column->pColumn->pColumn;
    }
}

void stRPHmm_forwardBackward(stRPHmm *hmm) {
    /*
     * Runs the forward and backward algorithms and sets the total column probabilities.
     *
     * (The backward algorithm on its own is now exposed as it requires running the forward
     * algorithm to set the total probabilities within the columns)
     *
     * This function must be run upon an HMM to calculate cell posterior probabilities.
     */
    stRPHmm_forward(hmm);
    stRPHmm_backward(hmm);
}

static int cellCmpFn(const void *a, const void *b, const void *extraArg) {
    /*
     * Sort cells by posterior probability in descending order.
     */
    stRPCell *cell1 = (stRPCell *)a, *cell2 = (stRPCell *)b;
    stRPColumn *column = (stRPColumn *)extraArg;
    double p1 = stRPCell_posteriorProb(cell1, column), p2 = stRPCell_posteriorProb(cell2, column);
    return p1 > p2 ? -1 : p1 < p2 ? 1 : 0;
}

static int mergeCellCmpFn(const void *a, const void *b, const void *extraArg) {
    /*
     * Sort merge cells by posterior probability in descending order.
     */
    stRPMergeCell *cell1 = (stRPMergeCell *)a, *cell2 = (stRPMergeCell *)b;
    stRPMergeColumn *column = (stRPMergeColumn *)extraArg;
    double p1 = stRPMergeCell_posteriorProb(cell1, column), p2 = stRPMergeCell_posteriorProb(cell2, column);
    return p1 > p2 ? -1 : p1 < p2 ? 1 : 0;
}

void filterMergeCells(stRPMergeColumn *mColumn, stSet *chosenMergeCellsSet) {
    /*
     * Removes merge cells from the column that are not in chosenMergeCellsSet
     */
    assert(stSet_size(chosenMergeCellsSet) > 0);
    stList *mergeCells = stHash_getValues(mColumn->mergeCellsFrom);
    for(int64_t i=0; i<stList_length(mergeCells); i++) {
        stRPMergeCell *mCell = stList_get(mergeCells, i);
        assert(mCell != NULL);
        if(stSet_search(chosenMergeCellsSet, mCell) == NULL) {
            // Remove the state from the merge column
            assert(stHash_search(mColumn->mergeCellsFrom, &(mCell->fromPartition)) == mCell);
            assert(stHash_search(mColumn->mergeCellsTo, &(mCell->toPartition)) == mCell);
            stHash_remove(mColumn->mergeCellsFrom, &(mCell->fromPartition));
            stHash_remove(mColumn->mergeCellsTo, &(mCell->toPartition));

            // Cleanup
            stRPMergeCell_destruct(mCell);
        }
    }
    stList_destruct(mergeCells);
    assert(stSet_size(chosenMergeCellsSet) == stHash_size(mColumn->mergeCellsFrom));
    assert(stSet_size(chosenMergeCellsSet) == stHash_size(mColumn->mergeCellsTo));
}

stSet *getLinkedMergeCells(stRPMergeColumn *mColumn,
        stRPMergeCell *(*getNCell)(stRPCell *, stRPMergeColumn *),
        stList *cells) {
    /*
     * Returns the set of merge cells in the column that are linked to a cell
     * in cells.
     */
    stSet *chosenMergeCellsSet = stSet_construct();
    for(int64_t i=0; i<stList_length(cells); i++) {
        stRPMergeCell *mCell = getNCell(stList_get(cells, i), mColumn);
        assert(mCell != NULL);
        stSet_insert(chosenMergeCellsSet, mCell);
    }
    assert(stSet_size(chosenMergeCellsSet) > 0);
    return chosenMergeCellsSet;
}

void relinkCells(stRPColumn *column, stList *cells) {
    /*
     * Re-links the cells in the list 'cells' to make up the list of cells in the column.
     */
    stRPCell **pCell = &column->head; // Pointer to previous cell, used to
    // remove cells from the linked list
    for(int64_t i=0; i<stList_length(cells); i++) {
        stRPCell *cell = stList_get(cells, i);
        *pCell = cell;
        pCell = &cell->nCell;
    }
    *pCell = NULL;
    assert(column->head != NULL);
}

stList *getLinkedCells(stRPColumn *column,
        stRPMergeCell *(*getPCell)(stRPCell *, stRPMergeColumn *),
        stRPMergeColumn *mColumn) {
    /*
     * Returns the set of cells in column that are linked to a cell in mColumn.
     */

    // Put cells into an array and sort by descending posterior prob
    // only keeping cells that still have a preceding merge cell

    stList *cells = stList_construct();
    stRPCell *cell = column->head;
    do {
        if(mColumn == NULL || getPCell(cell, mColumn) != NULL) {
            stList_append(cells, cell);
            cell = cell->nCell;
        }
        else {
            stRPCell *nCell = cell->nCell;
            stRPCell_destruct(cell);
            cell = nCell;
        }
    } while(cell != NULL);
    stList_sort2(cells, cellCmpFn, column);
    assert(stList_length(cells) > 0);

    return cells;
}

void stRPHmm_pruneForwards(stRPHmm *hmm) {
    /*
     * Remove cells from hmm whos posterior probability is below the given threshold
     */

    // For each column
    stRPColumn *column = hmm->firstColumn;
    stRPMergeColumn *mColumn = NULL;

    while(1) {
        assert(column->head != NULL);

        // Get cells that have a valid previous cell
        stList *cells = getLinkedCells(column, stRPMergeColumn_getPreviousMergeCell, mColumn);

        // Get rid of the excess cells
        while(stList_length(cells) > hmm->parameters->maxPartitionsInAColumn) {
            stRPCell_destruct(stList_pop(cells));
        }

        // Relink the cells (from most probable to least probable)
        relinkCells(column, cells);

        // Move on to the next merge column
        mColumn = column->nColumn;

        if(mColumn == NULL) {
            assert(column == hmm->lastColumn);
            stList_destruct(cells);
            break;
        }

        //  Get merge cells that are connected to a cell in the previous column
        stSet *chosenMergeCellsSet = getLinkedMergeCells(mColumn,
                stRPMergeColumn_getNextMergeCell, cells);

        // Shrink the the number of chosen cells to less than equal to the desired number
        stList *chosenMergeCellsList = stSet_getList(chosenMergeCellsSet);
        stList_sort2(chosenMergeCellsList, mergeCellCmpFn, mColumn);
        while(stList_length(chosenMergeCellsList) > hmm->parameters->maxPartitionsInAColumn) {
            stSet_remove(chosenMergeCellsSet, stList_pop(chosenMergeCellsList));
        }
        assert(stList_length(chosenMergeCellsList) == stSet_size(chosenMergeCellsSet));
        stList_destruct(chosenMergeCellsList);

        // Get rid of merge cells we don't need
        filterMergeCells(mColumn, chosenMergeCellsSet);

        // Cleanup
        stList_destruct(cells);
        stSet_destruct(chosenMergeCellsSet);

        column = mColumn->nColumn;
    }
}

void stRPHmm_pruneBackwards(stRPHmm *hmm) {
    /*
     * Remove cells from hmm whos posterior probability is below the given threshold
     */

    // For each column
    stRPColumn *column = hmm->lastColumn;
    stRPMergeColumn *mColumn = NULL;

    while(1) {
        assert(column->head != NULL);

        // Get cells that have a valid previous cell
        stList *cells = getLinkedCells(column, stRPMergeColumn_getNextMergeCell, mColumn);

        // This must be true because the forward pass has already winnowed the number below the
        // threshold
        assert(stList_length(cells) <= hmm->parameters->maxPartitionsInAColumn);

        // Relink the cells (from most probable to least probable)
        relinkCells(column, cells);

        // Move on to the next merge column
        mColumn = column->pColumn;

        if(mColumn == NULL) {
            assert(column == hmm->firstColumn);
            stList_destruct(cells);
            break;
        }

        //  Get merge cells that are connected to a cell in the previous column
        stSet *chosenMergeCellsSet = getLinkedMergeCells(mColumn,
                stRPMergeColumn_getPreviousMergeCell, cells);

        // By the same logic, this number if pruned on the forwards pass
        assert(stSet_size(chosenMergeCellsSet) <= hmm->parameters->maxPartitionsInAColumn);

        // Get rid of merge cells we don't need
        filterMergeCells(mColumn, chosenMergeCellsSet);

        // Cleanup
        stList_destruct(cells);
        stSet_destruct(chosenMergeCellsSet);

        column = mColumn->pColumn;
    }
}

void stRPHmm_prune(stRPHmm *hmm) {
    stRPHmm_pruneForwards(hmm);
    stRPHmm_pruneBackwards(hmm);
}

bool stRPHmm_overlapOnReference(stRPHmm *hmm1, stRPHmm *hmm2) {
    /*
     * Return non-zero iff hmm1 and hmm2 have the same reference sequence and overlapping
     * coordinates intervals on that reference sequence.
     */

    // If either interval is zero length this is not a well defined comparison
    if(hmm1->refLength <= 0 || hmm2->refLength <= 0) {
        st_errAbort("Trying to compare HMMs with a zero length coordinate interval");
    }

    // Check if on the same reference sequence
    if(!stString_eq(hmm1->referenceName, hmm2->referenceName)) {
        return 0;
    }

    // Check if intervals overlap

    // If hmm1 starts after hmm2's start coordinate then switch hmm1 for hmm2
    if(hmm1->refStart > hmm2->refStart) {
        return stRPHmm_overlapOnReference(hmm2, hmm1);
    }

    // The coordinates of the first interval overlap the second
    return hmm1->refStart + hmm1->refLength > hmm2->refStart;
}

/*
 * Read partitioning hmm column (stRPColumn) functions
 */

stRPColumn *stRPColumn_construct(int64_t refStart, int64_t length, int64_t depth,
        stProfileSeq **seqHeaders, uint8_t **seqs) {
    stRPColumn *column = st_calloc(1, sizeof(stRPColumn));

    // Reference coordinates
    column->refStart = refStart;
    column->length = length;

    // Sequences
    column->depth = depth;
    column->seqHeaders = seqHeaders;
    column->seqs = seqs;

    // Initially contains not states
    column->head = NULL;

    return column;
}

void stRPColumn_destruct(stRPColumn *column) {
    // Clean up the contained cells
    stRPCell *cell = column->head;
    while(cell != NULL) {
        stRPCell *pCell = cell;
        cell = cell->nCell;
        stRPCell_destruct(pCell);
    }

    free(column->seqHeaders);
    free(column->seqs);
    free(column);
}

void stRPColumn_print(stRPColumn *column, FILE *fileHandle, bool includeCells) {
    /*
     * Print a description of the column. If includeCells is true then print the
     * state of the cells too.
     */
    fprintf(fileHandle, "\tCOLUMN: REF_START: %" PRIi64
            " REF_LENGTH: %" PRIi64 " DEPTH: %" PRIi64
            " TOTAL_PROB: %f\n",
            column->refStart, column->length, column->depth,
            (float)column->totalLogProb);
    for(int64_t i=0; i<column->depth; i++) {
        stProfileSeq_print(column->seqHeaders[i], fileHandle, 0);
    }
    if(includeCells) {
        stRPCell *cell = column->head;
        do {
            fprintf(fileHandle, "\t\t");
            stRPCell_print(cell, fileHandle);
        } while((cell = cell->nCell) != NULL);
    }
}

void stRPColumn_split(stRPColumn *column, int64_t firstHalfLength, stRPHmm *hmm) {
    /*
     * Split the column into two to divide into two smaller reference intervals
     */

    // Create column
    stProfileSeq **seqHeaders = st_malloc(sizeof(stProfileSeq *) * column->depth);
    memcpy(seqHeaders, column->seqHeaders, sizeof(stProfileSeq *) * column->depth);
    uint8_t **seqs = st_malloc(sizeof(uint8_t *) * column->depth);
    // Update the pointers to the seqs
    for(int64_t i=0; i<column->depth; i++) {
        seqs[i] = &(column->seqs[i][firstHalfLength * seqHeaders[i]->alphabetSize]);
    }
    assert(firstHalfLength > 0); // Non-zero length for first half
    assert(column->length-firstHalfLength > 0); // Non-zero length for second half
    stRPColumn *rColumn = stRPColumn_construct(column->refStart+firstHalfLength,
            column->length-firstHalfLength, column->depth, seqHeaders, seqs);

    // Create merge column
    uint64_t acceptMask = makeAcceptMask(column->depth);
    stRPMergeColumn *mColumn = stRPMergeColumn_construct(acceptMask, acceptMask);

    // Copy cells
    stRPCell *cell = column->head;
    stRPCell **pCell = &(rColumn->head);
    do {
        *pCell = stRPCell_construct(cell->partition);
        stRPMergeCell_construct(cell->partition, cell->partition, mColumn);
        pCell = &((*pCell)->nCell);
    } while((cell = cell->nCell) != NULL);

    // Create links
    rColumn->pColumn = mColumn;
    mColumn->nColumn = rColumn;

    // If is the last column
    if(column->nColumn == NULL) {
       assert(hmm->lastColumn == column);
       hmm->lastColumn = rColumn;
       assert(rColumn->nColumn == NULL);
    }
    else {
        column->nColumn->pColumn = rColumn;
        rColumn->nColumn = column->nColumn;
    }
    column->nColumn = mColumn;
    mColumn->pColumn = column;

    // Increase column number
    hmm->columnNumber++;

    // Adjust length of previous column
    column->length = firstHalfLength;
}

/*
 * Read partitioning hmm state (stRPCell) functions
 */

stRPCell *stRPCell_construct(int64_t partition) {
    stRPCell *cell = st_calloc(1, sizeof(stRPCell));
    cell->partition = partition;
    return cell;
}

void stRPCell_destruct(stRPCell *cell) {
    free(cell);
}

void stRPCell_print(stRPCell *cell, FILE *fileHandle) {
    /*
     * Prints a debug representation of the cell.
     */
    char *partitionString = intToBinaryString(cell->partition);
    fprintf(fileHandle, "CELL PARTITION: %s FORWARD_PROB: %f BACKWARD_PROB: %f\n",
            partitionString, (float)cell->forwardLogProb, (float)cell->backwardLogProb);
    free(partitionString);
}

double stRPCell_posteriorProb(stRPCell *cell, stRPColumn *column) {
    /*
     * Calculate the posterior probability of visiting the given cell. Requires that the
     * forward and backward algorithms have been run.
     */
    double p = exp(cell->forwardLogProb + cell->backwardLogProb - column->totalLogProb);
    assert(p <= 1.1);
    assert(p >= 0.0);
    return p > 1.0 ? 1.0 : p;
}

/*
 * Read partitioning hmm merge column (stRPMergeColumn) functions
 */

static uint64_t intHashFn(const void *a) {
    return *(uint64_t *)a;
}

static int intEqualsFn(const void *key1, const void *key2) {
    return *(uint64_t *)key1 == *(uint64_t *)key2;
}

stRPMergeColumn *stRPMergeColumn_construct(uint64_t maskFrom, uint64_t maskTo) {
    stRPMergeColumn *mColumn = st_calloc(1, sizeof(stRPMergeColumn));
    mColumn->maskFrom = maskFrom;
    mColumn->maskTo = maskTo;

    // Maps between partitions and cells
    mColumn->mergeCellsFrom = stHash_construct3(intHashFn, intEqualsFn, NULL, (void (*)(void *))stRPMergeCell_destruct);
    mColumn->mergeCellsTo = stHash_construct3(intHashFn, intEqualsFn, NULL, NULL);

    return mColumn;
}

void stRPMergeColumn_destruct(stRPMergeColumn *mColumn) {
    stHash_destruct(mColumn->mergeCellsFrom);
    stHash_destruct(mColumn->mergeCellsTo);
    free(mColumn);
}

void stRPMergeColumn_print(stRPMergeColumn *mColumn, FILE *fileHandle, bool includeCells) {
    /*
     * Print a debug representation of the merge column.
     */
    char *maskFromString = intToBinaryString(mColumn->maskFrom);
    char *maskToString = intToBinaryString(mColumn->maskTo);
    fprintf(fileHandle, "\tMERGE_COLUMN MASK_FROM: %s MASK_TO: %s"
            " DEPTH: %" PRIi64 "\n", maskFromString, maskToString,
            stHash_size(mColumn->mergeCellsFrom));
    assert(stHash_size(mColumn->mergeCellsFrom) == stHash_size(mColumn->mergeCellsTo));
    free(maskFromString);
    free(maskToString);
    if(includeCells) {
        stHashIterator *it = stHash_getIterator(mColumn->mergeCellsFrom);
        stRPMergeCell *mCell;
        while((mCell = stHash_getNext(it)) != NULL) {
            fprintf(fileHandle, "\t\t");
            stRPMergeCell_print(mCell, fileHandle);
        }
        stHash_destructIterator(it);
    }
}

stRPMergeCell *stRPMergeColumn_getNextMergeCell(stRPCell *cell, stRPMergeColumn *mergeColumn) {
    /*
     * Get the merge cell that this cell feeds into.
     */
    uint64_t i = maskPartition(cell->partition, mergeColumn->maskFrom);
    stRPMergeCell *mCell = stHash_search(mergeColumn->mergeCellsFrom, &i);
    return mCell;
}

stRPMergeCell *stRPMergeColumn_getPreviousMergeCell(stRPCell *cell, stRPMergeColumn *mergeColumn) {
    /*
     * Get the merge cell that this cell feeds from.
     */
    uint64_t i = maskPartition(cell->partition,  mergeColumn->maskTo);
    stRPMergeCell *mCell = stHash_search(mergeColumn->mergeCellsTo, &i);
    return mCell;
}

int64_t stRPMergeColumn_numberOfPartitions(stRPMergeColumn *mColumn) {
    /*
     * Returns the number of cells in the column.
     */
    return stHash_size(mColumn->mergeCellsFrom);
}

/*
 * Read partitioning hmm merge cell (stRPMergeCell) functions
 */

stRPMergeCell *stRPMergeCell_construct(uint64_t fromPartition, uint64_t toPartition,
        stRPMergeColumn *mColumn) {
    /*
     * Create a merge cell, adding it to the merge column mColumn.
     */
    stRPMergeCell *mCell = st_calloc(1, sizeof(stRPMergeCell));
    mCell->fromPartition = fromPartition;
    mCell->toPartition = toPartition;
    stHash_insert(mColumn->mergeCellsFrom, &mCell->fromPartition, mCell);
    stHash_insert(mColumn->mergeCellsTo, &mCell->toPartition, mCell);
    return mCell;
}

void stRPMergeCell_destruct(stRPMergeCell *mCell) {
    free(mCell);
}

void stRPMergeCell_print(stRPMergeCell *mCell, FILE *fileHandle) {
    /*
     * Prints a debug representation of the cell.
     */
    char *fromPartitionString = intToBinaryString(mCell->fromPartition);
    char *toPartitionString = intToBinaryString(mCell->toPartition);
    fprintf(fileHandle, "MERGE_CELL FROM_PARTITION: %s TO_PARTITION: %s "
            "FORWARD_PROB: %f BACKWARD_PROB: %f\n",
            fromPartitionString, toPartitionString,
            (float)mCell->forwardLogProb, (float)mCell->backwardLogProb);
    free(fromPartitionString);
    free(toPartitionString);
}

double stRPMergeCell_posteriorProb(stRPMergeCell *mCell, stRPMergeColumn *mColumn) {
    /*
     * Calculate the posterior probability of visiting the given cell. Requires that the
     * forward and backward algorithms have been run.
     */
    double p = exp(mCell->forwardLogProb + mCell->backwardLogProb -
           mColumn->nColumn->totalLogProb);
    assert(p <= 1.001);
    assert(p >= 0.0);
    return p > 1.0 ? 1.0 : p;
}
