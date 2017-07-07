/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "stRPHmm.h"

/*
 * Functions for manipulating read partitions described in binary
 */

inline uint64_t makeAcceptMask(uint64_t depth) {
    /*
     * Returns a mask to the given sequence depth that includes all the sequences
     */
    assert(depth <= MAX_READ_PARTITIONING_DEPTH);
    return depth < 64 ? ~(0xFFFFFFFFFFFFFFFF << depth) : 0xFFFFFFFFFFFFFFFF;
}

inline uint64_t mergePartitionsOrMasks(uint64_t partition1, uint64_t partition2,
        uint64_t depthOfPartition1, uint64_t depthOfPartition2) {
    /*
     * Take two read partitions or masks and merge them together
     */
    assert(depthOfPartition1 + depthOfPartition2 <= MAX_READ_PARTITIONING_DEPTH);
    return (partition2 << depthOfPartition1) | partition1;
}

inline uint64_t maskPartition(uint64_t partition, uint64_t mask) {
    /*
     * Mask a read partition
     */
    return partition & mask;
}

inline uint64_t invertPartition(uint64_t partition, uint64_t depth) {
    /*
     * Invert a partition
     */
    return makeAcceptMask(depth) & ~partition;
}

inline bool seqInHap1(uint64_t partition, int64_t seqIndex) {
    /*
     * Returns non-zero if the sequence indexed by seqIndex is in the first haplotype,
     * rather than the second, according to the given partition.
     */
    assert(seqIndex < MAX_READ_PARTITIONING_DEPTH);
    return (partition >> seqIndex) & 1;
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
