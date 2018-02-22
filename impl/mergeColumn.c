/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "stRPHmm.h"

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
    assert(popcount64(fromPartition) == popcount64(toPartition));
    assert(popcount64(mColumn->maskFrom) == popcount64(mColumn->maskTo));
    assert(popcount64(fromPartition) <= popcount64(mColumn->maskFrom));

    stRPMergeCell *mCell = st_calloc(1, sizeof(stRPMergeCell));
    mCell->fromPartition = fromPartition;
    mCell->toPartition = toPartition;
    assert(stHash_search(mColumn->mergeCellsFrom, &mCell->fromPartition) == NULL);
    stHash_insert(mColumn->mergeCellsFrom, &mCell->fromPartition, mCell);
    assert(stHash_search(mColumn->mergeCellsTo, &mCell->toPartition) == NULL);
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

    // for debugging/breakpointing purposes
    if (p > 1.001)
        st_errAbort("\nERROR: invalid prob %f", p);
    if (p < 0)
        st_errAbort("\nERROR: invalid prob %f", p);

    assert(p <= 1.001);
    assert(p >= 0.0);
    return p > 1.0 ? 1.0 : p;
}
