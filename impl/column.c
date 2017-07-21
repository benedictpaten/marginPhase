/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "stRPHmm.h"

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

    // TODO: can refCoords be shared by hmm or do they need to be calloc'ed separately?
    column->refCoords = st_calloc(length, sizeof(int64_t));
    column->refCoordMap = stHash_construct3(stHash_stringKey, stHash_intPtrEqualKey, NULL, NULL);

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
//    free(column->refCoords);
    //    stHash_destruct(column->refCoordMap);
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
        seqs[i] = &(column->seqs[i][firstHalfLength * ALPHABET_SIZE]);
    }
    assert(firstHalfLength > 0); // Non-zero length for first half
    assert(column->length-firstHalfLength > 0); // Non-zero length for second half

    stRPColumn *rColumn = stRPColumn_construct(column->refCoords[firstHalfLength],
            column->length-firstHalfLength, column->depth, seqHeaders, seqs);

    // Set ref coords for new column
    int64_t *indexes = st_calloc(rColumn->length, sizeof(int64_t));
    for (int64_t i = 0; i < rColumn->length; i++) {
        rColumn->refCoords[i] = column->refCoords[i + firstHalfLength];
        indexes[i] = i;
        if (stHash_search(rColumn->refCoordMap, &rColumn->refCoords[i]) == NULL) {
            stHash_insert(rColumn->refCoordMap, &rColumn->refCoords[i], &indexes[i]);
        }
    }
    rColumn->refEnd = rColumn->refCoords[rColumn->length-1];
    assert(rColumn->refEnd == column->refEnd);

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
    column->refEnd = column->refCoords[column->length-1];
}

stSet *stRPColumn_getColumnSequencesAsSet(stRPColumn *column) {
    /*
     * Get profile sequences in the column as a set.
     */
    stSet *seqSet = stSet_construct();
    for(int64_t i=0; i<column->depth; i++) {
        stSet_insert(seqSet, column->seqHeaders[i]);
    }
    return seqSet;
}

stSet *stRPColumn_getSequencesInCommon(stRPColumn *column1, stRPColumn *column2) {
    /*
     * Returns set (as stSet) of profile sequences shared by both columns.
     */
    stSet *seqSet1 = stRPColumn_getColumnSequencesAsSet(column1);
    stSet *seqSet2 = stRPColumn_getColumnSequencesAsSet(column2);

    stSet *seqsShared = stSet_getIntersection(seqSet1, seqSet2);

    // Cleanup
    stSet_destruct(seqSet1);
    stSet_destruct(seqSet2);

    return seqsShared;
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

