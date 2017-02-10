/*
 * stRPHmm.c
 *
 *  Created on: Feb 4, 2017
 *      Author: benedictpaten
 */

#include "stRPHmm.h"


/* Algorithm overview
 *
 * Create a read partitioning HMM for every sequence.
 * Sort HMMs by start coordinate.
 * Merge HMMs into "tiling paths" consisting of sequences with minimal overlaps
 * While set of tiled HMMs is too large.
 *     Recursivelly split in two
 *     Merge together tiling HMMs
 *     Prune (if necessary)
 */
stList *getRPHmms(stList *profileSeqs) {
}

stRPHmm *stRPHmm_construct(stProfileSeq *profileSeq) {
    //
    stRPHmm *hmm = st_malloc(sizeof(stRPHmm));

    //
    hmm->referenceName = stString_copy(profileSeq->referenceName);
    hmm->refStart = profileSeq->refStart;
    hmm->refLength = profileSeq->length;

    //
    hmm->profileSeqs = stList_construct();
    stList_append(hmm->profileSeqs, profileSeq);
    hmm->profileSeqs[0] = profileSeq;

    //
    hmm->columnNumber = 1;
    hmm->maxDepth = 1;

    //
    stProfileProb **seqs = st_malloc(sizeof(stProfileProb *));
    seqs[0] = profileSeq->profileProbs;
    stRPColumn *column = stRPColumn_construct(hmm->refStart, hmm->refLength, 1, seqs);;
    hmm->firstColumn = column;
    hmm->lastColumn = column;

    return hmm;
}

void stRPHmm_destruct(stRPHmm *hmm) {
    //
    free(hmm->referenceName);
    free(hmm->profileSeqs);

    //
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

    //
    free(hmm);
}

void stRPHmm_alignColumns(stRPHmm *hmm1, stRPHmm *hmm2) {
    // If the two hmms don't overlap in reference space then complain
    if(!stRPHmm_overlapOnReference(stRPHmm *hmm1, stRPHmm *hmm2)) {
        // Raise an exception
    }

    // If hmm1 starts after hmm2 then call the other way around
    if(hmm1->refStart > hmm2->refStart) {
        return stRPHmm_alignColumns(hmm2, hmm1);
    }

    // If hmm1 starts before hmm2 add an emoty prefix interval to hmm2
    // so they have the same start coordinate
    if(hmm1->refStart < hmm2->refStart) {
        // Create column
        stRPColumn *column = stRPColumn_construct(hmm1->refStart, hmm2->refStart - hmm1->refStart, 0, NULL);
        // Add cell
        column->head = stRPCell_construct(0);
        // Create merge column
        stRPMergeColumn *mColumn = stRPMergeColumn_construct(0, INT32_MAX);
        // Add merge cell
        stRPMergeColumn_addCell(mColumn, 0, INT32_MAX);
        // Create links
        hmm2->firstColumn->pColumn = mColumn;
        mColumn->nColumn = hmm2->firstColumn;
        mColumn->pColumn = column;
        column->nColumn = mColumn;
        hmm2->firstColumn = column;
        //Adjust start and length of hmm2 interval
        hmm2->refLength += hmm2->refStart - hmm1->refStart;
        hmm2->refStart = hmm1->refStart;
    }

    // If hmm1 has a shorter reference interval length than hmm2 then call the function
    // with the hmms reversed.
    if(hmm1->refLength < hmm2->refLength) {
        return stRPHmm_alignColumns(hmm2, hmm1);
    }

    // If hmm1 has a longer reference interval than hmm2 append an empty suffix
    // interval to hmm2 to make them the same length.
    if(hmm1->refLength > hmm2->refLength) {
        // Create column
        stRPColumn *column = stRPColumn_construct(hmm1->lastColumn->refStart + hmm1->lastColumn->length, hmm1->refLength - hmm2->refLength, 0, NULL);
        // Add cell
        column->head = stRPCell_construct(0);
        // Create merge column
        stRPMergeColumn *mColumn = stRPMergeColumn_construct(INT32_MAX, 0);
        // Add merge cell
        stRPMergeColumn_addCell(mColumn, INT32_MAX, 0);
        // Create links
        hmm2->lastColumn->nColumn = mColumn;
        mColumn->pColumn = hmm2->lastColumn;
        mColumn->nColumn = column;
        column->pColumn = mColumn;
        hmm2->lastColumn = column;
        //Adjust start and length of hmm2 interval
        hmm2->refLength = hmm1->refLength;
    }

    // At this point both hmms have the same reference interval

    // While one hmm has a shorter reference interval than the other split the other interval
    // otherwise move on to the next
    stRPColumn *column1 = hmm1->firstColumn;
    stRPColumn *column2 = hmm2->firstColumn;
    while(1) {
        assert(column1->refStart == column2->refStart);

        if(column1->length > column2->length) {
            stRPColumn_split(column1, column2->length, hmm1);
        }
        else if(column1->length < column2->length) {
            stRPColumn_split(column2, column1->length, hmm2);
        }

        if(column1->nColumn == NULL) {
            assert(column2->nColumn == NULL);
            break;
        }

        column1 = column1->nColumn->nColumn;
        assert(column2->nColumn != NULL);
        column2 = column2->nColumn->nColumn;
        assert(column1 != NULL);
        assert(column2 != NULL);
    }
}

stRPHmm *stRPHmm_createCrossProductOfTwoAlignedHmm(stRPHmm *hmm1, stRPHmm *hmm2) {
    // If the two hmms don't have the same reference interval/columns
    if() {
        // Raise an exception
    }

    // Create a new empty hmm
    stRPHmm *hmm = st_malloc(sizeof(stRPHmm));
    // Set the reference interval
    hmm->referenceName = stString_copy(hmm1->referenceName);
    hmm->refStart = hmm1->refStart;
    hmm->refLength = hmm1->refLength;
    // Create the combined list of profile seqs
    hmm->profileSeqs = stList_copy(hmm1->profileSeqs, NULL);
    stList_appendAll(hmm->profileSeqs, hmm2->profileSeqs);
    // Set column number
    assert(hmm1->columnNumber == hmm2->columnNumber);
    hmm->columnNumber = hmm1->columnNumber;

    // For each pair of corresponding columns
    stRPColumn *column1 = hmm1->firstColumn;
    stRPColumn *column2 = hmm2->firstColumn;
    assert(column1 != NULL);
    assert(column2 != NULL);
    stRPMergeColumn *mColumn = NULL;

    while(1) {
        // Create the new column
        stProfileProb **seqs = st_malloc(sizeof(stProfileProb *) * column1->depth+column2->depth);
        memcpy(seqs, column1->seqs, sizeof(stProfileProb *) * column1->depth);
        memcpy(&seqs[column1->depth], column2->seqs, sizeof(stProfileProb *) * column2->depth);
        stRPColumn column = stRPColumn_construct(column1->refStart, column1->length,
                column1->depth+column2->depth, seqs);

        // If the there is a previous column
        if(mColumn != NULL) {
            mColumn->nColumn = column;
            column->pColumn = mColumn;
        }
        else {
            hmm->firstColumn = column;
        }

        // Create cross product of columns
        stRPCell **pCell = &(column->head);
        stRPCell *cell1 = column1->head;
        do {
            stRPCell *cell2 = column2->head;
            do {
                stRPCell *cell = stRPCell_construct((cell1->partition < column2->depth) | cell2->partition);
                // Link cells
                *pCell = cell;
                pCell = &cell;
            } while((cell2 = cell2->nCell) != NULL);
        } while((cell1 = cell1->nCell) != NULL);

        // Get the next merged column
        stRPMergeColumn *mColumn1 = column1->nColumn;
        stRPMergeColumn *mColumn2 = column2->nColumn;

        // If column is NULL, we have reached the last column
        // and we can exit
        if(mColumn1 == NULL) {
            assert(mColumn2 == NULL);

            // Set the last column pointer
            hmm->lastColumn = column;
            break;
        }

        // Create merged column
        mColumn = stRPMergeColumn_construct(mColumn1->maskFrom <
                mColumn2->pColumn->depth | mColumn2->maskFrom,
                mColumn1->maskTo < mColumn2->pColumn->depth | mColumn2->maskTo);
        mColumn->pColumn = column;

        // Create cross product of merged columns
        stHashIterator *cellIt1 = stHash_getIterator(mColumn1->mergeCellsFrom);
        stRPMergeCell *mCell1;
        while((mCell1 = stHash_getNext(cellIt1)) != NULL) {
            stHashIterator *cellIt2 = stHash_getIterator(mColumn2->mergeCellsFrom);
            stRPMergeCell *mCell2;
            while((mCell2 = stHash_getNext(cellIt2)) != NULL) {
                stRPMergeColumn_addCell(mColumn,
       mCell1->fromPartition < mColumn2->pColumn->depth | mCell2->fromPartition,
       mCell1->toPartition < mColumn2->nColumn->depth | mCell2->toPartition);
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

void stRPHmm_forward(stRPHmm *hmm) {
    stRPColumn *column = hmm->firstColumn;
    double totalProb = 0.0;

    while(1) {
        //
        stRPMergeColumn *pColumn = column->nColumn;
        stRPMergeColumn *mColumn = column->nColumn;

        //
        stRPCell *cell = column->head;
        do {
            //
            if(pColumn != NULL) {
                stRPMergeCell *mCell = stRPMergeColumn_getPreviousMergeCell(cell, mColumn);
                cell->forwardProb = mCell->forwardProb;
            }
            //
            else {
                cell->forwardProb = 0.0;
            }

            // Emission prob
            // TODO

            //
            if(mColumn != NULL) {
                // Add to the next merge cell
                stRPMergeCell *mCell = stRPMergeColumn_getNextMergeCell(cell, mColumn);

                //
                mCell->p = logAdd(cell->p, mCell->p);
            }
            else {
                totalProb = logAdd(totalProb, cell->p);
            }
        }
        while((cell = cell->nCell) != NULL);

        if(mColumn == NULL) {
            break;
        }
        column = mColumn->nColumn;
    }
}

void stRPHmm_backward(stRPHmm *hmm) {
    stRPColumn *column = hmm->firstColumn;
    double totalProb = 0.0;

    while(1) {
        //
        stRPMergeColumn *pColumn = column->nColumn;
        stRPMergeColumn *mColumn = column->nColumn;

        //
        stRPCell *cell = column->head;
        do {
            //
            if(pColumn != NULL) {
                stRPMergeCell *mCell = stRPMergeColumn_getPreviousMergeCell(cell, mColumn);
                cell->p = mCell->p;
            }
            //
            else {
                cell->p = 0.0;
            }

            // Emission prob
            // TODO

            //
            if(mColumn != NULL) {
                // Add to the next merge cell
                stRPMergeCell *mCell = stRPMergeColumn_getNextMergeCell(cell, mColumn);

                //
                mCell->p = logAdd(cell->p, mCell->p);
            }
            else {
                totalProb = logAdd(totalProb, cell->p);
            }
        }
        while((cell = cell->nCell) != NULL);

        if(mColumn == NULL) {
            break;
        }
        column = mColumn->nColumn;
    }
}

void stRPHmm_prune(stRPHmm *hmm, int64_t samples) {
    // Do the forward and backward matrix calculations.
    stRPHmm_forward(hmm);
    stRPHmm_backward(hmm);

    //
    stSet *visited = stSet_construct();

    //
    for(int64_t i=0; i<samples; i++) {
        stRPHMM_sample(hmm, visited);
    }

    //
    stRPColumn *column = hmm->firstColumn;
    while(1) {

        //
        stRPCell *cell = column->head;
        stRPCell **pCell = &(cell);
        do {
            //
            if(!stSet_search(visited, cell)) {
                //
                *pCell = cell->nCell;

                //
                stRPCell_destruct(cell, 0);
            }
        } while((cell = cell->nCell) != NULL);

        //
        stRPMergeColumn *mColumn = column->nColumn;

        //
        if(mColumn == NULL) {
            break;
        }

        //
        stList *mergeCells = stHash_getValues(mColumn->mergeCellsFrom);
        for(int64_t i=0; i<stList_length(mergeCells); i++) {
            stRPMergeCell *mCell = stList_get(mergeCells, i);

            //
            if(!stSet_search(visited, mCell)) {
                stHash_remove(mColumn->mergeCellsFrom, &(mCell->fromPartition));
                stHash_remove(mColumn->mergeCellsTo, &(mCell->toPartition));
                stRPMergeCell_destruct(mCell);
            }
        }
        stList_destruct(mergeCells);

        //
        column = mColumn->nColumn;
    }

    //
    stSet_destruct(visited);
}

bool stRPHmm_overlapOnReference(stRPHmm *hmm1, stRPHmm *hmm2) {
    // Check if on the same reference sequence
    if(!stString_eq(hmm1->referenceName, hmm2->referenceName)) {
        return 0;
    }

    // Check if intervals overlap

    //
    if(hmm1->refStart > hmm2->refStart) {
        return stRPHmm_overlapOnReference(hmm2, hmm1);
    }

    //
    assert(hmm1->refLength > 0);
    assert(hmm2->refLength > 0);

    //
    return hmm1->refStart + hmm1->refLength >= hmm2->refStart;
}


void stRPColumn_split(stRPColumn *column, int64_t firstHalfLength, stRPHmm *hmm) {
    // Create column
    stProfileProb **seqs = st_malloc(sizeof(stProfileProb *) * column->depth);
    memcpy(seqs, column->seqs, sizeof(stProfileProb *) * column->depth);
    stRPColumn *rColumn = stRPColumn_construct(column->refStart+firstHalfLength, column->length-firstHalfLength, seqs);

    // Create merge column
    stRPMergeColumn *mColumn = stRPMergeColumn_construct(0, 0);

    // Copy cells
    stRPCell *cell = column->head;
    stRPCell **pCell = &(rColumn->head);
    do {
        *pCell = stRPCell_construct(cell->partition);
        stRPMergeColumn_addCell(mColumn, cell->partition, cell->partition);
        pCell = &((*pCell)->nCell);
    } while((cell = cell->nCell) != NULL);

    // Create links
    rColumn->pColumn = mColumn;
    mColumn->nColumn = rColumn;
    // If is the last column
    if(column->nColumn == NULL) {
       hmm->lastColumn = rColumn;
    }
    else {
        column->nColumn->pColumn = rColumn;
        rColumn->nColumn = column->nColumn;
    }
    column->nColumn = mColumn;
    mColumn->pColumn = column;
}

stRPMergeCell *stRPMergeCell_construct(int32_t fromPartition, int32_t toPartition, stRPMergeColumn *mColumn) {
    stRPMergeCell *mCell = st_malloc(sizeof(stRPMergeCell));
    mCell->fromPartition = fromPartition;
    mCell->toPartition = toPartition;
    stHash_insert(mColumn->mergeCellsFrom, &mCell->fromPartition, mCell);
    stHash_insert(mColumn->mergeCellsTo, &mCell->toPartition, mCell);
    return mCell;
}
