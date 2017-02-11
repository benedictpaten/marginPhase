/*
 * stRPHmm.c
 *
 *  Created on: Feb 4, 2017
 *      Author: benedictpaten
 */

#include "stRPHmm.h"
#include "sonLib.h"

/*
 * Math functions, putting it here for now for purpose of fiddling.
 */

#define logUnderflowThreshold 7.5
#define posteriorMatchThreshold 0.01

static inline double lookup(double x) {
    //return log (exp (x) + 1);
    assert(x >= 0.00f);
    assert(x <= logUnderflowThreshold);
    if (x <= 1.00f)
        return ((-0.009350833524763f * x + 0.130659527668286f) * x + 0.498799810682272f) * x + 0.693203116424741f;
    if (x <= 2.50f)
        return ((-0.014532321752540f * x + 0.139942324101744f) * x + 0.495635523139337f) * x + 0.692140569840976f;
    if (x <= 4.50f)
        return ((-0.004605031767994f * x + 0.063427417320019f) * x + 0.695956496475118f) * x + 0.514272634594009f;
    return ((-0.000458661602210f * x + 0.009695946122598f) * x + 0.930734667215156f) * x + 0.168037164329057f;
}

double logAdd(double x, double y) {
    /*
     * Returns a reasonable approximation to log(exp(x) + exp(y)).
     */
    if (x < y)
        return (x == LOG_ZERO || y - x >= logUnderflowThreshold) ? y : lookup(y - x) + x;
    return (y == LOG_ZERO || x - y >= logUnderflowThreshold) ? x : lookup(x - y) + y;
}

/*
 * Functions to create a set of read partitioning HMMs that include a given input set of reads.
 */

int cmpint(int64_t i, int64_t j) {
    return i > j ? 1 : i < j ? -1 : 0;
}

int readHMMCmpFn(const void *a, const void *b) {
    /*
     * Compares two read partitioning HMMs by coordinate on the reference.
     */
    stRPHmm *hmm1 = (stRPHmm *)a, *hmm2 = (stRPHmm *)b;
    int i = strcmp(hmm1->referenceName, hmm2->referenceName);
    if(i == 0) {
        i = cmpint(hmm1->refStart,  hmm1->refStart);
        if(i == 0) {
            i = cmpint(hmm1->refLength,  hmm1->refLength);
        }
    }
    return i;
}

stRPHmm *getNextClosestNonoverlappingHmm(stRPHmm *hmm1, stSortedSet *readHmms) {
    /*
     * Returns the HMM from the set readHmms that does not overlap hmm1 but whose start coordinate is closest to
     * the end coordinate of hmm1. If does not exist returns NULL.
     */

    // Iterator in the set starting from hmm1
    stSortedSetIterator *it = stSortedSet_getIteratorFrom(readHmms, hmm1);
    stRPHmm *hmm2 = stSortedSet_getNext(it);
    assert(hmm2 == hmm1);

    // For each hmm in readHmms whose coordinate is >= than hmm1's
    while((hmm2 = stSortedSet_getNext(it)) != NULL) {
        // Compare the hmms coordinates just to check that hmm2 has a coordinate >= to hmm1s
        int i = readHMMCmpFn(hmm1, hmm2);
        assert(i >= 0);

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
    return hmm2;
}

stSet *makeComponent(stRPHmm *hmm, stSet *components) {
    /*
     * Create a component containing hmm and add the component to components.
     */
    stSortedSet *component = stSortedSet_construct3(readHMMCmpFn, NULL);
    stSortedSet_insert(component, hmm);
    stSet_insert(components, component);
    return component;
}

stSet *getOverlappingComponents(stList *tilingPath1, stList *tilingPath2) {
    /*
     * Two hmms overlap if their reference coordinate intervals overlaps. The transitive closure of the overlap relation
     * partitions a set of hmms into connected components. This function returns this partition for the hmms in tilingPath1
     * and tilingPath2, each of which is a set of hmms sorted by reference coordinate.
     */

    // A map of hmms to components
    stHash *componentsHash = stHash_construct();

    // The set of components
    stSet *components = stSet_construct();

    //
    int64_t j=0;
    for(int64_t i=0; i<stList_length(tilingPath1); i++) {
        stRPHmm *hmm1 = stList_get(tilingPath1, i);
        stSortedSet *component = NULL;

        //
        int64_t k = 0;
        while(j<stList_length(tilingPath2)) {
            stRPHmm *hmm2 = stList_get(tilingPath2, j+k);

            //
            if(stRPHmm_overlapOnReference(hmm1, hmm2)) {
                //
                k++;

                //
                if(component == NULL) {
                    component = stHash_search(components, hmm2);
                    if(component == NULL) {
                        component = makeComponent(hmm2, components);
                        stHash_insert(components, hmm2, component);
                    }
                    stSortedSet_insert(component, hmm1);
                    stHash_insert(components, hmm1, component);
                }
                //
                else {
                    assert(stHash_search(components, hmm2) == NULL);
                    stSortedSet_insert(component, hmm2);
                    stHash_insert(componentsHash, hmm2, component);
                }
            }
            else {
                //
                assert(readHMMCmpFn(hmm1, hmm2) != 0);

                //
                if(readHMMCmpFn(hmm1, hmm2) < 0) {
                    //
                    if(component == NULL) {
                        makeComponent(hmm1, components);
                    }

                    //
                    break;
                }
                //
                else {
                    // Add hmm2 to a trivial component if it does not overlap an HMM in tiling path1
                    if(stHash_search(componentsHash, hmm2) == NULL) {
                        makeComponent(hmm2, components);
                    }

                    //
                    j++;
                }

            }
        }
    }

    //
    while(j < stList_length(tilingPath2)) {
        stRPHmm *hmm2 = stList_get(tilingPath2, j++);
        if(stHash_search(componentsHash, hmm2) == NULL) {
            makeComponent(hmm2, components);
        }
    }

    //
    stHash_destruct(componentsHash);

    return components;
}

stList *mergeTwoTilingPaths(stList *tilingPath1, stList *tilingPath2) {
    //
    stSet *components = getOverlappingComponents(tilingPath1, tilingPath2);

    //
    stList_destruct(tilingPath1);
    stList_destruct(tilingPath2);

    //
    stList *newTilingPath = stList_construct();

    //
    stSetIterator *componentsIt = stSet_getIterator(components);
    stSortedSet *component;
    while((component = stSet_getNext(componentsIt)) != NULL) {

        //
        stRPHmm *hmm1 = stSortedSet_getFirst(component);
        stSortedSet_remove(component, hmm1);

        //
        while(stSortedSet_size(component) > 0) {
            stRPHmm *hmm2 = stSortedSet_getFirst(component);
            stSortedSet_remove(component, hmm2);

            // Align

            // Merge

            // Prune
        }
    }
    stSet_destructIterator(componentsIt);

    //
    return newTilingPath;
}

stList *mergeTilingPaths(stList *tilingPaths) {
    //
    if(stList_length(tilingPaths) == 0) {
        return NULL;
    }

    //
    if(stList_length(tilingPaths) == 1) {
        return stList_get(tilingPaths, 0);
    }

    //
    stList *tilingPath1;
    stList *tilingPath2;
    if(stList_length(tilingPaths) == 2) {

        //
        tilingPath1 = stList_get(tilingPaths, 0);
        tilingPath2 = stList_get(tilingPaths, 1);
    }
    else {
        //
        stList *tilingPaths1 = stList_construct();
        for(int64_t i=0; i<stList_length(tilingPaths)/2; i++) {
            stList_append(tilingPaths1, stList_get(tilingPaths, i));
        }
        tilingPath1 = mergeTilingPaths(tilingPaths1);
        stList_destruct(tilingPaths1);

        //
        stList *tilingPaths2 = stList_construct();
        for(int64_t i=stList_length(tilingPaths)/2; i < stList_length(tilingPaths); i++) {
            stList_append(tilingPaths2, stList_get(tilingPaths, i));
        }
        tilingPath2 = mergeTilingPaths(tilingPaths2);
        stList_destruct(tilingPaths2);
    }

    //
    assert(tilingPath1 != NULL);
    assert(tilingPath2 != NULL);
    stList *tilingPath = mergeTwoTilingPaths(tilingPath1, tilingPath2);
    stList_destruct(tilingPath1);
    stList_destruct(tilingPath2);
    return tilingPath;
}

stList *getRPHmms(stList *profileSeqs) {
    /*
     *
     */

    // Create a read partitioning HMM for every sequence and put in ordered set, ordered by reference coordinate
    stSortedSet *readHmms = stSortedSet_construct3(readHMMCmpFn, NULL);
    for(int64_t i=0; i<stList_length(profileSeqs); i++) {
        stSortedSet_insert(readHmms, stRPHmm_construct(stList_get(profileSeqs, i)));
    }

    // Organise HMMs into "tiling paths" consisting of sequences of hmms that do not overlap
    stList *tilingPaths = stList_construct();
    while(stSortedSet_size(readHmms) > 0) {

        //
        stList *tilingPath = stList_construct();
        stList_append(tilingPath);

        //
        stRPHmm *hmm = stSortedSet_getFirst(readHmms);
        stList_append(tilingPath, hmm);
        stSortedSet_remove(readHmms, hmm);

        //
        while((hmm = getNextClosestNonoverlappingHmm(hmm, readHmms)) != NULL) {
            stList_append(tilingPath, hmm);
        }
    }
    stSortedSet_destruct(readHmms);

    //
    stList *finalTilingPath = mergeTilingPaths(tilingPaths);

    //
    stList_destruct(tilingPaths);

    return finalTilingPath;
}

/*
 *
 */

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
    if(!stRPHmm_overlapOnReference(hmm1, hmm2)) {
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
