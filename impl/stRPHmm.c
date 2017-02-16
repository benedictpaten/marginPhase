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

int stRPHmm_cmpFn(const void *a, const void *b) {
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
        int i = stRPHmm_cmpFn(hmm1, hmm2);
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

stSortedSet *makeComponent(stRPHmm *hmm, stSet *components) {
    /*
     * Create a component containing hmm and add the component to components.
     */
    stSortedSet *component = stSortedSet_construct3(stRPHmm_cmpFn, NULL);
    stSortedSet_insert(component, hmm);
    stSet_insert(components, component);
    return component;
}

stSet *getOverlappingComponents(stList *tilingPath1, stList *tilingPath2) {
    /*
     * Two hmms overlap if their reference coordinate intervals overlaps. The transitive closure of the overlap relation
     * partitions a set of hmms into connected components. This function returns this partition for the hmms in tilingPath1
     * and tilingPath2, each of which is a set of hmms sorted by reference coordinate and which do not overlap in reference
     * coordinates. Each component is a stSortedSet.
     */

    // A map of hmms to components
    stHash *componentsHash = stHash_construct();

    // The set of components
    stSet *components = stSet_construct2((void (*)(void *))stSortedSet_destruct);

    // The "lagging" index of the hmm in tilingPath2 that could possibly overlap hmm1
    int64_t j=0;

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
                    component = stSet_search(components, hmm2);

                    // If hmm2 has no component make one
                    if(component == NULL) {
                        component = makeComponent(hmm2, components);
                        stHash_insert(componentsHash, hmm2, component);
                    }

                    // Add hmm1 to the component
                    stSortedSet_insert(component, hmm1);
                    stHash_insert(componentsHash, hmm1, component);
                }
                // Otherwise component is defined
                else {
                    // Add hmm2 to the component
                    assert(stHash_search(componentsHash, hmm2) == NULL); // Impossible to be defined, as implies that two
                    // hmms in tilingPath2 each both overlap two hmms in tilingPath1.
                    stSortedSet_insert(component, hmm2);
                    stHash_insert(componentsHash, hmm2, component);
                }
            }
            // Else hmm1 and hmm2 do not overlap
            else {
                assert(stRPHmm_cmpFn(hmm1, hmm2) != 0); // Not equal, obviously

                // If hmm1 occurs before hmm2 in the reference ordering
                if(stRPHmm_cmpFn(hmm1, hmm2) < 0) {

                    // If has no component, make a trivial component containing just hmm1 (it doesn't overlap with
                    // any other hmm)
                    if(component == NULL) {
                        makeComponent(hmm1, components);
                    }

                    // Done with hmm1
                    break;
                }
                // else hmm2 occurs after hmm1 in the reference ordering
                else {

                    // Add hmm2 to a trivial component if it does not overlap an HMM in tiling path1
                    if(stHash_search(componentsHash, hmm2) == NULL) {
                        makeComponent(hmm2, components);
                    }

                    // Increase the lagging index as hmm1 and proceding hmms can not overlap with hmm2
                    j++;
                }

            }
        }
    }

    // For any remaining hmms in tilingPath2 that have not been placed in a component
    // put them in a component
    while(j < stList_length(tilingPath2)) {
        stRPHmm *hmm2 = stList_get(tilingPath2, j++);
        if(stHash_search(componentsHash, hmm2) == NULL) {
            makeComponent(hmm2, components);
        }
    }

    // Cleanup
    stHash_destruct(componentsHash);

    return components;
}

stList *getTilingPaths(stSortedSet *hmms) {
    /*
     * Takes set of hmms ordered by reference coordinate (see stRPHmm_cmpFn) and returns a list of
     * tiling paths. Each tiling path consisting of maximal sequences of hmms that do not overlap.
     * Destroys sortedSet in the process.
     */
    stList *tilingPaths = stList_construct();
    while(stSortedSet_size(hmms) > 0) {

        // Make an empty tiling path and add to set of tiling paths built so far
        stList *tilingPath = stList_construct();
        stList_append(tilingPaths, tilingPath);

        // Get the hmm with lowest reference coordinate and add to the tiling path
        stRPHmm *hmm = stSortedSet_getFirst(hmms);
        stList_append(tilingPath, hmm);
        stSortedSet_remove(hmms, hmm);

        // While it exists, get the next closest non-overlapping hmm and add to the tiling path progressively
        while((hmm = getNextClosestNonoverlappingHmm(hmm, hmms)) != NULL) {
            stList_append(tilingPath, hmm);
        }
    }
    stSortedSet_destruct(hmms);

    return tilingPaths;
}

stRPHmm *fuseTilingPath(stList *tilingPath) {
    /*
     * Fuse together the hmms in the tiling path into one hmm. Destroys the tiling path and cleans it up.
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

stList *mergeTwoTilingPaths(stList *tilingPath1, stList *tilingPath2, double posteriorProbabilityThreshold) {
    /*
     *  Takes two lists tilingPath1 and tilingPath2, each of which is a set of hmms ordered by reference coordinates and
     *  non-overlapping in reference coordinates.
     *  Merges together the hmms and returns a single tiling path as a result in the same format as the input lists.
     *  Destroys the input tilingPaths in the process and cleans it up.
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

        // Make two sub tiling paths (there can only be two maximal paths, by definition)
        stList *tilingPaths = getTilingPaths(component);

        assert(stList_length(tilingPaths) == 2);

        stList *subTilingPath1 = stList_get(tilingPaths, 0);
        stList *subTilingPath2 = stList_get(tilingPaths, 1);

        // Fuse the hmms in each sub tiling path
        stRPHmm *hmm1 = fuseTilingPath(subTilingPath1);
        stRPHmm *hmm2 = fuseTilingPath(subTilingPath2);

        // Align
        stRPHmm_alignColumns(hmm1, hmm2);

        // Merge
        stRPHmm *hmm = stRPHmm_createCrossProductHmm(hmm1, hmm2);

        // Prune
        stRPHmm_prune(hmm, posteriorProbabilityThreshold);

        // Add to output tiling path
        stList_append(newTilingPath, hmm);
    }

    //Cleanup
    stList_destruct(componentsList);
    stSet_destruct(components);

    return newTilingPath;
}

stList *mergeTilingPaths(stList *tilingPaths, double posteriorProbabilityThreshold) {
    /*
     * Like mergeTwoTilingPaths(), except instead of just two tiling paths it takes a list.
     */

    // If no tiling paths in input warn and return an empty tiling path
    if(stList_length(tilingPaths) == 0) {
        st_logCritical("WARNING: Zero tiling paths to merge\n");
        return stList_construct();
    }

    // If only one tiling path in the input, the output is just the single input tiling path
    if(stList_length(tilingPaths) == 1) {
        return stList_get(tilingPaths, 0);
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
        tilingPath1 = mergeTilingPaths(tilingPaths1, posteriorProbabilityThreshold);
        stList_destruct(tilingPaths1);

        // Recursively turn the other half of the tiling paths into the other tiling path
        stList *tilingPaths2 = stList_construct();
        for(int64_t i=stList_length(tilingPaths)/2; i < stList_length(tilingPaths); i++) {
            stList_append(tilingPaths2, stList_get(tilingPaths, i));
        }
        tilingPath2 = mergeTilingPaths(tilingPaths2, posteriorProbabilityThreshold);
        stList_destruct(tilingPaths2);
    }
    // Otherwise the number of tiling paths is two
    else {
        tilingPath1 = stList_get(tilingPaths, 0);
        tilingPath2 = stList_get(tilingPaths, 1);
    }

    // Merge together the two tiling paths and return result
    assert(tilingPath1 != NULL);
    assert(tilingPath2 != NULL);
    return mergeTwoTilingPaths(tilingPath1, tilingPath2, posteriorProbabilityThreshold);
}

stList *getRPHmms(stList *profileSeqs, double posteriorProbabilityThreshold) {
    /*
     * Takes a set of profile sequences (stProfileSeq) and returns a list of read partitioning
     * hmms (stRPHmm) ordered and non-overlapping in reference coordinates
     */

    // Create a read partitioning HMM for every sequence and put in ordered set, ordered by reference coordinate
    stSortedSet *readHmms = stSortedSet_construct3(stRPHmm_cmpFn, NULL);
    for(int64_t i=0; i<stList_length(profileSeqs); i++) {
        stSortedSet_insert(readHmms, stRPHmm_construct(stList_get(profileSeqs, i)));
    }

    // Organise HMMs into "tiling paths" consisting of sequences of hmms that do not overlap
    stList *tilingPaths = getTilingPaths(readHmms);

    // Merge together the tiling paths into one merged tiling path, merging the individual hmms when
    // they overlap on the reference
    stList *finalTilingPath = mergeTilingPaths(tilingPaths, posteriorProbabilityThreshold);

    // Cleanup
    stList_destruct(tilingPaths);

    return finalTilingPath;
}

/*
 * Functions for the read partitioning hmm object stRPHmm.
 */

stRPHmm *stRPHmm_construct(stProfileSeq *profileSeq) {
    /*
     * Create a read partitioning HMM representing the single sequence profile.
     */

    stRPHmm *hmm = st_malloc(sizeof(stRPHmm));

    //  Set reference coordinates
    hmm->referenceName = stString_copy(profileSeq->referenceName);
    hmm->refStart = profileSeq->refStart;
    hmm->refLength = profileSeq->length;

    // Add the single profile sequence to the list of the hmm's sequences
    hmm->profileSeqs = stList_construct();
    stList_append(hmm->profileSeqs, profileSeq);

    hmm->columnNumber = 1; // The number of columns in the model, initially just 1
    hmm->maxDepth = 1; // The maximum number of states in a column, initially just 1

    // Create the first column of the model
    stProfileProb **seqs = st_malloc(sizeof(stProfileProb *));
    seqs[0] = profileSeq->profileProbs;
    stRPColumn *column = stRPColumn_construct(hmm->refStart, hmm->refLength, 1, seqs);;
    hmm->firstColumn = column;
    hmm->lastColumn = column;

    return hmm;
}

void stRPHmm_destruct(stRPHmm *hmm) {
    /*
     * Free memory owned by the hmm, including columns.
     */
    free(hmm->referenceName);
    stList_destruct(hmm->profileSeqs);

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

    free(hmm);
}

stRPHmm *stRPHmm_fuse(stRPHmm *leftHmm, stRPHmm *rightHmm) {
    /*
     * Fuses together two hmms, such that leftHmm and rightHMM are on the same reference sequence and non-overlapping and
     * left hmm preceds right hmm on the reference sequence.
     * Returns fused hmm, destroys input hmms in the process.
     */

    // Checks
    if(!stString_eq(leftHmm->referenceName, rightHmm->referenceName)) {
        st_errAbort("Attemping to fuse two hmms not on the same reference sequence");
    }

    if(stRPHmm_overlapOnReference(leftHmm, rightHmm)) {
        st_errAbort("Attemping to fuse two hmms that overlap in reference coordinates");
    }

    if(leftHmm->refStart >= rightHmm->refStart) {
        st_errAbort("Left hmm does not precede right hmm in reference coordinates");
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
    hmm->columnNumber = leftHmm->columnNumber + rightHmm->columnNumber + 1;

    // Make columns to fuse left hmm and right hmm's columns
    stRPMergeColumn *mColumn = stRPMergeColumn_construct(0, 0);
    leftHmm->lastColumn->nColumn = mColumn;
    mColumn->pColumn = leftHmm->lastColumn;
    int64_t gapLength = rightHmm->refStart - leftHmm->refStart + leftHmm->refLength;
    assert(gapLength >= 0);
    if(gapLength > 0) {
        stRPColumn *column = stRPColumn_construct(leftHmm->refStart + leftHmm->refLength, gapLength, 0, NULL);
        mColumn->nColumn = column;
        column->pColumn = mColumn;
        mColumn = stRPMergeColumn_construct(0, 0);
        column->nColumn = mColumn;
        mColumn->pColumn = column;
    }
    mColumn->nColumn = rightHmm->firstColumn;
    rightHmm->firstColumn->pColumn = mColumn;

    // Initialise first/last columns of fused hmm
    hmm->firstColumn = leftHmm->firstColumn;
    hmm->lastColumn = rightHmm->lastColumn;

    // Cleanup
    stRPHmm_destruct(leftHmm);
    stRPHmm_destruct(rightHmm);

    return hmm;
}

void stRPHmm_alignColumns(stRPHmm *hmm1, stRPHmm *hmm2) {
    /*
     * Align the input hmms, modifying them in place, so that they each
     *  (1) span the same reference interval,
     *  (2) have the same number of columns, and (3) so that for all i column i in each
     *  model span the same interval.
     */

    // If the two hmms don't overlap in reference space then complain
    if(!stRPHmm_overlapOnReference(hmm1, hmm2)) {
        st_errAbort("Attempting to align two HMMs that do not overlap in reference coordinate space");
    }

    // If hmm1 starts after hmm2 then call the other way around
    if(hmm1->refStart > hmm2->refStart) {
        stRPHmm_alignColumns(hmm2, hmm1);
        return;
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
        stRPHmm_alignColumns(hmm2, hmm1);
        return;
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

        // There are no more columns, so break
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
    /*
     *  For two aligned hmms (see stRPHmm_alignColumns) returns a new hmm that represents the
     *  cross product of all the states of the two input hmms.
     */

    // If the two hmms have not been previously aligned
    if(stRPHmm_cmpFn(hmm1, hmm2) != 0 || hmm1->columnNumber != hmm2->columnNumber) {
        st_errAbort("Trying to create cross product of two unalignd HMMs");
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
    hmm->columnNumber = hmm1->columnNumber;

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
        stProfileProb **seqs = st_malloc(sizeof(stProfileProb *) * column1->depth+column2->depth);
        memcpy(seqs, column1->seqs, sizeof(stProfileProb *) * column1->depth);
        memcpy(&seqs[column1->depth], column2->seqs, sizeof(stProfileProb *) * column2->depth);
        stRPColumn *column = stRPColumn_construct(column1->refStart, column1->length,
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
        mColumn = stRPMergeColumn_construct((mColumn1->maskFrom <
                mColumn2->pColumn->depth) | mColumn2->maskFrom,
                (mColumn1->maskTo < mColumn2->pColumn->depth) | mColumn2->maskTo);
        mColumn->pColumn = column;

        // Create cross product of merged columns
        stHashIterator *cellIt1 = stHash_getIterator(mColumn1->mergeCellsFrom);
        stRPMergeCell *mCell1;
        while((mCell1 = stHash_getNext(cellIt1)) != NULL) {
            stHashIterator *cellIt2 = stHash_getIterator(mColumn2->mergeCellsFrom);
            stRPMergeCell *mCell2;
            while((mCell2 = stHash_getNext(cellIt2)) != NULL) {
                stRPMergeColumn_addCell(mColumn,
       (mCell1->fromPartition < mColumn2->pColumn->depth) | mCell2->fromPartition,
       (mCell1->toPartition < mColumn2->nColumn->depth) | mCell2->toPartition);
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
    hmm->forwardProbability = LOG_ZERO;

    // Iterate through columns from first to last
    stRPColumn *column = hmm->firstColumn;
    while(1) {
        // Initialize column probability
        column->forwardProb = LOG_ZERO;

        // Initialise cells in the column
        stRPCell *cell = column->head;
        do {
            cell->forwardProb = LOG_ZER0;
        } while((cell = cell->nCell) != NULL);

        if(column->nColumn != NULL) {

            // Initialise cells in the next merge column
            stList *mergeCells = stHash_getValues(column->nColumn->mergeCellsFrom);
            for(int64_t i=0; i<stList_length(mergeCells); i++) {
                stRPMergeCell *mergeCell = stList_get(mergeCells, i);
                mergeCell->forwardProb = LOG_ZERO;
            }
            stList_destruct(mergeCells);

            column = column->nColumn->nColumn;
        }
        else {
            break;
        }
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
        // Iterate through states in column
        stRPCell *cell = column->head;
        do {
            // If the previous merge column exists then propagate forward probability from merge state
            if(column->pColumn != NULL) {
                stRPMergeCell *mCell = stRPMergeColumn_getPreviousMergeCell(cell, column->pColumn);
                cell->forwardProb = mCell->forwardProb;
            }
            // Otherwise initialize probability with log(1.0)
            else {
                cell->forwardProb = 0.0;
            }

            // Emission prob
            cell->forwardProb += hmm->emissionProbability(column, cell);

            // If the next merge column exists then propagate forward probability to the merge state
            if(column->nColumn != NULL) {
                // Add to the next merge cell
                stRPMergeCell *mCell = stRPMergeColumn_getNextMergeCell(cell, column->nColumn);
                mCell->forwardProb = logAdd(cell->forwardProb, mCell->forwardProb);
            }
            else {
                hmm->forwardProbability = logAdd(hmm->forwardProbability, cell->forwardProb);
            }

            // Add to column forward probability
            column->forwardProb = logAdd(column->forwardProb, cell->forwardProb);
        }
        while((cell = cell->nCell) != NULL);

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
    // Initialise total backward probability
    hmm->backwardProbability = LOG_ZERO;

    // Iterate through columns
    stRPColumn *column = hmm->firstColumn;
    while(1) {
        // Initialize column probability
        column->backwardProb = LOG_ZERO;

        // Initialise cells in the column
        stRPCell *cell = column->head;
        do {
            cell->backwardProb = LOG_ZER0;
        } while((cell = cell->nCell) != NULL);

        if(column->nColumn != NULL) {

            // Initialise cells in the next merge column
            stList *mergeCells = stHash_getValues(column->nColumn->mergeCellsFrom);
            for(int64_t i=0; i<stList_length(mergeCells); i++) {
                stRPMergeCell *mergeCell = stList_get(mergeCells, i);
                mergeCell->backwardProb = LOG_ZERO;
            }
            stList_destruct(mergeCells);

            column = column->nColumn->nColumn;
        }
        else {
            break;
        }
    }
}

void stRPHmm_backward(stRPHmm *hmm) {
    /*
     * Backward algorithm for hmm.
     */
    stRPColumn *column = hmm->lastColumn;

    // Initialise backward probabilites
    stRPHmm_initialiseBackwardProbs(hmm);

    // Iterate through columns from last to first
    while(1) {
        // Iterate through states in column
        stRPCell *cell = column->head;
        do {
            // If the next merge column exists then propagate backward probability from merge state
            if(column->nColumn != NULL) {
                stRPMergeCell *mCell = stRPMergeColumn_getNextMergeCell(cell, column->nColumn);
                cell->backwardProb = mCell->backwardProb;
            }
            // Otherwise initialize probability with log(1.0)
            else {
                cell->backwardProb = 0.0;
            }

            // Total backward prob to propagate
            double backwardProb = cell->backwardProb + hmm->emissionProbability(column, cell);

            // If the previous merge column exists then propagate backward probability to the merge state
            if(column->pColumn != NULL) {
                // Add to the previous merge cell
                stRPMergeCell *mCell = stRPMergeColumn_getPreviousMergeCell(cell, column->pColumn);
                mCell->backwardProb = logAdd(cell->backwardProb, mCell->backwardProb);
            }
            else {
                hmm->backwardProbability = logAdd(hmm->backwardProbability, backwardProb);
            }

            // Add to column backward probability
            column->backwardProb = logAdd(column->forwardProb, backwardProb);
        }
        while((cell = cell->nCell) != NULL);

        if(column->pColumn == NULL) {
            break;
        }
        column = column->pColumn->pColumn;
    }
}

void stRPHmm_prune(stRPHmm *hmm, double posteriorProbabilityThreshold) {
    /*
     * Remove cells from hmm whos posterior probability is below the given threshold
     */

    // For each column
    stRPColumn *column = hmm->firstColumn;
    while(1) {

        // For each state
        stRPCell *cell = column->head;
        stRPCell **pCell = &(cell); // Pointer to previous cell, used to remove cells from the linked list
        do {
            // If the posterior probability is below the given threshold
            if(stRPCell_posteriorProb(cell, column) < posteriorProbabilityThreshold) {
                // Remove the state from the linked list of states
                *pCell = cell->nCell;

                // Cleanup
                stRPCell_destruct(cell, 0);
            }
        } while((cell = cell->nCell) != NULL);

        // Move on to the next merge column
        stRPMergeColumn *mColumn = column->nColumn;

        if(mColumn == NULL) {
            break;
        }

        //  For each merge state
        stList *mergeCells = stHash_getValues(mColumn->mergeCellsFrom);
        for(int64_t i=0; i<stList_length(mergeCells); i++) {
            stRPMergeCell *mCell = stList_get(mergeCells, i);

            // If the merge state has posterior probability below the given threshold
            if(stRPMergeCell_posteriorProb(mCell, mColumn) < posteriorProbabilityThreshold) {
                // Remove the state from the merge column
                stHash_remove(mColumn->mergeCellsFrom, &(mCell->fromPartition));
                stHash_remove(mColumn->mergeCellsTo, &(mCell->toPartition));

                // Cleanup
                stRPMergeCell_destruct(mCell);
            }
        }
        stList_destruct(mergeCells);

        column = mColumn->nColumn;
    }
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

stRPColumn *stRPColumn_construct(int64_t refStart, int64_t length, int64_t depth, stProfileProb **seqs) {

}

void stRPColumn_destruct(stRPColumn *column) {

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

/*
 * Read partitioning hmm state (stRPCell) functions
 */

stRPCell *stRPCell_construct(int64_t partition) {

}

void stRPCell_destruct(stRPCell *cell, bool) {

}

double stRPCell_posteriorProb(stRPCell *cell, stRPColumn *column) {
    /*
     *
     */
}

/*
 * Read partitioning hmm merge column (stRPMergeColumn) functions
 */

stRPMergeCell *stRPMergeCell_construct(int32_t fromPartition, int32_t toPartition, stRPMergeColumn *mColumn) {
    stRPMergeCell *mCell = st_malloc(sizeof(stRPMergeCell));
    mCell->fromPartition = fromPartition;
    mCell->toPartition = toPartition;
    stHash_insert(mColumn->mergeCellsFrom, &mCell->fromPartition, mCell);
    stHash_insert(mColumn->mergeCellsTo, &mCell->toPartition, mCell);
    return mCell;
}

void stRPMergeColumn_destruct(stRPMergeColumn *mColumn) {

}

stRPMergeCell *stRPMergeColumn_getNextMergeCell(stRPCell *cell, stRPMergeColumn *mergeColumn) {
    /*
     *
     */
}

stRPMergeCell *stRPMergeColumn_getPreviousMergeCell(stRPCell *cell, stRPMergeColumn *mergeColumn) {
    /*
     *
     */
}

