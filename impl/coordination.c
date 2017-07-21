/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "stRPHmm.h"

// OpenMP
#if defined(_OPENMP)
#include <omp.h>
#define CELL_BUFFER_SIZE 1000
#endif

/*
 * Functions to create a set of read partitioning HMMs that include a given input set of reads.
 */

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
        if(hmm1->refCoords[hmm1->length-1] < hmm2->refStart) {
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
        // TODO delete
//        st_logInfo("Merging %d\t", i);

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

        // Recursively turn the other half of the tiling paths into the other tiling path
        stList *tilingPaths2 = stList_construct();
        for(int64_t i=stList_length(tilingPaths)/2; i < stList_length(tilingPaths); i++) {
            stList_append(tilingPaths2, stList_get(tilingPaths, i));
        }

#if defined(_OPENMP)
#pragma omp parallel
{
#pragma omp sections nowait
{
#pragma omp section
        tilingPath1 = mergeTilingPaths(tilingPaths1);

#pragma omp section
        tilingPath2 = mergeTilingPaths(tilingPaths2);

}
}
#else
        tilingPath1 = mergeTilingPaths(tilingPaths1);
        tilingPath2 = mergeTilingPaths(tilingPaths2);
#endif
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

stList *getRPHmms(stList *profileSeqs, stHash *referenceNamesToReferencePriors, stRPHmmParameters *params) {
    /*
     * Takes a set of profile sequences (stProfileSeq) and returns a list of read partitioning
     * hmms (stRPHmm) ordered and non-overlapping in reference coordinates.
     * referenceNamesToReferencePriors is a map from reference sequence names to corresponding
     * stReferencePriorProbs objects.
     */
    // Create a read partitioning HMM for every sequence and put in ordered set, ordered by reference coordinate
    stSortedSet *readHmms = stSortedSet_construct3(stRPHmm_cmpFn, (void (*)(void *))stRPHmm_destruct2);
    for(int64_t i=0; i<stList_length(profileSeqs); i++) {
        stProfileSeq *pSeq = stList_get(profileSeqs, i);
        stReferencePriorProbs *referencePriorProbs = stHash_search(referenceNamesToReferencePriors, pSeq->referenceName);
        assert(referencePriorProbs != NULL);
        stRPHmm *hmm = stRPHmm_construct(pSeq, referencePriorProbs, params);
        stSortedSet_insert(readHmms, hmm);
    }
    assert(stSortedSet_size(readHmms) == stList_length(profileSeqs));

    // Organise HMMs into "tiling paths" consisting of sequences of hmms that do not overlap
    stList *tilingPaths = getTilingPaths(readHmms);

    if (stList_length(tilingPaths) > params->maxCoverageDepth)
        st_logInfo("  > Eliminating %d extra tiling paths to reduce max depth\n", stList_length(tilingPaths) - params->maxCoverageDepth);

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
