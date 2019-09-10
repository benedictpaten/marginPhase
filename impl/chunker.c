/*
 * Copyright (C) 2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <omp.h>
#include "margin.h"

BamChunkRead *bamChunkRead_construct() {
    return bamChunkRead_construct2(NULL, NULL, NULL, TRUE, NULL);
}
BamChunkRead *bamChunkRead_construct2(char *readName, char *nucleotides, uint8_t *qualities, bool forwardStrand,
                                      BamChunk *parent) {
    BamChunkRead *r = malloc(sizeof(BamChunkRead));
    r->readName = readName;
    r->nucleotides = nucleotides;
    r->readLength = (nucleotides == NULL ? 0 : strlen(nucleotides));
    r->qualities = qualities;
    r->forwardStrand = forwardStrand;
    r->parent = parent;

    return r;
}
BamChunkRead *bamChunkRead_constructRLECopy(BamChunkRead  *read, RleString *rle) {
    BamChunkRead *r = st_calloc(1, sizeof(BamChunkRead));
    r->readName = read->readName ==  NULL ? NULL : stString_copy(read->readName);
    r->nucleotides = stString_copy(rle->rleString);
    r->readLength = rle->length;
    r->forwardStrand = read->forwardStrand;
    r->parent = read->parent;
    r->qualities = NULL;

    // calculate read qualities (if set)
    //TODO unit test this
    if (read->qualities != NULL) {
        r->qualities = st_calloc(rle->length, sizeof(uint8_t));
        int64_t rawPos = 0;
        for (int64_t rlePos = 0; rlePos < rle->length; rlePos++) {
            uint8_t min = UINT8_MAX;
            uint8_t max = 0;
            int64_t mean = 0;
            for (int64_t repeatIdx = 0; repeatIdx < rle->repeatCounts[rlePos]; repeatIdx++) {
                uint8_t q = read->qualities[rawPos];
                min = (q < min ? q : min);
                max = (q > max ? q : max);
                mean += q;

                rawPos++;
            }
            mean = mean / rle->repeatCounts[rlePos];
            assert(mean <= UINT8_MAX);
            // pick your favorite metric
            //r->qualities[rlePos] = min;
            //r->qualities[rlePos] = max;
            r->qualities[rlePos] = (uint8_t) mean;
        }
    }

    return r;
}
void bamChunkRead_destruct(BamChunkRead *r) {
    if (r->readName != NULL) free(r->readName);
    if (r->nucleotides != NULL) free(r->nucleotides);
    if (r->qualities != NULL) free(r->qualities);
    free(r);
}

char *mergeContigChunks(char **chunks, int64_t startIdx, int64_t endIdxExclusive, int64_t overlap, Params *params,
        char *missingChunkSpacer) {

    // merge chunks
    stList *polishedReferenceStrings = stList_construct3(0, free); // The polished reference strings, one for each chunk

    for (int64_t chunkIdx = startIdx; chunkIdx < endIdxExclusive; chunkIdx++) {
        // Get chunk and polished
        char* currentChunk = chunks[chunkIdx];

        if (stList_length(polishedReferenceStrings) == 0) {
            // special case for first chunk
            // we must copy the first one because the original and merged strings are freed separately
            currentChunk = stString_copy(currentChunk);
        } else {
            char *previousChunk = stList_peek(polishedReferenceStrings);

            // Trim the currrent and previous polished reference strings to remove overlap
            int64_t prefixStringCropEnd, suffixStringCropStart;
            int64_t overlapMatchWeight = removeOverlap(previousChunk, currentChunk,
                                                       overlap, params->polishParams,
                                                       &prefixStringCropEnd, &suffixStringCropStart);

            // we have an overlap
            if (overlapMatchWeight > 0) {
                st_logInfo(
                        "  Removed overlap between neighbouring chunks. Approx overlap size: %i, overlap-match weight: %f, "
                        "left-trim: %i, right-trim: %i:\n", overlap,
                        (float) overlapMatchWeight / PAIR_ALIGNMENT_PROB_1,
                        strlen(previousChunk) - prefixStringCropEnd, suffixStringCropStart);

                // Crop the suffix of the previous chunk's polished reference string
                previousChunk[prefixStringCropEnd] = '\0';

                // Crop the the prefix of the current chunk's polished reference string
                currentChunk = stString_copy(&(currentChunk[suffixStringCropStart]));

            // no good alignment, likely missing chunks but have to be able to handle freaky situations also
            } else {
                if (strlen(currentChunk) == 0) {
                    // missing chunk
                    st_logInfo("  No overlap found. Filling empty chunk with Ns.\n");
                    currentChunk = stString_copy(missingChunkSpacer);
                } else if (overlap == 0) {
                    // poorly configured but could be done (freaky)
                    currentChunk = stString_copy(currentChunk);
                } else {
                    // couldn't find an overlap (freaky)
                    st_logInfo("  No overlap found. Filling Ns in stitch position.\n");
                    stList_append(polishedReferenceStrings, stString_copy("NNNNNNNNNN"));
                    currentChunk = stString_copy(currentChunk);
                }
            }
        }

        // Add the polished sequence to the list of polished reference sequence chunks
        stList_append(polishedReferenceStrings, currentChunk);
    }

    // finish
    char *merged = stString_join2("", polishedReferenceStrings);
    stList_destruct(polishedReferenceStrings);
    return merged;
}

char *mergeContigChunksThreaded(char **chunks, int64_t startIdx, int64_t endIdxExclusive, int64_t numThreads,
        int64_t overlap, Params *params, char *missingChunkSpacer) {

    // special unthreaded case
    if (numThreads == 1) return mergeContigChunks(chunks, startIdx, endIdxExclusive, overlap, params, missingChunkSpacer);

    // divide into chunks
    int64_t totalChunks = endIdxExclusive - startIdx;
    int64_t chunksPerThread = (int64_t) ceil(1.0 * totalChunks / numThreads);
    if (chunksPerThread * (numThreads - 1) >= endIdxExclusive) numThreads--;
    char **outputChunks = st_calloc(numThreads, sizeof(char*));

    // multithread loop
    #pragma omp parallel for schedule(static,1)
    for (int64_t thread = 0; thread < numThreads; thread++) {
        int64_t threadedStartIdx = startIdx + chunksPerThread * thread;
        int64_t threadedEndIdxExclusive = threadedStartIdx + chunksPerThread;
        if (endIdxExclusive < threadedEndIdxExclusive) threadedEndIdxExclusive = endIdxExclusive;

        outputChunks[thread] = mergeContigChunks(chunks, threadedStartIdx, threadedEndIdxExclusive, overlap,
                params, missingChunkSpacer);
    }

    // finish
    char *contig = mergeContigChunks(outputChunks, 0, numThreads, overlap, params, missingChunkSpacer);
    for (int64_t i = 0; i < numThreads; i++) {
        free(outputChunks[i]);
    }
    free(outputChunks);
    return contig;
}