/*
 * Copyright (C) 2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "margin.h"
#include "bedidx.h"

#define DEFAULT_ALIGNMENT_SCORE 10

BamChunkRead *bamChunkRead_construct() {
    return bamChunkRead_construct2(NULL, NULL, NULL, TRUE, NULL);
}

BamChunkRead *bamChunkRead_construct2(char *readName, char *nucleotides,
                                      uint8_t *qualities, bool forwardStrand, bool useRunLengthEncoding) {
    BamChunkRead *r = calloc(1, sizeof(BamChunkRead));
    r->readName = stString_copy(readName);
    r->forwardStrand = forwardStrand;
    assert(nucleotides != NULL);
    r->rleRead = useRunLengthEncoding ? rleString_construct(nucleotides) : rleString_construct_no_rle(nucleotides);
    if(qualities != NULL) {
        r->qualities = rleString_rleQualities(r->rleRead, qualities);
    }

    return r;
}
BamChunkRead *bamChunkRead_constructCopy(BamChunkRead *copy) {
    BamChunkRead *r = calloc(1, sizeof(BamChunkRead));
    r->readName = stString_copy(copy->readName);
    r->forwardStrand = copy->forwardStrand;
    r->rleRead = rleString_copy(copy->rleRead);
    if (copy->qualities != NULL) {
        r->qualities = st_calloc(r->rleRead->length, sizeof(uint8_t));
        for (int64_t i = 0; i < r->rleRead->length; i++) {
            r->qualities[i] = copy->qualities[i];
        }
    }

    return r;
}

void bamChunkRead_destruct(BamChunkRead *r) {
    if (r->readName != NULL) free(r->readName);
    if (r->rleRead != NULL) rleString_destruct(r->rleRead);
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
                        "    Removed overlap between neighbouring chunks at %"PRId64". "
                        "Approx overlap size: %i, overlap-match weight: %f, "
                        "left-trim: %i, right-trim: %i:\n", chunkIdx, overlap,
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
                    st_logInfo("    No overlap found for empty chunk at %"PRId64". Filling empty chunk with Ns.\n", chunkIdx);
                    currentChunk = stString_copy(missingChunkSpacer);
                } else if (overlap == 0) {
                    // poorly configured but could be done (freaky)
                    st_logInfo("    No overlap configured with non-empty (len %"PRId64") chunk at %"PRId64". \n",
                               strlen(currentChunk), chunkIdx);
                    currentChunk = stString_copy(currentChunk);
                } else {
                    // couldn't find an overlap (freaky)
                    st_logInfo("    No overlap found at %"PRId64". Filling Ns in stitch position.\n", chunkIdx);
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
                                int64_t overlap, Params *params, char *missingChunkSpacer, char *referenceSequenceName) {

    // special unthreaded case
    if (numThreads == 1) return mergeContigChunks(chunks, startIdx, endIdxExclusive, overlap, params, missingChunkSpacer);

    // divide into chunks
    int64_t totalChunks = endIdxExclusive - startIdx;
    int64_t chunksPerThread = (int64_t) ceil(1.0 * totalChunks / numThreads);
    while (startIdx + chunksPerThread * (numThreads - 1) >= endIdxExclusive) {numThreads--;}
    char **outputChunks = st_calloc(numThreads, sizeof(char*));

    // multithread loop
    st_logInfo("  Merging chunks for %s from (%"PRId64", %"PRId64"] with %"PRId64" chunks per thread on %"PRId64" threads \n",
               referenceSequenceName, startIdx, endIdxExclusive, chunksPerThread, numThreads);
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
