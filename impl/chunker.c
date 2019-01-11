/*
 * Copyright (C) 2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

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

