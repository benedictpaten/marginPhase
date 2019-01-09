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
    BamChunkRead *r = malloc(sizeof(BamChunkRead));
    r->readName = read->readName ==  NULL ? NULL : stString_copy(read->readName);
    r->nucleotides = stString_copy(rle->rleString);
    r->readLength = rle->length;
    r->qualities = NULL; //TODO maybe put some exciting logic here?
    r->forwardStrand = read->forwardStrand;
    r->parent = read->parent;

    return r;
}
void bamChunkRead_destruct(BamChunkRead *r) {
    if (r->readName != NULL) free(r->readName);
    if (r->nucleotides != NULL) free(r->nucleotides);
    if (r->qualities != NULL) free(r->qualities);
    free(r);
}

