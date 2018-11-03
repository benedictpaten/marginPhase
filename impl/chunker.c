/*
 * Copyright (C) 2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "stRPHmm.h"
#include <stPolish.h>

//TODO put these in a parameters file
#define CHUNK_SIZE_DEFAULT 100000
#define CHUNK_MARGIN_DEFAULT 5000


int64_t saveContigChunks(stList *dest, BamChunker *parent, char *contig, int64_t contigStartPos, int64_t contigEndPos,
                         uint64_t chunkSize, uint64_t chunkMargin) {
    //TODO actually use chunkSize and chunkMargin

    int64_t chunkCount = 0;

    BamChunk *chunk = bamChunk_construct2(contig, contigStartPos, contigStartPos, contigEndPos, contigEndPos, parent);
    stList_append(dest, chunk);
    chunkCount++;

    return chunkCount;
}


BamChunker *bamChunker_construct2(char *bamFile, uint64_t chunkSize, uint64_t chunkMargin) {
    // the chunker we're building
    BamChunker *chunker = malloc(sizeof(BamChunker));
    chunker->bamFile = stString_copy(bamFile);
    chunker->chunkSize = chunkSize;
    chunker->chunkMargin = chunkMargin;
    chunker->chunks = stList_construct();
    chunker->chunkCount = 0;
    chunker->itorIdx = -1;

    // open bamfile
    samFile *in = hts_open(bamFile, "r");
    if (in == NULL) {
        st_errAbort("ERROR: Cannot open bam file %s\n", bamFile);
        return NULL;
    }
    bam_hdr_t *bamHdr = sam_hdr_read(in);
    bam1_t *aln = bam_init1();

    // list of chunk boundaries
    char *currentContig = NULL;
    int64_t contigStartPos = 0;
    int64_t contigEndPos = 0;

    // get all reads
    // there is probably a better way (bai?) to find min and max aligned positions (which we need for chunk divisions)
    while(sam_read1(in,bamHdr,aln) > 0) {

        // get aligned read length (and sanity check)
        int64_t unused = 0;
        int64_t readLength = getAlignedReadLength(aln, &unused);
        if (readLength <= 0) {
            continue;
        }

        // get start and stop position
        int64_t readStartPos = aln->core.pos + 1;           // Left most position of alignment
        int64_t readEndPos = readStartPos + readLength;

        // get contig
        char *contig = bamHdr->target_name[aln->core.tid];     // Contig name

        if (currentContig == NULL) {
            // first read
            currentContig = stString_copy(contig);
            contigStartPos = readStartPos;
            contigEndPos = readEndPos;
        } else if (stString_eq(currentContig, contig)) {
            // continue this contig's reads
            contigStartPos = readStartPos < contigStartPos ? readStartPos : contigStartPos;
            contigEndPos = readEndPos > contigEndPos ? readEndPos : contigEndPos;
        } else {
            // new contig
            int64_t savedChunkCount = saveContigChunks(chunker->chunks, chunker, currentContig,
                                                       contigStartPos, contigEndPos, chunkSize, chunkMargin);
            chunker->chunkCount += savedChunkCount;
            free(currentContig);
            currentContig = stString_copy(contig);
            contigStartPos = readStartPos;
            contigEndPos = readEndPos;
        }
    }
    // save last contig's chunks
    if (currentContig != NULL) {
        int64_t savedChunkCount = saveContigChunks(chunker->chunks, chunker, currentContig,
                                                   contigStartPos, contigEndPos, chunkSize, chunkMargin);
        chunker->chunkCount += savedChunkCount;
        free(currentContig);
    }

    // sanity check
    assert(stList_length(chunker->chunks) == chunker->chunkCount);

    // shut everything down
    bam_hdr_destroy(bamHdr);
    bam_destroy1(aln);
    sam_close(in);

    return chunker;
}

BamChunker *bamChunker_construct(char *bamFile) {
    return bamChunker_construct2(bamFile, CHUNK_SIZE_DEFAULT, CHUNK_MARGIN_DEFAULT);
}

void bamChunker_destruct(BamChunker *bamChunker) {
    //TODO
}

BamChunk *bamChunker_getNext(BamChunker *bamChunker) {
    // init if first invocation of getNext
    if (bamChunker->itorIdx < 0) {
        bamChunker->itorIdx = 0;
    }
    // handle end of list case
    if (bamChunker->itorIdx == bamChunker->chunkCount) {
        bamChunker->itorIdx = -1;
        return NULL;
    }

    // get chunk, increment, return
    BamChunk *chunk = stList_get(bamChunker->chunks, bamChunker->itorIdx);
    bamChunker->itorIdx += 1;
    return chunk;
}


BamChunk *bamChunk_construct() {
    return bamChunk_construct2(NULL, 0, 0, 0, 0, NULL);
}

BamChunk *bamChunk_construct2(char *refSeqName, int64_t chunkMarginStart, int64_t chunkStart, int64_t chunkEnd,
                              int64_t chunkMarginEnd, BamChunker *parent) {
    BamChunk *c = malloc(sizeof(BamChunk));
    c->refSeqName = stString_copy(refSeqName);
    c->chunkMarginStart = chunkMarginStart;
    c->chunkStart = chunkStart;
    c->chunkEnd = chunkEnd;
    c->chunkMarginEnd = chunkMarginEnd;
    c->parent = parent;
    return c;
}

void bamChunk_destruct(BamChunk *bamChunk) {
    //TODO
}