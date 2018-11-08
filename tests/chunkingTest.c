/*
 * Copyright (C) 2009-2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "stPolish.h"

#define INPUT_BAM "../tests/chunkingTest.bam"

static void test_getChunksByChrom(CuTest *testCase) {
    BamChunker *chunker = bamChunker_construct(INPUT_BAM);
    CuAssertTrue(testCase, chunker->chunkCount == 2);
    int itorCount = 0;
    while (bamChunker_getNext(chunker) != NULL) {
        itorCount++;
    }
    CuAssertTrue(testCase, itorCount == 2);
    bamChunker_destruct(chunker);
}

static void test_getChunksBy100kb(CuTest *testCase) {
    BamChunker *chunker = bamChunker_construct2(INPUT_BAM, 100000, 0, FALSE);

    // contig_1 alignments start at 100 000 and go to 2 100 008 (21 @ 100k size)
    // contig_2 alignments start at 100 000 and go to   100 032 ( 1 @ 100k size)
    CuAssertTrue(testCase, chunker->chunkCount == 22);
    int itorCount = 0;
    while (bamChunker_getNext(chunker) != NULL) {
        itorCount++;
    }
    CuAssertTrue(testCase, itorCount == 22);
    bamChunker_destruct(chunker);
}

static void test_getChunksWithBoundary(CuTest *testCase) {
    BamChunker *chunker = bamChunker_construct2(INPUT_BAM, 8, 4, FALSE);

    // has 9 reads of 8 characters aligned to position 100 000 and every 4 bases after (last read aligned to 100 032)
    int contig2ChunkCount = 0;
    BamChunk *chunk = NULL;
    while((chunk = bamChunker_getNext(chunker)) != NULL) {
        // margin testing is on contig_2
        if (!stString_eq(chunk->refSeqName, "contig_2")) continue;

        // see how many reads there are
        stList *reads = stList_construct();
        stList *alignments = stList_construct();
        uint32_t readCount = convertToReadsAndAlignments(chunk, reads, alignments);
        CuAssertTrue(testCase, readCount == stList_length(reads));
        CuAssertTrue(testCase, readCount == stList_length(alignments));

        // all these should cover any read whose alignment overlaps (inclusive) with chunkMarginStart and (exclusive)
        // with chunkMarginEnd.  exclusive as in, if an alignment starts at 0 and has length 100, it should be included
        // if chunkMarginEnd is 99 but not 100
        if (contig2ChunkCount == 0) {
            //0: cbs:99996 cs:100000 ce:100008 cbe:100012 -> reads 0, 4, 8 ; !12+
            CuAssertTrue(testCase, readCount == 3);
        } else if (contig2ChunkCount == 1) {
            //1: cbs:100004 cs:100008 ce:1000016 cbe:100020 -> reads 0, 4, 8, 12, 16 ; !20+
            CuAssertTrue(testCase, readCount == 5);
        } else if (contig2ChunkCount == 2) {
            //2: cbs:100012 cs:100016 ce:1000024 cbe:100028 -> reads !4- ; 8, 12, 16, 20, 24 ; !28+
            CuAssertTrue(testCase, readCount == 5);
        } else if (contig2ChunkCount == 3) {
            //3: cbs:100020 cs:100024 ce:1000032 cbe:100036 -> reads !12- ; 16, 20, 24, 28, 32
            CuAssertTrue(testCase, readCount == 5);
        } else if (contig2ChunkCount == 4) {
            //4: cbs:100028 cs:100032 ce:1000040 cbe:100044 -> reads !20- ; 24, 28, 32
            CuAssertTrue(testCase, readCount == 3);
        }

        // increment
        contig2ChunkCount++;
        stList_destruct(reads);
        stList_destruct(alignments);
    }
    CuAssertTrue(testCase, contig2ChunkCount == 5);
    bamChunker_destruct(chunker);
}


static void test_getChunksWithoutBoundary(CuTest *testCase) {
    BamChunker *chunker = bamChunker_construct2(INPUT_BAM, 8, 0, FALSE);

    // has 9 reads of 8 characters aligned to position 100 000 and every 4 bases after (last read aligned to 100 032)
    int contig2ChunkCount = 0;
    BamChunk *chunk = NULL;
    while((chunk = bamChunker_getNext(chunker)) != NULL) {
        // margin testing is on contig_2
        if (!stString_eq(chunk->refSeqName, "contig_2")) continue;

        // see how many reads there are
        stList *reads = stList_construct();
        stList *alignments = stList_construct();
        uint32_t readCount = convertToReadsAndAlignments(chunk, reads, alignments);
        CuAssertTrue(testCase, readCount == stList_length(reads));
        CuAssertTrue(testCase, readCount == stList_length(alignments));

        // all these should cover any read whose alignment overlaps (inclusive) with chunkMarginStart and (exclusive)
        // with chunkMarginEnd.  exclusive as in, if an alignment starts at 0 and has length 100, it should be included
        // if chunkMarginEnd is 99 but not 100
        if (contig2ChunkCount == 0) {
            //0: cbs:100000 cs:100000 ce:100008 cbe:100008 -> reads 0, 4 ; !8+
            CuAssertTrue(testCase, readCount == 2);
        } else if (contig2ChunkCount == 1) {
            //1: cbs:100008 cs:100008 ce:1000016 cbe:1000016 -> reads !0 ; 4, 8, 12 ; !16+
            CuAssertTrue(testCase, readCount == 3);
        } else if (contig2ChunkCount == 2) {
            //2: cbs:100016 cs:100016 ce:1000024 cbe:1000024 -> reads !8- ; 12, 16, 20 ; !24+
            CuAssertTrue(testCase, readCount == 3);
        } else if (contig2ChunkCount == 3) {
            //3: cbs:100024 cs:100024 ce:1000032 cbe:1000032 -> reads !16- ; 20, 24, 28 ; !32
            CuAssertTrue(testCase, readCount == 3);
        } else if (contig2ChunkCount == 4) {
            //4: cbs:100032 cs:100032 ce:1000040 cbe:1000040 -> reads !24- ; 28, 32
            CuAssertTrue(testCase, readCount == 2);
        }

        // increment
        contig2ChunkCount++;
        stList_destruct(reads);
        stList_destruct(alignments);
    }
    CuAssertTrue(testCase, contig2ChunkCount == 5);
    bamChunker_destruct(chunker);
}

void assertClippingAlignmentMatchCount(CuTest *testCase, int64_t idx, stList *alignment) {
    switch (idx) {
        case 0: //8S8M
        case 1: //8M8S
        case 2: //4S8M4S
        case 4: //4S4M2D4M4S
        case 6: //4S1M1D6M1D1M4S
        case 7: //4H8S8M
        case 8: //8M8S4H
        case 9: //4H4S8M4S4H
            CuAssertTrue(testCase, stList_length(alignment) == 8);
            break;
        case 3: //4S2M4I2M4S
            CuAssertTrue(testCase, stList_length(alignment) == 4);
            break;
        case 5: //4S1M1I4M1I1M4S
            CuAssertTrue(testCase, stList_length(alignment) == 6);
            break;
        default:
            CuAssertTrue(testCase, FALSE);
    }
}

static void test_getReadsWithoutSoftClipping(CuTest *testCase) {
    BamChunker *chunker = bamChunker_construct2(INPUT_BAM, 100000, 0, FALSE);

    // have reads aligned to a region wholly between 200 000 and 299 999 for soft clip testing
    BamChunk *chunk = NULL;
    bool foundChunk = FALSE;
    while((chunk = bamChunker_getNext(chunker)) != NULL) {
        if (!stString_eq(chunk->refSeqName, "contig_1") || (chunk->chunkBoundaryStart != 200000 &&
                chunk->chunkBoundaryEnd != 4000000)) continue;
        CuAssertTrue(testCase, !foundChunk);
        foundChunk = TRUE;

        // analyze reads and alignments
        stList *reads = stList_construct3(0, free);
        stList *alignments = stList_construct3(0, (void*)stIntTuple_destruct);
        uint32_t readCount = convertToReadsAndAlignments(chunk, reads, alignments);
        CuAssertTrue(testCase, readCount == 10);
        for (int64_t i = 0; i < 10; i++) {
            // by design, each read is 8 characters without softclipped bases
            CuAssertTrue(testCase, strlen((char*) stList_get(reads, i)) == 8);
            // check the length of the cigar strings
            assertClippingAlignmentMatchCount(testCase, i, stList_get(alignments, i));
        }
        stList_destruct(reads);
        stList_destruct(alignments);
    }
    CuAssertTrue(testCase, foundChunk);
    bamChunker_destruct(chunker);
}

static void test_getReadsWithSoftClipping(CuTest *testCase) {
    BamChunker *chunker = bamChunker_construct2(INPUT_BAM, 100000, 0, TRUE);

    // have reads aligned to a region wholly between 200 000 and 299 999 for soft clip testing
    BamChunk *chunk = NULL;
    bool foundChunk = FALSE;
    while((chunk = bamChunker_getNext(chunker)) != NULL) {
        if (!stString_eq(chunk->refSeqName, "contig_1") || (chunk->chunkBoundaryStart != 200000 &&
            chunk->chunkBoundaryEnd != 4000000)) continue;
        CuAssertTrue(testCase, !foundChunk);
        foundChunk = TRUE;

        // analyze reads and alignments
        stList *reads = stList_construct();
        stList *alignments = stList_construct();
        uint32_t readCount = convertToReadsAndAlignments(chunk, reads, alignments);
        CuAssertTrue(testCase, readCount == 10);
        for (int64_t i = 0; i < 10; i++) {
            // by design, each read is 16 characters with softclipped bases
            CuAssertTrue(testCase, strlen((char*) stList_get(reads, i)) == 16);
            // check the length of the cigar strings
            assertClippingAlignmentMatchCount(testCase, i, stList_get(alignments, i));
        }
        stList_destruct(reads);
        stList_destruct(alignments);
    }
    CuAssertTrue(testCase, foundChunk);
    bamChunker_destruct(chunker);
}


CuSuite* chunkingTestSuite(void) {
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, test_getChunksByChrom);
    SUITE_ADD_TEST(suite, test_getChunksBy100kb);
    SUITE_ADD_TEST(suite, test_getChunksWithBoundary);
    SUITE_ADD_TEST(suite, test_getChunksWithoutBoundary);
    SUITE_ADD_TEST(suite, test_getReadsWithSoftClipping);
    SUITE_ADD_TEST(suite, test_getReadsWithoutSoftClipping);

    return suite;
}
