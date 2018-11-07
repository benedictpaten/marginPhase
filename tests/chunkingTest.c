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
    }
    CuAssertTrue(testCase, contig2ChunkCount == 5);
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
    }
    CuAssertTrue(testCase, contig2ChunkCount == 5);
}


CuSuite* chunkingTestSuite(void) {
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, test_getChunksByChrom);
    SUITE_ADD_TEST(suite, test_getChunksBy100kb);
    SUITE_ADD_TEST(suite, test_getChunksWithBoundary);
    SUITE_ADD_TEST(suite, test_getChunksWithoutBoundary);

    return suite;
}
