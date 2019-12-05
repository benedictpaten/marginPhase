/*
 * Copyright (C) 2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <omp.h>
#include "margin.h"
#include "bedidx.h"

#define DEFAULT_ALIGNMENT_SCORE 10

//todo removed in THE_MERGE (lives in htsIntegration)
//int64_t saveContigChunks(stList *dest, BamChunker *parent, char *contig, int64_t contigStartPos, int64_t contigEndPos,
//                         uint64_t chunkSize, uint64_t chunkMargin) {
//
//    // whole contig case
//    if (chunkSize == 0) {
//        BamChunk *chunk = bamChunk_construct2(contig, contigStartPos, contigStartPos, contigEndPos, contigEndPos,
//                                              parent);
//        stList_append(dest, chunk);
//        return 1;
//    }
//
//    // specific chunk size
//    int64_t chunkCount = 0;
//    for (int64_t i = contigStartPos; i < contigEndPos; i += chunkSize) {
//        int64_t chunkEndPos = i + chunkSize;
//        chunkEndPos = (chunkEndPos > contigEndPos ? contigEndPos : chunkEndPos);
//        int64_t chunkMarginEndPos = chunkEndPos + chunkMargin;
//        chunkMarginEndPos = (chunkMarginEndPos > contigEndPos ? contigEndPos : chunkMarginEndPos);
//
//        BamChunk *chunk = bamChunk_construct2(contig, (chunkMargin > i ? 0 : i - chunkMargin), i, chunkEndPos,
//                                              chunkMarginEndPos, parent);
//        stList_append(dest, chunk);
//        chunkCount++;
//    }
//    return chunkCount;
//}
//
//
//BamChunker *bamChunker_construct(char *bamFile, PolishParams *params) {
//    return bamChunker_construct2(bamFile, NULL, params);
//}
//
//
//BamChunker *bamChunker_construct2(char *bamFile, char *region, PolishParams *params) {
//
//    // are we doing region filtering?
//    bool filterByRegion = false;
//    char regionContig[128] = "";
//    int regionStart = 0;
//    int regionEnd = 0;
//    if (region != NULL) {
//        int scanRet = sscanf(region, "%[^:]:%d-%d", regionContig, &regionStart, &regionEnd);
//        if (scanRet != 3 || strlen(regionContig) == 0) {
//            st_errAbort("Region in unexpected format (expected %%s:%%d-%%d): %s", region);
//        } else if (regionStart < 0 || regionEnd <= 0 || regionEnd <= regionStart) {
//            st_errAbort("Start and end locations in region must be positive, start must be less than end: %s", region);
//        }
//        filterByRegion = true;
//    }
//
//    // standard parameters
//    uint64_t chunkSize = params->chunkSize;
//    uint64_t chunkBoundary = params->chunkBoundary;
//    bool includeSoftClip = params->includeSoftClipping;
//
//    // the chunker we're building
//    BamChunker *chunker = malloc(sizeof(BamChunker));
//    chunker->bamFile = stString_copy(bamFile);
//    chunker->chunkSize = chunkSize;
//    chunker->chunkBoundary = chunkBoundary;
//    chunker->includeSoftClip = includeSoftClip;
//    chunker->params = params;
//    chunker->chunks = stList_construct3(0,(void*)bamChunk_destruct);
//    chunker->chunkCount = 0;
//    chunker->itorIdx = -1;
//
//    // open bamfile
//    samFile *in = hts_open(bamFile, "r");
//    if (in == NULL)
//        st_errAbort("ERROR: Cannot open bam file %s\n", bamFile);
//    hts_idx_t *idx = sam_index_load(in, bamFile); // load index (just to verify early that it exists)
//    if (idx == NULL)
//        st_errAbort("ERROR: Missing index for bam file %s\n", bamFile);
//    hts_idx_destroy(idx);
//    bam_hdr_t *bamHdr = sam_hdr_read(in);
//    bam1_t *aln = bam_init1();
//
//    // list of chunk boundaries
//    char *currentContig = NULL;
//    int64_t contigStartPos = 0;
//    int64_t contigEndPos = 0;
//
//    // get all reads
//    // there is probably a better way (bai?) to find min and max aligned positions (which we need for chunk divisions)
//    while(sam_read1(in,bamHdr,aln) > 0) {
//
//        // basic filtering (no read length, no cigar)
//        if (aln->core.l_qseq <= 0) continue;
//        if (aln->core.n_cigar == 0) continue;
//
//        //data
//        char *chr = bamHdr->target_name[aln->core.tid];
//        int64_t start_softclip = 0;
//        int64_t end_softclip = 0;
//        int64_t alnReadLength = getAlignedReadLength3(aln, &start_softclip, &end_softclip, FALSE);
//        int64_t alnStartPos = aln->core.pos + 1;
//        int64_t alnEndPos = alnStartPos + alnReadLength;
//
//        // does this belong in our chunk?
//        if (alnReadLength <= 0) continue;
//        if (filterByRegion && (!stString_eq(regionContig, chr) || alnStartPos >= regionEnd || alnEndPos <= regionStart))
//            continue;
//
//        // get start and stop position
//        int64_t readStartPos = aln->core.pos + 1;           // Left most position of alignment
//        int64_t readEndPos = readStartPos + alnReadLength;
//
//        // get contig
//        char *contig = bamHdr->target_name[aln->core.tid];     // Contig name
//
//        if (currentContig == NULL) {
//            // first read
//            currentContig = stString_copy(contig);
//            contigStartPos = readStartPos;
//            contigEndPos = readEndPos;
//        } else if (stString_eq(currentContig, contig)) {
//            // continue this contig's reads
//            contigStartPos = readStartPos < contigStartPos ? readStartPos : contigStartPos;
//            contigEndPos = readEndPos > contigEndPos ? readEndPos : contigEndPos;
//        } else {
//            // new contig (this should never happen if we're filtering by region)
//            assert(!filterByRegion);
//            int64_t savedChunkCount = saveContigChunks(chunker->chunks, chunker, currentContig,
//                                                       contigStartPos, contigEndPos, chunkSize, chunkBoundary);
//            chunker->chunkCount += savedChunkCount;
//            free(currentContig);
//            currentContig = stString_copy(contig);
//            contigStartPos = readStartPos;
//            contigEndPos = readEndPos;
//        }
//    }
//    // save last contig's chunks
//    if (currentContig != NULL) {
//        if (filterByRegion) {
//            contigStartPos = (contigStartPos < regionStart ? regionStart : contigStartPos);
//            contigEndPos = (contigEndPos > regionEnd ? regionEnd : contigEndPos);
//        }
//        int64_t savedChunkCount = saveContigChunks(chunker->chunks, chunker, currentContig,
//                                                   contigStartPos, contigEndPos, chunkSize, chunkBoundary);
//        chunker->chunkCount += savedChunkCount;
//        free(currentContig);
//    }
//
//    // sanity check
//    assert(stList_length(chunker->chunks) == chunker->chunkCount);
//
//    // shut everything down
//    bam_hdr_destroy(bamHdr);
//    bam_destroy1(aln);
//    sam_close(in);
//
//    return chunker;
//}
//
//void bamChunker_destruct(BamChunker *bamChunker) {
//    free(bamChunker->bamFile);
//    stList_destruct(bamChunker->chunks);
//    free(bamChunker);
//}
//
//BamChunk *bamChunker_getNext(BamChunker *bamChunker) {
//    // init if first invocation of getNext
//    if (bamChunker->itorIdx < 0) {
//        bamChunker->itorIdx = 0;
//    }
//    // handle end of list case
//    if (bamChunker->itorIdx == bamChunker->chunkCount) {
//        bamChunker->itorIdx = -1;
//        return NULL;
//    }
//
//    // get chunk, increment, return
//    BamChunk *chunk = stList_get(bamChunker->chunks, bamChunker->itorIdx);
//    bamChunker->itorIdx += 1;
//    return chunk;
//}
//
//BamChunk *bamChunk_construct() {
//    return bamChunk_construct2(NULL, 0, 0, 0, 0, NULL);
//}
//
//BamChunk *bamChunk_construct2(char *refSeqName, int64_t chunkBoundaryStart, int64_t chunkStart, int64_t chunkEnd,
//                              int64_t chunkBoundaryEnd, BamChunker *parent) {
//    BamChunk *c = malloc(sizeof(BamChunk));
//    c->refSeqName = stString_copy(refSeqName);
//    c->chunkBoundaryStart = chunkBoundaryStart;
//    c->chunkStart = chunkStart;
//    c->chunkEnd = chunkEnd;
//    c->chunkBoundaryEnd = chunkBoundaryEnd;
//    c->parent = parent;
//    return c;
//}
//
//void bamChunk_destruct(BamChunk *bamChunk) {
//    free(bamChunk->refSeqName);
//    free(bamChunk);
//}

BamChunkRead *bamChunkRead_construct() {
    return bamChunkRead_construct2(NULL, NULL, NULL, TRUE, NULL);
}

//todo
//BamChunkRead *bamChunkRead_construct2(char *readName, char *nucleotides, uint8_t *qualities, bool forwardStrand,
//                                      BamChunk *parent) {
//    BamChunkRead *r = malloc(sizeof(BamChunkRead));

BamChunkRead *bamChunkRead_construct2(char *readName, char *nucleotides,
                                      uint8_t *qualities, bool forwardStrand, bool useRunLengthEncoding) {
    BamChunkRead *r = calloc(1, sizeof(BamChunkRead));
    r->readName = readName;
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

// this all was moved to htsIntegration
//// This structure holds the bed information
//typedef struct samview_settings {
//    void* bed;
//} samview_settings_t;
//
//
//uint32_t convertToReadsAndAlignments(BamChunk *bamChunk, RleString *reference, stList *reads, stList *alignments) {
//    // sanity check
//    assert(stList_length(reads) == 0);
//    assert(stList_length(alignments) == 0);
//
//    uint64_t *ref_nonRleToRleCoordinateMap = reference == NULL ? NULL : rleString_getNonRleToRleCoordinateMap(reference);
//
//    // prep
//    int64_t chunkStart = bamChunk->chunkBoundaryStart;
//    int64_t chunkEnd = bamChunk->chunkBoundaryEnd;
//    bool includeSoftClip = bamChunk->parent->params->includeSoftClipping;
//    char *bamFile = bamChunk->parent->bamFile;
//    char *contig = bamChunk->refSeqName;
//    uint32_t savedAlignments = 0;
//
//    // prep for index (not entirely sure what all this does.  see samtools/sam_view.c
//    int filter_state = ALL, filter_op = 0;
//    int result;
//    samview_settings_t settings = { .bed = NULL };
//    char* region[1] = {};
//    region[0] = stString_print("%s:%d-%d", bamChunk->refSeqName, bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
//    settings.bed = bed_hash_regions(settings.bed, region, 0, 1, &filter_op); //insert(1) or filter out(0) the regions from the command line in the same hash table as the bed file
//    if (!filter_op) filter_state = FILTERED;
//    int regcount = 0;
//    hts_reglist_t *reglist = bed_reglist(settings.bed, filter_state, &regcount);
//    if(!reglist) {
//        st_errAbort("ERROR: Could not create list of regions for read conversion");
//    }
//
//    // file initialization
//    samFile *in = NULL;
//    hts_idx_t *idx = NULL;
//    hts_itr_multi_t *iter = NULL;
//    // bam file
//    if ((in = hts_open(bamFile, "r")) == 0) {
//        st_errAbort("ERROR: Cannot open bam file %s\n", bamFile);
//    }
//    // bam index
//    if ((idx = sam_index_load(in, bamFile)) == 0) {
//        st_errAbort("ERROR: Cannot open index for bam file %s\n", bamFile);
//    }
//    // header  //todo samFile *in = hts_open(bamFile, "r");
//    bam_hdr_t *bamHdr = sam_hdr_read(in);
//    // read object
//    bam1_t *aln = bam_init1();
//    // iterator for region
//    if ((iter = sam_itr_regions(idx, bamHdr, reglist, regcount)) == 0) {
//        st_errAbort("ERROR: Cannot open iterator for region %s for bam file %s\n", region[0], bamFile);
//    }
//
//    // fetch alignments //todo while(sam_read1(in,bamHdr,aln) > 0) {
//    while ((result = sam_itr_multi_next(in, iter, aln)) >= 0) {
//        // basic filtering (no read length, no cigar)
//        if (aln->core.l_qseq <= 0) continue;
//        if (aln->core.n_cigar == 0) continue;
//
//        //data
//        char *chr = bamHdr->target_name[aln->core.tid];
//        int64_t start_softclip = 0;
//        int64_t end_softclip = 0;
//        int64_t alnReadLength = getAlignedReadLength3(aln, &start_softclip, &end_softclip, FALSE);
//        if (alnReadLength <= 0) continue;
//        int64_t alnStartPos = aln->core.pos + 1;
//        int64_t alnEndPos = alnStartPos + alnReadLength;
//
//        // does this belong in our chunk?
//        if (!stString_eq(contig, chr)) continue;
//        if (alnStartPos >= chunkEnd) continue;
//        if (alnEndPos <= chunkStart) continue;
//
//        // get cigar and rep
//        uint32_t *cigar = bam_get_cigar(aln);
//        stList *cigRepr = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
//        // get sequence - all data we need is encoded in readStartIdxInChunk (start), readEnd idx, and seqLen
//        char *seq = st_calloc(seqLen + 1, sizeof(char));
//        uint8_t *seqBits = bam_get_seq(aln);
//        int64_t idxInOutputSeq = 0;
//        int64_t idxInBamRead = readStartIdxInChunk;
//        while (idxInBamRead < readEndIdxInChunk) {
//            seq[idxInOutputSeq] = seq_nt16_str[bam_seqi(seqBits, idxInBamRead)];
//            idxInBamRead++;
//            idxInOutputSeq++;
//        }
//        seq[seqLen] = '\0';
//
//        // get sequence qualities (if exists)
//        char *readName = stString_copy(bam_get_qname(aln));
//        uint8_t *qualBits = bam_get_qual(aln);
//        uint8_t *qual = NULL;
//        if (qualBits[0] != 0xff) { //inital score of 255 means qual scores are unavailable
//            idxInOutputSeq = 0;
//            idxInBamRead = readStartIdxInChunk;
//            qual = st_calloc(seqLen, sizeof(uint8_t));
//            while (idxInBamRead < readEndIdxInChunk) {
//                qual[idxInOutputSeq] = qualBits[idxInBamRead];
//                idxInBamRead++;
//                idxInOutputSeq++;
//
//            }
//            //assert(idxInBamRead == strlen(seq));
//        };
//
//        // save to read
//        bool forwardStrand = !bam_is_rev(aln);
//        BamChunkRead *chunkRead = bamChunkRead_construct2(readName, seq, qual, forwardStrand,
//        		bamChunk->parent->params->useRunLengthEncoding);
//        stList_append(reads, chunkRead);
//
//        if(bamChunk->parent->params->useRunLengthEncoding) {
//			// rle the alignment and save it
//			uint64_t *read_nonRleToRleCoordinateMap = rleString_getNonRleToRleCoordinateMap(chunkRead->rleRead);
//			stList_append(alignments, runLengthEncodeAlignment(cigRepr, ref_nonRleToRleCoordinateMap, read_nonRleToRleCoordinateMap));
//			stList_destruct(cigRepr);
//			free(read_nonRleToRleCoordinateMap);
//        }
//        else {
//        	stList_append(alignments, cigRepr);
//        }
//
//        savedAlignments++;
//    }
//    // the status from "get reads from iterator"
//    if (result < -1) {
//        st_errAbort("ERROR: Retrieval of region %d failed due to truncated file or corrupt BAM index file\n", iter->curr_tid);
//    }
//
//    // close it all down
//    hts_itr_multi_destroy(iter);
//    hts_idx_destroy(idx);
//    free(region[0]);
//    bed_destroy(settings.bed);
//    bam_hdr_destroy(bamHdr);
//    bam_destroy1(aln);
//    sam_close(in);
//    if(reference != NULL) {
//    	free(ref_nonRleToRleCoordinateMap);
//    }
//
//    return savedAlignments;
//}





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
    while (chunksPerThread * (numThreads - 1) >= endIdxExclusive) {numThreads--;}
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