/*
 * Copyright (C) 2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "stRPHmm.h"
#include <stPolish.h>

int64_t saveContigChunks(stList *dest, BamChunker *parent, char *contig, int64_t contigStartPos, int64_t contigEndPos,
                         uint64_t chunkSize, uint64_t chunkMargin) {

    // whole contig case
    if (chunkSize == 0) {
        BamChunk *chunk = bamChunk_construct2(contig, contigStartPos, contigStartPos, contigEndPos, contigEndPos,
                                              parent);
        stList_append(dest, chunk);
        return 1;
    }

    // specific chunk size
    int64_t chunkCount = 0;
    for (int64_t i = contigStartPos; i < contigEndPos; i += chunkSize) {
        BamChunk *chunk = bamChunk_construct2(contig, (chunkMargin > i ? 0 : i - chunkMargin), i, i + chunkSize,
                                              i + chunkSize + chunkMargin, parent);
        stList_append(dest, chunk);
        chunkCount++;
    }
    return chunkCount;
}


BamChunker *bamChunker_construct2(char *bamFile, uint64_t chunkSize, uint64_t chunkBoundary, bool includeSoftClip) {
    // the chunker we're building
    BamChunker *chunker = malloc(sizeof(BamChunker));
    chunker->bamFile = stString_copy(bamFile);
    chunker->chunkSize = chunkSize;
    chunker->chunkBoundary = chunkBoundary;
    chunker->includeSoftClip = includeSoftClip;
    chunker->chunks = stList_construct3(0,(void*)bamChunk_destruct);
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

        int64_t start_softclip = 0;
        int64_t end_softclip = 0;
        int64_t readLength =  getAlignedReadLength3(aln, &start_softclip, &end_softclip, FALSE);
        if (readLength <= 0) {
            continue;
        }

        // get start and stop position
        int64_t readStartPos = aln->core.pos + 1;           // Left most position of alignment
        int64_t readEndPos = readStartPos + readLength;
        //TODO extend by softclipping?

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
                                                       contigStartPos, contigEndPos, chunkSize, chunkBoundary);
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
                                                   contigStartPos, contigEndPos, chunkSize, chunkBoundary);
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
    return bamChunker_construct2(bamFile, UINT32_MAX, 0, FALSE);
}

void bamChunker_destruct(BamChunker *bamChunker) {
    free(bamChunker->bamFile);
    stList_destruct(bamChunker->chunks);
    free(bamChunker);
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

BamChunk *bamChunk_construct2(char *refSeqName, int64_t chunkBoundaryStart, int64_t chunkStart, int64_t chunkEnd,
                              int64_t chunkBoundaryEnd, BamChunker *parent) {
    BamChunk *c = malloc(sizeof(BamChunk));
    c->refSeqName = stString_copy(refSeqName);
    c->chunkBoundaryStart = chunkBoundaryStart;
    c->chunkStart = chunkStart;
    c->chunkEnd = chunkEnd;
    c->chunkBoundaryEnd = chunkBoundaryEnd;
    c->parent = parent;
    return c;
}

void bamChunk_destruct(BamChunk *bamChunk) {
    free(bamChunk->refSeqName);
    free(bamChunk);
}


uint32_t convertToReadsAndAlignments(BamChunk *bamChunk, stList *reads, stList *alignments) {
    // sanity check
    assert(stList_length(reads) == 0);
    assert(stList_length(alignments) == 0);

    // prep
    int64_t chunkStart = bamChunk->chunkBoundaryStart;
    int64_t chunkEnd = bamChunk->chunkBoundaryEnd;
    bool includeSoftClip = bamChunk->parent->includeSoftClip;
    char *bamFile = bamChunk->parent->bamFile;
    char *contig = bamChunk->refSeqName;
    uint32_t savedAlignments = 0;

    // get header, init align object
    samFile *in = hts_open(bamFile, "r");
    if (in == NULL) {
        st_errAbort("ERROR: Cannot open bam file %s\n", bamFile);
        return 0; //errAbort dies, this is not really a return value
    }
    bam_hdr_t *bamHdr = sam_hdr_read(in);
    bam1_t *aln = bam_init1();

    // read in alignments
    while(sam_read1(in,bamHdr,aln) > 0) {
        // basic filtering (no read length, no cigar)
        if (aln->core.l_qseq <= 0) continue;
        if (aln->core.n_cigar == 0) continue;

        //data
        char *chr = bamHdr->target_name[aln->core.tid];
        int64_t start_softclip = 0;
        int64_t end_softclip = 0;
        int64_t alnReadLength = getAlignedReadLength3(aln, &start_softclip, &end_softclip, FALSE);
        if (alnReadLength <= 0) continue;
        int64_t alnStartPos = aln->core.pos + 1;
        int64_t alnEndPos = alnStartPos + alnReadLength;

        // does this belong in our chunk?
        if (!stString_eq(contig, chr)) continue;
        if (alnStartPos >= chunkEnd) continue;
        if (alnEndPos <= chunkStart) continue;

        // get cigar and rep
        uint32_t *cigar = bam_get_cigar(aln);
        stList *cigRepr = stList_construct();

        // Variables to keep track of position in sequence / cigar operations
        int64_t cig_idx = 0;
        int64_t currPosInOp = 0;
        int64_t cigarOp = -1;
        int64_t cigarNum = -1;
        int64_t cigarIdxInSeq = 0;
        int64_t cigarIdxInRef = alnStartPos;

        // positional modifications
        int64_t refCigarModification = -1 * chunkStart;

        // we need to calculate:
        //  a. where in the (potentially softclipped read) to start storing characters
        //  b. what the alignments are wrt those characters
        // so we track the first aligned character in the read (for a.) and what alignment modification to make (for b.)
        int64_t seqCigarModification;
        int64_t firstNonSoftclipAlignedReadIdxInChunk;

        // the handling changes based on softclip inclusion and where the chunk boundaries are
        if (includeSoftClip) {
            if (alnStartPos < chunkStart) {
                // alignment spans chunkStart (this will not be affected by softclipping)
                firstNonSoftclipAlignedReadIdxInChunk = -1; //need to find position of first alignment
                seqCigarModification = 0;
            } else if (alnStartPos - start_softclip <= chunkStart) {
                // softclipped bases span chunkStart
                firstNonSoftclipAlignedReadIdxInChunk = 0;
                int64_t includedSoftclippedBases = alnStartPos - chunkStart;
                seqCigarModification = includedSoftclippedBases;
                assert(includedSoftclippedBases >= 0);
                assert(start_softclip - includedSoftclippedBases >= 0);
            } else {
                // softclipped bases are after chunkStart
                firstNonSoftclipAlignedReadIdxInChunk = 0;
                seqCigarModification = start_softclip;
            }
        } else {
            if (alnStartPos < chunkStart) {
                // alignment spans chunkStart
                firstNonSoftclipAlignedReadIdxInChunk = -1;
                seqCigarModification = 0;
            } else {
                // alignment starts after chunkStart
                firstNonSoftclipAlignedReadIdxInChunk = 0;
                seqCigarModification = 0;
            }
        }

        // track number of characters in aligned portion (will inform softclipping at end of read)
        int64_t alignedReadLength = 0;

        // iterate over cigar operations
        for (uint32_t i = 0; i < alnReadLength; i++) {

            // do we need the next cigar operation?
            if (currPosInOp == 0) {
                cigarOp = cigar[cig_idx] & BAM_CIGAR_MASK;
                cigarNum = cigar[cig_idx] >> BAM_CIGAR_SHIFT;
            }

            // handle current character
            if (cigarOp == BAM_CMATCH || cigarOp == BAM_CEQUAL || cigarOp == BAM_CDIFF) {
                if (cigarIdxInRef >= chunkStart && cigarIdxInRef < chunkEnd) {
                    stList_append(cigRepr, stIntTuple_construct2(cigarIdxInSeq + seqCigarModification,
                                                                 cigarIdxInRef + refCigarModification));
                    alignedReadLength++;
                }
                cigarIdxInSeq++;
                cigarIdxInRef++;
            } else if (cigarOp == BAM_CDEL || cigarOp == BAM_CREF_SKIP) {
                //delete
                cigarIdxInRef++;
            } else if (cigarOp == BAM_CINS) {
                //insert
                cigarIdxInSeq++;
                if (cigarIdxInRef >= chunkStart && cigarIdxInRef < chunkEnd) {
                    alignedReadLength++;
                }
                i--;
            } else if (cigarOp == BAM_CSOFT_CLIP || cigarOp == BAM_CHARD_CLIP || cigarOp == BAM_CPAD) {
                // nothing to do here. skip to next cigar operation
                currPosInOp = cigarNum - 1;
                i--;
            } else {
                st_logCritical("Unidentifiable cigar operation!\n");
            }

            // document read index in the chunk (for reads that span chunk boundary, used in read construction)
            if (firstNonSoftclipAlignedReadIdxInChunk < 0 && cigarIdxInRef >= chunkStart) {
                firstNonSoftclipAlignedReadIdxInChunk = cigarIdxInSeq;
                seqCigarModification = -1 * (firstNonSoftclipAlignedReadIdxInChunk + seqCigarModification);

            }

            // have we finished this last cigar
            currPosInOp++;
            if (currPosInOp == cigarNum) {
                cig_idx++;
                currPosInOp = 0;
            }
        }
        // sanity checks
        //TODO these may fail because of the existance of non-match end cigar operations
//        assert(cigarIdxInRef == alnEndPos);  //does not include soft clip
//        assert(cigarIdxInSeq == readEndIdx - (includeSoftClip ? end_softclip : 0));


        // get sequence positions
        int64_t seqLen = alignedReadLength;

        // modify start indices
        int64_t readCurrIdx = firstNonSoftclipAlignedReadIdxInChunk;
        if (firstNonSoftclipAlignedReadIdxInChunk != 0) {
            // the aligned portion spans chunkStart, so no softclipped bases are included
            readCurrIdx += start_softclip;
        } else if (!includeSoftClip) {
            // configured to not handle softclipped bases
            readCurrIdx += start_softclip;
        } else if (alnStartPos - start_softclip <= chunkStart) {
            // configured to handle softclipped bases; softclipped bases span chunkStart
            int64_t includedSoftclippedBases = alnStartPos - chunkStart;
            seqLen += includedSoftclippedBases;
            readCurrIdx += (start_softclip - includedSoftclippedBases);
        } else {
            // configured to handle softclipped bases; softclipped bases all occur after chunkStart
            seqLen += start_softclip;
            readCurrIdx = 0;
        }

        // modify end indices
        int64_t readEndIdx = readCurrIdx + seqLen;

        //TODO fix end-softclipped reads spanning end of chunk

        // get sequence
        char *seq = malloc((seqLen + 1) * sizeof(char));
        uint8_t *seqBits = bam_get_seq(aln);
        int64_t seqIdx = 0;
        while (readCurrIdx < readEndIdx) {
            seq[seqIdx] = seq_nt16_str[bam_seqi(seqBits, readCurrIdx)];
            readCurrIdx++;
            seqIdx++;
        }
        seq[seqLen] = '\0';


        // save
        stList_append(reads, seq);
        stList_append(alignments, cigRepr);
        savedAlignments++;
    }

    // close it all down
    bam_hdr_destroy(bamHdr);
    bam_destroy1(aln);
    sam_close(in);

    return savedAlignments;
}