//
// Created by tpesout on 1/8/19.
//

#include "htsIntegration.h"
#include "margin.h"

/*
 * getAlignedReadLength computes the length of the read sequence which is aligned to the reference.  Hard-clipped bases
 * are never included in this calculation.  Soft-clipped bases are similarly not included, but will be returned via
 * the two parameters if the 2nd or 3rd version is invoked.  The boundaryAtMatch parameter handles cases where an
 * insertion or deletion is the first or last operation after or before clipping.  If set, these will be treated as
 * soft-clipping; otherwise they will included in the final return value.
 */
int64_t getAlignedReadLength(bam1_t *aln) {
    int64_t start_softclip = 0;
    int64_t end_softclip = 0;
    return getAlignedReadLength3(aln, &start_softclip, &end_softclip, TRUE);
}
int64_t getAlignedReadLength2(bam1_t *aln, int64_t *start_softclip, int64_t *end_softclip) {
    return getAlignedReadLength3(aln, start_softclip, end_softclip, TRUE);
}
int64_t getAlignedReadLength3(bam1_t *aln, int64_t *start_softclip, int64_t *end_softclip, bool boundaryAtMatch) {
    // start read needs to be init'd to 0 (mostly this is to avoid misuse)
    if (*start_softclip != 0) st_errAbort("getAlignedReadLength2 invoked with improper start_softclip parameter");
    if (*end_softclip != 0) st_errAbort("getAlignedReadLength2 invoked with improper end_softclip parameter");

    // get relevant cigar info
    int64_t len = aln->core.l_qseq;
    uint32_t *cigar = bam_get_cigar(aln);

    // data for tracking
    int64_t start_ref = 0;
    int64_t cig_idx = 0;

    // Find the correct starting locations on the read and reference sequence,
    // to deal with things like inserts / deletions / soft clipping
    while (cig_idx < aln->core.n_cigar) {
        int cigarOp = cigar[cig_idx] & BAM_CIGAR_MASK;
        int cigarNum = cigar[cig_idx] >> BAM_CIGAR_SHIFT;

        if (cigarOp == BAM_CMATCH || cigarOp == BAM_CEQUAL || cigarOp == BAM_CDIFF) {
            break;
        } else if (cigarOp == BAM_CDEL || cigarOp == BAM_CREF_SKIP) {
            if (!boundaryAtMatch) break;
            cig_idx++;
        } else if (cigarOp == BAM_CINS) {
            if (!boundaryAtMatch) break;
            *start_softclip += cigarNum;
            cig_idx++;
        } else if (cigarOp == BAM_CSOFT_CLIP) {
            *start_softclip += cigarNum;
            cig_idx++;
        } else if (cigarOp == BAM_CHARD_CLIP || cigarOp == BAM_CPAD) {
            cig_idx++;
        } else {
            st_errAbort("Unidentifiable cigar operation\n");
        }
    }

    // Check for soft clipping at the end
    cig_idx = aln->core.n_cigar - 1;
    while (cig_idx > 0) {
        int cigarOp = cigar[cig_idx] & BAM_CIGAR_MASK;
        int cigarNum = cigar[cig_idx] >> BAM_CIGAR_SHIFT;

        if (cigarOp == BAM_CMATCH || cigarOp == BAM_CEQUAL || cigarOp == BAM_CDIFF) {
            break;
        } else if (cigarOp == BAM_CDEL || cigarOp == BAM_CREF_SKIP) {
            if (!boundaryAtMatch) break;
            cig_idx--;
        } else if (cigarOp == BAM_CINS) {
            if (!boundaryAtMatch) break;
            *end_softclip += cigarNum;
            cig_idx--;
        } else if (cigarOp == BAM_CSOFT_CLIP) {
            *end_softclip += cigarNum;
            cig_idx--;
        } else if (cigarOp == BAM_CHARD_CLIP || cigarOp == BAM_CPAD) {
            cig_idx--;
        } else {
            st_errAbort("Unidentifiable cigar operation\n");
        }
    }

    // Count number of insertions & deletions in sequence
    int64_t numInsertions = 0;
    int64_t numDeletions = 0;
    countIndels(cigar, aln->core.n_cigar, &numInsertions, &numDeletions);
    int64_t trueLength = len - *start_softclip - *end_softclip + numDeletions - numInsertions;

    return trueLength;
}


/*
 * Counts the number of insertions and deletions in a read, given its cigar string.
 */
void countIndels(uint32_t *cigar, uint32_t ncigar, int64_t *numInsertions, int64_t *numDeletions) {
    for (uint32_t i = 0; i < ncigar; i++) {
        int cigarOp = cigar[i] & BAM_CIGAR_MASK;
        int cigarNum = cigar[i] >> BAM_CIGAR_SHIFT;
        if (cigarOp == BAM_CINS) *numInsertions += cigarNum;
        if (cigarOp == BAM_CDEL) *numDeletions += cigarNum;
    }
}

/*
 * Utility function for BamChunk constructor
 */
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
        int64_t chunkEndPos = i + chunkSize;
        chunkEndPos = (chunkEndPos > contigEndPos ? contigEndPos : chunkEndPos);
        int64_t chunkMarginEndPos = chunkEndPos + chunkMargin;
        chunkMarginEndPos = (chunkMarginEndPos > contigEndPos ? contigEndPos : chunkMarginEndPos);

        BamChunk *chunk = bamChunk_construct2(contig, (chunkMargin > i ? 0 : i - chunkMargin), i, chunkEndPos,
                                              chunkMarginEndPos, parent);
        stList_append(dest, chunk);
        chunkCount++;
    }
    return chunkCount;
}


/*
 * These handle construction of the BamChunk object, by iterating through the bam (must be sorted), and finds the
 * first and last aligned location on each contig.  Then it generates a list of chunks based off of these positions,
 * with sizes determined by the parameters.
 */
BamChunker *bamChunker_construct(char *bamFile, PolishParams *params) {
    return bamChunker_construct2(bamFile, NULL, params);
}
BamChunker *bamChunker_construct2(char *bamFile, char *region, PolishParams *params) {

    // are we doing region filtering?
    bool filterByRegion = false;
    char regionContig[128] = "";
    int regionStart = 0;
    int regionEnd = 0;
    if (region != NULL) {
        int scanRet = sscanf(region, "%[^:]:%d-%d", regionContig, &regionStart, &regionEnd);
        if (scanRet != 3 || strlen(regionContig) == 0) {
            st_errAbort("Region in unexpected format (expected %%s:%%d-%%d): %s", region);
        } else if (regionStart < 0 || regionEnd <= 0 || regionEnd <= regionStart) {
            st_errAbort("Start and end locations in region must be positive, start must be less than end: %s", region);
        }
        filterByRegion = true;
    }

    // standard parameters
    uint64_t chunkSize = params->chunkSize;
    uint64_t chunkBoundary = params->chunkBoundary;
    bool includeSoftClip = params->includeSoftClipping;

    // the chunker we're building
    BamChunker *chunker = malloc(sizeof(BamChunker));
    chunker->bamFile = stString_copy(bamFile);
    chunker->chunkSize = chunkSize;
    chunker->chunkBoundary = chunkBoundary;
    chunker->includeSoftClip = includeSoftClip;
    chunker->params = params;
    chunker->chunks = stList_construct3(0,(void*)bamChunk_destruct);
    chunker->chunkCount = 0;
    chunker->itorIdx = -1;

    // open bamfile
    samFile *in = hts_open(bamFile, "r");
    if (in == NULL)
        st_errAbort("ERROR: Cannot open bam file %s\n", bamFile);
    hts_idx_t *idx = sam_index_load(in, bamFile); // load index (just to verify early that it exists)
    if (idx == NULL)
        st_errAbort("ERROR: Missing index for bam file %s\n", bamFile);
    hts_idx_destroy(idx);
    bam_hdr_t *bamHdr = sam_hdr_read(in);
    bam1_t *aln = bam_init1();

    // list of chunk boundaries
    char *currentContig = NULL;
    int64_t contigStartPos = 0;
    int64_t contigEndPos = 0;

    // get all reads
    // there is probably a better way (bai?) to find min and max aligned positions (which we need for chunk divisions)
    while(sam_read1(in,bamHdr,aln) > 0) {

        // basic filtering (no read length, no cigar)
        if (aln->core.l_qseq <= 0) continue;
        if (aln->core.n_cigar == 0) continue;

        //data
        char *chr = bamHdr->target_name[aln->core.tid];
        int64_t start_softclip = 0;
        int64_t end_softclip = 0;
        int64_t alnReadLength = getAlignedReadLength3(aln, &start_softclip, &end_softclip, FALSE);
        int64_t alnStartPos = aln->core.pos + 1;
        int64_t alnEndPos = alnStartPos + alnReadLength;

        // does this belong in our chunk?
        if (alnReadLength <= 0) continue;
        if (filterByRegion && (!stString_eq(regionContig, chr) || alnStartPos >= regionEnd || alnEndPos <= regionStart))
            continue;

        // get start and stop position
        int64_t readStartPos = aln->core.pos + 1;           // Left most position of alignment
        int64_t readEndPos = readStartPos + alnReadLength;

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
            // new contig (this should never happen if we're filtering by region)
            assert(!filterByRegion);
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
        if (filterByRegion) {
            contigStartPos = (contigStartPos < regionStart ? regionStart : contigStartPos);
            contigEndPos = (contigEndPos > regionEnd ? regionEnd : contigEndPos);
        }
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


// This structure holds the bed information
// TODO rewrite the code to just use a void*
typedef struct samview_settings {
    void* bed;
} samview_settings_t;


#define DEFAULT_ALIGNMENT_SCORE 10

/*
 * This generates a set of BamChunkReads (and alignments to the reference) from a BamChunk.  The BamChunk describes
 * positional information within the bam, from which the reads should be extracted.  The bam must be indexed.  Reads
 * which overlap the ends of the chunk are truncated.  A parameter in the BamChunk's parameters determines whether
 * softclipped portions of the reads should be included.
 */
uint32_t convertToReadsAndAlignments(BamChunk *bamChunk, stList *reads, stList *alignments) {
    // sanity check
    assert(stList_length(reads) == 0);
    assert(stList_length(alignments) == 0);

    // prep
    int64_t chunkStart = bamChunk->chunkBoundaryStart;
    int64_t chunkEnd = bamChunk->chunkBoundaryEnd;
    bool includeSoftClip = bamChunk->parent->params->includeSoftClipping;
    char *bamFile = bamChunk->parent->bamFile;
    char *contig = bamChunk->refSeqName;
    uint32_t savedAlignments = 0;

    // prep for index (not entirely sure what all this does.  see samtools/sam_view.c
    int filter_state = ALL, filter_op = 0;
    int result;
    samview_settings_t settings = { .bed = NULL };
    char* region[1] = {};
    region[0] = stString_print("%s:%d-%d", bamChunk->refSeqName, bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
    settings.bed = bed_hash_regions(settings.bed, region, 0, 1, &filter_op); //insert(1) or filter out(0) the regions from the command line in the same hash table as the bed file
    if (!filter_op) filter_state = FILTERED;
    int regcount = 0;
    hts_reglist_t *reglist = bed_reglist(settings.bed, filter_state, &regcount);
    if(!reglist) {
        st_errAbort("ERROR: Could not create list of regions for read conversion");
    }

    // file initialization
    samFile *in = NULL;
    hts_idx_t *idx = NULL;
    hts_itr_multi_t *iter = NULL;
    // bam file
    if ((in = hts_open(bamFile, "r")) == 0) {
        st_errAbort("ERROR: Cannot open bam file %s\n", bamFile);
    }
    // bam index
    if ((idx = sam_index_load(in, bamFile)) == 0) {
        st_errAbort("ERROR: Cannot open index for bam file %s\n", bamFile);
    }
    // header  //todo samFile *in = hts_open(bamFile, "r");
    bam_hdr_t *bamHdr = sam_hdr_read(in);
    // read object
    bam1_t *aln = bam_init1();
    // iterator for region
    if ((iter = sam_itr_regions(idx, bamHdr, reglist, regcount)) == 0) {
        st_errAbort("ERROR: Cannot open iterator for region %s for bam file %s\n", region[0], bamFile);
    }

    // fetch alignments //todo while(sam_read1(in,bamHdr,aln) > 0) {
    while ((result = sam_itr_multi_next(in, iter, aln)) >= 0) {
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
        stList *cigRepr = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);

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
        for (uint32_t i = 0; i <= alnReadLength; i++) {
            // handles cases where last alignment is an insert or last is match
            if (cig_idx == aln->core.n_cigar) break;

            // do we need the next cigar operation?
            if (currPosInOp == 0) {
                cigarOp = cigar[cig_idx] & BAM_CIGAR_MASK;
                cigarNum = cigar[cig_idx] >> BAM_CIGAR_SHIFT;
            }

            // handle current character
            if (cigarOp == BAM_CMATCH || cigarOp == BAM_CEQUAL || cigarOp == BAM_CDIFF) {
                if (cigarIdxInRef >= chunkStart && cigarIdxInRef < chunkEnd) {
                    stList_append(cigRepr, stIntTuple_construct3(cigarIdxInRef + refCigarModification,
                                                                 cigarIdxInSeq + seqCigarModification,
                                                                 DEFAULT_ALIGNMENT_SCORE));
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
        //assert(cigarIdxInRef == alnEndPos);  //does not include soft clip
        //assert(cigarIdxInSeq == readEndIdxInChunk - (includeSoftClip ? end_softclip : 0));

        // get sequence positions
        int64_t seqLen = alignedReadLength;

        // modify start indices
        int64_t readStartIdxInChunk = firstNonSoftclipAlignedReadIdxInChunk;
        if (firstNonSoftclipAlignedReadIdxInChunk != 0) {
            // the aligned portion spans chunkStart, so no softclipped bases are included
            readStartIdxInChunk += start_softclip;
        } else if (!includeSoftClip) {
            // configured to not handle softclipped bases
            readStartIdxInChunk += start_softclip;
        } else if (alnStartPos - start_softclip <= chunkStart) {
            // configured to handle softclipped bases; softclipped bases span chunkStart
            int64_t includedSoftclippedBases = alnStartPos - chunkStart;
            seqLen += includedSoftclippedBases;
            readStartIdxInChunk += (start_softclip - includedSoftclippedBases);
        } else {
            // configured to handle softclipped bases; softclipped bases all occur after chunkStart
            seqLen += start_softclip;
            readStartIdxInChunk = 0;
        }

        // modify end indices
        int64_t readEndIdxInChunk = readStartIdxInChunk + seqLen;
        if (alnEndPos < chunkEnd && includeSoftClip) {
            // all other cases mean we don't need to handle softclip (by config or aln extends past chunk end)
            if (alnEndPos + end_softclip <= chunkEnd) {
                // all softclipped bases fit in chunk
                readEndIdxInChunk += end_softclip;
                seqLen += end_softclip;
            } else {
                // softclipping spands chunkEnd
                int64_t includedSoftclippedBases = chunkEnd - alnEndPos;
                seqLen += includedSoftclippedBases;
                readEndIdxInChunk += includedSoftclippedBases;
            }
        }

        // get sequence - all data we need is encoded in readStartIdxInChunk (start), readEnd idx, and seqLen
        char *seq = st_calloc(seqLen + 1, sizeof(char));
        uint8_t *seqBits = bam_get_seq(aln);
        int64_t idxInOutputSeq = 0;
        int64_t idxInBamRead = readStartIdxInChunk;
        while (idxInBamRead < readEndIdxInChunk) {
            seq[idxInOutputSeq] = seq_nt16_str[bam_seqi(seqBits, idxInBamRead)];
            idxInBamRead++;
            idxInOutputSeq++;
        }
        seq[seqLen] = '\0';

        // get sequence qualities (if exists)
        char *readName = stString_copy(bam_get_qname(aln));
        uint8_t *qualBits = bam_get_qual(aln);
        uint8_t *qual = NULL;
        if (qualBits[0] != 0xff) { //inital score of 255 means qual scores are unavailable
            idxInOutputSeq = 0;
            idxInBamRead = readStartIdxInChunk;
            qual = st_calloc(seqLen, sizeof(uint8_t));
            while (idxInBamRead < readEndIdxInChunk) {
                qual[idxInOutputSeq] = qualBits[idxInBamRead];
                idxInBamRead++;
                idxInOutputSeq++;

            }
            assert(idxInOutputSeq == strlen(seq));
        };

        // sanity check
        assert(stIntTuple_get((stIntTuple *)stList_peek(cigRepr), 1) < strlen(seq));

        // save to read
        bool forwardStrand = !bam_is_rev(aln);
        BamChunkRead *chunkRead = bamChunkRead_construct2(readName, seq, qual, forwardStrand, bamChunk);
        stList_append(reads, chunkRead);
        stList_append(alignments, cigRepr);
        savedAlignments++;
    }
    // the status from "get reads from iterator"
    if (result < -1) {
        st_errAbort("ERROR: Retrieval of region %d failed due to truncated file or corrupt BAM index file\n", iter->curr_tid);
    }

    // close it all down
    hts_itr_multi_destroy(iter);
    hts_idx_destroy(idx);
    free(region[0]);
    bed_destroy(settings.bed);
    bam_hdr_destroy(bamHdr);
    bam_destroy1(aln);
    sam_close(in);

    return savedAlignments;
}

/*
 * Helper functions for reading alignment likelihood files.
 */


void appendProbsToList(stList *probabilityList, uint8_t pA, uint8_t pC, uint8_t pG, uint8_t pT, uint8_t pGap) {
    uint8_t *aPtr = calloc(1, sizeof(uint8_t));
    uint8_t *cPtr = calloc(1, sizeof(uint8_t));
    uint8_t *gPtr = calloc(1, sizeof(uint8_t));
    uint8_t *tPtr = calloc(1, sizeof(uint8_t));
    uint8_t *gapPtr = calloc(1, sizeof(uint8_t));
    *aPtr = pA;
    *cPtr = pC;
    *gPtr = pG;
    *tPtr = pT;
    *gapPtr = pGap;
    stList_append(probabilityList, aPtr);
    stList_append(probabilityList, cPtr);
    stList_append(probabilityList, gPtr);
    stList_append(probabilityList, tPtr);
    stList_append(probabilityList, gapPtr);
}

stProfileSeq* getProfileSequenceFromSingleNuclProbFile(char *signalAlignReadLocation, char *readName,
                                                       stBaseMapper *baseMapper, stRPHmmParameters *params) {
    // get signalAlign file
    FILE *fp = fopen(signalAlignReadLocation,"r");

    // for scanning the file
    int fieldSize = 2048;
    char *line = calloc(2048, sizeof(char));
    char *chromStr = calloc(fieldSize, sizeof(char));
    char *refPosStr = calloc(fieldSize, sizeof(char));
    char *pAStr = calloc(fieldSize, sizeof(char));
    char *pCStr = calloc(fieldSize, sizeof(char));
    char *pGStr = calloc(fieldSize, sizeof(char));
    char *pTStr = calloc(fieldSize, sizeof(char));
    char *pGapStr = calloc(fieldSize, sizeof(char));

    // for handling the data
    int64_t refPos;
    uint8_t pA;
    uint8_t pC;
    uint8_t pG;
    uint8_t pT;
    uint8_t pGap;
    uint8_t *aPtr;
    uint8_t *cPtr;
    uint8_t *gPtr;
    uint8_t *tPtr;
    uint8_t *gapPtr;

    // parse header
    while(!feof(fp)) {
        fscanf( fp, "%[^\n]\n", line);
        if (line[0] == '#') {
            if (line[1] == '#') continue;
            if (strcmp(line, "#CHROM\tPOS\tpA\tpC\tpG\tpT\tp_") != 0) {
                st_errAbort("SignalAlign output file %s has unexpected header format: %s",
                            signalAlignReadLocation, line);
            } else {
                break;
            }
        }
    }

    // get probabilities
    stList* probabilityList = stList_construct3(0, free);
    uint64_t firstReadPos = 0;
    uint64_t lastReadPos = 0;
    int64_t randomSeed = st_randomInt64(0,3);
    while(!feof(fp)) {
        // Scan
        fscanf( fp, "%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\n]\n",
                chromStr, refPosStr, pAStr, pCStr, pGStr, pTStr, pGapStr);

        // Get reference position
        refPos = atoi(refPosStr);

        // Check for gaps todo this might actually be a bug or something in signalAlign
        while (firstReadPos != 0 && refPos > lastReadPos + 1) {
            appendProbsToList(probabilityList, ALPHABET_MIN_PROB, ALPHABET_MIN_PROB, ALPHABET_MIN_PROB,
                              ALPHABET_MIN_PROB, ALPHABET_MAX_PROB);
            lastReadPos++;
        }

        // Check for inserts
        if (firstReadPos != 0 && refPos < lastReadPos + 1) continue;

        // Update first and last read position
        if (firstReadPos == 0) firstReadPos = refPos;
        lastReadPos = refPos;

        // Get probabilities and save
        pA = (uint8_t) (ALPHABET_MAX_PROB * atof(pAStr));
        pC = (uint8_t) (ALPHABET_MAX_PROB * atof(pCStr));
        pG = (uint8_t) (ALPHABET_MAX_PROB * atof(pGStr));
        pT = (uint8_t) (ALPHABET_MAX_PROB * atof(pTStr));
        pGap = (uint8_t) (ALPHABET_MAX_PROB * atof(pGapStr));

        // We need all integer probs to sum to MAX_PROB todo is there a way to do this better?
        while ((pA + pC + pG + pT + pGap) > ALPHABET_MAX_PROB) {
            if ((pA + pC + pG + pT + pGap) == 0) {
                break;
            }
            switch (randomSeed++ % 5) {
                case 0: if (pA != 0) pA--; break;
                case 1: if (pC != 0) pC--; break;
                case 2: if (pG != 0) pG--; break;
                case 3: if (pT != 0) pT--; break;
                case 4: if (pGap != 0) pGap--; break;
                default: assert(FALSE);
            }
        }
        while ((pA + pC + pG + pT + pGap) < ALPHABET_MAX_PROB) {
            if ((pA + pC + pG + pT + pGap) == 0) {
                break;
            }
            switch (randomSeed++ % 5) {
                case 0: if (pA != 0) pA++; break;
                case 1: if (pC != 0) pC++; break;
                case 2: if (pG != 0) pG++; break;
                case 3: if (pT != 0) pT++; break;
                case 4: if (pGap != 0) pGap++; break;
                default: assert(FALSE);
            }
        }

        // Save the values in a list todo this is not particularly efficient
        appendProbsToList(probabilityList, pA, pC, pG, pT, pGap);
    }
    // Now we're done with the file
    fclose(fp);

    // Create empty profile sequence
    uint64_t readLength = lastReadPos - firstReadPos + 1;
    stProfileSeq *pSeq = stProfileSeq_constructEmptyProfile(chromStr, readName, firstReadPos + 1, readLength);

    // Copy probabilities over
    uint64_t position = 0;
    stListIterator *itor = stList_getIterator(probabilityList);
    while (position < readLength) {

        // Get the locations of the probabilities
        aPtr = stList_getNext(itor);
        cPtr = stList_getNext(itor);
        gPtr = stList_getNext(itor);
        tPtr = stList_getNext(itor);
        gapPtr = stList_getNext(itor);

        // Assign the probabilities
        pSeq->profileProbs[position * ALPHABET_SIZE + stBaseMapper_getValueForChar(baseMapper, 'A')] =  *aPtr;
        pSeq->profileProbs[position * ALPHABET_SIZE + stBaseMapper_getValueForChar(baseMapper, 'C')] =  *cPtr;
        pSeq->profileProbs[position * ALPHABET_SIZE + stBaseMapper_getValueForChar(baseMapper, 'G')] =  *gPtr;
        pSeq->profileProbs[position * ALPHABET_SIZE + stBaseMapper_getValueForChar(baseMapper, 'T')] =  *tPtr;

        if (params->gapCharactersForDeletions) {
            // This assumes gap character is the last character in the alphabet given
            pSeq->profileProbs[position * ALPHABET_SIZE + (ALPHABET_SIZE - 1)] = *gapPtr;
        } else {
            pSeq->profileProbs[position * ALPHABET_SIZE + (ALPHABET_SIZE - 1)] = ALPHABET_MIN_PROB;
        }

        position++;
    }
    // We should have nothing left over in the list
    if (stList_getNext(itor) != NULL) {
        st_errAbort("Probability list has %d extra elements, with read length %d for file %s",
                    stList_length(probabilityList) - 5 * position, readLength, signalAlignReadLocation);
    }
    // Sanity check on the number of modifications to the probabilities
    // We only modify probability of bases with some probability, so to fix a rounding error, we should at worst have
    //  to make 4 modifications per location
    if (randomSeed > (4 * readLength)) {
        st_logDebug("\t\tNeeded average of %f modifications to base probs to ensure proper total probability for %s\n",
                    (1.0 * randomSeed / readLength), readName);
    }

    stList_destruct(probabilityList);
    free(line);
    free(chromStr);
    free(refPosStr);
    free(pAStr);
    free(pCStr);
    free(pGStr);
    free(pTStr);

    return pSeq;
}



/* Parse reads within an input interval of a reference sequence of a bam file
 * and create a list of profile sequences by turning characters into profile probabilities.
 *
 * In future, maybe use mapq scores to adjust profile (or posterior probabilities for
 * signal level alignments).
 * */
int64_t parseReads(stList *profileSequences, char *bamFile, stBaseMapper *baseMapper, stRPHmmParameters *params) {
    return parseReadsWithSingleNucleotideProbs(profileSequences, bamFile, baseMapper, params, NULL, false);
}

int64_t parseReadsWithSingleNucleotideProbs(stList *profileSequences, char *bamFile, stBaseMapper *baseMapper,
                                            stRPHmmParameters *params, char *singleNuclProbDirectory,
                                            bool onlySingleNuclProb) {
    if (singleNuclProbDirectory != NULL) {
        st_logInfo("\tModifying probabilities from single nucleotide probability files in %s\n",
                   singleNuclProbDirectory);
    }

    samFile *in = hts_open(bamFile, "r");
    if (in == NULL) {
        st_errAbort("ERROR: Cannot open bam file %s\n", bamFile);
        return -1;
    }
    bam_hdr_t *bamHdr = sam_hdr_read(in);
    bam1_t *aln = bam_init1();

    int64_t readCount = 0;
    int64_t singleNuclProbReadCount = 0;
    int64_t bamReadCount = 0;
    int64_t profileCount = 0;
    int64_t missingSingleNuclProbReads = 0;
    int64_t filteredReads = 0;
    int64_t filteredReads_flag = 0;
    int64_t filteredReads_mapq = 0;

    while(sam_read1(in,bamHdr,aln) > 0) {
        stProfileSeq *pSeq = NULL;

        int64_t pos = aln->core.pos+1;                      // Left most position of alignment
        char *chr = bamHdr->target_name[aln->core.tid] ;    // Contig name (chromosome)
        int64_t len = aln->core.l_qseq;                     // Length of the read.
        uint8_t *seq = bam_get_seq(aln);                    // DNA sequence
        char *readName = bam_get_qname(aln);
        uint32_t *cigar = bam_get_cigar(aln);

        if (aln->core.l_qseq <= 0) {
            filteredReads++;
            continue;
        }

        // Filter out any reads with specified flags
        if((aln->core.flag & params->filterAReadWithAnyOneOfTheseSamFlagsSet) > 0) {
            filteredReads++;
            filteredReads_flag++;
            continue;
        }

        // If there isn't a cigar string, don't bother including the read, since we don't
        // know how it aligns
        if (aln->core.n_cigar == 0) {
            filteredReads++;
            continue;
        }

        // If the mapq score is less than the given threshold, filter it out
        if (aln->core.qual < params->mapqFilter) {
            filteredReads++;
            filteredReads_mapq++;
            continue;
        }

        // Tracks how many reads there were
        readCount++;

        // Should we read from the signalAlign directory?
        if (singleNuclProbDirectory != NULL) {

            // Get signalAlign file (if exists)
            char *singleNuclProbReadLocation = stString_print("%s/%s.tsv", singleNuclProbDirectory, readName);
            if (access(singleNuclProbReadLocation, F_OK) == -1) {
                // Could not find the read file
                missingSingleNuclProbReads++;
            } else {
                // Found the read file
                pSeq = getProfileSequenceFromSingleNuclProbFile(singleNuclProbReadLocation, readName, baseMapper, params);
                singleNuclProbReadCount++;

                // We have a profile, so save it
                stList_append(profileSequences, pSeq);
                profileCount++;
            }
            free(singleNuclProbReadLocation);

            // If we found a SA file or if we don't want missing reads
            if (pSeq != NULL || onlySingleNuclProb) {
                continue;
            }
        }

        int64_t start_read = 0;
        int64_t end_read = 0;
        int64_t trueLength = getAlignedReadLength2(aln, &start_read, &end_read);

        if (trueLength <= 0) {
            filteredReads++;
            continue;
        }

        // Create empty profile sequence
        pSeq = stProfileSeq_constructEmptyProfile(chr, readName, pos, trueLength);

        // Variables to keep track of position in sequence / cigar operations
        int64_t cig_idx = 0;
        int64_t currPosInOp = 0;
        int64_t cigarOp = -1;
        int64_t cigarNum = -1;
        int64_t idxInSeq = start_read;

        // For each position turn character into profile probability
        // As is, this makes the probability 1 for the base read in, and 0 otherwise
        for (uint32_t i = 0; i < trueLength; i++) {

            if (currPosInOp == 0) {
                cigarOp = cigar[cig_idx] & BAM_CIGAR_MASK;
                cigarNum = cigar[cig_idx] >> BAM_CIGAR_SHIFT;
            }
            if (cigarOp == BAM_CMATCH || cigarOp == BAM_CEQUAL || cigarOp == BAM_CDIFF) {
                int64_t b = stBaseMapper_getValueForChar(baseMapper, seq_nt16_str[bam_seqi(seq, idxInSeq)]);
                pSeq->profileProbs[i * ALPHABET_SIZE + b] = ALPHABET_MAX_PROB;
                idxInSeq++;
            } else if (cigarOp == BAM_CDEL || cigarOp == BAM_CREF_SKIP) {
                // Only add a gap character when that param is on
                if (params->gapCharactersForDeletions) {
                    // This assumes gap character is the last character in the alphabet given
                    pSeq->profileProbs[i * ALPHABET_SIZE + (ALPHABET_SIZE - 1)] = ALPHABET_MAX_PROB;
                }
            } else if (cigarOp == BAM_CINS) {
                // Currently, ignore insertions
                idxInSeq++;
                i--;
            } else if (cigarOp == BAM_CSOFT_CLIP || cigarOp == BAM_CHARD_CLIP || cigarOp == BAM_CPAD) {
                // Nothing really to do here. skip to next cigar operation
                currPosInOp = cigarNum - 1;
                i--;
            } else {
                st_logCritical("Unidentifiable cigar operation\n");
            }

            currPosInOp++;
            if (currPosInOp == cigarNum) {
                cig_idx++;
                currPosInOp = 0;
            }
        }
        bamReadCount++;

        // Save profile seq
        if (pSeq->length > 0) {
            profileCount++;
            stList_append(profileSequences, pSeq);
        }

    }

    // Log signal align usage
    if (singleNuclProbDirectory != NULL) {
        if (missingSingleNuclProbReads > 0) {
            st_logInfo("\t%d/%d reads were missing single nucleotide probability file\n", missingSingleNuclProbReads, readCount);
        }
        st_logInfo("\tOf %d total reads: %d were loaded from single nucleotide probability data, and %d were from the bam\n",
                   profileCount, singleNuclProbReadCount, bamReadCount);

    }

    // Log filtering actions
    if(st_getLogLevel() == debug) {
        char *samFlagBitString = intToBinaryString(params->filterAReadWithAnyOneOfTheseSamFlagsSet);
        st_logDebug("\tFiltered %" PRIi64 " reads with either missing cigar lines, "
                                          "\n\t\tlow mapq scores (filtered %d reads with scores less than %d), "
                                          "\n\t\tand undesired sam flags "
                                          "(filtered %d reads with sam flags being filtered on: %s)\n",
                filteredReads, filteredReads_mapq, params->mapqFilter, filteredReads_flag, samFlagBitString);
        free(samFlagBitString);
    }

    // Sanity check (did we accidentally save profile sequences twice?)
    assert(stList_length(profileSequences) <= readCount);

    bam_hdr_destroy(bamHdr);
    bam_destroy1(aln);
    sam_close(in);

    return profileCount;
}


void writeHaplotypedSam(char *bamInFile, char *bamOutBase, stReadHaplotypePartitionTable *readHaplotypePartitions,
                        char *marginPhaseTag) {
    /*
     * Write out haplotyped sam file
     */

    // Prep
    char *haplotypedSamFile = stString_print("%s.sam", bamOutBase);

    // File management
    samFile *in = hts_open(bamInFile, "r");
    if (in == NULL) {
        st_errAbort("ERROR: Cannot open bam file %s\n", bamInFile);
    }
    bam_hdr_t *bamHdr = bam_hdr_read(in->fp.bgzf);
    bam1_t *aln = bam_init1();

    int r;
    st_logDebug("\tWriting haplotype output to: %s \n", haplotypedSamFile);

    samFile *out = hts_open(haplotypedSamFile, "w");
    r = sam_hdr_write(out, bamHdr);

    // Read in input file, write out each read to one sam file
    int32_t readCountH1 = 0;
    int32_t readCountH2 = 0;
    int32_t readCountFiltered = 0;
    char *haplotypeString;
    while(sam_read1(in,bamHdr,aln) > 0) {

        char *readName = bam_get_qname(aln);
        if (marginPhaseTag != NULL) {
            bam_aux_append(aln, MARGIN_PHASE_TAG, 'Z', (int)strlen(marginPhaseTag) + 1, (uint8_t*)marginPhaseTag);
        }

        stReadHaplotypeSequence *readHaplotypes = hashtable_search(readHaplotypePartitions, readName);
        if (readHaplotypes == NULL) {
            haplotypeString = stReadHaplotypeSequence_toStringEmpty();
            bam_aux_append(aln, HAPLOTYPE_TAG, 'Z', (int)strlen(haplotypeString) + 1, (uint8_t*)haplotypeString);
            r = sam_write1(out, bamHdr, aln);
            readCountFiltered++;
        } else {
            haplotypeString = stReadHaplotypeSequence_toString(readHaplotypes);
            bam_aux_append(aln, HAPLOTYPE_TAG, 'Z', (int)strlen(haplotypeString) + 1, (uint8_t*)haplotypeString);
            r = sam_write1(out, bamHdr, aln);

            // Document based on last recorded haplotype
            while (readHaplotypes->next != NULL) {readHaplotypes = readHaplotypes->next;}
            if (readHaplotypes->haplotype == 1)
                readCountH1++;
            else
                readCountH2++;
        }
        free(haplotypeString);
    }
    st_logInfo("\tSAM read counts:\n\t\thap1: %d\thap2: %d\tfiltered out: %d \n", readCountH1, readCountH2, readCountFiltered);

    // Cleanup
    bam_destroy1(aln);
    bam_hdr_destroy(bamHdr);
    sam_close(in);
    sam_close(out);
    free(haplotypedSamFile);
}

void writeSplitSams(char *bamInFile, char *bamOutBase, stReadHaplotypePartitionTable *readHaplotypePartitions,
                    char *marginPhaseTag) {
    /*
     * Write out sam files with reads in each split based on which haplotype partition they are in.
     */

    // Prep
    char *haplotype1SamOutFile = stString_print("%s.1.sam", bamOutBase);
    char *haplotype2SamOutFile = stString_print("%s.2.sam", bamOutBase);
    char *unmatchedSamOutFile = stString_print("%s.0.sam", bamOutBase);

    // File management
    samFile *in = hts_open(bamInFile, "r");
    if (in == NULL) {
        st_errAbort("ERROR: Cannot open bam file %s\n", bamInFile);
    }
    bam_hdr_t *bamHdr = sam_hdr_read(in);
    bam1_t *aln = bam_init1();

    int r;
    st_logDebug("\tWriting haplotype output to: %s, %s, and %s \n", haplotype1SamOutFile,
                haplotype2SamOutFile, unmatchedSamOutFile);
    samFile *out1 = hts_open(haplotype1SamOutFile, "w");
    r = sam_hdr_write(out1, bamHdr);

    samFile *out2 = hts_open(haplotype2SamOutFile, "w");
    r = sam_hdr_write(out2, bamHdr);

    samFile *outUnmatched = hts_open(unmatchedSamOutFile, "w");
    r = sam_hdr_write(outUnmatched, bamHdr);


    // Read in input file, write out each read to one sam file
    int32_t readCountH1 = 0;
    int32_t readCountH2 = 0;
    int32_t readCountFiltered = 0;
    char *haplotypeString;
    while(sam_read1(in,bamHdr,aln) > 0) {

        char *readName = bam_get_qname(aln);
        if (marginPhaseTag != NULL) {
            bam_aux_append(aln, MARGIN_PHASE_TAG, 'Z', (int)strlen(marginPhaseTag) + 1, (uint8_t*)marginPhaseTag);
        }

        stReadHaplotypeSequence *readHaplotypes = hashtable_search(readHaplotypePartitions, readName);
        if (readHaplotypes == NULL) {
            haplotypeString = stReadHaplotypeSequence_toStringEmpty();
            bam_aux_append(aln, HAPLOTYPE_TAG, 'Z', (int)strlen(haplotypeString) + 1, (uint8_t*)haplotypeString);
            r = sam_write1(outUnmatched, bamHdr, aln);
            readCountFiltered++;
        } else {
            haplotypeString = stReadHaplotypeSequence_toString(readHaplotypes);
            bam_aux_append(aln, HAPLOTYPE_TAG, 'Z', (int)strlen(haplotypeString) + 1, (uint8_t*)haplotypeString);

            // Document based on last recorded haplotype
            while (readHaplotypes->next != NULL) {readHaplotypes = readHaplotypes->next;}
            if (readHaplotypes->haplotype == 1) {
                r = sam_write1(out1, bamHdr, aln);
                readCountH1++;
            } else {
                r = sam_write1(out2, bamHdr, aln);
                readCountH2++;
            }
        }
        free(haplotypeString);

    }
    st_logInfo("\tSAM read counts:\n\t\thap1: %d\thap2: %d\tfiltered out: %d \n", readCountH1, readCountH2, readCountFiltered);

    // Cleanup
    bam_destroy1(aln);
    bam_hdr_destroy(bamHdr);
    sam_close(in);
    sam_close(out1);
    sam_close(out2);
    sam_close(outUnmatched);
    free(haplotype1SamOutFile);
    free(haplotype2SamOutFile);
    free(unmatchedSamOutFile);
}

bcf_hdr_t* writeVcfHeader(vcfFile *out, stList *genomeFragments, char *referenceName) {
    /*
     * Write the header of a vcf file.
     */

    bcf_hdr_t *hdr = bcf_hdr_init("w");
    kstring_t str = {0,0,NULL};

    // Generic info
    str.l = 0;
    ksprintf(&str, "##marginPhase=htslib-%s\n", hts_version());
    bcf_hdr_append(hdr, str.s);

    // Reference file used
    str.l = 0;
    ksprintf(&str, "##reference=file://%s\n", referenceName);
    bcf_hdr_append(hdr, str.s);

    // Contigs
    // TODO: assert unique fragments, get full chrom length
    for(int64_t i=0; i<stList_length(genomeFragments); i++) {
        stRPHmm *hmm = stList_get(genomeFragments, i);
        str.l = 0;
        ksprintf(&str, "##contig=<ID=%s>\n", hmm->referenceName); //hmm->referenceName is the chrom
        bcf_hdr_append(hdr, str.s);
    }

    // INFO fields
    str.l = 0;
    ksprintf(&str, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">");
    bcf_hdr_append(hdr, str.s);
    str.l = 0;
    ksprintf(&str, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate allele counts\">");
    bcf_hdr_append(hdr, str.s);

    // FORMAT fields
    str.l = 0;
    ksprintf(&str, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_append(hdr, str.s);
    str.l = 0;
    ksprintf(&str, "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype Likelihoods\">");
    bcf_hdr_append(hdr, str.s);
    str.l = 0;
    ksprintf(&str, "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set for GT\">");
    bcf_hdr_append(hdr, str.s);
    str.l = 0;
    ksprintf(&str, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">");
    bcf_hdr_append(hdr, str.s);
    str.l = 0;
    ksprintf(&str, "##FORMAT=<ID=MPI,Number=1,Type=Integer,Description=\"MarginPhase\tIdentifier\">");
    bcf_hdr_append(hdr, str.s);

    // Samples
    bcf_hdr_add_sample(hdr, "SMPL1"); //todo change the sample name
    bcf_hdr_add_sample(hdr, NULL);

    // Write header
    bcf_hdr_write(out, hdr);

    // Cleanup
    free(str.s);
    return hdr;
}


void writeIndelVariant(int32_t *gt_info, bcf_hdr_t *bcf_hdr, bcf1_t *bcf_rec, stGenomeFragment *gF,
                       stBaseMapper *baseMapper, char *referenceSeq, char refChar,
                       char h1AlphChar, char h2AlphChar, vcfFile *out, int64_t *index,
                       int32_t  *phaseSet, int32_t  *ps_info, int32_t *ac_info, float *gl_info,
                       bool *firstVariantInPhaseBlock, bool gvcf) {
    /*
     * Write a vcf record for a variant involving an insertion or deletion.
     */

    // Initialize strings
    kstring_t refstr = {0, 0, NULL};
    kputc(refChar, &refstr);
    kstring_t hap1str = {0, 0, NULL};
    kputc(h1AlphChar, &hap1str);
    kstring_t hap2str = {0, 0, NULL};
    kputc(h2AlphChar, &hap2str);

    uint64_t refCharVal = stBaseMapper_getValueForChar(baseMapper, refChar);
    uint64_t h1AlphVal = stBaseMapper_getValueForChar(baseMapper, h1AlphChar);
    uint64_t h2AlphVal = stBaseMapper_getValueForChar(baseMapper, h2AlphChar);
    uint64_t secondRefVal = stBaseMapper_getValueForChar(baseMapper, toupper(referenceSeq[*index + gF->refStart]));

    // Determine the sequence of the indel variant & reference sequence
    int64_t j = 1;
    int64_t i = *index;
    while (i + j < gF->length &&
           (stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString1[i + j]) == '-' ||
            stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[i + j]) == '-')) {

        int nextRefChar = toupper(referenceSeq[i + j + gF->refStart - 1]);
        kputc(nextRefChar, &refstr);
        if (stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString1[i + j]) != '-') {
            kputc(stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString1[i + j]), &hap1str);
        }
        if (stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[i + j]) != '-') {
            kputc(stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[i + j]), &hap2str);
        }
        j++;
    }

    if (strcmp(hap1str.s, hap2str.s) == 0) {
        // Homozygous alleles 1/1
        // Ref allele will be the reference string constructed
        if (*firstVariantInPhaseBlock) {
            gt_info[0] = bcf_gt_unphased(1);
            gt_info[1] = bcf_gt_unphased(1);
            *firstVariantInPhaseBlock = false;
            *phaseSet = bcf_rec->pos+1;
        } else {
            gt_info[0] = bcf_gt_phased(1);
            gt_info[1] = bcf_gt_phased(1);
        }
        kputc(',', &refstr);
        kputs(hap1str.s, &refstr);
        bcf_update_alleles_str(bcf_hdr, bcf_rec, refstr.s);

        // Genotype likelihoods (AA, AB, BB, AC, BC, CC)
        gl_info[0] = gF->genotypeLikelihoods[i][secondRefVal*ALPHABET_SIZE+secondRefVal];
        gl_info[1] = gF->genotypeLikelihoods[i][secondRefVal*ALPHABET_SIZE+(ALPHABET_SIZE-1)];
        gl_info[2] = gF->genotypeLikelihoods[i][(ALPHABET_SIZE-1)*ALPHABET_SIZE+(ALPHABET_SIZE-1)];
        // Update genotype likelihoods
        bcf_update_format(bcf_hdr, bcf_rec, "GL", gl_info, bcf_hdr_nsamples(bcf_hdr)*3, BCF_HT_REAL);

        // Update allele counts
        ac_info[0] = gF->alleleCountsHap1[i] + gF->alleleCountsHap2[i];
        bcf_update_info(bcf_hdr, bcf_rec, "AC", ac_info, bcf_hdr_nsamples(bcf_hdr), BCF_HT_INT);

    } else if (strcmp(hap1str.s, refstr.s) == 0) {
        // Het 0/1
        if (*firstVariantInPhaseBlock) {
            gt_info[0] = bcf_gt_unphased(0);
            gt_info[1] = bcf_gt_unphased(1);
            *firstVariantInPhaseBlock = false;
            *phaseSet = bcf_rec->pos+1;
        } else {
            gt_info[0] = bcf_gt_phased(0);
            gt_info[1] = bcf_gt_phased(1);
        }
        kputc(',', &hap1str);
        kputs(hap2str.s, &hap1str);
        bcf_update_alleles_str(bcf_hdr, bcf_rec, hap1str.s);

        // Genotype likelihoods (AA, AB, BB, AC, BC, CC)
        gl_info[0] = gF->genotypeLikelihoods[i][secondRefVal*ALPHABET_SIZE+secondRefVal];
        gl_info[1] = gF->genotypeLikelihoods[i][(ALPHABET_SIZE-1)*ALPHABET_SIZE+secondRefVal];
        gl_info[2] = gF->genotypeLikelihoods[i][(ALPHABET_SIZE-1)*ALPHABET_SIZE+(ALPHABET_SIZE-1)];
        // Update genotype likelihoods
        bcf_update_format(bcf_hdr, bcf_rec, "GL", gl_info, bcf_hdr_nsamples(bcf_hdr)*3, BCF_HT_REAL);

        // Update allele counts
        ac_info[0] = gF->alleleCountsHap1[i] + gF->alleleCountsHap2[i];
        bcf_update_info(bcf_hdr, bcf_rec, "AC", ac_info, bcf_hdr_nsamples(bcf_hdr), BCF_HT_INT);

    } else if (strcmp(hap2str.s, refstr.s) == 0){
        // Het 1/0
        if (*firstVariantInPhaseBlock) {
            gt_info[0] = bcf_gt_unphased(1);
            gt_info[1] = bcf_gt_unphased(0);
            *firstVariantInPhaseBlock = false;
            *phaseSet = bcf_rec->pos+1;
        } else {
            gt_info[0] = bcf_gt_phased(1);
            gt_info[1] = bcf_gt_phased(0);
        }
        kputc(',', &hap2str);
        kputs(hap1str.s, &hap2str);
        bcf_update_alleles_str(bcf_hdr, bcf_rec, hap2str.s);

        // Genotype likelihoods (AA, AB, BB, AC, BC, CC)
        gl_info[0] = gF->genotypeLikelihoods[i][secondRefVal*ALPHABET_SIZE+secondRefVal];
        gl_info[1] = gF->genotypeLikelihoods[i][(ALPHABET_SIZE-1)*ALPHABET_SIZE+secondRefVal];
        gl_info[2] = gF->genotypeLikelihoods[i][(ALPHABET_SIZE-1)*ALPHABET_SIZE+(ALPHABET_SIZE-1)];
        // Update genotype likelihoods
        bcf_update_format(bcf_hdr, bcf_rec, "GL", gl_info, bcf_hdr_nsamples(bcf_hdr)*3, BCF_HT_REAL);

        // Update allele counts
        ac_info[0] = gF->allele2CountsHap1[i] + gF->allele2CountsHap2[i];
        bcf_update_info(bcf_hdr, bcf_rec, "AC", ac_info, bcf_hdr_nsamples(bcf_hdr), BCF_HT_INT);

    } else {
        // Het 1/2
        // Neither matched the reference
        if (*firstVariantInPhaseBlock) {
            gt_info[0] = bcf_gt_unphased(1);
            gt_info[1] = bcf_gt_unphased(2);
            *firstVariantInPhaseBlock = false;
            *phaseSet = bcf_rec->pos+1;
        } else {
            gt_info[0] = bcf_gt_phased(1);
            gt_info[1] = bcf_gt_phased(2);
        }
        kputc(',', &refstr);
        kputs(hap1str.s, &refstr);
        kputc(',', &refstr);
        kputs(hap2str.s, &refstr);
        bcf_update_alleles_str(bcf_hdr, bcf_rec, refstr.s);

        // Genotype likelihoods (AA, AB, BB, AC, BC, CC)
        gl_info[0] = gF->genotypeLikelihoods[i][refCharVal*ALPHABET_SIZE+refCharVal];
        gl_info[1] = gF->genotypeLikelihoods[i][refCharVal*ALPHABET_SIZE+h1AlphVal];
        gl_info[2] = gF->genotypeLikelihoods[i][h1AlphVal*ALPHABET_SIZE+h1AlphVal];
        gl_info[3] = gF->genotypeLikelihoods[i][refCharVal*ALPHABET_SIZE+h2AlphVal];
        gl_info[4] = gF->genotypeLikelihoods[i][h1AlphVal*ALPHABET_SIZE+h2AlphVal];
        gl_info[5] = gF->genotypeLikelihoods[i][h2AlphVal*ALPHABET_SIZE+h2AlphVal];
        // Update genotype likelihoods
        bcf_update_format(bcf_hdr, bcf_rec, "GL", gl_info, bcf_hdr_nsamples(bcf_hdr)*6, BCF_HT_REAL);

        // Update allele counts
        ac_info[0] = gF->alleleCountsHap1[i] + gF->alleleCountsHap2[i];
        ac_info[1] = gF->allele2CountsHap1[i] + gF->allele2CountsHap2[i];
        bcf_update_info(bcf_hdr, bcf_rec, "AC", ac_info, bcf_hdr_nsamples(bcf_hdr), BCF_HT_INT);
    }

    // Write out the record
    if (gvcf) {
        // Only write out all extra positions inside indel in gvcf
        bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr) * 2);
        ps_info[0] = *phaseSet;
        bcf_update_format(bcf_hdr, bcf_rec, "PS", ps_info, bcf_hdr_nsamples(bcf_hdr), BCF_HT_INT);
        bcf_write1(out, bcf_hdr, bcf_rec);
        for (int64_t k = 1; k < j; k++) {
            bcf_rec->rid = bcf_hdr_name2id(bcf_hdr, gF->referenceName);
            bcf_rec->pos  = i + gF->refStart - 1 + k;
            // bcf_rec->qual = '.';
            kstring_t str = {0, 0, NULL};
            kputs(".,.", &str);
            bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
            bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr) * 2);
            ps_info[0] = *phaseSet;
            bcf_update_format(bcf_hdr, bcf_rec, "PS", ps_info, bcf_hdr_nsamples(bcf_hdr), BCF_HT_INT);
            bcf_write1(out, bcf_hdr, bcf_rec);
        }
    }
    else {
        bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr) * 2);
        ps_info[0] = *phaseSet;
        bcf_update_format(bcf_hdr, bcf_rec, "PS", ps_info, bcf_hdr_nsamples(bcf_hdr), BCF_HT_INT);
        bcf_write1(out, bcf_hdr, bcf_rec);
    }

    *index += (j - 1);

    free(refstr.s);
    free(hap1str.s);
    free(hap2str.s);
}

void writeHetSite(char h1AlphChar, char h2AlphChar, char refChar,
                  int32_t *phaseSet, bool *firstVariantInPhaseBlock,
                  int32_t *gt_info, bcf1_t *bcf_rec, kstring_t *str,
                  int32_t *ac_info, float *gl_info, stGenomeFragment *gF,
                  stBaseMapper *baseMapper, bcf_hdr_t *bcf_hdr, int64_t i) {
    /*
     * Write out a het site record.
     */
    int refCharVal = stBaseMapper_getValueForChar(baseMapper, refChar);
    uint64_t h1AlphVal = gF->haplotypeString1[i];
    uint64_t h2AlphVal = gF->haplotypeString2[i];

    if (h1AlphChar == refChar) {
        // 0|1
        if (*firstVariantInPhaseBlock) {
            gt_info[0] = bcf_gt_unphased(0);
            gt_info[1] = bcf_gt_unphased(1);
            *firstVariantInPhaseBlock = false;
            *phaseSet = bcf_rec->pos+1;
        } else {
            gt_info[0] = bcf_gt_phased(0);
            gt_info[1] = bcf_gt_phased(1);
        }
        kputc(h1AlphChar, str);
        kputc(',', str);
        kputc(h2AlphChar, str);
        // Allele counts - hap1
        ac_info[0] = gF->alleleCountsHap1[i] + gF->alleleCountsHap2[i];
        // Genotype likelihoods (AA, AB, BB, AC, BC, CC)
        gl_info[0] = gF->genotypeLikelihoods[i][refCharVal*ALPHABET_SIZE+refCharVal];
        gl_info[1] = gF->genotypeLikelihoods[i][refCharVal*ALPHABET_SIZE+h2AlphVal];
        gl_info[2] = gF->genotypeLikelihoods[i][h2AlphVal*ALPHABET_SIZE+h2AlphVal];
        // Update genotype likelihoods
        bcf_update_format(bcf_hdr, bcf_rec, "GL", gl_info, bcf_hdr_nsamples(bcf_hdr)*3, BCF_HT_REAL);

    } else if (h2AlphChar == refChar) {
        // 1|0
        if (*firstVariantInPhaseBlock) {
            gt_info[0] = bcf_gt_unphased(1);
            gt_info[1] = bcf_gt_unphased(0);
            *firstVariantInPhaseBlock = false;
            *phaseSet = bcf_rec->pos+1;
        } else {
            gt_info[0] = bcf_gt_phased(1);
            gt_info[1] = bcf_gt_phased(0);
        }
        kputc(h2AlphChar, str);
        kputc(',', str);
        kputc(h1AlphChar, str);
        // Allele counts - hap2
        ac_info[0] = gF->allele2CountsHap1[i] + gF->allele2CountsHap2[i];
        // Genotype likelihoods (AA, AB, BB, AC, BC, CC)
        gl_info[0] = gF->genotypeLikelihoods[i][refCharVal*ALPHABET_SIZE+refCharVal];
        gl_info[1] = gF->genotypeLikelihoods[i][h1AlphVal*ALPHABET_SIZE+refCharVal];
        gl_info[2] = gF->genotypeLikelihoods[i][h1AlphVal*ALPHABET_SIZE+h1AlphVal];
        // Update genotype likelihoods
        bcf_update_format(bcf_hdr, bcf_rec, "GL", gl_info, bcf_hdr_nsamples(bcf_hdr)*3, BCF_HT_REAL);
    } else {
        // 1|2
        if (*firstVariantInPhaseBlock) {
            gt_info[0] = bcf_gt_unphased(1);
            gt_info[1] = bcf_gt_unphased(2);
            *firstVariantInPhaseBlock = false;
            *phaseSet = bcf_rec->pos+1;
        } else {
            gt_info[0] = bcf_gt_phased(1);
            gt_info[1] = bcf_gt_phased(2);
        }
        kputc(refChar, str);
        kputc(',', str);
        kputc(h1AlphChar, str);
        kputc(',', str);
        kputc(h2AlphChar, str);
        // Allele counts - both hap1 and hap2
        ac_info[0] = gF->alleleCountsHap1[i] + gF->alleleCountsHap2[i];
        ac_info[1] = gF->allele2CountsHap1[i] + gF->allele2CountsHap2[i];
        // Genotype likelihoods (AA, AB, BB, AC, BC, CC)
        gl_info[0] = gF->genotypeLikelihoods[i][refCharVal*ALPHABET_SIZE+refCharVal];
        gl_info[1] = gF->genotypeLikelihoods[i][refCharVal*ALPHABET_SIZE+h1AlphVal];
        gl_info[2] = gF->genotypeLikelihoods[i][h1AlphVal*ALPHABET_SIZE+h1AlphVal];
        gl_info[3] = gF->genotypeLikelihoods[i][refCharVal*ALPHABET_SIZE+h2AlphVal];
        gl_info[4] = gF->genotypeLikelihoods[i][h1AlphVal*ALPHABET_SIZE+h2AlphVal];
        gl_info[5] = gF->genotypeLikelihoods[i][h2AlphVal*ALPHABET_SIZE+h2AlphVal];
        // Update genotype likelihoods
        bcf_update_format(bcf_hdr, bcf_rec, "GL", gl_info, bcf_hdr_nsamples(bcf_hdr)*6, BCF_HT_REAL);
    }

    // Update allele counts
    bcf_update_info(bcf_hdr, bcf_rec, "AC", ac_info, bcf_hdr_nsamples(bcf_hdr), BCF_HT_INT);


}

// This function writes out a vcf for the two haplotypes
// It optionally writes it relative to a reference fasta file or
// writes it for one of the haplotypes relative to the other
void writeVcfFragment(vcfFile *out, bcf_hdr_t *bcf_hdr, stGenomeFragment *gF,
                      char *referenceName, stBaseMapper *baseMapper, bool gvcf) {

    char *referenceSeq;
    // Get reference (needed for VCF generation)
    faidx_t *fai = fai_load(referenceName);
    if ( !fai ) {
        st_logCritical("Could not load fai index of %s.  Maybe you should run 'samtools faidx %s'\n",
                       referenceName, referenceName);
        return;
    }
    int seq_len;
    referenceSeq = fai_fetch(fai, gF->referenceName, &seq_len);
    if ( seq_len < 0 ) {
        st_logCritical("Failed to fetch reference sequence %s\n", referenceName);
        return;
    }
    fai_destroy(fai);

    // intialization
    bcf1_t *bcf_rec = bcf_init1();
    int32_t filter_info = bcf_hdr_id2int(bcf_hdr, BCF_DT_ID, "PASS"); //currently: everything passes
    int32_t *gt_info = (int*)malloc(bcf_hdr_nsamples(bcf_hdr)*2*sizeof(int)); //array specifying phasing
    int32_t *ps_info = (int*)malloc(bcf_hdr_nsamples(bcf_hdr)*sizeof(int)); //array specifying phase sets
    int32_t *dp_info = (int*)malloc(bcf_hdr_nsamples(bcf_hdr)*sizeof(int)); //array specifying read depths
    int32_t *ac_info = (int*)malloc(bcf_hdr_nsamples(bcf_hdr)*2*sizeof(int)); //array specifying allele counts
    float *gl_info = (float*)malloc(bcf_hdr_nsamples(bcf_hdr)*6*sizeof(float)); // array specifying genotype likelihoods
    kstring_t str = {0,0,NULL};
    bool firstVariantInPhaseBlock = true;
    int32_t phaseSet = gF->refStart - 1;

    // iterate over all positions
    for (int64_t i = 0; i < gF->length-1; i++) {

        uint64_t h1AlphVal = gF->haplotypeString1[i];
        uint64_t h2AlphVal = gF->haplotypeString2[i];
        char h1AlphChar = stBaseMapper_getCharForValue(baseMapper, h1AlphVal);
        char h2AlphChar = stBaseMapper_getCharForValue(baseMapper, h2AlphVal);


        uint64_t next_h1AlphVal = gF->haplotypeString1[i + 1];
        uint64_t next_h2AlphVal = gF->haplotypeString2[i + 1];
        char next_h1AlphChar = stBaseMapper_getCharForValue(baseMapper, next_h1AlphVal);
        char next_h2AlphChar = stBaseMapper_getCharForValue(baseMapper, next_h2AlphVal);
//        char nextRefChar = toupper(referenceSeq[i + gF->refStart]); // i + 1 + gF->refStart - 1

        //prep
        bcf_clear1(bcf_rec);
        str.l = 0;

        bcf_rec->rid = bcf_hdr_name2id(bcf_hdr, gF->referenceName); //defined in a contig in the top
        bcf_rec->pos  = i + gF->refStart - 1; // off by one?
        char refChar = toupper(referenceSeq[i + gF->refStart - 1]);
        int refCharVal = stBaseMapper_getValueForChar(baseMapper, refChar);

        // ID - skip
        // QUAL - currently writing out the genotype probability
        float genotypeQuality = -10 * log10f(1 - gF->genotypeProbs[i]);
        // Some programs restrict the maximum genotype quality to be 100.
        if (genotypeQuality > 100) genotypeQuality = 100;
        bcf_rec->qual = (int) genotypeQuality;

        // Get phasing info
        gt_info[0] = bcf_gt_phased(0);
        gt_info[1] = bcf_gt_phased(1);
        dp_info[0] = gF->hap1Depth[i] + gF->hap2Depth[i];
        bcf_update_info(bcf_hdr, bcf_rec, "DP", dp_info, bcf_hdr_nsamples(bcf_hdr), BCF_HT_INT);
        bcf_update_format(bcf_hdr, bcf_rec, "DP", dp_info, bcf_hdr_nsamples(bcf_hdr), BCF_HT_INT);

        if (i + 1 >= gF->length) break;

        if (next_h1AlphChar == '-' || next_h2AlphChar == '-'
            || h1AlphChar == '-' || h2AlphChar == '-') {
            // Insertion or deletion happening here
            writeIndelVariant(gt_info, bcf_hdr, bcf_rec, gF, baseMapper, referenceSeq, refChar,
                              h1AlphChar, h2AlphChar, out, &i, &phaseSet, ps_info, ac_info, gl_info,
                              &firstVariantInPhaseBlock, gvcf);
        }
        else if (h1AlphChar != h2AlphChar) {
            writeHetSite(h1AlphChar, h2AlphChar, refChar,
                         &phaseSet, &firstVariantInPhaseBlock,
                         gt_info, bcf_rec, &str, ac_info, gl_info, gF, baseMapper, bcf_hdr, i);

            // Update genotypes
            bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
            bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);
            // Update phase set
            ps_info[0] = phaseSet;
            bcf_update_format(bcf_hdr, bcf_rec, "PS", ps_info, bcf_hdr_nsamples(bcf_hdr), BCF_HT_INT);

            // Write record
            bcf_write1(out, bcf_hdr, bcf_rec);
        }
        else if ((h1AlphChar != refChar || h2AlphChar != refChar) && h1AlphChar == h2AlphChar) {
            // Doesn't match the reference
            if (firstVariantInPhaseBlock) {
                gt_info[0] = bcf_gt_unphased(1);
                gt_info[1] = bcf_gt_unphased(1);
                firstVariantInPhaseBlock = false;
                phaseSet = bcf_rec->pos+1;
            } else {
                gt_info[0] = bcf_gt_phased(1);
                gt_info[1] = bcf_gt_phased(1);
            }
            // Genotype likelihoods (AA, AB, BB, AC, BC, CC)
            gl_info[0] = gF->genotypeLikelihoods[i][refCharVal*ALPHABET_SIZE+refCharVal];
            gl_info[1] = gF->genotypeLikelihoods[i][refCharVal*ALPHABET_SIZE+h1AlphVal];
            gl_info[2] = gF->genotypeLikelihoods[i][h1AlphVal*ALPHABET_SIZE+h1AlphVal];

            kputc(refChar, &str);
            kputc(',', &str);
            kputc(h2AlphChar, &str);
            bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
            bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);
            // Update genotype likelihoods
            bcf_update_format(bcf_hdr, bcf_rec, "GL", gl_info, bcf_hdr_nsamples(bcf_hdr)*3, BCF_HT_REAL);
            // Update allele counts
            bcf_update_info(bcf_hdr, bcf_rec, "AC", ac_info, bcf_hdr_nsamples(bcf_hdr), BCF_HT_INT);
            ps_info[0] = phaseSet;
            bcf_update_format(bcf_hdr, bcf_rec, "PS", ps_info, bcf_hdr_nsamples(bcf_hdr), BCF_HT_INT);

            // Write record
            bcf_write1(out, bcf_hdr, bcf_rec);
        } else if (gvcf) {
            // Homozygous reference
            kputc(refChar, &str); // REF
            kputc(',', &str);
            kputc(h1AlphChar, &str);
            gt_info[0] = bcf_gt_phased(0);
            gt_info[1] = bcf_gt_phased(0);

            // TODO how to handle genotype likelihoods for this situation in a gvcf?

            bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
            // FORMAT / $SMPL1
            bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);
            ps_info[0] = phaseSet;
            bcf_update_format(bcf_hdr, bcf_rec, "PS", ps_info, bcf_hdr_nsamples(bcf_hdr), BCF_HT_INT);
            bcf_write1(out, bcf_hdr, bcf_rec);
        }
    }
    // Last position
    int h1AlphVal = gF->haplotypeString1[gF->length - 1];
    int h2AlphVal = gF->haplotypeString2[gF->length - 1];
    char h1AlphChar = stBaseMapper_getCharForValue(baseMapper, h1AlphVal);
    char h2AlphChar = stBaseMapper_getCharForValue(baseMapper, h2AlphVal);
    char refChar = toupper(referenceSeq[gF->length - 1 + gF->refStart - 1]);

    // prep
    bcf_clear1(bcf_rec);
    str.l = 0;
    bcf_rec->rid = bcf_hdr_name2id(bcf_hdr, gF->referenceName);
    bcf_rec->pos  = gF->length - 1 + gF->refStart - 1; // off by one?
    gt_info[0] = bcf_gt_phased(0);
    gt_info[1] = bcf_gt_phased(1);
    dp_info[0] = gF->hap1Depth[gF->length - 1] + gF->hap2Depth[gF->length - 1];
    bcf_update_info(bcf_hdr, bcf_rec, "DP", dp_info, bcf_hdr_nsamples(bcf_hdr), BCF_HT_INT);
    bcf_update_format(bcf_hdr, bcf_rec, "DP", dp_info, bcf_hdr_nsamples(bcf_hdr), BCF_HT_INT);
    // TODO: add phasing info, and allele counts, for last record

    if (h1AlphChar != '-' && h2AlphChar != '-') {
        if (gvcf || (h1AlphChar != h2AlphChar || h1AlphChar != refChar || h2AlphChar != refChar)) {
            kputc(refChar, &str); // REF
            kputc(',', &str);
            kputc(h1AlphChar, &str);
            kputc(',', &str);
            kputc(h2AlphChar, &str);

            bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
            bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);
            ps_info[0] = phaseSet;
            bcf_update_format(bcf_hdr, bcf_rec, "PS", ps_info, bcf_hdr_nsamples(bcf_hdr), BCF_HT_INT);
            bcf_write1(out, bcf_hdr, bcf_rec);
        }
    }

    // cleanup
    free(str.s);
    free(gt_info);
    if (referenceSeq) free(referenceSeq);
    bcf_destroy(bcf_rec);
}



stHash *createReferencePriorProbabilities(char *referenceFastaFile, stList *profileSequences,
                                          stBaseMapper *baseMapper, stRPHmmParameters *params) {
    /*
     * Create a set of stReferencePriorProbs that cover the reference intervals included in the profile sequences.
     * The return value is encoded as a map from the reference sequence name (as a string)
     * to the stReferencePriorProbs.
     */

    // Make map from reference sequence names to reference priors
    stHash *referenceNamesToReferencePriors = createEmptyReferencePriorProbabilities(profileSequences);

    // Load reference fasta index
    faidx_t *fai = fai_load(referenceFastaFile);
    if ( !fai ) {
        st_errAbort("Could not load fai index of %s.  Maybe you should run 'samtools faidx %s'\n",
                    referenceFastaFile, referenceFastaFile);
    }

    stHashIterator *hashIt = stHash_getIterator(referenceNamesToReferencePriors);
    char *referenceName;
    while((referenceName = stHash_getNext(hashIt)) != NULL) {
        stReferencePriorProbs *rProbs = stHash_search(referenceNamesToReferencePriors, referenceName);

        // Now get the corresponding reference sequence
        int seqLen;
        char *referenceSeq = fai_fetch(fai, rProbs->referenceName, &seqLen);
        if ( seqLen < 0 ) {
            st_errAbort("Failed to fetch reference sequence %s\n", rProbs->referenceName);
        }

        // Build probability profile
        assert(seqLen >= rProbs->length + rProbs->refStart);
        for(int64_t i=0; i<rProbs->length; i++) {
            uint8_t refChar = stBaseMapper_getValueForChar(baseMapper, referenceSeq[i+rProbs->refStart-1]);
            assert(refChar >= 0 && refChar < ALPHABET_SIZE);
            rProbs->referenceSequence[i] = refChar;
            for(int64_t j=0; j<ALPHABET_SIZE; j++) {
                rProbs->profileProbs[i*ALPHABET_SIZE + j] = *getSubstitutionProb(params->hetSubModel, refChar, j);
            }
        }
    }

    // Cleanup
    fai_destroy(fai);
    stHash_destructIterator(hashIt);

    return referenceNamesToReferencePriors;
}
