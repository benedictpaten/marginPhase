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
        int64_t chunkMarginStartPos = i - chunkMargin;
        chunkMarginStartPos = (chunkMarginStartPos < contigStartPos ? contigStartPos : chunkMarginStartPos);
        int64_t chunkMarginEndPos = chunkEndPos + chunkMargin;
        chunkMarginEndPos = (chunkMarginEndPos > contigEndPos ? contigEndPos : chunkMarginEndPos);

        BamChunk *chunk = bamChunk_construct2(contig, chunkMarginStartPos, i, chunkEndPos,
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
        if ((aln->core.flag & (uint16_t) 0x4) != 0)
            continue; //unaligned
        if (!params->includeSecondaryAlignments && (aln->core.flag & (uint16_t) 0x100) != 0)
            continue; //secondary
        if (!params->includeSupplementaryAlignments && (aln->core.flag & (uint16_t) 0x800) != 0)
            continue; //supplementary

        //data
        char *chr = bamHdr->target_name[aln->core.tid];
        int64_t start_softclip = 0;
        int64_t end_softclip = 0;
        int64_t alnReadLength = getAlignedReadLength3(aln, &start_softclip, &end_softclip, FALSE);
        int64_t alnStartPos = aln->core.pos;
        int64_t alnEndPos = alnStartPos + alnReadLength;

        // does this belong in our chunk?
        if (alnReadLength <= 0) continue;
        if (filterByRegion && (!stString_eq(regionContig, chr) || alnStartPos >= regionEnd || alnEndPos <= regionStart))
            continue;

        // get start and stop position
        int64_t readStartPos = aln->core.pos;           // Left most position of alignment
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

BamChunker *bamChunker_copyConstruct(BamChunker *toCopy) {
    BamChunker *chunker = malloc(sizeof(BamChunker));
    chunker->bamFile = stString_copy(toCopy->bamFile);
    chunker->chunkSize = toCopy->chunkSize;
    chunker->chunkBoundary = toCopy->chunkBoundary;
    chunker->includeSoftClip = toCopy->includeSoftClip;
    chunker->params = toCopy->params;
    chunker->chunks = stList_construct3(0,(void*)bamChunk_destruct);
    chunker->chunkCount = 0;
    return chunker;
}

void bamChunker_destruct(BamChunker *bamChunker) {
    free(bamChunker->bamFile);
    stList_destruct(bamChunker->chunks);
    free(bamChunker);
}

BamChunk *bamChunker_getChunk(BamChunker *bamChunker, int64_t chunkIdx) {
    BamChunk *chunk = stList_get(bamChunker->chunks, chunkIdx);
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

BamChunk *bamChunk_copyConstruct(BamChunk *toCopy) {
    BamChunk *c = malloc(sizeof(BamChunk));
    c->refSeqName = stString_copy(toCopy->refSeqName);
    c->chunkBoundaryStart = toCopy->chunkBoundaryStart;
    c->chunkStart = toCopy->chunkStart;
    c->chunkEnd = toCopy->chunkEnd;
    c->chunkBoundaryEnd = toCopy->chunkBoundaryEnd;
    c->parent = toCopy->parent;
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

uint32_t convertToReadsAndAlignments(BamChunk *bamChunk, RleString *reference,  stList *reads, stList *alignments) {

    // sanity check
    assert(stList_length(reads) == 0);
    assert(stList_length(alignments) == 0);

    uint64_t *ref_nonRleToRleCoordinateMap = reference == NULL ? NULL : rleString_getNonRleToRleCoordinateMap(reference);

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
        if ((aln->core.flag & (uint16_t) 0x4) != 0)
            continue; //unaligned
        if (!bamChunk->parent->params->includeSecondaryAlignments && (aln->core.flag & (uint16_t) 0x100) != 0)
            continue; //secondary
        if (!bamChunk->parent->params->includeSupplementaryAlignments && (aln->core.flag & (uint16_t) 0x800) != 0)
            continue; //supplementary
        if(aln->core.qual < bamChunk->parent->params->filterAlignmentsWithMapQBelowThisThreshold)
        	continue; //low mapping quality

        //data
        char *chr = bamHdr->target_name[aln->core.tid];
        int64_t start_softclip = 0;
        int64_t end_softclip = 0;
        int64_t alnReadLength = getAlignedReadLength3(aln, &start_softclip, &end_softclip, FALSE);
        if (alnReadLength <= 0) continue;
        int64_t alnStartPos = aln->core.pos;
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

        // failure case
        if (stList_length(cigRepr) == 0 || strlen(seq) == 0) {
            stList_destruct(cigRepr);
            free(readName);
            free(seq);
            if (qual != NULL) free(qual);
            continue;
        }

        // sanity check
        assert(stIntTuple_get((stIntTuple *)stList_peek(cigRepr), 1) < strlen(seq));

        // save to read
        bool forwardStrand = !bam_is_rev(aln);
        BamChunkRead *chunkRead = bamChunkRead_construct2(readName, seq, qual, forwardStrand,
                                                          bamChunk->parent->params->useRunLengthEncoding);
        stList_append(reads, chunkRead);

        // save alignment
        if(bamChunk->parent->params->useRunLengthEncoding) {
            // ref_nonRleToRleCoordinateMap should only be null w/ RLE in tests
            if (ref_nonRleToRleCoordinateMap != NULL) {
                // rle the alignment and save it
                uint64_t *read_nonRleToRleCoordinateMap = rleString_getNonRleToRleCoordinateMap(chunkRead->rleRead);
                stList_append(alignments, runLengthEncodeAlignment(cigRepr, ref_nonRleToRleCoordinateMap, read_nonRleToRleCoordinateMap));
                stList_destruct(cigRepr);
                free(read_nonRleToRleCoordinateMap);
            }
        }
        else {
            stList_append(alignments, cigRepr);
        }
        savedAlignments++;

        // cleanup
        free(readName);
        free(seq);
        if (qual != NULL) free(qual);
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
    if (ref_nonRleToRleCoordinateMap != NULL)
        free(ref_nonRleToRleCoordinateMap);
    return savedAlignments;
}

bool poorMansDownsample(int64_t intendedDepth, BamChunk *bamChunk, stList *reads, stList *alignments,
        stList *filteredReads, stList *filteredAlignments, stList *discardedReads, stList *discardedAlignments) {

    // calculate depth
    int64_t totalNucleotides = 0;
    for (int64_t i = 0; i < stList_length(reads); i++) {
        BamChunkRead *bcr = stList_get(reads,i);
        totalNucleotides += bcr->rleRead->length;
    }
    double averageDepth = 1.0 * totalNucleotides / (bamChunk->chunkBoundaryEnd - bamChunk->chunkBoundaryStart);

    // do we need to downsample?
    if (averageDepth < intendedDepth) {
        return FALSE;
    }

    // we do need to downsample
    char *logIdentifier = getLogIdentifier();
    st_logInfo(" %s Downsampling chunk with average depth %.2fx to %dx \n", logIdentifier, averageDepth, intendedDepth);
    free(logIdentifier);

    // keep some ratio of reads
    double ratioToKeep = intendedDepth / averageDepth;
    for (int64_t i = 0; i < stList_length(reads); i++) {
        if (st_random() < ratioToKeep) {
            stList_append(filteredReads, stList_get(reads, i));
            stList_append(filteredAlignments, stList_get(alignments, i));
        } else {
            stList_append(discardedReads, stList_get(reads, i));
            stList_append(discardedAlignments, stList_get(alignments, i));
        }
    }

    return TRUE;
}
