//
// Created by tpesout on 1/8/19.
//

#ifndef MARGINPHASE_EXTERNAL_INTEGRATION_H

#include <htslib/vcf.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/bgzf.h>
#include <htslib/hts.h>
#include "bedidx.h"
#include "margin.h"


BamChunker *bamChunker_construct(char *bamFile, PolishParams *params);
BamChunker *bamChunker_construct2(char *bamFile, char *region, PolishParams *params);
BamChunker *bamChunker_copyConstruct(BamChunker *toCopy);
void bamChunker_destruct(BamChunker *bamChunker);
BamChunk *bamChunker_getChunk(BamChunker *bamChunker, int64_t chunkIdx);

BamChunk *bamChunk_construct();
BamChunk *bamChunk_construct2(char *refSeqName, int64_t chunkBoundaryStart, int64_t chunkStart, int64_t chunkEnd,
                              int64_t chunkBoundaryEnd, BamChunker *parent);
BamChunk *bamChunk_copyConstruct(BamChunk *toCopy);
void bamChunk_destruct(BamChunk *bamChunk);

uint32_t convertToReadsAndAlignments(BamChunk *bamChunk, stList *reads, stList *alignments);
bool poorMansDownsample(int64_t intendedDepth, BamChunk *bamChunk, stList *reads, stList *alignments,
                        stList *filteredReads, stList *filteredAlignments, stList *discardedReads, stList *discardedAlignments);

int64_t getAlignedReadLength(bam1_t *aln);
int64_t getAlignedReadLength2(bam1_t *aln, int64_t *start_softclip, int64_t *end_softclip);
int64_t getAlignedReadLength3(bam1_t *aln, int64_t *start_softclip, int64_t *end_softclip, bool boundaryAtMatch);


int64_t parseReads(stList *profileSequences, char *bamFile, stBaseMapper *baseMapper, stRPHmmParameters *params);

int64_t parseReadsWithSingleNucleotideProbs(stList *profileSequences, char *bamFile, stBaseMapper *baseMapper,
                                            stRPHmmParameters *params, char *signalAlignDirectory, bool onlySignalAlign);

void countIndels(uint32_t *cigar, uint32_t ncigar, int64_t *numInsertions, int64_t *numDeletions);


// Output file writing methods

void writeVcfFragment(vcfFile *out, bcf_hdr_t *bcf_hdr, stGenomeFragment *gF, char *referenceName,
                      stBaseMapper *baseMapper, bool gvcf);

bcf_hdr_t* writeVcfHeader(vcfFile *out, stList *genomeFragments, char *referenceName);


void writeHaplotypedSam(char *bamInFile, char *bamOutBase, stReadHaplotypePartitionTable *readHaplotypePartitions,
                        char *marginPhaseTag);

void writeSplitSams(char *bamInFile, char *bamOutBase, stReadHaplotypePartitionTable *readHaplotypePartitions,
                    char *marginPhaseTag);

stHash *createReferencePriorProbabilities(char *referenceFastaFile, stList *profileSequences,
                                          stBaseMapper *baseMapper, stRPHmmParameters *params);

#define MARGINPHASE_EXTERNAL_INTEGRATION_H
#endif //MARGINPHASE_EXTERNAL_INTEGRATION_H
