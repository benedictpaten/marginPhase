//
// Created by tpesout on 1/8/19.
//

#ifndef MARGINPHASE_HTS_INTEGRATION_H

#include <htslib/vcf.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/bgzf.h>
#include <htslib/hts.h>
#include "bedidx.h"
#include "margin.h"


BamChunker *bamChunker_construct(char *bamFile, PolishParams *params);
BamChunker *bamChunker_construct2(char *bamFile, char *region, PolishParams *params);
void bamChunker_destruct(BamChunker *bamChunker);
BamChunk *bamChunker_getNext(BamChunker *bamChunker);

BamChunk *bamChunk_construct();
BamChunk *bamChunk_construct2(char *refSeqName, int64_t chunkBoundaryStart, int64_t chunkStart, int64_t chunkEnd,
                              int64_t chunkBoundaryEnd, BamChunker *parent);
void bamChunk_destruct(BamChunk *bamChunk);

uint32_t convertToReadsAndAlignments(BamChunk *bamChunk, stList *reads, stList *alignments);


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



stReferencePriorProbs *stReferencePriorProbs_constructEmptyProfile(char *referenceName, int64_t referenceStart, int64_t length);

void stReferencePriorProbs_destruct(stReferencePriorProbs *seq);

stHash *createEmptyReferencePriorProbabilities(stList *profileSequences);

stHash *createReferencePriorProbabilities(char *referenceFastaFile, stList *profileSequences,
                                          stBaseMapper *baseMapper, stRPHmmParameters *params);

int64_t filterHomozygousReferencePositions(stHash *referenceNamesToReferencePriors, stRPHmmParameters *params, int64_t *totalPositions);

double *stReferencePriorProbs_estimateReadErrorProbs(stHash *referenceNamesToReferencePriors, stRPHmmParameters *params);


#define MARGINPHASE_HTS_INTEGRATION_H
#endif //MARGINPHASE_HTS_INTEGRATION_H
