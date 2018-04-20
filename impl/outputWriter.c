/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com) & Arthur Rand (arand@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <ctype.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/hts.h>
#include <hashTableC.h>
#include "stRPHmm.h"



void writeHaplotypedBam(char *bamInFile, char *bamOutBase, stReadHaplotypePartitionTable *readHaplotypePartitions) {
    // prep
    char *haplotypedBamFile = stString_print("%s.haplotyped.bam", bamOutBase);

    // file management
    samFile *in = hts_open(bamInFile, "r");
    if (in == NULL) {
        st_errAbort("ERROR: Cannot open bam file %s\n", bamInFile);
    }
    bam_hdr_t *bamHdr = bam_hdr_read(in->fp.bgzf);
    bam1_t *aln = bam_init1();

    int r;
    st_logDebug("\tWriting haplotype output to: %s \n", haplotypedBamFile);
    BGZF *out = bgzf_open(haplotypedBamFile, "w");
    r = bam_hdr_write(out, bamHdr);

    // read in input file, write out each read to one sam file
    int32_t readCountH1 = 0;
    int32_t readCountH2 = 0;
    int32_t readCountFiltered = 0;
    char *haplotypeString;
    while(sam_read1(in,bamHdr,aln) > 0) {

        char *readName = bam_get_qname(aln);
        stReadHaplotypeSequence *readHaplotypes = hashtable_search(readHaplotypePartitions, readName);
        if (readHaplotypes == NULL) {
            haplotypeString = stReadHaplotypeSequence_toStringEmpty();
            bam_aux_append(aln, "MP", 'Z', (int)strlen(haplotypeString) + 1, (uint8_t*)haplotypeString);
            r = bam_write1(out, aln);
            readCountFiltered++;
        } else {
            haplotypeString = stReadHaplotypeSequence_toString(readHaplotypes);
            bam_aux_append(aln, "MP", 'Z', (int)strlen(haplotypeString) + 1, (uint8_t*)haplotypeString);
            r = bam_write1(out, aln);
            // document based on last recorded haplotype
            while (readHaplotypes->next != NULL) {readHaplotypes = readHaplotypes->next;}
            if (readHaplotypes->haplotype == 1)
                readCountH1++;
            else
                readCountH2++;
        }
        free(haplotypeString);
    }
    st_logInfo("\tBAM read counts:\n\t\thap1: %d\thap2: %d\tfiltered out: %d \n", readCountH1, readCountH2, readCountFiltered);

    // Cleanup
    bam_destroy1(aln);
    bam_hdr_destroy(bamHdr);
    sam_close(in);
    bgzf_close(out);
    free(haplotypedBamFile);
}

void writeSplitSams(char *bamInFile, char *bamOutBase,
                    stReadHaplotypePartitionTable *readHaplotypePartitions) {
    // prep
    char *haplotype1SamOutFile = stString_print("%s.1.sam", bamOutBase);
    char *haplotype2SamOutFile = stString_print("%s.2.sam", bamOutBase);
    char *unmatchedSamOutFile = stString_print("%s.filtered.sam", bamOutBase);

    // file management
    samFile *in = hts_open(bamInFile, "r");
    if (in == NULL) {
        st_errAbort("ERROR: Cannot open bam file %s\n", bamInFile);
    }
    bam_hdr_t *bamHdr = sam_hdr_read(in);
    bam1_t *aln = bam_init1();

    int r;
    st_logDebug("\tWriting haplotype output to: %s and %s \n", haplotype1SamOutFile, haplotype2SamOutFile);
    samFile *out1 = hts_open(haplotype1SamOutFile, "w");
    r = sam_hdr_write(out1, bamHdr);

    samFile *out2 = hts_open(haplotype2SamOutFile, "w");
    r = sam_hdr_write(out2, bamHdr);

    samFile *outUnmatched = hts_open(unmatchedSamOutFile, "w");
    r = sam_hdr_write(outUnmatched, bamHdr);


    // read in input file, write out each read to one sam file
    int32_t readCountH1 = 0;
    int32_t readCountH2 = 0;
    int32_t readCountFiltered = 0;
    char *haplotypeString;
    while(sam_read1(in,bamHdr,aln) > 0) {

        char *readName = bam_get_qname(aln);
        stReadHaplotypeSequence *readHaplotypes = hashtable_search(readHaplotypePartitions, readName);
        if (readHaplotypes == NULL) {
            haplotypeString = stReadHaplotypeSequence_toStringEmpty();
            bam_aux_append(aln, "MP", 'Z', (int)strlen(haplotypeString) + 1, (uint8_t*)haplotypeString);
            r = sam_write1(outUnmatched, bamHdr, aln);
            readCountFiltered++;
        } else {
            haplotypeString = stReadHaplotypeSequence_toString(readHaplotypes);
            bam_aux_append(aln, "MP", 'Z', (int)strlen(haplotypeString) + 1, (uint8_t*)haplotypeString);
            // document based on last recorded haplotype
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
    bcf_hdr_t *hdr = bcf_hdr_init("w");
    kstring_t str = {0,0,NULL};

    // generic info
    str.l = 0;
    ksprintf(&str, "##marginPhase=htslib-%s\n",hts_version());
    bcf_hdr_append(hdr, str.s);

    // reference file used
    str.l = 0;
    ksprintf(&str, "##reference=file://%s\n", referenceName);
    bcf_hdr_append(hdr, str.s);

    // contigs
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
    ksprintf(&str, "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set for GT\">");
    bcf_hdr_append(hdr, str.s);
    str.l = 0;
    ksprintf(&str, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">");
    bcf_hdr_append(hdr, str.s);


    // samples
    bcf_hdr_add_sample(hdr, "SMPL1"); //todo
    bcf_hdr_add_sample(hdr, NULL);

    // write header
    bcf_hdr_write(out, hdr);

    // cleanup
    free(str.s);
    return hdr;
}


void writeIndelVariant(int32_t *gt_info, bcf_hdr_t *bcf_hdr, bcf1_t *bcf_rec, stGenomeFragment *gF,
                       stBaseMapper *baseMapper, char *referenceSeq, int refChar,
                       int h1AlphChar, int h2AlphChar, vcfFile *out, int64_t *index,
                       int32_t  *phaseSet, int32_t  *ps_info, int32_t *ac_info,
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
                  int32_t *ac_info, stGenomeFragment *gF, int64_t i) {
    /*
     * Write out a het site record.
     */
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
    }

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
    kstring_t str = {0,0,NULL};
    bool firstVariantInPhaseBlock = true;
    int32_t phaseSet = gF->refStart - 1;

    // iterate over all positions
    for (int64_t i = 0; i < gF->length-1; i++) {

        int h1AlphVal = gF->haplotypeString1[i];
        int h2AlphVal = gF->haplotypeString2[i];
        char h1AlphChar = stBaseMapper_getCharForValue(baseMapper, h1AlphVal);
        char h2AlphChar = stBaseMapper_getCharForValue(baseMapper, h2AlphVal);


        int next_h1AlphVal = gF->haplotypeString1[i + 1];
        int next_h2AlphVal = gF->haplotypeString2[i + 1];
        char next_h1AlphChar = stBaseMapper_getCharForValue(baseMapper, next_h1AlphVal);
        char next_h2AlphChar = stBaseMapper_getCharForValue(baseMapper, next_h2AlphVal);
        char nextRefChar = toupper(referenceSeq[i + gF->refStart]); // i + 1 + gF->refStart - 1

        //prep
        bcf_clear1(bcf_rec);
        str.l = 0;

        bcf_rec->rid = bcf_hdr_name2id(bcf_hdr, gF->referenceName); //defined in a contig in the top
        bcf_rec->pos  = i + gF->refStart - 1; // off by one?

        // ID - skip
        // QUAL - currently writing out the genotype probability
        float genotypeQuality = -10 * log10f(1 - gF->genotypeProbs[i]);
        if (genotypeQuality > 100) genotypeQuality = 100;
        bcf_rec->qual = (int) genotypeQuality;

        // Get phasing info
        gt_info[0] = bcf_gt_phased(0);
        gt_info[1] = bcf_gt_phased(1);
        dp_info[0] = gF->hap1Depth[i] + gF->hap2Depth[i];
        bcf_update_info(bcf_hdr, bcf_rec, "DP", dp_info, bcf_hdr_nsamples(bcf_hdr), BCF_HT_INT);
        bcf_update_format(bcf_hdr, bcf_rec, "DP", dp_info, bcf_hdr_nsamples(bcf_hdr), BCF_HT_INT);

        char refChar = toupper(referenceSeq[i + gF->refStart - 1]);
        if (gvcf) {
            if (next_h1AlphChar == '-' || next_h2AlphChar == '-'
                    || h1AlphChar == '-' || h2AlphChar == '-') {
                // Insertion or deletion happening here
                writeIndelVariant(gt_info, bcf_hdr, bcf_rec, gF, baseMapper, referenceSeq, refChar,
                                  h1AlphChar, h2AlphChar, out, &i, &phaseSet, ps_info, ac_info,
                                  &firstVariantInPhaseBlock, gvcf);
            }
            else if (h1AlphChar != h2AlphChar) {
                writeHetSite(h1AlphChar, h2AlphChar, refChar,
                &phaseSet, &firstVariantInPhaseBlock,
                gt_info, bcf_rec, &str, ac_info, gF, i);
                // Update genotypes
                bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
                bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);
                ps_info[0] = phaseSet;
                bcf_update_format(bcf_hdr, bcf_rec, "PS", ps_info, bcf_hdr_nsamples(bcf_hdr), BCF_HT_INT);
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
                kputc(refChar, &str);
                kputc(',', &str);
                kputc(h2AlphChar, &str);
                bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
                bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);
                ps_info[0] = phaseSet;
                bcf_update_format(bcf_hdr, bcf_rec, "PS", ps_info, bcf_hdr_nsamples(bcf_hdr), BCF_HT_INT);
                bcf_write1(out, bcf_hdr, bcf_rec);
            }
            // TODO only difference is that this last 'else' exists to write all the other records...
            else {
                // Homozygous reference
                kputc(refChar, &str); // REF
                kputc(',', &str);
                kputc(h1AlphChar, &str);
                gt_info[0] = bcf_gt_phased(0);
                gt_info[1] = bcf_gt_phased(0);

                bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
                // FORMAT / $SMPL1
                bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);
                ps_info[0] = phaseSet;
                bcf_update_format(bcf_hdr, bcf_rec, "PS", ps_info, bcf_hdr_nsamples(bcf_hdr), BCF_HT_INT);
                bcf_write1(out, bcf_hdr, bcf_rec);
            }
        }
        // Regular vcf
        else {
            if (i + 1 >= gF->length) break;

            if (next_h1AlphChar == '-' || next_h2AlphChar == '-'
                    || h1AlphChar == '-' || h2AlphChar == '-') {
                // Insertion or deletion happening here
                writeIndelVariant(gt_info, bcf_hdr, bcf_rec, gF, baseMapper, referenceSeq, refChar,
                                  h1AlphChar, h2AlphChar, out, &i, &phaseSet, ps_info, ac_info,
                                  &firstVariantInPhaseBlock, gvcf);
            }
            else if (h1AlphChar != h2AlphChar) {
                writeHetSite(h1AlphChar, h2AlphChar, refChar,
                             &phaseSet, &firstVariantInPhaseBlock,
                             gt_info, bcf_rec, &str, ac_info, gF, i);
                // Update genotypes
                bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
                bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);
                // Update phase set
                ps_info[0] = phaseSet;
                bcf_update_format(bcf_hdr, bcf_rec, "PS", ps_info, bcf_hdr_nsamples(bcf_hdr), BCF_HT_INT);
                // Update allele counts
                bcf_update_info(bcf_hdr, bcf_rec, "AC", ac_info, bcf_hdr_nsamples(bcf_hdr), BCF_HT_INT);
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
                kputc(refChar, &str);
                kputc(',', &str);
                kputc(h2AlphChar, &str);
                // Update genotypes
                bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
                bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);
                // Update phase set
                ps_info[0] = phaseSet;
                bcf_update_format(bcf_hdr, bcf_rec, "PS", ps_info, bcf_hdr_nsamples(bcf_hdr), BCF_HT_INT);
                // Update allele counts
                ac_info[0] = gF->alleleCountsHap1[i] + gF->alleleCountsHap2[i];
                bcf_update_info(bcf_hdr, bcf_rec, "AC", ac_info, bcf_hdr_nsamples(bcf_hdr), BCF_HT_INT);
                // Write record
                bcf_write1(out, bcf_hdr, bcf_rec);
            }
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




/*
 * Prints information contained in genotypeResults struct.
 */
void printGenotypeResults(stGenotypeResults *results) {
    // Variables
    int64_t hetsInRefIndels = results->hetsInRef_Insertions + results->hetsInRef_Deletions;
    int64_t homozygousVariantsInRefIndels = results->homozygousVariantsInRef_Insertions + results->homozygousVariantsInRef_Deletions;

    st_logInfo("\n----- RESULTS -----\n");

    // Sensitivity
    st_logInfo("\nSensitivity / recall (fraction of true positives compared to reference): \n");
    st_logInfo("\tOverall: %.4f \t\t\t(%d out of %d)\n",
               (float)results->truePositives/results->positives,
               results->truePositives,
               results->positives);
    st_logInfo("\tSNVs: %.4f \t\t\t\t(%d out of %d)\n",
               (float) (results->truePositives-results->truePositiveIndels-results->truePositiveHomozygousIndels)/
                       (results->positives-hetsInRefIndels-homozygousVariantsInRefIndels),
               results->truePositives-results->truePositiveIndels-results->truePositiveHomozygousIndels,
               results->positives-hetsInRefIndels-homozygousVariantsInRefIndels);
    st_logInfo("\t\tHets: %.4f \t\t\t(%d out of %d)\n",
               (float) (results->truePositiveHet - results->truePositiveHetIndels)
               / (results->hetsInRef-hetsInRefIndels),
               results->truePositiveHet - results->truePositiveHetIndels,
               results->hetsInRef-hetsInRefIndels);
    st_logInfo("\t\tHomozygous alts: %.4f \t(%d out of %d)\n",
               (float) (results->truePositiveHomozygous - results->truePositiveHomozygousIndels)
               /(results->homozygousVariantsInRef-homozygousVariantsInRefIndels),
               results->truePositiveHomozygous - results->truePositiveHomozygousIndels,
               results->homozygousVariantsInRef - homozygousVariantsInRefIndels);
    st_logInfo("\tIndels: %.4f \t\t\t\t(%d out of %d)\n",
               (float)results->truePositiveIndels/
                       (homozygousVariantsInRefIndels + hetsInRefIndels),
               results->truePositiveIndels,
               hetsInRefIndels + homozygousVariantsInRefIndels);
    st_logInfo("\t\tHets: %.4f \t\t\t(%d out of %d)\n",
               (float) results->truePositiveHetIndels / hetsInRefIndels,
               results->truePositiveHetIndels,
               hetsInRefIndels);
    st_logInfo("\t\tHomozygous alts: %.4f \t(%d out of %d)\n",
               (float) results->truePositiveHomozygousIndels/homozygousVariantsInRefIndels,
               results->truePositiveHomozygousIndels,
               homozygousVariantsInRefIndels);

    // Specificity
    st_logInfo("\nSpecificity (fraction of true negatives compared to reference): \n");
    st_logInfo("\tOverall: %.4f \t\t(%" PRIi64 " out of %"PRIi64 ")\n",
               (float)results->trueNegatives/results->negatives,
               results->trueNegatives, results->negatives);

    // Precision
    st_logInfo("\nPrecision (fraction of true positives compared to all calls):\n");
    st_logInfo("\tOverall: %.4f \t(%d out of %d)\n",
               (float)results->truePositives/(results->truePositives+results->falsePositives),
               results->truePositives, results->truePositives+results->falsePositives);
    st_logInfo("\tSNVs: %.4f \t\t(%d out of %d)\n",
               (float) (results->truePositives-results->truePositiveIndels) /
                       (results->truePositives + results->falsePositives - results->truePositiveIndels - results->falsePositiveIndels),
               (results->truePositives - results->truePositiveIndels),
               (results->truePositives + results->falsePositives - results->truePositiveIndels - results->falsePositiveIndels));
    st_logInfo("\tIndels: %.4f \t\t(%d out of %d)\n",
               (float)results->truePositiveIndels / (results->truePositiveIndels + results->falsePositiveIndels),
               results->truePositiveIndels, results->truePositiveIndels + results->falsePositiveIndels);

    // False positives
    st_logInfo("\nFalse positives:\n");
    st_logInfo("\tOverall: %" PRIi64 "\n", results->falsePositives);
    st_logInfo("\tSNVs: %" PRIi64 "\n", results->falsePositives - results->falsePositiveIndels);
    st_logInfo("\tIndels: %" PRIi64 "\n", results->falsePositiveIndels);


    // More detailed numbers about errors
    st_logInfo("\nFalse negatives:\n");
    st_logInfo("\tOverall: %" PRIi64 "\n", results->falseNegatives);
    st_logInfo("\tSNV - Missed hets: %" PRIi64 "\n", results->error_missedHet -
                       results->error_missedHet_Deletions - results->error_missedHet_Insertions);
    st_logInfo("\tSNV - Incorrect homozygous variants: %" PRIi64 "\n",
               results->error_homozygousInRef -
               results->error_homozygous_Insertions - results->error_homozygous_Deletions);

    st_logInfo("\tIndels - Missed hets: %" PRIi64 "\n",
               results->error_missedHet_Insertions + results->error_missedHet_Deletions);
    st_logInfo("\t\tInsertions: %d \t(found %d out of %d, %.3f)\n",
               results->error_missedHet_Insertions,
               results->hetsInRef_Insertions - results->error_missedHet_Insertions,
               results->hetsInRef_Insertions,
               (float)(results->hetsInRef_Insertions - results->error_missedHet_Insertions) /
                       results->hetsInRef_Insertions);
    st_logInfo("\t\tDeletions: %d \t(found %d out of %d, %.3f)\n",
               results->error_missedHet_Deletions,
               results->hetsInRef_Deletions - results->error_missedHet_Deletions,
               results->hetsInRef_Deletions,
               (float)(results->hetsInRef_Deletions - results->error_missedHet_Deletions) /
                       results->hetsInRef_Deletions);

    st_logInfo("\tIndels - Incorrect homozygous variants: %" PRIi64 "\n",
               results->error_homozygous_Insertions + results->error_homozygous_Deletions);
    st_logInfo("\t\tInsertions: %d \t(found %d out of %d, %.3f)\n",
               results->error_homozygous_Insertions,
               results->homozygousVariantsInRef_Insertions - results->error_homozygous_Insertions,
               results->homozygousVariantsInRef_Insertions,
               (float)(results->homozygousVariantsInRef_Insertions -
                       results->error_homozygous_Insertions) /
                       results->homozygousVariantsInRef_Insertions);
    st_logInfo("\t\tDeletions: %d \t(found %d out of %d, %.3f)\n",
               results->error_homozygous_Deletions,
               results->homozygousVariantsInRef_Deletions - results->error_homozygous_Deletions,
               results->homozygousVariantsInRef_Deletions,
               (float)(results->homozygousVariantsInRef_Deletions -
                       results->error_homozygous_Deletions) /
                       results->homozygousVariantsInRef_Deletions);

    // Phasing
    st_logInfo("\nPhasing:\n");
    st_logInfo("\tSwitch error rate: %.4f \t (%" PRIi64 " out of %"PRIi64 ", ", (float)results->switchErrors/(results->truePositives-results->uncertainPhasing), results->switchErrors, results->truePositives-results->uncertainPhasing);
    st_logInfo("fraction correct: %.4f)\n", 1.0 - (float)results->switchErrors/(results->truePositives-results->uncertainPhasing));
    st_logInfo("\tAverage distance between switch errors: %.2f\n", results->switchErrorDistance);
    st_logInfo("\tUncertain phasing spots: %" PRIi64 "\n\n", results->uncertainPhasing);
}


void writeParamFile(char *outputFilename, stRPHmmParameters *params) {
    // get file
    FILE *fd = fopen(outputFilename, "w");
    if (fd == NULL) {
        st_logCritical("Failed to open output param file '%s'. No file will be written\n", outputFilename);
        return;
    }
    //for whether to print the last comma
    int64_t noCommaIdx = (ALPHABET_SIZE) * (ALPHABET_SIZE) - 1;

    fprintf(fd, "{\n");
    fprintf(fd, "  \"alphabet\" : [ \"Aa\", \"Cc\", \"Gg\", \"Tt\", \"-\" ],\n");
    fprintf(fd, "    \n");
    fprintf(fd, "  \"wildcard\" : \"Nn\",\n");
    fprintf(fd, "    \n");
    fprintf(fd, "  \"haplotypeSubstitutionModel\" : [ \n");
    for (int64_t i = 0; i < ALPHABET_SIZE; i++) {
        fprintf(fd, "    ");
        for (int64_t j = 0; j < ALPHABET_SIZE; j++) {
            int64_t idx = i * ALPHABET_SIZE + j;
            fprintf(fd, "%8f", exp(params->hetSubModelSlow[idx]));
            if (idx != noCommaIdx) fprintf(fd, ", ");
        }
        fprintf(fd, "\n");
    }
    fprintf(fd, "   ],\n");
    fprintf(fd, "    \n");
    fprintf(fd, "  \"readErrorModel\" : [ \n");
    for (int64_t i = 0; i < ALPHABET_SIZE; i++) {
        fprintf(fd, "    ");
        for (int64_t j = 0; j < ALPHABET_SIZE; j++) {
            int64_t idx = i * ALPHABET_SIZE + j;
            fprintf(fd, "%8f", exp(params->readErrorSubModelSlow[idx]));
            if (idx != noCommaIdx) fprintf(fd, ", ");
        }
        fprintf(fd, "\n");
    }
    fprintf(fd, "   ],\n");
    fprintf(fd, "    \n");
    fprintf(fd, "  \"maxNotSumTransitions\" : %s,\n", params->maxNotSumTransitions ? "true" : "false");
    fprintf(fd, "    \n");
    fprintf(fd, "  \"maxPartitionsInAColumn\" : %" PRIi64 ",\n", params->maxPartitionsInAColumn);
    fprintf(fd, "\n");
    fprintf(fd, "  \"maxCoverageDepth\" : %" PRIi64 ",\n", params->maxCoverageDepth);
    fprintf(fd, "\n");
    fprintf(fd, "  \"minReadCoverageToSupportPhasingBetweenHeterozygousSites\" : %" PRIi64 ",\n", params->minReadCoverageToSupportPhasingBetweenHeterozygousSites);
    fprintf(fd, "  \n");
    fprintf(fd, "  \"onDiagonalReadErrorPseudoCount\" : %f,\n", params->onDiagonalReadErrorPseudoCount);
    fprintf(fd, "  \n");
    fprintf(fd, "  \"offDiagonalReadErrorPseudoCount\" : %f,\n", params->offDiagonalReadErrorPseudoCount);
    fprintf(fd, "  \n");
    fprintf(fd, "  \"trainingIterations\" : %" PRIi64 ",\n", params->trainingIterations);
    fprintf(fd, "  \n");
    int64_t verbosityBitstring = 0
        | (params->verboseTruePositives ? LOG_TRUE_POSITIVES : 0);
    fprintf(fd, "  \"verbose\" : %" PRIi64 ",\n", verbosityBitstring);
    fprintf(fd, "}");

    if (fclose(fd) != 0) st_logCritical("Failed to close output param file: %s\n", outputFilename);
}

//
//  Read Haplotype Sequence functions
//

stReadHaplotypeSequence *stReadHaplotypeSequence_construct(
        int64_t readStart, int64_t phaseBlock, int64_t length, int8_t haplotype) {
    stReadHaplotypeSequence *s = malloc(sizeof(stReadHaplotypeSequence));
    s->readStart = readStart;
    s->phaseBlock = phaseBlock;
    s->length = length;
    s->haplotype = haplotype;
    s->next = NULL;
    return s;
}

char *stReadHaplotypeSequence_toString(stReadHaplotypeSequence *rhs) {
    char *curr = stString_print("h%"PRIi8",p%"PRIi64",r%"PRIi64",l%"PRIi64, rhs->haplotype, rhs->phaseBlock,
                                   rhs->readStart, rhs->length);
    if (rhs->next != NULL) {
        char *next = stReadHaplotypeSequence_toString(rhs->next);
        char *tmp = curr;
        curr = stString_print("%s;%s", tmp, next);
        free(tmp);
        free(next);
    }

    return curr;
}

char *stReadHaplotypeSequence_toStringEmpty() {
    return stString_print("h%"PRIi8",p%"PRIi64",r%"PRIi64",l%"PRIi64, 0, -1, -1, -1);
}

void stReadHaplotypeSequence_destruct(stReadHaplotypeSequence *rhs) {
    if (rhs->next != NULL) {
        stReadHaplotypeSequence_destruct(rhs->next);
    }
    free(rhs);
}


//
//  HaplotypePartitionTable
//  Tracks ReadHaplotypeSequence information
//

//todo destroy keys (char*)
stReadHaplotypePartitionTable *stReadHaplotypePartitionTable_construct(int64_t initialSize) {
    return create_hashtable((uint64_t) initialSize, stHash_stringKey, stHash_stringEqualKey,
                            NULL, (void *)stReadHaplotypeSequence_destruct);
}
void stReadHaplotypePartitionTable_add(stReadHaplotypePartitionTable *hpt, char *readName, int64_t readStart,
                                       int64_t phaseBlock, int64_t length, int8_t haplotype) {

    stReadHaplotypeSequence *new = stReadHaplotypeSequence_construct(readStart, phaseBlock, length, haplotype);
    stReadHaplotypeSequence *curr = hashtable_search(hpt, readName);
    if (NULL == curr) {
        hashtable_insert(hpt, readName, new);
    } else {
        stReadHaplotypeSequence *prev;
        do {
            prev = curr;
            if (curr->phaseBlock == phaseBlock) {
                if (haplotype != curr->haplotype) {
                    st_logCritical("\tRead %s found in both haplotypes in phase block %"PRIi64"\n", readName, phaseBlock);
                }

                // don't need to record
                stReadHaplotypeSequence_destruct(new);
                return;
            }
            curr = curr->next;
        } while (curr != NULL);
        prev->next = new;
    }
}
void stReadHaplotypePartitionTable_destruct(stReadHaplotypePartitionTable *hpt) {
    hashtable_destroy(hpt, true, false);
}


//
//  Populating HaplotypePartitionTable from GenomeFragment and HMM
//

void populateReadHaplotypePartitionTable(stReadHaplotypePartitionTable *hpt, stGenomeFragment *gF, stRPHmm *hmm,
                                         stList *path) {
    //todo track all partitioned reads and quit early if examined
    // same for whole GenomeFragment
    int64_t phaseBlock = gF->refStart - 1;
    int64_t phaseBlockEnd = phaseBlock + gF->length;

    // variables for partitions
    char *readName;
    int64_t readStart;
    int64_t length;
    int8_t haplotype;

    // For each cell/column pair
    stRPColumn *column = hmm->firstColumn;
    for(int64_t i=0; i<stList_length(path); i++) {
        stRPCell *cell = stList_get(path, i);

        // Get reads in partitions
        for(int64_t j=0; j<column->depth; j++) {
            stProfileSeq *read = column->seqHeaders[j];
            readName = read->readId;
            readStart = (read->refStart < phaseBlock ? phaseBlock - read->refStart : 0);
            length = (read->refStart + read->length > phaseBlockEnd ?
                      phaseBlockEnd - read->refStart :
                      read->length - readStart);
            haplotype = (int8_t) (seqInHap1(cell->partition, j) ? 1 : 2);

            //todo more sanity check
            assert(length >= 0);
            assert(length <= read->length);
            assert(readStart >= 0);
            assert(readStart <= read->length);

            // save to hpt
            stReadHaplotypePartitionTable_add(hpt, readName, readStart, phaseBlock, length, haplotype);
        }

        // iterate
        if(column->nColumn != NULL) {
            column = column->nColumn->nColumn;
        }
    }

}