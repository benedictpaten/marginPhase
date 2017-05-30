/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com) & Arthur Rand (arand@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <ctype.h>
#include <htslib/vcf.h>
#include "stRPHmm.h"

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

    // formatting
    str.l = 0;
    ksprintf(&str, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
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

// This function writes out a vcf for the two haplotypes
// It optionally writes it relative to a reference fasta file or
// writes it for one of the haplotypes relative to the other
void writeVcfFragment(vcfFile *out, bcf_hdr_t *bcf_hdr, stGenomeFragment *gF, char *referenceName, stBaseMapper *baseMapper, bool differencesOnly) {

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
    kstring_t str = {0,0,NULL};

    int totalLocs = 0;

    // iterate over all positions
    for (int64_t i = 0; i < gF->length; i++) {

        int h1AlphVal = gF->haplotypeString1[i];
        int h2AlphVal = gF->haplotypeString2[i];
        char h1AlphChar = stBaseMapper_getCharForValue(baseMapper, h1AlphVal);
        char h2AlphChar = stBaseMapper_getCharForValue(baseMapper, h2AlphVal);
        float h1Prob = gF->haplotypeProbs1[i];
        float h2Prob = gF->haplotypeProbs2[i];
        float genotypeProb = gF->genotypeProbs[i];
        float genotype = (float)gF->genotypeString[i];

        totalLocs++;

        //prep
        bcf_clear1(bcf_rec);
        str.l = 0;

        // contig (CHROM)
        bcf_rec->rid = bcf_hdr_name2id(bcf_hdr, gF->referenceName); //defined in a contig in the top
        // POS
        bcf_rec->pos  = i + gF->refStart; // off by one?
        // ID - skip
        // QUAL - currently writing out the genotype probability
        bcf_rec->qual = genotypeProb;

        // Get phasing info
        gt_info[0] = bcf_gt_phased(0);
        gt_info[1] = bcf_gt_phased(1);

        char refChar = toupper(referenceSeq[i + gF->refStart]);
        if (!differencesOnly) {
            kputc(refChar, &str); // REF
            kputc(',', &str);
            kputc(h1AlphChar, &str);
            kputc(',', &str);
            kputc(h2AlphChar, &str);

            bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
            // FORMAT / $SMPL1
            bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);
            bcf_write1(out, bcf_hdr, bcf_rec);

        } else {
            if (i + 1 >= gF->length) break;

            int next_h1AlphVal = gF->haplotypeString1[i+1];
            int next_h2AlphVal = gF->haplotypeString2[i+1];
            char next_h1AlphChar = stBaseMapper_getCharForValue(baseMapper, next_h1AlphVal);
            char next_h2AlphChar = stBaseMapper_getCharForValue(baseMapper, next_h2AlphVal);

//            fprintf(stderr, "pos: %d\n", bcf_rec->pos);
//            fprintf(stderr, "h1: %c \t h2: %c \n", h1AlphChar, h2AlphChar);
//            fprintf(stderr, "n h1: %c \t n h2: %c \n", next_h1AlphChar, next_h2AlphChar);

            if (next_h1AlphChar != next_h2AlphChar) {
                if (h1AlphChar != '-' && h2AlphChar != '-') {
                    // Check to see if there was an insertion or deletion in the next spot
                    if (next_h1AlphChar == '-' && next_h2AlphChar != '-') {
                        kputc(h1AlphChar, &str);
                        kputc(',', &str);
                        kputc(h2AlphChar, &str);
                        kputc(next_h2AlphChar, &str);
                        int j = 2;
                        while(i+j < gF->length && stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString1[i+j]) == '-' && stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[i+j]) != '-') {
                            kputc(stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[i+j]), &str);
                            j++;
                        }
                        bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
                        bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);
                        bcf_write1(out, bcf_hdr, bcf_rec);
                    } else if (next_h2AlphChar == '-' && next_h1AlphChar != '-') {
                        kputc(h1AlphChar, &str);
                        kputc(next_h1AlphChar, &str);
                        int j = 2;
                        while(i+j < gF->length && stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[i+j]) == '-' && stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString1[i+j]) != '-') {
                            kputc(stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString1[i+j]), &str);
                            j++;
                        }
                        kputc(',', &str);
                        kputc(h2AlphChar, &str);
                        bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
                        bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);
                        bcf_write1(out, bcf_hdr, bcf_rec);

                    } else if (h1AlphChar != h2AlphChar){
                        kputc(h1AlphChar, &str);
                        kputc(',', &str);
                        kputc(h2AlphChar, &str);
                        bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
                        bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);
                        bcf_write1(out, bcf_hdr, bcf_rec);
                    }

                }
            }
            else if (h1AlphChar != h2AlphChar && (h1AlphChar != '-' && h2AlphChar != '-')) {
                // Could also list things that don't match the reference if
                // h1AlphChar != refChar || h2AlphChar != refChar
                kputc(h1AlphChar, &str);
                kputc(',', &str);
                kputc(h2AlphChar, &str);
                bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
                bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);
                bcf_write1(out, bcf_hdr, bcf_rec);
            }

        }
    }

    // cleanup
    free(str.s);
    free(gt_info);
    if (referenceSeq) free(referenceSeq);
    bcf_destroy(bcf_rec);
}
