/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com) & Arthur Rand (arand@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <htslib/vcf.h>
#include "stRPHmm.h"

bcf_hdr_t* writeVcfHeader(vcfFile *out, stList *genomeFragments) {
    bcf_hdr_t *hdr = bcf_hdr_init("w");
    kstring_t str = {0,0,NULL};

    // generic info
    str.l = 0;
    ksprintf(&str, "##marginPhase=htslib-%s\n",hts_version());
    bcf_hdr_append(hdr, str.s);

    // reference file used
    str.l = 0;
    ksprintf(&str, "##reference=file://%s\n", "TODO"); //todo
    bcf_hdr_append(hdr, str.s);

    // contigs
    // TODO: assert unique fragments, get full chrom length
    for(int64_t i=0; i<stList_length(genomeFragments); i++) {
        stRPHmm *hmm = stList_get(genomeFragments, i);
        str.l = 0;
        //ksprintf(&str, "##contig=<ID=%s,length=%d>\n", "chr1", 249250621);
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

// This function writes a vcf using one of the haplotypes as the reference
// (not using the actual reference file the reads were originally aligned to)

void writeVcfNoReference(vcfFile *out, bcf_hdr_t *bcf_hdr, stGenomeFragment *gF, stBaseMapper *baseMapper) {

    // intialization
    bcf1_t *bcf_rec = bcf_init1();
    int32_t filter_info = bcf_hdr_id2int(bcf_hdr, BCF_DT_ID, "PASS"); //currently: everything passes
    int32_t *gt_info = (int*)malloc(bcf_hdr_nsamples(bcf_hdr)*2*sizeof(int)); //array specifying phasing
    kstring_t str = {0,0,NULL};

    st_logDebug("Writing out positions with differences between haplotypes\n");
    int numDifferences = 0;
    int totalLocs = 0;

    // iterate over all positions
    for (int64_t i = 0; i < gF->length; ++i) {

        int h1AlphVal = gF->haplotypeString1[i];
        int h2AlphVal = gF->haplotypeString2[i];
        char h1AlphChar = stBaseMapper_getBaseForValue(baseMapper, h1AlphVal);
        char h2AlphChar = stBaseMapper_getBaseForValue(baseMapper, h2AlphVal);

        // Only write characters that differ to vcf
        if (h2AlphChar != h1AlphChar) {
            numDifferences++;
            //prep
            bcf_clear1(bcf_rec);
            str.l = 0;
            int64_t pos = gF->refStart + i;

            // contig (CHROM)
            bcf_rec->rid = bcf_hdr_name2id(bcf_hdr, gF->referenceName); //defined in a contig in the top
            // POS
            bcf_rec->pos  = i + gF->refStart; // off by one?
            // ID - skip
            // QUAL - skip (TODO for now?)
            kputc(h1AlphChar, &str); // REF
            kputc(',', &str);
            kputc(h2AlphChar, &str);
            //st_logDebug("pos: %d, %c-%c\t", pos, h1AlphChar, h2AlphChar);


            bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
            // FORMAT / $SMPL1
            gt_info[0] = bcf_gt_phased(0);
            gt_info[1] = bcf_gt_phased(1);
            bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);

            // save it
            bcf_write1(out, bcf_hdr, bcf_rec);
        }
        totalLocs++;

    }
    st_logDebug("\nNumber of differences between the records: %d\n", numDifferences);
    st_logDebug("Which is %f percent\n", 100 * (float)numDifferences/(float)gF->length);

    // cleanup
    free(str.s); bcf_destroy(bcf_rec);
}



void writeVcfFragment(vcfFile *out, bcf_hdr_t *bcf_hdr, stGenomeFragment *gF, char *referenceSeq, char *referenceName, stBaseMapper *baseMapper) {
    // intialization
    bcf1_t *bcf_rec = bcf_init1();
    int32_t filter_info = bcf_hdr_id2int(bcf_hdr, BCF_DT_ID, "PASS"); //currently: everything passes
    int32_t *gt_info = (int*)malloc(bcf_hdr_nsamples(bcf_hdr)*2*sizeof(int)); //array specifying phasing
    kstring_t str = {0,0,NULL};

    st_logInfo("Reference file name: %s\n", referenceName);
    faidx_t *fai = fai_load(referenceName);
    if ( !fai ) {
        st_logCritical("Could not load fai index of %s.  Maybe you should run 'samtools faidx %s'\n",
                       referenceName, referenceName);
    }
    int seq_len;
    char *ref = fai_fetch(fai, "chr3", &seq_len); //TODO: make this not generic
    if ( seq_len < 0 ) {
        st_logCritical("Failed to fetch reference sequence %s in %s\n", "chr3", referenceName);
    }

    // iterate over all positions
    for (int64_t i = 0; i < gF->length; ++i) {

        int h1AlphVal = gF->haplotypeString1[i];
        int h2AlphVal = gF->haplotypeString2[i];
        char h1AlphChar = stBaseMapper_getBaseForValue(baseMapper, h1AlphVal);
        char h2AlphChar = stBaseMapper_getBaseForValue(baseMapper, h2AlphVal);

        char refChar = referenceSeq[i + gF->refStart];
        int refAlphVal = stBaseMapper_getValueForBase(baseMapper, refChar);

        //prep
        bcf_clear1(bcf_rec);
        str.l = 0;

        // contig (CHROM)
        bcf_rec->rid = bcf_hdr_name2id(bcf_hdr, gF->referenceName); //defined in a contig in the top
        // POS
        bcf_rec->pos  = i + gF->refStart; // off by one?
        // ID - skip
        // QUAL - skip (TODO for now?)

        // Check for where there are actual variations relative to the reference
        kputc(refChar, &str); // REF
        kputc(',', &str);
        kputc(h1AlphChar, &str);
        kputc(',', &str);
        kputc(h2AlphChar, &str);
        if (refChar != h1AlphChar) {
            bcf_rec->qual = 'x';
        }
        if (h1AlphChar != h2AlphChar) {
            bcf_rec->qual = '*';
        }
        bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
        // FORMAT / $SMPL1
        gt_info[0] = bcf_gt_phased(0);
        gt_info[1] = bcf_gt_phased(1);
        bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);

        // save it
        bcf_write1(out, bcf_hdr, bcf_rec);
    }
    // cleanup
    free(str.s); bcf_destroy(bcf_rec);
}
