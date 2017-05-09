/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "stRPHmm.h"


void test_100kbGenotyping(CuTest *testCase) {

    fprintf(stderr, "Testing haplotype inference on NA12878.pb.chr3.100kb.bam\n");

    char *paramsFile = "../tests/params.json";
    char *bamFile = "../tests/NA12878.pb.chr3.100kb.bam";
    char *vcfOutFile = "test_100kb.vcf";
    char *vcfOutFileDiff = "test_100kb_diff.vcf";

    fprintf(stderr, "Parsing parameters\n");
    stBaseMapper *baseMapper = stBaseMapper_construct();
    stRPHmmParameters *params = parseParameters(paramsFile, baseMapper);

    fprintf(stderr, "Creating profile sequences\n");
    stList *profileSequences = stList_construct();
    parseReads(profileSequences, bamFile, baseMapper);

    fprintf(stderr, "Building hmms\n");
    stList *hmms = getRPHmms(profileSequences, params);
    stList *l = stList_construct3(0, (void (*)(void *))stRPHmm_destruct2);
    while(stList_length(hmms) > 0) {
        stList_appendAll(l, stRPHMM_splitWherePhasingIsUncertain(stList_pop(hmms)));
    }
    hmms = l;

    fprintf(stderr, "Writing vcf files\n");
    vcfFile *vcfOutFP = vcf_open(vcfOutFile, "w");
    bcf_hdr_t *bcf_hdr = writeVcfHeader(vcfOutFP, l);

    vcfFile *vcfOutFP_diff = vcf_open(vcfOutFileDiff, "w");
    bcf_hdr_t *bcf_hdr_diff = writeVcfHeader(vcfOutFP_diff, l);

    kstring_t str = {0,0,NULL};

    for(int64_t i=0; i<stList_length(hmms); i++) {
        stRPHmm *hmm = stList_get(hmms, i);

        // Run the forward-backward algorithm
        stRPHmm_forwardBackward(hmm);

        // Now compute a high probability path through the hmm
        stList *path = stRPHmm_forwardTraceBack(hmm);

        // Compute the genome fragment
        stGenomeFragment *gF = stGenomeFragment_construct(hmm, path);

        int numDifferences = 0;

        bcf1_t *bcf_rec = bcf_init1();
        int32_t filter_info = bcf_hdr_id2int(bcf_hdr, BCF_DT_ID, "PASS"); //currently: everything passes
        int32_t *gt_info = (int*)malloc(bcf_hdr_nsamples(bcf_hdr)*2*sizeof(int)); //array specifying phasing
        kstring_t str = {0,0,NULL};

        // iterate over all positions
        for (int64_t i = 0; i < gF->length; ++i) {

            int h1AlphVal = gF->haplotypeString1[i];
            int h2AlphVal = gF->haplotypeString2[i];
            char h1AlphChar = stBaseMapper_getBaseForValue(baseMapper, h1AlphVal);
            char h2AlphChar = stBaseMapper_getBaseForValue(baseMapper, h2AlphVal);

            //prep
            bcf_clear1(bcf_rec);
            str.l = 0;
            int64_t pos = gF->refStart + i;

            // contig (CHROM)
            bcf_rec->rid = bcf_hdr_name2id(bcf_hdr, gF->referenceName); //defined in a contig in the top
            // POS
            bcf_rec->pos  = i + gF->refStart; // off by one?
            // ID - skip

            kputc(h1AlphChar, &str); // REF
            kputc(',', &str);
            kputc(h2AlphChar, &str);

            bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
            // FORMAT / $SMPL1
            gt_info[0] = bcf_gt_phased(0);
            gt_info[1] = bcf_gt_phased(1);
            bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);

            if (h2AlphChar != h1AlphChar) {
                numDifferences++;
                // save it
                bcf_write1(vcfOutFP, bcf_hdr, bcf_rec);
                bcf_write1(vcfOutFP_diff, bcf_hdr, bcf_rec);
            } else {
                // save it
                bcf_rec->qual = 0;
                bcf_write1(vcfOutFP, bcf_hdr, bcf_rec);
            }
        }
        fprintf(stderr, "\nNumber of differences between the records: %d\n", numDifferences);
        fprintf(stderr, "Which is %f percent\n", 100 * (float)numDifferences/(float)gF->length);
    }
}

CuSuite *marginPhaseTestSuite(void) {
    CuSuite* suite = CuSuiteNew();

    // System level tests
    SUITE_ADD_TEST(suite, test_100kbGenotyping);

    return suite;
}