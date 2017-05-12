/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "stRPHmm.h"


void genotypingTest(char *paramsFile, char *bamFile, char *vcfOutFile, char *vcfOutFileDiff, char *referenceFile) {
    fprintf(stderr, "Parsing parameters\n");
    stBaseMapper *baseMapper = stBaseMapper_construct();
    stRPHmmParameters *params = parseParameters(paramsFile, baseMapper);
    // Print a report of the parsed parameters
    stRPHmmParameters_printParameters(params, stderr);

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


    for(int64_t i=0; i<stList_length(hmms); i++) {
        stRPHmm *hmm = stList_get(hmms, i);

        // Run the forward-backward algorithm
        stRPHmm_forwardBackward(hmm);

        // Now compute a high probability path through the hmm
        stList *path = stRPHmm_forwardTraceBack(hmm);

        // Compute the genome fragment
        stGenomeFragment *gF = stGenomeFragment_construct(hmm, path);

        writeVcfFragment(vcfOutFP, bcf_hdr, gF, referenceFile, baseMapper, true);
        writeVcfFragment(vcfOutFP_diff, bcf_hdr_diff, gF, NULL, baseMapper, false);
    }

    vcf_close(vcfOutFP);
    vcf_close(vcfOutFP_diff);
}

void test_5kbGenotyping(CuTest *testCase) {

    fprintf(stderr, "Testing haplotype inference on NA12878.pb.chr3.5kb.bam\n");

    char *paramsFile = "../tests/params.json";
    char *bamFile = "../tests/NA12878.pb.chr3.5kb.bam";
    char *vcfOutFile = "test_5kb.vcf";
    char *vcfOutFileDiff = "test_5kb_diff.vcf";
    char *referenceFile = "../tests/hg19.chr3.9mb.fa";

    genotypingTest(paramsFile, bamFile, vcfOutFile, vcfOutFileDiff, referenceFile);
}

void test_100kbGenotyping(CuTest *testCase) {

    fprintf(stderr, "Testing haplotype inference on NA12878.pb.chr3.100kb.0.bam\n");

    char *paramsFile = "../tests/params.json";
    char *bamFile = "../tests/NA12878.pb.chr3.100kb.0.bam";
    char *vcfOutFile = "test_100kb.vcf";
    char *vcfOutFileDiff = "test_100kb_diff.vcf";
    char *referenceFile = "../tests/hg19.chr3.9mb.fa";

    genotypingTest(paramsFile, bamFile, vcfOutFile, vcfOutFileDiff, referenceFile);
}

CuSuite *marginPhaseTestSuite(void) {
    CuSuite* suite = CuSuiteNew();

    // System level tests
    SUITE_ADD_TEST(suite, test_5kbGenotyping);
    SUITE_ADD_TEST(suite, test_100kbGenotyping);

    return suite;
}
