/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */


#include <htslib/vcf.h>
#include "CuTest.h"
#include "sonLib.h"


int64_t genotypingTest2(char *paramsFile, char *bamFile, char *outputBase, char *referenceFile, char *vcfReference,
                        char* singleNuclProbDir, bool verbose) {

    // Run margin phase
    char *logString = verbose ? "--logLevel DEBUG" : "--logLevel INFO";
    char *singleNuclProbCmd = singleNuclProbDir == NULL ? "" :
                              stString_print("--singleNuclProbDir %s", singleNuclProbDir);
    char *command = stString_print("./marginPhase %s %s %s %s --outputBase %s --referenceVcf %s %s",
                                   bamFile, referenceFile, paramsFile, logString, outputBase, vcfReference,
                                   singleNuclProbCmd);
    st_logInfo("> Running command: %s\n", command);

    int64_t i = st_system(command);
    free(command);
    if (singleNuclProbDir != NULL) free(singleNuclProbCmd);
    return i;
}

int64_t genotypingTest(char *paramsFile, char *bamFile, char *outputBase,
                       char *referenceFile, char *vcfReference, bool verbose) {

    return genotypingTest2(paramsFile, bamFile, outputBase, referenceFile, vcfReference, NULL, verbose);
}


void test_5kbGenotyping(CuTest *testCase) {

    char *paramsFile = "../params/params.pacbio.json";
    char *referenceFile = "../tests/hg19.chr3.9mb.fa";
    char *outputBase = "test_5kb";
    bool verbose = true;

    char *bamFile = "../tests/NA12878.pb.chr3.5kb.bam";
    // TODO: create vcf specifically for this 5 kb region

    char *vcfReference = "../tests/NA12878.PG.chr3.100kb.0.vcf";


    st_logInfo("\n\nTesting haplotype inference on %s\n", bamFile);
    int64_t i = genotypingTest(paramsFile, bamFile, outputBase,
                               referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);
}

void test_5kbGenotyping_singleNuclProb(CuTest *testCase) {

    char *paramsFile = "../params/params.pacbio.json";
    char *referenceFile = "../tests/hg19.chr3.9mb.fa";
    char *outputBase = "test_5kb_singleNuclProb";
    char *bamFile = "../tests/NA12878.pb.chr3.5kb.bam";
    char *vcfReference = "../tests/NA12878.PG.chr3.100kb.0.vcf";
    char *singleNuclProbDir = "../tests/NA12878.pb.chr3.5kb";
    bool verbose = true;

    st_logInfo("\n\nTesting haplotype inference on %s with singleNuclProbs\n", bamFile);

    char *command = "tar xvf ../tests/NA12878.pb.chr3.5kb.singleNuclProb.zip --directory ../tests";
    int64_t i = st_system(command);
    CuAssertTrue(testCase, i == 0);

    st_logInfo("\nRunning command: %s\n\n", command);
    i = genotypingTest2(paramsFile, bamFile, outputBase, referenceFile, vcfReference, singleNuclProbDir, verbose);
    CuAssertTrue(testCase, i == 0);
}

void test_100kbGenotyping_pacbio(CuTest *testCase) {

    char *paramsFile = "../params/params.pacbio.json";
    char *referenceFile = "../tests/hg19.chr3.9mb.fa";
    char *outputBase = "test_100kb_pb";
    char *bamFile = "../tests/NA12878.pb.chr3.100kb.4.bam";
    char *vcfReference = "../tests/NA12878.PG.chr3.100kb.4.vcf";
    bool verbose = false;

    st_logInfo("\n\nTesting haplotype inference on %s\n", bamFile);

    int64_t i = genotypingTest(paramsFile, bamFile, outputBase,
                               referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);

}

void test_100kbGenotyping_nanopore(CuTest *testCase) {

    char *paramsFile = "../params/params.nanopore.json";
    char *referenceFile = "../tests/hg19.chr3.9mb.fa";
    char *outputBase = "test_100kb_np";
    char *bamFile = "../tests/NA12878.np.chr3.100kb.4.bam";
    char *vcfReference = "../tests/NA12878.PG.chr3.100kb.4.vcf";
    bool verbose = false;

    st_logInfo("\n\nTesting haplotype inference on %s\n", bamFile);

    int64_t i = genotypingTest(paramsFile, bamFile, outputBase,
                               referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);

}

void test_multiple100kbGenotyping_pacbio(CuTest *testCase) {

    st_logInfo("\n\nTesting all PacBio regions\n");

    char *paramsFile = "../params/params.pacbio.json";
    char *referenceFile = "../tests/hg19.chr3.9mb.fa";
    bool verbose = false;

    char *bamFile = "../tests/NA12878.pb.chr3.100kb.0.bam";
    char *vcfReference = "../tests/NA12878.PG.chr3.100kb.0.vcf";
    char *outputBase = "test_100kb_pb_0";
    st_logInfo("\nTesting haplotype inference on %s\n", bamFile);

    int64_t i = genotypingTest(paramsFile, bamFile, outputBase,
                           referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);

    bamFile = "../tests/NA12878.pb.chr3.100kb.1.bam";
    vcfReference = "../tests/NA12878.PG.chr3.100kb.1.vcf";
    outputBase = "test_100kb_pb_1";
    st_logInfo("\nTesting haplotype inference on %s\n", bamFile);

    i = genotypingTest(paramsFile, bamFile, outputBase,
                       referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);

    bamFile = "../tests/NA12878.pb.chr3.100kb.2.bam";
    vcfReference = "../tests/NA12878.PG.chr3.100kb.2.vcf";
    outputBase = "test_100kb_pb_2";
    st_logInfo("\nTesting haplotype inference on %s\n", bamFile);

    i = genotypingTest(paramsFile, bamFile, outputBase,
                       referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);


    bamFile = "../tests/NA12878.pb.chr3.100kb.3.bam";
    vcfReference = "../tests/NA12878.PG.chr3.100kb.3.vcf";
    outputBase = "test_100kb_pb_3";
    st_logInfo("\nTesting haplotype inference on %s\n", bamFile);

    i = genotypingTest(paramsFile, bamFile, outputBase,
                       referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);

    bamFile = "../tests/NA12878.pb.chr3.100kb.4.bam";
    vcfReference = "../tests/NA12878.PG.chr3.100kb.4.vcf";
    outputBase = "test_100kb_pb_4";
    st_logInfo("\nTesting haplotype inference on %s\n", bamFile);

    i = genotypingTest(paramsFile, bamFile, outputBase,
                       referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);
}

void test_multiple100kbGenotyping_nanopore(CuTest *testCase) {

    st_logInfo("\n\nTesting all nanopore regions\n");

    char *paramsFile = "../params/params.nanopore.json";
    char *referenceFile = "../tests/hg19.chr3.9mb.fa";
    bool verbose = false;

    char *bamFile = "../tests/NA12878.np.chr3.100kb.0.bam";
    char *vcfReference = "../tests/NA12878.PG.chr3.100kb.0.vcf";
    char *outputBase = "test_100kb_np_0";
    st_logInfo("\nTesting haplotype inference on %s\n", bamFile);

    int64_t i = genotypingTest(paramsFile, bamFile, outputBase,
                           referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);

    bamFile = "../tests/NA12878.np.chr3.100kb.1.bam";
    vcfReference = "../tests/NA12878.PG.chr3.100kb.1.vcf";
    outputBase = "test_100kb_np_1";
    st_logInfo("\nTesting haplotype inference on %s\n", bamFile);

    i = genotypingTest(paramsFile, bamFile, outputBase,
                       referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);

    bamFile = "../tests/NA12878.np.chr3.100kb.2.bam";
    vcfReference = "../tests/NA12878.PG.chr3.100kb.2.vcf";
    outputBase = "test_100kb_np_2";
    st_logInfo("\nTesting haplotype inference on %s\n", bamFile);

    i = genotypingTest(paramsFile, bamFile, outputBase,
                       referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);


    bamFile = "../tests/NA12878.np.chr3.100kb.3.bam";
    vcfReference = "../tests/NA12878.PG.chr3.100kb.3.vcf";
    outputBase = "test_100kb_np_3";
    st_logInfo("\nTesting haplotype inference on %s\n", bamFile);

    i = genotypingTest(paramsFile, bamFile, outputBase,
                       referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);

    bamFile = "../tests/NA12878.np.chr3.100kb.4.bam";
    vcfReference = "../tests/NA12878.PG.chr3.100kb.4.vcf";
    outputBase = "test_100kb_np_4";
    st_logInfo("\nTesting haplotype inference on %s\n", bamFile);

    i = genotypingTest(paramsFile, bamFile, outputBase,
                       referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);
}

void test_multiple100kbGenotyping_illuminaHiSeq(CuTest *testCase) {

    st_logInfo("\n\nTesting all Illumina HiSeq regions\n");

    char *paramsFile = "../params/params.pacbio.json";
    char *referenceFile = "../tests/hg19.chr3.9mb.fa";
    bool verbose = false;

    char *bamFile = "../tests/NA12878.ihs.chr3.100kb.0.bam";
    char *vcfReference = "../tests/NA12878.PG.chr3.100kb.0.vcf";
    char *outputBase = "test_100kb_ihs_0";
    st_logInfo("\nTesting haplotype inference on %s\n", bamFile);

    int64_t i = genotypingTest(paramsFile, bamFile, outputBase,
                           referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);

    bamFile = "../tests/NA12878.ihs.chr3.100kb.1.bam";
    vcfReference = "../tests/NA12878.PG.chr3.100kb.1.vcf";
    outputBase = "test_100kb_ihs_1";
    st_logInfo("\nTesting haplotype inference on %s\n", bamFile);

    i = genotypingTest(paramsFile, bamFile, outputBase,
                       referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);

    bamFile = "../tests/NA12878.ihs.chr3.100kb.2.bam";
    vcfReference = "../tests/NA12878.PG.chr3.100kb.2.vcf";
    outputBase = "test_100kb_ihs_2";
    st_logInfo("\nTesting haplotype inference on %s\n", bamFile);

    i = genotypingTest(paramsFile, bamFile, outputBase,
                       referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);

    bamFile = "../tests/NA12878.ihs.chr3.100kb.3.bam";
    vcfReference = "../tests/NA12878.PG.chr3.100kb.3.vcf";
    outputBase = "test_100kb_ihs_3";
    st_logInfo("\nTesting haplotype inference on %s\n", bamFile);

    i = genotypingTest(paramsFile, bamFile, outputBase,
                       referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);

    bamFile = "../tests/NA12878.ihs.chr3.100kb.4.bam";
    vcfReference = "../tests/NA12878.PG.chr3.100kb.4.vcf";
    outputBase = "test_100kb_ihs_4";
    st_logInfo("\nTesting haplotype inference on %s\n", bamFile);

    i = genotypingTest(paramsFile, bamFile, outputBase,
                       referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);
}

CuSuite *marginPhaseTestSuite(void) {

    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, test_5kbGenotyping);
    SUITE_ADD_TEST(suite, test_5kbGenotyping_singleNuclProb);
    SUITE_ADD_TEST(suite, test_100kbGenotyping_pacbio);
    SUITE_ADD_TEST(suite, test_100kbGenotyping_nanopore);

//    SUITE_ADD_TEST(suite, test_multiple100kbGenotyping_pacbio);
//    SUITE_ADD_TEST(suite, test_multiple100kbGenotyping_nanopore);
//    SUITE_ADD_TEST(suite, test_multiple100kbGenotyping_illuminaHiSeq);

    return suite;
}
