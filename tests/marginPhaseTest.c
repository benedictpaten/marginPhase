/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "margin.h"

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

/*
 * Test for a 5kb region
 */
void test_5kbGenotyping(CuTest *testCase) {

    char *paramsFile = "../params/params.pacbio.json";
    char *referenceFile = "../tests/data/hg19.chr3.9mb.fa";
    char *outputBase = "test_5kb";
    bool verbose = true;

    char *bamFile = "../tests/data/NA12878.pb.chr3.5kb.bam";
    // TODO: create vcf specifically for this 5 kb region

    char *vcfReference = "../tests/data/NA12878.PG.chr3.100kb.0.vcf";


    st_logInfo("\n\nTesting haplotype inference on %s\n", bamFile);
    int64_t i = genotypingTest(paramsFile, bamFile, outputBase,
                               referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);
}

/*
 * Test for a 5kb region with single nucleotide probabilities
 */
void test_5kbGenotyping_singleNuclProb(CuTest *testCase) {

    char *paramsFile = "../params/params.pacbio.json";
    char *referenceFile = "../tests/data/hg19.chr3.9mb.fa";
    char *outputBase = "test_5kb_singleNuclProb";
    char *bamFile = "../tests/data/NA12878.pb.chr3.5kb.bam";
    char *vcfReference = "../tests/data/NA12878.PG.chr3.100kb.0.vcf";
    char *singleNuclProbDir = "../tests/data/NA12878.pb.chr3.5kb";
    bool verbose = true;

    st_logInfo("\n\nTesting haplotype inference on %s with singleNuclProbs\n", bamFile);

    char *command = "tar xvf ../tests/data/NA12878.pb.chr3.5kb.singleNuclProb.zip --directory ../tests/data";
    int64_t i = st_system(command);
    CuAssertTrue(testCase, i == 0);

    st_logInfo("\nRunning command: %s\n\n", command);
    i = genotypingTest2(paramsFile, bamFile, outputBase, referenceFile, vcfReference, singleNuclProbDir, verbose);
    CuAssertTrue(testCase, i == 0);
}

/*
 * Test for a 100 kb region with PacBio reads
 */
void test_100kbGenotyping_pacbio(CuTest *testCase) {

    char *paramsFile = "../params/params.pacbio.json";
    char *referenceFile = "../tests/data/hg19.chr3.9mb.fa";
    char *outputBase = "test_100kb_pb";
    char *bamFile = "../tests/data/NA12878.pb.chr3.100kb.4.bam";
    char *vcfReference = "../tests/data/NA12878.PG.chr3.100kb.4.vcf";
    bool verbose = false;

    st_logInfo("\n\nTesting haplotype inference on %s\n", bamFile);

    int64_t i = genotypingTest(paramsFile, bamFile, outputBase,
                               referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);
}

/*
 * Test for a 100kb region with Oxford Nanopore reads
 */
void test_100kbGenotyping_nanopore(CuTest *testCase) {

    char *paramsFile = "../params/params.nanopore.json";
    char *referenceFile = "../tests/data/hg19.chr3.9mb.fa";
    char *outputBase = "test_100kb_np";
    char *bamFile = "../tests/data/NA12878.np.chr3.100kb.4.bam";
    char *vcfReference = "../tests/data/NA12878.PG.chr3.100kb.4.vcf";
    bool verbose = false;

    st_logInfo("\n\nTesting haplotype inference on %s\n", bamFile);

    int64_t i = genotypingTest(paramsFile, bamFile, outputBase,
                               referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);
}

/*
 * Test to run on five 100kb regions for PacBio
 */
void test_multiple100kbGenotyping_pacbio(CuTest *testCase) {

    st_logInfo("\n\nTesting all PacBio regions\n");

    // Params for test
    char *paramsFile = "../params/params.pacbio.json";
    char *referenceFile = "../tests/data/hg19.chr3.9mb.fa";
    bool verbose = false;

    // Test 1
    char *bamFile = "../tests/data/NA12878.pb.chr3.100kb.0.bam";
    char *vcfReference = "../tests/data/NA12878.PG.chr3.100kb.0.vcf";
    char *outputBase = "test_100kb_pb_0";
    st_logInfo("\nTesting haplotype inference on %s\n", bamFile);

    int64_t i = genotypingTest(paramsFile, bamFile, outputBase,
                           referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);

    // Test 2
    bamFile = "../tests/data/NA12878.pb.chr3.100kb.1.bam";
    vcfReference = "../tests/data/NA12878.PG.chr3.100kb.1.vcf";
    outputBase = "test_100kb_pb_1";
    st_logInfo("\nTesting haplotype inference on %s\n", bamFile);

    i = genotypingTest(paramsFile, bamFile, outputBase,
                       referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);

    // Test 3
    bamFile = "../tests/data/NA12878.pb.chr3.100kb.2.bam";
    vcfReference = "../tests/data/NA12878.PG.chr3.100kb.2.vcf";
    outputBase = "test_100kb_pb_2";
    st_logInfo("\nTesting haplotype inference on %s\n", bamFile);

    i = genotypingTest(paramsFile, bamFile, outputBase,
                       referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);

    // Test 4
    bamFile = "../tests/data/NA12878.pb.chr3.100kb.3.bam";
    vcfReference = "../tests/data/NA12878.PG.chr3.100kb.3.vcf";
    outputBase = "test_100kb_pb_3";
    st_logInfo("\nTesting haplotype inference on %s\n", bamFile);

    i = genotypingTest(paramsFile, bamFile, outputBase,
                       referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);

    // Test 5
    bamFile = "../tests/data/NA12878.pb.chr3.100kb.4.bam";
    vcfReference = "../tests/data/NA12878.PG.chr3.100kb.4.vcf";
    outputBase = "test_100kb_pb_4";
    st_logInfo("\nTesting haplotype inference on %s\n", bamFile);

    i = genotypingTest(paramsFile, bamFile, outputBase,
                       referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);
}

/*
 * Test to run on five 100kb regions for Oxford Nanopore
 */
void test_multiple100kbGenotyping_nanopore(CuTest *testCase) {

    st_logInfo("\n\nTesting all nanopore regions\n");

    // Params for test
    char *paramsFile = "../params/params.nanopore.json";
    char *referenceFile = "../tests/data/hg19.chr3.9mb.fa";
    bool verbose = false;

    // Test 1
    char *bamFile = "../tests/data/NA12878.np.chr3.100kb.0.bam";
    char *vcfReference = "../tests/data/NA12878.PG.chr3.100kb.0.vcf";
    char *outputBase = "test_100kb_np_0";
    st_logInfo("\nTesting haplotype inference on %s\n", bamFile);

    int64_t i = genotypingTest(paramsFile, bamFile, outputBase,
                           referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);

    // Test 2
    bamFile = "../tests/data/NA12878.np.chr3.100kb.1.bam";
    vcfReference = "../tests/data/NA12878.PG.chr3.100kb.1.vcf";
    outputBase = "test_100kb_np_1";
    st_logInfo("\nTesting haplotype inference on %s\n", bamFile);

    i = genotypingTest(paramsFile, bamFile, outputBase,
                       referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);

    // Test 3
    bamFile = "../tests/data/NA12878.np.chr3.100kb.2.bam";
    vcfReference = "../tests/data/NA12878.PG.chr3.100kb.2.vcf";
    outputBase = "test_100kb_np_2";
    st_logInfo("\nTesting haplotype inference on %s\n", bamFile);

    i = genotypingTest(paramsFile, bamFile, outputBase,
                       referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);

    // Test 4
    bamFile = "../tests/data/NA12878.np.chr3.100kb.3.bam";
    vcfReference = "../tests/data/NA12878.PG.chr3.100kb.3.vcf";
    outputBase = "test_100kb_np_3";
    st_logInfo("\nTesting haplotype inference on %s\n", bamFile);

    i = genotypingTest(paramsFile, bamFile, outputBase,
                       referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);

    // Test 5
    bamFile = "../tests/data/NA12878.np.chr3.100kb.4.bam";
    vcfReference = "../tests/data/NA12878.PG.chr3.100kb.4.vcf";
    outputBase = "test_100kb_np_4";
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

    return suite;
}
