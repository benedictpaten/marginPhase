/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */


#include <htslib/vcf.h>
#include "CuTest.h"
#include "sonLib.h"


int64_t genotypingTest(char *paramsFile, char *bamFile, char *outputBase,
                    char *referenceFile, char *vcfReference, bool verbose) {

    // Run margin phase
    char *logString = verbose ? "--logLevel DEBUG" : "--logLevel INFO";
    char *command = stString_print("./marginPhase %s %s %s --params %s --outputBase %s "
                                           " --referenceVcf %s",
                                   bamFile, referenceFile, logString,
                                   paramsFile, outputBase, vcfReference);
    st_logInfo("Running command: %s\n", command);
    st_logInfo("> Running margin phase on %s\n", bamFile);

    int64_t i = st_system(command);
    free(command);
    return i;
    // TODO : Do VCF comparison using VCF eval
}


void test_5kbGenotyping(CuTest *testCase) {

    char *paramsFile = "../params_pacbio_gaps.json";
    char *referenceFile = "../tests/hg19.chr3.9mb.fa";
    char *outputBase = "test_5kb";
    bool verbose = true;

    char *bamFile = "../tests/NA12878.pb.chr3.5kb.bam";
    // TODO: create vcf specifically for this 5 kb region

    char *vcfReference = "../tests/NA12878.PG.chr3.100kb.0.vcf";


    st_logInfo("Testing haplotype inference on %s\n", bamFile);
    int64_t i = genotypingTest(paramsFile, bamFile, outputBase,
                            referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);
}

void test_100kbGenotyping(CuTest *testCase) {

//    char *paramsFile = "../params_pacbio_gaps.json";
//    char *paramsFile = "../params_pacbio_currentbest.json";
//    char *paramsFile = "../params_nanopore_currentbest.json";
//    char *paramsFile = "../../params/params.np.averaged.9995.rp.json";
    char *paramsFile = "../../params/params.np.transitions.9995.rp.json";
    char *referenceFile = "../tests/hg19.chr3.9mb.fa";
    char *outputBase = "test_100kb";
    bool verbose = false;

//    char *bamFile = "../tests/NA12878.pb.chr3.100kb.4.bam";
    char *bamFile = "../tests/NA12878.np.chr3.100kb.4.bam";
//    char *bamFile = "../tests/NA12878.ihs.chr3.100kb.4.bam";
    char *vcfReference = "../tests/NA12878.PG.chr3.100kb.4.vcf";
//    char *vcfReference = "../tests/HG001.GRCh37.chr3.100kb.vcf";

    st_logInfo("\nTesting haplotype inference on %s\n", bamFile);

    int64_t i = genotypingTest(paramsFile, bamFile, outputBase,
                        referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);

}

void test_multiple100kbGenotyping_pacbio(CuTest *testCase) {

    st_logInfo("Testing all PacBio regions\n");

//    char *paramsFile = "../../params/params.pb.transitions.998.plain.json";
//    char *paramsFile = "../../params/params.pb.averaged.9995.plain.testing.json";
    char *paramsFile = "../params_pacbio_currentbest.json";
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

    st_logInfo("Testing all nanopore regions\n");

//    char *paramsFile = "../params_nanopore_currentbest.json";
//        char *paramsFile = "../../params/params.np.averaged.9995.rp.json";
    char *paramsFile = "../../params/params.np.transitions.9995.rp.json";
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

    st_logInfo("Testing all Illumina HiSeq regions\n");

    char *paramsFile = "../params.json";
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

//    SUITE_ADD_TEST(suite, test_5kbGenotyping);
//    SUITE_ADD_TEST(suite, test_100kbGenotyping);
//    SUITE_ADD_TEST(suite, test_multiple100kbGenotyping_pacbio);
//    SUITE_ADD_TEST(suite, test_multiple100kbGenotyping_nanopore);
//    SUITE_ADD_TEST(suite, test_multiple100kbGenotyping_illuminaHiSeq);

    return suite;
}
