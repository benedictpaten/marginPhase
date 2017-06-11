/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */


#include <htslib/vcf.h>
#include <htslib/sam.h>
#include "CuTest.h"
#include "sonLib.h"
#include "stRPHmm.h"
#include "../externalTools/sonLib/C/impl/sonLibListPrivate.h"

int genotypingTest(char *paramsFile, char *bamFile, char *vcfOutFile,
        char *referenceFile, char *vcfReference, bool verbose) {

    // Run margin phase
    char *logString = verbose ? "--logLevel DEBUG" : "--logLevel INFO";
    char *command = stString_print("./marginPhase --bamFile %s --referenceFasta %s %s --params %s --vcfFile %s "
            " --referenceVcf %s",
            bamFile, referenceFile, logString,
            paramsFile, vcfOutFile, vcfReference);
    st_logInfo("> Running margin phase on %s\n", bamFile);
    return st_system(command);

    // TODO : Do VCF comparison using VCF eval
}

void test_5kbGenotyping(CuTest *testCase) {

    char *paramsFile = "../params.json";
    char *referenceFile = "../tests/hg19.chr3.9mb.fa";
    char *vcfOutFile = "test_5kb.vcf";
    char *vcfOutFileDiff = "test_5kb_diff.vcf";
    char *samOutBase = "test_100kb";
    bool verbose = true;

    char *bamFile = "../tests/NA12878.pb.chr3.5kb.bam";
    // TODO: create vcf specifically for this 5 kb region

    char *vcfReference = "../tests/NA12878.PG.chr3.100kb.0.vcf";


    st_logInfo("Testing haplotype inference on %s\n", bamFile);
    int i = genotypingTest(paramsFile, bamFile, vcfOutFile,
                            referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);
}

void test_100kbGenotyping(CuTest *testCase) {

    char *paramsFile = "../params.json";
    char *referenceFile = "../tests/hg19.chr3.9mb.fa";
    char *vcfOutFile = "test_100kb.vcf";
    char *vcfOutFileDiff = "test_100kb_diff.vcf";
    char *samOutBase = "test_100kb";
    bool verbose = true;

    char *bamFile = "../tests/NA12878.pb.chr3.100kb.4.bam";
//    char *bamFile = "../tests/NA12878.np.chr3.100kb.4.bam";
//    char *bamFile = "../tests/NA12878.ihs.chr3.100kb.3.bam";
    char *vcfReference = "../tests/NA12878.PG.chr3.100kb.4.vcf";
//    char *vcfReference = "../tests/HG001.GRCh37.chr3.100kb.vcf";
//    char *bamFile = "../tests/NA12878.pb.chr3.2mb.bam";
//    char *vcfReference = "../tests/HG001.GRCh37.chr3.2mb.vcf";

    st_logInfo("Testing haplotype inference on %s\n", bamFile);

    int i = genotypingTest(paramsFile, bamFile, vcfOutFile,
                        referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);

}

void test_multiple100kbGenotyping(CuTest *testCase) {

    st_setLogLevelFromString("info");

    char *paramsFile = "../params.json";
    char *referenceFile = "../tests/hg19.chr3.9mb.fa";
    char *vcfOutFile = "test_100kb.vcf";
    char *vcfOutFileDiff = "test_100kb_diff.vcf";
    char *samOutBase = "test_100kb_0";
    bool verbose = true;

    char *bamFile = "../tests/NA12878.pb.chr3.100kb.0.bam";
    char *vcfReference = "../tests/NA12878.PG.chr3.100kb.0.vcf";
    st_logInfo("Testing haplotype inference on %s\n", bamFile);

    int i = genotypingTest(paramsFile, bamFile, vcfOutFile,
                           referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);

    bamFile = "../tests/NA12878.pb.chr3.100kb.1.bam";
    vcfReference = "../tests/NA12878.PG.chr3.100kb.1.vcf";
    samOutBase = "test_100kb_1";
    st_logInfo("Testing haplotype inference on %s\n", bamFile);

    i = genotypingTest(paramsFile, bamFile, vcfOutFile,
                       referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);

    bamFile = "../tests/NA12878.pb.chr3.100kb.2.bam";
    vcfReference = "../tests/NA12878.PG.chr3.100kb.2.vcf";
    samOutBase = "test_100kb_2";
    st_logInfo("Testing haplotype inference on %s\n", bamFile);

    i = genotypingTest(paramsFile, bamFile, vcfOutFile,
                       referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);


    bamFile = "../tests/NA12878.pb.chr3.100kb.3.bam";
    vcfReference = "../tests/NA12878.PG.chr3.100kb.3.vcf";
    samOutBase = "test_100kb_3";
    st_logInfo("Testing haplotype inference on %s\n", bamFile);

    i = genotypingTest(paramsFile, bamFile, vcfOutFile,
                       referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);

    bamFile = "../tests/NA12878.pb.chr3.100kb.4.bam";
    vcfReference = "../tests/NA12878.PG.chr3.100kb.4.vcf";
    samOutBase = "test_100kb_4";
    st_logInfo("Testing haplotype inference on %s\n", bamFile);

    i = genotypingTest(paramsFile, bamFile, vcfOutFile,
                       referenceFile, vcfReference, verbose);
    CuAssertTrue(testCase, i == 0);
}

CuSuite *marginPhaseTestSuite(void) {

    CuSuite* suite = CuSuiteNew();

//    st_setLogLevelFromString("debug");
    st_setLogLevelFromString("info");

//    SUITE_ADD_TEST(suite, test_5kbGenotyping);
    SUITE_ADD_TEST(suite, test_100kbGenotyping);
//    SUITE_ADD_TEST(suite, test_multiple100kbGenotyping);

    return suite;
}
