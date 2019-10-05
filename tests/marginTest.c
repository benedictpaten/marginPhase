/*
 * Copyright (C) 2009-2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "margin.h"

static char *polishParamsFile = "../params/allParams.np.json";

int64_t marginIntegrationTest(char *bamFile, char *referenceFile, char *paramsFile, char *region, bool verbose, bool diploid) {

    // Run margin phase
    char *logString = verbose ? "--logLevel DEBUG" : "--logLevel INFO";
    char *regionStr = region == NULL ? stString_print("") : stString_print("--region %s", region);
    char *diploidString = diploid ? "--diploid" : "";
    char *command = stString_print("./margin %s %s %s %s %s %s", bamFile, referenceFile, paramsFile, regionStr, logString, diploidString);
    st_logInfo("> Running command: %s\n", command);

    int64_t i = st_system(command);
    free(regionStr);
    free(command);
    return i;
}

void test_marginIntegration(CuTest *testCase) {
    char *referenceFile = "../tests/shasta_diploid/shasta_phasing_test.ref.fasta";
    bool verbose = false;
    char *bamFile = "../tests/shasta_diploid/shasta_phasing_test.align.bam";
    char *region = NULL;

    st_logInfo("\tTesting diploid polishing on %s\n", bamFile);
    int i = marginIntegrationTest(bamFile, referenceFile, polishParamsFile, region, verbose, 1);
    CuAssertTrue(testCase, i == 0);
}

CuSuite* marginIntegrationTestSuite(void) {
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, test_marginIntegration);


    return suite;
}
