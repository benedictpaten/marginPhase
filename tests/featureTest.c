/*
 * Copyright (C) 2009-2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <htsIntegration.h>
#include "CuTest.h"
#include "margin.h"

static char *INPUT_PARAMS = "../params/allParams.np.json";
static char *INPUT_BAM = "../tests/data/featureTest.bam";
static char *INPUT_FA = "../tests/data/featureTest.fa";
static char *INPUT_TRUTH_BAM = "../tests/data/featureTestTruth.bam";


int64_t polishingFeatureTest(char *bamFile, char *referenceFile, char *paramsFile, char *featureType,
        char *featureTruth, char *outputName, bool verbose) {

    // Run margin phase
    char *logString = verbose ? "--logLevel DEBUG" : "--logLevel INFO";
    char *featureTruthCmd = featureTruth == NULL ? stString_copy("") : stString_print("--trueReferenceBam %s", featureTruth) ;
    char *command = stString_print("./marginPolish %s %s %s --outputBase %s %s --outputFeatureType %s %s",
            bamFile, referenceFile, paramsFile, outputName, logString, featureType, featureTruthCmd);
    st_logInfo("> Running command: %s\n", command);

    int64_t i = st_system(command);
    free(featureTruthCmd);
    free(command);
    return i;
}


void test_simpleWeightFeatureGeneration(CuTest *testCase) {

    char *featureType = "simpleWeight";
    char *outputName = "test.simpleWeightFeature";
    st_logInfo("\n\nTesting %s feature polishing on %s and %s with %s\n", featureType, INPUT_BAM, INPUT_FA, INPUT_TRUTH_BAM);
    int64_t i = polishingFeatureTest(INPUT_BAM, INPUT_FA, INPUT_PARAMS, featureType, INPUT_TRUTH_BAM, outputName, FALSE);

    char *expectedOutputFa = stString_print("%s.fa", outputName);
    char *expectedOutputFeature = stString_print("%s.simpleWeight.C00000.feature_contig-0-51.tsv", outputName);
    CuAssertTrue(testCase, i == 0);
    CuAssertTrue(testCase, access(expectedOutputFa, F_OK ) == 0);
    CuAssertTrue(testCase, access(expectedOutputFeature, F_OK ) == 0);
}

CuSuite* featureTestSuite(void) {
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, test_simpleWeightFeatureGeneration);

    return suite;
}
