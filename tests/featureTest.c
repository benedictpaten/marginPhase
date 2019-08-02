/*
 * Copyright (C) 2009-2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <stdio.h>
#include <htsIntegration.h>
#include "CuTest.h"
#include "margin.h"
#include "helenFeatures.h"

static char *FEATURE_TEST_PARAMS = "../params/allParams.np.json";
static char *FEATURE_TEST_NO_RLE_PARAMS = "../params/allParams.np.no_rle.json";
static char *FEATURE_TEST_BAM = "../tests/data/featureTest/featureTest.bam";
static char *FEATURE_TEST_FA = "../tests/data/featureTest/featureTest.fa";
static char *FEATURE_TEST_TRUTH_BAM = "../tests/data/featureTest/featureTestTruth.bam";
static char *FEATURE_TEST_TRUTH_SEQ = "ACGATAACCGGTTAAACCCCGGGTTTCAAACCCCGGGGTTGATTACACAT";


int64_t polishingFeatureTest(char *bamFile, char *referenceFile, char *paramsFile, char *featureType,
        char *featureTruth, char *outputName, bool verbose) {

    // Run margin phase
    char *logString = verbose ? "--logLevel DEBUG" : "--logLevel INFO";
    char *featureTruthCmd = featureTruth == NULL ? stString_copy("") : stString_print("--trueReferenceBam %s", featureTruth) ;
    char *command = stString_print("./marginPolish %s %s %s --outputBase %s %s -F %s %s",
            bamFile, referenceFile, paramsFile, outputName, logString, featureType, featureTruthCmd);
    st_logInfo("> Running command: %s\n", command);
    fprintf(stderr, "> Running command: %s\n", command);

    int64_t i = st_system(command);
    free(featureTruthCmd);
    free(command);
    return i;
}


void test_simpleWeightFeatureGeneration(CuTest *testCase) {

    char *featureType = "simpleWeight";
    char *outputName = "test.simpleWeightFeature";
    st_logInfo("\n\nTesting %s feature polishing on %s and %s\n", featureType, FEATURE_TEST_BAM, FEATURE_TEST_FA);
    int64_t retCode = polishingFeatureTest(FEATURE_TEST_BAM, FEATURE_TEST_FA, FEATURE_TEST_NO_RLE_PARAMS,
                                           featureType, NULL, outputName, FALSE);

    char *expectedOutputFa = stString_print("%s.fa", outputName);
    char *expectedOutputFeature = stString_print("%s.T00.h5", outputName);
    CuAssertTrue(testCase, retCode == 0);
    CuAssertTrue(testCase, access(expectedOutputFa, F_OK ) == 0);
    CuAssertTrue(testCase, access(expectedOutputFeature, F_OK ) == 0);
    remove(expectedOutputFa);
    remove(expectedOutputFeature);
}


void test_splitRleWeightFeatureGeneration(CuTest *testCase) {

    char *featureType = "splitRleWeight";
    char *outputName = "test.splitRleWeightFeature";
    st_logInfo("\n\nTesting %s feature polishing on %s and %s\n", featureType, FEATURE_TEST_BAM, FEATURE_TEST_FA);
    int64_t retCode = polishingFeatureTest(FEATURE_TEST_BAM, FEATURE_TEST_FA, FEATURE_TEST_PARAMS, featureType,
            NULL, outputName, FALSE);

    char *expectedOutputFa = stString_print("%s.fa", outputName);
    char *expectedOutputFeature = stString_print("%s.T00.h5", outputName);
    CuAssertTrue(testCase, retCode == 0);
    CuAssertTrue(testCase, access(expectedOutputFa, F_OK ) == 0);
    CuAssertTrue(testCase, access(expectedOutputFeature, F_OK ) == 0);
    remove(expectedOutputFa);
    remove(expectedOutputFeature);
}

void test_simpleWeightIndex(CuTest *testCase) {

    int idx;
    PoaFeatureSimpleWeight *feature = PoaFeature_SimpleWeight_construct(0, 0);

    for (int64_t c = 0; c < SYMBOL_NUMBER; c++) {
        idx = PoaFeature_SimpleWeight_charIndex((Symbol) c, TRUE);
        CuAssertTrue(testCase, idx < POAFEATURE_SIMPLE_WEIGHT_TOTAL_SIZE);
        feature->weights[idx] += 1;

        idx = PoaFeature_SimpleWeight_charIndex((Symbol) c, FALSE);
        CuAssertTrue(testCase, idx < POAFEATURE_SIMPLE_WEIGHT_TOTAL_SIZE);
        feature->weights[idx] += 1;
    }

    idx = PoaFeature_SimpleWeight_gapIndex(TRUE);
    CuAssertTrue(testCase, idx < POAFEATURE_SIMPLE_WEIGHT_TOTAL_SIZE);
    feature->weights[idx] += 1;

    idx = PoaFeature_SimpleWeight_gapIndex(FALSE);
    CuAssertTrue(testCase, idx < POAFEATURE_SIMPLE_WEIGHT_TOTAL_SIZE);
    feature->weights[idx] += 1;

    for (int64_t i = 0; i < POAFEATURE_SIMPLE_WEIGHT_TOTAL_SIZE; i++) {
        CuAssertTrue(testCase, feature->weights[i] == 1);
    }

    PoaFeature_SimpleWeight_destruct(feature);

}


void test_RleWeightIndex(CuTest *testCase) {

    int idx;
    PoaFeatureRleWeight *feature = PoaFeature_RleWeight_construct(0, 0);

    for (int64_t c = 0; c < SYMBOL_NUMBER; c++) {
        for (int64_t l = 1; l <= POAFEATURE_MAX_RUN_LENGTH; l++) {
            idx = PoaFeature_RleWeight_charIndex((Symbol) c, l, TRUE);
            CuAssertTrue(testCase, idx < POAFEATURE_RLE_WEIGHT_TOTAL_SIZE);
            feature->weights[idx] += 1;

            idx = PoaFeature_RleWeight_charIndex((Symbol) c, l, FALSE);
            CuAssertTrue(testCase, idx < POAFEATURE_RLE_WEIGHT_TOTAL_SIZE);
            feature->weights[idx] += 1;
        }
    }

    idx = PoaFeature_RleWeight_gapIndex(TRUE);
    CuAssertTrue(testCase, idx < POAFEATURE_RLE_WEIGHT_TOTAL_SIZE);
    feature->weights[idx] += 1;

    idx = PoaFeature_RleWeight_gapIndex(FALSE);
    CuAssertTrue(testCase, idx < POAFEATURE_RLE_WEIGHT_TOTAL_SIZE);
    feature->weights[idx] += 1;

    for (int64_t i = 0; i < POAFEATURE_RLE_WEIGHT_TOTAL_SIZE; i++) {
        if (feature->weights[i] != 1) {
            CuAssertTrue(testCase, feature->weights[i] == 1);
        }
    }

    PoaFeature_RleWeight_destruct(feature);
}

CuSuite* featureTestSuite(void) {
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, test_simpleWeightIndex);
    SUITE_ADD_TEST(suite, test_RleWeightIndex);
    SUITE_ADD_TEST(suite, test_simpleWeightFeatureGeneration);
    SUITE_ADD_TEST(suite, test_splitRleWeightFeatureGeneration);

    return suite;
}
