/*
 * Copyright (C) 2009-2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <htsIntegration.h>
#include "CuTest.h"
#include "margin.h"
#include "helenFeatures.h"
#include <stdio.h>

static char *FEATURE_TEST_PARAMS = "../params/allParams.np.json";
static char *FEATURE_TEST_BAM = "../tests/data/featureTest.bam";
static char *FEATURE_TEST_FA = "../tests/data/featureTest.fa";
static char *FEATURE_TEST_TRUTH_BAM = "../tests/data/featureTestTruth.bam";
static char *FEATURE_TEST_TRUTH_SEQ = "ACGATAACCGGTTAAACATCCCGGGTTTCAAACCCCGGGGTTGATTACACAT";

stList *getSimpleWeightFeatureFromTSV(CuTest *testCase, char *tsvFile, bool includesLabels) {
    stList *features = stList_construct3(0, (void*)PoaFeature_SimpleCharacterCount_destruct);

    FILE *fp = fopen(tsvFile, "r");
    char * line = NULL;
    size_t len = 0;

    while (getline(&line, &len, fp) != -1) {
        if (line[0] == '#') {
            continue;
        }
        stList *parts = stString_splitByString(line, "\t");
        int64_t totalSize = SYMBOL_NUMBER * 2 + 2 + (includesLabels ? 1 : 0);
        CuAssertTrue(testCase, stList_length(parts) == totalSize);
        int64_t refPos = stSafeStrToInt64(stList_get(parts, 0));
        int64_t insPos = stSafeStrToInt64(stList_get(parts, 1));

        PoaFeatureSimpleWeight *feature = PoaFeature_SimpleWeight_construct(refPos, insPos);
        stList_append(features, feature);
        if (includesLabels) feature->label = ((char*)stList_get(parts, 2))[0];

        // get the character spots
        int64_t featureWeightPos = 0;
        int64_t i = (includesLabels ? 3 : 2);
        for (; i < totalSize - 2; i++) {
            double weight = atof(stList_get(parts, i));
            feature->weights[featureWeightPos] = weight;
            featureWeightPos++;
        }
        feature->weights[featureWeightPos+2] = atof(stList_get(parts, i));
        feature->weights[featureWeightPos+3] = atof(stList_get(parts, i+1));

        stList_destruct(parts);
    }

    free(line);
    fclose(fp);
    return features;
}


int64_t polishingFeatureTest(char *bamFile, char *referenceFile, char *paramsFile, char *featureType,
        char *featureTruth, char *outputName, bool verbose) {

    // Run margin phase
    char *logString = verbose ? "--logLevel DEBUG" : "--logLevel INFO";
    char *featureTruthCmd = featureTruth == NULL ? stString_copy("") : stString_print("--trueReferenceBam %s", featureTruth) ;
    char *command = stString_print("./marginPolish %s %s %s --outputBase %s %s --outputFeatureType %s %s",
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
    st_logInfo("\n\nTesting %s feature polishing on %s and %s with %s\n", featureType, FEATURE_TEST_BAM, FEATURE_TEST_FA, FEATURE_TEST_TRUTH_BAM);
    int64_t retCode = polishingFeatureTest(FEATURE_TEST_BAM, FEATURE_TEST_FA, FEATURE_TEST_PARAMS, featureType, FEATURE_TEST_TRUTH_BAM, outputName, FALSE);

    char *expectedOutputFa = stString_print("%s.fa", outputName);
    char *expectedOutputFeature = stString_print("%s.simpleWeight.C00000.feature_contig-0-50.tsv", outputName);
    CuAssertTrue(testCase, retCode == 0);
    CuAssertTrue(testCase, access(expectedOutputFa, F_OK ) == 0);
    CuAssertTrue(testCase, access(expectedOutputFeature, F_OK ) == 0);

    stList *features = getSimpleWeightFeatureFromTSV(testCase, expectedOutputFeature, TRUE);
    CuAssertTrue(testCase, stList_length(features) > 0);

    char *truthSeq = st_calloc(strlen(FEATURE_TEST_TRUTH_SEQ) + 1, sizeof(char));
    int64_t truthIdx = 0;
    double totalFwdWeight = 0.0;
    double totalBkwdWeight = 0.0;
    for (int64_t i = 0; i < stList_length(features); i++) {
        PoaFeatureSimpleWeight *feature = stList_get(features, i);

        // assertions about truth sequence presence
        if (feature->label != '_') {
            truthSeq[truthIdx] = feature->label;
            truthIdx++;
        }

        // total weights
        double totalPosFwd = 0.0;
        double totalPosBkwd = 0.0;
        for (int64_t j = 0; j < SYMBOL_NUMBER * 2; j++) {
            if (j % 2 == 0) totalPosFwd += feature->weights[j];
            else totalPosBkwd += feature->weights[j];
        }
        double posRatio = totalPosFwd / totalPosBkwd;
        double posTotal = totalPosFwd + totalPosBkwd;
        totalFwdWeight += totalPosFwd;
        totalBkwdWeight += totalPosBkwd;

        if (feature->insertPosition == 0) {
            CuAssertTrue(testCase, posRatio > 1.5 && posRatio < 3.0);
            CuAssertTrue(testCase, posTotal > 4.9 && posTotal < 7.0);
        }
    }

    double totalRatio = totalFwdWeight / totalBkwdWeight;
    double totalWeightPerPos = (totalFwdWeight + totalBkwdWeight) / 50;

    CuAssertTrue(testCase, stString_eq(truthSeq, FEATURE_TEST_TRUTH_SEQ));
    CuAssertTrue(testCase, totalRatio > 4.9 / 2.0 && totalRatio < 5.1 / 2.0);
    CuAssertTrue(testCase, totalWeightPerPos > 6.9 && totalWeightPerPos < 7.1);
}

CuSuite* featureTestSuite(void) {
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, test_simpleWeightFeatureGeneration);

    return suite;
}
