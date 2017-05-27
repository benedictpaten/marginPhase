/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "stRPHmm.h"

void test_jsmnParsing(CuTest *testCase) {

    char *paramsFile = "../tests/parsingTest.json";

    stBaseMapper *baseMapper = stBaseMapper_construct();
    stRPHmmParameters *params = parseParameters(paramsFile, baseMapper);

    // Check that alphabet was parsed as expected, and that conversions
    // between types of bases work properlu
    CuAssertIntEquals(testCase, baseMapper->size, 5);
    CuAssertStrEquals(testCase, baseMapper->wildcard, "Nn");

    // Check that numerical bases are mapped to characters correctly
    CuAssertIntEquals(testCase, stBaseMapper_getCharForValue(baseMapper, 0), 'A');
    CuAssertIntEquals(testCase, stBaseMapper_getCharForValue(baseMapper, 1), 'C');
    CuAssertIntEquals(testCase, stBaseMapper_getCharForValue(baseMapper, 2), 'G');
    CuAssertIntEquals(testCase, stBaseMapper_getCharForValue(baseMapper, 3), 'T');
    CuAssertIntEquals(testCase, stBaseMapper_getCharForValue(baseMapper, 4), '-');


    // Check that character bases are mapped to numbers correctly
    CuAssertIntEquals(testCase, stBaseMapper_getValueForChar(baseMapper, 'A'), 0);
    CuAssertIntEquals(testCase, stBaseMapper_getValueForChar(baseMapper, 'a'), 0);
    CuAssertIntEquals(testCase, stBaseMapper_getValueForChar(baseMapper, 'C'), 1);
    CuAssertIntEquals(testCase, stBaseMapper_getValueForChar(baseMapper, 'c'), 1);
    CuAssertIntEquals(testCase, stBaseMapper_getValueForChar(baseMapper, 'G'), 2);
    CuAssertIntEquals(testCase, stBaseMapper_getValueForChar(baseMapper, 'g'), 2);
    CuAssertIntEquals(testCase, stBaseMapper_getValueForChar(baseMapper, 'T'), 3);
    CuAssertIntEquals(testCase, stBaseMapper_getValueForChar(baseMapper, 't'), 3);
    CuAssertIntEquals(testCase, stBaseMapper_getValueForChar(baseMapper, '-'), 4);

    // Check stRPHmmParameters

    // Check haplotype substitution model and error model parsed correctly
    // and that the proper values were set in the model parameters
    double delta = 0.0001;
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            if (i < 4 && j < 4) {
                if (i == j) {
                    CuAssertDblEquals(testCase, params->hetSubModelSlow[i*5+j], log(0.998), delta);
                    CuAssertDblEquals(testCase, params->hetSubModel[i*5+j], scaleToLogIntegerSubMatrix(log(0.998)), delta);
                    CuAssertDblEquals(testCase, params->readErrorSubModelSlow[i*5+j], log(0.9), delta);
                    CuAssertDblEquals(testCase, params->readErrorSubModel[i*5+j], scaleToLogIntegerSubMatrix(log(0.9)), delta);
                } else {
                    CuAssertDblEquals(testCase, params->hetSubModelSlow[i*5+j], log(0.000333), delta);
                    CuAssertDblEquals(testCase, params->hetSubModel[i*5+j], scaleToLogIntegerSubMatrix(log(0.000333)), delta);
                    CuAssertDblEquals(testCase, params->readErrorSubModelSlow[i*5+j], log(0.01), delta);
                    CuAssertDblEquals(testCase, params->readErrorSubModel[i*5+j], scaleToLogIntegerSubMatrix(log(0.01)), delta);
                }
            }
            else if (i != j) {
                CuAssertDblEquals(testCase, params->hetSubModelSlow[i*5+j], log(0.001), delta);
                CuAssertDblEquals(testCase, params->hetSubModel[i*5+j], scaleToLogIntegerSubMatrix(log(0.001)), delta);
                CuAssertDblEquals(testCase, params->readErrorSubModelSlow[i*5+j], log(0.07), delta);
                CuAssertDblEquals(testCase, params->readErrorSubModel[i*5+j], scaleToLogIntegerSubMatrix(log(0.07)), delta);
            } else {
                CuAssertDblEquals(testCase, params->hetSubModelSlow[i*5+j], log(0.996), delta);
                CuAssertDblEquals(testCase, params->hetSubModel[i*5+j], scaleToLogIntegerSubMatrix(log(0.996)), delta);
                CuAssertDblEquals(testCase, params->readErrorSubModelSlow[i*5+j], log(0.72), delta);
                CuAssertDblEquals(testCase, params->readErrorSubModel[i*5+j], scaleToLogIntegerSubMatrix(log(0.72)), delta);
            }
        }
    }
    // Check remaining parameters parsed correctly
    CuAssertTrue(testCase, params->maxNotSumTransitions);
    CuAssertIntEquals(testCase, params->maxPartitionsInAColumn, 50);
    CuAssertIntEquals(testCase, params->maxCoverageDepth, 64);
    CuAssertIntEquals(testCase, params->minReadCoverageToSupportPhasingBetweenHeterozygousSites, 4);

    // cleanup
    stBaseMapper_destruct(baseMapper);
    stRPHmmParameters_destruct(params);
}

CuSuite *marginPhaseParserTestSuite(void) {
    st_setLogLevelFromString("debug");
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, test_jsmnParsing);

    return suite;
}
