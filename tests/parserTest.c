/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "margin.h"


void test_jsmnParsing(CuTest *testCase) {

    char *paramsFile = "../tests/parsingTest.json";

    stRPHmmParameters *params = parseParameters(paramsFile);

    // Check stRPHmmParameters

    // Check haplotype substitution model and error model parsed correctly
    // and that the proper values were set in the model parameters
    double delta = 0.0001;

    // Check remaining parameters parsed correctly
    CuAssertTrue(testCase, params->maxNotSumTransitions);
    CuAssertIntEquals(testCase, params->maxPartitionsInAColumn, 50);
    CuAssertIntEquals(testCase, params->maxCoverageDepth, 64);
    CuAssertIntEquals(testCase, params->minReadCoverageToSupportPhasingBetweenHeterozygousSites, 4);

    // cleanup
    stRPHmmParameters_destruct(params);
}

CuSuite *marginPhaseParserTestSuite(void) {
    st_setLogLevelFromString("debug");
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, test_jsmnParsing);

    return suite;
}
