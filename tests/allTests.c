/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "margin.h"

CuSuite *stRPHmmTestSuite(void);
CuSuite* polisherTestSuite(void);
CuSuite* parserTestSuite(void);
CuSuite* viewTestSuite(void);
CuSuite* chunkingTestSuite(void);
CuSuite* callConsensusTestSuite(void);
//CuSuite* featureTestSuite(void);
CuSuite* marginIntegrationTestSuite(void);
CuSuite* pairwiseAlignmentTestSuite(void);

// New tests for marginPhase interface
int marginPhaseTests(void) {

	st_setLogLevel(info);

	CuString *output = CuStringNew();
	CuSuite* suite = CuSuiteNew();

//	CuSuiteAddSuite(suite, stRPHmmTestSuite());
//	CuSuiteAddSuite(suite, parserTestSuite());
    CuSuiteAddSuite(suite, polisherTestSuite()); //todo fails
    CuSuiteAddSuite(suite, marginIntegrationTestSuite()); //todo fails
//	CuSuiteAddSuite(suite, viewTestSuite());
//	CuSuiteAddSuite(suite, pairwiseAlignmentTestSuite());
//	CuSuiteAddSuite(suite, stRPHmmTestSuite());
//	CuSuiteAddSuite(suite, chunkingTestSuite());
//	CuSuiteAddSuite(suite, callConsensusTestSuite()); //todo fails
//    #ifdef _HDF5
//    CuSuiteAddSuite(suite, featureTestSuite()); //todo fails
//    #endif

	CuSuiteRun(suite);
	CuSuiteSummary(suite, output);
	CuSuiteDetails(suite, output);
	printf("%s\n", output->buffer);
	int i = suite->failCount > 0;
	CuSuiteDelete(suite);
	CuStringDelete(output);
	return i;
}

int main(int argc, char *argv[]) {
    if(argc == 2) {
        st_setLogLevelFromString(argv[1]);
    }
	int i = marginPhaseTests();

	//while(1);

	return i;
}
