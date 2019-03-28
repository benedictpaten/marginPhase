/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "margin.h"

CuSuite *stRPHmmTestSuite(void);
CuSuite *marginPhaseTestSuite(void);
CuSuite* polisherTestSuite(void);
CuSuite* parserTestSuite(void);
CuSuite* viewTestSuite(void);
CuSuite* chunkingTestSuite(void);
CuSuite* callConsensusTestSuite(void);
CuSuite* featureTestSuite(void);


// New tests for marginPhase interface
int marginPhaseTests(void) {
	CuString *output = CuStringNew();
	CuSuite* suite = CuSuiteNew();

//	CuSuiteAddSuite(suite, stRPHmmTestSuite());
//	CuSuiteAddSuite(suite, parserTestSuite());
//	CuSuiteAddSuite(suite, marginPhaseTestSuite());
	CuSuiteAddSuite(suite, polisherTestSuite());
	CuSuiteAddSuite(suite, viewTestSuite());
	CuSuiteAddSuite(suite, chunkingTestSuite());
	CuSuiteAddSuite(suite, callConsensusTestSuite());
	CuSuiteAddSuite(suite, featureTestSuite());

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
