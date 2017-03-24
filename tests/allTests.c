/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sonLib.h"

CuSuite *stRPHmmTestSuite(void);

int stMarginPhaseTests(void) {
	CuString *output = CuStringNew();
	CuSuite* suite = CuSuiteNew();
	CuSuiteAddSuite(suite, stRPHmmTestSuite());
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
	int i = stMarginPhaseTests();

	//st_uglyf("Done\n");
	//while(1);

	return i;
}
