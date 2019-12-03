//#include <multipleAligner.h>
#include "margin.h"


// Make RLEStrings representing reads and list of the RLE strings
char* callConsensus(int readNo, char *readArray[], char *reference, char *paramsPath) {
	// Load parameters / models
	FILE *paramsFile = fopen(paramsPath, "r");
	if (paramsFile == NULL) {
		printf("Cannot open file '%s'\n", paramsPath);
		return "";
	}

	Params *p = params_readParams(paramsFile);


	stList *reads = stList_construct3(0, (void (*)(void*)) bamChunkRead_destruct);
    bool *readStrandArray = st_calloc(readNo, sizeof(bool));

    for (int64_t i = 0; i < readNo; i++) {
    	stList_append(reads, bamChunkRead_construct2(stString_print("read_%d", i), readArray[i], NULL, TRUE,
    			p->polishParams->useRunLengthEncoding));
    }

    // RLE reference (reference could be randomly chosen read)
    RleString *rleReference = rleString_construct(reference);

	// Build the POA
    Poa *poaRefined = poa_realignAll(reads, NULL, rleReference, p->polishParams);

    // Restimate the repeat counts of the backbone bases in the POA
    poa_estimateRepeatCountsUsingBayesianModel(poaRefined, reads, p->polishParams->repeatSubMatrix);

    char* nonRleString = rleString_expand(poaRefined->refString);

    //cleanup
    stList_destruct(reads);
    rleString_destruct(rleReference);
    poa_destruct(poaRefined);
    // TODO: Cleanup memory!

    return nonRleString;
}

