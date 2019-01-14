#include <multipleAligner.h>
#include "margin.h"


// Make RLEStrings representing reads and list of the RLE strings
PolishParams* getConsensusParameters(char *paramsPath) {

    // Load parameters / models
    FILE *paramsFile = fopen(paramsPath, "r");
    if (paramsFile == NULL) {
        fprintf(stderr, "[getConsensusParameters] cannot open file '%s'\n", paramsPath);
        return NULL;
    }

    Params *p = params_readParams(paramsFile);
    PolishParams *polish = p->polishParams;
    p->polishParams = NULL;
    params_destruct(p);
    return polish;
}

void destroyConsensusParameters(PolishParams *params) {
    polishParams_destruct(params);
}

char* callConsensus(int64_t readCount, char *nucleotides[], uint8_t *runLengths[], uint8_t strands[], PolishParams *params) {
    stList *rleReads = stList_construct3(0, (void (*)(void*)) bamChunkRead_destruct);
    stList *rleStrings = stList_construct3(0, (void (*)(void *)) rleString_destruct);

    for (int64_t i = 0; i < readCount; i++) {
        RleString *rleString = rleString_construct2(stString_copy(nucleotides[i]), runLengths[i]);
        stList_append(rleStrings, rleString);
        stList_append(rleReads, bamChunkRead_construct2(stString_print("read_%d", i), stString_copy(rleString->rleString),
                NULL, (strands[i] == 0 ? TRUE : FALSE), NULL)); // strands defined as 0 -> forward, 1 -> backward
    }

    // RLE reference starts as one of the input string
    RleString *rleReference = stList_get(rleStrings, 0);

    // run poa
    Poa *poa = poa_realignAll(rleReads, NULL, rleReference->rleString, params);

    // get consensus
    RleString *consensusRleString = expandRLEConsensus(poa, rleStrings, rleReads, params->repeatSubMatrix);
    char *consensus = stString_copy(consensusRleString->rleString);

    //cleanup
    stList_destruct(rleStrings);
    stList_destruct(rleReads);
    poa_destruct(poa);
    rleString_destruct(consensusRleString);

    return consensus;
}

