//#include <multipleAligner.h>
#include "margin.h"


// Make RLEStrings representing reads and list of the RLE strings
PolishParams* getConsensusParameters(char *paramsPath) {

    // Load parameters / models
    Params *p = params_readParams(paramsPath);
    PolishParams *polish = p->polishParams;
    p->polishParams = NULL;
    params_destruct(p);
    return polish;
}


void destroyConsensusParameters(PolishParams *params) {
    polishParams_destruct(params);
}

RleString* callConsensus(int64_t readCount, char *nucleotides[], uint8_t *runLengths[], uint8_t strands[], PolishParams *params) {
    stList *reads = stList_construct3(0, (void (*)(void*)) bamChunkRead_destruct);

    for (int64_t i = 0; i < readCount; i++) {
        RleString *rleString = rleString_constructPreComputed(stString_copy(nucleotides[i]), runLengths[i]);
        char *rawString = rleString_expand(rleString);
        stList_append(reads, bamChunkRead_construct2(stString_print("read_%d", i), rawString,
                NULL, (strands[i] == 0 ? TRUE : FALSE), NULL)); // strands defined as 0 -> forward, 1 -> backward
        free(rawString);
        rleString_destruct(rleString);
    }

    // RLE reference starts as one of the input string
    RleString *rleReference = ((BamChunkRead*)stList_get(reads, 0))->rleRead;

    // run poa
    Poa *poa = poa_realignAll(reads, NULL, rleReference, params);

    // get consensus
    poa_estimateRepeatCountsUsingBayesianModel(poa, reads, params->repeatSubMatrix);
    RleString *polishedRleConsensus = poa->refString;

    //cleanup
    stList_destruct(reads);
    poa_destruct(poa);

    return polishedRleConsensus;
}

void destroyRleString(RleString *r) {
    rleString_destruct(r);
}

