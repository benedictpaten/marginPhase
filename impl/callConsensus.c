#include <multipleAligner.h>
#include "margin.h"


// Make RLEStrings representing reads and list of the RLE strings
char* callConsensus(int readNo, char *readArray[], char *reference, char *paramsPath) {
    stList *reads = stList_construct3(0, (void (*)(void*)) bamChunkRead_destruct);
    stList *rleStrings = stList_construct3(0, (void (*)(void *)) rleString_destruct);
    bool *readStrandArray = st_calloc(readNo, sizeof(bool));

    for (int64_t i = 0; i < readNo; i++) {
        RleString *rleString = rleString_construct((char *) stString_copy(readArray[i]));
        stList_append(rleStrings, rleString);
        stList_append(reads, bamChunkRead_construct2(stString_print("read_%d", i), stString_copy(rleString->rleString),
                NULL, TRUE, NULL));
    }

    // RLE reference (reference could be randomly chosen read)
    RleString *rleReference = rleString_construct(stString_copy(reference));

    // Load parameters / models
    Params *p = params_readParams(paramsPath);

    Poa *poaRefined = poa_realignIterative(reads, NULL, rleReference->rleString, p->polishParams);

    // Now get a non-RLE (expanded) string
    RleString* rleConsensus = expandRLEConsensus(poaRefined, rleStrings, reads, p->polishParams->repeatSubMatrix);

    char* nonRleString = rleString_expand(rleConsensus);
    char *nonRLEConsensusString = rleString_expand(rleConsensus);

    //cleanup
    stList_destruct(rleStrings);
    stList_destruct(reads);
    rleString_destruct(rleReference);
    // TODO: Cleanup memory!

    return nonRleString;
}

