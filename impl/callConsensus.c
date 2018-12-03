//
// Created by ryan on 11/20/18.
//

//
// Created by ryan on 10/30/18.
//

#include "margin.h"

static char* paramsPath = "../params/polish/polishParams.json";
#define TEST_POLISH_FILES_DIR "../tests/polishTestExamples/"

void testPrint(){
    printf("working\n");
    printf("%s\n", paramsPath);
    printf("\n");
}

// Make RLEStrings representing reads and list of the RLE strings
char* callConsensus(int readNo, char *readArray[], char *reference) {
    stList *reads = stList_construct();
    stList *rleStrings = stList_construct3(0, (void (*)(void *)) rleString_destruct);
    bool *readStrandArray = st_calloc(readNo, sizeof(bool));

    for (int64_t i = 0; i < readNo; i++) {
        RleString *rleString = rleString_construct((char *) readArray[i]);
        stList_append(rleStrings, rleString);
        stList_append(reads, rleString->rleString);
    }

    // RLE reference (reference could be randomly chosen read)
    RleString *rleReference = rleString_construct(reference);

    // Load parameters / models
    FILE *paramsFile = fopen(paramsPath, "r");
    if (paramsFile == NULL) {
        printf("Cannot open file '%s'\n", paramsPath);
     return "";
    }

    Params *p = params_readParams(paramsFile);

    Poa *poaRefined = poa_realignIterative(reads, readStrandArray, NULL, rleReference->rleString, p->polishParams);    // Now get a non-RLE (expanded) string

    RleString *rleConsensusString = expandRLEConsensus(poaRefined, rleStrings, readStrandArray, p->polishParams->repeatSubMatrix);
    char *nonRLEConsensusString = rleString_expand(rleConsensusString);

    // TODO: Cleanup memory!

    return nonRLEConsensusString;
}


//int main(){
//    int readNo;
//
//    readNo = 45;
//    callConsensus(readNo, readArrayExample1, trueReferenceExample1);
//
//    readNo = 42;
//    callConsensus(readNo, readArrayExample2, trueReferenceExample2);
//}
