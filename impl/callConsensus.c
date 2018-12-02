#include <multipleAligner.h>
#include "margin.h"
#include "randomSequences.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sonLib.h"

static char* paramsPath = "/home/ryan/code/polarbear_assembly/external/polisher_cpp/external/src/project_marginPolish/params/allParams.np.json";
//#define TEST_POLISH_FILES_DIR "../tests/polishTestExamples/"

void testPrint(){
    printf("working\n");
    printf(paramsPath);
    printf("\n");
}

// Make RLEStrings representing reads and list of the RLE strings
char* callConsensus(int readNo, char *readArray[], char *reference) {
    stList *reads = stList_construct();
    stList *rleStrings = stList_construct3(0, (void (*)(void *)) rleString_destruct);

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

    Poa *poaRefined = poa_realignIterative2(reads, rleReference->rleString, p->polishParams, FALSE);

    // Now get a non-RLE (expanded) string
    RleString* rleConsensus = expandRLEConsensus(poaRefined, rleStrings, p->polishParams->repeatSubMatrix);

    char* nonRleString = rleString_expand(rleConsensus);

    return nonRleString;
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
