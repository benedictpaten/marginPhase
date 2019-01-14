//
// Created by ryan on 11/20/18.
//


#ifndef MARGINPHASE_CALLCONSENSUS_H
#define MARGINPHASE_CALLCONSENSUS_H

#include "margin.h"

#ifdef __cplusplus
extern "C" {
#endif

    // polish parameters functions
    PolishParams* getConsensusParameters(char *paramsPath);
    void destroyConsensusParameters(PolishParams *params);

    // consensus calling function
    char* callConsensus(int64_t readCount, char *nucleotides[], uint8_t *runLengths[], uint8_t strands[], PolishParams *params);

#ifdef __cplusplus
}
#endif

#endif //MARGINPHASE_CALLCONSENSUS_H
