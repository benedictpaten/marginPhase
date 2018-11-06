//
// Created by ryan on 11/2/18.
//

#ifndef MARGINPHASE_CALLCONSENSUS_H
#define MARGINPHASE_CALLCONSENSUS_H

#ifdef __cplusplus
extern "C" {
#endif

void testPrint();

// iteratively construct a consensus sequence using profile HMM/POA
char* callConsensus(int readNo, char *readArray[], char *reference);

#ifdef __cplusplus
}
#endif

#endif //MARGINPHASE_CALLCONSENSUS_H
