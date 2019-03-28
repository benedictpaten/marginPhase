/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com) & Arthur Rand (arand@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "margin.h"
#include "externalIntegration.h"



void writeParamFile(char *outputFilename, stRPHmmParameters *params) {
    // get file
    FILE *fd = fopen(outputFilename, "w");
    if (fd == NULL) {
        st_logCritical("Failed to open output param file '%s'. No file will be written\n", outputFilename);
        return;
    }
    //for whether to print the last comma
    int64_t noCommaIdx = (ALPHABET_SIZE) * (ALPHABET_SIZE) - 1;

    fprintf(fd, "{\n");
    fprintf(fd, "  \"alphabet\" : [ \"Aa\", \"Cc\", \"Gg\", \"Tt\", \"-\" ],\n");
    fprintf(fd, "    \n");
    fprintf(fd, "  \"wildcard\" : \"Nn\",\n");
    fprintf(fd, "    \n");
    fprintf(fd, "  \"haplotypeSubstitutionModel\" : [ \n");
    for (int64_t i = 0; i < ALPHABET_SIZE; i++) {
        fprintf(fd, "    ");
        for (int64_t j = 0; j < ALPHABET_SIZE; j++) {
            int64_t idx = i * ALPHABET_SIZE + j;
            fprintf(fd, "%8f", exp(params->hetSubModelSlow[idx]));
            if (idx != noCommaIdx) fprintf(fd, ", ");
        }
        fprintf(fd, "\n");
    }
    fprintf(fd, "   ],\n");
    fprintf(fd, "    \n");
    fprintf(fd, "  \"readErrorModel\" : [ \n");
    for (int64_t i = 0; i < ALPHABET_SIZE; i++) {
        fprintf(fd, "    ");
        for (int64_t j = 0; j < ALPHABET_SIZE; j++) {
            int64_t idx = i * ALPHABET_SIZE + j;
            fprintf(fd, "%8f", exp(params->readErrorSubModelSlow[idx]));
            if (idx != noCommaIdx) fprintf(fd, ", ");
        }
        fprintf(fd, "\n");
    }
    fprintf(fd, "   ],\n");
    fprintf(fd, "    \n");
    fprintf(fd, "  \"maxNotSumTransitions\" : %s,\n", params->maxNotSumTransitions ? "true" : "false");
    fprintf(fd, "    \n");
    fprintf(fd, "  \"maxPartitionsInAColumn\" : %" PRIi64 ",\n", params->maxPartitionsInAColumn);
    fprintf(fd, "\n");
    fprintf(fd, "  \"maxCoverageDepth\" : %" PRIi64 ",\n", params->maxCoverageDepth);
    fprintf(fd, "\n");
    fprintf(fd, "  \"minReadCoverageToSupportPhasingBetweenHeterozygousSites\" : %" PRIi64 ",\n", params->minReadCoverageToSupportPhasingBetweenHeterozygousSites);
    fprintf(fd, "  \n");
    fprintf(fd, "  \"onDiagonalReadErrorPseudoCount\" : %f,\n", params->onDiagonalReadErrorPseudoCount);
    fprintf(fd, "  \n");
    fprintf(fd, "  \"offDiagonalReadErrorPseudoCount\" : %f,\n", params->offDiagonalReadErrorPseudoCount);
    fprintf(fd, "  \n");
    fprintf(fd, "  \"trainingIterations\" : %" PRIi64 ",\n", params->trainingIterations);
    fprintf(fd, "  \n");
    int64_t verbosityBitstring = 0
        | (params->verboseTruePositives ? LOG_TRUE_POSITIVES : 0);
    fprintf(fd, "  \"verbose\" : %" PRIi64 ",\n", verbosityBitstring);
    fprintf(fd, "}");

    if (fclose(fd) != 0) st_logCritical("Failed to close output param file: %s\n", outputFilename);
}

//
//  Read Haplotype Sequence functions
//

stReadHaplotypeSequence *stReadHaplotypeSequence_construct(
        int64_t readStart, int64_t phaseBlock, int64_t length, int8_t haplotype) {
    stReadHaplotypeSequence *s = malloc(sizeof(stReadHaplotypeSequence));
    s->readStart = readStart;
    s->phaseBlock = phaseBlock;
    s->length = length;
    s->haplotype = haplotype;
    s->next = NULL;
    return s;
}

char *stReadHaplotypeSequence_toString(stReadHaplotypeSequence *rhs) {
    char *curr = stString_print("h%"PRIi8",p%"PRIi64",r%"PRIi64",l%"PRIi64, rhs->haplotype, rhs->phaseBlock,
                                   rhs->readStart, rhs->length);
    if (rhs->next != NULL) {
        char *next = stReadHaplotypeSequence_toString(rhs->next);
        char *tmp = curr;
        curr = stString_print("%s;%s", tmp, next);
        free(tmp);
        free(next);
    }

    return curr;
}

char *stReadHaplotypeSequence_toStringEmpty() {
    return stString_print("h0,p0,r0,l0");
}

void stReadHaplotypeSequence_destruct(stReadHaplotypeSequence *rhs) {
    if (rhs->next != NULL) {
        stReadHaplotypeSequence_destruct(rhs->next);
    }
    free(rhs);
}


//
//  HaplotypePartitionTable
//  Tracks ReadHaplotypeSequence information
//

//todo destroy keys (char*)
stReadHaplotypePartitionTable *stReadHaplotypePartitionTable_construct(int64_t initialSize) {
    return create_hashtable((uint64_t) initialSize, stHash_stringKey, stHash_stringEqualKey,
                            NULL, (void *)stReadHaplotypeSequence_destruct);
}
void stReadHaplotypePartitionTable_add(stReadHaplotypePartitionTable *hpt, char *readName, int64_t readStart,
                                       int64_t phaseBlock, int64_t length, int8_t haplotype) {

    stReadHaplotypeSequence *new = stReadHaplotypeSequence_construct(readStart, phaseBlock, length, haplotype);
    stReadHaplotypeSequence *curr = hashtable_search(hpt, readName);
    if (NULL == curr) {
        hashtable_insert(hpt, readName, new);
    } else {
        stReadHaplotypeSequence *prev;
        do {
            prev = curr;
            if (curr->phaseBlock == phaseBlock) {
                if (haplotype != curr->haplotype) {
                    st_logCritical("\tRead %s found in both haplotypes in phase block %"PRIi64"\n", readName, phaseBlock);
                }

                // don't need to record
                stReadHaplotypeSequence_destruct(new);
                return;
            }
            curr = curr->next;
        } while (curr != NULL);
        prev->next = new;
    }
}
void stReadHaplotypePartitionTable_destruct(stReadHaplotypePartitionTable *hpt) {
    hashtable_destroy(hpt, true, false);
}


//
//  Populating HaplotypePartitionTable from GenomeFragment and HMM
//

void populateReadHaplotypePartitionTable(stReadHaplotypePartitionTable *hpt, stGenomeFragment *gF, stRPHmm *hmm,
                                         stList *path) {
    //todo track all partitioned reads and quit early if examined
    // same for whole GenomeFragment
    int64_t phaseBlock = gF->refStart - 1;
    int64_t phaseBlockEnd = phaseBlock + gF->length;

    // variables for partitions
    char *readName;
    int64_t readStart;
    int64_t length;
    int8_t haplotype;

    // For each cell/column pair
    stRPColumn *column = hmm->firstColumn;
    for(int64_t i=0; i<stList_length(path); i++) {
        stRPCell *cell = stList_get(path, i);

        // Get reads in partitions
        for(int64_t j=0; j<column->depth; j++) {
            stProfileSeq *read = column->seqHeaders[j];
            readName = read->readId;
            readStart = (read->refStart < phaseBlock ? phaseBlock - read->refStart : 0);
            length = (read->refStart + read->length > phaseBlockEnd ?
                      phaseBlockEnd - read->refStart :
                      read->length - readStart);
            haplotype = (int8_t) (seqInHap1(cell->partition, j) ? 1 : 2);

            //todo more sanity check
            assert(length >= 0);
            assert(length <= read->length);
            assert(readStart >= 0);
            assert(readStart <= read->length);

            // save to hpt
            stReadHaplotypePartitionTable_add(hpt, readName, readStart, phaseBlock, length, haplotype);
        }

        // iterate
        if(column->nColumn != NULL) {
            column = column->nColumn->nColumn;
        }
    }

}
