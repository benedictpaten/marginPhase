/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "stRPHmm.h"

// OpenMP
#if defined(_OPENMP)
#include <omp.h>
#define CELL_BUFFER_SIZE 1000
#endif

inline double logAddP(double a, double b, bool maxNotSum) {
    /*
     * Local function for doing addition of logs or (if doing Viterbi style calculation), to take the max.
     */
    return maxNotSum ? (a > b ? a : b) : stMath_logAdd(a, b);
}

/*
 * Functions for the read partitioning hmm object stRPHmm.
 */

void stRPHmmParameters_destruct(stRPHmmParameters *params) {
    free(params->hetSubModel);
    free(params->hetSubModelSlow);
    free(params->readErrorSubModel);
    free(params->readErrorSubModelSlow);
    free(params);
}

static void printMatrix(FILE *fH, double *matrixSlow, uint16_t *matrixFast) {
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        fprintf(fH, "\t\t");
        for(int64_t j=0; j<ALPHABET_SIZE; j++) {
            fprintf(fH, "(FROM %" PRIi64 ", TO: %" PRIi64 "): %8f (%6i); ", i, j, exp(matrixSlow[i*ALPHABET_SIZE + j]), matrixFast[i*ALPHABET_SIZE + j]);
        }
        fprintf(fH, "\n");
    }
}

double *getColumnBaseComposition(stRPColumn *column, int64_t pos) {
    /*
     * Get the observed counts for each base seen at a particular position in a column
     */
    double *baseCounts = st_calloc(ALPHABET_SIZE, sizeof(double));
    for (int64_t i=0; i<column->depth; i++) {
        stProfileSeq *seq = column->seqHeaders[i];
        int64_t p = pos - seq->refStart;
        if (seq->refCoords[p] >= seq->refStart &&
                seq->refCoords[p] < seq->refEnd) {
            for(int64_t j=0; j<ALPHABET_SIZE; j++) {
                baseCounts[j] += getProb(&(seq->profileProbs[(p) * ALPHABET_SIZE]), j);
            }
        }
    }
    return baseCounts;
}

void printBaseComposition2(double *baseCounts) {
    /*
     * Print the counts/fraction of each alphabet character in a slightly more compressed form.
     */
    double totalCount = 0;
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        totalCount += baseCounts[i];
    }
    st_logDebug("\t\t0 (A)\t1 (C)\t2 (G)\t3 (T)\t4 (-) \n");
    st_logDebug("    Counts:");
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        st_logDebug("\t%0.1f", baseCounts[i]);
    }
    st_logDebug("\n    Fraction:");
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        st_logDebug("\t%0.3f", baseCounts[i]/totalCount);
    }
    st_logDebug( "\n");
}

void printColumnAtPosition(stRPHmm *hmm, int64_t pos) {
    /*
     * Print out the columns of the hmm at a specific position
     */
    stRPColumn *column = hmm->firstColumn;
    while(1) {
        if (pos >= column->refStart && pos <= column->refCoords[column->length-1]) {
            double *columnBaseCounts = getColumnBaseComposition(column, pos);
            printBaseComposition2(columnBaseCounts);
            free(columnBaseCounts);
        }
        if (column->nColumn == NULL) {
            break;
        }
        column = column->nColumn->nColumn;
    }
}

double *getProfileSequenceBaseCompositionAtPosition(stSet *profileSeqs, int64_t pos) {
    /*
     * Get the expected count of each alphabet character in the profile sequences, returned
     * as an array.
     */
    double *baseCounts = st_calloc(ALPHABET_SIZE, sizeof(double));
    stSetIterator *it = stSet_getIterator(profileSeqs);
    stProfileSeq *pSeq;
    while((pSeq = stSet_getNext(it)) != NULL) {
        if (pos > pSeq->refStart && pos < pSeq->refEnd) {
            for(int64_t j=0; j<ALPHABET_SIZE; j++) {
                baseCounts[j] += getProb(&(pSeq->profileProbs[(pos - pSeq->refStart)*ALPHABET_SIZE]), j);
            }
        }
    }
    return baseCounts;
}

void stRPHmmParameters_printParameters(stRPHmmParameters *params, FILE *fH) {
    /*
     * Print the parameters in the parameters object in a human readable form.
     */
    fprintf(fH, "Read Partitioning HMM Parameters\n");
    fprintf(fH, "\tAlphabet_size: %i\n"
            "\tMax_read coverage_depth: %" PRIi64 "\n"
            "\tMax_not sum transitions?: %i\n"
            "\tMax_partitions in a column of an HMM: %" PRIi64 "\n"
            "\tMin read coverage to support phasing between heterozygous sites: %" PRIi64 "\n",
            ALPHABET_SIZE, params->maxCoverageDepth,
            (int)params->maxNotSumTransitions, params->maxPartitionsInAColumn,
            params->minReadCoverageToSupportPhasingBetweenHeterozygousSites);

    fprintf(fH, "\tHeterozygous substitution rates:\n");
    printMatrix(fH, params->hetSubModelSlow, params->hetSubModel);

    fprintf(fH, "\tRead error substitution rates:\n");
    printMatrix(fH, params->readErrorSubModelSlow, params->readErrorSubModel);

    fprintf(fH, "\tIterations of parameter learning: %" PRIi64 "\n", params->trainingIterations);
    fprintf(fH, "\tFilter bad reads?: %i\n", (int)params->filterBadReads);
    fprintf(fH, "\tFilter match threshold: %f\n", params->filterMatchThreshold);
    fprintf(fH, "\tVerbose Attributes:\n");
    if (params->verboseTruePositives) fprintf(fH, "\t\tTRUE_POSITIVES\n");
    if (params->verboseFalsePositives) fprintf(fH, "\t\tFALSE_POSITIVES\n");
}

static void calculateReadErrorSubModel(double *readErrorSubModel, stGenomeFragment *gF, uint64_t *haplotypeSeq, stSet *reads) {
    /*
     * Returns a normalized substitution matrix estimating the probability of read error substitutions by ML.
     */
    stSetIterator *readIt = stSet_getIterator(reads);
    stProfileSeq *pSeq;
    int64_t end = gF->refCoords[gF->length - 1] + 1;
    while((pSeq = stSet_getNext(readIt)) != NULL) {
        // Get the overlapping interval
        int64_t i = gF->refStart > pSeq->refStart ? gF->refStart : pSeq->refStart;
        int64_t j = end < pSeq->refEnd ? end : pSeq->refEnd;
        // For each pair of read and haplotype characters
        for(;i<j;i++) {
            // Check coordinates in bounds
            // FIXME
            int64_t p = findCorrespondingRefCoordIndex(i, gF->refCoords, pSeq->refCoordMap);
//            int64_t p = i - pSeq->refStart;
            assert(i - gF->refStart >= 0 && gF->refCoords[i] < gF->length);
            assert(p >= 0 && p < pSeq->length);
            // Does that index work?
            int64_t hapChar = haplotypeSeq[i - gF->refStart];
            for(int64_t readChar=0; readChar<ALPHABET_SIZE; readChar++) {
                double probOfReadChar = getProb(&(pSeq->profileProbs[(p) * ALPHABET_SIZE]), readChar);
                *getSubstitutionProbSlow(readErrorSubModel, hapChar, readChar) += probOfReadChar;
            }
        }
    }
    stSet_destructIterator(readIt);
}

static void normaliseSubstitutionMatrix(double *subMatrix) {
    /*
     * Normalise matrix so that counts are converted to conditional probabilities of observing
     * derived character given source character.
     */
    for(int64_t fromChar=0; fromChar<ALPHABET_SIZE; fromChar++) {
        double totalSubCount = 0.0;
        for(int64_t toChar=0; toChar<ALPHABET_SIZE; toChar++) {
            totalSubCount += *getSubstitutionProbSlow(subMatrix, fromChar, toChar);
        }
        for(int64_t toChar=0; toChar<ALPHABET_SIZE; toChar++) {
            double p = *getSubstitutionProbSlow(subMatrix, fromChar, toChar) / totalSubCount;
            *getSubstitutionProbSlow(subMatrix, fromChar, toChar) = p <= 0.0001 ? 0.0001 : p;
        }
    }
}

void stRPHmmParameters_learnParameters(stRPHmmParameters *params, stList *profileSequences,
        stHash *referenceNamesToReferencePriors) {
    /*
     * Learn the substitution matrices iteratively, updating the params object in place. Iterations is the number of cycles
     * of stochastic parameter search to do.
     */

    // For each iteration construct a set of HMMs and estimate the parameters from it.
    for(int64_t i=0; i<params->trainingIterations; i++) {
        st_logDebug("\tStarting training iteration %" PRIi64 "\n", i);
        // Substitution model for haplotypes to reads
        double *readErrorSubModel = st_calloc(ALPHABET_SIZE * ALPHABET_SIZE, sizeof(double));
        for(int64_t j=0; j<ALPHABET_SIZE*ALPHABET_SIZE; j++) {
            readErrorSubModel[j] = params->offDiagonalReadErrorPseudoCount;
        }
        for(int64_t j=0; j<ALPHABET_SIZE; j++) {
            readErrorSubModel[j*ALPHABET_SIZE + j] = params->onDiagonalReadErrorPseudoCount;
        }

        stList *hmms = getRPHmms(profileSequences, referenceNamesToReferencePriors, params);

        st_logDebug("\t\thandling %d HMMs\n", stList_length(hmms));
        for(int64_t j=0; j<stList_length(hmms); j++) {
            st_logDebug("\t\t\titer %" PRIi64 "\n", j);
            stRPHmm *hmm = stList_get(hmms, j);

            // Run the forward-backward algorithm
            st_logDebug("\t\t\tforward backward\n", j);
            stRPHmm_forwardBackward(hmm);

            // Now compute a high probability path through the hmm
            st_logDebug("\t\t\tpath\n", j);
            stList *path = stRPHmm_forwardTraceBack(hmm);

            // Compute the genome fragment
            st_logDebug("\t\t\tget fragment\n", j);
            stGenomeFragment *gF = stGenomeFragment_construct(hmm, path);

            // Get partitioned sequences
            st_logDebug("\t\t\tpartition\n", j);
            stSet *reads1 = stRPHmm_partitionSequencesByStatePath(hmm, path, 1);
            stSet *reads2 = stRPHmm_partitionSequencesByStatePath(hmm, path, 0);

            // Estimate the read error substitution parameters
            st_logDebug("\t\t\tread error sub model\n", j);
            calculateReadErrorSubModel(readErrorSubModel, gF, gF->haplotypeString1, reads1);
            calculateReadErrorSubModel(readErrorSubModel, gF, gF->haplotypeString2, reads2);

            // Cleanup
            st_logDebug("\t\t\tcleanup\n", j);
            stSet_destruct(reads1);
            stSet_destruct(reads2);
            stGenomeFragment_destruct(gF);
            stList_destruct(path);
        }


        // Cleanup
        st_logDebug("\t\tcleaning\n");
        //TODO I think we need to destruct each hmm in here too
        stList_destruct(hmms);

        // Normalise the probabilities
        st_logDebug("\t\tnormalizing\n");
        normaliseSubstitutionMatrix(readErrorSubModel);

        // Update the read error substitution parameters of the parameters object
        st_logDebug("\t\tupdating sub prob\n");
        for(int64_t j=0; j<ALPHABET_SIZE; j++) {
            for(int64_t k=0; k<ALPHABET_SIZE; k++) {
                setSubstitutionProb(params->readErrorSubModel, params->readErrorSubModelSlow, j, k,
                        *getSubstitutionProbSlow(readErrorSubModel, j, k));
            }
        }

        // Cleanup
        st_logDebug("\t\tcleaning\n");
        free(readErrorSubModel);

        //Log the parameters info
        if(st_getLogLevel() == debug) {
            st_logDebug("\tParameters learned after iteration %" PRIi64 " of training:\n", i);
            stRPHmmParameters_printParameters(params, stderr);
        }
    }
}

static int cmpint64(int64_t i, int64_t j) {
    return i > j ? 1 : i < j ? -1 : 0;
}

inline int stRPHmm_cmpFn(const void *a, const void *b) {
    /*
     * Compares two read partitioning HMMs by coordinate on the reference.
     * Will return equal only if they are the same HMM, with the same memory
     * address, otherwise compares pointers for equal HMMs.
     */
    stRPHmm *hmm1 = (stRPHmm *)a, *hmm2 = (stRPHmm *)b;
    int i = strcmp(hmm1->referenceName, hmm2->referenceName);
    if(i == 0) {
        i = cmpint64(hmm1->refStart,  hmm2->refStart);
        if(i == 0) {
            i = cmpint64(hmm1->refLength,  hmm2->refLength);
            if(i == 0) {
                i = hmm1 > hmm2 ? 1 : (hmm1 < hmm2 ? -1 : 0);
            }
        }
    }
    return i;
}

stRPHmm *stRPHmm_construct(stProfileSeq *profileSeq, stReferencePriorProbs *referencePriorProbs, stRPHmmParameters *params) {
    /*
     * Create a read partitioning HMM representing the single sequence profile.
     */

    stRPHmm *hmm = st_calloc(1, sizeof(stRPHmm));

    //  Set reference coordinates
    hmm->referenceName = stString_copy(profileSeq->referenceName);
    hmm->refStart = profileSeq->refStart;
    hmm->refLength = profileSeq->length;

    // Add the single profile sequence to the list of the hmm's sequences
    hmm->profileSeqs = stList_construct();
    stList_append(hmm->profileSeqs, profileSeq);

    hmm->parameters = params; // Parameters for the model for computation, this is shared by different HMMs

    hmm->referencePriorProbs = referencePriorProbs;
    assert(stString_eq(hmm->referenceName, referencePriorProbs->referenceName));
    assert(hmm->refStart >= referencePriorProbs->refStart);


    // Set reference coordinates from those in the profile sequence
    hmm->refCoords = st_calloc(hmm->refLength, sizeof(int64_t));
    hmm->refCoordMap = stHash_construct3(stHash_stringKey, stHash_intPtrEqualKey, NULL, NULL);
    int64_t *indexes = st_calloc(hmm->refLength, sizeof(int64_t));
    for (int64_t i = 0; i < profileSeq->length; i++) {
        hmm->refCoords[i] = profileSeq->refCoords[i];
        indexes[i] = i;
        if (stHash_search(hmm->refCoordMap, &hmm->refCoords[i]) == NULL) {
            stHash_insert(hmm->refCoordMap, &hmm->refCoords[i], &indexes[i]);
        }
    }


    assert(hmm->refCoords[hmm->refLength - 1] <= referencePriorProbs->refCoords[referencePriorProbs->length - 1]);

    hmm->columnNumber = 1; // The number of columns in the model, initially just 1
    hmm->maxDepth = 1; // The maximum number of states in a column, initially just 1

    // Create the first column of the model
    stProfileSeq **seqHeaders = st_malloc(sizeof(stProfileSeq *));
    seqHeaders[0] = profileSeq;
    uint8_t **seqs = st_malloc(sizeof(uint8_t *));
    seqs[0] = profileSeq->profileProbs;
    stRPColumn *column = stRPColumn_construct(hmm->refStart, hmm->refLength, 1, seqHeaders, seqs);
    int64_t *colIndexes = st_calloc(column->length, sizeof(int64_t));
    for (int64_t i = 0; i < hmm->refLength; i++) {
        column->refCoords[i] = hmm->refCoords[i];
        colIndexes[i] = i;
        if (stHash_search(column->refCoordMap, &column->refCoords[i]) == NULL) {
            stHash_insert(column->refCoordMap, &column->refCoords[i], &colIndexes[i]);
        }
    }
    hmm->firstColumn = column;
    hmm->lastColumn = column;

    // Add two cells to the column to represent the two possible partitions of the single profile sequence
    stRPCell *cell = stRPCell_construct(1);
    column->head = cell;
    cell->nCell = stRPCell_construct(0);

    return hmm;
}

void stRPHmm_destruct(stRPHmm *hmm, bool destructColumns) {
    /*
     * Free memory owned by the hmm, including columns.
     */
    free(hmm->referenceName);
    stList_destruct(hmm->profileSeqs);

    if(destructColumns) {
        // Cleanup the columns of the hmm
        stRPColumn *column = hmm->firstColumn;
        while(1) {
            stRPMergeColumn *mColumn = column->nColumn;
            stRPColumn_destruct(column);
            if(mColumn == NULL) {
                break;
            }
            column = mColumn->nColumn;
            stRPMergeColumn_destruct(mColumn);
        }
        // Only free if columns are being destroyed?
//        free(hmm->refCoords);
    }
    // FIXME actually destroy
    //    stHash_destruct(hmm->refCoordMap);

    free(hmm);
}

void stRPHmm_destruct2(stRPHmm *hmm) {
    /*
     * Cleans up hmm and columns
     */
    stRPHmm_destruct(hmm, 1);
}

stList *stRPHmm_forwardTraceBack(stRPHmm *hmm) {
    /*
     * Traces back through the forward matrix picking the most probable path.
     * (yes, this is non-symmetric)
     * Returns the result as a list of cells, one from each column.
     */
    stList *path = stList_construct();

    stRPColumn *column = hmm->lastColumn;

    // Pick cell in the last column with highest probability
    stRPCell *cell = column->head;
    double maxProb = cell->forwardLogProb;
    stRPCell *maxCell = cell;
    while((cell = cell->nCell) != NULL) {
        if(cell->forwardLogProb > maxProb) {
            maxProb = cell->forwardLogProb;
            maxCell = cell;
        }
    }

    stList_append(path, maxCell); // Add chosen cell to output

    // Walk back through previous columns
    while(column->pColumn != NULL) {
        // Get previous merge cell
        stRPMergeCell *mCell = stRPMergeColumn_getPreviousMergeCell(maxCell, column->pColumn);
        assert(mCell != NULL);

        // Switch to previous column
        column = column->pColumn->pColumn;

        // Walk through cells in the previous column to find the one with the
        // highest forward probability that transitions to maxCell
        cell = column->head;
        maxCell = NULL;
        maxProb = ST_MATH_LOG_ZERO;
        do {
            // If compatible and has greater probability
            if(stRPMergeColumn_getNextMergeCell(cell, column->nColumn) == mCell && cell->forwardLogProb > maxProb) {
                maxProb = cell->forwardLogProb;
                maxCell = cell;
            }
        } while((cell = cell->nCell) != NULL);

        assert(maxCell != NULL);
        stList_append(path, maxCell);
    }

    stList_reverse(path); // So cells go in order

    return path;
}

stSet *stRPHmm_partitionSequencesByStatePath(stRPHmm *hmm, stList *path, bool partition1) {
    /*
     * For an hmm and path through the hmm (e.g. computed with stRPHmm_forwardTraceBack) returns the
     * set of sequences in the hmm that are predicted to come from one given haplotype.
     */

    stSet *seqsInHap1 = stSet_construct();

    // For each cell/column pair
    stRPColumn *column = hmm->firstColumn;
    for(int64_t i=0; i<stList_length(path); i++) {
        stRPCell *cell = stList_get(path, i);

        // Get sequences in first or second partition
        for(int64_t j=0; j<column->depth; j++) {
            if((seqInHap1(cell->partition, j) && partition1) ||
                    (!seqInHap1(cell->partition, j) && !partition1)) {
                stSet_insert(seqsInHap1, column->seqHeaders[j]);
            }
        }

        if(column->nColumn != NULL) {
            column = column->nColumn->nColumn;
        }
    }

    return seqsInHap1;
}

void stRPHmm_print(stRPHmm *hmm, FILE *fileHandle, bool includeColumns, bool includeCells) {
    /*
     * Prints a debug friendly representation of the state of an hmm.
     */
    //Header line
    fprintf(fileHandle, "HMM REF_NAME: %s REF_START: %" PRIi64 " REF_LENGTH %" PRIi64
            " COLUMN_NUMBER %" PRIi64 " MAX_DEPTH: %" PRIi64 " FORWARD_PROB: %f BACKWARD_PROB: %f\n",
            hmm->referenceName, hmm->refStart, hmm->refLength,
            hmm->columnNumber, hmm->maxDepth,
            (float)hmm->forwardLogProb, (float)hmm->backwardLogProb);

    if(includeColumns) {
        stRPColumn *column = hmm->firstColumn;
        int64_t i=0;
        while(1) {
            fprintf(fileHandle, "Column %" PRIi64 "\n", i++);

            // Print the column
            stRPColumn_print(column, fileHandle, includeCells);

            if(column->nColumn == NULL) {
                break;
            }

            // Print the merge column
            stRPMergeColumn_print(column->nColumn, fileHandle, includeCells);

            column = column->nColumn->nColumn;
        }
    }
}

stRPHmm *stRPHmm_fuse(stRPHmm *leftHmm, stRPHmm *rightHmm) {
    /*
     * Fuses together two hmms, such that leftHmm and rightHMM are on the same reference sequence and non-overlapping and
     * left hmm precedes right hmm on the reference sequence.
     * Returns fused hmm, destroys input hmms in the process.
     */

    // Checks
    if(!stString_eq(leftHmm->referenceName, rightHmm->referenceName)) {
        st_errAbort("Attempting to fuse two hmms not on the same reference sequence");
    }

    if(stRPHmm_overlapOnReference(leftHmm, rightHmm)) {
        st_errAbort("Attemping to fuse two hmms that overlap in reference coordinates");
    }

    if(leftHmm->refStart >= rightHmm->refStart) {
        st_errAbort("Left hmm does not precede right hmm in reference coordinates for merge");
    }

    // Create a new empty hmm
    stRPHmm *hmm = st_malloc(sizeof(stRPHmm));
    // Set the reference interval
    hmm->referenceName = stString_copy(leftHmm->referenceName);
    hmm->refStart = leftHmm->refStart;

    // TODO what about insertions in the gap?
    int64_t gapLeftIndex = findCorrespondingRefCoordIndex(leftHmm->refLength - 1, leftHmm->refCoords, leftHmm->referencePriorProbs->refCoordMap);
    int64_t gapRightIndex = findCorrespondingRefCoordIndex(0, rightHmm->refCoords, rightHmm->referencePriorProbs->refCoordMap);
    int64_t gapLength = gapRightIndex - gapLeftIndex - 1;
//    int64_t gapLength = rightHmm->refStart - (leftHmm->refCoords[leftHmm->refLength - 1]) - 1;
    if (gapLength < 0) {
        st_logInfo("gapLeftIndex: %d   gapRightIndex: %d   gapLength: %d \n", gapLeftIndex, gapRightIndex, gapLength);
        st_logInfo("leftHmm: %d - %d\tgap: %d - %d (%d - %d, %d)\trightHmm: %d - %d\n", leftHmm->refStart, leftHmm->refCoords[leftHmm->refLength-1], leftHmm->referencePriorProbs->refCoords[gapLeftIndex], leftHmm->referencePriorProbs->refCoords[gapRightIndex], gapLeftIndex, gapRightIndex, gapLength, rightHmm->refStart, rightHmm->refCoords[rightHmm->refLength-1]);
    }
    assert(gapLength >= 0);
//    hmm->refLength = rightHmm->refStart + rightHmm->refLength - leftHmm->refStart;
    hmm->refLength = leftHmm->refLength + rightHmm->refLength + gapLength;

    // Create the combined list of profile seqs
    hmm->profileSeqs = stList_copy(leftHmm->profileSeqs, NULL);
    stList_appendAll(hmm->profileSeqs, rightHmm->profileSeqs);
    // Set column number
    hmm->columnNumber = leftHmm->columnNumber + rightHmm->columnNumber;
    // Max depth
    hmm->maxDepth = leftHmm->maxDepth > rightHmm->maxDepth ? leftHmm->maxDepth : rightHmm->maxDepth;
    // Parameters
    if(leftHmm->parameters != rightHmm->parameters) {
        st_errAbort("HMM parameters differ in fuse function, panic.");
    }
    hmm->parameters = leftHmm->parameters;

    // Set reference position prior probabilities
    if(leftHmm->referencePriorProbs != rightHmm->referencePriorProbs) {
        st_errAbort("Hmm reference prior probs differ in fuse function, panic.");
    }
    hmm->referencePriorProbs = leftHmm->referencePriorProbs;

    // Make reference coords from those in left & right hmms
    // TODO is this indexing correctly?
    hmm->refCoords = st_calloc(hmm->refLength, sizeof(int64_t));
    hmm->refCoordMap = stHash_construct3(stHash_stringKey, stHash_intPtrEqualKey, NULL, NULL);
    int64_t *indexes = st_calloc(hmm->refLength, sizeof(int64_t));
    for (int64_t i = 0; i < leftHmm->refLength; i++) {
        hmm->refCoords[i] = leftHmm->refCoords[i];
        indexes[i] = i;
//        if (i == 8953) {
//            st_logInfo("Left: index = %d  refCoord = %d\n", i, hmm->refCoords[i]);
//        }
        if (stHash_search(hmm->refCoordMap, &hmm->refCoords[i]) == NULL) {
            stHash_insert(hmm->refCoordMap, &hmm->refCoords[i], &indexes[i]);
        }
//        if (leftHmm->refCoords[i] == 8117648 || leftHmm->refCoords[i] == 8176169) {
//            st_logInfo("\tin left hmm, refCoords[%d + stuff] = %d\n", i, leftHmm->refCoords[i]);
//        }
    }
    for (int64_t i = 0; i < rightHmm->refLength; i++) {
        int64_t rightIndex = i + leftHmm->refLength + gapLength;
        hmm->refCoords[rightIndex] = rightHmm->refCoords[i];
//        if (rightIndex == 8953) {
//            st_logInfo("Right: index = %d  refCoord = %d\n", rightIndex, hmm->refCoords[rightIndex]);
//        }
        indexes[rightIndex] = rightIndex;
        if (stHash_search(hmm->refCoordMap, &hmm->refCoords[rightIndex]) == NULL) {
            stHash_insert(hmm->refCoordMap, &hmm->refCoords[rightIndex], &indexes[rightIndex]);
        }
//        if (rightHmm->refCoords[i] == 8117648 || rightHmm->refCoords[i] == 8176169) {
//            st_logInfo("\tin right hmm, refCoords[%d + stuff] = %d\n", i, rightHmm->refCoords[i]);
//        }
    }
    // Make columns to fuse left hmm and right hmm's columns
    stRPMergeColumn *mColumn = stRPMergeColumn_construct(0, 0);

    // Links
    leftHmm->lastColumn->nColumn = mColumn;
    mColumn->pColumn = leftHmm->lastColumn;

    // Add merge cell to connect the cells in the two columns
    stRPMergeCell_construct(0, 0, mColumn);

    if(gapLength > 0) {
        // Make column in the gap
        stRPColumn *column = stRPColumn_construct(leftHmm->refCoords[leftHmm->refLength - 1] + 1,
                gapLength, 0, NULL, NULL);
        int64_t  *colIndexes = st_calloc(column->length, sizeof(int64_t));
        // TODO: what about refCoords in the gap?
        for (int64_t i = 0; i < gapLength; i++) {
            int64_t rProbsIndex = i + gapLeftIndex;
            int64_t refCoord = leftHmm->referencePriorProbs->refCoords[rProbsIndex]+1;
//            int64_t refCoord = leftHmm->referencePriorProbs->refCoords[i + gapLeftIndex+1];
            int64_t hmmIndex = i + leftHmm->refLength;
            hmm->refCoords[hmmIndex] = refCoord;
//            if (hmmIndex == 8953) {
//                st_logInfo("Gap: index = %d  refCoord = %d\n", hmmIndex, refCoord);
//            }
            indexes[hmmIndex] = hmmIndex;
            if (stHash_search(hmm->refCoordMap, &hmm->refCoords[hmmIndex]) == NULL) {
                stHash_insert(hmm->refCoordMap, &hmm->refCoords[hmmIndex], &indexes[hmmIndex]);
            }
//            hmm->refCoords[i + leftHmm->refLength] = i + leftHmm->refCoords[leftHmm->refLength - 1] + 1;
            column->refCoords[i] = refCoord;
            colIndexes[i] = i;
            if (stHash_search(column->refCoordMap, &column->refCoords[i]) == NULL) {
                stHash_insert(column->refCoordMap, &column->refCoords[i], &colIndexes[i]);
            }
//            if (refCoord == 8117648 || refCoord == 8176169) {
//                st_logInfo("\tin gap, refCoord = %d at index %d  col->refStart: %d (%d)\n", refCoord, i, column->refStart, column->refCoords[0]);
//                st_logInfo("rProbs[i-1] = %d, rProbs[i] = %d, rProbs[i+1] = %d, rProbs[i+2] = %d\n", leftHmm->referencePriorProbs->refCoords[i + gapLeftIndex], refCoord, leftHmm->referencePriorProbs->refCoords[i + gapLeftIndex + 2], leftHmm->referencePriorProbs->refCoords[i + gapLeftIndex + 3]);
//            }

        }

        // Links
        mColumn->nColumn = column;
        column->pColumn = mColumn;
        // Make cell for empty column
        column->head = stRPCell_construct(0);

        // Add right merge column
        mColumn = stRPMergeColumn_construct(0, 0);

        // Add merge cell to connect the cells in the two columns
        stRPMergeCell_construct(0, 0, mColumn);

        // Links
        column->nColumn = mColumn;
        mColumn->pColumn = column;
        // Increase the column number to account for the introduced gap column
        hmm->columnNumber += 1;
    }
    mColumn->nColumn = rightHmm->firstColumn;
    rightHmm->firstColumn->pColumn = mColumn;

    // Initialise first/last columns of fused hmm
    hmm->firstColumn = leftHmm->firstColumn;
    hmm->lastColumn = rightHmm->lastColumn;

//    if ((hmm->refStart <= 8117648 && hmm->refCoords[hmm->refLength-1] >= 8117648) ||
//            (hmm->refStart <= 8176169 && hmm->refCoords[hmm->refLength-1] >= 8176169)){
//        st_logInfo("leftHmm: %d - %d\tgap: %d - %d (%d - %d, %d)\trightHmm: %d - %d\n", leftHmm->refStart, leftHmm->refCoords[leftHmm->refLength-1], leftHmm->referencePriorProbs->refCoords[gapLeftIndex], leftHmm->referencePriorProbs->refCoords[gapRightIndex], gapLeftIndex, gapRightIndex, gapLength, rightHmm->refStart, rightHmm->refCoords[rightHmm->refLength-1]);
//        st_logInfo("rProbs[leftGap-1] = %d   leftHmm[leftGap-1] = %d \n", leftHmm->referencePriorProbs->refCoords[gapLeftIndex - 1], leftHmm->refCoords[leftHmm->refLength-2]);
//    }
    for (int64_t i = 0; i < hmm->refLength; i++) {
//        if (hmm->refCoords[i] == 8130004 || hmm->refCoords[i] == 8130005) {
//            st_logInfo("i = %d   hmm->refCoords[i] = %d \n", i, hmm->refCoords[i]);
//        }
        if (hmm->refCoords[i] == 0) {
            st_logInfo("! hmm->refCoords[%d] = 0   leftGap: %d   rightGap: %d   (gaplen: %d  hmmlen: %d)\n", i, gapLeftIndex, gapRightIndex, gapLength, hmm->refLength);
        }
    }
//    st_logInfo("\n");

    // Cleanup
    stRPHmm_destruct(leftHmm, 0);
    stRPHmm_destruct(rightHmm, 0);

    assert(hmm->lastColumn->refCoords[hmm->lastColumn->length - 1] == hmm->refCoords[hmm->refLength - 1]);

    return hmm;
}

void stRPHmm_alignColumns(stRPHmm *hmm1, stRPHmm *hmm2) {
    /*
     * Align the input hmms, modifying them in place, so that they each
     *  (1) span the same reference interval,
     *  (2) have the same number of columns, and
     *  (3) so that for all i, column i in each model span the same interval.
     */
    assert(hmm1 != hmm2);

    // If the two hmms don't overlap in reference space then complain
    if(!stRPHmm_overlapOnReference(hmm1, hmm2)) {
        st_errAbort("Attempting to align two HMMs that do not overlap in reference coordinate space");
    }

    // If hmm1 starts after hmm2 then call the other way around
    if(hmm1->refStart > hmm2->refStart) {
        stRPHmm_alignColumns(hmm2, hmm1);
        return;
    }

    // If hmm1 starts before hmm2 add an empty prefix interval to hmm2
    // so they have the same start coordinate
    if(hmm1->refStart < hmm2->refStart) {
        // Create column
        int64_t gapLeftIndex = findCorrespondingRefCoordIndex(0, hmm1->refCoords, hmm1->referencePriorProbs->refCoordMap);
        int64_t gapRightIndex = findCorrespondingRefCoordIndex(0, hmm2->refCoords, hmm2->referencePriorProbs->refCoordMap);
        int64_t gapLength = gapRightIndex - gapLeftIndex;
//        int64_t numInsertionsInGap= 0;
//        for (int64_t i = 1; hmm1->refCoords[i] <= hmm2->refStart; i++) {
//            if (hmm1->refCoords[i] == hmm1->refCoords[i-1]) numInsertionsInGap++;
//        }
        stRPColumn *column = stRPColumn_construct(hmm1->refStart, gapLength,
                0, NULL, NULL);
        int64_t *colIndexes = st_calloc(column->length, sizeof(int64_t));
        for (int64_t i = 0; i < column->length; i++) {
            // Column will have same reference coordinates as hmm1
            column->refCoords[i] = hmm1->refCoords[i];
            colIndexes[i] = i;
            if (stHash_search(column->refCoordMap, &column->refCoords[i]) == NULL) {
                stHash_insert(column->refCoordMap, &column->refCoords[i], &colIndexes[i]);
            }
        }

        // Add cell
        column->head = stRPCell_construct(0);
        // Create merge column
        stRPMergeColumn *mColumn = stRPMergeColumn_construct(0,0);
        // Add merge cell
        stRPMergeCell_construct(0, 0, mColumn);
        // Create links
        hmm2->firstColumn->pColumn = mColumn;
        mColumn->nColumn = hmm2->firstColumn;
        mColumn->pColumn = column;
        column->nColumn = mColumn;
        assert(column->pColumn == NULL);
        hmm2->firstColumn = column;
        //Adjust start and length of hmm2 interval
        hmm2->refLength += gapLength;
        hmm2->refStart = hmm1->refStart;
        // Increase column number
        hmm2->columnNumber++;

        // Fix ref coords
        int64_t *newRefCoords = st_calloc(hmm2->refLength, sizeof(int64_t));
        int64_t *indexes = st_calloc(hmm2->refLength, sizeof(int64_t));


        // Fill in ref coords at the prefix
        for (int64_t i = 0; i < gapLength; i++) {
            newRefCoords[i] = hmm1->refCoords[i];
            indexes[i] = i;
            if (stHash_search(hmm2->refCoordMap, &newRefCoords[i]) == NULL) {
                stHash_insert(hmm2->refCoordMap, &newRefCoords[i], &indexes[i]);
            }
        }
        // Fill in the rest of the ref coords from hmm2
        for (int64_t i = 0; i < hmm2->refLength-gapLength; i++) {
            int64_t hmm2Index = i + gapLength;
            newRefCoords[hmm2Index] = hmm2->refCoords[i];
            indexes[hmm2Index] = hmm2Index;
            if (stHash_search(hmm2->refCoordMap, &newRefCoords[hmm2Index]) == NULL) {
                stHash_insert(hmm2->refCoordMap, &newRefCoords[hmm2Index], &indexes[hmm2Index]);
            }
        }
        free(hmm2->refCoords);
        hmm2->refCoords = newRefCoords;
        for (int64_t j = 0; j < hmm2->refLength && j < hmm1->refLength; j++) {
            if (hmm2->refCoords[j] == 0) {
                st_logInfo("RefCoords[%d] = 0\n", j);
            }
            if (hmm1->refCoords[j] != hmm2->refCoords[j] ){
                st_logInfo("j: %d   hmm1: %d !=  hmm2: %d  \n", j, hmm1->refCoords[j], hmm2->refCoords[j]);
                int64_t gapLeftIndex = findCorrespondingRefCoordIndex(0, hmm1->refCoords, hmm1->referencePriorProbs->refCoordMap);
                int64_t gapRightIndex = findCorrespondingRefCoordIndex(0, hmm2->refCoords, hmm2->referencePriorProbs->refCoordMap);
                int64_t gapLength = gapRightIndex - gapLeftIndex;
            }
            assert(hmm1->refCoords[j] == hmm2->refCoords[j]);
        }
        assert(hmm2->lastColumn->refCoords[hmm2->lastColumn->length - 1] == hmm2->refCoords[hmm2->refLength - 1]);

    }

    // If hmm1 has a shorter reference interval length than hmm2 then call the function
    // with the hmms reversed.
    if(hmm1->refLength < hmm2->refLength) {
        stRPHmm_alignColumns(hmm2, hmm1);
        return;
    }

    // If hmm1 has a longer reference interval than hmm2 append an empty suffix
    // interval to hmm2 to make them the same length.
    if(hmm1->refLength > hmm2->refLength) {
        // Create column
        int64_t gapRightIndex = findCorrespondingRefCoordIndex(hmm1->refLength-1, hmm1->refCoords, hmm1->referencePriorProbs->refCoordMap);
        int64_t gapLeftIndex = findCorrespondingRefCoordIndex(hmm2->refLength-1, hmm2->refCoords, hmm2->referencePriorProbs->refCoordMap);
        int64_t gapLength = gapRightIndex - gapLeftIndex;
//        int64_t numInsertionsInGap= 0;
//        for (int64_t i = hmm2->refLength; hmm1->refCoords[i] <= hmm1->refCoords[hmm1->refLength - 1]; i++) {
//            if (hmm1->refCoords[i] == hmm1->refCoords[i-1]) numInsertionsInGap++;
//        }
//        int64_t columnLength = hmm1->refCoords[hmm1->refLength - 1] - hmm2->refCoords[hmm2->refLength - 1] + numInsertionsInGap;
//        int64_t columnLength = hmm1->refLength - hmm2->refLength + numInsertionsInGap
        stRPColumn *column = stRPColumn_construct(hmm2->refCoords[hmm2->refLength - 1] + 1,
                gapLength , 0, NULL, NULL);
        int64_t *colIndexes = st_calloc(column->length, sizeof(int64_t));
        for (int64_t j = 0; j < hmm2->refLength; j++) {
            assert(hmm1->refCoords[j] == hmm2->refCoords[j]);
        }
        for (int64_t i = 0; i < column->length; i++) {
            // RefStarts of hmm1 and hmm2 should be the same
            assert(hmm1->refStart == hmm2->refStart);
            column->refCoords[i] = hmm1->refCoords[i + hmm2->refLength];
            colIndexes[i] = i;
            if (stHash_search(column->refCoordMap, &column->refCoords[i]) == NULL) {
                stHash_insert(column->refCoordMap, &column->refCoords[i], &colIndexes[i]);
            }
        }

        // Add cell
        column->head = stRPCell_construct(0);
        // Create merge column
        stRPMergeColumn *mColumn = stRPMergeColumn_construct(0, 0);
        // Add merge cell
        stRPMergeCell_construct(0, 0, mColumn);
        // Create links
        hmm2->lastColumn->nColumn = mColumn;
        mColumn->pColumn = hmm2->lastColumn;
        mColumn->nColumn = column;
        column->pColumn = mColumn;
        assert(column->nColumn == NULL);
        hmm2->lastColumn = column;
        //Adjust start and length of hmm2 interval
        hmm2->refLength = hmm1->refLength;
        // Increase column number
        hmm2->columnNumber++;

        // Fix ref coords
        free(hmm2->refCoords);
        hmm2->refCoords = st_calloc(hmm2->refLength, sizeof(int64_t));
        int64_t *indexes = st_calloc(hmm2->refLength, sizeof(int64_t));
        for (int64_t i = 0; i < hmm2->refLength; i++) {
            hmm2->refCoords[i] = hmm1->refCoords[i];
            indexes[i] = i;
            if (stHash_search(hmm2->refCoordMap, &hmm2->refCoords[i]) == NULL) {
                stHash_insert(hmm2->refCoordMap, &hmm2->refCoords[i], &indexes[i]);
            }
        }
        for (int64_t j = 0; j < hmm2->refLength && j < hmm1->refLength; j++) {
            if (hmm2->refCoords[j] == 0) {
                st_logInfo("RefCoords[%d] = 0\n", j);
            }
            if (hmm1->refCoords[j] != hmm2->refCoords[j] ){
                st_logInfo("j: %d   hmm1: %d !=  hmm2: %d  \n", j, hmm1->refCoords[j], hmm2->refCoords[j]);
                int64_t gapLeftIndex = findCorrespondingRefCoordIndex(0, hmm1->refCoords, hmm1->referencePriorProbs->refCoordMap);
                int64_t gapRightIndex = findCorrespondingRefCoordIndex(0, hmm2->refCoords, hmm2->referencePriorProbs->refCoordMap);
                int64_t gapLength = gapRightIndex - gapLeftIndex;
            }
            assert(hmm1->refCoords[j] == hmm2->refCoords[j]);
        }
        assert(hmm2->lastColumn->refCoords[hmm2->lastColumn->length - 1] == hmm2->refCoords[hmm2->refLength - 1]);
    }

    // Quick coordinate checks
    assert(hmm1->refStart == hmm2->refStart);
    assert(hmm1->refLength == hmm2->refLength);
    assert(hmm1->firstColumn->refStart == hmm1->refStart);
    assert(hmm2->firstColumn->refStart == hmm2->refStart);
    assert(hmm1->lastColumn->refCoords[hmm1->lastColumn->length - 1] == hmm1->refCoords[hmm1->refLength - 1]);
    assert(hmm2->lastColumn->refCoords[hmm2->lastColumn->length - 1] == hmm2->refCoords[hmm2->refLength - 1]);
//    assert(hmm1->lastColumn->refStart + hmm1->lastColumn->length == hmm1->refStart + hmm1->refLength);
//    assert(hmm2->lastColumn->refStart + hmm2->lastColumn->length == hmm2->refStart + hmm2->refLength);

    // At this point both hmms have the same reference interval

    // While one hmm has a shorter reference interval than the other split the other interval
    // otherwise move on to the next
    stRPColumn *column1 = hmm1->firstColumn;
    stRPColumn *column2 = hmm2->firstColumn;
    while(1) {
        assert(column1->refStart == column2->refStart);

        if(column1->length > column2->length) {
            stRPColumn_split(column1, column2->length, hmm1);
//            assert(column1->nColumn->nColumn->refStart == column1->refStart + column2->length);
        }
        else if(column1->length < column2->length) {
            stRPColumn_split(column2, column1->length, hmm2);
        }

        assert(column1->refStart == column2->refStart);
        assert(column1->length == column2->length); // Now have equal length/start

        // There are no more columns, so break
        if(column1->nColumn == NULL) {
            assert(hmm1->lastColumn == column1);
            assert(column2->nColumn == NULL);
            assert(hmm2->lastColumn == column2);
            break;
        }

        column1 = column1->nColumn->nColumn;
        assert(column2->nColumn != NULL);
        column2 = column2->nColumn->nColumn;
        assert(column1 != NULL);
        assert(column2 != NULL);
    }

    assert(hmm1->columnNumber == hmm2->columnNumber);
}

stRPHmm *stRPHmm_createCrossProductOfTwoAlignedHmm(stRPHmm *hmm1, stRPHmm *hmm2) {
    /*
     *  For two aligned hmms (see stRPHmm_alignColumns) returns a new hmm that represents the
     *  cross product of all the states of the two input hmms.
     */

    // Do sanity checks that the two hmms have been aligned
    if(!stString_eq(hmm1->referenceName, hmm2->referenceName)) {
        st_errAbort("Trying to create cross product of two HMMs "
                "on different reference sequences");
    }
    if(hmm1->refStart != hmm2->refStart) {
        st_errAbort("Trying to create cross product of two HMMs "
                "with different reference interval starts");
    }
    if(hmm1->refLength != hmm2->refLength) {
        st_errAbort("Trying to create cross product of two HMMs "
                "with different reference interval length");
    }
    if(hmm1->columnNumber != hmm2->columnNumber) {
        st_errAbort("Trying to create cross product of two HMMs "
                "with different column numbers");
    }

    // Create a new empty hmm
    stRPHmm *hmm = st_calloc(1, sizeof(stRPHmm));
    // Set the reference interval
    hmm->referenceName = stString_copy(hmm1->referenceName);
    hmm->refStart = hmm1->refStart;
    hmm->refLength = hmm1->refLength;
    // Create the combined list of profile seqs
    hmm->profileSeqs = stList_copy(hmm1->profileSeqs, NULL);
    stList_appendAll(hmm->profileSeqs, hmm2->profileSeqs);
    // Set column number
    hmm->columnNumber = hmm1->columnNumber;
    // Set substitution matrices
    if(hmm1->parameters != hmm2->parameters) {
        st_errAbort("Hmm parameters differ in fuse function, panic.");
    }
    hmm->parameters = hmm1->parameters;
    // Set reference position prior probabilities
    if(hmm1->referencePriorProbs != hmm2->referencePriorProbs) {
        st_errAbort("Hmm reference prior probs differ in hmm cross product function, panic.");
    }
    hmm->referencePriorProbs = hmm1->referencePriorProbs;

    // Set reference coordinates
    hmm->refCoords = st_calloc(hmm->refLength, sizeof(int64_t));
    hmm->refCoordMap = stHash_construct3(stHash_stringKey, stHash_intPtrEqualKey, NULL, NULL);
    int64_t *indexes = st_calloc(hmm->refLength, sizeof(int64_t));
    for (int64_t i = 0; i < hmm->refLength; i++) {
        hmm->refCoords[i] = hmm1->refCoords[i];
        indexes[i] = i;
        if (stHash_search(hmm->refCoordMap, &hmm->refCoords[i]) == NULL) {
            stHash_insert(hmm->refCoordMap, &hmm->refCoords[i], &indexes[i]);
        }
    }

    // For each pair of corresponding columns
    stRPColumn *column1 = hmm1->firstColumn;
    stRPColumn *column2 = hmm2->firstColumn;
    assert(column1 != NULL);
    assert(column2 != NULL);
    stRPMergeColumn *mColumn = NULL;

    while(1) {
        // Check columns aligned
        assert(column1->refStart == column2->refStart);
        assert(column1->length == column2->length);

        // Create the new column

        // Depth
        int64_t newColumnDepth = column1->depth+column2->depth;
        if(newColumnDepth > hmm->maxDepth) {
            hmm->maxDepth = newColumnDepth;
        }

        // Seq headers
        stProfileSeq **seqHeaders = st_malloc(sizeof(stProfileSeq *) * newColumnDepth);
        memcpy(seqHeaders, column1->seqHeaders, sizeof(stProfileSeq *) * column1->depth);
        memcpy(&seqHeaders[column1->depth], column2->seqHeaders, sizeof(stProfileSeq *) * column2->depth);

        // Profiles
        uint8_t **seqs = st_malloc(sizeof(uint8_t *) * newColumnDepth);
        memcpy(seqs, column1->seqs, sizeof(uint8_t *) * column1->depth);
        memcpy(&seqs[column1->depth], column2->seqs, sizeof(uint8_t *) * column2->depth);

        stRPColumn *column = stRPColumn_construct(column1->refStart, column1->length,
                newColumnDepth, seqHeaders, seqs);
        int64_t *indexes = st_calloc(column->length, sizeof(int64_t));
        for (int64_t i = 0; i < column->length; i++) {
            column->refCoords[i] = column1->refCoords[i];
            indexes[i] = i;
            if (stHash_search(column->refCoordMap, &column->refCoords[i]) == NULL) {
                stHash_insert(column->refCoordMap, &column->refCoords[i], &indexes[i]);
            }
        }

        // If the there is a previous column
        if(mColumn != NULL) {
            mColumn->nColumn = column;
            column->pColumn = mColumn;
        }
        else {
            hmm->firstColumn = column;
            assert(column->pColumn == NULL);
        }

        // Create cross product of columns
        stRPCell **pCell = &column->head;
        stRPCell *cell1 = column1->head;
        do {
            stRPCell *cell2 = column2->head;
            do {
                stRPCell *cell = stRPCell_construct(mergePartitionsOrMasks(cell1->partition, cell2->partition,
                        column1->depth, column2->depth));
                // Link cells
                *pCell = cell;
                pCell = &cell->nCell;
            } while((cell2 = cell2->nCell) != NULL);
        } while((cell1 = cell1->nCell) != NULL);

        // Get the next merged column
        stRPMergeColumn *mColumn1 = column1->nColumn;
        stRPMergeColumn *mColumn2 = column2->nColumn;

        // If column is NULL, we have reached the last column
        // and we can exit
        if(mColumn1 == NULL) {
            assert(mColumn2 == NULL);
            assert(hmm1->lastColumn == column1);
            assert(hmm2->lastColumn == column2);

            // Set the last column pointer
            hmm->lastColumn = column;
            break;
        }

        // Create new merged column
        uint64_t fromMask = mergePartitionsOrMasks(mColumn1->maskFrom, mColumn2->maskFrom,
                mColumn1->pColumn->depth, mColumn2->pColumn->depth);
        uint64_t toMask = mergePartitionsOrMasks(mColumn1->maskTo, mColumn2->maskTo,
                        mColumn1->nColumn->depth, mColumn2->nColumn->depth);
        mColumn = stRPMergeColumn_construct(fromMask, toMask);

        // Connect links
        mColumn->pColumn = column;
        column->nColumn = mColumn;

        // Create cross product of merged columns
        stHashIterator *cellIt1 = stHash_getIterator(mColumn1->mergeCellsFrom);
        stRPMergeCell *mCell1;
        while((mCell1 = stHash_getNext(cellIt1)) != NULL) {
            stHashIterator *cellIt2 = stHash_getIterator(mColumn2->mergeCellsFrom);
            stRPMergeCell *mCell2;
            while((mCell2 = stHash_getNext(cellIt2)) != NULL) {
                uint64_t fromPartition = mergePartitionsOrMasks(mCell1->fromPartition,
                        mCell2->fromPartition,
                        mColumn1->pColumn->depth, mColumn2->pColumn->depth);

                uint64_t toPartition = mergePartitionsOrMasks(mCell1->toPartition,
                        mCell2->toPartition,
                        mColumn1->nColumn->depth, mColumn2->nColumn->depth);

                stRPMergeCell_construct(fromPartition, toPartition, mColumn);
            }
            stHash_destructIterator(cellIt2);
        }
        stHash_destructIterator(cellIt1);

        // Get next column
        column1 = mColumn1->nColumn;
        column2 = mColumn2->nColumn;
        assert(column1 != NULL);
        assert(column2 != NULL);
    }

    return hmm;
}

static void stRPHmm_initialiseProbs(stRPHmm *hmm) {
    /*
     * Initialize the forward and backward matrices.
     */
    // Initialize total forward and backward probabilities
    hmm->forwardLogProb = ST_MATH_LOG_ZERO;
    hmm->backwardLogProb = ST_MATH_LOG_ZERO;

    // Iterate through columns from first to last
    stRPColumn *column = hmm->firstColumn;
    while(1) {
        // Set total log prob
        column->totalLogProb = ST_MATH_LOG_ZERO;

        // Initialise cells in the column
        stRPCell *cell = column->head;
        do {
            cell->forwardLogProb = ST_MATH_LOG_ZERO;
            cell->backwardLogProb = ST_MATH_LOG_ZERO;
        } while((cell = cell->nCell) != NULL);

        if(column->nColumn == NULL) {
            break;
        }

        // Initialise cells in the next merge column
        stList *mergeCells = stHash_getValues(column->nColumn->mergeCellsFrom);
        for(int64_t i=0; i<stList_length(mergeCells); i++) {
            stRPMergeCell *mergeCell = stList_get(mergeCells, i);
            mergeCell->forwardLogProb = ST_MATH_LOG_ZERO;
            mergeCell->backwardLogProb = ST_MATH_LOG_ZERO;
        }
        stList_destruct(mergeCells);

        column = column->nColumn->nColumn;
    }
}

static inline void forwardCellCalc1(stRPHmm *hmm, stRPColumn *column, stRPCell *cell, uint64_t *bitCountVectors) {
    // If the previous merge column exists then propagate forward probability from merge state
    if(column->pColumn != NULL) {
        stRPMergeCell *mCell = stRPMergeColumn_getPreviousMergeCell(cell, column->pColumn);
        cell->forwardLogProb = mCell->forwardLogProb;
    }
    // Otherwise initialize probability with log(1.0)
    else {
        cell->forwardLogProb = ST_MATH_LOG_ONE;
    }

    // Calculate the emission prob
    double emissionProb = emissionLogProbability(column, cell, bitCountVectors,
            hmm->referencePriorProbs, (stRPHmmParameters *)hmm->parameters);

    // Add emission prob to forward log prob
    cell->forwardLogProb += emissionProb;

    // Store the emission probability for the cell in the backwardLogProb field temporarily
    // (is corrected during the backward pass)
    cell->backwardLogProb = emissionProb;
}

static inline void forwardCellCalc2(stRPHmm *hmm, stRPColumn *column, stRPCell *cell) {
    // If the next merge column exists then propagate forward probability to the merge state
    if (column->nColumn != NULL) {
        // Add to the next merge cell
        stRPMergeCell *mCell = stRPMergeColumn_getNextMergeCell(cell, column->nColumn);
        mCell->forwardLogProb = logAddP(mCell->forwardLogProb, cell->forwardLogProb,
                hmm->parameters->maxNotSumTransitions);
    } else {
        // Else propagate probability to total forward probability of model
        hmm->forwardLogProb = logAddP(hmm->forwardLogProb, cell->forwardLogProb, hmm->parameters->maxNotSumTransitions);
    }
}

static void stRPHmm_forward(stRPHmm *hmm) {
    /*
     * Forward algorithm for hmm.
     */
    stRPColumn *column = hmm->firstColumn;

    // Iterate through columns from first to last
    while(1) {
        // Get the bit count vectors for the column
        uint64_t *bitCountVectors = calculateCountBitVectors(column->seqs, column->depth,
                column->length);

        // Iterate through states in column
        stRPCell *cell = column->head;

        // If OpenMP is available then parallelize the calculation of the emission calcs
#if defined(_OPENMP)
        stRPCell *cells[CELL_BUFFER_SIZE];
        do {
            // Get as many cells as the buffer will fit / there are cells
            int64_t cellsInBuffer=0;
            do {
                cells[cellsInBuffer++] = cell;
            } while((cell = cell->nCell) != NULL && cellsInBuffer < CELL_BUFFER_SIZE);

#pragma omp parallel
{
#pragma omp for
            for(int64_t i=0; i<cellsInBuffer; i++) {
                forwardCellCalc1(hmm, column, cells[i], bitCountVectors);
            }
}
            for(int64_t i=0; i<cellsInBuffer; i++) {
                forwardCellCalc2(hmm, column, cells[i]);
            }
        } while(cell != NULL);
#else
        // Otherwise do it without the need for the cell buffer
        do {
            forwardCellCalc1(hmm, column, cell, bitCountVectors);
            forwardCellCalc2(hmm, column, cell);
        }
        while((cell = cell->nCell) != NULL);
#endif

        // Cleanup the bit count vectors
        free(bitCountVectors);

        if(column->nColumn == NULL) {
            break;
        }
        column = column->nColumn->nColumn;
    }
}

static inline void backwardCellCalc(stRPHmm *hmm, stRPColumn *column, stRPCell *cell) {
    // Retrieve the emission probability that was stored by the forward pass
    double probabilityToPropagateLogProb = cell->backwardLogProb;

    // If the next merge column exists then propagate backward probability from merge state
    if(column->nColumn != NULL) {
        stRPMergeCell *mCell = stRPMergeColumn_getNextMergeCell(cell, column->nColumn);
        cell->backwardLogProb = mCell->backwardLogProb;
        probabilityToPropagateLogProb += mCell->backwardLogProb;
    }
    else { // Else set the backward prob to log(1)
        cell->backwardLogProb = ST_MATH_LOG_ONE;
    }

    // If the previous merge column exists then propagate backward probability to the merge state
    if(column->pColumn != NULL) {
        // Add to the previous merge cell
        stRPMergeCell *mCell = stRPMergeColumn_getPreviousMergeCell(cell, column->pColumn);
        mCell->backwardLogProb = logAddP(mCell->backwardLogProb, probabilityToPropagateLogProb,
                hmm->parameters->maxNotSumTransitions);
    }
    else {
        hmm->backwardLogProb = logAddP(hmm->backwardLogProb, probabilityToPropagateLogProb,
                hmm->parameters->maxNotSumTransitions);
    }

    // Add to column total probability
    column->totalLogProb = logAddP(column->totalLogProb,
                 cell->forwardLogProb + cell->backwardLogProb, hmm->parameters->maxNotSumTransitions);
}

static void stRPHmm_backward(stRPHmm *hmm) {
    /*
     * Backward algorithm for hmm.
     */
    stRPColumn *column = hmm->lastColumn;

    // Iterate through columns from last to first
    while(1) {
        // Iterate through states in column
        stRPCell *cell = column->head;

        // Otherwise do it without the need for the cell buffer
        do {
            backwardCellCalc(hmm, column, cell);
        }
        while((cell = cell->nCell) != NULL);

        if(column->pColumn == NULL) {
            break;
        }
        column = column->pColumn->pColumn;
    }
}

void stRPHmm_forwardBackward(stRPHmm *hmm) {
    /*
     * Runs the forward and backward algorithms and sets the total column probabilities.
     *
     * This function must be run upon an HMM to calculate cell posterior probabilities.
     */
    // Initialise state values
    stRPHmm_initialiseProbs(hmm);
    // Run the forward and backward passes
    stRPHmm_forward(hmm);
    stRPHmm_backward(hmm);
}

static int cellCmpFn(const void *a, const void *b, const void *extraArg) {
    /*
     * Sort cells by posterior probability in descending order.
     */
    stRPCell *cell1 = (stRPCell *)a, *cell2 = (stRPCell *)b;
    stRPColumn *column = (stRPColumn *)extraArg;
    double p1 = stRPCell_posteriorProb(cell1, column), p2 = stRPCell_posteriorProb(cell2, column);
    return p1 > p2 ? -1 : p1 < p2 ? 1 : 0;
}

static int mergeCellCmpFn(const void *a, const void *b, const void *extraArg) {
    /*
     * Sort merge cells by posterior probability in descending order.
     */
    stRPMergeCell *cell1 = (stRPMergeCell *)a, *cell2 = (stRPMergeCell *)b;
    stRPMergeColumn *column = (stRPMergeColumn *)extraArg;
    double p1 = stRPMergeCell_posteriorProb(cell1, column), p2 = stRPMergeCell_posteriorProb(cell2, column);
    return p1 > p2 ? -1 : p1 < p2 ? 1 : 0;
}

void filterMergeCells(stRPMergeColumn *mColumn, stSet *chosenMergeCellsSet) {
    /*
     * Removes merge cells from the column that are not in chosenMergeCellsSet
     */
    assert(stSet_size(chosenMergeCellsSet) > 0);
    stList *mergeCells = stHash_getValues(mColumn->mergeCellsFrom);
    for(int64_t i=0; i<stList_length(mergeCells); i++) {
        stRPMergeCell *mCell = stList_get(mergeCells, i);
        assert(mCell != NULL);
        if(stSet_search(chosenMergeCellsSet, mCell) == NULL) {
            // Remove the state from the merge column
            assert(stHash_search(mColumn->mergeCellsFrom, &(mCell->fromPartition)) == mCell);
            assert(stHash_search(mColumn->mergeCellsTo, &(mCell->toPartition)) == mCell);
            stHash_remove(mColumn->mergeCellsFrom, &(mCell->fromPartition));
            stHash_remove(mColumn->mergeCellsTo, &(mCell->toPartition));

            // Cleanup
            stRPMergeCell_destruct(mCell);
        }
    }
    stList_destruct(mergeCells);
    assert(stSet_size(chosenMergeCellsSet) == stHash_size(mColumn->mergeCellsFrom));
    assert(stSet_size(chosenMergeCellsSet) == stHash_size(mColumn->mergeCellsTo));
}

stSet *getLinkedMergeCells(stRPMergeColumn *mColumn,
        stRPMergeCell *(*getNCell)(stRPCell *, stRPMergeColumn *),
        stList *cells) {
    /*
     * Returns the set of merge cells in the column that are linked to a cell
     * in cells.
     */
    stSet *chosenMergeCellsSet = stSet_construct();
    for(int64_t i=0; i<stList_length(cells); i++) {
        stRPMergeCell *mCell = getNCell(stList_get(cells, i), mColumn);
        assert(mCell != NULL);
        stSet_insert(chosenMergeCellsSet, mCell);
    }
    assert(stSet_size(chosenMergeCellsSet) > 0);
    return chosenMergeCellsSet;
}

void relinkCells(stRPColumn *column, stList *cells) {
    /*
     * Re-links the cells in the list 'cells' to make up the list of cells in the column.
     */
    stRPCell **pCell = &column->head; // Pointer to previous cell, used to
    // remove cells from the linked list
    for(int64_t i=0; i<stList_length(cells); i++) {
        stRPCell *cell = stList_get(cells, i);
        *pCell = cell;
        pCell = &cell->nCell;
    }
    *pCell = NULL;
    assert(column->head != NULL);
}

stList *getLinkedCells(stRPColumn *column,
        stRPMergeCell *(*getPCell)(stRPCell *, stRPMergeColumn *),
        stRPMergeColumn *mColumn) {
    /*
     * Returns the set of cells in column that are linked to a cell in mColumn.
     */

    // Put cells into an array and sort by descending posterior prob
    // only keeping cells that still have a preceding merge cell

    stList *cells = stList_construct();
    stRPCell *cell = column->head;
    do {
        if(mColumn == NULL || getPCell(cell, mColumn) != NULL) {
            stList_append(cells, cell);
            cell = cell->nCell;
        }
        else {
            stRPCell *nCell = cell->nCell;
            stRPCell_destruct(cell);
            cell = nCell;
        }
    } while(cell != NULL);
    stList_sort2(cells, cellCmpFn, column);
    assert(stList_length(cells) > 0);

    return cells;
}

void stRPHmm_pruneForwards(stRPHmm *hmm) {
    /*
     * Remove cells from hmm whos posterior probability is below the given threshold
     */

    // For each column
    stRPColumn *column = hmm->firstColumn;
    stRPMergeColumn *mColumn = NULL;

    while(1) {
        assert(column->head != NULL);

        // Get cells that have a valid previous cell
        stList *cells = getLinkedCells(column, stRPMergeColumn_getPreviousMergeCell, mColumn);

        // Get rid of the excess cells
        while(stList_length(cells) > hmm->parameters->maxPartitionsInAColumn) {
            stRPCell_destruct(stList_pop(cells));
        }

        // Relink the cells (from most probable to least probable)
        relinkCells(column, cells);

        // Move on to the next merge column
        mColumn = column->nColumn;

        if(mColumn == NULL) {
            assert(column == hmm->lastColumn);
            stList_destruct(cells);
            break;
        }

        //  Get merge cells that are connected to a cell in the previous column
        stSet *chosenMergeCellsSet = getLinkedMergeCells(mColumn,
                stRPMergeColumn_getNextMergeCell, cells);

        // Shrink the the number of chosen cells to less than equal to the desired number
        stList *chosenMergeCellsList = stSet_getList(chosenMergeCellsSet);
        stList_sort2(chosenMergeCellsList, mergeCellCmpFn, mColumn);
        while(stList_length(chosenMergeCellsList) > hmm->parameters->maxPartitionsInAColumn) {
            stSet_remove(chosenMergeCellsSet, stList_pop(chosenMergeCellsList));
        }
        assert(stList_length(chosenMergeCellsList) == stSet_size(chosenMergeCellsSet));
        stList_destruct(chosenMergeCellsList);

        // Get rid of merge cells we don't need
        filterMergeCells(mColumn, chosenMergeCellsSet);

        // Cleanup
        stList_destruct(cells);
        stSet_destruct(chosenMergeCellsSet);

        column = mColumn->nColumn;
    }
}

void stRPHmm_pruneBackwards(stRPHmm *hmm) {
    /*
     * Remove cells from hmm whos posterior probability is below the given threshold
     */

    // For each column
    stRPColumn *column = hmm->lastColumn;
    stRPMergeColumn *mColumn = NULL;

    while(1) {
        assert(column->head != NULL);

        // Get cells that have a valid previous cell
        stList *cells = getLinkedCells(column, stRPMergeColumn_getNextMergeCell, mColumn);

        // This must be true because the forward pass has already winnowed the number below the
        // threshold
        assert(stList_length(cells) <= hmm->parameters->maxPartitionsInAColumn);

        // Relink the cells (from most probable to least probable)
        relinkCells(column, cells);

        // Move on to the next merge column
        mColumn = column->pColumn;

        if(mColumn == NULL) {
            assert(column == hmm->firstColumn);
            stList_destruct(cells);
            break;
        }

        //  Get merge cells that are connected to a cell in the previous column
        stSet *chosenMergeCellsSet = getLinkedMergeCells(mColumn,
                stRPMergeColumn_getPreviousMergeCell, cells);

        // By the same logic, this number if pruned on the forwards pass
        assert(stSet_size(chosenMergeCellsSet) <= hmm->parameters->maxPartitionsInAColumn);

        // Get rid of merge cells we don't need
        filterMergeCells(mColumn, chosenMergeCellsSet);

        // Cleanup
        stList_destruct(cells);
        stSet_destruct(chosenMergeCellsSet);

        column = mColumn->pColumn;
    }
}

void stRPHmm_prune(stRPHmm *hmm) {
    stRPHmm_pruneForwards(hmm);
    stRPHmm_pruneBackwards(hmm);
}

bool stRPHmm_overlapOnReference(stRPHmm *hmm1, stRPHmm *hmm2) {
    /*
     * Return non-zero iff hmm1 and hmm2 have the same reference sequence and overlapping
     * coordinates intervals on that reference sequence.
     */

    // If either interval is zero length this is not a well defined comparison
    if(hmm1->refLength <= 0 || hmm2->refLength <= 0) {
        st_errAbort("Trying to compare HMMs with a zero length coordinate interval");
    }

    // Check if on the same reference sequence
    if(!stString_eq(hmm1->referenceName, hmm2->referenceName)) {
        return 0;
    }

    // Check if intervals overlap

    // If hmm1 starts after hmm2's start coordinate then switch hmm1 for hmm2
    if(hmm1->refStart > hmm2->refStart) {
        return stRPHmm_overlapOnReference(hmm2, hmm1);
    }

    // The coordinates of the first interval overlap the second
    return hmm1->refCoords[hmm1->refLength-1] >= hmm2->refStart;
//    return hmm1->refStart + hmm1->refLength > hmm2->refStart;
}

static stRPColumn *getColumn(stRPColumn *column, int64_t site) {
    /*
     * Returns column containing the given reference position, starting from the linked, preceding column "column".
     */
    assert(column != NULL);
    while(1) {
        assert(site >= column->refStart);
//        if(site < column->refStart + column->length) {
        if (site <= column->refCoords[column->length-1]) {
            return column;
        }
        if(column->nColumn == NULL) {
            break;
        }
        column = column->nColumn->nColumn;
    }
    st_errAbort("Site: %" PRIi64 " not contained in hmm\n", site);
    return column;
}

void stRPHmm_resetColumnNumberAndDepth(stRPHmm *hmm) {
    /*
     * Walk through the hmm calculate and set the maxDepth and column number.
     */
    hmm->columnNumber = 0;
    hmm->maxDepth = 0;
    stRPColumn *column = hmm->firstColumn;
    while(1) {
        hmm->columnNumber++;
        if(hmm->maxDepth < column->depth) {
            hmm->maxDepth = column->depth;
        }
        if(column->nColumn == NULL) {
            break;
        }
        column = column->nColumn->nColumn;
    }
}

stRPHmm *stRPHmm_split(stRPHmm *hmm, int64_t splitPoint) {
    /*
     * Splits the hmm into two at the specified point, given by the reference coordinate splitPiunt. The return value
     * is the suffix of the split, whose reference start is splitPoint.
     * The prefix of the split is the input hmm, which has its suffix cleaved off. Its length is then splitPoint-hmm->refStart.
     */

//    if (splitPoint == 8092070) {
        st_logInfo("Hmm for splitPoint %d:  %d - %d (%d)\n", splitPoint, hmm->refStart, hmm->refCoords[hmm->refLength-1], hmm->refLength);
//    }
    if(splitPoint <= hmm->refStart) {
        st_errAbort("The split point is at or before the start of the reference interval\n");
    }
    if(splitPoint > hmm->refCoords[hmm->refLength -1]) {
        st_errAbort("The split point %" PRIi64 " is after the last position of the reference interval\n", splitPoint);
    }
    assert(splitPoint <= hmm->refCoords[hmm->refLength-1]); // necessary?



    stRPHmm *suffixHmm = st_calloc(1, sizeof(stRPHmm));

    // Set the reference interval for the two hmms
    suffixHmm->referenceName = stString_copy(hmm->referenceName);
    suffixHmm->refStart = splitPoint;
    int64_t *splitPointIndex = stHash_search(hmm->refCoordMap, &splitPoint);
    assert(splitPointIndex != NULL);
    suffixHmm->refLength = hmm->refLength + hmm->refStart - splitPoint;
//    hmm->refLength = splitPoint - hmm->refStart;
    hmm->refLength = *splitPointIndex;
    assert(hmm->refLength > 0);
    assert(suffixHmm->refLength > 0);
    if (splitPoint == 8092070) {
        st_logInfo("Suffix length: %d   hmm length: %d   total: %d \n", suffixHmm->refLength, hmm->refLength, suffixHmm->refLength + hmm->refLength);
    }

    // Parameters
    suffixHmm->parameters = hmm->parameters;

    // Reference prior probabilities
    suffixHmm->referencePriorProbs = hmm->referencePriorProbs;

    // Divide the profile sequences between the two hmms (some may end in both if they span the interval)
    suffixHmm->profileSeqs = stList_construct();
    stList *prefixProfileSeqs = stList_construct();
    for(int64_t i=0; i<stList_length(hmm->profileSeqs); i++) {
        stProfileSeq *pSeq = stList_get(hmm->profileSeqs, i);
        if(pSeq->refStart < splitPoint) {
            stList_append(prefixProfileSeqs, pSeq);
        }
        if(pSeq->refEnd > splitPoint) {
            stList_append(suffixHmm->profileSeqs, pSeq);
        }
    }
    stList_destruct(hmm->profileSeqs);
    hmm->profileSeqs = prefixProfileSeqs;

    // Set ref coords
    //TODO can this share with hmm? how to deal with deallocation?
    suffixHmm->refCoords = st_calloc(suffixHmm->refLength, sizeof(int64_t));
    suffixHmm->refCoordMap = stHash_construct3(stHash_stringKey, stHash_intPtrEqualKey, NULL, NULL);
    int64_t *indexes = st_calloc(suffixHmm->refLength, sizeof(int64_t));
    for (int64_t i = 0; i < suffixHmm->refLength; i++) {
        suffixHmm->refCoords[i] = hmm->refCoords[i];
        indexes[i] = i;
        if (stHash_search(suffixHmm->refCoordMap, &suffixHmm->refCoords[i]) == NULL) {
            stHash_insert(suffixHmm->refCoordMap, &suffixHmm->refCoords[i], &indexes[i]);
        }
    }

    // Get the column containing the split point
    stRPColumn *splitColumn = getColumn(hmm->firstColumn, splitPoint);
    assert(splitColumn != NULL);
    assert(splitColumn->refStart <= splitPoint);
    assert(splitPoint <= splitColumn->refCoords[splitColumn->length-1]);

    // If the split point is within the column, split the column
    if(splitPoint > splitColumn->refStart) {
        // TODO check this
        int64_t *colSplitPointIndex = stHash_search(splitColumn->refCoordMap, &splitPoint);
//        int64_t colSplitPointIndex = splitPoint-splitColumn->refStart;
        assert(colSplitPointIndex != NULL);
//        stRPColumn_split(splitColumn, *splitPointIndex, hmm);
        stRPColumn_split(splitColumn, *colSplitPointIndex, hmm);
        splitColumn = splitColumn->nColumn->nColumn;
        assert(splitPoint == splitColumn->refStart);
    }

    // Set links between columns
    suffixHmm->firstColumn = splitColumn;
    suffixHmm->lastColumn = hmm->lastColumn;
    hmm->lastColumn = splitColumn->pColumn->pColumn;
    hmm->lastColumn->nColumn = NULL;
    stRPMergeColumn_destruct(splitColumn->pColumn); // Cleanup the merge column that is deleted by this pointer setting
    splitColumn->pColumn = NULL;

    // Set depth and column numbers
    stRPHmm_resetColumnNumberAndDepth(hmm);
    stRPHmm_resetColumnNumberAndDepth(suffixHmm);

    return suffixHmm;
}

static bool sitesLinkageIsWellSupported(stRPHmm *hmm, int64_t leftSite, int64_t rightSite) {
    /*
     * Returns true if the two sites, specified by reference coordinates leftSite and rightSite, are linked by
     * hmm->parameters->minReadCoverageToSupportPhasingBetweenHeterozygousSites, otherwise false.
     */
    stRPColumn *leftColumn = getColumn(hmm->firstColumn, leftSite);
    stRPColumn *rightColumn = getColumn(leftColumn, rightSite);

    stSet *sequencesInCommon = stRPColumn_getSequencesInCommon(leftColumn, rightColumn);

    // Condition to determine if well supported by reads
    bool wellSupported = stSet_size(sequencesInCommon) >= hmm->parameters->minReadCoverageToSupportPhasingBetweenHeterozygousSites;

    // Cleanup
    stSet_destruct(sequencesInCommon);

    return wellSupported;
}

stList *stRPHMM_splitWherePhasingIsUncertain(stRPHmm *hmm) {
    /*
     * Takes the input hmm and splits into a sequence of contiguous fragments covering the same reference interval,
     * returned as an ordered list of hmm fragments. Hmms are split where there is insufficient support between heterozygous
     * sites to support phasing between the two haplotypes. See sitesLinkageIsWellSupported for details.
     */

    // Run the forward-backward algorithm
    stRPHmm_forwardBackward(hmm);

    // Now compute a high probability path through the hmm
    stList *path = stRPHmm_forwardTraceBack(hmm);

    // Get two haplotypes for the path through the HMM
    stGenomeFragment *gF = stGenomeFragment_construct(hmm, path);

    // Find high confidence heterozygous sites
    stList *hetSites = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    for(int64_t i=0; i<gF->length; i++) {
        // If heterozygous site
        if (gF->refCoords[i] == 8207224) {
            st_logInfo("POS: %d   h1: %d   h2: %d\n", gF->refCoords[i], gF->haplotypeString1[i], gF->haplotypeString2[i]);
        }
        if(gF->haplotypeString1[i] != gF->haplotypeString2[i]) {
//            stList_append(hetSites, stIntTuple_construct1(gF->refStart + i));
            stList_append(hetSites, stIntTuple_construct1(gF->refCoords[i]));
            st_logInfo("HET SITE: gF->RC = %d \n", gF->refCoords[i]);
        }
    }

    // Split hmms
    stList *splitHmms = stList_construct3(0, (void (*)(void *))stRPHmm_destruct2);

    // For each pair of contiguous het sites if not supported by sufficient reads split the hmm
    for(int64_t i=0; i<stList_length(hetSites)-1; i++) {
        int64_t j = stIntTuple_get(stList_get(hetSites, i), 0);
        int64_t k = stIntTuple_get(stList_get(hetSites, i+1), 0);
        assert(k > j);

        // If not well supported by reads
        if(!sitesLinkageIsWellSupported(hmm, j, k)) {
            // Split hmm
            int64_t splitPoint = j+(k-j+1)/2;
            st_logInfo("splitPoint: %d   hmm->RC: %d\n", splitPoint, hmm->refCoords[hmm->refLength -1]);
            st_logInfo("j = %d,  k = %d\n", j, k);
            stRPHmm *rightHmm = stRPHmm_split(hmm, splitPoint);
            assert(rightHmm->refStart == splitPoint);
            if (hmm->refCoords[hmm->refLength-1] +1 != splitPoint) {
                st_logInfo("splitPoint: %d   hmm->RC: %d\n", splitPoint, hmm->refCoords[hmm->refLength -1]);
                st_logInfo("j = %d,  k = %d\n", j, k);
            }
            assert(hmm->refCoords[hmm->refLength - 1]+1 == splitPoint);
//            assert(hmm->refStart + hmm->refLength == splitPoint);
            // Add prefix of hmm to list of split hmms
            stList_append(splitHmms, hmm);
            // Set hmm as right hmm
            hmm = rightHmm;
        }
    }

    // Add the remaining part of the hmm to split hmms
    stList_append(splitHmms, hmm);

    // Cleanup
    stList_destruct(hetSites);
    stList_destruct(path);
    stGenomeFragment_destruct(gF);

    return splitHmms;
}
