/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */


#include <htslib/vcf.h>
#include "CuTest.h"
#include "sonLib.h"
#include "stRPHmm.h"
#include "../externalTools/sonLib/C/impl/sonLibListPrivate.h"

void printSequenceStats(FILE *fH, stList *profileSequences) {
    /*
     * Print stats about the set of profile sequences.
     */
    int64_t totalLength=0;
    for(int64_t i=0; i<stList_length(profileSequences); i++) {
        stProfileSeq *profileSeq = stList_get(profileSequences, i);
        totalLength += profileSeq->length;
    }
    fprintf(fH, "Got %" PRIi64 " profile sequences, with total length: %" PRIi64 ", average length: %f\n",
            stList_length(profileSequences), totalLength, ((float)totalLength)/stList_length(profileSequences));
}

double getExpectedNumberOfMatches(uint64_t *haplotypeString, int64_t start, int64_t length, stProfileSeq *profileSeq) {
    /*
     * Returns the expected number of positions in the profile sequence that are identical to the given haplotype string.
     */
    double totalExpectedMatches = 0.0;

    for(int64_t i=0; i<profileSeq->length; i++) {
        // Get base in the haplotype sequence
        int64_t j = i + profileSeq->refStart - start;
        if(j >= 0 && j < length) {
            uint64_t hapBase = haplotypeString[j];
            assert(hapBase < ALPHABET_SIZE);

            // Expectation of a match
            totalExpectedMatches += getProb(&(profileSeq->profileProbs[i * ALPHABET_SIZE]), hapBase);
        }
    }
    return totalExpectedMatches;
}

double getExpectedIdentity(uint64_t *haplotypeString, int64_t start, int64_t length, stSet *profileSeqs) {
    /*
     * Returns the expected fraction of positions in the profile sequences that match their corresponding position in the
     * given haplotype string.
     */
    double totalExpectedNumberOfMatches = 0.0;
    int64_t totalLength = 0;
    stSetIterator *it = stSet_getIterator(profileSeqs);
    stProfileSeq *pSeq;
    while((pSeq = stSet_getNext(it)) != NULL) {
        totalExpectedNumberOfMatches += getExpectedNumberOfMatches(haplotypeString, start, length, pSeq);
        totalLength += pSeq->length;
    }
    return totalExpectedNumberOfMatches/totalLength;
}

double getIdentityBetweenHaplotypes(uint64_t *hap1String, uint64_t *hap2String, int64_t length) {
    /*
     * Returns the fraction of positions in two haplotypes that are identical.
     */
    int64_t totalMatches = 0;
    for(int64_t i=0; i<length; i++) {
        if(hap1String[i] == hap2String[i]) {
            totalMatches++;
        }
    }
    return ((double)totalMatches)/length;
}

double getIdentityBetweenHaplotypesExcludingIndels(uint64_t *hap1String, uint64_t *hap2String, int64_t length) {
    /*
     * Returns the fraction of positions in two haplotypes that are identical.
     */
    int64_t totalMatches = 0;
    int64_t numGaps = 0;
    for(int64_t i=0; i<length; i++) {
        if(hap1String[i] == hap2String[i]) {
            totalMatches++;
        } else if (hap1String[i] == ALPHABET_SIZE-1 || hap2String[i] == ALPHABET_SIZE-1) {
            numGaps++;
        }
    }
    return ((double)totalMatches)/(length - numGaps);
}

void getExpectedMatcheBetweenProfileSeqs(stProfileSeq *pSeq1, stProfileSeq *pSeq2, int64_t *totalAlignedPositions, double *totalExpectedMatches) {
    /*
     * Calculates the number of base overlaps and expected base matches between two profile sequences.
     */

    for(int64_t i=0; i<pSeq1->length; i++) {
        // Establish if the coordinate is in both sequences
        int64_t j = i + pSeq1->refStart - pSeq2->refStart;
        if(j >= 0 && j < pSeq2->length) {
            (*totalAlignedPositions)++;

            // Calculate expectation of match
            for(int64_t k=0; k<ALPHABET_SIZE; k++) {
                double e1 = getProb(&(pSeq1->profileProbs[i * ALPHABET_SIZE]), k);
                double e2 = getProb(&(pSeq2->profileProbs[j * ALPHABET_SIZE]), k);
                assert(e1 * e2 <= 1.0);
                *totalExpectedMatches += e1 * e2;
            }
        }
    }
}

void printAvgIdentityBetweenProfileSequences(FILE *fH, stList *profileSequences, int64_t maxSeqs) {
    /*
     * Prints the average base identity between pairwise base overlaps between the given profile sequences
     */

    double totalExpectedMatches = 0.0;
    int64_t totalAlignedPositions = 0;

    for(int64_t i=0; i<stList_length(profileSequences) && i<maxSeqs; i++) {
        for(int64_t j=i+1; j<stList_length(profileSequences) && j<maxSeqs; j++) {
            getExpectedMatcheBetweenProfileSeqs(stList_get(profileSequences, i), stList_get(profileSequences, j),
                                                &totalAlignedPositions, &totalExpectedMatches);

        }
    }

    fprintf(fH, "Avg. pairwise identity between profile sequences: %f measured at %" PRIi64 " overlapping sites\n",
            totalExpectedMatches/totalAlignedPositions, totalAlignedPositions);
}

double *getHaplotypeBaseComposition(uint64_t *hapString, int64_t length) {
    /*
     * Get the count of each alphabet character in the haplotype sequence, returned
     * as an array.
     */
    double *baseCounts = st_calloc(ALPHABET_SIZE, sizeof(double));
    for(int64_t i=0; i<length; i++) {
        baseCounts[hapString[i]] += 1;
    }
    return baseCounts;
}

double *getExpectedProfileSequenceBaseComposition(stSet *profileSeqs) {
    /*
     * Get the expected count of each alphabet character in the profile sequences, returned
     * as an array.
     */
    double *baseCounts = st_calloc(ALPHABET_SIZE, sizeof(double));
    stSetIterator *it = stSet_getIterator(profileSeqs);
    stProfileSeq *pSeq;
    while((pSeq = stSet_getNext(it)) != NULL) {
        for(int64_t i=0; i<pSeq->length; i++) {
            for(int64_t j=0; j<ALPHABET_SIZE; j++) {
                baseCounts[j] += getProb(&(pSeq->profileProbs[i*ALPHABET_SIZE]), j);
            }
        }
    }
    return baseCounts;
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
        if (pos > pSeq->refStart && pos < pSeq->refStart+pSeq->length) {
            for(int64_t j=0; j<ALPHABET_SIZE; j++) {
                baseCounts[j] += getProb(&(pSeq->profileProbs[(pos - pSeq->refStart)*ALPHABET_SIZE]), j);
            }
        }
    }
    return baseCounts;
}

double *getColumnBaseComposition(stRPColumn *column, int64_t pos) {
    /*
     * Get the observed counts for each base seen at a particular position in a column
     */
    double *baseCounts = st_calloc(ALPHABET_SIZE, sizeof(double));
    for (int64_t i=0; i<column->depth; i++) {
        stProfileSeq *seq = column->seqHeaders[i];

        if (pos >= seq->refStart && pos < seq->length+seq->refStart) {
            for(int64_t j=0; j<ALPHABET_SIZE; j++) {
                baseCounts[j] += getProb(&(seq->profileProbs[(pos - seq->refStart) * ALPHABET_SIZE]), j);
            }
        }
    }
    return baseCounts;

}

void printBaseComposition(FILE *fH, double *baseCounts) {
    /*
     * Print the counts/fraction of each alphabet character.
     */
    double totalCount = 0;
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        totalCount += baseCounts[i];
    }
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        fprintf(fH, "Base %" PRIi64 " count: %f fraction: %f\n", i, baseCounts[i], baseCounts[i]/totalCount);
    }
}

void printBaseComposition2(FILE *fH, double *baseCounts) {
    /*
     * Print the counts/fraction of each alphabet character in a slightly more compressed form.
     */
    double totalCount = 0;
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        totalCount += baseCounts[i];
    }
    fprintf(fH, "\t\t0 (A)\t1 (C)\t2 (G)\t3 (T)\t4 (-) \n");
    fprintf(fH, "    Counts:");
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        fprintf(fH, "\t%0.1f", baseCounts[i]);
    }
    fprintf(fH, "\n    Fraction:");
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        fprintf(fH, "\t%0.3f", baseCounts[i]/totalCount);
    }
    fprintf(fH, "\n");
}


void printColumnAtPosition(stRPHmm *hmm, int64_t pos) {
    /*
     * Print out the columns of the hmm at a specific position
     */

    //stRPHmm_print(hmm, stderr, 1, 0);
    stRPColumn *column = hmm->firstColumn;
    while(1) {
        if (pos >= column->refStart && pos < column->refStart+column->length) {
            double *columnBaseCounts = getColumnBaseComposition(column, pos);
            printBaseComposition2(stderr, columnBaseCounts);
            free(columnBaseCounts);
        }

        if (column->nColumn == NULL) {
            break;
        }
        column = column->nColumn->nColumn;
    }
}

void test_jsmnParsing(CuTest *testCase) {

    char *paramsFile = "../tests/parsingTest.json";

    stBaseMapper *baseMapper = stBaseMapper_construct();
    stRPHmmParameters *params = parseParameters(paramsFile, baseMapper);

    // Check that alphabet was parsed as expected, and that conversions
    // between types of bases work properlu
    CuAssertIntEquals(testCase, baseMapper->size, 5);
    CuAssertStrEquals(testCase, baseMapper->wildcard, "Nn");

    // Check that numerical bases are mapped to characters correctly
    CuAssertIntEquals(testCase, stBaseMapper_getCharForValue(baseMapper, 0), 'A');
    CuAssertIntEquals(testCase, stBaseMapper_getCharForValue(baseMapper, 1), 'C');
    CuAssertIntEquals(testCase, stBaseMapper_getCharForValue(baseMapper, 2), 'G');
    CuAssertIntEquals(testCase, stBaseMapper_getCharForValue(baseMapper, 3), 'T');
    CuAssertIntEquals(testCase, stBaseMapper_getCharForValue(baseMapper, 4), '-');


    // Check that character bases are mapped to numbers correctly
    CuAssertIntEquals(testCase, stBaseMapper_getValueForChar(baseMapper, 'A'), 0);
    CuAssertIntEquals(testCase, stBaseMapper_getValueForChar(baseMapper, 'a'), 0);
    CuAssertIntEquals(testCase, stBaseMapper_getValueForChar(baseMapper, 'C'), 1);
    CuAssertIntEquals(testCase, stBaseMapper_getValueForChar(baseMapper, 'c'), 1);
    CuAssertIntEquals(testCase, stBaseMapper_getValueForChar(baseMapper, 'G'), 2);
    CuAssertIntEquals(testCase, stBaseMapper_getValueForChar(baseMapper, 'g'), 2);
    CuAssertIntEquals(testCase, stBaseMapper_getValueForChar(baseMapper, 'T'), 3);
    CuAssertIntEquals(testCase, stBaseMapper_getValueForChar(baseMapper, 't'), 3);
    CuAssertIntEquals(testCase, stBaseMapper_getValueForChar(baseMapper, '-'), 4);

    // Check stRPHmmParameters

    // Check haplotype substitution model and error model parsed correctly
    // and that the proper values were set in the model parameters
    double delta = 0.0001;
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            if (i < 4 && j < 4) {
                if (i == j) {
                    CuAssertDblEquals(testCase, params->hetSubModelSlow[i*5+j], log(0.998), delta);
                    CuAssertDblEquals(testCase, params->hetSubModel[i*5+j], scaleToLogIntegerSubMatrix(log(0.998)), delta);
                    CuAssertDblEquals(testCase, params->readErrorSubModelSlow[i*5+j], log(0.9), delta);
                    CuAssertDblEquals(testCase, params->readErrorSubModel[i*5+j], scaleToLogIntegerSubMatrix(log(0.9)), delta);
                } else {
                    CuAssertDblEquals(testCase, params->hetSubModelSlow[i*5+j], log(0.000333), delta);
                    CuAssertDblEquals(testCase, params->hetSubModel[i*5+j], scaleToLogIntegerSubMatrix(log(0.000333)), delta);
                    CuAssertDblEquals(testCase, params->readErrorSubModelSlow[i*5+j], log(0.01), delta);
                    CuAssertDblEquals(testCase, params->readErrorSubModel[i*5+j], scaleToLogIntegerSubMatrix(log(0.01)), delta);
                }
            }
            else if (i != j) {
                CuAssertDblEquals(testCase, params->hetSubModelSlow[i*5+j], log(0.001), delta);
                CuAssertDblEquals(testCase, params->hetSubModel[i*5+j], scaleToLogIntegerSubMatrix(log(0.001)), delta);
                CuAssertDblEquals(testCase, params->readErrorSubModelSlow[i*5+j], log(0.07), delta);
                CuAssertDblEquals(testCase, params->readErrorSubModel[i*5+j], scaleToLogIntegerSubMatrix(log(0.07)), delta);
            } else {
                CuAssertDblEquals(testCase, params->hetSubModelSlow[i*5+j], log(0.996), delta);
                CuAssertDblEquals(testCase, params->hetSubModel[i*5+j], scaleToLogIntegerSubMatrix(log(0.996)), delta);
                CuAssertDblEquals(testCase, params->readErrorSubModelSlow[i*5+j], log(0.72), delta);
                CuAssertDblEquals(testCase, params->readErrorSubModel[i*5+j], scaleToLogIntegerSubMatrix(log(0.72)), delta);
            }
        }
    }
    // Check remaining parameters parsed correctly
    CuAssertTrue(testCase, params->maxNotSumTransitions);
    CuAssertIntEquals(testCase, params->maxPartitionsInAColumn, 50);
    CuAssertIntEquals(testCase, params->maxCoverageDepth, 64);
    CuAssertIntEquals(testCase, params->minReadCoverageToSupportPhasingBetweenHeterozygousSites, 4);

    // cleanup
    stBaseMapper_destruct(baseMapper);
    stRPHmmParameters_destruct(params);
}

/*
 * Stores information about relevant test results.
 */
struct genotypeResults {
    int64_t matches;
    int64_t totalRefRecords;
    int64_t negatives;
    int64_t matchedGaps;
    int64_t error_trueVariantWrong;
    int64_t error_badPartition;
    int64_t falsePositives;
    int64_t falsePositiveGaps;
};

void compareVCFs(FILE *fh, stList *hmms,
                 char *vcf_toEval, char *vcf_ref, int64_t refStart, int64_t refEnd,
                 stBaseMapper *baseMapper, struct genotypeResults *results, bool printVerbose) {
    /*
     * Test to compare a vcf to a truth vcf containing known variants for the region.
     * This test currently requires knowledge of the specific interval that should be looked at
     *
     * Test depends on the format of the vcf files written in vcfWriter.c
     * (Currently they don't follow a quite standard format)
     *
     * Information about some of the results saved in the genotypeResults struct
     */
    fprintf(fh, "VCF reference: %s \n", vcf_ref);
    fprintf(fh, "VCF being evaluated: %s \n", vcf_toEval);

    vcfFile *inRef = vcf_open(vcf_ref,"r"); //open vcf file
    if (inRef == NULL) {
        fprintf(fh, "ERROR: cannot open reference vcf, %s\n", vcf_ref);
        return;
    }
    bcf_hdr_t *hdrRef = bcf_hdr_read(inRef); //read header
    bcf1_t *refRecord = bcf_init1(); //initialize for reading

    vcfFile *inEval = vcf_open(vcf_toEval,"r"); //open vcf file
    if (inEval == NULL) {
        fprintf(fh, "ERROR: cannot open vcf to evaluate, %s\n", vcf_toEval);
        return;
    }
    bcf_hdr_t *hdrEval = bcf_hdr_read(inEval); //read header
    bcf1_t *evalRecord = bcf_init1(); //initialize for reading
    int64_t referencePos = 0;

    fprintf(fh, "Writing out false negatives \n");

    // Start by looking at the first hmm
    int64_t hmmIndex = 0;
    stRPHmm *hmm = stList_get(hmms, hmmIndex);
    // Inefficient, but recalculate the info relevant to the hmm to get bipartitions
    // TODO: how to save this info from earlier?
    stRPHmm_forwardBackward(hmm);
    stList *path = stRPHmm_forwardTraceBack(hmm);
    stGenomeFragment *gF = stGenomeFragment_construct(hmm, path);
    stSet *reads1 = stRPHmm_partitionSequencesByStatePath(hmm, path, 1);
    stSet *reads2 = stRPHmm_partitionSequencesByStatePath(hmm, path, 0);

    // Iterate through the vcf being checked until getting to the start of the specified interval
    // Don't bother analyzing these records
    int64_t evalPos =  0;
    bcf1_t *unpackedRecord;
    while (evalPos < refStart - 1 && (bcf_read(inEval, hdrEval, evalRecord) == 0)) {
        unpackedRecord = evalRecord;
        bcf_unpack(unpackedRecord, BCF_UN_INFO);
        evalPos = unpackedRecord->pos;
    }

    while(bcf_read(inRef, hdrRef, refRecord) == 0) {
        // Unpack reference record
        bcf1_t *unpackedRecordRef = refRecord;
        bcf_unpack(unpackedRecordRef, BCF_UN_INFO);
        referencePos = unpackedRecordRef->pos;

        // Make sure to only look at records in the specified interval
        if (referencePos < refStart || referencePos < hmm->refStart) continue;
        if (referencePos > refEnd) break;

        // If the position is beyond the end of this hmm, get the next one
        if ((hmm->refStart + hmm->refLength) < referencePos) {
            hmmIndex++;
            if (hmmIndex < stList_length(hmms)) {
                hmm = stList_get(hmms, hmmIndex);
                path = stRPHmm_forwardTraceBack(hmm);
                gF = stGenomeFragment_construct(hmm, path);
                reads1 = stRPHmm_partitionSequencesByStatePath(hmm, path, 1);
                reads2 = stRPHmm_partitionSequencesByStatePath(hmm, path, 0);
            } else {
                break;
            }
        }

        results->totalRefRecords++;
        char *refChar = unpackedRecordRef->d.als;
        char *altRefChar = unpackedRecordRef->d.allele[1];

        // Iterate through vcf until getting to the position of the variant
        // from the reference vcf currently being looked at
        while (evalPos < referencePos && (bcf_read(inEval, hdrEval, evalRecord) == 0)) {
            unpackedRecord = evalRecord;                // unpack record
            bcf_unpack(unpackedRecord, BCF_UN_INFO);
            evalPos = unpackedRecord->pos;
            if (evalPos < referencePos) {
                char *evalRefChar = unpackedRecord->d.als;
                char *altEvalChar1 = unpackedRecord->d.allele[1];
                char *altEvalChar2 = unpackedRecord->d.allele[2];

                if (strcmp(altEvalChar1, altEvalChar2) == 0 && strcmp(evalRefChar, altEvalChar1) == 0) results->negatives++;
                else if ((strcmp(altEvalChar1, "-") == 0) && (strcmp("-", altEvalChar2) == 0)) results->matchedGaps++;
                else {
                    // False positive
                    results->falsePositives++;
                    if ((strcmp(altEvalChar1, "-") == 0) || (strcmp("-", altEvalChar2) == 0)) {
                        results->falsePositiveGaps++;
                    }
                    else fprintf(fh, "*FP  pos: %" PRIi64 " ref:%s alt1: %s alt2: %s \n",
                            evalPos, evalRefChar, altEvalChar1, altEvalChar2);
                }
            }
        }
        // At locus of known variation
        if (evalPos == referencePos) {
            char *evalRefChar = unpackedRecord->d.als;
            char *altEvalChar1 = unpackedRecord->d.allele[1];
            char *altEvalChar2 = unpackedRecord->d.allele[2];
            bool missedVariant = false;
            char *missingChar;

            if (strcmp(altRefChar, altEvalChar1) == 0 && strcmp(altRefChar, altEvalChar2) == 0) {
                // Both match alternate - no variation was found
                missedVariant = true;
                missingChar = refChar;
            } else if (strcmp(refChar, altEvalChar1) == 0 && strcmp(refChar, altEvalChar2) == 0) {
                // Both match reference - no variation was found
                missedVariant = true;
                missingChar = altRefChar;
            } else if (strcmp(altRefChar, altEvalChar1) == 0 || strcmp(altRefChar, altEvalChar2) == 0) {
                // Only one matches - variation found that matches the reference!
                results->matches++;
            } else {
                // Leftovers - don't fit into the other categories
                // TODO: classify these
                missedVariant = true;
                missingChar = altRefChar;
            }
            if (missedVariant) {
                // False negative - no variation was found, but truth vcf has one
                // Print out info about the record
                fprintf(fh, "pos: %" PRIi64 "\n\tref: %s\talt: ", referencePos, refChar);
                for (int i = 1; i < unpackedRecordRef->n_allele; i++) {
                    if (i != 1) fprintf(fh, ",");
                    fprintf(fh, "%s", unpackedRecordRef->d.allele[i]);
                }
                fprintf(fh, "\n\toutput alleles: ");
                for (int i = 1; i < unpackedRecord->n_allele; i++) {
                    if (i != 1) fprintf(fh, ",");
                    fprintf(fh, "%s", unpackedRecord->d.allele[i]);
                }
                fprintf(fh, "\n");

                double *read1BaseCounts = getProfileSequenceBaseCompositionAtPosition(reads1, evalPos);
                double *read2BaseCounts = getProfileSequenceBaseCompositionAtPosition(reads2, evalPos);

                if (printVerbose) {
                    printColumnAtPosition(hmm, evalPos);
                    fprintf(fh, "\tPartition 1: \n");
                    printBaseComposition2(stderr, read1BaseCounts);
                    fprintf(fh, "\tPartition 2: \n");
                    printBaseComposition2(stderr, read2BaseCounts);
                }

                fprintf(fh, "\tposterior prob: %f\n", unpackedRecord->qual);

                // Quantify the type of false negative it was
                double totalCount = 0;
                for (int64_t i = 0; i < ALPHABET_SIZE; i++) {
                    totalCount += read1BaseCounts[i];
                    totalCount += read2BaseCounts[i];
                }
                int missingBase = stBaseMapper_getValueForChar(baseMapper, *missingChar);
                float fractionMissingBase =
                        (read1BaseCounts[missingBase] + read2BaseCounts[missingBase]) / totalCount;
                fprintf(fh, "\tfraction of 'missing' base seen in reads: %f\n", fractionMissingBase);
                if (fractionMissingBase < 0.10) {       // TODO: threshold arbitrarily chosen
                    results->error_trueVariantWrong++;
                    fprintf(fh, "\tTRUE VARIANT WRONG\n");
                } else {
                    results->error_badPartition++;
                    fprintf(fh, "\tBAD PARTITION\n");
                }
                free(read1BaseCounts);
                free(read2BaseCounts);
            }
        }
    }

    // Look through remaining positions after last known variation in reference
    while (evalPos < refEnd && evalPos > refStart) {

        if (evalPos > (hmm->refStart+hmm->refLength)) {
            break;      // If the hmm didn't cover this region don't look through the rest
        }
        bcf1_t *unpackedRecord;
        if (bcf_read(inEval, hdrEval, evalRecord) == 0) {
            unpackedRecord = evalRecord;
            bcf_unpack(unpackedRecord, BCF_UN_STR);
            evalPos = unpackedRecord->pos;
            char *evalRefChar = unpackedRecord->d.als;
            char *altEvalChar1 = unpackedRecord->d.allele[1];
            char *altEvalChar2 = unpackedRecord->d.allele[2];

            if (strcmp(altEvalChar1, altEvalChar2) == 0 && strcmp(evalRefChar, altEvalChar1) == 0) results->negatives++;
            else if ((strcmp(altEvalChar1, "-") == 0) && (strcmp("-", altEvalChar2) == 0)) results->matchedGaps++;
            else {
                // False positive
                results->falsePositives++;
                if ((strcmp(altEvalChar1, "-") == 0) || (strcmp("-", altEvalChar2) == 0)) {
                    results->falsePositiveGaps++;
                }
                else fprintf(fh, "*FP  pos: %" PRIi64 " ref:%s alt1: %s alt2: %s \n",
                        evalPos, evalRefChar, altEvalChar1, altEvalChar2);
            }
        } else {
            break; // Can't read the next line (probably reached end of what vcf has written so far)
        }
    }

    // cleanup
    vcf_close(inRef);
    vcf_close(inEval);
    bcf_hdr_destroy(hdrRef);
    bcf_hdr_destroy(hdrEval);
    bcf_destroy(refRecord);
    bcf_destroy(evalRecord);
    stSet_destruct(reads1);
    stSet_destruct(reads2);
    stGenomeFragment_destruct(gF);
    stList_destruct(path);
}


void genotypingTest(FILE *fh, char *paramsFile, char *bamFile, char *vcfOutFile, char *vcfOutFileDiff, char *referenceFile, char *vcfReference, int64_t refStart, int64_t refEnd, bool printVerbose) {
    fprintf(stderr, "> Parsing parameters\n");
    stBaseMapper *baseMapper = stBaseMapper_construct();
    stRPHmmParameters *params = parseParameters(paramsFile, baseMapper);
    // Print a report of the parsed parameters
    stRPHmmParameters_printParameters(params, stderr);

    fprintf(stderr, "> Creating profile sequences\n");
    stList *profileSequences = stList_construct();
    parseReads(profileSequences, bamFile, baseMapper);

    // Print some stats about the input sequences
    printSequenceStats(stderr, profileSequences);
    printAvgIdentityBetweenProfileSequences(stderr, profileSequences, 100);

    fprintf(stderr, "> Building hmms\n");
    stList *hmms = getRPHmms(profileSequences, params);
    stList *l = stList_construct3(0, (void (*)(void *))stRPHmm_destruct2);
    while(stList_length(hmms) > 0) {
        stList_appendAll(l, stRPHMM_splitWherePhasingIsUncertain(stList_pop(hmms)));
    }
    hmms = l;
    fprintf(stderr, "\tGot %" PRIi64 " hmms\n", stList_length(hmms));

    fprintf(stderr, "> Writing vcf files\n");
    vcfFile *vcfOutFP = vcf_open(vcfOutFile, "w");
    bcf_hdr_t *bcf_hdr = writeVcfHeader(vcfOutFP, l, referenceFile);

    vcfFile *vcfOutFP_diff = vcf_open(vcfOutFileDiff, "w");
    bcf_hdr_t *bcf_hdr_diff = writeVcfHeader(vcfOutFP_diff, l, referenceFile);

    // TODO: add better constructor
    struct genotypeResults *results = st_calloc(8, sizeof(int64_t));

    for(int64_t i=0; i<stList_length(hmms); i++) {
        stRPHmm *hmm = stList_get(hmms, i);

        // Print stats about the HMM
        fprintf(stderr, "The %" PRIi64 "th hmm\n", i);
        stRPHmm_print(hmm, stderr, 0, 0);

        // Run the forward-backward algorithm
        stRPHmm_forwardBackward(hmm);

        // Now compute a high probability path through the hmm
        stList *path = stRPHmm_forwardTraceBack(hmm);

        // Compute the genome fragment
        stGenomeFragment *gF = stGenomeFragment_construct(hmm, path);

        // Get partitioned sequences
        stSet *reads1 = stRPHmm_partitionSequencesByStatePath(hmm, path, 1);
        stSet *reads2 = stRPHmm_partitionSequencesByStatePath(hmm, path, 0);

        if (printVerbose) {
            fprintf(stderr, "\nThere are %" PRIi64 " reads covered by the hmm, bipartitioned into sets of %" PRIi64 " and %" PRIi64 " reads\n",
                    stList_length(hmm->profileSeqs), stSet_size(reads1), stSet_size(reads2));

            // Print the similarity between the two imputed haplotypes sequences
            fprintf(stderr, "The haplotypes have an %f identity\n", getIdentityBetweenHaplotypes(gF->haplotypeString1, gF->haplotypeString2, gF->length));
            fprintf(stderr, "\tIdentity excluding indels: %f \n\n", getIdentityBetweenHaplotypesExcludingIndels(gF->haplotypeString1, gF->haplotypeString2, gF->length));

            // Print the base composition of the haplotype sequences
            double *hap1BaseCounts = getHaplotypeBaseComposition(gF->haplotypeString1, gF->length);
            fprintf(stderr, "The base composition of haplotype 1:\n");
            printBaseComposition(stderr, hap1BaseCounts);
            free(hap1BaseCounts);

            double *hap2BaseCounts = getHaplotypeBaseComposition(gF->haplotypeString2, gF->length);
            fprintf(stderr, "The base composition of haplotype 2:\n");
            printBaseComposition(stderr, hap2BaseCounts);
            free(hap2BaseCounts);

            // Print the base composition of the reads
            double *reads1BaseCounts =getExpectedProfileSequenceBaseComposition(reads1);
            fprintf(stderr, "The base composition of reads1 set:\n");
            printBaseComposition(stderr, reads1BaseCounts);
            free(reads1BaseCounts);

            double *reads2BaseCounts =getExpectedProfileSequenceBaseComposition(reads2);
            fprintf(stderr, "The base composition of reads2 set:\n");
            printBaseComposition(stderr, reads2BaseCounts);
            free(reads2BaseCounts);

            fprintf(stderr, "Genome fragment info: refStart = %" PRIi64 ", length = %" PRIi64 "\n", gF->refStart, gF->length);

            // Print some summary stats about the differences between haplotype sequences and the bipartitioned reads
            fprintf(stderr, "hap1 vs. reads1 identity: %f\n",
                    getExpectedIdentity(gF->haplotypeString1, gF->refStart, gF->length, reads1));
            fprintf(stderr, "hap1 vs. reads2 identity: %f\n",
                    getExpectedIdentity(gF->haplotypeString1, gF->refStart, gF->length, reads2));

            fprintf(stderr, "hap2 vs. reads2 identity: %f\n",
                    getExpectedIdentity(gF->haplotypeString2, gF->refStart, gF->length, reads2));
            fprintf(stderr, "hap2 vs. reads1 identity: %f\n",
                    getExpectedIdentity(gF->haplotypeString2, gF->refStart, gF->length, reads1));
        }


        writeVcfFragment(vcfOutFP, bcf_hdr, gF, referenceFile, baseMapper, false);
        writeVcfFragment(vcfOutFP_diff, bcf_hdr_diff, gF, referenceFile, baseMapper, true);

        // cleanup for this hmm
        stSet_destruct(reads1);
        stSet_destruct(reads2);
        stGenomeFragment_destruct(gF);
        stList_destruct(path);
    }
    compareVCFs(stderr, hmms, vcfOutFile, vcfReference,
                refStart, refEnd, baseMapper, results, 1);

    fprintf(stderr, "\nStats for all %" PRIi64 " hmms:\n", stList_length(hmms));

    int64_t regionLength = refEnd - refStart + 1;
    fprintf(fh, "\nSensitivity: %f \n(= fraction of true positives compared to reference, \t%" PRIi64 " out of %"PRIi64 ")\n",
            (float)results->matches/results->totalRefRecords, results->matches, results->totalRefRecords) ;
    fprintf(fh, "\t \t(Number of false negatives: %" PRIi64 ")\n", results->totalRefRecords-results->matches);
    fprintf(fh, "Specificity: %f \n(= fraction of true negatives compared to reference, \t%" PRIi64 " out of % "PRIi64 ")\n",
            (float)results->negatives/(regionLength-results->totalRefRecords), results->negatives, regionLength-results->totalRefRecords);
    fprintf(fh, "Number of deletions matched in both haplotypes: %"PRIi64 "\n", results->matchedGaps);
    fprintf(fh, "False positives found: %" PRIi64 " \n", results->falsePositives);
    fprintf(fh, "False positives that included a gap: %" PRIi64 " \n", results->falsePositiveGaps);
    fprintf(fh, "False negatives where true variant appears to be wrong: %" PRIi64 " \t(%f)\n", results->error_trueVariantWrong, (float)results->error_trueVariantWrong/(results->error_trueVariantWrong+results->error_badPartition));
    fprintf(fh, "False negatives where partition appears to be bad: %" PRIi64 " \t\t(%f)\n", results->error_badPartition, (float)results->error_badPartition/(results->error_trueVariantWrong+results->error_badPartition));

    // cleanup
    free(results);

    stRPHmmParameters_destruct(params);
    stBaseMapper_destruct(baseMapper);
    stList_destruct(hmms);
    stList_destruct(profileSequences);

    vcf_close(vcfOutFP);
    vcf_close(vcfOutFP_diff);
    bcf_hdr_destroy(bcf_hdr);
    bcf_hdr_destroy(bcf_hdr_diff);
}



void test_5kbGenotyping(CuTest *testCase) {

    char *paramsFile = "../tests/params.json";
    char *referenceFile = "../tests/hg19.chr3.9mb.fa";
    char *vcfOutFile = "test_5kb.vcf";
    char *vcfOutFileDiff = "test_5kb_diff.vcf";

    char *bamFile = "../tests/NA12878.pb.chr3.5kb.bam";
    char *vcfReference = "../tests/HG001.GRCh37.chr3.100kb.vcf";

    fprintf(stderr, "Testing haplotype inference on %s\n", bamFile);

    genotypingTest(stderr, paramsFile, bamFile, vcfOutFile, vcfOutFileDiff, referenceFile, vcfReference, 190000, 195000, 1);

    // TODO: create vcf specifically for this 5 kb region
    //compareVCFs(vcfOutFile, vcfReference, 150000, 155000);

}

void test_100kbGenotyping(CuTest *testCase) {

    char *paramsFile = "../tests/params.json";
    char *referenceFile = "../tests/hg19.chr3.9mb.fa";
    char *vcfOutFile = "test_100kb.vcf";
    char *vcfOutFileDiff = "test_100kb_diff.vcf";

    char *bamFile = "../tests/NA12878.pb.chr3.100kb.0.bam";
    char *vcfReference = "../tests/NA12878.PG.chr3.100kb.0.vcf";

    //bamFile = "../examples/KRAS_chr3.bam";

    fprintf(stderr, "Testing haplotype inference on %s\n", bamFile);
    genotypingTest(stderr, paramsFile, bamFile, vcfOutFile, vcfOutFileDiff, referenceFile, vcfReference, 100000, 200000, 1);
    //compareVCFs(vcfOutFile, vcfReference, 100000, 200000);
}




CuSuite *marginPhaseTestSuite(void) {

    st_setLogLevelFromString("debug");

    CuSuite* suite = CuSuiteNew();

    //SUITE_ADD_TEST(suite, test_jsmnParsing);
    //SUITE_ADD_TEST(suite, test_5kbGenotyping);
    SUITE_ADD_TEST(suite, test_100kbGenotyping);

    return suite;
}
