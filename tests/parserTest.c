/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "margin.h"
#include "htsIntegration.h"

/*
 * Test that a bam file is read in correctly.
 * Checks:
 * - Number of reads is correct.
 * - Characters in a read are read correctly.
 * - Profile sequence probabilities are updated correctly.
 * - Soft clipping is handled
 * - Insertions and deletions are handled
 */
void test_bamReadParsing(CuTest *testCase) {

    char *bamFile = "../tests/data/NA12878.pb.chr3.5kb.bam";
    char *paramsFile = "../tests/data/parsingTest.json";

    FILE *fh = fopen(paramsFile, "rb");
    Params *fullParams = params_readParams2(fh, FALSE, TRUE);
    fclose(fh);
    stBaseMapper *baseMapper = fullParams->baseMapper;
    stRPHmmParameters *params = fullParams->phaseParams;

    // Parse reads for interval
    st_logInfo("> Parsing input reads\n");
    stList *profileSequences = stList_construct3(0, (void (*)(void *))stProfileSeq_destruct);
    int64_t readCount = parseReads(profileSequences, bamFile, baseMapper, params);

    // Test to see number of reads is correct - should be 149
    CuAssertIntEquals(testCase, readCount, 149);

    // Test to see first read is being parsed correctly
    samFile *in = hts_open(bamFile, "r");
    bam_hdr_t *bamHdr = sam_hdr_read(in);
    bam1_t *aln = bam_init1();

    // Get first read
    CuAssertTrue(testCase, sam_read1(in,bamHdr,aln) > 0);

    int64_t pos = aln->core.pos; //left most position of alignment
    char *chr = bamHdr->target_name[aln->core.tid] ; //contig name (chromosome)
    int64_t len = aln->core.l_qseq; //length of the read.
    uint8_t *seq = bam_get_seq(aln);  // DNA sequence
    char *readName = bam_get_qname(aln);
    uint32_t *cigar = bam_get_cigar(aln);


    int64_t start_read = 0;
    int64_t end_read = 0;
    int64_t start_ref = pos;
    int64_t cig_idx = 0;

    // Find the correct starting locations on the read and reference sequence,
    // to deal with things like inserts / deletions / soft clipping
    while(cig_idx < aln->core.n_cigar) {
        int cigarOp = cigar[cig_idx] & BAM_CIGAR_MASK;
        int cigarNum = cigar[cig_idx] >> BAM_CIGAR_SHIFT;

        if (cigarOp == BAM_CMATCH || cigarOp == BAM_CEQUAL || cigarOp==BAM_CDIFF) {
            break;
        }
        else if (cigarOp == BAM_CDEL || cigarOp == BAM_CREF_SKIP) {
            start_ref += cigarNum;
            cig_idx++;
        } else if (cigarOp == BAM_CINS || cigarOp == BAM_CSOFT_CLIP) {
            start_read += cigarNum;
            cig_idx++;
        } else if (cigarOp == BAM_CHARD_CLIP || cigarOp == BAM_CPAD) {
            cig_idx++;
        } else {
            st_logCritical("Unidentifiable cigar operation\n");
        }
    }
    // Check for soft clipping at the end
    int lastCigarOp = cigar[aln->core.n_cigar-1] & BAM_CIGAR_MASK;
    int lastCigarNum = cigar[aln->core.n_cigar-1] >> BAM_CIGAR_SHIFT;
    if (lastCigarOp == BAM_CSOFT_CLIP) {
        end_read += lastCigarNum;
    }
    // Count number of insertions & deletions in sequence
    int64_t numInsertions = 0;
    int64_t numDeletions = 0;
    countIndels(cigar, aln->core.n_cigar, &numInsertions, &numDeletions);
    int64_t trueLength = len-start_read-end_read+numDeletions-numInsertions;

    // Create empty profile sequence
    stProfileSeq *pSeq = stProfileSeq_constructEmptyProfile(chr, readName, pos, trueLength);

    // Variables to keep track of position in sequence / cigar operations
    cig_idx = 0;
    int64_t currPosInOp = 0;
    int64_t cigarOp = -1;
    int64_t cigarNum = -1;
    int64_t idxInSeq = start_read;
    char *firstMatches = st_calloc(18, sizeof(char));
    char *lastMatches = st_calloc(29, sizeof(char));
    // For each position turn character into profile probability
    // As is, this makes the probability 1 for the base read in, and 0 otherwise
    for (uint32_t i = 0; i < trueLength; i++) {

        if (currPosInOp == 0) {
            cigarOp = cigar[cig_idx] & BAM_CIGAR_MASK;
            cigarNum = cigar[cig_idx] >> BAM_CIGAR_SHIFT;
        }
        // Build first and last sequences to check
        if (cig_idx == 1) {
            char c = seq_nt16_str[bam_seqi(seq, idxInSeq)];
            firstMatches[currPosInOp] = c;
        }
        else if (cig_idx == aln->core.n_cigar-2) {
            char c = seq_nt16_str[bam_seqi(seq, idxInSeq)];
            lastMatches[currPosInOp] = c;
        }

        if (cigarOp == BAM_CMATCH || cigarOp == BAM_CEQUAL || cigarOp == BAM_CDIFF) {
            int64_t b = stBaseMapper_getValueForChar(baseMapper, seq_nt16_str[bam_seqi(seq, idxInSeq)]);
            pSeq->profileProbs[i * ALPHABET_SIZE + b] = ALPHABET_MAX_PROB;
            idxInSeq++;
        }
        else if (cigarOp == BAM_CDEL || cigarOp == BAM_CREF_SKIP) {
            int64_t b = stBaseMapper_getValueForChar(baseMapper, '-');
            pSeq->profileProbs[i * ALPHABET_SIZE + b] = ALPHABET_MAX_PROB;
        } else if (cigarOp == BAM_CINS) {
            // Currently, ignore insertions
            idxInSeq++;
            i--;
        } else if (cigarOp == BAM_CSOFT_CLIP || cigarOp == BAM_CHARD_CLIP || cigarOp == BAM_CPAD) {
            // nothing really to do here. skip to next cigar operation
            currPosInOp = cigarNum-1;
            i--;
        } else {
            st_logCritical("Unidentifiable cigar operation\n");
        }

        currPosInOp++;
        if (currPosInOp == cigarNum) {
            cig_idx++;
            currPosInOp = 0;
        }
    }

    // Cigar string
    CuAssertIntEquals(testCase, aln->core.n_cigar, 2887);
    CuAssertIntEquals(testCase, start_read, 76);
    CuAssertIntEquals(testCase, end_read, 22);
    CuAssertIntEquals(testCase, numInsertions, 1265);
    CuAssertIntEquals(testCase, numDeletions, 646);

    // Read info
    CuAssertIntEquals(testCase, trueLength, 14786);
    CuAssertStrEquals(testCase, firstMatches, "CCAGTGCTTGACTTACAT");
    CuAssertStrEquals(testCase, lastMatches, "CACTGCTCTTGTGTACCAGATAAGTAAAA");
    uint8_t *p = &pSeq->profileProbs[1 * ALPHABET_SIZE];
    double delta = 0.0001;
    CuAssertDblEquals(testCase, getProb(p, 0), 0.0, delta);
    CuAssertDblEquals(testCase, getProb(p, 1), 1.0, delta);
    CuAssertDblEquals(testCase, getProb(p, 2), 0.0, delta);
    CuAssertDblEquals(testCase, getProb(p, 3), 0.0, delta);
    CuAssertDblEquals(testCase, getProb(p, 4), 0.0, delta);

    // Cleanup
    params_destruct(fullParams);
    stList_destruct(profileSequences);
}


void test_jsmnParsing(CuTest *testCase) {

    char *paramsFiles[2] = {"../tests/data/parsingTest.json", "../tests/data/parsingTest.allParams.json"};

    for (int p = 0; p < 2; p++) {
        char *paramsFile = paramsFiles[p];

        FILE *fh = fopen(paramsFile, "rb");
        Params *fullParams = params_readParams2(fh, FALSE, TRUE);
        fclose(fh);
        stBaseMapper *baseMapper = fullParams->baseMapper;
        stRPHmmParameters *params = fullParams->phaseParams;

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
                        CuAssertDblEquals(testCase, params->hetSubModelSlow[i * 5 + j], log(0.998), delta);
                        CuAssertDblEquals(testCase, params->hetSubModel[i * 5 + j],
                                          scaleToLogIntegerSubMatrix(log(0.998)), delta);
                        CuAssertDblEquals(testCase, params->readErrorSubModelSlow[i * 5 + j], log(0.9), delta);
                        CuAssertDblEquals(testCase, params->readErrorSubModel[i * 5 + j],
                                          scaleToLogIntegerSubMatrix(log(0.9)), delta);
                    } else {
                        CuAssertDblEquals(testCase, params->hetSubModelSlow[i * 5 + j], log(0.000333), delta);
                        CuAssertDblEquals(testCase, params->hetSubModel[i * 5 + j],
                                          scaleToLogIntegerSubMatrix(log(0.000333)), delta);
                        CuAssertDblEquals(testCase, params->readErrorSubModelSlow[i * 5 + j], log(0.01), delta);
                        CuAssertDblEquals(testCase, params->readErrorSubModel[i * 5 + j],
                                          scaleToLogIntegerSubMatrix(log(0.01)), delta);
                    }
                } else if (i != j) {
                    CuAssertDblEquals(testCase, params->hetSubModelSlow[i * 5 + j], log(0.001), delta);
                    CuAssertDblEquals(testCase, params->hetSubModel[i * 5 + j], scaleToLogIntegerSubMatrix(log(0.001)),
                                      delta);
                    CuAssertDblEquals(testCase, params->readErrorSubModelSlow[i * 5 + j], log(0.07), delta);
                    CuAssertDblEquals(testCase, params->readErrorSubModel[i * 5 + j],
                                      scaleToLogIntegerSubMatrix(log(0.07)), delta);
                } else {
                    CuAssertDblEquals(testCase, params->hetSubModelSlow[i * 5 + j], log(0.996), delta);
                    CuAssertDblEquals(testCase, params->hetSubModel[i * 5 + j], scaleToLogIntegerSubMatrix(log(0.996)),
                                      delta);
                    CuAssertDblEquals(testCase, params->readErrorSubModelSlow[i * 5 + j], log(0.72), delta);
                    CuAssertDblEquals(testCase, params->readErrorSubModel[i * 5 + j],
                                      scaleToLogIntegerSubMatrix(log(0.72)), delta);
                }
            }
        }
        // Check remaining parameters parsed correctly
        CuAssertTrue(testCase, params->maxNotSumTransitions);
        CuAssertIntEquals(testCase, params->maxPartitionsInAColumn, 50);
        CuAssertIntEquals(testCase, params->maxCoverageDepth, 64);
        CuAssertIntEquals(testCase, params->minReadCoverageToSupportPhasingBetweenHeterozygousSites, 4);

        // cleanup
        params_destruct(fullParams);
    }
}

CuSuite *marginPhaseParserTestSuite(void) {
    st_setLogLevelFromString("debug");
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, test_jsmnParsing);
    SUITE_ADD_TEST(suite, test_bamReadParsing);

    return suite;
}
