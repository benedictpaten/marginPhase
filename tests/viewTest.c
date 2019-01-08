/*
 * Copyright (C) 2009-2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "margin.h"

int64_t calcSequenceMatches(char *seq1, char *seq2); // in polisherTest

static char *polishParamsFile = "../params/allParams.np.json";
#define TEST_POLISH_FILES_DIR "../tests/polishTestExamples/"

stList *makeAlignmentList(const int64_t *anchorPairs, int64_t anchorPairNo) {
    stList *alignmentPairs = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);

    for(int64_t i=0; i<anchorPairNo; i++) {
        stList_append(alignmentPairs, stIntTuple_construct3(1, anchorPairs[i*2], anchorPairs[i*2+1]));
    }

    return alignmentPairs;
}

void test_view(CuTest *testCase) {
    // REF: -GA-TT--ACA-
    // S1 : ----TT------
    // S2 : -GA-T-CCACAA
    // S3 : ---GTT--ACA-

    int64_t refLength = 7;
    char *refString = "GATTACA";
    char *refName = "ref";

    int64_t seqNo = 3;
    stList *seqs = stList_construct3(0, free);
    stList_append(seqs, stString_copy("TT"));
    stList_append(seqs, stString_copy("GATCCACAA"));
    stList_append(seqs, stString_copy("GTTACA"));

    stList *seqNames = stList_construct3(0, free);
    stList_append(seqNames, stString_copy("S1"));
    stList_append(seqNames, stString_copy("S2"));
    stList_append(seqNames, stString_copy("S3"));

    stList *refToSeqAlignments = stList_construct3(0, (void (*)(void *))stList_destruct);
    stList_append(refToSeqAlignments, makeAlignmentList((const int64_t[]){ 2, 0, 3, 1 }, 2));
    stList_append(refToSeqAlignments, makeAlignmentList((const int64_t[]){ 0, 0, 1, 1, 2, 2, 4, 5, 5, 6, 6, 7 }, 6));
    stList_append(refToSeqAlignments, makeAlignmentList((const int64_t[]){ 2, 1, 3, 2, 4, 3, 5, 4, 6, 5 }, 5));

    MsaView *view = msaView_construct(refString, refName, refToSeqAlignments, seqs, seqNames);

    // Do a bunch of checks

    // Check that the alignment info is right
    int64_t alignmentMatrix[] = { -1,  -1,  0,  1, -1, -1, -1,
                                  0,   1,  2, -1,  5,  6,  7,
                                  -1,  -1,  1,  2,  3,  4,  5 };

    for(int64_t i=0; i<refLength; i++) {
        for(int64_t j=0; j<seqNo; j++) {
            CuAssertIntEquals(testCase, alignmentMatrix[j*refLength + i], msaView_getSeqCoordinate(view, i, j));
        }
    }

    int64_t precedingInsertLengthMatrix[] = { 0,  0,  0,  0,  0, 0, 0, 0,
                                              0,  0,  0,  0,  2, 0, 0, 1,
                                              0,  0,  1,  0,  0, 0, 0, 0 };

    int64_t precedingInsertStartMatrix[] = { -1,  -1,  -1,  -1,  -1, -1, -1, -1,
                                             -1,  -1,  -1,  -1,   3, -1, -1,  8,
                                             -1,  -1,   0,  -1,  -1, -1, -1, -1 };

    for(int64_t j=0; j<seqNo; j++) {
        for(int64_t i=0; i<refLength+1; i++) {
            CuAssertIntEquals(testCase, precedingInsertLengthMatrix[(j*(refLength+1)) + i], msaView_getPrecedingInsertLength(view, i, j));
            CuAssertIntEquals(testCase, precedingInsertStartMatrix[(j*(refLength+1)) + i], msaView_getPrecedingInsertStart(view, i, j));
        }
    }

    int64_t maxIndelLengths[] = { 0, 0, 1, 0, 2, 0, 0, 1 };

    for(int64_t i=0; i<refLength+1; i++) {
        CuAssertIntEquals(testCase, maxIndelLengths[i], msaView_getMaxPrecedingInsertLength(view, i));
    }

    // Print it
    if (st_getLogLevel() >= debug) {
        msaView_print(view, 1, stderr);
    }

    // Cleanup
    stList_destruct(seqs);
    stList_destruct(seqNames);
    stList_destruct(refToSeqAlignments);
    msaView_destruct(view);
}

static struct List *readSequences(char *fastaFile) {
    struct List *seqs = constructEmptyList(0, free);
    struct List *seqLengths = constructEmptyList(0, free);
    struct List *headers = constructEmptyList(0, free);

    FILE *fH = fopen(fastaFile, "r");
    fastaRead(fH, seqs, seqLengths, headers);
    fclose(fH);

    destructList(seqLengths);
    destructList(headers);

    return seqs;
}

char *getString(char *string, bool rle) {
    if(rle) {
        RleString *r = rleString_construct(string);
        string = stString_copy(r->rleString);
        rleString_destruct(r);
    }
    else {
        string = stString_copy(string);
    }
    return string;
}

void test_viewExamples(CuTest *testCase) {
	char *path=TEST_POLISH_FILES_DIR"largeExamples";
	int64_t exampleNo = 1;
	bool rle = 1;

	for(int64_t example=0; example<exampleNo; example++) {
		char *readFile = stString_print("%s/%i.fasta", path, (int)example);
		char *trueRefFile = stString_print("%s/%i.ref.fasta", path, (int)example);

		st_logInfo("Doing view test with %s read files and %s true ref file\n", readFile, trueRefFile);

        // Parse reads & reference
        struct List *r = readSequences((char *)readFile);
        assert(r->length > 1);
        RleString *rleReference = rleString_construct(r->list[0]);
        char *reference = getString(r->list[0], rle);
        stList *nucleotides = stList_construct3(0, free);
        stList *bamChunkReads = stList_construct3(0, (void (*)(void*))bamChunkRead_destruct);
        stList *rleReads = stList_construct3(0, (void (*)(void *))rleString_destruct);
        // TODO: Get examples with strands specified
        for(int64_t i=1; i<r->length; i++) {
            bool forwardStrand = TRUE;
            char *nucl = getString(r->list[i], rle);
            BamChunkRead *bcr = bamChunkRead_construct2(stString_print("read_%d", i), stString_copy(nucl),
                    NULL, forwardStrand, NULL);
            stList_append(nucleotides, nucl);
            stList_append(bamChunkReads, bcr);
            stList_append(rleReads, rleString_construct(r->list[i]));
        }
        destructList(r);

		// Parse reference
		struct List *trueReferenceList = readSequences((char *)trueRefFile);
		assert(trueReferenceList->length == 1);
		char *trueReference = getString(trueReferenceList->list[0], rle);
		RleString *rleTrueReference = rleString_construct(trueReferenceList->list[0]);
		destructList(trueReferenceList);

		// Polish params
		FILE *fh = fopen(polishParamsFile, "r");
		Params *params = params_readParams(fh);
		fclose(fh);

		// Set parameters
		params->polishParams->maxPoaConsensusIterations = 100;
		params->polishParams->minPoaConsensusIterations = 0;
		params->polishParams->maxRealignmentPolishIterations = 3;
		params->polishParams->minRealignmentPolishIterations = 3;

		// Generate alignment
		Poa *poa = poa_realignAll(bamChunkReads, NULL, reference, params->polishParams);

		// Generate final MEA read alignments to POA
		stList *alignments = poa_getReadAlignmentsToConsensus(poa, bamChunkReads, params->polishParams);

		// Make seq names
		stList *seqNames = stList_construct3(0, free);
		for(int64_t i=0; i<stList_length(bamChunkReads); i++) {
			stList_append(seqNames, stString_print("SEQ:%i", (int)i));
		}
		stList_append(seqNames, stString_print("TRUE_REF"));

		// Get an alignment between the inferred reference and the true reference and add it

		double alignmentScore;
		stList *refToTrueRefAlignment = getShiftedMEAAlignment(poa->refString, trueReference, NULL,
															   params->polishParams->p, params->polishParams->sM, 0, 0, &alignmentScore);

		stList_append(alignments, refToTrueRefAlignment);
		stList_append(nucleotides, trueReference);
		stList_append(rleReads, rleTrueReference);

		// Print alignment
		//TODO msaView_construct takes in nucleotides, not BCRs
		MsaView *view = msaView_construct(poa->refString, NULL,
								   	      alignments, nucleotides, seqNames);

		if (st_getLogLevel() >= debug) {
			msaView_print(view, 2, stderr);

			if(rle) {
				// Expand the RLE string
				RleString *rleConsensusString = expandRLEConsensus(poa, rleReads, bamChunkReads, params->polishParams->repeatSubMatrix);
				CuAssertIntEquals(testCase, rleConsensusString->length, stList_length(poa->nodes)-1);

				//msaView_printRepeatCounts(view, 1,
				//		rleConsensusString, rleReads, stderr);

				rleString_destruct(rleConsensusString);
			}
		}

		int64_t indelLength = 0;
		for(int64_t i=0; i<view->refLength; i++) {
			indelLength += msaView_getMaxPrecedingInsertLength(view, i);
		}

		st_logInfo("Got %i indels\n", (int)indelLength);

		// Simple stats
        int64_t totalMatches = calcSequenceMatches(poa->refString, trueReference);
		st_logInfo("Got %f sequence identity between predicted and true reference.\n", 2.0*totalMatches/(strlen(poa->refString)+strlen(trueReference)));

		// Cleanup
		rleString_destruct(rleReference);
		stList_destruct(rleReads);
		free(readFile);
		free(trueRefFile);
		stList_destruct(bamChunkReads);
		free(reference);
		stList_destruct(alignments);
		poa_destruct(poa);
		stList_destruct(seqNames);
		msaView_destruct(view);
		params_destruct(params);
	}
}

CuSuite* viewTestSuite(void) {
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, test_view);
    SUITE_ADD_TEST(suite, test_viewExamples);

    return suite;
}
