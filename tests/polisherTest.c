/*
 * Copyright (C) 2009-2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "margin.h"

static char *polishParamsFile = "../params/allParams.np.json";
#define TEST_POLISH_FILES_DIR "../tests/polishTestExamples/"

static void test_poa_getReferenceGraph(CuTest *testCase) {
	/*
	 * Test building a trivial poa graph containing just a reference string.
	 */

	char *reference = "GATTACA";

	Alphabet *alphabet = alphabet_constructNucleotide();
	Poa *poa = poa_getReferenceGraph(reference, alphabet);

	CuAssertTrue(testCase, stList_length(poa->nodes) == strlen(reference) + 1);
	for(int64_t i=0; i<strlen(reference); i++) {
		PoaNode *node = stList_get(poa->nodes, i+1);

		CuAssertTrue(testCase, node->base == reference[i]);
		CuAssertTrue(testCase, stList_length(node->inserts) == 0);
		CuAssertTrue(testCase, stList_length(node->deletes) == 0);
	}

	PoaNode *node = stList_get(poa->nodes, 0);
	CuAssertTrue(testCase, node->base == 'N');
	CuAssertTrue(testCase, stList_length(node->inserts) == 0);
	CuAssertTrue(testCase, stList_length(node->deletes) == 0);

	poa_destruct(poa);
}

static char *makeShiftedString(char *str, char *insert, int64_t insertPoint) {
	char *suffix = stString_copy(&str[insertPoint]);
	char *prefix = stString_copy(str);
	prefix[insertPoint] = '\0';
	char *shiftedStr = stString_print("%s%s%s", prefix, insert, suffix);
	free(suffix);
	free(prefix);
	return shiftedStr;
}

static void test_getShift(CuTest *testCase) {
	/*
	 * Test left shifting code.
	 */

	for(int64_t test=0; test<10000; test++) {
		// Make random string
		int64_t length = st_randomInt(1, 20);
		char *str = getRandomACGTSequence(length);

		// Make random insert of length m
		int64_t m = st_randomInt(1, 4);
		char *insert = getRandomACGTSequence(m);

		// Run get shift
		int64_t i = getShift(str, length, insert, m);

		//if(i + 2 < length) {
		//	fprintf(stderr, "Str: %s, str-length:%" PRIi64 " insert: %s, insert:%" PRIi64 "\n", str, length, insert, i);
		//}

		// Test resulting transplanted string is same as concatenated str+insert
		char *shiftedStr = makeShiftedString(str, insert, i);
		char *concatenatedStr = stString_print("%s%s", str, insert);

		CuAssertStrEquals(testCase, concatenatedStr, shiftedStr);

		// Cleanup
		free(shiftedStr);

		// Test no further left shift would work
		for(int64_t j=0; j<i; j++) {
			shiftedStr = makeShiftedString(str, insert, j);
			CuAssertTrue(testCase, !stString_eq(shiftedStr, concatenatedStr));
			free(shiftedStr);
		}

		// Cleanup
		free(concatenatedStr);
		free(str);
		free(insert);
	}
}

static void checkInserts(CuTest *testCase, Poa *poa, int64_t nodeIndex,
					     int64_t insertNumber, const char **inserts, const double *insertWeights, bool divideWeights) {
	PoaNode *node = stList_get(poa->nodes, nodeIndex);

	CuAssertIntEquals(testCase, stList_length(node->inserts), insertNumber);

	for(int64_t i=0; i<insertNumber; i++) {
		PoaInsert *poaInsert = stList_get(node->inserts, i);
		CuAssertStrEquals(testCase, inserts[i], poaInsert->insert);
		CuAssertDblEquals(testCase, insertWeights[i], poaInsert_getWeight(poaInsert) / (divideWeights ? PAIR_ALIGNMENT_PROB_1 : 1.0), 0.001);
	}
}

static void checkDeletes(CuTest *testCase, Poa *poa, int64_t nodeIndex,
					     int64_t deleteNumber, const int64_t *deleteLengths, const double *deleteWeights, bool divideWeights) {
	PoaNode *node = stList_get(poa->nodes, nodeIndex);

	CuAssertIntEquals(testCase, stList_length(node->deletes), deleteNumber);

	for(int64_t i=0; i<deleteNumber; i++) {
		PoaDelete *poaDelete = stList_get(node->deletes, i);
		CuAssertIntEquals(testCase, deleteLengths[i], poaDelete->length);
		CuAssertDblEquals(testCase, deleteWeights[i], poaDelete_getWeight(poaDelete) / (divideWeights ? PAIR_ALIGNMENT_PROB_1 : 1.0), 0.001);
	}
}

static void checkNode(CuTest *testCase, Poa *poa, int64_t nodeIndex, char base, const double *baseWeights,
		int64_t insertNumber, const char **inserts, const double *insertWeights,
		int64_t deleteNumber, const int64_t *deleteLengths, const double *deleteWeights) {

	PoaNode *node = stList_get(poa->nodes, nodeIndex);
	CuAssertTrue(testCase, node->base == base);

	// Matches
	for(int64_t i=0; i<poa->alphabet->alphabetSize; i++) {
		CuAssertDblEquals(testCase, node->baseWeights[i], baseWeights[i], 0.0);
	}

	// Inserts
	checkInserts(testCase, poa, nodeIndex, insertNumber, inserts, insertWeights, 0);

	// Deletes
	checkDeletes(testCase, poa, nodeIndex, deleteNumber, deleteLengths, deleteWeights, 0);
}

static void test_poa_augment_example(CuTest *testCase) {
	/*
	 * Test poa_augment gives works as expected on a small example.
	 */

	char *reference = "GATTACA";

	Alphabet *alphabet = alphabet_constructNucleotide();
	Poa *poa = poa_getReferenceGraph(reference, alphabet);

	char *read = "GATACGGT";

	stList *matches = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
	stList *inserts = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
	stList *deletes = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);

	stList_append(matches, stIntTuple_construct3(100, 0, 0));
	stList_append(matches, stIntTuple_construct3(100, 1, 1));
	stList_append(matches, stIntTuple_construct3(50, 2, 2));
	stList_append(matches, stIntTuple_construct3(50, 3, 2));
	stList_append(matches, stIntTuple_construct3(100, 4, 3));
	stList_append(matches, stIntTuple_construct3(100, 5, 4));
	stList_append(matches, stIntTuple_construct3(50, 6, 5));
	stList_append(matches, stIntTuple_construct3(25, 6, 6));
	stList_append(matches, stIntTuple_construct3(25, 6, 7));

	stList_append(inserts, stIntTuple_construct3(50, 5, 5));
	stList_append(inserts, stIntTuple_construct3(25, 5, 6));
	stList_append(inserts, stIntTuple_construct3(50, 6, 6));
	stList_append(inserts, stIntTuple_construct3(75, 6, 7));

	stList_append(deletes, stIntTuple_construct3(50, 2, 1));
	stList_append(deletes, stIntTuple_construct3(50, 3, 2));

	poa_augment(poa, read, 1, 0, matches, inserts, deletes);

	// Check POA graph is what we expect

	CuAssertTrue(testCase, stList_length(poa->nodes) == 8); // Length + prefix node

	checkNode(testCase, poa, 0, 'N', (const double[]){ 0.0, 0.0, 0.0, 0.0, 0.0 },
			0, (const char *[]){ "" }, (const double[]){ 0.0 },
			0, (const int64_t[]){ 0 }, (const double[]){ 0.0 });

	checkNode(testCase, poa, 1, 'G', (const double[]){ 0.0, 0.0, 100.0, 0.0, 0.0 },
				0, (const char *[]){ "" }, (const double[]){ 0.0 },
				0, (const int64_t[]){ 0 }, (const double[]){ 0.0 });

	checkNode(testCase, poa, 2, 'A', (const double[]){ 100.0, 0.0, 0.0, 0.0, 0.0 },
					0, (const char *[]){ "" }, (const double[]){ 0.0 },
					1, (const int64_t[]){ 1 }, (const double[]){ 100.0 });

	checkNode(testCase, poa, 3, 'T', (const double[]){ 0.0, 0.0, 0.0, 50.0, 0.0 },
					0, (const char *[]){ "" }, (const double[]){ 0.0 },
					0, (const int64_t[]){ 0 }, (const double[]){ 0.0 });

	checkNode(testCase, poa, 4, 'T', (const double[]){ 0.0, 0.0, 0.0, 50.0, 0.0 },
					0, (const char *[]){ "" }, (const double[]){ 0.0 },
					0, (const int64_t[]){ 0 }, (const double[]){ 0.0 });

	checkNode(testCase, poa, 5, 'A', (const double[]){ 100.0, 0.0, 0.0, 0.0, 0.0 },
					0, (const char *[]){ "" }, (const double[]){ 0.0 },
					0, (const int64_t[]){ 0 }, (const double[]){ 0.0 });

	checkNode(testCase, poa, 6, 'C', (const double[]){ 0.0, 100.0, 0.0, 0.0, 0.0 },
					2, (const char *[]){ "G", "GG" }, (const double[]){ 50.0, 25.0 },
					0, (const int64_t[]){ 0 }, (const double[]){ 0.0 });

	checkNode(testCase, poa, 7, 'A', (const double[]){ 0.0, 0.0, 75.0, 25.0, 0.0 },
					2, (const char *[]){ "GT", "T" }, (const double[]){ 50.0, 75.0 },
					0, (const int64_t[]){ 0 }, (const double[]){ 0.0 });

	// Cleanup
	poa_destruct(poa);
	stList_destruct(matches);
	stList_destruct(inserts);
	stList_destruct(deletes);
}

static void test_poa_realign_tiny_example1(CuTest *testCase) {
	/*
	 * Tests that poa_realign builds the expected poa graph for a small example of input sequences
	 */

	char *reference = "GATACAGCGGG";
	BamChunkRead *read = bamChunkRead_construct2(stString_print("read"), stString_print("GATTACAGCG"), NULL, 1, 0);

	stList *reads = stList_construct3(0,(void (*)(void *))bamChunkRead_destruct);
	stList_append(reads, read);
	bool readStrand = 1;

	FILE *fh = fopen(polishParamsFile, "r");
	Params *params = params_readParams(fh);
	fclose(fh);
	PolishParams *polishParams = params->polishParams;
	
	// This test used the default state machine in cPecan
	stateMachine_destruct(polishParams->sM);
	polishParams->sM = stateMachine3_constructNucleotide(threeState);

	/*
	// Generate set of posterior probabilities for matches, deletes and inserts with respect to reference.
	stList *matches = NULL, *inserts = NULL, *deletes = NULL;
	getAlignedPairsWithIndels(sM, reference, read, p, &matches, &deletes, &inserts, 0, 0);

	for(int64_t i=0; i<stList_length(matches); i++) {
		stIntTuple *alignedPair = stList_get(matches, i);
		fprintf(stderr, "Match: (x:%i) (y:%i) (weight:%f)\n", stIntTuple_get(alignedPair, 1),
				stIntTuple_get(alignedPair, 2), ((float)stIntTuple_get(alignedPair, 0))/PAIR_ALIGNMENT_PROB_1);
	}

	for(int64_t i=0; i<stList_length(inserts); i++) {
		stIntTuple *alignedPair = stList_get(inserts, i);
		fprintf(stderr, "Insert: (x:%i) (y:%i) (weight:%f)\n", stIntTuple_get(alignedPair, 1),
				stIntTuple_get(alignedPair, 2), ((float)stIntTuple_get(alignedPair, 0))/PAIR_ALIGNMENT_PROB_1);
	}

	for(int64_t i=0; i<stList_length(deletes); i++) {
		stIntTuple *alignedPair = stList_get(deletes, i);
		fprintf(stderr, "Delete: (x:%i) (y:%i) (weight:%f)\n", stIntTuple_get(alignedPair, 1),
					stIntTuple_get(alignedPair, 2), ((float)stIntTuple_get(alignedPair, 0))/PAIR_ALIGNMENT_PROB_1);
	}*/

	Poa *poa = poa_realign(reads, NULL, reference, polishParams);
	// Check we get the set of inserts and deletes we expect

	// . = match
	// | = insert
	// - = delete
	// : = insert and match
	// % = delete and match

	//     Reference
	//     -1 0 1 2 3 4 5 6 7 8 9 10
	//      N G A T A C A G C G G G
	// -1 N .
	//  0 G   .
	//  1 A   | .
	//  2 T     : . -
	//  3 T       : . %
	//  4 A         :   .
	//  5 C           .   .
	//  6 A             . - %
	//  7 G               . - %
	//  8 C                 . - %
	//  9 G                   . - %

	st_logInfo("Read:%s\n", read);
	st_logInfo("Reference:%s\n", reference);
	if (st_getLogLevel() >= info) {
		poa_print(poa, stderr, reads, 0.0, 0.0);
	}

	// Check inserts

	// A after ref 0
	checkInserts(testCase, poa, 1, 1, (const char *[]){ "A" }, (const double[]){ 0.038656 }, 1);
	// T after ref 1
	checkInserts(testCase, poa, 2, 1, (const char *[]){ "T" }, (const double[]){ 0.91074 }, 1);
	// T after ref 2
	checkInserts(testCase, poa, 3, 1, (const char *[]){ "A" }, (const double[]){ 0.038831 }, 1);
	// A after ref 3
	checkInserts(testCase, poa, 4, 0, (const char *[]){ "" }, (const double[]){ 1 }, 1);

	checkInserts(testCase, poa, 0, 0, (const char *[]){ "" }, (const double[]){ 1 }, 1);
	checkInserts(testCase, poa, 5, 0, (const char *[]){ "" }, (const double[]){ 1 }, 1);
	checkInserts(testCase, poa, 6, 0, (const char *[]){ "" }, (const double[]){ 1 }, 1);
	checkInserts(testCase, poa, 7, 0, (const char *[]){ "" }, (const double[]){ 1 }, 1);
	checkInserts(testCase, poa, 8, 0, (const char *[]){ "" }, (const double[]){ 1 }, 1);
	checkInserts(testCase, poa, 9, 0, (const char *[]){ "" }, (const double[]){ 1 }, 1);
	checkInserts(testCase, poa, 10, 0, (const char *[]){ "" }, (const double[]){ 1 }, 1);

	// Check deletes

	/*
	Delete: (x:3) (y:2) (weight:0.021429)
	Delete: (x:4) (y:3) (weight:0.011958)

	Delete: (x:6) (y:6) (weight:0.039542)
	Delete: (x:7) (y:6) (weight:0.041150)

	Delete: (x:7) (y:7) (weight:0.039140)
	Delete: (x:8) (y:7) (weight:0.038841)

	Delete: (x:8) (y:8) (weight:0.440979)
	Delete: (x:9) (y:8) (weight:0.438735)

	Delete: (x:9) (y:9) (weight:0.437247)
	Delete: (x:10) (y:9) (weight:0.440302)*/

	// No deletes first three positions
	checkDeletes(testCase, poa, 0, 0, (const int64_t[]){ 1 }, (const double[]){ 1 }, 1);
	checkDeletes(testCase, poa, 1, 0, (const int64_t[]){ 1 }, (const double[]){ 1 }, 1);
	checkDeletes(testCase, poa, 2, 0, (const int64_t[]){ 1 }, (const double[]){ 1 }, 1);
	checkDeletes(testCase, poa, 5, 0, (const int64_t[]){ 1 }, (const double[]){ 1 }, 1);
	checkDeletes(testCase, poa, 7, 0, (const int64_t[]){ 1 }, (const double[]){ 1 }, 1);
	checkDeletes(testCase, poa, 9, 0, (const int64_t[]){ 1 }, (const double[]){ 1 }, 1);
	checkDeletes(testCase, poa, 10, 0, (const int64_t[]){ 1 }, (const double[]){ 1 }, 1);

	// L1 after ref 2
	checkDeletes(testCase, poa, 3, 0, (const int64_t[]){ 1 }, (const double[]){ 1 }, 1);
	// L1 after ref 3
	checkDeletes(testCase, poa, 4, 0, (const int64_t[]){ 1 }, (const double[]){ 1 }, 1);
	// L2 after ref 5
	checkDeletes(testCase, poa, 6, 1, (const int64_t[]){ 2 }, (const double[]){ 0.077178 }, 1);
	// L2 after ref 7
	checkDeletes(testCase, poa, 8, 2, (const int64_t[]){ 1, 2 }, (const double[]){ 0.935157, 0.885846 }, 1);

	params_destruct(params);
	poa_destruct(poa);
	stList_destruct(reads);
}

static void test_poa_realign(CuTest *testCase) {
	/*
	 * Test poa_realign by generating lots of random examples
	 */

	for (int64_t test = 0; test < 100; test++) {

		//Make true reference
		char *trueReference = getRandomSequence(st_randomInt(1, 100));

		// Make starting reference
		char *reference = evolveSequence(trueReference);

		FILE *fh = fopen(polishParamsFile, "r");
		Params *params = params_readParams(fh);
		fclose(fh);
		PolishParams *polishParams = params->polishParams;

		// Reads
		int64_t readNumber = st_randomInt(0, 20);
		stList *reads = stList_construct3(0, (void(*)(void*)) bamChunkRead_destruct);
		for(int64_t i=0; i<readNumber; i++) {
			stList_append(reads, bamChunkRead_construct2(stString_print("read_%d", i), evolveSequence(trueReference),NULL,TRUE,
					params->polishParams->useRunLengthEncoding));
		}

		Poa *poa = poa_realign(reads, NULL, reference, polishParams);

		// Generate the read alignments and check the matches
		// Currently don't check the insert and deletes

		double *baseWeights = st_calloc(poa->alphabet->alphabetSize*strlen(reference), sizeof(double));

		for(int64_t i=0; i<readNumber; i++) {
			char *read = ((BamChunkRead*)stList_get(reads, i))->rleRead->rleString;

			// Generate set of posterior probabilities for matches, deletes and inserts with respect to reference.
			stList *matches = NULL, *inserts = NULL, *deletes = NULL;
			getAlignedPairsWithIndels(polishParams->sMConditional, reference, read, polishParams->p, &matches, &deletes, &inserts, 0, 0);

			// Collate matches
			for(int64_t j=0; j<stList_length(matches); j++) {
				stIntTuple *match = stList_get(matches, j);
				baseWeights[stIntTuple_get(match, 1) * poa->alphabet->alphabetSize + poa->alphabet->convertCharToSymbol(read[stIntTuple_get(match, 2)])] += stIntTuple_get(match, 0);
			}

			// Cleanup
			stList_destruct(matches);
			stList_destruct(inserts);
			stList_destruct(deletes);
		}

		// Check match weights tally
		for(int64_t i=0; i<strlen(reference); i++) {
			PoaNode *poaNode = stList_get(poa->nodes, i+1);
			for(int64_t j=0; j<poa->alphabet->alphabetSize; j++) {
				CuAssertDblEquals(testCase, poaNode->baseWeights[j], baseWeights[i*poa->alphabet->alphabetSize + j], 0.0001);
			}
		}

		st_logInfo("True-reference:%s\n", trueReference);
		if (st_getLogLevel() >= info) {
			poa_print(poa, stderr, reads, 5, 5);
		}

		//Cleanup
		free(baseWeights);
		free(trueReference);
		free(reference);
		stList_destruct(reads);
		poa_destruct(poa);
		params_destruct(params);
	}
}

static void test_poa_realignIterative(CuTest *testCase) {
	/*
	 * Test random small examples against poa_realignIterative
	 */

	for (int64_t test = 0; test < 100; test++) {
		FILE *fh = fopen(polishParamsFile, "r");
		Params *params = params_readParams(fh);
		fclose(fh);
		PolishParams *polishParams = params->polishParams;

		//Make true reference
		char *trueReference = getRandomSequence(st_randomInt(1, 100));

		// Make starting reference
		char *reference = evolveSequence(trueReference);

		// Reads
		int64_t readNumber = st_randomInt(0, 20);
		stList *reads = stList_construct3(0, (void(*)(void*)) bamChunkRead_destruct);
		for(int64_t i=0; i<readNumber; i++) {
			stList_append(reads, bamChunkRead_construct2(stString_print("Read_%d", i), evolveSequence(trueReference),
					NULL, st_random() > 0.5, polishParams->useRunLengthEncoding));
		}

		Poa *poa = poa_realignIterative(reads, NULL, reference, polishParams);

		st_logInfo("True-reference:%s\n", trueReference);
		if (st_getLogLevel() >= info) {
			poa_print(poa, stderr, reads, 5, 0);
		}

		//Cleanup
		free(trueReference);
		free(reference);
		stList_destruct(reads);
		poa_destruct(poa);
		params_destruct(params);
	}
}

int64_t calcSequenceMatches(char *seq1, char *seq2) {
	FILE *fh = fopen(polishParamsFile, "r");
	Params *params = params_readParams(fh);
	fclose(fh);
	PolishParams *polishParams = params->polishParams;

	//Get identity
	stList *allAlignedPairs = getAlignedPairs(polishParams->sM, seq1, seq2, polishParams->p, 0, 0);
	stList *alignedPairs = filterPairwiseAlignmentToMakePairsOrdered(allAlignedPairs, seq1, seq2, params->polishParams->p);

	int64_t matches = getNumberOfMatchingAlignedPairs(seq1, seq2, alignedPairs);

	// Cleanup
	params_destruct(params);
	//stList_destruct(alignedPairs);

	return matches;
}

typedef struct _alignmentMetrics {
	int64_t totalConsensusMatches;
	int64_t totalReferenceMatches;
	int64_t totalConsensusLength;
	int64_t totalReferenceLength;
	int64_t totalTrueReferenceLength;
} AlignmentMetrics;

static void test_poa_realign_example(CuTest *testCase, RleString *rleTrueReference, RleString *rleReference,
		stList *reads, AlignmentMetrics *rleAlignmentMetrics, AlignmentMetrics *nonRleAlignmentMetrics) {
	FILE *fh = fopen(polishParamsFile, "r");
	Params *params = params_readParams(fh);
	fclose(fh);
	PolishParams *polishParams = params->polishParams;
	
	// Set parameters
	params->polishParams->maxPoaConsensusIterations = 100;
	params->polishParams->minPoaConsensusIterations = 0;
	params->polishParams->referenceBasePenalty = 0.6;
	params->polishParams->maxRealignmentPolishIterations = 3;
	params->polishParams->minRealignmentPolishIterations = 3;

	Poa *poa = poa_realign(reads, NULL, rleReference->rleString, polishParams);
	Poa *poaRefined = poa_realignAll(reads, NULL, rleReference->rleString, polishParams);

	Poa *poaTrue = poa_realign(reads, NULL, rleTrueReference->rleString, polishParams);

	// Run phasing
	//stList *anchorAlignments = poa_getAnchorAlignments(poaRefined, NULL, stList_length(reads), params->polishParams);
	//stList *reads1, *reads2;
	//phaseReads(poaRefined->refString, strlen(poaRefined->refString), reads, anchorAlignments, &reads1, &reads2, params);
	//Poa *poaReads1 = poa_realignIterative(reads1, NULL, poaRefined->refString, polishParams);
	//Poa *poaReads2 = poa_realignIterative(reads2, NULL, poaRefined->refString, polishParams);

	// Look at non-rle comparison
	RleString *consensusRleString = expandRLEConsensus(poaRefined, reads, polishParams->repeatSubMatrix);
	char *nonRLEConsensusString = rleString_expand(consensusRleString);
	rleString_destruct(consensusRleString);

	// Calculate alignments between true reference and consensus and starting reference sequences
	int64_t consensusMatches = calcSequenceMatches(rleTrueReference->rleString, poaRefined->refString);
	int64_t referenceMatches = calcSequenceMatches(rleTrueReference->rleString, rleReference->rleString);
	char *trueReference = rleString_expand(rleTrueReference);
	char *reference = rleString_expand(rleReference);
	int64_t nonRLEConsensusMatches = calcSequenceMatches(trueReference, nonRLEConsensusString);
	int64_t nonRLEReferenceMatches = calcSequenceMatches(trueReference, reference);
	//int64_t consensusMatchesReads1 = calcSequenceMatches(rleTrueReference->rleString, poaReads1->refString);
	//int64_t consensusMatchesReads2 = calcSequenceMatches(rleTrueReference->rleString, poaReads2->refString);

	// Update the running total alignment metrics
	if(rleAlignmentMetrics != NULL) {
		rleAlignmentMetrics->totalConsensusMatches += consensusMatches;
		rleAlignmentMetrics->totalReferenceMatches += referenceMatches;
		rleAlignmentMetrics->totalConsensusLength += strlen(poaRefined->refString);
		rleAlignmentMetrics->totalReferenceLength += rleReference->length;
		rleAlignmentMetrics->totalTrueReferenceLength += rleTrueReference->length;
	}

	// Update the running total alignment metrics
	if(nonRleAlignmentMetrics != NULL) {
		nonRleAlignmentMetrics->totalConsensusMatches += nonRLEConsensusMatches;
		nonRleAlignmentMetrics->totalReferenceMatches += nonRLEReferenceMatches;
		nonRleAlignmentMetrics->totalConsensusLength += strlen(nonRLEConsensusString);
		nonRleAlignmentMetrics->totalReferenceLength += strlen(reference);
		nonRleAlignmentMetrics->totalTrueReferenceLength += strlen(trueReference);
	}

	// Log some stuff
	if (st_getLogLevel() >= info) {
		st_logInfo("Reference:      \t\t%s\n", rleReference->rleString);
		st_logInfo("True-reference: \t\t%s\n", rleTrueReference->rleString);
		st_logInfo("Consensus:      \t\t%s\n", poaRefined->refString);
		//st_logInfo("Consensus Reads1:\t%s\n", poaReads1->refString);
		//st_logInfo("Consensus Reads2:\t%s\n", poaReads2->refString);
		st_logInfo("Reference stats:     \t");
		poa_printSummaryStats(poa, stderr);
		st_logInfo("Consensus stats:     \t");
		poa_printSummaryStats(poaRefined, stderr);
		//st_logInfo("Reads 1 stats\t");
		//poa_printSummaryStats(poaReads1, stderr);
		//st_logInfo("Reads 2 stats\t");
		//poa_printSummaryStats(poaReads2, stderr);
		st_logInfo("True-reference stats:\t");
		poa_printSummaryStats(poaTrue, stderr);
		st_logInfo("Consensus : true-ref identity: %f\n", 2.0*consensusMatches/(rleTrueReference->length + strlen(poaRefined->refString)));
		//st_logInfo("Reads 1 consensus : true-ref identity: %f\n", 2.0*consensusMatchesReads1/(rleTrueReference->length + strlen(poaReads1->refString)));
		//st_logInfo("Reads 2 consensus : true-ref identity: %f\n", 2.0*consensusMatchesReads2/(rleTrueReference->length + strlen(poaReads2->refString)));
		st_logInfo("Start-ref : true-ref identity: %f\n", 2.0*referenceMatches/(rleTrueReference->length + rleReference->length));
		//st_logInfo("Total reads: %i, # reads partition1: %i, # reads partition2: %i\n", (int)stList_length(reads), (int)stList_length(reads1), (int)stList_length(reads2));
		// Non-RLE stats
		st_logInfo("Non-RLE Reference:     \t\t%s\n", reference);
		st_logInfo("Non-RLE True-reference:\t\t%s\n", trueReference);
		st_logInfo("Non-RLE Consensus:     \t\t%s\n", nonRLEConsensusString);
		st_logInfo("Non-RLE Consensus : true-ref identity: %f\n", 2.0*nonRLEConsensusMatches/(strlen(trueReference) + strlen(nonRLEConsensusString)));
		st_logInfo("Non-RLE Start-ref : true-ref identity: %f\n\n", 2.0*nonRLEReferenceMatches/(strlen(trueReference) + strlen(reference)));
	}

	if (st_getLogLevel() >= debug && !stString_eq(rleTrueReference->rleString, poaRefined->refString)) {
		//poa_print(poa, stderr, 5);
		//poa_print(poaRefined, stderr, reads, 2, 0);

		//poa_printTSV(poa, stderr, reads, 2, 0);

		///poa_printRepeatCounts(poa, stderr, reads);
	}

	// Cleanup
	params_destruct(params);
	poa_destruct(poa);
	poa_destruct(poaRefined);
	poa_destruct(poaTrue);
	free(trueReference);
	free(reference);
	//rleString_destruct(rleTrueReference);
	//rleString_destruct(rleReference);
	free(nonRLEConsensusString);
	//poa_destruct(poaReads1);
	//poa_destruct(poaReads2);
	//stList_destruct(anchorAlignments);
	//stList_destruct(reads1);
	//stList_destruct(reads2);
}

static struct List *readSequences(char *fastaFile, struct List **headers) {
	struct List *seqs = constructEmptyList(0, free);
	struct List *seqLengths = constructEmptyList(0, free);
	*headers = constructEmptyList(0, free);

	FILE *fH = fopen(fastaFile, "r");
	fastaRead(fH, seqs, seqLengths, *headers);
	fclose(fH);

	destructList(seqLengths);

	return seqs;
}

static void test_poa_realign_examples(CuTest *testCase, const char **examples, int64_t exampleNo, bool rle) {
	AlignmentMetrics *alignmentMetrics = st_calloc(1, sizeof(AlignmentMetrics));
	AlignmentMetrics *rleAlignmentMetrics = st_calloc(1, sizeof(AlignmentMetrics));
	for(int64_t example=0; example<exampleNo; example++) {
		const char *readFile = examples[example*2];
		const char *trueRefFile = examples[example*2+1];

		st_logInfo("\nDoing polish test with %s read files and %s true ref file\n", readFile, trueRefFile);

		// Parse sequences
		struct List *readHeaders;
		struct List *nucleotides = readSequences((char *)readFile, &readHeaders);
		assert(nucleotides->length > 1);
		struct List *trueReferenceHeaders;
		struct List *trueReferenceList = readSequences((char *)trueRefFile, &trueReferenceHeaders);
		assert(trueReferenceList->length == 1);

		// Parse strands (and create reads)
		stList *reads = stList_construct3(0, (void (*)(void*))bamChunkRead_destruct);
		for(int64_t i=1; i<readHeaders->length; i++) {
			char *header = readHeaders->list[i];
			char strand = header[strlen(header)-1];
			CuAssertTrue(testCase, strand == 'F' || strand == 'R');
			stList_append(reads, bamChunkRead_construct2(stString_print("read_%d", i),
					stString_copy(nucleotides->list[i]), NULL, strand == 'F', rle));
		}

		//if(strlen(reads->list[0]) > strlen(trueReferenceList->list[0]) * 0.8 || reads->length < 30) {
		//	fprintf(stderr, "Got short input ref:\n%s\n%s\n", reads->list[0], trueReferenceList->list[0]);
		//	continue;
		//}

		RleString *trueReference = rle ? rleString_construct(trueReferenceList->list[0]) : rleString_construct_no_rle(trueReferenceList->list[0]);
		RleString *reference = rle ? rleString_construct(nucleotides->list[0]) : rleString_construct_no_rle(nucleotides->list[0]);

		// Run poa iterative realign
		test_poa_realign_example(testCase, trueReference, reference,
								reads, rleAlignmentMetrics, alignmentMetrics);

		// Cleanup
		destructList(nucleotides);
		destructList(readHeaders);
		destructList(trueReferenceList);
		destructList(trueReferenceHeaders);
		free(reads);
		rleString_destruct(trueReference);
		rleString_destruct(reference);
	}

	// Alignment metrics for set
	st_logInfo("Total consensus identity: %f for %i bases, vs. %f starting identity\n",
			2.0*alignmentMetrics->totalConsensusMatches/(alignmentMetrics->totalConsensusLength+alignmentMetrics->totalTrueReferenceLength),
			alignmentMetrics->totalConsensusLength,
			2.0*alignmentMetrics->totalReferenceMatches/(alignmentMetrics->totalReferenceLength+alignmentMetrics->totalTrueReferenceLength));
	if(rle) {
		st_logInfo("RLE Space: Total consensus identity: %f for %i bases, vs. %f starting identity\n",
					2.0*rleAlignmentMetrics->totalConsensusMatches/(rleAlignmentMetrics->totalConsensusLength+rleAlignmentMetrics->totalTrueReferenceLength),
					rleAlignmentMetrics->totalConsensusLength,
					2.0*rleAlignmentMetrics->totalReferenceMatches/(rleAlignmentMetrics->totalReferenceLength+rleAlignmentMetrics->totalTrueReferenceLength));
	}

	free(alignmentMetrics);
	free(rleAlignmentMetrics);
}

static void test_poa_realign_examples_large(CuTest *testCase, int64_t exampleNo, const char *path, bool rle) {
	// Build strings
	const char **examples = st_calloc(2*exampleNo, sizeof(char *));
	for(int64_t i=0; i<exampleNo; i++) {
		examples[2*i] = stString_print("%s/%i.fasta", path, (int)i);
		examples[2*i+1] = stString_print("%s/%i.ref.fasta", path, (int)i);
	}

	test_poa_realign_examples(testCase, examples, exampleNo, rle);

	// Cleanup
	for(int64_t i=0; i<exampleNo*2; i++) {
		free((char *)examples[i]);
	}
	free(examples);
}

void test_poa_realign_ecoli_examples_rle(CuTest *testCase) {
	test_poa_realign_examples_large(testCase, 20, TEST_POLISH_FILES_DIR"20_random_100bp_windows_directional_ecoli_guppy", 1);
}

void test_poa_realign_ecoli_examples_no_rle(CuTest *testCase) {
	test_poa_realign_examples_large(testCase, 20, TEST_POLISH_FILES_DIR"20_random_100bp_windows_directional_ecoli_guppy", 0);
}

void test_poa_realign_ecoli_many_examples_rle(CuTest *testCase) {
	test_poa_realign_examples_large(testCase, 100, TEST_POLISH_FILES_DIR"500_random_100bp_windows_directional_ecoli_guppy", 1);
}

void test_poa_realign_ecoli_many_examples_no_rle(CuTest *testCase) {
	test_poa_realign_examples_large(testCase, 100, TEST_POLISH_FILES_DIR"500_random_100bp_windows_directional_ecoli_guppy", 0);
}

static void test_rleString_example(CuTest *testCase, const char *testStr,
		int64_t rleLength, int64_t nonRleLength,
		const char *testStrRLE, const int64_t *repeatCounts,
		const int64_t *nonRleToRleCoordinateMap) {
	RleString *rleString = rleString_construct((char *)testStr);
	uint64_t *nonRleToRleCoordinateMap2 = rleString_getNonRleToRleCoordinateMap(rleString);

	CuAssertIntEquals(testCase, rleLength, rleString->length);
	CuAssertStrEquals(testCase, testStrRLE, rleString->rleString);
	for(int64_t i=0; i<rleLength; i++) {
		CuAssertIntEquals(testCase, repeatCounts[i], rleString->repeatCounts[i]);
	}

	CuAssertIntEquals(testCase, nonRleLength, rleString->nonRleLength);
	for(int64_t i=0; i<nonRleLength; i++) {
		CuAssertIntEquals(testCase, nonRleToRleCoordinateMap[i], nonRleToRleCoordinateMap2[i]);
	}

	char *expandedRleString = rleString_expand(rleString);
	CuAssertStrEquals(testCase, testStr, expandedRleString);

	free(expandedRleString);
	rleString_destruct(rleString);
	free(nonRleToRleCoordinateMap2);
}

static void test_rleString_examples(CuTest *testCase) {
	test_rleString_example(testCase, "GATTACAGGGGTT", 8, 13, "GATACAGT", (const int64_t[]){ 1,1,2,1,1,1,4,2 },
			(const int64_t[]){ 0,1,2,2,3,4,5,6,6,6,6,7,7 });

	test_rleString_example(testCase, "TTTTT", 1, 5, "T", (const int64_t[]){ 5 },
			(const int64_t[]){ 0,0,0,0,0 });

	test_rleString_example(testCase, "", 0, 0, "", (const int64_t[]){ 1 },
			(const int64_t[]){ 0 });

	test_rleString_example(testCase, "TTTTTCC", 2, 7, "TC", (const int64_t[]){ 5, 2 },
			(const int64_t[]){ 0,0,0,0,0,1,1 });
}

void checkStringsAndFree(CuTest *testCase, const char *expected, char *temp) {
	CuAssertStrEquals(testCase, expected, temp);
	free(temp);
}

void test_addInsert(CuTest *testCase) {
	checkStringsAndFree(testCase, "GATTACA", addInsert("GAACA", "TT", 2));
	checkStringsAndFree(testCase, "GATTACA", addInsert("", "GATTACA", 0));
	checkStringsAndFree(testCase, "GATTACA", addInsert("ATTACA", "G", 0));
	checkStringsAndFree(testCase, "GATTACA", addInsert("GATTAC", "A", 6));
	checkStringsAndFree(testCase, "GATTACA", addInsert("GATTACA", "", 6));
	checkStringsAndFree(testCase, "GATTACA", addInsert("GATTACA", "", 3));
}

void test_removeDelete(CuTest *testCase) {
	checkStringsAndFree(testCase, "GATTACA", removeDelete("GATTGGACA", 2, 4));
	checkStringsAndFree(testCase, "GATTACA", removeDelete("GATTACA", 0, 0));
	checkStringsAndFree(testCase, "GATTACA", removeDelete("GATTACATT", 2, 7));
	checkStringsAndFree(testCase, "GATTACA", removeDelete("AGATTACA", 1, 0));
}

void test_polishParams(CuTest *testCase) {
	FILE *fh = fopen(polishParamsFile, "r");
	Params *params = params_readParams(fh);
	fclose(fh);
	PolishParams *polishParams = params->polishParams;

	CuAssertTrue(testCase, polishParams->useRunLengthEncoding);
	CuAssertDblEquals(testCase, polishParams->referenceBasePenalty, 0.5, 0);
	CuAssertDblEquals(testCase, polishParams->minPosteriorProbForAlignmentAnchorsLength, 6, 0);

	CuAssertDblEquals(testCase, polishParams->minPosteriorProbForAlignmentAnchors[0], 0.9, 0);
	CuAssertDblEquals(testCase, polishParams->minPosteriorProbForAlignmentAnchors[1], 10, 0);
	CuAssertDblEquals(testCase, polishParams->minPosteriorProbForAlignmentAnchors[2], 0.95, 0);
	CuAssertDblEquals(testCase, polishParams->minPosteriorProbForAlignmentAnchors[3], 4, 0);
	CuAssertDblEquals(testCase, polishParams->minPosteriorProbForAlignmentAnchors[4], 0.99, 0);
	CuAssertDblEquals(testCase, polishParams->minPosteriorProbForAlignmentAnchors[5], 0, 0);

	CuAssertDblEquals(testCase, polishParams->p->threshold, 0.01, 0);
	CuAssertDblEquals(testCase, polishParams->p->minDiagsBetweenTraceBack, 10000, 0);
	CuAssertDblEquals(testCase, polishParams->p->traceBackDiagonals, 40, 0);
	CuAssertDblEquals(testCase, polishParams->p->diagonalExpansion, 10, 0);
	CuAssertDblEquals(testCase, polishParams->p->constraintDiagonalTrim, 0, 0);
	CuAssertDblEquals(testCase, polishParams->p->splitMatrixBiggerThanThis, 250000000, 0);
	CuAssertDblEquals(testCase, polishParams->p->gapGamma, 0.5, 0);
	CuAssertTrue(testCase, !polishParams->p->alignAmbiguityCharacters);

	/*CuAssertDblEquals(testCase,  repeatSubMatrix_getLogProb(polishParams->repeatSubMatrix, a, 0, 0, 0), -0.059686935, 0);
	CuAssertDblEquals(testCase,  repeatSubMatrix_getLogProb(polishParams->repeatSubMatrix, c, 0, 0, 0), -0.055418707, 0);
	CuAssertDblEquals(testCase,  repeatSubMatrix_getLogProb(polishParams->repeatSubMatrix, g, 0, 0, 0), -0.05438334, 0);
	CuAssertDblEquals(testCase,  repeatSubMatrix_getLogProb(polishParams->repeatSubMatrix, t, 0, 0, 0), -0.035762809, 0);
	CuAssertDblEquals(testCase,  repeatSubMatrix_getLogProb(polishParams->repeatSubMatrix, a, 1, 0, 0), -0.036856437, 0);
	CuAssertDblEquals(testCase,  repeatSubMatrix_getLogProb(polishParams->repeatSubMatrix, c, 1, 0, 0), -0.062816805, 0);
	CuAssertDblEquals(testCase,  repeatSubMatrix_getLogProb(polishParams->repeatSubMatrix, g, 1, 0, 0), -0.055853556, 0);
	CuAssertDblEquals(testCase,  repeatSubMatrix_getLogProb(polishParams->repeatSubMatrix, t, 1, 0, 0), -0.065273937, 0);*/

	params_destruct(params);
}

void test_removeOverlapExample(CuTest *testCase) {
	FILE *fh = fopen(polishParamsFile, "r");
	Params *params = params_readParams(fh);
	fclose(fh);
	PolishParams *polishParams = params->polishParams;

	//Make prefix
	char *prefixString = stString_copy("ACGTGATTTCA");

	// Make sufix
	char *suffixString = stString_copy("GATTTCAACGT");

	int64_t approxOverlap = 10;

	// Run overlap remover
	int64_t prefixStringCropEnd, suffixStringCropStart;
	double overlapWeight = removeOverlap(prefixString, suffixString, approxOverlap, polishParams,
				  	  	  &prefixStringCropEnd, &suffixStringCropStart);

	CuAssertIntEquals(testCase, 7, prefixStringCropEnd);
	CuAssertIntEquals(testCase, 3, suffixStringCropStart);

	// Cleanup
	params_destruct(params);
	free(prefixString);
	free(suffixString);
}

void test_removeOverlap_RandomExamples(CuTest *testCase) {
	FILE *fh = fopen(polishParamsFile, "r");
	Params *params = params_readParams(fh);
	fclose(fh);
	PolishParams *polishParams = params->polishParams;

	for (int64_t test = 0; test < 100; test++) {
		//Make prefix
		char *prefixString = getRandomSequence(st_randomInt(1, 100));

		// Make sufix
		char *suffixString = getRandomSequence(st_randomInt(1, 100));

		int64_t approxOverlap = st_randomInt(0, 100);

		// Run overlap remover
		int64_t prefixStringCropEnd, suffixStringCropStart;
		removeOverlap(prefixString, suffixString, approxOverlap, polishParams,
					  &prefixStringCropEnd, &suffixStringCropStart);

		CuAssertTrue(testCase, prefixStringCropEnd >= 0);
		CuAssertTrue(testCase, prefixStringCropEnd <= strlen(prefixString));

		CuAssertTrue(testCase, suffixStringCropStart >= 0);
		CuAssertTrue(testCase, suffixStringCropStart <= strlen(suffixString));

		free(prefixString);
		free(suffixString);
	}

	params_destruct(params);
}

int64_t polishingTest(char *bamFile, char *referenceFile, char *paramsFile, char *region, bool verbose, bool diploid) {

    // Run margin phase
    char *logString = verbose ? "--logLevel DEBUG" : "--logLevel INFO";
    char *regionStr = region == NULL ? stString_print("") : stString_print("--region %s", region);
    char *diploidString = diploid ? "--diploid" : "";
    char *command = stString_print("./margin %s %s %s %s %s %s", bamFile, referenceFile, paramsFile, regionStr, logString, diploidString);
    st_logInfo("> Running command: %s\n", command);

    int64_t i = st_system(command);
    free(regionStr);
    free(command);
    return i;
}

void test_polish5kb_rle(CuTest *testCase) {
    char *referenceFile = "../tests/hg19.chr3.9mb.fa";
    bool verbose = false;
    char *bamFile = "../tests/NA12878.np.chr3.5kb.bam";
    char *region = "chr3:2150000-2155000";

    st_logInfo("\n\nTesting polishing on %s\n", bamFile);
    int64_t i = polishingTest(bamFile, referenceFile, polishParamsFile, region, verbose, 0);
    CuAssertTrue(testCase, i == 0);

    st_logInfo("\n\nTesting diploid polishing on %s\n", bamFile);
    i = polishingTest(bamFile, referenceFile, polishParamsFile, region, verbose, 1);
    CuAssertTrue(testCase, i == 0);
}

void test_polish5kb_no_region(CuTest *testCase) {
    char *referenceFile = "../tests/hg19.chr3.9mb.fa";
    bool verbose = false;
    char *bamFile = "../tests/NA12878.np.chr3.5kb.bam";

    st_logInfo("\n\nTesting polishing on %s\n", bamFile);
    int64_t i = polishingTest(bamFile, referenceFile, polishParamsFile, NULL, verbose, 0);
    CuAssertTrue(testCase, i == 0);

    st_logInfo("\n\nTesting diploid polishing on %s\n", bamFile);
    i = polishingTest(bamFile, referenceFile, polishParamsFile, NULL, verbose, 1);
    CuAssertTrue(testCase, i == 0);
}

void test_polish100kb(CuTest *testCase) {
	char *referenceFile = "../tests/hg19.chr3.9mb.fa";
	bool verbose = false;
	char *bamFile = "../tests/NA12878.np.chr3.100kb.4.bam";
	char *region = "chr3:8100000-8200000";

	st_logInfo("\n\nTesting polishing on %s\n", bamFile);
	int64_t i = polishingTest(bamFile, referenceFile, polishParamsFile, region, verbose, 0);
	CuAssertTrue(testCase, i == 0);

	st_logInfo("\n\nTesting polishing on %s\n", bamFile);
	i = polishingTest(bamFile, referenceFile, polishParamsFile, region, verbose, 1);
	CuAssertTrue(testCase, i == 0);
}

void test_binomialPValue(CuTest *testCase) {
	CuAssertDblEquals(testCase, 252.0, bionomialCoefficient(10, 5), 0.001);
	CuAssertDblEquals(testCase, 15504.0, bionomialCoefficient(20, 15), 0.001);
	CuAssertDblEquals(testCase, 80347448443237920.0, bionomialCoefficient(64, 22), 0.001);
	CuAssertDblEquals(testCase, 151473214816.0, bionomialCoefficient(64, 10), 0.001);
	CuAssertDblEquals(testCase, 1832624140942590534.0, bionomialCoefficient(64, 32), 0.001);
}

CuSuite* polisherTestSuite(void) {
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, test_poa_getReferenceGraph);
    SUITE_ADD_TEST(suite, test_poa_augment_example);
    SUITE_ADD_TEST(suite, test_poa_realign_tiny_example1);
    SUITE_ADD_TEST(suite, test_poa_realign);
    SUITE_ADD_TEST(suite, test_poa_realignIterative);
    SUITE_ADD_TEST(suite, test_getShift);
    SUITE_ADD_TEST(suite, test_rleString_examples);
    SUITE_ADD_TEST(suite, test_addInsert);
    SUITE_ADD_TEST(suite, test_removeDelete);
    SUITE_ADD_TEST(suite, test_polishParams);
    SUITE_ADD_TEST(suite, test_removeOverlapExample);
    SUITE_ADD_TEST(suite, test_removeOverlap_RandomExamples);

    SUITE_ADD_TEST(suite, test_polish5kb_rle);
    SUITE_ADD_TEST(suite, test_polish5kb_no_region);
    SUITE_ADD_TEST(suite, test_polish100kb);
    SUITE_ADD_TEST(suite, test_binomialPValue);

    SUITE_ADD_TEST(suite, test_poa_realign_ecoli_examples_rle);
    SUITE_ADD_TEST(suite, test_poa_realign_ecoli_examples_no_rle);
    SUITE_ADD_TEST(suite, test_poa_realign_ecoli_many_examples_rle);
    SUITE_ADD_TEST(suite, test_poa_realign_ecoli_many_examples_no_rle);

    return suite;
}
