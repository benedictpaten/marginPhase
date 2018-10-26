/*
 * Copyright (C) 2009-2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <multipleAligner.h>
#include "CuTest.h"
#include "stPolish.h"
#include "randomSequences.h"

static char *nanoporeHmmFile = "../params/polish/threeStateNanopore.hmm";
static char *repeatCountsModelFile = "../params/polish/log_prob_matrices_fasta_one_liners_2x_pseudocounts.txt";
#define TEST_POLISH_FILES_DIR "../tests/polishTestExamples/"

static void test_poa_getReferenceGraph(CuTest *testCase) {
	/*
	 * Test building a trivial poa graph containing just a reference string.
	 */

	char *reference = "GATTACA";

	Poa *poa = poa_getReferenceGraph(reference);

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
		CuAssertDblEquals(testCase, insertWeights[i], poaInsert->weight / (divideWeights ? PAIR_ALIGNMENT_PROB_1 : 1.0), 0.001);
	}
}

static void checkDeletes(CuTest *testCase, Poa *poa, int64_t nodeIndex,
					     int64_t deleteNumber, const int64_t *deleteLengths, const double *deleteWeights, bool divideWeights) {
	PoaNode *node = stList_get(poa->nodes, nodeIndex);

	CuAssertIntEquals(testCase, stList_length(node->deletes), deleteNumber);

	for(int64_t i=0; i<deleteNumber; i++) {
		PoaDelete *poaDelete = stList_get(node->deletes, i);
		CuAssertIntEquals(testCase, deleteLengths[i], poaDelete->length);
		CuAssertDblEquals(testCase, deleteWeights[i], poaDelete->weight / (divideWeights ? PAIR_ALIGNMENT_PROB_1 : 1.0), 0.001);
	}
}

static void checkNode(CuTest *testCase, Poa *poa, int64_t nodeIndex, char base, const double *baseWeights,
		int64_t insertNumber, const char **inserts, const double *insertWeights,
		int64_t deleteNumber, const int64_t *deleteLengths, const double *deleteWeights) {

	PoaNode *node = stList_get(poa->nodes, nodeIndex);
	CuAssertTrue(testCase, node->base == base);

	// Matches
	for(int64_t i=0; i<SYMBOL_NUMBER; i++) {
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

	Poa *poa = poa_getReferenceGraph(reference);

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

	poa_augment(poa, read, 0, matches, inserts, deletes);

	if (st_getLogLevel() >= info) {
		poa_print(poa, stderr, 0.0);
	}

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
	char *read = "GATTACAGCG";

	stList *reads = stList_construct();
	stList_append(reads, read);

	PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
	StateMachine *sM = stateMachine3_construct(threeState);

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

	Poa *poa = poa_realign(reads, NULL, reference, sM, p);

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
		poa_print(poa, stderr, 0.0);
	}

	// Check inserts

	// A after ref 0
	checkInserts(testCase, poa, 1, 1, (const char *[]){ "A" }, (const double[]){ 0.038656 }, 1);
	// T after ref 1
	checkInserts(testCase, poa, 2, 1, (const char *[]){ "T" }, (const double[]){ 0.874535 }, 1);
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
	checkDeletes(testCase, poa, 3, 1, (const int64_t[]){ 1 }, (const double[]){ 0.021429 }, 1);
	// L1 after ref 3
	checkDeletes(testCase, poa, 4, 1, (const int64_t[]){ 1 }, (const double[]){ 0.011958 }, 1);
	// L2 after ref 5
	checkDeletes(testCase, poa, 6, 1, (const int64_t[]){ 2 }, (const double[]){ 0.078383 }, 1);
	// L2 after ref 7
	checkDeletes(testCase, poa, 8, 1, (const int64_t[]){ 2 }, (const double[]){ 0.87598 }, 1);

	stateMachine_destruct(sM);
	pairwiseAlignmentBandingParameters_destruct(p);
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

		// Reads
		int64_t readNumber = st_randomInt(0, 20);
		stList *reads = stList_construct3(0, free);
		for(int64_t i=0; i<readNumber; i++) {
			stList_append(reads, evolveSequence(trueReference));
		}

		PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();

		Hmm *hmm = hmm_loadFromFile(nanoporeHmmFile);

		StateMachine *sM = hmm_getStateMachine(hmm); //stateMachine3_construct(threeState);

		Poa *poa = poa_realign(reads, NULL, reference, sM, p);

		poa_normalize(poa); // Shift all the indels

		// Generate the read alignments and check the matches
		// Currently don't check the insert and deletes

		double *baseWeights = st_calloc(SYMBOL_NUMBER*strlen(reference), sizeof(double));

		for(int64_t i=0; i<readNumber; i++) {
			char *read = stList_get(reads, i);

			// Generate set of posterior probabilities for matches, deletes and inserts with respect to reference.
			stList *matches = NULL, *inserts = NULL, *deletes = NULL;
			getAlignedPairsWithIndels(sM, reference, read, p, &matches, &deletes, &inserts, 0, 0);

			// Collate matches
			for(int64_t j=0; j<stList_length(matches); j++) {
				stIntTuple *match = stList_get(matches, j);
				baseWeights[stIntTuple_get(match, 1) * SYMBOL_NUMBER + symbol_convertCharToSymbol(read[stIntTuple_get(match, 2)])] += stIntTuple_get(match, 0);
			}

			// Cleanup
			stList_destruct(matches);
			stList_destruct(inserts);
			stList_destruct(deletes);
		}

		// Check match weights tally
		for(int64_t i=0; i<strlen(reference); i++) {
			PoaNode *poaNode = stList_get(poa->nodes, i+1);
			for(int64_t j=0; j<SYMBOL_NUMBER; j++) {
				CuAssertDblEquals(testCase, poaNode->baseWeights[j], baseWeights[i*SYMBOL_NUMBER + j], 0.0001);
			}
		}

		st_logInfo("True-reference:%s\n", trueReference);
		if (st_getLogLevel() >= info) {
			poa_print(poa, stderr, 5);
		}

		//Cleanup
		free(baseWeights);
		stateMachine_destruct(sM);
		free(trueReference);
		free(reference);
		stList_destruct(reads);
		pairwiseAlignmentBandingParameters_destruct(p);
		poa_destruct(poa);
		hmm_destruct(hmm);
	}
}

static void test_poa_realignIterative(CuTest *testCase) {
	/*
	 * Test random small examples against poa_realignIterative
	 */

	for (int64_t test = 0; test < 100; test++) {

		//Make true reference
		char *trueReference = getRandomSequence(st_randomInt(1, 100));

		// Make starting reference
		char *reference = evolveSequence(trueReference);

		// Reads
		int64_t readNumber = st_randomInt(0, 20);
		stList *reads = stList_construct3(0, free);
		for(int64_t i=0; i<readNumber; i++) {
			stList_append(reads, evolveSequence(trueReference));
		}

		PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();

		Hmm *hmm = hmm_loadFromFile(nanoporeHmmFile);

		StateMachine *sM = hmm_getStateMachine(hmm); //stateMachine3_construct(threeState);

		Poa *poa = poa_realignIterative(reads, NULL, reference, sM, p);

		st_logInfo("True-reference:%s\n", trueReference);
		if (st_getLogLevel() >= info) {
			poa_print(poa, stderr, 5);
		}

		//Cleanup
		stateMachine_destruct(sM);
		free(trueReference);
		free(reference);
		stList_destruct(reads);
		pairwiseAlignmentBandingParameters_destruct(p);
		poa_destruct(poa);
		hmm_destruct(hmm);
	}
}

static const char *readArrayExample1[] = {
		"CATTTTTCTCCTCCACCTGCAACAGAAGATAAAAACGCGCATCACAAACTACTTTATTG",
		"CATTTTTCTCTCCGTCACGTAATAGGAAAACAGATGAAAATGTGCACCATAAAACGCATTTTTATTT",
		"CATTTTCTCTCTCCGTCACGACAGGAAACAGATGAAAATGGGCACAAGACCACAAACGCATTTTGAT",
		"CATTTTTCTCCGGTCATTTAATGAAAACAGATGGTACTGCGTATGTGACATAAACGCATTTTTATTT",
		"CATTTCCTCCGTCACTGCACAGGAAAACAGATGAAAATGCAAGTATGGACCCACAAAACGCATTTTATTT",
		"CATTTTTTCTCTCTCCGTCAGCTGCATTGAAAATGATGAAATGCGGGTATGACTATAAACGCATTTATTT",
		"CATTTTTTTTCTCTCCTCCACACACAGGAAACAGATGAAAAATGTATGTGACCATAAAACGCATTTTATTT",
		"TATTTTCTCCGTCATTGCAGGAAAACAGATGAAATGTAAAGTATGTGAATTACAAACGGTTTTTTTTATTT",
		"CATTTTTCTCCTCCGTCATTGCACAGGAGTCAGATGAAAATGCGCATGTGACCATAACGCATTTTTTTATTT",
		"CATTTTTCTCCTCCGTCATACCGTGAAACAGATGAAAAATGCGGGCATGGGACCATAAAACGCATTTTTATTT",
		"CATTTTTCTCCTCCGTCATTGCACAGGAAAACAGATGAAAACGTGGGGCATGTGACCATAAACGCATTTTTATT",
		"CATTTTCTCTCCTCGTGTTGCACAGGAAAACAGATGAAAAATGCGAGATATGTGATCCACAAACATTTTTATTT",
		"CATTTTTCTCCTCCGTCATTGCACAGGAAAATGATGAAAATGCGGGGCATGTGACCATAAAACGCATTTTTATTT",
		"CATTTTCTCTCTCCCTCGTCATTGCACAGGAAAACAGATGAAAATGCAGGGCATGTGACCATAAAACGCATTTTTT",
		"CATTTTCTCTCCTCCACATTGCACAGGAAAACAGATGAAAATGCGGCATGTGACCATAAAACGCATTTCTTTATTT",
		"CATTTTCTCCGTCAGTCAACAATATGAAAACAGATGAAACGCGGGCACGTGACCATAAAACGCATTTTTTTTATTT",
		"CATTTTTCTCCTCCGTCATTGCATTGTGGAACAGATGAAAATGCGGGGTATGTGAATCATAAAACGCATTTTATTT",
		"CATTTTTCTCTCCGTCATTGCATTAGAAAACAGGGATGAAAATGCGGGCATGTGACCATAAAAACGCATTTTTATTT",
		"CATTTTTCTCTCTCCTCCGTCATTGCACAGGAAAACAGATGAAAAATGCGCGTGACTATAAAACGCATTTTTATTTT",
		"CATTTTCCTCTCCCTCCGTCATTTGCACAGGAAAACAGATGAAAAAATGCGGAATGGCTATTATAAACATTTTTAACT",
		"CATTTTTTTCTCCTCTGTCATTGCACAGGAAAACAGATGAAAAATGCGTATGTGACCATAAAATCCATTTCTTTTATTT",
		"CATTTTTCTCCTCCGTCATTGCACAGGAAAATGATGAAAAAATGCGGGCATGTGACCATAAAACGTGCATTTTTTATTT",
		"CATTTTCTCTCTCCTCCGTGTTGCACAGGAAAACCAGATGAAAATGCGGAACATGTGTTCATAAAACGCATTTTTATTT",
		"CATTTTCTCTCCCTCCGTCATTGCACAGGAAAACAGATGAAAATGCAGGGCAATAATGACCATAAAACGCATTTTTATTT",
		"CATTTTCTCTCCTCTCGTCATTTGCACAGGAAGAGCAGATGAAAATGCAGGGCATGTGACCATAAAACGCATTTTTATTT",
		"CCATTTTCTCTCTCCCTCCGTCATTGCACAGGAAAAGCAGATGAAAAATGCGGGCATGTGACCATAAAACGCATTTTATTT",
		"CATTTTTTTCTCTCCTGTCATTGCACAGGAAACAAAGAGATGAAAAATGCGGGCATGTGACCATAAAACGCATTTTTATTT",
		"CATTTTTCTCTCCCTCCGTCATTGCATAGGAAAACAGATGAAAATGCGGGGTATGTGGACCATAAAACGCATTTTTTATTT",
		"CATTTTCTCTCTCCCTCCGTCATTGCACAGGGAAAACAGATGAAAATTGCGGGGCATGTGACCATAAAACGCATTTTTATTT",
		"CATTTTCTCTCTCCCTCCGTCATTTGCACAGGAAAACAGATGAAAAATGCGGGGCATGTGACCATAAAACGCATTTTTTATTT",
		"CATTTTTCTCTCCCTCCGTCACTGCACAGGAAAAACAGATGAAAATGCGGGGCATGCATCATAAAACGTATTTTTATTGAATTT",
		"CATTTTCTCTCTCCCTCCGTCATTGCACAGGAAAACAGATAAGAAAAATGCAGGGGCATGTGACCATAAAACGCATTTTTATTT",
		"CATTTTTTCACTACTCTCCCTCCGTCGTACTGGAAAACAAACAGATAAATGCAGGGCATGTGACCATAAAACATTTTTTTATTT",
		"CATTTTTCTCTCTCCCTCCGTCATTGCACAGGAAAACAGATAAAAAAAAATGCAGGGGCATGTGACCATAAAACATTTTTATTT",
		"CATTTTTTCTCTCTCTCGTGTTGCACACAGGAAAACAGATGAAAAATGCCGGGGCATCATGACCATAAAACGCGTTTTTTTATTT",
		"CATTTTCTCTCTCCCTCCGTCATTGCACAGGAAAACAGATGAAAAATGCAGGGGCGTAACTGACCATAAAACGCATTTTTTATTT",
		"CATTTTTTCTCTCCTCCGTCATTGCACAGGAAAAATGTGATGAAAATGCGGGGTATGTGACCATAAAACGCATTTTTATGCTTCT",
		"CATTTTCTCTCTCCCTCCGTCATTGCACAGGAAAACAGATGAAAAATGCGAGGACATGTGACCATAAAACGCATTTTTTTTATTT",
		"CATTTTCTCTCTCCCTCCGTCCATTGCACAGGAAAACAGATATAAAAAATGCAGGGCATCAAACCATAAAACATTTTTTTTATTT",
		"CATTTTTCTCTCCCTCCGTCATTGCAATAGGAAAACAGATATTTTGGTGTACCGCAAGTATGTGACCATAAAACGTATTTTTATTT",
		"CATTTTCTCTCTCCCTCCGTCATTGCACAGGAAAACCAGATAGAAAAAATACAGGGCATGTGTTCATAAAGCACGCATTTTTATTT",
		"CATTTTTTTCTCTCTCCTCCGCTTTCACACACAGGAGTAAACAGATGAAAAATGTGGGCATGTGACCATAAAACGCATTTTTTATTT",
		"CATTTTCTCTCTCCCTCAAAATCATTTGCACAGGAAAACAGATAGAAAAATGCAACGGGGCATGTGATATAAAACGCATTTTTTATTT",
		"CATTTTCTACTCTCTCCCTCCGTCATTGCAGGAAAACAGATGAAAATGCAGGGAACATATATGACCATAAAACGCATTTTTTTTTATTT",
		"CATTTTCTCTCTCCCTCCGTCATTGCACAGGAAAACAGATGAAAAAAGAGCTGGCATGCGGGGCATGTGACCATAAAACGCATTTTTTTGT" };
static char *referenceExample1 =     "CATTTTCTCTCCCTCCGTCATTGCACAGGAAAACAGATGAAAATGCAGGGCAATAATGACCATAAAACGCATTTTTATTT";
static char *trueReferenceExample1 = "CATTTTTCTCTCTCCCTCCGTCATTGCACAGGAAAACAGATGAAAAATGCGGGGCATGTGACCATAAAACGCATTTTTTATTT";

static const char *readArrayExample2[] = {
			"GATGTAAAAATGACTGAGTTAGAACAGGCATAAATACATCTGT",
			"GATGTAAAAAAAAATGACAGAGAATAAAACTATCCTTATCTATT",
			"GATGTAAAAAGAAGCGGAAGTTAGAACAGGCATAAATACATCTGT",
			"GATGTAAAAAGAAATGACGGAAGAACAGAGCATAACACACATCTGT",
			"GATGTAAAAAAAGAATGATTTAGTTGAACAGAGCATAAATATCTGT",
			"GATGTAAAAAAAAGAAATGACGGAAGAACAGAGCATAACACATCTGT",
			"GATGTAAAAGAAATGGAGGTTAGAACAGAGCATAAATACACATCTGT",
			"GATGTAAAAAAAGAAATGATTTGGAAGAACAGAGCATAAATATCTGT",
			"GATGTAAAAGAAATGACGGAAGTTAGAATATATATAACACACATCTGT",
			"GATGTAAAAAAAGAATGGACGGTTAGAACAGAGCACAACACACATCTGT",
			"GATGTAAAAAAGAATGATAAAGTTAGAATAGAGCATAAATAACATCTGT",
			"GATGTAAAAGAAATGTGGAAGTTAGAACAGAGCATAAATACACATCTAT",
			"GATGTAAAAAAAAAGAAATGAAGCTAGAACAGAGCATAAATACATCTGT",
			"GATGTAAAAAAAAAATGACCCGGAAGTTGAACAGAGCATAATACATCTGT",
			"GATGTAAAAAAGAAATGATTTAAAGGAACAGAGCATAAATACACATCTGT",
			"GATGTAAAAAAAAGTGACGGAAGTTAAGACAGGCATAAATACACATCTGT",
			"GATGTAAAAAAGAAATGATTTGCTAGAACAGAGCATAAATACACATCTGT",
			"GATGTAAAAAAAGAAATGATGGGTTAGAATAGAGCATAAATACACATCTGT",
			"GATGTAAAAAAAAGAAATGACGAGTTAGAACCAGAGCACCATCTACATCAT"
			"GATGTAAAAAAAAAAATGACGGAAGTTAGACAAGCATAAATACACATCTAT",
			"GATGTAAAAAAAAGAATGATTTGAAGTTAGAACAGAGCATAACACATCTGT",
			"GATGTAAAAGAAAATCGACTGAAGTTAGAACAGAGCATAAATACACATCTGT",
			"GATGTAAAAAAAGAATGACGGAAGTTAGAACAGAGCATAAATACACATCTGT",
			"GATGTAAAAAAAGAAATGACGGAGTTAGAACAGAGCATAAATACACATCTAT",
			"GATGTAAAAAAAAAGAAATGTGTGAGTTAAGACAGAGCATAAATACATCTAT",
			"GATGTAAAAAAGAAATGACGGAAGTTAAGACAGAGCATAAATACACATCTATT",
			"GATGTAAAAAAAAAAATGTGGAAGTTAAAACAGAGCATAAATACACATCTATT",
			"GATGCAAAAAAAAAGAAATGACGGAAGTTAAATTAGAGCATAAATACATCTGT",
			"GATGTAAAAAAAAGAAATGATTTGGAAGTTACAGAGCATAAATACACATCTGT",
			"GATGTAAAGAAAATGATTTTAGAAGTTAGAACAGAGCATAACACAATATCTGT",
			"GATAAAAAAAAAGGAATGATTGGAAGCTAGAACAGAGCATAAATACACATCTGT",
			"GATGTAAAAAAAAGAAATGACGGAAGTTAAGACAGAGCATAAATACACATCTGT",
			"GATGTAAAAAAAAATGACGGAAGTTCTGAAACAGGCATAAATACACATCTGTAT",
			"GATGTAAAAAGAATGATTTGAAGTTAGAACAGAGTATATTAAATACACATCTGT",
			"GATGTAAAAAAAAAGAAATGACGGAAGTTAAGACAGAGCATAAATACACATCTATT",
			"GATGTAAAAAAAAAGAAATGATTTGAAGCAGAACAGAGCATAAATACAAGATCTGT",
			"GATGTAAAAAAAAGAAGAAATGACGGAAGTTAAACAGAGCATAAATACACATCTGT",
			"GATGTAAAAAAAAAAGAAATGATTGAAGTTAGAAATATACATAAATACACATCTGT",
			"GATGTAAAAAAAAAGAAATGATTTTAAAGTGAACAGAGCATAAATACACACCTTGGT",
			"GATGTAAAAAAAAAAAAGAAATGACGGAAGTTGAACTAGGCTTATAAATACATCTGT",
			"GATGCCAAAAAAAAAAAGAAATGGCCAGAGTTAGAACAGAGCATAAATACACATCTGT",
			"GATGTAAAAAAAAAGAAATGCGGATTTGGAAGTTAGAACAGTATATAAAGCACACATCCGT" };
static char *referenceExample2 =     "GATGTAAAAAAGAAATGATTTGCTAGAACAGAGCATAAATACACATCTGT";
static char *trueReferenceExample2 = "GATGTAAAAAAAAAGAAATGACGGAAGTTAGAACAGAGCATAAATACACATCTGT";

static double calcSequenceMatches(char *seq1, char *seq2) {
	PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
	Hmm *hmm = hmm_loadFromFile(nanoporeHmmFile);
	StateMachine *sM = hmm_getStateMachine(hmm); //stateMachine3_construct(threeState);

	//Get identity
	stList *allAlignedPairs = getAlignedPairs(sM, seq1, seq2, p, 0, 0);
	stList *alignedPairs = filterPairwiseAlignmentToMakePairsOrdered(allAlignedPairs, seq1, seq2, 0.0);

	double matches = getNumberOfMatchingAlignedPairs(seq1, seq2, alignedPairs);

	// Cleanup
	stateMachine_destruct(sM);
	pairwiseAlignmentBandingParameters_destruct(p);
	hmm_destruct(hmm);
	//stList_destruct(allAlignedPairs);
	stList_destruct(alignedPairs);

	return matches;
}

typedef struct _alignmentMetrics {
	int64_t totalConsensusMatches;
	int64_t totalReferenceMatches;
	int64_t totalConsensusLength;
	int64_t totalReferenceLength;
	int64_t totalTrueReferenceLength;
} AlignmentMetrics;

static void test_poa_realign_example_rle(CuTest *testCase, char *trueReference, char *reference, const char **readArray,
										 int64_t readNo, AlignmentMetrics *rleAlignmentMetrics, AlignmentMetrics *nonRleAlignmentMetrics) {
	stList *reads = stList_construct();
	stList *rleStrings = stList_construct3(0, (void (*)(void *))rleString_destruct);
	for(int64_t i=0; i<readNo; i++) {
		RleString *rleString = rleString_construct((char *)readArray[i]);
		stList_append(rleStrings, rleString);
		stList_append(reads, rleString->rleString);
	}
	RleString *rleReference = rleString_construct(reference);
	RleString *rleTrueReference = rleString_construct(trueReference);

	PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();

	Hmm *hmm = hmm_loadFromFile(nanoporeHmmFile);

	StateMachine *sM = hmm_getStateMachine(hmm); //stateMachine3_construct(threeState);

	Poa *poa = poa_realign(reads, NULL, rleReference->rleString, sM, p);
	Poa *poaRefined = poa_realignIterative(reads, NULL, rleReference->rleString, sM, p);

	//poaRefined = poa_checkMajorIndelEditsGreedily(poaRefined, reads, sM, p);

	Poa *poaTrue = poa_realign(reads, NULL, rleTrueReference->rleString, sM, p);

	RepeatSubMatrix *repeatSubMatrix = repeatSubMatrix_parse(repeatCountsModelFile);

	// Look at non-rle comparison
	char *nonRLEConsensusString = expandRLEConsensus(poaRefined, rleStrings, repeatSubMatrix);

	// Calculate alignments between true reference and consensus and starting reference sequences
	int64_t consensusMatches = calcSequenceMatches(rleTrueReference->rleString, poaRefined->refString);
	int64_t referenceMatches = calcSequenceMatches(rleTrueReference->rleString, rleReference->rleString);
	int64_t nonRLEConsensusMatches = calcSequenceMatches(trueReference, nonRLEConsensusString);
	int64_t nonRLEReferenceMatches = calcSequenceMatches(trueReference, reference);

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
		st_logInfo("Reference:\t\t%s\n", rleReference->rleString);
		st_logInfo("True-reference:\t\t%s\n", rleTrueReference->rleString);
		st_logInfo("Consensus:\t\t%s\n", poaRefined->refString);
		st_logInfo("Reference stats\t");
		poa_printSummaryStats(poa, stderr);
		st_logInfo("Consensus stats\t");
		poa_printSummaryStats(poaRefined, stderr);
		st_logInfo("True-reference stats\t");
		poa_printSummaryStats(poaTrue, stderr);
		st_logInfo("Consensus : true-ref identity: %f\n", 2.0*consensusMatches/(rleTrueReference->length + strlen(poaRefined->refString)));
		st_logInfo("Start-ref : true-ref identity: %f\n", 2.0*referenceMatches/(rleTrueReference->length + rleReference->length));
		// Non-RLE stats
		st_logInfo("Non-RLE Reference:\t\t%s\n", reference);
		st_logInfo("Non-RLE True-reference:\t\t%s\n", trueReference);
		st_logInfo("Non-RLE Consensus:\t\t%s\n", nonRLEConsensusString);
		st_logInfo("Non-RLE Consensus : true-ref identity: %f\n", 2.0*nonRLEConsensusMatches/(strlen(trueReference) + strlen(nonRLEConsensusString)));
		st_logInfo("Non-RLE Start-ref : true-ref identity: %f\n", 2.0*nonRLEReferenceMatches/(strlen(trueReference) + strlen(reference)));
	}

	if (st_getLogLevel() >= debug && !stString_eq(rleTrueReference->rleString, poaRefined->refString)) {
		//poa_print(poa, stderr, 5);
		poa_print(poaRefined, stderr, 2);
	}

	// Cleanup
	repeatSubMatrix_destruct(repeatSubMatrix);
	stateMachine_destruct(sM);
	pairwiseAlignmentBandingParameters_destruct(p);
	poa_destruct(poa);
	poa_destruct(poaRefined);
	poa_destruct(poaTrue);
	stList_destruct(reads);
	hmm_destruct(hmm);
	rleString_destruct(rleTrueReference);
	rleString_destruct(rleReference);
	stList_destruct(rleStrings);
	free(nonRLEConsensusString);
}

static void test_poa_realign_example(CuTest *testCase, char *trueReference, char *reference, const char **readArray,
									int64_t readNo, AlignmentMetrics *alignmentMetrics) {
	stList *reads = stList_construct();
	for(int64_t i=0; i<readNo; i++) {
		stList_append(reads, (char *)readArray[i]);
	}

	PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();

	Hmm *hmm = hmm_loadFromFile(nanoporeHmmFile);

	StateMachine *sM = hmm_getStateMachine(hmm); //stateMachine3_construct(threeState);

	Poa *poa = poa_realign(reads, NULL, reference, sM, p);
	Poa *poaRefined = poa_realignIterative(reads, NULL, reference, sM, p);

	poa_normalize(poa); // Shift all the indels

	// Calculate alignments between true reference and consensus and starting reference sequences
	int64_t consensusMatches = calcSequenceMatches(trueReference, poaRefined->refString);
	int64_t referenceMatches = calcSequenceMatches(trueReference, reference);

	// Update the running total alignment metrics
	if(alignmentMetrics != NULL) {
		alignmentMetrics->totalConsensusMatches += consensusMatches;
		alignmentMetrics->totalReferenceMatches += referenceMatches;
		alignmentMetrics->totalConsensusLength += strlen(poaRefined->refString);
		alignmentMetrics->totalReferenceLength += strlen(reference);
		alignmentMetrics->totalTrueReferenceLength += strlen(trueReference);
	}

	// Log some stuff
	if (st_getLogLevel() >= info) {
		st_logInfo("Reference:\t\t%s\n", reference);
		st_logInfo("True-reference:\t\t%s\n", trueReference);
		st_logInfo("Consensus:\t\t%s\n", poaRefined->refString);
		st_logInfo("Reference stats\t");
		poa_printSummaryStats(poa, stderr);
		st_logInfo("Consensus stats\t");
		poa_printSummaryStats(poaRefined, stderr);
		st_logInfo("Consensus : true-ref identity: %f\n", 2.0*consensusMatches/(strlen(trueReference) + strlen(poaRefined->refString)));
		st_logInfo("Start-ref : true-ref identity: %f\n", 2.0*referenceMatches/(strlen(trueReference) + strlen(reference)));
	}

	if (st_getLogLevel() >= debug && !stString_eq(trueReference, poaRefined->refString)) {
		//poa_print(poa, stderr, 5);
		poa_print(poaRefined, stderr, 2);
	}

	// Cleanup
	stateMachine_destruct(sM);
	pairwiseAlignmentBandingParameters_destruct(p);
	poa_destruct(poa);
	poa_destruct(poaRefined);
	stList_destruct(reads);
	hmm_destruct(hmm);
}

static void test_poa_realign_example1(CuTest *testCase) {
	//                 **                              *    *
	//               00  00000000111111111122222222223333333333444444444455
	//               01  23456789012345678901234567890123456789012345678901
	//Reference:     CA  TCTCTCTCGTCATGCACAGACAGATGATGCAGCATATGACATACGCATAT
	//True-reference:CATCTCTCTCTCGTCATGCACAGACAGATGATGC GCATGTGACATACGCATAT

	test_poa_realign_example(testCase, trueReferenceExample1, referenceExample1, readArrayExample1, 45, NULL);
}

static void test_poa_realign_rle_example1(CuTest *testCase) {
	//                 000000   0000111111111122222222223333333333444 444444455555555556666666666777777 7777
	//                 012345   6789012345678901234567890123456789012 345678901234567890123456789012345 6789
	//reference =     "CATTTT   CTCTCCCTCCGTCATTGCACAGGAAAACAGATGAAAA TGCAGGGCAATAATGACCATAAAACGCATTTTT ATTT";
	//trueReference = "CATTTTTCTCTCTCCCTCCGTCATTGCACAGGAAAACAGATGAAAAATGCGGGGCATG  TGACCATAAAACGCATTTTTTATTT";

	test_poa_realign_example_rle(testCase, trueReferenceExample1, referenceExample1, readArrayExample1, 45, NULL, NULL);
}

static void test_poa_realign_example2(CuTest *testCase) {
	//Reference:		    GATGTAAAAAA   GAAATGATTTGCTAGAACAGAGCATAAATACACATCTGT
	//True-reference:		GATGTAAAAAAAAAGAAATGACGGAAGTTAGAACAGAGCATAAATACACATCTGT
	//Predicted reference:	GATGTAAAAAAAA GAAATGATGGAAGTTAGAACAGAGCATAAATACACATCTGT

	test_poa_realign_example(testCase, trueReferenceExample2, referenceExample2, readArrayExample2, 41, NULL);
}

static void test_poa_realign_rle_example2(CuTest *testCase) {
	//                           * **
	//                000000000011111111112222222222333333333
	//                012345678901234567890123456789012345678
	//True-reference: GATGTAGATGACGAGTAGACAGAGCATATACACATCTGT
	//Reference:      GATGTAGATGATG CTAGACAGAGCATATACACATCTGT

	test_poa_realign_example_rle(testCase, trueReferenceExample2, referenceExample2, readArrayExample2, 41, NULL, NULL);
}

static int64_t exampleNo = 20;
static const char *examples[] = {
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_11396386_11396426.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_11396386_11396426_ref.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_11614960_11614996.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_11614960_11614996_ref.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_12932558_12932578.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_12932558_12932578_ref.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_13007488_13007541.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_13007488_13007541_ref.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_13107718_13107792.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_13107718_13107792_ref.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_1343860_1343898.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_1343860_1343898_ref.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_13941170_13941228.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_13941170_13941228_ref.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_2318461_2318484.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_2318461_2318484_ref.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_2730008_2730037.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_2730008_2730037_ref.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_3485131_3485192.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_3485131_3485192_ref.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_3931232_3931306.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_3931232_3931306_ref.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_4333072_4333143.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_4333072_4333143_ref.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_5398020_5398073.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_5398020_5398073_ref.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_5525849_5525911.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_5525849_5525911_ref.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_6769182_6769243.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_6769182_6769243_ref.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_7599992_7600003.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_7599992_7600003_ref.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_8020725_8020765.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_8020725_8020765_ref.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_8072702_8072725.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_8072702_8072725_ref.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_8980246_8980274.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_8980246_8980274_ref.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_9712745_9712803.fasta",
		TEST_POLISH_FILES_DIR"random_chr1_windows/NC_003279.8_9712745_9712803_ref.fasta"
};

static int64_t messyExampleNo = 8;
static const char *messyExamples[] = {
		TEST_POLISH_FILES_DIR"messy_windows/NC_003279.8_10029532_10029615.fasta",
		TEST_POLISH_FILES_DIR"messy_windows/NC_003279.8_10029532_10029615_ref.fasta",
		TEST_POLISH_FILES_DIR"messy_windows/NC_003279.8_10031827_10031861.fasta",
		TEST_POLISH_FILES_DIR"messy_windows/NC_003279.8_10031827_10031861_ref.fasta",
		TEST_POLISH_FILES_DIR"messy_windows/NC_003279.8_10037167_10037249.fasta",
		TEST_POLISH_FILES_DIR"messy_windows/NC_003279.8_10037167_10037249_ref.fasta",
		TEST_POLISH_FILES_DIR"messy_windows/NC_003279.8_10039004_10039029.fasta",
		TEST_POLISH_FILES_DIR"messy_windows/NC_003279.8_10039004_10039029_ref.fasta",
		TEST_POLISH_FILES_DIR"messy_windows/NC_003279.8_10040234_10040295.fasta",
		TEST_POLISH_FILES_DIR"messy_windows/NC_003279.8_10040234_10040295_ref.fasta",
		TEST_POLISH_FILES_DIR"messy_windows/NC_003279.8_1004298_1004407.fasta",
		TEST_POLISH_FILES_DIR"messy_windows/NC_003279.8_1004298_1004407_ref.fasta",
		TEST_POLISH_FILES_DIR"messy_windows/NC_003279.8_10044514_10044568.fasta",
		TEST_POLISH_FILES_DIR"messy_windows/NC_003279.8_10044514_10044568_ref.fasta",
		TEST_POLISH_FILES_DIR"messy_windows/NC_003279.8_14952113_14952138.fasta",
		TEST_POLISH_FILES_DIR"messy_windows/NC_003279.8_14952113_14952138_ref.fasta"
};

struct List *readSequences(char *fastaFile) {
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

static void test_poa_realign_examples(CuTest *testCase, const char **examples, int64_t exampleNo, bool rle) {
	AlignmentMetrics *alignmentMetrics = st_calloc(1, sizeof(AlignmentMetrics));
	AlignmentMetrics *rleAlignmentMetrics = st_calloc(1, sizeof(AlignmentMetrics));
	for(int64_t example=0; example<exampleNo; example++) {
		const char *readFile = examples[example*2];
		const char *trueRefFile = examples[example*2+1];

		// Parse sequences
		struct List *reads = readSequences((char *)readFile);
		assert(reads->length > 1);
		struct List *trueReferenceList = readSequences((char *)trueRefFile);
		assert(trueReferenceList->length == 1);

		// Run poa iterative realign
		if(rle) {
			test_poa_realign_example_rle(testCase, trueReferenceList->list[0], reads->list[0],
							//(char **)(&reads->list[1]), reads->length-1 > 20 ? 20 : reads->length-1, rle, alignmentMetrics);
							(const char **)(&reads->list[1]), reads->length-1, rleAlignmentMetrics, alignmentMetrics);
		}
		else {
			test_poa_realign_example(testCase, trueReferenceList->list[0], reads->list[0],
					//(char **)(&reads->list[1]), reads->length-1 > 20 ? 20 : reads->length-1, rle, alignmentMetrics);
					(const char **)(&reads->list[1]), reads->length-1, alignmentMetrics);
		}

		// Cleanup
		destructList(reads);
		destructList(trueReferenceList);
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

static void test_poa_realign_examples_no_rle(CuTest *testCase) {
	test_poa_realign_examples(testCase, examples, exampleNo, 0);
}

static void test_poa_realign_examples_rle(CuTest *testCase) {
	test_poa_realign_examples(testCase, examples, exampleNo, 1);
}

static void test_poa_realign_messy_examples_no_rle(CuTest *testCase) {
	test_poa_realign_examples(testCase, messyExamples, messyExampleNo, 0);
}

static void test_poa_realign_messy_examples_rle(CuTest *testCase) {
	test_poa_realign_examples(testCase, messyExamples, messyExampleNo, 1);
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

//void test_poa_realign_examples_very_large_rle(CuTest *testCase) {
//	test_poa_realign_examples_large(testCase, 2000, TEST_POLISH_FILES_DIR"2000_random_windows_chr1_celegans_guppy", 1);
//}

//void test_poa_realign_examples_very_large_no_rle(CuTest *testCase) {
//	test_poa_realign_examples_large(testCase, 2000, TEST_POLISH_FILES_DIR"2000_random_windows_chr1_celegans_guppy", 0);
//}

void test_poa_realign_examples_large_rle(CuTest *testCase) {
	test_poa_realign_examples_large(testCase, 200, TEST_POLISH_FILES_DIR"200_random_windows_chr1_celegans_guppy", 1);
}

void test_poa_realign_examples_large_no_rle(CuTest *testCase) {
	test_poa_realign_examples_large(testCase, 200, TEST_POLISH_FILES_DIR"200_random_windows_chr1_celegans_guppy", 0);
}

void test_poa_realign_examples_long_rle(CuTest *testCase) {
	test_poa_realign_examples_large(testCase, 20, TEST_POLISH_FILES_DIR"largeExamples", 1);
}

void test_poa_realign_examples_long_no_rle(CuTest *testCase) {
	test_poa_realign_examples_large(testCase, 20, TEST_POLISH_FILES_DIR"largeExamples", 0);
}

static void test_rleString_example(CuTest *testCase, const char *testStr, int64_t rleLength, const char *testStrRLE, const int64_t *repeatCounts) {
	RleString *rleString = rleString_construct((char *)testStr);

	CuAssertIntEquals(testCase, rleLength, rleString->length);
	CuAssertStrEquals(testCase, testStrRLE, rleString->rleString);
	for(int64_t i=0; i<rleLength; i++) {
		CuAssertIntEquals(testCase, repeatCounts[i], rleString->repeatCounts[i]);
	}

	rleString_destruct(rleString);
}

static void test_rleString_examples(CuTest *testCase) {
	test_rleString_example(testCase, "GATTACAGGGGTT", 8, "GATACAGT", (const int64_t[]){ 1,1,2,1,1,1,4,2 });

	test_rleString_example(testCase, "TTTTT", 1, "T", (const int64_t[]){ 5 });

	test_rleString_example(testCase, "", 0, "", (const int64_t[]){ 1 });

	test_rleString_example(testCase, "TTTTTCC", 2, "TC", (const int64_t[]){ 5, 2 });
}

void test_addInsert(CuTest *testCase) {
	CuAssertStrEquals(testCase, "GATTACA", addInsert("GAACA", "TT", 2));
	CuAssertStrEquals(testCase, "GATTACA", addInsert("", "GATTACA", 0));
	CuAssertStrEquals(testCase, "GATTACA", addInsert("ATTACA", "G", 0));
	CuAssertStrEquals(testCase, "GATTACA", addInsert("GATTAC", "A", 6));
	CuAssertStrEquals(testCase, "GATTACA", addInsert("GATTACA", "", 6));
	CuAssertStrEquals(testCase, "GATTACA", addInsert("GATTACA", "", 3));
}

void test_removeDelete(CuTest *testCase) {
	CuAssertStrEquals(testCase, "GATTACA", removeDelete("GATTGGACA", 2, 4));
	CuAssertStrEquals(testCase, "GATTACA", removeDelete("GATTACA", 0, 0));
	CuAssertStrEquals(testCase, "GATTACA", removeDelete("GATTACATT", 2, 7));
	CuAssertStrEquals(testCase, "GATTACA", removeDelete("AGATTACA", 1, 0));
}

//char *removeDelete(char *string, int64_t deleteLength, int64_t editStart);

/*
 // Crufty code used to generate an initial model for nanopore alignment

typedef enum {
    match = 0, shortGapX = 1, shortGapY = 2, longGapX = 3, longGapY = 4
} State;

static inline double *hmm_getTransition2(Hmm *hmm, int64_t from, int64_t to) {
    return &(hmm->transitions[from * hmm->stateNumber + to]);
}

static inline double *hmm_getEmissionsExpectation2(Hmm *hmm, int64_t state, Symbol x, Symbol y) {
    return &(hmm->emissions[state * SYMBOL_NUMBER_NO_N * SYMBOL_NUMBER_NO_N + x * SYMBOL_NUMBER_NO_N + y]);
}

static void test_hmm(CuTest *testCase) {
	Hmm *hmm = hmm_constructEmpty(0.0, threeState);

	//                 000000   0000111111111122222222223333333333444 444444455555555556666666666777777 7777
	//                 012345   6789012345678901234567890123456789012 345678901234567890123456789012345 6789
	//reference =     "CATTTT   CTCTCCCTCCGTCATTGCACAGGAAAACAGATGAAAA TGCAGGGCAATAATGACCATAAAACGCATTTTT ATTT";
	//trueReference = "CATTTTTCTCTCTCCCTCCGTCATTGCACAGGAAAACAGATGAAAAATGCGGGGCATG  TGACCATAAAACGCATTTTTTATTT";

	//reference =      GATGTAAAAAA   GAAATGATTT     GCTAGAACAGAGCATAAATACACATCTGT
	//trueReference =  GATGTAAAAAAAAAGAAATGA   CGGAAGTTAGAACAGAGCATAAATACACATCTGT

	hmm_getTransition2(hmm, match, match)[0] = 0.9;
	hmm_getTransition2(hmm, match, shortGapX)[0] = 0.05;
	hmm_getTransition2(hmm, match, shortGapY)[0] = 0.05;

	hmm_getTransition2(hmm, shortGapX, match)[0] = 0.5;
	hmm_getTransition2(hmm, shortGapX, shortGapX)[0] = 0.5;
	hmm_getTransition2(hmm, shortGapX, shortGapY)[0] = 0.0;

	hmm_getTransition2(hmm, shortGapY, match)[0] = 0.5;
	hmm_getTransition2(hmm, shortGapY, shortGapX)[0] = 0.0;
	hmm_getTransition2(hmm, shortGapY, shortGapY)[0] = 0.5;

	for(int64_t i=0; i<SYMBOL_NUMBER_NO_N; i++) {
		hmm_getEmissionsExpectation2(hmm, match, i, i)[0] = 0.98;
		for(int64_t j=0; j<SYMBOL_NUMBER_NO_N; j++) {
			if(j != i) {
				hmm_getEmissionsExpectation2(hmm, match, i, j)[0] = 0.02/3;
			}
		}
	}

	// Gaps
	for(int64_t i=0; i<SYMBOL_NUMBER_NO_N; i++) {
		for(int64_t j=0; j<SYMBOL_NUMBER_NO_N; j++) {
			hmm_getEmissionsExpectation2(hmm, shortGapX, i, j)[0] = 0.025;
		}
	}

	for(int64_t i=0; i<SYMBOL_NUMBER_NO_N; i++) {
		for(int64_t j=0; j<SYMBOL_NUMBER_NO_N; j++) {
			hmm_getEmissionsExpectation2(hmm, shortGapY, i, j)[0] = 0.025;
		}
	}

	FILE *fH = fopen("./threeStateNanopore.hmm", "w");
	hmm_write(hmm, fH);
	fclose(fH);
}*/

CuSuite* realignmentTestSuite(void) {
    CuSuite* suite = CuSuiteNew();

    /*SUITE_ADD_TEST(suite, test_poa_getReferenceGraph);
    SUITE_ADD_TEST(suite, test_poa_augment_example);
    SUITE_ADD_TEST(suite, test_poa_realign_tiny_example1);
    SUITE_ADD_TEST(suite, test_poa_realign_example1);
    SUITE_ADD_TEST(suite, test_poa_realign_example2);
    SUITE_ADD_TEST(suite, test_poa_realign_rle_example1);
    SUITE_ADD_TEST(suite, test_poa_realign_rle_example2);
    SUITE_ADD_TEST(suite, test_poa_realign);
    SUITE_ADD_TEST(suite, test_poa_realignIterative);
    SUITE_ADD_TEST(suite, test_getShift);
    SUITE_ADD_TEST(suite, test_rleString_examples);
    SUITE_ADD_TEST(suite, test_addInsert);
    SUITE_ADD_TEST(suite, test_removeDelete);

    SUITE_ADD_TEST(suite, test_poa_realign_examples_no_rle);
    SUITE_ADD_TEST(suite, test_poa_realign_examples_rle);

    SUITE_ADD_TEST(suite, test_poa_realign_messy_examples_no_rle);
    SUITE_ADD_TEST(suite, test_poa_realign_messy_examples_rle);

    SUITE_ADD_TEST(suite, test_poa_realign_examples_large_rle);
    SUITE_ADD_TEST(suite, test_poa_realign_examples_large_no_rle);*/

    SUITE_ADD_TEST(suite, test_poa_realign_examples_long_rle);
    //SUITE_ADD_TEST(suite, test_poa_realign_examples_long_no_rle);

    //SUITE_ADD_TEST(suite, test_poa_realign_examples_very_large_rle);
    //SUITE_ADD_TEST(suite, test_poa_realign_examples_very_large_no_rle);

    return suite;
}
