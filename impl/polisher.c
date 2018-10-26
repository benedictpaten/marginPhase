/*
 * Copyright (C) 2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "stPolish.h"

PoaBaseObservation *poaBaseObservation_construct(int64_t readNo, int64_t offset, double weight) {
	PoaBaseObservation *poaBaseObservation = st_calloc(1, sizeof(PoaBaseObservation));

	poaBaseObservation->readNo = readNo;
	poaBaseObservation->offset = offset;
	poaBaseObservation->weight = weight;

	return poaBaseObservation;
}

void poaBaseObservation_destruct(PoaBaseObservation *poaBaseObservation) {
	free(poaBaseObservation);
}

PoaInsert *poaInsert_construct(char *insert, double weight) {
	PoaInsert *poaInsert = st_calloc(1, sizeof(PoaInsert));

	poaInsert->insert = insert;
	poaInsert->weight = weight;

	return poaInsert;
}

void poaInsert_destruct(PoaInsert *poaInsert) {
	free(poaInsert->insert);
	free(poaInsert);
}

PoaDelete *poaDelete_construct(int64_t length, double weight) {
	PoaDelete *poaDelete = st_calloc(1, sizeof(PoaDelete));

	poaDelete->length = length;
	poaDelete->weight = weight;

	return poaDelete;
}

void poaDelete_destruct(PoaDelete *poaDelete) {
	free(poaDelete);
}

PoaNode *poaNode_construct(char base) {
	PoaNode *poaNode = st_calloc(1, sizeof(PoaNode));

	poaNode->inserts = stList_construct3(0, (void(*)(void *)) poaInsert_destruct);
	poaNode->deletes = stList_construct3(0, (void(*)(void *)) poaDelete_destruct);
	poaNode->base = base;
	poaNode->baseWeights = st_calloc(SYMBOL_NUMBER, sizeof(double)); // Encoded using Symbol enum
	poaNode->observations = stList_construct3(0, (void (*)(void *))poaBaseObservation_destruct);

	return poaNode;
}

void poaNode_destruct(PoaNode *poaNode) {
	stList_destruct(poaNode->inserts);
	stList_destruct(poaNode->deletes);
	stList_destruct(poaNode->observations);
	free(poaNode->baseWeights);
	free(poaNode);
}

Poa *poa_getReferenceGraph(char *reference) {
	Poa *poa = st_calloc(1, sizeof(Poa));

	poa->nodes = stList_construct3(0, (void (*)(void *))poaNode_destruct);
	poa->refString = stString_copy(reference);

	int64_t refLength = strlen(reference);
	stList_append(poa->nodes, poaNode_construct('N')); // Add empty prefix node
	for(int64_t i=0; i<refLength; i++) {
		stList_append(poa->nodes, poaNode_construct(toupper(reference[i])));
	}

	return poa;
}

void poa_destruct(Poa *poa) {
	free(poa->refString);
	stList_destruct(poa->nodes);
	free(poa);
}

int cmpAlignedPairsByCoordinates(const void *a, const void *b) {
	/*
	 * Compares aligned pairs, represented as stIntTuples of the form (weight, x, y) first by
	 * x coordinate and then by y coordinate, ignoring the weight.
	 */
	stIntTuple *one = (stIntTuple *)a, *two = (stIntTuple *)b;

	int i = stIntTuple_get(one, 1) < stIntTuple_get(two, 1) ? -1 : stIntTuple_get(one, 1) > stIntTuple_get(two, 1) ? 1 : 0;
	if(i == 0) {
		i = stIntTuple_get(one, 2) < stIntTuple_get(two, 2) ? -1 : stIntTuple_get(one, 2) > stIntTuple_get(two, 2) ? 1 : 0;
	}

	return i;
}

int cmpAlignedPairsByInvertedCoordinates(const void *a, const void *b) {
	/*
	 * As cmpAlignedPairsByCoordinates, but compares by y coordinate and then x coordinate.
	 */
	stIntTuple *one = (stIntTuple *)a, *two = (stIntTuple *)b;

	int i = stIntTuple_get(one, 2) < stIntTuple_get(two, 2) ? -1 : stIntTuple_get(one, 2) > stIntTuple_get(two, 2) ? 1 : 0;
	if(i == 0) {
		i = stIntTuple_get(one, 1) < stIntTuple_get(two, 1) ? -1 : stIntTuple_get(one, 1) > stIntTuple_get(two, 1) ? 1 : 0;
	}

	return i;
}

bool isMatch(stSortedSet *matchesSet, int64_t x, int64_t y) {
	stIntTuple pair[4];
	pair[0] = 3;
	pair[2] = x;
	pair[3] = y;
	return stSortedSet_search(matchesSet, &pair) != NULL;
}

static void addToInserts(PoaNode *node, char *insert, double weight) {
	/*
	 * Add given insert to node.
	 */

	PoaInsert *poaInsert = NULL;
	// Check if the complete insert is already in the poa graph:
	for(int64_t m=0; m<stList_length(node->inserts); m++) {
		poaInsert = stList_get(node->inserts, m);
		if(stString_eq(poaInsert->insert, insert)) {
			poaInsert->weight += weight;
			return;
		}
	}
	// otherwise add the insert to the poa graph
	stList_append(node->inserts, poaInsert_construct(stString_copy(insert), weight));
}

static void addToDeletes(PoaNode *node, int64_t length, double weight) {
	/*
	 * Add given deletion to node.
	 */

	PoaDelete *poaDelete = NULL;
	// Check if the delete is already in the poa graph:
	for(int64_t m=0; m<stList_length(node->deletes); m++) {
		poaDelete = stList_get(node->deletes, m);
		if(poaDelete->length == length) {
			poaDelete->weight += weight;
			return;
		}
	}
	// otherwise add the delete to the poa graph
	stList_append(node->deletes, poaDelete_construct(length, weight));
}

char *poa_getReferenceSubstring(Poa *poa, int64_t startIndex, int64_t length) {
	/*
	 * Get substring of the reference, starting from given node.
	 */
	char *refSubString = st_malloc(sizeof(char) * (1 + length));
	refSubString[length] = '\0';
	for(int64_t j=0; j<length; j++) {
		PoaNode *node = stList_get(poa->nodes, startIndex+j);
		refSubString[j] = node->base;
	}
	return refSubString;
}

static bool matchesReferenceSubstring(char *refString, int64_t refStart, char *str, int64_t length) {
	/*
	 * Returns true if the given string str matches the given reference substring starting
	 * from the given reference position, refStart
	 */
	for(int64_t l=0; l<length; l++) {
		if(refString[refStart+l] != str[l]) {
			return 0;
		}
	}
	return 1;
}

static bool hasInternalRepeat(char *str, int64_t length, int64_t repeatLength) {
	/*
	 * Establishes if str has an internal repeat of length repeatLength.
	 * e.g. if ATATAT, internal repeat AT (length 2) is an internal repeat, but ATA is not.
	 */
	if(length % repeatLength != 0) { // If not divisible by repeatLength then can not be repeat
		return 0;
	}
	for(int64_t i=repeatLength; i<length; i+=repeatLength) {
		for(int64_t j=0; j<repeatLength; j++) {
			if(str[j] != str[j+i]) {
				return 0;
			}
		}
	}
	return 1;
}

int64_t getShift(char *refString, int64_t refStart, char *str, int64_t length) {
	// Walk back over reference sequence and see if indel can be shifted

	// Establish minimal internal repeat length
	// if ATATAT, minimal internal repeat is AT,
	// similarly if AAAAAAA then minimal internal repeat is A
	int64_t minRepeatLength = 0;
	while(minRepeatLength++ < length) {
		if(hasInternalRepeat(str, length, minRepeatLength)) {
			break;
		}
	}

	// Now walk back by multiples of minimal internal repeat length
	for(int64_t k=refStart-minRepeatLength; k>=0; k-=minRepeatLength) {
		if(!matchesReferenceSubstring(refString, k, str, minRepeatLength)) {
			break;
		}
		refStart = k;
	}

	return refStart;
}

int64_t getMaxCommonSuffixLength(char *str1, int64_t length1, char *str2, int64_t length2) {
	/*
	 * Returns the length of the maximum suffix of the reference string ending at refStart (inclusive)
	 * that is the same as a suffix of str.
	 */
	int64_t i=0;
	while(length1-i-1 >= 0 && length2-i+1 >= 0) {
		if(str1[length1-1-i] != str2[length2-1-i]) {
			break;
		}
		i++;
	}

	return i;
}

char *rotateString(char *str, int64_t length, int64_t rotationLength) {
	/*
	 * Cyclic rotates the string so that the reverse suffix of str of rotationLength is removed and made the prefix
	 * of the returned string.
	 */
	char *str2 = st_calloc(length, sizeof(char));
	for(int64_t i=0; i<rotationLength; i++) {
		str2[i] = str[length-1-i];
	}
	for(int64_t i=0; i<length-rotationLength; i++) {
		str2[i+rotationLength] = str[i];
	}

	return str2;
}

void poa_augment(Poa *poa, char *read, int64_t readNo, stList *matches, stList *inserts, stList *deletes) {
	// Add weights of matches to the POA graph

	// For each match in alignment subgraph identify its corresponding node in the POA graph
	// add the weight of the match to the POA node
	for(int64_t i=0; i<stList_length(matches); i++) {
		stIntTuple *match = stList_get(matches, i);

		PoaNode *node = stList_get(poa->nodes, stIntTuple_get(match, 1)+1); // Get corresponding POA node

		// Add base weight to POA node
		node->baseWeights[symbol_convertCharToSymbol(read[stIntTuple_get(match, 2)])] += stIntTuple_get(match, 0);

		// PoaObservation
		stList_append(node->observations, poaBaseObservation_construct(readNo, stIntTuple_get(match, 2), stIntTuple_get(match, 0)));
	}

	// Create a set of match coordinates

	stSortedSet *matchesSet = stSortedSet_construct3(cmpAlignedPairsByCoordinates, NULL);
	for(int64_t i=0; i<stList_length(matches); i++) {
		stIntTuple *match = stList_get(matches, i);

		assert(stSortedSet_search(matchesSet, match) == NULL);
		stSortedSet_insert(matchesSet, match); // Add to matches
	}
	
	// Add inserts to the POA graph

	// Sort the inserts first by the reference coordinate and then by read coordinate
	stList_sort(inserts, cmpAlignedPairsByCoordinates);

	// Let a complete-insert be a sequence of inserts with the same reference coordinate i
	// and consecutive read coordinates j, j+1, ..., j+n, such that the i,j-1 is a match or equal to (-1,-1) (the beginning)
	// in the alignment subgraph and i+1,j+n+1 is similarly a match or equal to (N, M) (the end of the alignment).

	// Enumerate set of complete inserts
	int64_t readLength = strlen(read), refLength = stList_length(poa->nodes)-1;
	for(int64_t i=0; i<stList_length(inserts);) {

		stIntTuple *insertStart = stList_get(inserts, i); // Start of putative complete-insert

		int64_t j=i+1;
		for(;j<stList_length(inserts); j++) {
			stIntTuple *insertEnd = stList_get(inserts, j); // End of putative complete-insert

			// If they don't have the same reference coordinate then not part of same complete-insert
			if(stIntTuple_get(insertStart, 1) != stIntTuple_get(insertEnd, 1)) {
				break;
			}

			// If they don't form a contiguous sequence of read coordinates then not part of same complete-insert
			if(stIntTuple_get(insertStart, 2) + j - i != stIntTuple_get(insertEnd, 2)) {
				break;
			}
		}

		// At this point i (inclusive) and j (exclusive) form the start and end of a maximal putative
		// complete insert sequence in inserts

		// Now enumerate complete-inserts in this interval

		for(int64_t k=i; k<j; k++) {

			// If k position is not flanked by a preceding match or the beginning then can not be a complete insert
			if(!isMatch(matchesSet, stIntTuple_get(insertStart, 1), stIntTuple_get(insertStart, 2) + k - i - 1) &&
					(stIntTuple_get(insertStart, 1) > -1 || stIntTuple_get(insertStart, 2) + k - i - 1 > -1)) {
				continue;
			}

			for(int64_t l=k; l<j; l++) {

				// If l position is not flanked by a proceeding match or the end then can not be a complete insert
				if(!isMatch(matchesSet, stIntTuple_get(insertStart, 1) + 1, stIntTuple_get(insertStart, 2) + l - i + 1) &&
						(stIntTuple_get(insertStart, 1) + 1 < refLength || stIntTuple_get(insertStart, 2) + l - i + 1 < readLength)) {
					continue;
				}

				// At this point k (inclusive) and l (inclusive) represent a complete-insert

				// Calculate weight and label
				double insertWeight = UINT_MAX;
				int64_t insertLength = l+1-k;
				char insertLabel[insertLength+1];
				insertLabel[insertLength] = '\0';
				for(int64_t m=k; m<l+1; m++) {
					stIntTuple *insert = stList_get(inserts, m);
					insertLabel[m-k] = toupper(read[stIntTuple_get(insert, 2)]);
					insertWeight = insertWeight < stIntTuple_get(insert, 0) ? insertWeight : stIntTuple_get(insert, 0);
				}

				// Get the leftmost node in the poa graph to which the insert will connect

				// First find the left point to which the insert will be connected
				assert(stIntTuple_get(insertStart, 1) >= -1);
				int64_t insertPosition = stIntTuple_get(insertStart, 1)+1;

				// Now walk back over reference sequence and see if insert can be left-shifted
				insertPosition = getShift(poa->refString, insertPosition, insertLabel, insertLength);

				// Finally see if can be shifted by common suffix
				int64_t commonSuffixLength = getMaxCommonSuffixLength(poa->refString, insertPosition, insertLabel, insertLength);
				if(commonSuffixLength > 0) {
					char *newInsertStr = rotateString(insertLabel, insertLength, commonSuffixLength);
					memcpy(insertLabel, newInsertStr, insertLength);
					free(newInsertStr);
					insertPosition -= commonSuffixLength;
				}

				// Add insert to graph at leftmost position
				addToInserts(stList_get(poa->nodes, insertPosition), insertLabel, insertWeight);
			}
		}

		// Increase i to start of next maximal complete-insert
		i = j;
	}

	// Add deletes to the POA graph

	// Sort the deletes first by the read coordinate and then by reference coordinate
	stList_sort(deletes, cmpAlignedPairsByInvertedCoordinates);

	// Analogous to a complete-insert, let a complete-delete be a sequence of deletes with the same read coordinate j
	// and consecutive reference coordinates i, i+1, ..., i+m, such that the i-1,j is a match or equal to (-1,-1) (the beginning)
	// in the alignment subgraph and i+m+1,j+1 is similarly a match or (N, M) (the alignment end).

	// Enumerate set of complete-deletes, adding them to the graph
	for(int64_t i=0; i<stList_length(deletes);) {
		stIntTuple *deleteStart = stList_get(deletes, i); // Start of putative complete-delete

		int64_t j=i+1;
		for(;j<stList_length(deletes); j++) {
			stIntTuple *deleteEnd = stList_get(deletes, j); // End of putative complete-delete

			// If they don't have the same read coordinate then not part of same complete-insert
			if(stIntTuple_get(deleteStart, 2) != stIntTuple_get(deleteEnd, 2)) {
				break;
			}

			// If they don't form a contiguous sequence of read coordinates then not part of same complete-insert
			if(stIntTuple_get(deleteStart, 1) + j - i != stIntTuple_get(deleteEnd, 1)) {
				break;
			}
		}

		// At this point i (inclusive) and j (exclusive) form the start and end of a putative maximal
		// complete-delete sequence in deletes

		// Now enumerate complete-deletes in this interval

		for(int64_t k=i; k<j; k++) {

			// If k position is not flanked by a preceding match or alignment beginning then can not be a complete-delete
			if(!isMatch(matchesSet, stIntTuple_get(deleteStart, 1) + k - i - 1, stIntTuple_get(deleteStart, 2)) &&
					((stIntTuple_get(deleteStart, 1) + k - i - 1 > -1 || stIntTuple_get(deleteStart, 2) > -1))) {
				continue;
			}

			for(int64_t l=k; l<j; l++) {

				// If l position is not flanked by a proceeding match or alignment end then can not be a complete-delete
				if(!isMatch(matchesSet, stIntTuple_get(deleteStart, 1) + l - i + 1, stIntTuple_get(deleteStart, 2) + 1) &&
					(stIntTuple_get(deleteStart, 1) + l - i + 1 < refLength || stIntTuple_get(deleteStart, 2) + 1 < readLength)) {
					continue;
				}

				// At this point k (inclusive) and l (inclusive) represent a complete-delete

				// Delete length
				int64_t deleteLength = l-k+1;

				// Calculate weight
				double deleteWeight = UINT_MAX;
				for(int64_t m=k; m<l+1; m++) {
					stIntTuple *delete = stList_get(deletes, m);
					deleteWeight = deleteWeight < stIntTuple_get(delete, 0) ? deleteWeight : stIntTuple_get(delete, 0);
				}

				// Get the leftmost node in the poa graph to which the delete will connect

				// First find the left point to which the delete would be connected
				assert(stIntTuple_get(deleteStart, 1) >= 0);
				int64_t insertPosition = stIntTuple_get(deleteStart, 1);

				// Get string being deleted
				char *deleteLabel = poa_getReferenceSubstring(poa, insertPosition+1, deleteLength);

				// Now walk back over reference sequence and see if insert can be left-shifted
				insertPosition = getShift(poa->refString, insertPosition, deleteLabel, deleteLength);

				// Finally see if can be shifted by common suffix
				insertPosition -= getMaxCommonSuffixLength(poa->refString, insertPosition, deleteLabel, deleteLength);
				free(deleteLabel);

				// Add delete to graph at leftmost position
				addToDeletes(stList_get(poa->nodes, insertPosition), deleteLength, deleteWeight);
			}
		}

		// Increase i to start of next maximal complete-delete
		i = j;
	}

	// Cleanup
	stSortedSet_destruct(matchesSet);
}

/*
 * Generates a set of anchor alignments for the reads aligned to it.
 * Is used to restrict subsequent alignments.
 */
stList *poa_getAnchorAlignments(Poa *poa, int64_t *poaToConsensusMap, int64_t noOfReads) {

	// Allocate anchor alignments
	stList *anchorAlignments = stList_construct3(0, (void (*)(void *))stList_destruct);
	for(int64_t i=0; i<noOfReads; i++) {
		stList_append(anchorAlignments, stList_construct3(0, (void (*)(void *))stIntTuple_destruct));
	}

	// Walk through the weights of the POA to construct the anchor alignments
	for(int64_t i=1; i<stList_length(poa->nodes); i++) {
		PoaNode *poaNode = stList_get(poa->nodes, i);
		int64_t consensusIndex = poaToConsensusMap[i-1];
		if(consensusIndex != -1) { // Poa reference position is aligned to the consensus
			for(int64_t j=0; j<stList_length(poqNode->observations); j++) {
				poaBaseObservation_*obs = stList_get(poaNode->observations, j);
				if(obs->weight/PAIR_ALIGNMENT_PROB_1 > 0.9) { // High confidence anchor pair
					stList *anchorPairs = stList_get(anchorAlignments, obs->readNo);
					stList_append(anchorPairs, stIntTuple_construct2(consensusIndex, obs->offset));
				}
			}
		}
	}

	return anchorAlignments;
}

Poa *poa_realign(stList *reads, stList *anchorAlignments, char *reference,
			  	 StateMachine *sM, PairwiseAlignmentParameters *p) {
	// Build a reference graph with zero weights
	Poa *poa = poa_getReferenceGraph(reference);

	// For each read
	for(int64_t i=0; i<stList_length(reads); i++) {
		char *read = stList_get(reads, i);

		// Generate set of posterior probabilities for matches, deletes and inserts with respect to reference.
		stList *matches = NULL, *inserts = NULL, *deletes = NULL;

		if(anchorAlignments == NULL) {
			getAlignedPairsWithIndels(sM, reference, read, p, &matches, &deletes, &inserts, 0, 0);
		}
		else {
			stList *anchorPairs = stList_get(anchorAlignments, i); // An alignment anchoring the read to the reference
			getAlignedPairsWithIndelsUsingAnchors(sM, reference, read, anchorPairs, p, &matches, &deletes, &inserts, 0, 0);
		}

		// Add weights, edges and nodes to the poa
		poa_augment(poa, read, i, matches, inserts, deletes);

		// Cleanup
		stList_destruct(matches);
		stList_destruct(inserts);
		stList_destruct(deletes);
	}

	return poa;
}

static int cmpInsertsBySequence(const void *a, const void *b) {
	/*
	 * Compares PoaInserts by weight in ascending order.
	 */
	PoaInsert *one = (PoaInsert *)a, *two = (PoaInsert *)b;
	return strcmp(one->insert, two->insert);
}

static int cmpDeletesByLength(const void *a, const void *b) {
	/*
	 * Compares PoaDelete by weight in ascending order.
	 */
	PoaDelete *one = (PoaDelete *)a, *two = (PoaDelete *)b;
	return one->length < two->length ? -1 : one->length > two->length ? 1 : 0;
}

void poa_sortIndels(Poa *poa) {
	for(int64_t i=0; i<stList_length(poa->nodes); i++) {
		PoaNode *node = stList_get(poa->nodes, i);
		stList_sort(node->inserts, cmpInsertsBySequence);
		stList_sort(node->deletes, cmpDeletesByLength);
	}
}

void poa_normalize(Poa *poa) {
	poa_sortIndels(poa);
}

double poa_getReferenceNodeTotalMatchWeight(Poa *poa) {
	double weight = 0.0;
	for(int64_t i=0; i<stList_length(poa->nodes); i++) {
		PoaNode *node = stList_get(poa->nodes, i);
		weight += node->baseWeights[symbol_convertCharToSymbol(node->base)];
		//for(int64_t j=0; j<SYMBOL_NUMBER; j++) {
		//	weight += node->baseWeights[j];
		//}
	}
	return weight;
}

double poa_getReferenceNodeTotalDisagreementWeight(Poa *poa) {
	double weight = 0.0;
	for(int64_t i=0; i<stList_length(poa->nodes); i++) {
		PoaNode *node = stList_get(poa->nodes, i);
		int64_t refSymbol = symbol_convertCharToSymbol(node->base);
		for(int64_t j=0; j<SYMBOL_NUMBER; j++) {
			if(j != refSymbol) {
				weight += node->baseWeights[j];
			}
		}
	}
	return weight;
}

double poa_getInsertTotalWeight(Poa *poa) {
	double weight = 0.0;
	for(int64_t i=0; i<stList_length(poa->nodes); i++) {
		PoaNode *node = stList_get(poa->nodes, i);
		for(int64_t j=0; j<stList_length(node->inserts); j++) {
			PoaInsert *insert = stList_get(node->inserts, j);
			weight += insert->weight * strlen(insert->insert);
		}
	}
	return weight;
}

double poa_getDeleteTotalWeight(Poa *poa) {
	double weight = 0.0;
	for(int64_t i=0; i<stList_length(poa->nodes); i++) {
		PoaNode *node = stList_get(poa->nodes, i);
		for(int64_t j=0; j<stList_length(node->deletes); j++) {
			PoaDelete *delete = stList_get(node->deletes, j);
			weight += delete->weight * delete->length;
		}
	}
	return weight;
}

double poa_getTotalErrorWeight(Poa *poa)  {
	return poa_getDeleteTotalWeight(poa) + poa_getInsertTotalWeight(poa) + poa_getReferenceNodeTotalDisagreementWeight(poa);
}

void poa_print(Poa *poa, FILE *fH, float indelSignificanceThreshold) {
	// Print info for each base in reference in turn
	for(int64_t i=0; i<stList_length(poa->nodes); i++) {
		PoaNode *node = stList_get(poa->nodes, i);
		fprintf(fH, "%" PRIi64 "\t%c", i, node->base);

		// Bases
		double totalWeight = 0.0;
		for(int64_t j=0; j<SYMBOL_NUMBER; j++) {
			totalWeight += node->baseWeights[j];
		}
		for(int64_t j=0; j<SYMBOL_NUMBER; j++) {
			if(node->baseWeights[j]/totalWeight > 0.25) {
				fprintf(fH, "\t%c:%f (%f)", symbol_convertSymbolToChar(j), (float)node->baseWeights[j]/PAIR_ALIGNMENT_PROB_1, node->baseWeights[j]/totalWeight);
			}
		}
		fprintf(fH, "\tTotal-weight:%f\n", (float)totalWeight/PAIR_ALIGNMENT_PROB_1);

		// Inserts
		for(int64_t j=0; j<stList_length(node->inserts); j++) {
			PoaInsert *insert = stList_get(node->inserts, j);
			if(insert->weight/PAIR_ALIGNMENT_PROB_1 >= indelSignificanceThreshold) {
				fprintf(fH, "Insert\tSeq:%s\tWeight:%f\n", insert->insert, (float)insert->weight/PAIR_ALIGNMENT_PROB_1);
			}
		}

		// Deletes
		for(int64_t j=0; j<stList_length(node->deletes); j++) {
			PoaDelete *delete = stList_get(node->deletes, j);
			if(delete->weight/PAIR_ALIGNMENT_PROB_1 >= indelSignificanceThreshold) {
				fprintf(fH, "Delete\tLength:%" PRIi64 "\tWeight:%f\n", delete->length, (float)delete->weight/PAIR_ALIGNMENT_PROB_1);
			}
		}
	}
}

void poa_printSummaryStats(Poa *poa, FILE *fH) {
	double totalReferenceMatchWeight = poa_getReferenceNodeTotalMatchWeight(poa)/PAIR_ALIGNMENT_PROB_1;
	double totalReferenceMismatchWeight = poa_getReferenceNodeTotalDisagreementWeight(poa)/PAIR_ALIGNMENT_PROB_1;
	double totalInsertWeight = poa_getInsertTotalWeight(poa)/PAIR_ALIGNMENT_PROB_1;
	double totalDeleteWeight = poa_getDeleteTotalWeight(poa)/PAIR_ALIGNMENT_PROB_1;

	fprintf(fH, "Totals, reference match weight: %f reference mismatch weight: %f insert weight: %f delete weight: %f indel weight: %f, sum error: %f\n",
			totalReferenceMatchWeight, totalReferenceMismatchWeight,
			totalInsertWeight, totalDeleteWeight, totalInsertWeight + totalDeleteWeight, totalInsertWeight + totalDeleteWeight + totalReferenceMismatchWeight);
}

char *poa_getConsensus(Poa *poa, int64_t **poaToConsensusMap) {
	// Cheesy profile HMM like algorithm
	// Calculates forward probabilities through model, then
	// traces back through max prob local path, greedily

	// Probabilities/weights we keep track of.

	// Total weight of outgoing transitions
	double *totalOutgoingWeights = st_calloc(stList_length(poa->nodes), sizeof(double));

	// Forward probabilities
	double *nodeForwardLogProbs = st_calloc(stList_length(poa->nodes)+1, sizeof(double));
	// Initialize, only start state has log(1) = 0 prob
	for(int64_t i=1; i<stList_length(poa->nodes)+1; i++) {
		nodeForwardLogProbs[i] = LOG_ZERO;
	}

	// Forward probabilities of transitioning from a node to the its successor without
	// an indel
	double *matchTransitionForwardLogProbs = st_calloc(stList_length(poa->nodes), sizeof(double));

	// Calculate incoming deletions for each node

	stList *incomingDeletions = stList_construct3(0, (void (*)(void *))stList_destruct);
	for(int64_t i=0; i<stList_length(poa->nodes) + 1; i++) {
		stList_append(incomingDeletions, stList_construct());
	}
	for(int64_t i=0; i<stList_length(poa->nodes); i++) {
		PoaNode *node = stList_get(poa->nodes, i);
		for(int64_t j=0; j<stList_length(node->deletes); j++) {
			PoaDelete *delete = stList_get(node->deletes, j);
			stList_append(stList_get(incomingDeletions, i + delete->length + 1), delete);
		}
	}

	// Walk through the graph left-to-right calculating forward probabilities

	for(int64_t i=0; i<stList_length(poa->nodes); i++) {
		PoaNode *node = stList_get(poa->nodes, i);

		// Calculate total weight of indels connecting from this node

		double totalIndelWeight = 0.0;

		for(int64_t j=0; j<stList_length(node->inserts); j++) {
			PoaInsert *insert = stList_get(node->inserts, j);
			totalIndelWeight += insert->weight;
		}
		for(int64_t j=0; j<stList_length(node->deletes); j++) {
			PoaDelete *delete = stList_get(node->deletes, j);
			totalIndelWeight += delete->weight;
		}

		// Calculate the match transition weight of node
		// that is, we make an estimate of the weight/expectation
		// of transitioning from this node to the next node without an indel

		double matchTransitionWeight = 0.0;
		if(i == 0) {
			// Set the initiation probability according to the average base weight
			for(int64_t j=1; j<stList_length(poa->nodes); j++) {
				PoaNode *nNode = stList_get(poa->nodes, j);
				for(int64_t k=0; k<SYMBOL_NUMBER_NO_N; k++) {
					matchTransitionWeight += nNode->baseWeights[k];
				}
			}
			matchTransitionWeight /= stList_length(poa->nodes)-1;
		}
		else {
			for(int64_t j=0; j<SYMBOL_NUMBER; j++) {
				matchTransitionWeight += node->baseWeights[j];
			}
			matchTransitionWeight -= totalIndelWeight;
		}

		// Hack to stop zero weights
		matchTransitionWeight = matchTransitionWeight <= 0 ? 0.0001 : matchTransitionWeight; // Make a small value

		// Calculate the total weight of outgoing transitions
		totalOutgoingWeights[i] = matchTransitionWeight + totalIndelWeight;

		// Update the probabilities of nodes that connect by to this node

		// Inserts
		for(int64_t j=0; j<stList_length(node->inserts); j++) {
			PoaInsert *insert = stList_get(node->inserts, j);
			nodeForwardLogProbs[i+1] = logAdd(nodeForwardLogProbs[i+1],
					nodeForwardLogProbs[i] + log(insert->weight/totalOutgoingWeights[i]));
		}

		// Deletes
		for(int64_t j=0; j<stList_length(node->deletes); j++) {
			PoaDelete *delete = stList_get(node->deletes, j);
			nodeForwardLogProbs[i+delete->length+1] = logAdd(nodeForwardLogProbs[i+delete->length+1],
					nodeForwardLogProbs[i] + log(delete->weight/totalOutgoingWeights[i]));
		}

		// Match
		matchTransitionForwardLogProbs[i] = nodeForwardLogProbs[i] + log(matchTransitionWeight/totalOutgoingWeights[i]);
		nodeForwardLogProbs[i+1] = logAdd(nodeForwardLogProbs[i+1], matchTransitionForwardLogProbs[i]);
	}

	// Now traceback picking consensus greedily

	// Allocate consensus map, setting the alignment of reference
	// string positions initially all to gaps.
	*poaToConsensusMap = st_malloc((stList->length(poa->nodes)-1) * sizeof(int64_t));
	for(int64_t i=0; i<stList_len(poa->node)-1; i++) {
		*poaToConsensusMap[i] = -1;
	}

	stList *consensusStrings = stList_construct3(0, free);
	int64_t runningConsensusLength = 0;

	for(int64_t i=stList_length(poa->nodes); i>0;) {

		//  Add base if not at end
		if(i < stList_length(poa->nodes)) {
			PoaNode *node = stList_get(poa->nodes, i);
			
			// Picks a base, giving a discount to the reference base,
			// because the alignment is biased towards it

			int64_t refBaseIndex = symbol_convertCharToSymbol(node->base);

			double maxBaseWeight = 0;
			int64_t maxBaseIndex = -1;
			for(int64_t j=0; j<SYMBOL_NUMBER; j++) {
				if(j != refBaseIndex && node->baseWeights[j] > maxBaseWeight) {
					maxBaseWeight = node->baseWeights[j];
					maxBaseIndex = j;
				}
			}

			double refBaseWeight = node->baseWeights[refBaseIndex];

			if(refBaseWeight * 0.5 > maxBaseWeight) {
				maxBaseIndex = refBaseIndex;
			}

			stList_append(consensusStrings, stString_print("%c", symbol_convertSymbolToChar(maxBaseIndex)));

			// Update poa to consensus map
			*poaToConsensusMap[i-1] = runningConsensusLength++;
		}

		// Get max insert
		double maxInsertProb = LOG_ZERO;
		double totalInsertProb = LOG_ZERO;
		PoaInsert *maxInsert = NULL;
		PoaNode *pNode = stList_get(poa->nodes, i-1);
		for(int64_t j=0; j<stList_length(pNode->inserts); j++) {
			PoaInsert *insert = stList_get(pNode->inserts, j);
			double p = log(insert->weight/totalOutgoingWeights[i-1]) + nodeForwardLogProbs[i-1];
			if(p > maxInsertProb) {
				maxInsertProb = p;
				maxInsert = insert;
			}
			totalInsertProb = logAdd(totalInsertProb, p);
		}

		// Get max delete
		double maxDeleteProb = LOG_ZERO;
		double totalDeleteProb = LOG_ZERO;
		PoaDelete *maxDelete = NULL;
		stList *incidentDeletes = stList_get(incomingDeletions, i);
		for(int64_t j=0; j<stList_length(incidentDeletes); j++) {
			PoaDelete *delete = stList_get(incidentDeletes, j);
			double p = log(delete->weight/totalOutgoingWeights[i-delete->length-1]) + nodeForwardLogProbs[i-delete->length-1];
			if(p > maxDeleteProb) {
				maxDeleteProb = p;
				maxDelete = delete;
			}
			totalDeleteProb = logAdd(totalDeleteProb, p);
		}

		//fprintf(stderr, "%i Maxinsert: %f Maxdelete %f Match %f\n", (int)i, (float)maxInsertProb, (float)maxDeleteProb, (float)matchTransitionForwardLogProbs[i-1]);

		if(matchTransitionForwardLogProbs[i-1] >= totalDeleteProb && matchTransitionForwardLogProbs[i-1] >= totalInsertProb) {
			// Is likely a match, move back to previous reference base
			i--;
		}
		else if(totalInsertProb >= totalDeleteProb) {
			// Is likely an insert, append insert to consensus string
			// and move to a previous reference base
			stList_append(consensusStrings, stString_copy(maxInsert->insert));
			runningConsensusLength += strlen(maxInsert->insert);
			i--;
		}
		else {
			// Is likely a delete, jump back to skip deleted bases
			i -= maxDelete->length+1;
		}
	}

	// Concatenate backwards to make consensus string
	stList_reverse(consensusStrings);
	char *consensusString = stString_join2("", consensusStrings);
	assert(runningConsensusLength == strlen(consensusString));

	// Now reverse the poaToConsensusMap, because offsets are from  end of string but need them to be from beginning
	for(int64_t i=0; i<stList_length(poa->nodes)-1; i++) {
		if(*poaToConsensusMap[i] != -1) {
			*poaToConsensusMap[i] = runningConsensusLength - 1 - *poaToConsensusMap[i];
		}
	}

	// Cleanup
	stList_destruct(consensusStrings);
	stList_destruct(incomingDeletions);
	free(nodeForwardLogProbs);
	free(matchTransitionForwardLogProbs);
	free(totalOutgoingWeights);

	return consensusString;
}

Poa *poa_realignIterative(stList *reads, stList *anchorAlignments, char *reference,
			  	 StateMachine *sM, PairwiseAlignmentParameters *p) {
	Poa *poa = poa_realign(reads, anchorAlignments, reference, sM, p);
	double score = poa_getReferenceNodeTotalMatchWeight(poa) - poa_getTotalErrorWeight(poa);

	int64_t i=0;
	while(1) {
		int64_t *poaToConsensusMap;
		reference = poa_getConsensus(poa, &poaToConsensusMap);

		// Stop in case consensus string is same as old reference (i.e. greedy convergence)
		if(stString_eq(reference, poa->refString)) {
			free(reference);
			free(poaToConsensusMap);
			break;
		}

		// Get anchor alignments
		stList *anchorAlignments = poa_getAnchorAlignments(poa, poaToConsensusMap, stList_length(reads));

		// Generated updated poa
		Poa *poa2 = poa_realign(reads, anchorAlignments, reference, sM, p);

		// Cleanup
		free(reference);
		free(poaToConsensusMap);
		stList_destruct(anchorAlignments);

		double score2 = poa_getReferenceNodeTotalMatchWeight(poa2) - poa_getTotalErrorWeight(poa2);

		// Stop if score decreases (greedy stopping)
		if(score2 <= score) {
			poa_destruct(poa2);
			break;
		}

		poa_destruct(poa);
		poa = poa2;
		score = score2;
	}

	return poa;
}

char *addInsert(char *string, char *insertString, int64_t editStart) {
	int64_t insertLength = strlen(insertString);
	int64_t stringLength = strlen(string);

	// Allocate new string
	char *editedString = st_malloc(sizeof(char)*(stringLength + insertLength + 1));

	// Make new reference string
	memcpy(editedString, string, sizeof(char) * editStart);
	memcpy(&(editedString[editStart]), insertString, sizeof(char) * insertLength);
	memcpy(&(editedString[editStart + insertLength]), &(string[editStart]),
			sizeof(char) * (1 + stringLength - editStart));

	return editedString;
}

char *removeDelete(char *string, int64_t deleteLength, int64_t editStart) {
	int64_t stringLength = strlen(string);

	// Allocate new string
	assert(stringLength >= deleteLength);
	char *editedString = st_malloc(sizeof(char)*(stringLength - deleteLength + 1));

	// Make new reference string
	memcpy(editedString, string, sizeof(char)*editStart);
	memcpy(&(editedString[editStart]), &(string[editStart+deleteLength]),
			sizeof(char)*(1 + stringLength - editStart - deleteLength));

	return editedString;
}

Poa *poa_checkMajorIndelEditsGreedily(Poa *poa, stList *reads, StateMachine *sM, PairwiseAlignmentParameters *p) {
	double score = poa_getReferenceNodeTotalMatchWeight(poa) - poa_getTotalErrorWeight(poa);

	while(1) {
		int64_t insertStart = 0;
		PoaInsert *maxInsert = NULL;
		int64_t deleteStart = 0;
		PoaDelete *maxDelete = NULL;

		// Get highest value indel
		for(int64_t i=0; i<stList_length(poa->nodes); i++) {
			PoaNode *node = stList_get(poa->nodes, i);

			// Get max insert
			for(int64_t j=0; j<stList_length(node->inserts); j++) {
				PoaInsert *insert = stList_get(node->inserts, j);
				if(maxInsert == NULL || insert->weight > maxInsert->weight) {
					maxInsert = insert;
					insertStart = i;
				}
			}

			// Get max delete
			for(int64_t j=0; j<stList_length(node->deletes); j++) {
				PoaDelete *delete = stList_get(node->deletes, j);
				if(maxDelete == NULL || delete->weight > maxDelete->weight) {
					maxDelete = delete;
					deleteStart = i;
				}
			}
		}

		if(maxInsert == NULL || maxDelete == NULL) {
			return poa;
		}

		// Create new graph with edit
		char *editRef = maxInsert->weight >= maxDelete->weight ? addInsert(poa->refString, maxInsert->insert, insertStart) :
				removeDelete(poa->refString, maxDelete->length, deleteStart);
		Poa *poa2 = poa_realign(reads, editRef, sM, p);
		free(editRef);
		double score2 = poa_getReferenceNodeTotalMatchWeight(poa2) - poa_getTotalErrorWeight(poa2);

		// If new graph has better score, keep it, otherwise return old graph
		if(score2 <= score) {
			poa_destruct(poa2);
			return poa;
		}

		// We got a better graph, so keep it
		poa_destruct(poa);
		poa = poa2;
		score = score2;
	}
}

/*
 * Functions for run-length encoding/decoding with POAs
 */

RleString *rleString_construct(char *str) {
	RleString *rleString = st_calloc(1, sizeof(RleString));

	int64_t strLength = strlen(str);

	// Calc length of rle'd str
	for(int64_t i=0; i<strLength; i++) {
		if(i+1 == strLength || str[i] != str[i+1]) {
			rleString->length++;
		}
	}

	// Allocate
	rleString->rleString = st_calloc(rleString->length+1, sizeof(char));
	rleString->repeatCounts = st_calloc(rleString->length, sizeof(int64_t));

	// Fill out
	int64_t j=0, k=1;
	for(int64_t i=0; i<strLength; i++) {
		if(i+1 == strLength || str[i] != str[i+1]) {
			rleString->rleString[j] = str[i];
			rleString->repeatCounts[j++] = k;
			k=1;
		}
		else {
			k++;
		}
	}
	assert(j == rleString->length);

	return rleString;
}

void rleString_destruct(RleString *rleString) {
	free(rleString->rleString);
	free(rleString->repeatCounts);
	free(rleString);
}

static char *expandRLEConsensus2(PoaNode *node, stList *rleReads, RepeatSubMatrix *repeatSubMatrix) {
	// Pick the base
	double maxBaseWeight = node->baseWeights[0];
	int64_t maxBaseIndex = 0;
	for(int64_t j=1; j<SYMBOL_NUMBER; j++) {
		if(node->baseWeights[j] > maxBaseWeight) {
			maxBaseWeight = node->baseWeights[j];
			maxBaseIndex = j;
		}
	}
	char base = symbol_convertSymbolToChar(maxBaseIndex);

	// Repeat count
	double logProbability;
	int64_t repeatCount = repeatSubMatrix_getMLRepeatCount(repeatSubMatrix, maxBaseIndex, node->observations,
			rleReads, &logProbability);

	// Make repeat string
	char *str = st_calloc(repeatCount+1, sizeof(char));
	for(int64_t j=0; j<repeatCount; j++) {
		str[j] = base;
	}
	str[repeatCount] = '\0';

	return str;
}

char *expandRLEConsensus(Poa *poa, stList *rleReads, RepeatSubMatrix *repeatSubMatrix) {
	stList *consensus = stList_construct3(0, free);
	for(int64_t i=1; i<stList_length(poa->nodes); i++) {
		stList_append(consensus, expandRLEConsensus2(stList_get(poa->nodes, i),
				rleReads, repeatSubMatrix));
	}
	char *consensusString = stString_join2("", consensus);
	stList_destruct(consensus);
	return consensusString;
}

/*
 * Functions for modeling repeat counts
 */

double *repeatSubMatrix_setLogProb(RepeatSubMatrix *repeatSubMatrix, Symbol base, int64_t observedRepeatCount, int64_t underlyingRepeatCount) {
	return &(repeatSubMatrix->logProbabilities[base * repeatSubMatrix->maximumRepeatLength * repeatSubMatrix->maximumRepeatLength +
											 underlyingRepeatCount * repeatSubMatrix->maximumRepeatLength + observedRepeatCount]);
}

double repeatSubMatrix_getLogProb(RepeatSubMatrix *repeatSubMatrix, Symbol base, int64_t observedRepeatCount, int64_t underlyingRepeatCount) {
	return *repeatSubMatrix_setLogProb(repeatSubMatrix, base, observedRepeatCount, underlyingRepeatCount);
}

RepeatSubMatrix *repeatSubMatrix_parse(char *fileName) {
	RepeatSubMatrix *repeatSubMatrix = st_calloc(1, sizeof(RepeatSubMatrix));

	repeatSubMatrix->maximumRepeatLength = 51;
	repeatSubMatrix->logProbabilities = st_calloc(SYMBOL_NUMBER_NO_N *
			repeatSubMatrix->maximumRepeatLength * repeatSubMatrix->maximumRepeatLength, sizeof(double));

	FILE *fh = fopen(fileName, "r");

	for(int64_t i=0; i<SYMBOL_NUMBER_NO_N; i++) {
		char *header = stFile_getLineFromFile(fh);
		Symbol baseIndex = symbol_convertCharToSymbol(header[1]);
		char *probsStr = stFile_getLineFromFile(fh);
		stList *tokens = stString_split(probsStr);
		for(int64_t j=0; j<repeatSubMatrix->maximumRepeatLength; j++) {
			for(int64_t k=0; k<repeatSubMatrix->maximumRepeatLength; k++) {
				char *logProbStr = stList_get(tokens, j*repeatSubMatrix->maximumRepeatLength + k);
				int l = sscanf(logProbStr, "%lf", repeatSubMatrix_setLogProb(repeatSubMatrix, baseIndex, k, j));
				(void)l;
				assert(l == 1);
			}
		}
		// Cleanup
		free(header);
		free(probsStr);
		stList_destruct(tokens);
	}

	fclose(fh);

	return repeatSubMatrix;
}

void repeatSubMatrix_destruct(RepeatSubMatrix *repeatSubMatrix) {
	free(repeatSubMatrix->logProbabilities);
	free(repeatSubMatrix);
}

double repeatSubMatrix_getLogProbForGivenRepeatCount(RepeatSubMatrix *repeatSubMatrix, Symbol base, stList *observations,
												     stList *rleReads, int64_t underlyingRepeatCount) {
	double logProb = LOG_ONE;
	for(int64_t i=0; i<stList_length(observations); i++) {
		PoaBaseObservation *observation = stList_get(observations, i);
		RleString *read = stList_get(rleReads, observation->readNo);
		int64_t observedRepeatCount = read->repeatCounts[observation->offset];
		assert(underlyingRepeatCount < repeatSubMatrix->maximumRepeatLength);
		logProb += repeatSubMatrix_getLogProb(repeatSubMatrix, base, observedRepeatCount, underlyingRepeatCount) * observation->weight;
	}

	return logProb;
}

int64_t repeatSubMatrix_getMLRepeatCount(RepeatSubMatrix *repeatSubMatrix, Symbol base, stList *observations,
		stList *rleReads, double *logProbability) {
	assert(stList_length(observations) > 0);

	// Get the range or repeat observations, used to avoid calculating all repeat lengths, heuristically
	int64_t minRepeatLength=repeatSubMatrix->maximumRepeatLength-1, maxRepeatLength=0; // Mins and maxs inclusive
	for(int64_t i=0; i<stList_length(observations); i++) {
		PoaBaseObservation *observation = stList_get(observations, i);
		RleString *read = stList_get(rleReads, observation->readNo);
		int64_t observedRepeatCount = read->repeatCounts[observation->offset];
		assert(observedRepeatCount < repeatSubMatrix->maximumRepeatLength);
		if(observedRepeatCount < minRepeatLength) {
			minRepeatLength = observedRepeatCount;
		}
		if(observedRepeatCount > maxRepeatLength) {
			maxRepeatLength = observedRepeatCount;
		}
	}

	// Calc the range of repeat observations
	double mlLogProb = repeatSubMatrix_getLogProbForGivenRepeatCount(repeatSubMatrix, base, observations, rleReads, minRepeatLength);
	int64_t mlRepeatLength = minRepeatLength;
	for(int64_t i=minRepeatLength+1; i<maxRepeatLength+1; i++) {
		double p = repeatSubMatrix_getLogProbForGivenRepeatCount(repeatSubMatrix, base, observations, rleReads, i);
		if(p > mlLogProb) {
			mlLogProb = p;
			mlRepeatLength = i;
		}
	}
	*logProbability = mlLogProb;
	return mlRepeatLength;
}

