/*
 * Copyright (C) 2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "margin.h"

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

PoaInsert *poaInsert_construct(char *insert, double weight, bool strand) {
	PoaInsert *poaInsert = st_calloc(1, sizeof(PoaInsert));

	poaInsert->insert = insert;
	if(strand) {
		poaInsert->weightForwardStrand = weight;
	}
	else {
		poaInsert->weightReverseStrand = weight;
	}

	return poaInsert;
}

void poaInsert_destruct(PoaInsert *poaInsert) {
	free(poaInsert->insert);
	free(poaInsert);
}

static bool isBalanced(double forwardStrandWeight, double reverseStrandWeight, double balanceRatio) {
	return balanceRatio == 0 || (forwardStrandWeight > reverseStrandWeight/balanceRatio && reverseStrandWeight > forwardStrandWeight/balanceRatio);
}

double poaInsert_getWeight(PoaInsert *insert) {
	return insert->weightForwardStrand + insert->weightReverseStrand;
	//return isBalanced(insert->weightForwardStrand, insert->weightReverseStrand, 100) ? insert->weightForwardStrand + insert->weightReverseStrand : 0.0;
}

PoaDelete *poaDelete_construct(int64_t length, double weight, bool strand) {
	PoaDelete *poaDelete = st_calloc(1, sizeof(PoaDelete));

	poaDelete->length = length;
	if(strand) {
		poaDelete->weightForwardStrand = weight;
	}
	else {
		poaDelete->weightReverseStrand = weight;
	}

	return poaDelete;
}

void poaDelete_destruct(PoaDelete *poaDelete) {
	free(poaDelete);
}

double poaDelete_getWeight(PoaDelete *delete) {
	return delete->weightForwardStrand + delete->weightReverseStrand;
	//return isBalanced(delete->weightForwardStrand, delete->weightReverseStrand, 100) ? delete->weightForwardStrand + delete->weightReverseStrand : 0.0;
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

static void addToInserts(PoaNode *node, char *insert, double weight, bool strand) {
	/*
	 * Add given insert to node.
	 */

	PoaInsert *poaInsert = NULL;
	// Check if the complete insert is already in the poa graph:
	for(int64_t m=0; m<stList_length(node->inserts); m++) {
		poaInsert = stList_get(node->inserts, m);
		if(stString_eq(poaInsert->insert, insert)) {
			if(strand) {
				poaInsert->weightForwardStrand += weight;
			}
			else {
				poaInsert->weightReverseStrand += weight;
			}
			return;
		}
	}
	// otherwise add the insert to the poa graph
	stList_append(node->inserts, poaInsert_construct(stString_copy(insert), weight, strand));
}

static void addToDeletes(PoaNode *node, int64_t length, double weight, bool strand) {
	/*
	 * Add given deletion to node.
	 */

	PoaDelete *poaDelete = NULL;
	// Check if the delete is already in the poa graph:
	for(int64_t m=0; m<stList_length(node->deletes); m++) {
		poaDelete = stList_get(node->deletes, m);
		if(poaDelete->length == length) {
			if(strand) {
				poaDelete->weightForwardStrand += weight;
			}
			else {
				poaDelete->weightReverseStrand += weight;
			}
			return;
		}
	}
	// otherwise add the delete to the poa graph
	stList_append(node->deletes, poaDelete_construct(length, weight, strand));
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
	while(length1-i-1 >= 0 && length2-i-1 >= 0) {
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

void poa_augment(Poa *poa, char *read, bool readStrand, int64_t readNo, stList *matches, stList *inserts, stList *deletes) {
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
				stIntTuple_get(insertStart, 2) + k - i - 1 > -1) {
				continue;
			}

			for(int64_t l=k; l<j; l++) {

				// If l position is not flanked by a proceeding match or the end then can not be a complete insert
				if(!isMatch(matchesSet, stIntTuple_get(insertStart, 1) + 1, stIntTuple_get(insertStart, 2) + l - i + 1) &&
					stIntTuple_get(insertStart, 2) + l - i + 1 < readLength) {
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
				addToInserts(stList_get(poa->nodes, insertPosition), insertLabel, insertWeight, readStrand);
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

			// If k position is not flanked by a preceding match or the alignment beginning then can not be a complete-delete
			if(!isMatch(matchesSet, stIntTuple_get(deleteStart, 1) + k - i - 1, stIntTuple_get(deleteStart, 2)) &&
					stIntTuple_get(deleteStart, 1) + k - i - 1 > -1) {
				continue;
			}

			for(int64_t l=k; l<j; l++) {

				// If l position is not flanked by a proceeding match or the alignment end then can not be a complete-delete
				if(!isMatch(matchesSet, stIntTuple_get(deleteStart, 1) + l - i + 1, stIntTuple_get(deleteStart, 2) + 1) &&
					stIntTuple_get(deleteStart, 1) + l - i + 1 < refLength) {
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
				assert(stIntTuple_get(deleteStart, 1) + k - i >= 0);
				int64_t insertPosition = stIntTuple_get(deleteStart, 1) + k - i;

				// Get string being deleted
				char *deleteLabel = poa_getReferenceSubstring(poa, insertPosition+1, deleteLength);

				// Now walk back over reference sequence and see if insert can be left-shifted
				insertPosition = getShift(poa->refString, insertPosition, deleteLabel, deleteLength);

				// Finally see if can be shifted by common suffix
				insertPosition -= getMaxCommonSuffixLength(poa->refString, insertPosition, deleteLabel, deleteLength);
				free(deleteLabel);

				// Add delete to graph at leftmost position
				addToDeletes(stList_get(poa->nodes, insertPosition), deleteLength, deleteWeight, readStrand);
			}
		}

		// Increase i to start of next maximal complete-delete
		i = j;
	}

	// Cleanup
	stSortedSet_destruct(matchesSet);
}

stList *poa_getAnchorAlignments(Poa *poa, int64_t *poaToConsensusMap, int64_t noOfReads, PolishParams *pp) {

	// Allocate anchor alignments
	stList *anchorAlignments = stList_construct3(0, (void (*)(void *))stList_destruct);
	for(int64_t i=0; i<noOfReads; i++) {
		stList_append(anchorAlignments, stList_construct3(0, (void (*)(void *))stIntTuple_destruct));
	}

	// Walk through the weights of the POA to construct the anchor alignments
	for(int64_t i=1; i<stList_length(poa->nodes); i++) {
		PoaNode *poaNode = stList_get(poa->nodes, i);
		// If consensus index is null, then just use the underlying reference sequence of the POA.
		int64_t consensusIndex = poaToConsensusMap == NULL ? i-1 : poaToConsensusMap[i-1];
		if(consensusIndex != -1) { // Poa reference position is aligned to the consensus
			for(int64_t j=0; j<stList_length(poaNode->observations); j++) {
				PoaBaseObservation *obs = stList_get(poaNode->observations, j);
				double normalizedObsWeight = obs->weight/PAIR_ALIGNMENT_PROB_1;

				if(normalizedObsWeight > pp->minPosteriorProbForAlignmentAnchors[0]) { // High confidence anchor pair
					stList *anchorPairs = stList_get(anchorAlignments, obs->readNo);

					// Figure out the exact diagonal expansion
					int64_t expansion = (int64_t) pp->minPosteriorProbForAlignmentAnchors[1];
					for(int64_t i=2; i<pp->minPosteriorProbForAlignmentAnchorsLength; i+=2) {
						if(normalizedObsWeight >= pp->minPosteriorProbForAlignmentAnchors[i]) {
							expansion = (int64_t) pp->minPosteriorProbForAlignmentAnchors[i+1];
						}
						else {
							break;
						}
					}

					// The following is masking an underlying bug that allows for multiple high confidence
					// alignments per position
					if(stList_length(anchorPairs) == 0) {
						stList_append(anchorPairs, stIntTuple_construct3(consensusIndex, obs->offset, expansion));
					}
					else {
						stIntTuple *pPair = stList_peek(anchorPairs);

						if(stIntTuple_get(pPair, 0) < consensusIndex && stIntTuple_get(pPair, 1) < obs->offset) {
							stList_append(anchorPairs, stIntTuple_construct3(consensusIndex, obs->offset, expansion));
						}
						//else {
						//	fprintf(stderr, "Ooops read: %i x1: %i y1: %i x2: %i y2: %i %f\n", (int)obs->readNo,
						//			(int)stIntTuple_get(pPair, 0), (int)stIntTuple_get(pPair, 1),
						//			(int)consensusIndex, (int)obs->offset, (float)obs->weight/PAIR_ALIGNMENT_PROB_1);
						//}
					}

				}
			}
		}
	}

	return anchorAlignments;
}

static void adjustAnchors(stList *anchorPairs, int64_t index, int64_t adjustment) {
	for(int64_t i=0; i<stList_length(anchorPairs); i++) {
		stIntTuple *pair = stList_get(anchorPairs, i);
		((int64_t *)pair)[index+1] += adjustment;
	}
}

/*
 * Generates aligned pairs and indel probs, but first crops reference to only include sequence from first
 * to last anchor position.
 */
void getAlignedPairsWithIndelsCroppingReference(char *reference, int64_t refLength,
		char *read, stList *anchorPairs,
		stList **matches, stList **inserts, stList **deletes, PolishParams *polishParams) {
	// Crop reference, to avoid long unaligned prefix and suffix
	// that generates a lot of delete pairs

	// Get cropping coordinates
	// TODO I think we may want to extend refStart and refEnd by the length of the read before and after the first and last aligned positions
	int64_t firstRefPosition, endRefPosition;
	if(stList_length(anchorPairs) > 0) {
		stIntTuple *fPair = stList_get(anchorPairs, 0);
		firstRefPosition = stIntTuple_get(fPair, 0) - stIntTuple_get(fPair, 1);
		firstRefPosition = firstRefPosition < 0 ? 0 : firstRefPosition;

		stIntTuple *lPair = stList_peek(anchorPairs);
		endRefPosition = 1 + stIntTuple_get(lPair, 0) + (strlen(read) - stIntTuple_get(lPair, 1));
		endRefPosition = endRefPosition > refLength ? refLength : endRefPosition;
	}
	else {
		firstRefPosition = 0;
		endRefPosition = refLength;
	}
	assert(firstRefPosition < refLength && firstRefPosition >= 0);
	assert(endRefPosition <= refLength && endRefPosition >= 0);

	// Adjust anchor positions
	adjustAnchors(anchorPairs, 0, -firstRefPosition);

	// Crop reference
	char c = reference[endRefPosition];
	reference[endRefPosition] = '\0';

	// Get alignment
	getAlignedPairsWithIndelsUsingAnchors(polishParams->sM, &(reference[firstRefPosition]), read,
										  anchorPairs, polishParams->p, matches, deletes, inserts, 0, 0);
	//TODO are the delete and insert lists inverted here?

	// De-crop reference
	reference[endRefPosition] = c;

	// Adjust back anchors
	adjustAnchors(anchorPairs, 0, firstRefPosition);

	// Shift matches/inserts/deletes
	adjustAnchors(*matches, 1, firstRefPosition);
	adjustAnchors(*inserts, 1, firstRefPosition);
	adjustAnchors(*deletes, 1, firstRefPosition);
}

Poa *poa_realign(stList *bamChunkReads, stList *anchorAlignments, char *reference, PolishParams *polishParams) {
	// Build a reference graph with zero weights
	Poa *poa = poa_getReferenceGraph(reference);
	int64_t refLength = stList_length(poa->nodes)-1;

	// For each read
	for(int64_t i=0; i<stList_length(bamChunkReads); i++) {
        BamChunkRead *chunkRead = stList_get(bamChunkReads, i);

		// Generate set of posterior probabilities for matches, deletes and inserts with respect to reference.
		stList *matches = NULL, *inserts = NULL, *deletes = NULL;

		if(anchorAlignments == NULL) {
			getAlignedPairsWithIndels(polishParams->sM, reference, chunkRead->nucleotides, polishParams->p,
                                      &matches, &deletes, &inserts, 0, 0);
		}
		else {
			getAlignedPairsWithIndelsCroppingReference(reference, refLength, chunkRead->nucleotides, stList_get(anchorAlignments, i),
                                                       &matches, &inserts, &deletes, polishParams);
		}

		// Add weights, edges and nodes to the poa
		poa_augment(poa, chunkRead->nucleotides, chunkRead->forwardStrand, i, matches, inserts, deletes);

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

double poa_getReferenceNodeTotalMatchWeight(Poa *poa) {
	double weight = 0.0;
	for(int64_t i=0; i<stList_length(poa->nodes); i++) {
		PoaNode *node = stList_get(poa->nodes, i);
		weight += node->baseWeights[symbol_convertCharToSymbol(node->base)];
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
			weight += poaInsert_getWeight(insert) * strlen(insert->insert);
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
			weight += poaDelete_getWeight(delete) * delete->length;
		}
	}
	return weight;
}

double poa_getTotalErrorWeight(Poa *poa)  {
	return poa_getDeleteTotalWeight(poa) + poa_getInsertTotalWeight(poa) + poa_getReferenceNodeTotalDisagreementWeight(poa);
}

double *poaNode_getStrandSpecificBaseWeights(PoaNode *node, stList *bamChunkReads,
		double *totalWeight, double *totalPositiveWeight, double *totalNegativeWeight) {
	/*
	 * Calculate strand specific base weights.
	 */
	*totalWeight = 0.0;
	*totalPositiveWeight = 0.0;
	*totalNegativeWeight = 0.0;
	double *baseWeights = st_calloc(SYMBOL_NUMBER*2, sizeof(double));
	for(int64_t i=0; i<stList_length(node->observations); i++) {
		PoaBaseObservation *baseObs = stList_get(node->observations, i);
		*totalWeight += baseObs->weight;
		BamChunkRead *read = stList_get(bamChunkReads, baseObs->readNo);
		char base = read->nucleotides[baseObs->offset];
		baseWeights[symbol_convertCharToSymbol(base) * 2 + (read->forwardStrand ? 1 : 0)] += baseObs->weight;
		if(read->forwardStrand) {
			*totalPositiveWeight += baseObs->weight;
		}
		else {
			*totalNegativeWeight += baseObs->weight;
		}
	}

	return baseWeights;
}

void poa_printRepeatCounts(Poa *poa, FILE *fH,
		stList *rleReads, stList *bamChunkReads) {
		fprintf(fH, "REF_INDEX\tREF_BASE");
		fprintf(fH, "\tREPEAT_COUNT_OBSxN(READ_BASE:READ_STRAND:REPEAT_COUNT,WEIGHT)\n");

		// Print info for each base in reference in turn
		for(int64_t i=0; i<stList_length(poa->nodes); i++) {
			PoaNode *node = stList_get(poa->nodes, i);

			fprintf(fH, "%" PRIi64 "\t%c", i, node->base);


			for(int64_t j=0; j<stList_length(node->observations); j++) {
				PoaBaseObservation *obs = stList_get(node->observations, j);
				RleString *rleRead = stList_get(rleReads, obs->readNo);
				BamChunkRead *bamChunkRead = stList_get(bamChunkReads, obs->readNo);
				int64_t repeatCount = rleRead->repeatCounts[obs->offset];
				char base = rleRead->rleString[obs->offset];
				fprintf(fH, "\t%c%c%" PRIi64 ",%.3f", base, bamChunkRead->forwardStrand ? '+' : '-', repeatCount, obs->weight/PAIR_ALIGNMENT_PROB_1);
			}

			fprintf(fH, "\n");
		}
}

void poa_printTSV(Poa *poa, FILE *fH,
		stList *bamChunkReads,
		float indelSignificanceThreshold, float strandBalanceRatio) {

	fprintf(fH, "REF_INDEX\tREF_BASE\tTOTAL_WEIGHT\tPOS_STRAND_WEIGHT\tNEG_STRAND_WEIGHT");
	for(int64_t j=0; j<SYMBOL_NUMBER; j++) {
		char c = symbol_convertSymbolToChar(j);
		fprintf(fH, "\tNORM_BASE_%c_WEIGHT\tNORM_POS_STRAND_BASE_%c_WEIGHT\tNORM_NEG_BASE_%c_WEIGHT", c, c, c);
	}

	fprintf(fH, "\tINSERTSxN(INSERT_SEQ, TOTAL_WEIGHT, TOTAL_POS_STRAND_WEIGHT, TOTAL_NEG_STRAND_WEIGHT)\tDELETESxN(DELETE_LENGTH, TOTAL_WEIGHT, TOTAL_POS_STRAND_WEIGHT, TOTAL_NEG_STRAND_WEIGHT)\n");

	// Print info for each base in reference in turn
	for(int64_t i=0; i<stList_length(poa->nodes); i++) {
		PoaNode *node = stList_get(poa->nodes, i);

		// Print base weights first

		// Calculate strand specific base weights
		double totalWeight, totalPositiveWeight, totalNegativeWeight;
		double *baseWeights = poaNode_getStrandSpecificBaseWeights(node, bamChunkReads,
				&totalWeight, &totalPositiveWeight, &totalNegativeWeight);
		
		fprintf(fH, "%" PRIi64 "\t%c\t%f\t%f\t%f", i, node->base,
				totalWeight/PAIR_ALIGNMENT_PROB_1, totalPositiveWeight/PAIR_ALIGNMENT_PROB_1,
				totalNegativeWeight/PAIR_ALIGNMENT_PROB_1);

		for(int64_t j=0; j<SYMBOL_NUMBER; j++) {
			double positiveStrandBaseWeight = baseWeights[j*2 + 1];
			double negativeStrandBaseWeight = baseWeights[j*2 + 0];
			double totalBaseWeight = positiveStrandBaseWeight + negativeStrandBaseWeight;

			fprintf(fH, "\t%f\t%f\t%f,", node->baseWeights[j]/totalWeight,
					positiveStrandBaseWeight/totalPositiveWeight, negativeStrandBaseWeight/totalNegativeWeight);
		}

		free(baseWeights);

		// Inserts
		for(int64_t j=0; j<stList_length(node->inserts); j++) {
			PoaInsert *insert = stList_get(node->inserts, j);
			if(poaInsert_getWeight(insert)/PAIR_ALIGNMENT_PROB_1 >= indelSignificanceThreshold &&
					isBalanced(insert->weightForwardStrand, insert->weightReverseStrand, strandBalanceRatio)) {
				fprintf(fH, "\tINSERT\t%s\t%f\t%f\t%f",
						insert->insert,
						(float)poaInsert_getWeight(insert)/PAIR_ALIGNMENT_PROB_1,
						(float)insert->weightForwardStrand/PAIR_ALIGNMENT_PROB_1,
						(float)insert->weightReverseStrand/PAIR_ALIGNMENT_PROB_1);
			}
		}

		// Deletes
		for(int64_t j=0; j<stList_length(node->deletes); j++) {
			PoaDelete *delete = stList_get(node->deletes, j);
			if(poaDelete_getWeight(delete)/PAIR_ALIGNMENT_PROB_1 >= indelSignificanceThreshold &&
					isBalanced(delete->weightForwardStrand, delete->weightReverseStrand, strandBalanceRatio)) {
				fprintf(fH, "\tDELETE\t%" PRIi64 "\t%f\t%f\t%f",
						delete->length,
						(float)poaDelete_getWeight(delete)/PAIR_ALIGNMENT_PROB_1,
						(float)delete->weightForwardStrand/PAIR_ALIGNMENT_PROB_1,
						(float)delete->weightReverseStrand/PAIR_ALIGNMENT_PROB_1);
			}
		}

		fprintf(fH, "\n");
	}
}

void poa_print(Poa *poa, FILE *fH,
		stList *bamChunkReads,
		float indelSignificanceThreshold, float strandBalanceRatio) {
	// Print info for each base in reference in turn
	for(int64_t i=0; i<stList_length(poa->nodes); i++) {
		PoaNode *node = stList_get(poa->nodes, i);

		// Calculate strand specific base weights
		double totalWeight, totalPositiveWeight, totalNegativeWeight;
		double *baseWeights = poaNode_getStrandSpecificBaseWeights(node, bamChunkReads,
											&totalWeight, &totalPositiveWeight, &totalNegativeWeight);

		fprintf(fH, "%" PRIi64 "\t%c total-weight:%f\ttotal-pos-weight:%f\ttotal-neg-weight:%f", i, node->base,
				totalWeight/PAIR_ALIGNMENT_PROB_1, totalPositiveWeight/PAIR_ALIGNMENT_PROB_1, totalNegativeWeight/PAIR_ALIGNMENT_PROB_1);

		for(int64_t j=0; j<SYMBOL_NUMBER; j++) {
			double positiveStrandBaseWeight = baseWeights[j*2 + 1];
			double negativeStrandBaseWeight = baseWeights[j*2 + 0];
			double totalBaseWeight = positiveStrandBaseWeight + negativeStrandBaseWeight;

			if(totalBaseWeight/totalWeight > 0.25) {
				fprintf(fH, "\t%c:%f (%f) +str:%f, -str:%f,", symbol_convertSymbolToChar(j),
						(float)node->baseWeights[j]/PAIR_ALIGNMENT_PROB_1, node->baseWeights[j]/totalWeight,
						positiveStrandBaseWeight/totalPositiveWeight, negativeStrandBaseWeight/totalNegativeWeight);
			}
		}
		fprintf(fH, "\tTotal-weight:%f\n", (float)totalWeight/PAIR_ALIGNMENT_PROB_1);

		free(baseWeights);

		// Inserts
		for(int64_t j=0; j<stList_length(node->inserts); j++) {
			PoaInsert *insert = stList_get(node->inserts, j);
			if(poaInsert_getWeight(insert)/PAIR_ALIGNMENT_PROB_1 >= indelSignificanceThreshold &&
					isBalanced(insert->weightForwardStrand, insert->weightReverseStrand, strandBalanceRatio)) {
				fprintf(fH, "Insert\tSeq:%s\tTotal weight:%f\tForward Strand Weight:%f\tReverse Strand Weight:%f\n", insert->insert,
						(float)poaInsert_getWeight(insert)/PAIR_ALIGNMENT_PROB_1,
						(float)insert->weightForwardStrand/PAIR_ALIGNMENT_PROB_1,
						(float)insert->weightReverseStrand/PAIR_ALIGNMENT_PROB_1);
			}
		}

		// Deletes
		for(int64_t j=0; j<stList_length(node->deletes); j++) {
			PoaDelete *delete = stList_get(node->deletes, j);
			if(poaDelete_getWeight(delete)/PAIR_ALIGNMENT_PROB_1 >= indelSignificanceThreshold &&
					isBalanced(delete->weightForwardStrand, delete->weightReverseStrand, strandBalanceRatio)) {
				fprintf(fH, "Delete\tLength:%" PRIi64 "\tTotal weight:%f\tForward Strand Weight:%f\tReverse Strand Weight:%f\n",
						delete->length,
						(float)poaDelete_getWeight(delete)/PAIR_ALIGNMENT_PROB_1,
						(float)delete->weightForwardStrand/PAIR_ALIGNMENT_PROB_1,
						(float)delete->weightReverseStrand/PAIR_ALIGNMENT_PROB_1);
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

char *poa_getConsensus(Poa *poa, int64_t **poaToConsensusMap, PolishParams *pp) {
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
			totalIndelWeight += poaInsert_getWeight(insert);
		}
		for(int64_t j=0; j<stList_length(node->deletes); j++) {
			PoaDelete *delete = stList_get(node->deletes, j);
			totalIndelWeight += poaDelete_getWeight(delete);
		}

		// Calculate the match transition weight of node
		// that is, we make an estimate of the weight/expectation
		// of transitioning from this node to the next node without an indel

		double matchTransitionWeight = 0.0;
		if(i == 0) {
			if(stList_length(poa->nodes) == 1) { // In case is zero length reference
				matchTransitionWeight = 1.0;
			}
			else {
				// Set the initiation probability according to the average base weight
				for(int64_t j=1; j<stList_length(poa->nodes); j++) {
					PoaNode *nNode = stList_get(poa->nodes, j);
					for(int64_t k=0; k<SYMBOL_NUMBER_NO_N; k++) {
						matchTransitionWeight += nNode->baseWeights[k];
					}
				}
				matchTransitionWeight /= stList_length(poa->nodes)-1;
				matchTransitionWeight -= totalIndelWeight;
			}
		}
		else {
			for(int64_t j=0; j<SYMBOL_NUMBER; j++) {
				matchTransitionWeight += node->baseWeights[j];
			}
			matchTransitionWeight -= totalIndelWeight;
		}

		// Hack to stop zero weights
		matchTransitionWeight = matchTransitionWeight <= 0.0 ? 0.0001 : matchTransitionWeight; // Make a small value

		// Calculate the total weight of outgoing transitions
		totalOutgoingWeights[i] = matchTransitionWeight + totalIndelWeight;

		// Update the probabilities of nodes that connect by to this node

		// Inserts
		for(int64_t j=0; j<stList_length(node->inserts); j++) {
			PoaInsert *insert = stList_get(node->inserts, j);
			nodeForwardLogProbs[i+1] = logAdd(nodeForwardLogProbs[i+1],
					nodeForwardLogProbs[i] + log(poaInsert_getWeight(insert)/totalOutgoingWeights[i]));
		}

		// Deletes
		for(int64_t j=0; j<stList_length(node->deletes); j++) {
			PoaDelete *delete = stList_get(node->deletes, j);
			nodeForwardLogProbs[i+delete->length+1] = logAdd(nodeForwardLogProbs[i+delete->length+1],
					nodeForwardLogProbs[i] + log(poaDelete_getWeight(delete)/totalOutgoingWeights[i]));
		}

		// Match
		matchTransitionForwardLogProbs[i] = nodeForwardLogProbs[i] + log(matchTransitionWeight/totalOutgoingWeights[i]);
		nodeForwardLogProbs[i+1] = logAdd(nodeForwardLogProbs[i+1], matchTransitionForwardLogProbs[i]);
	}

	// Now traceback picking consensus greedily

	// Allocate consensus map, setting the alignment of reference
	// string positions initially all to gaps.
	*poaToConsensusMap = st_malloc((stList_length(poa->nodes)-1) * sizeof(int64_t));
	for(int64_t i=0; i<stList_length(poa->nodes)-1; i++) {
		(*poaToConsensusMap)[i] = -1;
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

			if(refBaseWeight * pp->referenceBasePenalty > maxBaseWeight) {
				maxBaseIndex = refBaseIndex;
			}

			stList_append(consensusStrings, stString_print("%c", symbol_convertSymbolToChar(maxBaseIndex)));

			// Update poa to consensus map
			(*poaToConsensusMap)[i-1] = runningConsensusLength++;
		}

		// Get max insert
		double maxInsertProb = LOG_ZERO;
		double totalInsertProb = LOG_ZERO;
		PoaInsert *maxInsert = NULL;
		PoaNode *pNode = stList_get(poa->nodes, i-1);
		for(int64_t j=0; j<stList_length(pNode->inserts); j++) {
			PoaInsert *insert = stList_get(pNode->inserts, j);
			double p = log(poaInsert_getWeight(insert)/totalOutgoingWeights[i-1]) + nodeForwardLogProbs[i-1];
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
			double p = log(poaDelete_getWeight(delete)/totalOutgoingWeights[i-delete->length-1]) + nodeForwardLogProbs[i-delete->length-1];
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
		if((*poaToConsensusMap)[i] != -1) {
			(*poaToConsensusMap)[i] = runningConsensusLength - 1 - (*poaToConsensusMap)[i];
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

/*
Poa *poa_realignIterative(stList *reads, stList *alignments, char *reference, PolishParams *polishParams) {
	Poa *poa = poa_realign(reads, alignments, reference, polishParams);

	time_t startTime = time(NULL);

	double score = poa_getReferenceNodeTotalMatchWeight(poa) - poa_getTotalErrorWeight(poa);

	st_logInfo("Starting realignment with score: %f\n", score/PAIR_ALIGNMENT_PROB_1);

	int64_t i=0;
	while(1) {
		int64_t *poaToConsensusMap;
		reference = poa_getConsensus(poa, &poaToConsensusMap, polishParams);

		// Stop in case consensus string is same as old reference (i.e. greedy convergence)
		if(stString_eq(reference, poa->refString)) {
			free(reference);
			free(poaToConsensusMap);
			break;
		}


		// Get anchor alignments
		stList *anchorAlignments = poa_getAnchorAlignments(poa, poaToConsensusMap, stList_length(reads), polishParams);

		time_t realignStartTime = time(NULL);

		// Generated updated poa
		Poa *poa2 = poa_realign(reads, anchorAlignments, reference, polishParams);

		// Cleanup
		free(reference);
		free(poaToConsensusMap);
		stList_destruct(anchorAlignments);

		double score2 = poa_getReferenceNodeTotalMatchWeight(poa2) - poa_getTotalErrorWeight(poa2);

		st_logInfo("Took %f seconds to do round %" PRIi64 " of realignment, Have score: %f (%f score diff)\n",
					(float)(time(NULL) - realignStartTime), i+1, score2/PAIR_ALIGNMENT_PROB_1, (score2-score)/PAIR_ALIGNMENT_PROB_1);

		// Stop if score decreases (greedy stopping)
		if(score2 <= score) {
			poa_destruct(poa2);
			break;
		}

		poa_destruct(poa);
		poa = poa2;
		score = score2;

		//if(++i >= 2) {
		//	break;
		//}
		i++;
	}

	st_logInfo("Took %f seconds to realign iterative through %" PRIi64 " iterations, got final score : %f\n",
			(float)(time(NULL) - startTime), i+1, score/PAIR_ALIGNMENT_PROB_1);

	return poa;
}
*/

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

Poa *poa_checkMajorIndelEditsGreedily(Poa *poa, stList *bamChunkReads, PolishParams *polishParams) {
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
				if(maxInsert == NULL || poaInsert_getWeight(insert) > poaInsert_getWeight(maxInsert)) {
					maxInsert = insert;
					insertStart = i;
				}
			}

			// Get max delete
			for(int64_t j=0; j<stList_length(node->deletes); j++) {
				PoaDelete *delete = stList_get(node->deletes, j);
				if(maxDelete == NULL || poaDelete_getWeight(delete) > poaDelete_getWeight(maxDelete)) {
					maxDelete = delete;
					deleteStart = i;
				}
			}
		}

		if(maxInsert == NULL || maxDelete == NULL) {
			return poa;
		}

		// Create new graph with edit
		char *editRef = poaInsert_getWeight(maxInsert) >= poaDelete_getWeight(maxDelete) ? addInsert(poa->refString, maxInsert->insert, insertStart) :
				removeDelete(poa->refString, maxDelete->length, deleteStart);
		// TODO: Add anchor constraints
		Poa *poa2 = poa_realign(bamChunkReads, NULL, editRef, polishParams);
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

stList *poa_getReadAlignmentsToConsensus(Poa *poa, stList *bamChunkReads, PolishParams *polishParams) {
	// Generate anchor alignments
	stList *anchorAlignments = poa_getAnchorAlignments(poa, NULL, stList_length(bamChunkReads), polishParams);

	// Alignments
	stList *alignments = stList_construct3(0, (void (*)(void *))stList_destruct);

	// Make the MEA alignments
	int64_t refLength = stList_length(poa->nodes)-1;
	for(int64_t i=0; i<stList_length(bamChunkReads); i++) {
		BamChunkRead* read = stList_get(bamChunkReads, i);
		char *nucleotides  = read->nucleotides;
		stList *anchorAlignment = stList_get(anchorAlignments, i);

		// Generate the posterior alignment probabilities
		stList *matches, *inserts, *deletes;
		getAlignedPairsWithIndelsCroppingReference(poa->refString, stList_length(poa->nodes)-1, nucleotides,
				anchorAlignment, &matches, &inserts, &deletes, polishParams);

		// Get the MEA alignment
		double alignmentScore;
		stList *alignment = getMaximalExpectedAccuracyPairwiseAlignment(matches, deletes, inserts,
				refLength, strlen(nucleotides), &alignmentScore, polishParams->p);

		// Left shift the alignment
		stList *leftShiftedAlignment = leftShiftAlignment(alignment, poa->refString, nucleotides);

		// Cleanup
		stList_destruct(inserts);
		stList_destruct(deletes);
		stList_destruct(matches);
		stList_destruct(alignment);

		stList_append(alignments, leftShiftedAlignment);
	}

	// Cleanup
	stList_destruct(anchorAlignments);

	return alignments;
}

/*
 * Functions for run-length encoding/decoding with POAs
 */

RleString *rleString_construct(char *str) {
	RleString *rleString = st_calloc(1, sizeof(RleString));

	rleString->nonRleLength = strlen(str);

	// Calc length of rle'd str
	for(int64_t i=0; i<rleString->nonRleLength; i++) {
		if(i+1 == rleString->nonRleLength || str[i] != str[i+1]) {
			rleString->length++;
		}
	}

	// Allocate
	rleString->rleString = st_calloc(rleString->length+1, sizeof(char));
	rleString->repeatCounts = st_calloc(rleString->length, sizeof(int64_t));
	rleString->rleToNonRleCoordinateMap = st_calloc(rleString->length, sizeof(int64_t));
	rleString->nonRleToRleCoordinateMap = st_calloc(rleString->nonRleLength, sizeof(int64_t));

	// Fill out
	int64_t j=0, k=1;
	for(int64_t i=0; i<rleString->nonRleLength; i++) {
		rleString->nonRleToRleCoordinateMap[i] = j;
		if(i+1 == rleString->nonRleLength || str[i] != str[i+1]) {
			rleString->rleString[j] = str[i];
			rleString->repeatCounts[j] = k;
			rleString->rleToNonRleCoordinateMap[j++] = i - k + 1;
			k=1;
		}
		else {
			k++;
		}
	}
	assert(j == rleString->length);

	return rleString;
}


RleString *rleString_construct2(char *rleChars, uint8_t *rleCounts) {
	RleString *rleString = st_calloc(1, sizeof(RleString));

	rleString->length = strlen(rleChars);
	rleString->nonRleLength = 0;
	for (int64_t i = 0; i < rleString->length; i++) {
		rleString->nonRleLength += rleCounts[i];
	}

	// Allocate
	rleString->rleString = st_calloc(rleString->length+1, sizeof(char));
	rleString->repeatCounts = st_calloc(rleString->length, sizeof(int64_t));
	rleString->rleToNonRleCoordinateMap = st_calloc(rleString->length, sizeof(int64_t));
	rleString->nonRleToRleCoordinateMap = st_calloc(rleString->nonRleLength, sizeof(int64_t));

	// Fill out
    //TODO verify this logic
	int64_t n=1;
	for(int64_t r=0; r<rleString->length; r++) {
		// counts
		rleString->repeatCounts[r] = rleCounts[r];
		// coordinates
		rleString->rleToNonRleCoordinateMap[r] = n;
		for (uint8_t c = 0; c < rleCounts[r]; c++) {
			rleString->nonRleToRleCoordinateMap[n + c] = r;
		}
	}
	assert(n == rleString->nonRleLength);

	return rleString;
}

void rleString_destruct(RleString *rleString) {
	free(rleString->rleString);
	free(rleString->repeatCounts);
	free(rleString->rleToNonRleCoordinateMap);
	free(rleString->nonRleToRleCoordinateMap);
	free(rleString);
}

char *rleString_expand(RleString *rleString) {
	char *s = st_calloc(rleString->nonRleLength+1, sizeof(char));
	int64_t j=0;
	for(int64_t i=0; i<rleString->length; i++) {
		for(int64_t k=0; k<rleString->repeatCounts[i]; k++) {
			s[j++] = rleString->rleString[i];
		}
	}
	s[rleString->nonRleLength] = '\0';
	return s;
}

static int64_t expandRLEConsensus2(PoaNode *node, stList *rleReads, stList *bamChunkReads, RepeatSubMatrix *repeatSubMatrix) {
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
	return repeatSubMatrix_getMLRepeatCount(repeatSubMatrix, maxBaseIndex, node->observations,
			rleReads, bamChunkReads, &logProbability);
}

RleString *expandRLEConsensus(Poa *poa, stList *rleReads, stList *bamChunkReads, RepeatSubMatrix *repeatSubMatrix) {
	RleString *rleString = st_calloc(1, sizeof(RleString));

	rleString->length = stList_length(poa->nodes)-1;
	rleString->rleString = stString_copy(poa->refString);
	rleString->repeatCounts = st_calloc(rleString->length, sizeof(int64_t));
	rleString->rleToNonRleCoordinateMap = st_calloc(rleString->length, sizeof(int64_t));
	for(int64_t i=1; i<stList_length(poa->nodes); i++) {
		int64_t repeatCount = expandRLEConsensus2(stList_get(poa->nodes, i), rleReads, bamChunkReads, repeatSubMatrix);
		rleString->repeatCounts[i-1] = repeatCount;
		rleString->rleToNonRleCoordinateMap[i-1] = rleString->nonRleLength;
		rleString->nonRleLength += repeatCount;
	}
	rleString->nonRleToRleCoordinateMap = st_calloc(rleString->nonRleLength, sizeof(int64_t));
	int64_t j=0;
	for(int64_t i=0; i<rleString->length; i++) {
		for(int64_t k=0; k<rleString->repeatCounts[i]; k++) {
			rleString->nonRleToRleCoordinateMap[j++] = i;
		}
	}
	assert(rleString->nonRleLength == j);

	return rleString;
}

stList *runLengthEncodeAlignment(stList *alignment,
							     RleString *seqX, RleString *seqY) {
	stList *rleAlignment = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);

	int64_t x=-1, y=-1;
	for(int64_t i=0; i<stList_length(alignment); i++) {
		stIntTuple *alignedPair = stList_get(alignment, i);

		int64_t x2 = seqX->nonRleToRleCoordinateMap[stIntTuple_get(alignedPair, 0)];
		int64_t y2 = seqY->nonRleToRleCoordinateMap[stIntTuple_get(alignedPair, 1)];

		if(x2 > x && y2 > y) {
			stList_append(rleAlignment, stIntTuple_construct3(x2, y2, stIntTuple_get(alignedPair, 2)));
			x = x2; y = y2;
		}
	}

	return rleAlignment;
}

/*
 * Functions for modeling repeat counts
 */

double *repeatSubMatrix_setLogProb(RepeatSubMatrix *repeatSubMatrix, Symbol base, bool strand, int64_t observedRepeatCount, int64_t underlyingRepeatCount) {
    if (base == SYMBOL_NUMBER_NO_N) {
        st_errAbort("[repeatSubMatrix_setLogProb] base 'Nn' not supported for repeat estimation\n");
    }
    int64_t idx = (2 * base + (strand ? 1 : 0)) * repeatSubMatrix->maximumRepeatLength * repeatSubMatrix->maximumRepeatLength +
            underlyingRepeatCount * repeatSubMatrix->maximumRepeatLength +
            observedRepeatCount;
    assert(idx < repeatSubMatrix->maxEntry);
	return &(repeatSubMatrix->logProbabilities[idx]);
}

double repeatSubMatrix_getLogProb(RepeatSubMatrix *repeatSubMatrix, Symbol base, bool strand, int64_t observedRepeatCount, int64_t underlyingRepeatCount) {
	double *loc = repeatSubMatrix_setLogProb(repeatSubMatrix, base, strand, observedRepeatCount, underlyingRepeatCount);
//	printf("%p\n", loc);
	return *loc;
}

void repeatSubMatrix_destruct(RepeatSubMatrix *repeatSubMatrix) {
	free(repeatSubMatrix->logProbabilities);
	free(repeatSubMatrix);
}

RepeatSubMatrix *repeatSubMatrix_constructEmpty() {
	RepeatSubMatrix *repeatSubMatrix = st_calloc(1, sizeof(RepeatSubMatrix));
	repeatSubMatrix->maximumRepeatLength = 51;
	repeatSubMatrix->maxEntry = 2 * SYMBOL_NUMBER_NO_N * repeatSubMatrix->maximumRepeatLength * repeatSubMatrix->maximumRepeatLength;
	repeatSubMatrix->logProbabilities = st_calloc(repeatSubMatrix->maxEntry, sizeof(double));
	return repeatSubMatrix;
}

double repeatSubMatrix_getLogProbForGivenRepeatCount(RepeatSubMatrix *repeatSubMatrix, Symbol base, stList *observations,
												     stList *rleReads, stList *bamChunkReads, int64_t underlyingRepeatCount) {
	assert(underlyingRepeatCount < repeatSubMatrix->maximumRepeatLength);
	double logProb = LOG_ONE;
	for(int64_t i=0; i<stList_length(observations); i++) {
		PoaBaseObservation *observation = stList_get(observations, i);
		BamChunkRead *read = stList_get(bamChunkReads, observation->readNo);
		RleString *rleRead = stList_get(rleReads, observation->readNo);
		int64_t observedRepeatCount = rleRead->repeatCounts[observation->offset];

		// Be robust to over-long repeat count observations
		observedRepeatCount = observedRepeatCount >= repeatSubMatrix->maximumRepeatLength ?
		        repeatSubMatrix->maximumRepeatLength-1 : observedRepeatCount;

		logProb += repeatSubMatrix_getLogProb(repeatSubMatrix, base, read->forwardStrand,
		        observedRepeatCount, underlyingRepeatCount) * observation->weight;
	}

	return logProb;
}

int64_t repeatSubMatrix_getMLRepeatCount(RepeatSubMatrix *repeatSubMatrix, Symbol base, stList *observations,
		stList *rleReads, stList *bamChunkReads, double *logProbability) {
	if(stList_length(observations) == 0) {
		return 0; // The case that we have no alignments, we assume there is no sequence there/
	}

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
	if(maxRepeatLength >= repeatSubMatrix->maximumRepeatLength) {
		st_logCritical("Got overlong repeat observation: %" PRIi64 ", ignoring this and cutting off overlong repeat counts to max\n", maxRepeatLength);
		maxRepeatLength = repeatSubMatrix->maximumRepeatLength-1;
	}

	// Calc the range of repeat observations
	double mlLogProb = repeatSubMatrix_getLogProbForGivenRepeatCount(repeatSubMatrix, base, observations, rleReads,
	        bamChunkReads, minRepeatLength);
	int64_t mlRepeatLength = minRepeatLength;
	for(int64_t i=minRepeatLength+1; i<maxRepeatLength+1; i++) {
		double p = repeatSubMatrix_getLogProbForGivenRepeatCount(repeatSubMatrix, base, observations, rleReads,
		        bamChunkReads, i);
		if(p > mlLogProb) {
			mlLogProb = p;
			mlRepeatLength = i;
		}
	}
	*logProbability = mlLogProb;
	return mlRepeatLength;
}

int64_t removeOverlap(char *prefixString, char *suffixString, int64_t approxOverlap, PolishParams *polishParams,
				   int64_t *prefixStringCropEnd, int64_t *suffixStringCropStart) {

	// Align the overlapping suffix of the prefixString and the prefix of the suffix string
	int64_t prefixStringLength = strlen(prefixString);
	int64_t suffixStringLength = strlen(suffixString);

	// Get coordinates of substrings to be aligned
	int64_t i = (prefixStringLength - approxOverlap) < 0 ? 0 : prefixStringLength - approxOverlap;
	int64_t j = approxOverlap < suffixStringLength ? approxOverlap : suffixStringLength;

	// Crop suffix
	char c = suffixString[j];
	suffixString[j] = '\0';

	// Run the alignment
	stList *alignedPairs = getAlignedPairs(polishParams->sM, &(prefixString[i]), suffixString, polishParams->p, 1, 1);

	if(stList_length(alignedPairs) == 0 && st_getLogLevel() >= info) {
		st_logInfo("Failed to find good overlap. Suffix-string: %s\n", &(prefixString[i]));
		st_logInfo("Failed to find good overlap. Prefix-string: %s\n", suffixString);
	}

	// Remove the suffix crop
	suffixString[j] = c;

	// Pick the median point
	stIntTuple *maxPair = NULL;
	for(int64_t k=0; k<stList_length(alignedPairs); k++) {
		stIntTuple *aPair = stList_get(alignedPairs, k);
		if(maxPair == NULL || stIntTuple_get(aPair, 0) > stIntTuple_get(maxPair, 0)) {
			maxPair = aPair;
		}
	}
	if(maxPair == NULL) {
		st_logCritical("Failed to find any aligned pairs between overlapping strings, not "
				"doing any trimming (approx overlap: %i, len x: %i, len y: %i)\n", approxOverlap, prefixStringLength, suffixStringLength);
		*prefixStringCropEnd = prefixStringLength;
		*suffixStringCropStart = 0;
	}
	else {
		*prefixStringCropEnd = stIntTuple_get(maxPair, 1) + i; // Exclusive
		*suffixStringCropStart = stIntTuple_get(maxPair, 2);  // Inclusive
	}

	int64_t overlapWeight = maxPair == NULL ? 0 : stIntTuple_get(maxPair, 0);

	stList_destruct(alignedPairs);

	return overlapWeight;
}

// New polish algorithm

char *getExistingSubstring(Poa *poa, int64_t from, int64_t to) {
	/*
	 * Gets substring of the poa reference string from "from" (inclusive)
	 * to "to" (exclusive).
	 */
	char *s = st_malloc(sizeof(char) * (to-from+1));
	s[to-from] = '\0';
	for(int64_t i=from; i<to; i++) {
		PoaNode *node = stList_get(poa->nodes, i);
		s[i-from] = node->base;
	}
	return s;
}

double getTotalWeight(PoaNode *node) {
	/*
	 * Returns the total base weight of reads aligned to the given node.
	 */
	double totalWeight = 0.0;
	for(int64_t i=0; i<SYMBOL_NUMBER; i++) {
		totalWeight += node->baseWeights[i];
	}
	return totalWeight;
}

double getAvgCoverage(Poa *poa, int64_t from, int64_t to) {
	// Calculate average coverage, which is used to determine candidate variants
	double avgCoverage = 0.0;
	for(int64_t j=from; j<to; j++) {
		avgCoverage += getTotalWeight(stList_get(poa->nodes, j));
	}
	return avgCoverage / (to-from);
}

char getNextCandidateBase(PoaNode *node, int64_t *i, double candidateWeight) {
	/*
	 * Iterates through candidate bases for a reference position returning those with sufficient weight.
	 * Always returns the reference base
	 */
	double totalWeight = getTotalWeight(node);
	while(*i<SYMBOL_NUMBER) {
		char base = symbol_convertSymbolToChar(*i);
		if(node->baseWeights[(*i)++] > candidateWeight || toupper(node->base) == base) {
			return base;
		}
	}
	return '-';
}

bool hasCandidateSubstitution(PoaNode *node, double candidateWeight) {
	/*
	 * Returns non-zero if the node has a candidate base that is different to the
	 * current base.
	 */
	int64_t i=0;
	char base;
	while((base = getNextCandidateBase(node, &i, candidateWeight)) != '-') {
		if(base != node->base) {
			return 1;
		}
	}
	return 0;
}

char *getNextCandidateInsert(PoaNode *node, int64_t *i, double candidateWeight) {
	/*
	 * Iterates through candidate insertions for a reference position returning those with sufficient weight.
	 */
	while((*i)++ < stList_length(node->inserts)) {
		PoaInsert *insert = stList_get(node->inserts, (*i)-1);
		if(poaInsert_getWeight(insert) > candidateWeight) {
			return insert->insert;
		}
	}
	return NULL;
}

bool hasCandidateInsert(PoaNode *node, double candidateWeight) {
	/*
	 * Returns non-zero if the node has a candidate insert.
	 */
	int64_t i=0;
	return getNextCandidateInsert(node, &i, candidateWeight) != NULL;
}

int64_t getNextCandidateDelete(PoaNode *node, int64_t *i, double candidateWeight) {
	/*
	 * Iterates through candidate deletions for a reference position returning those with sufficient weight.
	 */
	while((*i)++ < stList_length(node->deletes)) {
		PoaDelete *delete = stList_get(node->deletes, (*i)-1);
		if(poaDelete_getWeight(delete) > candidateWeight) {
			return delete->length;
		}
	}
	return -1;
}

bool maxCandidateDeleteLength(PoaNode *node, double candidateWeight) {
	/*
	 * Returns maximum length of a candidate deletion starting after this position.
	 */
	int64_t i=0;
	int64_t deleteLength, maxDeleteLength = 0;
	while((deleteLength = getNextCandidateDelete(node, &i, candidateWeight)) != -1) {
		if(deleteLength > maxDeleteLength) {
			maxDeleteLength = deleteLength;
		}
	}
	return maxDeleteLength;
}

stList *getCandidateConsensusSubstrings(Poa *poa, int64_t from, int64_t to,
										double *candidateWeights, double weightAdjustment, int64_t maximumStringNumber) {
	/*
	 *  A candidate variant is an edit (either insert, delete or substitution) to the poa reference
	 *  string with "high" weight. This function returns all possible combinations of candidate variants,
	 *  each as a new consensus substring, for the interval of the reference string from "from" (inclusive)
	 *  to "to" (exclusive). Returned list of strings always contains the reference string without edits (the no
	 *  candidate variants string).
	 */

	// Function is recursive

	// First get suffix substrings
	stList *suffixes;
	if(from+1 < to) {
		suffixes = getCandidateConsensusSubstrings(poa, from+1, to, candidateWeights, weightAdjustment, maximumStringNumber);

		if(suffixes == NULL) { // If too many combinations, return null.
			return NULL;
		}
	}
	else {
		suffixes = stList_construct3(0, free);
		stList_append(suffixes, stString_copy("")); // Start with the empty string
	}

	// Now extend by adding on prefix variants
	stList *consensusSubstrings = stList_construct3(0, free);

	PoaNode *node = stList_get(poa->nodes, from);

	double candidateWeight = candidateWeights[from] * weightAdjustment;

	int64_t i=0;
	char base;
	while((base = getNextCandidateBase(node, &i, candidateWeight)) != '-') { // Enumerate the possible bases at the reference node.

		// Create the consensus substrings with no inserts or deletes starting at this node
		for(int64_t j=0; j<stList_length(suffixes); j++) {
			stList_append(consensusSubstrings, stString_print("%c%s", base, stList_get(suffixes, j)));
		}

		// Now add insert cases
		int64_t k=0;
		char *insert;
		while((insert = getNextCandidateInsert(node, &k, candidateWeight)) != NULL) {
			for(int64_t j=0; j<stList_length(suffixes); j++) {
				stList_append(consensusSubstrings, stString_print("%c%s%s", base, insert, stList_get(suffixes, j)));
			}
		}

		// Add then deletes
		k = 0;
		int64_t deleteLength;
		while((deleteLength = getNextCandidateDelete(node, &k, candidateWeight)) > 0) {
			for(int64_t j=0; j<stList_length(suffixes); j++) {
				char *suffixHaplotype = stList_get(suffixes, j);
				stList_append(consensusSubstrings, stString_print("%c%s", base,
							((int64_t)strlen(suffixHaplotype) - deleteLength >= 0) ? &(suffixHaplotype[deleteLength]) : ""));
			}
		}
	}

	// Clean up
	stList_destruct(suffixes);

	if(stList_length(consensusSubstrings) > maximumStringNumber) {
		// Clean up and return null (too many combinations)
		stList_destruct(consensusSubstrings);
		return NULL;
	}

	return consensusSubstrings;
}

double computeLogLikelihoodOfConsensusString(char *reference, stList *nucleotides, PolishParams *params) {
	/*
	 * Computes the log probability of the reference given the reads.
	 */
	double logProb = LOG_ONE;
	stList *anchorPairs = stList_construct(); // Currently empty
	for(int64_t i=0; i<stList_length(nucleotides); i++) {
		char *read = stList_get(nucleotides, i);
		logProb += computeForwardProbability(reference, read, anchorPairs, params->p, params->sM, 0, 0);
	}

	// Cleanup
	stList_destruct(anchorPairs);

	return logProb;
}

int poaBaseObservation_cmp(const void *a, const void *b) {
	PoaBaseObservation *obs1 = (PoaBaseObservation *)a;
	PoaBaseObservation *obs2 = (PoaBaseObservation *)b;
	if(obs1->readNo != obs2->readNo) { // Sort first is ascending read number order
		return obs1->readNo < obs2->readNo ? -1 : 1;
	}
	if(obs1->weight != obs2->weight) { // Sort second in descending weight order
		return obs1->weight > obs2->weight ? -1 : 1;
	}
	return 0;
}

void sortBaseObservations(Poa *poa) {
	/*
	 * Sort the POA base observations to make them appropriate for getReadSubstrings.
	 */
	for(int64_t i=0; i<stList_length(poa->nodes); i++) {
		PoaNode *node = stList_get(poa->nodes, i);
		stList_sort(node->observations, poaBaseObservation_cmp);
	}
}

int64_t skipDupes(PoaNode *node, int64_t i, int64_t readNo) {
	while(i < stList_length(node->observations)) {
		PoaBaseObservation *obs = stList_get(node->observations, i);
		if(obs->readNo != readNo) {
			break;
		}
		i++;
	}
	return i;
}

void getReadSubstring(BamChunkRead *bamChunkRead, int64_t start, int64_t length, stList *readSubstrings, stList *qualityValues) {

	// Calculate the qual value
	double qual = 0;
	if(bamChunkRead->qualities != NULL) {
		int64_t j = 0;
		for(int64_t i=0; i<length; i++) {
			j += (int64_t)bamChunkRead->qualities[i+start];
		}
		qual = (float)j / length;
		qual /= -10.0; // Quals are phred, qual = -10 * log_10(p), convert to
			// log probability
	}

	// Add read substring
	char *read = &(bamChunkRead->nucleotides[start]);
	char c = read[length];
	read[length] = '\0';

	if(qual > 20.0 || 1) { // This is a hack
		stList_append(readSubstrings, stString_copy(read));
		stList_append(qualityValues, stDoubleTuple_construct(1, qual));
	}

	read[length] = c;
}

void getReadSubstrings(stList *bamChunkReads, Poa *poa, int64_t from, int64_t to,
						  stList *readSubstrings, stList *qualityValues) {
	/*
	 * Get the substrings of reads aligned to the interval from (inclusive) to to
	 * (exclusive) and their qual values. Adds them to readSubstrings and qualValues, respectively.
	 */

	// Deal with boundary cases
	if(from == 0) {
		if(to == stList_length(poa->nodes)) {
			// If from and to reference positions that bound the complete alignment just
			// copy the complete reads
			for(int64_t i=0; i<stList_length(bamChunkReads); i++) {
			    BamChunkRead *bamChunkRead = stList_get(bamChunkReads, i);
				getReadSubstring(bamChunkRead, 0, bamChunkRead->readLength, readSubstrings, qualityValues);
			}
			return;
		}

		// Otherwise, include the read prefixes that end at to
		PoaNode *node = stList_get(poa->nodes, to);
		int64_t i=0;
		while(i<stList_length(node->observations)) {
			PoaBaseObservation *obs = stList_get(node->observations, i);
            BamChunkRead *bamChunkRead = stList_get(bamChunkReads, obs->readNo);
            // Trim the read substring, copy it and add to the substrings list
			getReadSubstring(bamChunkRead, 0, obs->offset, readSubstrings, qualityValues);
			i = skipDupes(node, ++i, obs->readNo);
		}
		return;
	}
	else if(to == stList_length(poa->nodes)) {
		// Finally, include the read suffixs that start at from
		PoaNode *node = stList_get(poa->nodes, from);
		int64_t i = 0;
		while (i < stList_length(node->observations)) {
			PoaBaseObservation *obs = stList_get(node->observations, i);
            BamChunkRead *bamChunkRead = stList_get(bamChunkReads, obs->readNo);
			// Trim the read substring, copy it and add to the substrings list
			getReadSubstring(bamChunkRead, obs->offset, bamChunkRead->readLength-obs->offset, readSubstrings, qualityValues);
			i = skipDupes(node, ++i, obs->readNo);
		}
		return;
	}

	PoaNode *fromNode = stList_get(poa->nodes, from);
	PoaNode *toNode = stList_get(poa->nodes, to);

	int64_t i=0, j=0;
	while(i < stList_length(fromNode->observations) && j < stList_length(toNode->observations)) {
		PoaBaseObservation *obsFrom = stList_get(fromNode->observations, i);
		PoaBaseObservation *obsTo = stList_get(toNode->observations, j);

		if(obsFrom->readNo == obsTo->readNo) {
            BamChunkRead *bamChunkRead = stList_get(bamChunkReads, obsFrom->readNo);
            getReadSubstring(bamChunkRead, obsFrom->offset, obsTo->offset-obsFrom->offset, readSubstrings, qualityValues);
			i = skipDupes(fromNode, ++i, obsFrom->readNo);
			j = skipDupes(toNode, ++j, obsTo->readNo);
		}
		else if(obsFrom->readNo < obsTo->readNo) {
			i = skipDupes(fromNode, ++i, obsFrom->readNo);
		}
		else {
			assert(obsFrom->readNo > obsTo->readNo);
			j = skipDupes(toNode, ++j, obsTo->readNo);
		}
	}
}

char *getBestConsensusSubstring(Poa *poa, stList *bamChunkReads, int64_t from, int64_t to, double *candidateWeights, PolishParams *params) {
	/*
	 * Heuristically searches for the best consensus substring between from (inclusive) and to (exclusive).
	 * Does so by enumerating candidate variants in interval and then building and test all possible resulting
	 * consensus substrings
	 */
	if(to-from > 1000) { // If region to be anchored is too long give up
		return getExistingSubstring(poa, from, to);
	}

	// Get consensus substrings, ensuring we use a candidate weight that does not cause us to evaluate more than
	// a maximum number of strings.
	stList *consensusSubstrings = NULL;
	double weightAdjustment = 1.0;
	do {
		consensusSubstrings = getCandidateConsensusSubstrings(poa, from, to, candidateWeights, weightAdjustment, params->maxConsensusStrings);
		weightAdjustment *= 1.5; // Increase the candidate weight by 50%
	} while(consensusSubstrings == NULL);

	char *consensusSubstring = stList_pop(consensusSubstrings);

	if(stList_length(consensusSubstrings) > 0) {
		// Get read substrings
		stList *readSubstrings = stList_construct3(0, free);
		stList *qualityValues = stList_construct3(0, (void (*)(void *))stDoubleTuple_destruct);
		getReadSubstrings(bamChunkReads, poa, from, to, readSubstrings, qualityValues);

		if(st_getLogLevel() >= debug) {
			st_logDebug("Got %" PRIi64 " consensus strings from: %" PRIi64 " to %" PRIi64 " with %" PRIi64 " reads\n",
						stList_length(consensusSubstrings)+1, from, to, stList_length(readSubstrings));
			for(int64_t i=0; i<stList_length(readSubstrings); i++) {
				st_logDebug("\tGot read substring: %s %f\n", stList_get(readSubstrings, i), stDoubleTuple_getPosition(stList_get(qualityValues, i), 0));
			}
		}

		// Assess likelihood of each remaining substring in turn,
		// keeping the most probable
		double maxLogProb = computeLogLikelihoodOfConsensusString(consensusSubstring, readSubstrings, params);

		st_logDebug("\tFor consensus-string %s got log-prob: %f\n", consensusSubstring, maxLogProb);

		while(stList_length(consensusSubstrings) > 0) {
			char *c = stList_pop(consensusSubstrings);
			double logProb = computeLogLikelihoodOfConsensusString(c, readSubstrings, params);
			st_logDebug("\tFor consensus-string %s got log-prob: %f\n", c, logProb);
			if(logProb > maxLogProb) {
				maxLogProb = logProb;
				free(consensusSubstring);
				consensusSubstring = c;
			}
			else {
				free(c);
			}
		}

		// Cleanup
		stList_destruct(readSubstrings);
		stList_destruct(qualityValues);
	}

	// Cleanup
	stList_destruct(consensusSubstrings);

	return consensusSubstring;
}

// Code to create anchors

double *getCandidateWeights(Poa *poa, PolishParams *params) {
	double *candidateWeights = st_calloc(stList_length(poa->nodes), sizeof(double));

	int64_t window = 100; // Size of window to average coverage over

	if(window >= stList_length(poa->nodes)) {
		double candidateWeight = getAvgCoverage(poa, 0, stList_length(poa->nodes)) * params->candidateVariantWeight;
		for(int64_t i=0; i<stList_length(poa->nodes); i++) {
			candidateWeights[i] = candidateWeight;
		}
		return candidateWeights;
	}

	double totalWeight = 0;
	for(int64_t i=0; i<stList_length(poa->nodes); i++) {
		totalWeight += getTotalWeight(stList_get(poa->nodes, i));
		if(i >= window) {
			totalWeight -= getTotalWeight(stList_get(poa->nodes, i-window));
			candidateWeights[i-window/2] = totalWeight/window * params->candidateVariantWeight;
		}
	}

	// Fill in bounding bases
	for(int64_t i=0; i<window/2; i++) {
		candidateWeights[i] = candidateWeights[window/2];
		candidateWeights[stList_length(poa->nodes)-1-i] = candidateWeights[stList_length(poa->nodes)-1-window/2];
	}

	return candidateWeights;
}

bool *getCandidateVariantOverlapPositions(Poa *poa, double *candidateWeights) {
	/*
	 * Return a boolean for each poaNode (as an array) indicating if the node is a candidate variant
	 * site or is included in a candidate deletion.
	 */

	bool *candidateVariantPositions = st_calloc(stList_length(poa->nodes), sizeof(bool));

	// Calculate positions that overlap candidate variants
	for(int64_t i=0; i<stList_length(poa->nodes); i++) {
		PoaNode *node = stList_get(poa->nodes, i);

		// Mark as variant if has a candidate substitution or an insert starts at this position
		if(hasCandidateSubstitution(node, candidateWeights[i]) || hasCandidateInsert(node, candidateWeights[i])) {
			candidateVariantPositions[i] = 1;
		}

		int64_t j = maxCandidateDeleteLength(node, candidateWeights[i]);
		if(j > 0) { // Mark as variant if precedes the start of a deletion
			candidateVariantPositions[i] = 1;
		}
		// Mark as variant if is included in candidate deletion
		while(j > 0) {
			assert(i+j < stList_length(poa->nodes));
			candidateVariantPositions[i+(j--)] = 1;
		}
	}

	return candidateVariantPositions;
}

bool *expand(bool *b, int64_t length, int64_t expansion) {
	/*
	 * Returns a bool array in which a position is non-zero if a position
	 * in b +/i expansion is non-zero.
	 */
	bool *b2 = st_calloc(length, sizeof(bool));
	for(int64_t i=0; i<length; i++) {
		if(b[i]) {
			for(int64_t j=i-expansion; j<i+expansion; j++) {
				if(j >= 0 && j < length) {
					b2[j] = 1;
				}
			}
		}
	}

	return b2;
}

bool *getFilteredAnchorPositions(Poa *poa, double *candidateWeights, PolishParams *params) {
	/*
	 * Create set of anchor positions, using positions not close to candidate variants
	 */
	// Identity sites that overlap candidate variants, and expand to surrounding positions
	bool *candidateVariantPositions = getCandidateVariantOverlapPositions(poa, candidateWeights);
	bool *expandedCandidateVariantPositions = expand(candidateVariantPositions, stList_length(poa->nodes), params->columnAnchorTrim);

	// Anchors are those that are not close to expanded candidate variant positions
	bool *anchors = st_calloc(stList_length(poa->nodes), sizeof(bool));
	for(int64_t i=0; i<stList_length(poa->nodes); i++) {
		anchors[i] = !expandedCandidateVariantPositions[i];
	}

	// Cleanup
	free(candidateVariantPositions);
	free(expandedCandidateVariantPositions);

	// Log some stuff about the anchors
	if(st_getLogLevel() >= debug) {
		int64_t totalAnchorNo = 0;
		for(int64_t i=0; i<stList_length(poa->nodes); i++) {
			totalAnchorNo += anchors[i] ? 1 : 0;
		}
		st_logDebug("In poa_polish got %" PRIi64 " anchors for ref seq of length %" PRIi64 ", that's one every: %f bases\n",
					totalAnchorNo, stList_length(poa->nodes), stList_length(poa->nodes)/(double)totalAnchorNo);
	}

	return anchors;
}

// Core polishing logic functions

void poa_polishP(Poa *poa, stList *bamChunkReads, int64_t from, int64_t to,
					stList *consensusSubstrings, int64_t *consensusStringLength,
					int64_t *poaToConsensusMap, double *candidateWeights, PolishParams *params) {
	// Get existing reference string
	char *existingConsensusSubstring = getExistingSubstring(poa, from, to);

	// Get best consensus substring
	char *consensusSubstring = getBestConsensusSubstring(poa, bamChunkReads, from, to, candidateWeights, params);

	// Add new best consensus substring to the growing new consensus string
	stList_append(consensusSubstrings, consensusSubstring);

	// Now update the alignment between the existing reference substring and the new consensus sequences

	if(stString_eq(existingConsensusSubstring, consensusSubstring)) {
		// If the new and old strings are the same then copy the alignment across
		for(int64_t i=0; i<to-from; i++) {
			if(from + i > 0 && *consensusStringLength + i > 0) { // The > 0 checks are to avoid including alignments to the "N" prefix
				poaToConsensusMap[from + i - 1] = *consensusStringLength + i - 1;
			}
		}
	}
	else {
		st_logDebug("In poa polish got anchors, from: %" PRIi64 " to: %" PRIi64
					", \nexisting string:\t%s\nnew string:\t\t%s\n", from, to,
					existingConsensusSubstring, consensusSubstring);

		// Create alignment between new and old consensus strings and update the map
		double alignmentScore;
		stList *l = stList_construct(); // Empty set of alignment anchors
		stList *alignedPairs = getShiftedMEAAlignment(existingConsensusSubstring, consensusSubstring, l, params->p, params->sM,
												   0, 0, &alignmentScore);
		stList_destruct(l);

		for(int64_t k=0; k<stList_length(alignedPairs); k++) {
			stIntTuple *alignedPair = stList_get(alignedPairs, k);
			// Only take high confidence aligned pairs in updated map
			if(((double)stIntTuple_get(alignedPair, 0))/PAIR_ALIGNMENT_PROB_1 > 0.99
					&& from + stIntTuple_get(alignedPair, 1) > 0 && *consensusStringLength + stIntTuple_get(alignedPair, 2) > 0) {
				// The > 0 checks are to avoid including alignments to the "N" prefix
				poaToConsensusMap[from + stIntTuple_get(alignedPair, 1) - 1] = *consensusStringLength + stIntTuple_get(alignedPair, 2) - 1;
			}
		}

		stList_destruct(alignedPairs);
	}

	// Cleanup
	free(existingConsensusSubstring);


	// Updare consensus string length
	*consensusStringLength += strlen(consensusSubstring);
}

char *poa_polish2(Poa *poa, stList *bamChunkReads, PolishParams *params,
				  int64_t **poaToConsensusMap) {
	/*
	 * As pao_polish, but returns the polished reference string and a
	 * map back to the input poa reference sequence.
	 */

	// Setup
	double *candidateWeights = getCandidateWeights(poa, params);

	if(st_getLogLevel() >= info) {
		double avgCoverage = getAvgCoverage(poa, 0, stList_length(poa->nodes));
		double totalCandidateWeight = 0.0;
		for(int64_t i=0; i<stList_length(poa->nodes); i++) {
			totalCandidateWeight += candidateWeights[i];
		}
		st_logDebug("Got avg. coverage: %f for region of length: %" PRIi64 " and avg. candidate weight of: %f\n",
				avgCoverage/PAIR_ALIGNMENT_PROB_1, stList_length(poa->nodes), totalCandidateWeight/(PAIR_ALIGNMENT_PROB_1*stList_length(poa->nodes)));
	}

	// Sort the base observations to make the getReadSubstrings function work
	sortBaseObservations(poa);

	// Identify anchor points, represented as a binary array, one bit for each POA node
	bool *anchors = getFilteredAnchorPositions(poa, candidateWeights, params);

	// Map to track alignment between the new consensus sequence and the poa reference sequence
	*poaToConsensusMap = st_malloc((stList_length(poa->nodes)-1) * sizeof(int64_t));
	for(int64_t i=0; i<stList_length(poa->nodes)-1; i++) {
		(*poaToConsensusMap)[i] = -1;
	}

	// Substrings of the consensus string that when concatenated form the overall consensus string
	stList *consensusSubstrings = stList_construct3(0, free);
	int64_t j=0; // Length of the growing consensus substring

	// Enumerate candidate variant combinations between anchors

	int64_t pAnchor = 0; // Previous anchor, starting from first position of POA, which is the prefix "N"
	for(int64_t i=1; i<stList_length(poa->nodes); i++) {
		if(anchors[i]) { // If position i is an anchor
			// In case anchors are trivially adjacent
			if(i-pAnchor == 1)  {
				stList_append(consensusSubstrings, stString_print("%c", ((PoaNode *)stList_get(poa->nodes, pAnchor))->base));
				if(i > 0 && j > 0) {
					(*poaToConsensusMap)[i-1] = j-1;
				}
				j++;
			}
			else {
				poa_polishP(poa, bamChunkReads, pAnchor, i,
							consensusSubstrings, &j, *poaToConsensusMap, candidateWeights, params);
			}

			// Update previous anchor
			pAnchor = i;
		}
	}

	// Deal with the suffix
	poa_polishP(poa, bamChunkReads, pAnchor, stList_length(poa->nodes),
				consensusSubstrings, &j, *poaToConsensusMap, candidateWeights, params);

	// Build the new consensus string by concatenating the constituent pieces
	char *newConsensusString = stString_join2("", consensusSubstrings);

	// Remove the prefix "N" character (which is the first position of any POA)
	assert(strlen(newConsensusString) >= 1);
	assert(newConsensusString[0] == 'N');
	char *c = newConsensusString;
	newConsensusString = stString_copy(&(newConsensusString[1]));
	free(c);

	// Cleanup
	free(candidateWeights);
	free(anchors);
	stList_destruct(consensusSubstrings);

	return newConsensusString;
}

Poa *poa_polish(Poa *poa, stList *bamChunkReads, PolishParams *params) {
	/*
	 * "Polishes" the given POA reference string to create a new consensus reference string.
	 * Algorithm starts by dividing the reference into anchor points - points where the majority
	 * of bamChunkReads are confidently aligned. Then, between each pair of consecutive anchors it looks
	 * for "candidate variants", variants (either substitutions, insertions or deletions) with
	 * significant weight. For every combination of candidate variants between the two anchor points
	 * a candidate string is constructed and all the read substrings mapping between the two anchor points
	 * are aligned to it. The candidate string, including the current reference substring,
	 * with the highest likelihood is then selected.
	 */

	time_t startTime = time(NULL);

	int64_t *poaToConsensusMap;
	char *newConsensusString = poa_polish2(poa, bamChunkReads, params, &poaToConsensusMap);

	// Get anchor alignments
	stList *anchorAlignments = poa_getAnchorAlignments(poa, poaToConsensusMap, stList_length(bamChunkReads), params);

	// Generated updated poa
	Poa *poa2 = poa_realign(bamChunkReads, anchorAlignments, newConsensusString, params);

	// Cleanup
	free(newConsensusString);
	stList_destruct(anchorAlignments);
	free(poaToConsensusMap);

	if(st_getLogLevel() >= info) {
		double score = poa_getReferenceNodeTotalMatchWeight(poa) - poa_getTotalErrorWeight(poa);
		double score2 = poa_getReferenceNodeTotalMatchWeight(poa2) - poa_getTotalErrorWeight(poa2);
		st_logInfo("Took %f seconds to do a round of polishing, got score: %f, a diff: %f \n",
					(float)(time(NULL) - startTime), score2/PAIR_ALIGNMENT_PROB_1, (score2-score)/PAIR_ALIGNMENT_PROB_1);
	}

	return poa2;
}

// Functions to iteratively polish a sequence
Poa *poa_realignIterative3(Poa *poa, stList *bamChunkReads,
						   PolishParams *polishParams, bool hmmMNotRealign,
						   int64_t minIterations, int64_t maxIterations) {
	assert(maxIterations >= 0);
	assert(minIterations <= maxIterations);

	time_t startTime = time(NULL);

	double score = poa_getReferenceNodeTotalMatchWeight(poa) - poa_getTotalErrorWeight(poa);

	st_logInfo("Starting realignment with score: %f\n", score/PAIR_ALIGNMENT_PROB_1);

	int64_t i=0;
	while(i < maxIterations) {
		i++;

		time_t consensusFindingStartTime = time(NULL);

		int64_t *poaToConsensusMap;
		char *reference = hmmMNotRealign ? poa_getConsensus(poa, &poaToConsensusMap, polishParams) :
				poa_polish2(poa, bamChunkReads, polishParams, &poaToConsensusMap);

		st_logInfo("Took %f seconds to do round %" PRIi64 " of consensus finding using algorithm %s\n",
				(float)(time(NULL) - consensusFindingStartTime), i, hmmMNotRealign ? "consensus" : "polish");

		// Stop in case consensus string is same as old reference (i.e. greedy convergence)
		if(stString_eq(reference, poa->refString)) {
			free(reference);
			free(poaToConsensusMap);
			break;
		}

		// Get anchor alignments
		stList *anchorAlignments = poa_getAnchorAlignments(poa, poaToConsensusMap, stList_length(bamChunkReads), polishParams);

		time_t realignStartTime = time(NULL);

		// Generated updated poa
		Poa *poa2 = poa_realign(bamChunkReads, anchorAlignments, reference, polishParams);

		// Cleanup
		free(reference);
		free(poaToConsensusMap);
		stList_destruct(anchorAlignments);

		double score2 = poa_getReferenceNodeTotalMatchWeight(poa2) - poa_getTotalErrorWeight(poa2);

		st_logInfo("Took %f seconds to do round %" PRIi64 " of realignment, Have score: %f (%f score diff)\n",
					(float)(time(NULL) - realignStartTime), i, score2/PAIR_ALIGNMENT_PROB_1, (score2-score)/PAIR_ALIGNMENT_PROB_1);

		// Stop if score decreases (greedy stopping)
		if(score2 <= score && i >= minIterations) {
			poa_destruct(poa2);
			break;
		}

		poa_destruct(poa);
		poa = poa2;
		score = score2;
	}

	st_logInfo("Took %f seconds to realign iterative using algorithm: %s through %" PRIi64 " iterations, got final score : %f\n",
			(float)(time(NULL) - startTime), hmmMNotRealign ? "consensus" : "polish", i, score/PAIR_ALIGNMENT_PROB_1);

	return poa;
}

Poa *poa_realignIterative2(stList *bamChunkReads,
						   stList *anchorAlignments, char *reference,
						   PolishParams *polishParams, bool hmmMNotRealign,
						   int64_t minIterations, int64_t maxIterations) {
	time_t startTime = time(NULL);
	Poa *poa = poa_realign(bamChunkReads, anchorAlignments, reference, polishParams);
	st_logInfo("Took %f seconds to generate initial POA\n", (float)(time(NULL) - startTime));
	return maxIterations == 0 ? poa : poa_realignIterative3(poa, bamChunkReads, polishParams, hmmMNotRealign, minIterations, maxIterations);
}

Poa *poa_realignIterative(stList *bamChunkReads, stList *anchorAlignments, char *reference,
						  PolishParams *polishParams) {
	return poa_realignIterative2(bamChunkReads, anchorAlignments, reference, polishParams, 1, 0, (int64_t)10000000);
}

Poa *poa_realignAll(stList *bamChunkReads, stList *anchorAlignments, char *reference,
						  PolishParams *polishParams) {
	Poa *poa = poa_realignIterative2(bamChunkReads, anchorAlignments, reference, polishParams, 1,
									 polishParams->minPoaConsensusIterations, polishParams->maxPoaConsensusIterations);
	poa = poa_realignIterative3(poa, bamChunkReads, polishParams, 0,
			polishParams->minRealignmentPolishIterations, polishParams->maxRealignmentPolishIterations);
	return poa;
}
