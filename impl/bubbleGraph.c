/*
 * Copyright (C) 2019 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "margin.h"

/*
 * Bubble graphs
 */

double bubble_getLogLikelihoodOfAllele(Bubble *b, int64_t allele) {
	double logLikelihood = 0.0;
	for(int64_t i=0; i<b->readNo; i++) {
		logLikelihood += b->alleleReadSupports[allele*b->readNo + i];
	}
	return logLikelihood;
}

int64_t bubble_getIndexOfHighestLikelihoodAllele(Bubble *b) {
	int64_t maxAllele = 0;
	assert(b->alleleNo > 0);
	double maxAlleleLikelihood = bubble_getLogLikelihoodOfAllele(b, 0);
	for(int64_t i=1; i<b->alleleNo; i++) {
		double alleleLikelihood = bubble_getLogLikelihoodOfAllele(b, i);
		if(alleleLikelihood > maxAlleleLikelihood) {
			maxAllele = i;
			maxAlleleLikelihood = alleleLikelihood;
		}
	}
	return maxAllele;
}

uint64_t *bubbleGraph_getConsensusPath(BubbleGraph *bg, PolishParams *polishParams) {
	uint64_t *consensusPath = st_calloc(bg->bubbleNo, sizeof(uint64_t));
	for(int64_t i=0; i<bg->bubbleNo; i++) {
		Bubble *b = &(bg->bubbles[i]);
		consensusPath[i] = bubble_getIndexOfHighestLikelihoodAllele(b);
	}
	return consensusPath;
}

char *bubbleGraph_getConsensusString(BubbleGraph *bg, uint64_t *consensusPath, int64_t **poaToConsensusMap, PolishParams *polishParams) {
	// Map to track alignment between the new consensus sequence and the current reference sequence
	*poaToConsensusMap = st_malloc(bg->refLength * sizeof(int64_t));
	for(int64_t i=0; i<bg->refLength; i++) {
		(*poaToConsensusMap)[i] = -1;
	}

	// Substrings of the consensus string that when concatenated form the overall consensus string
	stList *consensusSubstrings = stList_construct3(0, free);
	int64_t j=0; // Index in the consensus substring
	int64_t k=0; // Index in the reference string
	for(int64_t i=0; i<bg->bubbleNo; i++) {
		Bubble *b = &(bg->bubbles[i]);

		// Add prefix after the last bubble (or start) but before the new bubble start
		if(k < b->refStart) {
			stList_append(consensusSubstrings, stString_getSubString(bg->refString, k, b->refStart-k));
			do {
				(*poaToConsensusMap)[k++] = j++;
			} while(k < b->refStart);
		}

		// Add the bubble string itself
		// Noting, if there are not sufficient numbers of sequences to call the consensus
		// use the current reference sequence
		char *consensusSubstring = stString_copy(b->alleles[consensusPath[i]]);
		stList_append(consensusSubstrings, consensusSubstring);

		if(st_getLogLevel() >= debug && 0) {
			if(strcmp(consensusSubstring, b->refAllele) != 0) {
				st_logDebug("In bubbleGraph_getConsensus, from: %" PRIi64 " to: %" PRIi64
							", \nexisting string:\t%s\nnew string:\t\t%s\n", k, k+b->length,
							b->refAllele, consensusSubstring);

				for(int64_t j=0; j<b->alleleNo; j++) {
					st_logDebug("\tGot allele: \t%s with log-likelihood: %f\n",
							b->alleles[j], bubble_getLogLikelihoodOfAllele(b, j));
				}

				for(int64_t j=0; j<b->readNo; j++) {
					st_logDebug("\tGot read: \t%s, q-value: %f\n", b->reads[j]->readSubstring, b->reads[j]->qualValue);
				}
			}
		}

		// Check if the same as the existing reference
		// in which case we can maintain the alignment
		if(strcmp(b->refAllele, consensusSubstring) == 0) {
			do {
				(*poaToConsensusMap)[k++] = j++;
			} while(k < b->refStart + b->length);
		}
		else {
			// Otherwise just update coordinates
			k += b->length;
			j += strlen(consensusSubstring);
		}
	}

	// Add the suffix of the reference after the last bubble
	if(k < bg->refLength) {
		stList_append(consensusSubstrings, stString_getSubString(bg->refString, k, bg->refLength-k));
		do {
			(*poaToConsensusMap)[k++] = j++;
		} while(k < bg->refLength);
	}

	// Build the new consensus string by concatenating the constituent pieces
	char *newConsensusString = stString_join2("", consensusSubstrings);
	assert(strlen(newConsensusString) == j);

	// Cleanup
	stList_destruct(consensusSubstrings);

	return newConsensusString;
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

double getTotalWeight(Poa *poa, PoaNode *node) {
	/*
	 * Returns the total base weight of reads aligned to the given node.
	 */
	double totalWeight = 0.0;
	for(int64_t i=0; i<poa->alphabet->alphabetSize; i++) {
		totalWeight += node->baseWeights[i];
	}
	return totalWeight;
}

double getAvgCoverage(Poa *poa, int64_t from, int64_t to) {
	// Calculate average coverage, which is used to determine candidate variants
	double avgCoverage = 0.0;
	for(int64_t j=from; j<to; j++) {
		avgCoverage += getTotalWeight(poa, stList_get(poa->nodes, j));
	}
	return avgCoverage / (to-from);
}

char getNextCandidateBase(Poa *poa, PoaNode *node, int64_t *i, double candidateWeight) {
	/*
	 * Iterates through candidate bases for a reference position returning those with sufficient weight.
	 * Always returns the reference base
	 */
	double totalWeight = getTotalWeight(poa, node);
	while(*i<poa->alphabet->alphabetSize) {
		char base = poa->alphabet->convertSymbolToChar(*i);
		if(node->baseWeights[(*i)++] > candidateWeight || toupper(node->base) == base) {
			return base;
		}
	}
	return '-';
}

bool hasCandidateSubstitution(Poa *poa, PoaNode *node, double candidateWeight) {
	/*
	 * Returns non-zero if the node has a candidate base that is different to the
	 * current base.
	 */
	int64_t i=0;
	char base;
	while((base = getNextCandidateBase(poa, node, &i, candidateWeight)) != '-') {
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
	while((base = getNextCandidateBase(poa, node, &i, candidateWeight)) != '-') { // Enumerate the possible bases at the reference node.

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

BamChunkReadSubstring *getReadSubstring(BamChunkRead *bamChunkRead, int64_t start, int64_t length, PolishParams *params) {
	assert(length >= 0);

	BamChunkReadSubstring *rs = st_calloc(1, sizeof(BamChunkReadSubstring));

	// Basic attributes
	rs->read = bamChunkRead;
	rs->start = start;
	rs->length = length;

	// Calculate the qual value
	if(bamChunkRead->qualities != NULL) {
		int64_t j = 0;
		for(int64_t i=0; i<length; i++) {
			j += (int64_t)bamChunkRead->qualities[i+start];
		}
		rs->qualValue = (double)j / length; // Quals are phred, qual = -10 * log_10(p)
	}
	else {
		rs->qualValue = -1.0;
	}

	// Add read substring - TODO: fix not thread safe
	char *read = &(bamChunkRead->rleRead->rleString[start]);
	char c = read[length];
	read[length] = '\0';
	rs->readSubstring = stString_copy(read);
	read[length] = c;

	return rs;
}

void readSubstring_destruct(BamChunkReadSubstring *rs) {
	free(rs->readSubstring);
	free(rs);
}

double computeLogLikelihoodOfConsensusString(char *reference, stList *nucleotides, PolishParams *params) {
	/*
	 * Computes the log probability of the reference given the reads.
	 */
	double logProb = LOG_ONE;
	stList *anchorPairs = stList_construct(); // Currently empty
	for(int64_t i=0; i<stList_length(nucleotides); i++) {
		BamChunkReadSubstring *rs = stList_get(nucleotides, i);
		logProb += computeForwardProbability(reference, rs->readSubstring, anchorPairs, params->p, params->sM, 0, 0);
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

int readSubstrings_cmpByQual(const void *a, const void *b) {
	/*
	 * Compares read substrings by quality in descending order
	 */
	BamChunkReadSubstring *rs1 = (BamChunkReadSubstring *)a;
	BamChunkReadSubstring *rs2 = (BamChunkReadSubstring *)b;

	return rs1->qualValue < rs2->qualValue ? 1 : (rs1->qualValue > rs2->qualValue ? -1 : 0);
}

stList *filterReadSubstrings(stList *readSubstrings, PolishParams *params) {
	// Sort the substrings by descending qvalue
	stList_sort(readSubstrings, readSubstrings_cmpByQual);

	while(stList_length(readSubstrings) > params->filterReadsWhileHaveAtLeastThisCoverage) {
		BamChunkReadSubstring *rs = stList_peek(readSubstrings);
		if(rs->qualValue >= params->minAvgBaseQuality || rs->qualValue == -1) { //Filter by qvalue, but don't filter if some or all reads
			// don't have q-values
			break;
		}
		readSubstring_destruct(rs);
		stList_pop(readSubstrings);
	}

	return readSubstrings;
}

stList *getReadSubstrings(stList *bamChunkReads, Poa *poa, int64_t from, int64_t to, PolishParams *params) {
	/*
	 * Get the substrings of reads aligned to the interval from (inclusive) to to
	 * (exclusive) and their qual values. Adds them to readSubstrings and qualValues, respectively.
	 *
	 */
	stList *readSubstrings = stList_construct3(0, (void (*)(void *))readSubstring_destruct);

	// Deal with boundary cases
	if(from == 0) {
		if(to == stList_length(poa->nodes)) {
			// If from and to reference positions that bound the complete alignment just
			// copy the complete reads
			for(int64_t i=0; i<stList_length(bamChunkReads); i++) {
			    BamChunkRead *bamChunkRead = stList_get(bamChunkReads, i);
				stList_append(readSubstrings, getReadSubstring(bamChunkRead, 0, bamChunkRead->rleRead->length, params));
			}
			return filterReadSubstrings(readSubstrings, params);
		}

		// Otherwise, include the read prefixes that end at to
		PoaNode *node = stList_get(poa->nodes, to);
		int64_t i=0;
		while(i<stList_length(node->observations)) {
			PoaBaseObservation *obs = stList_get(node->observations, i);
            BamChunkRead *bamChunkRead = stList_get(bamChunkReads, obs->readNo);
            // Trim the read substring, copy it and add to the substrings list
            stList_append(readSubstrings, getReadSubstring(bamChunkRead, 0, obs->offset, params));
			i = skipDupes(node, ++i, obs->readNo);
		}
		return filterReadSubstrings(readSubstrings, params);
	}
	else if(to == stList_length(poa->nodes)) {
		// Finally, include the read suffixs that start at from
		PoaNode *node = stList_get(poa->nodes, from);
		int64_t i = 0;
		while (i < stList_length(node->observations)) {
			PoaBaseObservation *obs = stList_get(node->observations, i);
            BamChunkRead *bamChunkRead = stList_get(bamChunkReads, obs->readNo);
			// Trim the read substring, copy it and add to the substrings list
            stList_append(readSubstrings, getReadSubstring(bamChunkRead, obs->offset, bamChunkRead->rleRead->length-obs->offset, params));
			i = skipDupes(node, ++i, obs->readNo);
		}
		return filterReadSubstrings(readSubstrings, params);
	}

	PoaNode *fromNode = stList_get(poa->nodes, from);
	PoaNode *toNode = stList_get(poa->nodes, to);

	int64_t i=0, j=0;
	while(i < stList_length(fromNode->observations) && j < stList_length(toNode->observations)) {
		PoaBaseObservation *obsFrom = stList_get(fromNode->observations, i);
		PoaBaseObservation *obsTo = stList_get(toNode->observations, j);

		if(obsFrom->readNo == obsTo->readNo) {
            BamChunkRead *bamChunkRead = stList_get(bamChunkReads, obsFrom->readNo);
            if(obsTo->offset-obsFrom->offset > 0) { // If a non zero run of bases
            	stList_append(readSubstrings, getReadSubstring(bamChunkRead, obsFrom->offset, obsTo->offset-obsFrom->offset, params));
            }
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

	return filterReadSubstrings(readSubstrings, params);
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
		totalWeight += getTotalWeight(poa, stList_get(poa->nodes, i));
		if(i >= window) {
			totalWeight -= getTotalWeight(poa, stList_get(poa->nodes, i-window));
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
		if(hasCandidateSubstitution(poa, node, candidateWeights[i]) || hasCandidateInsert(node, candidateWeights[i])) {
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
		st_logDebug("Creating filtered anchor positions got: %" PRIi64 " anchors for ref seq of length: %" PRIi64 ", that's one every: %f bases\n",
					totalAnchorNo, stList_length(poa->nodes), stList_length(poa->nodes)/(double)totalAnchorNo);
	}

	return anchors;
}


BubbleGraph *bubbleGraph_constructFromPoa(Poa *poa, stList *bamChunkReads, PolishParams *params) {
	// Setup
	double *candidateWeights = getCandidateWeights(poa, params);

	// Log info about the alignment
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

	// Make a list of bubbles
	stList *bubbles = stList_construct3(0, free);
	int64_t pAnchor = 0; // Previous anchor, starting from first position of POA, which is the prefix "N"
	for(int64_t i=1; i<stList_length(poa->nodes); i++) {
		if(anchors[i]) { // If position i is an anchor
			assert(i > pAnchor);
			if(i-pAnchor != 1)  { // In case anchors are not trivially adjacent there exists a potential bubble
				// with start coordinate on the reference sequence of pAnchor and length pAnchor-i

				// Calculate the list of alleles
				stList *alleles = NULL;
				double weightAdjustment = 1.0;
				do {
					alleles = getCandidateConsensusSubstrings(poa, pAnchor+1, i,
							candidateWeights, weightAdjustment, params->maxConsensusStrings);
					weightAdjustment *= 1.5; // Increase the candidate weight by 50%
				} while(alleles == NULL);

				// Get existing reference string
				char *existingRefSubstring = getExistingSubstring(poa, pAnchor+1, i);
				assert(strlen(existingRefSubstring) == i-pAnchor-1);

				// If it is not trivial because it contains more than one allele, or an allele different
				// to the reference
				if(stList_length(alleles) > 1 || strcmp(stList_peek(alleles), existingRefSubstring) != 0) {

					Bubble *b = st_malloc(sizeof(Bubble)); // Make a bubble and add to list of bubbles
					stList_append(bubbles, b);

					// Set the coordinates
					b->refStart = pAnchor;
					b->length = i-pAnchor-1;
					assert(b->length > 0);

					// The reference allele
					b->refAllele = existingRefSubstring;

					// Now copy the alleles list to the bubble's array of alleles
					b->alleleNo = stList_length(alleles);
					b->alleles = st_malloc(sizeof(char *) * b->alleleNo);
					for(int64_t j=0; j<b->alleleNo; j++) {
						b->alleles[j] = stList_pop(alleles);
					}

					// Get read substrings
					stList *readSubstrings = getReadSubstrings(bamChunkReads, poa, pAnchor+1, i, params);
					b->readNo = stList_length(readSubstrings);
					b->reads = st_malloc(sizeof(BamChunkReadSubstring *) * b->readNo);
					for(int64_t j=0; j<b->readNo; j++) {
						b->reads[j] = stList_pop(readSubstrings);
					}
					stList_destruct(readSubstrings);

					// Get allele supports
					b->alleleReadSupports = st_calloc(b->readNo*b->alleleNo, sizeof(float));
					stList *anchorPairs = stList_construct(); // Currently empty
					for(int64_t j=0; j<b->alleleNo; j++) {
						for(int64_t k=0; k<b->readNo; k++) {
							b->alleleReadSupports[j*b->readNo + k] =
					computeForwardProbability(b->alleles[j], b->reads[k]->readSubstring, anchorPairs, params->p, params->sM, 0, 0);
						}
					}
					stList_destruct(anchorPairs);
				}
				// Cleanup
				else {
					free(existingRefSubstring);
				}
				stList_destruct(alleles);
			}
			// Update previous anchor
			pAnchor = i;
		}
	}

	// Build the the graph

	BubbleGraph *bg = st_malloc(sizeof(BubbleGraph));
	bg->refString = poa->refString;
	bg->refLength = stList_length(poa->nodes)-1;

	// Copy the bubbles
	bg->bubbleNo = stList_length(bubbles);
	bg->bubbles = st_calloc(bg->bubbleNo, sizeof(Bubble)); // allocate bubbles
	for(int64_t i=0; i<bg->bubbleNo; i++) {
		bg->bubbles[i] = *(Bubble *)stList_get(bubbles, i);
	}

	// Fill in the bubble allele offsets
	int64_t alleleOffset = 0;
	for(int64_t i=0; i<bg->bubbleNo; i++) {
		bg->bubbles[i].alleleOffset = alleleOffset;
		alleleOffset += bg->bubbles[i].alleleNo;
	}
	bg->totalAlleles = alleleOffset;

	// Cleanup
	free(anchors);
	free(candidateWeights);
	stList_destruct(bubbles);

	return bg;
}

void bubble_destruct(Bubble b) {
	// Cleanup the reads
	for(int64_t j=0; j<b.readNo; j++) {
		free(b.reads[j]);
	}
	free(b.reads);
	// Cleanup the alleles
	for(int64_t j=0; j<b.alleleNo; j++) {
		free(b.alleles[j]);
	}
	free(b.alleles);
	// Cleanup the allele supports
	free(b.alleleReadSupports);
	// Cleanup the reference allele
	free(b.refAllele);
}

void bubbleGraph_destruct(BubbleGraph *bg) {
	// Clean up the memory for each bubble
	for(int64_t i=0; i<bg->bubbleNo; i++) {
		bubble_destruct(bg->bubbles[i]);
	}
	free(bg->bubbles);
	free(bg);
}

void bubbleGraph_print(BubbleGraph *bg, FILE *fh) {

}

stHash *bubbleGraph_getProfileSeqs(BubbleGraph *bg, stReference *ref) {
	// First calculate the length of all the profile sequences

	stHash *readEnds = stHash_construct2(NULL, free); // The last bubble the read is observed to be part of

	// For each bubble in the bubble graph
	for(uint64_t i=0; i<bg->bubbleNo; i++) {
		Bubble *b = &(bg->bubbles[i]);

		// For each read aligned to this bubble
		for(uint64_t j=0; j<b->readNo; j++) {
			BamChunkReadSubstring *s = b->reads[j];
			assert(s->read != NULL);

			// Look up the corresponding profile sequence
			uint64_t *k = stHash_search(readEnds, s->read);

			if(k == NULL) { // We are starting a new read, so make a new
				// length entry
				k = st_calloc(1, sizeof(uint64_t));
				stHash_insert(readEnds, s->read, k);
			}

			k[0] = i; // Update last time read observed aligned to a bubble
		}
	}

	stHash *readsToPSeqs = stHash_construct();

	// Now build the profile sequences

	// For each bubble in the bubble graph
	for(uint64_t i=0; i<bg->bubbleNo; i++) {
		Bubble *b = &(bg->bubbles[i]);

		// For each read aligned to this bubble
		for(uint64_t j=0; j<b->readNo; j++) {
			BamChunkReadSubstring *s = b->reads[j];
			assert(stHash_search(readEnds, s->read) != NULL);

			// Look up the corresponding profile sequence
			stProfileSeq *pSeq = stHash_search(readsToPSeqs, s->read);

			if(pSeq == NULL) { // We are starting a new read, so make a new
				// profile sequence

				// Calculate the length in bubbles of the profile sequence
				uint64_t *k = stHash_search(readEnds, s->read);
				assert(k != NULL);
				assert(i <= k[0]); // The first bubble the read is part of must precede or be equal to the last
				uint64_t pSeqLength = k[0] - i + 1;
				assert(i+pSeqLength <= bg->bubbleNo);

				pSeq = stProfileSeq_constructEmptyProfile(ref, s->read->readName, i, pSeqLength);
				stHash_insert(readsToPSeqs, s->read, pSeq);
			}

			// Sanity check the pSeq
			assert(b->alleleOffset >= pSeq->alleleOffset);
			assert(i < pSeq->refStart + pSeq->length);

			// For each allele in bubble add the prob that the read was generated by
			// the read

			// Set prob as diff to most probable allele
			uint64_t alleleOffset = b->alleleOffset-pSeq->alleleOffset;
			for(uint64_t k=0; k<b->alleleNo; k++) {
				int64_t l = -roundf(b->alleleReadSupports[b->readNo * k + j]);
				//int64_t l = roundf(2.0 * (maxLogProb - b->alleleReadSupports[b->readNo * k + j]));
				assert(l >= 0);
				pSeq->profileProbs[alleleOffset+k] = l > 255 ? 255 : l;
			}
		}
	}

	// Cleanup
	stHash_destruct(readEnds);

	return readsToPSeqs;
}

stReference *bubbleGraph_getReference(BubbleGraph *bg, char *refName, Params *params) {
	stReference *ref = st_calloc(1, sizeof(stReference));

	ref->referenceName = stString_copy(refName);
	ref->length = bg->bubbleNo;
	ref->sites = st_calloc(bg->bubbleNo, sizeof(stSite));
	ref->totalAlleles = 0;

	stList *anchorPairs = stList_construct(); // Currently empty
	for(uint64_t i=0; i<bg->bubbleNo; i++) {
		Bubble *b = &bg->bubbles[i];
		stSite *s = &ref->sites[i];
		s->alleleNumber = b->alleleNo;
		s->alleleOffset = b->alleleOffset;
		ref->totalAlleles += b->alleleNo;
		//TODO: fix: make probs proper
		s->allelePriorLogProbs = st_calloc(b->alleleNo, sizeof(uint16_t));
		s->substitutionLogProbs = st_calloc(b->alleleNo*b->alleleNo, sizeof(uint16_t));

		// We set the probability of a substitution between two different alleles to -evoScale multiplied by the forward probability of the
		// alignment of the two alleles, rounded to the nearest integer

		for(uint64_t j=0; j<b->alleleNo; j++) {
			for(uint64_t k=j; k<b->alleleNo; k++) {

				float f = -computeForwardProbability(b->alleles[j], b->alleles[k], anchorPairs, params->polishParams->p, params->polishParams->sM, 0, 0);

				int64_t l = roundf(k == j ? 0 : f * params->polishParams->hetScalingParameter);
				l = l > 255 ? 255 : l;
				assert(l >= 0);
				assert(l <= 255);

				s->substitutionLogProbs[j * b->alleleNo + k] = l;
				s->substitutionLogProbs[k * b->alleleNo + j] = l;
			}
		}
	}

	return ref;
}

/*
 * Phasing of bubble graphs
 */

void bubbleGraph_logPhasedBubbleGraph(BubbleGraph *bg, stRPHmm *hmm, stList *path,
		stHash *readsToPSeqs, stList *profileSeqs, stGenomeFragment *gF) {
	/*
	 * Sanity checks / logging for phased bubble graph
	 */

	if(st_getLogLevel() == debug) {
		// Check read partition is complete
		assert(stSet_size(gF->reads1) + stSet_size(gF->reads2) == stList_length(profileSeqs));
		stSet *intersection = stSet_getIntersection(gF->reads1, gF->reads2);
		assert(stSet_size(intersection) == 0);
		stSet_destruct(intersection);

		stRPColumn *column = hmm->firstColumn;
		assert(column->length > 0);
		uint64_t colIndex = 0, colCo = 0;

		for(uint64_t i=0; i<gF->length; i++) {
			assert(column != NULL);
			Bubble *b = &bg->bubbles[gF->refStart+i];

			stSite *s = &(hmm->ref->sites[gF->refStart+i]);
			assert(s->alleleNumber == b->alleleNo);

			assert(gF->haplotypeString1[i] < b->alleleNo);
			assert(gF->haplotypeString2[i] < b->alleleNo);
			//assert(column->depth >= b->readNo);

			if(gF->haplotypeString1[i] != gF->haplotypeString2[i]) {
				stRPCell *cell = stList_get(path, colIndex);

				double strandSkew = bubble_phasedStrandSkew(b, readsToPSeqs, gF);

				st_logDebug(">>Phasing Bubble Graph: At site: %i (of %i) with %i potential alleles got %s (%i) (log-prob: %f) for hap1 with %i reads and %s (%i) (log-prob: %f) for hap2 with %i reads (total depth %i), and ancestral allele %s (%i), genotype prob: %f, strand-skew p-value: %f\n",
						(int)i, (int)gF->length, (int)b->alleleNo,
						b->alleles[gF->haplotypeString1[i]], (int)gF->haplotypeString1[i], gF->haplotypeProbs1[i], popcount64(cell->partition),
						b->alleles[gF->haplotypeString2[i]], (int)gF->haplotypeString2[i], gF->haplotypeProbs2[i], (int)(column->depth-popcount64(cell->partition)), (int)column->depth,
						b->alleles[gF->ancestorString[i]], (int)gF->ancestorString[i], gF->genotypeProbs[i], (float)strandSkew);

				double strandSkews[b->alleleNo];
				bubble_calculateStrandSkews(b, strandSkews);

				for(uint64_t j=0; j<b->alleleNo; j++) {
					st_logDebug("\t>>Allele %i\t strand-skew: %f \t%s\t ", (int)j, (float)strandSkews[j], b->alleles[j]);
					for(uint64_t k=0; k<b->alleleNo; k++) {
						st_logDebug("%i \t", (int)s->substitutionLogProbs[j * b->alleleNo + k]);
					}
					st_logDebug("\n");
				}


				for(uint64_t k=0; k<2; k++) {
					uint64_t l=0;
					float supports[b->alleleNo];
					for(uint64_t j=0; j<b->alleleNo; j++) {
						supports[j] = 0.0;
					}

					for(uint64_t j=0; j<b->readNo; j++) {
						BamChunkReadSubstring *s = b->reads[j];
						stProfileSeq *pSeq = stHash_search(readsToPSeqs, s->read);
						assert(pSeq != NULL);
						if(stSet_search(k == 0 ? gF->reads1 : gF->reads2, pSeq) != NULL) {
							st_logDebug("\t\t>>Partition %i, read %i:\t strand %i\t ", (int)k+1, (int)l++, (int)s->read->forwardStrand);

							for(uint64_t m=0; m<b->alleleNo; m++) {
								st_logDebug("%f\t", b->alleleReadSupports[m*b->readNo + j]);
								supports[m] += roundf(b->alleleReadSupports[m*b->readNo + j]);
							}

							st_logDebug("%s\n", s->readSubstring);
						}
					}

					st_logDebug("\t\tCombined allele partition supports:\n");
					st_logDebug("\t\t\t");
					for(uint64_t j=0; j<b->alleleNo; j++) {
						st_logDebug("%f\t", supports[j]);
					}
					st_logDebug("\n");

				}
			}

			if(++colCo >= column->length) {
				colCo = 0; colIndex++;
				column = colIndex < stList_length(path) ? column->nColumn->nColumn : NULL;
				assert(column == NULL || column->length > 0);
			}
		}
		assert(colIndex == stList_length(path));

		st_logDebug(">>Fraction of bubbles skewed %f (of %i total)\n", (float)bubbleGraph_skewedBubbles(bg, readsToPSeqs, gF), (int)bg->bubbleNo);
	}
}

stGenomeFragment *bubbleGraph_phaseBubbleGraph(BubbleGraph *bg, char *refSeqName, stList *reads, Params *params) {
	/*
	 * Runs phasing algorithm to split the reads embedded in the bubble graph into two partitions.
	 */

	// Generate profile sequences and reference
	stReference *ref = bubbleGraph_getReference(bg, refSeqName, params);
	assert(ref->length == bg->bubbleNo);
	stHash *readsToPSeqs = bubbleGraph_getProfileSeqs(bg, ref);
	stList *profileSeqs = stHash_getValues(readsToPSeqs);

	if(stList_length(profileSeqs) == 0) {
		stGenomeFragment *gf = stGenomeFragment_constructEmpty(ref, 0, 0, stSet_construct(), stSet_construct());
		stList_destruct(profileSeqs);
		stHash_destruct(readsToPSeqs);
		return gf;
	}

	assert(stList_length(reads) >= stList_length(profileSeqs));
	if(stList_length(reads) != stList_length(profileSeqs)) {
		st_logInfo("In converting from reads to profile sequences have %" PRIi64 " reads and %" PRIi64 " profile sequences\n",
				stList_length(reads), stList_length(profileSeqs));
	}

	// Filter reads so that the maximum coverage depth does not exceed params->maxCoverageDepth
	st_logInfo("> Filtering reads by coverage depth\n");
	stList *filteredProfileSeqs = stList_construct();
	stList *discardedProfileSeqs = stList_construct();
	filterReadsByCoverageDepth(profileSeqs, params->phaseParams, filteredProfileSeqs, discardedProfileSeqs);
	st_logInfo("\tFiltered %" PRIi64 " reads of %" PRIi64
			" to achieve maximum coverage depth of %" PRIi64 "\n",
			stList_length(discardedProfileSeqs), stList_length(profileSeqs),
			params->phaseParams->maxCoverageDepth);
	stList_setDestructor(profileSeqs, NULL);
	stList_destruct(profileSeqs);
	profileSeqs = filteredProfileSeqs;

	// Run phasing
	stList *hmms = getRPHmms(profileSeqs, params->phaseParams);
	stRPHmm *hmm = stList_pop(hmms);
	assert(stList_length(hmms) == 0);
	stList_destruct(hmms);

	// Run the forward-backward algorithm
	stRPHmm_forwardBackward(hmm);

	// Now compute a high probability path through the hmm
	stList *path = stRPHmm_forwardTraceBack(hmm);

	assert(hmm->refStart >= 0);
	assert(hmm->refStart + hmm->refLength <= bg->bubbleNo);

	// Compute the genome fragment
	stGenomeFragment *gF = stGenomeFragment_construct(hmm, path);

	// Refine the genome fragment by repartitoning the reads iteratively
	if(params->phaseParams->roundsOfIterativeRefinement > 0) {
		stGenomeFragment_refineGenomeFragment(gF, hmm, path, params->phaseParams->roundsOfIterativeRefinement);
	}

	// Sanity checks
	assert(gF->refStart >= 0);
	assert(gF->refStart + gF->length <= bg->bubbleNo);
	assert(gF->length == hmm->refLength);

	// Check / log the result
	bubbleGraph_logPhasedBubbleGraph(bg, hmm, path, readsToPSeqs, profileSeqs, gF);

	// For reads that exceeded the coverage depth, add them back to the haplotype they fit best
	while(stList_length(discardedProfileSeqs) > 0) {
		stProfileSeq *pSeq = stList_pop(discardedProfileSeqs);
		double i = getLogProbOfReadGivenHaplotype(gF->haplotypeString1, gF->refStart, gF->length, pSeq, gF->reference);
		double j = getLogProbOfReadGivenHaplotype(gF->haplotypeString2, gF->refStart, gF->length, pSeq, gF->reference);
        //TODO is this right?  tpesout changed from (i < j ? reads2 : reads2)
		stSet_insert(i < j ? gF->reads2 : gF->reads1, pSeq);
	}
	stList_destruct(discardedProfileSeqs);

	// Cleanup
	stRPHmm_destruct(hmm, 1);
	stList_destruct(path);
	stList_destruct(profileSeqs);
	stHash_destruct(readsToPSeqs);

	return gF;
}

Poa *bubbleGraph_getNewPoa(BubbleGraph *bg, uint64_t *consensusPath, Poa *poa, stList *reads, Params *params) {

	// Get new consensus string
	int64_t *poaToConsensusMap;
	char *newConsensusString = bubbleGraph_getConsensusString(bg, consensusPath, &poaToConsensusMap, params->polishParams);

	// Get anchor alignments
	stList *anchorAlignments = poa_getAnchorAlignments(poa, poaToConsensusMap, stList_length(reads), params->polishParams);

	// Generated updated poa
	Poa *poa2 = poa_realign(reads, anchorAlignments, newConsensusString, params->polishParams);

	// Cleanup
	free(poaToConsensusMap);
	free(newConsensusString);
	stList_destruct(anchorAlignments);

	return poa2;
}

/*
 * Stuff to manage allele-strand-skew
 */

void bubble_calculateStrandSkews(Bubble *b, double *skews) {
	// Calculate the strand specific read supports
	double forwardStrandSupports[b->alleleNo];
	double reverseStrandSupports[b->alleleNo];
	uint64_t totalForward = 0, totalBackward = 0;
	for(int64_t j=0; j<b->alleleNo; j++) {
		forwardStrandSupports[j] = 0.0;
		reverseStrandSupports[j] = 0.0;
	}
	for(int64_t i=0; i<b->readNo; i++) {
		BamChunkReadSubstring *r = b->reads[i];
		double *d;
		if(r->read->forwardStrand) {
			totalForward++;
			d = forwardStrandSupports;
		}
		else {
			totalBackward++;
			d = reverseStrandSupports;
		}
		for(int64_t j=0; j<b->alleleNo; j++) {
			d[j] += b->alleleReadSupports[j*b->readNo + i];
		}
	}

	// Calculate the average allele skew
	for(int64_t j=0; j<b->alleleNo; j++) {
		skews[j] = (forwardStrandSupports[j]/totalForward - reverseStrandSupports[j]/totalBackward) /
				(fabs(forwardStrandSupports[j] + reverseStrandSupports[j])/(totalForward+totalBackward));
	}
}

double *bubbleGraph_getAlleleStrandSkews(BubbleGraph *bg) {
	double *skews = st_calloc(bg->totalAlleles, sizeof(double));
	for(int64_t i=0; i<bg->bubbleNo; i++) {
		Bubble *b = &bg->bubbles[i];
		bubble_calculateStrandSkews(b, &skews[b->alleleOffset]);
	}

	return skews;
}

double bubbleGraph_getAlleleStrandSkewCutoff(BubbleGraph *bg, Params *p) {
	double *skews = bubbleGraph_getAlleleStrandSkews(bg);

	// Calc average skew
	double avgSkew = 0.0;
	for(int64_t i=0; i<bg->totalAlleles; i++) {
		avgSkew += fabs(skews[i]);
	}
	avgSkew /= bg->totalAlleles;

	// Calc variance
	double variance = 0.0;
	for(int64_t i=0; i<bg->totalAlleles; i++) {
		double d = skews[i] - avgSkew;
		variance += d * d;
	}
	variance /= bg->totalAlleles;

	double sd = sqrt(variance);
	double cutoff = avgSkew + p->polishParams->alleleStrandSkew * sd;

	st_logDebug("Got avg. allele skew: %f, variance: %f, standard-deviation: %f, cutoff: %f from %i alleles across %i bubbles (alleleStrandSkew parameter: %f)\n",
			avgSkew, variance, sd, cutoff, (int)bg->totalAlleles, (int)bg->bubbleNo, (float)p->polishParams->alleleStrandSkew);

	free(skews);

	return cutoff;
}

uint128_t bionomialCoefficient(int64_t n, int64_t k) {
	uint128_t ans=1;
    k=k>n-k?n-k:k;
    int64_t j=1;
    for(;j<=k;j++,n--) {
        if(n%j==0) {
            ans*=n/j;
        } else if(ans%j==0) {
            ans=ans/j*n;
        } else {
            ans=(ans*n)/j;
        }
    }
    return ans;
}

double binomialPValue(int64_t n, int64_t k) {
	uint128_t j=0.0;
	k = k < n/2 ? n-k : k;
	for(int64_t i=k; i<=n; i++) {
		j += bionomialCoefficient(n, i);
	}
	return j/pow(2.0, n);
}

double bubble_phasedStrandSkew(Bubble *b, stHash *readsToPSeqs, stGenomeFragment *gf) {
	int64_t reads = 0, positives = 0;
	for(int64_t i=0; i<b->readNo; i++) {
		stProfileSeq *pSeq = stHash_search(readsToPSeqs, b->reads[i]->read);
		assert(pSeq != NULL);
		if(stSet_search(gf->reads1, pSeq) != NULL) {
			reads++;
			if(b->reads[i]->read->forwardStrand) {
				positives++;
			}
		}
		else if(stSet_search(gf->reads2, pSeq) != NULL) {
			reads++;
			if(!b->reads[i]->read->forwardStrand) {
				positives++;
			}
		}
	}
	return binomialPValue(reads, positives);
}

double bubbleGraph_skewedBubbles(BubbleGraph *bg, stHash *readsToPSeqs, stGenomeFragment *gf) {
	int64_t skewedBubbles = 0;
	for(int64_t i=0; i<bg->bubbleNo; i++) {
		skewedBubbles += bubble_phasedStrandSkew(&bg->bubbles[i], readsToPSeqs, gf) < 0.05 ? 1 : 0;
	}
	return ((float)skewedBubbles)/bg->bubbleNo;
}

bool filterByBubbleStrandSkew(Bubble *b, double *alleleSkewCutoff) {
	double skews[b->alleleNo];
	bubble_calculateStrandSkews(b, skews);
	for(int64_t i=0; i<b->alleleNo; i++) {
		if(fabs(skews[i]) > *alleleSkewCutoff) {
			return 1;
		}
	}
	return 0;
}

uint64_t bubbleGraph_filterBubbles(BubbleGraph *bg,
		bool (*filterFn)(Bubble *, void *), void *extraArg) {
	uint64_t j=0, alleleOffset = 0;
	for(int64_t i=0; i<bg->bubbleNo; i++) {
		if(!filterFn(&bg->bubbles[i], extraArg)) {
			bg->bubbles[j] = bg->bubbles[i];
			bg->bubbles[j].alleleOffset = alleleOffset;
			alleleOffset += bg->bubbles[j++].alleleNo;
		}
		else {
			bubble_destruct(bg->bubbles[i]);
		}
	}
	uint64_t filteredBubbles = bg->bubbleNo - j;
	bg->bubbleNo = j;
	bg->totalAlleles = alleleOffset;
	bg->bubbles = st_realloc(bg->bubbles, sizeof(Bubble)*bg->bubbleNo); // Resize bubbles array

	return filteredBubbles;
}

void bubbleGraph_filterBubblesByAlleleStrandSkew(BubbleGraph *bg, Params *p) {
	double alleleStrandSkewCutoff = bubbleGraph_getAlleleStrandSkewCutoff(bg, p);

	uint64_t bubblesFiltered = bubbleGraph_filterBubbles(bg, (bool (*)(Bubble *, void *))filterByBubbleStrandSkew,
			&alleleStrandSkewCutoff);

	st_logDebug(">Allele strand skew filtering: Filtered %i bubbles of %i total bubbles using allele skew cut off of: %f\n",
			(int)bubblesFiltered, (int)(bg->bubbleNo+bubblesFiltered), alleleStrandSkewCutoff);
}

bool filterByBubbleIndel(Bubble *b, void *extraArg) {
	if(b->alleleNo != 2) {
		return 1;
	}
	return 0;
	for(int64_t i=1; i<b->alleleNo; i++) {
		if(strlen(b->alleles[i]) != strlen(b->alleles[i-1])) {
			return 1;
		}
	}
	return 0;
}

void bubbleGraph_filterBubblesToRemoveIndels(BubbleGraph *bg, Params *p) {
	uint64_t bubblesFiltered = bubbleGraph_filterBubbles(bg, (bool (*)(Bubble *, void *))filterByBubbleIndel,
				NULL);

	st_logDebug("> Indel filtering: Filtered %i bubbles of %i total bubbles using indel filter\n",
				(int)bubblesFiltered, (int)(bg->bubbleNo+bubblesFiltered));
}

