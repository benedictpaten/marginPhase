


static void adjustAnchors(stList *anchorPairs, int64_t index, int64_t adjustment) {
	for(int64_t i=0; i<stList_length(anchorPairs); i++) {
		stIntTuple *pair = stList_get(anchorPairs, i);
		((int64_t *)pair)[index+1] += adjustment;
	}
}

void getAlignedPairsWithIndelsCroppingReference(char *reference, int64_t refLength,
		char *read, stList *anchorPairs,
		stList **matches, stList **inserts, stList **deletes, PairwiseAlignmentParameters *p, StateMachine *sM) {
	// Crop reference, to avoid long unaligned prefix and suffix
	// that generates a lot of delete pairs

	// Get cropping coordinates
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
	getAlignedPairsWithIndelsUsingAnchors(sM, &(reference[firstRefPosition]), read,
										  anchorPairs, p, matches, deletes, inserts, 0, 0);

	// De-crop reference
	reference[endRefPosition] = c;

	// Adjust back anchors
	adjustAnchors(anchorPairs, 0, firstRefPosition);

	// Shift matches/inserts/deletes
	adjustAnchors(*matches, 1, firstRefPosition);
	adjustAnchors(*inserts, 1, firstRefPosition);
	adjustAnchors(*deletes, 1, firstRefPosition);
}
