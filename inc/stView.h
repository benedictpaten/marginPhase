/*
 * Copyright (C) 2009-2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef VIEW_H_
#define VIEW_H_

#include "sonLib.h"

typedef struct _refMsaView {
	int64_t refLength; // The length of the reference sequence
	char *refSeq; // The reference sequence - this is not copied by the constructor
	char *refSeqName; // The reference sequence name - this is not copied by the constructor
	int64_t seqNo; // The number of non-ref sequences aligned to the reference
	stList *seqs; // The non-ref sequences - - this is not copied by the constructor
	stList *seqNames; // The non-ref sequence names - this is not copied by the constructor
	int64_t *seqCoordinates; // A matrix giving the coordinates of the non-reference sequence
	// as aligned to the reference sequence
	int64_t *maxPrecedingInsertLengths; // The maximum length of an insert in
	// any of the sequences preceding the reference positions
} MsaView;

/*
 * Get the coordinate in the given sequence aligned to the given reference position. Returns -1 if
 * no sequence position is aligned to the reference position.
 */
int64_t msaView_getSeqCoordinate(MsaView *view, int64_t refCoordinate, int64_t seqIndex);

/*
 * Gets the length of any insert in the given sequence preceding the given reference position. If
 * the sequence is not aligned at the given reference position returns 0.
 */
int64_t msaView_getPrecedingInsertLength(MsaView *view, int64_t rightRefCoordinate, int64_t seqIndex);

/*
 * Gets the first position in the sequence of an insert preceding the given reference position. If there
 * is no such insert returns -1.
 */
int64_t msaView_getPrecedingIndelStart(MsaView *view, int64_t rightRefCoordinate, int64_t seqIndex);

/*
 * Gets the maximum length of an indel preceding the given reference position
 */
int64_t msaView_getMaxPrecedingIndelLength(MsaView *view, int64_t rightRefCoordinate);

/*
 * Builds an MSA view for the given reference and aligned sequences.
 * Does not copy the strings or string names, just holds references.
 */
MsaView *msaView_construct(char *refSeq, char *refName,
						   stList *refToSeqAlignments, stList *seqs, stList *seqNames);

void msaView_destruct(MsaView * view);

/*
 * Prints a quick view of the MSA for debugging/browsing.
 */
void msaView_print(MsaView *view, FILE *fh);

#endif /* VIEW_H_ */