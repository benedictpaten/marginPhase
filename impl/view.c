/*
 * Copyright (C) 2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "stPolish.h"
#include  "stView.h"
#include <time.h>
#include <htslib/sam.h>
#include <stRPHmm.h>

static int64_t *msaView_set(MsaView *view, int64_t refCoordinate, int64_t seqIndex) {
	return &(view->seqCoordinates[(view->refLength+1) * seqIndex + refCoordinate]);
}

int64_t msaView_getSeqCoordinate(MsaView *view, int64_t refCoordinate, int64_t seqIndex) {
	int64_t i = msaView_set(view, refCoordinate, seqIndex)[0];
	return i < 0 ? -1 : i-2;
}

int64_t msaView_getUpToSeqCoordinate(MsaView *view, int64_t refCoordinate, int64_t seqIndex) {
	int64_t i = msaView_set(view, refCoordinate, seqIndex)[0];
	return i < 0 ? -i-2 : i-2;
}

int64_t msaView_getPrecedingInsertLength(MsaView *view, int64_t rightRefCoordinate, int64_t seqIndex) {
	int64_t i = msaView_set(view, rightRefCoordinate, seqIndex)[0];
	if(i < 0) {
		return 0;
	}
	if(rightRefCoordinate == 0) {
		return i-2;
	}
	int64_t j = msaView_set(view, rightRefCoordinate-1, seqIndex)[0];
	if(j < 0) {
		return i + j - 1;
	}
	return i - j - 1;
}

int64_t msaView_getPrecedingInsertStart(MsaView *view, int64_t rightRefCoordinate, int64_t seqIndex) {
	int64_t indelLength = msaView_getPrecedingInsertLength(view, rightRefCoordinate, seqIndex);
	if(indelLength == 0) {
		return -1;
	}
	return msaView_getSeqCoordinate(view, rightRefCoordinate, seqIndex) - indelLength;
}

int64_t msaView_getMaxPrecedingInsertLength(MsaView *view, int64_t rightRefCoordinate) {
	return view->maxPrecedingInsertLengths[rightRefCoordinate];
}

int64_t msaView_getPrecedingCoverageDepth(MsaView *view, int64_t rightRefCoordinate, int64_t indelOffset) {
	return  view->precedingInsertCoverages[rightRefCoordinate][indelOffset];
}

int64_t msaView_getMaxPrecedingInsertLengthWithGivenCoverage(MsaView *view, int64_t rightRefCoordinate, int64_t minCoverage) {
	for(int64_t i=0; i<msaView_getMaxPrecedingInsertLength(view, rightRefCoordinate); i++) {
		if(msaView_getPrecedingCoverageDepth(view, rightRefCoordinate, i) < minCoverage) {
			return i;
		}
	}
	return msaView_getMaxPrecedingInsertLength(view, rightRefCoordinate);
}

MsaView *msaView_construct(char *refSeq, char *refName,
		stList *refToSeqAlignments, stList *seqs, stList *seqNames) {
	MsaView *view = st_malloc(sizeof(MsaView));

	view->refSeq = refSeq; // This is not copied
	view->refLength = strlen(refSeq);
	view->refSeqName = refName; // This is not copied
	view->seqNo = stList_length(refToSeqAlignments);
	view->seqs = seqs; // This is not copied
	view->seqNames = seqNames; // Ditto
	view->seqCoordinates = st_calloc(view->seqNo * (view->refLength+1), sizeof(int64_t)); // At each
	// reference position for each non-ref sequence stores the coordinate of the position + 1 in the non-ref sequence aligned
	// to the reference position, if non-ref sequence is aligned at that position stores then stores -1 times the index
	// of the rightmost position aligned to the prefix of the reference up to that position + 1. The plus ones are to avoid
	// dealing with difference between 0 and -0. This storage format is sufficient to represent the alignment in an easy
	// to access format

	for(int64_t i=0; i<view->seqNo; i++) {
		stList *alignment = stList_get(refToSeqAlignments, i);
		for(int64_t j=0; j<stList_length(alignment); j++) {
			stIntTuple *alignedPair = stList_get(alignment, j);
			msaView_set(view, stIntTuple_get(alignedPair, 1), i)[0] = stIntTuple_get(alignedPair, 2)+2;
		}
		msaView_set(view, view->refLength, i)[0] = strlen(stList_get(view->seqs, i)) + 2;
		int64_t k = 1;
		for(int64_t j=0; j<view->refLength; j++) {
			int64_t *l = msaView_set(view, j, i);
			if(l[0] == 0) {
				l[0] = -k;
			}
			else {
				k = l[0];
			}
		}
	}

	view->maxPrecedingInsertLengths = st_calloc(view->refLength+1, sizeof(int64_t));
	view->precedingInsertCoverages = st_calloc(view->refLength+1, sizeof(int64_t *));
	for(int64_t j=0; j<view->refLength+1; j++) {
		int64_t maxIndelLength=0;
		for(int64_t i=0; i<view->seqNo; i++) {
			int64_t k=msaView_getPrecedingInsertLength(view, j, i);
			if(k > maxIndelLength) {
				maxIndelLength = k;
			}
		}
		view->maxPrecedingInsertLengths[j] = maxIndelLength;
		view->precedingInsertCoverages[j] = st_calloc(maxIndelLength, sizeof(int64_t));
		for(int64_t i=0; i<view->seqNo; i++) {
			int64_t k=msaView_getPrecedingInsertLength(view, j, i);
			for(int64_t l=0; l<k; l++) {
				view->precedingInsertCoverages[j][l]++;
			}
		}
	}

	return view;
}

void msaView_destruct(MsaView * view) {
	free(view->maxPrecedingInsertLengths);
	free(view->seqCoordinates);
	free(view);
}

static void printRepeatChar(FILE *fh, char repeatChar, int64_t repeatCount) {
	for(int64_t i=0; i<repeatCount; i++) {
		fprintf(fh, "%c", repeatChar);
	}
}

static void printSeqName(FILE *fh, char *seqName, int64_t seqCoordinate) {
	int64_t j = strlen(seqName);
	for(int64_t i=0; i<10; i++) {
		fprintf(fh, "%c", i < j ? seqName[i] : ' ');
	}
	fprintf(fh, "\t%" PRIi64 "\t", seqCoordinate);
}

static void msaView_printP2(MsaView *view, int64_t refStart, int64_t length,
						   int64_t minInsertCoverage,
						   char (*refCharFn)(int64_t refCoordinate, void *extraArg),
						   char (*charFn)(int64_t seq, int64_t seqCoordinate, int64_t refCoordinate, void *extraArg),
						   void *extraArg,
						   FILE *fh) {
	// Calculate which indels to print
	int64_t indelLengthsToPrint[length];
	for(int64_t i=0; i<length; i++) {
		indelLengthsToPrint[i] = msaView_getMaxPrecedingInsertLengthWithGivenCoverage(view, i+refStart, minInsertCoverage);
	}

	// Print the reference
	printSeqName(fh, view->refSeqName == NULL ? "REF" : view->refSeqName, refStart);
	for(int64_t i=refStart; i<refStart+length; i++) {
		printRepeatChar(fh, '-', indelLengthsToPrint[i-refStart]);
		fprintf(fh, "%c", refCharFn(i, extraArg));
	}
	fprintf(fh, "\n");

	// Print the reads
	for(int64_t j=0; j<view->seqNo; j++) {
		if(view->seqNames == NULL) {
			char *seqName = stString_print("SEQ:%i", j);
			printSeqName(fh, seqName, msaView_getUpToSeqCoordinate(view, refStart, j));
			free(seqName);
		}
		else {
			printSeqName(fh, stList_get(view->seqNames, j), msaView_getUpToSeqCoordinate(view, refStart, j));
		}
		for(int64_t i=refStart; i<refStart+length; i++) {
			int64_t indelLength = msaView_getPrecedingInsertLength(view, i, j);
			if(indelLength > indelLengthsToPrint[i-refStart]) {
				indelLength = indelLengthsToPrint[i-refStart];
			}
			if(indelLength > 0) {
				int64_t indelStart = msaView_getPrecedingInsertStart(view, i, j);
				for(int64_t k=0; k<indelLength; k++) {
					fprintf(fh, "%c", charFn(j, indelStart+k, -1, extraArg));
				}
			}
			printRepeatChar(fh, '-', indelLengthsToPrint[i-refStart] - indelLength);

			int64_t seqCoordinate = msaView_getSeqCoordinate(view, i, j);
			if(seqCoordinate != -1) {
				fprintf(fh, "%c", charFn(j, seqCoordinate, i, extraArg));
			}
			else {
				fprintf(fh, "+");
			}
		}
		fprintf(fh, "\n");
	}
	fprintf(fh, "\n");
}

void msaView_printP(MsaView *view, int64_t minInsertCoverage,
				   char (*refCharFn)(int64_t refCoordinate, void *extraArg),
				   char (*charFn)(int64_t seq, int64_t seqCoordinate, int64_t refCoordinate, void *extraArg),
				   void *extraArg, FILE *fh) {
	int64_t width = 30;
	for(int64_t i=0; i<view->refLength; i+=width) {
		msaView_printP2(view, i, (i+width < view->refLength) ? width : view->refLength-i, minInsertCoverage,
					   refCharFn, charFn, extraArg, fh);
	}
}

char refCharBaseFn(int64_t refCoordinate, void *extraArg) {
	MsaView *view = extraArg;
	return  view->refSeq[refCoordinate];
}

char seqCharBaseFn(int64_t seq, int64_t seqCoordinate, int64_t refCoordinate, void *extraArg) {
	MsaView *view = extraArg;
	char *sequence = stList_get(view->seqs, seq);
	return view->refSeq[refCoordinate] == sequence[seqCoordinate] ? '*' : sequence[seqCoordinate];
}

void msaView_print(MsaView *view, int64_t minInsertCoverage, FILE *fh) {
	msaView_printP(view, minInsertCoverage, refCharBaseFn, seqCharBaseFn, view, fh);
}

char repeatCountToChar(int64_t repeatCount) {
	return (char)(repeatCount + 48);
}

char refCharRepeatCountFn(int64_t refCoordinate, void *extraArg) {
	RleString *refString = ((void **)extraArg)[0];
	return repeatCountToChar(refString->repeatCounts[refCoordinate]);
}

char seqCharRepeatCountFn(int64_t seq, int64_t seqCoordinate, int64_t refCoordinate, void *extraArg) {
	RleString *refString = ((void **)extraArg)[0];
	stList *rleStrings = ((void **)extraArg)[1];
	RleString *rleString = stList_get(rleStrings, seq);

	int64_t refRepeatCount = refCoordinate >= 0 ? refString->repeatCounts[refCoordinate] : -1;
	int64_t seqRepeatCount = rleString->repeatCounts[seqCoordinate];

	return refRepeatCount == seqRepeatCount ? '*' : repeatCountToChar(seqRepeatCount);
}

void msaView_printRepeatCounts(MsaView *view, int64_t minInsertCoverage,
							   RleString *refString, stList *rleStrings, FILE *fh) {
	msaView_printP(view, minInsertCoverage, refCharRepeatCountFn, seqCharRepeatCountFn, (const void *[]){ refString, rleStrings }, fh);
}

