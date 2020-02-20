/*
 * Copyright (C) 2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <getopt.h>
#include <stdio.h>
#include <ctype.h>
#include <memory.h>
#include <hashTableC.h>
#include <unistd.h>

#include "marginVersion.h"
#include "margin.h"
#include "htsIntegration.h"


/*
 * Main functions
 */

void usage() {
    fprintf(stderr, "usage: margin <BAM_FILE> <ASSEMBLY_FASTA> <PARAMS> [options]\n");
    fprintf(stderr, "Version: %s \n\n", MARGIN_POLISH_VERSION_H);
    fprintf(stderr, "Polishes an assembly using the reads in a BAM file and produces polished sequences using a haploid or diploid model:\n");
    fprintf(stderr, "    1) a fasta file giving an updated reference.\n");
    fprintf(stderr, "    2) and (optionally) a set of outputs useful further polishing algorithms\n");

    fprintf(stderr, "\nRequired arguments:\n");
    fprintf(stderr, "    BAM_FILE BAM_FILE is the alignment of reads to the assembly (or reference).\n");
    fprintf(stderr, "    ASSEMBLY_FASTA is the reference sequence BAM file in fasta format.\n");
    fprintf(stderr, "    PARAMS is the file with marginPolish parameters.\n");

    fprintf(stderr, "\nDefault options:\n");
    fprintf(stderr, "    -h --help              : Print this help screen\n");
    fprintf(stderr, "    -d --diploid           : Do diploid polishing, outputting two polished sequences per reference sequence\n");
    fprintf(stderr, "    -a --logLevel          : Set the log level [default = info]\n");
    fprintf(stderr, "    -o --outputBase        : Name to use for output files [default = output]\n");
    fprintf(stderr, "    -r --region            : If set, will only compute for given chromosomal region.\n");
    fprintf(stderr, "                               Format: chr:start_pos-end_pos (chr3:2000-3000).\n");
    fprintf(stderr, "    -i --outputRepeatCounts        : File to write out the repeat counts [default = NULL]\n");
    fprintf(stderr, "    -j --outputPoaTsv        : File to write out the poa as TSV file [default = NULL]\n");
}

stHash *parseReferenceSequences(char *referenceFastaFile) {
	/*
	 * Get hash of reference sequence names in fasta to their sequences, doing some munging on the sequence names.
	 */
	st_logInfo("> Parsing reference sequences from file: %s\n", referenceFastaFile);
	FILE *fh = fopen(referenceFastaFile, "r");
	stHash *referenceSequences = fastaReadToMap(fh);  //valgrind says blocks from this allocation are "still reachable"
	fclose(fh);
	// log names and transform (if necessary)
	stList *refSeqNames = stHash_getKeys(referenceSequences);
	int64_t origRefSeqLen = stList_length(refSeqNames);
	st_logDebug("\tReference contigs: \n");
	for (int64_t i = 0; i < origRefSeqLen; ++i) {
		char *fullRefSeqName = (char *) stList_get(refSeqNames, i);
		st_logDebug("\t\t%s\n", fullRefSeqName);
		char refSeqName[128] = "";
		if (sscanf(fullRefSeqName, "%s", refSeqName) == 1 && !stString_eq(fullRefSeqName, refSeqName)) {
			// this transformation is necessary for cases where the reference has metadata after the contig name:
			// >contig001 length=1000 date=1999-12-31
			char *newKey = stString_copy(refSeqName);
			char *refSeq = stHash_search(referenceSequences, fullRefSeqName);
			stHash_insert(referenceSequences, newKey, refSeq);
			stHash_removeAndFreeKey(referenceSequences, fullRefSeqName);
			st_logDebug("\t\t\t-> %s\n", newKey);
		}
	}
	stList_destruct(refSeqNames);

	return referenceSequences;
}

RleString *bamChunk_getReferenceSubstring(BamChunk *bamChunk, stHash *referenceSequences, Params *params) {
	/*
	 * Get corresponding substring of the reference for a given bamChunk.
	 */
	char *fullReferenceString = stHash_search(referenceSequences, bamChunk->refSeqName);
	if (fullReferenceString == NULL) {
		st_logCritical("> ERROR: Reference sequence missing from reference map: %s \n", bamChunk->refSeqName);
		return NULL;
	}
	int64_t refLen = strlen(fullReferenceString);
	char *referenceString = stString_getSubString(fullReferenceString, bamChunk->chunkBoundaryStart,
		(refLen < bamChunk->chunkBoundaryEnd ? refLen : bamChunk->chunkBoundaryEnd) - bamChunk->chunkBoundaryStart);

	RleString *rleRef = params->polishParams->useRunLengthEncoding ? rleString_construct(referenceString) : rleString_construct_no_rle(referenceString);
	free(referenceString);

	return rleRef;
}

typedef struct _polishedReferenceSequence {
	/*
	 * Object for managing the output of a polished reference sequence.
	 */
	char *referenceSequenceName;
	char *referenceSequenceNameSuffix; // Suffix appended to the name of each reference sequence when written
	// out, allows the distinction on the haplotypes
	char *referenceSequenceNameForPrinting;
	stList *polishedReferenceStrings;
	FILE *polishedReferenceFileHandle;
	FILE *outputPoaTsvFileHandle;
	FILE *outputRepeatCountFileHandle;
} PolishedReferenceSequence;

PolishedReferenceSequence *polishedReferenceSequence_construct(Params *params, char *referenceSequenceNameSuffix,
		FILE *polishedReferenceFileHandle, FILE *outputPoaTsvFileHandle, FILE *outputRepeatCountFileHandle) {
	PolishedReferenceSequence *rSeq = st_calloc(1, sizeof(PolishedReferenceSequence));

	rSeq->referenceSequenceNameSuffix = stString_copy(referenceSequenceNameSuffix);
	rSeq->polishedReferenceFileHandle = polishedReferenceFileHandle;
	rSeq->outputPoaTsvFileHandle = outputPoaTsvFileHandle;
	rSeq->outputRepeatCountFileHandle = outputRepeatCountFileHandle;

	return rSeq;
}

void polishedReferenceSequence_processChunkSequence(PolishedReferenceSequence *rSeq,
		BamChunk *bamChunk, Poa *poa, stList *reads, Params *params) {
	// Do run-length decoding
	//poa_estimateRepeatCountsUsingBayesianModel(poa, reads, params->polishParams->repeatSubMatrix);
	char *polishedReferenceString = rleString_expand(poa->refString);

	// Log info about the POA
	if (st_getLogLevel() >= info) {
		st_logInfo("Summary stats for POA:\t");
		poa_printSummaryStats(poa, stderr);
	}

	// Write any optional outputs about repeat count and POA, etc.
	if(rSeq->outputPoaTsvFileHandle != NULL) {
		poa_printTSV(poa, rSeq->outputPoaTsvFileHandle, reads, 5, 0);
	}
	if(rSeq->outputRepeatCountFileHandle != NULL) {
		poa_printRepeatCounts(poa, rSeq->outputRepeatCountFileHandle, reads);
	}

	// If there is no prior chunk
	if(rSeq->referenceSequenceName == NULL) {
		assert(rSeq->polishedReferenceStrings == NULL);
		rSeq->polishedReferenceStrings = stList_construct3(0, free);
		rSeq->referenceSequenceName = stString_copy(bamChunk->refSeqName);
		assert(rSeq->referenceSequenceNameForPrinting == NULL);
		rSeq->referenceSequenceNameForPrinting = stString_print("%s%s", rSeq->referenceSequenceName, rSeq->referenceSequenceNameSuffix);
	}
	// Else, print the prior reference sequence if current chunk not part of that sequence
	else if(!stString_eq(bamChunk->refSeqName, rSeq->referenceSequenceName)) {
		assert(stList_length(rSeq->polishedReferenceStrings) > 0);

		// Write the previous polished reference string out
		char *s = stString_join2("", rSeq->polishedReferenceStrings);
		fastaWrite(s, rSeq->referenceSequenceNameForPrinting, rSeq->polishedReferenceFileHandle);

		// Clean up
		free(s);
		stList_destruct(rSeq->polishedReferenceStrings);
		free(rSeq->referenceSequenceName);
		free(rSeq->referenceSequenceNameForPrinting);

		// Reset for next reference sequence
		rSeq->polishedReferenceStrings = stList_construct3(0, free);
		rSeq->referenceSequenceName = stString_copy(bamChunk->refSeqName);
		rSeq->referenceSequenceNameForPrinting = stString_print("%s%s", rSeq->referenceSequenceName, rSeq->referenceSequenceNameSuffix);
	}
	// If there was a previous chunk then trim it's polished reference sequence
	// to remove overlap with the current chunk's polished reference sequence
	else if(stList_length(rSeq->polishedReferenceStrings) > 0) {
		char *previousPolishedReferenceString = stList_peek(rSeq->polishedReferenceStrings);

		// Trim the currrent and previous polished reference strings to remove overlap
		int64_t prefixStringCropEnd, suffixStringCropStart;
		int64_t overlapMatchWeight = removeOverlap(previousPolishedReferenceString, polishedReferenceString,
												   bamChunk->parent->chunkBoundary * 2, params->polishParams,
												   &prefixStringCropEnd, &suffixStringCropStart);

		st_logInfo("Removed overlap between neighbouring chunks. Approx overlap size: %i, overlap-match weight: %f, "
				"left-trim: %i, right-trim: %i:\n", (int)bamChunk->parent->chunkBoundary * 2, (float)overlapMatchWeight/PAIR_ALIGNMENT_PROB_1,
				strlen(previousPolishedReferenceString) - prefixStringCropEnd, suffixStringCropStart);

		// Crop the suffix of the previous chunk's polished reference string
		previousPolishedReferenceString[prefixStringCropEnd] = '\0';

		// Crop the the prefix of the current chunk's polished reference string
		char *c = polishedReferenceString;
		polishedReferenceString = stString_copy(&(polishedReferenceString[suffixStringCropStart]));
		free(c);
	}

	// Add the polished sequence to the list of polished reference sequence chunks
	stList_append(rSeq->polishedReferenceStrings, polishedReferenceString);
}

void polishedReferenceSequence_flush(PolishedReferenceSequence *rSeq) {
	// Write out the last chunk
    if(rSeq->referenceSequenceName != NULL) {
    	// Write the previous polished reference string out
    	char *s = stString_join2("", rSeq->polishedReferenceStrings);
    	assert(rSeq->referenceSequenceNameForPrinting != NULL);
    	fastaWrite(s, rSeq->referenceSequenceNameForPrinting, rSeq->polishedReferenceFileHandle);

    	// Clean up
    	free(s);
    	stList_destruct(rSeq->polishedReferenceStrings);
    	free(rSeq->referenceSequenceName);
    	free(rSeq->referenceSequenceNameForPrinting);
    	rSeq->referenceSequenceName = NULL;
    }
}

void polishedReferenceSequence_destruct(PolishedReferenceSequence *rSeq) {
	polishedReferenceSequence_flush(rSeq);
	free(rSeq->referenceSequenceNameSuffix);
	free(rSeq);
}

uint64_t *getPaddedHaplotypeString(uint64_t *hap, stGenomeFragment *gf, BubbleGraph *bg, Params *params) {
	/*
	 * Pads a haplotype string from the genome fragment to account for any missing prefix or suffix.
	 */
	uint64_t *paddedHap = bubbleGraph_getConsensusPath(bg, params->polishParams);

	for(uint64_t i=0; i<gf->length; i++) {
		paddedHap[i+gf->refStart] = hap[i];
	}

	return paddedHap;
}

int main(int argc, char *argv[]) {
    // Parameters / arguments
    char *logLevelString = stString_copy("info");
    bool diploid = 0; // By default assuume a haploid model
    char *bamInFile = NULL;
    char *paramsFile = NULL;
    char *referenceFastaFile = NULL;
    char *outputBase = stString_copy("output");
    char *regionStr = NULL;
    int64_t verboseBitstring = -1;
    char *outputRepeatCountFile = NULL;
    char *outputPoaTsvFile = NULL;

    // TODO: When done testing, optionally set random seed using st_randomSeed();

    if(argc < 4) {
        usage();
        return 0;
    }

    bamInFile = stString_copy(argv[1]);
    referenceFastaFile = stString_copy(argv[2]);
    paramsFile = stString_copy(argv[3]);

    // Parse the options
    while (1) {
        static struct option long_options[] = {
                { "logLevel", required_argument, 0, 'a' },
                { "help", no_argument, 0, 'h' },
				{ "diploid", no_argument, 0, 'd'},
                { "outputBase", required_argument, 0, 'o'},
                { "region", required_argument, 0, 'r'},
                { "verbose", required_argument, 0, 'v'},
				{ "outputRepeatCounts", required_argument, 0, 'i'},
				{ "outputPoaTsv", required_argument, 0, 'j'},
                { 0, 0, 0, 0 } };

        int option_index = 0;
        int key = getopt_long(argc-2, &argv[2], "a:o:v:r:hdi:j:", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
        case 'a':
            free(logLevelString);
            logLevelString = stString_copy(optarg);
            break;
        case 'h':
            usage();
            return 0;
        case 'd':
            diploid = 1;
            break;
        case 'o':
            free(outputBase);
            outputBase = stString_copy(optarg);
            break;
        case 'r':
            regionStr = stString_copy(optarg);
            break;
        case 'v':
            verboseBitstring = atoi(optarg);
            break;
        case 'i':
        	outputRepeatCountFile = stString_copy(optarg);
        	break;
        case 'j':
        	outputPoaTsvFile = stString_copy(optarg);
        	break;
        default:
            usage();
            return 0;
        }
    }

    // Initialization from arguments
    st_setLogLevelFromString(logLevelString);
    free(logLevelString);

    // Parse parameters
    st_logInfo("> Using the diploid model: %s\n", diploid ? "True" : "False");
    st_logInfo("> Parsing model parameters from file: %s\n", paramsFile);
    Params *params = params_readParams(paramsFile);

    // Print a report of the parsed parameters
    if(st_getLogLevel() == debug) {
    	params_printParameters(params, stderr);
    }

    // Parse reference as map of header string to nucleotide sequences
    stHash *referenceSequences = parseReferenceSequences(referenceFastaFile);

    // Make output formatting object(s)
    char *polishedReferenceOutFile = stString_print("%s.fa", outputBase);
    st_logInfo("> Going to write polished reference in : %s\n", polishedReferenceOutFile);
    FILE *polishedReferenceFileHandle = fopen(polishedReferenceOutFile, "w");
    FILE *outputPoaTsvFileHandle = outputPoaTsvFile != NULL ? fopen(outputPoaTsvFile, "w") : NULL;
    FILE *outputRepeatCountFileHandle = outputRepeatCountFile != NULL ? fopen(outputRepeatCountFile, "w") : NULL;

    PolishedReferenceSequence *rSeq1 = polishedReferenceSequence_construct(params, diploid ? "_hap_1" : "",
    		polishedReferenceFileHandle, outputPoaTsvFileHandle, outputRepeatCountFileHandle), *rSeq2 = NULL;
    if(diploid) {
    	rSeq2 = polishedReferenceSequence_construct(params, "_hap_2",
    			polishedReferenceFileHandle, outputPoaTsvFileHandle, outputRepeatCountFileHandle);
    }
	free(polishedReferenceOutFile);

    // if regionStr is NULL, it will be ignored in construct2
    BamChunker *bamChunker = bamChunker_construct2(bamInFile, regionStr, params->polishParams);
    st_logInfo("> Set up bam chunker with chunk size: %i and overlap %i (for region=%s)\n",
    		   (int)bamChunker->chunkSize, (int)bamChunker->chunkBoundary, regionStr == NULL ? "all" : regionStr);

    // For each chunk of the BAM
	for (int64_t chunkIdx = 0; chunkIdx < bamChunker->chunkCount; chunkIdx++) {
		BamChunk *bamChunk = bamChunker_getChunk(bamChunker, chunkIdx);
    	RleString *reference = bamChunk_getReferenceSubstring(bamChunk, referenceSequences, params);

        st_logInfo("> Going to process a chunk for reference sequence: %s, starting at: %i and ending at: %i\n",
        		   bamChunk->refSeqName, (int)bamChunk->chunkBoundaryStart,
				   (int)bamChunk->chunkBoundaryEnd);

		// Convert bam lines into corresponding reads and alignments
		st_logInfo("> Parsing input reads from file: %s\n", bamInFile);
		stList *reads = stList_construct3(0, (void (*)(void *))bamChunkRead_destruct);
        stList *alignments = stList_construct3(0, (void (*)(void *))stList_destruct);
        convertToReadsAndAlignments(bamChunk, reference, reads, alignments);

		// Now run the polishing method

		// Generate the haploid partial order alignment (POA)
		Poa *poa = poa_realignAll(reads, alignments, reference, params->polishParams);

		// If diploid
		if(diploid) {
			// Get the bubble graph representation
			bool useReadAlleles = params->polishParams->useReadAlleles;
			params->polishParams->useReadAlleles = params->polishParams->useReadAllelesInPhasing;
			BubbleGraph *bg = bubbleGraph_constructFromPoa(poa, reads, params->polishParams);
			params->polishParams->useReadAlleles = useReadAlleles;

			// Now make a POA for each of the haplotypes
			stHash *readsToPSeqs;
			stGenomeFragment *gf = bubbleGraph_phaseBubbleGraph(bg, bamChunk->refSeqName, reads, params, &readsToPSeqs);

			stSet *readsBelongingToHap1, *readsBelongingToHap2;
			stGenomeFragment_phaseBamChunkReads(gf, readsToPSeqs, reads, &readsBelongingToHap1, &readsBelongingToHap2);
			st_logInfo("After phasing, of %i reads got %i reads partitioned into hap1 and %i reads partitioned into hap2 (%i unphased)\n",
			(int)stList_length(reads), (int)stSet_size(readsBelongingToHap1), (int)stSet_size(readsBelongingToHap2), 
			(int)(stList_length(reads) - stSet_size(readsBelongingToHap1)- stSet_size(readsBelongingToHap2)));

			// Debug report of hets
			uint64_t totalHets = 0;
			for(uint64_t i=0; i<gf->length; i++) {
				Bubble *b = &bg->bubbles[i+gf->refStart];
				if(gf->haplotypeString1[i] != gf->haplotypeString2[i]) {
					st_logDebug("Got predicted het at bubble %i %s %s\n", (int)i+gf->refStart, b->alleles[gf->haplotypeString1[i]]->rleString,
							b->alleles[gf->haplotypeString2[i]]->rleString);
					totalHets++;
				}
				else if(!rleString_eq(b->alleles[gf->haplotypeString1[i]], b->refAllele)) {
					st_logDebug("Got predicted hom alt at bubble %i %i\n", (int)i+gf->refStart, (int)gf->haplotypeString1[i]);
				}
			}
			st_logInfo("In phasing chunk, got: %i hets from: %i total sites (fraction: %f)\n", (int)totalHets, (int)gf->length, (float)totalHets/gf->length);

			st_logInfo("Building POA for each haplotype\n");
			uint64_t *hap1 = getPaddedHaplotypeString(gf->haplotypeString1, gf, bg, params);
			uint64_t *hap2 = getPaddedHaplotypeString(gf->haplotypeString2, gf, bg, params);

			Poa *poa_hap1 = bubbleGraph_getNewPoa(bg, hap1, poa, reads, params);
			Poa *poa_hap2 = bubbleGraph_getNewPoa(bg, hap2, poa, reads, params);

			/*st_logInfo("Using read phasing to reestimate bases in phased manner\n");
			poa_estimatePhasedBasesUsingBayesianModel(poa_hap1, reads,
					readsBelongingToHap1, readsBelongingToHap2, params->polishParams);

			poa_estimatePhasedBasesUsingBayesianModel(poa_hap2, reads,
								readsBelongingToHap2, readsBelongingToHap1, params->polishParams);*/

			if(params->polishParams->useRunLengthEncoding) {
				st_logInfo("Using read phasing to reestimate repeat counts in phased manner\n");
				poa_estimatePhasedRepeatCountsUsingBayesianModel(poa_hap1, reads,
						params->polishParams->repeatSubMatrix, readsBelongingToHap1, readsBelongingToHap2, params->polishParams);

				poa_estimatePhasedRepeatCountsUsingBayesianModel(poa_hap2, reads,
						params->polishParams->repeatSubMatrix, readsBelongingToHap2, readsBelongingToHap1, params->polishParams);
			}

			// Output
			polishedReferenceSequence_processChunkSequence(rSeq1, bamChunk, poa_hap1, reads, params);
			polishedReferenceSequence_processChunkSequence(rSeq2, bamChunk, poa_hap2, reads, params);

			// Cleanup
			free(hap1);
			free(hap2);
			bubbleGraph_destruct(bg);
			stGenomeFragment_destruct(gf);
			poa_destruct(poa_hap1);
			poa_destruct(poa_hap2);
			stSet_destruct(readsBelongingToHap1);
			stSet_destruct(readsBelongingToHap2);
			stHash_destruct(readsToPSeqs);
		}
		else {
			polishedReferenceSequence_processChunkSequence(rSeq1, bamChunk, poa, reads, params);
		}

		// Cleanup
		poa_destruct(poa);
		stList_destruct(reads);
        stList_destruct(alignments);
        rleString_destruct(reference);
    }

    polishedReferenceSequence_destruct(rSeq1);
    if(diploid) {
    	polishedReferenceSequence_destruct(rSeq2);
    }

    // Cleanup
    st_logInfo("> Finished polishing.\n");
    fclose(polishedReferenceFileHandle);
    if(outputPoaTsvFile != NULL) {
    	fclose(outputPoaTsvFileHandle);
    	free(outputPoaTsvFile);
    }
    if(outputPoaTsvFileHandle != NULL) {
    	fclose(outputRepeatCountFileHandle);
    	free(outputRepeatCountFile);
    }
    bamChunker_destruct(bamChunker);
    stHash_destruct(referenceSequences);
    params_destruct(params);

	free(bamInFile);
	free(referenceFastaFile);
	free(paramsFile);

	if (outputBase != NULL) free(outputBase);
	if (regionStr != NULL) free(regionStr);

    //while(1); // Use this for testing for memory leaks

    return 0;
}

