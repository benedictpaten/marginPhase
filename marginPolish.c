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

#include "margin_phase_version.h"
#include "margin.h"


/*
 * Main functions
 */

void usage() {
    fprintf(stderr, "usage: marginPolish <BAM_FILE> <ASSEMBLY_FASTA> <PARAMS> [options]\n");
    fprintf(stderr, "Version: %s \n\n", MARGINPHASE_MARGIN_PHASE_VERSION_H);
    fprintf(stderr, "Polishes an assembly using the reads in a BAM file and produces:\n");
    fprintf(stderr, "    1) a fasta file giving an updated reference.\n");
    fprintf(stderr, "    2) and (optionally) a SAM/BAM/CRAM file of the reads giving their alignment to the updated reference\n");

    fprintf(stderr, "\nRequired arguments:\n");
    fprintf(stderr, "    BAM_FILE BAM_FILE is the alignment of reads to the assembly (or reference).\n");
    fprintf(stderr, "    ASSEMBLY_FASTA is the reference sequence BAM file in fasta format.\n");
    fprintf(stderr, "    PARAMS is the file with marginPolish parameters.\n");

    fprintf(stderr, "\nDefault options:\n");
    fprintf(stderr, "    -h --help              : Print this help screen\n");
    fprintf(stderr, "    -a --logLevel          : Set the log level [default = info]\n");
    fprintf(stderr, "    -o --outputBase        : Name to use for output files [default = output]\n");
    fprintf(stderr, "    -r --region            : If set, will only compute for given chromosomal region.\n");
    fprintf(stderr, "                               Format: chr:start_pos-end_pos (chr3:2000-3000).\n");
}

int main(int argc, char *argv[]) {
    // Parameters / arguments
    char *logLevelString = stString_copy("info");
    char *bamInFile = NULL;
    char *paramsFile = NULL;
    char *referenceFastaFile = NULL;
    char *outputBase = stString_copy("output");
    char *regionStr = NULL;
    int64_t verboseBitstring = -1;

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
                { "outputBase", required_argument, 0, 'o'},
                { "region", required_argument, 0, 'r'},
                { "verbose", required_argument, 0, 'v'},
                { 0, 0, 0, 0 } };

        int option_index = 0;
        int key = getopt_long(argc-2, &argv[2], "a:o:v:h", long_options, &option_index);

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
        default:
            usage();
            return 0;
        }
    }

    // Initialization from arguments
    st_setLogLevelFromString(logLevelString);
    free(logLevelString);

    // Parse parameters
    st_logInfo("> Parsing model parameters from file: %s\n", paramsFile);
    FILE *fh = fopen(paramsFile, "r");
    Params *params = params_readParams(fh);
    fclose(fh);

    // Print a report of the parsed parameters
    if(st_getLogLevel() == debug) {
    	params_printParameters(params, stderr);
    }

    // Parse reference as map of header string to nucleotide sequences
    st_logInfo("> Parsing reference sequences from file: %s\n", referenceFastaFile);
    fh = fopen(referenceFastaFile, "r");
    stHash *referenceSequences = fastaReadToMap(fh);
    fclose(fh);
    // log names and provide transform
    stList *refSeqNames = stHash_getKeys(referenceSequences);
    int64_t origRefSeqLen = stList_length(refSeqNames);
    st_logDebug("\tReference contigs: \n");
    for (int64_t i = 0; i < origRefSeqLen; ++i) {
        char *fullRefSeqName = (char *) stList_get(refSeqNames, i);
        st_logDebug("\t\t%s\n", fullRefSeqName);
        char refSeqName[128] = "";
        if (sscanf(fullRefSeqName, "%s", refSeqName) == 1 && !stString_eq(fullRefSeqName, refSeqName)) {
            char *newKey = stString_copy(refSeqName);
            char *refSeq = stHash_search(referenceSequences, fullRefSeqName);
            stHash_insert(referenceSequences, newKey, refSeq);
            stHash_removeAndFreeKey(referenceSequences, fullRefSeqName);
            st_logDebug("\t\t\t-> %s\n", newKey);
        }
    }
    stList_destruct(refSeqNames);

    // Open output files
    char *polishedReferenceOutFile = stString_print("%s.fa", outputBase);
    st_logInfo("> Going to write polished reference in : %s\n", polishedReferenceOutFile);

    FILE *polishedReferenceOutFh = fopen(polishedReferenceOutFile, "w");
    free(polishedReferenceOutFile);

    // if regionStr is NULL, it will be ignored in construct2
    BamChunker *bamChunker = bamChunker_construct2(bamInFile, regionStr, params->polishParams);
    st_logInfo("> Set up bam chunker with chunk size: %i and overlap %i (for region=%s)\n",
    		   (int)bamChunker->chunkSize, (int)bamChunker->chunkBoundary, regionStr == NULL ? "all" : regionStr);

    stList *polishedReferenceStrings = NULL; // The polished reference strings, one
    // for each chunk
    char *referenceSequenceName = NULL;

    // For each chunk of the BAM
    BamChunk *bamChunk = NULL;
    while((bamChunk = bamChunker_getNext(bamChunker)) != NULL) {
    	// Get reference string for chunk of alignment
    	char *fullReferenceString = stHash_search(referenceSequences, bamChunk->refSeqName);
        if (fullReferenceString == NULL) {
            st_logCritical("> ERROR: Reference sequence missing from reference map: %s \n", bamChunk->refSeqName);
            continue;
        }
        int64_t refLen = strlen(fullReferenceString);
        char *referenceString = stString_getSubString(fullReferenceString, bamChunk->chunkBoundaryStart,
            (refLen < bamChunk->chunkBoundaryEnd ? refLen : bamChunk->chunkBoundaryEnd) - bamChunk->chunkBoundaryStart);

        st_logInfo("> Going to process a chunk for reference sequence: %s, starting at: %i and ending at: %i\n",
        		   bamChunk->refSeqName, (int)bamChunk->chunkBoundaryStart,
				   (int)(refLen < bamChunk->chunkBoundaryEnd ? refLen : bamChunk->chunkBoundaryEnd));

		// Convert bam lines into corresponding reads and alignments
		st_logInfo("> Parsing input reads from file: %s\n", bamInFile);
		stList *reads = stList_construct3(0, (void (*)(void *))bamChunkRead_destruct);
        stList *alignments = stList_construct3(0, (void (*)(void *))stList_destruct);
        convertToReadsAndAlignments(bamChunk, reads, alignments);
        // TODO: remove
//		bool *readStrandArray = st_calloc(stList_length(reads), sizeof(bool));

		Poa *poa = NULL; // The poa alignment
		char *polishedReferenceString = NULL; // The polished reference string

		// Now run the polishing method

		if(params->polishParams->useRunLengthEncoding) {
			st_logInfo("> Running polishing algorithm using run-length encoding\n");

			// Run-length encoded polishing

			// Do run length encoding (RLE)

			// First RLE the reference
			RleString *rleReference = rleString_construct(referenceString);

			// Now RLE the reads
			stList *rleNucleotides = stList_construct3(0, (void (*)(void *))rleString_destruct);
			stList *rleReads = stList_construct3(0, (void (*)(void *))bamChunkRead_destruct);
            stList *rleAlignments = stList_construct3(0, (void (*)(void*))stList_destruct);
			for(int64_t j=0; j<stList_length(reads); j++) {
				BamChunkRead *read = stList_get(reads, j);
				RleString *rleNucleotideString = rleString_construct(read->nucleotides);
				stList_append(rleNucleotides, rleNucleotideString);
				stList_append(rleReads, bamChunkRead_constructRLECopy(read, rleNucleotideString));
				stList_append(rleAlignments, runLengthEncodeAlignment(stList_get(alignments, j), rleReference, rleNucleotideString));
			}

			// Generate partial order alignment (POA)
			poa = poa_realignIterative(reads, rleAlignments, rleReference->rleString, params->polishParams);

			// Now optionally do phasing and haplotype specific polishing

			//stList *anchorAlignments = poa_getAnchorAlignments(poa, NULL, stList_length(reads), params->polishParams);
			//stList *reads1, *reads2;
			//phaseReads(poa->refString, stList_length(poa->nodes)-1, l, anchorAlignments, &reads1, &reads2, params);

			// Do run-length decoding
			RleString *polishedRLEReference = expandRLEConsensus(poa, rleReads, reads, params->polishParams->repeatSubMatrix);
			polishedReferenceString = rleString_expand(polishedRLEReference);

			// Now cleanup run-length stuff
			stList_destruct(rleNucleotides);
			stList_destruct(rleReads);
            stList_destruct(rleAlignments);
			rleString_destruct(rleReference);
			rleString_destruct(polishedRLEReference);
		}
		else { // Non-run-length encoded polishing
			st_logInfo("> Running polishing algorithm without using run-length encoding\n");

			// Generate partial order alignment (POA)
			poa = poa_realignIterative(reads, alignments, referenceString, params->polishParams);

			// Polished string is the final backbone of the POA
			polishedReferenceString = stString_copy(poa->refString);
		}

		// Log info about the POA
		if (st_getLogLevel() >= info) {
			st_logInfo("Summary stats for POA:\t");
			poa_printSummaryStats(poa, stderr);
		}
		if (st_getLogLevel() >= debug) {
			poa_print(poa, stderr, 5, 5);
		}

		// If there is no prior chunk
		if(referenceSequenceName == NULL) {
			polishedReferenceStrings = stList_construct3(0, free);
			referenceSequenceName = stString_copy(bamChunk->refSeqName);
		}
		// Else, print the prior reference sequence if current chunk not part of that sequence
		else if(!stString_eq(bamChunk->refSeqName, referenceSequenceName)) {
			assert(stList_length(polishedReferenceStrings) > 0);

			// Write the previous polished reference string out
			char *s = stString_join2("", polishedReferenceStrings);
			fastaWrite(s, referenceSequenceName, polishedReferenceOutFh);

			// Clean up
			free(s);
			stList_destruct(polishedReferenceStrings);
			free(referenceSequenceName);

			// Reset for next reference sequence
			polishedReferenceStrings = stList_construct3(0, free);
			referenceSequenceName = stString_copy(bamChunk->refSeqName);
		}
		// If there was a previous chunk then trim it's polished reference sequence
		// to remove overlap with the current chunk's polished reference sequence
		else if(stList_length(polishedReferenceStrings) > 0) {
			char *previousPolishedReferenceString = stList_peek(polishedReferenceStrings);

			// Trim the currrent and previous polished reference strings to remove overlap
			int64_t prefixStringCropEnd, suffixStringCropStart;
			int64_t overlapMatchWeight = removeOverlap(previousPolishedReferenceString, polishedReferenceString,
													   bamChunker->chunkBoundary * 2, params->polishParams,
													   &prefixStringCropEnd, &suffixStringCropStart);

			st_logInfo("Removed overlap between neighbouring chunks. Approx overlap size: %i, overlap-match weight: %f, "
					"left-trim: %i, right-trim: %i:\n", (int)bamChunker->chunkBoundary * 2, (float)overlapMatchWeight/PAIR_ALIGNMENT_PROB_1,
					strlen(previousPolishedReferenceString) - prefixStringCropEnd, suffixStringCropStart);

			// Crop the suffix of the previous chunk's polished reference string
			previousPolishedReferenceString[prefixStringCropEnd] = '\0';

			// Crop the the prefix of the current chunk's polished reference string
			char *c = polishedReferenceString;
			polishedReferenceString = stString_copy(&(polishedReferenceString[suffixStringCropStart]));
			free(c);
		}

		// Add the polished sequence to the list of polished reference sequence chunks
		stList_append(polishedReferenceStrings, polishedReferenceString);

		// Cleanup
		poa_destruct(poa);
		stList_destruct(reads);
        stList_destruct(alignments);
		free(referenceString);
    }

    // Write out the last chunk
    if(referenceSequenceName != NULL) {
    	// Write the previous polished reference string out
    	char *s = stString_join2("", polishedReferenceStrings);
    	fastaWrite(s, referenceSequenceName, polishedReferenceOutFh);

    	// Clean up
    	free(s);
    	stList_destruct(polishedReferenceStrings);
    	free(referenceSequenceName);
    }

    st_logInfo("> Finished polishing.\n");

    // Cleanup
    bamChunker_destruct(bamChunker);
    fclose(polishedReferenceOutFh);
    stHash_destruct(referenceSequences);
    params_destruct(params);

    //while(1); // Use this for testing for memory leaks

    return 0;
}

