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
#include "sonLib.h"
#include "stPolish.h"

/*
 * Main functions
 */

void usage() {
    fprintf(stderr, "usage: marginPolish <SAM/BAM/CRAM> <REFERENCE_FASTA> <PARAMS> [options]\n");
    fprintf(stderr, "Version: %s \n\n", MARGINPHASE_MARGIN_PHASE_VERSION_H);
    fprintf(stderr, "Polishes an assembly using the reads in a BAM file and produces:\n");
    fprintf(stderr, "    1) a fasta file giving an updated reference.\n");
    fprintf(stderr, "    2) and (optionally) a SAM/BAM/CRAM file of the reads giving their alignment to the updated reference\n");

    fprintf(stderr, "\nRequired arguments:\n");
    fprintf(stderr, "    SAM/BAM/CRAM is the alignment of reads.  All reads must be aligned to the same contig \n");
    fprintf(stderr, "        and be in sam/bam/cram format.\n");
    fprintf(stderr, "    REFERENCE_FASTA is the reference sequence for the SAM/BAM/CRAM's contig in fasta format.\n");
    fprintf(stderr, "    PARAMS is the file with marginPolish parameters.\n");

    fprintf(stderr, "\nDefault options:\n");
    fprintf(stderr, "    -h --help              : Print this help screen\n");
    fprintf(stderr, "    -a --logLevel          : Set the log level [default = info]\n");
    fprintf(stderr, "    -o --outputBase        : Name to use for output files [default = output]\n");
    fprintf(stderr, "    -r --region            : If set, will only compute for this region.\n");
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
    PolishParams *params = polishParams_readParams(fh);
    fclose(fh);

    // Print a report of the parsed parameters
    if(st_getLogLevel() == debug) {
    	polishParams_printParameters(params, stderr);
    }

    // Parse reference as map of header string to nucleotide sequences
    st_logInfo("> Parsing reference sequences from file: %s\n", referenceFastaFile);
    fh = fopen(referenceFastaFile, "r");
    stHash *referenceSequences = fastaReadToMap(fh);
    fclose(fh);

    // Open output files
    char *referenceOutFile = stString_print("%s.fa", outputBase);
    char *bamOutFile = stString_print("%s.bam", outputBase);
    st_logInfo("> Going to write polished reference in : %s\n", referenceOutFile);
    st_logInfo("> Going to write polished bam file in : %s\n", bamOutFile);
    FILE *referenceOutFh = fopen(referenceOutFile, "w");
    FILE *bamOutFh = fopen(bamOutFile, "w");
    free(referenceOutFile);
    free(bamOutFile);

    // if regionStr is NULL, it will be ignored in construct2
    BamChunker *bamChunker = bamChunker_construct2(bamInFile, regionStr, params);

    // For each chunk of the BAM
    BamChunk *bamChunk = NULL;
    while((bamChunk = bamChunker_getNext(bamChunker)) != NULL) {
    	// Get reference string for chunk of alignment
    	char *fullReferenceString = stHash_search(referenceSequences, bamChunk->refSeqName);
        int64_t refLen = strlen(fullReferenceString);
        char *referenceString = stString_getSubString(fullReferenceString, bamChunk->chunkBoundaryStart,
            (refLen < bamChunk->chunkBoundaryEnd ? refLen : bamChunk->chunkBoundaryEnd) - bamChunk->chunkBoundaryStart);

		// Convert bam lines into corresponding reads and alignments
		st_logInfo("> Parsing input reads from file: %s\n", bamInFile);
		stList *reads = stList_construct3(0, free);
		stList *alignments = stList_construct3(0, (void (*)(void *))stList_destruct);
		convertToReadsAndAlignments(bamChunk, reads, alignments);

		// Output data structures
		char *consensusReferenceString = NULL; // The polished reference string
		stList *updatedAlignments = NULL; // TODO
		Poa *poa = NULL; // The partial order alignment
		stList *readToConsensusAlignments; // Final alignments of reads to reference

		// Now run the polishing method

		if(params->useRunLengthEncoding) {
			st_logInfo("> Running polishing algorithm using run-length encoding\n");

			// Run-length encoded polishing

			// Do run length encoding (RLE)

			// First RLE the reference
			RleString *rleReference = rleString_construct(referenceString);

			// Now RLE the reads
			stList *rleReads = stList_construct3(0, (void (*)(void *))rleString_destruct);
			stList *l = stList_construct(); // Just the rle nucleotide strings
			stList *rleAlignments = stList_construct3(0, (void (*)(void *))stList_destruct);
			for(int64_t j=0; j<stList_length(reads); j++) {
				char *read = stList_get(reads, j);
				RleString *rleRead = rleString_construct(read);
				stList_append(rleReads, rleRead);
				stList_append(l, rleRead->rleString);
				stList_append(rleAlignments, runLengthEncodeAlignment(stList_get(alignments, j), rleReference, rleRead));
			}

			// Generate partial order alignment (POA)
			poa = poa_realignIterative(l, rleAlignments, rleReference->rleString, params);

			// Do run-length decoding
			RleString *consensusReference = expandRLEConsensus(poa, rleReads, params->repeatSubMatrix);
			consensusReferenceString = rleString_expand(consensusReference);

			// Generate final MEA alignments in RLE space
			//stList *rleMEAAlignments = poa_getReadAlignmentsToConsensus(poa, reads, params);

			// get anchor alignments
			// get pairwise alignments with mea function
			// do left shifts
			//TODO

			// Expand alignments into non-RLE space
			//TODO

			// Now cleanup run-length stuff
			stList_destruct(rleReads);
			stList_destruct(l);
			stList_destruct(rleAlignments);
			rleString_destruct(rleReference);
		}
		else { // Non-run-length encoded polishing
			st_logInfo("> Running polishing algorithm without using run-length encoding\n");

			// Generate partial order alignment (POA)
			poa = poa_realignIterative(reads, alignments, referenceString, params);

			// Consensus is the final backbone of the POA
			consensusReferenceString = stString_copy(poa->refString);

			// Generate updated alignments
			readToConsensusAlignments = poa_getReadAlignmentsToConsensus(poa, reads, params);

			// Clean up
			poa_destruct(poa);
		}

		// Log info about the POA
		if (st_getLogLevel() >= info) {
			st_logInfo("Summary stats for POA:\t");
			poa_printSummaryStats(poa, stderr);
		}
		if (st_getLogLevel() >= debug) {
			poa_print(poa, stderr, 5);
		}

		// Output the finished sequence
        char *chunkIdentifier = stString_print("%s:%"PRIu64"-%"PRIu64, bamChunk->refSeqName,
                                               bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
		fastaWrite(consensusReferenceString, chunkIdentifier, referenceOutFh);
        free(chunkIdentifier);

		//TODO Write bam reads

		// Cleanup
		poa_destruct(poa);
		free(consensusReferenceString);
    }

    // Cleanup
    bamChunker_destruct(bamChunker);
    fclose(referenceOutFh);
    fclose(bamOutFh);
    stHash_destruct(referenceSequences);
    polishParams_destruct(params);

    //while(1); // Use this for testing for memory leaks

    return 0;
}

