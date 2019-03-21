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
#include <omp.h>
#include <time.h>

#include "marginVersion.h"
#include "margin.h"
#include "htsIntegration.h"


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
    fprintf(stderr, "    -h --help                : Print this help screen\n");
    fprintf(stderr, "    -a --logLevel            : Set the log level [default = info]\n");
    # ifdef _OPENMP
    fprintf(stderr, "    -t --threads             : Set number of concurrent threads (default: 1)\n");
    #endif
    fprintf(stderr, "    -o --outputBase          : Name to use for output files [default = output]\n");
    fprintf(stderr, "    -r --region              : If set, will only compute for given chromosomal region.\n");
    fprintf(stderr, "                                 Format: chr:start_pos-end_pos (chr3:2000-3000).\n");

    fprintf(stderr, "    -i --outputRepeatCounts  : File to write out the repeat counts [default = NULL]\n");
    fprintf(stderr, "    -j --outputPoaTsv        : File to write out the poa as TSV file [default = NULL]\n");
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
    int numThreads = 0;
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
                # ifdef _OPENMP
                { "threads", required_argument, 0, 't'},
                #endif
                { "outputBase", required_argument, 0, 'o'},
                { "region", required_argument, 0, 'r'},
                { "verbose", required_argument, 0, 'v'},
				{ "outputRepeatCounts", required_argument, 0, 'i'},
				{ "outputPoaTsv", required_argument, 0, 'j'},
                { 0, 0, 0, 0 } };

        int option_index = 0;
        int key = getopt_long(argc-2, &argv[2], "a:o:v:r:hi:j:t:", long_options, &option_index);

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
        case 'i':
        	outputRepeatCountFile = stString_copy(optarg);
        	break;
        case 'j':
        	outputPoaTsvFile = stString_copy(optarg);
        	break;
        case 't':
            numThreads = atoi(optarg);
            if (numThreads <= 0) {
                st_errAbort("Invalid thread count: %d", numThreads);
            }
            break;
        default:
            usage();
            return 0;
        }
    }

    // Initialization from arguments
    st_setLogLevelFromString(logLevelString);
    free(logLevelString);
    # ifdef _OPENMP
    if (numThreads > 0) {
        omp_set_num_threads(numThreads);
    } else {
        omp_set_num_threads(1);
    }
    st_logInfo("Running OpenMP with %d threads.\n", omp_get_max_threads());
    # endif

    // Parse parameters
    st_logInfo("> Parsing model parameters from file: %s\n", paramsFile);
    Params *params = params_readParams(paramsFile);

    // Print a report of the parsed parameters
    if(st_getLogLevel() == debug) {
    	params_printParameters(params, stderr);
    }

    // Parse reference as map of header string to nucleotide sequences
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

    // Open output files
    char *polishedReferenceOutFile = stString_print("%s.fa", outputBase);
    st_logInfo("> Going to write polished reference in : %s\n", polishedReferenceOutFile);
    FILE *polishedReferenceOutFh = fopen(polishedReferenceOutFile, "w");
    free(polishedReferenceOutFile);

    // Open optional files for writing out repeat counts and such
    FILE *outputRepeatCountFileHandle = NULL;
    if(outputRepeatCountFile != NULL) {
    	outputRepeatCountFileHandle = fopen(outputRepeatCountFile, "w");
    }
    FILE *outputPoaTsvFileHandle = NULL;
    if(outputPoaTsvFile != NULL) {
    	outputPoaTsvFileHandle = fopen(outputPoaTsvFile, "w");
    }

    // if regionStr is NULL, it will be ignored in construct2
    BamChunker *bamChunker = bamChunker_construct2(bamInFile, regionStr, params->polishParams);
    st_logInfo("> Set up bam chunker with chunk size: %i and overlap %i (for region=%s)\n",
    		   (int)bamChunker->chunkSize, (int)bamChunker->chunkBoundary, regionStr == NULL ? "all" : regionStr);

    // Polish chunks
    // Each chunk produces a char* as output which is saved here
    char **chunkResults = st_calloc(bamChunker->chunkCount, sizeof(char*));

    // multiproccess the chunks, save to results
    int64_t chunkIdx;
    #pragma omp parallel for
    for (chunkIdx = 0; chunkIdx < bamChunker->chunkCount; chunkIdx++) {
        // Time all chunks
        clock_t start = clock();

        // Get chunk
        BamChunk *bamChunk = bamChunker_getChunk(bamChunker, chunkIdx);
        char *logIdentifier;
        # ifdef _OPENMP
        logIdentifier = stString_print("T%02d_C%05"PRId64, omp_get_thread_num(), chunkIdx);
        # else
        logIdentifier = stString_copy("");
        # endif

        // Get reference string for chunk of alignment
        //TODO maybe we don't need to stringcopy this
        char *fullReferenceString = stString_copy(stHash_search(referenceSequences, bamChunk->refSeqName));
        if (fullReferenceString == NULL) {
            st_logCritical("> ERROR: Reference sequence missing from reference map: %s \n", bamChunk->refSeqName);
            continue;
        }
        int64_t refLen = strlen(fullReferenceString);
        char *referenceString = stString_getSubString(fullReferenceString, bamChunk->chunkBoundaryStart,
                                                      (refLen < bamChunk->chunkBoundaryEnd ? refLen
                                                                                           : bamChunk->chunkBoundaryEnd) -
                                                      bamChunk->chunkBoundaryStart);

        st_logInfo("> %s: Going to process a chunk for reference sequence: %s, starting at: %i and ending at: %i\n",
                   logIdentifier, bamChunk->refSeqName, (int) bamChunk->chunkBoundaryStart,
                   (int) (refLen < bamChunk->chunkBoundaryEnd ? refLen : bamChunk->chunkBoundaryEnd));

        // Convert bam lines into corresponding reads and alignments
        st_logInfo("> %s Parsing input reads from file: %s\n", logIdentifier, bamInFile);
        stList *reads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
        stList *alignments = stList_construct3(0, (void (*)(void *)) stList_destruct);
        convertToReadsAndAlignments(bamChunk, reads, alignments);

        Poa *poa = NULL; // The poa alignment
        char *polishedReferenceString = NULL; // The polished reference string

        // Now run the polishing method

        if (params->polishParams->useRunLengthEncoding) {
            st_logInfo("> %s Running polishing algorithm using run-length encoding\n", logIdentifier);

            // Run-length encoded polishing

            // Do run length encoding (RLE)

            // First RLE the reference
            RleString *rleReference = rleString_construct(referenceString);

            // Now RLE the reads
            stList *rleNucleotides = stList_construct3(0, (void (*)(void *)) rleString_destruct);
            stList *rleReads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
            stList *rleAlignments = stList_construct3(0, (void (*)(void *)) stList_destruct);
            for (int64_t j = 0; j < stList_length(reads); j++) {
                BamChunkRead *read = stList_get(reads, j);
                RleString *rleNucleotideString = rleString_construct(read->nucleotides);
                stList_append(rleNucleotides, rleNucleotideString);
                stList_append(rleReads, bamChunkRead_constructRLECopy(read, rleNucleotideString));
                stList_append(rleAlignments,
                              runLengthEncodeAlignment(stList_get(alignments, j), rleReference, rleNucleotideString));
            }

            // Generate partial order alignment (POA) (destroys rleAlignments in the process)
            poa = poa_realignAll(rleReads, rleAlignments, rleReference->rleString, params->polishParams);

            // Now optionally do phasing and haplotype specific polishing

            //stList *anchorAlignments = poa_getAnchorAlignments(poa, NULL, stList_length(reads), params->polishParams);
            //stList *reads1, *reads2;
            //phaseReads(poa->refString, stList_length(poa->nodes)-1, l, anchorAlignments, &reads1, &reads2, params);

            // Do run-length decoding
            RleString *polishedRLEReference = expandRLEConsensus(poa, rleNucleotides, rleReads,
                                                                 params->polishParams->repeatSubMatrix);
            polishedReferenceString = rleString_expand(polishedRLEReference);

            // Log info about the POA
            if (st_getLogLevel() >= info) {
                st_logInfo("> %s Summary stats for POA:\t", logIdentifier);
                poa_printSummaryStats(poa, stderr);
            }
            if (st_getLogLevel() >= debug) {
                poa_print(poa, stderr, rleReads, 5, 5);
            }

            // Write any optional outputs about repeat count and POA, etc.
            /*TODO fix this after multithreading
            if(outputPoaTsvFileHandle != NULL) {
                poa_printTSV(poa, outputPoaTsvFileHandle, rleReads, 5, 0);
            }
            if(outputRepeatCountFileHandle != NULL) {
                poa_printRepeatCounts(poa, outputRepeatCountFileHandle, rleNucleotides, rleReads);
            } */

            // Now cleanup run-length stuff
            stList_destruct(rleNucleotides);
            stList_destruct(rleReads);
            stList_destruct(rleAlignments);
            rleString_destruct(rleReference);
            rleString_destruct(polishedRLEReference);
        } else { // Non-run-length encoded polishing
            st_logInfo("> %s Running polishing algorithm without using run-length encoding\n", logIdentifier);

            // Generate partial order alignment (POA)
            poa = poa_realignAll(reads, alignments, referenceString, params->polishParams);

            // Polished string is the final backbone of the POA
            polishedReferenceString = stString_copy(poa->refString);

            // Log info about the POA
            if (st_getLogLevel() >= info) {
                st_logInfo("> %s Summary stats for POA:\t", logIdentifier);
                poa_printSummaryStats(poa, stderr);
            }
            if (st_getLogLevel() >= debug) {
                poa_print(poa, stderr, reads, 5, 5);
            }
        }

        // save polished reference string to chunk output array
        chunkResults[chunkIdx] = polishedReferenceString;

        // report timing
        clock_t end = clock();
        st_logInfo("> %s: Chunk with %d reads processed in %d sec\n",
                   logIdentifier, stList_length(reads), (int) (end - start) / CLOCKS_PER_SEC);

        // cleanup
        poa_destruct(poa);
        stList_destruct(reads);
        stList_destruct(alignments);
        free(referenceString);
        free(fullReferenceString);
    }

    // merge chunks
    stList *polishedReferenceStrings = NULL; // The polished reference strings, one for each chunk
    char *referenceSequenceName = NULL;
    for (chunkIdx = 0; chunkIdx < bamChunker->chunkCount; chunkIdx++) {
        // Get chunk and polished
        BamChunk *bamChunk = bamChunker_getChunk(bamChunker, chunkIdx);
        char* polishedReferenceString = chunkResults[chunkIdx];

		// If there is no prior chunk for this contig
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
    fclose(polishedReferenceOutFh);

    // Cleanup
    st_logInfo("> Finished polishing.\n");
    bamChunker_destruct(bamChunker);
    stHash_destruct(referenceSequences);
    params_destruct(params);

    if(outputPoaTsvFileHandle != NULL) {
    	fclose(outputPoaTsvFileHandle);
    	free(outputPoaTsvFile);
    }
    if(outputRepeatCountFileHandle != NULL) {
    	fclose(outputRepeatCountFileHandle);
    	free(outputRepeatCountFile);
    }

    //while(1); // Use this for testing for memory leaks

    return 0;
}

