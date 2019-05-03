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
#include "helenFeatures.h"


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
    fprintf(stderr, "    -t --threads             : Set number of concurrent threads [default = 1]\n");
    #endif
    fprintf(stderr, "    -o --outputBase          : Name to use for output files [default = 'output']\n");
    fprintf(stderr, "    -r --region              : If set, will only compute for given chromosomal region.\n");
    fprintf(stderr, "                                 Format: chr:start_pos-end_pos (chr3:2000-3000).\n");

    fprintf(stderr, "\nHELEN feature generation options:\n");
    fprintf(stderr, "    -f --outputFeatureType   : output features of chunks for HELEN.  Valid types:\n");
    fprintf(stderr, "                                 simpleWeight:    weighted likelihood from POA nodes (non-RLE)\n");
    fprintf(stderr, "                                 rleWeight:       weighted likelihood from POA nodes (RLE)\n");
    fprintf(stderr, "                                 nuclAndRlWeight: weighted likelihood, split into nucleotide and run length\n");
    fprintf(stderr, "                                 splitRleWeight:  weighted likelihood, with run lengths split into chunks\n");
    fprintf(stderr, "    -L --splitRleWeightMaxRL : max run length (for 'splitRleWeight' type only) [default = %d]\n", POAFEATURE_SPLIT_MAX_RUN_LENGTH_DEFAULT);
    fprintf(stderr, "    -u --trueReferenceBam    : true reference aligned to ASSEMBLY_FASTA, for HELEN\n");
    fprintf(stderr, "                               features.  Setting this parameter will include labels\n");
    fprintf(stderr, "                               in output.\n");
    fprintf(stderr, "    -5 --hdf5Only            : only output H5 feature files.  Default behavior is to output\n");
    fprintf(stderr, "                               h5, tsv, and fa for each chunk.\n");

    fprintf(stderr, "\nMiscellaneous supplementary output options:\n");
    fprintf(stderr, "    -i --outputRepeatCounts  : Output base to write out the repeat counts [default = NULL]\n");
    fprintf(stderr, "    -j --outputPoaTsv        : Output base to write out the poa as TSV file [default = NULL]\n");
    fprintf(stderr, "\n");
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
    int numThreads = 1;
    char *outputRepeatCountBase = NULL;
    char *outputPoaTsvBase = NULL;

    // for feature generation
    HelenFeatureType helenFeatureType = HFEAT_NONE;
    char *trueReferenceBam = NULL;
    BamChunker *trueReferenceChunker = NULL;
    bool fullFeatureOutput = TRUE;
    int64_t splitWeightMaxRunLength = POAFEATURE_SPLIT_MAX_RUN_LENGTH_DEFAULT;
    void **splitWeightHDF5Files = NULL;

    // TODO: When done testing, optionally set random seed using st_randomSeed();

    if(argc < 4) {
        free(outputBase);
        free(logLevelString);
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
                { "outputFeatureType", required_argument, 0, 'f'},
                { "trueReferenceBam", required_argument, 0, 'u'},
                { "hdf5Only", no_argument, 0, '5'},
                { "splitRleWeightMaxRL", required_argument, 0, 'L'},
				{ "outputRepeatCounts", required_argument, 0, 'i'},
				{ "outputPoaTsv", required_argument, 0, 'j'},
                { 0, 0, 0, 0 } };

        int option_index = 0;
        int key = getopt_long(argc-2, &argv[2], "a:o:v:r:f:u:h5L:i:j:t:", long_options, &option_index);

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
        	outputRepeatCountBase = stString_copy(optarg);
        	break;
        case 'j':
            outputPoaTsvBase = stString_copy(optarg);
            break;
        case 'f':
            if (stString_eq(optarg, "simpleWeight")) {
                helenFeatureType = HFEAT_SIMPLE_WEIGHT;
            } else if (stString_eq(optarg, "rleWeight")) {
                helenFeatureType = HFEAT_RLE_WEIGHT;
            } else if (stString_eq(optarg, "nuclAndRlWeight")) {
                helenFeatureType = HFEAT_NUCL_AND_RL_WEIGHT;
            } else if (stString_eq(optarg, "splitRleWeight")) {
                helenFeatureType = HFEAT_SPLIT_RLE_WEIGHT;
            } else {
                fprintf(stderr, "Unrecognized featureType for HELEN: %s\n\n", optarg);
                usage();
                return 1;
            }
            break;
        case 'u':
            trueReferenceBam = stString_copy(optarg);
            break;
        case '5':
            fullFeatureOutput = FALSE;
            break;
        case 'L':
            splitWeightMaxRunLength = atoi(optarg);
            if (splitWeightMaxRunLength <= 0) {
                st_errAbort("Invalid splitRleWeightMaxRL: %d", splitWeightMaxRunLength);
            }
            break;
        case 't':
            numThreads = atoi(optarg);
            if (numThreads <= 0) {
                st_errAbort("Invalid thread count: %d", numThreads);
            }
            break;
        default:
            usage();
            free(outputBase);
            free(logLevelString);
            free(bamInFile);
            free(referenceFastaFile);
            free(paramsFile);
            if (trueReferenceBam != NULL) free(trueReferenceBam);
            return 0;
        }
    }
    // sanity check (verify files exist)
    if (access(bamInFile, R_OK ) != 0) {
        st_errAbort("Could not read from file: %s\n", bamInFile);
        char *idx = stString_print("%s.bai", bamInFile);
        if (access(idx, R_OK ) != 0 ) {
            st_errAbort("BAM does not appear to be indexed: %s\n", bamInFile);
        }
        free(idx);
    } else if (access(referenceFastaFile, R_OK ) != 0 ) {
        st_errAbort("Could not read from file: %s\n", referenceFastaFile);
    } else if (access(paramsFile, R_OK ) != 0 ) {
        st_errAbort("Could not read from file: %s\n", paramsFile);
    } else if (trueReferenceBam != NULL && access(trueReferenceBam, R_OK ) != 0 ) {
        st_errAbort("Could not read from file: %s\n", trueReferenceBam);
        char *idx = stString_print("%s.bai", trueReferenceBam);
        if (access(idx, R_OK ) != 0 ) {
            st_errAbort("BAM does not appear to be indexed: %s\n", trueReferenceBam);
        }
        free(idx);
    }

    // Initialization from arguments
    st_setLogLevelFromString(logLevelString);
    free(logLevelString);
    # ifdef _OPENMP
    if (numThreads <= 0) {
        numThreads = 1;
    }
    omp_set_num_threads(numThreads);
    st_logInfo("Running OpenMP with %d threads.\n", omp_get_max_threads());
    # endif

    // Parse parameters
    st_logInfo("> Parsing model parameters from file: %s\n", paramsFile);
    Params *params = params_readParams(paramsFile);

    // Set no RLE if appropriate feature type is set
    if (helenFeatureType == HFEAT_SIMPLE_WEIGHT) {
        if (params->polishParams->useRunLengthEncoding) {
            st_logInfo("> Changing runLengthEncoding parameter to FALSE because of HELEN feature type.\n");
            params->polishParams->useRunLengthEncoding = FALSE;
        }
    // everthing else requires RLE
    } else if (helenFeatureType != HFEAT_NONE) {
        if (!params->polishParams->useRunLengthEncoding) {
            st_logInfo("> Changing runLengthEncoding parameter to TRUE because of HELEN feature type.\n");
            params->polishParams->useRunLengthEncoding = TRUE;
        }
    }

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

    // get chunker for bam.  if regionStr is NULL, it will be ignored
    BamChunker *bamChunker = bamChunker_construct2(bamInFile, regionStr, params->polishParams);
    st_logInfo("> Set up bam chunker with chunk size %i and overlap %i (for region=%s), resulting in %i total chunks\n",
    		   (int)bamChunker->chunkSize, (int)bamChunker->chunkBoundary, regionStr == NULL ? "all" : regionStr,
    		   bamChunker->chunkCount);

    // for feature generation
    BamChunker *trueReferenceBamChunker = NULL;
    if (trueReferenceBam != NULL) {
        trueReferenceBamChunker = bamChunker_copyConstruct(bamChunker);
        free(trueReferenceBamChunker->bamFile);
        trueReferenceBamChunker->bamFile = stString_copy(trueReferenceBam);
    }
    #ifdef _HDF5
    if (helenFeatureType == HFEAT_SPLIT_RLE_WEIGHT) {
        splitWeightHDF5Files = (void**) openSplitRleFeatureHDF5FilesByThreadCount(outputBase, numThreads);
    }
    #endif


    // Polish chunks
    // Each chunk produces a char* as output which is saved here
    char **chunkResults = st_calloc(bamChunker->chunkCount, sizeof(char*));

    // multiproccess the chunks, save to results
    int64_t chunkIdx;
    #pragma omp parallel for schedule(dynamic,1)
    for (chunkIdx = 0; chunkIdx < bamChunker->chunkCount; chunkIdx++) {
        // Time all chunks
        time_t start = time(NULL);

        // Get chunk
        BamChunk *bamChunk = bamChunker_getChunk(bamChunker, chunkIdx);
        char *logIdentifier;
        # ifdef _OPENMP
        logIdentifier = stString_print(" T%02d_C%05"PRId64, omp_get_thread_num(), chunkIdx);
        # else
        logIdentifier = stString_copy("");
        # endif

        // Get reference string for chunk of alignment
        char *fullReferenceString = stHash_search(referenceSequences, bamChunk->refSeqName);
        if (fullReferenceString == NULL) {
            st_logCritical("> ERROR: Reference sequence missing from reference map: %s \n", bamChunk->refSeqName);
            continue;
        }
        int64_t fullRefLen = strlen(fullReferenceString);
        char *referenceString = stString_getSubString(fullReferenceString, bamChunk->chunkBoundaryStart,
                                                      (fullRefLen < bamChunk->chunkBoundaryEnd ? fullRefLen
                                                                                           : bamChunk->chunkBoundaryEnd) -
                                                      bamChunk->chunkBoundaryStart);


        st_logInfo(">%s Going to process a chunk for reference sequence: %s, starting at: %i and ending at: %i\n",
                   logIdentifier, bamChunk->refSeqName, (int) bamChunk->chunkBoundaryStart,
                   (int) (fullRefLen < bamChunk->chunkBoundaryEnd ? fullRefLen : bamChunk->chunkBoundaryEnd));

        // Convert bam lines into corresponding reads and alignments
        st_logInfo(">%s Parsing input reads from file: %s\n", logIdentifier, bamInFile);
        stList *reads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
        stList *alignments = stList_construct3(0, (void (*)(void *)) stList_destruct);
        convertToReadsAndAlignments(bamChunk, reads, alignments);

        // do downsampling if appropriate
        if (params->polishParams->maxDepth > 0) {
            // get downsampling structures
            stList *filteredReads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
            stList *discardedReads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
            stList *filteredAlignments = stList_construct3(0, (void (*)(void *)) stList_destruct);
            stList *discardedAlignments = stList_construct3(0, (void (*)(void *)) stList_destruct);

            bool didDownsample = poorMansDownsample(params->polishParams->maxDepth, bamChunk, reads, alignments,
                    filteredReads, filteredAlignments, discardedReads, discardedAlignments);

            // we need to destroy the discarded reads and structures
            if (didDownsample) {
                st_logInfo(" %s Downsampled from %"PRId64" to %"PRId64" reads\n", logIdentifier,
                        stList_length(reads), stList_length(filteredReads));
                // free all reads and alignments not used
                stList_destruct(discardedReads);
                stList_destruct(discardedAlignments);
                // still has all the old reads, need to not free these
                stList_setDestructor(reads, NULL);
                stList_setDestructor(alignments, NULL);
                stList_destruct(reads);
                stList_destruct(alignments);
                // and keep the filtered reads
                reads = filteredReads;
                alignments = filteredAlignments;
            }
            // no downsampling, we just need to free the (empty) objects
            else {
                stList_destruct(filteredReads);
                stList_destruct(filteredAlignments);
                stList_destruct(discardedReads);
                stList_destruct(discardedAlignments);
            }

        }

        Poa *poa = NULL; // The poa alignment
        char *polishedConsensusString = NULL; // The polished reference string


        // prep for RLE work
        RleString *rleReference = NULL;
        stList *rleNucleotides = stList_construct3(0, (void (*)(void *)) rleString_destruct);
        stList *rleReads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
        stList *rleAlignments = stList_construct3(0, (void (*)(void *)) stList_destruct);
        uint64_t totalNucleotides = 0;

        // Note RLE status (and handle reference)
        if (params->polishParams->useRunLengthEncoding) {
            st_logInfo(">%s Applying RLE\n", logIdentifier);
            rleReference = rleString_construct(referenceString);
        } else {
            st_logInfo(">%s Skipping RLE\n", logIdentifier);
            rleReference = rleString_constructNoRLE(referenceString);
        }

        // RLE the reads
        for (int64_t j = 0; j < stList_length(reads); j++) {
            BamChunkRead *read = stList_get(reads, j);
            stList *alignment = stList_get(alignments, j);
            RleString *rleNucleotideString = NULL;

            // Perform or skip RLE
            if (params->polishParams->useRunLengthEncoding) {
                rleNucleotideString = rleString_construct(read->nucleotides);
            } else {
                rleNucleotideString = rleString_constructNoRLE(read->nucleotides);
            }
            totalNucleotides += rleNucleotideString->length;

            // Do RLE follow up regardless of whether RLE is applied
            stList_append(rleNucleotides, rleNucleotideString);
            stList_append(rleReads, bamChunkRead_constructRLECopy(read, rleNucleotideString));
            stList_append(rleAlignments, runLengthEncodeAlignment(alignment, rleReference, rleNucleotideString));
        }


        // Run the polishing method
        st_logInfo(">%s Running polishing algorithm with %"PRId64" reads and %"PRIu64"K nucleotides\n",
                logIdentifier, stList_length(reads), totalNucleotides >> 10);

        // Generate partial order alignment (POA) (destroys rleAlignments in the process)
        poa = poa_realignAll(rleReads, rleAlignments, rleReference->rleString, params->polishParams);

        // Now optionally do phasing and haplotype specific polishing

        /*TODO needs to be implemented
        stList *anchorAlignments = poa_getAnchorAlignments(poa, NULL, stList_length(reads), params->polishParams);
        stList *reads1, *reads2;
        phaseReads(poa->refString, stList_length(poa->nodes)-1, l, anchorAlignments, &reads1, &reads2, params);
        */

        // get polished reference string and expand RLE (regardless of whether RLE was applied)
        RleString *polishedRleConsensus = expandRLEConsensus(poa, rleNucleotides, rleReads,
                                                             params->polishParams->repeatSubMatrix);
        polishedConsensusString = rleString_expand(polishedRleConsensus);

        // Log info about the POA
        if (st_getLogLevel() >= info) {
            st_logInfo(">%s Summary stats for POA:\t", logIdentifier);
            poa_printSummaryStats(poa, stderr);
        }
        if (st_getLogLevel() >= debug) {
            poa_print(poa, stderr, rleReads, 5, 5);
        }

        // Write any optional outputs about repeat count and POA, etc.
        if(outputPoaTsvBase != NULL) {
            char *outputPoaTsvFilename = stString_print("%s.poa.C%05"PRId64".%s-%"PRId64"-%"PRId64".tsv",
                                                        outputPoaTsvBase, chunkIdx, bamChunk->refSeqName,
                                                        bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
            FILE *outputPoaTsvFileHandle = fopen(outputPoaTsvFilename, "w");
            poa_printTSV(poa, outputPoaTsvFileHandle, rleReads, 5, 0);
            fclose(outputPoaTsvFileHandle);
            free(outputPoaTsvFilename);
        }
        if(outputRepeatCountBase != NULL) {
            char *outputRepeatCountFilename = stString_print("%s.repeatCount.C%05"PRId64".%s-%"PRId64"-%"PRId64".tsv",
                                                             outputRepeatCountBase, chunkIdx, bamChunk->refSeqName,
                                                             bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
            FILE *outputRepeatCountFileHandle = fopen(outputRepeatCountFilename, "w");
            poa_printRepeatCounts(poa, outputRepeatCountFileHandle, rleNucleotides, rleReads);
            fclose(outputRepeatCountFileHandle);
            free(outputRepeatCountFilename);
        }


        // save polished reference string to chunk output array
        chunkResults[chunkIdx] = polishedConsensusString;

        // HELEN feature outputs
        if (helenFeatureType != HFEAT_NONE) {
            st_logInfo(">%s Performing feature generation for chunk.\n", logIdentifier);
            // get filename
            char *helenFeatureOutfileBase = NULL;
            switch (helenFeatureType) {
                case HFEAT_SIMPLE_WEIGHT:
                    helenFeatureOutfileBase = stString_print("%s.simpleWeight.C%05"PRId64".%s-%"PRId64"-%"PRId64,
                                                             outputBase, chunkIdx, bamChunk->refSeqName,
                                                             bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
                    break;
                case HFEAT_RLE_WEIGHT:
                    helenFeatureOutfileBase = stString_print("%s.rleWeight.C%05"PRId64".%s-%"PRId64"-%"PRId64,
                                                             outputBase, chunkIdx, bamChunk->refSeqName,
                                                             bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
                    break;
                case HFEAT_NUCL_AND_RL_WEIGHT:
                    helenFeatureOutfileBase = stString_print("%s.nuclAndRlWeight.C%05"PRId64".%s-%"PRId64"-%"PRId64,
                                                             outputBase, chunkIdx, bamChunk->refSeqName,
                                                             bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
                    break;
                case HFEAT_SPLIT_RLE_WEIGHT:
                    // name of folder, not of file
                    helenFeatureOutfileBase = stString_print("splitRleWeight.C%05"PRId64".%s-%"PRId64"-%"PRId64,
                                                             chunkIdx, bamChunk->refSeqName,
                                                             bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
                    break;
                default:
                    st_errAbort("Unhandled HELEN feature type!\n");
            }

            // necessary to annotate poa with truth (if true reference BAM has been specified)
            stList *trueRefAlignment = NULL;
            RleString *trueRefRleString = NULL;
            bool validReferenceAlignment = FALSE;

            // get reference chunk
            if (trueReferenceBam != NULL) {
                // get alignment of true ref to assembly
                stList *trueRefReads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
                stList *unused = stList_construct3(0, (void (*)(void *)) stList_destruct);
                // construct new chunk
                BamChunk *trueRefBamChunk = bamChunk_copyConstruct(bamChunk);
                trueRefBamChunk->parent = trueReferenceBamChunker;
                // get true ref as "read"
                uint32_t trueAlignmentCount = convertToReadsAndAlignments(trueRefBamChunk, trueRefReads, unused);

                // poor man's "do we have a unique alignment"
                if (trueAlignmentCount == 1) {
                    BamChunkRead *trueRefRead = stList_get(trueRefReads, 0);

                    stList *trueRefAlignmentRawSpace = alignConsensusAndTruth(polishedConsensusString, trueRefRead->nucleotides);
                    if (st_getLogLevel() == debug) {
                        printMEAAlignment(polishedConsensusString, trueRefRead->nucleotides,
                                          strlen(polishedConsensusString), strlen(trueRefRead->nucleotides),
                                          trueRefAlignmentRawSpace, NULL, NULL);
                    }


                    // convert to rleSpace if appropriate
                    if (params->polishParams->useRunLengthEncoding) {
                        trueRefRleString = rleString_construct(trueRefRead->nucleotides);
                        trueRefAlignment = runLengthEncodeAlignment2(trueRefAlignmentRawSpace, polishedRleConsensus,
                                trueRefRleString, 1, 2, 0);
                        if (st_getLogLevel() == debug) {
                            printMEAAlignment(polishedRleConsensus->rleString, trueRefRleString->rleString,
                                              strlen(polishedRleConsensus->rleString),
                                              strlen(trueRefRleString->rleString),
                                              trueRefAlignment, polishedRleConsensus->repeatCounts,
                                              trueRefRleString->repeatCounts);
                        }
                        stList_destruct(trueRefAlignmentRawSpace);
                    } else {
                        trueRefRleString = rleString_constructNoRLE(trueRefRead->nucleotides);
                        trueRefAlignment = trueRefAlignmentRawSpace;
                    }


                    // we found a single alignment of reference
                    double refLengthRatio = 1.0 * trueRefRleString->length / polishedRleConsensus->length;
                    double alnLengthRatio = 1.0 * stList_length(trueRefAlignment) / polishedRleConsensus->length;
                    int refLengthRatioHundredthsOffOne = abs((int) (100 * (1.0 - refLengthRatio)));
                    int alnLengthRatioHundredthsOffOne = abs((int) (100 * (1.0 - alnLengthRatio)));
                    if (stList_length(trueRefAlignment) > 0 && refLengthRatioHundredthsOffOne < 10 &&
                            alnLengthRatioHundredthsOffOne < 10) {
                        validReferenceAlignment = TRUE;
                    } else {
                        st_logInfo(" %s True reference alignment QC failed:  polished length %"PRId64", true ref length"
                                   " ratio (true/polished) %f, aligned pairs length ratio (true/polished): %f\n",
                                   logIdentifier, polishedRleConsensus->length, refLengthRatio, alnLengthRatio);
                    }
                }

                stList_destruct(trueRefReads);
                stList_destruct(unused);
                bamChunk_destruct(trueRefBamChunk);
            }

            // either write it, or note that we failed to find a valid reference alignment
            if (trueReferenceBam != NULL && !validReferenceAlignment) {
                st_logInfo(" %s No valid reference alignment was found, skipping HELEN feature output.\n", logIdentifier);
            } else {
                st_logInfo(" %s Writing HELEN features with filename base: %s\n", logIdentifier, helenFeatureOutfileBase);

                // write the actual features (type dependent)
                poa_writeHelenFeatures(helenFeatureType, poa, rleReads, rleNucleotides, helenFeatureOutfileBase,
                        bamChunk, trueRefAlignment, polishedRleConsensus, trueRefRleString, fullFeatureOutput,
                        splitWeightMaxRunLength, splitWeightHDF5Files);

                // write the polished chunk in fasta format
                if (fullFeatureOutput) {
                    char *chunkPolishedRefFilename = stString_print("%s.fa", helenFeatureOutfileBase);
                    char *chunkPolishedRefContigName = stString_print("%s\t%"PRId64"\t%"PRId64"\t%s",
                                                                      bamChunk->refSeqName,
                                                                      bamChunk->chunkBoundaryStart,
                                                                      bamChunk->chunkBoundaryEnd,
                                                                      helenFeatureOutfileBase);
                    FILE *chunkPolishedRefOutFh = fopen(chunkPolishedRefFilename, "w");
                    fastaWrite(polishedConsensusString, chunkPolishedRefContigName, chunkPolishedRefOutFh);
                    fclose(chunkPolishedRefOutFh);
                    free(chunkPolishedRefFilename);
                    free(chunkPolishedRefContigName);
                }
            }

            // cleanup
            free(helenFeatureOutfileBase);
            if (trueRefAlignment != NULL) stList_destruct(trueRefAlignment);
            if (trueRefRleString != NULL) rleString_destruct(trueRefRleString);
        }

        // report timing
        st_logInfo(">%s Chunk with %"PRId64" reads and %"PRIu64"K nucleotides processed in %d sec\n",
                   logIdentifier, stList_length(reads), totalNucleotides >> 10, (int) (time(NULL) - start));

        // Cleanup
        stList_destruct(rleNucleotides);
        stList_destruct(rleReads);
        stList_destruct(rleAlignments);
        rleString_destruct(rleReference);
        rleString_destruct(polishedRleConsensus);
        poa_destruct(poa);
        stList_destruct(reads);
        stList_destruct(alignments);
        free(referenceString);
        free(logIdentifier);
    }

    // merge chunks
    st_logInfo("> Merging polished reference strings from %"PRIu64" chunks.\n", bamChunker->chunkCount);
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

			st_logInfo("  Removed overlap between neighbouring chunks. Approx overlap size: %i, overlap-match weight: %f, "
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

    if (trueReferenceBam != NULL) free(trueReferenceBam);
    if (trueReferenceBamChunker != NULL) bamChunker_destruct(trueReferenceBamChunker);

    if (regionStr != NULL) free(regionStr);
    #ifdef _HDF5
    if (splitWeightHDF5Files != NULL) {
        for (int64_t i = 0; i < numThreads; i++) {
            splitRleFeatureHDF5FileInfo_destruct((SplitRleFeatureHDF5FileInfo*) splitWeightHDF5Files[i]);
        }
    }
    #endif
    free(chunkResults);
    free(outputBase);
    free(bamInFile);
    free(referenceFastaFile);
    free(paramsFile);


//    while(1); // Use this for testing for memory leaks

    return 0;
}

