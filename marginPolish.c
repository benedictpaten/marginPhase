a new line
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
#include <time.h>
#include <sys/stat.h>
#include "marginVersion.h"

#include "margin.h"
#include "htsIntegration.h"
#include "helenFeatures.h"


//TODO move these to a better spot
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

    RleString *rleRef = params->polishParams->useRunLengthEncoding ?
            rleString_construct(referenceString) : rleString_construct_no_rle(referenceString);
    free(referenceString);

    return rleRef;
}


/*
 * Main functions
 */

void usage() {
    fprintf(stderr, "usage: marginPolish <BAM_FILE> <ASSEMBLY_FASTA> <PARAMS> [options]\n");
    fprintf(stderr, "Version: %s \n\n", MARGIN_POLISH_VERSION_H);
    fprintf(stderr, "Polishes the ASSEMBLY_FASTA using alignments in BAM_FILE.\n");

    fprintf(stderr, "\nRequired arguments:\n");
    fprintf(stderr, "    BAM_FILE is the alignment of reads to the assembly (or reference).\n");
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

    # ifdef _HDF5
    fprintf(stderr, "\nHELEN feature generation options:\n");
    fprintf(stderr, "    -f --produceFeatures     : output features for HELEN.\n");
    fprintf(stderr, "    -F --featureType         : output features of chunks for HELEN.  Valid types:\n");
    fprintf(stderr, "                                 splitRleWeight:   [default] run lengths split into chunks\n");
    fprintf(stderr, "                                 channelRleWeight: run lengths split into per-nucleotide channels\n");
    fprintf(stderr, "                                 simpleWeight:     weighted likelihood from POA nodes (non-RLE)\n");
    fprintf(stderr, "    -L --splitRleWeightMaxRL : max run length (for 'splitRleWeight' and 'channelRleWeight' types) \n");
    fprintf(stderr, "                                 [splitRleWeight default = %d, channelRleWeight default = %d]\n",
            POAFEATURE_SPLIT_MAX_RUN_LENGTH_DEFAULT, POAFEATURE_CHANNEL_MAX_RUN_LENGTH_DEFAULT);
    fprintf(stderr, "    -u --trueReferenceBam    : true reference aligned to ASSEMBLY_FASTA, for HELEN\n");
    fprintf(stderr, "                               features.  Setting this parameter will include labels\n");
    fprintf(stderr, "                               in output.\n");
    # endif

    fprintf(stderr, "\nMiscellaneous supplementary output options:\n");
    fprintf(stderr, "    -i --outputRepeatCounts  : Output base to write out the repeat counts [default = NULL]\n");
    fprintf(stderr, "    -j --outputPoaTsv        : Output base to write out the poa as TSV file [default = NULL]\n");
    fprintf(stderr, "    -d --outputPoaDot        : Output base to write out the poa as DOT file [default = NULL]\n");
    fprintf(stderr, "\n");
}

char *getTimeDescriptorFromSeconds(int64_t seconds) {
    int64_t minutes = (int64_t) (seconds / 60);
    int64_t hours = (int64_t) (minutes / 60);
    char *timeDescriptor;

    if (hours > 0) {
        timeDescriptor = stString_print("%"PRId64"h %"PRId64"m", hours,
                                        minutes - (hours * 60));
    } else if (minutes > 0) {
        timeDescriptor = stString_print("%"PRId64"m %"PRId64"s", minutes,
                                        seconds - (minutes * 60));
    } else {
        timeDescriptor = stString_print("%"PRId64"s", seconds);
    }
    return timeDescriptor;
}

char *getFileBase(char *base, char *defawlt) {
    struct stat fileStat;
    int64_t rc = stat(base, &fileStat);
    if (S_ISDIR(fileStat.st_mode)) {
        if (optarg[strlen(base) - 1] == '/') optarg[strlen(base) - 1] = '\0';
        return stString_print("%s/%s", base, defawlt);
    } else {
        return stString_copy(base);
    }
}

int main(int argc, char *argv[]) {

    // Parameters / arguments
    char *logLevelString = stString_copy("critical");
    char *bamInFile = NULL;
    char *paramsFile = NULL;
    char *referenceFastaFile = NULL;
    char *outputBase = stString_copy("output");
    char *regionStr = NULL;
    int numThreads = 1;
    char *outputRepeatCountBase = NULL;
    char *outputPoaTsvBase = NULL;
    char *outputPoaDotBase = NULL;

    // for feature generation
    HelenFeatureType helenFeatureType = HFEAT_NONE;
    char *trueReferenceBam = NULL;
    BamChunker *trueReferenceChunker = NULL;
    bool fullFeatureOutput = FALSE;
    int64_t splitWeightMaxRunLength = 0;
    void **helenHDF5Files = NULL;

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
                { "produceFeatures", no_argument, 0, 'f'},
                { "featureType", required_argument, 0, 'F'},
                { "trueReferenceBam", required_argument, 0, 'u'},
                { "splitRleWeightMaxRL", required_argument, 0, 'L'},
				{ "outputRepeatCounts", required_argument, 0, 'i'},
				{ "outputPoaTsv", required_argument, 0, 'j'},
				{ "outputPoaDot", required_argument, 0, 'd'},
                { 0, 0, 0, 0 } };

        int option_index = 0;
        int key = getopt_long(argc-2, &argv[2], "a:o:v:r:fF:u:hL:i:j:d:t:", long_options, &option_index);

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
            outputBase = getFileBase(optarg, "output");
            break;
        case 'r':
            regionStr = stString_copy(optarg);
            break;
        case 'i':
            outputRepeatCountBase = getFileBase(optarg, "repeatCount");
            break;
        case 'j':
            outputPoaTsvBase = getFileBase(optarg, "poa");
            break;
        case 'd':
            outputPoaDotBase = getFileBase(optarg, "poa");
            break;
        case 'F':
            if (stString_eq(optarg, "simpleWeight")) {
                helenFeatureType = HFEAT_SIMPLE_WEIGHT;
            } else if (stString_eq(optarg, "splitRleWeight")) {
                helenFeatureType = HFEAT_SPLIT_RLE_WEIGHT;
            } else if (stString_eq(optarg, "channelRleWeight")) {
                helenFeatureType = HFEAT_CHANNEL_RLE_WEIGHT;
            } else {
                fprintf(stderr, "Unrecognized featureType for HELEN: %s\n\n", optarg);
                usage();
                return 1;
            }
            break;
        case 'u':
            trueReferenceBam = stString_copy(optarg);
            break;
        case 'f':
            if (helenFeatureType == HFEAT_NONE) helenFeatureType = HFEAT_SPLIT_RLE_WEIGHT;
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

    // sanitiy check for poa plotting
    if ((outputPoaTsvBase != NULL || outputPoaDotBase != NULL) && regionStr == NULL) {
        st_logCritical("--outputPoaTsv and --outputPoaDot options should only be used for a specific region!\n");
    }

    // Initialization from arguments
    time_t startTime = time(NULL);
    st_setLogLevelFromString(logLevelString);
    free(logLevelString);
    # ifdef _OPENMP
    if (numThreads <= 0) {
        numThreads = 1;
    }
    omp_set_num_threads(numThreads);
    st_logCritical("Running OpenMP with %d threads.\n", omp_get_max_threads());
    # endif
    if (helenFeatureType != HFEAT_NONE && splitWeightMaxRunLength == 0) {
        switch (helenFeatureType) {
            case HFEAT_SPLIT_RLE_WEIGHT:
                splitWeightMaxRunLength = POAFEATURE_SPLIT_MAX_RUN_LENGTH_DEFAULT;
                break;
            case HFEAT_CHANNEL_RLE_WEIGHT:
                splitWeightMaxRunLength = POAFEATURE_CHANNEL_MAX_RUN_LENGTH_DEFAULT;
                break;
            default:
                break;
        }
    }

    // Parse parameters
    st_logCritical("> Parsing model parameters from file: %s\n", paramsFile);
    Params *params = params_readParams(paramsFile);

    // Set no RLE if appropriate feature type is set
    if (helenFeatureType == HFEAT_SIMPLE_WEIGHT) {
        if (params->polishParams->useRunLengthEncoding) {
            st_errAbort("Invalid runLengthEncoding parameter because of HELEN feature type.\n");
        }
    // everthing else requires RLE
    } else if (helenFeatureType != HFEAT_NONE) {
        if (!params->polishParams->useRunLengthEncoding) {
            st_errAbort("Invalid runLengthEncoding parameter because of HELEN feature type.\n");
        }
    }

    // Print a report of the parsed parameters
    if(st_getLogLevel() == debug) {
    	params_printParameters(params, stderr);
    }

    // Parse reference as map of header string to nucleotide sequences
    st_logCritical("> Parsing reference sequences from file: %s\n", referenceFastaFile);
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
    st_logCritical("> Going to write polished reference in : %s\n", polishedReferenceOutFile);
    FILE *polishedReferenceOutFh = fopen(polishedReferenceOutFile, "w");
    if (polishedReferenceOutFh == NULL) {
        st_errAbort("Could not open %s for writing!\n", polishedReferenceOutFile);
    }
    free(polishedReferenceOutFile);

    // get chunker for bam.  if regionStr is NULL, it will be ignored
    BamChunker *bamChunker = bamChunker_construct2(bamInFile, regionStr, params->polishParams);
    st_logCritical("> Set up bam chunker with chunk size %i and overlap %i (for region=%s), resulting in %i total chunks\n",
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
    if (helenFeatureType != HFEAT_NONE) {
        helenHDF5Files = (void**) openHelenFeatureHDF5FilesByThreadCount(outputBase, numThreads);
    }
    #endif

    // Polish chunks
    // Each chunk produces a char* as output which is saved here
    char **chunkResults = st_calloc(bamChunker->chunkCount, sizeof(char*));

    // (may) need to shuffle chunks
    stList *chunkOrder = stList_construct3(0, (void (*)(void*))stIntTuple_destruct);
    for (int64_t i = 0; i < bamChunker->chunkCount; i++) {
        stList_append(chunkOrder, stIntTuple_construct1(i));
    }
    if (params->polishParams->shuffleChunks) {
        stList_shuffle(chunkOrder);
    }

    // multiproccess the chunks, save to results
    int64_t lastReportedPercentage = 0;
    time_t polishStartTime = time(NULL);

    # ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,1)
    # endif
    for (int64_t i = 0; i < bamChunker->chunkCount; i++) {
        int64_t chunkIdx = stIntTuple_get(stList_get(chunkOrder, i), 0);
        // Time all chunks
        time_t chunkStartTime = time(NULL);

        // Get chunk
        BamChunk *bamChunk = bamChunker_getChunk(bamChunker, chunkIdx);

        // logging
        char *logIdentifier;
        bool logProgress = FALSE;
        int64_t currentPercentage = (int64_t) (100 * i / bamChunker->chunkCount);
        # ifdef _OPENMP
        logIdentifier = stString_print(" T%02d_C%05"PRId64, omp_get_thread_num(), chunkIdx);
        if (omp_get_thread_num() == 0) {
            if (currentPercentage != lastReportedPercentage) {
                logProgress = TRUE;
                lastReportedPercentage = currentPercentage;
            }
        }
        # else
        logIdentifier = stString_copy("");
        if (currentPercentage != lastReportedPercentage) {
            logProgress = TRUE;
            lastReportedPercentage = currentPercentage;
        }
        # endif

        if (logProgress) {
            // log progress
            int64_t timeTaken = (int64_t) (time(NULL) - polishStartTime);
            int64_t secondsRemaining = (int64_t) floor(1.0 * timeTaken / currentPercentage * (100 - currentPercentage));
            char *timeDescriptor = (secondsRemaining == 0 && currentPercentage <= 50 ?
                    stString_print("unknown") : getTimeDescriptorFromSeconds(secondsRemaining));
            st_logCritical("> Polishing %2"PRId64"%% complete (%"PRId64"/%"PRId64").  Estimated time remaining: %s\n",
                    currentPercentage, i, bamChunker->chunkCount, timeDescriptor);
            free(timeDescriptor);
        }

        // Get reference string for chunk of alignment
        char *fullReferenceString = stHash_search(referenceSequences, bamChunk->refSeqName);
        if (fullReferenceString == NULL) {
            st_errAbort("ERROR: Reference sequence missing from reference map: %s. Perhaps the BAM and REF are mismatched?",
                    bamChunk->refSeqName);
        }
        int64_t fullRefLen = strlen(fullReferenceString);
        if (bamChunk->chunkBoundaryStart > fullRefLen) {
            st_errAbort("ERROR: Reference sequence %s has length %"PRId64", chunk %"PRId64" has start position %"
            PRId64". Perhaps the BAM and REF are mismatched?",
                    bamChunk->refSeqName, fullRefLen, chunkIdx, bamChunk->chunkBoundaryStart);
        }

        RleString *rleReference = bamChunk_getReferenceSubstring(bamChunk, referenceSequences, params);

        st_logInfo(">%s Going to process a chunk for reference sequence: %s, starting at: %i and ending at: %i\n",
                   logIdentifier, bamChunk->refSeqName, (int) bamChunk->chunkBoundaryStart,
                   (int) (fullRefLen < bamChunk->chunkBoundaryEnd ? fullRefLen : bamChunk->chunkBoundaryEnd));

        // Convert bam lines into corresponding reads and alignments
        st_logInfo(">%s Parsing input reads from file: %s\n", logIdentifier, bamInFile);
        stList *reads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
        stList *alignments = stList_construct3(0, (void (*)(void *)) stList_destruct);
        convertToReadsAndAlignments(bamChunk, rleReference, reads, alignments);

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

        //todo rle-ing moved to convertToReadsAndAlignmentsOrWhatever
        // prep for RLE work
//        stList *rleNucleotides = stList_construct3(0, (void (*)(void *)) rleString_destruct);
//        stList *rleReads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
//        stList *rleAlignments = stList_construct3(0, (void (*)(void *)) stList_destruct);
//        uint64_t totalNucleotides = 0;
//
//        // RLE the reads
//        for (int64_t j = 0; j < stList_length(reads); j++) {
//            BamChunkRead *read = stList_get(reads, j);
//            stList *alignment = stList_get(alignments, j);
//            RleString *rleNucleotideString = NULL;
//
//            // Perform or skip RLE
//            if (params->polishParams->useRunLengthEncoding) {
//                rleNucleotideString = rleString_construct(read->nucleotides);
//            } else {
//                rleNucleotideString = rleString_constructNoRLE(read->nucleotides);
//            }
//            totalNucleotides += rleNucleotideString->length;
//
//            // Do RLE follow up regardless of whether RLE is applied
//            stList_append(rleNucleotides, rleNucleotideString);
//            stList_append(rleReads, bamChunkRead_constructRLECopy(read, rleNucleotideString));
//            stList_append(rleAlignments, runLengthEncodeAlignment(alignment, rleReference, rleNucleotideString));
//        }


        // Run the polishing method
        int64_t totalNucleotides = 0;
        if (st_getLogLevel() >= info) {
            for (int64_t u = 0 ; u < stList_length(reads); u++) {
                totalNucleotides += strlen(((BamChunkRead*)stList_get(reads, u))->rleRead->rleString);
            }
            st_logInfo(">%s Running polishing algorithm with %"PRId64" reads and %"PRIu64"K nucleotides\n",
                       logIdentifier, stList_length(reads), totalNucleotides >> 10);
        }

        // Generate partial order alignment (POA) (destroys rleAlignments in the process)
        poa = poa_realignAll(reads, alignments, rleReference, params->polishParams);


        // get polished reference string and expand RLE (regardless of whether RLE was applied)
        poa_estimateRepeatCountsUsingBayesianModel(poa, reads, params->polishParams->repeatSubMatrix);
        RleString *polishedRleConsensus = poa->refString;
        polishedConsensusString = rleString_expand(polishedRleConsensus);


        // Log info about the POA
        if (st_getLogLevel() >= info) {
            st_logInfo(">%s Summary stats for POA:\t", logIdentifier);
            poa_printSummaryStats(poa, stderr);
        }
        if (st_getLogLevel() >= debug) {
            poa_print(poa, stderr, reads, 5, 5);
        }

        // Write any optional outputs about repeat count and POA, etc.
        if(outputPoaDotBase != NULL) {
            char *outputPoaDotFilename = stString_print("%s.poa.C%05"PRId64".%s-%"PRId64"-%"PRId64".dot",
                                                        outputPoaDotBase, chunkIdx, bamChunk->refSeqName,
                                                        bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
            FILE *outputPoaTsvFileHandle = fopen(outputPoaDotFilename, "w");
            poa_printDOT(poa, outputPoaTsvFileHandle, reads);
            fclose(outputPoaTsvFileHandle);
            free(outputPoaDotFilename);
        }
        if(outputPoaTsvBase != NULL) {
            char *outputPoaTsvFilename = stString_print("%s.poa.C%05"PRId64".%s-%"PRId64"-%"PRId64".tsv",
                                                        outputPoaTsvBase, chunkIdx, bamChunk->refSeqName,
                                                        bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
            FILE *outputPoaTsvFileHandle = fopen(outputPoaTsvFilename, "w");
            poa_printTSV(poa, outputPoaTsvFileHandle, reads, 5, 0);
            fclose(outputPoaTsvFileHandle);
            free(outputPoaTsvFilename);
        }
        if(outputRepeatCountBase != NULL) {
            char *outputRepeatCountFilename = stString_print("%s.repeatCount.C%05"PRId64".%s-%"PRId64"-%"PRId64".tsv",
                                                             outputRepeatCountBase, chunkIdx, bamChunk->refSeqName,
                                                             bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
            FILE *outputRepeatCountFileHandle = fopen(outputRepeatCountFilename, "w");
            poa_printRepeatCounts(poa, outputRepeatCountFileHandle, reads);
            fclose(outputRepeatCountFileHandle);
            free(outputRepeatCountFilename);
        }


        // save polished reference string to chunk output array
        chunkResults[chunkIdx] = polishedConsensusString;

        // HELEN feature outputs

        #ifdef _HDF5
        if (helenFeatureType != HFEAT_NONE) {
            handleHelenFeatures(helenFeatureType, trueReferenceBamChunker, splitWeightMaxRunLength,
                    helenHDF5Files, fullFeatureOutput, trueReferenceBam, params, logIdentifier, chunkIdx,
                    bamChunk, poa, reads, polishedConsensusString, polishedRleConsensus);

        }
        #endif

        // report timing
        if (st_getLogLevel() >= info) {
            st_logInfo(">%s Chunk with %"PRId64" reads and %"PRIu64"K nucleotides processed in %d sec\n",
                       logIdentifier, stList_length(reads), totalNucleotides >> 10,
                       (int) (time(NULL) - chunkStartTime));
        }

        // Cleanup
        rleString_destruct(rleReference);
        poa_destruct(poa);
        stList_destruct(reads);
        stList_destruct(alignments);
        free(logIdentifier);
    }

    // prep for merge
    assert(bamChunker->chunkCount > 0);
    int64_t contigStartIdx = 0;
    char *referenceSequenceName = stString_copy(bamChunker_getChunk(bamChunker, 0)->refSeqName);
    lastReportedPercentage = 0;
    time_t mergeStartTime = time(NULL);

    // for filling missing chunks with N's
    int64_t spacerSize = (bamChunker->chunkBoundary == 0 ? 50 : bamChunker->chunkBoundary * 3);
    char *missingChunkSpacer = st_calloc(spacerSize + 1, sizeof(char));
    for (int64_t i = 0; i < spacerSize; i++) {
        missingChunkSpacer[i] = 'N';
    }
    missingChunkSpacer[spacerSize] = '\0';

    // merge chunks
    st_logCritical("> Merging polished reference strings from %"PRIu64" chunks.\n", bamChunker->chunkCount);

    // find which chunks belong to each contig, merge each contig threaded, write out
    for (int64_t chunkIdx = 1; chunkIdx <= bamChunker->chunkCount; chunkIdx++) {
        
        // we encountered the last chunk in the contig (end of list or new refSeqName)
        if (chunkIdx == bamChunker->chunkCount || !stString_eq(referenceSequenceName,
                bamChunker_getChunk(bamChunker, chunkIdx)->refSeqName)) {

            // generate and save sequence
            char *contigSequence = mergeContigChunksThreaded(chunkResults, contigStartIdx, chunkIdx, numThreads, 
                    bamChunker->chunkBoundary * 2, params, missingChunkSpacer, referenceSequenceName);
            fastaWrite(contigSequence, referenceSequenceName, polishedReferenceOutFh);

            // log progress
            int64_t currentPercentage = (int64_t) (100 * chunkIdx / bamChunker->chunkCount);
            if (currentPercentage != lastReportedPercentage) {
                lastReportedPercentage = currentPercentage;
                int64_t timeTaken = (int64_t) (time(NULL) - mergeStartTime);
                int64_t secondsRemaining = (int64_t) floor(1.0 * timeTaken / currentPercentage * (100 - currentPercentage));
                char *timeDescriptor = (secondsRemaining == 0 && currentPercentage <= 50 ?
                                        stString_print("unknown") : getTimeDescriptorFromSeconds(secondsRemaining));
                st_logCritical("> Merging %2"PRId64"%% complete (%"PRId64"/%"PRId64").  Estimated time remaining: %s\n",
                        currentPercentage, chunkIdx, bamChunker->chunkCount, timeDescriptor);
                free(timeDescriptor);
            }

            // Clean up
            free(contigSequence);
            free(referenceSequenceName);

            // Reset for next reference sequence
            if (chunkIdx != bamChunker->chunkCount) {
                contigStartIdx = chunkIdx;
                referenceSequenceName = stString_copy(bamChunker_getChunk(bamChunker, chunkIdx)->refSeqName);
            }
        }
        // nothing to do otherwise, just wait until end or new contig
    }

    // everything has been written, cleanup merging infrastructure
    fclose(polishedReferenceOutFh);
    free(missingChunkSpacer);
    for (int64_t chunkIdx = 0; chunkIdx < bamChunker->chunkCount; chunkIdx++) {
        free(chunkResults[chunkIdx]);
    }

    // Cleanup
    bamChunker_destruct(bamChunker);
    stHash_destruct(referenceSequences);
    params_destruct(params);
    if (trueReferenceBam != NULL) free(trueReferenceBam);
    if (trueReferenceBamChunker != NULL) bamChunker_destruct(trueReferenceBamChunker);
    if (regionStr != NULL) free(regionStr);
    #ifdef _HDF5
    if (helenHDF5Files != NULL) {
        for (int64_t i = 0; i < numThreads; i++) {
            HelenFeatureHDF5FileInfo_destruct((HelenFeatureHDF5FileInfo *) helenHDF5Files[i]);
        }
        free(helenHDF5Files);
    }
    #endif
    free(chunkResults);
    free(outputBase);
    free(bamInFile);
    free(referenceFastaFile);
    free(paramsFile);

    // log completion
    char *timeDescriptor = getTimeDescriptorFromSeconds(time(NULL) - startTime);
    st_logCritical("> Finished polishing in %s.\n", timeDescriptor);
    free(timeDescriptor);

//    while(1); // Use this for testing for memory leaks

    return 0;
}

