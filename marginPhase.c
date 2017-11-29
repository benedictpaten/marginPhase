/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com) & Arthur Rand (arand@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <getopt.h>
#include <stdio.h>
#include <ctype.h>
#include <memory.h>

#include "stRPHmm.h"
#include "sonLib.h"
#include "externalTools/sonLib/C/impl/sonLibListPrivate.h"

/*
 * Functions used to print debug output.
 */

void printSequenceStats(FILE *fH, stList *profileSequences) {
    /*
     * Print stats about the set of profile sequences.
     */
    int64_t totalLength=0;
    for(int64_t i=0; i<stList_length(profileSequences); i++) {
        stProfileSeq *profileSeq = stList_get(profileSequences, i);
        totalLength += profileSeq->length;
    }
    fprintf(fH, "\tGot %" PRIi64 " profile sequences, with total length: %" PRIi64 ", average length: %f\n",
            stList_length(profileSequences), totalLength, ((float)totalLength)/stList_length(profileSequences));
}

double getExpectedNumberOfMatches(uint64_t *haplotypeString, int64_t start, int64_t length, stProfileSeq *profileSeq) {
    /*
     * Returns the expected number of positions in the profile sequence that are identical to the given haplotype string.
     */
    double totalExpectedMatches = 0.0;

    for(int64_t i=0; i<profileSeq->length; i++) {
        // Get base in the haplotype sequence
        int64_t j = i + profileSeq->refStart - start;
        if(j >= 0 && j < length) {
            uint64_t hapBase = haplotypeString[j];
            assert(hapBase < ALPHABET_SIZE);

            // Expectation of a match
            totalExpectedMatches += getProb(&(profileSeq->profileProbs[i * ALPHABET_SIZE]), hapBase);
        }
    }
    return totalExpectedMatches;
}

double getExpectedIdentity(uint64_t *haplotypeString, int64_t start, int64_t length, stSet *profileSeqs) {
    /*
     * Returns the expected fraction of positions in the profile sequences that match their corresponding position in the
     * given haplotype string.
     */
    double totalExpectedNumberOfMatches = 0.0;
    int64_t totalLength = 0;
    stSetIterator *it = stSet_getIterator(profileSeqs);
    stProfileSeq *pSeq;
    while((pSeq = stSet_getNext(it)) != NULL) {
        totalExpectedNumberOfMatches += getExpectedNumberOfMatches(haplotypeString, start, length, pSeq);
        totalLength += pSeq->length;
    }
    stSet_destructIterator(it);
    return totalExpectedNumberOfMatches/totalLength;
}

double getIdentityBetweenHaplotypes(uint64_t *hap1String, uint64_t *hap2String, int64_t length) {
    /*
     * Returns the fraction of positions in two haplotypes that are identical.
     */
    int64_t totalMatches = 0;
    for(int64_t i=0; i<length; i++) {
        if(hap1String[i] == hap2String[i]) {
            totalMatches++;
        }
    }
    return ((double)totalMatches)/length;
}

double getIdentityBetweenHaplotypesExcludingIndels(uint64_t *hap1String, uint64_t *hap2String, int64_t length) {
    /*
     * Returns the fraction of positions in two haplotypes that are identical.
     */
    int64_t totalMatches = 0;
    int64_t numGaps = 0;
    for(int64_t i=0; i<length; i++) {
        if(hap1String[i] == hap2String[i]) {
            totalMatches++;
        } else if (hap1String[i] == ALPHABET_SIZE-1 || hap2String[i] == ALPHABET_SIZE-1) {
            numGaps++;
        }
    }
    return ((double)totalMatches)/(length - numGaps);
}

void getExpectedMatchesBetweenProfileSeqs(stProfileSeq *pSeq1, stProfileSeq *pSeq2, int64_t *totalAlignedPositions, double *totalExpectedMatches) {
    /*
     * Calculates the number of base overlaps and expected base matches between two profile sequences.
     */
    for(int64_t i=0; i<pSeq1->length; i++) {
        // Establish if the coordinate is in both sequences
        int64_t j = i + pSeq1->refStart - pSeq2->refStart;
        if(j >= 0 && j < pSeq2->length) {
            (*totalAlignedPositions)++;

            // Calculate expectation of match
            for(int64_t k=0; k<ALPHABET_SIZE; k++) {
                double e1 = getProb(&(pSeq1->profileProbs[i * ALPHABET_SIZE]), k);
                double e2 = getProb(&(pSeq2->profileProbs[j * ALPHABET_SIZE]), k);
                assert(e1 * e2 <= 1.0);
                *totalExpectedMatches += e1 * e2;
            }
        }
    }
}

void printAvgIdentityBetweenProfileSequences(FILE *fH, stList *profileSequences, int64_t maxSeqs) {
    /*
     * Prints the average base identity between pairwise base overlaps between the given profile sequences
     */
    double totalExpectedMatches = 0.0;
    int64_t totalAlignedPositions = 0;

    for(int64_t i=0; i<stList_length(profileSequences) && i<maxSeqs; i++) {
        for(int64_t j=i+1; j<stList_length(profileSequences) && j<maxSeqs; j++) {
            getExpectedMatchesBetweenProfileSeqs(stList_get(profileSequences, i), stList_get(profileSequences, j),
                                                &totalAlignedPositions, &totalExpectedMatches);

        }
    }

    fprintf(fH, "\tAvg. pairwise identity between profile sequences: %f measured at %" PRIi64 " overlapping sites\n",
            totalExpectedMatches/totalAlignedPositions, totalAlignedPositions);
}
double *getHaplotypeBaseComposition(uint64_t *hapString, int64_t length) {
    /*
     * Get the count of each alphabet character in the haplotype sequence, returned
     * as an array.
     */
    double *baseCounts = st_calloc(ALPHABET_SIZE, sizeof(double));
    for(int64_t i=0; i<length; i++) {
        baseCounts[hapString[i]] += 1;
    }
    return baseCounts;
}

double *getExpectedProfileSequenceBaseComposition(stSet *profileSeqs) {
    /*
     * Get the expected count of each alphabet character in the profile sequences, returned
     * as an array.
     */
    double *baseCounts = st_calloc(ALPHABET_SIZE, sizeof(double));
    stSetIterator *it = stSet_getIterator(profileSeqs);
    stProfileSeq *pSeq;
    while((pSeq = stSet_getNext(it)) != NULL) {
        for(int64_t i=0; i<pSeq->length; i++) {
            for(int64_t j=0; j<ALPHABET_SIZE; j++) {
                baseCounts[j] += getProb(&(pSeq->profileProbs[i*ALPHABET_SIZE]), j);
            }
        }
    }
    return baseCounts;
}

void printBaseComposition(FILE *fH, double *baseCounts) {
    /*
     * Print the counts/fraction of each alphabet character.
     */
    double totalCount = 0;
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        totalCount += baseCounts[i];
    }
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        fprintf(fH, "\tBase %" PRIi64 " count: %f fraction: %f\n", i, baseCounts[i], baseCounts[i]/totalCount);
    }
}

double matchesTopTwoBases(int64_t pos, stProfileSeq *profileSeq, stReferencePriorProbs *rProbs) {
    /*
     * Returns the probability that the profile sequence matched the most common base at
     * position pos, or that it matched either of the top two if the locus looks possibly heterozygous
     */
    uint8_t base1, base2;
    double max1 = 0.0;
    double max2 = 0.0;
    for (uint8_t i = 0; i < ALPHABET_SIZE; i++) {
        double p = rProbs->baseCounts[pos * ALPHABET_SIZE + i];
        if (p > max2) {
            if (p > max1) {
                max2 = max1;
                max1 = p;
                base2 = base1;
                base1 = i;
            } else {
                max2 = p;
                base2 = i;
            }
        }
    }
    int64_t pPos = pos - profileSeq->refStart + rProbs->refStart;
    double p = getProb(&(profileSeq->profileProbs[pPos * ALPHABET_SIZE]), base1);
    // TODO check these parameters
    if (max2 > (max1 / ALPHABET_SIZE) && max2 > 3) {
        p += getProb(&(profileSeq->profileProbs[pPos * ALPHABET_SIZE]), base2);
    }
    return p;
}

double getExpectedNumberOfConsensusMatches(stProfileSeq *profileSeq, stReferencePriorProbs *rProbs) {
    /*
     * Returns the expected number of positions in the profile sequence that are identical
     * to the given haplotype string.
     */
    double totalExpectedMatches = 0.0;

    for(int64_t i=0; i<profileSeq->length; i++) {
        // Get base in the haplotype sequence
        int64_t j = i + profileSeq->refStart - rProbs->refStart;
        if(j >= 0 && j < rProbs->length) {
            // Add a match if either of the top two bases are matched
            totalExpectedMatches +=matchesTopTwoBases(j, profileSeq, rProbs);
        }
    }
    return totalExpectedMatches;
}

stList *prefilterReads(stList *profileSequences, int64_t *misses, stHash *referenceNamesToReferencePriors, stRPHmmParameters *params) {
    /*
     * Filter out profile sequences that don't have identity with the consensus reference sequence
     * above a specified threshold.
     */
    stList *filteredProfileSequences = stList_construct3(0, (void (*)(void *))stProfileSeq_destruct);

    for(int64_t i=0; i<stList_length(profileSequences); i++) {
        stProfileSeq *pSeq = stList_get(profileSequences, i);
        stReferencePriorProbs *rProbs = stHash_search(referenceNamesToReferencePriors, pSeq->referenceName);
        double identity = getExpectedNumberOfConsensusMatches(pSeq, rProbs) / pSeq->length;
        if (identity < params->filterMatchThreshold) {
            (*misses)++;
            stProfileSeq_destruct(pSeq);
        } else {
            stList_append(filteredProfileSequences, pSeq);
        }
    }
    stList_setDestructor(profileSequences, NULL);
    stList_destruct(profileSequences);
    return filteredProfileSequences;
}

stList *createHMMs(stList *profileSequences, stHash *referenceNamesToReferencePriors, stRPHmmParameters *params) {
    /*
     * Create the set of hmms that the forward-backward algorithm will eventually be run on.
     */
    // Create the initial list of HMMs
    st_logInfo("> Creating read partitioning HMMs\n");
    stList *hmms = getRPHmms(profileSequences, referenceNamesToReferencePriors, params);
    st_logInfo("Got %" PRIi64 " hmms before splitting\n", stList_length(hmms));

    // Break up the hmms where the phasing is uncertain
    st_logInfo("> Breaking apart HMMs where the phasing is uncertain\n");

    stList *l = stList_construct3(0, (void (*)(void *))stRPHmm_destruct2);
    int64_t initialHmmListSize = stList_length(hmms);

    // Reverse to make sure hmm list stays in correct order
    stList_reverse(hmms);
    while(stList_length(hmms) > 0) {
        stList_appendAll(l, stRPHMM_splitWherePhasingIsUncertain(stList_pop(hmms)));
    }
    hmms = l;
    st_logInfo("Created %d hmms after splitting at uncertain regions of phasing (previously %d)\n",
                stList_length(hmms), initialHmmListSize);

    int64_t idx = 0;
    if(st_getLogLevel() == debug && stList_length(hmms) != initialHmmListSize) {
        stListIterator *itor = stList_getIterator(hmms);
        stRPHmm *hmm = NULL;
        while ((hmm = stList_getNext(itor)) != NULL) {
            st_logDebug("\thmm %3d: \tstart pos: %8d \tend pos: %8d\n", idx, hmm->refStart, (hmm->refStart + hmm->refLength));
            idx++;
        }
    }
    return hmms;
}

void logHmm(stRPHmm *hmm, stSet *reads1, stSet *reads2, stGenomeFragment *gF) {
    /*
     * Print debug-level logging information about an HMM and associated genome fragment.
     */
    if(st_getLogLevel() == debug) {
        st_logDebug("> Creating genome fragment for reference sequence: %s, start: %" PRIi64 ", length: %" PRIi64 "\n",
                    hmm->referenceName, hmm->refStart, hmm->refLength);
        st_logDebug("\nThere are %" PRIi64 " reads covered by the hmm, bipartitioned into sets of %" PRIi64 " and %" PRIi64 " reads\n",
                    stList_length(hmm->profileSeqs), stSet_size(reads1), stSet_size(reads2));

        // Print the similarity between the two imputed haplotypes sequences
        st_logDebug("The haplotypes have an %f identity\n", getIdentityBetweenHaplotypes(gF->haplotypeString1, gF->haplotypeString2, gF->length));
        st_logDebug("\tIdentity excluding indels: %f \n\n", getIdentityBetweenHaplotypesExcludingIndels(gF->haplotypeString1, gF->haplotypeString2, gF->length));

        // Print the base composition of the haplotype sequences
        double *hap1BaseCounts = getHaplotypeBaseComposition(gF->haplotypeString1, gF->length);
        st_logDebug("The base composition of haplotype 1:\n");
        printBaseComposition(stderr, hap1BaseCounts);
        free(hap1BaseCounts);

        double *hap2BaseCounts = getHaplotypeBaseComposition(gF->haplotypeString2, gF->length);
        st_logDebug("The base composition of haplotype 2:\n");
        printBaseComposition(stderr, hap2BaseCounts);
        free(hap2BaseCounts);

        // Print the base composition of the reads
        double *reads1BaseCounts =getExpectedProfileSequenceBaseComposition(reads1);
        st_logDebug("The base composition of reads1 set:\n");
        printBaseComposition(stderr, reads1BaseCounts);
        free(reads1BaseCounts);

        double *reads2BaseCounts =getExpectedProfileSequenceBaseComposition(reads2);
        st_logDebug("The base composition of reads2 set:\n");
        printBaseComposition(stderr, reads2BaseCounts);
        free(reads2BaseCounts);

        // Print some summary stats about the differences between haplotype sequences and the bipartitioned reads
        st_logDebug("hap1 vs. reads1 identity: %f\n",
                    getExpectedIdentity(gF->haplotypeString1, gF->refStart, gF->length, reads1));
        st_logDebug("hap1 vs. reads2 identity: %f\n",
                    getExpectedIdentity(gF->haplotypeString1, gF->refStart, gF->length, reads2));
        st_logDebug("hap2 vs. reads2 identity: %f\n",
                    getExpectedIdentity(gF->haplotypeString2, gF->refStart, gF->length, reads2));
        st_logDebug("hap2 vs. reads1 identity: %f\n",
                    getExpectedIdentity(gF->haplotypeString2, gF->refStart, gF->length, reads1));
    }
}


/*
 * Main functions
 */

void usage() {
fprintf(stderr, "marginPhase <BAM_FILE> <REFERENCE_FASTA> [options]\n");
    fprintf(stderr,
            "Phases the reads in an interval of a BAM file (BAM_FILE) reporting a gVCF file "
            "giving genotypes and haplotypes for region.\n"
            "REFERENCE_FASTA is the reference sequence for the region in fasta format.\n");
    fprintf(stderr, "-h --help            : Print this help screen\n");
    fprintf(stderr, "-a --logLevel        : Set the log level [default = info]\n");
    fprintf(stderr, "-o --outputBase      : Output Base (\"example\" -> \"example1.sam\", \"example2.sam\", \"example.vcf\")\n");
    fprintf(stderr, "-p --params          : Input params file\n");
    fprintf(stderr, "-r --referenceVCF    : Reference vcf file, to compare output to\n");
    fprintf(stderr, "-v --verbose         : Bitmask controlling outputs\n");
    fprintf(stderr, "                     \t%3d - LOG_TRUE_POSITIVES\n", LOG_TRUE_POSITIVES);
}

int main(int argc, char *argv[]) {

    // Parameters / arguments
    char *logLevelString = stString_copy("info");
    char *bamInFile = NULL;
    char *referenceFastaFile = NULL;
    char *referenceVCF = NULL;

    char *outputBase = "output";
    char *paramsFile = "params.json";
    int64_t verboseBitstring = -1;

    // TODO: When done testing, optionally set random seed using st_randomSeed();

    if(argc < 3) {
        usage();
        return 0;
    }

    bamInFile = stString_copy(argv[1]);
    referenceFastaFile = stString_copy(argv[2]);

    // Parse the options
    while (1) {
        static struct option long_options[] = {
                { "logLevel", required_argument, 0, 'a' },
                { "help", no_argument, 0, 'h' },
                { "outputBase", required_argument, 0, 'o'},
                { "params", required_argument, 0, 'p'},
                { "referenceVcf", required_argument, 0, 'r'},
                { "verbose", optional_argument, 0, 'v'},
                { 0, 0, 0, 0 } };

        int option_index = 0;
        int key = getopt_long(argc-2, &argv[2], "a:o:v::p:r:h", long_options, &option_index);

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
            outputBase = stString_copy(optarg);
            break;
        case 'p':
            paramsFile = stString_copy(optarg);
            break;
        case 'r':
            referenceVCF = stString_copy(optarg);
            break;
        case 'v':
            if (optarg == NULL) verboseBitstring = (1 << 16) - 1;
            else verboseBitstring = atoi(optarg);
            break;
        default:
            usage();
            return 0;
        }
    }

    // Initialization from arguments
    st_setLogLevelFromString(logLevelString);
    free(logLevelString);

    // Output file names
    char *vcfOutFile = stString_print("%s.vcf", outputBase);
    char *paramsOutFile = stString_print("%s.params.json", outputBase);
    char *vcfOutFile_all = stString_print("%s.gvcf", outputBase);

    // Parse any model parameters
    st_logInfo("> Parsing model parameters from file: %s\n", paramsFile);
    stBaseMapper *baseMapper = stBaseMapper_construct();
    stRPHmmParameters *params = parseParameters(paramsFile, baseMapper);
    if (verboseBitstring >= 0) setVerbosity(params, verboseBitstring); //run this AFTER parameters, so CL args overwrite

    // Print a report of the parsed parameters
    if(st_getLogLevel() == debug) {
        stRPHmmParameters_printParameters(params, stderr);
    }

    // Parse reads for interval
    st_logInfo("> Parsing input reads from file: %s\n", bamInFile);
    stList *profileSequences = stList_construct3(0, (void (*)(void *))stProfileSeq_destruct);
    int64_t readCount = parseReads(profileSequences, bamInFile, baseMapper, params);
    st_logInfo("\tCreated %d profile sequences\n", readCount);

    // Print some stats about the input sequences
    if(st_getLogLevel() == debug) {
        printSequenceStats(stderr, profileSequences);
        printAvgIdentityBetweenProfileSequences(stderr, profileSequences, 100);
    }

    // Getting reference sequence priors
    stHash *referenceNamesToReferencePriors;
    if(!params->useReferencePrior) {
        st_logInfo("> Using a flat prior over reference positions\n");
        referenceNamesToReferencePriors = createEmptyReferencePriorProbabilities(profileSequences);
    }
    else {
        st_logInfo("> Parsing prior probabilities on positions from reference sequences: %s\n", referenceFastaFile);
        referenceNamesToReferencePriors = createReferencePriorProbabilities(referenceFastaFile, profileSequences,
                baseMapper, params);
	}

    // Filter reads that are too divergent from the consensus sequence
    // (taking into account spots that are potentially heterozygous)
    if (params->filterBadReads) {
        int64_t initialSize = stList_length(profileSequences);
        int64_t misses = 0;
        st_logInfo("> Pre-filtering reads to remove reads with less than %f identity to the consensus sequence\n", params->filterMatchThreshold);
        profileSequences = prefilterReads(profileSequences, &misses, referenceNamesToReferencePriors, params);
        st_logInfo("\tFiltered %d profile sequences (%f percent)\n", misses, (float)misses*100/initialSize);
    }

    // Setup a filter to ignore likely homozygous reference positions
    if(params->filterLikelyHomozygousSites) {
        int64_t totalPositions;
        st_logInfo("> Filtering likely homozygous positions\n");
        int64_t filteredPositions = filterHomozygousReferencePositions(referenceNamesToReferencePriors, params, &totalPositions);
        st_logInfo("\tFiltered %" PRIi64 " (%f) likely homozygous positions, \n\teach with fewer than %" PRIi64
                " aligned occurrences of any second most frequent base, \n\tleaving only %" PRIi64 " (%f) positions of %" PRIi64
                " total positions\n", filteredPositions, (double)filteredPositions/totalPositions,
                (int64_t)params->minSecondMostFrequentBaseFilter, totalPositions - filteredPositions,
                (double)(totalPositions - filteredPositions)/totalPositions, totalPositions);
    }

    // Filter reads so that the maximum coverage depth does not exceed params->maxCoverageDepth
    st_logInfo("> Filtering reads by coverage depth\n");
    stList *filteredProfileSeqs = stList_construct3(0, (void (*)(void *))stProfileSeq_destruct);
    stList *discardedProfileSeqs = stList_construct3(0, (void (*)(void *))stProfileSeq_destruct);
    filterReadsByCoverageDepth(profileSequences, params, filteredProfileSeqs, discardedProfileSeqs, referenceNamesToReferencePriors);
    st_logInfo("\tFiltered %" PRIi64 " reads of %" PRIi64
                       " to achieve maximum coverage depth of %" PRIi64 "\n",
               stList_length(discardedProfileSeqs), stList_length(profileSequences),
               params->maxCoverageDepth);
    stList_destruct(discardedProfileSeqs);
    stList_setDestructor(profileSequences, NULL);
    stList_destruct(profileSequences);
    profileSequences = filteredProfileSeqs;

    // Estimate the read error substitution matrix from the alignment of the reads to the reference and set the read error
    // substitution probs
    if(params->estimateReadErrorProbsEmpirically) {
        st_logInfo("> Estimating read errors from alignment (quick and dirty)\n");

        // Make read error substitution matrix
        double *readErrorSubModel = stReferencePriorProbs_estimateReadErrorProbs(referenceNamesToReferencePriors, params);
        // Set substitution probabilities
        stRPHmmParameters_setReadErrorSubstitutionParameters(params, readErrorSubModel);

        // Cleanup
        free(readErrorSubModel);
        if(st_getLogLevel() == debug) {
            stRPHmmParameters_printParameters(params, stderr);
        }
    }

    // Learn the parameters for the input data
    if(params->trainingIterations > 0) {
        st_logInfo("> Learning parameters for HMM model (%" PRIi64 " iterations)\n", params->trainingIterations);
        stRPHmmParameters_learnParameters(params, profileSequences, referenceNamesToReferencePriors);

        // Print a report of the parsed parameters
        if(st_getLogLevel() == debug) {
            st_logDebug("> Learned parameters:\n");
            stRPHmmParameters_printParameters(params, stderr);
        }

        st_logInfo("\tWriting learned parameters to file: %s\n", paramsOutFile);
        writeParamFile(paramsOutFile, params);
    }

    // Get the final list of hmms
    stList *hmms = createHMMs(profileSequences, referenceNamesToReferencePriors, params);

    // Prep for BAM outputs
    stSet *read1Ids = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, NULL);
    stSet *read2Ids = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, NULL);
    int64_t totalGFlength = 0;

    // Start VCF generation
    vcfFile *vcfOutFP = vcf_open(vcfOutFile, "w");
    bcf_hdr_t *hdr = writeVcfHeader(vcfOutFP, hmms, referenceFastaFile);
    vcfFile *vcfOutFP_all;
    bcf_hdr_t *hdr2;
    if (params->writeGVCF) {
        // Write gVCF if specified to
        vcfOutFP_all = vcf_open(vcfOutFile_all, "w");
        hdr2 = writeVcfHeader(vcfOutFP_all, hmms, referenceFastaFile);
    }

    // For each read partitioning HMM
    for(int64_t i=0; i<stList_length(hmms); i++) {
        stRPHmm *hmm = stList_get(hmms, i);

        // Run the forward-backward algorithm
        stRPHmm_forwardBackward(hmm);

        // Now compute a high probability path through the hmm
        stList *path = stRPHmm_forwardTraceBack(hmm);

        // Compute the genome fragment
        stGenomeFragment *gF = stGenomeFragment_construct(hmm, path);
        totalGFlength += gF->length;

        // Get the reads which mapped to each path
        stSet *reads1 = stRPHmm_partitionSequencesByStatePath(hmm, path, true);
        stSet *reads2 = stRPHmm_partitionSequencesByStatePath(hmm, path, false);

        if (stSet_size(reads1) < 1 || stSet_size(reads2) < 1) {
            continue;
        }

        // Refine the genome fragment by repartitoning the reads iteratively
        if(params->roundsOfIterativeRefinement > 0) {
            stGenomeFragment_refineGenomeFragment(gF, reads1, reads2, hmm, path, params->roundsOfIterativeRefinement);
        }

        addProfileSeqIdsToSet(reads1, read1Ids);
        addProfileSeqIdsToSet(reads2, read2Ids);

        // Log information about the hmm
        logHmm(hmm, reads1, reads2, gF);

        // Write two vcfs, one using the reference fasta file and one not
        writeVcfFragment(vcfOutFP, hdr, gF, referenceFastaFile, baseMapper, false);
        if (params->writeGVCF) {
            writeVcfFragment(vcfOutFP_all, hdr2, gF, referenceFastaFile, baseMapper, true);
        }


        // Cleanup
        stGenomeFragment_destruct(gF);
        stSet_destruct(reads1);
        stSet_destruct(reads2);
        stList_destruct(path);
    }
    // Cleanup vcf
    vcf_close(vcfOutFP);
    bcf_hdr_destroy(hdr);
    if (params->writeGVCF) {
        vcf_close(vcfOutFP_all);
        bcf_hdr_destroy(hdr2);
    }

    // Write out VCF
    st_logInfo("Finished writing out VCF into file: %s\n", vcfOutFile);

    st_logInfo("\nThere were a total of %d genome fragments. Average length = %f\n", stList_length(hmms),
               (float) totalGFlength / stList_length(hmms));


    if (params->compareVCFs) {
        // Compare the output vcf with the reference vcf
        stGenotypeResults *results = st_calloc(1, sizeof(stGenotypeResults));
        compareVCFs(stderr, hmms, vcfOutFile, referenceVCF, baseMapper, results, params);
        printGenotypeResults(results);
        free(results);
    }

    // Write out two BAMs, one for each read partition
    if (params->writeSplitSams) {
        st_logInfo("\n> Writing out SAM files for each partition into files: %s.0.sam and %s.1.sam\n", outputBase,
                   outputBase);
        writeSplitSams(bamInFile, outputBase, read1Ids, read2Ids);
    }
    if (params->writeSplitBams) {
        st_logInfo("\n> Writing out BAM files for each partition into files: %s.0.bam and %s.1.bam\n", outputBase,
                   outputBase);
        writeSplitBams(bamInFile, outputBase, read1Ids, read2Ids);
    }

    stList_destruct(profileSequences);
    stSet_destruct(read1Ids);
    stSet_destruct(read2Ids);
    stList_destruct(hmms);

    stBaseMapper_destruct(baseMapper);
    stRPHmmParameters_destruct(params);
    stHash_destruct(referenceNamesToReferencePriors);

    // TODO: only free these if they need to be
//    free(paramsFile);
//    free(bamInFile);
//    free(samOutBase);
//    free(vcfOutFile);
//    free(referenceName);

    //while(1); // Use this for testing for memory leaks

    return 0;
}

