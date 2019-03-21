/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com) & Arthur Rand (arand@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <getopt.h>
#include <stdio.h>
#include <ctype.h>
#include <memory.h>
#include <hashTableC.h>
#include <unistd.h>
#include "vcf.h"

#include "margin.h"
#include "marginVersion.h"
#include "vcfComparison.h"
#include "htsIntegration.h"
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

void getExpectedMatchesBetweenProfileSeqs(stProfileSeq *pSeq1, stProfileSeq *pSeq2,
                                          int64_t *totalAlignedPositions, double *totalExpectedMatches) {
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

stList *prefilterReads(stList *profileSequences, int64_t *misses,
                       stHash *referenceNamesToReferencePriors, stRPHmmParameters *params) {
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
    st_logInfo("\tGot %" PRIi64 " hmms before splitting\n", stList_length(hmms));

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
    st_logInfo("\tCreated %d hmms after splitting at uncertain regions of phasing (previously %d)\n",
                stList_length(hmms), initialHmmListSize);

    int64_t idx = 0;
    if(st_getLogLevel() == debug && stList_length(hmms) != initialHmmListSize) {
        stListIterator *itor = stList_getIterator(hmms);
        stRPHmm *hmm = NULL;
        while ((hmm = stList_getNext(itor)) != NULL) {
            st_logDebug("\thmm %3d: \tstart pos: %8d \tend pos: %8d\n",
                        idx, hmm->refStart, (hmm->refStart + hmm->refLength));
            idx++;
        }
    }
    return hmms;
}

/*
 * Main functions
 */

void usage() {
    fprintf(stderr, "usage: marginPhase <BAM> <REFERENCE_FASTA> <PARAMS> [options]\n");
    fprintf(stderr, "Version: %s \n\n", MARGINPHASE_MARGIN_PHASE_VERSION_H);
    fprintf(stderr, "Phases the reads in a BAM file and produces:\n");
    fprintf(stderr, "    1) a VCF file giving genotypes and haplotypes.\n");
    fprintf(stderr, "    2) a SAM file where each read is annotated with haplotype information\n");

    fprintf(stderr, "\nRequired arguments:\n");
    fprintf(stderr, "    BAM is the alignment of reads.  All reads must be aligned to the same contig \n");
    fprintf(stderr, "        and be in bam format.\n");
    fprintf(stderr, "    REFERENCE_FASTA is the reference sequence for the BAM's contig in fasta format.\n");
    fprintf(stderr, "    PARAMS is the file with marginPhase parameters.\n");

    fprintf(stderr, "\nDefault options:\n");
    fprintf(stderr, "    -h --help              : Print this help screen\n");
    fprintf(stderr, "    -o --outputBase        : Base output identifier [default = \"output\"]\n");
    fprintf(stderr, "                               \"example\" -> \"example.sam\", \"example.vcf\"\n");
    fprintf(stderr, "    -a --logLevel          : Set the log level [default = info]\n");
    fprintf(stderr, "    -t --tag               : Annotate all output reads with this value for the \n");
    fprintf(stderr, "                               '"MARGIN_PHASE_TAG"' tag\n");

    fprintf(stderr, "\nNucleotide probabilities options:\n");
    fprintf(stderr, "    -s --singleNuclProbDir : Directory of single nucleotide probabilities files\n");
    fprintf(stderr, "    -S --onlySNP           : Use only single nucleotide probabilities information,\n");
    fprintf(stderr, "                               so reads that aren't in SNP dir are discarded\n");

    fprintf(stderr, "\nVCF Comparison options:\n");
    fprintf(stderr, "    -r --referenceVCF      : Reference vcf file for output comparison\n");
    fprintf(stderr, "    -v --verbose           : Bitmask controlling outputs \n");
    fprintf(stderr, "                               %3d - LOG_TRUE_POSITIVES\n", LOG_TRUE_POSITIVES);
    fprintf(stderr, "                               %3d - LOG_FALSE_POSITIVES\n", LOG_FALSE_POSITIVES);
    fprintf(stderr, "                               %3d - LOG_FALSE_NEGATIVES\n", LOG_FALSE_NEGATIVES);
    fprintf(stderr, "                               example: 0 -> N/A; 2 -> LFP; 7 -> LTP,LFP,LFN)\n");

}

int main(int argc, char *argv[]) {

    // Parameters / arguments
    char *logLevelString = stString_copy("info");
    char *bamInFile = NULL;
    char *paramsFile = NULL;
    char *referenceFastaFile = NULL;
    char *marginPhaseTag = NULL;
    char *referenceVCF = NULL;
    char *singleNucleotideProbabilityDirectory = NULL;
    char *outputBase = "output";
    int64_t verboseBitstring = -1;
    bool onlySNP = false;

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
                { "referenceVcf", required_argument, 0, 'r'},
                { "tag", required_argument, 0, 't'},
                { "singleNuclProbDir", required_argument, 0, 's'},
                { "onlySNP", no_argument, 0, 'S'},
                { "verbose", required_argument, 0, 'v'},
                { 0, 0, 0, 0 } };

        int option_index = 0;
        int key = getopt_long(argc-2, &argv[2], "a:o:v:r:s:hS", long_options, &option_index);

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
        case 't':
            marginPhaseTag = stString_copy(optarg);
            break;
        case 'r':
            referenceVCF = stString_copy(optarg);
            break;
        case 's':
            singleNucleotideProbabilityDirectory = stString_copy(optarg);
            if (singleNucleotideProbabilityDirectory[strlen(optarg) - 1] == '/') {
                singleNucleotideProbabilityDirectory[strlen(optarg) - 1] = '\0';
            }
            break;
        case 'S':
            onlySNP = true;
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

    // Output file names
    char *vcfOutFile = stString_print("%s.vcf", outputBase);
    char *paramsOutFile = stString_print("%s.params.json", outputBase);
    char *vcfOutFile_all = stString_print("%s.gvcf", outputBase);

    // Parse any model parameters
    st_logInfo("> Parsing model parameters from file: %s\n", paramsFile);
    Params *fullParams = params_readParams2(paramsFile, FALSE, TRUE);
    stBaseMapper *baseMapper = fullParams->baseMapper;
    stRPHmmParameters *params = fullParams->phaseParams;
    assert(baseMapper != NULL);
    assert(params != NULL);
    if (verboseBitstring >= 0) setVerbosity(params, verboseBitstring); //run this AFTER parameters, so CL args overwrite

    // Print a report of the parsed parameters
    if(st_getLogLevel() == debug) {
        stRPHmmParameters_printParameters(params, stderr);
    }

    // Parse reads for interval
    st_logInfo("> Parsing input reads from file: %s\n", bamInFile);
    stList *profileSequences = stList_construct3(0, (void (*)(void *))stProfileSeq_destruct);
    int64_t readCount = 0;
    if (singleNucleotideProbabilityDirectory == NULL) {
        readCount = parseReads(profileSequences, bamInFile, baseMapper, params);
    } else {
        readCount = parseReadsWithSingleNucleotideProbs(profileSequences, bamInFile, baseMapper, params,
                                                        singleNucleotideProbabilityDirectory, onlySNP);
    }
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
        st_logInfo("> Pre-filtering reads to remove reads with less than %f identity to the consensus sequence\n",
                   params->filterMatchThreshold);
        profileSequences = prefilterReads(profileSequences, &misses, referenceNamesToReferencePriors, params);
        st_logInfo("\tFiltered %d profile sequences (%f percent)\n", misses, (float)misses*100/initialSize);
    }

    // Setup a filter to ignore likely homozygous reference positions
    if(params->filterLikelyHomozygousSites) {
        int64_t totalPositions;
        st_logInfo("> Filtering likely homozygous positions\n");
        int64_t filteredPositions =
                filterHomozygousReferencePositions(referenceNamesToReferencePriors, params, &totalPositions);
        st_logInfo("\tFiltered %" PRIi64 " (%f) likely homozygous positions, \n\teach with fewer than %" PRIi64
                " aligned occurrences of any second most frequent base, \n\tleaving only %" PRIi64
                           " (%f) positions of %" PRIi64
                " total positions\n", filteredPositions, (double)filteredPositions/totalPositions,
                (int64_t)params->minSecondMostFrequentBaseFilter, totalPositions - filteredPositions,
                (double)(totalPositions - filteredPositions)/totalPositions, totalPositions);
    }

    // Filter reads so that the maximum coverage depth does not exceed params->maxCoverageDepth
    st_logInfo("> Filtering reads by coverage depth\n");
    stList *filteredProfileSeqs = stList_construct3(0, (void (*)(void *))stProfileSeq_destruct);
    stList *discardedProfileSeqs = stList_construct3(0, (void (*)(void *))stProfileSeq_destruct);
    filterReadsByCoverageDepth(profileSequences, params, filteredProfileSeqs, discardedProfileSeqs,
                               referenceNamesToReferencePriors);
    st_logInfo("\tFiltered %" PRIi64 " reads of %" PRIi64
                       " to achieve maximum coverage depth of %" PRIi64 "\n",
               stList_length(discardedProfileSeqs), stList_length(profileSequences),
               params->maxCoverageDepth);
    stList_destruct(discardedProfileSeqs);
    stList_setDestructor(profileSequences, NULL);
    stList_destruct(profileSequences);
    profileSequences = filteredProfileSeqs;

    // Estimate the read error substitution matrix from the alignment of the reads
    // to the reference and set the read error substitution probs
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
    stReadHaplotypePartitionTable *readHaplotypePartitions = stReadHaplotypePartitionTable_construct(
            stList_length(profileSequences));
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

        // Only one haplotype found (likely a small set of reads)
        if (stSet_size(reads1) < 1 || stSet_size(reads2) < 1) {
            populateReadHaplotypePartitionTable(readHaplotypePartitions, gF, hmm, path);
            continue;
        }

        // Refine the genome fragment by repartitoning the reads iteratively
        if(params->roundsOfIterativeRefinement > 0) {
            stGenomeFragment_refineGenomeFragment(gF, reads1, reads2, hmm, path, params->roundsOfIterativeRefinement);
        }

        // save bipartition
        populateReadHaplotypePartitionTable(readHaplotypePartitions, gF, hmm, path);

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
    st_logInfo("\n\tFinished writing out VCF into file: %s\n", vcfOutFile);

    st_logInfo("\n> There were a total of %d genome fragments. Average length = %f\n", stList_length(hmms),
               (float) totalGFlength / stList_length(hmms));

    // do comparison if referenceVCF is specified
    if (referenceVCF != NULL) {
        st_logInfo("\n> Comparing VCFs\n");
        if (access(referenceVCF, F_OK) == -1) {
            st_logInfo("\tReference VCF file does not exist: %s\n", referenceVCF);
        } else {
            // Compare the output vcf with the reference vcf
            stGenotypeResults *results = st_calloc(1, sizeof(stGenotypeResults));
            compareVCFs(stderr, hmms, vcfOutFile, referenceVCF, baseMapper, results, params);
            printGenotypeResults(results);
            free(results);
        }
    }

    // Write out two BAMs, one for each read partition
    if (params->writeSplitSams) {
        st_logInfo("\tWriting out SAM files for each partition\n", outputBase,
                   outputBase);
        writeSplitSams(bamInFile, outputBase, readHaplotypePartitions, marginPhaseTag);
    }
    if (params->writeUnifiedSam) {
        st_logInfo("\tWriting out single SAM file with read partitioning\n", outputBase,
                   outputBase);
        writeHaplotypedSam(bamInFile, outputBase, readHaplotypePartitions, marginPhaseTag);
    }

    stList_destruct(profileSequences);
    stReadHaplotypePartitionTable_destruct(readHaplotypePartitions);
    stList_destruct(hmms);

    params_destruct(fullParams);
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

