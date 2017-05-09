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

void usage() {
    fprintf(stderr, "marginPhase [options] BAM_FILE\n");
    fprintf(stderr,
            "Phases the reads in an interval of a BAM file reporting a gVCF file "
            "giving genotypes and haplotypes for region.\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
    fprintf(stderr, "-b --bamFile : Input BAM file\n");
    fprintf(stderr, "-v --vcfFile : Output VCF file\n");
    fprintf(stderr, "-p --params : Input params file\n");
    fprintf(stderr, "-r --reference : Reference fasta file\n");
}

void printGenotypeFragment(stGenomeFragment *gF, bool printProbs) {
    st_logDebug("*** Genome Fragment Information ***\n");
    st_logDebug("Reference name: %s\n", gF->referenceName);
    st_logDebug("Ref start: %d \n", gF->refStart);
    st_logDebug("Length: %d \n", gF->length);

    st_logDebug("Genotype string: %u\n", gF->genotypeString);
    for (int64_t j = 0; j < gF->length; j++) {
        st_logDebug("%u\t", gF->genotypeString[j]);
    }
    st_logDebug("\nHaplotype 1 string: %u\n", gF->haplotypeString1);
    for (int64_t j = 0; j < gF->length; j++) {
        st_logDebug("%u\t", gF->haplotypeString1[j]);
    }
    st_logDebug("\nHaplotype 2 string: %u\n", gF->haplotypeString2);
    for (int64_t j = 0; j < gF->length; j++) {
        st_logDebug("%u\t", gF->haplotypeString2[j]);
    }

    if (printProbs) {
        st_logDebug("Genotype Probabilities: \n");
        for (int64_t j = 0; j < gF->length; j++) {
            st_logDebug("%f \t", gF->genotypeProbs[j]);
        }
        st_logDebug("\nHaplotype 1 Probabilities: \n");
        for (int64_t j = 0; j < gF->length; j++) {
            st_logDebug("%f \t", gF->haplotypeProbs1[j]);
        }
        st_logDebug("\nHaplotype 2 Probabilities: \n");
        for (int64_t j = 0; j < gF->length; j++) {
            st_logDebug("%f \t", gF->haplotypeProbs2[j]);
        }
    }
}


int main(int argc, char *argv[]) {
    // Parameters / arguments
    char * logLevelString = NULL;
    char *bamFile = NULL;
    char *bamInFile = NULL;
    char *vcfOutFile = "output.vcf";
    char *vcfOutFile2 = "output2.vcf";
    char *paramsFile = "params.json";
    char *referenceName = "hg19.chr3.fa";

    // TODO: When done testing, set random seed
    // st_randomSeed();

    // Parse the options
    while (1) {
        static struct option long_options[] = {
                { "logLevel", required_argument, 0, 'a' },
                { "help", no_argument, 0, 'h' },
                { "bamFile", required_argument, 0, 'b'},
                { "vcfOutFile", required_argument, 0, 'v'},
                { "params", required_argument, 0, 'p'},
                { "reference", required_argument, 0, 'r'},
                { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:v:p:r:h", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
        case 'a':
            logLevelString = stString_copy(optarg);
            st_setLogLevelFromString(logLevelString);
            break;
        case 'h':
            usage();
            return 0;
        case 'b':
            bamInFile = stString_copy(optarg);
            break;
        case 'v':
            vcfOutFile = stString_copy(optarg);
            break;
        case 'p':
            paramsFile = stString_copy(optarg);
            break;
        case 'r':
            referenceName = stString_copy(optarg);
            break;
        default:
            usage();
            return 1;
        }
    }

    // Parse any model parameters
    st_logInfo("Parsing model parameters\n");
    stBaseMapper *baseMapper = stBaseMapper_construct();
    stRPHmmParameters *params = parseParameters(paramsFile, baseMapper);

    // Parse reads for interval
    st_logInfo("Parsing input reads\n");
    stList *profileSequences = stList_construct();
    parseReads(profileSequences, bamInFile, baseMapper);

    // Create HMMs
    st_logInfo("Creating read partitioning HMMs\n");
    stList *hmms = getRPHmms(profileSequences, params);

    // Break up the hmms where the phasing is uncertain
    st_logInfo("Breaking apart HMMs where the phasing is uncertain\n");

    stList *l = stList_construct3(0, (void (*)(void *))stRPHmm_destruct2);
    while(stList_length(hmms) > 0) {
        stList_appendAll(l, stRPHMM_splitWherePhasingIsUncertain(stList_pop(hmms)));
    }
    hmms = l;

    // Get reference (needed for VCF generation)

    faidx_t *fai = fai_load(referenceName);
    if ( !fai ) {
        st_logCritical("Could not load fai index of %s.  Maybe you should run 'samtools faidx %s'\n",
                       referenceName, referenceName);
        return EXIT_FAILURE; //todo close things?
    }

    // Start VCF generation
    st_logDebug("Writing out VCF header %s\n", vcfOutFile);
    vcfFile *vcfOutFP = vcf_open(vcfOutFile, "w");
    bcf_hdr_t *hdr = writeVcfHeader(vcfOutFP, l);

    vcfFile *vcfOutFP2 = vcf_open(vcfOutFile2, "w");
    bcf_hdr_t *hdr2 = writeVcfHeader(vcfOutFP2, l);
    kstring_t str = {0,0,NULL};

    st_logDebug("Created %d hmms \n", stList_length(hmms));
    // For each read partitioning HMM
    for(int64_t i=0; i<stList_length(hmms); i++) {
        stRPHmm *hmm = stList_get(hmms, i);

        st_logInfo("Creating genome fragment for reference sequence: %s, start: %" PRIi64 ", length: %" PRIi64 "\n",
                    hmm->referenceName, hmm->refStart, hmm->refLength);

        // Run the forward-backward algorithm
        stRPHmm_forwardBackward(hmm);

        // Now compute a high probability path through the hmm
        stList *path = stRPHmm_forwardTraceBack(hmm);

        // Compute the genome fragment
        stGenomeFragment *gF = stGenomeFragment_construct(hmm, path);

        printGenotypeFragment(gF, false);


        // Write out VCF
        st_logInfo("\nWriting out VCF for fragment\n");

        // Write vcf relative to the reference fasta file
        // Get reference sequence
        str.l = 0; ksprintf(&str, "%s:%d-%d", gF->referenceName, gF->refStart, gF->refStart + gF->length + 1);
        int seq_len;
        char *seq = fai_fetch(fai, gF->referenceName, &seq_len);
        if ( seq_len < 0 ) {
            st_logCritical("Failed to fetch reference sequence %s in %s\n", str.s, referenceName);
            return EXIT_FAILURE; //todo close/free?
        }
        writeVcfFragment(vcfOutFP, hdr, gF, seq, referenceName, baseMapper);
        free(seq);


        // This one doesn't use the reference file when writing the vcf
        // It only writes the differences between the two haplotypes
         writeVcfNoReference(vcfOutFP2, hdr2, gF, baseMapper);


        // Optionally write out two BAMs, one for each read partition
        // st_logInfo("Writing out BAM partitions for fragment\n");
        /*
         * TODO: Optionally, write out a new BAM file expressing the partition (which reads belong in which partition)
         * Not sure if we need to write out multiple files or if we can add a per read flag to express this information.
         */
    }

    // Cleanup
    vcf_close(vcfOutFP);
    free(str.s);
    bcf_hdr_destroy(hdr);
    stList_destruct(hmms);

    //while(1); // Use this for testing for memory leaks

    return 0;
}

