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
    fprintf(stderr, "-a --logLevel   : Set the log level [default = info]\n");
    fprintf(stderr, "-h --help       : Print this help screen\n");
    fprintf(stderr, "-b --bamFile    : Input BAM file\n");
    fprintf(stderr, "-o --outSamBase : Output SAM Base (\"example\" -> \"example1.sam\", \"example2.sam\")\n");
    fprintf(stderr, "-v --vcfFile    : Output VCF file\n");
    fprintf(stderr, "-p --params     : Input params file\n");
    fprintf(stderr, "-r --reference  : Reference fasta file\n");
}

int main(int argc, char *argv[]) {
    // Parameters / arguments
    char * logLevelString = "info";
    char *bamInFile = NULL;
    char *samOutBase = "haplotype";
    char *vcfOutFile = "output.vcf";
    char *vcfOutFile2 = "output2.vcf";
    char *paramsFile = "params.json";
    char *referenceName = "hg19.chr3.fa";

    // TODO: When done testing, set random seed using st_randomSeed();

    // Parse the options
    while (1) {
        static struct option long_options[] = {
                { "logLevel", required_argument, 0, 'a' },
                { "help", no_argument, 0, 'h' },
                { "bamFile", required_argument, 0, 'b'},
                { "outSamBase", required_argument, 0, 'o'},
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
            break;
        case 'h':
            usage();
            return 0;
        case 'b':
            bamInFile = stString_copy(optarg);
            break;
        case 'o':
            samOutBase = stString_copy(optarg);
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
    st_setLogLevelFromString(logLevelString);

    // Parse any model parameters
    st_logInfo("> Parsing model parameters\n");
    stBaseMapper *baseMapper = stBaseMapper_construct();
    stRPHmmParameters *params = parseParameters(paramsFile, baseMapper);

    // Parse reads for interval
    st_logInfo("> Parsing input reads\n");
    stList *profileSequences = stList_construct();
    parseReads(profileSequences, bamInFile, baseMapper);

    // Create HMMs
    st_logInfo("> Creating read partitioning HMMs\n");
    stList *hmms = getRPHmms(profileSequences, params);

    // Break up the hmms where the phasing is uncertain
    st_logInfo("> Breaking apart HMMs where the phasing is uncertain\n");

    stList *l = stList_construct3(0, (void (*)(void *))stRPHmm_destruct2);
    while(stList_length(hmms) > 0) {
        stList_appendAll(l, stRPHMM_splitWherePhasingIsUncertain(stList_pop(hmms)));
    }
    hmms = l;
    st_logDebug("\tCreated %d hmms \n", stList_length(hmms));

    // Start VCF generation
    vcfFile *vcfOutFP = vcf_open(vcfOutFile, "w");
    bcf_hdr_t *hdr = writeVcfHeader(vcfOutFP, l);

    vcfFile *vcfOutFP2 = vcf_open(vcfOutFile2, "w");
    bcf_hdr_t *hdr2 = writeVcfHeader(vcfOutFP2, l);

    // Prep for BAM outputs
    stSet *haplotype1Ids = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, NULL);
    stSet *haplotype2Ids = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, NULL);

    // For each read partitioning HMM
    for(int64_t i=0; i<stList_length(hmms); i++) {
        stRPHmm *hmm = stList_get(hmms, i);

        st_logInfo("> Creating genome fragment for reference sequence: %s, start: %" PRIi64 ", length: %" PRIi64 "\n",
                    hmm->referenceName, hmm->refStart, hmm->refLength);

        // Run the forward-backward algorithm
        stRPHmm_forwardBackward(hmm);

        // Now compute a high probability path through the hmm
        stList *path = stRPHmm_forwardTraceBack(hmm);

        // Compute the genome fragment
        stGenomeFragment *gF = stGenomeFragment_construct(hmm, path);

        // Write out VCF
        st_logInfo("> Writing out VCF for fragment\n");

        // Write two vcfs, one using the reference fasta file and one not
        writeVcfFragment(vcfOutFP, hdr, gF, referenceName, baseMapper, true);
        writeVcfFragment(vcfOutFP2, hdr2, gF, NULL, baseMapper, false);

        // Get the reads which mapped to each path
        stList * haplotype1 = stSet_getList(stRPHmm_partitionSequencesByStatePath(hmm, path, true));
        for (int64_t j=0; j<stList_length(haplotype1); j++) {
            stSet_insert(haplotype1Ids, ((stProfileSeq *)stList_get(haplotype1, j))->readId);
        }

        stList * haplotype2 = stSet_getList(stRPHmm_partitionSequencesByStatePath(hmm, path, false));
        for (int64_t j=0; j<stList_length(haplotype2); j++) {
            stSet_insert(haplotype2Ids, ((stProfileSeq *)stList_get(haplotype2, j))->readId);
        }

    }

    // Write out two BAMs, one for each read partition
    st_logInfo("> Writing out SAM files for each partition\n");
    writeSplitSams(bamInFile, samOutBase, haplotype1Ids, haplotype2Ids);

    // Cleanup
    vcf_close(vcfOutFP);
    vcf_close(vcfOutFP2);
    bcf_hdr_destroy(hdr);
    stList_destruct(hmms);

    //while(1); // Use this for testing for memory leaks

    return 0;
}

