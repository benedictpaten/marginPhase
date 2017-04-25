/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com) & Arthur Rand (arand@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <getopt.h>
#include <stdio.h>
#include <ctype.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <memory.h>
#include <htslib/faidx.h>
#include "stRPHmm.h"
#include "jsmn.h"

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
    fprintf(stderr, "-n --refSeqName : Name of reference sequence\n");
    fprintf(stderr, "-s --intervalStart : Starting position of interval to read\n");
    fprintf(stderr, "-e --intervalEnd : Ending position of interval to read\n");
}


bcf_hdr_t* writeVcfHeader(vcfFile *out, stList *genomeFragments) {
    bcf_hdr_t *hdr = bcf_hdr_init("w");
    kstring_t str = {0,0,NULL};

    // generic info
    str.l = 0;
    ksprintf(&str, "##marginPhase=htslib-%s\n",hts_version());
    bcf_hdr_append(hdr, str.s);

    // reference file used
    str.l = 0;
    ksprintf(&str, "##reference=file://%s\n", "TODO"); //todo
    bcf_hdr_append(hdr, str.s);

    // contigs
    // TODO: assert unique fragments, get full chrom length
    for(int64_t i=0; i<stList_length(genomeFragments); i++) {
        stRPHmm *hmm = stList_get(genomeFragments, i);
        str.l = 0;
        //ksprintf(&str, "##contig=<ID=%s,length=%d>\n", "chr1", 249250621);
        ksprintf(&str, "##contig=<ID=%s>\n", hmm->referenceName); //hmm->referenceName is the chrom
        bcf_hdr_append(hdr, str.s);
    }

    // formatting
    str.l = 0;
    ksprintf(&str, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_append(hdr, str.s);

    // samples
    bcf_hdr_add_sample(hdr, "SMPL1"); //todo
    bcf_hdr_add_sample(hdr, NULL);

    // write header
    bcf_hdr_write(out, hdr);

    // cleanup
    free(str.s);
    return hdr;
}

void writeVcfFragment(vcfFile *out, bcf_hdr_t *bcf_hdr, stGenomeFragment *gF, char *referenceSeq,
                      char *alphabet, char* wildcard) {
     // intialization
    bcf1_t *bcf_rec = bcf_init1();
    int32_t filter_info = bcf_hdr_id2int(bcf_hdr, BCF_DT_ID, "PASS"); //currently: everything passes
    int32_t *gt_info = (int*)malloc(bcf_hdr_nsamples(bcf_hdr)*2*sizeof(int)); //array specifying phasing
    kstring_t str = {0,0,NULL};

    // iterate over all positions
    for (int64_t i = 0; i < gF->length; ++i) {
        char refChar = referenceSeq[i];
        char refAlphChar = baseToAlphabet(refChar, alphabet, wildcard) - '0';
        st_logDebug("%d\t%c/%d:ref\t%8d:gs\t%.8f:gp\t%8d:hs1\t%.8f:hp1\t%8d:hs2\t%.8f:hp2\n",
                    i, refChar, refAlphChar,
                    gF->genotypeString[i], gF->genotypeProbs[i],
                    gF->haplotypeString1[i], gF->haplotypeProbs1[i],
                    gF->haplotypeString2[i], gF->haplotypeProbs2[i]);
        int h1AlphChar = gF->haplotypeString1[i];
        int h2AlphChar = gF->haplotypeString2[i];

        //prep
        bcf_clear1(bcf_rec);
        str.l = 0;
        int64_t pos = gF->refStart + i;

        // contig (CHROM)
        bcf_rec->rid = bcf_hdr_name2id(bcf_hdr, gF->referenceName); //defined in a contig in the top
        // POS
        bcf_rec->pos  = i + gF->refStart; // off by one?
        // ID - skip
        // QUAL - skip (TODO for now?)
        // REF
        kputc(refChar, &str);
        // ALT
        kputc(',', &str);
        kputc(h1AlphChar + '0', &str);
        kputc(',', &str);
        kputc(h2AlphChar + '0', &str);
        bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
        // FILTER
        bcf_update_filter(bcf_hdr, bcf_rec, &filter_info, 1);
        // FORMAT / $SMPL1
        gt_info[0] = bcf_gt_phased(0);
        gt_info[1] = bcf_gt_phased(1);
        bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);

        // save it
        bcf_write1(out, bcf_hdr, bcf_rec);
    }

    // cleanup
    free(str.s); bcf_destroy(bcf_rec);
}


int main(int argc, char *argv[]) {
    // Parameters / arguments
    char * logLevelString = NULL;
    char *bamFile = NULL;
    char *bamInFile = NULL;
    char *vcfOutFile = NULL;
    char *paramsFile = "params.json";
    char *referenceName = "hg19.chr3.fa";

    char *refSeqName = NULL;
    int32_t intervalStart = -1;
    int32_t intervalEnd = -1;

    // When done testing, set random seed
    // st_randomSeed();

    // Parse the options
    while (1) {
        static struct option long_options[] = {
                { "logLevel", required_argument, 0, 'a' },
                { "help", no_argument, 0, 'h' },
                { "bamFile", required_argument, 0, 'b'},
                { "vcfOutFile", required_argument, 0, 'v'},
                { "params", required_argument, 0, 'p'},
                { "refSeqName", required_argument, 0, 'n'},
                { "intervalStart", required_argument, 0, 's'},
                { "intervalEnd", required_argument, 0, 'e'},
                { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:v:p:n:s:e:h", long_options, &option_index);

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
        case 'n':
            refSeqName = stString_copy(optarg);
            break;
        case 's':
            intervalStart = atoi(optarg);
            break;
        case 'e':
            intervalEnd = atoi(optarg);
            break;
        default:
            usage();
            return 1;
        }
    }

    // Parse any model parameters
    st_logInfo("Parsing model parameters\n");

    char *alphabet[ALPHABET_SIZE];
    char *wildcard;
    stRPHmmParameters *params = parseParameters(paramsFile, alphabet, &wildcard);

    // Parse reads for interval
    st_logInfo("Parsing input reads\n");

    stList *profileSequences = stList_construct();
    parseReads(profileSequences, bamInFile, alphabet, wildcard, refSeqName, intervalStart, intervalEnd);



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

//        st_logDebug("*** Genome Fragment Information ***\n");
//        st_logDebug("Reference name: %s\n", gF->referenceName);
//        st_logDebug("Ref start: %d \n", gF->refStart);
//        st_logDebug("Length: %d \n", gF->length);
//
//        st_logDebug("Genotype string: %u\n", gF->genotypeString);
//        for (int64_t j = 0; j < gF->length; j++) {
//            st_logDebug("%u\t", gF->genotypeString[j]);
//        }
//        st_logDebug("Genotype Probabilities: \n");
//        for (int64_t j = 0; j < gF->length; j++) {
//            st_logDebug("%f \t", gF->genotypeProbs[j]);
//        }
//        st_logDebug("\nHaplotype 1 string: %u\n", gF->haplotypeString1);
//        for (int64_t j = 0; j < gF->length; j++) {
//            st_logDebug("%u\t", gF->haplotypeString1[j]);
//        }
//        st_logDebug("\nHaplotype 1 Probabilities: \n");
//        for (int64_t j = 0; j < gF->length; j++) {
//            st_logDebug("%f \t", gF->haplotypeProbs1[j]);
//        }
//        st_logDebug("\nHaplotype 2 string: %u\n", gF->haplotypeString2);
//        for (int64_t j = 0; j < gF->length; j++) {
//            st_logDebug("%u\t", gF->haplotypeString2[j]);
//        }
//        st_logDebug("\nHaplotype 2 Probabilities: \n");
//        for (int64_t j = 0; j < gF->length; j++) {
//            st_logDebug("%f \t", gF->haplotypeProbs2[j]);
//        }

        // Write out VCF
        st_logInfo("\nWriting out VCF for fragment\n");

        // Get reference sequence
        str.l = 0; ksprintf(&str, "%s:%d-%d", gF->referenceName, gF->refStart, gF->refStart + gF->length + 1);
        int seq_len;
        char *seq = fai_fetch(fai, str.s, &seq_len);
        printf(seq);
        if ( seq_len < 0 ) {
            st_logCritical("Failed to fetch reference sequence %s in %s\n", str.s, referenceName);
            return EXIT_FAILURE; //todo close/free?
        }

        // Write out VCF
        st_logInfo("Writing out VCF for fragment\n");
        writeVcfFragment(vcfOutFP, hdr, gF, seq, alphabet, wildcard);
        free(seq);


        // Optionally write out two BAMs, one for each read partition
        st_logInfo("Writing out BAM partitions for fragment\n");
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

