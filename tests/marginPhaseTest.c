/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "stRPHmm.h"

void printSequenceStats(FILE *fH, stList *profileSequences) {
    int64_t totalLength=0;
    for(int64_t i=0; i<stList_length(profileSequences); i++) {
        stProfileSeq *profileSeq = stList_get(profileSequences, i);
        totalLength += profileSeq->length;
    }
    fprintf(fH, "Got %" PRIi64 " profile sequences, with total length: %" PRIi64 ", average length: %f\n",
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
        assert(j >= 0 && j < length);
        uint64_t hapBase = haplotypeString[j];
        assert(j < ALPHABET_SIZE);

        // Expectation of a match
        totalExpectedMatches += getProb(&(profileSeq->profileProbs[i * ALPHABET_SIZE]), hapBase);
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
        fprintf(fH, "Base %" PRIi64 " count: %f fraction: %f\n", i, baseCounts[i], baseCounts[i]/totalCount);
    }
}

void genotypingTest(char *paramsFile, char *bamFile, char *vcfOutFile, char *vcfOutFileDiff, char *referenceFile) {
    fprintf(stderr, "Parsing parameters\n");
    stBaseMapper *baseMapper = stBaseMapper_construct();
    stRPHmmParameters *params = parseParameters(paramsFile, baseMapper);
    // Print a report of the parsed parameters
    stRPHmmParameters_printParameters(params, stderr);

    fprintf(stderr, "Creating profile sequences\n");
    stList *profileSequences = stList_construct();
    parseReads(profileSequences, bamFile, baseMapper);
    // Print some stats about the input sequences
    printSequenceStats(stderr, profileSequences);

    fprintf(stderr, "Building hmms\n");
    stList *hmms = getRPHmms(profileSequences, params);
    stList *l = stList_construct3(0, (void (*)(void *))stRPHmm_destruct2);
    while(stList_length(hmms) > 0) {
        stList_appendAll(l, stRPHMM_splitWherePhasingIsUncertain(stList_pop(hmms)));
    }
    hmms = l;
    fprintf(stderr, "Got %" PRIi64 " hmms\n", stList_length(hmms));

    fprintf(stderr, "Writing vcf files\n");
    vcfFile *vcfOutFP = vcf_open(vcfOutFile, "w");
    bcf_hdr_t *bcf_hdr = writeVcfHeader(vcfOutFP, l);

    vcfFile *vcfOutFP_diff = vcf_open(vcfOutFileDiff, "w");
    bcf_hdr_t *bcf_hdr_diff = writeVcfHeader(vcfOutFP_diff, l);


    for(int64_t i=0; i<stList_length(hmms); i++) {
        stRPHmm *hmm = stList_get(hmms, i);

        // Print stats about the HMM
        fprintf(stderr, "The %" PRIi64 "th hmm\n", i);
        stRPHmm_print(hmm, stderr, 0, 0);

        // Run the forward-backward algorithm
        stRPHmm_forwardBackward(hmm);

        // Now compute a high probability path through the hmm
        stList *path = stRPHmm_forwardTraceBack(hmm);

        // Compute the genome fragment
        stGenomeFragment *gF = stGenomeFragment_construct(hmm, path);

        // Get partitioned sequences
        stSet *reads1 = stRPHmm_partitionSequencesByStatePath(hmm, path, 1);
        stSet *reads2 = stRPHmm_partitionSequencesByStatePath(hmm, path, 0);
        fprintf(stderr, "There are %" PRIi64 " reads covered by the hmm, bipartitioned into sets of %" PRIi64 " and %" PRIi64 " reads\n",
                stList_length(hmm->profileSeqs), stSet_size(reads1), stSet_size(reads2));

        // Print the similarity between the two imputed haplotypes sequences
        fprintf(stderr, "The haplotypes have an %f identity\n", getIdentityBetweenHaplotypes(gF->haplotypeString1, gF->haplotypeString2, gF->length));

        // Print the base composition of the haplotype sequences
        double *hap1BaseCounts = getHaplotypeBaseComposition(gF->haplotypeString1, gF->length);
        fprintf(stderr, "The base composition of haplotype 1:\n");
        printBaseComposition(stderr, hap1BaseCounts);
        free(hap1BaseCounts);

        double *hap2BaseCounts = getHaplotypeBaseComposition(gF->haplotypeString2, gF->length);
        fprintf(stderr, "The base composition of haplotype 1:\n");
        printBaseComposition(stderr, hap2BaseCounts);
        free(hap2BaseCounts);

        // Print the base composition of the reads
        double *reads1BaseCounts =getExpectedProfileSequenceBaseComposition(reads1);
        fprintf(stderr, "The base composition of reads1 set:\n");
        printBaseComposition(stderr, reads1BaseCounts);
        free(reads1BaseCounts);

        double *reads2BaseCounts =getExpectedProfileSequenceBaseComposition(reads2);
        fprintf(stderr, "The base composition of reads2 set:\n");
        printBaseComposition(stderr, reads2BaseCounts);
        free(reads2BaseCounts);

        // Print some summary stats about the differences between haplotype sequences and the bipartitioned reads
        fprintf(stderr, "hap1 vs. reads1 identity: %f\n", getExpectedIdentity(gF->haplotypeString1, gF->refStart, gF->length, reads1));
        fprintf(stderr, "hap1 vs. reads2 identity: %f\n", getExpectedIdentity(gF->haplotypeString1, gF->refStart, gF->length, reads2));

        fprintf(stderr, "hap2 vs. reads2 identity: %f\n", getExpectedIdentity(gF->haplotypeString2, gF->refStart, gF->length, reads2));
        fprintf(stderr, "hap2 vs. reads1 identity: %f\n", getExpectedIdentity(gF->haplotypeString2, gF->refStart, gF->length, reads1));


        writeVcfFragment(vcfOutFP, bcf_hdr, gF, referenceFile, baseMapper, true);
        writeVcfFragment(vcfOutFP_diff, bcf_hdr_diff, gF, NULL, baseMapper, false);
    }

    vcf_close(vcfOutFP);
    vcf_close(vcfOutFP_diff);
}

void test_5kbGenotyping(CuTest *testCase) {

    fprintf(stderr, "Testing haplotype inference on NA12878.pb.chr3.5kb.bam\n");

    char *paramsFile = "../tests/params.json";
    char *bamFile = "../tests/NA12878.pb.chr3.5kb.bam";
    char *vcfOutFile = "test_5kb.vcf";
    char *vcfOutFileDiff = "test_5kb_diff.vcf";
    char *referenceFile = "../tests/hg19.chr3.9mb.fa";

    genotypingTest(paramsFile, bamFile, vcfOutFile, vcfOutFileDiff, referenceFile);
}

void test_100kbGenotyping(CuTest *testCase) {

    fprintf(stderr, "Testing haplotype inference on NA12878.pb.chr3.100kb.0.bam\n");

    char *paramsFile = "../tests/params.json";
    char *bamFile = "../tests/NA12878.pb.chr3.100kb.0.bam";
    char *vcfOutFile = "test_100kb.vcf";
    char *vcfOutFileDiff = "test_100kb_diff.vcf";
    char *referenceFile = "../tests/hg19.chr3.9mb.fa";

    genotypingTest(paramsFile, bamFile, vcfOutFile, vcfOutFileDiff, referenceFile);
}

CuSuite *marginPhaseTestSuite(void) {
    CuSuite* suite = CuSuiteNew();

    // System level tests
    SUITE_ADD_TEST(suite, test_5kbGenotyping);
    SUITE_ADD_TEST(suite, test_100kbGenotyping);

    return suite;
}
