/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com) & Arthur Rand (arand@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <ctype.h>
#include <htslib/vcf.h>
#include "stRPHmm.h"



void writeSplitBams(char *bamInFile, char *bamOutBase,
                    stSet *haplotype1Ids, stSet *haplotype2Ids) {
    // prep
    char *haplotype1BamOutFile = stString_print("%s.1.sam", bamOutBase);
    char *haplotype2BamOutFile = stString_print("%s.2.sam", bamOutBase);
    char *unmatchedBamOutFile = stString_print("%s.filtered.sam", bamOutBase);

    // file management
    samFile *in = hts_open(bamInFile, "r");
    if (in == NULL) {
        st_errAbort("ERROR: Cannot open bam file %s\n", bamInFile);
    }
    bam_hdr_t *bamHdr = sam_hdr_read(in);
    bam1_t *aln = bam_init1();

    int r;
    st_logDebug("Writing haplotype output to: %s and %s \n", haplotype1BamOutFile, haplotype2BamOutFile);
    samFile *out1 = hts_open(haplotype1BamOutFile, "w");
    r = sam_hdr_write(out1, bamHdr);

    samFile *out2 = hts_open(haplotype2BamOutFile, "w");
    r = sam_hdr_write(out2, bamHdr);

    samFile *out3 = hts_open(unmatchedBamOutFile, "w");
    r = sam_hdr_write(out3, bamHdr);

    // read in input file, write out each read to one sam file
    int32_t readCountH1 = 0;
    int32_t readCountH2 = 0;
    int32_t readCountFiltered = 0;
    while(sam_read1(in,bamHdr,aln) > 0) {

        char *readName = bam_get_qname(aln);
        if (stSet_search(haplotype1Ids, readName) != NULL) {
            r = sam_write1(out1, bamHdr, aln);
            readCountH1++;
        } else if (stSet_search(haplotype2Ids, readName) != NULL) {
            r = sam_write1(out2, bamHdr, aln);
            readCountH2++;
        } else {
            r = sam_write1(out3, bamHdr, aln);
            readCountFiltered++;
        }
    }
    st_logInfo("Read counts:\n\thap1: %d\thap2: %d\tfiltered out: %d \n", readCountH1, readCountH2, readCountFiltered);

    // Cleanup
    bam_destroy1(aln);
    bam_hdr_destroy(bamHdr);
    sam_close(in);
    sam_close(out1);
    sam_close(out2);
    sam_close(out3);
    free(haplotype1BamOutFile);
    free(haplotype2BamOutFile);
    free(unmatchedBamOutFile);
}


bcf_hdr_t* writeVcfHeader(vcfFile *out, stList *genomeFragments, char *referenceName) {
    bcf_hdr_t *hdr = bcf_hdr_init("w");
    kstring_t str = {0,0,NULL};

    // generic info
    str.l = 0;
    ksprintf(&str, "##marginPhase=htslib-%s\n",hts_version());
    bcf_hdr_append(hdr, str.s);

    // reference file used
    str.l = 0;
    ksprintf(&str, "##reference=file://%s\n", referenceName);
    bcf_hdr_append(hdr, str.s);

    // contigs
    // TODO: assert unique fragments, get full chrom length
    for(int64_t i=0; i<stList_length(genomeFragments); i++) {
        stRPHmm *hmm = stList_get(genomeFragments, i);
        str.l = 0;
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

// This function writes out a vcf for the two haplotypes
// It optionally writes it relative to a reference fasta file or
// writes it for one of the haplotypes relative to the other
void writeVcfFragment(vcfFile *out, bcf_hdr_t *bcf_hdr, stGenomeFragment *gF,
                      char *referenceName, stBaseMapper *baseMapper, bool differencesOnly) {

    char *referenceSeq;
    // Get reference (needed for VCF generation)
    faidx_t *fai = fai_load(referenceName);
    if ( !fai ) {
        st_logCritical("Could not load fai index of %s.  Maybe you should run 'samtools faidx %s'\n",
                       referenceName, referenceName);
        return;
    }
    int seq_len;
    referenceSeq = fai_fetch(fai, gF->referenceName, &seq_len);
    if ( seq_len < 0 ) {
        st_logCritical("Failed to fetch reference sequence %s\n", referenceName);
        return;
    }
    fai_destroy(fai);


    // intialization
    bcf1_t *bcf_rec = bcf_init1();
    int32_t filter_info = bcf_hdr_id2int(bcf_hdr, BCF_DT_ID, "PASS"); //currently: everything passes
    int32_t *gt_info = (int*)malloc(bcf_hdr_nsamples(bcf_hdr)*2*sizeof(int)); //array specifying phasing
    kstring_t str = {0,0,NULL};

    // iterate over all positions
    for (int64_t i = 0; i < gF->length; i++) {

        int h1AlphVal = gF->haplotypeString1[i];
        int h2AlphVal = gF->haplotypeString2[i];
        char h1AlphChar = stBaseMapper_getCharForValue(baseMapper, h1AlphVal);
        char h2AlphChar = stBaseMapper_getCharForValue(baseMapper, h2AlphVal);
        float h1Prob = gF->haplotypeProbs1[i];
        float h2Prob = gF->haplotypeProbs2[i];
        float genotypeProb = gF->genotypeProbs[i];
        float genotype = (float)gF->genotypeString[i];

        //prep
        bcf_clear1(bcf_rec);
        str.l = 0;

        // contig (CHROM)
        bcf_rec->rid = bcf_hdr_name2id(bcf_hdr, gF->referenceName); //defined in a contig in the top
        // POS
        bcf_rec->pos  = i + gF->refStart - 1; // off by one?
        // ID - skip
        // QUAL - currently writing out the genotype probability
        bcf_rec->qual = genotypeProb;

        // Get phasing info
        gt_info[0] = bcf_gt_phased(0);
        gt_info[1] = bcf_gt_phased(1);

        char refChar = toupper(referenceSeq[i + gF->refStart-1]);
        if (!differencesOnly) {
            kputc(refChar, &str); // REF
            kputc(',', &str);
            kputc(h1AlphChar, &str);
            kputc(',', &str);
            kputc(h2AlphChar, &str);

            bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
            // FORMAT / $SMPL1
            bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);
            bcf_write1(out, bcf_hdr, bcf_rec);

        } else {
            if (i + 1 >= gF->length) break;

            int next_h1AlphVal = gF->haplotypeString1[i+1];
            int next_h2AlphVal = gF->haplotypeString2[i+1];
            char next_h1AlphChar = stBaseMapper_getCharForValue(baseMapper, next_h1AlphVal);
            char next_h2AlphChar = stBaseMapper_getCharForValue(baseMapper, next_h2AlphVal);
            char nextRefChar = toupper(referenceSeq[i + 1 + gF->refStart-1]);

            if (next_h1AlphChar != next_h2AlphChar && (h1AlphChar != '-' && h2AlphChar != '-')) {
                // Check to see if there was an insertion or deletion in the next spot
                if (next_h1AlphChar == '-' || next_h2AlphChar == '-') {
                    kstring_t refstr = {0,0,NULL};
                    kputc(refChar, &refstr);
                    kputc(nextRefChar, &refstr);
                    kputc(h1AlphChar, &str);
                    // Count how many positions are involved in indel
                    int j = 2;
                    int genotype1length = 1;
                    int genotype2length = 1;
                    while(i+j < gF->length &&
                            (stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString1[i+j]) == '-' ||
                            stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[i+j]) == '-')) {

                        nextRefChar = toupper(referenceSeq[i + j + gF->refStart-1]);
                        kputc(nextRefChar, &str);
                        if (stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString1[i+j]) != '-') {
                            genotype1length++;
                        }
                        if (stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[i+j]) != '-') {
                            genotype2length++;
                        }
                        j++;
                    }
                    char *genotype1 = st_calloc(genotype1length, sizeof(char));
                    char *genotype2 = st_calloc(genotype2length, sizeof(char));
                    genotype1[0] = h1AlphChar;
                    genotype2[0] = h2AlphChar;
                    int genotype1index = 1;
                    int genotype2index = 1;

                    // Add the bases from the first haplotype
                    for (int64_t k = 1; k < j; k++) {
                        if (stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString1[i+k]) != '-') {
                            kputc(stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString1[i+k]), &str);
                            genotype1[genotype1index] = stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString1[i+k]);
                            genotype1index++;
                        }
                    }

                    // Go through the second haplotype
                    for (int64_t k = 1; k < j; k++) {
                        if (stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[i+k]) != '-') {
                            genotype2[genotype2index] = stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[i+k]);
                            genotype2index++;
                        }
                    }
                    if (strcmp(genotype1, genotype2) == 0) {
                        // Homozygous alleles
                        // Ref allele will be the reference string constructed
                        gt_info[0] = bcf_gt_phased(1);
                        gt_info[1] = bcf_gt_phased(1);
                        // Add the bases from the second haplotype
                        kputc(',', &refstr);
                        kputc(h2AlphChar, &refstr);
                        for (int64_t k = 1; k < j; k++) {
                            if (stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[i+k]) != '-') {
                                kputc(stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[i+k]), &refstr);
                            }
                        }
                        bcf_update_alleles_str(bcf_hdr, bcf_rec, refstr.s);
                    } else {
                        // Add the bases from the second haplotype
                        kputc(',', &str);
                        kputc(h2AlphChar, &str);
                        for (int64_t k = 1; k < j; k++) {
                            if (stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[i+k]) != '-') {
                                kputc(stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[i+k]), &str);
                            }
                        }
                        bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
                    }
                    i += (j-1);
                    bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);
                    bcf_write1(out, bcf_hdr, bcf_rec);
                }
                else if (h1AlphChar != h2AlphChar){
                    kputc(h1AlphChar, &str);
                    kputc(',', &str);
                    kputc(h2AlphChar, &str);
                    bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
                    bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);
                    bcf_write1(out, bcf_hdr, bcf_rec);
                }
            }
            else if (h1AlphChar != h2AlphChar && (h1AlphChar != '-' && h2AlphChar != '-')) {
                kputc(h1AlphChar, &str);
                kputc(',', &str);
                kputc(h2AlphChar, &str);
                bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
                bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);
                bcf_write1(out, bcf_hdr, bcf_rec);

            }
            else if (next_h1AlphChar == '-' && next_h2AlphChar == '-') {
                // Cases where both strands might have gaps relative to the reference
                kstring_t refstr = {0,0,NULL};
                kputc(refChar, &refstr);
                kputc(nextRefChar, &refstr);
                kputc(h1AlphChar, &str);
                // Count how many positions are involved in indel
                int j = 2;
                int genotype1length = 1;
                int genotype2length = 1;
                while(i+j < gF->length &&
                      (stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString1[i+j]) == '-' ||
                       stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[i+j]) == '-')) {
                    nextRefChar = toupper(referenceSeq[i + j + gF->refStart-1]);
                    kputc(nextRefChar, &str);
                    if (stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString1[i+j]) != '-') {
                        genotype1length++;
                    }
                    if (stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[i+j]) != '-') {
                        genotype2length++;
                    }
                    j++;
                }
                char *genotype1 = st_calloc(genotype1length, sizeof(char));
                char *genotype2 = st_calloc(genotype2length, sizeof(char));
                genotype1[0] = h1AlphChar;
                genotype2[0] = h2AlphChar;
                int genotype1index = 1;
                int genotype2index = 1;
                // Add the bases from the first haplotype
                for (int64_t k = 1; k < j; k++) {
                    if (stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString1[i+k]) != '-') {
                        kputc(stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString1[i+k]), &str);
                        genotype1[genotype1index] = stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString1[i+k]);
                        genotype1index++;
                    }
                }
                // Go through the second haplotype
                for (int64_t k = 1; k < j; k++) {
                    if (stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[i+k]) != '-') {
                        genotype2[genotype2index] = stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[i+k]);
                        genotype2index++;
                    }
                }
//                if (h1AlphChar == h2AlphChar) {
                if (strcmp(genotype1, genotype2) == 0) {
                    // Homozygous alleles
                    // Ref allele will be the reference string constructed
                    gt_info[0] = bcf_gt_phased(1);
                    gt_info[1] = bcf_gt_phased(1);
                    // Add the bases from the second haplotype
                    kputc(',', &refstr);
                    kputc(h2AlphChar, &refstr);
                    for (int64_t k = 1; k < j; k++) {
                        if (stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[i+k]) != '-') {
                            kputc(stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[i+k]), &refstr);
                        }
                    }
                    bcf_update_alleles_str(bcf_hdr, bcf_rec, refstr.s);
                } else {
                    // Add the bases from the second haplotype
                    kputc(',', &str);
                    kputc(h2AlphChar, &str);
                    for (int64_t k = 1; k < j; k++) {
                        if (stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[i+k]) != '-') {
                            kputc(stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[i+k]), &str);
                        }
                    }
                    bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
                }
                i += (j-1);

                bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);
                bcf_write1(out, bcf_hdr, bcf_rec);

            } else if (h1AlphChar != refChar || h2AlphChar != refChar) {
                // Doesn't match the reference
                gt_info[0] = bcf_gt_phased(1);
                gt_info[1] = bcf_gt_phased(1);
                kputc(refChar, &str);
                kputc(',', &str);
                kputc(h2AlphChar, &str);
                bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
                bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);
                bcf_write1(out, bcf_hdr, bcf_rec);
            }
        }
    }

    // cleanup
    free(str.s);
    free(gt_info);
    if (referenceSeq) free(referenceSeq);
    bcf_destroy(bcf_rec);
}




/*
 * Prints information contained in genotypeResults struct.
 */
void printGenotypeResults(stGenotypeResults *results) {
    // Variables
    int64_t hetsInRefIndels = results->hetsInRef_Insertions + results->hetsInRef_Deletions;
    int64_t homozygousVariantsInRefIndels = results->homozygousVariantsInRef_Insertions + results->homozygousVariantsInRef_Deletions;

    // Sensitivity
    st_logInfo("\nSensitivity / recall (fraction of true positives compared to reference): \n");
    st_logInfo("\tOverall: %.4f \t\t\t(%d out of %d)\n",
               (float)results->truePositives/results->positives,
               results->truePositives,
               results->positives);
    st_logInfo("\tNo indels: %.4f \t\t\t(%d out of %d)\n",
               (float) (results->truePositives-results->truePositiveIndels-results->truePositiveHomozygousIndels)/
                       (results->positives-hetsInRefIndels-homozygousVariantsInRefIndels),
               results->truePositives-results->truePositiveIndels-results->truePositiveHomozygousIndels,
               results->positives-hetsInRefIndels-homozygousVariantsInRefIndels);
    st_logInfo("\t\tHets: %.4f \t\t\t(%d out of %d)\n",
               (float) (results->truePositiveHet - results->truePositiveHetIndels)
               / (results->hetsInRef-hetsInRefIndels),
               results->truePositiveHet - results->truePositiveHetIndels,
               results->hetsInRef-hetsInRefIndels);
    st_logInfo("\t\tHomozygous alts: %.4f \t(%d out of %d)\n",
               (float) (results->truePositiveHomozygous - results->truePositiveHomozygousIndels)
               /(results->homozygousVariantsInRef-homozygousVariantsInRefIndels),
               results->truePositiveHomozygous - results->truePositiveHomozygousIndels,
               results->homozygousVariantsInRef - homozygousVariantsInRefIndels);
    st_logInfo("\tIndels only: %.4f \t\t\t(%d out of %d)\n",
               (float)results->truePositiveIndels/
                       (homozygousVariantsInRefIndels + hetsInRefIndels),
               results->truePositiveIndels,
               hetsInRefIndels + homozygousVariantsInRefIndels);
    st_logInfo("\t\tHets: %.4f \t\t\t(%d out of %d)\n",
               (float) results->truePositiveHetIndels / hetsInRefIndels,
               results->truePositiveHetIndels,
               hetsInRefIndels);
    st_logInfo("\t\tHomozygous alts: %.4f \t(%d out of %d)\n",
               (float) results->truePositiveHomozygousIndels/homozygousVariantsInRefIndels,
               results->truePositiveHomozygousIndels,
               homozygousVariantsInRefIndels);

    // Specificity
    st_logInfo("\nSpecificity (fraction of true negatives compared to reference): \n");
    st_logInfo("\tOverall: %.4f \t\t(%" PRIi64 " out of %"PRIi64 ")\n",
               (float)results->trueNegatives/results->negatives,
               results->trueNegatives, results->negatives);

    // Precision
    st_logInfo("\nPrecision (fraction of true positives compared to all calls):\n");
    st_logInfo("\tOverall: %.4f \t\t(%d out of %d)\n",
               (float)results->truePositives/(results->truePositives+results->falsePositives),
               results->truePositives, results->truePositives+results->falsePositives);
    st_logInfo("\tNo indels: %.4f \t\t(%d out of %d)\n",
               (float) (results->truePositives-results->truePositiveIndels) /
                       (results->truePositives + results->falsePositives - results->truePositiveIndels - results->falsePositiveIndels),
               (results->truePositives - results->truePositiveIndels),
               (results->truePositives + results->falsePositives - results->truePositiveIndels - results->falsePositiveIndels));
    st_logInfo("\tIndels only: %.4f \t\t(%d out of %d)\n",
               (float)results->truePositiveIndels / (results->truePositiveIndels + results->falsePositiveIndels),
               results->truePositiveIndels, results->truePositiveIndels + results->falsePositiveIndels);

    // False positives
    st_logInfo("\nFalse positives:\n");
    st_logInfo("\tOverall: %" PRIi64 "\n", results->falsePositives);
    st_logInfo("\tNo indels: %" PRIi64 "\n", results->falsePositives - results->falsePositiveIndels);
    st_logInfo("\tIndels only: %" PRIi64 "\n", results->falsePositiveIndels);


    // More detailed numbers about errors
    st_logInfo("\nFalse negatives:\n");
    st_logInfo("\tOverall: %" PRIi64 "\n", results->falseNegatives);
    st_logInfo("\tMissed hets - SNV: %" PRIi64 "\n", results->error_missedHet -
                       results->error_missedHet_Deletions - results->error_missedHet_Insertions);
    st_logInfo("\tMissed hets - Indels: %" PRIi64 "\n",
               results->error_missedHet_Insertions + results->error_missedHet_Deletions);
    st_logInfo("\t\tInsertions: %d \t(out of %d)\n",
               results->error_missedHet_Insertions, results->hetsInRef_Insertions);
    st_logInfo("\t\tDeletions: %d \t(out of %d)\n",
               results->error_missedHet_Deletions, results->hetsInRef_Deletions);
    st_logInfo("\tIncorrect homozygous variants - SNV: %" PRIi64 "\n", results->error_homozygousInRef -
                       results->error_homozygous_Insertions - results->error_homozygous_Deletions);
    st_logInfo("\tIncorrect homozygous variants - Indels: %" PRIi64 "\n",
               results->error_homozygous_Insertions + results->error_homozygous_Deletions);
    st_logInfo("\t\tInsertions: %d \t(out of %d)\n",
               results->error_homozygous_Insertions, results->homozygousVariantsInRef_Insertions);
    st_logInfo("\t\tDeletions: %d \t(out of %d)\n",
               results->error_homozygous_Deletions, results->homozygousVariantsInRef_Deletions);

    // Phasing
    st_logInfo("\nPhasing:\n");
    st_logInfo("\tSwitch error rate: %.4f \t (%" PRIi64 " out of %"PRIi64 ", ", (float)results->switchErrors/(results->truePositives-results->uncertainPhasing), results->switchErrors, results->truePositives-results->uncertainPhasing);
    st_logInfo("fraction correct: %.2f)\n", 1.0 - (float)results->switchErrors/(results->truePositives-results->uncertainPhasing));
    st_logInfo("\tAverage distance between switch errors: %.3f\n", results->switchErrorDistance);
    st_logInfo("\tUncertain phasing spots: %" PRIi64 "\n\n", results->uncertainPhasing);
}


void writeParamFile(char *outputFilename, stRPHmmParameters *params) {
    // get file
    FILE *fd = fopen(outputFilename, "w");
    if (fd == NULL) {
        st_logCritical("Failed to open output param file '%s'. No file will be written\n", outputFilename);
        return;
    }
    //for whether to print the last comma
    int64_t noCommaIdx = (ALPHABET_SIZE) * (ALPHABET_SIZE) - 1;

    fprintf(fd, "{\n");
    fprintf(fd, "  \"alphabet\" : [ \"Aa\", \"Cc\", \"Gg\", \"Tt\", \"-\" ],\n");
    fprintf(fd, "    \n");
    fprintf(fd, "  \"wildcard\" : \"Nn\",\n");
    fprintf(fd, "    \n");
    fprintf(fd, "  \"haplotypeSubstitutionModel\" : [ \n");
    for (int64_t i = 0; i < ALPHABET_SIZE; i++) {
        fprintf(fd, "    ");
        for (int64_t j = 0; j < ALPHABET_SIZE; j++) {
            int64_t idx = i * ALPHABET_SIZE + j;
            fprintf(fd, "%8f", exp(params->hetSubModelSlow[idx]));
            if (idx != noCommaIdx) fprintf(fd, ", ");
        }
        fprintf(fd, "\n");
    }
    fprintf(fd, "   ],\n");
    fprintf(fd, "    \n");
    fprintf(fd, "  \"readErrorModel\" : [ \n");
    for (int64_t i = 0; i < ALPHABET_SIZE; i++) {
        fprintf(fd, "    ");
        for (int64_t j = 0; j < ALPHABET_SIZE; j++) {
            int64_t idx = i * ALPHABET_SIZE + j;
            fprintf(fd, "%8f", exp(params->readErrorSubModelSlow[idx]));
            if (idx != noCommaIdx) fprintf(fd, ", ");
        }
        fprintf(fd, "\n");
    }
    fprintf(fd, "   ],\n");
    fprintf(fd, "    \n");
    fprintf(fd, "  \"maxNotSumTransitions\" : %s,\n", params->maxNotSumTransitions ? "true" : "false");
    fprintf(fd, "    \n");
    fprintf(fd, "  \"maxPartitionsInAColumn\" : %" PRIi64 ",\n", params->maxPartitionsInAColumn);
    fprintf(fd, "\n");
    fprintf(fd, "  \"maxCoverageDepth\" : %" PRIi64 ",\n", params->maxCoverageDepth);
    fprintf(fd, "\n");
    fprintf(fd, "  \"minReadCoverageToSupportPhasingBetweenHeterozygousSites\" : %" PRIi64 ",\n", params->minReadCoverageToSupportPhasingBetweenHeterozygousSites);
    fprintf(fd, "  \n");
    fprintf(fd, "  \"onDiagonalReadErrorPseudoCount\" : %f,\n", params->onDiagonalReadErrorPseudoCount);
    fprintf(fd, "  \n");
    fprintf(fd, "  \"offDiagonalReadErrorPseudoCount\" : %f,\n", params->offDiagonalReadErrorPseudoCount);
    fprintf(fd, "  \n");
    fprintf(fd, "  \"trainingIterations\" : %" PRIi64 ",\n", params->trainingIterations);
    fprintf(fd, "  \n");
    int64_t verbosityBitstring = 0
        | (params->verboseTruePositives ? LOG_TRUE_POSITIVES : 0);
    fprintf(fd, "  \"verbose\" : %" PRIi64 ",\n", verbosityBitstring);
    fprintf(fd, "}");

    if (fclose(fd) != 0) st_logCritical("Failed to close output param file: %s\n", outputFilename);
}
