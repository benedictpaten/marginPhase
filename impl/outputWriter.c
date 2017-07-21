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
void writeVcfFragment(vcfFile *out, bcf_hdr_t *bcf_hdr, stGenomeFragment *gF, char *referenceName, stBaseMapper *baseMapper, bool differencesOnly) {

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

        int64_t h1AlphVal = gF->haplotypeString1[i];
        int64_t h2AlphVal = gF->haplotypeString2[i];
        char h1AlphChar = stBaseMapper_getCharForValue(baseMapper, h1AlphVal);
        char h2AlphChar = stBaseMapper_getCharForValue(baseMapper, h2AlphVal);
        float genotypeProb = gF->genotypeProbs[i];

        // TODO: not currently being used
        float h1Prob = gF->haplotypeProbs1[i];
        float h2Prob = gF->haplotypeProbs2[i];
        float genotype = (float)gF->genotypeString[i];

        //prep
        bcf_clear1(bcf_rec);
        str.l = 0;

        // contig (CHROM)
        bcf_rec->rid = bcf_hdr_name2id(bcf_hdr, gF->referenceName); //defined in a contig in the top
        // POS
//        bcf_rec->pos  = i + gF->refStart - 1; // off by one?
//        int64_t recordPos = gF->refCoords[i] - 1;
        int64_t recordPos = gF->refCoords[i]-1;
        int64_t gapSize = gapSizeAtIndex(gF->refCoords, i);

//        int64_t recordPos = i + gF->refStart - 1;
        bcf_rec->pos = recordPos;
        // ID - skip
        // QUAL - currently writing out the genotype probability
        bcf_rec->qual = genotypeProb;

        // Get phasing info
        gt_info[0] = bcf_gt_phased(0);
        gt_info[1] = bcf_gt_phased(1);



        char refChar = toupper(referenceSeq[recordPos]);
        if (!differencesOnly) {
            // TODO: What to do on this where there is a gap?
            if (gF->refCoords[i] != gF->refCoords[i-1]) {
                kputc(refChar, &str); // REF
                kputc(',', &str);
                kputc(h1AlphChar, &str);
                kputc(',', &str);
                kputc(h2AlphChar, &str);

                bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
                // FORMAT / $SMPL1
                bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);
                bcf_write1(out, bcf_hdr, bcf_rec);
            }

        } else {
            if (i + 1 >= gF->length) break;

            int64_t next_h1AlphVal = gF->haplotypeString1[i+1];
            int64_t next_h2AlphVal = gF->haplotypeString2[i+1];
            char next_h1AlphChar = stBaseMapper_getCharForValue(baseMapper, next_h1AlphVal);
            char next_h2AlphChar = stBaseMapper_getCharForValue(baseMapper, next_h2AlphVal);

            if (next_h1AlphChar != next_h2AlphChar) {
                if (h1AlphChar != '-' && h2AlphChar != '-') {
                    // Check to see if there was an insertion or deletion in the next spot
                    if (next_h1AlphChar == '-' && next_h2AlphChar != '-') {
                        kputc(h1AlphChar, &str);
                        kputc(',', &str);
                        kputc(h2AlphChar, &str);
                        kputc(next_h2AlphChar, &str);
                        int j = 2;
                        while(i+j < gF->length
                              && stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString1[i+j]) == '-'
                              && stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[i+j]) != '-') {
                            kputc(stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[i+j]), &str);
                            j++;
                        }
                        bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
                        bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);
                        bcf_write1(out, bcf_hdr, bcf_rec);
                    } else if (next_h2AlphChar == '-' && next_h1AlphChar != '-') {
                        kputc(h1AlphChar, &str);
                        kputc(next_h1AlphChar, &str);
                        int j = 2;
                        while(i+j < gF->length &&
                                stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[i+j]) == '-'
                              && stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString1[i+j]) != '-') {
                            kputc(stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString1[i+j]), &str);
                            j++;
                        }
                        kputc(',', &str);
                        kputc(h2AlphChar, &str);
                        bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
                        bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);
                        bcf_write1(out, bcf_hdr, bcf_rec);

                    } else if (h1AlphChar != h2AlphChar){
                        kputc(h1AlphChar, &str);
                        kputc(',', &str);
                        kputc(h2AlphChar, &str);
                        bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
                        bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);
                        bcf_write1(out, bcf_hdr, bcf_rec);
                    }
                }
            }
            else if (h1AlphChar != h2AlphChar && (h1AlphChar != '-' && h2AlphChar != '-')) {
                // Could also list things that don't match the reference if
                // h1AlphChar != refChar || h2AlphChar != refChar
                kputc(h1AlphChar, &str);
                kputc(',', &str);
                kputc(h2AlphChar, &str);
                bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
                bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);
                bcf_write1(out, bcf_hdr, bcf_rec);
            }
        }
//        i += gapSize;
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
    // Sensitivity
    st_logInfo("\nSensitivity: %f, \t without indels: %f \n\t(= fraction of true positives compared to reference, \t%"
                       PRIi64 " out of %"PRIi64 " / %" PRIi64 " out of %" PRIi64 ")\n",
               (float)results->truePositives/results->positives,
               (float) (results->truePositives-results->truePositiveGaps)/(results->positives-results->truePositiveGaps-results->error_missedIndels),
               results->truePositives, results->positives,
               results->truePositives-results->truePositiveGaps, results->positives-results->truePositiveGaps-results->error_missedIndels) ;
    st_logInfo("\tHomozygous variants in sample: %" PRIi64 "\n",
               results->error_homozygousInRef);
    st_logInfo("\tFalse negatives: %" PRIi64 "\n", results->falseNegatives);

    // Specificity
    st_logInfo("\nSpecificity: %f \n\t(= fraction of true negatives compared to reference, \t%"
                       PRIi64 " out of % "PRIi64 ")\n",
               (float)results->trueNegatives/results->negatives,
               results->trueNegatives, results->negatives);
    st_logInfo("\tIncorrect positives: %" PRIi64 "\n", results->error_incorrectVariant);
    st_logInfo("\tFalse positives: %" PRIi64 ",\twith gaps: %" PRIi64 "\n", results->falsePositives, results->falsePositiveGaps);


    // More detailed numbers about errors
    st_logInfo("\nFalse negatives:\n");
    st_logInfo("\tPartition bad: %" PRIi64 " \t\t(%f)\n",
               results->error_badPartition, (float)results->error_badPartition/results->falseNegatives);
    st_logInfo("\tIndel missed: %" PRIi64 " \t\t(%f)\n",
               results->error_missedIndels, (float)results->error_missedIndels/results->falseNegatives);

    // Phasing
    st_logInfo("\nPhasing:\n");
    st_logInfo("\tSwitch error rate: %f \t (%" PRIi64 " out of %"PRIi64 ", ", (float)results->switchErrors/(results->truePositives-results->uncertainPhasing), results->switchErrors, results->truePositives-results->uncertainPhasing);
    st_logInfo("fraction correct: %f)\n", 1.0 - (float)results->switchErrors/(results->truePositives-results->uncertainPhasing));
    st_logInfo("\tAverage distance between switch errors: %f\n", results->switchErrorDistance);
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

void printFalsePositive(bcf1_t *unpackedRecord, int64_t evalPos, stRPHmm *hmm) {

    char *evalRefChar = unpackedRecord->d.als;
    char *evalAltChar = unpackedRecord->d.allele[1];

    st_logDebug("  pos: %" PRIi64 " \n  ref:%s alt: %s \n",
                evalPos, evalRefChar, evalAltChar);
    printColumnAtPosition(hmm, evalPos);

    if (strlen(evalRefChar) > 1 || strlen(evalAltChar) > 1) {
        size_t indelLen = strlen(evalRefChar);
        if (strlen(evalAltChar) > indelLen) indelLen = strlen(evalAltChar);
        for (int64_t j = 1; j < indelLen; j++) {
            st_logDebug("\tNext pos: %" PRIi64 "\n", evalPos+j);
            printColumnAtPosition(hmm, evalPos+j);
        }
    }
}

void printAlleleInfo(bcf1_t *unpackedRecordRef, stRPHmm *hmm, int64_t referencePos, char *refChar, char *evalRefChar, char *evalAltChar) {
    st_logDebug("  pos: %" PRIi64 "\n  ref: %s   alt: ", referencePos, refChar);
    for (int i = 1; i < unpackedRecordRef->n_allele; i++) {
        if (i != 1) st_logDebug(",");
        st_logDebug("%s", unpackedRecordRef->d.allele[i]);
        st_logDebug("\n  output: %s, %s\n", evalRefChar, evalAltChar);
        printColumnAtPosition(hmm, referencePos);
    }
}

void printPartitionInfo(int64_t referencePos, int64_t evalPos, stSet *reads1, stSet *reads2, stGenomeFragment *gF, stRPHmm *hmm) {
    // print additional partition info
    double *read1BaseCounts = getProfileSequenceBaseCompositionAtPosition(reads1, referencePos);
    double *read2BaseCounts = getProfileSequenceBaseCompositionAtPosition(reads2, referencePos);
    st_logDebug("\tPartition 1: \n");
    printBaseComposition2(read1BaseCounts);
    st_logDebug("\tPartition 2: \n");
    printBaseComposition2(read2BaseCounts);
    int64_t *gFIndex = stHash_search(gF->refCoordMap, &evalPos);
    st_logDebug("\tposterior prob: %f\n", gF->genotypeProbs[*gFIndex]);
}

void recordHomozygousVariant(stGenotypeResults *results) {
    results->falsePositives++;
    results->error_homozygousInRef++;
    results->positives--;
    results->negatives++;
}

void recordTruePositive(stGenotypeResults *results, stRPHmmParameters *params, bcf1_t *unpackedRecordRef, stBaseMapper *baseMapper, stRPHmm *hmm, stGenomeFragment *gF, stSet *reads1, stSet *reads2) {
    // Get all of these again
    char *refChar = unpackedRecordRef->d.als;
    char *refAltChar = unpackedRecordRef->d.allele[1];
    int64_t referencePos = unpackedRecordRef->pos+1;
    int64_t *gFIndex = stHash_search(gF->refCoordMap, &referencePos);
    int h1AlphChar = stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString1[*gFIndex]);
    int h2AlphChar = stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[*gFIndex]);

    results->truePositives++;
    if (strlen(refChar) > 1 || strlen(refAltChar) > 1) results->truePositiveGaps++;
    if (params->verboseTruePositives) {
        st_logDebug("\nTRUE POSITIVE\n");
        printAlleleInfo(unpackedRecordRef, hmm, referencePos, refChar, (char *) &h1AlphChar, (char *) &h2AlphChar);
        printPartitionInfo(referencePos, referencePos, reads1, reads2, gF, hmm);
    }
}

/*
     * Test to compare a vcf to a truth vcf containing known variants for the region.
     *
     * Test depends on the format of the vcf files written in vcfWriter.c
     * (Currently they don't follow a quite standard format)
     *
     * Information about some of the results saved in the genotypeResults struct
     *
     */
void compareVCFs(FILE *fh, stList *hmms, char *vcf_toEval, char *vcf_ref,
                 stBaseMapper *baseMapper, stGenotypeResults *results, stRPHmmParameters *params) {

    st_logInfo("VCF reference: %s \n", vcf_ref);
    st_logInfo("VCF being evaluated: %s \n", vcf_toEval);

    vcfFile *inRef = vcf_open(vcf_ref,"r"); //open vcf file
    if (inRef == NULL) {
        st_logCritical("ERROR: cannot open reference vcf, %s\n", vcf_ref);
        return;
    }
    bcf_hdr_t *hdrRef = bcf_hdr_read(inRef); //read header
    bcf1_t *refRecord = bcf_init1(); //initialize for reading

    vcfFile *inEval = vcf_open(vcf_toEval,"r"); //open vcf file
    if (inEval == NULL) {
        st_logCritical("ERROR: cannot open vcf to evaluate, %s\n", vcf_toEval);
        return;
    }
    bcf_hdr_t *hdrEval = bcf_hdr_read(inEval); //read header
    bcf1_t *evalRecord = bcf_init1(); //initialize for reading
    int64_t referencePos = 0;

    // Start by looking at the first hmm
    int64_t hmmIndex = 0;
    stRPHmm *hmm = stList_get(hmms, hmmIndex);

    // Inefficient, but recalculate the info relevant to the hmm to get bipartitions
    stRPHmm_forwardBackward(hmm);
    stList *path = stRPHmm_forwardTraceBack(hmm);
    stGenomeFragment *gF = stGenomeFragment_construct(hmm, path);
    stSet *reads1 = stRPHmm_partitionSequencesByStatePath(hmm, path, 1);
    stSet *reads2 = stRPHmm_partitionSequencesByStatePath(hmm, path, 0);

    // Iterate through the vcf being checked until getting to the start of the specified interval
    // Don't bother analyzing these records
    int64_t refStart = 0;
    int64_t evalPos =  0;
    bcf1_t *unpackedRecord;

    // Variables for keeping track of phasing info
    bool phasingHap1 = false;
    bool phasingHap2 = false;
    float switchErrorDistance = 0;

//    int64_t curiousIndex = 8098619;
//    if (curiousIndex >= gF->refStart && curiousIndex <= gF->refCoords[gF->length - 1]) {
//        st_logInfo("Got to the curious index, %d\n", curiousIndex);
//        st_setLogLevelFromString("debug");
//        printPartitionInfo(curiousIndex, curiousIndex, reads1, reads2, gF, hmm);
//        st_setLogLevelFromString("info");
//    }


    while(bcf_read(inRef, hdrRef, refRecord) == 0) {

        // To take care of the case where a false positive may have been skipped
        // over if the previous eval location was a false negative
        bool maybeFalsePositive = false;
        if (referencePos < evalPos) {
            maybeFalsePositive = true;
        }
        // Unpack reference record
        bcf1_t *unpackedRecordRef = refRecord;
        bcf_unpack(unpackedRecordRef, BCF_UN_ALL);
        referencePos = unpackedRecordRef->pos+1;
        char *refChar = unpackedRecordRef->d.als;
        char *refAltChar = unpackedRecordRef->d.allele[1];

        // Skip to the first known location of variation in file being evaluated
        if (results->positives == 0) {
            refStart = referencePos;
        }

        // Make sure to only look at records in the specified interval
        if (referencePos < hmm->refStart) continue;

        // If the position is beyond the end of this hmm, get the next one
        while (hmm->refEnd < referencePos) {
            hmmIndex++;
            if (hmmIndex < stList_length(hmms)) {
                // Cleanup old stuff
                stList_destruct(path);
                stGenomeFragment_destruct(gF);
                stSet_destruct(reads1);
                stSet_destruct(reads2);

                // Get new stuff
                hmm = stList_get(hmms, hmmIndex);
                path = stRPHmm_forwardTraceBack(hmm);
                gF = stGenomeFragment_construct(hmm, path);
                reads1 = stRPHmm_partitionSequencesByStatePath(hmm, path, 1);
                reads2 = stRPHmm_partitionSequencesByStatePath(hmm, path, 0);
                phasingHap1 = false;
                phasingHap2 = false;

//                if (curiousIndex >= gF->refStart && curiousIndex <= gF->refCoords[gF->length - 1]) {
//                    st_logInfo("Got to the curious index, %d\n", curiousIndex);
//                    st_setLogLevelFromString("debug");
//                    printPartitionInfo(curiousIndex, curiousIndex, reads1, reads2, gF, hmm);
//                    st_setLogLevelFromString("info");
//                }
            } else {
                break;
            }
        }
        // No more fragments to look through
        if (hmmIndex == stList_length(hmms)) break;

        // Get genotype info
        int i, j, ngt, nsmpl = bcf_hdr_nsamples(hdrRef);
        int32_t *gt_arr = NULL, ngt_arr = 0;
        ngt = bcf_get_genotypes(hdrRef, unpackedRecordRef, &gt_arr, &ngt_arr);

        int allele1 = bcf_gt_allele(gt_arr[0]);
        int allele2 = bcf_gt_allele(gt_arr[1]);
        int64_t *gFIndex = stHash_search(gF->refCoordMap, &referencePos);
        int h1AlphChar = stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString1[*gFIndex]);
        int h2AlphChar = stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[*gFIndex]);
        results->positives++;

        if (maybeFalsePositive && evalPos < referencePos) {
            results->falsePositives++;
            char *evalRefChar = unpackedRecord->d.als;
            char *evalAltChar = unpackedRecord->d.allele[1];
            if (strlen(evalRefChar) > 1 || strlen(evalAltChar) > 1) results->falsePositiveGaps++;
            if (params->verboseFalsePositives) {
                st_logDebug("\nFALSE POSITIVE\n");
                printFalsePositive(unpackedRecord, evalPos, hmm);
                printPartitionInfo(referencePos, evalPos, reads1, reads2, gF, hmm);
            }
        }

        // Iterate through vcf until getting to the position of the variant
        // from the reference vcf currently being looked at
        while (evalPos < referencePos) {
            if (bcf_read(inEval, hdrEval, evalRecord) != 0) {
                break;  // can't read record - no more records in file to evaluate
            }
            unpackedRecord = evalRecord;                // unpack record
            bcf_unpack(unpackedRecord, BCF_UN_INFO);
            evalPos = unpackedRecord->pos+1;
            if (evalPos < refStart) continue;           // skip this record

            // Check for false positives - variations found not in reference
            if (evalPos < referencePos) {
                results->falsePositives++;
                char *evalRefChar = unpackedRecord->d.als;
                char *evalAltChar = unpackedRecord->d.allele[1];
                if (strlen(evalRefChar) > 1 || strlen(evalAltChar) > 1) results->falsePositiveGaps++;
                if (params->verboseFalsePositives) {
                    st_logDebug("\nFALSE POSITIVE \n");
                    printFalsePositive(unpackedRecord, evalPos, hmm);
                    printPartitionInfo(referencePos, evalPos, reads1, reads2, gF, hmm);
                }
            } else {
                break;
            }
        }

        // At locus of known variation
        if (evalPos == referencePos) {
            char *evalRefChar = unpackedRecord->d.als;
            char *evalAltChar = unpackedRecord->d.allele[1];

            if (allele1 == allele2) {
                if ((strcmp(refChar, evalRefChar) == 0 && strcmp(evalAltChar, refAltChar) == 0) || (strcmp(refChar, evalAltChar) == 0 && strcmp(evalRefChar, refAltChar) == 0)) {
                    recordHomozygousVariant(results);
                    if (params->verboseFalsePositives) {
                        st_logDebug("\nHOMOZYGOUS VARIANT IN REF\n");
                        printFalsePositive(unpackedRecord, evalPos, hmm);
                    }
                } else {
                    st_logDebug("\nUNKNOWN ERROR \n");
                    printAlleleInfo(unpackedRecordRef, hmm, referencePos, refChar, evalRefChar, evalAltChar);
                }
            } else if (!phasingHap1 && !phasingHap2) {
                // Beginning of genotype fragment, figure out which haplotype matches with ref
                results->uncertainPhasing++;
                if ((strcmp(refChar, evalRefChar) == 0 && strcmp(evalAltChar, refAltChar) == 0) || (strcmp(refChar, evalAltChar) == 0 && strcmp(evalRefChar, refAltChar) == 0)) {
                    results->truePositives++;
                } else results->falsePositives++;

                if (allele1 == 0 && allele2 == 1) {
                    if (strcmp(refChar, evalRefChar) == 0 && strcmp(evalAltChar, refAltChar) == 0)
                        phasingHap1 = true;
                    else phasingHap2 = true;
                } else if (allele1 == 1 && allele2 == 0) {
                    if (strcmp(refChar, evalAltChar) == 0 && strcmp(evalRefChar, refAltChar) == 0)
                        phasingHap1 = true;
                    else phasingHap2 = true;
                } else {
                    // Homozygous at start of fragment
                }
            } else if (phasingHap1) {
                if (allele1 == 0 && allele2 == 1) {
                    if (strcmp(refChar, evalRefChar) == 0 && strcmp(evalAltChar, refAltChar) == 0) {
                        switchErrorDistance++;
                        recordTruePositive(results, params, unpackedRecordRef, baseMapper, hmm, gF, reads1, reads2);
                    } else if (strcmp(refChar, evalAltChar) == 0 && strcmp(evalRefChar, refAltChar) == 0) {
                        results->switchErrors++;
                        results->switchErrorDistance += switchErrorDistance;
                        switchErrorDistance = 0;
                        phasingHap2 = true;
                        phasingHap1 = false;
                        recordTruePositive(results, params, unpackedRecordRef, baseMapper, hmm, gF, reads1, reads2);
                        st_logDebug("\nSwitch error @ position %" PRIi64 "\n", referencePos);
                    } else {
                        st_logDebug("\nINCORRECT POSITIVE\n");
                        results->falsePositives++;
                        results->error_incorrectVariant++;
                        printAlleleInfo(unpackedRecordRef, hmm, referencePos, refChar, evalRefChar, evalAltChar);
                        printPartitionInfo(referencePos, evalPos, reads1, reads2, gF, hmm);
                    }
                } else {
                    if (strcmp(refChar, evalRefChar) == 0 && strcmp(evalAltChar, refAltChar) == 0) {
                        results->switchErrors++;
                        results->switchErrorDistance += switchErrorDistance;
                        switchErrorDistance = 0;
                        phasingHap2 = true;
                        phasingHap1 = false;
                        recordTruePositive(results, params, unpackedRecordRef, baseMapper, hmm, gF, reads1, reads2);
                    } else if (strcmp(refChar, evalAltChar) == 0 && strcmp(evalRefChar, refAltChar) == 0) {
                        switchErrorDistance++;
                        recordTruePositive(results, params, unpackedRecordRef, baseMapper, hmm, gF, reads1, reads2);
                    } else {
                        st_logDebug("\nINCORRECT POSITIVE\n");
                        results->falsePositives++;
                        results->error_incorrectVariant++;
                        printAlleleInfo(unpackedRecordRef, hmm, referencePos, refChar, evalRefChar, evalAltChar);
                        printPartitionInfo(referencePos, evalPos, reads1, reads2, gF, hmm);
                    }
                }
            } else if (phasingHap2) {
                if (allele1 == 0 && allele2 == 1) {
                    if (strcmp(refChar, evalAltChar) == 0 && strcmp(evalRefChar, refAltChar) == 0) {
                        switchErrorDistance++;
                        recordTruePositive(results, params, unpackedRecordRef, baseMapper, hmm, gF, reads1, reads2);
                    } else if (strcmp(refChar, evalRefChar) == 0 && strcmp(evalAltChar, refAltChar) == 0) {
                        results->switchErrors++;
                        results->switchErrorDistance += switchErrorDistance;
                        switchErrorDistance = 0;
                        phasingHap1 = true;
                        phasingHap2 = false;
                        recordTruePositive(results, params, unpackedRecordRef, baseMapper, hmm, gF, reads1, reads2);
                    } else {
                        st_logDebug("\nINCORRECT POSITIVE\n");
                        printAlleleInfo(unpackedRecordRef, hmm, referencePos, refChar, evalRefChar, evalAltChar);
                        printPartitionInfo(referencePos, evalPos, reads1, reads2, gF, hmm);
                        results->falsePositives++;
                        results->error_incorrectVariant++;
                    }
                } else {
                    if (strcmp(refChar, evalAltChar) == 0 && strcmp(evalRefChar, refAltChar) == 0) {
                        results->switchErrors++;
                        results->switchErrorDistance += switchErrorDistance;
                        switchErrorDistance = 0;
                        phasingHap1 = true;
                        phasingHap2 = false;
                        recordTruePositive(results, params, unpackedRecordRef, baseMapper, hmm, gF, reads1, reads2);
                    } else if (strcmp(refChar, evalRefChar) == 0 && strcmp(evalAltChar, refAltChar) == 0) {
                        switchErrorDistance++;
                        recordTruePositive(results, params, unpackedRecordRef, baseMapper, hmm, gF, reads1, reads2);
                    } else {
                        st_logDebug("\nINCORRECT POSITIVE\n");
                        printAlleleInfo(unpackedRecordRef, hmm, referencePos, refChar, evalRefChar, evalAltChar);
                        printPartitionInfo(referencePos, evalPos, reads1, reads2, gF, hmm);
                        results->falsePositives++;
                        results->error_incorrectVariant++;
                    }
                }
            }

        } else if (evalPos > referencePos){
            // Missed the variant
            double *read1BaseCounts = getProfileSequenceBaseCompositionAtPosition(reads1, referencePos);
            double *read2BaseCounts = getProfileSequenceBaseCompositionAtPosition(reads2, referencePos);

            // Check to make sure that the reference actually had a variation there
            if (allele1 == allele2) {
                results->error_homozygousInRef++;
                results->positives--;
                results->negatives++;
            }
                // False negative - no variation was found, but truth vcf has one
            else if (allele1 != allele2){
                results->falseNegatives++;
                int64_t *gFIndex = stHash_search(gF->refCoordMap, &referencePos);

                // Check if record was an insertion or deletion
                if (strlen(refChar) > 1 || strlen(refAltChar) > 1) {
                    results->error_missedIndels++;
                    size_t indelLen = strlen(refChar);
                    if (strlen(refAltChar) > indelLen) indelLen = strlen(refAltChar);


                    st_logDebug("\nMISS: INDEL\n");
                    printAlleleInfo(unpackedRecordRef, hmm, referencePos, refChar,
                                    (char *) &h1AlphChar, (char *) &h2AlphChar);

                    for (int64_t j = 1; j < indelLen; j++) {
                        st_logDebug("\tNext pos: %" PRIi64 "\n", referencePos+j);
                        printColumnAtPosition(hmm, referencePos+j);
                    }

                    st_logDebug("\tposterior prob: %f\n", gF->genotypeProbs[*gFIndex]);
                } else {
                    results->error_badPartition++;
                    st_logDebug("\nMISS: SNV\n");
                    printAlleleInfo(unpackedRecordRef, hmm, referencePos, refChar, (char *) &h1AlphChar, (char *) &h2AlphChar);

                    st_logDebug("\tPartition 1: \n");
                    printBaseComposition2(read1BaseCounts);
                    st_logDebug("\tPartition 2: \n");
                    printBaseComposition2(read2BaseCounts);
                    st_logDebug("\tposterior prob: %f\n", gF->genotypeProbs[*gFIndex]);

                }
                free(read1BaseCounts);
                free(read2BaseCounts);
            }
            free(gt_arr);
        }
    }
    if (results->truePositives == 0) st_logInfo("No matches between vcfs found - did you compare against the correct vcf?\n");

    // Remaining positions after the last variant in the reference are not currently being looked through
    // False positives in this region could therefore be missed
    // (In addition to false positives after the first variant)
    results->negatives += (referencePos - refStart - results->positives);
    results->trueNegatives += (results->negatives - results->falsePositives);
    results->switchErrorDistance = results->switchErrorDistance/results->switchErrors;

    // cleanup
    vcf_close(inRef);
    vcf_close(inEval);
    bcf_hdr_destroy(hdrRef);
    bcf_hdr_destroy(hdrEval);
    bcf_destroy(refRecord);
    bcf_destroy(evalRecord);
    stSet_destruct(reads1);
    stSet_destruct(reads2);
    stGenomeFragment_destruct(gF);
    stList_destruct(path);
}