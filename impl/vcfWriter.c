/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com) & Arthur Rand (arand@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <ctype.h>
#include <htslib/vcf.h>
#include "stRPHmm.h"

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

//            fprintf(stderr, "pos: %d\n", bcf_rec->pos);
//            fprintf(stderr, "h1: %c \t h2: %c \n", h1AlphChar, h2AlphChar);
//            fprintf(stderr, "n h1: %c \t n h2: %c \n", next_h1AlphChar, next_h2AlphChar);

            if (next_h1AlphChar != next_h2AlphChar) {
                if (h1AlphChar != '-' && h2AlphChar != '-') {
                    // Check to see if there was an insertion or deletion in the next spot
                    if (next_h1AlphChar == '-' && next_h2AlphChar != '-') {
                        kputc(h1AlphChar, &str);
                        kputc(',', &str);
                        kputc(h2AlphChar, &str);
                        kputc(next_h2AlphChar, &str);
                        int j = 2;
                        while(i+j < gF->length && stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString1[i+j]) == '-' && stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[i+j]) != '-') {
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
                        while(i+j < gF->length && stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[i+j]) == '-' && stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString1[i+j]) != '-') {
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
    }

    // cleanup
    free(str.s);
    free(gt_info);
    if (referenceSeq) free(referenceSeq);
    bcf_destroy(bcf_rec);
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
void compareVCFs(FILE *fh, stList *hmms,
                 char *vcf_toEval, char *vcf_ref, double threshold,
                 stBaseMapper *baseMapper, stGenotypeResults *results) {

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

    st_logInfo("> Comparing vcf files \n");

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
    int recordCount = 0;

    while(bcf_read(inRef, hdrRef, refRecord) == 0) {
        // To take care of the case where a false positive may have been skipped
        // over if the previous eval location was a false negative
        bool maybeFalsePositive = false;
        if (referencePos < evalPos) {
            maybeFalsePositive = true;
        }
        // Unpack reference record
        bcf1_t *unpackedRecordRef = refRecord;
        bcf_unpack(unpackedRecordRef, BCF_UN_INFO);
        referencePos = unpackedRecordRef->pos+1;

        if (maybeFalsePositive && evalPos < referencePos) {
            results->falsePositives++;
            char *evalRefChar = unpackedRecord->d.als;
            char *evalAltChar = unpackedRecord->d.allele[1];
            st_logDebug("\nFALSE POSITIVE \n\t pos: %" PRIi64 " ref:%s alt: %s \n",
                        evalPos, evalRefChar, evalAltChar);
            printColumnAtPosition(hmm, evalPos);
            if (strlen(evalRefChar) > 1 || strlen(evalAltChar) > 1) {
                results->falsePositiveGaps++;

                size_t indelLen = strlen(evalRefChar);
                if (strlen(evalAltChar) > indelLen) indelLen = strlen(evalAltChar);
                for (int64_t j = 1; j < indelLen; j++) {
                    st_logDebug("\tNext pos: %" PRIi64 "\n", evalPos+j);
                    printColumnAtPosition(hmm, evalPos+j);
                }
            }
        }

        // Skip to the first known location of variation in file being evaluated
        if (results->positives == 0) refStart = referencePos;

        // Make sure to only look at records in the specified interval
        if (referencePos < hmm->refStart) continue;

        // If the position is beyond the end of this hmm, get the next one
        while ((hmm->refStart + hmm->refLength) < referencePos) {
            hmmIndex++;
            if (hmmIndex < stList_length(hmms)) {
                hmm = stList_get(hmms, hmmIndex);
                path = stRPHmm_forwardTraceBack(hmm);
                gF = stGenomeFragment_construct(hmm, path);
                reads1 = stRPHmm_partitionSequencesByStatePath(hmm, path, 1);
                reads2 = stRPHmm_partitionSequencesByStatePath(hmm, path, 0);
            } else {
                break;
            }
        }

        results->positives++;
        char *refChar = unpackedRecordRef->d.als;
        char *refAltChar = unpackedRecordRef->d.allele[1];


        // Iterate through vcf until getting to the position of the variant
        // from the reference vcf currently being looked at
        while (evalPos < referencePos) {

            if (bcf_read(inEval, hdrEval, evalRecord) != 0) {
                st_logDebug("Error: bcf_read\n");
                break;
            }

            unpackedRecord = evalRecord;                // unpack record
            bcf_unpack(unpackedRecord, BCF_UN_INFO);
            evalPos = unpackedRecord->pos+1;
            char *evalRefChar = unpackedRecord->d.als;
            char *evalAltChar = unpackedRecord->d.allele[1];
            if (evalPos < refStart) continue;           // skip this record
            recordCount++;


            // Check for false positives - variations found not in reference
            if (evalPos < referencePos) {
                results->falsePositives++;
                st_logDebug("\nFALSE POSITIVE \npos: %" PRIi64 " ref:%s alt: %s \n",
                            evalPos, evalRefChar, evalAltChar);
                printColumnAtPosition(hmm, evalPos);
                if (strlen(evalRefChar) > 1 || strlen(evalAltChar) > 1) {
                    results->falsePositiveGaps++;

                    size_t indelLen = strlen(evalRefChar);
                    if (strlen(evalAltChar) > indelLen) indelLen = strlen(evalAltChar);
                    for (int64_t j = 1; j < indelLen; j++) {
                        st_logDebug("\tNext pos: %" PRIi64 "\n", evalPos+j);
                        printColumnAtPosition(hmm, evalPos+j);
                    }
                }
            } else {
                break;
            }
        }
        // At locus of known variation
        if (evalPos == referencePos) {
            char *evalRefChar = unpackedRecord->d.als;
            char *evalAltChar = unpackedRecord->d.allele[1];

            if (strcmp(refChar, evalRefChar) == 0 && strcmp(evalAltChar, refAltChar) == 0) {
                results->truePositives++;
            } else if (strcmp(refChar, evalAltChar) == 0 && strcmp(evalRefChar, refAltChar) == 0) {
                results->truePositives++;
            } else {
                results->falsePositives++;
                st_logDebug("\nINCORRECT POSITIVE \n");
                st_logDebug("pos: %" PRIi64 "\nref: %s\talt: %s\n", referencePos, evalRefChar, evalAltChar);
                st_logDebug("reference: ref: %s\talt: %s\n", refChar, refAltChar);
            }

        } else {
            // Missed the variant
            // False negative - no variation was found, but truth vcf has one
            results->falseNegatives++;

            double *read1BaseCounts = getProfileSequenceBaseCompositionAtPosition(reads1, referencePos);
            double *read2BaseCounts = getProfileSequenceBaseCompositionAtPosition(reads2, referencePos);

            char h1AlphChar = stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString1[referencePos-gF->refStart]);
            char h2AlphChar = stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[referencePos-gF->refStart]);

            // Check if record was an insertion or deletion
            if (strlen(refChar) > 1 || strlen(refAltChar) > 1) {
                results->error_missedIndels++;
                size_t indelLen = strlen(refChar);
                if (strlen(refAltChar) > indelLen) indelLen = strlen(refAltChar);

                st_logDebug("\nMISS: INSERTION / DELETION\n");
                st_logDebug("pos: %" PRIi64 "\nref: %s\talt: ", referencePos, refChar);
                for (int i = 1; i < unpackedRecordRef->n_allele; i++) {
                    if (i != 1) st_logDebug(",");
                    st_logDebug("%s", unpackedRecordRef->d.allele[i]);
                    st_logDebug("\n\toutput alleles: %c, %c\n", h1AlphChar, h2AlphChar);
                    printColumnAtPosition(hmm, referencePos);
                }

                for (int64_t j = 1; j < indelLen; j++) {
                    st_logDebug("\tNext pos: %" PRIi64 "\n", referencePos+j);
                    printColumnAtPosition(hmm, referencePos+j);
                }
            } else {
                // Quantify the type of false negative it was
                double totalCount = 0;
                for (int64_t i = 0; i < ALPHABET_SIZE; i++) {
                    totalCount += read1BaseCounts[i];
                    totalCount += read2BaseCounts[i];
                }
                int refBase = stBaseMapper_getValueForChar(baseMapper, *refChar);
                int altBase = stBaseMapper_getValueForChar(baseMapper, *refAltChar);
                float fractionRefBase = (read1BaseCounts[refBase] + read2BaseCounts[refBase]) / totalCount;
                float fractionAltBase = (read1BaseCounts[altBase] + read2BaseCounts[altBase]) / totalCount;

                if (fractionRefBase < threshold || fractionAltBase < threshold) {
                    results->error_trueVariantWrong++;
                    st_logDebug("\nMISS: TRUE VARIANT WRONG\n");
                } else {
                    results->error_badPartition++;
                    st_logDebug("\nMISS: BAD PARTITION\n");
                }

                st_logDebug("pos: %" PRIi64 "\nref: %s\talt: ", referencePos, refChar);
                for (int i = 1; i < unpackedRecordRef->n_allele; i++) {
                    if (i != 1) st_logDebug(",");
                    st_logDebug("%s", unpackedRecordRef->d.allele[i]);
                    st_logDebug("\n\toutput alleles: %c, %c\n", h1AlphChar, h2AlphChar);
                    printColumnAtPosition(hmm, referencePos);
                }
                st_logDebug("\tPartition 1: \n");
                printBaseComposition2(read1BaseCounts);
                st_logDebug("\tPartition 2: \n");
                printBaseComposition2(read2BaseCounts);

                st_logDebug("\tfraction of reference base seen in reads: %f\n", fractionRefBase);
                st_logDebug("\tfraction of alternate base seen in reads: %f\n", fractionAltBase);
            }
            st_logDebug("\tposterior prob: %f\n", gF->genotypeProbs[referencePos-gF->refStart]);

            free(read1BaseCounts);
            free(read2BaseCounts);
        }
    }

    // Remaining positions after the last variant in the reference are not currently being looked through
    // False positives in this region could therefore be missed
    // (In addition to false positives after the first variant)
    results->negatives = referencePos - refStart - results->positives;
    results->trueNegatives = results->negatives - results->falsePositives;


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


/*
 * Prints information contained in genotypeResults struct.
 */
void printGenotypeResults(stGenotypeResults *results) {
    // Sensitivity
    st_logInfo("\nSensitivity: %f \n(= fraction of true positives compared to reference, \t%"
                       PRIi64 " out of %"PRIi64 ")\n",
               (float)results->truePositives/results->positives,
               results->truePositives, results->positives) ;
    st_logInfo("\t \t(Number of false negatives: %" PRIi64 ")\n", results->falseNegatives);
    st_logInfo("Sensitivity ignoring false negatives not supported by reads:\t %f \n",
               (float)results->truePositives/(results->positives-results->error_trueVariantWrong)) ;

    // Specificity
    st_logInfo("\nSpecificity: %f \n(= fraction of true negatives compared to reference, \t%"
                       PRIi64 " out of % "PRIi64 ")\n",
               (float)results->trueNegatives/results->negatives,
               results->trueNegatives, results->negatives);
    st_logInfo("\t \t(Number of false positives: %" PRIi64 ",\twithout gaps: %" PRIi64 ")\n", results->falsePositives, results->falsePositives-results->falsePositiveGaps);

    // More detailed numbers about errors
    st_logInfo("\nFalse negatives:\n");
    st_logInfo("Partition bad: %" PRIi64 " \t\t\t\t(%f)\n",
               results->error_badPartition, (float)results->error_badPartition/results->falseNegatives);

    st_logInfo("Involve indels in the ref vcf: %" PRIi64 " \t\t(%f)\n",
               results->error_missedIndels, (float)results->error_missedIndels/results->falseNegatives);

    st_logInfo("True variant not supported by reads: %" PRIi64 " \t(%f)\n",
               results->error_trueVariantWrong,
               (float)results->error_trueVariantWrong/results->falseNegatives);
}
