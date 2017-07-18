/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com) & Arthur Rand (arand@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <getopt.h>
#include <stdio.h>
#include <ctype.h>
#include <htslib/vcf.h>
#include <htslib/sam.h>
#include "stRPHmm.h"
#include "sam.h"
#include "bgzf.h"

void usage() {
    fprintf(stderr, "vcfCompare [options] -r VCF_REFERENCE -e VCF_EVALUATED\n");
    fprintf(stderr,
            "Compares VCF_EVALUATED to VCF_REFERENCE on the following metrics:\n"
            "\t1. The concordance of all calls in VCF_REFERENCE compared to VCF_EVALUATED\n"
            "\t2. The concordance of all reference positions in VCF_REFERENCE compared to VCF_EVALUATED\n"
            "\t3. The phasing in VCF_REFERENCE compared to VCF_EVALUATED\n"
    );
    fprintf(stderr, "-h --help : Print this help screen\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-r --vcfReference : Specify reference VCF\n");
    fprintf(stderr, "-e --vcfEvaluated : Specify evaluated VCF\n");
}

void create_vcf() {
    // HEADER // (source:bam_plcmd.mpileup)

    // header objects
    bcf_hdr_t *bcf_hdr = bcf_hdr_init("w");
    kstring_t str = {0,0,NULL};

    // generic info
    ksprintf(&str, "##vcfCompare=htslib-%s\n",hts_version());
    bcf_hdr_append(bcf_hdr, str.s);

    // reference file used
    str.l = 0;
    ksprintf(&str, "##reference=file://%s\n", "TODO");
    bcf_hdr_append(bcf_hdr, str.s);

    // contig
    str.l = 0;
    //ksprintf(&str, "##contig=<ID=%s,length=%d>\n", "chr1", 249250621);
    ksprintf(&str, "##contig=<ID=%s>\n", "chr1");
    bcf_hdr_append(bcf_hdr, str.s);

    // formatting
    bcf_hdr_append(bcf_hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");

    // samples
    bcf_hdr_add_sample(bcf_hdr, "SMPL1");
    bcf_hdr_add_sample(bcf_hdr, NULL);

    // write header
    vcfFile *out = vcf_open("out.vcf", "w");
    bcf_hdr_write(out, bcf_hdr);
    vcf_close(out);


    // RECORD // (source:test-vcf-api.c)
    out = vcf_open("out.vcf", "a");
    bcf1_t *bcf_rec = bcf_init1();
    int32_t filter_info = bcf_hdr_id2int(bcf_hdr, BCF_DT_ID, "PASS"); //how to specify not pass? do we care?
    int32_t *gt_info = (int*)malloc(bcf_hdr_nsamples(bcf_hdr)*2*sizeof(int)); //array specifying phasing

    //prep
    bcf_clear1(bcf_rec);
    str.l = 0;

    // contig (CHROM)
    bcf_rec->rid = bcf_hdr_name2id(bcf_hdr, "chr1"); //defined in a contig in the top
    // position (POS) - this is 0-based indexing, printed as 1-based
    bcf_rec->pos  = 1;
    // identifier (ID)
    bcf_update_id(bcf_hdr, bcf_rec, "ex12345");
    // quality (QUAL)
    bcf_rec->qual = 0;
    // reference (REF) - comma separated list, first is REF column
    kputc('A', &str);
    // alternate (ALT) - ..rest are in ALT
    kputc(',C,T', &str);
    bcf_update_alleles_str(bcf_hdr, bcf_rec, str.s);
    // filtering (FILTER)
    bcf_update_filter(bcf_hdr, bcf_rec, &filter_info, 1);
    // genotype (FORMAT / $SMPL1 [$SMPL2..])
    gt_info[0] = bcf_gt_phased(0);
    gt_info[1] = bcf_gt_phased(1);
    bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);
    // save it
    bcf_write1(out, bcf_hdr, bcf_rec);



    // prep
    bcf_clear1(bcf_rec);
    // CHROM
    bcf_rec->rid = bcf_hdr_name2id(bcf_hdr, "chr1");
    // POS
    bcf_rec->pos  = 2;
    // REF, ALT
    bcf_update_alleles_str(bcf_hdr, bcf_rec, "A");
    // FILTER
    filter_info = bcf_hdr_id2int(bcf_hdr, BCF_DT_ID, "PASS");
    bcf_update_filter(bcf_hdr, bcf_rec, &filter_info, 1);
    // FORMAT - GT
    gt_info = (int*)malloc(bcf_hdr_nsamples(bcf_hdr)*2*sizeof(int));
    gt_info[0] = bcf_gt_phased(0);
    gt_info[1] = bcf_gt_unphased(1);
    bcf_update_genotypes(bcf_hdr, bcf_rec, gt_info, bcf_hdr_nsamples(bcf_hdr)*2);
    // save it
    bcf_write1(out, bcf_hdr, bcf_rec);


    free(str.s);
    vcf_close(out);
}

void fill_hdr(bcf_hdr_t **hdr) {
    *hdr = bcf_hdr_init("w");
    if (hdr == NULL) {
        printf("Null");
    } else {
        printf("not null");
    }
}

void validate_hdr(bcf_hdr_t *hdr) {
    if (hdr == NULL) {
        printf("Null");
    } else {
        printf("not null");
    }
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

void printAlleleInfo(bcf1_t *unpackedRecordRef, stRPHmm *hmm, int64_t referencePos, char *refChar, int h1AlphChar, int h2AlphChar) {
    st_logDebug("  pos: %" PRIi64 "\n  ref: %s   alt: ", referencePos, refChar);
    for (int i = 1; i < unpackedRecordRef->n_allele; i++) {
        if (i != 1) st_logDebug(",");
        st_logDebug("%s", unpackedRecordRef->d.allele[i]);
        st_logDebug("\n  output: %c, %c\n", h1AlphChar, h2AlphChar);
        printColumnAtPosition(hmm, referencePos);
    }
}

void printPartitionInfo(int64_t referencePos, int64_t evalPos, stSet *reads1, stSet *reads2, stGenomeFragment *gF) {
    // print additional partition info
    double *read1BaseCounts = getProfileSequenceBaseCompositionAtPosition(reads1, referencePos);
    double *read2BaseCounts = getProfileSequenceBaseCompositionAtPosition(reads2, referencePos);
    st_logDebug("\tPartition 1: \n");
    printBaseComposition2(read1BaseCounts);
    st_logDebug("\tPartition 2: \n");
    printBaseComposition2(read2BaseCounts);
    st_logDebug("\tposterior prob: %f\n", gF->genotypeProbs[evalPos-gF->refStart]);
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
    int h1AlphChar = stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString1[referencePos-gF->refStart]);
    int h2AlphChar = stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[referencePos-gF->refStart]);

    results->truePositives++;
    if (strlen(refChar) > 1 || strlen(refAltChar) > 1) results->truePositiveGaps++;
    if (params->verboseTruePositives) {
        st_logDebug("\nTRUE POSITIVE\n");
        printAlleleInfo(unpackedRecordRef, hmm, referencePos, refChar, h1AlphChar, h2AlphChar);
        printPartitionInfo(referencePos, referencePos, reads1, reads2, gF);
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
     * This is a more limited version of the compareVCFs found in outputWriter.c,
     * as it does not have the knowledge of the hmms used in marginPhase
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
    stGenomeFragment_setInsertionCounts(gF);
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

    int64_t curiousIndex = 8098619;
    if (curiousIndex >= gF->refStart && curiousIndex < gF->refCoords[gF->length - 1]) {
        st_logInfo("Got to the curious index, %d\n", curiousIndex);
        st_setLogLevelFromString("debug");
        printPartitionInfo(curiousIndex, curiousIndex, reads1, reads2, gF);
        st_setLogLevelFromString("info");
    }


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
        while ((hmm->refStart + hmm->refLength) < referencePos) {
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
                stGenomeFragment_setInsertionCounts(gF);
                reads1 = stRPHmm_partitionSequencesByStatePath(hmm, path, 1);
                reads2 = stRPHmm_partitionSequencesByStatePath(hmm, path, 0);
                phasingHap1 = false;
                phasingHap2 = false;

                if (curiousIndex >= gF->refStart && curiousIndex < gF->refCoords[gF->length - 1]) {
                    st_logInfo("Got to the curious index, %d\n", curiousIndex);
                    st_setLogLevelFromString("debug");
                    printPartitionInfo(curiousIndex, curiousIndex, reads1, reads2, gF);
                    st_setLogLevelFromString("info");
                }
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
        int h1AlphChar = stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString1[referencePos-gF->refStart]);
        int h2AlphChar = stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[referencePos-gF->refStart]);
        results->positives++;

        if (maybeFalsePositive && evalPos < referencePos) {
            results->falsePositives++;
            char *evalRefChar = unpackedRecord->d.als;
            char *evalAltChar = unpackedRecord->d.allele[1];
            if (strlen(evalRefChar) > 1 || strlen(evalAltChar) > 1) results->falsePositiveGaps++;
            if (params->verboseFalsePositives) {
                st_logDebug("\nFALSE POSITIVE\n");
                printFalsePositive(unpackedRecord, evalPos, hmm);
                printPartitionInfo(referencePos, evalPos, reads1, reads2, gF);
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

//            st_logInfo("Reference pos: %d \t eval pos: %d \n", referencePos, evalPos);

            // Check for false positives - variations found not in reference
            if (evalPos < referencePos) {
                results->falsePositives++;
                char *evalRefChar = unpackedRecord->d.als;
                char *evalAltChar = unpackedRecord->d.allele[1];
                if (strlen(evalRefChar) > 1 || strlen(evalAltChar) > 1) results->falsePositiveGaps++;
                if (params->verboseFalsePositives) {
                    st_logDebug("\nFALSE POSITIVE \n");
                    printFalsePositive(unpackedRecord, evalPos, hmm);
                    printPartitionInfo(referencePos, evalPos, reads1, reads2, gF);
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
                    printAlleleInfo(unpackedRecordRef, hmm, referencePos, refChar, h1AlphChar, h2AlphChar);
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
                        printAlleleInfo(unpackedRecordRef, hmm, referencePos, refChar, h1AlphChar, h2AlphChar);
                        printPartitionInfo(referencePos, evalPos, reads1, reads2, gF);
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
                        printAlleleInfo(unpackedRecordRef, hmm, referencePos, refChar, h1AlphChar, h2AlphChar);
                        printPartitionInfo(referencePos, evalPos, reads1, reads2, gF);
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
                        printAlleleInfo(unpackedRecordRef, hmm, referencePos, refChar, h1AlphChar, h2AlphChar);
                        printPartitionInfo(referencePos, evalPos, reads1, reads2, gF);
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
                        printAlleleInfo(unpackedRecordRef, hmm, referencePos, refChar, h1AlphChar, h2AlphChar);
                        printPartitionInfo(referencePos, evalPos, reads1, reads2, gF);
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
                // Check if record was an insertion or deletion
                if (strlen(refChar) > 1 || strlen(refAltChar) > 1) {
                    results->error_missedIndels++;
                    size_t indelLen = strlen(refChar);
                    if (strlen(refAltChar) > indelLen) indelLen = strlen(refAltChar);

                    st_logDebug("\nMISS: INDEL\n");
                    printAlleleInfo(unpackedRecordRef, hmm, referencePos, refChar, h1AlphChar, h2AlphChar);

                    for (int64_t j = 1; j < indelLen; j++) {
                        st_logDebug("\tNext pos: %" PRIi64 "\n", referencePos+j);
                        printColumnAtPosition(hmm, referencePos+j);
                    }
                    st_logDebug("\tposterior prob: %f\n", gF->genotypeProbs[referencePos-gF->refStart]);
                } else {
                    results->error_badPartition++;
                    st_logDebug("\nMISS: SNV\n");
                    printAlleleInfo(unpackedRecordRef, hmm, referencePos, refChar, h1AlphChar, h2AlphChar);

                    st_logDebug("\tPartition 1: \n");
                    printBaseComposition2(read1BaseCounts);
                    st_logDebug("\tPartition 2: \n");
                    printBaseComposition2(read2BaseCounts);
                    st_logDebug("\tposterior prob: %f\n", gF->genotypeProbs[referencePos-gF->refStart]);

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


int main(int argc, char *argv[]) {
    // Parameters / arguments

    char * logLevelString = NULL;
    char *vcfReference = NULL;
    char *vcfEvaluated = NULL;

    // Parse the options
    while (1) {
        static struct option long_options[] = {
                { "logLevel", required_argument, 0, 'a' },
                { "help", no_argument, 0, 'h' },
                { "vcfReference", required_argument, 0, 'r'},
                { "vcfEvaluated", required_argument, 0, 'e'},
                { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:r:e:h", long_options, &option_index);

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
        case 'r':
            vcfReference = stString_copy(optarg);
            break;
        case 'e':
            vcfEvaluated = stString_copy(optarg);
            break;
        default:
            usage();
            return 1;
        }
    }

    bcf_hdr_t *hdr = NULL;
    fill_hdr(&hdr);
    if (hdr == NULL) {
        printf("Null");
    } else {
        printf("not null");
    }
    validate_hdr(hdr);

    if (true) return 0;
////    create_vcf();
////    fai_build()
//
//    st_logInfo("VCF Reference: %s \n", vcfReference);
//    st_logInfo("VCF Evaluated: %s \n", vcfEvaluated);
//
//    vcfFile *inRef = vcf_open(vcfReference,"r"); //open vcf file
//    bcf_hdr_t *hdrRef = bcf_hdr_read(inRef); //read header
//    bcf_t *record = bcf_init1(); //initialize for reading
//
//
//    vcfFile *out = vcf_open("out.vcf", "w");
//    vcfFile *outOrig = vcf_open("outOrig.vcf", "w");
////    bcf_hdr_write(out, hdrRef);
//
//    int lineNr = 0;
//    while(bcf_read(inRef,hdrRef,record) == 0){
//
//        st_logDebug("%d:\n", lineNr);
//        st_logDebug("\ttid:%d pos:%d bin:%d qual:%d l_qname:%d flag:%d unused1:%d l_extranul:%d n_cigar:%d "
//                            "l_qseq:%d mtid:%d mpos:%d isize:%d \n", record->core.tid, record->core.pos,
//                    record->core.bin, record->core.qual, record->core.l_qname, record->core.flag, record->core.unused1,
//                    record->core.l_extranul, record->core.n_cigar, record->core.l_qseq, record->core.mtid,
//                    record->core.mpos, record->core.isize);
//
//        bcf1_t *unpackedRecord = record;
//        bcf_unpack(unpackedRecord, BCF_UN_ALL);
//
//        st_logDebug("\tid:%s ref:%s alleles:", unpackedRecord->d.id, unpackedRecord->d.als);
//        for (int i = 1; i < unpackedRecord->n_allele; i++) {
//            if (i!=1) st_logDebug(",");
//            st_logDebug(unpackedRecord->d.allele[i]);
//        }
//
//        int nsmpl = bcf_hdr_nsamples(hdrRef);
//        int32_t *gt_arr = NULL, ngt_arr = 0; //locations for get_genotype types and count
//
//        int ngt = bcf_get_genotypes(hdrRef, record, &gt_arr, &ngt_arr);
//        if ( ngt > 0 ) { // GT is present
//
//            st_logDebug(" phasing:");
//            int max_ploidy = ngt / nsmpl;
//            for (int i = 0; i < nsmpl; i++) {
//                int32_t *ptr = gt_arr + i * max_ploidy;
//                for (int j = 0; j < max_ploidy; j++) {
//                    if (ptr[j] == bcf_int32_vector_end) break;// if true, the sample has smaller ploidy
//                    if (bcf_gt_is_missing(ptr[j])) continue;// missing allele
//                    int allele_index = bcf_gt_allele(ptr[j]); // the VCF 0-based allele index
//                    int is_phased = bcf_gt_is_phased(ptr[j]);
//
//                    if (j!=0) {
//                        if (is_phased) st_logDebug("|");
//                        else st_logDebug("/");
//                    }
//                    st_logDebug("%d", allele_index);
//                }
//            }
//            free(gt_arr);
//            st_logDebug("\n");
//        }
//
////        //write original record again
////        bcf_write(outOrig, hdrRef, record);
////
////        //modify, write, reset
////        record->id += 1;
////        record->l_data += 1;
////        record->m_data += 1;
////        bcf_write(out, hdrRef, record);
////        record = bcf_init1();
//
//        lineNr++;
//        if (lineNr > 8) {
//            break;
//        }
//    }
//
//    vcf_close(inRef);
//    vcf_close(out);
//    vcf_close(outOrig);
//
//
//
//    //while(1); // Use this for testing for memory leaks
//
//    return 0;
}

