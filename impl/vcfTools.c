/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com) & Arthur Rand (arand@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <getopt.h>
#include <stdio.h>
#include <htslib/vcf.h>
#include "stRPHmm.h"

typedef struct _vcfRecordComparisonInfo vcfRecordComparisonInfo;
struct _vcfRecordComparisonInfo {
    // The vcf records
    bcf1_t *unpackedRecordRef;
    bcf1_t *unpackedRecordEval;

    // The positions
    int64_t referencePos;
    int64_t evalPos;

    // The genotypes in the vcf
    char *refAllele;
    char *altAllele1;
    char *evalRefAllele;
    char *evalAltAllele;
    char *gt_ref_hap1;
    char *gt_ref_hap2;
    char *gt_eval_hap1;
    char *gt_eval_hap2;
    char h1AlphChar;
    char h2AlphChar;

    // The phasing info
    int refPhasing1;
    int refPhasing2;
    int evalPhasing1;
    int evalPhasing2;

    // Phasing stuff
    bool phasingHap1;
    bool phasingHap2;
};

void printBaseComposition2(double *baseCounts) {
    /*
     * Print the counts/fraction of each alphabet character in a slightly more compressed form.
     */
    double totalCount = 0;
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        totalCount += baseCounts[i];
    }
    st_logDebug("\t\t0 (A)\t1 (C)\t2 (G)\t3 (T)\t4 (-) \n");
    st_logDebug("    Counts:");
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        st_logDebug("\t%0.1f", baseCounts[i]);
    }
    st_logDebug("\n    Fraction:");
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        st_logDebug("\t%0.3f", baseCounts[i]/totalCount);
    }
    st_logDebug( "\n");
}


void printColumnAtPosition(stRPHmm *hmm, int64_t pos) {
    /*
     * Print out the columns of the hmm at a specific positionk
     */
    stRPColumn *column = hmm->firstColumn;
    while(1) {
        if (pos >= column->refStart && pos < column->refStart+column->length) {
            double *columnBaseCounts = getColumnBaseComposition(column, pos);
            st_logDebug("\tTotal counts: \n");
            printBaseComposition2(columnBaseCounts);
            free(columnBaseCounts);
        }
        if (column->nColumn == NULL) {
            break;
        }
        column = column->nColumn->nColumn;
    }
}

void printAlleleInfo(vcfRecordComparisonInfo *vcfInfo) {
    /*
     * Print information about the alleles found in the reference vcf and in the vcf being evaluated.
     */
    st_logDebug("  pos: %" PRIi64 "\n", vcfInfo->referencePos);
    st_logDebug("  baseline  ref: %s\t alt: ", vcfInfo->refAllele);
    for (int i = 1; i < vcfInfo->unpackedRecordRef->n_allele; i++) {
        if (i != 1) st_logDebug(",");
        st_logDebug("%s", vcfInfo->unpackedRecordRef->d.allele[i]);
    }
    st_logDebug(" \tphasing: %d | %d", vcfInfo->refPhasing1, vcfInfo->refPhasing2);
    st_logDebug("\n  output    ref: %s\t alt: %s   \tphasing: %d | %d\n",
                vcfInfo->evalRefAllele, vcfInfo->evalAltAllele, vcfInfo->evalPhasing1, vcfInfo->evalPhasing2);
}

void printPartitionInfo(int64_t pos, stSet *reads1, stSet *reads2, stGenomeFragment *gF) {
    /*
     * Print information about the base composition of the partitions.
     */
    double *read1BaseCounts = getProfileSequenceBaseCompositionAtPosition(reads1, pos);
    double *read2BaseCounts = getProfileSequenceBaseCompositionAtPosition(reads2, pos);
    st_logDebug("\tPartition 1: \n");
    printBaseComposition2(read1BaseCounts);
    st_logDebug("\tPartition 2: \n");
    printBaseComposition2(read2BaseCounts);
    st_logDebug("\t\t\tposterior prob: %f\n", gF->genotypeProbs[pos-gF->refStart]);
}

void printBaseCompositionAtPosition(int64_t pos, stReferencePriorProbs *rProbs) {
    /*
     * Print the counts/fraction of each alphabet character in a slightly more compressed form.
     */
    double totalCount = 0;
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        totalCount += rProbs->baseCounts[(pos - rProbs->refStart)*ALPHABET_SIZE + i];
    }
    st_logDebug("\t\t0 (A)\t1 (C)\t2 (G)\t3 (T)\t4 (-) \n");
    st_logDebug("    Counts:");
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        st_logDebug("\t%0.1f", rProbs->baseCounts[(pos - rProbs->refStart)*ALPHABET_SIZE + i]);
    }
    st_logDebug("\n    Fraction:");
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        st_logDebug("\t%0.3f", rProbs->baseCounts[(pos - rProbs->refStart)*ALPHABET_SIZE + i]/totalCount);
    }
    st_logDebug( "\n");
}

void printTotalBaseCompositionAtPosition(int64_t pos, stReferencePriorProbs *rProbs1, stReferencePriorProbs *rProbs2) {
    /*
     * Print the counts/fraction of each alphabet character in a slightly more compressed form.
     */
    double totalCount = 0;
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        totalCount += rProbs1->baseCounts[(pos - rProbs1->refStart)*ALPHABET_SIZE + i];
        totalCount += rProbs2->baseCounts[(pos - rProbs2->refStart)*ALPHABET_SIZE + i];
    }
    st_logDebug("\t\t0 (A)\t1 (C)\t2 (G)\t3 (T)\t4 (-) \n");
    st_logDebug("    Counts:");
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        st_logDebug("\t%0.1f", rProbs1->baseCounts[(pos - rProbs1->refStart)*ALPHABET_SIZE + i] +
                rProbs2->baseCounts[(pos - rProbs2->refStart)*ALPHABET_SIZE + i]);
    }
    st_logDebug("\n    Fraction:");
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        st_logDebug("\t%0.3f", (rProbs1->baseCounts[(pos - rProbs1->refStart)*ALPHABET_SIZE + i] +
                rProbs2->baseCounts[(pos - rProbs2->refStart)*ALPHABET_SIZE + i]) /totalCount);
    }
    st_logDebug( "\n");
}

void printPartitionInfo2(int64_t pos, stReferencePriorProbs *rProbs1, stReferencePriorProbs *rProbs2) {
    /*
     * Print information about the base composition of the partitions.
     */
    st_logDebug("\tTotal counts: \n");
    printTotalBaseCompositionAtPosition(pos, rProbs1, rProbs2);
    st_logDebug("\tPartition 1: \n");
    printBaseCompositionAtPosition(pos, rProbs1);
    st_logDebug("\tPartition 2: \n");
    printBaseCompositionAtPosition(pos, rProbs2);
}

static void recordHomozygousVariantBasic(stGenotypeResults *results, vcfRecordComparisonInfo *vcfInfo) {
    /*
     * Record a variant that was homozygous in the reference
     */
    // VCF eval seems to count these both as false positives and as false negatives
    results->falsePositives++;
    results->falseNegatives++;
    results->error_homozygousInRef++;
    if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1)) {
        results->error_homozygous_Deletions++;
    }
    else if (strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) {
        results->error_homozygous_Insertions++;
    }
}

static void recordHomozygousVariant(stGenotypeResults *results, vcfRecordComparisonInfo *vcfInfo,
stBaseMapper *baseMapper, stRPHmm *hmm, stGenomeFragment *gF, stSet *reads1, stSet *reads2) {
    /*
     * Record a variant that was homozygous in the reference (phasing 1/1)
     */
    recordHomozygousVariantBasic(results, vcfInfo);


    if (hmm->parameters->verboseFalseNegatives || hmm->parameters->verboseFalsePositives) {

        size_t indelLen = strlen(vcfInfo->refAllele) > strlen(vcfInfo->gt_ref_hap1) ? strlen(vcfInfo->refAllele) : strlen(vcfInfo->gt_ref_hap1);

        char *refSeq = st_calloc(indelLen, sizeof(char));
        for (int64_t i = 0; i < indelLen; i++) {
            refSeq[i] = stBaseMapper_getCharForValue(baseMapper, gF->referenceSequence[vcfInfo->referencePos - gF->refStart + i]);
        }
        st_logDebug("  pos: %" PRIi64 "\n", vcfInfo->referencePos);
        st_logDebug("  baseline  ref: %s\t alt: %s \tphasing: %d | %d\n", refSeq,
                    vcfInfo->gt_ref_hap1, vcfInfo->refPhasing1, vcfInfo->refPhasing2);
        st_logDebug("  output    ref: %s\t alt: %s \tphasing: %d | %d\n",
                    vcfInfo->evalRefAllele, vcfInfo->evalAltAllele, vcfInfo->evalPhasing1, vcfInfo->evalPhasing2);
        printColumnAtPosition(hmm, vcfInfo->referencePos);
        printPartitionInfo(vcfInfo->referencePos, reads1, reads2, gF);

        for (int64_t j = 1; j < indelLen; j++) {
            st_logDebug(" -> pos: %" PRIi64 "\n", vcfInfo->referencePos+j);
            printColumnAtPosition(hmm, vcfInfo->referencePos+j);
            printPartitionInfo(vcfInfo->referencePos+j, reads1, reads2, gF);
        }
        free(refSeq);
    }
}

static void recordHomozygousVariant2(stGenotypeResults *results, vcfRecordComparisonInfo *vcfInfo,
                                    stBaseMapper *baseMapper, stRPHmmParameters *params,
                                     stReferencePriorProbs *rProbs1, stReferencePriorProbs *rProbs2) {
    /*
     * Record a variant that was homozygous in the reference (phasing 1/1)
     */
    recordHomozygousVariantBasic(results, vcfInfo);

    if (params->verboseFalseNegatives || params->verboseFalsePositives) {

        size_t indelLen = strlen(vcfInfo->refAllele) > strlen(vcfInfo->gt_ref_hap1) ? strlen(vcfInfo->refAllele) : strlen(vcfInfo->gt_ref_hap1);

        char *refSeq = st_calloc(indelLen, sizeof(char));
        for (int64_t i = 0; i < indelLen; i++) {
            refSeq[i] = stBaseMapper_getCharForValue(baseMapper, rProbs1->referenceSequence[vcfInfo->referencePos - rProbs1->refStart + i]);
        }
        st_logDebug("  pos: %" PRIi64 "\n", vcfInfo->referencePos);
        st_logDebug("  baseline  ref: %s\t alt: %s \tphasing: %d | %d\n", refSeq,
                    vcfInfo->gt_ref_hap1, vcfInfo->refPhasing1, vcfInfo->refPhasing2);
        st_logDebug("  output    ref: %s\t alt: %s \tphasing: %d | %d\n",
                    vcfInfo->evalRefAllele, vcfInfo->evalAltAllele, vcfInfo->evalPhasing1, vcfInfo->evalPhasing2);
        printPartitionInfo2(vcfInfo->referencePos, rProbs1, rProbs2);

        for (int64_t j = 1; j < indelLen; j++) {
            st_logDebug(" -> pos: %" PRIi64 "\n", vcfInfo->referencePos+j);
            printPartitionInfo2(vcfInfo->referencePos+j, rProbs1, rProbs2);
        }
        free(refSeq);
    }
}

static void recordHomozygousVariantMissBasic(stGenotypeResults *results, vcfRecordComparisonInfo *vcfInfo) {
    /*
     * Record a variant that was homozygous in the reference
     */
    // VCF eval counts these as only false negative
    results->falseNegatives++;
    results->error_homozygousInRef++;
    if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1)) {
        results->error_homozygous_Deletions++;
    }
    else if (strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) {
        results->error_homozygous_Insertions++;
    }
}

static void recordHomozygousVariantMiss(stGenotypeResults *results, vcfRecordComparisonInfo *vcfInfo,
                                    stBaseMapper *baseMapper, stRPHmm *hmm, stGenomeFragment *gF, stSet *reads1, stSet *reads2) {
    /*
     * Record a variant that was homozygous in the reference
     */
    recordHomozygousVariantMissBasic(results, vcfInfo);

    if (hmm->parameters->verboseFalseNegatives) {
        if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1)) st_logDebug("\nMISS - HOMOZYGOUS ALT VARIANT - DELETION  (FN)\n");
        else if (strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) st_logDebug("\nMISS - HOMOZYGOUS ALT VARIANT - INSERTION  (FN)\n");
        else st_logDebug("\nMISS - HOMOZYGOUS ALT VARIANT - SNV  (FN)\n");

        size_t indelLen = strlen(vcfInfo->refAllele) > strlen(vcfInfo->gt_ref_hap1) ? strlen(vcfInfo->refAllele) : strlen(vcfInfo->gt_ref_hap1);

        char *refSeq = st_calloc(indelLen, sizeof(char));
        for (int64_t i = 0; i < indelLen; i++) {
            refSeq[i] = stBaseMapper_getCharForValue(baseMapper, gF->referenceSequence[vcfInfo->referencePos - gF->refStart + i]);
        }
        char refChar = stBaseMapper_getCharForValue(baseMapper, gF->referenceSequence[vcfInfo->referencePos - gF->refStart]);

        st_logDebug("  pos: %" PRIi64 "\n", vcfInfo->referencePos);
        st_logDebug("  baseline  ref: %s\t alt: %s \tphasing: %d | %d\n", vcfInfo->refAllele,
                    vcfInfo->gt_ref_hap1, vcfInfo->refPhasing1, vcfInfo->refPhasing2);
        if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1)) {
            st_logDebug("  output    ref: %s\t alt: %s \tphasing: %d | %d\n", refSeq, refSeq, 0, 0);
        } else {
            st_logDebug("  output    ref: %c\t alt: %c \tphasing: %d | %d\n", refChar, refChar, 0, 0);
        }
        printColumnAtPosition(hmm, vcfInfo->referencePos);
        printPartitionInfo(vcfInfo->referencePos, reads1, reads2, gF);

        for (int64_t j = 1; j < indelLen; j++) {
            st_logDebug(" -> pos: %" PRIi64 "\n", vcfInfo->referencePos+j);
            printColumnAtPosition(hmm, vcfInfo->referencePos+j);
            printPartitionInfo(vcfInfo->referencePos+j, reads1, reads2, gF);
        }
        free(refSeq);
    }
}

static void recordHomozygousVariantMiss2(stGenotypeResults *results, vcfRecordComparisonInfo *vcfInfo,
                                        stBaseMapper *baseMapper, stRPHmmParameters *params,
                                         stReferencePriorProbs *rProbs1, stReferencePriorProbs *rProbs2) {
    /*
     * Record a variant that was homozygous in the reference
     */
    recordHomozygousVariantMissBasic(results, vcfInfo);

    if (params->verboseFalseNegatives) {
        if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1)) st_logDebug("\nMISS - HOMOZYGOUS ALT VARIANT - DELETION  (FN)\n");
        else if (strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) st_logDebug("\nMISS - HOMOZYGOUS ALT VARIANT - INSERTION  (FN)\n");
        else st_logDebug("\nMISS - HOMOZYGOUS ALT VARIANT - SNV  (FN)\n");


        size_t indelLen = strlen(vcfInfo->refAllele) > strlen(vcfInfo->gt_ref_hap1) ? strlen(vcfInfo->refAllele) : strlen(vcfInfo->gt_ref_hap1);

        char *refSeq = st_calloc(indelLen, sizeof(char));
        for (int64_t i = 0; i < indelLen; i++) {
            refSeq[i] = stBaseMapper_getCharForValue(baseMapper, rProbs1->referenceSequence[vcfInfo->referencePos - rProbs1->refStart + i]);
        }
        char refChar = stBaseMapper_getCharForValue(baseMapper, rProbs1->referenceSequence[vcfInfo->referencePos - rProbs1->refStart]);

        st_logDebug("  pos: %" PRIi64 "\n", vcfInfo->referencePos);
        st_logDebug("  baseline  ref: %s\t alt: %s \tphasing: %d | %d\n", vcfInfo->refAllele,
                    vcfInfo->gt_ref_hap1, vcfInfo->refPhasing1, vcfInfo->refPhasing2);
        if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1)) {
            st_logDebug("  output    ref: %s\t alt: %s \tphasing: %d | %d\n", refSeq, refSeq, 0, 0);
        } else {
            st_logDebug("  output    ref: %c\t alt: %c \tphasing: %d | %d\n", refChar, refChar, 0, 0);
        }

        printPartitionInfo2(vcfInfo->referencePos, rProbs1, rProbs2);

        for (int64_t j = 1; j < indelLen; j++) {
            st_logDebug(" -> pos: %" PRIi64 "\n", vcfInfo->referencePos+j);
            printPartitionInfo2(vcfInfo->referencePos+j, rProbs1, rProbs2);
        }
        free(refSeq);
    }
}






static void recordTruePositiveBasic(stGenotypeResults *results, vcfRecordComparisonInfo *vcfInfo) {
    /*
     * Record stats and print info about a true positive result.
     */
    results->truePositives++;
    if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1)  ||
            strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) {
        results->truePositiveIndels++;
    }
}

static void recordTruePositive(stGenotypeResults *results, stRPHmmParameters *params, vcfRecordComparisonInfo *vcfInfo, stRPHmm *hmm, stGenomeFragment *gF, stSet *reads1, stSet *reads2) {
    /*
     * Record stats and print info about a true positive result.
     */
    recordTruePositiveBasic(results, vcfInfo);

    if (params->verboseTruePositives) {
        if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1)  ||
            strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) {
            st_logDebug("\nTRUE POSITIVE - INDEL\n");
        } else {
            st_logDebug("\nTRUE POSITIVE\n");
        }
        printAlleleInfo(vcfInfo);
        printColumnAtPosition(hmm, vcfInfo->referencePos);
        printPartitionInfo(vcfInfo->referencePos, reads1, reads2, gF);
    }
}

static void recordTruePositive2(stGenotypeResults *results, stRPHmmParameters *params, vcfRecordComparisonInfo *vcfInfo, stReferencePriorProbs *rProbs1, stReferencePriorProbs *rProbs2) {
    /*
     * Record stats and print info about a true positive result.
     */
    recordTruePositiveBasic(results, vcfInfo);

    if (params->verboseTruePositives) {
        if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1)  ||
            strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) {
            st_logDebug("\nTRUE POSITIVE - INDEL\n");
        } else {
            st_logDebug("\nTRUE POSITIVE\n");
        }
        printAlleleInfo(vcfInfo);
        printPartitionInfo2(vcfInfo->referencePos, rProbs1, rProbs2);
    }
}



static void recordIncorrectVariantBasic(vcfRecordComparisonInfo *vcfInfo, stGenotypeResults *results) {
    /*
     * Records an error where the bases predicted are not consistent with the ref
     */
    // in vcfEval, these appear as both FP and FN
    results->falseNegatives++;
    results->falsePositives++;

    if (vcfInfo->refPhasing1 == vcfInfo->refPhasing2) {
        results->error_homozygousInRef++;
        if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1)) results->error_homozygous_Deletions++;
        else if (strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) results->error_homozygous_Insertions++;
    } else {
        results->error_missedHet++;
        if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1)) {
            results->error_missedHet_Deletions++;
            results->falsePositiveIndels++;
        } else if (strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) {
            results->error_missedHet_Insertions++;
            results->falsePositiveIndels++;
        }
    }
}


static void printIncorrectVariantDescription(vcfRecordComparisonInfo *vcfInfo, const stRPHmmParameters *params) {
    /*
     * Prints the description of the type of incorrect variant found.
     */
    if (vcfInfo->refPhasing1 == vcfInfo->refPhasing2) {
        if (params->verboseFalseNegatives || params->verboseFalsePositives) {
            st_logDebug("\nERROR: PREDICTED HET FOR HOMOZYGOUS VARIANT  (FN, FP)\n");
        }
    } else {
        if (params->verboseFalseNegatives || params->verboseFalsePositives) {
            if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1)) {
                st_logDebug("\nERROR: INCORRECT VARIANT - DELETION  (FN, FP)\n");
            } else if (strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) {
                st_logDebug("\nERROR: INCORRECT VARIANT - INSERTION  (FN, FP)\n");
            } else {
                st_logDebug("\nERROR: INCORRECT VARIANT  (FN, FP)\n");
            }
        }
    }
}


static void recordIncorrectVariant(vcfRecordComparisonInfo *vcfInfo, stGenotypeResults *results, stRPHmm *hmm, stGenomeFragment *gF, stSet *reads1, stSet *reads2) {
    /*
     * Records an error where the bases predicted are not consistent with the ref
     */
    recordIncorrectVariantBasic(vcfInfo, results);
    printIncorrectVariantDescription(vcfInfo, hmm->parameters);

    if (hmm->parameters->verboseFalseNegatives || hmm->parameters->verboseFalsePositives) {
        printAlleleInfo(vcfInfo);
        printColumnAtPosition(hmm, vcfInfo->referencePos);
        printPartitionInfo(vcfInfo->referencePos, reads1, reads2, gF);

        size_t indelLen = strlen(vcfInfo->gt_eval_hap1) > strlen(vcfInfo->gt_eval_hap2) ?
                          strlen(vcfInfo->gt_eval_hap1) : strlen(vcfInfo->gt_eval_hap2);

        for (int64_t j = 1; j < indelLen; j++) {
            st_logDebug(" -> pos: %" PRIi64 "\n", vcfInfo->evalPos+j);
            printColumnAtPosition(hmm, vcfInfo->evalPos+j);
            printPartitionInfo(vcfInfo->evalPos +j, reads1, reads2, gF);
        }
    }
}

static void recordIncorrectVariant2(vcfRecordComparisonInfo *vcfInfo, stGenotypeResults *results, stRPHmmParameters *params, stReferencePriorProbs *rProbs1, stReferencePriorProbs *rProbs2) {
    /*
     * Records an error where the bases predicted are not consistent with the ref
     */
    recordIncorrectVariantBasic(vcfInfo, results);
    printIncorrectVariantDescription(vcfInfo, params);

    if (params->verboseFalseNegatives || params->verboseFalsePositives) {
        printAlleleInfo(vcfInfo);
        printPartitionInfo2(vcfInfo->referencePos, rProbs1, rProbs2);

        size_t indelLen = strlen(vcfInfo->gt_eval_hap1) > strlen(vcfInfo->gt_eval_hap2) ?
                          strlen(vcfInfo->gt_eval_hap1) : strlen(vcfInfo->gt_eval_hap2);

        for (int64_t j = 1; j < indelLen; j++) {
            st_logDebug(" -> pos: %" PRIi64 "\n", vcfInfo->evalPos+j);
            printPartitionInfo2(vcfInfo->evalPos+j, rProbs1, rProbs2);
        }
    }
}

static void determinePhasingConsistencyBasic(vcfRecordComparisonInfo *vcfInfo,
                                             float *switchErrorDistance, stGenotypeResults *results) {
    /*
     * Determines error and phasing information for a record in the evaluated vcf
     */
    char *gt_eval_1 = vcfInfo->phasingHap1 ? vcfInfo->gt_eval_hap1 : vcfInfo->gt_eval_hap2;
    char *gt_eval_2 = vcfInfo->phasingHap1 ? vcfInfo->gt_eval_hap2 : vcfInfo->gt_eval_hap1;

    if (strcmp(vcfInfo->gt_ref_hap1, gt_eval_1) == 0 && strcmp(vcfInfo->gt_ref_hap2, gt_eval_2) == 0) {
        (*switchErrorDistance)++;
        results->truePositiveHet++;
        if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1) ||
                strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) {
            results->truePositiveHetIndels++;
        }
        recordTruePositiveBasic(results, vcfInfo);
    }
    else if (strcmp(vcfInfo->gt_ref_hap1, gt_eval_2) == 0 && strcmp(vcfInfo->gt_ref_hap2, gt_eval_1) == 0) {
        results->switchErrors++;
        results->switchErrorDistance += (*switchErrorDistance);
        (*switchErrorDistance) = 0;
        vcfInfo->phasingHap1 = !vcfInfo->phasingHap1;
        vcfInfo->phasingHap2 = !vcfInfo->phasingHap2;
        results->truePositiveHet++;
        if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1)||
                strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) {
            results->truePositiveHetIndels++;
        }
        recordTruePositiveBasic(results, vcfInfo);
    }
    else {
        recordIncorrectVariantBasic(vcfInfo, results);
    }
}

static void determinePhasingConsistency(vcfRecordComparisonInfo *vcfInfo,
                                        float *switchErrorDistance, stRPHmm *hmm, stGenomeFragment *gF,
                                        stSet *reads1, stSet *reads2,
                                        stRPHmmParameters *params, stGenotypeResults *results) {
    /*
     * Determines error and phasing information for a record in the evaluated vcf
     */
    char *gt_eval_1 = vcfInfo->phasingHap1 ? vcfInfo->gt_eval_hap1 : vcfInfo->gt_eval_hap2;
    char *gt_eval_2 = vcfInfo->phasingHap1 ? vcfInfo->gt_eval_hap2 : vcfInfo->gt_eval_hap1;

    if (strcmp(vcfInfo->gt_ref_hap1, gt_eval_1) == 0 && strcmp(vcfInfo->gt_ref_hap2, gt_eval_2) == 0) {
        (*switchErrorDistance)++;
        results->truePositiveHet++;
        if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1) ||
                strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) {
            results->truePositiveHetIndels++;
        }
        recordTruePositive(results, params, vcfInfo, hmm, gF, reads1, reads2);
    }
    else if (strcmp(vcfInfo->gt_ref_hap1, gt_eval_2) == 0 && strcmp(vcfInfo->gt_ref_hap2, gt_eval_1) == 0) {
        results->switchErrors++;
        results->switchErrorDistance += (*switchErrorDistance);
        (*switchErrorDistance) = 0;
        vcfInfo->phasingHap1 = !vcfInfo->phasingHap1;
        vcfInfo->phasingHap2 = !vcfInfo->phasingHap2;
        results->truePositiveHet++;
        if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1) ||
                strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) {
            results->truePositiveHetIndels++;
        }
        recordTruePositive(results, params, vcfInfo, hmm, gF, reads1, reads2);
        st_logDebug("\nSwitch error @ position %" PRIi64 "\n", vcfInfo->referencePos);

    }
    else {
        recordIncorrectVariant(vcfInfo, results, hmm, gF, reads1, reads2);
    }
}

static void determinePhasingConsistency2(vcfRecordComparisonInfo *vcfInfo,
                                        float *switchErrorDistance,
                                         stReferencePriorProbs *rProbs1, stReferencePriorProbs *rProbs2,
                                        stRPHmmParameters *params, stGenotypeResults *results) {
    /*
     * Determines error and phasing information for a record in the evaluated vcf
     */
    char *gt_eval_1 = vcfInfo->phasingHap1 ? vcfInfo->gt_eval_hap1 : vcfInfo->gt_eval_hap2;
    char *gt_eval_2 = vcfInfo->phasingHap1 ? vcfInfo->gt_eval_hap2 : vcfInfo->gt_eval_hap1;

    if (strcmp(vcfInfo->gt_ref_hap1, gt_eval_1) == 0 && strcmp(vcfInfo->gt_ref_hap2, gt_eval_2) == 0) {
        (*switchErrorDistance)++;
        results->truePositiveHet++;
        if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1) ||
                strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) {
            results->truePositiveHetIndels++;
        }
        recordTruePositive2(results, params, vcfInfo, rProbs1, rProbs2);
    }
    else if (strcmp(vcfInfo->gt_ref_hap1, gt_eval_2) == 0 && strcmp(vcfInfo->gt_ref_hap2, gt_eval_1) == 0) {
        results->switchErrors++;
        results->switchErrorDistance += (*switchErrorDistance);
        (*switchErrorDistance) = 0;
        vcfInfo->phasingHap1 = !vcfInfo->phasingHap1;
        vcfInfo->phasingHap2 = !vcfInfo->phasingHap2;
        results->truePositiveHet++;
        if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1) ||
                strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) {
            results->truePositiveHetIndels++;
        }
        recordTruePositive2(results, params, vcfInfo, rProbs1, rProbs2);
        st_logDebug("\nSwitch error @ position %" PRIi64 "\n", vcfInfo->referencePos);

    }
    else {
        recordIncorrectVariant2(vcfInfo, results, params, rProbs1, rProbs2);
    }
}



void printFalsePositive(vcfRecordComparisonInfo *vcfInfo, stRPHmm *hmm,
                        stSet *reads1, stSet *reads2, stGenomeFragment *gF, stBaseMapper *baseMapper) {
    /*
     * Prints out alleles and base compositions for a false positive.
     */
    st_logDebug("\nFALSE POSITIVE");

    size_t indelLen = strlen(vcfInfo->evalRefAllele) > strlen(vcfInfo->evalAltAllele) ?
                      strlen(vcfInfo->evalRefAllele) : strlen(vcfInfo->evalAltAllele);
    char *refSeq = st_calloc(indelLen, sizeof(char));
    for (int64_t i = 0; i < indelLen; i++) {
        refSeq[i] = stBaseMapper_getCharForValue(baseMapper, gF->referenceSequence[vcfInfo->evalPos - gF->refStart + i]);
    }
    // TODO this assumes that the output matches the reference sequence in the case of a FP
    // example: position 4197709,  homozygous variant T->A at position 4197710
    // should output GA, GA (not GT, GT)

    if (indelLen > 1) st_logDebug("  -  INDEL  (FP)");
    else st_logDebug("  -  SNV  (FP)");
    st_logDebug("\n  pos: %" PRIi64 "\n", vcfInfo->evalPos);
    st_logDebug("  baseline  ref: %s\t alt: %s \tphasing: %d | %d\n", refSeq, refSeq,
                0, 0);
    st_logDebug("  output    ref: %s\t alt: %s \tphasing: %d | %d\n",
                vcfInfo->evalRefAllele, vcfInfo->evalAltAllele, vcfInfo->evalPhasing1, vcfInfo->evalPhasing2);

    printColumnAtPosition(hmm, vcfInfo->evalPos);
    printPartitionInfo(vcfInfo->evalPos, reads1, reads2, gF);

    // TODO take out limit on indel size to print
    for (int64_t j = 1; (j < indelLen && j < 5); j++) {
        st_logDebug(" -> pos: %" PRIi64 "\n", vcfInfo->evalPos+j);
        printColumnAtPosition(hmm, vcfInfo->evalPos+j);
        printPartitionInfo(vcfInfo->evalPos+j, reads1, reads2, gF);
    }
    free(refSeq);
}

void printFalsePositive2(vcfRecordComparisonInfo *vcfInfo, stReferencePriorProbs *rProbs1, stReferencePriorProbs *rProbs2, stBaseMapper *baseMapper) {
    /*
     * Prints out alleles and base compositions for a false positive.
     */
    st_logDebug("\nFALSE POSITIVE");

    size_t indelLen = strlen(vcfInfo->evalRefAllele) > strlen(vcfInfo->evalAltAllele) ?
                      strlen(vcfInfo->evalRefAllele) : strlen(vcfInfo->evalAltAllele);
    char *refSeq = st_calloc(indelLen, sizeof(char));
    for (int64_t i = 0; i < indelLen; i++) {
        refSeq[i] = stBaseMapper_getCharForValue(baseMapper, rProbs1->referenceSequence[vcfInfo->evalPos - rProbs1->refStart + i]);
    }
    // TODO this assumes that the output matches the reference sequence in the case of a FP
    // example: position 4197709,  homozygous variant T->A at position 4197710
    // should output GA, GA (not GT, GT)

    if (indelLen > 1) st_logDebug("  -  INDEL  (FP)");
    else st_logDebug("  -  SNV  (FP)s");
    st_logDebug("\n  pos: %" PRIi64 "\n", vcfInfo->evalPos);
    st_logDebug("  baseline  ref: %s\t alt: %s \tphasing: %d | %d\n", refSeq, refSeq,
                0, 0);
    st_logDebug("  output    ref: %s\t alt: %s \tphasing: %d | %d\n",
                vcfInfo->evalRefAllele, vcfInfo->evalAltAllele, vcfInfo->evalPhasing1, vcfInfo->evalPhasing2);
    printPartitionInfo2(vcfInfo->evalPos, rProbs1, rProbs2);

    // TODO take out limit on indel size to print
    for (int64_t j = 1; (j < indelLen && j < 5); j++) {
        st_logDebug(" -> pos: %" PRIi64 "\n", vcfInfo->evalPos+j);
        printPartitionInfo2(vcfInfo->evalPos+j, rProbs1, rProbs2);
    }
    free(refSeq);
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
void compareVCFsBasic(FILE *fh, char *vcf_toEval, char *vcf_ref, stGenotypeResults *results) {

    st_logInfo("> Comparing vcf files \n");
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
    int evalRecordPhased;

    // Iterate through the vcf being checked until getting to the start of the specified interval
    // Don't bother analyzing these records
    int64_t refStart = 0;

    // Variables for keeping track of phasing info
    float switchErrorDistance = 0;

    vcfRecordComparisonInfo *vcfInfo = st_calloc(1, sizeof(vcfRecordComparisonInfo));

    while(bcf_read(inRef, hdrRef, refRecord) == 0) {

        // To take care of the case where a false positive may have been skipped
        // over if the previous eval location was a false negative
        bool maybeFalsePositive = false;
        if (vcfInfo->referencePos < vcfInfo->evalPos) {
            maybeFalsePositive = true;
        }

        // Unpack reference record
        bcf_unpack(refRecord, BCF_UN_ALL);
        vcfInfo->referencePos = refRecord->pos+1;
        vcfInfo->unpackedRecordRef = refRecord;

        // Get genotype info
        int i, j, ngt, nsmpl = bcf_hdr_nsamples(hdrRef);
        int32_t *gt_arr = NULL, ngt_arr = 0;
        ngt = bcf_get_genotypes(hdrRef, vcfInfo->unpackedRecordRef, &gt_arr, &ngt_arr);
        vcfInfo->refPhasing1 = bcf_gt_allele(gt_arr[0]);
        vcfInfo->refPhasing2 = bcf_gt_allele(gt_arr[1]);
        vcfInfo->refAllele = refRecord->d.als;
        vcfInfo->altAllele1 = refRecord->d.allele[1];
        vcfInfo->gt_ref_hap1 = vcfInfo->refPhasing1 == 0 ? refRecord->d.als
                                                         : refRecord->d.allele[vcfInfo->refPhasing1];
        vcfInfo->gt_ref_hap2 = vcfInfo->refPhasing2 == 0 ? refRecord->d.als
                                                         : refRecord->d.allele[vcfInfo->refPhasing2];


        // Skip to the first known location of variation in file being evaluated
        if (results->positives == 0) {
            refStart = vcfInfo->referencePos;
        }

        results->positives++;
        if (vcfInfo->refPhasing1 != vcfInfo->refPhasing2) {
            results->hetsInRef++;
            if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1)) {
                results->hetsInRef_Deletions++;
            }
            else if (strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) {
                results->hetsInRef_Insertions++;
            }
        } else {
            results->homozygousVariantsInRef++;
            if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1)) {
                results->homozygousVariantsInRef_Deletions++;
            }
            else if (strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) {
                results->homozygousVariantsInRef_Insertions++;
            }
        }

        if (maybeFalsePositive && vcfInfo->evalPos < vcfInfo->referencePos) {
            results->falsePositives++;
            if (vcfInfo->evalPhasing1 < 0 || vcfInfo->evalPhasing2 < 0) {
                continue;
            }
            if (strlen(vcfInfo->gt_eval_hap1) > 1 || strlen(vcfInfo->gt_eval_hap2) > 1) {
                results->falsePositiveIndels++;
            }
            if (!vcfInfo->phasingHap1 && !vcfInfo->phasingHap2) {
                if (vcfInfo->evalPhasing1 == 0 && vcfInfo->evalPhasing2 == 1) {
                    vcfInfo->phasingHap1 = true;
                } else if (vcfInfo->evalPhasing1 == 1 && vcfInfo->evalPhasing2 == 0) {
                    vcfInfo->phasingHap2 = true;
                }
            }
        }

        // Iterate through vcf until getting to the position of the variant
        // from the reference vcf currently being looked at
        while (vcfInfo->evalPos < vcfInfo->referencePos) {
            if (bcf_read(inEval, hdrEval, evalRecord) != 0) {
                break;  // can't read record - no more records in file to evaluate
            }
            // Unpack record
            bcf_unpack(evalRecord, BCF_UN_INFO);
            vcfInfo->evalPos = evalRecord->pos+1;
            vcfInfo->unpackedRecordEval = evalRecord;

            if (strcmp(bcf_hdr_id2name(hdrRef, refRecord->rid),
                       bcf_hdr_id2name(hdrEval, evalRecord->rid)) != 0) {
                st_errAbort("ERROR: contig names of records being compared not equal. "
                                    "(contig: %s, pos: %d vs. contig: %s, pos: %d)\n",
                            bcf_hdr_id2name(hdrRef, refRecord->rid), vcfInfo->referencePos,
                            bcf_hdr_id2name(hdrEval, evalRecord->rid), vcfInfo->evalPos);
            }

            // Get genotype info
            int i, j, ngt, nsmpl = bcf_hdr_nsamples(hdrEval);
            int32_t *eval_gt_arr = NULL, eval_ngt_arr = 0;
            ngt = bcf_get_genotypes(hdrEval, vcfInfo->unpackedRecordEval, &eval_gt_arr, &eval_ngt_arr);
            evalRecordPhased = bcf_gt_is_phased(eval_gt_arr[1]);
            if (evalRecordPhased == 0) {
                vcfInfo->phasingHap1 = false;
                vcfInfo->phasingHap2 = false;
            }
            vcfInfo->evalPhasing1 = bcf_gt_allele(eval_gt_arr[0]);
            vcfInfo->evalPhasing2 = bcf_gt_allele(eval_gt_arr[1]);
            vcfInfo->evalRefAllele = evalRecord->d.als;
            vcfInfo->gt_eval_hap1 = vcfInfo->evalPhasing1 == 0 ? evalRecord->d.als
                                                               : evalRecord->d.allele[vcfInfo->evalPhasing1];
            vcfInfo->gt_eval_hap2 = vcfInfo->evalPhasing2 == 0 ? evalRecord->d.als
                                                               : evalRecord->d.allele[vcfInfo->evalPhasing2];


            if (vcfInfo->evalPos < refStart) continue;           // skip this record

            // Check for false positives - variations found not in reference
            if (vcfInfo->evalPos < vcfInfo->referencePos) {
                results->falsePositives++;
                if (vcfInfo->evalPhasing1 < 0 || vcfInfo->evalPhasing2 < 0) {
                    continue;
                }
                if (strlen(vcfInfo->gt_eval_hap1) > 1 || strlen(vcfInfo->gt_eval_hap2) > 1) {
                    results->falsePositiveIndels++;
                }
                if (!vcfInfo->phasingHap1 && !vcfInfo->phasingHap2) {
                    if (vcfInfo->evalPhasing1 == 0 && vcfInfo->evalPhasing2 == 1) {
                        vcfInfo->phasingHap1 = true;
                    } else if (vcfInfo->evalPhasing1 == 1 && vcfInfo->evalPhasing2 == 0) {
                        vcfInfo->phasingHap2 = true;
                    }
                }
            }
        }

        // At locus of known variation
        if (vcfInfo->evalPos == vcfInfo->referencePos) {
            if (vcfInfo->evalPhasing1 < 0 || vcfInfo->evalPhasing2 < 0) {
                // Genotype unknown
                if (vcfInfo->evalPos < vcfInfo->referencePos) {
                    results->falsePositives++;
                } else {
                    results->falseNegatives++;
                }
            }
            // Reference variant 1/1
            else if (vcfInfo->refPhasing1 == vcfInfo->refPhasing2) {
                // Variant being evaluated also 1/1 - homozygous
                if (vcfInfo->evalPhasing1 == vcfInfo->evalPhasing2) {
                    if (strcmp(vcfInfo->gt_ref_hap1, vcfInfo->gt_eval_hap1) == 0
                        && strcmp(vcfInfo->gt_ref_hap2, vcfInfo->gt_eval_hap2) == 0) {
                        // True positive homozygous variant
                        results->truePositiveHomozygous++;
                        if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1) ||
                                strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) {
                            results->truePositiveHomozygousIndels++;
                        }
                        recordTruePositiveBasic(results, vcfInfo);
                    }
                    else {
                        // Predicted homozygous, but doesn't match reference
                        // This could also be an "incorrect" variant...
                        recordHomozygousVariantBasic(results, vcfInfo);
                    }
                } else {
                    // Heterozygous variant in evaluated vcf :/
                    if ((strcmp(vcfInfo->gt_ref_hap1, vcfInfo->gt_eval_hap1) == 0
                         && strcmp(vcfInfo->gt_ref_hap2, vcfInfo->gt_eval_hap2) == 0)
                        || (strcmp(vcfInfo->gt_ref_hap1, vcfInfo->gt_eval_hap2) == 0
                            && strcmp(vcfInfo->gt_eval_hap2, vcfInfo->gt_eval_hap1) == 0)) {

                        // Correct bases, but should have been homozygous
                        recordHomozygousVariantBasic(results, vcfInfo);

                    } else {
                        // Incorrect prediction
                        recordIncorrectVariantBasic(vcfInfo, results);
                    }
                }

            } else if (!vcfInfo->phasingHap1 && !vcfInfo->phasingHap2) {
                // Beginning of genotype fragment, figure out which haplotype matches with ref
                results->uncertainPhasing++;
                if ((strcmp(vcfInfo->gt_ref_hap1,vcfInfo-> gt_eval_hap1) == 0 &&
                     strcmp(vcfInfo->gt_ref_hap2, vcfInfo->gt_eval_hap2) == 0) ||
                    (strcmp(vcfInfo->gt_ref_hap1, vcfInfo->gt_eval_hap2) == 0 &&
                     strcmp(vcfInfo->gt_ref_hap2, vcfInfo->gt_eval_hap1) == 0)) {
                    results->truePositiveHet++;
                    if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1) ||
                            strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) {
                        results->truePositiveHetIndels++;
                    }
                    recordTruePositiveBasic(results, vcfInfo);
                } else {
                    results->falsePositives++;
                    if (strlen(vcfInfo->gt_eval_hap1) > 1 || strlen(vcfInfo->gt_eval_hap2) > 1) {
                        results->falsePositiveIndels++;
                    }
                }
                if (vcfInfo->evalPhasing1 < vcfInfo->evalPhasing2) {
                    vcfInfo->phasingHap1 = true;
                } else if (vcfInfo->evalPhasing1 > vcfInfo->evalPhasing2) {
                    vcfInfo->phasingHap2 = true;
                }
            } else {
                determinePhasingConsistencyBasic(vcfInfo, &switchErrorDistance, results);
            }
        } else if (vcfInfo->evalPos > vcfInfo->referencePos) {
            // False negative - no variation was found, but truth vcf has one
            if (vcfInfo->refPhasing1 != vcfInfo->refPhasing2){
                results->error_missedHet++;
                results->falseNegatives++;

                if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1)) {
                    results->error_missedHet_Deletions++;
                } else if (strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) {
                    results->error_missedHet_Insertions++;
                }
            } else {
                recordHomozygousVariantMissBasic(results, vcfInfo);
            }
            free(gt_arr);
        }
    }
    if (results->truePositives == 0) st_logInfo("No matches between vcfs found - did you compare against the correct vcf?\n");

    // Remaining positions after the last variant in the reference are not currently being looked through
    // False positives in this region could therefore be missed
    // (In addition to false positives after the first variant)
    results->negatives += (vcfInfo->referencePos - refStart - results->positives);
    results->trueNegatives += (results->negatives - results->falsePositives);
    results->switchErrorDistance = results->switchErrorDistance/results->switchErrors;

    // cleanup
    vcf_close(inRef);
    vcf_close(inEval);
    bcf_hdr_destroy(hdrRef);
    bcf_hdr_destroy(hdrEval);
    bcf_destroy(refRecord);
    bcf_destroy(evalRecord);
    free(vcfInfo);

}



/*
     * Test to compare a vcf to a truth vcf containing known variants for the region.
     *
     * Information about some of the results saved in the genotypeResults struct
     *
     */
void compareVCFs(FILE *fh, stList *hmms, char *vcf_toEval, char *vcf_ref,
                 stBaseMapper *baseMapper, stGenotypeResults *results, stRPHmmParameters *params) {

    st_logInfo("> Comparing vcf files \n");
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
    int evalRecordPhased;

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

    // Variables for keeping track of phasing info
    float switchErrorDistance = 0;

    vcfRecordComparisonInfo *vcfInfo = st_calloc(1, sizeof(vcfRecordComparisonInfo));

    while(bcf_read(inRef, hdrRef, refRecord) == 0) {

        // To take care of the case where a false positive may have been skipped
        // over if the previous eval location was a false negative
        bool maybeFalsePositive = false;
        if (vcfInfo->referencePos < vcfInfo->evalPos) {
            maybeFalsePositive = true;
        }

        // Unpack reference record
        bcf_unpack(refRecord, BCF_UN_ALL);
        vcfInfo->referencePos = refRecord->pos+1;
        vcfInfo->unpackedRecordRef = refRecord;

        // Get genotype info
        int i, j, ngt, nsmpl = bcf_hdr_nsamples(hdrRef);
        int32_t *gt_arr = NULL, ngt_arr = 0;
        ngt = bcf_get_genotypes(hdrRef, vcfInfo->unpackedRecordRef, &gt_arr, &ngt_arr);
        vcfInfo->refPhasing1 = bcf_gt_allele(gt_arr[0]);
        vcfInfo->refPhasing2 = bcf_gt_allele(gt_arr[1]);
        vcfInfo->refAllele = refRecord->d.als;
        vcfInfo->altAllele1 = refRecord->d.allele[1];
        vcfInfo->gt_ref_hap1 = vcfInfo->refPhasing1 == 0 ? refRecord->d.als
                                                        : refRecord->d.allele[vcfInfo->refPhasing1];
        vcfInfo->gt_ref_hap2 = vcfInfo->refPhasing2 == 0 ? refRecord->d.als
                                                        : refRecord->d.allele[vcfInfo->refPhasing2];

        // Skip to the first known location of variation in file being evaluated
        if (results->positives == 0) {
            refStart = vcfInfo->referencePos;
        }

        // Make sure to only look at records in the specified interval
        if (vcfInfo->referencePos < hmm->refStart) continue;

        // If the position is beyond the end of this hmm, get the next one
        while ((hmm->refStart + hmm->refLength) < vcfInfo->referencePos) {
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
                vcfInfo->phasingHap1 = false;
                vcfInfo->phasingHap2 = false;
            } else {
                break;
            }
        }
        // No more fragments to look through
        if (hmmIndex == stList_length(hmms)) break;

        vcfInfo->h1AlphChar = stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString1[vcfInfo->referencePos - gF->refStart]);
        vcfInfo->h2AlphChar = stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[vcfInfo->referencePos - gF->refStart]);

        results->positives++;
        if (vcfInfo->refPhasing1 != vcfInfo->refPhasing2) {
            results->hetsInRef++;
            if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1)) {
                results->hetsInRef_Deletions++;
            }
            else if (strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) {
                results->hetsInRef_Insertions++;
            }
        } else {
            results->homozygousVariantsInRef++;
            if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1)) {
                results->homozygousVariantsInRef_Deletions++;
            }
            else if (strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) {
                results->homozygousVariantsInRef_Insertions++;
            }
        }

        if (maybeFalsePositive && vcfInfo->evalPos < vcfInfo->referencePos) {
            results->falsePositives++;
            if (params->verboseFalsePositives) {
                printFalsePositive(vcfInfo, hmm, reads1, reads2, gF, baseMapper);
            }
            if (vcfInfo->evalPhasing1 < 0 || vcfInfo->evalPhasing2 < 0) {
                continue;
            }
            if (strlen(vcfInfo->evalRefAllele) > 1 || strlen(vcfInfo->evalAltAllele) > 1) {
                results->falsePositiveIndels++;
            }
            if (!vcfInfo->phasingHap1 && !vcfInfo->phasingHap2) {
                if (vcfInfo->evalPhasing1 == 0 && vcfInfo->evalPhasing2 == 1) {
                    vcfInfo->phasingHap1 = true;
                } else if (vcfInfo->evalPhasing1 == 1 && vcfInfo->evalPhasing2 == 0) {
                    vcfInfo->phasingHap2 = true;
                }
            }
        }

        // Iterate through vcf until getting to the position of the variant
        // from the reference vcf currently being looked at
        while (vcfInfo->evalPos < vcfInfo->referencePos) {
            if (bcf_read(inEval, hdrEval, evalRecord) != 0) {
                break;  // can't read record - no more records in file to evaluate
            }
            // Unpack record
            bcf_unpack(evalRecord, BCF_UN_INFO);
            vcfInfo->evalPos = evalRecord->pos+1;
            vcfInfo->unpackedRecordEval = evalRecord;

            if (strcmp(bcf_hdr_id2name(hdrRef, refRecord->rid),
                       bcf_hdr_id2name(hdrEval, evalRecord->rid)) != 0) {
                st_errAbort("ERROR: contig names of records being compared not equal. "
                                    "(contig: %s, pos: %d vs. contig: %s, pos: %d)\n",
                            bcf_hdr_id2name(hdrRef, refRecord->rid), vcfInfo->referencePos,
                            bcf_hdr_id2name(hdrEval, evalRecord->rid), vcfInfo->evalPos);
            }

            // Get genotype info
            int i, j, ngt, nsmpl = bcf_hdr_nsamples(hdrEval);
            int32_t *eval_gt_arr = NULL, eval_ngt_arr = 0;
            ngt = bcf_get_genotypes(hdrEval, vcfInfo->unpackedRecordEval, &eval_gt_arr, &eval_ngt_arr);
            vcfInfo->evalPhasing1 = bcf_gt_allele(eval_gt_arr[0]);
            vcfInfo->evalPhasing2 = bcf_gt_allele(eval_gt_arr[1]);
            evalRecordPhased = bcf_gt_is_phased(eval_gt_arr[1]);

            vcfInfo->evalRefAllele = evalRecord->d.als;
            vcfInfo->evalAltAllele = evalRecord->d.allele[1];
            vcfInfo->gt_eval_hap1 = vcfInfo->evalPhasing1 == 0 ? evalRecord->d.als
                                                              : evalRecord->d.allele[vcfInfo->evalPhasing1];
            vcfInfo->gt_eval_hap2 = vcfInfo->evalPhasing2 == 0 ? evalRecord->d.als
                                                              : evalRecord->d.allele[vcfInfo->evalPhasing2];

            if (vcfInfo->evalPos < refStart) continue;           // skip this record

            // Check for false positives - variations found not in reference
            if (vcfInfo->evalPos < vcfInfo->referencePos) {
                results->falsePositives++;
                if (params->verboseFalsePositives) {
                    printFalsePositive(vcfInfo, hmm, reads1, reads2, gF, baseMapper);
                }
                if (vcfInfo->evalPhasing1 < 0 || vcfInfo->evalPhasing2 < 0) {
                    continue;
                }
                if (strlen(vcfInfo->evalRefAllele) > 1 || strlen(vcfInfo->evalAltAllele) > 1) {
                    results->falsePositiveIndels++;
                }
                if (!vcfInfo->phasingHap1 && !vcfInfo->phasingHap2) {
                    if (vcfInfo->evalPhasing1 == 0 && vcfInfo->evalPhasing2 == 1) {
                        vcfInfo->phasingHap1 = true;
                    } else if (vcfInfo->evalPhasing1 == 1 && vcfInfo->evalPhasing2 == 0) {
                        vcfInfo->phasingHap2 = true;
                    }
                }
            } else {
                break;
            }
        }

        // At locus of known variation
        if (vcfInfo->evalPos == vcfInfo->referencePos) {
            if (vcfInfo->evalPhasing1 < 0 || vcfInfo->evalPhasing2 < 0) {
                // Genotype unknown
                if (vcfInfo->evalPos < vcfInfo->referencePos) {
                    results->falsePositives++;
                } else {
                    results->falseNegatives++;
                }
            }
            // Reference variant 1/1
            else if (vcfInfo->refPhasing1 == vcfInfo->refPhasing2) {
                // Variant being evaluated also 1/1 - homozygous
                if (vcfInfo->evalPhasing1 == vcfInfo->evalPhasing2) {
                    if (strcmp(vcfInfo->gt_ref_hap1, vcfInfo->gt_eval_hap1) == 0
                        && strcmp(vcfInfo->gt_ref_hap2, vcfInfo->gt_eval_hap2) == 0) {
                        // True positive homozygous variant
                        results->truePositiveHomozygous++;
                        if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1) ||
                                strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) {
                            results->truePositiveHomozygousIndels++;
                        }
                        recordTruePositive(results, params, vcfInfo, hmm, gF, reads1, reads2);
                    } else {
                        // Predicted homozygous, but doesn't match reference
                        // This could also be an "incorrect" variant...
                        if (params->verboseFalseNegatives || params->verboseFalsePositives) {
                            st_logDebug("\nMISS - INCORRECT VARIANT, HOMOZYGOUS  (FN, FP)\n");
                        }
                        recordHomozygousVariant(results, vcfInfo, baseMapper, hmm, gF, reads1, reads2);
                    }
                } else {
                    // Heterozygous variant in evaluated vcf :/
                    if ((strcmp(vcfInfo->gt_ref_hap1, vcfInfo->gt_eval_hap1) == 0
                         && strcmp(vcfInfo->gt_ref_hap2, vcfInfo->gt_eval_hap2) == 0)
                        || (strcmp(vcfInfo->gt_ref_hap1, vcfInfo->gt_eval_hap2) == 0
                            && strcmp(vcfInfo->gt_eval_hap2, vcfInfo->gt_eval_hap1) == 0)) {

                        // Correct bases, but should have been homozygous
                        if (params->verboseFalseNegatives || params->verboseFalsePositives) {
                            st_logDebug("\nMISS - PREDICTED HET VARIANT, HOMOZYGOUS IN REFERENCE  (FN, FP)\n");
                        }
                        recordHomozygousVariant(results, vcfInfo, baseMapper, hmm, gF, reads1, reads2);

                    } else {
                        // Incorrect prediction
                        recordIncorrectVariant(vcfInfo, results, hmm, gF, reads1, reads2);
                    }
                }

            } else if (!vcfInfo->phasingHap1 && !vcfInfo->phasingHap2) {

                // Beginning of genotype fragment, figure out which haplotype matches with ref
                results->uncertainPhasing++;
                if ((strcmp(vcfInfo->gt_ref_hap1,vcfInfo-> gt_eval_hap1) == 0 &&
                        strcmp(vcfInfo->gt_ref_hap2, vcfInfo->gt_eval_hap2) == 0) ||
                        (strcmp(vcfInfo->gt_ref_hap1, vcfInfo->gt_eval_hap2) == 0 &&
                                strcmp(vcfInfo->gt_ref_hap2, vcfInfo->gt_eval_hap1) == 0)) {
                    results->truePositiveHet++;
                    if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1) ||
                            strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) {
                        results->truePositiveHetIndels++;
                    }
                    recordTruePositive(results, params, vcfInfo, hmm, gF, reads1, reads2);
                } else {
                    results->falsePositives++;
                    if (strlen(vcfInfo->evalRefAllele) > 1 || strlen(vcfInfo->evalAltAllele) > 1) {
                        results->falsePositiveIndels++;
                    }
                    if (params->verboseFalsePositives) {
                        printFalsePositive(vcfInfo, hmm, reads1, reads2, gF, baseMapper);
                    }
                }
                if (vcfInfo->evalPhasing1 < vcfInfo->evalPhasing2) {
                    vcfInfo->phasingHap1 = true;
                } else if (vcfInfo->evalPhasing1 > vcfInfo->evalPhasing2) {
                    vcfInfo->phasingHap2 = true;
                }
            } else {
                determinePhasingConsistency(vcfInfo, &switchErrorDistance,
                                   hmm, gF, reads1, reads2, params, results);
            }
        } else if (vcfInfo->evalPos > vcfInfo->referencePos){
            // Missed the variant - False negative

            double *read1BaseCounts = getProfileSequenceBaseCompositionAtPosition(reads1, vcfInfo->referencePos);
            double *read2BaseCounts = getProfileSequenceBaseCompositionAtPosition(reads2, vcfInfo->referencePos);

            // False negative - no variation was found, but truth vcf has one
            if (vcfInfo->refPhasing1 != vcfInfo->refPhasing2){
                results->error_missedHet++;
                results->falseNegatives++;

                if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1) ||
                        strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) {
                    // Indel
                    if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1)) {
                        results->error_missedHet_Deletions++;
                    } else {
                        results->error_missedHet_Insertions++;
                    }

                    if (params->verboseFalseNegatives) {
                        size_t indelLen = strlen(vcfInfo->gt_ref_hap1) > strlen(vcfInfo->gt_ref_hap2)
                                          ? strlen(vcfInfo->gt_ref_hap1) : strlen(vcfInfo->gt_ref_hap2);
                        char *refSeq = st_calloc(indelLen, sizeof(char));
                        for (int64_t i = 0; i < indelLen; i++) {
                            refSeq[i] = stBaseMapper_getCharForValue(baseMapper, gF->referenceSequence[vcfInfo->referencePos - gF->refStart + i]);
                        }
                        char refChar = stBaseMapper_getCharForValue(baseMapper, gF->referenceSequence[vcfInfo->referencePos - gF->refStart]);

                        if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1)) {
                            st_logDebug("\nMISS  -  DELETION  (FN)\n");
                        } else {
                            st_logDebug("\nMISS  -  INSERTION  (FN)\n");
                        }
                        st_logDebug("  pos: %" PRIi64 "\n", vcfInfo->referencePos);
                        st_logDebug("  baseline  ref: %s\t alt: ", vcfInfo->refAllele);
                        for (int i = 1; i < vcfInfo->unpackedRecordRef->n_allele; i++) {
                            if (i != 1) st_logDebug(",");
                            st_logDebug("%s", vcfInfo->unpackedRecordRef->d.allele[i]);
                        }
                        st_logDebug(" \tphasing: %d | %d", vcfInfo->refPhasing1, vcfInfo->refPhasing2);
                        if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1)) {
                            st_logDebug("\n  output    ref: %s\t alt: %s \tphasing: %d | %d\n",
                                        refSeq, refSeq, 0, 0);
                        } else {
                            st_logDebug("\n  output    ref: %c\t alt: %c \tphasing: %d | %d\n",
                                        refChar, refChar, 0, 0);
                        }
                        printColumnAtPosition(hmm, vcfInfo->referencePos);
                        printPartitionInfo(vcfInfo->referencePos, reads1, reads2, gF);

                        for (int64_t j = 1; j < indelLen; j++) {
                            st_logDebug(" -> pos: %" PRIi64 "\n", vcfInfo->referencePos+j);
                            printColumnAtPosition(hmm, vcfInfo->referencePos+j);
                            printPartitionInfo(vcfInfo->referencePos+j, reads1, reads2, gF);
                        }
                        free(refSeq);
                    }
                } else {
                    // SNV
                    if (params->verboseFalseNegatives) {
                        st_logDebug("\nMISS  -  SNV  (FN)\n");
                        st_logDebug("  pos: %" PRIi64 "\n", vcfInfo->referencePos);
                        st_logDebug("  baseline  ref: %s\t alt: ", vcfInfo->refAllele);

                        for (int i = 1; i < vcfInfo->unpackedRecordRef->n_allele; i++) {
                            if (i != 1) st_logDebug(",");
                            st_logDebug("%s", vcfInfo->unpackedRecordRef->d.allele[i]);
                        }
                        st_logDebug(" \tphasing: %d | %d", vcfInfo->refPhasing1, vcfInfo->refPhasing2);
                        st_logDebug("\n  output    ref: %c\t alt: %c \tphasing: %d | %d\n",
                                    vcfInfo->h1AlphChar, vcfInfo->h2AlphChar, 0, 0);
                        printColumnAtPosition(hmm, vcfInfo->referencePos);

                        st_logDebug("\tPartition 1: \n");
                        printBaseComposition2(read1BaseCounts);
                        st_logDebug("\tPartition 2: \n");
                        printBaseComposition2(read2BaseCounts);
                        st_logDebug("\t\t\tposterior prob: %f\n",
                                    gF->genotypeProbs[vcfInfo->referencePos-gF->refStart]);
                    }
                }
            } else {
                recordHomozygousVariantMiss(results, vcfInfo, baseMapper, hmm, gF, reads1, reads2);
            }

            free(read1BaseCounts);
            free(read2BaseCounts);
            free(gt_arr);
        }
    }
    if (results->truePositives == 0) st_logInfo("No matches between vcfs found - did you compare against the correct vcf?\n");

    // Remaining positions after the last variant in the reference are not currently being looked through
    // False positives in this region could therefore be missed
    // (In addition to false positives after the first variant)
    results->negatives += (vcfInfo->referencePos - refStart - results->positives);
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
    free(vcfInfo);
}


/*
     * Test to compare a vcf to a truth vcf containing known variants for the region.
     *
     * Information about some of the results saved in the genotypeResults struct.
     *
     * This version will do output debugging using the bam files output from marginPhase.
     *
     */
void compareVCFs_debugWithBams(char *vcf_toEval, char *vcf_ref, char *bamFile1, char *bamFile2, char *referenceFasta, stBaseMapper *baseMapper, stGenotypeResults *results, stRPHmmParameters *params) {

    st_logInfo("> Comparing vcf files \n");
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

    // Parse reads
    st_logInfo("> Parsing input reads from file: %s\n", bamFile1);
    stList *profileSequences1 = stList_construct3(0, (void (*)(void *))stProfileSeq_destruct);
    int64_t readCount1 = parseReads(profileSequences1, bamFile1, baseMapper, params);
    st_logInfo("\tCreated %d profile sequences\n", readCount1);
    stHash *referenceNamesToReferencePriors1 = createReferencePriorProbabilities(referenceFasta, profileSequences1, baseMapper, params);
    stHashIterator *hashIt1 = stHash_getIterator(referenceNamesToReferencePriors1);
    stReferencePriorProbs *rProbs1;
    char *referenceName1;
    if ((referenceName1 = stHash_getNext(hashIt1)) != NULL) {
        rProbs1 = stHash_search(referenceNamesToReferencePriors1, referenceName1);
    }

    st_logInfo("> Parsing input reads from file: %s\n", bamFile2);
    stList *profileSequences2 = stList_construct3(0, (void (*)(void *))stProfileSeq_destruct);
    int64_t readCount2 = parseReads(profileSequences2, bamFile2, baseMapper, params);
    st_logInfo("\tCreated %d profile sequences\n", readCount2);
    stHash *referenceNamesToReferencePriors2 = createReferencePriorProbabilities(referenceFasta, profileSequences2, baseMapper, params);
    stHashIterator *hashIt2 = stHash_getIterator(referenceNamesToReferencePriors2);
    stReferencePriorProbs *rProbs2;
    char *referenceName2;
    if ((referenceName2 = stHash_getNext(hashIt2)) != NULL) {
        rProbs2 = stHash_search(referenceNamesToReferencePriors2, referenceName2);
    }

    // Variables to keep track of some stats
    int64_t refStart = 0;
    float switchErrorDistance = 0;
    vcfRecordComparisonInfo *vcfInfo = st_calloc(1, sizeof(vcfRecordComparisonInfo));

    while(bcf_read(inRef, hdrRef, refRecord) == 0) {

        // To take care of the case where a false positive may have been skipped
        // over if the previous eval location was a false negative
        bool maybeFalsePositive = false;
        if (vcfInfo->referencePos < vcfInfo->evalPos) {
            maybeFalsePositive = true;
        }

        // Unpack reference record
        bcf_unpack(refRecord, BCF_UN_ALL);
        vcfInfo->referencePos = refRecord->pos+1;
        vcfInfo->unpackedRecordRef = refRecord;

        // Get genotype info
        int i, j, ngt, nsmpl = bcf_hdr_nsamples(hdrRef);
        int32_t *gt_arr = NULL, ngt_arr = 0;
        ngt = bcf_get_genotypes(hdrRef, vcfInfo->unpackedRecordRef, &gt_arr, &ngt_arr);
        vcfInfo->refPhasing1 = bcf_gt_allele(gt_arr[0]);
        vcfInfo->refPhasing2 = bcf_gt_allele(gt_arr[1]);
        vcfInfo->refAllele = refRecord->d.als;
        vcfInfo->altAllele1 = refRecord->d.allele[1];
        vcfInfo->gt_ref_hap1 = vcfInfo->refPhasing1 == 0 ? refRecord->d.als
                                                         : refRecord->d.allele[vcfInfo->refPhasing1];
        vcfInfo->gt_ref_hap2 = vcfInfo->refPhasing2 == 0 ? refRecord->d.als
                                                         : refRecord->d.allele[vcfInfo->refPhasing2];

        // Skip to the first known location of variation in file being evaluated
        if (results->positives == 0) {
            refStart = vcfInfo->referencePos;
        }

        results->positives++;
        if (vcfInfo->refPhasing1 != vcfInfo->refPhasing2) {
            results->hetsInRef++;
            if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1)) {
                results->hetsInRef_Deletions++;
            }
            else if (strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) {
                results->hetsInRef_Insertions++;
            }
        } else {
            results->homozygousVariantsInRef++;
            if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1)) {
                results->homozygousVariantsInRef_Deletions++;
            }
            else if (strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) {
                results->homozygousVariantsInRef_Insertions++;
            }
        }

        if (maybeFalsePositive && vcfInfo->evalPos < vcfInfo->referencePos) {
            results->falsePositives++;
            if (params->verboseFalsePositives) {
                printFalsePositive2(vcfInfo, rProbs1, rProbs2, baseMapper);
            }
            if (vcfInfo->evalPhasing1 < 0 || vcfInfo->evalPhasing2 < 0) {
                continue;
            }
            if (strlen(vcfInfo->evalRefAllele) > 1 || strlen(vcfInfo->evalAltAllele) > 1) {
                results->falsePositiveIndels++;
            }
            if (!vcfInfo->phasingHap1 && !vcfInfo->phasingHap2) {
                if (vcfInfo->evalPhasing1 == 0 && vcfInfo->evalPhasing2 == 1) {
                    vcfInfo->phasingHap1 = true;
                } else if (vcfInfo->evalPhasing1 == 1 && vcfInfo->evalPhasing2 == 0) {
                    vcfInfo->phasingHap2 = true;
                }
            }
        }

        // Iterate through vcf until getting to the position of the variant
        // from the reference vcf currently being looked at
        while (vcfInfo->evalPos < vcfInfo->referencePos) {
            if (bcf_read(inEval, hdrEval, evalRecord) != 0) {
                break;  // can't read record - no more records in file to evaluate
            }
            // Unpack record
            bcf_unpack(evalRecord, BCF_UN_INFO);
            vcfInfo->evalPos = evalRecord->pos+1;
            vcfInfo->unpackedRecordEval = evalRecord;

            if (strcmp(bcf_hdr_id2name(hdrRef, refRecord->rid),
                       bcf_hdr_id2name(hdrEval, evalRecord->rid)) != 0) {
                st_errAbort("ERROR: contig names of records being compared not equal. "
                                    "(contig: %s, pos: %d vs. contig: %s, pos: %d)\n",
                            bcf_hdr_id2name(hdrRef, refRecord->rid), vcfInfo->referencePos,
                            bcf_hdr_id2name(hdrEval, evalRecord->rid), vcfInfo->evalPos);
            }

            // Get genotype info
            int i, j, ngt, nsmpl = bcf_hdr_nsamples(hdrEval);
            int32_t *eval_gt_arr = NULL, eval_ngt_arr = 0;
            ngt = bcf_get_genotypes(hdrEval, vcfInfo->unpackedRecordEval, &eval_gt_arr, &eval_ngt_arr);
            vcfInfo->evalPhasing1 = bcf_gt_allele(eval_gt_arr[0]);
            vcfInfo->evalPhasing2 = bcf_gt_allele(eval_gt_arr[1]);

            vcfInfo->evalRefAllele = evalRecord->d.als;
            vcfInfo->evalAltAllele = evalRecord->d.allele[1];
            vcfInfo->gt_eval_hap1 = vcfInfo->evalPhasing1 == 0 ? evalRecord->d.als
                                                               : evalRecord->d.allele[vcfInfo->evalPhasing1];
            vcfInfo->gt_eval_hap2 = vcfInfo->evalPhasing2 == 0 ? evalRecord->d.als
                                                               : evalRecord->d.allele[vcfInfo->evalPhasing2];

            if (vcfInfo->evalPos < refStart) continue;           // skip this record

            // Check for false positives - variations found not in reference
            if (vcfInfo->evalPos < vcfInfo->referencePos) {
                results->falsePositives++;
                if (params->verboseFalsePositives) {
                    printFalsePositive2(vcfInfo, rProbs1, rProbs2, baseMapper);
                }
                if (vcfInfo->evalPhasing1 < 0 || vcfInfo->evalPhasing2 < 0) {
                    continue;
                }
                if (strlen(vcfInfo->evalRefAllele) > 1 || strlen(vcfInfo->evalAltAllele) > 1) {
                    results->falsePositiveIndels++;
                }
                if (!vcfInfo->phasingHap1 && !vcfInfo->phasingHap2) {
                    if (vcfInfo->evalPhasing1 == 0 && vcfInfo->evalPhasing2 == 1) {
                        vcfInfo->phasingHap1 = true;
                    } else if (vcfInfo->evalPhasing1 == 1 && vcfInfo->evalPhasing2 == 0) {
                        vcfInfo->phasingHap2 = true;
                    }
                }
            } else {
                break;
            }
        }

        // At locus of known variation
        if (vcfInfo->evalPos == vcfInfo->referencePos) {
            if (vcfInfo->evalPhasing1 < 0 || vcfInfo->evalPhasing2 < 0) {
                // Genotype unknown
                if (vcfInfo->evalPos < vcfInfo->referencePos) {
                    results->falsePositives++;
                } else {
                    results->falseNegatives++;
                }
            }
            // Reference variant 1/1
            else if (vcfInfo->refPhasing1 == vcfInfo->refPhasing2) {
                // Variant being evaluated also 1/1 - homozygous
                if (vcfInfo->evalPhasing1 == vcfInfo->evalPhasing2) {
                    if (strcmp(vcfInfo->gt_ref_hap1, vcfInfo->gt_eval_hap1) == 0
                        && strcmp(vcfInfo->gt_ref_hap2, vcfInfo->gt_eval_hap2) == 0) {
                        // True positive homozygous variant
                        results->truePositiveHomozygous++;
                        if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1) ||
                                strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) {
                            results->truePositiveHomozygousIndels++;
                        }
                        recordTruePositive2(results, params, vcfInfo, rProbs1, rProbs2);
                    } else {
                        // Predicted homozygous, but doesn't match reference
                        // This could also be an "incorrect" variant...
                        if (params->verboseFalseNegatives || params->verboseFalsePositives) {
                            st_logDebug("\nMISS - INCORRECT VARIANT, HOMOZYGOUS  (FN, FP)\n");
                        }
                        recordHomozygousVariant2(results, vcfInfo, baseMapper, params, rProbs1, rProbs2);
                    }
                } else {
                    // Heterozygous variant in evaluated vcf :/
                    if ((strcmp(vcfInfo->gt_ref_hap1, vcfInfo->gt_eval_hap1) == 0
                         && strcmp(vcfInfo->gt_ref_hap2, vcfInfo->gt_eval_hap2) == 0)
                        || (strcmp(vcfInfo->gt_ref_hap1, vcfInfo->gt_eval_hap2) == 0
                            && strcmp(vcfInfo->gt_eval_hap2, vcfInfo->gt_eval_hap1) == 0)) {

                        // Correct bases, but should have been homozygous
                        if (params->verboseFalseNegatives || params->verboseFalsePositives) {
                            st_logDebug("\nMISS - PREDICTED HET VARIANT, HOMOZYGOUS IN REFERENCE  (FN, FP)\n");
                        }
                        recordHomozygousVariant2(results, vcfInfo, baseMapper, params, rProbs1, rProbs2);

                    } else {
                        // Incorrect prediction
                        recordIncorrectVariant2(vcfInfo, results, params, rProbs1, rProbs2);
                    }
                }

            } else if (!vcfInfo->phasingHap1 && !vcfInfo->phasingHap2) {

                // Beginning of genotype fragment, figure out which haplotype matches with ref
                results->uncertainPhasing++;
                if ((strcmp(vcfInfo->gt_ref_hap1,vcfInfo-> gt_eval_hap1) == 0 &&
                     strcmp(vcfInfo->gt_ref_hap2, vcfInfo->gt_eval_hap2) == 0) ||
                    (strcmp(vcfInfo->gt_ref_hap1, vcfInfo->gt_eval_hap2) == 0 &&
                     strcmp(vcfInfo->gt_ref_hap2, vcfInfo->gt_eval_hap1) == 0)) {
                    results->truePositiveHet++;
                    if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1) ||
                            strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) {
                        results->truePositiveHetIndels++;
                    }
                    recordTruePositive2(results, params, vcfInfo, rProbs1, rProbs2);
                } else {
                    results->falsePositives++;
                    if (strlen(vcfInfo->evalRefAllele) > 1 || strlen(vcfInfo->evalAltAllele) > 1) {
                        results->falsePositiveIndels++;
                    }
                    if (params->verboseFalsePositives) {
                        printFalsePositive2(vcfInfo, rProbs1, rProbs2, baseMapper);
                    }
                }
                if (vcfInfo->evalPhasing1 < vcfInfo->evalPhasing2) {
                    vcfInfo->phasingHap1 = true;
                } else if (vcfInfo->evalPhasing1 > vcfInfo->evalPhasing2) {
                    vcfInfo->phasingHap2 = true;
                }
            } else {
                determinePhasingConsistency2(vcfInfo, &switchErrorDistance, rProbs1, rProbs2, params, results);
            }
        } else if (vcfInfo->evalPos > vcfInfo->referencePos){
            // False negative - no variation was found, but truth vcf has one
            if (vcfInfo->refPhasing1 != vcfInfo->refPhasing2){
                results->error_missedHet++;
                results->falseNegatives++;

                if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1) ||
                        strlen(vcfInfo->refAllele) < strlen(vcfInfo->altAllele1)) {
                    // Indel
                    if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1)) {
                        results->error_missedHet_Deletions++;
                    } else {
                        results->error_missedHet_Insertions++;
                    }

                    if (params->verboseFalseNegatives) {
                        size_t indelLen = strlen(vcfInfo->gt_ref_hap1) > strlen(vcfInfo->gt_ref_hap2)
                                          ? strlen(vcfInfo->gt_ref_hap1) : strlen(vcfInfo->gt_ref_hap2);
                        char *refSeq = st_calloc(indelLen, sizeof(char));
                        for (int64_t i = 0; i < indelLen; i++) {
                            refSeq[i] = stBaseMapper_getCharForValue(baseMapper, rProbs1->referenceSequence[vcfInfo->referencePos - rProbs1->refStart + i]);
                        }
                        char refChar = stBaseMapper_getCharForValue(baseMapper, rProbs1->referenceSequence[vcfInfo->referencePos - rProbs1->refStart]);

                        if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1)) {
                            st_logDebug("\nMISS  -  DELETION  (FN)\n");
                        } else {
                            st_logDebug("\nMISS  -  INSERTION  (FN)\n");
                        }

                        st_logDebug("  pos: %" PRIi64 "\n", vcfInfo->referencePos);
                        st_logDebug("  baseline  ref: %s\t alt: ", vcfInfo->refAllele);
                        for (int i = 1; i < vcfInfo->unpackedRecordRef->n_allele; i++) {
                            if (i != 1) st_logDebug(",");
                            st_logDebug("%s", vcfInfo->unpackedRecordRef->d.allele[i]);
                        }
                        st_logDebug(" \tphasing: %d | %d", vcfInfo->refPhasing1, vcfInfo->refPhasing2);
                        if (strlen(vcfInfo->refAllele) > strlen(vcfInfo->altAllele1)) {
                            st_logDebug("\n  output    ref: %s\t alt: %s \tphasing: %d | %d\n",
                                        refSeq, refSeq, 0, 0);
                        } else {
                            st_logDebug("\n  output    ref: %c\t alt: %c \tphasing: %d | %d\n",
                                        refChar, refChar, 0, 0);
                        }


                        printPartitionInfo2(vcfInfo->referencePos, rProbs1, rProbs2);

                        for (int64_t j = 1; j < indelLen; j++) {
                            st_logDebug(" -> pos: %" PRIi64 "\n", vcfInfo->referencePos+j);
                            printPartitionInfo2(vcfInfo->referencePos+j, rProbs1, rProbs2);
                        }
                        free(refSeq);
                    }
                } else {
                    // SNV
                    if (params->verboseFalseNegatives) {
                        st_logDebug("\nMISS  -  SNV  (FN)\n");
                        st_logDebug("  pos: %" PRIi64 "\n", vcfInfo->referencePos);
                        st_logDebug("  baseline  ref: %s\t alt: ", vcfInfo->refAllele);

                        for (int i = 1; i < vcfInfo->unpackedRecordRef->n_allele; i++) {
                            if (i != 1) st_logDebug(",");
                            st_logDebug("%s", vcfInfo->unpackedRecordRef->d.allele[i]);
                        }
                        st_logDebug(" \tphasing: %d | %d", vcfInfo->refPhasing1, vcfInfo->refPhasing2);
                        st_logDebug("\n  output    ref: %s\t alt: %s \tphasing: %d | %d\n",
                                    vcfInfo->refAllele, vcfInfo->refAllele, 0, 0);
                        printPartitionInfo2(vcfInfo->referencePos, rProbs1, rProbs2);
                    }
                }
            } else {
                recordHomozygousVariantMiss2(results, vcfInfo, baseMapper, params, rProbs1, rProbs2);
            }
            free(gt_arr);
        }
    }
    if (results->truePositives == 0) st_logInfo("No matches between vcfs found - did you compare against the correct vcf?\n");

    // Remaining positions after the last variant in the reference are not currently being looked through
    // False positives in this region could therefore be missed
    // (In addition to false positives after the first variant)
    results->negatives += (vcfInfo->referencePos - refStart - results->positives);
    results->trueNegatives += (results->negatives - results->falsePositives);
    results->switchErrorDistance = results->switchErrorDistance/results->switchErrors;

    // cleanup
    vcf_close(inRef);
    vcf_close(inEval);
    bcf_hdr_destroy(hdrRef);
    bcf_hdr_destroy(hdrEval);
    bcf_destroy(refRecord);
    bcf_destroy(evalRecord);
    free(vcfInfo);
    stHash_destruct(referenceNamesToReferencePriors1);
    stHash_destruct(referenceNamesToReferencePriors2);
    stList_destruct(profileSequences1);
    stList_destruct(profileSequences2);
}