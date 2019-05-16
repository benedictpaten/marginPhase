//
// Created by tpesout on 1/8/19.
//

#ifndef MARGINPHASE_VCF_COMPARISON_H

#include "margin.h"

/*
 * _stGenotypeResults
 * Struct which stores information about relevant test results.
 */
struct _stGenotypeResults {

    // Variants in reference
    int64_t negatives;
    int64_t positives;
    int64_t homozygousVariantsInRef;
    int64_t homozygousVariantsInRef_Insertions;
    int64_t homozygousVariantsInRef_Deletions;
    int64_t hetsInRef;
    int64_t hetsInRef_Insertions;
    int64_t hetsInRef_Deletions;

    // Variants in evaluated vcf
    int64_t truePositives;
    int64_t falsePositives;
    int64_t trueNegatives;
    int64_t falseNegatives;

    // Stats for specific types of variants
    int64_t truePositiveIndels;
    int64_t falsePositiveIndels;
    int64_t truePositiveHomozygous;
    int64_t truePositiveHet;
    int64_t truePositiveHomozygousIndels;
    int64_t truePositiveHetIndels;

    // Types of errors
    int64_t error_missedHet;
    int64_t error_missedHet_Insertions;
    int64_t error_missedHet_Deletions;
    int64_t error_homozygousInRef;
    int64_t error_homozygous_Insertions;
    int64_t error_homozygous_Deletions;

    // Phasing
    int64_t switchErrors;
    float switchErrorDistance;
    int64_t uncertainPhasing;
};
typedef struct _stGenotypeResults stGenotypeResults;
void printGenotypeResults(stGenotypeResults *results);

/*
 * VCF comparison methods
 */

void compareVCFs(FILE *fh, stList *hmms, char *vcf_toEval, char *vcf_ref,
                 stBaseMapper *baseMapper, stGenotypeResults *results, stRPHmmParameters *params);

void compareVCFsBasic(FILE *fh, char *vcf_toEval, char *vcf_ref, stGenotypeResults *results);

void compareVCFs_debugWithBams(char *vcf_toEval, char *vcf_ref, char *bamFile1, char *bamFile2, char *referenceFasta,
                               stBaseMapper *baseMapper, stGenotypeResults *results, stRPHmmParameters *params);


#define MARGINPHASE_VCF_COMPARISON_H
#endif //MARGINPHASE_VCF_COMPARISON_H
