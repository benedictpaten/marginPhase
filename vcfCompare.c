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

static void recordHomozygousVariant(stGenotypeResults *results) {
    results->falsePositives++;
    results->error_homozygousInRef++;
    results->positives--;
    results->negatives++;
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

    // Iterate through the vcf being checked until getting to the start of the specified interval
    // Don't bother analyzing these records
    int64_t refStart = 0;
    int64_t evalPos =  0;
    bcf1_t *unpackedRecord;

    // Variables for keeping track of phasing info
    bool phasingHap1 = false;
    bool phasingHap2 = false;
    float switchErrorDistance = 0;

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

        if (strlen(refChar) > 1 || strlen(refAltChar) > 1) {
            results->indels++;
        }

        // Skip to the first known location of variation in file being evaluated
        if (results->positives == 0) {
            refStart = referencePos;
        }

        // Get genotype info
        int i, j, ngt, nsmpl = bcf_hdr_nsamples(hdrRef);
        int32_t *gt_arr = NULL, ngt_arr = 0;
        ngt = bcf_get_genotypes(hdrRef, unpackedRecordRef, &gt_arr, &ngt_arr);

        int allele1 = bcf_gt_allele(gt_arr[0]);
        int allele2 = bcf_gt_allele(gt_arr[1]);
        results->positives++;

        if (maybeFalsePositive && evalPos < referencePos) {
            results->falsePositives++;
            char *evalRefChar = unpackedRecord->d.als;
            char *evalAltChar = unpackedRecord->d.allele[1];
            if (strlen(evalRefChar) > 1 || strlen(evalAltChar) > 1) {
                results->falsePositiveGaps++;
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
                if (strlen(evalRefChar) > 1 || strlen(evalAltChar) > 1){
                    results->falsePositiveGaps++;
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
                } else {
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
                        results->truePositives++;
                        if (strlen(refChar) > 1 || strlen(refAltChar) > 1){
                            results->truePositiveGaps++;
                        }
                    } else if (strcmp(refChar, evalAltChar) == 0
                               && strcmp(evalRefChar, refAltChar) == 0) {
                        results->switchErrors++;
                        results->switchErrorDistance += switchErrorDistance;
                        switchErrorDistance = 0;
                        phasingHap2 = true;
                        phasingHap1 = false;
                        results->truePositives++;
                        if (strlen(refChar) > 1 || strlen(refAltChar) > 1){
                            results->truePositiveGaps++;
                        }
                    } else {
                        results->falsePositives++;
                        results->error_incorrectVariant++;
                    }
                } else {
                    if (strcmp(refChar, evalRefChar) == 0 && strcmp(evalAltChar, refAltChar) == 0) {
                        results->switchErrors++;
                        results->switchErrorDistance += switchErrorDistance;
                        switchErrorDistance = 0;
                        phasingHap2 = true;
                        phasingHap1 = false;
                        results->truePositives++;
                        if (strlen(refChar) > 1 || strlen(refAltChar) > 1){
                            results->truePositiveGaps++;
                        }
                    } else if (strcmp(refChar, evalAltChar) == 0
                               && strcmp(evalRefChar, refAltChar) == 0) {
                        switchErrorDistance++;
                        results->truePositives++;
                        if (strlen(refChar) > 1 || strlen(refAltChar) > 1){
                            results->truePositiveGaps++;
                        }
                    } else {
                        results->falsePositives++;
                        results->error_incorrectVariant++;
                    }
                }
            } else if (phasingHap2) {
                if (allele1 == 0 && allele2 == 1) {
                    if (strcmp(refChar, evalAltChar) == 0 && strcmp(evalRefChar, refAltChar) == 0) {
                        switchErrorDistance++;
                        results->truePositives++;
                        if (strlen(refChar) > 1 || strlen(refAltChar) > 1){
                            results->truePositiveGaps++;
                        }
                    } else if (strcmp(refChar, evalRefChar) == 0
                               && strcmp(evalAltChar, refAltChar) == 0) {
                        results->switchErrors++;
                        results->switchErrorDistance += switchErrorDistance;
                        switchErrorDistance = 0;
                        phasingHap1 = true;
                        phasingHap2 = false;
                        results->truePositives++;
                        if (strlen(refChar) > 1 || strlen(refAltChar) > 1){
                            results->truePositiveGaps++;
                        }
                    } else {
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
                        results->truePositives++;
                        if (strlen(refChar) > 1 || strlen(refAltChar) > 1){
                            results->truePositiveGaps++;
                        }
                    } else if (strcmp(refChar, evalRefChar) == 0 && strcmp(evalAltChar, refAltChar) == 0) {
                        switchErrorDistance++;
                        results->truePositives++;
                        if (strlen(refChar) > 1 || strlen(refAltChar) > 1){
                            results->truePositiveGaps++;
                        }
                    } else {
                        results->falsePositives++;
                        results->error_incorrectVariant++;
                    }
                }
            }

        } else if (evalPos > referencePos){
            // Missed the variant
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
                } else {
                    results->error_badPartition++;
                }
            }
            free(gt_arr);
        }
    }
    if (results->truePositives == 0) {
        st_logInfo("No matches between vcfs found - did you compare against the correct vcf?\n");
    }

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
}


int main(int argc, char *argv[]) {
    // Parameters / arguments

    char *logLevelString = "info";
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
    st_setLogLevelFromString(logLevelString);

    stGenotypeResults *results = st_calloc(1, sizeof(stGenotypeResults));
    compareVCFsBasic(stderr, vcfEvaluated, vcfReference, results);
    printGenotypeResults(results);

    if (true) return 0;
}

