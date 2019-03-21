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
#include "vcfComparison.h"
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



int main(int argc, char *argv[]) {
    /*
     * Run the basic vf comparison program.
     */
    // Parameters / arguments
    char *logLevelString = "info";
    char *vcfReference = NULL;
    char *vcfEvaluated = NULL;
    bool debug = false;
    char *bamFile1 = NULL;
    char *bamFile2 = NULL;
    char *paramsFile = NULL;
    char *referenceFasta = NULL;
    char *inputBase = NULL;

    // Parse the options
    while (1) {
        static struct option long_options[] = {
                { "logLevel", required_argument, 0, 'a' },
                { "help", no_argument, 0, 'h' },
                { "vcfReference", required_argument, 0, 'r'},
                { "vcfEvaluated", required_argument, 0, 'e'},
                { "debug", no_argument, 0, 'd'},
                { "bamFile1", required_argument, 0, '1'},
                { "bamFile2", required_argument, 0, '2'},
                { "params", required_argument, 0, 'p'},
                { "fasta", required_argument, 0, 'f'},
                { "inputBase", required_argument, 0, 'i'},
                { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:r:e:p:1:2:f:i:hd", long_options, &option_index);

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
            case 'd':
                debug = true;
                logLevelString = "debug";
                break;
            case '1':
                bamFile1 = stString_copy(optarg);
                break;
            case '2':
                bamFile2 = stString_copy(optarg);
                break;
            case 'p':
                paramsFile = stString_copy(optarg);
                break;
            case 'i':
                inputBase = stString_copy(optarg);
                break;
            case 'f':
                referenceFasta = stString_copy(optarg);
                break;
            default:
                usage();
                return 1;
        }
    }
    st_setLogLevelFromString(logLevelString);
    if (inputBase != NULL) {
        if (bamFile1 == NULL) bamFile1 = stString_print("%s.1.bam", inputBase);
        if (bamFile2 == NULL) bamFile2 = stString_print("%s.2.bam", inputBase);
        if (vcfEvaluated == NULL) vcfEvaluated = stString_print("%s.vcf", inputBase);
    }

    stGenotypeResults *results = st_calloc(1, sizeof(stGenotypeResults));

    if (vcfReference == NULL) {
        st_errAbort("ERROR: must specify reference vcf file. Use flag -r <FILE>\n");
    }
    if (vcfEvaluated == NULL) {
        st_errAbort("ERROR: must specify vcf file with calls to evaluate. Use flag -e <FILE>\n");
    }

    if (debug) {
        // Error checking
        if (paramsFile == NULL) {
            st_errAbort("ERROR: must specify parameters file when in debug mode. Use flag -p <FILE>\n");
        }
        if (bamFile1 == NULL || bamFile2 == NULL) {
            st_errAbort("ERROR: must specify bam files for each partition when in debug mode. Use flags -1 <FILE> -2 <FILE>\n");
        }
        if (referenceFasta == NULL) {
            st_errAbort("ERROR: must specify reference fasta while when in debug mode. Use flag -f <FILE>\n");
        }
        // Parse files
        st_logInfo("> Parsing input files. params file: %s  bam1file: %s  bam2file: %s\n",
                   paramsFile, bamFile1, bamFile2);
        FILE *fh = fopen(paramsFile, "rb");
        if (fh == NULL) {
            st_errAbort("ERROR: Cannot open parameters file %s\n", paramsFile);
        }
        Params *fullParams = params_readParams2(fh, FALSE, TRUE);
        fclose(fh);
        stBaseMapper *baseMapper = fullParams->baseMapper;
        stRPHmmParameters *params = fullParams->phaseParams;

        // Compare vcfs
        compareVCFs_debugWithBams(vcfEvaluated, vcfReference, bamFile1, bamFile2,
                                  referenceFasta, baseMapper, results, params);
        params_destruct(fullParams);
    } else {
        compareVCFsBasic(stderr, vcfEvaluated, vcfReference, results);
    }
    printGenotypeResults(results);

    free(results);

    return 0;
}
