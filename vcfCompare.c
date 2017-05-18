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

