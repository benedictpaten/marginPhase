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

    // BCF header creation bam_plcmd.mpileup
    bcf_hdr_t *bcf_hdr = bcf_hdr_init("w");

    kstring_t str = {0,0,NULL};

    ksprintf(&str, "##vcfCompare=htslib-%s\n",hts_version());
    bcf_hdr_append(bcf_hdr, str.s);

    str.l = 0;
    ksprintf(&str, "##reference=file://%s\n", "TODO");
    bcf_hdr_append(bcf_hdr, str.s);

    str.l = 0;
    ksprintf(&str, "##contig=<ID=%s,length=%d>\n", "chr1", 249250621);
    bcf_hdr_append(bcf_hdr, str.s);

//    bcf_hdr_append(bcf_hdr,"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">");
//    bcf_hdr_append(bcf_hdr,"##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Total allelic depths\">");


    bcf_hdr_add_sample(bcf_hdr, "SMPL1");
    bcf_hdr_add_sample(bcf_hdr, NULL);

    vcfFile *out = vcf_open("out.vcf", "w");
    bcf_hdr_write(out, bcf_hdr);


    bcf1_t *bcf_rec = bcf_init1();
    bcf_clear1(bcf_rec);
    bcf_rec->rid  = 0; //this is the index of the contig specified in the header
    bcf_rec->pos  = 1; //this is the position in the reference
    bcf_rec->qual = 0; //the quality score something something?
    bcf_rec->n_sample = 1;
    int i, j, nals = 1;
    str.l = 0;
    kstring_t tmp = {0,0,NULL};

    // we put comma-separated values into tmp.  first is REF column, rest are in ALT
    kputc('A', &tmp);
    kputc(',', &tmp);
    kputc('C', &tmp);
    kputc(',', &tmp);
    kputc('T', &tmp);
//    kputc("ACGTN"[bc->ori_ref], &bc->tmp);
//    for (i=1; i<5; i++)
//    {
//        if (bc->a[i] < 0) break;
//        kputc(',', &bc->tmp);
//        if ( bc->unseen==i ) kputs("<*>", &bc->tmp);
//        else kputc("ACGT"[bc->a[i]], &bc->tmp);
//        nals++;
//    }
    bcf_update_alleles_str(bcf_hdr, bcf_rec, tmp.s);

    //todo: this doesn't work but I think it's close
//    char *tmp_str[] = {"0|1"};
//    bcf_update_format_string(bcf_hdr, bcf_rec, "GT", (char **)tmp_str, 1);

    bcf_write1(out, bcf_hdr, bcf_rec);




    //writing variants: call2bcf
//    bcf_hdr_t *hdr = bc->bcf_hdr;
//    rec->rid  = bc->tid;
//    rec->pos  = bc->pos;
//    rec->qual = 0;
//
//
//    bc->tmp.l = 0;
//    if (bc->ori_ref < 0)    // indel
//    {
//        // REF
//        kputc(ref[bc->pos], &bc->tmp);
//        for (j = 0; j < bca->indelreg; ++j) kputc(ref[bc->pos+1+j], &bc->tmp);
//
//        // ALT
//        for (i=1; i<4; i++)
//        {
//            if (bc->a[i] < 0) break;
//            kputc(',', &bc->tmp); kputc(ref[bc->pos], &bc->tmp);
//
//            if (bca->indel_types[bc->a[i]] < 0) { // deletion
//                for (j = -bca->indel_types[bc->a[i]]; j < bca->indelreg; ++j)
//                    kputc(ref[bc->pos+1+j], &bc->tmp);
//            } else { // insertion; cannot be a reference unless a bug
//                char *inscns = &bca->inscns[bc->a[i] * bca->maxins];
//                for (j = 0; j < bca->indel_types[bc->a[i]]; ++j)
//                    kputc("ACGTN"[(int)inscns[j]], &bc->tmp);
//                for (j = 0; j < bca->indelreg; ++j) kputc(ref[bc->pos+1+j], &bc->tmp);
//            }
//            nals++;
//        }
//    }
//    else    // SNP
//    {
//        kputc("ACGTN"[bc->ori_ref], &bc->tmp);
//        for (i=1; i<5; i++)
//        {
//            if (bc->a[i] < 0) break;
//            kputc(',', &bc->tmp);
//            if ( bc->unseen==i ) kputs("<*>", &bc->tmp);
//            else kputc("ACGT"[bc->a[i]], &bc->tmp);
//            nals++;
//        }
//    }
//    char *alleles_str = 0;
//    bcf_update_alleles_str(hdr, bcf_rec, alleles_str);




    free(str.s);
    free(tmp.s);
    vcf_close(out);
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


    create_vcf();
    if (true) return 0;


    st_logInfo("VCF Reference: %s \n", vcfReference);
    st_logInfo("VCF Evaluated: %s \n", vcfEvaluated);

    vcfFile *inRef = vcf_open(vcfReference,"r"); //open vcf file
    bcf_hdr_t *hdrRef = bcf_hdr_read(inRef); //read header
    bam1_t *record = bcf_init1(); //initialize for reading


    vcfFile *out = vcf_open("out.vcf", "w");
    vcfFile *outOrig = vcf_open("outOrig.vcf", "w");
//    bcf_hdr_write(out, hdrRef);

    int lineNr = 0;
    while(bcf_read(inRef,hdrRef,record) == 0){

        st_logDebug("%d:\n", lineNr);
        st_logDebug("\ttid:%d pos:%d bin:%d qual:%d l_qname:%d flag:%d unused1:%d l_extranul:%d n_cigar:%d "
                            "l_qseq:%d mtid:%d mpos:%d isize:%d \n", record->core.tid, record->core.pos,
                    record->core.bin, record->core.qual, record->core.l_qname, record->core.flag, record->core.unused1,
                    record->core.l_extranul, record->core.n_cigar, record->core.l_qseq, record->core.mtid,
                    record->core.mpos, record->core.isize);

        bcf1_t *unpackedRecord = record;
        bcf_unpack(unpackedRecord, BCF_UN_ALL);

        st_logDebug("\tid:%s ref:%s alleles:", unpackedRecord->d.id, unpackedRecord->d.als);
        for (int i = 1; i < unpackedRecord->n_allele; i++) {
            if (i!=1) st_logDebug(",");
            st_logDebug(unpackedRecord->d.allele[i]);
        }

        int nsmpl = bcf_hdr_nsamples(hdrRef);
        int32_t *gt_arr = NULL, ngt_arr = 0; //locations for get_genotype types and count

        int ngt = bcf_get_genotypes(hdrRef, record, &gt_arr, &ngt_arr);
        if ( ngt > 0 ) { // GT is present

            st_logDebug(" phasing:");
            int max_ploidy = ngt / nsmpl;
            for (int i = 0; i < nsmpl; i++) {
                int32_t *ptr = gt_arr + i * max_ploidy;
                for (int j = 0; j < max_ploidy; j++) {
                    if (ptr[j] == bcf_int32_vector_end) break;// if true, the sample has smaller ploidy
                    if (bcf_gt_is_missing(ptr[j])) continue;// missing allele
                    int allele_index = bcf_gt_allele(ptr[j]); // the VCF 0-based allele index
                    int is_phased = bcf_gt_is_phased(ptr[j]);

                    if (j!=0) {
                        if (is_phased) st_logDebug("|");
                        else st_logDebug("/");
                    }
                    st_logDebug("%d", allele_index);
                }
            }
            free(gt_arr);
            st_logDebug("\n");
        }

//        //write original record again
//        bcf_write(outOrig, hdrRef, record);
//
//        //modify, write, reset
//        record->id += 1;
//        record->l_data += 1;
//        record->m_data += 1;
//        bcf_write(out, hdrRef, record);
//        record = bcf_init1();

        lineNr++;
        if (lineNr > 8) {
            break;
        }
    }

    vcf_close(inRef);
    vcf_close(out);
    vcf_close(outOrig);



    //while(1); // Use this for testing for memory leaks

    return 0;
}

