/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com) & Arthur Rand (arand@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <htslib/hts.h>
#include "stRPHmm.h"


void writeSplitSams(char *bamInFile, char *bamOutBase,
                          stSet *haplotype1Ids, stSet *haplotype2Ids) {
    // prep
    char haplotype1BamOutFile[strlen(bamOutBase) + 7];
    strcpy(haplotype1BamOutFile, bamOutBase);
    strcat(haplotype1BamOutFile, ".1.bam");
    char haplotype2BamOutFile[strlen(bamOutBase) + 7];
    strcpy(haplotype2BamOutFile, bamOutBase);
    strcat(haplotype2BamOutFile, ".2.bam");
    char unmatchedBamOutFile[strlen(bamOutBase) + 15];
    strcpy(unmatchedBamOutFile, bamOutBase);
    strcat(unmatchedBamOutFile, ".unmatched.bam");

    // file management
    samFile *in = hts_open(bamInFile, "r");
    if (in == NULL) {
        st_errAbort("ERROR: Cannot open bam file %s\n", bamInFile);
    }
    bam_hdr_t *bamHdr = sam_hdr_read(in);
    bam1_t *aln = bam_init1();

    int r;
    st_logDebug("Writing haplotype output to: %s and %s \n", haplotype1BamOutFile, haplotype2BamOutFile);
    samFile *out1 = hts_open_format(haplotype1BamOutFile, "w", &(in->format));
    r = sam_hdr_write(out1, bamHdr);
    if (r != 0) st_logDebug("Couldn't write header to %s\n", haplotype1BamOutFile);

    samFile *out2 = hts_open_format(haplotype2BamOutFile, "w", &(in->format));
    r = sam_hdr_write(out2, bamHdr);
    if (r != 0) st_logDebug("Couldn't write header to %s\n", haplotype2BamOutFile);

    samFile *out3 = hts_open_format(unmatchedBamOutFile, "w", &(in->format));
    r = sam_hdr_write(out3, bamHdr);
    if (r != 0) st_logDebug("Couldn't write header to %s\n", unmatchedBamOutFile);

    // read in input file, write out each read to one sam file
    int64_t readCountH1 = 0;
    int64_t readCountH2 = 0;
    int64_t readCountNeither = 0;
    while(bam_read1(in->fp.bgzf,aln) > 0) {

        char *readName = bam_get_qname(aln);
        if (stSet_search(haplotype1Ids, readName) != NULL) {
            r = bam_write1(out1->fp.bgzf, aln);
            readCountH1++;
        } else if (stSet_search(haplotype2Ids, readName) != NULL) {
            r = bam_write1(out2->fp.bgzf, aln);
            readCountH2++;
        } else {
            r = bam_write1(out3->fp.bgzf, aln);
            readCountNeither++;
        }
    }
    st_logDebug("Read counts:\n\thap1:%" PRIi64 "\thap2:%" PRIi64 "\tneither:%" PRIi64 "\n", readCountH1, readCountH2, readCountNeither);

    sam_close(in);
    sam_close(out1);
    sam_close(out2);
    sam_close(out3);
    bam_destroy1(aln);
    bam_hdr_destroy(bamHdr);
}

