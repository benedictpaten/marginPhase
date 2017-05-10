/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com) & Arthur Rand (arand@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <htslib/sam.h>
#include "stRPHmm.h"
#include "jsmn.h"

stBaseMapper* stBaseMapper_construct() {
    stBaseMapper *bm = (stBaseMapper*)calloc(1, sizeof(stBaseMapper));
    bm->baseToNum = calloc(256, sizeof(char));
    bm->numToBase = calloc(256, sizeof(int));
    bm->wildcard = "";
    bm->size = 0;
    for (int i = 0; i < 256; i++) {
        bm->baseToNum[i] = -1;
        bm->numToBase[i] = -1;
    }
    return bm;
}

void stBaseMapper_addBases(stBaseMapper *bm, char *bases) {
    for (int i = 0; i < strlen(bases); i++) {
        char base = bases[i];
        if (bm->numToBase[bm->size] < 0) bm->numToBase[bm->size] = base;
        bm->baseToNum[base] = bm->size;
    }
    bm->size++;
    if (bm->size > ALPHABET_SIZE) {
        st_errAbort("BaseMapper size has exceeded ALPHABET_SIZE parameter (%d)", ALPHABET_SIZE);
    }
}

void stBaseMapper_setWildcard(stBaseMapper* bm, char *wildcard) {
    bm->wildcard = wildcard;
}

int stBaseMapper_getValueForBase(stBaseMapper *bm, char base) {
    int value = bm->baseToNum[base];
    if (value >= 0) return value;
    for (int i = 0; i < strlen(bm->wildcard); i++) {
        if (bm->wildcard[i] == base) {
            return st_randomInt(0, bm->size-1);
        }
    }
    st_errAbort("Base '%c' (%d) not in alphabet", base, base);
    return -1;
}

int stBaseMapper_getBaseForValue(stBaseMapper *bm, int value) {
    char base = bm->numToBase[value];
    if (base >= 0) return base;
    st_errAbort("Value '%d' not specified in alphabet", value);
    return -1;
}



char *json_token_tostr(char *js, jsmntok_t *t)
{
    js[t->end] = '\0';
    return js + t->start;
}

/*
     * Get model parameters from params file.
     *
     * Minimally we need a heterozygozity rate (the fraction of reference positions that are different between the
     * haplotypes (excluding gaps)
     * We also need to figure out what to do with gap positions (which are just treated as an additional character)
     * Once these are read in we need to construct (as shown in the tests) the different matrices.

     * Notes:
     *      "alphabet" is an array specifying the conversion of symbols from the read alignments into the non-negative integer space
     *      used by the program.  e.g. "alphabet" : [ "Aa", "Cc", "Gg", "Tt", "-" ] specifies an alphabet of cardinality 5,
     *      with each string in the array specifying which characters map to which integer, starting from 0. In the example,
     *      "C" or a "c" character is mapped to 1 while "-" is mapped to 4.
     *
     *      The wildcard symbols are treated as mapping to each possible integer symbol with equal probability
     *      Any other symbol encountered by the parsing of reads should be treated as an error
     *
     *      If ALPHABET_SIZE does not equal the cardinality of the input alphabet then an error should be thrown.
     *
     *      The "haplotypeSubstitutionModel" gives probabilities of substitutions between haplotypes in the model, the "readErrorModel"
     *      gives the probabilities of errors in the reads.
     *      Each is a square matrix of size alphabet size * alphabet size
     *      Each should be
     *      converted to log space for the model. Each row of each matrix should sum to 1.0 (roughly) and be normalised to 1.0
     *
     *
*/
stRPHmmParameters *parseParameters(char *paramsFile, stBaseMapper *baseMapper) {

    // Variables for parsing
    int numTokens = 256;
    jsmn_parser parser;
    jsmn_init(&parser);
    jsmntok_t tokens[numTokens];
    char *js = NULL;
    size_t jslen = 0;
    char buf[BUFSIZ];
    int r;

    // Variables for hmm parameters (initialize & set default values)
    uint16_t  *hetSubModel = st_calloc(ALPHABET_SIZE*ALPHABET_SIZE, sizeof(uint16_t));
    double *hetSubModelSlow = st_calloc(ALPHABET_SIZE*ALPHABET_SIZE, sizeof(double));
    uint16_t  *readErrorSubModel = st_calloc(ALPHABET_SIZE*ALPHABET_SIZE, sizeof(uint16_t));
    double *readErrorSubModelSlow = st_calloc(ALPHABET_SIZE*ALPHABET_SIZE, sizeof(double));
    bool maxNotSumTransitions = true;
    int64_t  maxPartitionsInAColumn = 50;
    int64_t  maxCoverageDepth = MAX_READ_PARTITIONING_DEPTH;
    int64_t minReadCoverageToSupportPhasingBetweenHeterozygousSites = 0;

    FILE *fp;
    fp = fopen(paramsFile, "rb");
    if (fp == NULL) {
        st_errAbort("ERROR: Cannot open parameters file %s\n", paramsFile);
    }
    r = fread(buf, 1, sizeof(buf), fp);
    if (r > 0) {
        js = realloc(js, jslen + r + 1);
    }
    fclose(fp);
    if (js != NULL) {
        strncpy(js + jslen, buf, r);
        jslen = jslen + r;
        r = jsmn_parse(&parser, js, jslen, tokens, numTokens);
    }

    if (r == JSMN_ERROR_NOMEM) {
        st_logDebug("Error when parsing json: not enough tokens allocated. Is the JSON file too big? %d\n", r);
        return NULL;
    }

    // Parse tokens, starting at token 1
    // (token 0 is entire object)
    for (int64_t i = 1; i < r; i++) {
        jsmntok_t key = tokens[i];
        char *keyString = json_token_tostr(js, &key);

        if (strcmp(keyString, "alphabet") == 0) {
            jsmntok_t alphabetTok = tokens[i+1];
            if (alphabetTok.size != ALPHABET_SIZE) {
                st_errAbort("Alphabet size in JSON does not match constant ALPHABET_SIZE \n");
            }
            for (int j = 0; j < ALPHABET_SIZE; j++) {
                jsmntok_t tok = tokens[i+j+2];
                char *tokStr = json_token_tostr(js, &tok);
                stBaseMapper_addBases(baseMapper, tokStr);
            }
            i += ALPHABET_SIZE + 1;
        }
        if (strcmp(keyString, "wildcard") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            stBaseMapper_setWildcard(baseMapper, tokStr);
            i++;
        }
        if (strcmp(keyString, "haplotypeSubstitutionModel") == 0) {
            jsmntok_t hapSubTok = tokens[i+1];
            if (hapSubTok.size != ALPHABET_SIZE * ALPHABET_SIZE) {
                st_errAbort("ERROR: Size of haplotype substitution model in JSON does not match ALPHABET_SIZE * ALPHABET_SIZE \n");
            }
            for (int j = 0; j < ALPHABET_SIZE * ALPHABET_SIZE; j++) {
                jsmntok_t tok = tokens[i+j+2];
                char *tokStr = json_token_tostr(js, &tok);
                setSubstitutionProb(hetSubModel, hetSubModelSlow, j/ALPHABET_SIZE, j%ALPHABET_SIZE, atof(tokStr));
            }
            i += hapSubTok.size + 1;
        }
        if (strcmp(keyString, "readErrorModel") == 0) {
            jsmntok_t readErrTok = tokens[i+1];
            if (readErrTok.size != ALPHABET_SIZE * ALPHABET_SIZE) {
                st_errAbort("ERROR: Size of read error model in JSON does not match ALPHABET_SIZE * ALPHABET_SIZE \n");
            }
            for (int j = 0; j < ALPHABET_SIZE * ALPHABET_SIZE; j++) {
                jsmntok_t tok = tokens[i+j+2];
                char *tokStr = json_token_tostr(js, &tok);
                setSubstitutionProb(readErrorSubModel, readErrorSubModelSlow, j/ALPHABET_SIZE, j%ALPHABET_SIZE, atof(tokStr));

            }
            i += readErrTok.size + 1;
        }
        if (strcmp(keyString, "maxNotSumTransitions") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            if (strcmp(tokStr, "false") == 0) maxNotSumTransitions = false;
            i++;
        }
        if (strcmp(keyString, "maximumPartitionsInAColumn") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            maxPartitionsInAColumn = atoi(tokStr);
            i++;
        }
        if (strcmp(keyString, "maxCoverageDepth") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            maxCoverageDepth = atoi(tokStr);
            i++;
        }
        if (strcmp(keyString, "minReadCoverageToSupportPhasingBetweenHeterozygousSites") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            minReadCoverageToSupportPhasingBetweenHeterozygousSites = atoi(tokStr);
            i++;
        }
    }

    stRPHmmParameters *params = stRPHmmParameters_construct(
            hetSubModel, hetSubModelSlow, readErrorSubModel, readErrorSubModelSlow,
            maxNotSumTransitions, maxPartitionsInAColumn, maxCoverageDepth,
            minReadCoverageToSupportPhasingBetweenHeterozygousSites);
    return params;
}


/* Parse reads within an input interval of a reference sequence of a bam file
 * and create a list of profile sequences by turning characters into profile probabilities.
 *
 * In future, maybe use mapq scores to adjust profile (or posterior probabilities for
 * signal level alignments).
 * */
void parseReads(stList *profileSequences, char *bamFile, stBaseMapper *baseMapper) {

    samFile *in = hts_open(bamFile, "r");
    if (in == NULL) {
        st_errAbort("ERROR: Cannot open bam file %s\n", bamFile);
    }
    bam_hdr_t *bamHdr = sam_hdr_read(in);
    bam1_t *aln = bam_init1();

    int32_t readCount = 0;
    int32_t numInsertions = 0;
    int32_t numDeletions = 0;

    while(sam_read1(in,bamHdr,aln) > 0){

        int32_t pos = aln->core.pos; //left most position of alignment
        char *chr = bamHdr->target_name[aln->core.tid] ; //contig name (chromosome)
        uint32_t len = aln->core.l_qseq; //length of the read.
        uint8_t *seq = bam_get_seq(aln);  // DNA sequence
        uint32_t *cigar = bam_get_cigar(aln);

        readCount++;
        uint32_t start_read = 0;
        uint32_t end_read = 0;
        uint32_t start_ref = pos;
        uint32_t cig_idx = 0;

        // Find the correct starting locations on the read and reference sequence,
        // to deal with things like inserts / deletions / soft clipping
        while(cig_idx < aln->core.n_cigar) {
            int cigarOp = cigar[cig_idx] & BAM_CIGAR_MASK;
            int cigarNum = cigar[cig_idx] >> BAM_CIGAR_SHIFT;

            if (cigarOp == BAM_CMATCH || cigarOp == BAM_CEQUAL || cigarOp==BAM_CDIFF) {
                break;
            }
            else if (cigarOp == BAM_CDEL || cigarOp == BAM_CREF_SKIP) {
                start_ref += cigarNum;
                cig_idx++;
                numDeletions++;
            } else if (cigarOp == BAM_CINS || cigarOp == BAM_CSOFT_CLIP) {
                start_read += cigarNum;
                cig_idx++;
                if (cigarOp == BAM_CINS) numInsertions++;
            } else if (cigarOp == BAM_CHARD_CLIP || cigarOp == BAM_CPAD) {
                cig_idx++;
            } else {
                st_logCritical("Unidentifiable cigar operation\n");
            }
        }
        int lastCigarOp = cigar[aln->core.n_cigar-1] & BAM_CIGAR_MASK;
        int lastCigarNum = cigar[aln->core.n_cigar-1] >> BAM_CIGAR_SHIFT;
        if (lastCigarOp == BAM_CSOFT_CLIP) {
            end_read += lastCigarNum;
        }

        char *readName = bam_get_qname(aln);
        stProfileSeq *pSeq = stProfileSeq_constructEmptyProfile(chr, readName, pos, len-start_read-end_read);

        for (uint32_t i = 0; i < len-start_read-end_read; i++) {
            // For each position turn character into profile probability
            // As is, this makes the probability 1 for the base read in, and 0 otherwise
            int b = stBaseMapper_getValueForBase(baseMapper, seq_nt16_str[bam_seqi(seq, start_read+i)]);

            pSeq->profileProbs[i * ALPHABET_SIZE + b] = ALPHABET_MAX_PROB;
        }
        stList_append(profileSequences, pSeq);
    }
    st_logDebug("\tCreated %d profile sequences\n", readCount);

    bam_destroy1(aln);
    sam_close(in);
}

