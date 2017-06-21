/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com) & Arthur Rand (arand@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <htslib/sam.h>
#include "stRPHmm.h"
#include "jsmn.h"

stBaseMapper* stBaseMapper_construct() {
    stBaseMapper *bm = (stBaseMapper*)st_malloc(sizeof(stBaseMapper));
    bm->charToNum = st_calloc(256, sizeof(int8_t));
    bm->numToChar = st_calloc(ALPHABET_SIZE, sizeof(char));
    bm->wildcard = "";
    bm->size = 0;

    return bm;
}

void stBaseMapper_destruct(stBaseMapper *bm) {
    free(bm->charToNum);
    free(bm->numToChar);
    free(bm);
}

void stBaseMapper_addBases(stBaseMapper *bm, char *bases) {
    for (int i = 0; i < strlen(bases); i++) {
        char base = bases[i];
        if (bm->numToChar[bm->size] == 0) bm->numToChar[bm->size] = base;
        bm->charToNum[base] = bm->size;
    }
    bm->size++;
    if (bm->size > ALPHABET_SIZE) {
        st_errAbort("BaseMapper size has exceeded ALPHABET_SIZE parameter (%d)", ALPHABET_SIZE);
    }
}

void stBaseMapper_setWildcard(stBaseMapper* bm, char *wildcard) {
    bm->wildcard = wildcard;
}

int stBaseMapper_getValueForChar(stBaseMapper *bm, char base) {
    int value = bm->charToNum[base];
    if (value >= 0) return value;
    for (int i = 0; i < strlen(bm->wildcard); i++) {
        if (bm->wildcard[i] == base) {
            return st_randomInt(0, bm->size-1);
        }
    }
    st_errAbort("Base '%c' (%d) not in alphabet", base, base);
    return -1;
}

int stBaseMapper_getCharForValue(stBaseMapper *bm, int value) {
    char base = bm->numToChar[value];
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
     * Set hmm parameters.
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

    // Params object
    stRPHmmParameters *params = st_calloc(1, sizeof(stRPHmmParameters));

    // Variables for hmm parameters (initialize & set default values)
    params->hetSubModel = st_calloc(ALPHABET_SIZE*ALPHABET_SIZE, sizeof(uint16_t));
    params->hetSubModelSlow = st_calloc(ALPHABET_SIZE*ALPHABET_SIZE, sizeof(double));
    params->readErrorSubModel = st_calloc(ALPHABET_SIZE*ALPHABET_SIZE, sizeof(uint16_t));
    params->readErrorSubModelSlow = st_calloc(ALPHABET_SIZE*ALPHABET_SIZE, sizeof(double));
    params->maxNotSumTransitions = true;
    params->maxPartitionsInAColumn = 50;
    params->maxCoverageDepth = MAX_READ_PARTITIONING_DEPTH;
    params->minReadCoverageToSupportPhasingBetweenHeterozygousSites = 0;
    params->offDiagonalReadErrorPseudoCount = 1;
    params->onDiagonalReadErrorPseudoCount = 1;
    params->trainingIterations = 0;
    params->filterBadReads = false;
    params->filterMatchThreshold = 0.90;

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
                setSubstitutionProb(params->hetSubModel, params->hetSubModelSlow, j/ALPHABET_SIZE, j%ALPHABET_SIZE, atof(tokStr));
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
                setSubstitutionProb(params->readErrorSubModel, params->readErrorSubModelSlow, j/ALPHABET_SIZE, j%ALPHABET_SIZE, atof(tokStr));

            }
            i += readErrTok.size + 1;
        }
        if (strcmp(keyString, "maxNotSumTransitions") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            if (strcmp(tokStr, "false") == 0) params->maxNotSumTransitions = false;
            i++;
        }
        if (strcmp(keyString, "maxPartitionsInAColumn") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            params->maxPartitionsInAColumn = atoi(tokStr);
            i++;
        }
        if (strcmp(keyString, "maxCoverageDepth") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            params->maxCoverageDepth = atoi(tokStr);
            i++;
        }
        if (strcmp(keyString, "minReadCoverageToSupportPhasingBetweenHeterozygousSites") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            params->minReadCoverageToSupportPhasingBetweenHeterozygousSites = atoi(tokStr);
            i++;
        }
        if (strcmp(keyString, "onDiagonalReadErrorPseudoCount") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            params->onDiagonalReadErrorPseudoCount = atof(tokStr);
            i++;
        }
        if (strcmp(keyString, "offDiagonalReadErrorPseudoCount") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            params->offDiagonalReadErrorPseudoCount = atof(tokStr);
            i++;
        }
        if (strcmp(keyString, "trainingIterations") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            params->trainingIterations = atoi(tokStr);
            i++;
        }
        if (strcmp(keyString, "filterBadReads") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            if (strcmp(tokStr, "true") == 0) params->filterBadReads = true;
            i++;
        }
        if (strcmp(keyString, "filterMatchThreshold") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            params->filterMatchThreshold = atof(tokStr);
            i++;
        }
    }

    free(js);
    return params;
}

void countIndels(uint32_t *cigar, uint32_t ncigar, int64_t *numInsertions, int64_t *numDeletions) {
    for (uint32_t i = 0; i < ncigar; i++) {
        int cigarOp = cigar[i] & BAM_CIGAR_MASK;
        int cigarNum = cigar[i] >> BAM_CIGAR_SHIFT;

        if (cigarOp == BAM_CINS) *numInsertions += cigarNum;
        if (cigarOp == BAM_CDEL) *numDeletions += cigarNum;
    }
}


/* Parse reads within an input interval of a reference sequence of a bam file
 * and create a list of profile sequences by turning characters into profile probabilities.
 *
 * In future, maybe use mapq scores to adjust profile (or posterior probabilities for
 * signal level alignments).
 * */
int64_t parseReads(stList *profileSequences, char *bamFile, stBaseMapper *baseMapper) {

    samFile *in = hts_open(bamFile, "r");
    if (in == NULL) {
        st_errAbort("ERROR: Cannot open bam file %s\n", bamFile);
        return -1;
    }
    bam_hdr_t *bamHdr = sam_hdr_read(in);
    bam1_t *aln = bam_init1();

    int64_t readCount = 0;

    while(sam_read1(in,bamHdr,aln) > 0){

        int64_t pos = aln->core.pos+1; //left most position of alignment
        char *chr = bamHdr->target_name[aln->core.tid] ; //contig name (chromosome)
        int64_t len = aln->core.l_qseq; //length of the read.
        uint8_t *seq = bam_get_seq(aln);  // DNA sequence
        char *readName = bam_get_qname(aln);
        uint32_t *cigar = bam_get_cigar(aln);

        // If there isn't a cigar string, don't bother including the read, since we don't
        // know how it aligns
        if (aln->core.n_cigar == 0) {
            continue;
        }

        readCount++;
        int64_t start_read = 0;
        int64_t end_read = 0;
        int64_t start_ref = pos;
        int64_t cig_idx = 0;

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
            } else if (cigarOp == BAM_CINS || cigarOp == BAM_CSOFT_CLIP) {
                start_read += cigarNum;
                cig_idx++;
            } else if (cigarOp == BAM_CHARD_CLIP || cigarOp == BAM_CPAD) {
                cig_idx++;
            } else {
                st_errAbort("Unidentifiable cigar operation\n");
            }
        }

        // Check for soft clipping at the end
        int lastCigarOp = cigar[aln->core.n_cigar-1] & BAM_CIGAR_MASK;
        int lastCigarNum = cigar[aln->core.n_cigar-1] >> BAM_CIGAR_SHIFT;
        if (lastCigarOp == BAM_CSOFT_CLIP) {
            end_read += lastCigarNum;
        }

        // Count number of insertions & deletions in sequence
        int64_t numInsertions = 0;
        int64_t numDeletions = 0;
        countIndels(cigar, aln->core.n_cigar, &numInsertions, &numDeletions);
        int64_t trueLength = len-start_read-end_read+numDeletions-numInsertions;

        // Create empty profile sequence
        stProfileSeq *pSeq = stProfileSeq_constructEmptyProfile(chr, readName, pos, trueLength);

        // Variables to keep track of position in sequence / cigar operations
        cig_idx = 0;
        int64_t currPosInOp = 0;
        int64_t cigarOp = -1;
        int64_t cigarNum = -1;
        int64_t idxInSeq = start_read;

        // For each position turn character into profile probability
        // As is, this makes the probability 1 for the base read in, and 0 otherwise
        for (uint32_t i = 0; i < trueLength; i++) {

            if (currPosInOp == 0) {
                cigarOp = cigar[cig_idx] & BAM_CIGAR_MASK;
                cigarNum = cigar[cig_idx] >> BAM_CIGAR_SHIFT;
            }
            if (cigarOp == BAM_CMATCH || cigarOp == BAM_CEQUAL || cigarOp == BAM_CDIFF) {
                int64_t b = stBaseMapper_getValueForChar(baseMapper, seq_nt16_str[bam_seqi(seq, idxInSeq)]);
                pSeq->profileProbs[i * ALPHABET_SIZE + b] = ALPHABET_MAX_PROB;
                idxInSeq++;
            }
            else if (cigarOp == BAM_CDEL || cigarOp == BAM_CREF_SKIP) {
                // Currently, all profile probabilities are set to 0 here
                // This worked better than adding in a '-' character or having a flat distribution of probabilities
            } else if (cigarOp == BAM_CINS) {
                // Currently, ignore insertions
                idxInSeq++;
                i--;
            } else if (cigarOp == BAM_CSOFT_CLIP || cigarOp == BAM_CHARD_CLIP || cigarOp == BAM_CPAD) {
                // nothing really to do here. skip to next cigar operation
                currPosInOp = cigarNum-1;
                i--;
            } else {
                st_logCritical("Unidentifiable cigar operation\n");
            }

            currPosInOp++;
            if (currPosInOp == cigarNum) {
                cig_idx++;
                currPosInOp = 0;
            }
        }
        stList_append(profileSequences, pSeq);
    }


    bam_hdr_destroy(bamHdr);
    bam_destroy1(aln);
    sam_close(in);

    return readCount;
}

