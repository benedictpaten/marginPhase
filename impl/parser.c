//
// Created by Marina Haukness on 4/24/17.
//

#include "stRPHmm.h"
#include "jsmn.h"
#include <htslib/sam.h>
#include <htslib/faidx.h>


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
}

int stBaseMapper_getBaseForValue(stBaseMapper *bm, int value) {
    char base = bm->numToBase[value];
    if (base >= 0) return base;
    st_errAbort("Value '%d' not specified in alphabet", value);
}



char *json_token_tostr(char *js, jsmntok_t *t)
{
    js[t->end] = '\0';
    return js + t->start;
}

/*
     * TODO: Get model parameters from params file.
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
*/
stRPHmmParameters *parseParameters(char *paramsFile, stBaseMapper *baseMapper) {

    // Variables for parsing
    st_logDebug("Parsing json file: %s \n", paramsFile);
    int numTokens = 256;
    jsmn_parser parser;
    jsmn_init(&parser);
    jsmntok_t tokens[numTokens];
    char *js = NULL;
    size_t jslen = 0;
    char buf[BUFSIZ];
    int r;

    // Variables for hmm parameters
    uint16_t  *hetSubModel = st_calloc(ALPHABET_SIZE*ALPHABET_SIZE, sizeof(uint16_t));
    double *hetSubModelSlow = st_calloc(ALPHABET_SIZE*ALPHABET_SIZE, sizeof(double));
    uint16_t  *readErrorSubModel = st_calloc(ALPHABET_SIZE*ALPHABET_SIZE, sizeof(uint16_t));
    double *readErrorSubModelSlow = st_calloc(ALPHABET_SIZE*ALPHABET_SIZE, sizeof(double));
    bool maxNotSumTransitions = true;
    int64_t  maxPartitionsInAColumn = 50;
    int64_t  maxCoverageDepth = MAX_READ_PARTITIONING_DEPTH;
    int64_t minReadCoverageToSupportPhasingBetweenHeterozygousSites = 0;

    // TODO: Error checking
    // Make sure params file exists
    // In case where buffer too small / not enough tokens, reallocate
    FILE *fp;
    fp = fopen(paramsFile, "rb");
    r = fread(buf, 1, sizeof(buf), fp);
    if (r > 0) {
        js = realloc(js, jslen + r + 1);
    }
    fclose(fp);
    if (js != NULL) {
        strncpy(js + jslen, buf, r);
        jslen = jslen + r;
        r = jsmn_parse(&parser, js, jslen, tokens, numTokens);
        st_logDebug("Number of tokens needed to parse: %d\n", r);
    }

    if (r < 0) {
        st_logDebug("Error when parsing json: %d\n", r);
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
//                alphabet[j] = tokStr;
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
                st_errAbort("Size of haplotype substitution model in JSON does not match ALPHABET_SIZE * ALPHABET_SIZE \n");
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
                st_errAbort("Size of read error model in JSON does not match ALPHABET_SIZE * ALPHABET_SIZE \n");
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
    // Construct actual hmm parameters

//    st_logDebug("HAPLOTYPE SUBSTITUTION MODEL: \n");
//    for (int64_t i = 0; i < ALPHABET_SIZE * ALPHABET_SIZE; i++) {
//        st_logDebug(" %f \t", hetSubModelSlow[i]);
//        if ((i+1) % ALPHABET_SIZE == 0) {
//            st_logDebug("\n");
//        }
//    }
//    st_logDebug("READ ERROR MODEL: \n");
//    for (int64_t i = 0; i < ALPHABET_SIZE * ALPHABET_SIZE; i++) {
//        st_logDebug(" %f \t", readErrorSubModelSlow[i]);
//        if ((i+1) % ALPHABET_SIZE == 0) {
//            st_logDebug("\n");
//        }
//    }


    stRPHmmParameters *params = stRPHmmParameters_construct(
            hetSubModel, hetSubModelSlow, readErrorSubModel, readErrorSubModelSlow,
            maxNotSumTransitions, maxPartitionsInAColumn, maxCoverageDepth,
            minReadCoverageToSupportPhasingBetweenHeterozygousSites);
    return params;
}

//char baseToAlphabet(char b, char **alphabet, char *wildcard) {
//    for (size_t i = 0; i < ALPHABET_SIZE; i++) {
//        char *bases = alphabet[i];
//        size_t len = strlen(bases);
//        for (size_t j = 0; j < len; j++) {
//            if (b == bases[j]) return FIRST_ALPHABET_CHAR + i;
//        }
//    }
//    // Wildcard becomes random base (equal probability of any)
//    for (size_t i = 0; i < strlen(wildcard); i++) {
//        if (b == wildcard[i]) return st_randomInt(FIRST_ALPHABET_CHAR, FIRST_ALPHABET_CHAR+ALPHABET_SIZE-1);
//    }
//    st_logInfo("ERROR: Character in sequence not in alphabet\n");
//    return  FIRST_ALPHABET_CHAR - 1;
//}


/*
     * TODO: Use htslib to parse the reads within an input interval of a reference sequence of a bam file
     * and create a list of profile sequences
       where for each position you turn the character into a profile probability, as shown in the tests

       In future we can use the mapq scores for reads to adjust these profiles, or for signal level alignments
       use the posterior probabilities.
     */
void parseReads(stList *profileSequences, char *bamFile, stBaseMapper *baseMapper, char *refSeqName, int32_t intervalStart, int32_t intervalEnd) {

    st_logDebug("Bam file: %s \n", bamFile);
    samFile *in = hts_open(bamFile, "r");
    bam_hdr_t *bamHdr = sam_hdr_read(in);
    bam1_t *aln = bam_init1();

    int32_t readCount = 0;

    // TODO: add implementation to read from specific intervals in bam file (?)
    bool readWholeFile = true;
//    bool readWholeFile = false;
//    if (!refSeqName) {
//        st_logInfo("No reference sequence name given, reading whole file.\n");
//        readWholeFile = true;
//    }
//    if (intervalStart < 0) {
//        intervalStart = 0;
//        st_logInfo("Reading from start of chromosome.\n");
//    }
//    if (intervalEnd < 0) {
//        st_logInfo("Reading until end of chromosome.\n", intervalEnd);
//    }
//
//    hts_idx_t *hts_idx_init(int n, int fmt, uint64_t offset0, int min_shift, int n_lvls);
//    hts_itr_t *sam_itr_queryi(const hts_idx_t *idx, int tid, int beg, int end);
//    hts_itr_t *sam_itr_querys(const hts_idx_t *idx, bam_hdr_t *hdr, const char *region);

    char *referenceName = "hg19.chr3.fa"; //http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr3.fa.gz 
    faidx_t *fai = fai_load(referenceName);
    if ( !fai ) {
        st_logCritical("Could not load fai index of %s.  Maybe you should run 'samtools faidx %s'\n",
                       referenceName, referenceName);
    }
    int seq_len;
    char *ref = fai_fetch(fai, "chr3", &seq_len);
    if ( seq_len < 0 ) {
        st_logCritical("Failed to fetch reference sequence %s in %s\n", "chr3", referenceName);
    }



    while(sam_read1(in,bamHdr,aln) > 0){

        int32_t pos = aln->core.pos +1; //left most position of alignment
        char *chr = bamHdr->target_name[aln->core.tid] ; //contig name (chromosome)
        uint32_t len = aln->core.l_qseq; //length of the read.
        uint8_t *seq = bam_get_seq(aln);  // DNA sequence

        float matched = 0.0;
        float unmatched = 0.0;
        if(readWholeFile || strcmp(chr, refSeqName) == 0) {
            // todo Should reads that cross boundaries be counted?
            if (pos >= intervalStart && (intervalEnd < 0 || pos + len <= intervalEnd)) {
                readCount++;

                stProfileSeq *pSeq = stProfileSeq_constructEmptyProfile(chr, pos, len);

                for (int64_t i = 0; i < len; i++) {
                    // For each position turn character into profile probability
                    // As is, this makes the probability 1 for the base read in, and 0 otherwise
                    // Should this be modified to take into account error rates?
                    // What about coverage from other profile sequences?
                    int b = stBaseMapper_getValueForBase(baseMapper, seq_nt16_str[bam_seqi(seq, i)]);
                    int r = stBaseMapper_getValueForBase(baseMapper, ref[pos + i +1]);
                    if (b == r) matched++;
                    else unmatched++;

                    pSeq->profileProbs[i * ALPHABET_SIZE + b] = ALPHABET_MAX_PROB;
                }
                stList_append(profileSequences, pSeq);
//                if (readCount % 1000 == 0) {
//                    stProfileSeq_print(pSeq, stderr, true);
//                }
            }
        }
        st_logInfo("Read had %2.f%% (out of %5.f) matched chars\n", (100*matched / (matched + unmatched)), (matched + unmatched));
    }
    st_logDebug("Number of profile sequences created: %d\n", readCount);

    bam_destroy1(aln);
    sam_close(in);
}
