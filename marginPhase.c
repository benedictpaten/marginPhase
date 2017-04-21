/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com) & Arthur Rand (arand@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <getopt.h>
#include <stdio.h>
#include <ctype.h>
#include <htslib/sam.h>
#include "stRPHmm.h"
#include "jsmn.h"



void usage() {
    fprintf(stderr, "marginPhase [options] BAM_FILE\n");
    fprintf(stderr,
            "Phases the reads in an interval of a BAM file reporting a gVCF file "
            "giving genotypes and haplotypes for region.\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
    fprintf(stderr, "-b --bamFile : Input BAM file\n");
    fprintf(stderr, "-v --vcfFile : Output VCF file\n");
    fprintf(stderr, "-p --params : Input params file\n");
    fprintf(stderr, "-n --refSeqName : Name of reference sequence\n");
    fprintf(stderr, "-s --intervalStart : Starting position of interval to read\n");
    fprintf(stderr, "-e --intervalEnd : Ending position of interval to read\n");
}

char baseToAlphabet(char b, char **alphabet, char *wildcard) {
    // TODO: use actual alphabet & wildcard that was input
    for (size_t i = 0; i < ALPHABET_SIZE; i++) {
        char *bases = alphabet[i];
        size_t len = strlen(bases);
        for (size_t j = 0; j < len; j++) {
            if (b == bases[j]) return FIRST_ALPHABET_CHAR + i;
        }
    }
    for (size_t i = 0; i < strlen(wildcard); i++) {
        if (b == wildcard[i]) return st_randomInt(FIRST_ALPHABET_CHAR, FIRST_ALPHABET_CHAR+ALPHABET_SIZE-1);
    }
    st_logInfo("ERROR: Character in sequence not in alphabet\n");
    return  FIRST_ALPHABET_CHAR - 1;
}

char *json_token_tostr(char *js, jsmntok_t *t)
{
    js[t->end] = '\0';
    return js + t->start;
}

void parseReads(stList *profileSequences, char *bamFile, char **alphabet, char *wildcard, char *refSeqName, int32_t intervalStart, int32_t intervalEnd) {
    /*
     * TODO: Use htslib to parse the reads within an input interval of a reference sequence of a bam file
     * and create a list of profile sequences using
     * stProfileSeq *stProfileSeq_constructEmptyProfile(char *referenceName,
            int64_t referenceStart, int64_t length);
       where for each position you turn the character into a profile probability, as shown in the tests

       In future we can use the mapq scores for reads to adjust these profiles, or for signal level alignments
       use the posterior probabilities.
     */

    st_logDebug("Bam file: %s \n", bamFile);
    samFile *in = hts_open(bamFile, "r");
    bam_hdr_t *bamHdr = sam_hdr_read(in);
//    BGZF *in = bgzf_open(bamFile,"r"); //open bam file
//    bam_hdr_t *bamHdr = bam_hdr_read(in); //read header
    bam1_t *aln = bam_init1(); //initialize an alignment

    bool readWholeFile = false;
    if (!refSeqName) {
        st_logInfo("No reference sequence name given, reading whole file.\n");
        readWholeFile = true;
    }
    if (intervalStart < 0) {
        intervalStart = 0;
        st_logInfo("Reading from start of chromosome.\n");
    }
    if (intervalEnd < 0) {
        st_logInfo("Reading until end of chromosome.\n", intervalEnd);
    }
    int32_t readCount = 0;

    // Figure out these parameters....
//    int min_shift = 14;
//    int n_lvls = 5;
//    //hts_idx_t *idx = hts_idx_init(bamHdr->n_targets, HTS_FMT_BAI, 0, min_shift, n_lvls);
//    hts_idx_t *idx;
//    char *region;
//    hts_itr_t *itr = sam_itr_querys(idx, bamHdr, region);


//    hts_idx_t *hts_idx_init(int n, int fmt, uint64_t offset0, int min_shift, int n_lvls);
//    hts_itr_t *sam_itr_queryi(const hts_idx_t *idx, int tid, int beg, int end);
//    hts_itr_t *sam_itr_querys(const hts_idx_t *idx, bam_hdr_t *hdr, const char *region);

    while(sam_read1(in,bamHdr,aln) > 0){

        int32_t pos = aln->core.pos +1; //left most position of alignment in zero based coordianate (+1)
        char *chr = bamHdr->target_name[aln->core.tid] ; //contig name (chromosome)
        uint32_t len = aln->core.l_qseq; //length of the read.
        uint8_t *seq = bam_get_seq(aln);  // DNA sequence

        // TODO: is this the right way to deal with the reference sequence / input interval?
        if(readWholeFile || strcmp(chr, refSeqName) == 0) {
            // Should reads that cross boundaries be counted?
            if (pos >= intervalStart && (intervalEnd < 0 || pos + len <= intervalEnd)) {
                readCount++;

                stProfileSeq *pSeq = stProfileSeq_constructEmptyProfile(chr, pos, len);

                for (int64_t i = 0; i < len; i++) {
                    // For each position turn character into profile probability
                    // As is, this makes the probability 1 for the base read by the file, and 0 otherwise
                    // Should this be modified to take into account error rates?
                    char b = baseToAlphabet(seq_nt16_str[bam_seqi(seq, i)], alphabet, wildcard);
                    pSeq->profileProbs[i * ALPHABET_SIZE + b - FIRST_ALPHABET_CHAR] = ALPHABET_MAX_PROB;
                }
                stList_append(profileSequences, pSeq);
                if (readCount % 1000 == 0) {
                    stProfileSeq_print(pSeq, stderr, true);
                }
            }
        }
    }
    st_logDebug("Number of reads found with name: %s in interval %d - %d : %d\n", refSeqName, intervalStart, intervalEnd, readCount);

    bam_destroy1(aln);
//    bgzf_close(in);
    sam_close(in);
}




stRPHmmParameters *parseParameters(char *paramsFile, char **alphabet, char **wildcard) {
    /*
     * TODO: Read model parameters from params file
     * See params.json.
     * Suggest writing basic parser / tests and putting in impl/parser.c (or something similar)
     * Suggest using http://zserge.com/jsmn.html for parsing the json
     *
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

    st_logDebug("Printing tokens... \n");
    for (int64_t i = 0; i < r; i++) {
        jsmntok_t key = tokens[i];
        char *keyString = json_token_tostr(js, &key);

        if (i == 0) {
            st_logDebug("Root Key: %s\n", keyString);
        }
        if (strcmp(keyString, "alphabet") == 0) {
            jsmntok_t alphabetTok = tokens[i+1];
            if (alphabetTok.size != ALPHABET_SIZE) {
                st_errAbort("Alphabet size in JSON does not match constant ALPHABET_SIZE \n");
            }
            for (int j = 0; j < ALPHABET_SIZE; j++) {
                jsmntok_t tok = tokens[i+j+2];
                char *tokStr = json_token_tostr(js, &tok);
                alphabet[j] = tokStr;
            }
            st_logDebug("\n");
            i += ALPHABET_SIZE + 1;

        }
        if (strcmp(keyString, "wildcard") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            *wildcard = tokStr;
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
            st_logDebug("\n");
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
            st_logDebug("\n");
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

    st_logDebug("HAPLOTYPE SUBSTITUTION MODEL: \n");
    for (int64_t i = 0; i < ALPHABET_SIZE * ALPHABET_SIZE; i++) {
        st_logDebug(" %f \t", hetSubModelSlow[i]);
        if ((i+1) % ALPHABET_SIZE == 0) {
            st_logDebug("\n");
        }
    }
    st_logDebug("READ ERROR MODEL: \n");
    for (int64_t i = 0; i < ALPHABET_SIZE * ALPHABET_SIZE; i++) {
        st_logDebug(" %f \t", readErrorSubModelSlow[i]);
        if ((i+1) % ALPHABET_SIZE == 0) {
            st_logDebug("\n");
        }
    }


    stRPHmmParameters *params = stRPHmmParameters_construct(hetSubModel, hetSubModelSlow, readErrorSubModel,
                                                            readErrorSubModelSlow, maxNotSumTransitions,
                                                            maxPartitionsInAColumn, maxCoverageDepth,
                                                            minReadCoverageToSupportPhasingBetweenHeterozygousSites);
    return params;
}




int main(int argc, char *argv[]) {
    // Parameters / arguments
    char * logLevelString = NULL;
    char *bamFile = NULL;
    char *vcfFile = NULL;
    char *paramsFile = "params.json";

    char *refSeqName = NULL;
    int32_t intervalStart = -1;
    int32_t intervalEnd = -1;

    // When done testing, set random seed
    // st_randomSeed(0);

    // Parse the options
    while (1) {
        static struct option long_options[] = {
                { "logLevel", required_argument, 0, 'a' },
                { "help", no_argument, 0, 'h' },
                { "bamFile", required_argument, 0, 'b'},
                { "vcfFile", required_argument, 0, 'v'},
                { "params", required_argument, 0, 'p'},
                { "refSeqName", required_argument, 0, 'n'},
                { "intervalStart", required_argument, 0, 's'},
                { "intervalEnd", required_argument, 0, 'e'},
                { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:v:p:n:s:e:h", long_options, &option_index);

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
        case 'b':
            bamFile = stString_copy(optarg);
            break;
        case 'v':
            vcfFile = stString_copy(optarg);
            break;
        case 'p':
            paramsFile = stString_copy(optarg);
            break;
        case 'n':
            refSeqName = stString_copy(optarg);
            break;
        case 's':
            intervalStart = atoi(optarg);
            break;
        case 'e':
            intervalEnd = atoi(optarg);
            break;
        default:
            usage();
            return 1;
        }
    }

    // Parse any model parameters
    /*
     * TODO: Get model parameters. I suggest we make a simple json or yaml file to hold these parameters.
     * Minimally we need a heterozygozity rate (the fraction of reference positions that are different between the
     * haplotypes (excluding gaps)
     * We also need to figure out what to do with gap positions (which are just treated as an additional character)
     * Once these are read in we need to construct (as shown in the tests) the different matrices.
     * We will need:
     * stRPHmmParameters *stRPHmmParameters_construct
       And:
       void setSubstitutionProb
     */
    st_logInfo("Parsing model parameters\n");

    char *alphabet[ALPHABET_SIZE];
    char *wildcard;
    stRPHmmParameters *params = parseParameters(paramsFile, alphabet, &wildcard);

    st_logDebug("Alphabet array: \n");
    for (int j = 0; j < ALPHABET_SIZE; j++ ) {
        st_logDebug("%d: %s \t%d\n", j, alphabet[j], strlen(alphabet[j]));
    }
    st_logDebug("Wildcard: %s\n", wildcard);

    char c = baseToAlphabet('G', alphabet, wildcard);
    st_logDebug("G was translated to %c\n", c);
    char c2 = baseToAlphabet('-', alphabet, wildcard);
    st_logDebug("- was translated to %c\n", c2);
    char c3 = baseToAlphabet('n', alphabet, wildcard);
    st_logDebug("n was translated to %c\n", c3);

    st_logDebug("Random chars... \n");
    for (int i = 0;i < 10; i++) {
        st_logDebug("%c \t", baseToAlphabet('N', alphabet, wildcard));
        st_logDebug("%c \t", baseToAlphabet('n', alphabet, wildcard));
    }


    // Parse reads for interval
    st_logInfo("Parsing input reads\n");

    stList *profileSequences = stList_construct();
    parseReads(profileSequences, bamFile, alphabet, wildcard, refSeqName, intervalStart, intervalEnd);

    

    // Create HMMs
    st_logInfo("Creating read partitioning HMMs\n");
    stList *hmms = getRPHmms(profileSequences, params);

    // Break up the hmms where the phasing is uncertain
    st_logInfo("Breaking apart HMMs where the phasing is uncertain\n");

    stList *l = stList_construct3(0, (void (*)(void *))stRPHmm_destruct2);
    while(stList_length(hmms) > 0) {
        stList_appendAll(l, stRPHMM_splitWherePhasingIsUncertain(stList_pop(hmms)));
    }
    hmms = l;

    // Create HMM traceback and genome fragments

    // For each read partitioning HMM
    for(int64_t i=0; i<stList_length(hmms); i++) {
        stRPHmm *hmm = stList_get(hmms, i);

        st_logInfo("Creating genome fragment for reference sequence: %s, start: %" PRIi64 ", length: %" PRIi64 "\n",
                    hmm->referenceName, hmm->refStart, hmm->refLength);

        // Run the forward-backward algorithm
        stRPHmm_forwardBackward(hmm);

        // Now compute a high probability path through the hmm
        stList *path = stRPHmm_forwardTraceBack(hmm);

        // Compute the genome fragment
        stGenomeFragment *gF = stGenomeFragment_construct(hmm, path);

        // Write out VCF
        st_logInfo("Writing out VCF for fragment\n");
        /*
         * TODO: Convert the genome fragment into a portion of a VCF file (we'll need to write the header out earlier)
         * We can express the genotypes and (using phase sets) the phasing relationships.
         */

        // Optionally write out two BAMs, one for each read partition
        st_logInfo("Writing out BAM partitions for fragment\n");
        /*
         * TODO: Optionally, write out a new BAM file expressing the partition (which reads belong in which partition)
         * Not sure if we need to write out multiple files or if we can add a per read flag to express this information.
         */
    }

    // Cleanup
    stList_destruct(hmms);

    //while(1); // Use this for testing for memory leaks

    return 0;
}

