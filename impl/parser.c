/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com) & Arthur Rand (arand@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <unistd.h>
#include <htslib/sam.h>
#include <util.h>
#include "stRPHmm.h"
#include "jsmn.h"

/*
 * stBaseMapper constructor
 */
stBaseMapper* stBaseMapper_construct() {
    stBaseMapper *bm = (stBaseMapper*)st_malloc(sizeof(stBaseMapper));
    bm->charToNum = st_calloc(256, sizeof(uint8_t));
    bm->numToChar = st_calloc(ALPHABET_SIZE, sizeof(uint8_t));
    bm->wildcard = "";
    bm->size = 0;
    return bm;
}

/*
 * stBaseMapper destructor
 */
void stBaseMapper_destruct(stBaseMapper *bm) {
    free(bm->charToNum);
    free(bm->numToChar);
    free(bm);
}

/*
 * Add bases into the baseMapper object.
 */
void stBaseMapper_addBases(stBaseMapper *bm, char *bases) {
    for (uint8_t i = 0; i < strlen(bases); i++) {
        char base = bases[i];
        if (bm->numToChar[bm->size] == 0) bm->numToChar[bm->size] = base;
        bm->charToNum[base] = bm->size;
    }
    bm->size++;
    if (bm->size > ALPHABET_SIZE) {
        st_errAbort("BaseMapper size has exceeded ALPHABET_SIZE parameter (%d)", ALPHABET_SIZE);
    }
}

/*
 * Set the baseMapper wildcard.
 */
void stBaseMapper_setWildcard(stBaseMapper* bm, char *wildcard) {
    bm->wildcard = wildcard;
}

/*
 * Given a character for a base, return the numeric value.
 */
uint8_t stBaseMapper_getValueForChar(stBaseMapper *bm, char base) {
    uint8_t value = bm->charToNum[base];
    if (value >= 0) return value;
    for (int i = 0; i < strlen(bm->wildcard); i++) {
        if (bm->wildcard[i] == base) {
            assert(bm->size-1 < UINT8_MAX);
            return st_randomInt(0, bm->size-1);
        }
    }
    st_errAbort("Base '%c' (%d) not in alphabet", base, base);
    return UINT8_MAX;
}

/*
 * Given the numeric value for a base, return the char.
 */
char stBaseMapper_getCharForValue(stBaseMapper *bm, int value) {
    char base = bm->numToChar[value];
    if (base >= 0) return base;
    st_errAbort("Value '%d' not specified in alphabet", value);
    return -1;
}

/*
 * Convert json token into a string.
 */
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
    char buf[BUFSIZ * 10]; // TODO: FIX, This is terrible code, we should not assume the size of the file is less than this
    int r;

    // Params object
    stRPHmmParameters *params = st_calloc(1, sizeof(stRPHmmParameters));

    // Variables for hmm parameters (initialize & set default values)
    params->hetSubModel = st_calloc(ALPHABET_SIZE*ALPHABET_SIZE, sizeof(uint16_t));
    params->hetSubModelSlow = st_calloc(ALPHABET_SIZE*ALPHABET_SIZE, sizeof(double));
    params->readErrorSubModel = st_calloc(ALPHABET_SIZE*ALPHABET_SIZE, sizeof(uint16_t));
    params->readErrorSubModelSlow = st_calloc(ALPHABET_SIZE*ALPHABET_SIZE, sizeof(double));

    // More variables for hmm stuff
    params->maxCoverageDepth = MAX_READ_PARTITIONING_DEPTH;
    params->maxNotSumTransitions = true;
    params->minPartitionsInAColumn = 50;
    params->maxPartitionsInAColumn = 200;
    params->minPosteriorProbabilityForPartition = 0.001;
    params->minReadCoverageToSupportPhasingBetweenHeterozygousSites = 0;

    // Hmm training options
    params->trainingIterations = 0;
    params->offDiagonalReadErrorPseudoCount = 1;
    params->onDiagonalReadErrorPseudoCount = 1;

    // Read filtering options
    params->minSecondMostFrequentBaseFilter = 2;
    params->minSecondMostFrequentBaseLogProbFilter = 0;
    params->filterAReadWithAnyOneOfTheseSamFlagsSet = 0;
    params->filterBadReads = false;
    params->filterMatchThreshold = 0.90;
    params->filterLikelyHomozygousSites = false;
    params->mapqFilter = 0;

    // Other marginPhase program options
    params->useReferencePrior = false;
    params->gapCharactersForDeletions = true;
    params->estimateReadErrorProbsEmpirically = false;
    params->roundsOfIterativeRefinement = 0;
    params->includeInvertedPartitions = true;
    params->writeGVCF = false;
    params->writeSplitSams = true;
    params->writeUnifiedSam = true;

    setVerbosity(params, 0);

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
        st_errAbort("Error when parsing json: not enough tokens allocated. Is the JSON file too big? %d\n", r);
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
        else if (strcmp(keyString, "wildcard") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            stBaseMapper_setWildcard(baseMapper, tokStr);
            i++;
        }
        else if (strcmp(keyString, "haplotypeSubstitutionModel") == 0) {
            jsmntok_t hapSubTok = tokens[i+1];
            if (hapSubTok.size != ALPHABET_SIZE * ALPHABET_SIZE) {
                st_errAbort("ERROR: Size of haplotype substitution model in JSON "
                                    "does not match ALPHABET_SIZE * ALPHABET_SIZE \n");
            }
            for (int j = 0; j < ALPHABET_SIZE * ALPHABET_SIZE; j++) {
                jsmntok_t tok = tokens[i+j+2];
                char *tokStr = json_token_tostr(js, &tok);
                setSubstitutionProb(params->hetSubModel, params->hetSubModelSlow,
                                    j/ALPHABET_SIZE, j%ALPHABET_SIZE, atof(tokStr));
            }
            i += hapSubTok.size + 1;
        }
        else if (strcmp(keyString, "readErrorModel") == 0) {
            jsmntok_t readErrTok = tokens[i+1];
            if (readErrTok.size != ALPHABET_SIZE * ALPHABET_SIZE) {
                st_errAbort("ERROR: Size of read error model in JSON "
                                    "does not match ALPHABET_SIZE * ALPHABET_SIZE \n");
            }
            for (int j = 0; j < ALPHABET_SIZE * ALPHABET_SIZE; j++) {
                jsmntok_t tok = tokens[i+j+2];
                char *tokStr = json_token_tostr(js, &tok);
                setSubstitutionProb(params->readErrorSubModel, params->readErrorSubModelSlow,
                                    j/ALPHABET_SIZE, j%ALPHABET_SIZE, atof(tokStr));

            }
            i += readErrTok.size + 1;
        }
        else if (strcmp(keyString, "maxNotSumTransitions") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            assert(strcmp(tokStr, "true") || strcmp(tokStr, "false"));
            params->maxNotSumTransitions = strcmp(tokStr, "true") == 0;
            i++;
        }
        else if (strcmp(keyString, "minPartitionsInAColumn") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            params->minPartitionsInAColumn = atoi(tokStr);
            i++;
        }
        else if (strcmp(keyString, "maxPartitionsInAColumn") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            params->maxPartitionsInAColumn = atoi(tokStr);
            i++;
        }
        else if (strcmp(keyString, "minPosteriorProbabilityForPartition") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            params->minPosteriorProbabilityForPartition = atof(tokStr);
            i++;
        }
        else if (strcmp(keyString, "maxCoverageDepth") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            params->maxCoverageDepth = atoi(tokStr);
            i++;
        }
        else if (strcmp(keyString, "minReadCoverageToSupportPhasingBetweenHeterozygousSites") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            params->minReadCoverageToSupportPhasingBetweenHeterozygousSites = atoi(tokStr);
            i++;
        }
        else if (strcmp(keyString, "onDiagonalReadErrorPseudoCount") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            params->onDiagonalReadErrorPseudoCount = atof(tokStr);
            i++;
        }
        else if (strcmp(keyString, "offDiagonalReadErrorPseudoCount") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            params->offDiagonalReadErrorPseudoCount = atof(tokStr);
            i++;
        }
        else if (strcmp(keyString, "trainingIterations") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            params->trainingIterations = atoi(tokStr);
            i++;
        }
        else if (strcmp(keyString, "filterBadReads") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            if (strcmp(tokStr, "true") == 0) params->filterBadReads = true;
            i++;
        }
        else if (strcmp(keyString, "filterMatchThreshold") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            params->filterMatchThreshold = atof(tokStr);
            i++;
        }
        else if (strcmp(keyString, "useReferencePrior") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            assert(strcmp(tokStr, "true") || strcmp(tokStr, "false"));
            params->useReferencePrior = strcmp(tokStr, "true") == 0;
            i++;
        }
        else if (strcmp(keyString, "includeInvertedPartitions") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            assert(strcmp(tokStr, "true") || strcmp(tokStr, "false"));
            params->includeInvertedPartitions = strcmp(tokStr, "true") == 0;
            i++;
        }
        else if (strcmp(keyString, "verbose") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            int64_t bitString = atoi(tokStr);
            setVerbosity(params, bitString);
            i++;
        }
        else if(strcmp(keyString, "filterLikelyHomozygousSites") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            assert(strcmp(tokStr, "true") || strcmp(tokStr, "false"));
            params->filterLikelyHomozygousSites = strcmp(tokStr, "true") == 0;
            i++;
        }
        else if (strcmp(keyString, "minSecondMostFrequentBaseFilter") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            params->minSecondMostFrequentBaseFilter = atof(tokStr);
            i++;
        }
        else if (strcmp(keyString, "minSecondMostFrequentBaseLogProbFilter") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            params->minSecondMostFrequentBaseLogProbFilter = atof(tokStr);
            i++;
        }
        else if (strcmp(keyString, "gapCharactersForDeletions") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            assert(strcmp(tokStr, "true") || strcmp(tokStr, "false"));
            params->gapCharactersForDeletions = strcmp(tokStr, "true") == 0;
            i++;
        }
        else if (strcmp(keyString, "filterAReadWithAnyOneOfTheseSamFlagsSet") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            int64_t bitString = atoi(tokStr);
            if(bitString < 0 || bitString > UINT16_MAX) {
                st_errAbort("ERROR: Attempting to set 16-bit string with invalid argument: %s", tokStr);
            }
            params->filterAReadWithAnyOneOfTheseSamFlagsSet = bitString;
            i++;
        }
        else if (strcmp(keyString, "estimateReadErrorProbsEmpirically") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            assert(strcmp(tokStr, "true") || strcmp(tokStr, "false"));
            params->estimateReadErrorProbsEmpirically = strcmp(tokStr, "true") == 0;
            i++;
        }
        else if (strcmp(keyString, "roundsOfIterativeRefinement") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            params->roundsOfIterativeRefinement = atoi(tokStr);
            i++;
        }
        else if (strcmp(keyString, "compareVCFs") == 0) {
            //removed, but legacy params may still contain this
            i++;
        }
        else if (strcmp(keyString, "writeGVCF") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            assert(strcmp(tokStr, "true") || strcmp(tokStr, "false"));
            params->writeGVCF = strcmp(tokStr, "true") == 0;
            i++;
        } else if (strcmp(keyString, "writeSplitSams") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            assert(strcmp(tokStr, "true") || strcmp(tokStr, "false"));
            params->writeSplitSams = strcmp(tokStr, "true") == 0;
            i++;
        }else if (strcmp(keyString, "writeUnifiedSam") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            assert(strcmp(tokStr, "true") || strcmp(tokStr, "false"));
            params->writeUnifiedSam = strcmp(tokStr, "true") == 0;
            i++;
        }
        else if (strcmp(keyString, "mapqFilter") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = json_token_tostr(js, &tok);
            params->mapqFilter = atoi(tokStr);
            i++;
        }
        else {
            st_errAbort("ERROR: Unrecognised key in params file: %s\n", keyString);
        }

    }

    free(js);
    return params;
}

/*
 * Sets the level of verbosity for vcf comparison.
 */
void setVerbosity(stRPHmmParameters *params, int64_t bitstring) {
    params->verboseTruePositives = (bitstring & LOG_TRUE_POSITIVES) != 0;
    params->verboseFalsePositives = (bitstring & LOG_FALSE_POSITIVES) != 0;
    params->verboseFalseNegatives = (bitstring & LOG_FALSE_NEGATIVES) != 0;
}

/*
 * Counts the number of insertions and deletions in a read, given its cigar string.
 */
void countIndels(uint32_t *cigar, uint32_t ncigar, int64_t *numInsertions, int64_t *numDeletions) {
    for (uint32_t i = 0; i < ncigar; i++) {
        int cigarOp = cigar[i] & BAM_CIGAR_MASK;
        int cigarNum = cigar[i] >> BAM_CIGAR_SHIFT;
        if (cigarOp == BAM_CINS) *numInsertions += cigarNum;
        if (cigarOp == BAM_CDEL) *numDeletions += cigarNum;
    }
}

void appendProbsToList(stList *probabilityList, uint8_t pA, uint8_t pC, uint8_t pG, uint8_t pT, uint8_t pGap) {
    uint8_t *aPtr = calloc(1, sizeof(uint8_t));
    uint8_t *cPtr = calloc(1, sizeof(uint8_t));
    uint8_t *gPtr = calloc(1, sizeof(uint8_t));
    uint8_t *tPtr = calloc(1, sizeof(uint8_t));
    uint8_t *gapPtr = calloc(1, sizeof(uint8_t));
    *aPtr = pA;
    *cPtr = pC;
    *gPtr = pG;
    *tPtr = pT;
    *gapPtr = pGap;
    stList_append(probabilityList, aPtr);
    stList_append(probabilityList, cPtr);
    stList_append(probabilityList, gPtr);
    stList_append(probabilityList, tPtr);
    stList_append(probabilityList, gapPtr);
}

stProfileSeq* getProfileSequenceFromSingleNuclProbFile(char *signalAlignReadLocation, char *readName,
                                                       stBaseMapper *baseMapper, stRPHmmParameters *params) {
    // get signalAlign file
    FILE *fp = fopen(signalAlignReadLocation,"r");

    // for scanning the file
    int fieldSize = 2048;
    char *line = calloc(2048, sizeof(char));
    char *chromStr = calloc(fieldSize, sizeof(char));
    char *refPosStr = calloc(fieldSize, sizeof(char));
    char *pAStr = calloc(fieldSize, sizeof(char));
    char *pCStr = calloc(fieldSize, sizeof(char));
    char *pGStr = calloc(fieldSize, sizeof(char));
    char *pTStr = calloc(fieldSize, sizeof(char));
    char *pGapStr = calloc(fieldSize, sizeof(char));

    // for handling the data
    int64_t refPos;
    uint8_t pA;
    uint8_t pC;
    uint8_t pG;
    uint8_t pT;
    uint8_t pGap;
    uint8_t *aPtr;
    uint8_t *cPtr;
    uint8_t *gPtr;
    uint8_t *tPtr;
    uint8_t *gapPtr;

    // parse header
    while(!feof(fp)) {
        fscanf( fp, "%[^\n]\n", line);
        if (line[0] == '#') {
            if (line[1] == '#') continue;
            if (strcmp(line, "#CHROM\tPOS\tpA\tpC\tpG\tpT\tp_") != 0) {
                st_errAbort("SignalAlign output file %s has unexpected header format: %s",
                            signalAlignReadLocation, line);
            } else {
                break;
            }
        }
    }

    // get probabilities
    stList* probabilityList = stList_construct3(0, free);
    uint64_t firstReadPos = 0;
    uint64_t lastReadPos = 0;
    int64_t randomSeed = st_randomInt64(0,3);
    while(!feof(fp)) {
        // Scan
        fscanf( fp, "%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\n]\n",
                chromStr, refPosStr, pAStr, pCStr, pGStr, pTStr, pGapStr);

        // Get reference position
        refPos = atoi(refPosStr);

        // Check for gaps todo this might actually be a bug or something in signalAlign
        while (firstReadPos != 0 && refPos > lastReadPos + 1) {
            appendProbsToList(probabilityList, ALPHABET_MIN_PROB, ALPHABET_MIN_PROB, ALPHABET_MIN_PROB,
                              ALPHABET_MIN_PROB, ALPHABET_MAX_PROB);
            lastReadPos++;
        }

        // Check for inserts
        if (firstReadPos != 0 && refPos < lastReadPos + 1) continue;

        // Update first and last read position
        if (firstReadPos == 0) firstReadPos = refPos;
        lastReadPos = refPos;

        // Get probabilities and save
        pA = (uint8_t) (ALPHABET_MAX_PROB * atof(pAStr));
        pC = (uint8_t) (ALPHABET_MAX_PROB * atof(pCStr));
        pG = (uint8_t) (ALPHABET_MAX_PROB * atof(pGStr));
        pT = (uint8_t) (ALPHABET_MAX_PROB * atof(pTStr));
        pGap = (uint8_t) (ALPHABET_MAX_PROB * atof(pGapStr));

        // We need all integer probs to sum to MAX_PROB todo is there a way to do this better?
        while ((pA + pC + pG + pT + pGap) > ALPHABET_MAX_PROB) {
            if ((pA + pC + pG + pT + pGap) == 0) {
                break;
            }
            switch (randomSeed++ % 5) {
                case 0: if (pA != 0) pA--; break;
                case 1: if (pC != 0) pC--; break;
                case 2: if (pG != 0) pG--; break;
                case 3: if (pT != 0) pT--; break;
                case 4: if (pGap != 0) pGap--; break;
                default: assert(FALSE);
            }
        }
        while ((pA + pC + pG + pT + pGap) < ALPHABET_MAX_PROB) {
            if ((pA + pC + pG + pT + pGap) == 0) {
                break;
            }
            switch (randomSeed++ % 5) {
                case 0: if (pA != 0) pA++; break;
                case 1: if (pC != 0) pC++; break;
                case 2: if (pG != 0) pG++; break;
                case 3: if (pT != 0) pT++; break;
                case 4: if (pGap != 0) pGap++; break;
                default: assert(FALSE);
            }
        }

        // Save the values in a list todo this is not particularly efficient
        appendProbsToList(probabilityList, pA, pC, pG, pT, pGap);
    }
    // Now we're done with the file
    fclose(fp);

    // Create empty profile sequence
    uint64_t readLength = lastReadPos - firstReadPos + 1;
    stProfileSeq *pSeq = stProfileSeq_constructEmptyProfile(chromStr, readName, firstReadPos + 1, readLength);

    // Copy probabilities over
    uint64_t position = 0;
    stListIterator *itor = stList_getIterator(probabilityList);
    while (position < readLength) {

        // Get the locations of the probabilities
        aPtr = stList_getNext(itor);
        cPtr = stList_getNext(itor);
        gPtr = stList_getNext(itor);
        tPtr = stList_getNext(itor);
        gapPtr = stList_getNext(itor);

        // Assign the probabilities
        pSeq->profileProbs[position * ALPHABET_SIZE + stBaseMapper_getValueForChar(baseMapper, 'A')] =  *aPtr;
        pSeq->profileProbs[position * ALPHABET_SIZE + stBaseMapper_getValueForChar(baseMapper, 'C')] =  *cPtr;
        pSeq->profileProbs[position * ALPHABET_SIZE + stBaseMapper_getValueForChar(baseMapper, 'G')] =  *gPtr;
        pSeq->profileProbs[position * ALPHABET_SIZE + stBaseMapper_getValueForChar(baseMapper, 'T')] =  *tPtr;

        if (params->gapCharactersForDeletions) {
            // This assumes gap character is the last character in the alphabet given
            pSeq->profileProbs[position * ALPHABET_SIZE + (ALPHABET_SIZE - 1)] = *gapPtr;
        } else {
            pSeq->profileProbs[position * ALPHABET_SIZE + (ALPHABET_SIZE - 1)] = ALPHABET_MIN_PROB;
        }

        position++;
    }
    // We should have nothing left over in the list
    if (stList_getNext(itor) != NULL) {
        st_errAbort("Probability list has %d extra elements, with read length %d for file %s",
                    stList_length(probabilityList) - 5 * position, readLength, signalAlignReadLocation);
    }
    // Sanity check on the number of modifications to the probabilities
    // We only modify probability of bases with some probability, so to fix a rounding error, we should at worst have
    //  to make 4 modifications per location
    if (randomSeed > (4 * readLength)) {
        st_logDebug("\t\tNeeded average of %f modifications to base probs to ensure proper total probability for %s\n",
                    (1.0 * randomSeed / readLength), readName);
    }

    stList_destruct(probabilityList);
    free(line);
    free(chromStr);
    free(refPosStr);
    free(pAStr);
    free(pCStr);
    free(pGStr);
    free(pTStr);

    return pSeq;
}



/* Parse reads within an input interval of a reference sequence of a bam file
 * and create a list of profile sequences by turning characters into profile probabilities.
 *
 * In future, maybe use mapq scores to adjust profile (or posterior probabilities for
 * signal level alignments).
 * */
int64_t parseReads(stList *profileSequences, char *bamFile, stBaseMapper *baseMapper, stRPHmmParameters *params) {
    return parseReadsWithSingleNucleotideProbs(profileSequences, bamFile, baseMapper, params, NULL, false);
}

int64_t parseReadsWithSingleNucleotideProbs(stList *profileSequences, char *bamFile, stBaseMapper *baseMapper,
                                            stRPHmmParameters *params, char *singleNuclProbDirectory,
                                            bool onlySingleNuclProb) {
    if (singleNuclProbDirectory != NULL) {
        st_logInfo("\tModifying probabilities from single nucleotide probability files in %s\n",
                   singleNuclProbDirectory);
    }

    samFile *in = hts_open(bamFile, "r");
    if (in == NULL) {
        st_errAbort("ERROR: Cannot open bam file %s\n", bamFile);
        return -1;
    }
    bam_hdr_t *bamHdr = sam_hdr_read(in);
    bam1_t *aln = bam_init1();

    int64_t readCount = 0;
    int64_t singleNuclProbReadCount = 0;
    int64_t bamReadCount = 0;
    int64_t profileCount = 0;
    int64_t missingSingleNuclProbReads = 0;
    int64_t filteredReads = 0;
    int64_t filteredReads_flag = 0;
    int64_t filteredReads_mapq = 0;

    while(sam_read1(in,bamHdr,aln) > 0) {
        stProfileSeq *pSeq = NULL;

        int64_t pos = aln->core.pos+1;                      // Left most position of alignment
        char *chr = bamHdr->target_name[aln->core.tid] ;    // Contig name (chromosome)
        int64_t len = aln->core.l_qseq;                     // Length of the read.
        uint8_t *seq = bam_get_seq(aln);                    // DNA sequence
        char *readName = bam_get_qname(aln);
        uint32_t *cigar = bam_get_cigar(aln);

        if (aln->core.l_qseq <= 0) {
            filteredReads++;
            continue;
        }

        // Filter out any reads with specified flags
        if((aln->core.flag & params->filterAReadWithAnyOneOfTheseSamFlagsSet) > 0) {
            filteredReads++;
            filteredReads_flag++;
            continue;
        }

        // If there isn't a cigar string, don't bother including the read, since we don't
        // know how it aligns
        if (aln->core.n_cigar == 0) {
            filteredReads++;
            continue;
        }

        // If the mapq score is less than the given threshold, filter it out
        if (aln->core.qual < params->mapqFilter) {
            filteredReads++;
            filteredReads_mapq++;
            continue;
        }

        // Tracks how many reads there were
        readCount++;

        // Should we read from the signalAlign directory?
        if (singleNuclProbDirectory != NULL) {

            // Get signalAlign file (if exists)
            char *singleNuclProbReadLocation = stString_print("%s/%s.tsv", singleNuclProbDirectory, readName);
            if (access(singleNuclProbReadLocation, F_OK) == -1) {
                // Could not find the read file
                missingSingleNuclProbReads++;
            } else {
                // Found the read file
                pSeq = getProfileSequenceFromSingleNuclProbFile(singleNuclProbReadLocation, readName, baseMapper, params);
                singleNuclProbReadCount++;

                // We have a profile, so save it
                stList_append(profileSequences, pSeq);
                profileCount++;
            }
            free(singleNuclProbReadLocation);

            // If we found a SA file or if we don't want missing reads
            if (pSeq != NULL || onlySingleNuclProb) {
                continue;
            }
        }

        int64_t start_read = 0;
        int64_t end_read = 0;
        int64_t start_ref = pos;
        int64_t cig_idx = 0;

        // Find the correct starting locations on the read and reference sequence,
        // to deal with things like inserts / deletions / soft clipping
        while (cig_idx < aln->core.n_cigar) {
            int cigarOp = cigar[cig_idx] & BAM_CIGAR_MASK;
            int cigarNum = cigar[cig_idx] >> BAM_CIGAR_SHIFT;

            if (cigarOp == BAM_CMATCH || cigarOp == BAM_CEQUAL || cigarOp == BAM_CDIFF) {
                break;
            } else if (cigarOp == BAM_CDEL || cigarOp == BAM_CREF_SKIP) {
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
        if (aln->core.n_cigar > 1) {
            int lastCigarOp = cigar[aln->core.n_cigar-1] & BAM_CIGAR_MASK;
            int lastCigarNum = cigar[aln->core.n_cigar-1] >> BAM_CIGAR_SHIFT;
            if (lastCigarOp == BAM_CSOFT_CLIP) {
                end_read += lastCigarNum;
            }
        }

        // Count number of insertions & deletions in sequence
        int64_t numInsertions = 0;
        int64_t numDeletions = 0;
        countIndels(cigar, aln->core.n_cigar, &numInsertions, &numDeletions);
        int64_t trueLength = len - start_read - end_read + numDeletions - numInsertions;

        if (trueLength <= 0) {
            filteredReads++;
            continue;
        }

        // Create empty profile sequence
        pSeq = stProfileSeq_constructEmptyProfile(chr, readName, pos, trueLength);

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
            } else if (cigarOp == BAM_CDEL || cigarOp == BAM_CREF_SKIP) {
                // Only add a gap character when that param is on
                if (params->gapCharactersForDeletions) {
                    // This assumes gap character is the last character in the alphabet given
                    pSeq->profileProbs[i * ALPHABET_SIZE + (ALPHABET_SIZE - 1)] = ALPHABET_MAX_PROB;
                }
            } else if (cigarOp == BAM_CINS) {
                // Currently, ignore insertions
                idxInSeq++;
                i--;
            } else if (cigarOp == BAM_CSOFT_CLIP || cigarOp == BAM_CHARD_CLIP || cigarOp == BAM_CPAD) {
                // Nothing really to do here. skip to next cigar operation
                currPosInOp = cigarNum - 1;
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
        bamReadCount++;

        // Save profile seq
        if (pSeq->length > 0) {
            profileCount++;
            stList_append(profileSequences, pSeq);
        }

    }

    // Log signal align usage
    if (singleNuclProbDirectory != NULL) {
        if (missingSingleNuclProbReads > 0) {
            st_logInfo("\t%d/%d reads were missing single nucleotide probability file\n", missingSingleNuclProbReads, readCount);
        }
        st_logInfo("\tOf %d total reads: %d were loaded from single nucleotide probability data, and %d were from the bam\n",
                   profileCount, singleNuclProbReadCount, bamReadCount);

    }

    // Log filtering actions
    if(st_getLogLevel() == debug) {
        char *samFlagBitString = intToBinaryString(params->filterAReadWithAnyOneOfTheseSamFlagsSet);
        st_logDebug("\tFiltered %" PRIi64 " reads with either missing cigar lines, "
                            "\n\t\tlow mapq scores (filtered %d reads with scores less than %d), "
                            "\n\t\tand undesired sam flags "
                            "(filtered %d reads with sam flags being filtered on: %s)\n",
                filteredReads, filteredReads_mapq, params->mapqFilter, filteredReads_flag, samFlagBitString);
        free(samFlagBitString);
    }

    // Sanity check (did we accidentally save profile sequences twice?)
    assert(stList_length(profileSequences) <= readCount);

    bam_hdr_destroy(bamHdr);
    bam_destroy1(aln);
    sam_close(in);

    return profileCount;
}

