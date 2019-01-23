/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com) & Arthur Rand (arand@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "margin.h"

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
char stBaseMapper_getCharForValue(stBaseMapper *bm, uint64_t value) {
    char base = bm->numToChar[value];
    if (base >= 0) return base;
    st_errAbort("Value '%d' not specified in alphabet", value);
    return -1;
}

/*
 * Get model parameters from params file.
 * Set hmm parameters.
*/

stRPHmmParameters *stRPHmmParameters_construct() {
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

	return params;
}

stRPHmmParameters *phaseParams_fromJson(char *buf, size_t r, stBaseMapper *baseMapper) {
	// Setup parser
	jsmntok_t *tokens;
	char *js;
	int64_t tokenNumber = stJson_setupParser(buf, r, &tokens, &js);

	stRPHmmParameters *params = stRPHmmParameters_construct();

    setVerbosity(params, 0);

    //TODO: refactor the following to use the json parsing functions

    // Parse tokens, starting at token 1
    // (token 0 is entire object)
    for (int64_t i = 1; i < tokenNumber; i++) {
        jsmntok_t key = tokens[i];
        char *keyString = stJson_token_tostr(js, &key);

        if (strcmp(keyString, "alphabet") == 0) {
            jsmntok_t alphabetTok = tokens[i+1];
            if (alphabetTok.size != ALPHABET_SIZE) {
                st_errAbort("Alphabet size in JSON does not match constant ALPHABET_SIZE \n");
            }
            for (int j = 0; j < ALPHABET_SIZE; j++) {
                jsmntok_t tok = tokens[i+j+2];
                char *tokStr = stJson_token_tostr(js, &tok);
                stBaseMapper_addBases(baseMapper, tokStr);
            }
            i += ALPHABET_SIZE + 1;
        }
        else if (strcmp(keyString, "wildcard") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = stJson_token_tostr(js, &tok);
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
                char *tokStr = stJson_token_tostr(js, &tok);
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
                char *tokStr = stJson_token_tostr(js, &tok);
                setSubstitutionProb(params->readErrorSubModel, params->readErrorSubModelSlow,
                                    j/ALPHABET_SIZE, j%ALPHABET_SIZE, atof(tokStr));

            }
            i += readErrTok.size + 1;
        }
        else if (strcmp(keyString, "maxNotSumTransitions") == 0) {
        	params->maxNotSumTransitions = stJson_parseBool(js, tokens, ++i);
        }
        else if (strcmp(keyString, "minPartitionsInAColumn") == 0) {
        	params->minPartitionsInAColumn = stJson_parseInt(js, tokens, ++i);
        }
        else if (strcmp(keyString, "maxPartitionsInAColumn") == 0) {
        	params->maxPartitionsInAColumn = stJson_parseInt(js, tokens, ++i);
        }
        else if (strcmp(keyString, "minPosteriorProbabilityForPartition") == 0) {
        	params->minPosteriorProbabilityForPartition = stJson_parseFloat(js, tokens, ++i);
        }
        else if (strcmp(keyString, "maxCoverageDepth") == 0) {
        	params->maxCoverageDepth = stJson_parseInt(js, tokens, ++i);
        }
        else if (strcmp(keyString, "minReadCoverageToSupportPhasingBetweenHeterozygousSites") == 0) {
        	params->minReadCoverageToSupportPhasingBetweenHeterozygousSites = stJson_parseInt(js, tokens, ++i);
        }
        else if (strcmp(keyString, "onDiagonalReadErrorPseudoCount") == 0) {
        	params->onDiagonalReadErrorPseudoCount = stJson_parseFloat(js, tokens, ++i);
        }
        else if (strcmp(keyString, "offDiagonalReadErrorPseudoCount") == 0) {
        	params->offDiagonalReadErrorPseudoCount = stJson_parseFloat(js, tokens, ++i);
        }
        else if (strcmp(keyString, "trainingIterations") == 0) {
        	params->trainingIterations = stJson_parseInt(js, tokens, ++i);
        }
        else if (strcmp(keyString, "filterBadReads") == 0) {
        	params->filterBadReads = stJson_parseBool(js, tokens, ++i);
        }
        else if (strcmp(keyString, "filterMatchThreshold") == 0) {
        	params->filterMatchThreshold = stJson_parseFloat(js, tokens, ++i);
        }
        else if (strcmp(keyString, "useReferencePrior") == 0) {
        	params->useReferencePrior = stJson_parseBool(js, tokens, ++i);
        }
        else if (strcmp(keyString, "includeInvertedPartitions") == 0) {
        	params->includeInvertedPartitions = stJson_parseBool(js, tokens, ++i);
        }
        else if (strcmp(keyString, "verbose") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = stJson_token_tostr(js, &tok);
            int64_t bitString = atoi(tokStr);
            setVerbosity(params, bitString);
            i++;
        }
        else if(strcmp(keyString, "filterLikelyHomozygousSites") == 0) {
        	params->filterLikelyHomozygousSites = stJson_parseBool(js, tokens, ++i);
        }
        else if (strcmp(keyString, "minSecondMostFrequentBaseFilter") == 0) {
        	params->minSecondMostFrequentBaseFilter = stJson_parseFloat(js, tokens, ++i);
        }
        else if (strcmp(keyString, "minSecondMostFrequentBaseLogProbFilter") == 0) {
        	params->minSecondMostFrequentBaseLogProbFilter = stJson_parseFloat(js, tokens, ++i);
        }
        else if (strcmp(keyString, "gapCharactersForDeletions") == 0) {
        	params->gapCharactersForDeletions = stJson_parseBool(js, tokens, ++i);
        }
        else if (strcmp(keyString, "filterAReadWithAnyOneOfTheseSamFlagsSet") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = stJson_token_tostr(js, &tok);
            int64_t bitString = atoi(tokStr);
            if(bitString < 0 || bitString > UINT16_MAX) {
                st_errAbort("ERROR: Attempting to set 16-bit string with invalid argument: %s", tokStr);
            }
            params->filterAReadWithAnyOneOfTheseSamFlagsSet = bitString;
            i++;
        }
        else if (strcmp(keyString, "estimateReadErrorProbsEmpirically") == 0) {
        	params->estimateReadErrorProbsEmpirically = stJson_parseBool(js, tokens, ++i);
        }
        else if (strcmp(keyString, "roundsOfIterativeRefinement") == 0) {
        	params->roundsOfIterativeRefinement = stJson_parseInt(js, tokens, ++i);
        }
        else if (strcmp(keyString, "compareVCFs") == 0) {
            //removed, but legacy params may still contain this
            i++;
        }
        else if (strcmp(keyString, "writeGVCF") == 0) {
        	params->writeGVCF = stJson_parseBool(js, tokens, ++i);
        } else if (strcmp(keyString, "writeSplitSams") == 0) {
        	params->writeSplitSams = stJson_parseBool(js, tokens, ++i);
        }else if (strcmp(keyString, "writeUnifiedSam") == 0) {
        	params->writeUnifiedSam = stJson_parseBool(js, tokens, ++i);
        }
        else if (strcmp(keyString, "mapqFilter") == 0) {
        	params->mapqFilter = stJson_parseInt(js, tokens, ++i);
        }
        else {
            st_errAbort("ERROR: Unrecognised key in params file: %s\n", keyString);
        }
    }

    // Cleanup
    free(js);
    free(tokens);

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
 * Params object for polisher
 */

int64_t repeatSubMatrix_parseLogProbabilities(RepeatSubMatrix *repeatSubMatrix, Symbol base, bool strand, char *js, jsmntok_t *tokens, int64_t tokenIndex) {
	int64_t maxRepeatCount = repeatSubMatrix->maximumRepeatLength;
	int64_t i = stJson_parseFloatArray(repeatSubMatrix_setLogProb(repeatSubMatrix, base, strand, 0, 0), maxRepeatCount*maxRepeatCount, js, tokens, tokenIndex);
	return i;
}

RepeatSubMatrix *repeatSubMatrix_jsonParse(char *buf, size_t r) {
	// Setup parser
	jsmntok_t *tokens;
	char *js;
	int64_t tokenNumber = stJson_setupParser(buf, r, &tokens, &js);

	RepeatSubMatrix *repeatSubMatrix = repeatSubMatrix_constructEmpty();

	for(int64_t tokenIndex=1; tokenIndex < tokenNumber; tokenIndex++) {
		jsmntok_t key = tokens[tokenIndex];
		char *keyString = stJson_token_tostr(js, &key);
		if(strlen(keyString) != 31) {
			st_errAbort("ERROR: Unrecognised key in repeat sub matrix json: %s\n", keyString);
		}
		char base = keyString[28];
		if(base != 'A' && base != 'C' && base != 'G' && base != 'T') {
			st_errAbort("ERROR: Unrecognised base in repeat sub matrix json: %s, base=%c\n", keyString, base);
		}
		if(keyString[30] != 'F' && keyString[30] != 'R') {
			st_errAbort("ERROR: Unrecognised strand in repeat sub matrix json: %s, strand:%c\n", keyString, keyString[30]);
		}
		bool strand = keyString[30] == 'F';
		tokenIndex = repeatSubMatrix_parseLogProbabilities(repeatSubMatrix, symbol_convertCharToSymbol(base), strand, js, tokens, tokenIndex+1);
	}

	// Cleanup
	free(js);
	free(tokens);

	return repeatSubMatrix;
}

PolishParams *polishParams_jsonParse(char *buf, size_t r) {
	// Setup parser
	jsmntok_t *tokens;
	char *js;
	int64_t tokenNumber = stJson_setupParser(buf, r, &tokens, &js);

	// Make empty params object
	PolishParams *params = st_calloc(1, sizeof(PolishParams));

	// Intelligent defaults
	params->useRunLengthEncoding = 1;
	params->referenceBasePenalty = 0.5;
	params->minPosteriorProbForAlignmentAnchors = st_calloc(2, sizeof(double));
	params->minPosteriorProbForAlignmentAnchors[0] = 0.9;
	params->minPosteriorProbForAlignmentAnchors[0] = 10;
	params->minPosteriorProbForAlignmentAnchorsLength = 2;
    params->includeSoftClipping = FALSE; //todo add this in
    params->chunkSize = 0;
    params->chunkBoundary = 0;
    params->candidateVariantWeight = 0.2;
    params->columnAnchorTrim = 5;
    params->maxConsensusStrings = 100;

	// Parse tokens, starting at token 1
    // (token 0 is entire object)
	bool gotHmm = 0, gotPairwiseAlignmentParameters = 0, gotRepeatCountMatrix = 0;
    for (int64_t tokenIndex=1; tokenIndex < tokenNumber; tokenIndex++) {
        jsmntok_t key = tokens[tokenIndex];
        char *keyString = stJson_token_tostr(js, &key);

        if (strcmp(keyString, "useRunLengthEncoding") == 0) {
        	params->useRunLengthEncoding = stJson_parseBool(js, tokens, ++tokenIndex);
        }
        else if (strcmp(keyString, "referenceBasePenalty") == 0) {
        	params->referenceBasePenalty = stJson_parseFloat(js, tokens, ++tokenIndex);
        }
        else if (strcmp(keyString, "minPosteriorProbForAlignmentAnchors") == 0) {
        	free(params->minPosteriorProbForAlignmentAnchors); //Cleanup the old one
        	int64_t arraySize = tokens[tokenIndex+1].size;
        	if(arraySize % 2 != 0 && arraySize > 0) {
        		st_errAbort("ERROR: length of minPosteriorProbForAlignmentAnchors must be even and greater than zero\n");
        	}
        	params->minPosteriorProbForAlignmentAnchors = st_calloc(arraySize, sizeof(double));
        	params->minPosteriorProbForAlignmentAnchorsLength = arraySize;
        	tokenIndex = stJson_parseFloatArray(params->minPosteriorProbForAlignmentAnchors,
        			 	 	 	 	 	 	 	 arraySize, js, tokens, ++tokenIndex);
        	double pValue = 0.0;
        	for(int64_t i=0; i<params->minPosteriorProbForAlignmentAnchorsLength; i+=2) {
        		if(params->minPosteriorProbForAlignmentAnchors[i] < pValue || params->minPosteriorProbForAlignmentAnchors[i] > 1) {
        			st_errAbort("ERROR: minPosteriorProbForAlignmentAnchors must be even, greater than zero and increasing\n");
        		}
        		pValue = params->minPosteriorProbForAlignmentAnchors[i];
        		if(params->minPosteriorProbForAlignmentAnchors[i+1] < 0 || ((int64_t)params->minPosteriorProbForAlignmentAnchors[i+1]) % 2 != 0) {
        			st_errAbort("ERROR: minPosteriorProbForAlignmentAnchors diagonal expansion must be, greater than zero and even\n");
        		}
        	}
        }
        else if (strcmp(keyString, "repeatCountSubstitutionMatrix") == 0) {
        	jsmntok_t tok = tokens[tokenIndex+1];
        	char *tokStr = stJson_token_tostr(js, &tok);
        	params->repeatSubMatrix = repeatSubMatrix_jsonParse(tokStr, strlen(tokStr));
        	tokenIndex += stJson_getNestedTokenCount(tokens, tokenIndex+1);
        	gotRepeatCountMatrix = 1;
        }
        else if (strcmp(keyString, "hmm") == 0) {
        	jsmntok_t tok = tokens[tokenIndex+1];
        	char *tokStr = stJson_token_tostr(js, &tok);
        	params->hmm = hmm_jsonParse(tokStr, strlen(tokStr));
        	params->sM = hmm_getStateMachine(params->hmm);
        	tokenIndex += stJson_getNestedTokenCount(tokens, tokenIndex+1);
        	gotHmm = 1;
        }
        else if (strcmp(keyString, "pairwiseAlignmentParameters") == 0) {
        	jsmntok_t tok = tokens[tokenIndex+1];
        	char *tokStr = stJson_token_tostr(js, &tok);
        	params->p = pairwiseAlignmentParameters_jsonParse(tokStr, strlen(tokStr));
            if (params->p->diagonalExpansion % 2 != 0) {
                st_errAbort("ERROR: pairwiseAlignmentParameters.diagonalExpansion must be even\n");
            }
        	tokenIndex += stJson_getNestedTokenCount(tokens, tokenIndex+1);
        	gotPairwiseAlignmentParameters = 1;
        }
        else if (strcmp(keyString, "includeSoftClipping") == 0) {
            params->includeSoftClipping = stJson_parseBool(js, tokens, ++tokenIndex);
        }
        else if (strcmp(keyString, "chunkSize") == 0) {
            if (stJson_parseInt(js, tokens, ++tokenIndex) < 0) {
                st_errAbort("ERROR: chunkSize parameter must zero or greater\n");
            }
            params->chunkSize = (uint64_t) stJson_parseInt(js, tokens, tokenIndex);
        }
        else if (strcmp(keyString, "chunkBoundary") == 0) {
            if (stJson_parseInt(js, tokens, ++tokenIndex) < 0) {
                st_errAbort("ERROR: chunkBoundary parameter must zero or greater\n");
            }
            params->chunkBoundary = (uint64_t) stJson_parseInt(js, tokens, tokenIndex);
        }
        else if (strcmp(keyString, "candidateVariantWeight") == 0) {
			if (stJson_parseFloat(js, tokens, ++tokenIndex) < 0) {
				st_errAbort("ERROR: candidateVariantWeight parameter must zero or greater\n");
			}
			params->candidateVariantWeight = stJson_parseFloat(js, tokens, tokenIndex);
		} else if (strcmp(keyString, "columnAnchorTrim") == 0) {
			if (stJson_parseInt(js, tokens, ++tokenIndex) < 0) {
				st_errAbort("ERROR: columnAnchorTrim parameter must zero or greater\n");
			}
			params->columnAnchorTrim = (uint64_t) stJson_parseInt(js, tokens, tokenIndex);
		} else if (strcmp(keyString, "maxConsensusStrings") == 0) {
			if (stJson_parseInt(js, tokens, ++tokenIndex) < 0) {
				st_errAbort("ERROR: maxConsensusStrings parameter must zero or greater\n");
			}
			params->maxConsensusStrings = (uint64_t) stJson_parseInt(js, tokens, tokenIndex);
		}
		else if (strcmp(keyString, "maxPoaConsensusIterations") == 0) {
			if (stJson_parseInt(js, tokens, ++tokenIndex) < 0) {
				st_errAbort("ERROR: maxPoaConsensusIterations parameter must zero or greater\n");
			}
			params->maxPoaConsensusIterations = (uint64_t) stJson_parseInt(js, tokens, tokenIndex);
		} else if (strcmp(keyString, "minPoaConsensusIterations") == 0) {
			if (stJson_parseInt(js, tokens, ++tokenIndex) < 0) {
				st_errAbort("ERROR: minPoaConsensusIterations parameter must zero or greater\n");
			}
			params->minPoaConsensusIterations = (uint64_t) stJson_parseInt(js, tokens, tokenIndex);
		} else if (strcmp(keyString, "maxRealignmentPolishIterations") == 0) {
			if (stJson_parseInt(js, tokens, ++tokenIndex) < 0) {
				st_errAbort("ERROR: maxRealignmentPolishIterations parameter must zero or greater\n");
			}
			params->maxRealignmentPolishIterations = (uint64_t) stJson_parseInt(js, tokens, tokenIndex);
		} else if (strcmp(keyString, "minRealignmentPolishIterations") == 0) {
			if (stJson_parseInt(js, tokens, ++tokenIndex) < 0) {
				st_errAbort("ERROR: minRealignmentPolishIterations parameter must zero or greater\n");
			}
			params->minRealignmentPolishIterations = (uint64_t) stJson_parseInt(js, tokens, tokenIndex);
		}
        else {
            st_errAbort("ERROR: Unrecognised key in polish params json: %s\n", keyString);
        }
    }

    if(!gotRepeatCountMatrix) {
    	st_errAbort("ERROR: Did not find repeat counts specified in json polish params\n");
    }
    if(!gotHmm) {
    	st_errAbort("ERROR: Did not find HMM specified in json polish params\n");
    }
    if(!gotPairwiseAlignmentParameters) {
    	st_errAbort("ERROR: Did not find pairwise alignment params specified in json polish params\n");
    }

    // Cleanup
    free(js);
    free(tokens);

    return params;
}

void polishParams_printParameters(PolishParams *polishParams, FILE *fh) {
    //TODO
    st_logCritical("Need to implement polishParams_printParameters\n");
}

void polishParams_destruct(PolishParams *params) {
	repeatSubMatrix_destruct(params->repeatSubMatrix);
	stateMachine_destruct(params->sM);
	hmm_destruct(params->hmm);
	pairwiseAlignmentBandingParameters_destruct(params->p);
	free(params);
}

/*
 * Global params objects
 */
Params *params_jsonParse(char *buf, size_t r, bool requirePolish, bool requirePhase) {
	// Setup parser
	jsmntok_t *tokens;
	char *js;
	int64_t tokenNumber = stJson_setupParser(buf, r, &tokens, &js);

	// Make empty params object
	Params *params = st_calloc(1, sizeof(Params));

	// Parse tokens, starting at token 1
	// (token 0 is entire object)
	bool gotPolish = 0, gotPhase = 0;
	for (int64_t tokenIndex=1; tokenIndex < tokenNumber; tokenIndex++) {
		jsmntok_t key = tokens[tokenIndex];
		char *keyString = stJson_token_tostr(js, &key);

		if (strcmp(keyString, "polish") == 0) {
		    if (requirePolish) {
                jsmntok_t tok = tokens[tokenIndex + 1];
                char *tokStr = stJson_token_tostr(js, &tok);
                params->polishParams = polishParams_jsonParse(tokStr, strlen(tokStr));
            }
			tokenIndex += stJson_getNestedTokenCount(tokens, tokenIndex+1);
			gotPolish = 1;
		}
		else if (strcmp(keyString, "phase") == 0) {
            if (requirePhase) {
                jsmntok_t tok = tokens[tokenIndex + 1];
                char *tokStr = stJson_token_tostr(js, &tok);
                params->baseMapper = stBaseMapper_construct();
                params->phaseParams = phaseParams_fromJson(tokStr, strlen(tokStr), params->baseMapper);
            }
			tokenIndex += stJson_getNestedTokenCount(tokens, tokenIndex+1);
			gotPhase = 1;
		}
		else {
		    // maybe we only need one of these?
		    if (requirePolish && !requirePhase) {
		        st_logInfo("WARN: parameters file missing 'polish' and 'phase' top-level entries.  Interpreting as 'polish'.\n");
                params->polishParams = polishParams_jsonParse(buf, r);
                tokenIndex = tokenNumber;
                gotPolish = 1;
		    } else if (requirePhase && ! requirePolish) {
                st_logInfo("WARN: parameters file missing 'polish' and 'phase' top-level entries.  Interpreting as 'phase'.\n");
                params->baseMapper = stBaseMapper_construct();
                params->phaseParams = phaseParams_fromJson(buf, r, params->baseMapper);
                tokenIndex = tokenNumber;
		        gotPhase = 1;
		    } else {
                st_errAbort("ERROR: Unrecognised key in params json: %s\n", keyString);
            }
		}
	}

	if(!gotPolish && requirePolish) {
		st_errAbort("ERROR: Did not find polish parameters in json params\n");
	}
	if(!gotPhase && requirePhase) {
		st_errAbort("ERROR: Did not find phase parameters in json params\n");
	}

	// Cleanup
	free(js);
	free(tokens);

	return params;
}

Params *params_readParams(FILE *fp) {
    return params_readParams2(fp, TRUE, TRUE);
}
Params *params_readParams2(FILE *fp, bool requirePolish, bool requirePhase) {
    char buf[BUFSIZ * 300]; // TODO: FIX, This is terrible code, we should not assume the size of the file is less than this
    return params_jsonParse(buf, fread(buf, sizeof(char), sizeof(buf), fp), requirePolish, requirePhase);
}

void params_destruct(Params *params) {
	if (params->phaseParams != NULL) stRPHmmParameters_destruct(params->phaseParams);
    if (params->polishParams != NULL) polishParams_destruct(params->polishParams);
    if (params->baseMapper != NULL) stBaseMapper_destruct(params->baseMapper);
	free(params);
}

void params_printParameters(Params *params, FILE *fh) {
	fprintf(fh, "Polish parameters:\n");
	polishParams_printParameters(params->polishParams, fh);
	fprintf(fh, "Phase parameters:\n");
	stRPHmmParameters_printParameters(params->phaseParams, fh);
}
