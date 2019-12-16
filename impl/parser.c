/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com) & Arthur Rand (arand@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <sys/stat.h>
#include "margin.h"

/*
 * Get model parameters from params file.
 * Set hmm parameters.
*/

stRPHmmParameters *stRPHmmParameters_construct() {
	// Params object
	stRPHmmParameters *params = st_calloc(1, sizeof(stRPHmmParameters));

	// More variables for hmm stuff
	params->maxCoverageDepth = MAX_READ_PARTITIONING_DEPTH;
	params->maxNotSumTransitions = true;
	params->minPartitionsInAColumn = 50;
	params->maxPartitionsInAColumn = 200;
	params->minPosteriorProbabilityForPartition = 0.001;
	params->minReadCoverageToSupportPhasingBetweenHeterozygousSites = 0;

	// Other marginPhase program options
	params->roundsOfIterativeRefinement = 0;
	params->includeInvertedPartitions = true;

	return params;
}

stRPHmmParameters *parseParameters_fromJson(char *buf, size_t r) {
	// Setup parser
	jsmntok_t *tokens;
	char *js;
	int64_t tokenNumber = stJson_setupParser(buf, r, &tokens, &js);

	stRPHmmParameters *params = stRPHmmParameters_construct();

    //TODO: refactor the following to use the json parsing functions

    // Parse tokens, starting at token 1
    // (token 0 is entire object)
    for (int64_t i = 1; i < tokenNumber; i++) {
        jsmntok_t key = tokens[i];
        char *keyString = stJson_token_tostr(js, &key);

        if (strcmp(keyString, "maxNotSumTransitions") == 0) {
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
        else if (strcmp(keyString, "includeInvertedPartitions") == 0) {
        	params->includeInvertedPartitions = stJson_parseBool(js, tokens, ++i);
        }
        else if (strcmp(keyString, "verbose") == 0) {
            jsmntok_t tok = tokens[i+1];
            char *tokStr = stJson_token_tostr(js, &tok);
            int64_t bitString = atoi(tokStr);
            // TODO - currently does nothing
            i++;
        }
        else if (strcmp(keyString, "roundsOfIterativeRefinement") == 0) {
        	params->roundsOfIterativeRefinement = stJson_parseInt(js, tokens, ++i);
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
//void setVerbosity(stRPHmmParameters *params, int64_t bitstring) {
//    params->verboseTruePositives = (bitstring & LOG_TRUE_POSITIVES) != 0;
//    params->verboseFalsePositives = (bitstring & LOG_FALSE_POSITIVES) != 0;
//    params->verboseFalseNegatives = (bitstring & LOG_FALSE_NEGATIVES) != 0;
//}
//
//stRPHmmParameters *parseParameters(char *paramsFile) {
//	char buf[BUFSIZ * 3000]; // TODO: FIX, This is terrible code, we should not assume the size of the file is less than this
//	FILE *fh = fopen(paramsFile, "rb");
//	if (fh == NULL) {
//		st_errAbort("ERROR: Cannot open parameters file %s\n", paramsFile);
//	}
//	stRPHmmParameters *p =  parseParameters_fromJson(buf, fread(buf, sizeof(char), sizeof(buf), fh));
//	fclose(fh);
//	return p;
//}
//
///*
// * Counts the number of insertions and deletions in a read, given its cigar string.
// */
//void countIndels(uint32_t *cigar, uint32_t ncigar, int64_t *numInsertions, int64_t *numDeletions) {
//    for (uint32_t i = 0; i < ncigar; i++) {
//        int cigarOp = cigar[i] & BAM_CIGAR_MASK;
//        int cigarNum = cigar[i] >> BAM_CIGAR_SHIFT;
//        if (cigarOp == BAM_CINS) *numInsertions += cigarNum;
//        if (cigarOp == BAM_CDEL) *numDeletions += cigarNum;
//    }
//}
//
//int64_t getAlignedReadLength3(bam1_t *aln, int64_t *start_softclip, int64_t *end_softclip, bool boundaryAtMatch) {
//    // start read needs to be init'd to 0 (mostly this is to avoid misuse)
//    if (*start_softclip != 0) st_errAbort("getAlignedReadLength2 invoked with improper start_softclip parameter");
//    if (*end_softclip != 0) st_errAbort("getAlignedReadLength2 invoked with improper end_softclip parameter");
//
//    // get relevant cigar info
//    int64_t len = aln->core.l_qseq;
//    uint32_t *cigar = bam_get_cigar(aln);
//
//    // data for tracking
//    int64_t start_ref = 0;
//    int64_t cig_idx = 0;
//
//    // Find the correct starting locations on the read and reference sequence,
//    // to deal with things like inserts / deletions / soft clipping
//    while (cig_idx < aln->core.n_cigar) {
//        int cigarOp = cigar[cig_idx] & BAM_CIGAR_MASK;
//        int cigarNum = cigar[cig_idx] >> BAM_CIGAR_SHIFT;
//
//        if (cigarOp == BAM_CMATCH || cigarOp == BAM_CEQUAL || cigarOp == BAM_CDIFF) {
//            break;
//        } else if (cigarOp == BAM_CDEL || cigarOp == BAM_CREF_SKIP) {
//            if (!boundaryAtMatch) break;
//            cig_idx++;
//        } else if (cigarOp == BAM_CINS) {
//            if (!boundaryAtMatch) break;
//            *start_softclip += cigarNum;
//            cig_idx++;
//        } else if (cigarOp == BAM_CSOFT_CLIP) {
//            *start_softclip += cigarNum;
//            cig_idx++;
//        } else if (cigarOp == BAM_CHARD_CLIP || cigarOp == BAM_CPAD) {
//            cig_idx++;
//        } else {
//            st_errAbort("Unidentifiable cigar operation\n");
//        }
//    }
//
//    // Check for soft clipping at the end
//    cig_idx = aln->core.n_cigar - 1;
//    while (cig_idx > 0) {
//        int cigarOp = cigar[cig_idx] & BAM_CIGAR_MASK;
//        int cigarNum = cigar[cig_idx] >> BAM_CIGAR_SHIFT;
//
//        if (cigarOp == BAM_CMATCH || cigarOp == BAM_CEQUAL || cigarOp == BAM_CDIFF) {
//            break;
//        } else if (cigarOp == BAM_CDEL || cigarOp == BAM_CREF_SKIP) {
//            if (!boundaryAtMatch) break;
//            cig_idx--;
//        } else if (cigarOp == BAM_CINS) {
//            if (!boundaryAtMatch) break;
//            *end_softclip += cigarNum;
//            cig_idx--;
//        } else if (cigarOp == BAM_CSOFT_CLIP) {
//            *end_softclip += cigarNum;
//            cig_idx--;
//        } else if (cigarOp == BAM_CHARD_CLIP || cigarOp == BAM_CPAD) {
//            cig_idx--;
//        } else {
//            st_errAbort("Unidentifiable cigar operation\n");
//        }
//    }
//
//    // Count number of insertions & deletions in sequence
//    int64_t numInsertions = 0;
//    int64_t numDeletions = 0;
//    countIndels(cigar, aln->core.n_cigar, &numInsertions, &numDeletions);
//    int64_t trueLength = len - *start_softclip - *end_softclip + numDeletions - numInsertions;
//
//    return trueLength;
//>>>>>>> 38cbd8720b51472c90061b42658b6e5665bd1106
//}

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

	// TODO: Generalize to not be nucleotide only
	RepeatSubMatrix *repeatSubMatrix = repeatSubMatrix_constructEmpty(alphabet_constructNucleotide());
	// TODO: ADD support for Ns
	for(int64_t tokenIndex=1; tokenIndex < tokenNumber; tokenIndex++) {
		jsmntok_t key = tokens[tokenIndex];
		char *keyString = stJson_token_tostr(js, &key);


		if(stString_eq("baseLogRepeatCounts_AT", keyString)) {
			tokenIndex = stJson_parseFloatArray(repeatSubMatrix->baseLogProbs_AT,
					repeatSubMatrix->maximumRepeatLength, js, tokens, tokenIndex+1);
			continue;
		}

		if(stString_eq("baseLogRepeatCounts_GC", keyString)) {
			tokenIndex = stJson_parseFloatArray(repeatSubMatrix->baseLogProbs_GC,
					repeatSubMatrix->maximumRepeatLength, js, tokens, tokenIndex+1);
			continue;
		}

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
		tokenIndex = repeatSubMatrix_parseLogProbabilities(repeatSubMatrix,
				repeatSubMatrix->alphabet->convertCharToSymbol(base), strand, js, tokens, tokenIndex+1);
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
    params->includeSoftClipping = FALSE;
    params->shuffleChunks = TRUE;
    params->chunkSize = 0;
    params->chunkBoundary = 0;
    params->maxDepth = 0;
    params->includeSecondaryAlignments=FALSE;
    params->includeSupplementaryAlignments=TRUE;
    params->candidateVariantWeight = 0.2;
    params->columnAnchorTrim = 5;
    params->maxConsensusStrings = 100;
    params->repeatSubMatrix = NULL;

	// Parse tokens, starting at token 1
    // (token 0 is entire object)
	bool gotHmm = 0, gotHmmConditional = 0, gotPairwiseAlignmentParameters = 0, gotRepeatCountMatrix = 0, gotAlphabet = 0;
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
        else if (strcmp(keyString, "poaConstructCompareRepeatCounts") == 0) {
        	params->poaConstructCompareRepeatCounts = stJson_parseBool(js, tokens, ++tokenIndex);
        }
        else if (strcmp(keyString, "hmm") == 0) {
        	jsmntok_t tok = tokens[tokenIndex+1];
        	char *tokStr = stJson_token_tostr(js, &tok);
        	params->hmm = hmm_jsonParse(tokStr, strlen(tokStr));
        	params->sM = hmm_getStateMachine(params->hmm);
        	tokenIndex += stJson_getNestedTokenCount(tokens, tokenIndex+1);
        	gotHmm = 1;
        }
        else if (strcmp(keyString, "hmmConditional") == 0) {
			jsmntok_t tok = tokens[tokenIndex + 1];
			char *tokStr = stJson_token_tostr(js, &tok);
			params->hmmConditional = hmm_jsonParse(tokStr, strlen(tokStr));
			params->sMConditional = hmm_getStateMachine(params->hmmConditional);
			tokenIndex += stJson_getNestedTokenCount(tokens, tokenIndex + 1);
			gotHmmConditional = 1;
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
        else if (strcmp(keyString, "shuffleChunks") == 0) {
            params->shuffleChunks = stJson_parseBool(js, tokens, ++tokenIndex);
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
        else if (strcmp(keyString, "maxDepth") == 0) {
            if (stJson_parseInt(js, tokens, ++tokenIndex) < 0) {
                st_errAbort("ERROR: maxDepth parameter must zero or greater\n");
            }
            params->maxDepth = (uint64_t) stJson_parseInt(js, tokens, tokenIndex);
        }
        else if (strcmp(keyString, "includeSecondaryAlignments") == 0) {
            params->includeSecondaryAlignments = stJson_parseBool(js, tokens, ++tokenIndex);
        }
        else if (strcmp(keyString, "includeSupplementaryAlignments") == 0) {
            params->includeSupplementaryAlignments = stJson_parseBool(js, tokens, ++tokenIndex);
        }
        else if (strcmp(keyString, "candidateVariantWeight") == 0) {
			if (stJson_parseFloat(js, tokens, ++tokenIndex) < 0) {
				st_errAbort("ERROR: candidateVariantWeight parameter must zero or greater\n");
			}
			params->candidateVariantWeight = stJson_parseFloat(js, tokens, tokenIndex);
		}
        else if (strcmp(keyString, "columnAnchorTrim") == 0) {
			if (stJson_parseInt(js, tokens, ++tokenIndex) < 0) {
				st_errAbort("ERROR: columnAnchorTrim parameter must zero or greater\n");
			}
			params->columnAnchorTrim = (uint64_t) stJson_parseInt(js, tokens, tokenIndex);
		} else if (strcmp(keyString, "maxConsensusStrings") == 0) {
			if (stJson_parseInt(js, tokens, ++tokenIndex) < 0) {
				st_errAbort("ERROR: maxConsensusStrings parameter must zero or greater\n");
			}
			params->maxConsensusStrings = (uint64_t) stJson_parseInt(js, tokens, tokenIndex);
		} else if (strcmp(keyString, "maxPoaConsensusIterations") == 0) {
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
		} else if (strcmp(keyString, "minReadsToCallConsensus") == 0) {
			if (stJson_parseInt(js, tokens, ++tokenIndex) < 0) {
				st_errAbort("ERROR: minReadsToCallConsensus parameter must zero or greater\n");
			}
			params->minReadsToCallConsensus = (uint64_t) stJson_parseInt(js, tokens, tokenIndex);
		} else if (strcmp(keyString, "filterReadsWhileHaveAtLeastThisCoverage") == 0) {
			if (stJson_parseInt(js, tokens, ++tokenIndex) < 0) {
				st_errAbort("ERROR: filterReadsWhileHaveAtLeastThisCoverage parameter must zero or greater\n");
			}
			params->filterReadsWhileHaveAtLeastThisCoverage = (uint64_t) stJson_parseInt(js, tokens, tokenIndex);
		} else if (strcmp(keyString, "minAvgBaseQuality") == 0) {
			if (stJson_parseFloat(js, tokens, ++tokenIndex) < 0) {
				st_errAbort("ERROR: minAvgBaseQuality parameter must zero or greater\n");
			}
			params->minAvgBaseQuality = stJson_parseFloat(js, tokens, tokenIndex);
		} else if (strcmp(keyString, "hetScalingParameter") == 0) {
			if (stJson_parseFloat(js, tokens, ++tokenIndex) < 0) {
				st_errAbort("ERROR: hetScalingParameter parameter must zero or greater\n");
			}
			params->hetScalingParameter = stJson_parseFloat(js, tokens, tokenIndex);
		} else if (strcmp(keyString, "alleleStrandSkew") == 0) {
			if (stJson_parseFloat(js, tokens, ++tokenIndex) < 0) {
						st_errAbort("ERROR: alleleStrandSkew parameter must zero or greater\n");
			}
			params->alleleStrandSkew = stJson_parseFloat(js, tokens, tokenIndex);
		} else if (strcmp(keyString, "useOnlySubstitutionsForPhasing") == 0) {
		            params->useOnlySubstitutionsForPhasing = stJson_parseBool(js, tokens, ++tokenIndex);
		} else if (strcmp(keyString, "alphabet") == 0) {
			jsmntok_t tok = tokens[++tokenIndex];
			char *tokStr = stJson_token_tostr(js, &tok);
			if(stString_eq(tokStr, "nucleotide")) {
				params->alphabet = alphabet_constructNucleotide();
			}
			else {
				st_errAbort("ERROR: Unrecognised alphabet type json: %s\n", tokStr);
			}
			gotAlphabet = TRUE;
		} else {
            st_errAbort("ERROR: Unrecognised key in polish params json: %s\n", keyString);
        }
    }

    if(!gotRepeatCountMatrix && params->useRunLengthEncoding) {
    	st_logCritical("  ERROR: Did not find repeat counts specified in json polish params! Will default to MODE estimation\n");
    }
    if(!gotHmm) {
    	st_errAbort("ERROR: Did not find HMM specified in json polish params\n");
    }
    if(!gotHmmConditional) {
    	st_errAbort("ERROR: Did not find HMM conditional specified in json polish params\n");
    }
	if(!gotPairwiseAlignmentParameters) {
		st_errAbort("ERROR: Did not find pairwise alignment params specified in json polish params\n");
	}
	if(!gotAlphabet) {
		st_errAbort("ERROR: Did not find alphabet params specified in json polish params\n");
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

Hmm *polishParams_convertHmmEmissionsNuclToNuclRL(Hmm *hmmNucl, RepeatSubMatrix *rlMatrix) {
    assert(hmmNucl->emissionsType == nucleotideEmissions);
    Hmm *hmmNuclRL = hmm_constructEmpty(0, hmmNucl->type, runlengthNucleotideEmissions);




    return hmmNuclRL;
}

void polishParams_destruct(PolishParams *params) {
	if (params->repeatSubMatrix != NULL) repeatSubMatrix_destruct(params->repeatSubMatrix);
	stateMachine_destruct(params->sM);
    stateMachine_destruct(params->sMConditional);
    hmm_destruct(params->hmm);
    hmm_destruct(params->hmmConditional);
	pairwiseAlignmentBandingParameters_destruct(params->p);
    free(params->minPosteriorProbForAlignmentAnchors);
	alphabet_destruct(params->alphabet);
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
                params->phaseParams = parseParameters_fromJson(tokStr, strlen(tokStr));
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
                params->phaseParams = parseParameters_fromJson(buf, r);
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

Params *params_readParams(char *paramsFile) {
    return params_readParams2(paramsFile, TRUE, TRUE);
}
Params *params_readParams2(char *paramsFile, bool requirePolish, bool requirePhase) {
    // open file and check for existence
    FILE *fh = fopen(paramsFile, "rb");
    if (fh == NULL) {
        st_errAbort("ERROR: Cannot open parameters file %s\n", paramsFile);
    }

    // get file length
    struct stat st;
    fstat(fileno(fh), &st);

    // read file and parse json
    char *buf = st_calloc(st.st_size + 1, sizeof(char));
    size_t readLen = fread(buf, sizeof(char), st.st_size , fh);
    buf[st.st_size] = '\0';
    Params *params = params_jsonParse(buf, readLen, requirePolish, requirePhase);

    // close
    free(buf);
    fclose(fh);

    return params;
}

void params_destruct(Params *params) {
	if (params->phaseParams != NULL) stRPHmmParameters_destruct(params->phaseParams);
    if (params->polishParams != NULL) polishParams_destruct(params->polishParams);
	free(params);
}

void params_printParameters(Params *params, FILE *fh) {
	fprintf(fh, "Polish parameters:\n");
	polishParams_printParameters(params->polishParams, fh);
	fprintf(fh, "Phase parameters:\n");
	stRPHmmParameters_printParameters(params->phaseParams, fh);
}
