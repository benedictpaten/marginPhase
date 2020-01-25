/*
 * stateMachine.c
 *
 *  Created on: 1 Aug 2014
 *      Author: benedictpaten
 */

#include "margin.h"

#define RLENE_MAX_REPEAT_LENGTH 51

typedef enum {
    match = 0, shortGapX = 1, shortGapY = 2, longGapX = 3, longGapY = 4
} State;

static void state_check(StateMachine *sM, State s) {
    assert(s >= 0 && s < sM->stateNumber);
}

///////////////////////////////////
///////////////////////////////////
//Alphabet
//
///////////////////////////////////
///////////////////////////////////

static Symbol convertNucleotideCharToSymbol(char i) {
    switch (i) {
		case 'A':
		case 'a':
			return 0;
		case 'C':
		case 'c':
			return 1;
		case 'G':
		case 'g':
			return 2;
		case 'T':
		case 't':
			return 3;
		default:
			return 4;
    }
}

static char convertNucleotideSymbolToChar(Symbol i) {
    switch (i) {
		case 0:
			return 'A';
		case 1:
			return 'C';
		case 2:
			return 'G';
		case 3:
			return 'T';
		default:
			return 'N';
    }
}

bool nucleotidesEqual(Symbol x, Symbol y) {
	return x == y;
}

Alphabet *alphabet_constructNucleotide() {
	Alphabet *a = st_calloc(1, sizeof(Alphabet));
	a->alphabetSize = 5; // Fifth character represents "N"
	a->convertCharToSymbol = convertNucleotideCharToSymbol;
	a->convertSymbolToChar = convertNucleotideSymbolToChar;
	a->symbolsEqual = nucleotidesEqual;

	return a;
}

void alphabet_destruct(Alphabet *a) {
	free(a);
}

///////////////////////////////////
///////////////////////////////////
//SymbolStrings
//
///////////////////////////////////
///////////////////////////////////

Symbol *symbol_convertStringToSymbols(const char *s, int64_t start, int64_t length, Alphabet *a) {
    assert(length >= 0 && start >= 0);
    assert(strlen(s) >= start+length);
    Symbol *cS = st_malloc(length * sizeof(Symbol));
    for (int64_t i = 0; i < length; i++) {
        cS[i] = a->convertCharToSymbol(s[i+start]);
    }
    return cS;
}

SymbolString symbolString_construct(const char *sequence, int64_t start, int64_t length, Alphabet *a) {
    SymbolString symbolString;
    symbolString.alphabet = a;
    symbolString.sequence = symbol_convertStringToSymbols(sequence, start, length, a);
    symbolString.length = length;
    return symbolString;
}

SymbolString symbolString_getSubString(SymbolString s, uint64_t start, uint64_t length) {
	SymbolString s2 = s;
	s2.sequence = st_calloc(length, sizeof(Symbol));
	assert(start+length <= s.length); // Check bounds
	for(int64_t i=0; i<length; i++) {
		s2.sequence[i] = s.sequence[start+i];
	}
	return s2;
}

void symbolString_destruct(SymbolString s) {
    free(s.sequence);
}

Symbol symbol_getRepeatLength(Symbol s) {
	return s >> 8;  // Last 13 bits are length
}

Symbol symbol_stripRepeatCount(Symbol s) {
	return s & 255; // First eight bits encode symbol
}

Symbol symbol_addRepeatCount(Symbol character, uint64_t runLength) {
	assert(character <= 255);
	assert(runLength <= 255);
	return (runLength << 8) | character;
}

SymbolString rleString_constructSymbolString(RleString *s, int64_t start, int64_t length, Alphabet *a) {
	assert(start+length <= s->length);
	SymbolString symbolString;
	symbolString.alphabet = a;
	symbolString.sequence = symbol_convertStringToSymbols(s->rleString, start, length, a);
	for (int64_t i = 0; i < length; i++) {
		symbolString.sequence[i] = symbol_addRepeatCount(symbolString.sequence[i], s->repeatCounts[i+start]);
	}
	symbolString.length = length;
	return symbolString;
}

///////////////////////////////////
///////////////////////////////////
//Hmms
///////////////////////////////////
///////////////////////////////////

Hmm *hmm_constructEmpty(double pseudoExpectation, StateMachineType type, EmissionType emissionType) {
    Hmm *hmm = st_calloc(1, sizeof(Hmm));
    hmm->type = type;
    switch (type) {
    case threeState:
    case threeStateAsymmetric:
        hmm->stateNumber = 3;
        break;
    default:
        st_errAbort("Unrecognised state type: %i\n", type);
    }
    hmm->transitions = st_malloc(hmm->stateNumber * hmm->stateNumber * sizeof(double));
    for (int64_t i = 0; i < hmm->stateNumber * hmm->stateNumber; i++) {
        hmm->transitions[i] = pseudoExpectation;
    }

    hmm->emissionsType = emissionType;
    hmm->emissionNoPerState = st_calloc(hmm->stateNumber, sizeof(int64_t));
    hmm->emissionOffsetPerState = st_calloc(hmm->stateNumber, sizeof(int64_t));

    switch (emissionType) {
    case nucleotideEmissions:
    	for(int64_t i=0; i<hmm->stateNumber; i++) {
    		hmm->emissionNoPerState[i] = 16;
    	}
    	break;
    case runlengthNucleotideEmissions:
    	for(int64_t i=0; i<hmm->stateNumber; i++) {
    		hmm->emissionNoPerState[i] = 4 + RLENE_MAX_REPEAT_LENGTH;
    	}
    	hmm->emissionNoPerState[match] = 16 + RLENE_MAX_REPEAT_LENGTH * RLENE_MAX_REPEAT_LENGTH;
		break;
    default:
         st_errAbort("Unrecognised emission type: %i\n", type);
    }
    int64_t j=0;
    for(int64_t i=0; i<hmm->stateNumber; i++) {
    	hmm->emissionOffsetPerState[i] = j;
        j += hmm->emissionNoPerState[i];
    }
    hmm->totalEmissions = j;
    hmm->emissions = st_malloc(hmm->totalEmissions * sizeof(double));
    for (int64_t i = 0; i < hmm->totalEmissions; i++) {
        hmm->emissions[i] = pseudoExpectation;
    }
    hmm->likelihood = 0.0;
    return hmm;
}

void hmm_destruct(Hmm *hmm) {
    free(hmm->transitions);
    free(hmm->emissions);
    free(hmm->emissionNoPerState);
    free(hmm->emissionOffsetPerState);
    free(hmm);
}

Hmm *hmm_jsonParse(char *buf, size_t r) {
	// Setup parser
	jsmntok_t *tokens;
	char *js;
	int64_t tokenNumber = stJson_setupParser(buf, r, &tokens, &js);
	if(tokenNumber < 4) {
		st_errAbort("ERROR: too few tokens to parse in hmm json: %i\n", (int)tokenNumber);
	}

	// Find out type of hmm
	char *keyString = stJson_token_tostr(js, &(tokens[1]));
	if(strcmp(keyString, "type") != 0) {
		st_errAbort("ERROR: Unrecognised or unordered ('type', 'emissionsType', 'transitions', 'emissions') key in polish params json: %s\n", keyString);
	}
	int64_t stateMachineType = stJson_parseInt(js, tokens, 2);

	// Find out emission type of hmm
	keyString = stJson_token_tostr(js, &(tokens[3]));
	if(strcmp(keyString, "emissionsType") != 0) {
		st_errAbort("ERROR: Unrecognised or unordered ('type', 'emissionsType', 'transitions', 'emissions') key in polish params json: %s\n", keyString);
	}
	int64_t emissionsType = stJson_parseInt(js, tokens, 4);

	// Make empty hmm object
	Hmm *hmm = hmm_constructEmpty(0, stateMachineType, emissionsType);

	// Parse tokens, starting at token 1
    // (token 0 is entire object)
	bool gotEmissions = 0, gotTransitions = 0;
    for (int64_t tokenIndex=5; tokenIndex < tokenNumber; tokenIndex++) {
        keyString = stJson_token_tostr(js, &(tokens[tokenIndex]));
        if(strcmp(keyString, "transitions") == 0) {
        	tokenIndex = stJson_parseFloatArray(hmm->transitions, hmm->stateNumber * hmm->stateNumber, js, tokens, ++tokenIndex);
        	gotTransitions = 1;
        }
        else if(strcmp(keyString, "emissions") == 0) {
        	tokenIndex = stJson_parseFloatArray(hmm->emissions, hmm->totalEmissions,
        			js, tokens, ++tokenIndex);
        	gotEmissions = 1;
        }
        else if(strcmp(keyString, "likelihood") == 0) {
        	hmm->likelihood = stJson_parseFloat(js, tokens, ++tokenIndex);
        }
        else {
        	st_errAbort("ERROR: Unrecognised key in hmm json: %s\n", keyString);
        }
    }

    if(!gotEmissions) {
    	st_errAbort("ERROR: Did not find emissions specified in json HMM\n");
    }
    if(!gotTransitions) {
        st_errAbort("ERROR: Did not find transitions specified in json HMM\n");
    }

    // Cleanup
    free(js);
    free(tokens);

    return hmm;
}

// Transitions

static inline double *hmm_getTransition2(Hmm *hmm, int64_t from, int64_t to) {
    return &(hmm->transitions[from * hmm->stateNumber + to]);
}

double hmm_getTransition(Hmm *hmm, int64_t from, int64_t to) {
    return *hmm_getTransition2(hmm, from, to);
}

void hmm_addToTransitionExpectation(Hmm *hmm, int64_t from, int64_t to, double p) {
    *hmm_getTransition2(hmm, from, to) += p;
}

void hmm_setTransition(Hmm *hmm, int64_t from, int64_t to, double p) {
    *hmm_getTransition2(hmm, from, to) = p;
}

// Emissions

static inline double *hmm_getEmissionsExpectation2(Hmm *hmm, int64_t state, int64_t emissionNo) {
    assert(state >= 0 && state < hmm->stateNumber);
    assert(emissionNo >= 0 && emissionNo < hmm->emissionNoPerState[state]);
    assert(hmm->emissionOffsetPerState[state] + emissionNo < hmm->totalEmissions);
	return &(hmm->emissions[hmm->emissionOffsetPerState[state] + emissionNo]);
}

double hmm_getEmissionsExpectation(Hmm *hmm, int64_t state, int64_t emissionNo) {
    return *hmm_getEmissionsExpectation2(hmm, state, emissionNo);
}

void hmm_addToEmissionsExpectation(Hmm *hmm, int64_t state, int64_t emissionNo, double p) {
    *hmm_getEmissionsExpectation2(hmm, state, emissionNo) += p;
}

void hmm_setEmissionsExpectation(Hmm *hmm, int64_t state, int64_t emissionNo, double p) {
    *hmm_getEmissionsExpectation2(hmm, state, emissionNo) = p;
}

static void hmm_emissions_loadProbs(Hmm *hmm, double *emissionProbs, int64_t state, int64_t start, int64_t emissionNumber) {
    for(int64_t x=0; x<emissionNumber; x++) {
    	emissionProbs[x] = log(hmm_getEmissionsExpectation(hmm, state, start+x));
    }
}

// Em stuff

void hmm_normalise(Hmm *hmm) {
    for (int64_t from = 0; from < hmm->stateNumber; from++) {
        double total = 0.0;
        for (int64_t to = 0; to < hmm->stateNumber; to++) {
            total += hmm_getTransition(hmm, from, to);
        }
        for (int64_t to = 0; to < hmm->stateNumber; to++) {
            hmm_setTransition(hmm, from, to, hmm_getTransition(hmm, from, to) / total);
        }
    }
    //Normalise the emissions
    for (int64_t state = 0; state < hmm->stateNumber; state++) {
        double total = 0.0;
        for (int64_t x = 0; x < hmm->emissionNoPerState[state]; x++) {
        	total += hmm_getEmissionsExpectation(hmm, state, x);
        }
        for (int64_t x = 0; x < hmm->emissionNoPerState[state]; x++) {
        	hmm_setEmissionsExpectation(hmm, state, x, hmm_getEmissionsExpectation(hmm, state, x) / total);
        }
    }
}

void hmm_randomise(Hmm *hmm) {
    //Transitions
    for (int64_t from = 0; from < hmm->stateNumber; from++) {
        for (int64_t to = 0; to < hmm->stateNumber; to++) {
            hmm_setTransition(hmm, from, to, st_random());
        }
    }
    //Emissions
    for (int64_t state = 0; state < hmm->stateNumber; state++) {
        for (int64_t x = 0; x < hmm->emissionNoPerState[state]; x++) {
        	hmm_setEmissionsExpectation(hmm, state, x, st_random());
        }
    }

    hmm_normalise(hmm);
}

///////////////////////////////////
///////////////////////////////////
//Emissions
///////////////////////////////////
///////////////////////////////////

typedef struct _nucleotideEmissions {
	Emissions e;
	double EMISSION_MATCH_PROBS[16]; //Match emission probs
	double EMISSION_GAP_X_PROBS[4]; //Gap X emission probs
	double EMISSION_GAP_Y_PROBS[4]; //Gap Y emission probs
} NucleotideEmissions;

static inline double getNucleotideGapProb(const double *emissionGapProbs, Symbol i) {
    if(i >= 4) {
        return -1.386294361; // log(0.25) ambiguous character
    }
    return emissionGapProbs[i];
}

static inline double getNucleotideGapProbX(NucleotideEmissions *ne, Symbol i) {
    return getNucleotideGapProb(ne->EMISSION_GAP_X_PROBS, i);
}

static inline double getNucleotideGapProbY(NucleotideEmissions *ne, Symbol i) {
    return getNucleotideGapProb(ne->EMISSION_GAP_Y_PROBS, i);
}

static inline double getNucleotideMatchProb(NucleotideEmissions *e, Symbol x, Symbol y) {
    if(x >= 4 || y >= 4) {
        return -2.772588722; //log(0.25**2)
    }
    return e->EMISSION_MATCH_PROBS[x * 4 + y];
}

static void setNucleotideEmissionMatchProbsToDefaults(double *emissionMatchProbs) {
    /*
     * This is used to set the emissions to reasonable values.
     */
    const double EMISSION_MATCH=-1.8917761142; //ln(0.1508037262);
    const double EMISSION_TRANSVERSION=-4.3459578861; //ln(0.01295908897);
    const double EMISSION_TRANSITION=-3.760242452; //ln(0.02327809587);
    //Symmetric matrix of transition probabilities.
    const double i[16] = {
            EMISSION_MATCH, EMISSION_TRANSVERSION, EMISSION_TRANSITION, EMISSION_TRANSVERSION,
            EMISSION_TRANSVERSION, EMISSION_MATCH, EMISSION_TRANSVERSION, EMISSION_TRANSITION,
            EMISSION_TRANSITION, EMISSION_TRANSVERSION, EMISSION_MATCH, EMISSION_TRANSVERSION,
            EMISSION_TRANSVERSION, EMISSION_TRANSITION, EMISSION_TRANSVERSION, EMISSION_MATCH };
    memcpy(emissionMatchProbs, i, sizeof(double)*16);
}

static void setNucleotideEmissionGapProbsToDefaults(double *emissionGapProbs) {
    /*
     * This is used to set the emissions to reasonable values.
     */
    const double EMISSION_GAP = -1.3862943611; //log(0.25)
    const double i[4] = { EMISSION_GAP, EMISSION_GAP, EMISSION_GAP, EMISSION_GAP };
    memcpy(emissionGapProbs, i, sizeof(double)*4);
}

Emissions *nucleotideEmissions_construct() {
	NucleotideEmissions *ne = st_calloc(1, sizeof(NucleotideEmissions));

	ne->e.alphabet = alphabet_constructNucleotide();
	ne->e.emission = (double (*)(Emissions *, Symbol, Symbol)) getNucleotideMatchProb;
	ne->e.gapEmissionX = (double (*)(Emissions *, Symbol)) getNucleotideGapProbX;
	ne->e.gapEmissionY = (double (*)(Emissions *, Symbol)) getNucleotideGapProbY;

	setNucleotideEmissionMatchProbsToDefaults(ne->EMISSION_MATCH_PROBS);
	setNucleotideEmissionGapProbsToDefaults(ne->EMISSION_GAP_X_PROBS);
	setNucleotideEmissionGapProbsToDefaults(ne->EMISSION_GAP_Y_PROBS);

	return (Emissions *)ne;
}

///////////////////////////////////
///////////////////////////////////
//RLE Emissions
///////////////////////////////////
///////////////////////////////////

typedef struct _rleNucleotideEmissions {
	NucleotideEmissions ne;
	double repeatLengthSubstitutionProbs[RLENE_MAX_REPEAT_LENGTH*RLENE_MAX_REPEAT_LENGTH];
	double repeatLengthGapXProbs[RLENE_MAX_REPEAT_LENGTH];
	double repeatLengthGapYProbs[RLENE_MAX_REPEAT_LENGTH];
	RepeatSubMatrix *repeatSubMatrix;
} RleNucleotideEmissions;

static inline double getRLENucleotideGapProb(const double *nucleotideGapProbs,
		const double *repeatLengthGapProbs, Symbol i) {
    return getNucleotideGapProb(nucleotideGapProbs, symbol_stripRepeatCount(i)); // + 1.0 * repeatLengthGapProbs[symbol_getRepeatLength(i)];
}

static inline double getRleNucleotideGapProbX(RleNucleotideEmissions *rlene, Symbol x) {
    return getRLENucleotideGapProb(rlene->ne.EMISSION_GAP_X_PROBS, rlene->repeatLengthGapXProbs, x);
}

static inline double getRleNucleotideGapProbY(RleNucleotideEmissions *rlene, Symbol y) {
    return getRLENucleotideGapProb(rlene->ne.EMISSION_GAP_Y_PROBS, rlene->repeatLengthGapYProbs, y);// + 2.3025 * rlene->repeatSubMatrix->baseLogProbs_AT[symbol_getRepeatLength(y)];
}

static inline double getRleNucleotideMatchProb(RleNucleotideEmissions *rlene, Symbol x, Symbol y) {
	Symbol xBase = symbol_stripRepeatCount(x);
    return getNucleotideMatchProb(&(rlene->ne), xBase, symbol_stripRepeatCount(y)) +
    		2.3025 * repeatSubMatrix_getLogProb(rlene->repeatSubMatrix, xBase, 1, symbol_getRepeatLength(y), symbol_getRepeatLength(x));
    		//0.0;
    		//(symbol_getRepeatLength(x) == symbol_getRepeatLength(y) ? 0.0 : -4.0);
    		//-1.0 * abs(symbol_getRepeatLength(x) - symbol_getRepeatLength(y));
    		//1.0 * rlene->repeatLengthSubstitutionProbs[symbol_getRepeatLength(x) * RLENE_MAX_REPEAT_LENGTH + symbol_getRepeatLength(y)];
}

Emissions *rleNucleotideEmissions_construct() {
	RleNucleotideEmissions *rlene = st_calloc(1, sizeof(RleNucleotideEmissions));
	NucleotideEmissions *ne = (NucleotideEmissions *)nucleotideEmissions_construct();
	rlene->ne = *ne;
	free(ne);
	rlene->ne.e.emission = (double (*)(Emissions *, Symbol, Symbol)) getRleNucleotideMatchProb;
	rlene->ne.e.gapEmissionX = (double (*)(Emissions *, Symbol)) getRleNucleotideGapProbX;
	rlene->ne.e.gapEmissionY = (double (*)(Emissions *, Symbol)) getRleNucleotideGapProbY;

	// Default to making all probs 1.
	for(int64_t i=0; i<RLENE_MAX_REPEAT_LENGTH; i++) {
		rlene->repeatLengthGapXProbs[i] = LOG_ONE;
		rlene->repeatLengthGapYProbs[i] = LOG_ONE;
		for(int64_t j=0; j<RLENE_MAX_REPEAT_LENGTH; j++) {
			rlene->repeatLengthSubstitutionProbs[i*RLENE_MAX_REPEAT_LENGTH+j] = LOG_ONE;
		}
	}

	return (Emissions *)rlene;
}

///////////////////////////////////
///////////////////////////////////
//Emissions constructor
///////////////////////////////////
///////////////////////////////////

Emissions *emissions_construct(Hmm *hmm) {
	if (hmm->emissionsType == nucleotideEmissions) {
		NucleotideEmissions *ne = (NucleotideEmissions *)nucleotideEmissions_construct();
		hmm_emissions_loadProbs(hmm, ne->EMISSION_MATCH_PROBS, 0, 0, 16);
		hmm_emissions_loadProbs(hmm, ne->EMISSION_GAP_X_PROBS, 1, 0, 4);
		hmm_emissions_loadProbs(hmm, ne->EMISSION_GAP_Y_PROBS, 2, 0, 4);
		return (Emissions *)ne;
	}
	else if(hmm->emissionsType == runlengthNucleotideEmissions) {
		RleNucleotideEmissions *rlene = (RleNucleotideEmissions *)rleNucleotideEmissions_construct();

		hmm_emissions_loadProbs(hmm, rlene->ne.EMISSION_MATCH_PROBS, 0, 0, 16);
		hmm_emissions_loadProbs(hmm, rlene->repeatLengthSubstitutionProbs, 0, 16,
				RLENE_MAX_REPEAT_LENGTH*RLENE_MAX_REPEAT_LENGTH);

		hmm_emissions_loadProbs(hmm, rlene->ne.EMISSION_GAP_X_PROBS, 1, 0, 4);
		hmm_emissions_loadProbs(hmm, rlene->repeatLengthGapXProbs, 1, 4,
				RLENE_MAX_REPEAT_LENGTH);

		hmm_emissions_loadProbs(hmm, rlene->ne.EMISSION_GAP_Y_PROBS, 2, 0, 4);
		hmm_emissions_loadProbs(hmm, rlene->repeatLengthGapYProbs, 2, 4,
				RLENE_MAX_REPEAT_LENGTH);

		return (Emissions *)rlene;
	}
	st_errAbort("Load from hmm: unrecognized emission type");
	return NULL;
}

void emissions_destruct(Emissions *e) {
	alphabet_destruct(e->alphabet);
	free(e);
}

///////////////////////////////////
///////////////////////////////////
//Three state state-machine
///////////////////////////////////
///////////////////////////////////

//Transitions
typedef struct _StateMachine3 StateMachine3;

struct _StateMachine3 {
    //3 state state machine, allowing for symmetry in x and y.
    StateMachine model;
    double TRANSITION_MATCH_CONTINUE;
    double TRANSITION_MATCH_FROM_GAP_X;
    double TRANSITION_MATCH_FROM_GAP_Y;
    double TRANSITION_GAP_OPEN_X;
    double TRANSITION_GAP_OPEN_Y;
    double TRANSITION_GAP_EXTEND_X;
    double TRANSITION_GAP_EXTEND_Y;
    double TRANSITION_GAP_SWITCH_TO_X;
    double TRANSITION_GAP_SWITCH_TO_Y;
};

static double stateMachine3_startStateProb(StateMachine *sM, int64_t state) {
    //Match state is like going to a match.
    state_check(sM, state);
    return state == match ? 0 : LOG_ZERO;
}

static double stateMachine3_raggedStartStateProb(StateMachine *sM, int64_t state) {
    state_check(sM, state);
    return (state == shortGapX || state == shortGapY) ? 0 : LOG_ZERO;
}

static double stateMachine3_endStateProb(StateMachine *sM, int64_t state) {
    //End state is like to going to a match
    StateMachine3 *sM3 = (StateMachine3 *) sM;
    state_check(sM, state);
    switch (state) {
    case match:
        return sM3->TRANSITION_MATCH_CONTINUE;
    case shortGapX:
        return sM3->TRANSITION_MATCH_FROM_GAP_X;
    case shortGapY:
        return sM3->TRANSITION_MATCH_FROM_GAP_Y;
    }
    return 0.0;
}

static double stateMachine3_raggedEndStateProb(StateMachine *sM, int64_t state) {
    //End state is like to going to a gap
    StateMachine3 *sM3 = (StateMachine3 *) sM;
    state_check(sM, state);
    switch (state) {
    case match:
        return (sM3->TRANSITION_GAP_OPEN_X + sM3->TRANSITION_GAP_OPEN_Y) / 2.0;
    case shortGapX:
        return sM3->TRANSITION_GAP_EXTEND_X;
    case shortGapY:
        return sM3->TRANSITION_GAP_EXTEND_Y;
    }
    return 0.0;
}

static void stateMachine3_cellCalculate(StateMachine *sM, double *current, double *lower, double *middle, double *upper,
        Symbol cX, Symbol cY, void (*doTransition)(double *, double *, int64_t, int64_t, double, double, void *),
        void *extraArgs) {
    StateMachine3 *sM3 = (StateMachine3 *) sM;
    if (lower != NULL) {
        double eP = sM->emissions->gapEmissionX(sM->emissions, cX);
        doTransition(lower, current, match, shortGapX, eP, sM3->TRANSITION_GAP_OPEN_X, extraArgs);
        doTransition(lower, current, shortGapX, shortGapX, eP, sM3->TRANSITION_GAP_EXTEND_X, extraArgs);
        doTransition(lower, current, shortGapY, shortGapX, eP, sM3->TRANSITION_GAP_SWITCH_TO_X, extraArgs);
    }
    if (middle != NULL) {
        double eP = sM->emissions->emission(sM->emissions, cX, cY);
        doTransition(middle, current, match, match, eP, sM3->TRANSITION_MATCH_CONTINUE, extraArgs);
        doTransition(middle, current, shortGapX, match, eP, sM3->TRANSITION_MATCH_FROM_GAP_X, extraArgs);
        doTransition(middle, current, shortGapY, match, eP, sM3->TRANSITION_MATCH_FROM_GAP_Y, extraArgs);
    }
    if (upper != NULL) {
        double eP = sM->emissions->gapEmissionY(sM->emissions, cY);
        doTransition(upper, current, match, shortGapY, eP, sM3->TRANSITION_GAP_OPEN_Y, extraArgs);
        doTransition(upper, current, shortGapY, shortGapY, eP, sM3->TRANSITION_GAP_EXTEND_Y, extraArgs);
        doTransition(upper, current, shortGapX, shortGapY, eP, sM3->TRANSITION_GAP_SWITCH_TO_Y, extraArgs);
    }
}

StateMachine *stateMachine3_construct(StateMachineType type, Emissions *e) {
    StateMachine3 *sM3 = st_malloc(sizeof(StateMachine3));
    sM3->TRANSITION_MATCH_CONTINUE = -0.030064059121770816; //0.9703833696510062f
    sM3->TRANSITION_MATCH_FROM_GAP_X = -1.272871422049609; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    sM3->TRANSITION_MATCH_FROM_GAP_Y = -1.272871422049609; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    sM3->TRANSITION_GAP_OPEN_X = -4.21256642; //0.0129868352330243
    sM3->TRANSITION_GAP_OPEN_Y = -4.21256642; //0.0129868352330243
    sM3->TRANSITION_GAP_EXTEND_X = -0.3388262689231553; //0.7126062401851738f;
    sM3->TRANSITION_GAP_EXTEND_Y = -0.3388262689231553; //0.7126062401851738f;
    sM3->TRANSITION_GAP_SWITCH_TO_X = -4.910694825551255; //0.0073673675173412815f;
    sM3->TRANSITION_GAP_SWITCH_TO_Y = -4.910694825551255; //0.0073673675173412815f;
    if (type != threeState && type != threeStateAsymmetric) {
        st_errAbort("Tried to create a three state state-machine with the wrong type");
    }
    sM3->model.type = type;
    sM3->model.stateNumber = 3;
    sM3->model.matchState = match;
    sM3->model.gapXState = shortGapX;
    sM3->model.gapYState = shortGapY;
    sM3->model.startStateProb = stateMachine3_startStateProb;
    sM3->model.endStateProb = stateMachine3_endStateProb;
    sM3->model.raggedStartStateProb = stateMachine3_raggedStartStateProb;
    sM3->model.raggedEndStateProb = stateMachine3_raggedEndStateProb;
    sM3->model.cellCalculate = stateMachine3_cellCalculate;
    sM3->model.emissions = e;

    return (StateMachine *) sM3;
}

StateMachine *stateMachine3_constructNucleotide(StateMachineType type) {
	return stateMachine3_construct(type, nucleotideEmissions_construct());
}

static void stateMachine3_loadAsymmetric(StateMachine3 *sM3, Hmm *hmm) {
    if (hmm->type != threeStateAsymmetric) {
        st_errAbort("Wrong hmm type");
    }
    sM3->TRANSITION_MATCH_CONTINUE = log(hmm_getTransition(hmm, match, match));
    sM3->TRANSITION_MATCH_FROM_GAP_X = log(hmm_getTransition(hmm, shortGapX, match));
    sM3->TRANSITION_MATCH_FROM_GAP_Y = log(hmm_getTransition(hmm, shortGapY, match));
    sM3->TRANSITION_GAP_OPEN_X = log(hmm_getTransition(hmm, match, shortGapX));
    sM3->TRANSITION_GAP_OPEN_Y = log(hmm_getTransition(hmm, match, shortGapY));
    sM3->TRANSITION_GAP_EXTEND_X = log(hmm_getTransition(hmm, shortGapX, shortGapX));
    sM3->TRANSITION_GAP_EXTEND_Y = log(hmm_getTransition(hmm, shortGapY, shortGapY));
    sM3->TRANSITION_GAP_SWITCH_TO_X = log(hmm_getTransition(hmm, shortGapY, shortGapX));
    sM3->TRANSITION_GAP_SWITCH_TO_Y = log(hmm_getTransition(hmm, shortGapX, shortGapY));
    int64_t xGapStates[1] = { shortGapX };
    int64_t yGapStates[1] = { shortGapY };
}

static void stateMachine3_loadSymmetric(StateMachine3 *sM3, Hmm *hmm) {
    if (hmm->type != threeState) {
        st_errAbort("Wrong hmm type");
    }
    sM3->TRANSITION_MATCH_CONTINUE = log(hmm_getTransition(hmm, match, match));
    sM3->TRANSITION_MATCH_FROM_GAP_X = log(
            (hmm_getTransition(hmm, shortGapX, match) + hmm_getTransition(hmm, shortGapY, match)) / 2.0);
    sM3->TRANSITION_MATCH_FROM_GAP_Y = sM3->TRANSITION_MATCH_FROM_GAP_X;
    sM3->TRANSITION_GAP_OPEN_X = log(
            (hmm_getTransition(hmm, match, shortGapX) + hmm_getTransition(hmm, match, shortGapY)) / 2.0);
    sM3->TRANSITION_GAP_OPEN_Y = sM3->TRANSITION_GAP_OPEN_X;
    sM3->TRANSITION_GAP_EXTEND_X = log(
            (hmm_getTransition(hmm, shortGapX, shortGapX) + hmm_getTransition(hmm, shortGapY, shortGapY)) / 2.0);
    sM3->TRANSITION_GAP_EXTEND_Y = sM3->TRANSITION_GAP_EXTEND_X;
    sM3->TRANSITION_GAP_SWITCH_TO_X = log(
            (hmm_getTransition(hmm, shortGapY, shortGapX) + hmm_getTransition(hmm, shortGapX, shortGapY)) / 2.0);
    sM3->TRANSITION_GAP_SWITCH_TO_Y = sM3->TRANSITION_GAP_SWITCH_TO_X;
    int64_t xGapStates[2] = { shortGapX };
    int64_t yGapStates[2] = { shortGapY };
}

///////////////////////////////////
///////////////////////////////////
//Public functions
///////////////////////////////////
///////////////////////////////////

StateMachine *hmm_getStateMachine(Hmm *hmm) {
	Emissions *e = emissions_construct(hmm);
    if (hmm->type == threeStateAsymmetric) {
        StateMachine3 *sM3 = (StateMachine3 *) stateMachine3_construct(hmm->type, e);
        stateMachine3_loadAsymmetric(sM3, hmm);
        return (StateMachine *) sM3;
    }
    if (hmm->type == threeState) {
        StateMachine3 *sM3 = (StateMachine3 *) stateMachine3_construct(hmm->type, e);
        stateMachine3_loadSymmetric(sM3, hmm);
        return (StateMachine *) sM3;
    }
    return NULL;
}

void stateMachine_addRepeatSubMatrix(StateMachine *stateMachine, RepeatSubMatrix *repeatSubMatrix) {
	RleNucleotideEmissions *rleEmissions = (RleNucleotideEmissions *)stateMachine->emissions;
	rleEmissions->repeatSubMatrix = repeatSubMatrix;
}

void stateMachine_destruct(StateMachine *stateMachine) {
	emissions_destruct(stateMachine->emissions);
    free(stateMachine);
}
