/*
 * stateMachine.h
 *
 *  Created on: Aug 8, 2014
 *      Author: benedictpaten
 */

#ifndef STATEMACHINE_H_
#define STATEMACHINE_H_

#include "sonLib.h"

typedef struct _alphabet Alphabet;
typedef struct _emissions Emissions;
typedef uint16_t Symbol;
typedef struct _hmm Hmm;

typedef struct _symbolString {
	Alphabet *alphabet;
	Symbol *sequence;
	int64_t length;
} SymbolString;

SymbolString symbolString_construct(const char *sequence, int64_t length, Alphabet *a);

void symbolString_destruct(SymbolString s);

/*
 * Alphabet object
 */

struct _alphabet {
	uint64_t alphabetSize;

	Symbol (*convertCharToSymbol)(char i);

	char (*convertSymbolToChar)(Symbol i);
};

Alphabet *alphabet_constructNucleotide();

void alphabet_destruct(Alphabet *alphabet);

/*
 * Emissions object
 */

typedef enum {
    nucleotideEmissions=0,
	nucleotideEmissionsSymmetric=1,
	runlengthNucleotideEmissions=2
} EmissionType;

struct _emissions {
	Alphabet *alphabet;

    double (*emission)(Emissions *e, Symbol cX, Symbol cY);

    double (*gapEmissionX)(Emissions *e, Symbol cX);

    double (*gapEmissionY)(Emissions *e, Symbol cY);
};

Emissions *nucleotideEmissions_construct();

Emissions *emissions_construct(Hmm *hmm);

void emissions_destruct(Emissions *e);

/*
 * The statemachine object for computing pairwise alignments with.
 */

typedef enum {
    threeState=2,
    threeStateAsymmetric=3
} StateMachineType;

typedef struct _stateMachine StateMachine;

struct _stateMachine {
    StateMachineType type;
    int64_t stateNumber;
    int64_t matchState;
    int64_t gapXState; // This is the "primary" gap x state (e.g. short)
    int64_t gapYState;  // This is the "primary" gap y state (e.g. short)
    Emissions *emissions; // Used to calculate emissions probabilities

    double (*startStateProb)(StateMachine *sM, int64_t state);

    double (*endStateProb)(StateMachine *sM, int64_t state);

    double (*raggedEndStateProb)(StateMachine *sM, int64_t state);

    double (*raggedStartStateProb)(StateMachine *sM, int64_t state);

    //Cells (states at a given coordinate(
    void (*cellCalculate)(StateMachine *sM, double *current, double *lower, double *middle, double *upper, Symbol cX, Symbol cY,
            void(*doTransition)(double *, double *, int64_t, int64_t, double, double, void *), void *extraArgs);
};

StateMachine *hmm_getStateMachine(Hmm *hmm);

StateMachine *stateMachine3_construct(StateMachineType type, Emissions *e); //the type is to specify symmetric/asymmetric

StateMachine *stateMachine3_constructNucleotide(StateMachineType type); // Construct with default nucleotide emission model

void stateMachine_destruct(StateMachine *stateMachine);

/*
 * Hmm for loading/unloading HMMs and storing expectations.
 */

struct _hmm {
    StateMachineType type;
    EmissionType emissionsType;
    int64_t emissionNoPerState;
    double *transitions;
    double *emissions;
    double likelihood;
    int64_t stateNumber;
};

void hmm_destruct(Hmm *hmm);

Hmm *hmm_jsonParse(char *buf, size_t r);

Hmm *hmm_constructEmpty(double pseudoExpectation, StateMachineType type, EmissionType emissionType);

// Stuff related to EM - legacy right now

void hmm_normalise(Hmm *hmm);

void hmm_randomise(Hmm *hmm); //Creates normalised HMM with parameters set to small random values.

void hmm_addToTransitionExpectation(Hmm *hmmExpectations, int64_t from, int64_t to, double p);

double hmm_getTransition(Hmm *hmmExpectations, int64_t from, int64_t to);

void hmm_setTransition(Hmm *hmm, int64_t from, int64_t to, double p);

void hmm_addToEmissionsExpectation(Hmm *hmmExpectations, int64_t state, int64_t emissionNo, double p);

double hmm_getEmissionsExpectation(Hmm *hmm, int64_t state, int64_t emissionNo);

void hmm_setEmissionsExpectation(Hmm *hmm, int64_t state, int64_t emissionNo, double p);

#endif /* STATEMACHINE_H_ */
