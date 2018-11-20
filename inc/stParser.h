/*
 * Copyright (C) 2009-2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef ST_PARSER_H_
#define ST_PARSER_H_

#include "sonLib.h"
#include "stPolish.h"
#include "stRPHmm.h"

/*
 * Combined parameter object for phase, polish, view, etc.
 */

typedef struct _params {
	PolishParams *polishParams;
	stRPHmmParameters *phaseParams;
	stBaseMapper *baseMapper;
} Params;

Params *params_readParams(FILE *fp);

void params_destruct(Params *params);

void params_printParameters(Params *params, FILE *fh);

#endif /* ST_PARSER_H_ */
