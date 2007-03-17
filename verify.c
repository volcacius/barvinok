#include <stdlib.h>
#include "verify.h"

/* RANGE : normal range for evalutations (-RANGE -> RANGE) */
#define RANGE 50

/* SRANGE : small range for evalutations */
#define SRANGE 15

/* if dimension >= BIDDIM, use SRANGE */
#define BIGDIM 5

/* VSRANGE : very small range for evalutations */
#define VSRANGE 5

/* if dimension >= VBIDDIM, use VSRANGE */
#define VBIGDIM 8

static struct argp_option argp_options[] = {
    { "verify",     	    'T',    0,	    0 },
    { "print-all",  	    'A',    0,	    0 },
    { "continue-on-error",  'C',    0,	    0 },
    { "min",   	    	    'm',    "int",  0 },
    { "max",   	    	    'M',    "int",  0 },
    { "range",      	    'r',    "int",  0 },
    { 0 }
};

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
    struct verify_options *options = state->input;

    switch (key) {
    case ARGP_KEY_INIT:
	options->verify = 0;
	options->print_all = 0;
	options->continue_on_error = 0;
	options->m = INT_MAX;
	options->M = INT_MIN;
	break;
    case ARGP_KEY_FINI:
	break;
    case 'T':
	options->verify = 1;
	break;
    case 'A':
	options->print_all = 1;
	break;
    case 'C':
	options->continue_on_error = 1;
	break;
    case 'm':
	options->m = atoi(arg);
	options->verify = 1;
	break;
    case 'M':
	options->M = atoi(arg);
	options->verify = 1;
	break;
    case 'r':
	options->M = atoi(arg);
	options->m = -options->M;
	options->verify = 1;
	break;
    default:
	return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

void verify_options_set_range(struct verify_options *options, Polyhedron *P)
{
    int r;

    if (P->Dimension >= VBIGDIM)
	r = VSRANGE;
    else if (P->Dimension >= BIGDIM)
	r = SRANGE;
    else
	r = RANGE;
    if (options->M == INT_MIN)
	options->M = r;
    if (options->m == INT_MAX)
	options->m = -r;

    if (options->verify && options->m > options->M) {
	fprintf(stderr,"Nothing to do: min > max !\n");
	exit(0);
    }
}

struct argp verify_argp = {
    argp_options, parse_opt, 0, 0
};
