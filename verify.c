#include <stdlib.h>
#include <barvinok/options.h>
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
    { "exact",     	    'E',    0,	    0 },
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
	options->exact = 0;
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
    case 'E':
	options->exact = 1;
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

void verify_options_set_range(struct verify_options *options, int dim)
{
    int r;

    if (dim >= VBIGDIM)
	r = VSRANGE;
    else if (dim >= BIGDIM)
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

Polyhedron *check_poly_context_scan(Polyhedron *C,
				    const struct verify_options *options)
{
    int i;
    Matrix *MM;
    Polyhedron *CC, *CS, *U;
    unsigned MaxRays = options->barvinok->MaxRays;

    if (C->Dimension <= 0)
	return NULL;

    /* Intersect context with range */
    MM = Matrix_Alloc(2*C->Dimension, C->Dimension+2);
    for (i = 0; i < C->Dimension; ++i) {
	value_set_si(MM->p[2*i][0], 1);
	value_set_si(MM->p[2*i][1+i], 1);
	value_set_si(MM->p[2*i][1+C->Dimension], -options->m);
	value_set_si(MM->p[2*i+1][0], 1);
	value_set_si(MM->p[2*i+1][1+i], -1);
	value_set_si(MM->p[2*i+1][1+C->Dimension], options->M);
    }
    CC = AddConstraints(MM->p[0], 2*C->Dimension, C, options->barvinok->MaxRays);
    U = Universe_Polyhedron(0);
    CS = Polyhedron_Scan(CC, U, MaxRays & POL_NO_DUAL ? 0 : MaxRays);
    Polyhedron_Free(U);
    Polyhedron_Free(CC);
    Matrix_Free(MM);
    return CS;
}

void check_poly_init(Polyhedron *C, struct verify_options *options)
{
    int d, i;

    if (options->print_all)
	return;
    if (C->Dimension <= 0)
	return;

    d = options->M - options->m;
    if (d > 80)
	options->st = 1+d/80;
    else
	options->st = 1;
    for (i = options->m; i <= options->M; i += options->st)
	printf(".");
    printf( "\r" );
    fflush(stdout);
}

/****************************************************/
/* function check_poly :                            */
/* scans the parameter space from Min to Max (all   */
/* directions). Computes the number of points in    */
/* the polytope using both methods, and compare them*/
/* returns 1 on success                             */
/****************************************************/

int check_poly(Polyhedron *CS, const struct check_poly_data *data,
	       int nparam, int pos, Value *z,
	       const struct verify_options *options)
{
    if (pos == nparam) {
	if (!data->check(data, nparam, z, options) && !options->continue_on_error)
	    return 0;
    } else {
	Value LB, UB;
	int ok;
	value_init(LB);
	value_init(UB);
	ok = !(lower_upper_bounds(1+pos, CS, z-1, &LB, &UB));
	assert(ok);
	for (; value_le(LB, UB); value_increment(LB, LB)) {
	    if (!options->print_all) {
		int k = VALUE_TO_INT(LB);
		if (!pos && !(k % options->st)) {
		    printf("o");
		    fflush(stdout);
		}
	    }
	      
	    value_assign(z[pos], LB);
	    if (!check_poly(CS->next, data, nparam, pos+1, z, options)) {
		value_clear(LB);
		value_clear(UB);
		return 0;
	    }
	}
	value_set_si(z[pos], 0);
	value_clear(LB);
	value_clear(UB);
    }
    return 1;
} /* check_poly */
