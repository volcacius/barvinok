#include <assert.h>
#include <stdlib.h>
#include <barvinok/options.h>
#include <barvinok/util.h>
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
	options->r = -1;
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
	options->r = atoi(arg);
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
    /* If the user didn't set m or M, then we try to adjust the window
     * to the context in check_poly_context_scan.
     */
    if (options->m == INT_MAX && options->M == INT_MIN)
	options->r = r;
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

static Polyhedron *project_on(Polyhedron *P, int i)
{
    unsigned dim = P->Dimension;
    Matrix *T;
    Polyhedron *I;

    if (dim == 1)
	return Polyhedron_Copy(P);

    T = Matrix_Alloc(2, dim+1);
    value_set_si(T->p[0][i], 1);
    value_set_si(T->p[1][dim], 1);
    I = Polyhedron_Image(P, T, P->NbConstraints);
    Matrix_Free(T);
    return I;
}

static void set_bounds(Polyhedron *C, Value **rows, int i, unsigned nparam,
			const struct verify_options *options)
{
    Value min;
    Value max;

    value_init(min);
    value_init(max);
    value_set_si(min, options->m);
    value_set_si(max, options->M);

    if (options->r > 0) {
	Polyhedron *I = project_on(C, i);
	line_minmax(I, &min, &max);
	if (value_cmp_si(min, options->M) >= 0)
	    value_add_int(max, min, options->r);
	else if (value_cmp_si(max, options->m) <= 0)
	    value_sub_int(min, max, options->r);
	else {
	    value_set_si(min, options->m);
	    value_set_si(max, options->M);
	}
    }

    value_set_si(rows[0][0], 1);
    value_set_si(rows[0][1+i], 1);
    value_oppose(rows[0][1+nparam], min);
    value_set_si(rows[1][0], 1);
    value_set_si(rows[1][1+i], -1);
    value_assign(rows[1][1+nparam], max);

    value_clear(min);
    value_clear(max);
}

Polyhedron *check_poly_context_scan(Polyhedron *P, Polyhedron **C,
				    unsigned nparam,
				    const struct verify_options *options)
{
    int i;
    Matrix *MM;
    Polyhedron *CC, *CC2, *CS, *U;
    unsigned MaxRays = options->barvinok->MaxRays;

    if (nparam <= 0)
	return NULL;

    if (!P)
	CC = *C;
    else {
	CC = Polyhedron_Project(P, nparam);
	if (*C) {
	    CC2 = DomainIntersection(*C, CC, MaxRays);
	    Domain_Free(CC);
	    CC = CC2;
	}
    }

    /* Intersect context with range */
    MM = Matrix_Alloc(2*nparam, nparam+2);
    for (i = 0; i < nparam; ++i)
	set_bounds(CC, &MM->p[2*i], i, nparam, options);
    CC2 = AddConstraints(MM->p[0], 2*nparam, CC, options->barvinok->MaxRays);
    if (CC != *C)
	Domain_Free(CC);
    CC = CC2;
    U = Universe_Polyhedron(0);
    CS = Polyhedron_Scan(CC, U, MaxRays & POL_NO_DUAL ? 0 : MaxRays);
    Polyhedron_Free(U);
    *C = CC;
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

static void print_rational(FILE *out, Value n, Value d)
{
    value_print(out, VALUE_FMT, n);
    if (value_notone_p(d)) {
	fprintf(out, "/");
	value_print(out, VALUE_FMT, d);
    }
}

void check_poly_print(int ok, int nparam, Value *z,
		      Value want_n, Value want_d,
		      Value got_n, Value got_d,
		      const char *op, const char *check,
		      const char *long_op,
		      const struct verify_options *options)
{
    int k;

    if (options->print_all) {
	printf("%s(", op);
	value_print(stdout, VALUE_FMT, z[0]);
	for (k = 1; k < nparam; ++k) {
	    printf(", ");
	    value_print(stdout, VALUE_FMT, z[k]);
	}
	printf(") = ");
	print_rational(stdout, got_n, got_d);
	printf(", %s = ", check);
	print_rational(stdout, want_n, want_d);
	printf(". ");
    }

    if (!ok) {
	printf("\n");
	fflush(stdout);
	fprintf(stderr, "Error !\n");
	fprintf(stderr, "%s(", op);
	value_print(stderr, VALUE_FMT, z[0]);
	for (k = 1; k < nparam; ++k) {
	    fprintf(stderr,", ");
	    value_print(stderr, VALUE_FMT, z[k]);
	}
	fprintf(stderr, ") should be ");
	print_rational(stderr, want_n, want_d);
	fprintf(stderr, ", while %s gives ", long_op);
	print_rational(stderr, got_n, got_d);
	fprintf(stderr, ".\n");
    } else if (options->print_all)
	printf("OK.\n");
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
