#include <assert.h>
#include <iostream>
#include <barvinok/barvinok.h>
#include <barvinok/options.h>
#include <barvinok/util.h>
#include "argp.h"
#include "progname.h"
#include "evalue_convert.h"
#include "evalue_read.h"
#include "verify.h"

using std::cout;
using std::cerr;
using std::endl;

#define ALLOCN(type,n) (type*)malloc((n) * sizeof(type))

#define OPT_VARS  	    (BV_OPT_LAST+1)

struct argp_option argp_options[] = {
    { "variables",	    OPT_VARS,  	"list",	0,
	"comma separated list of variables over which to sum" },
    { "verbose",	    'v',  	0,	0, },
    { 0 }
};

struct options {
    struct convert_options   convert;
    struct verify_options    verify;
    char* var_list;
    int verbose;
};

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
    struct options *options = (struct options*) state->input;

    switch (key) {
    case ARGP_KEY_INIT:
	state->child_inputs[0] = &options->convert;
	state->child_inputs[1] = &options->verify;
	state->child_inputs[2] = options->verify.barvinok;
	options->var_list = NULL;
	options->verbose = 0;
	break;
    case 'v':
	options->verbose = 1;
	break;
    case OPT_VARS:
	options->var_list = strdup(arg);
	break;
    default:
	return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static int check_poly_sum(const struct check_poly_data *data,
			  int nparam, Value *z,
			  const struct verify_options *options);

struct check_poly_sum_data : public check_poly_data {
    Polyhedron	    	**S;
    evalue		 *EP;
    evalue		 *sum;

    check_poly_sum_data(Value *z, evalue *EP, evalue *sum) :
	    		EP(EP), sum(sum) {
	this->z = z;
	this->check = check_poly_sum;
    }
};

static void sum(Polyhedron *S, int pos, const check_poly_sum_data *data,
		evalue *s, const struct verify_options *options)
{
    if (!S) {
	evalue *e = evalue_eval(data->EP, data->z+1);
	eadd(e, s);
	free_evalue_refs(e);
	free(e);
    } else {
	Value LB, UB;
	int ok;
	value_init(LB);
	value_init(UB);
	ok = !(lower_upper_bounds(1+pos, S, data->z, &LB, &UB));
	assert(ok);
	for (; value_le(LB, UB); value_increment(LB, LB)) {
	    value_assign(data->z[1+pos], LB);
	    sum(S->next, pos+1, data, s, options);
	}
	value_set_si(data->z[1+pos], 0);
	value_clear(LB);
	value_clear(UB);
    }
}

static evalue *sum(const check_poly_sum_data *data,
			  const struct verify_options *options)
{
    evalue *s = evalue_zero();
    for (int i = 0; i < data->EP->x.p->size/2; ++i)
	if (!emptyQ2(data->S[i]))
	    sum(data->S[i], 0, data, s, options);
    return s;
}

static int check_poly_sum(const struct check_poly_data *data,
			  int nparam, Value *z,
			  const struct verify_options *options)
{
    const check_poly_sum_data *sum_data;
    sum_data = static_cast<const check_poly_sum_data *>(data);
    evalue *e, *s;
    int k;
    int ok;

    e = evalue_eval(sum_data->sum, z);
    if (options->print_all) {
	printf("sum(");
	value_print(stdout, VALUE_FMT, z[0]);
	for (k = 1; k < nparam; ++k) {
	    printf(", ");
	    value_print(stdout, VALUE_FMT, z[k]);
	}
	printf(") = ");
	value_print(stdout, VALUE_FMT, e->x.n);
	if (value_notone_p(e->d)) {
	    printf("/");
	    value_print(stdout, VALUE_FMT, e->d);
	}
    }

    s = sum(sum_data, options);

    if (options->print_all) {
	printf(", sum(EP) = ");
	value_print(stdout, VALUE_FMT, s->x.n);
	if (value_notone_p(s->d)) {
	    printf("/");
	    value_print(stdout, VALUE_FMT, s->d);
	}
	printf(". ");
    }

    ok = eequal(e, s);

    if (!ok) {
	printf("\n"); 
	fflush(stdout);
	fprintf(stderr,"Error !\n");
	fprintf(stderr,"sum( ");
	value_print(stderr, VALUE_FMT, z[0]);
	for (k = 1; k < nparam; ++k) {
	    fprintf(stderr, ", ");
	    value_print(stderr, VALUE_FMT, z[k]);
	}
	fprintf(stderr," ) should be ");
	value_print(stderr, VALUE_FMT, s->x.n);
	if (value_notone_p(s->d)) {
	    fprintf(stderr, "/");
	    value_print(stderr, VALUE_FMT, s->d);
	}
	fprintf(stderr,", while summation gives ");
	value_print(stderr, VALUE_FMT, e->x.n);
	if (value_notone_p(e->d)) {
	    fprintf(stderr, "/");
	    value_print(stderr, VALUE_FMT, e->d);
	}
	fprintf(stderr, ".\n");
    } else if (options->print_all)
	printf("OK.\n");

    free_evalue_refs(s);
    free(s);
    free_evalue_refs(e);
    free(e);

    return ok;
}

static int verify(Polyhedron *P, evalue *sum, evalue *EP,
		  unsigned nvar, unsigned nparam, Vector *p,
		  struct verify_options *options)
{
    Polyhedron *CS;
    unsigned MaxRays = options->barvinok->MaxRays;
    int error = 0;

    CS = check_poly_context_scan(NULL, &P, P->Dimension, options);

    check_poly_init(P, options);

    if (!(CS && emptyQ2(CS))) {
	check_poly_sum_data data(p->p, EP, sum);
	data.S = ALLOCN(Polyhedron *, EP->x.p->size/2);
	for (int i = 0; i < EP->x.p->size/2; ++i) {
	    Polyhedron *A = EVALUE_DOMAIN(EP->x.p->arr[2*i]);
	    data.S[i] = Polyhedron_Scan(A, P, MaxRays & POL_NO_DUAL ? 0 : MaxRays);
	}
	error = !check_poly(CS, &data, nparam, 0, p->p+1+nvar, options);
	for (int i = 0; i < EP->x.p->size/2; ++i)
	    Domain_Free(data.S[i]);
	free(data.S);
    }

    if (!options->print_all)
	printf("\n");

    if (CS) {
	Domain_Free(CS);
	Domain_Free(P);
    }

    return error;
}

/*
 * Project on final dim dimensions
 */
Polyhedron *DomainProject(Polyhedron *D, unsigned dim, unsigned MaxRays)
{
    Polyhedron *P;
    Polyhedron *R;

    R = Polyhedron_Project(D, dim);
    for (P = D->next; P; P = P->next) {
	Polyhedron *R2 = Polyhedron_Project(P, dim);
	Polyhedron *R3 = DomainUnion(R, R2, MaxRays);
	Polyhedron_Free(R2);
	Domain_Free(R);
	R = R3;
    }
    return R;
}

static int verify(evalue *EP, evalue *sum, unsigned nvar, unsigned nparam,
		  struct verify_options *options)
{
    Vector *p;

    p = Vector_Alloc(nvar+nparam+2);
    value_set_si(p->p[nvar+nparam+1], 1);

    assert(value_zero_p(EP->d));
    assert(EP->x.p->type == partition);

    Polyhedron *EP_D = EVALUE_DOMAIN(EP->x.p->arr[0]);
    Polyhedron *D = Polyhedron_Project(EP_D, nparam);

    for (int i = 1; i < EP->x.p->size/2; ++i) {
	Polyhedron *D2 = D;
	EP_D = DomainProject(EVALUE_DOMAIN(EP->x.p->arr[2*i]), nparam,
			     options->barvinok->MaxRays);
	D = DomainUnion(EP_D, D, options->barvinok->MaxRays);
	Domain_Free(D2);
    }

    int error = 0;

    for (Polyhedron *P = D; P; P = P->next) {
	error = verify(P, sum, EP, nvar, nparam, p, options);
	if (error && !options->continue_on_error)
	    break;
    }

    Domain_Free(D);
    Vector_Free(p);

    return error;
}

int main(int argc, char **argv)
{
    evalue *EP;
    char **all_vars = NULL;
    unsigned nvar;
    unsigned nparam;
    struct options options;
    struct barvinok_options *bv_options = barvinok_options_new_with_defaults();
    static struct argp_child argp_children[] = {
	{ &convert_argp,    	0,	"input conversion",	1 },
	{ &verify_argp,    	0,	"verification",		2 },
	{ &barvinok_argp,    	0,	"barvinok options",	3 },
	{ 0 }
    };
    static struct argp argp = { argp_options, parse_opt, 0, 0, argp_children };
    int result = 0;

    options.verify.barvinok = bv_options;
    set_program_name(argv[0]);
    argp_parse(&argp, argc, argv, 0, 0, &options);

    EP = evalue_read_from_file(stdin, options.var_list, &all_vars,
			       &nvar, &nparam, bv_options->MaxRays);
    assert(EP);

    if (options.verify.verify)
	verify_options_set_range(&options.verify, nvar+nparam);

    evalue_convert(EP, &options.convert, options.verbose, nparam, all_vars);

    if (EVALUE_IS_ZERO(*EP))
	print_evalue(stdout, EP, all_vars);
    else {
	evalue *sum = barvinok_summate(EP, nvar, bv_options);
	if (options.verify.verify)
	    result = verify(EP, sum, nvar, nparam, &options.verify);
	else
	    print_evalue(stdout, sum, all_vars+nvar);
	free_evalue_refs(sum);
	free(sum);
    }

    free_evalue_refs(EP);
    free(EP);

    if (options.var_list)
	free(options.var_list);
    Free_ParamNames(all_vars, nvar+nparam);
    barvinok_options_free(bv_options);
    return result;
}
