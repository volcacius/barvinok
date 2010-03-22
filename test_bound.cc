#include <assert.h>
#include <limits.h>
#include <math.h>
#include <bernstein/bernstein.h>
#include <barvinok/options.h>
#include <barvinok/bernstein.h>
#include <barvinok/util.h>
#include "argp.h"
#include "progname.h"
#include "evalue_read.h"
#include "verify.h"
#include "range.h"

#ifdef HAVE_SYS_TIMES_H

#include <sys/times.h>

typedef clock_t		my_clock_t;

static my_clock_t time_diff(struct tms *before, struct tms *after)
{
	return after->tms_utime - before->tms_utime;
}

#else

typedef int		my_clock_t;

struct tms {};
static void times(struct tms* time)
{
}
static my_clock_t time_diff(struct tms *before, struct tms *after)
{
	return 0;
}

#endif

using std::cout;
using std::cerr;
using std::endl;
using namespace GiNaC;
using namespace bernstein;
using namespace barvinok;

static struct {
    int	    method;
} methods[] = {
    { BV_BOUND_BERNSTEIN },
    { BV_BOUND_RANGE },
};

#define nr_methods (sizeof(methods)/sizeof(*methods))

struct argp_option argp_options[] = {
    { "quiet",	    'q' },
    { 0 }
};

struct options {
    struct verify_options    verify;
    int quiet;
};

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
    struct options *options = (struct options*) state->input;

    switch (key) {
    case ARGP_KEY_INIT:
	state->child_inputs[0] = options->verify.barvinok;
	state->child_inputs[1] = &options->verify;
	options->quiet = 0;
	break;
    case 'q':
	options->quiet = 1;
	break;
    default:
	return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

struct result_data {
    isl_int		    n;
    double		    RE_sum[nr_methods];

    my_clock_t		    ticks[nr_methods];
    size_t		    size[nr_methods];
};

void result_data_init(struct result_data *result)
{
    int i;
    for (i = 0; i < nr_methods; ++i) {
	result->RE_sum[i] = 0;
	result->ticks[i] = 0;
	result->size[i] = 0;
    }
    isl_int_init(result->n);
}

void result_data_clear(struct result_data *result)
{
    int i;
    isl_int_clear(result->n);
}

void result_data_print(struct result_data *result, int n)
{
    int i;

    fprintf(stderr, "%d", result->ticks[0]/n);
    for (i = 1; i < nr_methods; ++i)
	fprintf(stderr, ", %d", result->ticks[i]/n);
    fprintf(stderr, "\n");

    fprintf(stderr, "%zd/%d", result->size[0], n);
    for (i = 1; i < nr_methods; ++i)
	fprintf(stderr, ", %zd/%d", result->size[i], n);
    fprintf(stderr, "\n");

    fprintf(stderr, "%g\n", isl_int_get_d(result->n));
    fprintf(stderr, "%g", result->RE_sum[0]/isl_int_get_d(result->n));
    for (i = 1; i < nr_methods; ++i)
	fprintf(stderr, ", %g", result->RE_sum[i]/isl_int_get_d(result->n));
    fprintf(stderr, "\n");
}

struct verify_point_bound {
	struct verify_point_data vpd;
	struct result_data *result;
	isl_pw_qpolynomial *pwqp;
	isl_pw_qpolynomial_fold **pwf;
};

static int verify_point(__isl_take isl_point *pnt, void *user)
{
	struct verify_point_bound *vpb = (struct verify_point_bound *) user;
	const struct verify_options *options = vpb->vpd.options;
	int i;
	unsigned nparam;
	isl_int max, min, exact, approx;
	isl_int n, d;
	isl_qpolynomial *opt;
	isl_pw_qpolynomial *pwqp;
	int cst;

	vpb->vpd.n--;

	isl_int_init(max);
	isl_int_init(min);
	isl_int_init(exact);
	isl_int_init(approx);
	isl_int_init(n);
	isl_int_init(d);

	pwqp = isl_pw_qpolynomial_copy(vpb->pwqp);

	nparam = isl_pw_qpolynomial_dim(pwqp, isl_dim_param);
	for (i = 0; i < nparam; ++i) {
		isl_point_get_coordinate(pnt, isl_dim_param, i, &n);
		pwqp = isl_pw_qpolynomial_fix_dim(pwqp, isl_dim_param, i, n);
	}

	opt = isl_pw_qpolynomial_max(isl_pw_qpolynomial_copy(pwqp));
	cst = isl_qpolynomial_is_cst(opt, &n, &d);
	isl_qpolynomial_free(opt);
	assert(cst == 1);
	isl_int_fdiv_q(max, n, d);

	opt = isl_pw_qpolynomial_min(pwqp);
	cst = isl_qpolynomial_is_cst(opt, &n, &d);
	isl_qpolynomial_free(opt);
	assert(cst == 1);
	isl_int_cdiv_q(min, n, d);

	isl_int_sub(exact, max, min);
	isl_int_add_ui(exact, exact, 1);

	if (options->print_all) {
		fprintf(stderr, "max: "); isl_int_print(stderr, max, 0);
		fprintf(stderr, ", min: "); isl_int_print(stderr, min, 0);
		fprintf(stderr, ", range: "); isl_int_print(stderr, exact, 0);
	}

	for (int i = 0; i < nr_methods; ++i) {
		double error;

		opt = isl_pw_qpolynomial_fold_eval(
				isl_pw_qpolynomial_fold_copy(vpb->pwf[2 * i]),
				isl_point_copy(pnt));
		cst = isl_qpolynomial_is_cst(opt, &n, &d);
		isl_qpolynomial_free(opt);
		assert(cst == 1);
		isl_int_fdiv_q(max, n, d);
	
		opt = isl_pw_qpolynomial_fold_eval(
				isl_pw_qpolynomial_fold_copy(vpb->pwf[2 * i + 1]),
				isl_point_copy(pnt));
		cst = isl_qpolynomial_is_cst(opt, &n, &d);
		isl_qpolynomial_free(opt);
		assert(cst == 1);
		isl_int_cdiv_q(min, n, d);

		isl_int_sub(approx, max, min);
		isl_int_add_ui(approx, approx, 1);
		if (options->print_all) {
			fprintf(stderr, ", "); isl_int_print(stderr, approx, 0);
		}

		assert(isl_int_ge(approx, exact));
		isl_int_sub(approx, approx, exact);

		error = fabs(isl_int_get_d(approx)) / isl_int_get_d(exact);
		if (options->print_all)
			fprintf(stderr, " (%g)", error);
		vpb->result->RE_sum[i] += error;
	}

	if (options->print_all) {
		fprintf(stderr, "\n");
	} else if ((vpb->vpd.n % vpb->vpd.s) == 0) {
		printf("o");
		fflush(stdout);
	}

	isl_int_clear(max);
	isl_int_clear(min);
	isl_int_clear(exact);
	isl_int_clear(approx);
	isl_int_clear(n);
	isl_int_clear(d);

	isl_point_free(pnt);

	return 0;
}

static void test(evalue *EP, unsigned nvar,
		 const GiNaC::exvector &params,
		 piecewise_lst **pl, struct result_data *result,
		 struct verify_options *options)
{
	struct verify_point_bound vpb = { { options }, result };
	isl_ctx *ctx = isl_ctx_alloc();
	isl_dim *dim;
	isl_set *context;
	int r;
	int i;

	vpb.pwf = isl_alloc_array(ctx, isl_pw_qpolynomial_fold *, 2 * nr_methods);

	dim = isl_dim_set_alloc(ctx, params.size(), 0);
	for (i = 0; i < 2 * nr_methods; ++i)
		vpb.pwf[i] = isl_pw_qpolynomial_fold_from_ginac(isl_dim_copy(dim),
								pl[i], params);
	isl_dim_free(dim);
	dim = isl_dim_set_alloc(ctx, nvar + params.size(), 0);
	vpb.pwqp = isl_pw_qpolynomial_from_evalue(dim, EP);
	vpb.pwqp = isl_pw_qpolynomial_move(vpb.pwqp, isl_dim_set, 0,
						isl_dim_param, 0, nvar);
	context = isl_pw_qpolynomial_domain(isl_pw_qpolynomial_copy(vpb.pwqp));
	context = isl_set_remove(context, isl_dim_set, 0, nvar);
	context = verify_context_set_bounds(context, options);

	r = verify_point_data_init(&vpb.vpd, context);
	assert(r == 0);
	isl_int_set_si(result->n, vpb.vpd.n);

	isl_set_foreach_point(context, verify_point, &vpb);
	assert(!vpb.vpd.error);

	isl_set_free(context);
	isl_pw_qpolynomial_free(vpb.pwqp);
	for (i = 0; i < 2 * nr_methods; ++i)
		isl_pw_qpolynomial_fold_free(vpb.pwf[i]);

	isl_ctx_free(ctx);

	verify_point_data_fini(&vpb.vpd);
	free(vpb.pwf);
}

static int number_of_polynomials(piecewise_lst *pl)
{
    int n = 0;
    for (int i = 0; i < pl->list.size(); ++i)
	n += pl->list[i].second.nops();
    return n;
}

void handle(FILE *in, struct result_data *result, struct verify_options *options)
{
    evalue *EP, *upper, *lower;
    const char **all_vars = NULL;
    unsigned nvar;
    unsigned nparam;
    Polyhedron *U;
    exvector params;
    piecewise_lst *pl[2*nr_methods];

    EP = evalue_read_from_file(in, NULL, &all_vars,
			       &nvar, &nparam, options->barvinok->MaxRays);
    assert(EP);
    if (EVALUE_IS_ZERO(*EP)) {
	evalue_free(EP);
	Free_ParamNames(all_vars, nvar+nparam);
	return;
    }

    upper = evalue_dup(EP);
    lower = evalue_dup(EP);
    evalue_frac2polynomial(upper, 1, options->barvinok->MaxRays);
    evalue_frac2polynomial(lower, -1, options->barvinok->MaxRays);

    U = Universe_Polyhedron(nparam);
    params = constructParameterVector(all_vars+nvar, nparam);

    for (int i = 0; i < nr_methods; ++i) {
	struct tms st_cpu;
	struct tms en_cpu;

	times(&st_cpu);
	for (int j = 0; j < 2; ++j) {
	    evalue *poly = j == 0 ? upper : lower;
	    int sign = j == 0 ? BV_BERNSTEIN_MAX : BV_BERNSTEIN_MIN;
	    options->barvinok->bernstein_optimize = sign;
	    if (methods[i].method == BV_BOUND_BERNSTEIN) {
		pl[2*i+j] = evalue_bernstein_coefficients(NULL, poly, U, params,
						   options->barvinok);
		if (sign == BV_BERNSTEIN_MIN)
		    pl[2*i+j]->minimize();
		else
		    pl[2*i+j]->maximize();
	    } else {
		pl[2*i+j] = evalue_range_propagation(NULL, poly, params,
					      options->barvinok);
	    }
	}
	times(&en_cpu);
	result->ticks[i] = time_diff(&en_cpu, &st_cpu);
	if (options->barvinok->verbose)
	    for (int j = 0; j < 2; ++j)
		cerr << *pl[2*i+j] << endl;
	result->size[i] = number_of_polynomials(pl[2*i]);
	result->size[i] += number_of_polynomials(pl[2*i+1]);
    }
    test(EP, nvar, params, pl, result, options);

    for (int i = 0; i < 2*nr_methods; ++i)
	delete pl[i];
    Polyhedron_Free(U);
    evalue_free(EP);
    evalue_free(lower);
    evalue_free(upper);
    Free_ParamNames(all_vars, nvar+nparam);
}

int main(int argc, char **argv)
{
    struct barvinok_options *bv_options = barvinok_options_new_with_defaults();
    char path[PATH_MAX+1];
    struct result_data all_result;
    int n = 0;
    static struct argp_child argp_children[] = {
	{ &barvinok_argp,    	0,	0,  		0 },
	{ &verify_argp,    	0,	"verification",		BV_GRP_LAST+1 },
	{ 0 }
    };
    static struct argp argp = { argp_options, parse_opt, 0, 0, argp_children };
    struct options options;

    options.verify.barvinok = bv_options;
    set_program_name(argv[0]);
    argp_parse(&argp, argc, argv, 0, 0, &options);

    if (options.verify.M == INT_MIN)
	options.verify.M = 100;
    if (options.verify.m == INT_MAX)
	options.verify.m = -100;

    result_data_init(&all_result);

    while (fgets(path, sizeof(path), stdin)) {
	struct result_data result;
	FILE *in;
	int i;

	++n;
	result_data_init(&result);
	fprintf(stderr, "%s", path);
	*strchr(path, '\n') = '\0';
	in = fopen(path, "r");
	assert(in);
	handle(in, &result, &options.verify);
	fclose(in);

	if (!options.quiet)
	    result_data_print(&result, 1);

	isl_int_add(all_result.n, all_result.n, result.n);
	for (i = 0; i < nr_methods; ++i) {
	    all_result.RE_sum[i] += result.RE_sum[i];
	    all_result.ticks[i] += result.ticks[i];
	    all_result.size[i] += result.size[i];
	}

	result_data_clear(&result);

	if (!options.quiet) {
	    fprintf(stderr, "average\n");
	    result_data_print(&all_result, n);
	}
    }

    result_data_clear(&all_result);

    barvinok_options_free(bv_options);

    return 0;
}
