#include <assert.h>
#include <limits.h>
#include <math.h>
#include <sys/times.h>
#include <bernstein/bernstein.h>
#include <barvinok/options.h>
#include <barvinok/bernstein.h>
#include <barvinok/util.h>
#include "argp.h"
#include "progname.h"
#include "evalue_read.h"
#include "verify.h"
#include "range.h"

using std::cout;
using std::cerr;
using std::endl;
using namespace GiNaC;
using namespace bernstein;
using namespace barvinok;

#define	METHOD_BERNSTEIN	0
#define METHOD_PROPAGATION	1

static struct {
    int	    method;
} methods[] = {
    { METHOD_BERNSTEIN },
    { METHOD_PROPAGATION },
};

#define nr_methods (sizeof(methods)/sizeof(*methods))

struct argp_option argp_options[] = {
    { 0 }
};

struct options {
    struct verify_options    verify;
};

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
    struct options *options = (struct options*) state->input;

    switch (key) {
    case ARGP_KEY_INIT:
	state->child_inputs[0] = options->verify.barvinok;
	state->child_inputs[1] = &options->verify;
	break;
    default:
	return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

struct result_data {
    Value		    n;
    double		    RE_sum[nr_methods];

    clock_t		    ticks[nr_methods];
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
    value_init(result->n);
}

void result_data_clear(struct result_data *result)
{
    int i;
    value_clear(result->n);
}

void result_data_print(struct result_data *result, int n)
{
    int i;

    fprintf(stderr, "%d", result->ticks[0]/n);
    for (i = 1; i < nr_methods; ++i)
	fprintf(stderr, ", %d", result->ticks[i]/n);
    fprintf(stderr, "\n");

    fprintf(stderr, "%d/%d", result->size[0], n);
    for (i = 1; i < nr_methods; ++i)
	fprintf(stderr, ", %d/%d", result->size[i], n);
    fprintf(stderr, "\n");

    fprintf(stderr, "%g\n", VALUE_TO_DOUBLE(result->n));
    fprintf(stderr, "%g", result->RE_sum[0]/VALUE_TO_DOUBLE(result->n));
    for (i = 1; i < nr_methods; ++i)
	fprintf(stderr, ", %g", result->RE_sum[i]/VALUE_TO_DOUBLE(result->n));
    fprintf(stderr, "\n");
}

static int test_bound(const struct check_poly_data *data,
			  int nparam, Value *z,
			  const struct verify_options *options);

struct test_bound_data : public check_EP_data {
    piecewise_lst	**pl;
    struct result_data   *result;

    test_bound_data(evalue *EP, piecewise_lst **pl, result_data *result) :
			pl(pl), result(result) {
	this->EP = EP;
	this->cp.check = test_bound;
    }
};

static int test_bound(const struct check_poly_data *data,
			  int nparam, Value *z,
			  const struct verify_options *options)
{
    const test_bound_data *tb_data = (const test_bound_data *)data;
    Value max, min, exact, approx;
    Value n, d;

    value_init(exact);
    value_init(approx);
    value_init(max);
    value_init(min);
    value_init(n);
    value_init(d);

    evalue_optimum(tb_data, &max, 1);
    evalue_optimum(tb_data, &min, -1);
    value_assign(exact, max);
    value_subtract(exact, exact, min);
    value_add_int(exact, exact, 1);

    if (options->print_all) {
	value_print(stderr, "max: "VALUE_FMT, max);
	value_print(stderr, ", min: "VALUE_FMT, min);
	value_print(stderr, ", range: "VALUE_FMT, exact);
    }

    value_increment(tb_data->result->n, tb_data->result->n);
    for (int i = 0; i < nr_methods; ++i) {
	double error;

	tb_data->pl[2*i]->evaluate(nparam, z, &n, &d);
	mpz_fdiv_q(max, n, d);
	tb_data->pl[2*i+1]->evaluate(nparam, z, &n, &d);
	mpz_cdiv_q(min, n, d);

	value_assign(approx, max);
	value_subtract(approx, approx, min);
	value_add_int(approx, approx, 1);
	if (options->print_all)
	    value_print(stderr, ", "VALUE_FMT, approx);

	assert(value_ge(approx, exact));
	value_subtract(approx, approx, exact);

	error = ::abs(VALUE_TO_DOUBLE(approx)) / VALUE_TO_DOUBLE(exact);
	if (options->print_all)
	    fprintf(stderr, " (%g)", error);
	tb_data->result->RE_sum[i] += error;
    }

    if (options->print_all)
	fprintf(stderr, "\n");

    value_clear(n);
    value_clear(d);
    value_clear(max);
    value_clear(min);
    value_clear(exact);
    value_clear(approx);

    return 1;
}

static void test(evalue *EP, unsigned nvar, unsigned nparam,
		 piecewise_lst **pl, struct result_data *result,
		 struct verify_options *options)
{
    test_bound_data data(EP, pl, result);
    check_EP(&data, nvar, nparam, options);
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
    char **all_vars = NULL;
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
	    if (methods[i].method == METHOD_BERNSTEIN) {
		pl[2*i+j] = evalue_bernstein_coefficients(NULL, poly, U, params,
						   options->barvinok);
		if (sign == BV_BERNSTEIN_MIN)
		    pl[2*i+j]->minimize();
		else
		    pl[2*i+j]->maximize();
	    } else {
		pl[2*i+j] = evalue_range_propagation(NULL, poly, params,
					      options->barvinok);
		if (sign == BV_BERNSTEIN_MIN)
		    pl[2*i+j]->sign = -1;
		else
		    pl[2*i+j]->sign = 1;
	    }
	}
	times(&en_cpu);
	result->ticks[i] = en_cpu.tms_utime - st_cpu.tms_utime;
	if (options->barvinok->verbose)
	    for (int j = 0; j < 2; ++j)
		cerr << *pl[2*i+j] << endl;
	result->size[i] = number_of_polynomials(pl[2*i]);
	result->size[i] += number_of_polynomials(pl[2*i+1]);
    }
    test(EP, nvar, nparam, pl, result, options);

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

	result_data_print(&result, 1);

	value_addto(all_result.n, all_result.n, result.n);
	for (i = 0; i < nr_methods; ++i) {
	    all_result.RE_sum[i] += result.RE_sum[i];
	    all_result.ticks[i] += result.ticks[i];
	    all_result.size[i] += result.size[i];
	}

	result_data_clear(&result);

	fprintf(stderr, "average\n");
	result_data_print(&all_result, n);
    }

    result_data_clear(&all_result);

    barvinok_options_free(bv_options);

    return 0;
}
