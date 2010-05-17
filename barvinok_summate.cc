#include <assert.h>
#include <iostream>
#include <barvinok/barvinok.h>
#include <barvinok/options.h>
#include <barvinok/util.h>
#include "barvinok_summate_options.h"
#include "evalue_convert.h"
#include "evalue_read.h"
#include "verify.h"

using std::cout;
using std::cerr;
using std::endl;

#define ALLOCN(type,n) (type*)malloc((n) * sizeof(type))

static int check_poly_sum(const struct check_poly_data *data,
			  int nparam, Value *z,
			  const struct verify_options *options);

struct check_poly_sum_data : public check_EP_data {
    evalue		 *sum;

    check_poly_sum_data(const evalue *EP, evalue *sum) : sum(sum) {
	this->EP = EP;
	this->cp.check = check_poly_sum;
    }
};

static void sum(Polyhedron *S, int pos, const check_poly_sum_data *data,
		evalue *s, const struct verify_options *options)
{
    if (!S) {
	evalue *e = evalue_eval(data->EP, data->cp.z+1);
	eadd(e, s);
	evalue_free(e);
    } else {
	Value LB, UB;
	int ok;
	value_init(LB);
	value_init(UB);
	ok = !(lower_upper_bounds(1+pos, S, data->cp.z, &LB, &UB));
	assert(ok);
	for (; value_le(LB, UB); value_increment(LB, LB)) {
	    value_assign(data->cp.z[1+pos], LB);
	    sum(S->next, pos+1, data, s, options);
	}
	value_set_si(data->cp.z[1+pos], 0);
	value_clear(LB);
	value_clear(UB);
    }
}

static evalue *sum(const check_poly_sum_data *data,
			  const struct verify_options *options)
{
    evalue *s = evalue_zero();
    for (int i = 0; i < data->n_S; ++i)
	if (!emptyQ2(data->S[i]))
	    sum(data->S[i], 0, data, s, options);
    return s;
}

static int check_poly_sum(const struct check_poly_data *data,
			  int nparam, Value *z,
			  const struct verify_options *options)
{
    const check_poly_sum_data *sum_data;
    sum_data = (const check_poly_sum_data *)data;
    evalue *e, *s;
    int k;
    int ok;

    e = evalue_eval(sum_data->sum, z);
    s = sum(sum_data, options);

    ok = eequal(e, s);

    check_poly_print(ok, nparam, z, s->x.n, s->d, e->x.n, e->d,
		     "sum", "sum(EP)", "summation", options);

    evalue_free(s);
    evalue_free(e);

    return ok;
}

static int verify(evalue *EP, evalue *sum, unsigned nvar, unsigned nparam,
		  struct verify_options *options)
{
    check_poly_sum_data data(EP, sum);
    return !check_EP(&data, nvar, nparam, options);
}

int main(int argc, char **argv)
{
    evalue *EP;
    const char **all_vars = NULL;
    unsigned nvar;
    unsigned nparam;
    int result = 0;
    struct options *options = options_new_with_defaults();

    argc = options_parse(options, argc, argv, ISL_ARG_ALL);

    EP = evalue_read_from_file(stdin, options->var_list, &all_vars,
			       &nvar, &nparam, options->verify->barvinok->MaxRays);
    assert(EP);

    if (options->verify->verify)
	verify_options_set_range(options->verify, nvar+nparam);

    evalue_convert(EP, options->convert, options->verify->barvinok->verbose,
			nparam, all_vars);

    if (EVALUE_IS_ZERO(*EP))
	print_evalue(stdout, EP, all_vars);
    else {
	evalue *sum = barvinok_summate(EP, nvar, options->verify->barvinok);
	assert(sum);
	if (options->verify->verify)
	    result = verify(EP, sum, nvar, nparam, options->verify);
	else
	    print_evalue(stdout, sum, all_vars+nvar);
	evalue_free(sum);
    }

    evalue_free(EP);

    Free_ParamNames(all_vars, nvar+nparam);
    options_free(options);
    return result;
}
