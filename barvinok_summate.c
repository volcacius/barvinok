#include <assert.h>
#include <barvinok/barvinok.h>
#include <barvinok/options.h>
#include <barvinok/util.h>
#include "evalue_convert.h"
#include "evalue_read.h"
#include "verify.h"

struct options {
	struct convert_options   *convert;
	struct verify_options    *verify;
	char *var_list;
};

struct isl_arg options_arg[] = {
ISL_ARG_CHILD(struct options, verify, NULL,
	verify_options_arg, "verification")
ISL_ARG_CHILD(struct options, convert, NULL,
	convert_options_arg, "output conversion")
ISL_ARG_STR(struct options, var_list, 0, "variables", "list", NULL,
	"comma separated list of variables over which to sum")
ISL_ARG_END
};

ISL_ARG_DEF(options, struct options, options_arg)

struct verify_point_sum {
	struct verify_point_data vpd;
	isl_pw_qpolynomial *pwqp;
	isl_pw_qpolynomial *sum;

	isl_pw_qpolynomial *fixed;
	isl_qpolynomial *manual;
};

static int manual_sum(__isl_take isl_point *pnt, void *user)
{
	struct verify_point_sum *vps = (struct verify_point_sum *) user;
	isl_qpolynomial *qp;

	qp = isl_pw_qpolynomial_eval(isl_pw_qpolynomial_copy(vps->fixed), pnt);
	vps->manual = isl_qpolynomial_add(vps->manual, qp);

	return 0;
}

static int verify_point(__isl_take isl_point *pnt, void *user)
{
	struct verify_point_sum *vps = (struct verify_point_sum *) user;
	int i;
	int ok;
	unsigned nparam;
	isl_int v;
	isl_set *dom;
	isl_qpolynomial *eval;
	int r;
	FILE *out = vps->vpd.options->print_all ? stdout : stderr;

	vps->vpd.n--;

	isl_int_init(v);
	vps->fixed = isl_pw_qpolynomial_copy(vps->pwqp);
	nparam = isl_pw_qpolynomial_dim(vps->sum, isl_dim_param);
	for (i = 0; i < nparam; ++i) {
		isl_point_get_coordinate(pnt, isl_dim_param, i, &v);
		vps->fixed = isl_pw_qpolynomial_fix_dim(vps->fixed,
						    	isl_dim_param, i, v);
	}

	eval = isl_pw_qpolynomial_eval(isl_pw_qpolynomial_copy(vps->sum),
					isl_point_copy(pnt));

	vps->manual = isl_qpolynomial_zero(isl_pw_qpolynomial_get_dim(vps->pwqp));
	dom = isl_pw_qpolynomial_domain(isl_pw_qpolynomial_copy(vps->fixed));
	r = isl_set_foreach_point(dom, &manual_sum, user);
	isl_set_free(dom);
	if (r < 0)
		goto error;

	ok = isl_qpolynomial_is_equal(eval, vps->manual);

	if (vps->vpd.options->print_all || !ok) {
		isl_ctx *ctx = isl_pw_qpolynomial_get_ctx(vps->pwqp);
		isl_printer *p = isl_printer_to_file(ctx, out);
		fprintf(out, "sum(");
		for (i = 0; i < nparam; ++i) {
			if (i)
				fprintf(out, ", ");
			isl_point_get_coordinate(pnt, isl_dim_param, i, &v);
			isl_int_print(out, v, 0);
		}
		fprintf(out, ") = ");
		p = isl_printer_print_qpolynomial(p, eval);
		fprintf(out, ", sum(EP) = ");
		p = isl_printer_print_qpolynomial(p, vps->manual);
		if (ok)
			fprintf(out, ". OK\n");
		else
			fprintf(out, ". NOT OK\n");
		isl_printer_free(p);
	} else if ((vps->vpd.n % vps->vpd.s) == 0) {
		printf("o");
		fflush(stdout);
	}

	if (0) {
error:
		ok = 0;
	}
	isl_qpolynomial_free(vps->manual);
	isl_pw_qpolynomial_free(vps->fixed);
	isl_qpolynomial_free(eval);
	isl_int_clear(v);
	isl_point_free(pnt);

	if (!ok)
		vps->vpd.error = 1;

	if (vps->vpd.options->continue_on_error)
		ok = 1;

	return (vps->vpd.n >= 1 && ok) ? 0 : -1;
}

static int verify(__isl_keep isl_pw_qpolynomial *pwqp,
	__isl_take isl_pw_qpolynomial *sum, struct verify_options *options)
{
	struct verify_point_sum vps = { { options } };
	isl_set *context;
	int r;

	vps.pwqp = pwqp;
	vps.sum = sum;

	context = isl_pw_qpolynomial_domain(isl_pw_qpolynomial_copy(sum));
	context = verify_context_set_bounds(context, options);

	r = verify_point_data_init(&vps.vpd, context);

	if (r == 0)
		isl_set_foreach_point(context, verify_point, &vps);
	if (vps.vpd.error)
		r = -1;

	isl_set_free(context);

	verify_point_data_fini(&vps.vpd);

	return r;
}

int main(int argc, char **argv)
{
    int i;
    evalue *EP;
    const char **all_vars = NULL;
    unsigned nvar;
    unsigned nparam;
    int result = 0;
    isl_ctx *ctx;
    isl_dim *dim;
    isl_pw_qpolynomial *pwqp;
    isl_pw_qpolynomial *sum;
    struct options *options = options_new_with_defaults();

    argc = options_parse(options, argc, argv, ISL_ARG_ALL);
    ctx = isl_ctx_alloc_with_options(options_arg, options);

    EP = evalue_read_from_file(stdin, options->var_list, &all_vars,
			       &nvar, &nparam, options->verify->barvinok->MaxRays);
    assert(EP);

    if (options->verify->verify)
	verify_options_set_range(options->verify, nvar+nparam);

    evalue_convert(EP, options->convert, options->verify->barvinok->verbose,
			nparam, all_vars);

    dim = isl_dim_set_alloc(ctx, nparam, 0);
    for (i = 0; i < nparam; ++i)
	dim = isl_dim_set_name(dim, isl_dim_param, i, all_vars[nvar + i]);
    dim = isl_dim_insert(dim, isl_dim_param, 0, nvar);
    pwqp = isl_pw_qpolynomial_from_evalue(dim, EP);
    pwqp = isl_pw_qpolynomial_move_dims(pwqp, isl_dim_set, 0, isl_dim_param, 0, nvar);
    evalue_free(EP);

    sum = isl_pw_qpolynomial_sum(isl_pw_qpolynomial_copy(pwqp));
    if (options->verify->verify)
	result = verify(pwqp, sum, options->verify);
    else {
	isl_printer *p = isl_printer_to_file(ctx, stdout);
	p = isl_printer_print_pw_qpolynomial(p, sum);
	p = isl_printer_end_line(p);
	isl_printer_free(p);
    }
    isl_pw_qpolynomial_free(sum);
    isl_pw_qpolynomial_free(pwqp);

    Free_ParamNames(all_vars, nvar+nparam);
    isl_ctx_free(ctx);
    return result;
}
