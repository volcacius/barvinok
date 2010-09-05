#include <assert.h>
#include <iostream>
#include <barvinok/evalue.h>
#include <barvinok/options.h>
#include <barvinok/util.h>
#include <barvinok/barvinok.h>
#include "bound_options.h"
#include "evalue_convert.h"
#include "evalue_read.h"
#include "verify.h"

using std::cout;
using std::cerr;
using std::endl;

#define ALLOCN(type,n) (type*)malloc((n) * sizeof(type))

struct verify_point_bound {
	struct verify_point_data vpd;
	isl_pw_qpolynomial *pwqp;
	isl_pw_qpolynomial_fold *pwf;
	enum isl_fold type;
};

static int verify_point(__isl_take isl_point *pnt, void *user)
{
	int i;
	unsigned nparam;
	struct verify_point_bound *vpb = (struct verify_point_bound *) user;
	isl_int v, n, d, b, t;
	isl_pw_qpolynomial *pwqp;
	isl_qpolynomial *bound;
	isl_qpolynomial *opt;
	const char *minmax;
	int sign;
	int ok;
	int cst;
	FILE *out = vpb->vpd.options->print_all ? stdout : stderr;

	vpb->vpd.n--;

	if (vpb->type == isl_fold_max) {
		minmax = "max";
		sign = 1;
	} else {
		minmax = "min";
		sign = -1;
	}

	isl_int_init(b);
	isl_int_init(t);
	isl_int_init(v);
	isl_int_init(n);
	isl_int_init(d);

	pwqp = isl_pw_qpolynomial_copy(vpb->pwqp);

	nparam = isl_pw_qpolynomial_dim(pwqp, isl_dim_param);
	for (i = 0; i < nparam; ++i) {
		isl_point_get_coordinate(pnt, isl_dim_param, i, &v);
		pwqp = isl_pw_qpolynomial_fix_dim(pwqp, isl_dim_param, i, v);
	}

	bound = isl_pw_qpolynomial_fold_eval(isl_pw_qpolynomial_fold_copy(vpb->pwf),
						isl_point_copy(pnt));

	if (sign > 0)
		opt = isl_pw_qpolynomial_max(pwqp);
	else
		opt = isl_pw_qpolynomial_min(pwqp);

	cst = isl_qpolynomial_is_cst(opt, &n, &d);
	if (cst != 1)
		goto error;
	if (sign > 0)
		isl_int_fdiv_q(v, n, d);
	else
		isl_int_cdiv_q(v, n, d);

	cst = isl_qpolynomial_is_cst(bound, &n, &d);
	if (cst != 1)
		goto error;
	if (sign > 0)
		isl_int_fdiv_q(b, n, d);
	else
		isl_int_cdiv_q(b, n, d);

	if (sign > 0)
		ok = value_ge(b, v);
	else
		ok = value_le(b, v);

	if (vpb->vpd.options->print_all || !ok) {
		fprintf(out, "%s(", minmax);
		for (i = 0; i < nparam; ++i) {
			if (i)
				fprintf(out, ", ");
			isl_point_get_coordinate(pnt, isl_dim_param, i, &t);
			isl_int_print(out, t, 0);
		}
		fprintf(out, ") = ");
		isl_int_print(out, n, 0);
		if (!isl_int_is_one(d)) {
			fprintf(out, "/");
			isl_int_print(out, d, 0);
			fprintf(out, " (");
			isl_int_print(out, b, 0);
			fprintf(out, ")");
		}
		fprintf(out, ", %s(EP) = ", minmax);
		isl_int_print(out, v, 0);
		if (ok)
			fprintf(out, ". OK\n");
		else
			fprintf(out, ". NOT OK\n");
	} else if ((vpb->vpd.n % vpb->vpd.s) == 0) {
		printf("o");
		fflush(stdout);
	}

	if (0) {
error:
		ok = 0;
	}

	isl_qpolynomial_free(bound);
	isl_qpolynomial_free(opt);
	isl_point_free(pnt);

	isl_int_clear(t);
	isl_int_clear(d);
	isl_int_clear(n);
	isl_int_clear(v);
	isl_int_clear(b);

	if (!ok)
		vpb->vpd.error = 1;

	if (vpb->vpd.options->continue_on_error)
		ok = 1;

	return (vpb->vpd.n >= 1 && ok) ? 0 : -1;
}

static int verify(__isl_keep isl_pw_qpolynomial_fold *pwf,
	__isl_take isl_pw_qpolynomial *pwqp, enum isl_fold type,
	struct verify_options *options)
{
	struct verify_point_bound vpb = { { options } };
	isl_ctx *ctx;
	isl_set *context;
	int r;

	ctx = isl_pw_qpolynomial_fold_get_ctx(pwf);
	vpb.pwf = pwf;
	vpb.type = type;
	vpb.pwqp = pwqp;
	context = isl_pw_qpolynomial_fold_domain(
					isl_pw_qpolynomial_fold_copy(vpb.pwf));
	context = verify_context_set_bounds(context, options);

	r = verify_point_data_init(&vpb.vpd, context);

	if (r == 0)
		isl_set_foreach_point(context, verify_point, &vpb);
	if (vpb.vpd.error)
		r = -1;

	isl_set_free(context);
	isl_pw_qpolynomial_free(pwqp);

	verify_point_data_fini(&vpb.vpd);

	return r;
}

static __isl_give isl_pw_qpolynomial_fold *iterate(
	__isl_take isl_pw_qpolynomial *pwqp, enum isl_fold type)
{
	isl_dim *dim = isl_pw_qpolynomial_get_dim(pwqp);
	isl_set *set;
	isl_qpolynomial *qp;
	isl_qpolynomial_fold *fold;
	unsigned nvar;

	assert(isl_dim_size(dim, isl_dim_param) == 0);
	nvar = isl_dim_size(dim, isl_dim_set);

	if (type == isl_fold_min)
		qp = isl_pw_qpolynomial_min(pwqp);
	else
		qp = isl_pw_qpolynomial_max(pwqp);

	qp = isl_qpolynomial_drop_dims(qp, isl_dim_set, 0, nvar);
	fold = isl_qpolynomial_fold_alloc(type, qp);
	dim = isl_dim_drop(dim, isl_dim_set, 0, nvar);
	set = isl_set_universe(dim);

	return isl_pw_qpolynomial_fold_alloc(type, set, fold);
}

struct bv_split_data {
	int size;
	__isl_give isl_pw_qpolynomial *pwqp_less;
	__isl_give isl_pw_qpolynomial *pwqp_more;
};

static int split_on_size(__isl_take isl_set *set,
	__isl_take isl_qpolynomial *qp, void *user)
{
	struct bv_split_data *data = (struct bv_split_data *)user;
	int bounded;
	isl_set *set_np;
	isl_pw_qpolynomial *pwqp;
	int nparam;

	nparam = isl_set_dim(set, isl_dim_param);
	set_np = isl_set_move_dims(isl_set_copy(set), isl_dim_set, 0,
					isl_dim_param, 0, nparam);
	bounded = isl_set_is_bounded(set_np);
	assert(bounded >= 0);
	if (bounded) {
		isl_pw_qpolynomial *pwqp;
		isl_qpolynomial *cst;
		isl_int m;
		int is_cst;

		pwqp = isl_set_card(set_np);
		cst = isl_pw_qpolynomial_max(pwqp);
		isl_int_init(m);
		is_cst = isl_qpolynomial_is_cst(cst, &m, NULL);
		isl_qpolynomial_free(cst);
		assert(is_cst);
		bounded = isl_int_cmp_si(m, data->size) <= 0;
		isl_int_clear(m);
	} else
		isl_set_free(set_np);

	pwqp = isl_pw_qpolynomial_alloc(set, qp);
	if (bounded)
		data->pwqp_less = isl_pw_qpolynomial_add_disjoint(
						data->pwqp_less, pwqp);
	else
		data->pwqp_more = isl_pw_qpolynomial_add_disjoint(
						data->pwqp_more, pwqp);

	return 0;
}

/*
 * Split (partition) pwpq into a partition with (sub)domains containing
 * size integer points or less and a partition with (sub)domains
 * containing more integer points.
 */
static int split_on_domain_size(__isl_take isl_pw_qpolynomial *pwqp,
	__isl_give isl_pw_qpolynomial **pwqp_less,
	__isl_give isl_pw_qpolynomial **pwqp_more,
	int size)
{
	isl_dim *dim;
	struct bv_split_data data = { size };
	int r;

	dim = isl_pw_qpolynomial_get_dim(pwqp);
	data.pwqp_less = isl_pw_qpolynomial_zero(isl_dim_copy(dim));
	data.pwqp_more = isl_pw_qpolynomial_zero(dim);
	r = isl_pw_qpolynomial_foreach_piece(pwqp, &split_on_size, &data);
	*pwqp_less = data.pwqp_less;
	*pwqp_more = data.pwqp_more;
	isl_pw_qpolynomial_free(pwqp);

	return r;
}

static __isl_give isl_pw_qpolynomial_fold *optimize(
	__isl_take isl_pw_qpolynomial *pwqp, enum isl_fold type,
	struct options *options)
{
    isl_pw_qpolynomial_fold *pwf;

    if (options->iterate > 0) {
	isl_pw_qpolynomial *pwqp_less, *pwqp_more;
	isl_pw_qpolynomial_fold *pwf_less, *pwf_more;
	split_on_domain_size(pwqp, &pwqp_less, &pwqp_more, options->iterate);
	pwf_less = iterate(pwqp_less, type);
	pwf_more = isl_pw_qpolynomial_bound(pwqp_more, type, NULL);
	pwf = isl_pw_qpolynomial_fold_fold(pwf_less, pwf_more);
    } else if (options->iterate)
	pwf = iterate(pwqp, type);
    else
	pwf = isl_pw_qpolynomial_bound(pwqp, type, NULL);
    return pwf;
}

static int optimize(__isl_take isl_pw_qpolynomial *pwqp, struct options *options)
{
    int print_solution = 1;
    int result = 0;
    isl_pw_qpolynomial_fold *pwf;
    enum isl_fold type = options->lower ? isl_fold_min : isl_fold_max;

    if (options->verify->verify) {
	isl_dim *dim = isl_pw_qpolynomial_get_dim(pwqp);
	unsigned total = isl_dim_total(dim);
	isl_dim_free(dim);
	verify_options_set_range(options->verify, total);
	if (!options->verify->barvinok->verbose)
	    print_solution = 0;
    }

    pwf = optimize(isl_pw_qpolynomial_copy(pwqp), type, options);
    assert(pwf);
    if (print_solution) {
	isl_ctx *ctx = isl_pw_qpolynomial_get_ctx(pwqp);
	isl_printer *p = isl_printer_to_file(ctx, stdout);
	p = isl_printer_print_pw_qpolynomial_fold(p, pwf);
	p = isl_printer_end_line(p);
	isl_printer_free(p);
    }
    if (options->verify->verify) {
	enum isl_fold type = options->lower ? isl_fold_min : isl_fold_max;
	result = verify(pwf, pwqp, type, options->verify);
    } else
	isl_pw_qpolynomial_free(pwqp);
    isl_pw_qpolynomial_fold_free(pwf);

    return result;
}

int main(int argc, char **argv)
{
    evalue *EP;
    const char **all_vars = NULL;
    unsigned nvar;
    unsigned nparam;
    struct options *options = options_new_with_defaults();
    isl_ctx *ctx;
    isl_dim *dim;
    isl_pw_qpolynomial *pwqp;
    int result = 0;

    argc = options_parse(options, argc, argv, ISL_ARG_ALL);
    ctx = isl_ctx_alloc_with_options(options_arg, options);

    EP = evalue_read_from_file(stdin, options->var_list, &all_vars,
			       &nvar, &nparam, options->verify->barvinok->MaxRays);
    assert(EP);

    evalue_convert(EP, options->convert, options->verify->barvinok->verbose,
			nvar+nparam, all_vars);

    dim = isl_dim_set_alloc(ctx, nparam, 0);
    for (int i = 0; i < nparam; ++i)
	dim = isl_dim_set_name(dim, isl_dim_param, i, all_vars[nvar + i]);
    dim = isl_dim_insert(dim, isl_dim_param, 0, nvar);
    pwqp = isl_pw_qpolynomial_from_evalue(dim, EP);
    pwqp = isl_pw_qpolynomial_move_dims(pwqp, isl_dim_set, 0, isl_dim_param, 0, nvar);
    evalue_free(EP);

    if (options->split)
	pwqp = isl_pw_qpolynomial_split_periods(pwqp, options->split);

    result = optimize(pwqp, options);

    Free_ParamNames(all_vars, nvar+nparam);
    isl_ctx_free(ctx);
    return result;
}
