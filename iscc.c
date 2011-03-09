#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <isl/obj.h>
#include <isl/stream.h>
#include <isl/vertices.h>
#include <isl/flow.h>
#include <isl_obj_list.h>
#include <isl_obj_str.h>
#include <barvinok/isl.h>
#include <barvinok/options.h>

#include "config.h"

#ifdef HAVE_CLOOG
#include <cloog/isl/cloog.h>
#endif

static int isl_bool_false = 0;
static int isl_bool_true = 1;
static int isl_bool_error = -1;

enum iscc_op { ISCC_READ, ISCC_WRITE, ISCC_SOURCE, ISCC_VERTICES,
	       ISCC_LAST, ISCC_ANY, ISCC_BEFORE, ISCC_UNDER,
	       ISCC_TYPEOF,
	       ISCC_N_OP };
static const char *op_name[ISCC_N_OP] = {
	[ISCC_READ] = "read",
	[ISCC_WRITE] = "write",
	[ISCC_SOURCE] = "source",
	[ISCC_VERTICES] = "vertices",
	[ISCC_LAST] = "last",
	[ISCC_ANY] = "any",
	[ISCC_BEFORE] = "before",
	[ISCC_UNDER] = "under",
	[ISCC_TYPEOF] = "typeof"
};
static enum isl_token_type iscc_op[ISCC_N_OP];

struct isl_arg_choice iscc_format[] = {
	{"isl",		ISL_FORMAT_ISL},
	{"omega",	ISL_FORMAT_OMEGA},
	{"polylib",	ISL_FORMAT_POLYLIB},
	{"ext-polylib",	ISL_FORMAT_EXT_POLYLIB},
	{"latex",	ISL_FORMAT_LATEX},
	{"C",		ISL_FORMAT_C},
	{0}
};

struct iscc_options {
	struct barvinok_options	*barvinok;
	unsigned		 format;
	int			 io;
};

struct isl_arg iscc_options_arg[] = {
ISL_ARG_CHILD(struct iscc_options, barvinok, "barvinok", barvinok_options_arg,
	"barvinok options")
ISL_ARG_CHOICE(struct iscc_options, format, 0, "format", \
	iscc_format,	ISL_FORMAT_ISL, "output format")
ISL_ARG_BOOL(struct iscc_options, io, 0, "io", 1,
	"allow read and write operations")
ISL_ARG_END
};

ISL_ARG_DEF(iscc_options, struct iscc_options, iscc_options_arg)
ISL_ARG_CTX_DEF(iscc_options, struct iscc_options, iscc_options_arg)

static void *isl_obj_bool_copy(void *v)
{
	return v;
}

static void isl_obj_bool_free(void *v)
{
}

static __isl_give isl_printer *isl_obj_bool_print(__isl_take isl_printer *p,
	void *v)
{
	if (v == &isl_bool_true)
		return isl_printer_print_str(p, "True");
	else if (v == &isl_bool_false)
		return isl_printer_print_str(p, "False");
	else
		return isl_printer_print_str(p, "Error");
}

static void *isl_obj_bool_add(void *v1, void *v2)
{
	return v1;
}

struct isl_obj_vtable isl_obj_bool_vtable = {
	isl_obj_bool_copy,
	isl_obj_bool_add,
	isl_obj_bool_print,
	isl_obj_bool_free
};
#define isl_obj_bool		(&isl_obj_bool_vtable)

int *isl_bool_from_int(int res)
{
	return res < 0 ? &isl_bool_error : res ? &isl_bool_true : &isl_bool_false;
}

int *union_map_is_equal(__isl_take isl_union_map *map1,
	__isl_take isl_union_map *map2)
{
	int res = isl_union_map_is_equal(map1, map2);
	isl_union_map_free(map1);
	isl_union_map_free(map2);
	return isl_bool_from_int(res);
}
int *union_set_is_equal(__isl_take isl_union_set *set1,
	__isl_take isl_union_set *set2)
{
	return union_map_is_equal((isl_union_map *)set1, (isl_union_map *)set2);
}

int *union_map_is_subset(__isl_take isl_union_map *map1,
	__isl_take isl_union_map *map2)
{
	int res = isl_union_map_is_subset(map1, map2);
	isl_union_map_free(map1);
	isl_union_map_free(map2);
	return isl_bool_from_int(res);
}
int *union_set_is_subset(__isl_take isl_union_set *set1,
	__isl_take isl_union_set *set2)
{
	return union_map_is_subset((isl_union_map *)set1, (isl_union_map *)set2);
}

int *union_map_is_strict_subset(__isl_take isl_union_map *map1,
	__isl_take isl_union_map *map2)
{
	int res = isl_union_map_is_strict_subset(map1, map2);
	isl_union_map_free(map1);
	isl_union_map_free(map2);
	return isl_bool_from_int(res);
}
int *union_set_is_strict_subset(__isl_take isl_union_set *set1,
	__isl_take isl_union_set *set2)
{
	return union_map_is_strict_subset((isl_union_map *)set1,
					  (isl_union_map *)set2);
}

int *union_map_is_superset(__isl_take isl_union_map *map1,
	__isl_take isl_union_map *map2)
{
	return union_map_is_subset(map2, map1);
}
int *union_set_is_superset(__isl_take isl_union_set *set1,
	__isl_take isl_union_set *set2)
{
	return union_set_is_subset(set2, set1);
}

int *union_map_is_strict_superset(__isl_take isl_union_map *map1,
	__isl_take isl_union_map *map2)
{
	return union_map_is_strict_subset(map2, map1);
}
int *union_set_is_strict_superset(__isl_take isl_union_set *set1,
	__isl_take isl_union_set *set2)
{
	return union_set_is_strict_subset(set2, set1);
}

extern struct isl_obj_vtable isl_obj_list_vtable;
#define isl_obj_list		(&isl_obj_list_vtable)

typedef void *(*isc_bin_op_fn)(void *lhs, void *rhs);
struct isc_bin_op {
	enum isl_token_type	op;
	isl_obj_type		lhs;
	isl_obj_type		rhs;
	isl_obj_type		res;
	isc_bin_op_fn		fn;
};
struct isc_named_bin_op {
	char			*name;
	struct isc_bin_op	op;
};

struct iscc_at {
	isl_union_pw_qpolynomial *upwqp;
	isl_union_pw_qpolynomial *res;
};

static int eval_at(__isl_take isl_point *pnt, void *user)
{
	struct iscc_at *at = (struct iscc_at *) user;
	isl_qpolynomial *qp;
	isl_set *set;

	set = isl_set_from_point(isl_point_copy(pnt));
	qp = isl_union_pw_qpolynomial_eval(
				isl_union_pw_qpolynomial_copy(at->upwqp), pnt);

	at->res = isl_union_pw_qpolynomial_add(at->res,
			isl_union_pw_qpolynomial_from_pw_qpolynomial(
				isl_pw_qpolynomial_alloc(set, qp)));

	return 0;
}

__isl_give isl_union_pw_qpolynomial *isl_union_pw_qpolynomial_at(
	__isl_take isl_union_pw_qpolynomial *upwqp,
	__isl_take isl_union_set *uset)
{
	struct iscc_at at;

	at.upwqp = upwqp;
	at.res = isl_union_pw_qpolynomial_zero(isl_union_set_get_dim(uset));

	isl_union_set_foreach_point(uset, eval_at, &at);

	isl_union_pw_qpolynomial_free(upwqp);
	isl_union_set_free(uset);

	return at.res;
}

struct iscc_fold_at {
	isl_union_pw_qpolynomial_fold *upwf;
	isl_union_pw_qpolynomial *res;
};

static int eval_fold_at(__isl_take isl_point *pnt, void *user)
{
	struct iscc_fold_at *at = (struct iscc_fold_at *) user;
	isl_qpolynomial *qp;
	isl_set *set;

	set = isl_set_from_point(isl_point_copy(pnt));
	qp = isl_union_pw_qpolynomial_fold_eval(
			    isl_union_pw_qpolynomial_fold_copy(at->upwf), pnt);

	at->res = isl_union_pw_qpolynomial_add(at->res,
			isl_union_pw_qpolynomial_from_pw_qpolynomial(
				isl_pw_qpolynomial_alloc(set, qp)));

	return 0;
}

__isl_give isl_union_pw_qpolynomial *isl_union_pw_qpolynomial_fold_at(
	__isl_take isl_union_pw_qpolynomial_fold *upwf,
	__isl_take isl_union_set *uset)
{
	struct iscc_fold_at at;

	at.upwf = upwf;
	at.res = isl_union_pw_qpolynomial_zero(isl_union_set_get_dim(uset));

	isl_union_set_foreach_point(uset, eval_fold_at, &at);

	isl_union_pw_qpolynomial_fold_free(upwf);
	isl_union_set_free(uset);

	return at.res;
}

static __isl_give isl_union_pw_qpolynomial_fold *union_pw_qpolynomial_add_union_pw_qpolynomial_fold(
	__isl_take isl_union_pw_qpolynomial *upwqp,
	__isl_take isl_union_pw_qpolynomial_fold *upwf)
{
	return isl_union_pw_qpolynomial_fold_add_union_pw_qpolynomial(upwf,
									upwqp);
}

static __isl_give struct isl_list *union_map_apply_union_pw_qpolynomial_fold(
	__isl_take isl_union_map *umap,
	__isl_take isl_union_pw_qpolynomial_fold *upwf)
{
	isl_ctx *ctx;
	struct isl_list *list;
	int tight;

	ctx = isl_union_map_get_ctx(umap);
	list = isl_list_alloc(ctx, 2);
	if (!list)
		goto error2;

	list->obj[0].type = isl_obj_union_pw_qpolynomial_fold;
	list->obj[0].v = isl_union_map_apply_union_pw_qpolynomial_fold(umap,
							upwf, &tight);
	list->obj[1].type = isl_obj_bool;
	list->obj[1].v = tight ? &isl_bool_true : &isl_bool_false;
	if (tight < 0 || !list->obj[0].v)
		goto error;

	return list;
error2:
	isl_union_map_free(umap);
	isl_union_pw_qpolynomial_fold_free(upwf);
error:
	isl_list_free(list);
	return NULL;
}

static __isl_give isl_union_pw_qpolynomial *union_pw_qpolynomial_int_mul(
	__isl_take isl_union_pw_qpolynomial *upwqp, __isl_take isl_int_obj *i)
{
	isl_int v;

	if (!i)
		goto error;

	isl_int_init(v);
	isl_int_obj_get_int(i, &v);
	upwqp = isl_union_pw_qpolynomial_mul_isl_int(upwqp, v);
	isl_int_clear(v);

	isl_int_obj_free(i);

	return upwqp;
error:
	isl_union_pw_qpolynomial_free(upwqp);
	return NULL;
}

static __isl_give isl_union_pw_qpolynomial *int_union_pw_qpolynomial_mul(
	__isl_take isl_int_obj *i, __isl_take isl_union_pw_qpolynomial *upwqp)
{
	return union_pw_qpolynomial_int_mul(upwqp, i);
}

static __isl_give isl_union_pw_qpolynomial_fold *union_pw_qpolynomial_fold_int_mul(
	__isl_take isl_union_pw_qpolynomial_fold *upwf,
	__isl_take isl_int_obj *i)
{
	isl_int v;

	if (!i)
		goto error;

	isl_int_init(v);
	isl_int_obj_get_int(i, &v);
	upwf = isl_union_pw_qpolynomial_fold_mul_isl_int(upwf, v);
	isl_int_clear(v);

	isl_int_obj_free(i);

	return upwf;
error:
	isl_union_pw_qpolynomial_fold_free(upwf);
	return NULL;
}

static __isl_give isl_union_pw_qpolynomial_fold *int_union_pw_qpolynomial_fold_mul(
	__isl_take isl_int_obj *i,
	__isl_take isl_union_pw_qpolynomial_fold *upwf)
{
	return union_pw_qpolynomial_fold_int_mul(upwf, i);
}

struct isc_bin_op bin_ops[] = {
	{ '+',	isl_obj_int,	isl_obj_int, isl_obj_int,
		(isc_bin_op_fn) &isl_int_obj_add },
	{ '-',	isl_obj_int,	isl_obj_int, isl_obj_int,
		(isc_bin_op_fn) &isl_int_obj_sub },
	{ '*',	isl_obj_int,	isl_obj_int, isl_obj_int,
		(isc_bin_op_fn) &isl_int_obj_mul },
	{ '+',	isl_obj_union_set,	isl_obj_union_set,
		isl_obj_union_set,
		(isc_bin_op_fn) &isl_union_set_union },
	{ '+',	isl_obj_union_map,	isl_obj_union_map,
		isl_obj_union_map,
		(isc_bin_op_fn) &isl_union_map_union },
	{ '-',	isl_obj_union_set,	isl_obj_union_set,
		isl_obj_union_set,
		(isc_bin_op_fn) &isl_union_set_subtract },
	{ '-',	isl_obj_union_map,	isl_obj_union_map,
		isl_obj_union_map,
		(isc_bin_op_fn) &isl_union_map_subtract },
	{ '*',	isl_obj_union_set,	isl_obj_union_set,
		isl_obj_union_set,
		(isc_bin_op_fn) &isl_union_set_intersect },
	{ '*',	isl_obj_union_map,	isl_obj_union_map,
		isl_obj_union_map,
		(isc_bin_op_fn) &isl_union_map_intersect },
	{ '*',	isl_obj_union_map,	isl_obj_union_set,
		isl_obj_union_map,
		(isc_bin_op_fn) &isl_union_map_intersect_domain },
	{ '.',	isl_obj_union_map,	isl_obj_union_map,
		isl_obj_union_map,
		(isc_bin_op_fn) &isl_union_map_apply_range },
	{ '.',	isl_obj_union_map,	isl_obj_union_pw_qpolynomial,
		isl_obj_union_pw_qpolynomial,
		(isc_bin_op_fn) &isl_union_map_apply_union_pw_qpolynomial },
	{ '.',	isl_obj_union_map,	isl_obj_union_pw_qpolynomial_fold,
		isl_obj_list,
		(isc_bin_op_fn) &union_map_apply_union_pw_qpolynomial_fold },
	{ ISL_TOKEN_TO,	isl_obj_union_set,	isl_obj_union_set,
		isl_obj_union_map,
		(isc_bin_op_fn) &isl_union_map_from_domain_and_range },
	{ '=', isl_obj_union_set,	isl_obj_union_set,	isl_obj_bool,
		(isc_bin_op_fn) &union_set_is_equal },
	{ '=', isl_obj_union_map,	isl_obj_union_map,	isl_obj_bool,
		(isc_bin_op_fn) &union_map_is_equal },
	{ ISL_TOKEN_LE, isl_obj_union_set,	isl_obj_union_set,
		isl_obj_bool, (isc_bin_op_fn) &union_set_is_subset },
	{ ISL_TOKEN_LE, isl_obj_union_map,	isl_obj_union_map,
		isl_obj_bool, (isc_bin_op_fn) &union_map_is_subset },
	{ ISL_TOKEN_LT, isl_obj_union_set,	isl_obj_union_set,
		isl_obj_bool, (isc_bin_op_fn) &union_set_is_strict_subset },
	{ ISL_TOKEN_LT, isl_obj_union_map,	isl_obj_union_map,
		isl_obj_bool, (isc_bin_op_fn) &union_map_is_strict_subset },
	{ ISL_TOKEN_GE, isl_obj_union_set,	isl_obj_union_set,
		isl_obj_bool, (isc_bin_op_fn) &union_set_is_superset },
	{ ISL_TOKEN_GE, isl_obj_union_map,	isl_obj_union_map,
		isl_obj_bool, (isc_bin_op_fn) &union_map_is_superset },
	{ ISL_TOKEN_GT, isl_obj_union_set,	isl_obj_union_set,
		isl_obj_bool, (isc_bin_op_fn) &union_set_is_strict_superset },
	{ ISL_TOKEN_GT, isl_obj_union_map,	isl_obj_union_map,
		isl_obj_bool, (isc_bin_op_fn) &union_map_is_strict_superset },
	{ ISL_TOKEN_LEX_LE,	isl_obj_union_set,	isl_obj_union_set,
		isl_obj_union_map,
		(isc_bin_op_fn) &isl_union_set_lex_le_union_set },
	{ ISL_TOKEN_LEX_LT,	isl_obj_union_set,	isl_obj_union_set,
		isl_obj_union_map,
		(isc_bin_op_fn) &isl_union_set_lex_lt_union_set },
	{ ISL_TOKEN_LEX_GE,	isl_obj_union_set,	isl_obj_union_set,
		isl_obj_union_map,
		(isc_bin_op_fn) &isl_union_set_lex_ge_union_set },
	{ ISL_TOKEN_LEX_GT,	isl_obj_union_set,	isl_obj_union_set,
		isl_obj_union_map,
		(isc_bin_op_fn) &isl_union_set_lex_gt_union_set },
	{ ISL_TOKEN_LEX_LE,	isl_obj_union_map,	isl_obj_union_map,
		isl_obj_union_map,
		(isc_bin_op_fn) &isl_union_map_lex_le_union_map },
	{ ISL_TOKEN_LEX_LT,	isl_obj_union_map,	isl_obj_union_map,
		isl_obj_union_map,
		(isc_bin_op_fn) &isl_union_map_lex_lt_union_map },
	{ ISL_TOKEN_LEX_GE,	isl_obj_union_map,	isl_obj_union_map,
		isl_obj_union_map,
		(isc_bin_op_fn) &isl_union_map_lex_ge_union_map },
	{ ISL_TOKEN_LEX_GT,	isl_obj_union_map,	isl_obj_union_map,
		isl_obj_union_map,
		(isc_bin_op_fn) &isl_union_map_lex_gt_union_map },
	{ '.',	isl_obj_union_pw_qpolynomial_fold,
		isl_obj_union_pw_qpolynomial_fold,
		isl_obj_union_pw_qpolynomial_fold,
		(isc_bin_op_fn) &isl_union_pw_qpolynomial_fold_fold },
	{ '+',	isl_obj_union_pw_qpolynomial,	isl_obj_union_pw_qpolynomial,
		isl_obj_union_pw_qpolynomial,
		(isc_bin_op_fn) &isl_union_pw_qpolynomial_add },
	{ '+',	isl_obj_union_pw_qpolynomial,
		isl_obj_union_pw_qpolynomial_fold,
		isl_obj_union_pw_qpolynomial_fold,
		(isc_bin_op_fn) &union_pw_qpolynomial_add_union_pw_qpolynomial_fold },
	{ '+',	isl_obj_union_pw_qpolynomial_fold,
		isl_obj_union_pw_qpolynomial,
		isl_obj_union_pw_qpolynomial_fold,
		(isc_bin_op_fn) &isl_union_pw_qpolynomial_fold_add_union_pw_qpolynomial },
	{ '-',	isl_obj_union_pw_qpolynomial,	isl_obj_union_pw_qpolynomial,
		isl_obj_union_pw_qpolynomial,
		(isc_bin_op_fn) &isl_union_pw_qpolynomial_sub },
	{ '*',	isl_obj_int,	isl_obj_union_pw_qpolynomial,
		isl_obj_union_pw_qpolynomial,
		(isc_bin_op_fn) &int_union_pw_qpolynomial_mul },
	{ '*',	isl_obj_union_pw_qpolynomial,	isl_obj_int,
		isl_obj_union_pw_qpolynomial,
		(isc_bin_op_fn) &union_pw_qpolynomial_int_mul },
	{ '*',	isl_obj_int,	isl_obj_union_pw_qpolynomial_fold,
		isl_obj_union_pw_qpolynomial_fold,
		(isc_bin_op_fn) &int_union_pw_qpolynomial_fold_mul },
	{ '*',	isl_obj_union_pw_qpolynomial_fold,	isl_obj_int,
		isl_obj_union_pw_qpolynomial_fold,
		(isc_bin_op_fn) &union_pw_qpolynomial_fold_int_mul },
	{ '*',	isl_obj_union_pw_qpolynomial,	isl_obj_union_pw_qpolynomial,
		isl_obj_union_pw_qpolynomial,
		(isc_bin_op_fn) &isl_union_pw_qpolynomial_mul },
	{ '*',	isl_obj_union_pw_qpolynomial,	isl_obj_union_set,
		isl_obj_union_pw_qpolynomial,
		(isc_bin_op_fn) &isl_union_pw_qpolynomial_intersect_domain },
	{ '*',	isl_obj_union_pw_qpolynomial_fold,	isl_obj_union_set,
		isl_obj_union_pw_qpolynomial_fold,
		(isc_bin_op_fn) &isl_union_pw_qpolynomial_fold_intersect_domain },
	{ '@',	isl_obj_union_pw_qpolynomial, isl_obj_union_set,
		isl_obj_union_pw_qpolynomial,
		(isc_bin_op_fn) &isl_union_pw_qpolynomial_at },
	{ '@',	isl_obj_union_pw_qpolynomial_fold, isl_obj_union_set,
		isl_obj_union_pw_qpolynomial,
		(isc_bin_op_fn) &isl_union_pw_qpolynomial_fold_at },
	{ '%',	isl_obj_union_set,	isl_obj_union_set,
		isl_obj_union_set,
		(isc_bin_op_fn) &isl_union_set_gist },
	{ '%',	isl_obj_union_map,	isl_obj_union_map,
		isl_obj_union_map,
		(isc_bin_op_fn) &isl_union_map_gist },
	{ '%',	isl_obj_union_pw_qpolynomial,	isl_obj_union_set,
		isl_obj_union_pw_qpolynomial,
		(isc_bin_op_fn) &isl_union_pw_qpolynomial_gist },
	{ '%',	isl_obj_union_pw_qpolynomial_fold,	isl_obj_union_set,
		isl_obj_union_pw_qpolynomial_fold,
		(isc_bin_op_fn) &isl_union_pw_qpolynomial_fold_gist },
	{ '+',	isl_obj_str,	isl_obj_str,	isl_obj_str,
		(isc_bin_op_fn) &isl_str_concat },
	0
};

static __isl_give isl_union_map *map_after_map(__isl_take isl_union_map *umap1,
	__isl_take isl_union_map *umap2)
{
	return isl_union_map_apply_range(umap2, umap1);
}

static __isl_give isl_union_pw_qpolynomial *qpolynomial_after_map(
	__isl_take isl_union_pw_qpolynomial *upwqp,
	__isl_take isl_union_map *umap)
{
	return isl_union_map_apply_union_pw_qpolynomial(umap, upwqp);
}

static __isl_give struct isl_list *qpolynomial_fold_after_map(
	__isl_take isl_union_pw_qpolynomial_fold *upwf,
	__isl_take isl_union_map *umap)
{
	return union_map_apply_union_pw_qpolynomial_fold(umap, upwf);
}

struct isc_named_bin_op named_bin_ops[] = {
	{ "after",	{ -1, isl_obj_union_map,	isl_obj_union_map,
		isl_obj_union_map,
		(isc_bin_op_fn) &map_after_map } },
	{ "after",	{ -1, isl_obj_union_pw_qpolynomial,
		isl_obj_union_map, isl_obj_union_pw_qpolynomial,
		(isc_bin_op_fn) &qpolynomial_after_map } },
	{ "after",	{ -1, isl_obj_union_pw_qpolynomial_fold,
		isl_obj_union_map, isl_obj_list,
		(isc_bin_op_fn) &qpolynomial_fold_after_map } },
	{ "before",	{ -1, isl_obj_union_map,	isl_obj_union_map,
		isl_obj_union_map,
		(isc_bin_op_fn) &isl_union_map_apply_range } },
	{ "before",	{ -1, isl_obj_union_map,
		isl_obj_union_pw_qpolynomial, isl_obj_union_pw_qpolynomial,
		(isc_bin_op_fn) &isl_union_map_apply_union_pw_qpolynomial } },
	{ "before",	{ -1, isl_obj_union_map,
		isl_obj_union_pw_qpolynomial_fold, isl_obj_list,
		(isc_bin_op_fn) &union_map_apply_union_pw_qpolynomial_fold } },
	{ "cross",	{ -1,	isl_obj_union_set,	isl_obj_union_set,
		isl_obj_union_set,
		(isc_bin_op_fn) &isl_union_set_product } },
	{ "cross",	{ -1,	isl_obj_union_map,	isl_obj_union_map,
		isl_obj_union_map,
		(isc_bin_op_fn) &isl_union_map_product } },
	NULL
};

__isl_give isl_set *union_set_sample(__isl_take isl_union_set *uset)
{
	return isl_set_from_basic_set(isl_union_set_sample(uset));
}

__isl_give isl_map *union_map_sample(__isl_take isl_union_map *umap)
{
	return isl_map_from_basic_map(isl_union_map_sample(umap));
}

static __isl_give struct isl_list *union_pw_qpolynomial_upper_bound(
	__isl_take isl_union_pw_qpolynomial *upwqp)
{
	isl_ctx *ctx;
	struct isl_list *list;
	int tight;

	ctx = isl_union_pw_qpolynomial_get_ctx(upwqp);
	list = isl_list_alloc(ctx, 2);
	if (!list)
		goto error2;

	list->obj[0].type = isl_obj_union_pw_qpolynomial_fold;
	list->obj[0].v = isl_union_pw_qpolynomial_bound(upwqp,
							isl_fold_max, &tight);
	list->obj[1].type = isl_obj_bool;
	list->obj[1].v = tight ? &isl_bool_true : &isl_bool_false;
	if (tight < 0 || !list->obj[0].v)
		goto error;

	return list;
error2:
	isl_union_pw_qpolynomial_free(upwqp);
error:
	isl_list_free(list);
	return NULL;
}

#ifdef HAVE_CLOOG
void *map_codegen(void *arg)
{
	isl_dim *dim;
	isl_union_map *umap = (isl_union_map *)arg;
	isl_ctx *ctx = isl_union_map_get_ctx(umap);
	CloogState *state;
	CloogOptions *options;
	CloogDomain *context;
	CloogUnionDomain *ud;
	CloogInput *input;
	struct clast_stmt *stmt;

	state = cloog_isl_state_malloc(ctx);
	options = cloog_options_malloc(state);
	options->language = LANGUAGE_C;
	options->strides = 1;

	ud = cloog_union_domain_from_isl_union_map(isl_union_map_copy(umap));

	dim = isl_union_map_get_dim(umap);
	context = cloog_domain_from_isl_set(isl_set_universe(dim));

	input = cloog_input_alloc(context, ud);

	stmt = cloog_clast_create_from_input(input, options);
	clast_pprint(stdout, stmt, 0, options);
	cloog_clast_free(stmt);

error:
	cloog_options_free(options);
	cloog_state_free(state);
	isl_union_map_free(umap);
	return NULL;
}

void *set_codegen(void *arg)
{
	isl_dim *dim;
	isl_union_set *uset = (isl_union_set *)arg;
	isl_ctx *ctx = isl_union_set_get_ctx(uset);
	CloogState *state;
	CloogOptions *options;
	CloogDomain *context;
	CloogUnionDomain *ud;
	CloogInput *input;
	struct clast_stmt *stmt;

	if (isl_union_set_n_set(uset) > 1)
		isl_die(ctx, isl_error_invalid,
			"code generation for more than one domain "
			"requires a schedule", goto error);

	state = cloog_isl_state_malloc(ctx);
	options = cloog_options_malloc(state);
	options->language = LANGUAGE_C;
	options->strides = 1;

	ud = cloog_union_domain_from_isl_union_set(isl_union_set_copy(uset));

	dim = isl_union_set_get_dim(uset);
	context = cloog_domain_from_isl_set(isl_set_universe(dim));

	input = cloog_input_alloc(context, ud);

	stmt = cloog_clast_create_from_input(input, options);
	clast_pprint(stdout, stmt, 0, options);
	cloog_clast_free(stmt);

	cloog_options_free(options);
	cloog_state_free(state);
error:
	isl_union_set_free(uset);
	return NULL;
}
#endif

static int add_point(__isl_take isl_point *pnt, void *user)
{
	isl_union_set **scan = (isl_union_set **) user;

	*scan = isl_union_set_add_set(*scan, isl_set_from_point(pnt));

	return 0;
}

static __isl_give isl_union_set *union_set_scan(__isl_take isl_union_set *uset)
{
	isl_union_set *scan;

	scan = isl_union_set_empty(isl_union_set_get_dim(uset));

	if (isl_union_set_foreach_point(uset, add_point, &scan) < 0) {
		isl_union_set_free(scan);
		return uset;
	}

	isl_union_set_free(uset);
	return scan;
}

static __isl_give isl_union_map *union_map_scan(__isl_take isl_union_map *umap)
{
	return isl_union_set_unwrap(union_set_scan(isl_union_map_wrap(umap)));
}

static __isl_give isl_union_pw_qpolynomial *union_pw_qpolynomial_poly(
	__isl_take isl_union_pw_qpolynomial *upwqp)
{
	return isl_union_pw_qpolynomial_to_polynomial(upwqp, 0);
}

static __isl_give isl_union_pw_qpolynomial *union_pw_qpolynomial_lpoly(
	__isl_take isl_union_pw_qpolynomial *upwqp)
{
	return isl_union_pw_qpolynomial_to_polynomial(upwqp, -1);
}

static __isl_give isl_union_pw_qpolynomial *union_pw_qpolynomial_upoly(
	__isl_take isl_union_pw_qpolynomial *upwqp)
{
	return isl_union_pw_qpolynomial_to_polynomial(upwqp, 1);
}

typedef void *(*isc_un_op_fn)(void *arg);
struct isc_un_op {
	enum isl_token_type	op;
	isl_obj_type		arg;
	isl_obj_type		res;
	isc_un_op_fn		fn;
};
struct isc_named_un_op {
	char			*name;
	struct isc_un_op	op;
};
struct isc_named_un_op named_un_ops[] = {
	{"aff",	{ -1,	isl_obj_union_map,	isl_obj_union_map,
		(isc_un_op_fn) &isl_union_map_affine_hull } },
	{"aff",	{ -1,	isl_obj_union_set,	isl_obj_union_set,
		(isc_un_op_fn) &isl_union_set_affine_hull } },
	{"card",	{ -1,	isl_obj_union_set,
		isl_obj_union_pw_qpolynomial,
		(isc_un_op_fn) &isl_union_set_card } },
	{"card",	{ -1,	isl_obj_union_map,
		isl_obj_union_pw_qpolynomial,
		(isc_un_op_fn) &isl_union_map_card } },
	{"coalesce",	{ -1,	isl_obj_union_set,	isl_obj_union_set,
		(isc_un_op_fn) &isl_union_set_coalesce } },
	{"coalesce",	{ -1,	isl_obj_union_map,	isl_obj_union_map,
		(isc_un_op_fn) &isl_union_map_coalesce } },
	{"coalesce",	{ -1,	isl_obj_union_pw_qpolynomial,
		isl_obj_union_pw_qpolynomial,
		(isc_un_op_fn) &isl_union_pw_qpolynomial_coalesce } },
	{"coalesce",	{ -1,	isl_obj_union_pw_qpolynomial_fold,
		isl_obj_union_pw_qpolynomial_fold,
		(isc_un_op_fn) &isl_union_pw_qpolynomial_fold_coalesce } },
#ifdef HAVE_CLOOG
	{"codegen",	{ -1,	isl_obj_union_set, isl_obj_none,
		&set_codegen } },
	{"codegen",	{ -1,	isl_obj_union_map, isl_obj_none,
		&map_codegen } },
#endif
	{"deltas",	{ -1,	isl_obj_union_map,	isl_obj_union_set,
		(isc_un_op_fn) &isl_union_map_deltas } },
	{"deltas_map",	{ -1,	isl_obj_union_map,	isl_obj_union_map,
		(isc_un_op_fn) &isl_union_map_deltas_map } },
	{"dom",	{ -1,	isl_obj_union_map,	isl_obj_union_set,
		(isc_un_op_fn) &isl_union_map_domain } },
	{"dom",	{ -1,	isl_obj_union_pw_qpolynomial,	isl_obj_union_set,
		(isc_un_op_fn) &isl_union_pw_qpolynomial_domain } },
	{"dom",	{ -1,	isl_obj_union_pw_qpolynomial_fold,
		isl_obj_union_set,
		(isc_un_op_fn) &isl_union_pw_qpolynomial_fold_domain } },
	{"domain",	{ -1,	isl_obj_union_map,	isl_obj_union_set,
		(isc_un_op_fn) &isl_union_map_domain } },
	{"domain",	{ -1,	isl_obj_union_pw_qpolynomial,
		isl_obj_union_set,
		(isc_un_op_fn) &isl_union_pw_qpolynomial_domain } },
	{"domain",	{ -1,	isl_obj_union_pw_qpolynomial_fold,
		isl_obj_union_set,
		(isc_un_op_fn) &isl_union_pw_qpolynomial_fold_domain } },
	{"domain_map",	{ -1,	isl_obj_union_map,	isl_obj_union_map,
		(isc_un_op_fn) &isl_union_map_domain_map } },
	{"ran",	{ -1,	isl_obj_union_map,	isl_obj_union_set,
		(isc_un_op_fn) &isl_union_map_range } },
	{"range",	{ -1,	isl_obj_union_map,	isl_obj_union_set,
		(isc_un_op_fn) &isl_union_map_range } },
	{"range_map",	{ -1,	isl_obj_union_map,	isl_obj_union_map,
		(isc_un_op_fn) &isl_union_map_range_map } },
	{"identity",	{ -1,	isl_obj_union_set,	isl_obj_union_map,
		(isc_un_op_fn) &isl_union_set_identity } },
	{"lexmin",	{ -1,	isl_obj_union_map,	isl_obj_union_map,
		(isc_un_op_fn) &isl_union_map_lexmin } },
	{"lexmax",	{ -1,	isl_obj_union_map,	isl_obj_union_map,
		(isc_un_op_fn) &isl_union_map_lexmax } },
	{"lexmin",	{ -1,	isl_obj_union_set,	isl_obj_union_set,
		(isc_un_op_fn) &isl_union_set_lexmin } },
	{"lexmax",	{ -1,	isl_obj_union_set,	isl_obj_union_set,
		(isc_un_op_fn) &isl_union_set_lexmax } },
	{"poly",	{ -1,	isl_obj_union_map,	isl_obj_union_map,
		(isc_un_op_fn) &isl_union_map_polyhedral_hull } },
	{"poly",	{ -1,	isl_obj_union_set,	isl_obj_union_set,
		(isc_un_op_fn) &isl_union_set_polyhedral_hull } },
	{"poly",	{ -1,	isl_obj_union_pw_qpolynomial,
		isl_obj_union_pw_qpolynomial,
		(isc_un_op_fn) &union_pw_qpolynomial_poly } },
	{"lpoly",	{ -1,	isl_obj_union_pw_qpolynomial,
		isl_obj_union_pw_qpolynomial,
		(isc_un_op_fn) &union_pw_qpolynomial_lpoly } },
	{"upoly",	{ -1,	isl_obj_union_pw_qpolynomial,
		isl_obj_union_pw_qpolynomial,
		(isc_un_op_fn) &union_pw_qpolynomial_upoly } },
	{"sample",	{ -1,	isl_obj_union_set,	isl_obj_set,
		(isc_un_op_fn) &union_set_sample } },
	{"sample",	{ -1,	isl_obj_union_map,	isl_obj_map,
		(isc_un_op_fn) &union_map_sample } },
	{"scan",	{ -1,	isl_obj_union_set,	isl_obj_union_set,
		(isc_un_op_fn) &union_set_scan } },
	{"scan",	{ -1,	isl_obj_union_map,	isl_obj_union_map,
		(isc_un_op_fn) &union_map_scan } },
	{"sum",		{ -1,	isl_obj_union_pw_qpolynomial,
		isl_obj_union_pw_qpolynomial,
		(isc_un_op_fn) &isl_union_pw_qpolynomial_sum } },
	{"ub",		{ -1,	isl_obj_union_pw_qpolynomial, isl_obj_list,
		(isc_un_op_fn) &union_pw_qpolynomial_upper_bound } },
	{"unwrap",	{ -1,	isl_obj_union_set,	isl_obj_union_map,
		(isc_un_op_fn) &isl_union_set_unwrap } },
	{"wrap",	{ -1,	isl_obj_union_map,	isl_obj_union_set,
		(isc_un_op_fn) &isl_union_map_wrap } },
	NULL
};

struct isl_named_obj {
	char		*name;
	struct isl_obj	obj;
};

static void free_obj(struct isl_obj obj)
{
	obj.type->free(obj.v);
}

static int same_name(const void *entry, const void *val)
{
	const struct isl_named_obj *named = (const struct isl_named_obj *)entry;

	return !strcmp(named->name, val);
}

static int do_assign(struct isl_ctx *ctx, struct isl_hash_table *table,
	char *name, struct isl_obj obj)
{
	struct isl_hash_table_entry *entry;
	uint32_t name_hash;
	struct isl_named_obj *named;

	name_hash = isl_hash_string(isl_hash_init(), name);
	entry = isl_hash_table_find(ctx, table, name_hash, same_name, name, 1);
	if (!entry)
		goto error;
	if (entry->data) {
		named = entry->data;
		free_obj(named->obj);
		free(name);
	} else {
		named = isl_alloc_type(ctx, struct isl_named_obj);
		if (!named)
			goto error;
		named->name = name;
		entry->data = named;
	}
	named->obj = obj;

	return 0;
error:
	free_obj(obj);
	free(name);
	return -1;
}

static struct isl_obj stored_obj(struct isl_ctx *ctx,
	struct isl_hash_table *table, char *name)
{
	struct isl_obj obj = { isl_obj_none, NULL };
	struct isl_hash_table_entry *entry;
	uint32_t name_hash;

	name_hash = isl_hash_string(isl_hash_init(), name);
	entry = isl_hash_table_find(ctx, table, name_hash, same_name, name, 0);
	if (entry) {
		struct isl_named_obj *named;
		named = entry->data;
		obj = named->obj;
	} else
		fprintf(stderr, "unknown identifier '%s'\n", name);

	free(name);
	obj.v = obj.type->copy(obj.v);
	return obj;
}

static int is_subtype(struct isl_obj obj, isl_obj_type super)
{
	if (obj.type == super)
		return 1;
	if (obj.type == isl_obj_map && super == isl_obj_union_map)
		return 1;
	if (obj.type == isl_obj_set && super == isl_obj_union_set)
		return 1;
	if (obj.type == isl_obj_pw_qpolynomial &&
	    super == isl_obj_union_pw_qpolynomial)
		return 1;
	if (obj.type == isl_obj_pw_qpolynomial_fold &&
	    super == isl_obj_union_pw_qpolynomial_fold)
		return 1;
	if (obj.type == isl_obj_union_set && isl_union_set_is_empty(obj.v))
		return 1;
	if (obj.type == isl_obj_list) {
		struct isl_list *list = obj.v;
		if (list->n == 2 && list->obj[1].type == isl_obj_bool)
			return is_subtype(list->obj[0], super);
	}
	if (super == isl_obj_str)
		return 1;
	return 0;
}

static struct isl_obj obj_at(struct isl_obj obj, int i)
{
	struct isl_list *list = obj.v;

	obj = list->obj[i];
	obj.v = obj.type->copy(obj.v);

	isl_list_free(list);

	return obj;
}

static struct isl_obj convert(isl_ctx *ctx, struct isl_obj obj,
	isl_obj_type type)
{
	if (obj.type == type)
		return obj;
	if (obj.type == isl_obj_map && type == isl_obj_union_map) {
		obj.type = isl_obj_union_map;
		obj.v = isl_union_map_from_map(obj.v);
		return obj;
	}
	if (obj.type == isl_obj_set && type == isl_obj_union_set) {
		obj.type = isl_obj_union_set;
		obj.v = isl_union_set_from_set(obj.v);
		return obj;
	}
	if (obj.type == isl_obj_pw_qpolynomial &&
	    type == isl_obj_union_pw_qpolynomial) {
		obj.type = isl_obj_union_pw_qpolynomial;
		obj.v = isl_union_pw_qpolynomial_from_pw_qpolynomial(obj.v);
		return obj;
	}
	if (obj.type == isl_obj_pw_qpolynomial_fold &&
	    type == isl_obj_union_pw_qpolynomial_fold) {
		obj.type = isl_obj_union_pw_qpolynomial_fold;
		obj.v = isl_union_pw_qpolynomial_fold_from_pw_qpolynomial_fold(obj.v);
		return obj;
	}
	if (obj.type == isl_obj_union_set && isl_union_set_is_empty(obj.v)) {
		if (type == isl_obj_union_map) {
			obj.type = isl_obj_union_map;
			return obj;
		}
		if (type == isl_obj_union_pw_qpolynomial) {
			isl_dim *dim = isl_union_set_get_dim(obj.v);
			isl_union_set_free(obj.v);
			obj.v = isl_union_pw_qpolynomial_zero(dim);
			obj.type = isl_obj_union_pw_qpolynomial;
			return obj;
		}
		if (type == isl_obj_union_pw_qpolynomial_fold) {
			isl_dim *dim = isl_union_set_get_dim(obj.v);
			isl_union_set_free(obj.v);
			obj.v = isl_union_pw_qpolynomial_fold_zero(dim,
								isl_fold_list);
			obj.type = isl_obj_union_pw_qpolynomial_fold;
			return obj;
		}
	}
	if (obj.type == isl_obj_list) {
		struct isl_list *list = obj.v;
		if (list->n == 2 && list->obj[1].type == isl_obj_bool)
			return convert(ctx, obj_at(obj, 0), type);
	}
	if (type == isl_obj_str) {
		isl_str *str;
		isl_printer *p;
		char *s;

		p = isl_printer_to_str(ctx);
		if (!p)
			goto error;
		p = obj.type->print(p, obj.v);
		s = isl_printer_get_str(p);
		isl_printer_free(p);

		str = isl_str_from_string(ctx, s);
		if (!str)
			goto error;
		free_obj(obj);
		obj.v = str;
		obj.type = isl_obj_str;
		return obj;

	}
error:
	free_obj(obj);
	obj.type = isl_obj_none;
	obj.v = NULL;
	return obj;
}

static struct isc_bin_op *read_bin_op_if_available(struct isl_stream *s,
	struct isl_obj lhs)
{
	int i;
	struct isl_token *tok;

	tok = isl_stream_next_token(s);
	if (!tok)
		return NULL;

	for (i = 0; ; ++i) {
		if (!bin_ops[i].op)
			break;
		if (bin_ops[i].op != tok->type)
			continue;
		if (!is_subtype(lhs, bin_ops[i].lhs))
			continue;

		isl_token_free(tok);
		return &bin_ops[i];
	}

	for (i = 0; ; ++i) {
		if (!named_bin_ops[i].name)
			break;
		if (named_bin_ops[i].op.op != tok->type)
			continue;
		if (!is_subtype(lhs, named_bin_ops[i].op.lhs))
			continue;

		isl_token_free(tok);
		return &named_bin_ops[i].op;
	}

	isl_stream_push_token(s, tok);

	return NULL;
}

static struct isc_un_op *read_prefix_un_op_if_available(struct isl_stream *s)
{
	int i;
	struct isl_token *tok;

	tok = isl_stream_next_token(s);
	if (!tok)
		return NULL;

	for (i = 0; ; ++i) {
		if (!named_un_ops[i].name)
			break;
		if (named_un_ops[i].op.op != tok->type)
			continue;

		isl_token_free(tok);
		return &named_un_ops[i].op;
	}

	isl_stream_push_token(s, tok);

	return NULL;
}

static struct isc_un_op *find_matching_un_op(struct isc_un_op *like,
	struct isl_obj arg)
{
	int i;

	for (i = 0; ; ++i) {
		if (!named_un_ops[i].name)
			break;
		if (named_un_ops[i].op.op != like->op)
			continue;
		if (!is_subtype(arg, named_un_ops[i].op.arg))
			continue;

		return &named_un_ops[i].op;
	}

	return NULL;
}

static int is_assign(struct isl_stream *s)
{
	struct isl_token *tok;
	struct isl_token *tok2;
	int assign;

	tok = isl_stream_next_token(s);
	if (!tok)
		return 0;
	if (tok->type != ISL_TOKEN_IDENT) {
		isl_stream_push_token(s, tok);
		return 0;
	}

	tok2 = isl_stream_next_token(s);
	if (!tok2) {
		isl_stream_push_token(s, tok);
		return 0;
	}
	assign = tok2->type == ISL_TOKEN_DEF;
	isl_stream_push_token(s, tok2);
	isl_stream_push_token(s, tok);

	return assign;
}

static struct isl_obj read_obj(struct isl_stream *s,
	struct isl_hash_table *table);
static struct isl_obj read_expr(struct isl_stream *s,
	struct isl_hash_table *table);

static struct isl_obj read_un_op_expr(struct isl_stream *s,
	struct isl_hash_table *table, struct isc_un_op *op)
{
	struct isl_obj obj = { isl_obj_none, NULL };

	obj = read_obj(s, table);
	if (!obj.v)
		goto error;

	op = find_matching_un_op(op, obj);

	if (!op)
		isl_die(s->ctx, isl_error_invalid,
			"no such unary operator defined on given operand",
			goto error);

	obj = convert(s->ctx, obj, op->arg);
	obj.v = op->fn(obj.v);
	obj.type = op->res;

	return obj;
error:
	free_obj(obj);
	obj.type = isl_obj_none;
	obj.v = NULL;
	return obj;
}

static struct isl_obj transitive_closure(struct isl_ctx *ctx, struct isl_obj obj)
{
	struct isl_list *list;
	int exact;

	if (obj.type != isl_obj_union_map)
		obj = convert(ctx, obj, isl_obj_union_map);
	isl_assert(ctx, obj.type == isl_obj_union_map, goto error);
	list = isl_list_alloc(ctx, 2);
	if (!list)
		goto error;

	list->obj[0].type = isl_obj_union_map;
	list->obj[0].v = isl_union_map_transitive_closure(obj.v, &exact);
	list->obj[1].type = isl_obj_bool;
	list->obj[1].v = exact ? &isl_bool_true : &isl_bool_false;
	obj.v = list;
	obj.type = isl_obj_list;
	if (exact < 0 || !list->obj[0].v)
		goto error;

	return obj;
error:
	free_obj(obj);
	obj.type = isl_obj_none;
	obj.v = NULL;
	return obj;
}

static struct isl_obj obj_at_index(struct isl_stream *s, struct isl_obj obj)
{
	struct isl_list *list = obj.v;
	struct isl_token *tok;
	int i;

	tok = isl_stream_next_token(s);
	if (!tok || tok->type != ISL_TOKEN_VALUE) {
		isl_stream_error(s, tok, "expecting index");
		if (tok)
			isl_stream_push_token(s, tok);
		goto error;
	}
	i = isl_int_get_si(tok->u.v);
	isl_token_free(tok);
	isl_assert(s->ctx, i < list->n, goto error);
	if (isl_stream_eat(s, ']'))
		goto error;

	return obj_at(obj, i);
error:
	free_obj(obj);
	obj.type = isl_obj_none;
	obj.v = NULL;
	return obj;
}

static struct isl_obj apply(struct isl_stream *s, __isl_take isl_union_map *umap,
	struct isl_hash_table *table)
{
	struct isl_obj obj;

	obj = read_expr(s, table);
	isl_assert(s->ctx, is_subtype(obj, isl_obj_union_set) ||
			   is_subtype(obj, isl_obj_union_map), goto error);

	if (obj.type == isl_obj_list) {
		struct isl_list *list = obj.v;
		if (list->n == 2 && list->obj[1].type == isl_obj_bool)
			obj = obj_at(obj, 0);
	}
	if (obj.type == isl_obj_set)
		obj = convert(s->ctx, obj, isl_obj_union_set);
	else if (obj.type == isl_obj_map)
		obj = convert(s->ctx, obj, isl_obj_union_map);
	if (obj.type == isl_obj_union_set) {
		obj.v = isl_union_set_apply(obj.v, umap);
	} else
		obj.v = isl_union_map_apply_range(obj.v, umap);
	if (!obj.v)
		goto error2;

	if (isl_stream_eat(s, ')'))
		goto error2;

	return obj;
error:
	isl_union_map_free(umap);
error2:
	free_obj(obj);
	obj.type = isl_obj_none;
	obj.v = NULL;
	return obj;
}

static struct isl_obj apply_fun(struct isl_stream *s,
	struct isl_obj obj, struct isl_hash_table *table)
{
	struct isl_obj arg;

	arg = read_expr(s, table);
	if (!is_subtype(arg, isl_obj_union_map) &&
	    !is_subtype(arg, isl_obj_union_set))
		isl_die(s->ctx, isl_error_invalid,
			"expecting set of map argument", goto error);

	if (arg.type == isl_obj_list) {
		struct isl_list *list = arg.v;
		if (list->n == 2 && list->obj[1].type == isl_obj_bool)
			arg = obj_at(arg, 0);
	}
	if (arg.type == isl_obj_set)
		arg = convert(s->ctx, arg, isl_obj_union_set);
	else if (arg.type == isl_obj_map)
		arg = convert(s->ctx, arg, isl_obj_union_map);
	if (arg.type == isl_obj_union_set) {
		arg.v = isl_union_map_from_range(arg.v);
		arg.type = isl_obj_union_map;
	}
	if (obj.type == isl_obj_union_pw_qpolynomial) {
		obj.v = isl_union_map_apply_union_pw_qpolynomial(arg.v, obj.v);
	} else {
		obj.type = isl_obj_list;
		obj.v = union_map_apply_union_pw_qpolynomial_fold(arg.v, obj.v);
	}
	if (!obj.v)
		goto error2;

	if (isl_stream_eat(s, ')'))
		goto error2;

	return obj;
error:
	free_obj(arg);
error2:
	free_obj(obj);
	obj.type = isl_obj_none;
	obj.v = NULL;
	return obj;
}

struct add_vertex_data {
	struct isl_list *list;
	int i;
};

static int add_vertex(__isl_take isl_vertex *vertex, void *user)
{
	struct add_vertex_data *data = (struct add_vertex_data *)user;
	isl_basic_set *expr;

	expr = isl_vertex_get_expr(vertex);

	data->list->obj[data->i].type = isl_obj_set;
	data->list->obj[data->i].v = isl_set_from_basic_set(expr);
	data->i++;

	isl_vertex_free(vertex);

	return 0;
}

static int set_vertices(__isl_take isl_set *set, void *user)
{
	isl_ctx *ctx;
	isl_basic_set *hull;
	isl_vertices *vertices = NULL;
	struct isl_list *list = NULL;
	int r;
	struct add_vertex_data *data = (struct add_vertex_data *)user;

	set = isl_set_remove_divs(set);
	hull = isl_set_convex_hull(set);
	vertices = isl_basic_set_compute_vertices(hull);
	isl_basic_set_free(hull);

	list = data->list;

	ctx = isl_vertices_get_ctx(vertices);
	data->list = isl_list_alloc(ctx, isl_vertices_get_n_vertices(vertices));
	if (!data->list)
		goto error;

	data->i = 0;
	r = isl_vertices_foreach_vertex(vertices, &add_vertex, user);

	data->list = isl_list_concat(list, data->list);

	isl_vertices_free(vertices);

	return r;
error:
	data->list = list;
	isl_vertices_free(vertices);
	return -1;
}

static struct isl_obj vertices(struct isl_stream *s,
	struct isl_hash_table *table)
{
	isl_ctx *ctx;
	struct isl_obj obj;
	struct isl_list *list = NULL;
	isl_union_set *uset;
	struct add_vertex_data data = { NULL };

	obj = read_expr(s, table);
	obj = convert(s->ctx, obj, isl_obj_union_set);
	isl_assert(s->ctx, obj.type == isl_obj_union_set, goto error);
	uset = obj.v;
	obj.v = NULL;

	ctx = isl_union_set_get_ctx(uset);
	list = isl_list_alloc(ctx, 0);
	if (!list)
		goto error;

	data.list = list;

	if (isl_union_set_foreach_set(uset, &set_vertices, &data) < 0)
		goto error;

	isl_union_set_free(uset);

	obj.type = isl_obj_list;
	obj.v = data.list;

	return obj;
error:
	isl_union_set_free(uset);
	isl_list_free(data.list);
	free_obj(obj);
	obj.type = isl_obj_none;
	obj.v = NULL;
	return obj;
}

static struct isl_obj type_of(struct isl_stream *s,
	struct isl_hash_table *table)
{
	isl_ctx *ctx;
	struct isl_obj obj;
	const char *type = "unknown";

	obj = read_expr(s, table);

	if (obj.type == isl_obj_map ||
	    obj.type == isl_obj_union_map)
		type = "map";
	if (obj.type == isl_obj_set ||
	    obj.type == isl_obj_union_set)
		type = "set";
	if (obj.type == isl_obj_pw_qpolynomial ||
	    obj.type == isl_obj_union_pw_qpolynomial)
		type = "piecewise quasipolynomial";
	if (obj.type == isl_obj_pw_qpolynomial_fold ||
	    obj.type == isl_obj_union_pw_qpolynomial_fold)
		type = "piecewise quasipolynomial fold";
	if (obj.type == isl_obj_list)
		type = "list";
	if (obj.type == isl_obj_bool)
		type = "boolean";
	if (obj.type == isl_obj_str)
		type = "string";
	if (obj.type == isl_obj_int)
		type = "int";

	free_obj(obj);
	obj.type = isl_obj_str;
	obj.v = isl_str_from_string(s->ctx, strdup(type));

	return obj;
}

static __isl_give isl_union_map *read_map(struct isl_stream *s,
	struct isl_hash_table *table)
{
	struct isl_obj obj;

	obj = read_obj(s, table);
	obj = convert(s->ctx, obj, isl_obj_union_map);
	isl_assert(s->ctx, obj.type == isl_obj_union_map, goto error);
	return obj.v;
error:
	free_obj(obj);
	return NULL;
}

static struct isl_obj last_any(struct isl_stream *s,
	struct isl_hash_table *table, __isl_take isl_union_map *must_source,
	__isl_take isl_union_map *may_source)
{
	struct isl_obj obj = { isl_obj_none, NULL };
	isl_union_map *sink = NULL;
	isl_union_map *schedule = NULL;
	isl_union_map *may_dep;
	isl_union_map *must_dep;

	if (isl_stream_eat(s, iscc_op[ISCC_BEFORE]))
		goto error;

	sink = read_map(s, table);
	if (!sink)
		goto error;

	if (isl_stream_eat(s, iscc_op[ISCC_UNDER]))
		goto error;

	schedule = read_map(s, table);
	if (!schedule)
		goto error;

	if (isl_union_map_compute_flow(sink, must_source, may_source,
				       schedule, &must_dep, &may_dep,
				       NULL, NULL) < 0)
		return obj;

	obj.type = isl_obj_union_map;
	obj.v = isl_union_map_union(must_dep, may_dep);

	return obj;
error:
	isl_union_map_free(may_source);
	isl_union_map_free(must_source);
	isl_union_map_free(sink);
	isl_union_map_free(schedule);
	free_obj(obj);
	obj.type = isl_obj_none;
	obj.v = NULL;
	return obj;
}

static struct isl_obj any(struct isl_stream *s, struct isl_hash_table *table)
{
	struct isl_obj obj = { isl_obj_none, NULL };
	isl_union_map *must_source = NULL;
	isl_union_map *may_source = NULL;
	isl_union_map *sink = NULL;
	isl_union_map *schedule = NULL;
	isl_union_map *may_dep;

	may_source = read_map(s, table);
	if (!may_source)
		goto error;

	if (isl_stream_eat_if_available(s, iscc_op[ISCC_LAST])) {
		must_source = read_map(s, table);
		if (!must_source)
			goto error;
		return last_any(s, table, must_source, may_source);
	}

	if (isl_stream_eat(s, iscc_op[ISCC_BEFORE]))
		goto error;

	sink = read_map(s, table);
	if (!sink)
		goto error;

	if (isl_stream_eat(s, iscc_op[ISCC_UNDER]))
		goto error;

	schedule = read_map(s, table);
	if (!schedule)
		goto error;

	must_source = isl_union_map_empty(isl_union_map_get_dim(sink));
	if (isl_union_map_compute_flow(sink, must_source, may_source,
				       schedule, NULL, &may_dep,
				       NULL, NULL) < 0)
		return obj;

	obj.type = isl_obj_union_map;
	obj.v = may_dep;

	return obj;
error:
	isl_union_map_free(may_source);
	isl_union_map_free(must_source);
	isl_union_map_free(sink);
	isl_union_map_free(schedule);
	free_obj(obj);
	obj.type = isl_obj_none;
	obj.v = NULL;
	return obj;
}

static struct isl_obj last(struct isl_stream *s, struct isl_hash_table *table)
{
	struct isl_obj obj = { isl_obj_none, NULL };
	struct isl_list *list = NULL;
	isl_union_map *must_source = NULL;
	isl_union_map *may_source = NULL;
	isl_union_map *sink = NULL;
	isl_union_map *schedule = NULL;
	isl_union_map *must_dep;
	isl_union_map *must_no_source;

	must_source = read_map(s, table);
	if (!must_source)
		goto error;

	if (isl_stream_eat_if_available(s, iscc_op[ISCC_ANY])) {
		may_source = read_map(s, table);
		if (!may_source)
			goto error;
		return last_any(s, table, must_source, may_source);
	}

	list = isl_list_alloc(s->ctx, 2);
	if (!list)
		goto error;

	if (isl_stream_eat(s, iscc_op[ISCC_BEFORE]))
		goto error;

	sink = read_map(s, table);
	if (!sink)
		goto error;

	if (isl_stream_eat(s, iscc_op[ISCC_UNDER]))
		goto error;

	schedule = read_map(s, table);
	if (!schedule)
		goto error;

	may_source = isl_union_map_empty(isl_union_map_get_dim(sink));
	if (isl_union_map_compute_flow(sink, must_source, may_source,
				       schedule, &must_dep, NULL,
				       &must_no_source, NULL) < 0)
		return obj;

	list->obj[0].type = isl_obj_union_map;
	list->obj[0].v = must_dep;
	list->obj[1].type = isl_obj_union_map;
	list->obj[1].v = must_no_source;

	obj.v = list;
	obj.type = isl_obj_list;

	return obj;
error:
	isl_list_free(list);
	isl_union_map_free(may_source);
	isl_union_map_free(must_source);
	isl_union_map_free(sink);
	isl_union_map_free(schedule);
	free_obj(obj);
	obj.type = isl_obj_none;
	obj.v = NULL;
	return obj;
}

static struct isl_obj power(struct isl_stream *s, struct isl_obj obj)
{
	struct isl_token *tok;

	if (isl_stream_eat_if_available(s, '+'))
		return transitive_closure(s->ctx, obj);

	tok = isl_stream_next_token(s);
	if (!tok || tok->type != ISL_TOKEN_VALUE || isl_int_cmp_si(tok->u.v, -1)) {
		isl_stream_error(s, tok, "expecting -1");
		if (tok)
			isl_stream_push_token(s, tok);
		goto error;
	}
	isl_token_free(tok);
	isl_assert(s->ctx, is_subtype(obj, isl_obj_union_map), goto error);
	if (obj.type != isl_obj_union_map)
		obj = convert(s->ctx, obj, isl_obj_union_map);

	obj.v = isl_union_map_reverse(obj.v);
	if (!obj.v)
		goto error;

	return obj;
error:
	free_obj(obj);
	obj.type = isl_obj_none;
	obj.v = NULL;
	return obj;
}

static struct isl_obj read_from_file(struct isl_stream *s)
{
	struct isl_obj obj;
	struct isl_token *tok;
	struct isl_stream *s_file;
	struct iscc_options *options;
	FILE *file;

	tok = isl_stream_next_token(s);
	if (!tok || tok->type != ISL_TOKEN_STRING) {
		isl_stream_error(s, tok, "expecting filename");
		isl_token_free(tok);
		goto error;
	}

	options = isl_ctx_peek_iscc_options(s->ctx);
	if (!options || !options->io) {
		isl_token_free(tok);
		isl_die(s->ctx, isl_error_invalid,
			"read operation not allowed", goto error);
	}

	file = fopen(tok->u.s, "r");
	isl_token_free(tok);
	isl_assert(s->ctx, file, goto error);

	s_file = isl_stream_new_file(s->ctx, file);
	if (!s_file) {
		fclose(file);
		goto error;
	}

	obj = isl_stream_read_obj(s_file);

	isl_stream_free(s_file);
	fclose(file);

	return obj;
error:
	obj.type = isl_obj_none;
	obj.v = NULL;
	return obj;
}

static struct isl_obj write_to_file(struct isl_stream *s,
	struct isl_hash_table *table)
{
	struct isl_obj obj;
	struct isl_token *tok;
	struct isl_stream *s_file;
	struct iscc_options *options;
	FILE *file;
	isl_printer *p;

	tok = isl_stream_next_token(s);
	if (!tok || tok->type != ISL_TOKEN_STRING) {
		isl_stream_error(s, tok, "expecting filename");
		isl_token_free(tok);
		goto error;
	}

	obj = read_expr(s, table);

	options = isl_ctx_peek_iscc_options(s->ctx);
	if (!options || !options->io) {
		isl_token_free(tok);
		isl_die(s->ctx, isl_error_invalid,
			"write operation not allowed", goto error);
	}

	file = fopen(tok->u.s, "w");
	isl_token_free(tok);
	if (!file)
		isl_die(s->ctx, isl_error_unknown,
			"could not open file for writing", goto error);

	p = isl_printer_to_file(s->ctx, file);
	p = isl_printer_set_output_format(p, options->format);
	p = obj.type->print(p, obj.v);
	p = isl_printer_end_line(p);
	isl_printer_free(p);

	fclose(file);
error:
	free_obj(obj);
	obj.type = isl_obj_none;
	obj.v = NULL;
	return obj;
}

static struct isl_obj read_string_if_available(struct isl_stream *s)
{
	struct isl_token *tok;
	struct isl_obj obj = { isl_obj_none, NULL };

	tok = isl_stream_next_token(s);
	if (!tok)
		return obj;
	if (tok->type == ISL_TOKEN_STRING) {
		isl_str *str;
		str = isl_str_alloc(s->ctx);
		if (!str)
			goto error;
		str->s = strdup(tok->u.s);
		isl_token_free(tok);
		obj.v = str;
		obj.type = isl_obj_str;
	} else
		isl_stream_push_token(s, tok);
	return obj;
error:
	isl_token_free(tok);
	return obj;
}

static struct isl_obj read_obj(struct isl_stream *s,
	struct isl_hash_table *table)
{
	struct isl_obj obj = { isl_obj_none, NULL };
	char *name = NULL;
	struct isc_un_op *op = NULL;

	obj = read_string_if_available(s);
	if (obj.v)
		return obj;
	if (isl_stream_eat_if_available(s, '(')) {
		obj = read_expr(s, table);
		if (!obj.v || isl_stream_eat(s, ')'))
			goto error;
	} else {
		op = read_prefix_un_op_if_available(s);
		if (op)
			return read_un_op_expr(s, table, op);

		if (isl_stream_eat_if_available(s, iscc_op[ISCC_READ]))
			return read_from_file(s);
		if (isl_stream_eat_if_available(s, iscc_op[ISCC_WRITE]))
			return write_to_file(s, table);
		if (isl_stream_eat_if_available(s, iscc_op[ISCC_VERTICES]))
			return vertices(s, table);
		if (isl_stream_eat_if_available(s, iscc_op[ISCC_ANY]))
			return any(s, table);
		if (isl_stream_eat_if_available(s, iscc_op[ISCC_LAST]))
			return last(s, table);
		if (isl_stream_eat_if_available(s, iscc_op[ISCC_TYPEOF]))
			return type_of(s, table);

		name = isl_stream_read_ident_if_available(s);
		if (name)
			obj = stored_obj(s->ctx, table, name);
		else
			obj = isl_stream_read_obj(s);
		if (!obj.v)
			goto error;
	}

	if (isl_stream_eat_if_available(s, '^'))
		obj = power(s, obj);
	else if (obj.type == isl_obj_list && isl_stream_eat_if_available(s, '['))
		obj = obj_at_index(s, obj);
	else if (is_subtype(obj, isl_obj_union_map) &&
		 isl_stream_eat_if_available(s, '(')) {
		obj = convert(s->ctx, obj, isl_obj_union_map);
		obj = apply(s, obj.v, table);
	} else if (is_subtype(obj, isl_obj_union_pw_qpolynomial) &&
		   isl_stream_eat_if_available(s, '(')) {
		obj = convert(s->ctx, obj, isl_obj_union_pw_qpolynomial);
		obj = apply_fun(s, obj, table);
	} else if (is_subtype(obj, isl_obj_union_pw_qpolynomial_fold) &&
		   isl_stream_eat_if_available(s, '(')) {
		obj = convert(s->ctx, obj, isl_obj_union_pw_qpolynomial_fold);
		obj = apply_fun(s, obj, table);
	}

	return obj;
error:
	free_obj(obj);
	obj.type = isl_obj_none;
	obj.v = NULL;
	return obj;
}

static struct isc_bin_op *find_matching_bin_op(struct isc_bin_op *like,
	struct isl_obj lhs, struct isl_obj rhs)
{
	int i;

	for (i = 0; ; ++i) {
		if (!bin_ops[i].op)
			break;
		if (bin_ops[i].op != like->op)
			continue;
		if (!is_subtype(lhs, bin_ops[i].lhs))
			continue;
		if (!is_subtype(rhs, bin_ops[i].rhs))
			continue;

		return &bin_ops[i];
	}

	for (i = 0; ; ++i) {
		if (!named_bin_ops[i].name)
			break;
		if (named_bin_ops[i].op.op != like->op)
			continue;
		if (!is_subtype(lhs, named_bin_ops[i].op.lhs))
			continue;
		if (!is_subtype(rhs, named_bin_ops[i].op.rhs))
			continue;

		return &named_bin_ops[i].op;
	}

	return NULL;
}

static int next_is_neg_int(struct isl_stream *s)
{
	struct isl_token *tok;
	int ret;

	tok = isl_stream_next_token(s);
	ret = tok && tok->type == ISL_TOKEN_VALUE && isl_int_is_neg(tok->u.v);
	isl_stream_push_token(s, tok);

	return ret;
}

static struct isl_obj read_expr(struct isl_stream *s,
	struct isl_hash_table *table)
{
	struct isl_obj obj = { isl_obj_none, NULL };
	struct isl_obj right_obj = { isl_obj_none, NULL };

	obj = read_obj(s, table);
	for (; obj.v;) {
		struct isc_bin_op *op = NULL;

		op = read_bin_op_if_available(s, obj);
		if (!op)
			break;

		right_obj = read_obj(s, table);

		op = find_matching_bin_op(op, obj, right_obj);

		if (!op)
			isl_die(s->ctx, isl_error_invalid,
			    "no such binary operator defined on given operands",
			    goto error);

		obj = convert(s->ctx, obj, op->lhs);
		right_obj = convert(s->ctx, right_obj, op->rhs);
		obj.v = op->fn(obj.v, right_obj.v);
		obj.type = op->res;
	}

	if (obj.type == isl_obj_int && next_is_neg_int(s)) {
		right_obj = read_obj(s, table);
		obj.v = isl_int_obj_add(obj.v, right_obj.v);
	}

	return obj;
error:
	free_obj(right_obj);
	free_obj(obj);
	obj.type = isl_obj_none;
	obj.v = NULL;
	return obj;
}

static __isl_give isl_printer *source_file(struct isl_stream *s,
	struct isl_hash_table *table, __isl_take isl_printer *p);

static __isl_give isl_printer *read_line(struct isl_stream *s,
	struct isl_hash_table *table, __isl_take isl_printer *p)
{
	struct isl_obj obj = { isl_obj_none, NULL };
	char *lhs = NULL;
	int assign = 0;
	struct isc_bin_op *op = NULL;

	if (!p)
		return NULL;
	if (isl_stream_is_empty(s))
		return p;

	if (isl_stream_eat_if_available(s, iscc_op[ISCC_SOURCE]))
		return source_file(s, table, p);

	assign = is_assign(s);
	if (assign) {
		lhs = isl_stream_read_ident_if_available(s);
		if (isl_stream_eat(s, ISL_TOKEN_DEF))
			goto error;
	}

	obj = read_expr(s, table);
	if (obj.type == isl_obj_none || obj.v == NULL)
		goto error;
	if (isl_stream_eat(s, ';'))
		goto error;

	if (assign) {
		if (do_assign(s->ctx, table, lhs, obj))
			return p;
	} else {
		p = obj.type->print(p, obj.v);
		p = isl_printer_end_line(p);
		free_obj(obj);
	}

	return p;
error:
	isl_stream_flush_tokens(s);
	isl_stream_skip_line(s);
	free(lhs);
	free_obj(obj);
	return p;
}

int free_cb(void **entry, void *user)
{
	struct isl_named_obj *named = *entry;

	free_obj(named->obj);
	free(named->name);
	free(named);

	return 0;
}

static void register_named_ops(struct isl_stream *s)
{
	int i;

	for (i = 0; i < ISCC_N_OP; ++i) {
		iscc_op[i] = isl_stream_register_keyword(s, op_name[i]);
		assert(iscc_op[i] != ISL_TOKEN_ERROR);
	}

	for (i = 0; ; ++i) {
		if (!named_un_ops[i].name)
			break;
		named_un_ops[i].op.op = isl_stream_register_keyword(s,
							named_un_ops[i].name);
		assert(named_un_ops[i].op.op != ISL_TOKEN_ERROR);
	}

	for (i = 0; ; ++i) {
		if (!named_bin_ops[i].name)
			break;
		named_bin_ops[i].op.op = isl_stream_register_keyword(s,
							named_bin_ops[i].name);
		assert(named_bin_ops[i].op.op != ISL_TOKEN_ERROR);
	}
}

static __isl_give isl_printer *source_file(struct isl_stream *s,
	struct isl_hash_table *table, __isl_take isl_printer *p)
{
	struct isl_token *tok;
	struct isl_stream *s_file;
	FILE *file;

	tok = isl_stream_next_token(s);
	if (!tok || tok->type != ISL_TOKEN_STRING) {
		isl_stream_error(s, tok, "expecting filename");
		isl_token_free(tok);
		return p;
	}

	file = fopen(tok->u.s, "r");
	isl_token_free(tok);
	isl_assert(s->ctx, file, return p);

	s_file = isl_stream_new_file(s->ctx, file);
	if (!s_file) {
		fclose(file);
		return p;
	}

	register_named_ops(s_file);

	while (!s_file->eof)
		p = read_line(s_file, table, p);

	isl_stream_free(s_file);
	fclose(file);

	isl_stream_eat(s, ';');

	return p;
}

int main(int argc, char **argv)
{
	struct isl_ctx *ctx;
	struct isl_stream *s;
	struct isl_hash_table *table;
	struct iscc_options *options;
	isl_printer *p;

	options = iscc_options_new_with_defaults();
	assert(options);
	argc = iscc_options_parse(options, argc, argv, ISL_ARG_ALL);

	ctx = isl_ctx_alloc_with_options(iscc_options_arg, options);
	s = isl_stream_new_file(ctx, stdin);
	assert(s);
	table = isl_hash_table_alloc(ctx, 10);
	assert(table);
	p = isl_printer_to_file(ctx, stdout);
	p = isl_printer_set_output_format(p, options->format);
	assert(p);

	register_named_ops(s);

	while (p && !s->eof) {
		p = read_line(s, table, p);
	}

	isl_printer_free(p);
	isl_hash_table_foreach(ctx, table, free_cb, NULL);
	isl_hash_table_free(ctx, table);
	isl_stream_free(s);
	isl_ctx_free(ctx);

	return 0;
}
