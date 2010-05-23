#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <isl_obj.h>
#include <isl_stream.h>
#include <isl_obj_list.h>
#include <barvinok/barvinok.h>
#include <bound_common.h>

#include "config.h"

static int isl_bool_false = 0;
static int isl_bool_true = 1;
static int isl_bool_error = -1;

static enum isl_token_type read_op;

struct isl_arg_choice iscc_format[] = {
	{"isl",		ISL_FORMAT_ISL},
	{"omega",	ISL_FORMAT_OMEGA},
	{"polylib",	ISL_FORMAT_POLYLIB},
	{"latex",	ISL_FORMAT_LATEX},
	{"C",		ISL_FORMAT_C},
	{0}
};

struct iscc_options {
	struct barvinok_options	*barvinok;
	unsigned		 format;
};

struct isl_arg iscc_options_arg[] = {
ISL_ARG_CHILD(struct iscc_options, barvinok, "barvinok", barvinok_options_arg,
	"barvinok options")
ISL_ARG_CHOICE(struct iscc_options, format, 0, "format", \
	iscc_format,	ISL_FORMAT_ISL, "output format")
ISL_ARG_END
};

ISL_ARG_DEF(iscc_options, struct iscc_options, iscc_options_arg)

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

int *map_is_equal(__isl_take isl_map *map1, __isl_take isl_map *map2)
{
	int res = isl_map_is_equal(map1, map2);
	isl_map_free(map1);
	isl_map_free(map2);
	return isl_bool_from_int(res);
}
int *set_is_equal(__isl_take isl_set *set1, __isl_take isl_set *set2)
{
	return map_is_equal((isl_map *)set1, (isl_map *)set2);
}

int *map_is_subset(__isl_take isl_map *map1, __isl_take isl_map *map2)
{
	int res = isl_map_is_subset(map1, map2);
	isl_map_free(map1);
	isl_map_free(map2);
	return isl_bool_from_int(res);
}
int *set_is_subset(__isl_take isl_set *set1, __isl_take isl_set *set2)
{
	return map_is_subset((isl_map *)set1, (isl_map *)set2);
}

int *map_is_strict_subset(__isl_take isl_map *map1, __isl_take isl_map *map2)
{
	int res = isl_map_is_strict_subset(map1, map2);
	isl_map_free(map1);
	isl_map_free(map2);
	return isl_bool_from_int(res);
}
int *set_is_strict_subset(__isl_take isl_set *set1, __isl_take isl_set *set2)
{
	return map_is_strict_subset((isl_map *)set1, (isl_map *)set2);
}

int *map_is_superset(__isl_take isl_map *map1, __isl_take isl_map *map2)
{
	return map_is_subset(map2, map1);
}
int *set_is_superset(__isl_take isl_set *set1, __isl_take isl_set *set2)
{
	return set_is_subset(set2, set1);
}

int *map_is_strict_superset(__isl_take isl_map *map1, __isl_take isl_map *map2)
{
	return map_is_strict_subset(map2, map1);
}
int *set_is_strict_superset(__isl_take isl_set *set1, __isl_take isl_set *set2)
{
	return set_is_strict_subset(set2, set1);
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
	isl_pw_qpolynomial *pwqp;
	isl_pw_qpolynomial *res;
};

static int eval_at(__isl_take isl_point *pnt, void *user)
{
	struct iscc_at *at = (struct iscc_at *) user;
	isl_qpolynomial *qp;
	isl_set *set;

	set = isl_set_from_point(isl_point_copy(pnt));
	qp = isl_pw_qpolynomial_eval(isl_pw_qpolynomial_copy(at->pwqp), pnt);

	at->res = isl_pw_qpolynomial_add_disjoint(at->res,
			isl_pw_qpolynomial_alloc(set, qp));

	return 0;
}

__isl_give isl_pw_qpolynomial *isl_pw_qpolynomial_at(
	__isl_take isl_pw_qpolynomial *pwqp, __isl_take isl_set *set)
{
	struct iscc_at at;

	at.pwqp = pwqp;
	at.res = isl_pw_qpolynomial_zero(isl_set_get_dim(set));

	isl_set_foreach_point(set, eval_at, &at);

	isl_pw_qpolynomial_free(pwqp);
	isl_set_free(set);

	return at.res;
}

struct iscc_fold_at {
	isl_pw_qpolynomial_fold *pwf;
	isl_pw_qpolynomial *res;
};

static int eval_fold_at(__isl_take isl_point *pnt, void *user)
{
	struct iscc_fold_at *at = (struct iscc_fold_at *) user;
	isl_qpolynomial *qp;
	isl_set *set;

	set = isl_set_from_point(isl_point_copy(pnt));
	qp = isl_pw_qpolynomial_fold_eval(isl_pw_qpolynomial_fold_copy(at->pwf),
						pnt);

	at->res = isl_pw_qpolynomial_add_disjoint(at->res,
			isl_pw_qpolynomial_alloc(set, qp));

	return 0;
}

__isl_give isl_pw_qpolynomial *isl_pw_qpolynomial_fold_at(
	__isl_take isl_pw_qpolynomial_fold *pwf, __isl_take isl_set *set)
{
	struct iscc_fold_at at;

	at.pwf = pwf;
	at.res = isl_pw_qpolynomial_zero(isl_set_get_dim(set));

	isl_set_foreach_point(set, eval_fold_at, &at);

	isl_pw_qpolynomial_fold_free(pwf);
	isl_set_free(set);

	return at.res;
}

struct isc_bin_op bin_ops[] = {
	{ '+',	isl_obj_set,	isl_obj_set,
		isl_obj_set,
		(isc_bin_op_fn) &isl_set_union },
	{ '+',	isl_obj_map,	isl_obj_map,
		isl_obj_map,
		(isc_bin_op_fn) &isl_map_union },
	{ '-',	isl_obj_set,	isl_obj_set,
		isl_obj_set,
		(isc_bin_op_fn) &isl_set_subtract },
	{ '-',	isl_obj_map,	isl_obj_map,
		isl_obj_map,
		(isc_bin_op_fn) &isl_map_subtract },
	{ '*',	isl_obj_set,	isl_obj_set,
		isl_obj_set,
		(isc_bin_op_fn) &isl_set_intersect },
	{ '*',	isl_obj_map,	isl_obj_map,
		isl_obj_map,
		(isc_bin_op_fn) &isl_map_intersect },
	{ '*',	isl_obj_map,	isl_obj_set,
		isl_obj_map,
		(isc_bin_op_fn) &isl_map_intersect_domain },
	{ '.',	isl_obj_map,	isl_obj_map,
		isl_obj_map,
		(isc_bin_op_fn) &isl_map_apply_range },
	{ ISL_TOKEN_TO,	isl_obj_set,	isl_obj_set,	isl_obj_map,
		(isc_bin_op_fn) &isl_map_from_domain_and_range },
	{ '=', isl_obj_set,	isl_obj_set,	isl_obj_bool,
		(isc_bin_op_fn) &set_is_equal },
	{ '=', isl_obj_map,	isl_obj_map,	isl_obj_bool,
		(isc_bin_op_fn) &map_is_equal },
	{ ISL_TOKEN_LE, isl_obj_set,	isl_obj_set,	isl_obj_bool,
		(isc_bin_op_fn) &set_is_subset },
	{ ISL_TOKEN_LE, isl_obj_map,	isl_obj_map,	isl_obj_bool,
		(isc_bin_op_fn) &map_is_subset },
	{ ISL_TOKEN_LT, isl_obj_set,	isl_obj_set,	isl_obj_bool,
		(isc_bin_op_fn) &set_is_strict_subset },
	{ ISL_TOKEN_LT, isl_obj_map,	isl_obj_map,	isl_obj_bool,
		(isc_bin_op_fn) &map_is_strict_subset },
	{ ISL_TOKEN_GE, isl_obj_set,	isl_obj_set,	isl_obj_bool,
		(isc_bin_op_fn) &set_is_superset },
	{ ISL_TOKEN_GE, isl_obj_map,	isl_obj_map,	isl_obj_bool,
		(isc_bin_op_fn) &map_is_superset },
	{ ISL_TOKEN_GT, isl_obj_set,	isl_obj_set,	isl_obj_bool,
		(isc_bin_op_fn) &set_is_strict_superset },
	{ ISL_TOKEN_GT, isl_obj_map,	isl_obj_map,	isl_obj_bool,
		(isc_bin_op_fn) &map_is_strict_superset },
	{ '+',	isl_obj_pw_qpolynomial,	isl_obj_pw_qpolynomial,
		isl_obj_pw_qpolynomial,
		(isc_bin_op_fn) &isl_pw_qpolynomial_add },
	{ '-',	isl_obj_pw_qpolynomial,	isl_obj_pw_qpolynomial,
		isl_obj_pw_qpolynomial,
		(isc_bin_op_fn) &isl_pw_qpolynomial_sub },
	{ '*',	isl_obj_pw_qpolynomial,	isl_obj_pw_qpolynomial,
		isl_obj_pw_qpolynomial,
		(isc_bin_op_fn) &isl_pw_qpolynomial_mul },
	{ '*',	isl_obj_pw_qpolynomial,	isl_obj_set,
		isl_obj_pw_qpolynomial,
		(isc_bin_op_fn) &isl_pw_qpolynomial_intersect_domain },
	{ '*',	isl_obj_pw_qpolynomial_fold,	isl_obj_set,
		isl_obj_pw_qpolynomial_fold,
		(isc_bin_op_fn) &isl_pw_qpolynomial_fold_intersect_domain },
	{ '@',	isl_obj_pw_qpolynomial, isl_obj_set,
		isl_obj_pw_qpolynomial,
		(isc_bin_op_fn) &isl_pw_qpolynomial_at },
	{ '@',	isl_obj_pw_qpolynomial_fold, isl_obj_set,
		isl_obj_pw_qpolynomial,
		(isc_bin_op_fn) &isl_pw_qpolynomial_fold_at },
	{ '%',	isl_obj_set,	isl_obj_set,
		isl_obj_set,
		(isc_bin_op_fn) &isl_set_gist },
	{ '%',	isl_obj_map,	isl_obj_map,
		isl_obj_map,
		(isc_bin_op_fn) &isl_map_gist },
	{ '%',	isl_obj_pw_qpolynomial,	isl_obj_set,
		isl_obj_pw_qpolynomial,
		(isc_bin_op_fn) &isl_pw_qpolynomial_gist },
	{ '%',	isl_obj_pw_qpolynomial_fold,	isl_obj_set,
		isl_obj_pw_qpolynomial_fold,
		(isc_bin_op_fn) &isl_pw_qpolynomial_fold_gist },
	0
};
struct isc_named_bin_op named_bin_ops[] = {
	{ "cross",	{ -1,	isl_obj_set,	isl_obj_set,	isl_obj_set,
		(isc_bin_op_fn) &isl_set_product } },
	{ "cross",	{ -1,	isl_obj_map,	isl_obj_map,	isl_obj_map,
		(isc_bin_op_fn) &isl_map_product } },
	NULL
};

__isl_give isl_set *set_sample(__isl_take isl_set *set)
{
	return isl_set_from_basic_set(isl_set_sample(set));
}

__isl_give isl_map *map_sample(__isl_take isl_map *map)
{
	return isl_map_from_basic_map(isl_map_sample(map));
}

static __isl_give isl_set *set_affine_hull(__isl_take isl_set *set)
{
	return isl_set_from_basic_set(isl_set_affine_hull(set));
}

static __isl_give isl_map *map_affine_hull(__isl_take isl_map *map)
{
	return isl_map_from_basic_map(isl_map_affine_hull(map));
}

static __isl_give isl_pw_qpolynomial_fold *pw_qpolynomial_upper_bound(
	__isl_take isl_pw_qpolynomial *pwqp)
{
#ifdef HAVE_GINAC
	return isl_pw_qpolynomial_bound(pwqp, isl_fold_max, BV_BOUND_BERNSTEIN);
#else
	return isl_pw_qpolynomial_bound(pwqp, isl_fold_max, BV_BOUND_RANGE);
#endif
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
	{"aff",	{ -1,	isl_obj_map,	isl_obj_map,
		(isc_un_op_fn) &map_affine_hull } },
	{"aff",	{ -1,	isl_obj_set,	isl_obj_set,
		(isc_un_op_fn) &set_affine_hull } },
	{"card",	{ -1,	isl_obj_set,	isl_obj_pw_qpolynomial,
		(isc_un_op_fn) &isl_set_card } },
	{"card",	{ -1,	isl_obj_map,	isl_obj_pw_qpolynomial,
		(isc_un_op_fn) &isl_map_card } },
	{"coalesce",	{ -1,	isl_obj_set,	isl_obj_set,
		(isc_un_op_fn) &isl_set_coalesce } },
	{"coalesce",	{ -1,	isl_obj_map,	isl_obj_map,
		(isc_un_op_fn) &isl_map_coalesce } },
	{"coalesce",	{ -1,	isl_obj_pw_qpolynomial,	isl_obj_pw_qpolynomial,
		(isc_un_op_fn) &isl_pw_qpolynomial_coalesce } },
	{"coalesce",	{ -1,	isl_obj_pw_qpolynomial_fold,
		isl_obj_pw_qpolynomial_fold,
		(isc_un_op_fn) &isl_pw_qpolynomial_fold_coalesce } },
	{"deltas",	{ -1,	isl_obj_map,	isl_obj_set,
		(isc_un_op_fn) &isl_map_deltas } },
	{"dom",	{ -1,	isl_obj_map,	isl_obj_set,
		(isc_un_op_fn) &isl_map_domain } },
	{"dom",	{ -1,	isl_obj_pw_qpolynomial,	isl_obj_set,
		(isc_un_op_fn) &isl_pw_qpolynomial_domain } },
	{"dom",	{ -1,	isl_obj_pw_qpolynomial_fold,	isl_obj_set,
		(isc_un_op_fn) &isl_pw_qpolynomial_fold_domain } },
	{"ran",	{ -1,	isl_obj_map,	isl_obj_set,
		(isc_un_op_fn) &isl_map_range } },
	{"lexmin",	{ -1,	isl_obj_map,	isl_obj_map,
		(isc_un_op_fn) &isl_map_lexmin } },
	{"lexmax",	{ -1,	isl_obj_map,	isl_obj_map,
		(isc_un_op_fn) &isl_map_lexmax } },
	{"lexmin",	{ -1,	isl_obj_set,	isl_obj_set,
		(isc_un_op_fn) &isl_set_lexmin } },
	{"lexmax",	{ -1,	isl_obj_set,	isl_obj_set,
		(isc_un_op_fn) &isl_set_lexmax } },
	{"sample",	{ -1,	isl_obj_set,	isl_obj_set,
		(isc_un_op_fn) &set_sample } },
	{"sample",	{ -1,	isl_obj_map,	isl_obj_map,
		(isc_un_op_fn) &map_sample } },
	{"sum",		{ -1,	isl_obj_pw_qpolynomial,	isl_obj_pw_qpolynomial,
		(isc_un_op_fn) &isl_pw_qpolynomial_sum } },
	{"ub",		{ -1,	isl_obj_pw_qpolynomial, isl_obj_pw_qpolynomial_fold,
		(isc_un_op_fn) &pw_qpolynomial_upper_bound } },
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
	}

	free(name);
	obj.v = obj.type->copy(obj.v);
	return obj;
}

static struct isc_bin_op *read_bin_op_if_available(struct isl_stream *s,
	isl_obj_type lhs)
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
		if (bin_ops[i].lhs != lhs)
			continue;

		isl_token_free(tok);
		return &bin_ops[i];
	}

	for (i = 0; ; ++i) {
		if (!named_bin_ops[i].name)
			break;
		if (named_bin_ops[i].op.op != tok->type)
			continue;
		if (named_bin_ops[i].op.lhs != lhs)
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
	isl_obj_type arg)
{
	int i;

	for (i = 0; ; ++i) {
		if (!named_un_ops[i].name)
			break;
		if (named_un_ops[i].op.op != like->op)
			continue;
		if (named_un_ops[i].op.arg != arg)
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

	op = find_matching_un_op(op, obj.type);

	isl_assert(s->ctx, op, goto error);
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

	isl_assert(ctx, obj.type == isl_obj_map, goto error);
	list = isl_list_alloc(ctx, 2);
	if (!list)
		goto error;

	list->obj[0].type = isl_obj_map;
	list->obj[0].v = isl_map_transitive_closure(obj.v, &exact);
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
	isl_assert(s, i < list->n, goto error);
	if (isl_stream_eat(s, ']'))
		goto error;

	obj = list->obj[i];
	obj.v = obj.type->copy(obj.v);

	isl_list_free(list);

	return obj;
error:
	free_obj(obj);
	obj.type = isl_obj_none;
	obj.v = NULL;
	return obj;
}

static struct isl_obj apply(struct isl_stream *s, __isl_take isl_map *map,
	struct isl_hash_table *table)
{
	struct isl_obj obj;

	obj = read_expr(s, table);
	isl_assert(s->ctx, obj.type == isl_obj_set || obj.type == isl_obj_map,
		goto error);

	if (obj.type == isl_obj_set) {
		obj.v = isl_set_apply(obj.v, map);
	} else
		obj.v = isl_map_apply_range(obj.v, map);
	if (!obj.v)
		goto error;

	if (isl_stream_eat(s, ')'))
		goto error;

	return obj;
error:
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
	isl_assert(s->ctx, obj.type == isl_obj_map, goto error);

	obj.v = isl_map_reverse(obj.v);
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
	FILE *file;

	tok = isl_stream_next_token(s);
	if (!tok || tok->type != ISL_TOKEN_STRING) {
		isl_stream_error(s, tok, "expecting filename");
		isl_token_free(tok);
		goto error;
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

static struct isl_obj read_obj(struct isl_stream *s,
	struct isl_hash_table *table)
{
	struct isl_obj obj = { isl_obj_none, NULL };
	char *name = NULL;
	struct isc_un_op *op = NULL;

	if (isl_stream_eat_if_available(s, '(')) {
		obj = read_expr(s, table);
		if (isl_stream_eat(s, ')'))
			goto error;
	} else {
		op = read_prefix_un_op_if_available(s);
		if (op)
			return read_un_op_expr(s, table, op);

		if (isl_stream_eat_if_available(s, read_op))
			return read_from_file(s);

		name = isl_stream_read_ident_if_available(s);
		if (name) {
			obj = stored_obj(s->ctx, table, name);
		} else {
			obj = isl_stream_read_obj(s);
			assert(obj.v);
		}
	}

	if (isl_stream_eat_if_available(s, '^'))
		obj = power(s, obj);
	else if (obj.type == isl_obj_list && isl_stream_eat_if_available(s, '['))
		obj = obj_at_index(s, obj);
	else if (obj.type == isl_obj_map && isl_stream_eat_if_available(s, '('))
		obj = apply(s, obj.v, table);

	return obj;
error:
	free_obj(obj);
	obj.type = isl_obj_none;
	obj.v = NULL;
	return obj;
}

static struct isc_bin_op *find_matching_bin_op(struct isc_bin_op *like,
	isl_obj_type lhs, isl_obj_type rhs)
{
	int i;

	for (i = 0; ; ++i) {
		if (!bin_ops[i].op)
			break;
		if (bin_ops[i].op != like->op)
			continue;
		if (bin_ops[i].lhs != lhs)
			continue;
		if (bin_ops[i].rhs != rhs)
			continue;

		return &bin_ops[i];
	}

	for (i = 0; ; ++i) {
		if (!named_bin_ops[i].name)
			break;
		if (named_bin_ops[i].op.op != like->op)
			continue;
		if (named_bin_ops[i].op.lhs != lhs)
			continue;
		if (named_bin_ops[i].op.rhs != rhs)
			continue;

		return &named_bin_ops[i].op;
	}

	return NULL;
}

static struct isl_obj read_expr(struct isl_stream *s,
	struct isl_hash_table *table)
{
	struct isl_obj obj = { isl_obj_none, NULL };
	struct isl_obj right_obj = { isl_obj_none, NULL };

	obj = read_obj(s, table);
	for (;;) {
		struct isc_bin_op *op = NULL;

		op = read_bin_op_if_available(s, obj.type);
		if (!op)
			break;

		right_obj = read_obj(s, table);

		op = find_matching_bin_op(op, obj.type, right_obj.type);

		isl_assert(s->ctx, op, goto error);
		obj.v = op->fn(obj.v, right_obj.v);
		obj.type = op->res;
	}

	return obj;
error:
	free_obj(right_obj);
	free_obj(obj);
	obj.type = isl_obj_none;
	obj.v = NULL;
	return obj;
}

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
			return;
	} else {
		p = obj.type->print(p, obj.v);
		p = isl_printer_end_line(p);
		free_obj(obj);
	}

	return p;
error:
	free(lhs);
	free_obj(obj);
	return p;
}

int free_cb(void *entry)
{
	struct isl_named_obj *named = entry;

	free_obj(named->obj);
	free(named->name);
	free(named);

	return 0;
}

static void register_named_ops(struct isl_stream *s)
{
	int i;

	read_op = isl_stream_register_keyword(s, "read");
	assert(read_op != ISL_TOKEN_ERROR);

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

	while (!s->eof) {
		p = read_line(s, table, p);
	}

	isl_printer_free(p);
	isl_hash_table_foreach(ctx, table, free_cb);
	isl_hash_table_free(ctx, table);
	isl_stream_free(s);
	isl_ctx_free(ctx);

	return 0;
}
