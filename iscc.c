#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <isl_obj.h>
#include <isl_stream.h>
#include <isl_obj_list.h>
#include <barvinok/barvinok.h>

#include "config.h"
#ifdef HAVE_GINAC
#include <barvinok/bernstein.h>
#endif

static int isl_bool_false = 0;
static int isl_bool_true = 1;
static int isl_bool_error = -1;

static void *isl_obj_bool_copy(void *v)
{
	return v;
}

static void isl_obj_bool_free(void *v)
{
}

static void isl_obj_bool_print(void *v, FILE *out)
{
	if (v == &isl_bool_true)
		fprintf(out, "True");
	else if (v == &isl_bool_false)
		fprintf(out, "False");
	else
		fprintf(out, "Error");
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
	0
};

__isl_give isl_set *set_sample(__isl_take isl_set *set)
{
	return isl_set_from_basic_set(isl_set_sample(set));
}

__isl_give isl_map *map_sample(__isl_take isl_map *map)
{
	return isl_map_from_basic_map(isl_map_sample(map));
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
	{"card",	{ -1,	isl_obj_set,	isl_obj_pw_qpolynomial,
		(isc_un_op_fn) &isl_set_card } },
	{"card",	{ -1,	isl_obj_map,	isl_obj_pw_qpolynomial,
		(isc_un_op_fn) &isl_map_card } },
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
#ifdef HAVE_GINAC
	{"ub",		{ -1,	isl_obj_pw_qpolynomial, isl_obj_pw_qpolynomial_fold,
		(isc_un_op_fn) &isl_pw_qpolynomial_upper_bound } },
#endif
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
		if (!named_un_ops[i].op.op)
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
		if (!named_un_ops[i].op.op)
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

		name = isl_stream_read_ident_if_available(s);
		if (name) {
			obj = stored_obj(s->ctx, table, name);
		} else {
			obj = isl_stream_read_obj(s);
			assert(obj.v);
		}
	}

	if (isl_stream_eat_if_available(s, '^')) {
		if (isl_stream_eat(s, '+'))
			goto error;
		obj = transitive_closure(s->ctx, obj);
	} else if (obj.type == isl_obj_list && isl_stream_eat_if_available(s, '['))
		obj = obj_at_index(s, obj);

	return obj;
error:
	free_obj(obj);
	obj.type = isl_obj_none;
	obj.v = NULL;
	return obj;
}

static struct isc_un_op *find_matching_bin_op(struct isc_un_op *like,
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

		return &bin_ops[i].op;
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

static void read_line(struct isl_stream *s, struct isl_hash_table *table)
{
	struct isl_obj obj = { isl_obj_none, NULL };
	char *lhs = NULL;
	int assign = 0;
	struct isc_bin_op *op = NULL;

	if (isl_stream_is_empty(s))
		return;

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
		obj.type->print(obj.v, stdout);
		printf("\n");
		free_obj(obj);
	}

	return;
error:
	free(lhs);
	free_obj(obj);
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

	for (i = 0; ; ++i) {
		if (!named_un_ops[i].name)
			break;
		named_un_ops[i].op.op = isl_stream_register_keyword(s,
							named_un_ops[i].name);
		assert(named_un_ops[i].op.op != ISL_TOKEN_ERROR);
	}
}

int main(int argc, char **argv)
{
	struct isl_ctx *ctx;
	struct isl_stream *s;
	struct isl_hash_table *table;

	ctx = isl_ctx_alloc();
	s = isl_stream_new_file(ctx, stdin);
	assert(s);
	table = isl_hash_table_alloc(ctx, 10);
	assert(table);

	register_named_ops(s);

	while (!feof(stdin)) {
		read_line(s, table);
	}

	isl_hash_table_foreach(ctx, table, free_cb);
	isl_hash_table_free(ctx, table);
	isl_stream_free(s);
	isl_ctx_free(ctx);

	return 0;
}
