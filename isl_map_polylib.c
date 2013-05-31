/*
 * Copyright 2008-2009 Katholieke Universiteit Leuven
 *
 * Use of this software is governed by the GNU GPLv2 license
 *
 * Written by Sven Verdoolaege, K.U.Leuven, Departement
 * Computerwetenschappen, Celestijnenlaan 200A, B-3001 Leuven, Belgium
 */

#include <isl/val_gmp.h>
#include <isl/set.h>
#include <isl/map.h>
#include <isl/constraint.h>
#include "isl_set_polylib.h"
#include "isl_map_polylib.h"

static __isl_give isl_constraint *copy_constraint_from(
	__isl_take isl_constraint *dst, Value *src)
{
	int i, j, k;
	isl_ctx *ctx = isl_constraint_get_ctx(dst);
	isl_val *v;
	enum isl_dim_type types[] = { isl_dim_in, isl_dim_out, isl_dim_param };

	k = 1;
	for (i = 0; i < 3; ++i) {
		int n = isl_constraint_dim(dst, types[i]);
		for (j = 0; j < n; ++j, ++k) {
			v = isl_val_int_from_gmp(ctx, src[k]);
			dst = isl_constraint_set_coefficient_val(dst, types[i],
								j, v);
		}
	}

	v = isl_val_int_from_gmp(ctx, src[k]);
	dst = isl_constraint_set_constant_val(dst, v);

	return dst;
}

static __isl_give isl_basic_map *add_equality(__isl_take isl_basic_map *bmap,
			 Value *constraint)
{
	isl_constraint *c;

	c = isl_equality_alloc(isl_basic_map_get_local_space(bmap));

	c = copy_constraint_from(c, constraint);

	bmap = isl_basic_map_add_constraint(bmap, c);

	return bmap;
}

static __isl_give isl_basic_map *add_inequality(__isl_take isl_basic_map *bmap,
			 Value *constraint)
{
	isl_constraint *c;

	c = isl_inequality_alloc(isl_basic_map_get_local_space(bmap));

	copy_constraint_from(c, constraint);

	bmap = isl_basic_map_add_constraint(bmap, c);

	return bmap;
}

static __isl_give isl_basic_map *copy_constraints(
			__isl_take isl_basic_map *bmap, Polyhedron *P)
{
	int i;

	for (i = 0; i < P->NbConstraints; ++i) {
		if (value_zero_p(P->Constraint[i][0]))
			bmap = add_equality(bmap, P->Constraint[i]);
		else
			bmap = add_inequality(bmap, P->Constraint[i]);
	}

	return bmap;
}

struct isl_basic_set *isl_basic_set_new_from_polylib(Polyhedron *P,
			struct isl_space *dim)
{
	isl_ctx *ctx;

	if (!dim)
		return NULL;
	ctx = isl_space_get_ctx(dim);
	isl_assert(ctx, isl_space_dim(dim, isl_dim_in) == 0, return NULL);

	return (struct isl_basic_set *)
		isl_basic_map_new_from_polylib(P, dim);
}

struct isl_basic_map *isl_basic_map_new_from_polylib(Polyhedron *P,
			struct isl_space *dim)
{
	isl_ctx *ctx;
	struct isl_basic_map *bmap;
	unsigned n_out;
	unsigned extra;

	if (!dim)
		return NULL;

	ctx = isl_space_get_ctx(dim);
	isl_assert(ctx, P, goto error);
	isl_assert(ctx, P->Dimension >= isl_space_dim(dim, isl_dim_all),
		    goto error);

	n_out = isl_space_dim(dim, isl_dim_out);
	extra = P->Dimension - isl_space_dim(dim, isl_dim_all);
	dim = isl_space_from_domain(isl_space_wrap(dim));
	dim = isl_space_add_dims(dim, isl_dim_out, extra);
	bmap = isl_basic_map_universe(dim);
	if (!bmap)
		return NULL;

	bmap = copy_constraints(bmap, P);
	bmap = isl_basic_set_unwrap(isl_basic_map_domain(bmap));

	return bmap;
error:
	isl_space_free(dim);
	return NULL;
}

struct isl_set *isl_set_new_from_polylib(Polyhedron *D, struct isl_space *dim)
{
	isl_ctx *ctx;
	struct isl_set *set = NULL;
	Polyhedron *P;

	if (!dim)
		return NULL;
	ctx = isl_space_get_ctx(dim);
	isl_assert(ctx, isl_space_dim(dim, isl_dim_in) == 0, return NULL);

	set = isl_set_empty(isl_space_copy(dim));
	if (!set)
		goto error;

	for (P = D; P; P = P->next)
		set = isl_set_union_disjoint(set,
		    isl_set_from_basic_set(
		    isl_basic_set_new_from_polylib(P, isl_space_copy(dim))));
	isl_space_free(dim);
	return set;
error:
	isl_space_free(dim);
	return NULL;
}

struct isl_map *isl_map_new_from_polylib(Polyhedron *D, struct isl_space *dim)
{
	struct isl_map *map = NULL;
	Polyhedron *P;

	if (!dim)
		return NULL;

	map = isl_map_empty(isl_space_copy(dim));
	if (!map)
		goto error;

	for (P = D; P; P = P->next)
		map = isl_map_union_disjoint(map,
		    isl_map_from_basic_map(
		    isl_basic_map_new_from_polylib(P, isl_space_copy(dim))));
	isl_space_free(dim);
	return map;
error:
	isl_space_free(dim);
	return NULL;
}

static int count_constraints(__isl_take isl_constraint *c, void *user)
{
	int *n = (int *)user;
	(*n)++;
	isl_constraint_free(c);
	return 0;
}

struct isl_poly_copy {
	int n;
	Matrix *M;
};

static int copy_constraint_to(__isl_take isl_constraint *c, void *user)
{
	int i, j, k;
	enum isl_dim_type types[] = { isl_dim_in, isl_dim_out,
					isl_dim_div, isl_dim_param };
	struct isl_poly_copy *data = (struct isl_poly_copy *)user;
	isl_val *v;

	if (isl_constraint_is_equality(c))
		value_set_si(data->M->p[data->n][0], 0);
	else
		value_set_si(data->M->p[data->n][0], 1);
	k = 1;
	for (i = 0; i < 4; ++i) {
		int n = isl_constraint_dim(c, types[i]);
		for (j = 0; j < n; ++j, ++k) {
			v = isl_constraint_get_coefficient_val(c, types[i], j);
			isl_val_get_num_gmp(v, data->M->p[data->n][k]);
			isl_val_free(v);
		}
	}
	v = isl_constraint_get_constant_val(c);
	isl_val_get_num_gmp(v, data->M->p[data->n][k]);
	isl_val_free(v);
	isl_constraint_free(c);
	data->n++;
	return 0;
}

Polyhedron *isl_basic_map_to_polylib(struct isl_basic_map *bmap)
{
	Polyhedron *P;
	unsigned off;
	unsigned nparam;
	unsigned n_in;
	unsigned n_out;
	unsigned max_rays;
	unsigned n_div;
	int n = 0;
	struct isl_poly_copy data;

	if (!bmap)
		return NULL;

	if (isl_basic_map_is_rational(bmap))
		max_rays = POL_NO_DUAL;
	else
		max_rays = POL_NO_DUAL | POL_INTEGER;

	if (isl_basic_map_foreach_constraint(bmap, &count_constraints, &n) < 0)
		return NULL;

	nparam = isl_basic_map_n_param(bmap);
	n_in = isl_basic_map_n_in(bmap);
	n_out = isl_basic_map_n_out(bmap);
	n_div = isl_basic_map_dim(bmap, isl_dim_div);
	data.M = Matrix_Alloc(n, 1 + n_in + n_out + n_div + nparam + 1);
	data.n = 0;
	if (isl_basic_map_foreach_constraint(bmap,
					    &copy_constraint_to, &data) < 0) {
		Matrix_Free(data.M);
		return NULL;
	}
	P = Constraints2Polyhedron(data.M, max_rays);
	Matrix_Free(data.M);

	return P;
}

static int add_basic_map(__isl_take isl_basic_map *bmap, void *user)
{
	Polyhedron ***next = user;

	**next = isl_basic_map_to_polylib(bmap);
	*next = &(**next)->next;

	isl_basic_map_free(bmap);
	return 0;
}

Polyhedron *isl_map_to_polylib(struct isl_map *map)
{
	int i;
	Polyhedron *R = NULL;
	Polyhedron **next = &R;

	if (!map)
		return NULL;

	if (isl_map_foreach_basic_map(map, &add_basic_map, &next) < 0)
		goto error;

	return R ? R : Empty_Polyhedron(isl_map_dim(map, isl_dim_all));
error:
	Domain_Free(R);
	return NULL;
}

Polyhedron *isl_basic_set_to_polylib(struct isl_basic_set *bset)
{
	return isl_basic_map_to_polylib((struct isl_basic_map *)bset);
}

Polyhedron *isl_set_to_polylib(struct isl_set *set)
{
	return isl_map_to_polylib((struct isl_map *)set);
}
