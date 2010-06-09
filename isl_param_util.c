#include <isl_set_polylib.h>
#include <isl_vertices.h>
#include "isl_param_util.h"

static Matrix *expr2vertex(Polyhedron *E, unsigned nvar)
{
	int i;
	Matrix *M;
	unsigned nparam = E->Dimension - nvar;
	Value mone;

	value_init(mone);
	value_set_si(mone, -1);
	M = Matrix_Alloc(nvar, nparam + 1 + 1);
	for (i = 0; i < nvar; ++i) {
		Vector_Scale(E->Constraint[i] + 1 + nvar, M->p[i],
				mone, nparam + 1);
		value_assign(M->p[i][nparam + 1], E->Constraint[i][1 + i]);
	}
	value_clear(mone);

	return M;
}

#define INT_BITS (sizeof(unsigned) * 8)

static int add_vertex(__isl_take isl_vertex *vertex, void *user)
{
	Param_Vertices ***next_V = (Param_Vertices ***) user;
	Param_Vertices *V;
	Polyhedron *D, *E;
	isl_basic_set *dom;
	isl_basic_set *expr;
	isl_ctx *ctx;
	unsigned nvar;

	ctx = isl_vertex_get_ctx(vertex);

	dom = isl_vertex_get_domain(vertex);
	D = isl_basic_set_to_polylib(dom);
	isl_basic_set_free(dom);

	expr = isl_vertex_get_expr(vertex);
	nvar = isl_basic_set_dim(expr, isl_dim_set);
	E = isl_basic_set_to_polylib(expr);
	isl_basic_set_free(expr);

	V = isl_alloc_type(ctx, Param_Vertices);
	V->Vertex = expr2vertex(E, nvar);
	V->Domain = Polyhedron2Constraints(D);
	V->Facets = NULL;
	V->next = NULL;

	Polyhedron_Free(D);
	Polyhedron_Free(E);

	**next_V = V;
	*next_V = &V->next;

	isl_vertex_free(vertex);

	return 0;
}

struct bv_add_chamber_data {
	Param_Domain **next_D;
	int vertex_len;
	Param_Domain *dom;
};

static int add_chamber_vertex(__isl_take isl_vertex *vertex, void *user)
{
	int j;
	struct bv_add_chamber_data *data = (struct bv_add_chamber_data *)user;
	unsigned v;

	v = isl_vertex_get_id(vertex);
	data->dom->F[v / INT_BITS] |= 1u << (INT_BITS - (v % INT_BITS) - 1);

	isl_vertex_free(vertex);

	return 0;
}

static int add_chamber(__isl_take isl_cell *cell, void *user)
{
	struct bv_add_chamber_data *data = (struct bv_add_chamber_data *)user;
	isl_ctx *ctx;
	isl_basic_set *domain;

	ctx = isl_cell_get_ctx(cell);

	domain = isl_cell_get_domain(cell);

	data->dom = isl_alloc_type(ctx, Param_Domain);
	data->dom->Domain = isl_basic_set_to_polylib(domain);
	data->dom->F = isl_calloc_array(ctx, unsigned, data->vertex_len);
	data->dom->next = NULL;

	isl_basic_set_free(domain);

	*data->next_D = data->dom;
	data->next_D = &data->dom->next;

	isl_cell_foreach_vertex(cell, &add_chamber_vertex, data);

	isl_cell_free(cell);

	return 0;
}

Param_Polyhedron *ISL_P2PP(Polyhedron *P, Polyhedron *C,
			  struct barvinok_options *options)
{
	int i, j;
	isl_ctx *ctx = isl_ctx_alloc();
	isl_dim *dim;
	isl_basic_set *bset, *context;
	isl_vertices *vertices;
	unsigned nparam = C->Dimension;
	unsigned nvar = P->Dimension - nparam;
	Param_Polyhedron *PP = isl_calloc_type(ctx, Param_Polyhedron);
	Param_Vertices **next_V;
	struct bv_add_chamber_data data;

	dim = isl_dim_set_alloc(ctx, nparam, nvar);
	bset = isl_basic_set_new_from_polylib(P, dim);
	dim = isl_dim_set_alloc(ctx, nparam, 0);
	context = isl_basic_set_new_from_polylib(C, dim);

	bset = isl_basic_set_intersect(bset, context);

	vertices = isl_basic_set_compute_vertices(bset);
	isl_basic_set_free(bset);

	PP->Rays = NULL;
	PP->nbV = isl_vertices_get_n_vertices(vertices);
	PP->Constraints = Polyhedron2Constraints(P);

	next_V = &PP->V;
	isl_vertices_foreach_vertex(vertices, &add_vertex, &next_V);

	data.next_D = &PP->D;
	data.vertex_len = (PP->nbV + INT_BITS - 1)/INT_BITS;
	isl_vertices_foreach_cell(vertices, &add_chamber, &data);

	isl_vertices_free(vertices);

	isl_ctx_free(ctx);

	return PP;
}
