#include <assert.h>
#include <barvinok/options.h>
#include <barvinok/sample.h>
#include <barvinok/util.h>
#include "hilbert.h"
#include "hull.h"
#include "ilp.h"
#include "polysign.h"

struct integer_hull {
    Polyhedron	*P;	/* Original polyhedron or cone */
    Polyhedron  *init;  /* Initial (under)approximation of integer hull */

    Matrix 	*F;	/* Set of already computed facets */
    int		n_F;	/* The number of computed facets  */

    /* Computes a lower bound for the objective function obj over
     * the integer hull, possibly exploiting the already computed
     * facets of the integer hull given in hull->F.
     */
    void (*set_lower_bound)(struct integer_hull *hull, Value *obj,
			    Value *lower, struct barvinok_options *options);
};

static int matrix_add(Matrix *M, int n, Value *v)
{
    if (n >= M->NbRows)
	Matrix_Extend(M, 3*(M->NbRows+10)/2);
    Vector_Copy(v, M->p[n], M->NbColumns);
    return n+1;
}

static int select_best(struct integer_hull *hull,
		       Polyhedron *P, Matrix *candidates,
		       int *n_candidates, int *candidate,
		       Value *min, Value *max,
		       struct barvinok_options *options)
{
    int i, j, k;
    Matrix *M = candidates;
    Vector *lower = Vector_Alloc(P->NbConstraints);
    Vector *upper = Vector_Alloc(P->NbConstraints+2);
    Vector *distances = Vector_Alloc(P->NbConstraints);
    int i_min;
    int non_zero = 0;

    i_min = -1;
    for (i = 0; i < P->NbConstraints; ++i) {
	if (First_Non_Zero(P->Constraint[i]+1, P->Dimension) == -1)
	    continue;
	for (j = 0; j < hull->n_F; ++j)
	    if (Vector_Equal(hull->F->p[j]+1, P->Constraint[i]+1, P->Dimension))
		break;
	if (j < hull->n_F)
	    continue;
	hull->set_lower_bound(hull, P->Constraint[i]+1, &lower->p[i], options);
	value_oppose(upper->p[i], P->Constraint[i][1+P->Dimension]);
	value_subtract(distances->p[i], upper->p[i], lower->p[i]);
	if (value_zero_p(distances->p[i]))
	    continue;
	if (i_min == -1 || value_lt(distances->p[i], distances->p[i_min]))
	    i_min = i;
	non_zero++;
    }
    if (i_min == -1) {
	Vector_Free(lower);
	Vector_Free(upper);
	Vector_Free(distances);
	return -1;
    }

    *candidate = -1;
    for (j = 0, k = 0; j < *n_candidates; ++j) {
	int keep = 0;
	for (i = 0; i < P->NbConstraints; ++i) {
	    if (value_zero_p(distances->p[i]))
		continue;
	    Inner_Product(M->p[j], P->Constraint[i]+1, P->Dimension,
			  &upper->p[P->NbConstraints]);
	    value_addto(upper->p[P->NbConstraints+1],
			upper->p[P->NbConstraints],
			P->Constraint[i][1+P->Dimension]);
	    if (value_neg_p(upper->p[P->NbConstraints+1]))
		keep = 1;
	    if (value_lt(upper->p[P->NbConstraints], upper->p[i])) {
		value_assign(upper->p[i], upper->p[P->NbConstraints]);
		value_subtract(distances->p[i], upper->p[i], lower->p[i]);
		if (i_min == i)
		    *candidate = k;
		else if (value_lt(distances->p[i], distances->p[i_min])) {
		    i_min = i;
		    *candidate = k;
		}
	    }
	}
	if (keep) {
	    if (k != j)
		Vector_Exchange(M->p[j], M->p[k], M->NbColumns);
	    ++k;
	}
    }
    *n_candidates = k;

    value_decrement(*max, upper->p[i_min]);
    value_assign(*min, lower->p[i_min]);

    Vector_Free(lower);
    Vector_Free(upper);
    Vector_Free(distances);

    return i_min;
}

/* Return the (integer) vertices of P */
static Matrix *Polyhedron_Vertices(Polyhedron *P)
{
    int i, j;
    Matrix *M = Matrix_Alloc(P->NbRays, P->Dimension+1);

    for (i = 0, j = 0; i < P->NbRays; ++i) {
	if (value_zero_p(P->Ray[i][1+P->Dimension]))
	    continue;
	assert(value_one_p(P->Ray[i][1+P->Dimension]));
	value_set_si(M->p[j][P->Dimension], 1);
	Vector_Copy(P->Ray[i]+1, M->p[j++], P->Dimension);
    }
    M->NbRows = j;

    return M;
}

/* Extends an initial approximation hull->init of the integer hull
 * of hull->P using generalized basis reduction.
 * In particular, it considers each facet of the current
 * approximation add optimizes in the direction of the normal,
 * adding the optimal point to the approximation as well as
 * all other points that may have been found along the way,
 * until no facet yields any more new points.
 * Returns a list of the vertices of this integer hull.
 */
static Matrix *gbr_hull_extend(struct integer_hull *hull,
				struct barvinok_options *options)
{
    Value min, max;
    Polyhedron *Q = hull->init;
    Vector *ray = Vector_Alloc(2+Q->Dimension);
    Matrix *candidates;
    int n_candidates = 0;
    Matrix *vertices;

    hull->F= Matrix_Alloc(Q->NbConstraints, 2+Q->Dimension);
    hull->n_F = 0;

    candidates = Matrix_Alloc(0, Q->Dimension);

    value_init(min);
    value_init(max);

    value_set_si(ray->p[0], 1);
    value_set_si(ray->p[1+Q->Dimension], 1);

    for (;;) {
	Polyhedron *R;
	int i, i_min, candidate;
	Vector *opt;
	Value *vertex = NULL;

	i_min = select_best(hull, Q, candidates, &n_candidates, &candidate,
				&min, &max, options);
	if (i_min == -1)
	    break;

	opt = Polyhedron_Integer_Minimum(hull->P, Q->Constraint[i_min]+1,
					 min, max, candidates, &n_candidates,
					 options);
	candidates->NbRows = n_candidates;

	hull->n_F = matrix_add(hull->F, hull->n_F, Q->Constraint[i_min]);

	if (opt)
	    vertex = opt->p;
	else if (candidate != -1)
	    vertex = candidates->p[candidate];

	if (!vertex)
	    continue;

	Inner_Product(hull->F->p[hull->n_F-1]+1, vertex, Q->Dimension,
		      &hull->F->p[hull->n_F-1][1+Q->Dimension]);
	value_oppose(hull->F->p[hull->n_F-1][1+Q->Dimension],
		     hull->F->p[hull->n_F-1][1+Q->Dimension]);

	Vector_Copy(vertex, ray->p+1, Q->Dimension);
	R = AddRays(ray->p, 1, Q, options->MaxRays);
	Polyhedron_Free(Q);
	Q = R;

	if (opt)
	    Vector_Free(opt);
    }

    vertices = Polyhedron_Vertices(Q);

    Polyhedron_Free(Q);
    hull->init = NULL;

    value_clear(min);
    value_clear(max);

    Vector_Free(ray);
    Matrix_Free(hull->F);
    hull->F = NULL;
    Matrix_Free(candidates);

    return vertices;
}

/* Returns the Minkowski sum of the cone and the polytope
 * that is the convex hull of its (integer) extremal rays.
 */
static Polyhedron *truncate_cone(Polyhedron *C,
				 struct barvinok_options *options)
{
    int i;
    Matrix *V;
    Polyhedron *P;

    V = Matrix_Alloc(2*C->NbRays, 2+C->Dimension);
    for (i = 0; i < C->NbRays; ++i) {
	if (value_notzero_p(C->Ray[i][1+C->Dimension]))
	    continue;
	Vector_Copy(C->Ray[i], V->p[C->NbRays+i], 1+C->Dimension+1);
	Vector_Copy(C->Ray[i], V->p[i], 1+C->Dimension);
	value_set_si(V->p[i][1+C->Dimension], 1);
    }
    P = Rays2Polyhedron(V, options->MaxRays);
    Matrix_Free(V);
    return P;
}

/* Frees original list and transformation matrix CV */
static Matrix *transform_list_of_vertices(Matrix *list, Matrix *CV)
{
    Matrix *T, *M;
    T = Transpose(CV);
    M = Matrix_Alloc(list->NbRows, T->NbColumns);
    Matrix_Product(list, T, M);
    Matrix_Free(list);
    Matrix_Free(CV);
    Matrix_Free(T);

    return M;
}

/* Set the lower bound for optimizing along the normal of a facet
 * in case of computing the integer hull of a truncated cone
 * (i.e., a cone with the origin removed).
 * In particular, if the constraint is one bounding the original
 * cone, then the lower bound is zero; if it is a new constraint,
 * then the lower bound is one.
 * A more accurate bound can be obtained by looking at the
 * already computed facets of the integer hull.
 */
static void set_to_one(struct integer_hull *hull, Value *obj,
		       Value *lower, struct barvinok_options *options)
{
    if (value_zero_p(obj[hull->P->Dimension]))
	value_set_si(*lower, 0);
    else
	value_set_si(*lower, 1);
}

static Matrix *gbr_hull(Polyhedron *C, struct barvinok_options *options)
{
    Matrix *vertices;
    Matrix *CV = NULL;
    struct integer_hull hull;

    POL_ENSURE_VERTICES(C);
    remove_all_equalities(&C, NULL, NULL, &CV, 0, options->MaxRays);

    POL_ENSURE_VERTICES(C);

    hull.P = C;
    hull.init = truncate_cone(C, options);
    hull.set_lower_bound = set_to_one;
    vertices = gbr_hull_extend(&hull, options);

    if (CV) {
	vertices = transform_list_of_vertices(vertices, CV);
	Polyhedron_Free(C);
    }

    return vertices;
}

Matrix *Cone_Integer_Hull(Polyhedron *C, struct barvinok_options *options)
{
    switch(options->integer_hull) {
    case BV_HULL_GBR:
	return gbr_hull(C, options);
    case BV_HULL_HILBERT:
	return Cone_Hilbert_Integer_Hull(C, options);
    default:
	assert(0);
    }
}

/* Computes initial full-dimensional approximation of the integer hull of P,
 * or discovers an implicit equality of P.
 * We start off with the integer vertices of P itself, if any.
 * If there are no such vertices (available), then we simply
 * take an arbitrary integer point in the polytope.
 * Then we keep optimizing over normals of equalities in the description
 * of the current approximation, until we find an equality that holds
 * for the entire integer hull of P or until we have obtained a
 * full-dimensional approximation.
 */
static Polyhedron *internal_polytope(Polyhedron *P,
				     struct barvinok_options *options)
{
    int i, j;
    Polyhedron *init;
    Matrix *vertices = Matrix_Alloc(P->NbRays, P->Dimension+2);

    for (i = 0, j = 0; i < P->NbRays; ++i) {
	if (value_notone_p(P->Ray[i][1+P->Dimension]))
	    continue;
	Vector_Copy(P->Ray[i], vertices->p[j++], 2+P->Dimension);
    }
    vertices->NbRows = j;
    if (!j) {
	Vector *sample = Polyhedron_Sample(P, options);
	if (sample) {
	    value_set_si(vertices->p[0][1], 1);
	    Vector_Copy(sample->p, vertices->p[0]+1, P->Dimension+1);
	    Vector_Free(sample);
	    vertices->NbRows = 1;
	}
    }
    init = Rays2Polyhedron(vertices, options->MaxRays);
    Matrix_Free(vertices);

    while (!emptyQ(init) && init->NbEq) {
	Vector *sample;
	Vector *v = Vector_Alloc(init->Dimension+2);
	Polyhedron *Q;

	value_set_si(v->p[0], 1);
	Vector_Copy(init->Constraint[0]+1, v->p+1, init->Dimension+1);
	value_decrement(v->p[1+init->Dimension],v->p[1+init->Dimension]);

	Q = AddConstraints(v->p, 1, P, options->MaxRays);
	sample = Polyhedron_Sample(Q, options);
	Polyhedron_Free(Q);
	if (!sample) {
	    Vector_Oppose(init->Constraint[0]+1, v->p+1, init->Dimension+1);
	    value_decrement(v->p[1+init->Dimension],v->p[1+init->Dimension]);

	    Q = AddConstraints(v->p, 1, P, options->MaxRays);
	    sample = Polyhedron_Sample(Q, options);
	    Polyhedron_Free(Q);
	}
	if (!sample) {
	    P = AddConstraints(init->Constraint[0], 1, P, options->MaxRays);
	    Polyhedron_Free(init);
	    Vector_Free(v);
	    return P;
	}
	assert(sample);

	Vector_Copy(sample->p, v->p+1, init->Dimension+1);
	Q = AddRays(v->p, 1, init, options->MaxRays);
	Polyhedron_Free(init);
	init = Q;

	Vector_Free(sample);
	Vector_Free(v);
    }

    return init;
}

/* Set the lower bound for optimizing along the normal of a facet
 * in case of computing the integer hull of a polyhedron.
 * Although obj contains the full affine objective function,
 * the calling function expect a lower bound on the linear part only.
 * We therefore need to subtract the constant from the computed
 * optimum of the linear relaxation.
 */
static void set_lower(struct integer_hull *hull, Value *obj,
		      Value *lower, struct barvinok_options *options)
{
    enum lp_result res;
    Value one;

    value_init(one);
    value_set_si(one, 1);

    res = polyhedron_opt(hull->P, obj, one, lp_min, lower, options);
    value_subtract(*lower, *lower, obj[hull->P->Dimension]);
    assert(res == lp_ok);

    value_clear(one);
}

/* Computes the vertices of the integer hull of P
 */
Matrix *Polyhedron_Integer_Hull(Polyhedron *P,
				struct barvinok_options *options)
{
    Matrix *vertices;
    Matrix *CV = NULL;
    Polyhedron *init;

    POL_ENSURE_VERTICES(P);
    remove_all_equalities(&P, NULL, NULL, &CV, 0, options->MaxRays);

    POL_ENSURE_VERTICES(P);

    if (P->Dimension == 0) {
	vertices = Matrix_Alloc(1, 1);
	value_set_si(vertices->p[0][0], 1);
    } else {
	init = internal_polytope(P, options);
	if (emptyQ(init)) {
	    vertices = Matrix_Alloc(0, P->Dimension+1);
	    Polyhedron_Free(init);
	} else if (init->NbEq) {
	    vertices = Polyhedron_Integer_Hull(init, options);
	    Polyhedron_Free(init);
	} else {
	    struct integer_hull hull;
	    hull.P = P;
	    hull.init = init;
	    hull.set_lower_bound = set_lower;
	    vertices = gbr_hull_extend(&hull, options);
	}
    }

    if (CV) {
	vertices = transform_list_of_vertices(vertices, CV);
	Polyhedron_Free(P);
    }

    return vertices;
}
