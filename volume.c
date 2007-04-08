#include <barvinok/polylib.h>
#include <barvinok/barvinok.h>
#include <barvinok/options.h>
#include <barvinok/util.h>
#include "reduce_domain.h"
#include "scale.h"
#include "volume.h"

#define ALLOC(type) (type*)malloc(sizeof(type))
#define ALLOCN(type,n) (type*)malloc((n) * sizeof(type))

/* returns an evalue that corresponds to
 *
 * c/den x_param
 */
static evalue *term(int param, Value c, Value den)
{
    evalue *EP = ALLOC(evalue);
    value_init(EP->d);
    value_set_si(EP->d,0);
    EP->x.p = new_enode(polynomial, 2, param + 1);
    evalue_set_si(&EP->x.p->arr[0], 0, 1);
    value_init(EP->x.p->arr[1].x.n);
    value_assign(EP->x.p->arr[1].d, den);
    value_assign(EP->x.p->arr[1].x.n, c);
    return EP;
}

/* Computes an evalue representation of a coordinate
 * in a ParamVertices.
 */
static evalue *vertex2evalue(Value *vertex, int nparam)
{
    int i;
    evalue *E = ALLOC(evalue);
    value_init(E->d);
    evalue_set(E, vertex[nparam], vertex[nparam+1]);
    for (i = 0; i < nparam; ++i) {
	evalue *t = term(i, vertex[i], vertex[nparam+1]);
	eadd(t, E);
	free_evalue_refs(t);
	free(t);
    }
    return E;
}

static void matrix_print(evalue ***matrix, int dim, int *cols,
			 char **param_names)
{
    int i, j;

    for (i = 0; i < dim; ++i)
	for (j = 0; j < dim; ++j) {
	    int k = cols ? cols[j] : j;
	    fprintf(stderr, "%d %d: ", i, j);
	    print_evalue(stderr, matrix[i][k], param_names);
	    fprintf(stderr, "\n");
	}
}

/* Compute determinant using Laplace's formula.
 * In particular, the determinant is expanded along the last row.
 * The cols array is a list of columns that remain in the currect submatrix.
 */
static evalue *determinant_cols(evalue ***matrix, int dim, int *cols)
{
    evalue *det, *tmp;
    evalue mone;

    if (dim == 1) {
	det = ALLOC(evalue);
	value_init(det->d);
	evalue_copy(det, matrix[0][cols[0]]);
	return det;
    }

    value_init(mone.d);
    evalue_set_si(&mone, -1, 1);
    int j;
    det = NULL;
    int *newcols = ALLOCN(int, dim-1);
    for (j = 1; j < dim; ++j)
	newcols[j-1] = cols[j];
    for (j = 0; j < dim; ++j) {
	if (j != 0)
	    newcols[j-1] = cols[j-1];
	tmp = determinant_cols(matrix, dim-1, newcols);
	emul(matrix[dim-1][cols[j]], tmp);
	if ((dim+j)%2 == 0)
	    emul(&mone, tmp);
	if (!det)
	    det = tmp;
	else {
	    eadd(tmp, det);
	    free_evalue_refs(tmp);
	    free(tmp);
	}
    }
    free(newcols);
    free_evalue_refs(&mone);

    return det;
}

static evalue *determinant(evalue ***matrix, int dim)
{
    int i;
    int *cols = ALLOCN(int, dim);
    evalue *det;

    for (i = 0; i < dim; ++i)
	cols[i] = i;

    det = determinant_cols(matrix, dim, cols);

    free(cols);

    return det;
}

/* Compute the facet of P that saturates constraint c.
 */
static Polyhedron *facet(Polyhedron *P, int c, unsigned MaxRays)
{
    Polyhedron *F;
    Vector *row = Vector_Alloc(1+P->Dimension+1);
    Vector_Copy(P->Constraint[c]+1, row->p+1, P->Dimension+1);
    F = AddConstraints(row->p, 1, P, MaxRays);
    Vector_Free(row);
    return F;
}

/* Compute a dummy Param_Domain that contains all vertices of Param_Domain D
 * (which contains the vertices of P) that lie on the facet obtain by
 * saturating constraint c of P
 */
static Param_Domain *face_vertices(Param_Polyhedron *PP, Param_Domain *D,
				   Polyhedron *P, int c)
{
    int nv;
    Param_Vertices *V;
    Param_Domain *FD = ALLOC(Param_Domain);
    FD->Domain = 0;
    FD->next = 0;

    nv = (PP->nbV - 1)/(8*sizeof(int)) + 1;
    FD->F = ALLOCN(unsigned, nv);
    memset(FD->F, 0, nv * sizeof(unsigned));

    FORALL_PVertex_in_ParamPolyhedron(V, D, PP) /* _i, _ix, _bx internal counters */
	int n;
	unsigned char *supporting = supporting_constraints(P, V, &n);
	if (supporting[c])
	    FD->F[_ix] |= _bx;
	free(supporting);
    END_FORALL_PVertex_in_ParamPolyhedron;

    return FD;
}

/* Substitute parameters by the corresponding element in subs
 */
static evalue *evalue_substitute(evalue *e, evalue **subs)
{
    evalue *res = NULL;
    evalue *c;
    int i;

    if (value_notzero_p(e->d)) {
	res = ALLOC(evalue);
	value_init(res->d);
	evalue_copy(res, e);
	return res;
    }
    assert(e->x.p->type == polynomial);

    res = evalue_zero();
    for (i = e->x.p->size-1; i > 0; --i) {
	c = evalue_substitute(&e->x.p->arr[i], subs);
	eadd(c, res);
	free_evalue_refs(c);
	free(c);
	emul(subs[e->x.p->pos-1], res);
    }
    c = evalue_substitute(&e->x.p->arr[0], subs);
    eadd(c, res);
    free_evalue_refs(c);
    free(c);

    return res;
}

/* Plug in the parametric vertex V in the constraint constraint.
 * The result is stored in row, with the denominator in position 0.
 */
static void Param_Inner_Product(Value *constraint, Matrix *Vertex,
				Value *row)
{
    unsigned nparam = Vertex->NbColumns - 2;
    unsigned nvar = Vertex->NbRows;
    int j;
    Value tmp, tmp2;

    value_set_si(row[0], 1);
    Vector_Set(row+1, 0, nparam+1);

    value_init(tmp);
    value_init(tmp2);

    for (j = 0 ; j < nvar; ++j) {
	value_set_si(tmp, 1);
	value_assign(tmp2,  constraint[1+j]);
	if (value_ne(row[0], Vertex->p[j][nparam+1])) {
	    value_assign(tmp, row[0]);
	    value_lcm(row[0], Vertex->p[j][nparam+1], &row[0]);
	    value_division(tmp, row[0], tmp);
	    value_multiply(tmp2, tmp2, row[0]);
	    value_division(tmp2, tmp2, Vertex->p[j][nparam+1]);
	}
	Vector_Combine(row+1, Vertex->p[j], row+1, tmp, tmp2, nparam+1);
    }
    value_set_si(tmp, 1);
    Vector_Combine(row+1, constraint+1+nvar, row+1, tmp, row[0], nparam+1);

    value_clear(tmp);
    value_clear(tmp2);
}

struct parameter_point {
    Vector *coord;
    evalue **e;
};

struct parameter_point *parameter_point_new(unsigned nparam)
{
    struct parameter_point *point = ALLOC(struct parameter_point);
    point->coord = Vector_Alloc(nparam+1);
    point->e = NULL;
    return point;
}

evalue **parameter_point_evalue(struct parameter_point *point)
{
    int j;
    unsigned nparam = point->coord->Size-1;

    if (point->e)
	return point->e;

    point->e = ALLOCN(evalue *, nparam);
    for (j = 0; j < nparam; ++j) {
	point->e[j] = ALLOC(evalue);
	value_init(point->e[j]->d);
	evalue_set(point->e[j], point->coord->p[j], point->coord->p[nparam]);
    }

    return point->e;
}

void parameter_point_free(struct parameter_point *point)
{
    int i;
    unsigned nparam = point->coord->Size-1;

    Vector_Free(point->coord);

    if (point->e) {
	for (i = 0; i < nparam; ++i) {
	    free_evalue_refs(point->e[i]);
	    free(point->e[i]);
	}
	free(point->e);
    }
    free(point);
}

/* Computes point in pameter space where polyhedron is non-empty.
 * For each of the parametric vertices, and each of the facets
 * not (always) containing the vertex, we remove the parameter
 * values for which the facet does contain the vertex.
 */
static struct parameter_point *non_empty_point(Param_Polyhedron *PP,
	Param_Domain *D, Polyhedron *P, Polyhedron *C, unsigned MaxRays)
{
    Param_Vertices *V;
    unsigned dim = P->Dimension;
    unsigned nparam = C->Dimension;
    unsigned nvar = dim - nparam;
    Polyhedron *RD, *cut, *tmp;
    Matrix *M;
    struct parameter_point *point;
    int i, j;
    unsigned cut_MaxRays = MaxRays;
    int nv;

    nv = (PP->nbV - 1)/(8*sizeof(int)) + 1;

    POL_UNSET(cut_MaxRays, POL_INTEGER);

    M = Matrix_Alloc(1, nparam+2);
    RD = C;
    FORALL_PVertex_in_ParamPolyhedron(V, D, PP) /* _ix, _bx internal counters */
	for (i = P->NbEq; i < P->NbConstraints; ++i) {
	    if (First_Non_Zero(P->Constraint[i]+1, nvar) == -1)
		continue;
	    Param_Inner_Product(P->Constraint[i], V->Vertex, M->p[0]);
	    if (First_Non_Zero(M->p[0]+1, nparam) == -1)
		/* supporting facet,
		 * or non-supporting facet independent of params
		 */
		continue;
	    value_set_si(M->p[0][0], 0);
	    cut = Constraints2Polyhedron(M, cut_MaxRays);
	    tmp = DomainDifference(RD, cut, MaxRays);
	    if (RD != C)
		Domain_Free(RD);
	    RD = tmp;
	    Polyhedron_Free(cut);
	}
	if (emptyQ2(RD))
	    break;
    END_FORALL_PVertex_in_ParamPolyhedron;
    Matrix_Free(M);

    POL_ENSURE_VERTICES(RD);
    if (emptyQ(RD))
	point = NULL;
    else {
	point = parameter_point_new(nparam);
	for (i = 0; i < RD->NbRays; ++i)
	    if (value_notzero_p(RD->Ray[i][1+nparam]))
		break;
	assert(i < RD->NbRays);
	Vector_Copy(RD->Ray[i]+1, point->coord->p, nparam+1);
    }

    if (RD != C)
	Domain_Free(RD);

    return point;
}

static Matrix *barycenter(Param_Polyhedron *PP, Param_Domain *D)
{
    int nbV;
    Matrix *center = NULL;
    Value denom;
    Value fc, fv;
    unsigned nparam;
    int i;
    Param_Vertices *V;

    value_init(fc);
    value_init(fv);
    nbV = 0;
    FORALL_PVertex_in_ParamPolyhedron(V, D, PP)
	++nbV;
	if (!center) {
	    center = Matrix_Copy(V->Vertex);
	    nparam = center->NbColumns - 2;
	} else {
	    for (i = 0; i < center->NbRows; ++i) {
		value_assign(fc, center->p[i][nparam+1]);
		value_lcm(fc, V->Vertex->p[i][nparam+1],
			    &center->p[i][nparam+1]);
		value_division(fc, center->p[i][nparam+1], fc);
		value_division(fv, center->p[i][nparam+1],
				V->Vertex->p[i][nparam+1]);
		Vector_Combine(center->p[i], V->Vertex->p[i], center->p[i],
			       fc, fv, nparam+1);
	    }
	}
    END_FORALL_PVertex_in_ParamPolyhedron;
    value_clear(fc);
    value_clear(fv);

    value_init(denom);
    value_set_si(denom, nbV);
    for (i = 0; i < center->NbRows; ++i) {
	value_multiply(center->p[i][nparam+1], center->p[i][nparam+1], denom);
	Vector_Normalize(center->p[i], nparam+2);
    }
    value_clear(denom);

    return center;
}

/* Compute dim! times the volume of polyhedron F in Param_Domain D.
 * If F is a simplex, then the volume is computed of a recursive pyramid
 * over F with the points already in matrix.
 * Otherwise, the barycenter of F is added to matrix and the function
 * is called recursively on the facets of F.
 *
 * The first row of matrix contain the _negative_ of the first point.
 * The remaining rows of matrix contain the distance of the corresponding
 * point to the first point.
 */
static evalue *volume_in_domain(Param_Polyhedron *PP, Param_Domain *D,
				unsigned dim, evalue ***matrix,
				struct parameter_point *point, Polyhedron *C,
				int row, Polyhedron *F,
				struct barvinok_options *options);

static evalue *volume_triangulate(Param_Polyhedron *PP, Param_Domain *D,
				  unsigned dim, evalue ***matrix,
				  struct parameter_point *point, Polyhedron *C,
				  int row, Polyhedron *F,
				  struct barvinok_options *options)
{
    int j;
    evalue *tmp;
    evalue *vol;
    evalue mone;
    Matrix *center;
    unsigned cut_MaxRays = options->MaxRays;
    unsigned nparam = C->Dimension;
    Matrix *M = NULL;

    POL_UNSET(cut_MaxRays, POL_INTEGER);

    value_init(mone.d);
    evalue_set_si(&mone, -1, 1);

    center = barycenter(PP, D);
    for (j = 0; j < dim; ++j)
	matrix[row][j] = vertex2evalue(center->p[j], center->NbColumns - 2);

    if (row == 0) {
	for (j = 0; j < dim; ++j)
	    emul(&mone, matrix[row][j]);
    } else {
	for (j = 0; j < dim; ++j)
	    eadd(matrix[0][j], matrix[row][j]);
    }

    if (!point)
	M = Matrix_Alloc(1, nparam+2);

    vol = NULL;
    POL_ENSURE_FACETS(F);
    for (j = F->NbEq; j < F->NbConstraints; ++j) {
	Polyhedron *FC;
	Polyhedron *FF;
	Param_Domain *FD;
	if (First_Non_Zero(F->Constraint[j]+1, dim) == -1)
	    continue;
	if (point)
	    FC = C;
	else {
	    Polyhedron *cut;
	    int pos;
	    Param_Inner_Product(F->Constraint[j], center, M->p[0]);
	    pos = First_Non_Zero(M->p[0]+1, nparam+1);
	    if (pos == -1)
		/* barycenter always lies on facet */
		continue;
	    if (pos == nparam)
		FC = C;
	    else {
		value_set_si(M->p[0][0], 0);
		cut = Constraints2Polyhedron(M, cut_MaxRays);
		FC = DomainDifference(C, cut, options->MaxRays);
		Polyhedron_Free(cut);
		POL_ENSURE_VERTICES(FC);
		if (emptyQ(FC)) {
		    /* barycenter lies on facet for all parameters in C */
		    Polyhedron_Free(FC);
		    continue;
		}
	    }
	}
	FF = facet(F, j, options->MaxRays);
	FD = face_vertices(PP, D, F, j);
	tmp = volume_in_domain(PP, FD, dim, matrix, point, FC,
			       row+1, FF, options);
	if (FC != C)
	    Domain_Free(FC);
	if (!vol)
	    vol = tmp;
	else {
	    eadd(tmp, vol);
	    free_evalue_refs(tmp);
	    free(tmp);
	}
	Polyhedron_Free(FF);
	Param_Domain_Free(FD);
    }

    Matrix_Free(center);
    if (!point)
	Matrix_Free(M);

    for (j = 0; j < dim; ++j) {
	free_evalue_refs(matrix[row][j]);
	free(matrix[row][j]);
    }

    free_evalue_refs(&mone);
    return vol;
}

static evalue *volume_simplex(Param_Polyhedron *PP, Param_Domain *D,
				unsigned dim, evalue ***matrix,
				struct parameter_point *point,
				int row, unsigned MaxRays)
{
    evalue mone;
    Param_Vertices *V;
    evalue *vol, *val;
    int i, j;

    if (!point)
	return evalue_zero();

    value_init(mone.d);
    evalue_set_si(&mone, -1, 1);

    i = row;
    FORALL_PVertex_in_ParamPolyhedron(V, D, PP) /* _ix, _bx internal counters */
	for (j = 0; j < dim; ++j) {
	    matrix[i][j] = vertex2evalue(V->Vertex->p[j],
					   V->Vertex->NbColumns - 2);
	    if (i == 0)
		emul(&mone, matrix[i][j]);
	    else
		eadd(matrix[0][j], matrix[i][j]);
	}
	++i;
    END_FORALL_PVertex_in_ParamPolyhedron;

    vol = determinant(matrix+1, dim);

    val = evalue_substitute(vol, parameter_point_evalue(point));

    assert(value_notzero_p(val->d));
    assert(value_notzero_p(val->x.n));
    if (value_neg_p(val->x.n))
	emul(&mone, vol);

    free_evalue_refs(val);
    free(val);

    for (i = row; i < dim+1; ++i) {
	for (j = 0; j < dim; ++j) {
	    free_evalue_refs(matrix[i][j]);
	    free(matrix[i][j]);
	}
    }

    free_evalue_refs(&mone);

    return vol;
}

static evalue *volume_triangulate_lift(Param_Polyhedron *PP, Param_Domain *D,
					unsigned dim, evalue ***matrix,
					struct parameter_point *point,
					int row, unsigned MaxRays)
{
    const static int MAX_TRY=10;
    Param_Vertices *V;
    int nbV, nv;
    int i;
    int t = 0;
    Matrix *FixedRays, *M;
    Polyhedron *L;
    Param_Domain SD;
    Value tmp;
    evalue *s, *vol;

    SD.Domain = 0;
    SD.next = 0;
    nv = (PP->nbV - 1)/(8*sizeof(int)) + 1;
    SD.F = ALLOCN(unsigned, nv);

    FixedRays = Matrix_Alloc(PP->nbV+1, 1+dim+2);
    nbV = 0;
    FORALL_PVertex_in_ParamPolyhedron(V, D, PP)
	unsigned nparam = V->Vertex->NbColumns - 2;
	Param_Vertex_Common_Denominator(V);
	for (i = 0; i < V->Vertex->NbRows; ++i) {
	    value_multiply(FixedRays->p[nbV][1+i], V->Vertex->p[i][nparam],
			   point->coord->p[nparam]);
	    Inner_Product(V->Vertex->p[i], point->coord->p, nparam,
			  &FixedRays->p[nbV][1+dim]);
	    value_addto(FixedRays->p[nbV][1+i], FixedRays->p[nbV][1+i],
			FixedRays->p[nbV][1+dim]);
	}
	value_multiply(FixedRays->p[nbV][1+dim+1], V->Vertex->p[0][nparam+1],
		       point->coord->p[nparam]);
	value_set_si(FixedRays->p[nbV][0], 1);
	++nbV;
    END_FORALL_PVertex_in_ParamPolyhedron;
    value_set_si(FixedRays->p[nbV][0], 1);
    value_set_si(FixedRays->p[nbV][1+dim], 1);
    FixedRays->NbRows = nbV+1;

    value_init(tmp);
    if (0) {
try_again:
	/* Usually vol should still be NULL */
	if (vol) {
	    free_evalue_refs(vol);
	    free(vol);
	}
	Polyhedron_Free(L);
	++t;
    }
    assert(t <= MAX_TRY);
    vol = NULL;

    for (i = 0; i < nbV; ++i)
	value_set_si(FixedRays->p[i][1+dim], random_int((t+1)*dim*nbV)+1);

    M = Matrix_Copy(FixedRays);
    L = Rays2Polyhedron(M, MaxRays);
    Matrix_Free(M);

    POL_ENSURE_FACETS(L);
    for (i = 0; i < L->NbConstraints; ++i) {
	int r;
	/* Ignore perpendicular facets, i.e., facets with 0 z-coordinate */
	if (value_negz_p(L->Constraint[i][1+dim]))
	    continue;

	memset(SD.F, 0, nv * sizeof(unsigned));
	nbV = 0;
	r = 0;
	FORALL_PVertex_in_ParamPolyhedron(V, D, PP) /* _ix, _bx internal */
	    Inner_Product(FixedRays->p[nbV]+1, L->Constraint[i]+1, dim+2, &tmp);
	    if (value_zero_p(tmp)) {
		if (r > dim-row)
		    goto try_again;
		SD.F[_ix] |= _bx;
		++r;
	    }
	    ++nbV;
	END_FORALL_PVertex_in_ParamPolyhedron;
	assert(r == (dim-row)+1);

	s = volume_simplex(PP, &SD, dim, matrix, point, row, MaxRays);
	if (!vol)
	    vol = s;
	else {
	    eadd(s, vol);
	    free_evalue_refs(s);
	    free(s);
	}
    }
    Polyhedron_Free(L);
    Matrix_Free(FixedRays);
    value_clear(tmp);
    free(SD.F);

    return vol;
}

static evalue *volume_in_domain(Param_Polyhedron *PP, Param_Domain *D,
				unsigned dim, evalue ***matrix,
				struct parameter_point *point, Polyhedron *C,
				int row, Polyhedron *F,
				struct barvinok_options *options)
{
    int nbV;
    Param_Vertices *V;
    evalue *vol;
    int point_computed = 0;

    if (!point) {
	point = non_empty_point(PP, D, F, C, options->MaxRays);
	if (point)
	    point_computed = 1;
    }

    nbV = 0;
    FORALL_PVertex_in_ParamPolyhedron(V, D, PP)
	++nbV;
    END_FORALL_PVertex_in_ParamPolyhedron;

    if (nbV > (dim-row) + 1) {
	if (point && options->volume_triangulate_lift)
	    vol = volume_triangulate_lift(PP, D, dim, matrix, point,
					  row, options->MaxRays);
	else
	    vol = volume_triangulate(PP, D, dim, matrix, point, C,
				     row, F, options);
    } else {
	assert(nbV == (dim-row) + 1);
	vol = volume_simplex(PP, D, dim, matrix, point, row, options->MaxRays);
    }

    if (point_computed)
	parameter_point_free(point);

    return vol;
}

evalue* Param_Polyhedron_Volume(Polyhedron *P, Polyhedron* C,
				struct barvinok_options *options)
{
    evalue ***matrix;
    unsigned nparam = C->Dimension;
    unsigned nvar = P->Dimension - C->Dimension;
    Param_Polyhedron *PP;
    unsigned PP_MaxRays = options->MaxRays;
    unsigned MaxRays;
    int i, j;
    Value fact;
    evalue *vol;
    int nd;
    struct section { Polyhedron *D; evalue *E; } *s;
    Polyhedron **fVD;
    Param_Domain *D, *next;
    Polyhedron *CA, *F;

    if (options->polynomial_approximation == BV_APPROX_SIGN_NONE)
	options->polynomial_approximation = BV_APPROX_SIGN_APPROX;

    if (options->polynomial_approximation != BV_APPROX_SIGN_APPROX) {
	int pa = options->polynomial_approximation;
	assert(pa == BV_APPROX_SIGN_UPPER || pa == BV_APPROX_SIGN_LOWER);

	P = Polyhedron_Flate(P, nparam, pa == BV_APPROX_SIGN_UPPER,
			     options->MaxRays);

	/* Don't deflate/inflate again (on this polytope) */
	options->polynomial_approximation = BV_APPROX_SIGN_APPROX;
	vol = barvinok_enumerate_with_options(P, C, options);
	options->polynomial_approximation = pa;

	Polyhedron_Free(P);
	return vol;
    }

    if (PP_MaxRays & POL_NO_DUAL)
	PP_MaxRays = 0;

    MaxRays = options->MaxRays;
    POL_UNSET(options->MaxRays, POL_INTEGER);

    value_init(fact);
    Factorial(nvar, &fact);

    PP = Polyhedron2Param_Domain(P, C, PP_MaxRays);

    for (nd = 0, D = PP->D; D; ++nd, D = D->next);
    s = ALLOCN(struct section, nd);
    fVD = ALLOCN(Polyhedron *, nd);

    matrix = ALLOCN(evalue **, nvar+1);
    for (i = 0; i < nvar+1; ++i)
	matrix[i] = ALLOCN(evalue *, nvar);

    for (nd = 0, D = PP->D; D; D = next) {
	Polyhedron *rVD = reduce_domain(D->Domain, NULL, NULL, fVD, nd, options);

	next = D->next;

	if (!rVD)
	    continue;

	CA = align_context(D->Domain, P->Dimension, MaxRays);
	F = DomainIntersection(P, CA, options->MaxRays);
	Domain_Free(CA);

	s[nd].D = rVD;
	s[nd].E = volume_in_domain(PP, D, nvar, matrix, NULL, rVD,
				   0, F, options);
	Domain_Free(F);
	evalue_div(s[nd].E, fact);

	++nd;
    }
    options->MaxRays = MaxRays;

    vol = ALLOC(evalue);
    value_init(vol->d);
    value_set_si(vol->d, 0);

    if (nd == 0)
	evalue_set_si(vol, 0, 1);
    else {
	vol->x.p = new_enode(partition, 2*nd, C->Dimension);
	for (i = 0; i < nd; ++i) {
	    EVALUE_SET_DOMAIN(vol->x.p->arr[2*i], s[i].D);
	    value_clear(vol->x.p->arr[2*i+1].d);
	    vol->x.p->arr[2*i+1] = *s[i].E;
	    free(s[i].E);
	    Domain_Free(fVD[i]);
	}
    }
    free(s);
    free(fVD);

    for (i = 0; i < nvar+1; ++i)
	free(matrix[i]);
    free(matrix);
    Param_Polyhedron_Free(PP);
    value_clear(fact);

    return vol;
}
