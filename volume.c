#include <barvinok/polylib.h>
#include <barvinok/options.h>
#include <barvinok/util.h>
#include "reduce_domain.h"
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
	    fprintf(stderr, "%d %d: ", i, j);
	    print_evalue(stderr, matrix[i][cols[j]], param_names);
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
				evalue **point, int *removed,
				int row, Polyhedron *F, unsigned MaxRays);

static evalue *volume_triangulate(Param_Polyhedron *PP, Param_Domain *D,
				  unsigned dim, evalue ***matrix,
				  evalue **point, int *removed,
				  int row, Polyhedron *F, unsigned MaxRays)
{
    int nbV, j;
    Param_Vertices *V;
    Value denom;
    evalue *tmp;
    evalue *vol;
    evalue mone;

    value_init(mone.d);
    evalue_set_si(&mone, -1, 1);

    nbV = 0;
    FORALL_PVertex_in_ParamPolyhedron(V, D, PP)
	for (j = 0; j < dim; ++j) {
	    tmp = vertex2evalue(V->Vertex->p[j], V->Vertex->NbColumns - 2);
	    if (nbV == 0)
		matrix[row][j] = tmp;
	    else {
		eadd(tmp, matrix[row][j]);
		free_evalue_refs(tmp);
		free(tmp);
	    }
	}
	++nbV;
    END_FORALL_PVertex_in_ParamPolyhedron;
    value_init(denom);
    value_set_si(denom, nbV);
    for (j = 0; j < dim; ++j)
	evalue_div(matrix[row][j], denom);
    value_clear(denom);

    if (row == 0) {
	for (j = 0; j < dim; ++j)
	    emul(&mone, matrix[row][j]);
    } else {
	for (j = 0; j < dim; ++j)
	    eadd(matrix[0][j], matrix[row][j]);
    }

    vol = NULL;
    POL_ENSURE_FACETS(F);
    for (j = F->NbEq; j < F->NbConstraints; ++j) {
	Polyhedron *FF;
	Param_Domain *FD;
	if (First_Non_Zero(F->Constraint[j]+1, dim) == -1)
	    continue;
	FF = facet(F, j, MaxRays);
	FD = face_vertices(PP, D, F, j);
	tmp = volume_in_domain(PP, FD, dim, matrix, point, removed,
			       row+1, FF, MaxRays);
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

    for (j = 0; j < dim; ++j) {
	free_evalue_refs(matrix[row][j]);
	free(matrix[row][j]);
    }

    free_evalue_refs(&mone);
    return vol;
}

static evalue *volume_in_domain(Param_Polyhedron *PP, Param_Domain *D,
				unsigned dim, evalue ***matrix,
				evalue **point, int *removed,
				int row, Polyhedron *F, unsigned MaxRays)
{
    int nbV;
    Param_Vertices *V;
    evalue *vol, *val;
    int i, j;
    evalue mone;
    int empty = 0;

    nbV = 0;
    FORALL_PVertex_in_ParamPolyhedron(V, D, PP)
	++nbV;
    END_FORALL_PVertex_in_ParamPolyhedron;

    if (nbV > (dim-row) + 1)
	return volume_triangulate(PP, D, dim, matrix, point, removed,
				  row, F, MaxRays);

    assert(nbV == (dim-row) + 1);

    value_init(mone.d);
    evalue_set_si(&mone, -1, 1);

    i = row;
    FORALL_PVertex_in_ParamPolyhedron(V, D, PP) /* _ix, _bx internal counters */
	if (removed[_ix] & _bx)
	    empty = 1;
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

    val = evalue_substitute(vol, point);

    assert(value_notzero_p(val->d));
    if (value_zero_p(val->x.n)) {
	evalue *tmp;
	assert(empty);
	tmp = val;
	val = vol;
	vol = tmp;
    } else if (value_neg_p(val->x.n))
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

/* Plug in the parametric vertex V in the constraint constraint.
 * The result is stored in row, with the denominator in position 0.
 */
static void Param_Inner_Product(Value *constraint, Param_Vertices *V,
				Value *row)
{
    unsigned nparam = V->Vertex->NbColumns - 2;
    unsigned nvar = V->Vertex->NbRows;
    int j;
    Value tmp, tmp2;

    value_set_si(row[0], 1);
    Vector_Set(row+1, 0, nparam+1);

    value_init(tmp);
    value_init(tmp2);

    for (j = 0 ; j < nvar; ++j) {
	value_set_si(tmp, 1);
	value_assign(tmp2,  constraint[1+j]);
	if (value_ne(row[0], V->Vertex->p[j][nparam+1])) {
	    value_assign(tmp, row[0]);
	    value_lcm(row[0], V->Vertex->p[j][nparam+1], &row[0]);
	    value_division(tmp, row[0], tmp);
	    value_multiply(tmp2, tmp2, row[0]);
	    value_division(tmp2, tmp2, V->Vertex->p[j][nparam+1]);
	}
	Vector_Combine(row+1, V->Vertex->p[j], row+1, tmp, tmp2, nparam+1);
    }
    value_set_si(tmp, 1);
    Vector_Combine(row+1, constraint+1+nvar, row+1, tmp, row[0], nparam+1);

    value_clear(tmp);
    value_clear(tmp2);
}

/* Computes point in pameter space where polyhedron is non-empty.
 * For each of the parametric vertices, and each of the facets
 * not (always) containing the vertex, we remove the parameter
 * values for which the facet does contain the vertex.
 */
static evalue **non_empty_point(Param_Polyhedron *PP, Param_Domain *D,
				Polyhedron *P, int **removed, unsigned MaxRays)
{
    Param_Vertices *V;
    unsigned dim = P->Dimension;
    unsigned nparam = D->Domain->Dimension;
    unsigned nvar = dim - nparam;
    Polyhedron *RD, *cut, *tmp;
    Matrix *M;
    evalue **point;
    int i, j;
    unsigned cut_MaxRays = MaxRays;
    int nv;

    nv = (PP->nbV - 1)/(8*sizeof(int)) + 1;
    removed[0] = ALLOCN(unsigned, nv);
    memset(removed[0], 0, nv * sizeof(unsigned));

    POL_UNSET(cut_MaxRays, POL_INTEGER);

    M = Matrix_Alloc(1, nparam+2);
    RD = NULL;
    FORALL_PVertex_in_ParamPolyhedron(V, D, PP) /* _ix, _bx internal counters */
	Polyhedron *VD = D->Domain;
	for (i = P->NbEq; i < P->NbConstraints; ++i) {
	    if (First_Non_Zero(P->Constraint[i]+1, nvar) == -1)
		continue;
	    Param_Inner_Product(P->Constraint[i], V, M->p[0]);
	    if (First_Non_Zero(M->p[0]+1, nparam) == -1)
		/* supporting facet,
		 * or non-supporting facet independent of params
		 */
		continue;
	    value_set_si(M->p[0][0], 0);
	    cut = Constraints2Polyhedron(M, cut_MaxRays);
	    tmp = DomainDifference(VD, cut, MaxRays);
	    if (VD != D->Domain)
		Domain_Free(VD);
	    VD = tmp;
	    Polyhedron_Free(cut);
	}
	if (VD != D->Domain)
	    VD = DomainConstraintSimplify(VD, MaxRays);
	POL_ENSURE_FACETS(VD);
	if (emptyQ(VD)) {
	    Domain_Free(VD);
	    removed[0][_ix] |= _bx;
	} else {
	    if (!RD)
		RD = VD;
	    else {
		tmp = DomainIntersection(RD, VD, MaxRays);
		if (VD != D->Domain)
		    Domain_Free(VD);
		if (RD != D->Domain)
		    Domain_Free(RD);
		RD = tmp;
	    }
	}
    END_FORALL_PVertex_in_ParamPolyhedron;
    Matrix_Free(M);

    if (!RD)
	point = NULL;
    else {
	POL_ENSURE_VERTICES(RD);
	point = ALLOCN(evalue *, nvar);
	for (i = 0; i < RD->NbRays; ++i)
	    if (value_notzero_p(RD->Ray[i][1+nparam]))
		break;
	assert(i < RD->NbRays);
	for (j = 0; j < nparam; ++j) {
	    point[j] = ALLOC(evalue);
	    value_init(point[j]->d);
	    evalue_set(point[j], RD->Ray[i][1+j], RD->Ray[i][1+nparam]);
	}
    }

    if (RD != D->Domain)
	Domain_Free(RD);

    return point;
}

evalue* Param_Polyhedron_Volume(Polyhedron *P, Polyhedron* C,
				struct barvinok_options *options)
{
    evalue ***matrix;
    unsigned nparam = C->Dimension;
    unsigned nvar = P->Dimension - C->Dimension;
    Param_Polyhedron *PP;
    unsigned PP_MaxRays = options->MaxRays;
    unsigned rat_MaxRays = options->MaxRays;
    int i, j;
    Value fact;
    evalue *vol;
    int nd;
    struct section { Polyhedron *D; evalue *E; } *s;
    Polyhedron **fVD;
    Param_Domain *D, *next;
    Polyhedron *CA, *F;

    if (PP_MaxRays & POL_NO_DUAL)
	PP_MaxRays = 0;

    POL_UNSET(rat_MaxRays, POL_INTEGER);

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
	evalue **point;
	int *removed;
	Polyhedron *rVD = reduce_domain(D->Domain, NULL, NULL, fVD, nd, options);

	next = D->next;

	if (!rVD)
	    continue;

	CA = align_context(D->Domain, P->Dimension, options->MaxRays);
	F = DomainIntersection(P, CA, rat_MaxRays);
	Domain_Free(CA);

	point = non_empty_point(PP, D, F, &removed, options->MaxRays);
	if (!point) {
	    free(removed);
	    Domain_Free(fVD[nd]);
	    Domain_Free(rVD);
	    Domain_Free(F);
	    continue;
	}

	s[nd].D = rVD;
	s[nd].E = volume_in_domain(PP, D, nvar, matrix, point, removed,
				   0, F, rat_MaxRays);
	Domain_Free(F);
	evalue_div(s[nd].E, fact);

	for (i = 0; i < nparam; ++i) {
	    free_evalue_refs(point[i]);
	    free(point[i]);
	}
	free(point);
	free(removed);

	++nd;
    }

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
