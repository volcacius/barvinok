#include <barvinok/util.h>
#include <barvinok/options.h>
#include "scale.h"

/* If a vertex is described by A x + B p + c = 0, then
 * M = [A B] and we want to compute a linear transformation L such
 * that H L = A and H \Z contains both A \Z and B \Z.
 * We compute
 *             [ A B ] = [ H 0 ] [ U_11  U_12 ]
 *                               [ U_21  U_22 ]
 *
 * U_11 is the required linear transformation.
 * Note that this also works if M has more rows than there are variables,
 * i.e., if some rows in M are linear combinations of other rows.
 * These extra rows only affect and H and not U.
 */
static Lattice *extract_lattice(Matrix *M, unsigned nvar)
{
    int row;
    Matrix *H, *Q, *U, *Li;
    Lattice *L;
    int ok;

    left_hermite(M, &H, &Q, &U);
    Matrix_Free(U);

    Li = Matrix_Alloc(nvar+1, nvar+1);
    L = Matrix_Alloc(nvar+1, nvar+1);
    value_set_si(Li->p[nvar][nvar], 1);

    for (row = 0; row < nvar; ++row)
	Vector_Copy(Q->p[row], Li->p[row], nvar);
    Matrix_Free(H);
    Matrix_Free(Q);

    ok = Matrix_Inverse(Li, L);
    assert(ok);
    Matrix_Free(Li);

    return L;
}

/* Returns the smallest (wrt inclusion) lattice that contains both X and Y */
static Lattice *LatticeJoin(Lattice *X, Lattice *Y)
{
    int i;
    int dim = X->NbRows-1;
    Value lcm;
    Value tmp;
    Lattice *L;
    Matrix *M, *H, *U, *Q;

    assert(X->NbColumns-1 == dim);
    assert(Y->NbRows-1 == dim);
    assert(Y->NbColumns-1 == dim);

    value_init(lcm);
    value_init(tmp);

    M = Matrix_Alloc(dim, 2*dim);
    value_lcm(X->p[dim][dim], Y->p[dim][dim], &lcm);

    value_division(tmp, lcm, X->p[dim][dim]);
    for (i = 0; i < dim; ++i)
	Vector_Scale(X->p[i], M->p[i], tmp, dim);
    value_division(tmp, lcm, Y->p[dim][dim]);
    for (i = 0; i < dim; ++i)
	Vector_Scale(Y->p[i], M->p[i]+dim, tmp, dim);

    left_hermite(M, &H, &Q, &U);
    Matrix_Free(M);
    Matrix_Free(Q);
    Matrix_Free(U);

    L = Matrix_Alloc(dim+1, dim+1);
    value_assign(L->p[dim][dim], lcm);
    for (i = 0; i < dim; ++i)
	Vector_Copy(H->p[i], L->p[i], dim);
    Matrix_Free(H);

    value_clear(tmp);
    value_clear(lcm);
    return L;
}

static void Param_Vertex_Common_Denominator(Param_Vertices *V)
{
    unsigned dim;
    Value lcm;
    int i;

    assert(V->Vertex->NbRows > 0);
    dim = V->Vertex->NbColumns-2;

    value_init(lcm);

    value_assign(lcm, V->Vertex->p[0][dim+1]);
    for (i = 1; i < V->Vertex->NbRows; ++i)
	value_lcm(V->Vertex->p[i][dim+1], lcm, &lcm);

    for (i = 0; i < V->Vertex->NbRows; ++i) {
	if (value_eq(V->Vertex->p[i][dim+1], lcm))
	    continue;
	value_division(V->Vertex->p[i][dim+1], lcm, V->Vertex->p[i][dim+1]);
	Vector_Scale(V->Vertex->p[i], V->Vertex->p[i],
		     V->Vertex->p[i][dim+1], dim+1);
	value_assign(V->Vertex->p[i][dim+1], lcm);
    }

    value_clear(lcm);
}

static void Param_Vertex_Image(Param_Vertices *V, Matrix *T)
{
    unsigned nvar  = V->Vertex->NbRows;
    unsigned nparam = V->Vertex->NbColumns - 2;
    Matrix *Vertex;
    int i;

    Param_Vertex_Common_Denominator(V);
    Vertex = Matrix_Alloc(V->Vertex->NbRows, V->Vertex->NbColumns);
    Matrix_Product(T, V->Vertex, Vertex);
    for (i = 0; i < nvar; ++i) {
	value_assign(Vertex->p[i][nparam+1], V->Vertex->p[i][nparam+1]);
	Vector_Normalize(Vertex->p[i], nparam+2);
    }
    Matrix_Free(V->Vertex);
    V->Vertex = Vertex;
}

/* Scales the parametric polyhedron with constraints *P and vertices PP
 * such that the number of integer points can be represented by a polynomial.
 * Both *P and P->Vertex are adapted according to the scaling.
 * The scaling factor is returned in *det.
 * The transformation that maps the new coordinates to the original
 * coordinates is returned in *Lat (if Lat != NULL).
 * The enumerator of the scaled parametric polyhedron should be divided
 * by this number to obtain an approximation of the enumerator of the
 * original parametric polyhedron.
 *
 * The algorithm is described in "Approximating Ehrhart Polynomials using
 * affine transformations" by B. Meister.
 */
void Param_Polyhedron_Scale_Integer_Slow(Param_Polyhedron *PP, Polyhedron **P,
				         Lattice **Lat,
					 Value *det, unsigned MaxRays)
{
    Param_Vertices *V;
    unsigned dim = (*P)->Dimension;
    unsigned nparam;
    unsigned nvar;
    Lattice *L = NULL, *Li;
    Matrix *T;
    Matrix *expansion;
    int i;
    int ok;

    if (!PP->nbV)
	return;

    nparam = PP->V->Vertex->NbColumns - 2;
    nvar = dim - nparam;

    for (V = PP->V; V; V = V->next) {
	Lattice *L2;
	Matrix *M;
	int i, j, n;
	unsigned char *supporting;

	supporting = supporting_constraints(*P, V, &n);
	M = Matrix_Alloc(n, (*P)->Dimension);
	for (i = 0, j = 0; i < (*P)->NbConstraints; ++i)
	    if (supporting[i])
		Vector_Copy((*P)->Constraint[i]+1, M->p[j++], (*P)->Dimension);
	free(supporting);
	L2 = extract_lattice(M, nvar);
	Matrix_Free(M);

	if (!L)
	    L = L2;
	else {
	    Lattice *L3 = LatticeJoin(L, L2);
	    Matrix_Free(L);
	    Matrix_Free(L2);
	    L = L3;
	}
    }

    if (Lat)
	*Lat = Matrix_Copy(L);

    /* apply the variable expansion to the polyhedron (constraints) */
    expansion = Matrix_Alloc(nvar + nparam + 1,  nvar + nparam + 1);
    for (i = 0; i < nvar; ++i)
	Vector_Copy(L->p[i], expansion->p[i], nvar);
    for (i = nvar; i < nvar+nparam+1; ++i)
	value_assign(expansion->p[i][i], L->p[nvar][nvar]);

    *P = Polyhedron_Preimage(*P, expansion, MaxRays);
    Matrix_Free(expansion);

    /* apply the variable expansion to the parametric vertices */
    Li = Matrix_Alloc(nvar+1, nvar+1);
    ok = Matrix_Inverse(L, Li);
    assert(ok);
    Matrix_Free(L);
    assert(value_one_p(Li->p[nvar][nvar]));
    T = Matrix_Alloc(nvar, nvar);
    value_set_si(*det, 1);
    for (i = 0; i < nvar; ++i) {
	value_multiply(*det, *det, Li->p[i][i]);
	Vector_Copy(Li->p[i], T->p[i], nvar);
    }
    Matrix_Free(Li);
    for (V = PP->V; V; V = V->next)
	Param_Vertex_Image(V, T);
    Matrix_Free(T);
}

/* Scales the parametric polyhedron with constraints *P and vertices PP
 * such that the number of integer points can be represented by a polynomial.
 * Both *P and P->Vertex are adapted according to the scaling.
 * The scaling factor is returned in *det.
 * The transformation that maps the new coordinates to the original
 * coordinates is returned in *Lat (if Lat != NULL).
 * The enumerator of the scaled parametric polyhedron should be divided
 * by this number to obtain an approximation of the enumerator of the
 * original parametric polyhedron.
 *
 * The algorithm is described in "Approximating Ehrhart Polynomials using
 * affine transformations" by B. Meister.
 */
void Param_Polyhedron_Scale_Integer_Fast(Param_Polyhedron *PP, Polyhedron **P,
				         Lattice **Lat,
					 Value *det, unsigned MaxRays)
{
  int i;
  int nb_param, nb_vars;
  Vector *denoms;
  Param_Vertices *V;
  Value global_var_lcm;
  Value tmp;
  Matrix *expansion;

  value_set_si(*det, 1);
  if (!PP->nbV)
    return;

  nb_param = PP->D->Domain->Dimension;
  nb_vars = PP->V->Vertex->NbRows;

  /* Scan the vertices and make an orthogonal expansion of the variable
     space */
  /* a- prepare the array of common denominators */
  denoms = Vector_Alloc(nb_vars);
  value_init(global_var_lcm);

  value_init(tmp);
  /* b- scan the vertices and compute the variables' global lcms */
  for (V = PP->V; V; V = V->next) {
    for (i = 0; i < nb_vars; i++) {
      Vector_Gcd(V->Vertex->p[i], nb_param, &tmp);
      Gcd(tmp, V->Vertex->p[i][nb_param+1], &tmp);
      value_division(tmp, V->Vertex->p[i][nb_param+1], tmp);
      Lcm3(denoms->p[i], tmp, &denoms->p[i]);
    }
  }
  value_clear(tmp);

  value_set_si(global_var_lcm, 1);
  for (i = 0; i < nb_vars; i++) {
    value_multiply(*det, *det, denoms->p[i]);
    Lcm3(global_var_lcm, denoms->p[i], &global_var_lcm);
  }

  /* scale vertices */
  for (V = PP->V; V; V = V->next)
    for (i = 0; i < nb_vars; i++) {
      Vector_Scale(V->Vertex->p[i], V->Vertex->p[i], denoms->p[i], nb_param+1);
      Vector_Normalize(V->Vertex->p[i], nb_param+2);
    }

  /* the expansion can be actually writen as global_var_lcm.L^{-1} */
  /* this is equivalent to multiply the rows of P by denoms_det */
  for (i = 0; i < nb_vars; i++)
    value_division(denoms->p[i], global_var_lcm, denoms->p[i]);

  /* OPT : we could use a vector instead of a diagonal matrix here (c- and d-).*/
  /* c- make the quick expansion matrix */
  expansion = Matrix_Alloc(nb_vars+nb_param+1, nb_vars+nb_param+1);
  for (i = 0; i < nb_vars; i++)
    value_assign(expansion->p[i][i], denoms->p[i]);
  for (i = nb_vars; i < nb_vars+nb_param+1; i++)
    value_assign(expansion->p[i][i], global_var_lcm);

  /* d- apply the variable expansion to the polyhedron */
  if (P)
    *P = Polyhedron_Preimage(*P, expansion, MaxRays);

  if (Lat) {
    Lattice *L = Matrix_Alloc(nb_vars+1, nb_vars+1);
    for (i = 0; i < nb_vars; ++i)
	value_assign(L->p[i][i], denoms->p[i]);
    value_assign(L->p[nb_vars][nb_vars], global_var_lcm);
    *Lat = L;
  }

  Matrix_Free(expansion);
  value_clear(global_var_lcm);
  Vector_Free(denoms);
}

/* adapted from mpolyhedron_inflate in PolyLib */
static Polyhedron *Polyhedron_Inflate(Polyhedron *P, unsigned nparam,
				      unsigned MaxRays)
{
    Value sum;
    int nvar = P->Dimension - nparam;
    Matrix *C = Polyhedron2Constraints(P);
    Polyhedron *P2;
    int i, j;

    value_init(sum);
    /* subtract the sum of the negative coefficients of each inequality */
    for (i = 0; i < C->NbRows; ++i) {
	value_set_si(sum, 0);
	for (j = 0; j < nvar; ++j)
	    if (value_neg_p(C->p[i][1+j]))
		value_addto(sum, sum, C->p[i][1+j]);
	value_subtract(C->p[i][1+P->Dimension], C->p[i][1+P->Dimension], sum);
    }
    value_clear(sum);
    P2 = Constraints2Polyhedron(C, MaxRays);
    Matrix_Free(C);
    return P2;
}

static void linear_min(Polyhedron *D, Value *obj, Value *min)
{
    int i;
    Value tmp;
    value_init(tmp);
    POL_ENSURE_VERTICES(D);
    for (i = 0; i < D->NbRays; ++i) {
	Inner_Product(obj, D->Ray[i]+1, D->Dimension, &tmp);
	mpz_cdiv_q(tmp, tmp, D->Ray[i][1+D->Dimension]);
	if (!i || value_lt(tmp, *min))
	    value_assign(*min, tmp);
    }
    value_clear(tmp);
}

static void Vector_Oppose(Value *p1, Value *p2, unsigned len)
{
    unsigned i;

    for (i = 0; i < len; ++i)
	value_oppose(p2[i], p1[i]);
}

static Polyhedron *inflate_deflate_domain(Lattice *L, unsigned MaxRays)
{
    unsigned nvar = L->NbRows-1;
    int i;
    Matrix *M;
    Polyhedron *D;

    M = Matrix_Alloc(2*nvar, 1+nvar+1);
    for (i = 0; i < nvar; ++i) {
	value_set_si(M->p[2*i][0], 1);
	Vector_Copy(L->p[i], M->p[2*i]+1, nvar);
	Vector_Normalize(M->p[2*i]+1, nvar);

	value_set_si(M->p[2*i+1][0], 1);
	Vector_Oppose(L->p[i], M->p[2*i+1]+1, nvar);
	value_assign(M->p[2*i+1][1+nvar], L->p[nvar][nvar]);
	Vector_Normalize(M->p[2*i+1]+1, nvar+1);
	value_decrement(M->p[2*i+1][1+nvar], M->p[2*i+1][1+nvar]);
    }
    D = Constraints2Polyhedron(M, MaxRays);
    Matrix_Free(M);

    return D;
}

static Polyhedron *Polyhedron_Inflate4(Polyhedron *P, Lattice *L,
				       unsigned nparam, unsigned MaxRays)
{
    int i;
    unsigned nvar = P->Dimension - nparam;
    Vector *obj;
    Value min;
    Matrix *C;
    Polyhedron *D;
    Polyhedron *P2;

    if (!L)
	return Polyhedron_Inflate(P, nparam, MaxRays);

    D = inflate_deflate_domain(L, MaxRays);
    value_init(min);
    obj = Vector_Alloc(nvar);
    C = Polyhedron2Constraints(P);

    for (i = 0; i < C->NbRows; ++i) {
	Vector_Copy(C->p[i]+1, obj->p, nvar);
	linear_min(D, obj->p, &min);
	value_subtract(C->p[i][1+P->Dimension], C->p[i][1+P->Dimension], min);
    }

    Polyhedron_Free(D);
    P2 = Constraints2Polyhedron(C, MaxRays);
    Matrix_Free(C);
    Vector_Free(obj);
    value_clear(min);
    return P2;
}

/* adapted from mpolyhedron_deflate in PolyLib */
static Polyhedron *Polyhedron_Deflate(Polyhedron *P, unsigned nparam,
				      unsigned MaxRays)
{
    Value sum;
    int nvar = P->Dimension - nparam;
    Matrix *C = Polyhedron2Constraints(P);
    Polyhedron *P2;
    int i, j;

    value_init(sum);
    /* subtract the sum of the positive coefficients of each inequality */
    for (i = 0; i < C->NbRows; ++i) {
	value_set_si(sum, 0);
	for (j = 0; j < nvar; ++j)
	    if (value_pos_p(C->p[i][1+j]))
		value_addto(sum, sum, C->p[i][1+j]);
	value_subtract(C->p[i][1+P->Dimension], C->p[i][1+P->Dimension], sum);
    }
    value_clear(sum);
    P2 = Constraints2Polyhedron(C, MaxRays);
    Matrix_Free(C);
    return P2;
}

static Polyhedron *Polyhedron_Deflate4(Polyhedron *P, Lattice *L,
				       unsigned nparam, unsigned MaxRays)
{
    unsigned nvar = P->Dimension - nparam;
    Vector *obj;
    Value min;
    Matrix *C;
    Polyhedron *D;
    Polyhedron *P2;
    int i;

    if (!L)
	return Polyhedron_Deflate(P, nparam, MaxRays);

    D = inflate_deflate_domain(L, MaxRays);
    value_init(min);
    obj = Vector_Alloc(nvar);
    C = Polyhedron2Constraints(P);

    for (i = 0; i < C->NbRows; ++i) {
	Vector_Oppose(C->p[i]+1, obj->p, nvar);
	linear_min(D, obj->p, &min);
	value_addto(C->p[i][1+P->Dimension], C->p[i][1+P->Dimension], min);
    }

    Polyhedron_Free(D);
    P2 = Constraints2Polyhedron(C, MaxRays);
    Matrix_Free(C);
    Vector_Free(obj);
    value_clear(min);
    return P2;
}

static void Param_Polyhedron_Scale(Param_Polyhedron *PP, Polyhedron **P,
				   Lattice **L,
				   Value *det, struct barvinok_options *options)
{
    if (options->scale_flags & BV_APPROX_SCALE_FAST)
	Param_Polyhedron_Scale_Integer_Fast(PP, P, L, det, options->MaxRays);
    else
	Param_Polyhedron_Scale_Integer_Slow(PP, P, L, det, options->MaxRays);
}

Polyhedron *scale_init(Polyhedron *P, Polyhedron *C, struct scale_data *scaling,
		       struct barvinok_options *options)
{
    unsigned nparam = C->Dimension;
    Polyhedron *Porig = P;
    Polyhedron *T;
    int scale_narrow = options->scale_flags & BV_APPROX_SCALE_NARROW;
    Lattice *L = NULL;

    value_init(scaling->det);
    value_set_si(scaling->det, 1);
    scaling->save_approximation = options->polynomial_approximation;

    if (options->polynomial_approximation == BV_APPROX_SIGN_NONE ||
        options->polynomial_approximation == BV_APPROX_SIGN_APPROX)
	return P;

    if (scale_narrow) {
	Param_Polyhedron *PP;
	unsigned PP_MaxRays = options->MaxRays;
	if (PP_MaxRays & POL_NO_DUAL)
	    PP_MaxRays = 0;
	PP = Polyhedron2Param_Domain(P, C, PP_MaxRays);
	Param_Polyhedron_Scale(PP, &P, &L, &scaling->det, options);
	Param_Polyhedron_Free(PP);
	/* Don't scale again (on this polytope) */
	options->approximation_method = BV_APPROX_NONE;
    }
    T = P;
    if (options->polynomial_approximation == BV_APPROX_SIGN_UPPER)
	P = Polyhedron_Inflate4(P, L, nparam, options->MaxRays);
    if (options->polynomial_approximation == BV_APPROX_SIGN_LOWER)
	P = Polyhedron_Deflate4(P, L, nparam, options->MaxRays);

    /* Don't deflate/inflate again (on this polytope) */
    options->polynomial_approximation = BV_APPROX_SIGN_NONE;

    if (T != Porig)
	Polyhedron_Free(T);
    if (L)
	Matrix_Free(L);
    return P;
}

Polyhedron *scale(Param_Polyhedron *PP, Polyhedron *P,
		  struct scale_data *scaling, int free_P,
		  struct barvinok_options *options)
{
    Polyhedron *T = P;
    unsigned MaxRays;

    MaxRays = options->MaxRays;
    POL_UNSET(options->MaxRays, POL_INTEGER);

    Param_Polyhedron_Scale(PP, &T, NULL, &scaling->det, options);
    if (free_P)
	Polyhedron_Free(P);

    options->MaxRays = MaxRays;

    return T;
}

void scale_finish(evalue *e, struct scale_data *scaling,
		  struct barvinok_options *options)
{
    if (value_notone_p(scaling->det))
	evalue_div(e, scaling->det);
    value_clear(scaling->det);
    /* reset options that may have been changed */
    options->approximation_method = BV_APPROX_SCALE;
    options->polynomial_approximation = scaling->save_approximation;
}
