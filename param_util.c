#include <barvinok/options.h>
#include "param_util.h"

#define ALLOC(type) (type*)malloc(sizeof(type))
#define ALLOCN(type,n) (type*)malloc((n) * sizeof(type))

void Param_Vertex_Common_Denominator(Param_Vertices *V)
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

/* Plug in the parametric vertex Vertex (nvar x (nparam + 2))
 * in the constraint constraint (1 + nvar + nparam + 1).
 * The result is stored in row (1 + nparam + 1),
 * with the denominator in position 0.
 */
void Param_Inner_Product(Value *constraint, Matrix *Vertex, Value *row)
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

Param_Polyhedron *Polyhedron2Param_Polyhedron(Polyhedron *Din, Polyhedron *Cin,
					      struct barvinok_options *options)
{
    unsigned MaxRays = options->MaxRays;
    if (MaxRays & POL_NO_DUAL)
	MaxRays = 0;
    return Polyhedron2Param_Domain(Din, Cin, MaxRays);
}

/* Compute a dummy Param_Domain that contains all vertices of Param_Domain D
 * (which contains the vertices of P) that lie on the facet obtained by
 * saturating constraint c of P
 */
Param_Domain *Param_Polyhedron_Facet(Param_Polyhedron *PP, Param_Domain *D,
				     Polyhedron *P, int c)
{
    int nv;
    Param_Vertices *V;
    unsigned nparam = PP->V->Vertex->NbColumns-2;
    Vector *row = Vector_Alloc(1+nparam+1);
    Param_Domain *FD = ALLOC(Param_Domain);
    FD->Domain = 0;
    FD->next = 0;

    nv = (PP->nbV - 1)/(8*sizeof(int)) + 1;
    FD->F = ALLOCN(unsigned, nv);
    memset(FD->F, 0, nv * sizeof(unsigned));

    FORALL_PVertex_in_ParamPolyhedron(V, D, PP) /* _i, _ix, _bx internal counters */
	int n;
	Param_Inner_Product(P->Constraint[c], V->Vertex, row->p);
	if (First_Non_Zero(row->p+1, nparam+1) == -1)
	    FD->F[_ix] |= _bx;
    END_FORALL_PVertex_in_ParamPolyhedron;

    Vector_Free(row);

    return FD;
}
