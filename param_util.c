#include <barvinok/options.h>
#include "param_util.h"
#include "topcom.h"
#include "config.h"

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

static Param_Polyhedron *PL_P2PP(Polyhedron *Din, Polyhedron *Cin,
					      struct barvinok_options *options)
{
    unsigned MaxRays = options->MaxRays;
    if (MaxRays & (POL_NO_DUAL | POL_INTEGER))
	MaxRays = 0;
    return Polyhedron2Param_Domain(Din, Cin, MaxRays);
}


#ifndef POINTS2TRIANGS_PATH
Param_Polyhedron *TC_P2PP(Polyhedron *P, Polyhedron *C,
			  struct barvinok_options *options)
{
    return NULL;
}
#endif

Param_Polyhedron *Polyhedron2Param_Polyhedron(Polyhedron *P, Polyhedron *C,
					      struct barvinok_options *options)
{
    switch(options->chambers) {
    case BV_CHAMBERS_POLYLIB:
	return PL_P2PP(P, C, options);
	break;
    case BV_CHAMBERS_TOPCOM:
	return TC_P2PP(P, C, options);
	break;
    default:
	assert(0);
    }
}
