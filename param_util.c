#include <barvinok/options.h>
#include "param_util.h"

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

Param_Polyhedron *Polyhedron2Param_Polyhedron(Polyhedron *Din, Polyhedron *Cin,
					      struct barvinok_options *options)
{
    unsigned MaxRays = options->MaxRays;
    if (MaxRays & POL_NO_DUAL)
	MaxRays = 0;
    return Polyhedron2Param_Domain(Din, Cin, MaxRays);
}
