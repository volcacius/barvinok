#include <assert.h>
#include <barvinok/options.h>
#include <barvinok/util.h>
#include "argp.h"
#include "progname.h"
#include "ilp.h"
#include "polysign.h"

int main(int argc, char **argv)
{
    struct barvinok_options *options = barvinok_options_new_with_defaults();
    Polyhedron *P;
    Vector *obj, *affine;
    Vector *minmax;
    Vector *opt;
    enum lp_result res;
    Value one;

    set_program_name(argv[0]);
    argp_parse(&barvinok_argp, argc, argv, 0, 0, options);

    P = Polyhedron_Read(options->MaxRays);
    obj = Vector_Read();
    minmax = Vector_Alloc(2);

    if (obj->Size == P->Dimension) {
	affine = Vector_Alloc(P->Dimension+1);
	Vector_Copy(obj->p, affine->p, P->Dimension);
	value_set_si(affine->p[P->Dimension], 1);
    } else if (obj->Size == P->Dimension+1)
	affine = obj;
    else
	assert(0);

    value_init(one);
    value_set_si(one, 1);
    res = polyhedron_range(P, affine->p, one, &minmax->p[0], &minmax->p[1],
			   options);
    assert(res == lp_ok);
    value_clear(one);

    opt = Polyhedron_Integer_Minimum(P, obj->p, minmax->p[0], minmax->p[1],
			     NULL, NULL, options);
    Vector_Print(stdout, P_VALUE_FMT, opt);
    Vector_Free(opt);

    if (obj != affine)
	Vector_Free(affine);
    Vector_Free(obj);
    Vector_Free(minmax);
    Polyhedron_Free(P);

    barvinok_options_free(options);
}
