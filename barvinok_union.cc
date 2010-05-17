#include <barvinok/barvinok.h>
#include <barvinok/util.h>
#include "barvinok_union_options.h"

/* The input of this example program is similar to that of ehrhart_union
 * in the PolyLib distribution, the difference being that the number of
 * polytopes in the union needs to be specified explicitly.
 * The input starts with this number, followed by this number of
 * polytopes in combined data and parameter space, a context polytope
 * in parameter space and (optionally) the names of the parameters.
 * All polytopes are in PolyLib notation.
 */


int main(int argc, char **argv)
{
    Matrix *M;
    Polyhedron *C, *D = NULL;
    int i, npol;
    const char **param_name;
    char s[128];
    int c, ind = 0;
    struct union_options *options = union_options_new_with_defaults();

    argc = union_options_parse(options, argc, argv, ISL_ARG_ALL);

    fgets(s, 128, stdin);
    while ((*s=='#') || (sscanf(s, "%d", &npol)<1))
	fgets(s, 128, stdin);

    for (i = 0; i < npol; ++i) {
	Polyhedron *P;
	M = Matrix_Read();
	P = Constraints2Polyhedron(M, options->barvinok->MaxRays);
	Matrix_Free(M);
	D = DomainConcat(P, D);
    }
    M = Matrix_Read();
    C = Constraints2Polyhedron(M, options->barvinok->MaxRays);
    Matrix_Free(M);
    Polyhedron_Print(stdout, P_VALUE_FMT, D);
    Polyhedron_Print(stdout, P_VALUE_FMT, C);
    param_name = Read_ParamNames(stdin, C->Dimension);
    if (options->series) {
	gen_fun *gf;
	gf = barvinok_enumerate_union_series(D, C, options->barvinok->MaxRays);
	gf->print(std::cout, C->Dimension, param_name);
	puts("");
	delete gf;
    } else {
	evalue *EP;
	EP = barvinok_enumerate_union(D, C, options->barvinok->MaxRays);
	print_evalue(stdout, EP, param_name);
	evalue_free(EP);
    }
    Free_ParamNames(param_name, C->Dimension);
    Domain_Free(D);
    Polyhedron_Free(C);
    union_options_free(options);
    return 0;
}
