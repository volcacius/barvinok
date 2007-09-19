#include <barvinok/barvinok.h>
#include <barvinok/util.h>
#include "argp.h"
#include "progname.h"

/* The input of this example program is similar to that of ehrhart_union
 * in the PolyLib distribution, the difference being that the number of
 * polytopes in the union needs to be specified explicitly.
 * The input starts with this number, followed by this number of
 * polytopes in combined data and parameter space, a context polytope
 * in parameter space and (optionally) the names of the parameters.
 * All polytopes are in PolyLib notation.
 */

struct argp_option argp_options[] = {
    { "series",    's', 0, 0, "compute rational generating function" },
    { 0 }
};

struct arguments {
    int series;
    struct barvinok_options	*barvinok;
};

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
    struct arguments *options = (struct arguments*) state->input;

    switch (key) {
    case ARGP_KEY_INIT:
	state->child_inputs[0] = options->barvinok;
	options->series = 0;
	break;
    case 's':
	options->series = 1;
	break;
    default:
	return ARGP_ERR_UNKNOWN;
    }
    return 0;
}


int main(int argc, char **argv)
{
    Matrix *M;
    Polyhedron *C, *D = NULL;
    int i, npol;
    char **param_name;
    char s[128];
    int c, ind = 0;
    int series = 0;
    struct arguments options;
    static struct argp_child argp_children[] = {
	{ &barvinok_argp,    	0,	0,  		0 },
	{ 0 }
    };
    static struct argp argp = { argp_options, parse_opt, 0, 0, argp_children };
    struct barvinok_options *bv_options = barvinok_options_new_with_defaults();

    options.barvinok = bv_options;
    set_program_name(argv[0]);
    argp_parse(&argp, argc, argv, 0, 0, &options);

    fgets(s, 128, stdin);
    while ((*s=='#') || (sscanf(s, "%d", &npol)<1))
	fgets(s, 128, stdin);

    for (i = 0; i < npol; ++i) {
	Polyhedron *P;
	M = Matrix_Read();
	P = Constraints2Polyhedron(M, bv_options->MaxRays);
	Matrix_Free(M);
	D = DomainConcat(P, D);
    }
    M = Matrix_Read();
    C = Constraints2Polyhedron(M, bv_options->MaxRays);
    Matrix_Free(M);
    Polyhedron_Print(stdout, P_VALUE_FMT, D);
    Polyhedron_Print(stdout, P_VALUE_FMT, C);
    param_name = Read_ParamNames(stdin, C->Dimension);
    if (series) {
	gen_fun *gf;
	gf = barvinok_enumerate_union_series(D, C, bv_options->MaxRays);
	gf->print(std::cout, C->Dimension, param_name);
	puts("");
	delete gf;
    } else {
	evalue *EP;
	EP = barvinok_enumerate_union(D, C, bv_options->MaxRays);
	print_evalue(stdout, EP, param_name);
	evalue_free(EP);
    }
    Free_ParamNames(param_name, C->Dimension);
    Domain_Free(D);
    Polyhedron_Free(C);
    barvinok_options_free(bv_options);
    return 0;
}
