#include <unistd.h>
#include <stdlib.h>
#include <barvinok/evalue.h>
#include <barvinok/util.h>
#include <barvinok/barvinok.h>
#include "argp.h"

/* The input of this example program is the same as that of testehrhart
 * in the PolyLib distribution, i.e., a polytope in combined
 * data and parameter space, a context polytope in parameter space
 * and (optionally) the names of the parameters.
 * Both polytopes are in PolyLib notation.
 */

struct argp_option argp_options[] = {
    { "convert",   'c', 0, 0, "convert fractionals to periodics" },
    { "floor",     'f', 0, 0, "convert fractionals to floorings" },
    { "size",      'S' },
    { "series",    's', 0, 0, "compute rational generating function" },
    { "explicit",  'e', 0, 0, "convert rgf to psp" },
    { 0 }
};

struct arguments {
    struct barvinok_options *barvinok;
    int convert;
    int floor;
    int size;
    int series;
    int function;
};

error_t parse_opt(int key, char *arg, struct argp_state *state)
{
    struct arguments *options = (struct arguments*) state->input;

    switch (key) {
    case ARGP_KEY_INIT:
	state->child_inputs[0] = options->barvinok;
	options->convert = 0;
	options->floor = 0;
	options->size = 0;
	options->series = 0;
	options->function = 0;
	break;
    case 'c':
	options->convert = 1;
	break;
    case 'f':
	options->floor = 1;
	break;
    case 'S':
	options->size = 1;
	break;
    case 'e':
	options->function = 1;
	/* fall through */
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
    Polyhedron *A, *C;
    Matrix *M;
    evalue *EP;
    char **param_name;
    struct arguments options;
    static struct argp_child argp_children[] = {
	{ &barvinok_argp,    0,	0,  0 },
	{ 0 }
    };
    static struct argp argp = { argp_options, parse_opt, 0, 0, argp_children };
    struct barvinok_options *bv_options = barvinok_options_new_with_defaults();

    options.barvinok = bv_options;
    argp_parse(&argp, argc, argv, 0, 0, &options);

    M = Matrix_Read();
    A = Constraints2Polyhedron(M, bv_options->MaxRays);
    Matrix_Free(M);
    M = Matrix_Read();
    C = Constraints2Polyhedron(M, bv_options->MaxRays);
    Matrix_Free(M);
    Polyhedron_Print(stdout, P_VALUE_FMT, A);
    Polyhedron_Print(stdout, P_VALUE_FMT, C);
    param_name = Read_ParamNames(stdin, C->Dimension);

    if (options.series) {
	gen_fun *gf;
	gf = barvinok_series_with_options(A, C, bv_options);
	gf->print(std::cout, C->Dimension, param_name);
	puts("");
	if (options.function) {
	    EP = *gf;
	    print_evalue(stdout, EP, param_name);
	    free_evalue_refs(EP);
	    free(EP);
	}
	delete gf;
    } else {
	EP = barvinok_enumerate_with_options(A, C, bv_options);
	print_evalue(stdout, EP, param_name);
	if (options.size)
	    printf("\nSize: %d\n", evalue_size(EP));
	if (options.floor) {
	    fprintf(stderr, "WARNING: floor conversion not supported\n");
	    evalue_frac2floor2(EP, 0);
	    print_evalue(stdout, EP, param_name);
	} else if (options.convert) {
	    evalue_mod2table(EP, C->Dimension);
	    print_evalue(stdout, EP, param_name);
	    if (options.size)
		printf("\nSize: %d\n", evalue_size(EP));
	}
	free_evalue_refs(EP);
	free(EP);
    }

    Free_ParamNames(param_name, C->Dimension);
    Polyhedron_Free(A);
    Polyhedron_Free(C);
    free(bv_options);
    return 0;
}
