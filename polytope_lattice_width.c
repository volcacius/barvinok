#include <assert.h>
#include <barvinok/polylib.h>
#include <barvinok/evalue.h>
#include <barvinok/options.h>
#include "lattice_width.h"
#include "argp.h"
#include "progname.h"

struct argp_option argp_options[] = {
    { "direction",	    'd',  0,	0,  "print width directions" },
    { 0 }
};

struct arguments {
    struct barvinok_options *options;
    int direction;
};

error_t parse_opt(int key, char *arg, struct argp_state *state)
{
    struct arguments *arguments = state->input;

    switch (key) {
    case ARGP_KEY_INIT:
	state->child_inputs[0] = arguments->options;
	arguments->direction = 0;
	break;
    case 'd':
	arguments->direction = 1;
	break;
    default:
	return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

int main(int argc, char **argv)
{
    Polyhedron *P, *C;
    Matrix *M;
    char **param_name;
    struct arguments arguments;
    static struct argp_child argp_children[] = {
	{ &barvinok_argp,    0,	0,  0 },
	{ 0 }
    };
    static struct argp argp = { argp_options, parse_opt, 0, 0, argp_children };
    struct barvinok_options *options = barvinok_options_new_with_defaults();
    struct width_direction_array *dirs;
    int i;

    arguments.options = options;

    set_program_name(argv[0]);
    argp_parse(&argp, argc, argv, 0, 0, &arguments);

    M = Matrix_Read();
    assert(M);
    P = Constraints2Polyhedron(M, options->MaxRays);
    Matrix_Free(M);
    M = Matrix_Read();
    assert(M);
    C = Constraints2Polyhedron(M, options->MaxRays);
    Matrix_Free(M);
    param_name = Read_ParamNames(stdin, C->Dimension);

    dirs = Polyhedron_Lattice_Width_Directions(P, C, options);
    for (i = 0; i < dirs->n; ++i) {
	evalue *E;

	Print_Domain(stdout, dirs->wd[i].domain, param_name);
	if (arguments.direction)
	    Vector_Print(stdout, P_VALUE_FMT, dirs->wd[i].dir);
	E = affine2evalue(dirs->wd[i].width->p,
			  dirs->wd[i].width->p[C->Dimension+1],
			  C->Dimension);
	print_evalue(stdout, E, param_name);
	evalue_free(E);
    }
    free_width_direction_array(dirs);

    Free_ParamNames(param_name, C->Dimension);
    Polyhedron_Free(P);
    Polyhedron_Free(C);
    barvinok_options_free(options);

    return 0;
}
