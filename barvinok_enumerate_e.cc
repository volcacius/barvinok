#include <ctype.h>
#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include <barvinok/util.h>
#include <barvinok/barvinok.h>
#include "argp.h"
#include "progname.h"
#include "error.h"
#include "config.h"
#ifdef HAVE_OMEGA
#include "omega_interface/convert.h"
#include "omega_interface/count.h"
#endif
#include "skewed_genfun.h"
#include "verify.h"
#include "verif_ehrhart.h"
#include "verify_series.h"
#include "evalue_convert.h"

/* The input of this example program is a polytope in combined
 * data and parameter space followed by two lines indicating
 * the number of existential variables and parameters respectively.
 * The first lines starts with "E ", followed by a number.
 * The second lines starts with "P ", followed by a number.
 * These two lines are (optionally) followed by the names of the parameters.
 * The polytope is in PolyLib notation.
 */

struct argp_option argp_options[] = {
    { "omega",      	    'o',    0,      0 },
    { "parker", 	    'P',    0,      0 },
    { "pip",   	    	    'p',    0,      0 },
    { "series",     	    's',    0,	    0 },
    { "series",		    's', 0, 0, "compute rational generating function" },
    { "explicit",	    'e', 0, 0, "convert rgf to psp" },
    { "scarf",      	    'S',    0,	    0 },
    { 0 }
};

struct arguments {
    struct verify_options    verify;
    struct convert_options   convert;
    int omega;
    int parker;
    int pip;
    int scarf;
    int series;
    int function;
};

error_t parse_opt(int key, char *arg, struct argp_state *state)
{
    struct arguments *arguments = (struct arguments *)(state->input);

    switch (key) {
    case ARGP_KEY_INIT:
	state->child_inputs[0] = arguments->verify.barvinok;
	state->child_inputs[1] = &arguments->verify;
	state->child_inputs[2] = &arguments->convert;
	break;
    case 'e':
	arguments->function = 1;
	/* fall through */
    case 's':
	arguments->series = 1;
	break;
    case 'S':
	arguments->scarf = 1;
	break;
    case 'o':
#ifdef HAVE_OMEGA
	arguments->omega = 1;
#else
	error(0, 0, "--omega option not supported");
#endif
	break;
    case 'P':
#ifdef USE_PARKER
	arguments->parker = 1;
#else
	error(0, 0, "--parker option not supported");
#endif
	break;
    case 'p':
	arguments->pip = 1;
	break;
    default:
	return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

#ifdef HAVE_OMEGA

Polyhedron *Omega_simplify(Polyhedron *P, 
			    unsigned exist, unsigned nparam, const char **parms,
			    unsigned MaxRays)
{
    varvector varv;
    varvector paramv;
    Relation r = Polyhedron2relation(P, exist, nparam, parms);
    Polyhedron_Free(P);
    return relation2Domain(r, varv, paramv, MaxRays);
}
#else
Polyhedron *Omega_simplify(Polyhedron *P, 
			    unsigned exist, unsigned nparam, const char **parms,
			    unsigned MaxRays)
{
    return P;
}

evalue *barvinok_enumerate_parker(Polyhedron *P,
					unsigned exist, unsigned nparam,
					unsigned MaxRays)
{
    assert(0);
    return NULL;
}
#endif

static void verify_results(Polyhedron *P, evalue *EP, gen_fun *gf,
			   int exist, int nparam,
			   arguments *options);

static char *next_line(FILE *input, char *line, unsigned len)
{
	char *p;

	do {
		if (!(p = fgets(line, len, input)))
			return NULL;
		while (isspace(*p) && *p != '\n')
			++p;
	} while (*p == '#' || *p == '\n');

	return p;
}

int main(int argc, char **argv)
{
    Polyhedron *A;
    Matrix *MA;
    const char **param_name;
    int exist, nparam, nvar;
    char s[128];
    evalue *EP = NULL;
    gen_fun *gf = NULL;
    int print_solution = 1;
    struct arguments arguments;
    static struct argp_child argp_children[] = {
	{ &barvinok_argp,    	0,	0,  			0 },
	{ &verify_argp,    	0,	"verification",		BV_GRP_LAST+1 },
	{ &convert_argp,    	0,	"output conversion",    BV_GRP_LAST+2 },
	{ 0 }
    };
    static struct argp argp = { argp_options, parse_opt, 0, 0, argp_children };
    struct barvinok_options *options = barvinok_options_new_with_defaults();

    arguments.verify.barvinok = options;
    arguments.omega = 0;
    arguments.parker = 0;
    arguments.pip = 0;
    arguments.scarf = 0;
    arguments.series = 0;
    arguments.function = 0;

    set_program_name(argv[0]);
    argp_parse(&argp, argc, argv, 0, 0, &arguments);

    MA = Matrix_Read();
    A = Constraints2Polyhedron(MA, options->MaxRays);
    Matrix_Free(MA);

    exist = -1;
    while (next_line(stdin, s, sizeof(s)))
	if (sscanf(s, "E %d", &exist) == 1)
	    break;
    assert(exist >= 0);

    nparam = -1;
    while (next_line(stdin, s, sizeof(s)))
	if (sscanf(s, "P %d", &nparam) == 1)
	    break;
    assert(nparam >= 0);

    /******* Read the options: initialize Min and Max ********/

    if (arguments.verify.verify) {
	verify_options_set_range(&arguments.verify, A->Dimension);
	if (!options->verbose)
	    print_solution = 0;
    }

    if (print_solution && options->verbose) {
	Polyhedron_Print(stdout, P_VALUE_FMT, A);
	printf("exist: %d, nparam: %d\n", exist, nparam);
    }
    param_name = Read_ParamNames(stdin, nparam);
    nvar = A->Dimension - exist - nparam;
    if (arguments.omega) {
	A = Omega_simplify(A, exist, nparam, param_name, options->MaxRays);
	assert(!A->next);
	exist = A->Dimension - nvar - nparam;
    }
    if (arguments.series) {
	if (arguments.scarf)
	    gf = barvinok_enumerate_scarf_series(A, exist, nparam, options);
	else
	    gf = barvinok_enumerate_e_series(A, exist, nparam, options);
	if (print_solution) {
	    gf->print(std::cout, nparam, param_name);
	    puts("");
	}
	if (arguments.function) {
	    EP = *gf;
	    if (print_solution)
		print_evalue(stdout, EP, param_name);
	}
    } else {
	if (arguments.parker)
	    EP = barvinok_enumerate_parker(A, A->Dimension-nparam-exist,
						nparam, options->MaxRays);
	else if (arguments.scarf)
	    EP = barvinok_enumerate_scarf(A, exist, nparam, options);
	else if (arguments.pip && exist > 0)
	    EP = barvinok_enumerate_pip_with_options(A, exist, nparam, options);
	else
	    EP = barvinok_enumerate_e_with_options(A, exist, nparam, options);
	reduce_evalue(EP);
	if (evalue_convert(EP, &arguments.convert, options->verbose, nparam,
			   param_name))
	    print_solution = 0;
	if (print_solution)
	    print_evalue(stdout, EP, param_name);
    }
    if (arguments.verify.verify) {
	arguments.verify.params = param_name;
	verify_results(A, EP, gf, exist, nparam, &arguments);
    }
    if (gf)
	delete gf;
    if (EP)
	evalue_free(EP);

    if (options->print_stats)
	barvinok_stats_print(options->stats, stdout);

    Free_ParamNames(param_name, nparam);
    Polyhedron_Free(A);
    barvinok_options_free(options);
    return 0;
}

void verify_results(Polyhedron *P, evalue *EP, gen_fun *gf,
		       int exist, int nparam,
		       arguments *options)
{
    int i;
    int res = 0;
    Vector *p;
    Value tmp;
    Polyhedron *S, *CS;
    unsigned MaxRays = options->verify.barvinok->MaxRays;
    Polyhedron *C = NULL;
    value_init(tmp);

    p = Vector_Alloc(P->Dimension+2);
    value_set_si(p->p[P->Dimension+1], 1);

    CS = check_poly_context_scan(P, &C, nparam, &options->verify);
    if (!C)
	C = Universe_Polyhedron(nparam);

    /* S = scanning list of polyhedra */
    S = Polyhedron_Scan(P, C, MaxRays & POL_NO_DUAL ? 0 : MaxRays);

    check_poly_init(C, &options->verify);

    /******* CHECK NOW *********/
    if (S) {
	if (!options->series || options->function) {
	    if (!check_poly_EP(S, CS, EP, exist, nparam, 0, p->p,
				&options->verify))
		res = -1;
	} else {
	    skewed_gen_fun *sgf = new skewed_gen_fun(new gen_fun(gf));
	    if (!check_poly_gf(S, CS, sgf, exist, nparam, 0, p->p,
				&options->verify))
		res = -1;
	    delete sgf;
	}
    }
      
    if (!options->verify.print_all)
	printf( "\n" );
    
    Vector_Free(p);
    value_clear(tmp);
    Domain_Free(S);
    Polyhedron_Free(C);
    if (CS)
	Domain_Free(CS);
}
