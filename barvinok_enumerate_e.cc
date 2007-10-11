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
#include "omega/convert.h"
#endif
#include "verify.h"
#include "verif_ehrhart.h"
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
    { "pip",   	    	    'p',    0,      0 },
    { "series",     	    's',    0,	    0 },
    { "series",		    's', 0, 0, "compute rational generating function" },
    { "explicit",	    'e', 0, 0, "convert rgf to psp" },
    { "scarf",      	    'S',    0,	    0 },
    { "verbose",    	    'v' },
    { 0 }
};

struct arguments {
    struct verify_options    verify;
    struct convert_options   convert;
    int omega;
    int pip;
    int scarf;
    int series;
    int verbose;
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
    case 'p':
	arguments->pip = 1;
	break;
    case 'v':
	arguments->verbose = 1;
	break;
    default:
	return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

#ifdef HAVE_OMEGA

Polyhedron *Omega_simplify(Polyhedron *P, 
			    unsigned exist, unsigned nparam, char **parms)
{
    varvector varv;
    varvector paramv;
    Relation r = Polyhedron2relation(P, exist, nparam, parms);
    Polyhedron_Free(P);
    return relation2Domain(r, varv, paramv);
}
#else
Polyhedron *Omega_simplify(Polyhedron *P, 
			    unsigned exist, unsigned nparam, char **parms)
{
    return P;
}
#endif

static void verify_results(Polyhedron *P, evalue *EP, int exist, int nparam,
			   verify_options *options);

int main(int argc, char **argv)
{
    Polyhedron *A;
    Matrix *MA;
    char **param_name;
    int exist, nparam, nvar;
    char s[128];
    evalue *EP = NULL;
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
    arguments.pip = 0;
    arguments.scarf = 0;
    arguments.series = 0;
    arguments.function = 0;
    arguments.verbose = 0;

    set_program_name(argv[0]);
    argp_parse(&argp, argc, argv, 0, 0, &arguments);

    if (arguments.series && !arguments.scarf) {
	fprintf(stderr, 
		"--series currently only available if --scarf is specified\n");
	exit(1);
    }

    MA = Matrix_Read();
    A = Constraints2Polyhedron(MA, options->MaxRays);
    Matrix_Free(MA);

    fgets(s, 128, stdin);
    while ((*s=='#') || (sscanf(s, "E %d", &exist)<1))
	fgets(s, 128, stdin);

    fgets(s, 128, stdin);
    while ((*s=='#') || (sscanf(s, "P %d", &nparam)<1))
	fgets(s, 128, stdin);

    /******* Read the options: initialize Min and Max ********/

    if (arguments.verify.verify) {
	verify_options_set_range(&arguments.verify, A->Dimension);
	if (!arguments.verbose)
	    print_solution = 0;
    }

    if (print_solution && arguments.verbose) {
	Polyhedron_Print(stdout, P_VALUE_FMT, A);
	printf("exist: %d, nparam: %d\n", exist, nparam);
    }
    param_name = Read_ParamNames(stdin, nparam);
    nvar = A->Dimension - exist - nparam;
    if (arguments.omega) {
	A = Omega_simplify(A, exist, nparam, param_name);
	assert(!A->next);
	exist = A->Dimension - nvar - nparam;
    }
    if (arguments.series) {
	gen_fun *gf;
	assert(arguments.scarf);
	gf = barvinok_enumerate_scarf_series(A, exist, nparam, options);
	if (print_solution) {
	    gf->print(std::cout, nparam, param_name);
	    puts("");
	}
	if (arguments.function) {
	    EP = *gf;
	    if (print_solution)
		print_evalue(stdout, EP, param_name);
	}
	delete gf;
    } else {
	if (arguments.scarf)
	    EP = barvinok_enumerate_scarf(A, exist, nparam, options);
	else if (arguments.pip && exist > 0)
	    EP = barvinok_enumerate_pip_with_options(A, exist, nparam, options);
	else
	    EP = barvinok_enumerate_e_with_options(A, exist, nparam, options);
	reduce_evalue(EP);
	if (evalue_convert(EP, &arguments.convert, arguments.verbose, nparam,
			   param_name))
	    print_solution = 0;
	if (print_solution)
	    print_evalue(stdout, EP, param_name);
    }
    if (EP && arguments.verify.verify) {
	arguments.verify.params = param_name;
	verify_results(A, EP, exist, nparam, &arguments.verify);
    }
    if (EP)
	evalue_free(EP);
    Free_ParamNames(param_name, nparam);
    Polyhedron_Free(A);
    barvinok_options_free(options);
    return 0;
}

void verify_results(Polyhedron *P, evalue *EP, int exist, int nparam,
		    verify_options *options)
{
    int i;
    int res;
    Value *p, tmp;
    Polyhedron *S, *CS;
    unsigned MaxRays = options->barvinok->MaxRays;
    Polyhedron *C = NULL;
    value_init(tmp);

    p = (Value *)malloc(sizeof(Value) * (P->Dimension+2));
    for(i=0;i<=P->Dimension;i++) {
      value_init(p[i]);
      value_set_si(p[i],0);
    }
    value_init(p[i]);
    value_set_si(p[i],1);

    CS = check_poly_context_scan(P, &C, nparam, options);

    /* S = scanning list of polyhedra */
    S = Polyhedron_Scan(P, C, MaxRays & POL_NO_DUAL ? 0 : MaxRays);

    check_poly_init(C, options);

    /******* CHECK NOW *********/
    res = 0;
    if (S && !check_poly_EP(S, CS, EP, exist, nparam, 0, p, options)) {
      fprintf(stderr,"Check failed !\n");
      res = -1;
    }
      
    if (!options->print_all)
	printf( "\n" );
    
    for(i=0;i<=(P->Dimension+1);i++) 
      value_clear(p[i]);
    free(p);
    value_clear(tmp);
    Domain_Free(S);
    Polyhedron_Free(C);
    if (CS)
	Domain_Free(CS);
}
