#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include <barvinok/util.h>
#include <barvinok/barvinok.h>
#include "argp.h"
#include "error.h"
#include "config.h"
#ifdef HAVE_OMEGA
#include "omega/convert.h"
#endif
#include "verify.h"

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
    { "scarf",      	    'S',    0,	    0 },
    { "convert",    	    'c',    0,	    0 },
    { "floor",      	    'f',    0,	    0 },
    { "range-reduction",    'R',    0,	    0 },
    { 0 }
};

struct arguments {
    struct barvinok_options *options;
    struct verify_options    verify;
    int range;
    int convert;
    int omega;
    int pip;
    int scarf;
    int series;
    int floor;
};

error_t parse_opt(int key, char *arg, struct argp_state *state)
{
    struct arguments *arguments = (struct arguments *)(state->input);

    switch (key) {
    case ARGP_KEY_INIT:
	state->child_inputs[0] = arguments->options;
	state->child_inputs[1] = &arguments->verify;
	break;
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
#ifdef HAVE_PIPLIB
	arguments->pip = 1;
#else
	error(0, 0, "--pip option not supported");
#endif
	break;
    case 'f':
	arguments->floor = 1;
	break;
    case 'c':
	arguments->convert = 1;
	break;
    case 'R':
	arguments->range = 1;
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

/* define this to continue the test after first error found */
/* #define DONT_BREAK_ON_ERROR */

static Value min_val, max_val;

static char **params;

static int st;

static int check_poly(Polyhedron *S, Polyhedron *C, evalue *EP,
		      int exist, int nparam, int pos, Value *z, int print_all);
static void verify_results(Polyhedron *P, evalue *EP, int exist, int nparam, 
			   int m, int M, int print_all, unsigned MaxRays);

int main(int argc, char **argv)
{
    Polyhedron *A;
    Matrix *MA;
    char **param_name;
    int exist, nparam, nvar;
    char s[128];
    evalue *EP;
    int print_solution = 1;
    struct arguments arguments;
    static struct argp_child argp_children[] = {
	{ &barvinok_argp,    	0,	0,  		0 },
	{ &verify_argp,    	0,	"verification",	1 },
	{ 0 }
    };
    static struct argp argp = { argp_options, parse_opt, 0, 0, argp_children };
    struct barvinok_options *options = barvinok_options_new_with_defaults();

    arguments.options = options;
    arguments.range = 0;
    arguments.convert = 0;
    arguments.omega = 0;
    arguments.pip = 0;
    arguments.scarf = 0;
    arguments.series = 0;
    arguments.floor = 0;

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
    verify_options_set_range(&arguments.verify, A);

    if (arguments.verify.verify)
	print_solution = 0;

    if (print_solution) {
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
	barvinok_options *options = barvinok_options_new_with_defaults();
	assert(arguments.scarf);
	gf = barvinok_enumerate_scarf_series(A, exist, nparam, options);
	if (print_solution) {
	    gf->print(std::cout, nparam, param_name);
	    puts("");
	}
	delete gf;
	barvinok_options_free(options);
    } else {
	if (arguments.scarf) {
	    barvinok_options *options = barvinok_options_new_with_defaults();
	    EP = barvinok_enumerate_scarf(A, exist, nparam, options);
	    barvinok_options_free(options);
	} else if (arguments.pip && exist > 0)
	    EP = barvinok_enumerate_pip_with_options(A, exist, nparam, options);
	else
	    EP = barvinok_enumerate_e_with_options(A, exist, nparam, options);
	reduce_evalue(EP);
	evalue_combine(EP);
	if (arguments.range)
	    evalue_range_reduction(EP);
	if (print_solution)
	    print_evalue(stdout, EP, param_name);
	if (arguments.floor) {
	    fprintf(stderr, "WARNING: floor conversion not supported\n");
	    evalue_frac2floor2(EP, 0);
	    if (print_solution)
		print_evalue(stdout, EP, param_name);
	} else if (arguments.convert) {
	    evalue_mod2table(EP, nparam);
	    if (print_solution)
		print_evalue(stdout, EP, param_name);
	}
	if (arguments.verify.verify)
	    verify_results(A, EP, exist, nparam, arguments.verify.m,
			    arguments.verify.M, arguments.verify.print_all,
			    options->MaxRays);
	free_evalue_refs(EP);
	free(EP);
    }
    Free_ParamNames(param_name, nparam);
    Polyhedron_Free(A);
    return 0;
}

void verify_results(Polyhedron *P, evalue *EP, int exist, int nparam, int m, int M,
		    int print_all, unsigned MaxRays)
{
    int i;
    int res;
    Value *p, tmp;
    Polyhedron *S;
    Polyhedron *C = Polyhedron_Project(P, nparam);
    value_init(min_val);
    value_init(max_val);
    value_set_si(min_val,m);
    value_set_si(max_val,M);
    value_init(tmp);

    p = (Value *)malloc(sizeof(Value) * (P->Dimension+2));
    for(i=0;i<=P->Dimension;i++) {
      value_init(p[i]);
      value_set_si(p[i],0);
    }
    value_init(p[i]);
    value_set_si(p[i],1);

    /* S = scanning list of polyhedra */
    S = Polyhedron_Scan(P, C, MaxRays & POL_NO_DUAL ? 0 : MaxRays);

    if (!print_all) {
	if (C->Dimension > 0) {
	  value_subtract(tmp,max_val,min_val);
	  if (VALUE_TO_INT(tmp) > 80)
	    st = 1+(VALUE_TO_INT(tmp))/80;
	  else
	    st=1;
	  for(i=VALUE_TO_INT(min_val);i<=VALUE_TO_INT(max_val);i+=st)
	    printf(".");
	  printf( "\r" );
	  fflush(stdout);
	}
    }

    /******* CHECK NOW *********/
    res = 0;
    if(S && !check_poly(S, C, EP, exist, nparam, 0, p, print_all)) {
      fprintf(stderr,"Check failed !\n");
      res = -1;
    }
      
    if (!print_all)
	printf( "\n" );
    
    for(i=0;i<=(P->Dimension+1);i++) 
      value_clear(p[i]);
    free(p);
    value_clear(tmp);
    Domain_Free(S);
    Polyhedron_Free(C);
}

/****************************************************/
/* function check_poly :                            */
/* scans the parameter space from min to max (all   */
/* directions). Computes the number of points in    */
/* the polytope using both methods, and compare them*/
/* returns 1 on success                             */
/****************************************************/

int check_poly(Polyhedron *S, Polyhedron *C, evalue *EP,
	       int exist, int nparam, int pos, Value *z, int print_all)
{  
  int k;
  Value c,tmp;
  
  value_init(c); value_init(tmp);
  
  if(pos == nparam) {
    
    /* Computes the ehrhart polynomial */
    value_set_double(c, compute_evalue(EP,&z[S->Dimension-nparam+1])+.25);
    /* if c=0 we may be out of context. */
    /* scanning is useless in this case*/
    if(!in_domain(C,&z[S->Dimension-nparam+1])) {
   
      /* ok */ ;
    }
    else {
      
    if (print_all) {
      printf("EP( ");
      value_print(stdout,VALUE_FMT,z[S->Dimension-nparam+1]);
      for(k=S->Dimension-nparam+2;k<=S->Dimension;++k) {
	printf(", ");
	value_print(stdout,VALUE_FMT,z[k]);
      }
      printf(" ) = ");
      value_print(stdout,VALUE_FMT,c);
      printf(" ");
    }

      /* Manually count the number of points */
      count_points_e(1, S, exist, nparam, z, &tmp);
    if (print_all) {
	printf(", count = ");
	value_print(stdout, P_VALUE_FMT, tmp);
	printf(". ");
    }

      if(value_ne(tmp,c)) {
        printf("\n"); 
        fflush(stdout);
        fprintf(stderr,"Error !\n");
        fprintf(stderr,"EP( ");
        value_print(stderr,VALUE_FMT,z[S->Dimension-nparam+1]);
        for(k=S->Dimension-nparam+2;k<=S->Dimension;++k) {
          fprintf(stderr,", ");
          value_print(stderr,VALUE_FMT,z[k]);
        }
        fprintf(stderr," ) should be ");
        value_print(stderr,VALUE_FMT,tmp);
        fprintf(stderr,", while EP eval gives ");
        value_print(stderr,VALUE_FMT,c);
        fprintf(stderr,".\n");
	print_evalue(stderr, EP, params);
#ifndef DONT_BREAK_ON_ERROR
	value_clear(c); value_clear(tmp);
	return(0);
#endif
      }
      else if (print_all)
	printf("OK.\n");
    }
  }
  else
    for(value_assign(tmp,min_val); value_le(tmp,max_val); value_increment(tmp,tmp)) {
      if (!print_all) {
	  k = VALUE_TO_INT(tmp);
	  if(!pos && !(k%st)) {
	    printf("o");
	    fflush(stdout);
	  }
       }
      
      value_assign(z[pos+S->Dimension-nparam+1],tmp);
      if(!check_poly(S, C, EP, exist, nparam, pos+1, z, print_all)) {
	value_clear(c); value_clear(tmp);
	return(0);
      }
    }
  value_clear(c); value_clear(tmp);
  return(1);
} /* check_poly */
