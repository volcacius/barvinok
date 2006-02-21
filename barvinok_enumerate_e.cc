#include <unistd.h>
#include <stdlib.h>
#include <polylib/polylibgmp.h>
#include <barvinok/util.h>
#include <barvinok/barvinok.h>
#include "config.h"
#ifdef HAVE_OMEGA
#include "omega/convert.h"
#endif

/* The input of this example program is a polytope in combined
 * data and parameter space followed by two lines indicating
 * the number of existential variables and parameters respectively.
 * The first lines starts with "E ", followed by a number.
 * The second lines starts with "P ", followed by a number.
 * These two lines are (optionally) followed by the names of the parameters.
 * The polytope is in PolyLib notation.
 */

#ifdef HAVE_GROWING_CHERNIKOVA
#define MAXRAYS    POL_NO_DUAL
#else
#define MAXRAYS  600
#endif

#ifndef HAVE_GETOPT_H
#define getopt_long(a,b,c,d,e) getopt(a,b,c)
#else
#include <getopt.h>
struct option options[] = {
#ifdef HAVE_OMEGA
    { "omega",      no_argument,  0,  'o' },
#endif
#ifdef HAVE_PIPLIB
    { "pip",   	    no_argument,  0,  'p' },
#endif
    { "convert",    no_argument,  0,  'c' },
    { "floor",      no_argument,  0,  'f' },
    { "range-reduction",	no_argument,	0,  'R' },
    { "verify",     no_argument,  0,  'T' },
    { "min",   	    required_argument,  0,  'm' },
    { "max",   	    required_argument,  0,  'M' },
    { "range",      required_argument,  0,  'r' },
    { "version",    no_argument,  0,  'V' },
    { 0, 0, 0, 0 }
};
#endif

#ifdef HAVE_PIPLIB
#define PIPLIB_OPT "p"
#else
#define PIPLIB_OPT ""
#endif

#ifdef HAVE_OMEGA
#define OMEGA_OPT "o"

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
#define OMEGA_OPT ""
Polyhedron *Omega_simplify(Polyhedron *P, 
			    unsigned exist, unsigned nparam, char **parms)
{
    return P;
}
#endif

/* define this to print all the results */
/* else, only a progress bar is printed */
/* #define PRINT_ALL_RESULTS	 
/* define this to continue the test after first error found */
/* #define DONT_BREAK_ON_ERROR */

/* RANGE : normal range for evalutations (-RANGE -> RANGE) */
#define RANGE 50

/* SRANGE : small range for evalutations */
#define SRANGE 15

/* if dimension >= BIDDIM, use SRANGE */
#define BIGDIM 5

/* VSRANGE : very small range for evalutations */
#define VSRANGE 5

/* if dimension >= VBIDDIM, use VSRANGE */
#define VBIGDIM 8

static Value min_val, max_val;

static char **params;

#ifdef DONT_BREAK_ON_ERROR
#define PRINT_ALL_RESULTS
#endif

#ifndef PRINT_ALL_RESULTS
static int st;
#endif

static int check_poly(Polyhedron *S, Polyhedron *C, evalue *EP,
		      int exist, int nparam, int pos, Value *z);
static void verify_results(Polyhedron *P, evalue *EP, int exist, int nparam, 
			   int m, int M);

int main(int argc, char **argv)
{
    Polyhedron *A;
    Matrix *MA;
    char **param_name;
    int exist, nparam, nvar;
    char s[128];
    evalue *EP;
    int c, ind = 0;
    int range = 0;
    int convert = 0;
    int omega = 0;
    int pip = 0;
    int floor = 0;
    int verify = 0;
    int m = INT_MAX, M = INT_MIN, r;
    int print_solution = 1;

    while ((c = getopt_long(argc, argv, 
		    OMEGA_OPT PIPLIB_OPT "fcRTm:M:r:V", options, &ind)) != -1) {
	switch (c) {
	case 'o':
	    omega = 1;
	    break;
	case 'p':
	    pip = 1;
	    break;
	case 'f':
	    floor = 1;
	    break;
	case 'c':
	    convert = 1;
	    break;
	case 'R':
	    range = 1;
	    break;
	case 'T':
	    verify = 1;
	    break;
	case 'm':
	    m = atoi(optarg);
	    verify = 1;
	    break;
	case 'M':
	    M = atoi(optarg);
	    verify = 1;
	    break;
	case 'r':
	    M = atoi(optarg);
	    m = -M;
	    verify = 1;
	    break;
	case 'V':
	    printf(barvinok_version());
	    exit(0);
	    break;
	}
    }

    MA = Matrix_Read();
    A = Constraints2Polyhedron(MA, MAXRAYS);
    Matrix_Free(MA);

    fgets(s, 128, stdin);
    while ((*s=='#') || (sscanf(s, "E %d", &exist)<1))
	fgets(s, 128, stdin);

    fgets(s, 128, stdin);
    while ((*s=='#') || (sscanf(s, "P %d", &nparam)<1))
	fgets(s, 128, stdin);

    /******* Read the options: initialize Min and Max ********/
    if (A->Dimension >= VBIGDIM)
	r = VSRANGE;
    else if (A->Dimension >= BIGDIM)
	r = SRANGE;
    else
	r = RANGE;
    if (M == INT_MIN)
	M = r;
    if (m == INT_MAX)
	m = -r;

    if (verify && m > M) {
	fprintf(stderr,"Nothing to do: min > max !\n");
	return(0);
    }
    if (verify)
	print_solution = 0;

    if (print_solution) {
	Polyhedron_Print(stdout, P_VALUE_FMT, A);
	printf("exist: %d, nparam: %d\n", exist, nparam);
    }
    param_name = Read_ParamNames(stdin, nparam);
    nvar = A->Dimension - exist - nparam;
    if (omega) {
	A = Omega_simplify(A, exist, nparam, param_name);
	assert(!A->next);
	exist = A->Dimension - nvar - nparam;
    }
    if (pip && exist > 0)
	EP = barvinok_enumerate_pip(A, exist, nparam, MAXRAYS);
    else
	EP = barvinok_enumerate_e(A, exist, nparam, MAXRAYS);
    reduce_evalue(EP);
    evalue_combine(EP);
    if (range)
	evalue_range_reduction(EP);
    if (print_solution)
	print_evalue(stdout, EP, param_name);
    if (floor) {
	fprintf(stderr, "WARNING: floor conversion not supported\n");
	evalue_frac2floor(EP);
	if (print_solution)
	    print_evalue(stdout, EP, param_name);
    } else if (convert) {
	evalue_mod2table(EP, nparam);
	if (print_solution)
	    print_evalue(stdout, EP, param_name);
    }
    if (verify)
	verify_results(A, EP, exist, nparam, m, M);
    free_evalue_refs(EP);
    free(EP);
    Free_ParamNames(param_name, nparam);
    Polyhedron_Free(A);
    return 0;
}

void verify_results(Polyhedron *P, evalue *EP, int exist, int nparam, int m, int M)
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
    S = Polyhedron_Scan(P, C, MAXRAYS & POL_NO_DUAL ? 0 : MAXRAYS);

#ifndef PRINT_ALL_RESULTS
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
#endif

    /******* CHECK NOW *********/
    res = 0;
    if(S && !check_poly(S, C, EP, exist, nparam, 0, p)) {
      fprintf(stderr,"Check failed !\n");
      res = -1;
    }
      
#ifndef PRINT_ALL_RESULTS
    printf( "\n" );
#endif
    
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
	       int exist, int nparam, int pos, Value *z)
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
      
#ifdef PRINT_ALL_RESULTS
      printf("EP( ");
      value_print(stdout,VALUE_FMT,z[S->Dimension-nparam+1]);
      for(k=S->Dimension-nparam+2;k<=S->Dimension;++k) {
	printf(", ");
	value_print(stdout,VALUE_FMT,z[k]);
      }
      printf(" ) = ");
      value_print(stdout,VALUE_FMT,c);
      printf(" ");
#endif

      /* Manually count the number of points */
      count_points_e(1, S, exist, nparam, z, &tmp);
#ifdef PRINT_ALL_RESULTS
	printf(", count = ");
	value_print(stdout, P_VALUE_FMT, tmp);
	printf(". ");
#endif

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

#ifdef PRINT_ALL_RESULTS
      else
	printf("OK.\n");
#endif
    }
  }
  else
    for(value_assign(tmp,min_val); value_le(tmp,max_val); value_increment(tmp,tmp)) {

#ifndef PRINT_ALL_RESULTS
      k = VALUE_TO_INT(tmp);
      if(!pos && !(k%st)) {
	printf("o");
	fflush(stdout);
      }
#endif
      
      value_assign(z[pos+S->Dimension-nparam+1],tmp);
      if(!check_poly(S, C, EP, exist, nparam, pos+1, z)) {
	value_clear(c); value_clear(tmp);
	return(0);
      }
    }
  value_clear(c); value_clear(tmp);
  return(1);
} /* check_poly */
