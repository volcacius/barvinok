/*************************************************/
/*     verif_ehrhart.c                           */
/* program to compare effective number of points */
/* in a polytope with the corresponding          */
/* evaluation of the Ehrhart polynomial.         */
/* Parameters vary in range -RANGE to RANGE      */
/* (define below) by default.                    */
/* Can be overridden by specifying               */
/* -r<RANGE>, or -m<min> and -M<max>             */
/*                                               */
/* written by Vincent Loechner (c) 2000.         */
/*  loechner@icps.u-strasbg.fr                   */
/*************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include <polylib/polylibgmp.h>
#include <barvinok/evalue.h>
#include <barvinok/util.h>
#include <barvinok/barvinok.h>
#include "config.h"

#ifdef HAVE_GROWING_CHERNIKOVA
#define MAXRAYS     0
#else
#define MAXRAYS  1024
#endif

#ifndef HAVE_GETOPT_H
#define getopt_long(a,b,c,d,e) getopt(a,b,c)
#else
#include <getopt.h>
struct option options[] = {
    { "min",   no_argument,  0,  'm' },
    { "max",   no_argument,  0,  'M' },
    { "range",   no_argument,  0,  'r' },
    { "pip",   no_argument,  0,  'p' },
    { "version",   no_argument,  0,  'V' },
    { 0, 0, 0, 0 }
};
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

Value min, max;

char **params;

#ifdef DONT_BREAK_ON_ERROR
#define PRINT_ALL_RESULTS
#endif

#ifndef PRINT_ALL_RESULTS
int st;
#endif

/****************************************************/
/* function check_poly :                            */
/* scans the parameter space from min to max (all   */
/* directions). Computes the number of points in    */
/* the polytope using both methods, and compare them*/
/* returns 1 on success                             */
/****************************************************/

int check_poly(Polyhedron *S, Polyhedron *C, evalue *EP,
	       int exist, int nparam, int pos, Value *z) {
  
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
    for(value_assign(tmp,min); value_le(tmp,max); value_increment(tmp,tmp)) {

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

int main(int argc,char *argv[])
{
  Matrix *C1, *P1;
  Polyhedron *C, *P, *S;
  Value *p, tmp;
  int i,j;
    int m = INT_MAX, M = INT_MIN, r;
  int exist, nparam;
  char s[128];
  evalue *EP;
  int res;
  int c, ind = 0;
  int pip = 0;

    while ((c = getopt_long(argc, argv, "pm:M:r:V", options, &ind)) != -1) {
	switch (c) {
	case 'p':
	    pip = 1;
	    break;
	case 'm':
	    m = atoi(optarg);
	    break;
	case 'M':
	    M = atoi(optarg);
	    break;
	case 'r':
	    M = atoi(optarg);
	    m = -M;
	    break;
	case 'V':
	    printf(barvinok_version());
	    exit(0);
	    break;
	}
    }
  
/******* Read the input *********/
  P1 = Matrix_Read();

    fgets(s, 128, stdin);
    while ((*s=='#') || (sscanf(s, "E %d", &exist)<1))
	fgets(s, 128, stdin);

    fgets(s, 128, stdin);
    while ((*s=='#') || (sscanf(s, "P %d", &nparam)<1))
	fgets(s, 128, stdin);

  P = Constraints2Polyhedron(P1,MAXRAYS);
  params = Read_ParamNames(stdin, nparam);
  Matrix_Free(P1);

  /******* Read the options: initialize Min and Max ********/
  if(P->Dimension >= VBIGDIM)
    r = VSRANGE;
  else if(P->Dimension >= BIGDIM)
    r = SRANGE;
  else
    r = RANGE;
  if (M == INT_MIN)
    M = r;
  if (m == INT_MAX)
    m = -r;

  if(m > M) {
    fprintf(stderr,"Nothing to do: min > max !\n");
    return(0);
  }
  value_init(min);
  value_init(max);
  value_set_si(min,m);
  value_set_si(max,M);
  value_init(tmp);

  /******* Compute true context *******/
  C1 = Matrix_Alloc(nparam+1,P->Dimension+1);

  for(i=0;i<C1->NbRows;i++)
    for(j=0;j<C1->NbColumns;j++)
      if(i==j-P->Dimension+nparam)
	value_set_si(C1->p[i][j],1);
      else
	value_set_si(C1->p[i][j],0);
  C = Polyhedron_Image(P,C1,MAXRAYS);
  Matrix_Free(C1);

  /******* Compute EP *********/
  if (pip)
    EP = barvinok_enumerate_pip(P, exist, nparam, MAXRAYS);
  else
    EP = barvinok_enumerate_e(P, exist, nparam, MAXRAYS);
  
  /******* Initializations for check *********/
  p = (Value *)malloc(sizeof(Value) * (P->Dimension+2));
  for(i=0;i<=P->Dimension;i++) {
    value_init(p[i]);
    value_set_si(p[i],0);
  }
  value_init(p[i]);
  value_set_si(p[i],1);

  /* S = scanning list of polyhedra */
  S = Polyhedron_Scan(P,C,MAXRAYS);

#ifndef PRINT_ALL_RESULTS
  if(C->Dimension > 0) {
    value_subtract(tmp,max,min);
    if (VALUE_TO_INT(tmp) > 80)
      st = 1+(VALUE_TO_INT(tmp))/80;
    else
      st=1;
    for(i=VALUE_TO_INT(min);i<=VALUE_TO_INT(max);i+=st)
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
  Free_ParamNames(params, C->Dimension);
  Polyhedron_Free(S);
  Polyhedron_Free(C);
  Polyhedron_Free(P);
  free_evalue_refs(EP);
  free(EP);
  return res;
} /* main */




