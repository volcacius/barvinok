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

#include <barvinok/evalue.h>
#include <barvinok/barvinok.h>
#include "verif_ehrhart.h"

#undef CS   /* for Solaris 10 */

#include "config.h"
#ifndef HAVE_COUNT_POINTS4
#define count_points(a,b,c,d) {				\
	    int cc = count_points(a,b,c);		\
	    value_set_si(*d,cc);			\
	}
#endif

Polyhedron *check_poly_context_scan(Polyhedron *C, struct verify_options *options)
{
    int i;
    Matrix *MM;
    Polyhedron *CC, *CS, *U;

    if (C->Dimension <= 0)
	return NULL;

    /* Intersect context with range */
    MM = Matrix_Alloc(2*C->Dimension, C->Dimension+2);
    for (i = 0; i < C->Dimension; ++i) {
	value_set_si(MM->p[2*i][0], 1);
	value_set_si(MM->p[2*i][1+i], 1);
	value_set_si(MM->p[2*i][1+C->Dimension], -options->m);
	value_set_si(MM->p[2*i+1][0], 1);
	value_set_si(MM->p[2*i+1][1+i], -1);
	value_set_si(MM->p[2*i+1][1+C->Dimension], options->M);
    }
    CC = AddConstraints(MM->p[0], 2*C->Dimension, C, options->barvinok->MaxRays);
    U = Universe_Polyhedron(0);
    CS = Polyhedron_Scan(CC, U, options->barvinok->MaxRays);
    Polyhedron_Free(U);
    Polyhedron_Free(CC);
    Matrix_Free(MM);
    return CS;
}

void check_poly_init(Polyhedron *C, struct verify_options *options)
{
    int d, i;

    if (options->print_all)
	return;
    if (C->Dimension <= 0)
	return;

    d = options->M - options->m;
    if (d > 80)
	options->st = 1+d/80;
    else
	options->st = 1;
    for (i = options->m; i <= options->M; i += options->st)
	printf(".");
    printf( "\r" );
    fflush(stdout);
}

/****************************************************/
/* function check_poly :                            */
/* scans the parameter space from Min to Max (all   */
/* directions). Computes the number of points in    */
/* the polytope using both methods, and compare them*/
/* returns 1 on success                             */
/****************************************************/

int check_poly(Polyhedron *S, Polyhedron *CS, evalue *EP, int exist,
	       int nparam, int pos, Value *z, const struct verify_options *options)
{
  int k;
  Value c, tmp;
  int ok;
  int pa = options->barvinok->polynomial_approximation;
  
  value_init(c); value_init(tmp);
  
  if(pos == nparam) {
    
    /* Computes the ehrhart polynomial */
    value_set_double(c, compute_evalue(EP,&z[S->Dimension-nparam+1])+.25);

    if (options->print_all) {
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
	if (exist)
	    count_points_e(1, S, exist, nparam, z, &tmp);
	else
	    count_points(1, S, z, &tmp);

    if (options->print_all) {
	printf(", count = ");
	value_print(stdout, P_VALUE_FMT, tmp);
	printf(". ");
    }

	if (pa == BV_POLAPPROX_PRE_APPROX)
	    /* just accept everything */
	    ok = 1;
	else if (pa == BV_POLAPPROX_PRE_LOWER || pa == BV_POLAPPROX_LOWER)
	    ok = value_le(c, tmp);
	else if (pa == BV_POLAPPROX_PRE_UPPER || pa == BV_POLAPPROX_UPPER)
	    ok = value_ge(c, tmp);
	else
	    ok = value_eq(c, tmp);

      if (!ok) {
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
	    print_evalue(stderr, EP, options->params);
	    if (value_zero_p(EP->d) && EP->x.p->type == partition)
		for (k = 0; k < EP->x.p->size/2; ++k) {
		    Polyhedron *D = EVALUE_DOMAIN(EP->x.p->arr[2*k]);
		    if (in_domain(D, &z[S->Dimension-nparam+1])) {
			Print_Domain(stderr, D, options->params);
			print_evalue(stderr, &EP->x.p->arr[2*k+1], options->params);
		    }
	    }
	if (!options->continue_on_error) {
	    value_clear(c); value_clear(tmp);
	    return 0;
	}
      } else if (options->print_all)
	printf("OK.\n");
  }
  else {
	Value LB, UB;
	int ok;
	value_init(LB);
	value_init(UB);
	ok = !(lower_upper_bounds(1+pos, CS, &z[S->Dimension-nparam], &LB, &UB));
	assert(ok);
	for(value_assign(tmp,LB); value_le(tmp,UB); value_increment(tmp,tmp)) {
	    if (!options->print_all) {
		k = VALUE_TO_INT(tmp);
		if (!pos && !(k % options->st)) {
		    printf("o");
		    fflush(stdout);
		}
	    }
	      
	    value_assign(z[pos+S->Dimension-nparam+1],tmp);
	    if (!check_poly(S, CS->next, EP, exist, nparam, pos+1, z, options)) {
		value_clear(c); value_clear(tmp);
		value_clear(LB);
		value_clear(UB);
		return 0;
	    }
	}
	value_set_si(z[pos+S->Dimension-nparam+1],0);
	value_clear(LB);
	value_clear(UB);
  }
  value_clear(c); value_clear(tmp);
  return(1);
} /* check_poly */

