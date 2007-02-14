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

/****************************************************/
/* function check_poly :                            */
/* scans the parameter space from Min to Max (all   */
/* directions). Computes the number of points in    */
/* the polytope using both methods, and compare them*/
/* returns 1 on success                             */
/****************************************************/

int check_poly(Polyhedron *S,Polyhedron *CS,Enumeration *en,
	       int nparam, int pos, Value *z, const struct verify_options *options)
{
  int k;
  Value c,tmp,*ctmp;
  Value LB, UB;
  
  value_init(c); value_init(tmp);
  value_init(LB);
  value_init(UB);
  
  if(pos == nparam) {
    
    /* Computes the ehrhart polynomial */
    value_assign(c,*(ctmp=compute_poly(en,&z[S->Dimension-nparam+1])));
    value_clear(*ctmp);
    free(ctmp);
    /* if c=0 we may be out of context. */
    /* scanning is useless in this case*/

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
      count_points(1,S,z,&tmp);
    if (options->print_all) {
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
        {
        	 Enumeration *ee;
        	 Enumeration_Print(stderr, en, options->params);
        	 ee = en;
        	 while (ee) {
		    if (in_domain(ee->ValidityDomain,&z[S->Dimension-nparam+1])) {
                	 Print_Domain(stderr, ee->ValidityDomain, options->params);
                	 print_evalue(stderr, &ee->EP, options->params);
		    }
		    ee = ee->next;
		 }
        }
#ifndef DONT_BREAK_ON_ERROR
	value_clear(c); value_clear(tmp);
	value_clear(LB);
	value_clear(UB);
	return(0);
#endif
      } else if (options->print_all)
	printf("OK.\n");
  }
  else {
    int ok = 
	!(lower_upper_bounds(1+pos, CS, &z[S->Dimension-nparam], &LB, &UB));
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
      if (!check_poly(S, CS->next, en, nparam, pos+1, z, options)) {
	value_clear(c); value_clear(tmp);
	value_clear(LB);
	value_clear(UB);
	return(0);
      }
    }
    value_set_si(z[pos+S->Dimension-nparam+1],0);
  }
  value_clear(c); value_clear(tmp);
  value_clear(LB);
  value_clear(UB);
  return(1);
} /* check_poly */

