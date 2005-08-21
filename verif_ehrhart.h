#include <gmp.h>
#include <NTL/mat_ZZ.h>
extern "C" {
#include <polylib/polylibgmp.h>
#include <barvinok/evalue.h>

/* define this to print all the results */
/* else, only a progress bar is printed */
/* #define PRINT_ALL_RESULTS	 */
/* define this to continue the test after first error found */
/* #define DONT_BREAK_ON_ERROR */

extern Value Min, Max;

extern char **params;

#ifdef DONT_BREAK_ON_ERROR
#define PRINT_ALL_RESULTS
#endif

#ifndef PRINT_ALL_RESULTS
extern int st;
#endif

int check_poly(Polyhedron *S,Polyhedron *C,Enumeration *en,
	       int nparam,int pos,Value *z);
}
