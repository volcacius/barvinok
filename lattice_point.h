#ifndef LATTICE_POINT_H
#define LATTICE_POINT_H

#include <gmp.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <barvinok/polylib.h>

#ifdef NTL_STD_CXX
using namespace NTL;
#endif

struct barvinok_options;

evalue *multi_monom(vec_ZZ& p);
int normal_mod(Value *coef, int len, Value *m);
evalue *fractional_part(Value *coef, Value denom, int nvar, Polyhedron *PD);
void lattice_point_fixed(Value *vertex, Value *vertex_res,
			 Matrix *Rays, Matrix *Rays_res,
			 Value *point);
void lattice_points_fixed(Value *vertex, Value *vertex_res,
			  Matrix *Rays, Matrix *Rays_res, Matrix *points,
			  unsigned long det);
void lattice_point(Param_Vertices *V, const mat_ZZ& rays, vec_ZZ& num, 
		   evalue **E_vertex, barvinok_options *options);

/* This structure encodes the power of the term in a rational generating function.
 * 
 * Either E == NULL or constant = 0
 * If E != NULL, then the power is 	    E
 * If E == NULL, then the power is 	    constant
 */
struct term_info {
    evalue	  **E;
    vec_ZZ	    constant;
};

void lattice_point(Param_Vertices* V, const mat_ZZ& rays, vec_ZZ& lambda,
    term_info* term, unsigned long det,
    barvinok_options *options);

#endif
