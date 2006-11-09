#ifndef BARVINOK_BASIS_REDUCTION_H
#define BARVINOK_BASIS_REDUCTION_H

#include <gmp.h>

#if defined(__cplusplus)
extern "C" {
#endif

#include <polylib/polylibgmp.h>

struct barvinok_options;

Matrix *Polyhedron_Reduced_Basis(Polyhedron *P, struct barvinok_options *options);

Matrix *glpk_Polyhedron_Reduced_Basis(Polyhedron *P);
Matrix *cdd_Polyhedron_Reduced_Basis(Polyhedron *P);

#if defined(__cplusplus)
}
#endif

#endif
