#ifndef BARVINOK_BASIS_REDUCTION_H
#define BARVINOK_BASIS_REDUCTION_H

#include <gmp.h>

#if defined(__cplusplus)
extern "C" {
#endif

#include <polylib/polylibgmp.h>

Matrix *Polyhedron_Reduced_Basis(Polyhedron *P);

#if defined(__cplusplus)
}
#endif

#endif
