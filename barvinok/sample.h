#ifndef BARVINOK_SAMPLE_H
#define BARVINOK_SAMPLE_H

#include <gmp.h>

#if defined(__cplusplus)
extern "C" {
#endif

#include <polylib/polylibgmp.h>

Vector *Polyhedron_Sample(Polyhedron *P, unsigned MaxRays);

#if defined(__cplusplus)
}
#endif

#endif
