#ifndef BARVINOK_H
#define BARVINOK_H

#include <gmp.h>

#if defined(__cplusplus)
extern "C" {
#endif

#include <polylib/polylibgmp.h>

Polyhedron *decompose(Polyhedron *C);

#if defined(__cplusplus)
}
#endif

#endif
