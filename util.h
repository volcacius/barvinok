#ifndef UTIL_H
#define UTIL_H

#if defined(__cplusplus)
extern "C" {
#endif

#include <polylib/polylibgmp.h>

Polyhedron* Polyhedron_Polar(Polyhedron *P, unsigned NbMaxRays);
Polyhedron* supporting_cone(Polyhedron *P, int v, unsigned NbMaxRays);
Polyhedron* triangularize_cone(Polyhedron *P, unsigned NbMaxCons);
Polyhedron *remove_equalities(Polyhedron *P);
void manual_count(Polyhedron *P, Value* result);
Polyhedron* reduce(Polyhedron *P, Value* factor);

#if defined(__cplusplus)
}
#endif

#endif
