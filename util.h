#ifndef UTIL_H
#define UTIL_H

#if defined(__cplusplus)
extern "C" {
#endif

#include <polylib/polylibgmp.h>

int random_int(int max);
Polyhedron* Polyhedron_Polar(Polyhedron *P, unsigned NbMaxRays);
Polyhedron* supporting_cone(Polyhedron *P, int v, unsigned NbMaxRays);
Polyhedron* triangularize_cone(Polyhedron *P, unsigned NbMaxCons);
void check_triangulization(Polyhedron *P, Polyhedron *T);
Polyhedron *remove_equalities(Polyhedron *P);
void manual_count(Polyhedron *P, Value* result);
Polyhedron* Polyhedron_Reduce(Polyhedron *P, Value* factor);

#if defined(__cplusplus)
}
#endif

#endif
