#ifndef UTIL_H
#define UTIL_H

#include <polylib/polylibgmp.h>

Polyhedron* Polyhedron_Polar(Polyhedron *P, unsigned NbMaxRays);
Polyhedron* supporting_cone(Polyhedron *P, int v, unsigned NbMaxRays);

#endif
