#ifndef BARVINOK_H
#define BARVINOK_H

#include <polylib/polylibgmp.h>

Polyhedron* Polyhedron_Polar(Polyhedron *P, unsigned NbMaxRays);
Polyhedron* supporting_cone(Polyhedron *P, int v, unsigned NbMaxRays);

#endif
