#ifndef SCARF_H
#define SCARF_H

#include <barvinok/genfun.h>

evalue *barvinok_enumerate_scarf(Polyhedron *P,
			  unsigned exist, unsigned nparam, unsigned MaxRays);
gen_fun *barvinok_enumerate_scarf_series(Polyhedron *P,
			  unsigned exist, unsigned nparam, unsigned MaxRays);

#endif
