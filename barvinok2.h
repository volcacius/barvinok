#ifndef BARVINOK2_H
#define BARVINOK2_H

#include <gmp.h>
extern "C" {
#include <polylib/polylibgmp.h>
}
#include <genfun.h>

gen_fun * barvinok_series(Polyhedron *P, Polyhedron* C, unsigned MaxRays);

#endif
