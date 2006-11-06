/* Don't make this header public without removing the include of config.h */
#if defined(__cplusplus)
extern "C" {
#endif

#include <polylib/polylibgmp.h>
#include "config.h"

#ifndef HAVE_COMPRESS_PARMS
static int remove_all_equalities(Polyhedron **P, Polyhedron **C, Matrix **CPP, 
				 Matrix **CVP, unsigned nparam, unsigned MaxRays)
{
    assert(0);
}
#else
int remove_all_equalities(Polyhedron **P, Polyhedron **C, Matrix **CPP, Matrix **CVP,
			  unsigned nparam, unsigned MaxRays);
#endif

#if defined(__cplusplus)
}
#endif
