#ifndef BERNSTEIN_MAXIMIZE_H
#define BERNSTEIN_MAXIMIZE_H

#include <gmp.h>
#include <ginac/ginac.h>
extern "C" {
#include <polylib/polylibgmp.h>
}

namespace bernstein {

GiNaC::lst maximize(Polyhedron *domain, GiNaC::lst coeffs,
		    const GiNaC::exvector& vars);
GiNaC::lst minimize(Polyhedron *domain, GiNaC::lst coeffs,
		    const GiNaC::exvector& vars);
GiNaC::lst remove_redundants(Polyhedron *domain, GiNaC::lst list1, GiNaC::lst list2,
			     const GiNaC::exvector& vars, int sign);

}

#endif
