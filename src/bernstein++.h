/* 
 *	Bernstein Expansion
 *
 *	- c to c++ functions
 */

#include <string.h>
#include <ginac/ginac.h>
#include <gmp.h>
#include "polylib++.h"

GiNaC::ex convertPolynomial(long long *m, unsigned int nbRows, unsigned int nbColumns,
			    const GiNaC::exvector& params);

GiNaC::ex polyConvertParameters(long long *m, unsigned int nbRows, 
			  unsigned int nbColumns, 
			  long long **llPolynomialCoefficients, unsigned int *llRows, 
			  unsigned int *llColumns, const GiNaC::exvector& vars,
			  const GiNaC::exvector& params);
