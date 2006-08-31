#include <ginac/ginac.h>
#include <bernstein/piecewise_lst.h>
#include <barvinok/evalue.h>

namespace barvinok {

bernstein::piecewise_lst_s *evalue_bernstein_coefficients(
	    bernstein::piecewise_lst_s *pl_all, evalue *e, 
	    Polyhedron *ctx, const GiNaC::exvector& params);
GiNaC::ex evalue2ex(evalue *e, const GiNaC::exvector& vars);

}
