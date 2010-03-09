#include <barvinok/evalue.h>

#if defined(__cplusplus)
extern "C" {
#endif

__isl_give isl_pw_qpolynomial_fold *isl_pw_qpolynomial_upper_bound(
	__isl_take isl_pw_qpolynomial *pwqp);

#if defined(__cplusplus)
}

#include <ginac/ginac.h>
#include <bernstein/piecewise_lst.h>

struct barvinok_options;

namespace barvinok {

bernstein::piecewise_lst *evalue_bernstein_coefficients(
	    bernstein::piecewise_lst *pl_all, evalue *e, 
	    Polyhedron *ctx, const GiNaC::exvector& params);
bernstein::piecewise_lst *evalue_bernstein_coefficients(
	    bernstein::piecewise_lst *pl_all, evalue *e, 
	    Polyhedron *ctx, const GiNaC::exvector& params,
	    barvinok_options *options);
GiNaC::ex evalue2ex(evalue *e, const GiNaC::exvector& vars);
GiNaC::ex evalue2ex(const evalue *e, const GiNaC::exvector& vars,
		     GiNaC::exvector& floorvar, Matrix **C);

__isl_give isl_qpolynomial *isl_qpolynomial_from_ginac(__isl_take isl_dim *dim,
	const GiNaC::ex &ex, const GiNaC::exvector &params);
__isl_give isl_qpolynomial_fold *isl_qpolynomial_fold_from_ginac(
	__isl_take isl_dim *dim, enum isl_fold type, const GiNaC::lst &lst,
	const GiNaC::exvector &params);
__isl_give isl_pw_qpolynomial_fold *isl_pw_qpolynomial_fold_from_ginac(
	__isl_take isl_dim *dim, bernstein::piecewise_lst *pl,
	const GiNaC::exvector &params);

}

#endif
