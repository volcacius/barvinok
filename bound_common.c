#include <assert.h>
#include <barvinok/options.h>
#include <bound_common.h>

#include "config.h"

#ifndef HAVE_GINAC
__isl_give isl_pw_qpolynomial_fold *isl_pw_qpolynomial_bound_bernstein(
	__isl_take isl_pw_qpolynomial *pwqp, enum isl_fold type)
{
	assert(0);
}
#endif

__isl_give isl_pw_qpolynomial_fold *isl_pw_qpolynomial_upper_bound(
	__isl_take isl_pw_qpolynomial *pwqp)
{
#ifdef HAVE_GINAC
	return isl_pw_qpolynomial_bound(pwqp, isl_fold_max, BV_BOUND_BERNSTEIN);
#else
	return isl_pw_qpolynomial_bound(pwqp, isl_fold_max, BV_BOUND_RANGE);
#endif
}

__isl_give isl_pw_qpolynomial_fold *isl_pw_qpolynomial_bound(
	__isl_take isl_pw_qpolynomial *pwqp, enum isl_fold type, int method)
{
	if (method == BV_BOUND_BERNSTEIN)
		return isl_pw_qpolynomial_bound_bernstein(pwqp, type);
	else
		return isl_pw_qpolynomial_bound_range(pwqp, type, NULL);
}
