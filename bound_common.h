#include <barvinok/bernstein.h>

#if defined(__cplusplus)
extern "C" {
#endif

__isl_give isl_pw_qpolynomial_fold *isl_pw_qpolynomial_bound_bernstein(
	__isl_take isl_pw_qpolynomial *pwqp, enum isl_fold type);
__isl_give isl_pw_qpolynomial_fold *isl_pw_qpolynomial_bound(
	__isl_take isl_pw_qpolynomial *pwqp, enum isl_fold type, int method);

#if defined(__cplusplus)
}
#endif
