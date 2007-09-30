#ifndef BARVINOK_POLYLIB_H
#define BARVINOK_POLYLIB_H

#include <gmp.h>

#if defined(__cplusplus)
extern "C" {
#endif

#include <polylib/polylibgmp.h>

#define	value_gcd(ref,val1,val2)	mpz_gcd(ref,val1,val2)
#define	value_lcm(ref,val1,val2)	mpz_lcm(ref,val1,val2)

#undef divide
#undef value_compare

#if defined(__cplusplus)
}
#endif

#endif
