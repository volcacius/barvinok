#ifndef DPOLY_H
#define DPOLY_H

#include <assert.h>
#include <vector>
#include <gmp.h>
#include <NTL/vec_ZZ.h>
extern "C" {
#include <polylib/polylibgmp.h>
}
#include "conversion.h"

#ifdef NTL_STD_CXX
using namespace NTL;
#endif

class dpoly {
public:
    vec_ZZ coeff;
    dpoly(int d, ZZ& degree, int offset = 0);
    void operator *= (dpoly& f);
    void div(dpoly& d, mpq_t count, ZZ& sign);
};

struct dpoly_r_term {
    int	    *powers;
    ZZ	    coeff;
};

/* len: number of elements in c
 * each element in c is the coefficient of a power of t
 * in the MacLaurin expansion
 */
struct dpoly_r {
    std::vector< dpoly_r_term * >	*c;
    int len;
    int dim;
    ZZ denom;

    void add_term(int i, int * powers, ZZ& coeff);
    dpoly_r(int len, int dim);
    dpoly_r(dpoly& num, int dim);
    dpoly_r(dpoly& num, dpoly& den, int pos, int dim);
    dpoly_r(dpoly_r* num, dpoly& den, int pos, int dim);
    ~dpoly_r();
    dpoly_r *div(dpoly& d);
    void dump(void);
};

#endif
