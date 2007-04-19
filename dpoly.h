#ifndef DPOLY_H
#define DPOLY_H

#include <assert.h>
#include <barvinok/set.h>
#include <vector>
#include <gmp.h>
#include <NTL/vec_ZZ.h>
#include <barvinok/polylib.h>
#include "conversion.h"

#ifdef NTL_STD_CXX
using namespace NTL;
#endif

class dpoly {
public:
    vec_ZZ coeff;
    dpoly(int d, ZZ& degree, int offset = 0);
    void operator += (const dpoly& t);
    void operator *= (const ZZ& f);
    void operator *= (dpoly& f);
    void div(dpoly& d, mpq_t count, ZZ& sign);
    void div(dpoly& d, mpq_t *count, const mpq_t& factor);
private:
    mpq_t *div(dpoly &d) const;
    void   clear_div(mpq_t *c) const;
};

/* Each element in powers corresponds to a factor of the form (1 - z^b)
 * and indicates the exponent of this factor in the denominator.
 * The constants b are stored elsewhere (den_r in reducer::reducer).
 */
struct dpoly_r_term {
    std::vector<int>    powers;
    ZZ	    	    	coeff;
};

struct dpoly_r_term_lex_smaller {
    bool operator()(const dpoly_r_term* t1, const dpoly_r_term* t2) const {
	return t1->powers < t2->powers;
    }
};

typedef std::set<dpoly_r_term*, dpoly_r_term_lex_smaller> dpoly_r_term_list;

/* len: number of elements in c
 * each element in c is the coefficient of a power of t
 * in the MacLaurin expansion
 */
struct dpoly_r {
    dpoly_r_term_list	*c;
    int len;
    int dim;
    ZZ denom;

    void add_term(int i, const std::vector<int>& powers, const ZZ& coeff);
    dpoly_r(int len, int dim);
    dpoly_r(dpoly& num, int dim);
    dpoly_r(dpoly& num, dpoly& den, int pos, int dim);
    dpoly_r(const dpoly_r* num, dpoly& den, int pos, int dim);
    ~dpoly_r();
    dpoly_r *div(const dpoly& d) const;
    void dump(void);
};

#endif
