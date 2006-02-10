#ifndef GENFUN_H
#define GENFUN_H

#include <vector>
#include <gmp.h>
#include <NTL/mat_ZZ.h>
extern "C" {
#include <polylib/polylibgmp.h>
#include <barvinok/evalue.h>
}

#ifdef NTL_STD_CXX
using namespace NTL;
#endif

struct short_rat {
    struct __short_rat_n {
	/* rows: terms in numerator */
	/* coeff has two columns: the numerator and the denominator */
	mat_ZZ	coeff;
	mat_ZZ	power;
    } n;
    struct __short_rat_d {
	/* rows: factors in denominator */
	mat_ZZ	power;
    } d;
};

struct gen_fun {
    std::vector< short_rat * > term;
    Polyhedron *context;

    void add(const ZZ& cn, const ZZ& cd, const vec_ZZ& num, 
	     const mat_ZZ& den);
    void substitute(Matrix *CP, const mat_ZZ& map, const vec_ZZ& offset);
    void print(unsigned int nparam, char **param_name) const;
    operator evalue *() const;
    void coefficient(Value* params, Value* c) const;

    gen_fun(Polyhedron *C = NULL) : context(C) {}
    ~gen_fun() {
	if (context)
	    Polyhedron_Free(context);
	for (int i = 0; i < term.size(); ++i)
	    delete term[i];
    }
};

#endif
