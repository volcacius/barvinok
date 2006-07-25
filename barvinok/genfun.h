#ifndef GENFUN_H
#define GENFUN_H

#include <vector>
#include <iostream>
#include <gmp.h>
#include <NTL/mat_ZZ.h>
extern "C" {
#include <polylib/polylibgmp.h>
#include <barvinok/evalue.h>
}
#include <barvinok/NTL_QQ.h>

#ifdef NTL_STD_CXX
using namespace NTL;
#endif

struct short_rat {
    struct __short_rat_n {
	/* rows: terms in numerator */
	/* coeff has two columns: the numerator and the denominator */
	vec_QQ	coeff;
	mat_ZZ	power;
    } n;
    struct __short_rat_d {
	/* rows: factors in denominator */
	mat_ZZ	power;
    } d;
    void add(short_rat *rat);
    bool reduced();
};

struct gen_fun {
    std::vector< short_rat * > term;
    Polyhedron *context;

    void add(const QQ& c, const vec_ZZ& num, const mat_ZZ& den);
    /* add cn/cd times gf */
    void add(const QQ& c, const gen_fun *gf);
    void substitute(Matrix *CP, const mat_ZZ& map, const vec_ZZ& offset);
    gen_fun *Hadamard_product(const gen_fun *gf, unsigned MaxRays);
    void add_union(gen_fun *gf, unsigned MaxRays);
    void shift(const vec_ZZ& offset);
    void print(std::ostream& os, unsigned int nparam, char **param_name) const;
    operator evalue *() const;
    void coefficient(Value* params, Value* c) const;
    gen_fun *summate(int nvar) const;

    gen_fun(const gen_fun *gf) {
	QQ one(1, 1);
	context = Polyhedron_Copy(gf->context);
	add(one, gf);
    }
    gen_fun(Polyhedron *C = NULL) : context(C) {}
    ~gen_fun() {
	if (context)
	    Polyhedron_Free(context);
	for (int i = 0; i < term.size(); ++i)
	    delete term[i];
    }
};

#endif
