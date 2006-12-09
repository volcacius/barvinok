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
#include <barvinok/options.h>

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
    short_rat(Value c);
    short_rat(const QQ& c, const vec_ZZ& num, const mat_ZZ& den);
};

struct gen_fun {
    std::vector< short_rat * > term;
    Polyhedron *context;

    void add(const QQ& c, const vec_ZZ& num, const mat_ZZ& den);
    /* add c times gf */
    void add(const QQ& c, const gen_fun *gf);
    void substitute(Matrix *CP);
    gen_fun *Hadamard_product(const gen_fun *gf, barvinok_options *options);
    void add_union(gen_fun *gf, barvinok_options *options);
    void shift(const vec_ZZ& offset);
    void divide(const vec_ZZ& power);
    void print(std::ostream& os, unsigned int nparam, char **param_name) const;
    operator evalue *() const;
    void coefficient(Value* params, Value* c) const;
    gen_fun *summate(int nvar, barvinok_options *options) const;
    bool summate(Value *sum) const;

    gen_fun(const gen_fun *gf) {
	QQ one(1, 1);
	context = Polyhedron_Copy(gf->context);
	add(one, gf);
    }
    gen_fun(Value c);
    gen_fun(Polyhedron *C = NULL) : context(C) {}
    ~gen_fun() {
	if (context)
	    Polyhedron_Free(context);
	for (int i = 0; i < term.size(); ++i)
	    delete term[i];
    }
};

#endif
