#ifndef GENFUN_H
#define GENFUN_H

#include <vector>
#include <gmp.h>
#include <NTL/mat_ZZ.h>
extern "C" {
#include <polylib/polylibgmp.h>
#include "ev_operations.h"
}

#ifdef NTL_STD_CXX
using namespace NTL;
#endif

struct short_rat {
    struct {
	/* rows: terms in numerator */
	mat_ZZ	coeff;
	mat_ZZ	power;
    } n;
    struct {
	/* rows: factors in denominator */
	mat_ZZ	power;
    } d;
};

struct gen_fun {
    std::vector< short_rat * > term;

    void add(ZZ& cn, ZZ& cd, vec_ZZ& num, mat_ZZ& den);
    void print(unsigned int nparam, char **param_name);
    operator evalue *();
};

#endif
