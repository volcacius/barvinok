#ifndef GENFUN_H
#define GENFUN_H

#include <vector>
#include <gmp.h>
#include <NTL/mat_ZZ.h>
extern "C" {
#include <polylib/polylibgmp.h>
}

#ifdef NTL_STD_CXX
using namespace NTL;
#endif

struct short_rat {
    struct {
	mat_ZZ	coeff;
	mat_ZZ	power;
    } n;
    struct {
	mat_ZZ	power;
    } d;
};

struct gen_fun {
    std::vector< short_rat * > term;

    void add(ZZ& cn, ZZ& cd, vec_ZZ& num, mat_ZZ& den);
    void print(void);
};

#endif
