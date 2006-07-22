#ifndef REDUCER_H
#define REDUCER_H

#include <NTL/mat_ZZ.h>
#include "decomposer.h"
#include "dpoly.h"

#ifdef NTL_STD_CXX
using namespace NTL;
#endif

/* base for non-parametric counting */
struct np_base : public polar_decomposer {
    unsigned dim;
    ZZ one;

    np_base(unsigned dim) {
	this->dim = dim;
	one = 1;
    }

    virtual void handle_polar(Polyhedron *C, Value *vertex, ZZ n, ZZ d) = 0;
    virtual void handle_polar(Polyhedron *C, int s);
    void start(Polyhedron *P, unsigned MaxRays);
    virtual void init(Polyhedron *P) {
    }

private:
    Value *current_vertex;
    ZZ sign;
};

struct reducer : public np_base {
    vec_ZZ vertex;
    //vec_ZZ den;
    ZZ num;
    mpq_t tcount;
    mpz_t tn;
    mpz_t td;
    int lower;	    // call base when only this many variables is left

    reducer(unsigned dim) : np_base(dim) {
	//den.SetLength(dim);
	mpq_init(tcount);
	mpz_init(tn);
	mpz_init(td);
    }

    ~reducer() {
	mpq_clear(tcount);
	mpz_clear(tn);
	mpz_clear(td);
    }

    virtual void handle_polar(Polyhedron *C, Value *vertex, ZZ n, ZZ d);
    void reduce(ZZ c, ZZ cd, vec_ZZ& num, mat_ZZ& den_f);
    virtual void base(ZZ& c, ZZ& cd, vec_ZZ& num, mat_ZZ& den_f) = 0;
    virtual void split(vec_ZZ& num, ZZ& num_s, vec_ZZ& num_p,
		       mat_ZZ& den_f, vec_ZZ& den_s, mat_ZZ& den_r) = 0;
};

struct ireducer : public reducer {
    ireducer(unsigned dim) : reducer(dim) {}

    virtual void split(vec_ZZ& num, ZZ& num_s, vec_ZZ& num_p,
		       mat_ZZ& den_f, vec_ZZ& den_s, mat_ZZ& den_r);
};

void normalize(ZZ& sign, ZZ& num_s, vec_ZZ& num_p, vec_ZZ& den_s, vec_ZZ& den_p,
	       mat_ZZ& f);

#endif
