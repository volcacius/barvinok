#ifndef REDUCER_H
#define REDUCER_H

#include <NTL/mat_ZZ.h>
#include <barvinok/NTL_QQ.h>
#include <barvinok/options.h>
#include "decomposer.h"
#include "dpoly.h"

#ifdef NTL_STD_CXX
using namespace NTL;
#endif

struct gen_fun;

/* base for non-parametric counting */
struct np_base : public polar_decomposer {
    unsigned dim;
    ZZ one;

    np_base(unsigned dim) {
	this->dim = dim;
	one = 1;
    }

    virtual void handle_polar(Polyhedron *C, Value *vertex, QQ c) = 0;
    virtual void handle_polar(Polyhedron *C, int s);
    virtual void start(Polyhedron *P, barvinok_options *options);
    void do_vertex_cone(const QQ& factor, Polyhedron *Cone, 
			Value *vertex, barvinok_options *options) {
	current_vertex = vertex;
	this->factor = factor;
	decompose(Cone, options);
    }
    virtual void init(Polyhedron *P) {
    }
    virtual void get_count(Value *result) {
	assert(0);
    }
    virtual ~np_base() {
    }

private:
    QQ factor;
    Value *current_vertex;
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

    virtual void handle_polar(Polyhedron *C, Value *vertex, QQ c);
    void reduce(QQ c, vec_ZZ& num, mat_ZZ& den_f);
    virtual void base(QQ& c, const vec_ZZ& num, const mat_ZZ& den_f) = 0;
    virtual void split(vec_ZZ& num, ZZ& num_s, vec_ZZ& num_p,
		       mat_ZZ& den_f, vec_ZZ& den_s, mat_ZZ& den_r) = 0;
    virtual gen_fun *get_gf() {
	assert(0);
    }
};

struct ireducer : public reducer {
    ireducer(unsigned dim) : reducer(dim) {}

    virtual void split(vec_ZZ& num, ZZ& num_s, vec_ZZ& num_p,
		       mat_ZZ& den_f, vec_ZZ& den_s, mat_ZZ& den_r);
};

void normalize(ZZ& sign, ZZ& num_s, vec_ZZ& num_p, vec_ZZ& den_s, vec_ZZ& den_p,
	       mat_ZZ& f);

// incremental counter
struct icounter : public ireducer {
    mpq_t count;

    icounter(unsigned dim) : ireducer(dim) {
	mpq_init(count);
	lower = 1;
    }
    ~icounter() {
	mpq_clear(count);
    }
    virtual void base(QQ& c, const vec_ZZ& num, const mat_ZZ& den_f);
    virtual void get_count(Value *result) {
	assert(value_one_p(&count[0]._mp_den));
	value_assign(*result, &count[0]._mp_num);
    }
};

void normalize(ZZ& sign, ZZ& num, vec_ZZ& den);

/* An incremental counter for possibly infinite sets.
 * Rather than just keeping track of the constant term
 * of the Laurent expansions, we also keep track of the
 * coefficients of negative powers.
 * If any of these is non-zero, then the counted set is infinite.
 */
struct infinite_icounter : public ireducer {
    /* an array of coefficients; count[i] is the coeffient of
     * the term with power -i.
     */
    mpq_t *count;
    unsigned len;

    infinite_icounter(unsigned dim, unsigned maxlen) : ireducer(dim), len(maxlen+1) {
	/* Not sure whether it works for dim != 1 */
	assert(dim == 1);
	count = new mpq_t[len];
	for (int i = 0; i < len; ++i)
	    mpq_init(count[i]);
	lower = 1;
    }
    ~infinite_icounter() {
	for (int i = 0; i < len; ++i)
	    mpq_clear(count[i]);
	delete [] count;
    }
    virtual void base(QQ& c, const vec_ZZ& num, const mat_ZZ& den_f);
};

#endif
