#ifndef REDUCER_H
#define REDUCER_H

#include <assert.h>
#include <NTL/mat_ZZ.h>
#include <barvinok/NTL_QQ.h>
#include <barvinok/options.h>
#include "decomposer.h"
#include "dpoly.h"

using namespace NTL;

struct gen_fun;

extern struct OrthogonalException {} Orthogonal;

/* base for non-parametric counting */
struct np_base : public signed_cone_consumer {
    unsigned dim;
    ZZ one;

    np_base(unsigned dim) {
	assert(dim > 0);
	this->dim = dim;
	one = 1;
    }

    virtual void handle(const mat_ZZ& rays, Value *vertex, const QQ& c,
			unsigned long det,
			barvinok_options *options) = 0;
    virtual void handle(const signed_cone& sc, barvinok_options *options);
    virtual void start(Polyhedron *P, barvinok_options *options);
    void do_vertex_cone(const QQ& factor, Polyhedron *Cone, 
			Value *vertex, barvinok_options *options) {
	current_vertex = vertex;
	this->factor = factor;
	barvinok_decompose(Cone, *this, options);
    }
    virtual void init(Polyhedron *P, int n_try) {
    }
    virtual void reset() {
	assert(0);
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
    mat_ZZ vertex;
    //vec_ZZ den;
    ZZ num;
    mpq_t tcount;
    mpz_t tn;
    mpz_t td;
    int lower;	    // call base when only this many variables is left
    Value tz;

    reducer(unsigned dim) : np_base(dim) {
	vertex.SetDims(1, dim);
	//den.SetLength(dim);
	mpq_init(tcount);
	mpz_init(tn);
	mpz_init(td);
	value_init(tz);
    }

    ~reducer() {
	value_clear(tz);
	mpq_clear(tcount);
	mpz_clear(tn);
	mpz_clear(td);
    }

    virtual void handle(const mat_ZZ& rays, Value *vertex, const QQ& c,
			unsigned long det, barvinok_options *options);
    void reduce(const vec_QQ& c, const mat_ZZ& num, const mat_ZZ& den_f);
    void reduce_canonical(const vec_QQ& c, const mat_ZZ& num,
			    const mat_ZZ& den_f);
    virtual void base(const QQ& c, const vec_ZZ& num, const mat_ZZ& den_f) = 0;
    virtual void base(const vec_QQ& c, const mat_ZZ& num, const mat_ZZ& den_f);
    virtual void split(const mat_ZZ& num, vec_ZZ& num_s, mat_ZZ& num_p,
		       const mat_ZZ& den_f, vec_ZZ& den_s, mat_ZZ& den_r) = 0;
    virtual gen_fun *get_gf() {
	assert(0);
	return NULL;
    }
};

void split_one(const mat_ZZ& num, vec_ZZ& num_s, mat_ZZ& num_p,
	       const mat_ZZ& den_f, vec_ZZ& den_s, mat_ZZ& den_r);

struct ireducer : public reducer {
    ireducer(unsigned dim) : reducer(dim) {}

    virtual void split(const mat_ZZ& num, vec_ZZ& num_s, mat_ZZ& num_p,
		       const mat_ZZ& den_f, vec_ZZ& den_s, mat_ZZ& den_r) {
	split_one(num, num_s, num_p, den_f, den_s, den_r);
    }
};

void normalize(ZZ& sign, vec_ZZ& num_s, mat_ZZ& num_p, vec_ZZ& den_s, vec_ZZ& den_p,
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
    virtual void base(const QQ& c, const vec_ZZ& num, const mat_ZZ& den_f);
    virtual void get_count(Value *result) {
	assert(value_one_p(&count[0]._mp_den));
	value_assign(*result, &count[0]._mp_num);
    }
};

#endif
