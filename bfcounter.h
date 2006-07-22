#include <vector>
#include <gmp.h>
extern "C" {
#include <polylib/polylibgmp.h>
}
#include <barvinok/NTL_QQ.h>
#include "reducer.h"

struct bfc_term_base {
    // the number of times a given factor appears in the denominator
    int	    *powers;
    mat_ZZ   terms;

    bfc_term_base(int len) {
	powers = new int[len];
    }

    virtual ~bfc_term_base() {
	delete [] powers;
    }
};

struct bfc_term : public bfc_term_base {
    vec_QQ   c;

    bfc_term(int len) : bfc_term_base(len) {}
};

typedef std::vector< bfc_term_base * > bfc_vec;

struct bf_reducer;

struct bf_base : public np_base {
    ZZ one;
    mpq_t tcount;
    mpz_t tn;
    mpz_t td;
    int lower;	    // call base when only this many variables is left

    bf_base(unsigned dim) : np_base(dim) {
	mpq_init(tcount);
	mpz_init(tn);
	mpz_init(td);
	one = 1;
    }

    ~bf_base() {
	mpq_clear(tcount);
	mpz_clear(tn);
	mpz_clear(td);
    }

    virtual void handle_polar(Polyhedron *C, Value *vertex, QQ c);
    int setup_factors(Polyhedron *P, mat_ZZ& factors, bfc_term_base* t, int s);

    bfc_term_base* find_bfc_term(bfc_vec& v, int *powers, int len);
    void add_term(bfc_term_base *t, vec_ZZ& num1, vec_ZZ& num);
    void add_term(bfc_term_base *t, vec_ZZ& num);

    void reduce(mat_ZZ& factors, bfc_vec& v);
    virtual void base(mat_ZZ& factors, bfc_vec& v) = 0;

    virtual bfc_term_base* new_bf_term(int len) = 0;
    virtual void set_factor(bfc_term_base *t, int k, int change) = 0;
    virtual void set_factor(bfc_term_base *t, int k, mpq_t &f, int change) = 0;
    virtual void set_factor(bfc_term_base *t, int k, const QQ& c, int change) = 0;
    virtual void update_term(bfc_term_base *t, int i) = 0;
    virtual void insert_term(bfc_term_base *t, int i) = 0;
    virtual bool constant_vertex(int dim) = 0;
    virtual void cum(bf_reducer *bfr, bfc_term_base *t, int k, 
		     dpoly_r *r) = 0;
};

struct bfcounter_base : public bf_base {
    QQ c;

    bfcounter_base(unsigned dim) : bf_base(dim) {
    }

    bfc_term_base* new_bf_term(int len);

    virtual void set_factor(bfc_term_base *t, int k, int change);
    virtual void set_factor(bfc_term_base *t, int k, mpq_t &f, int change);
    virtual void set_factor(bfc_term_base *t, int k, const QQ& c_factor, int change);
    virtual void insert_term(bfc_term_base *t, int i);
    virtual void update_term(bfc_term_base *t, int i);

    virtual bool constant_vertex(int dim) { return true; }
    virtual void cum(bf_reducer *bfr, bfc_term_base *t, int k, dpoly_r *r) {
	assert(0);
    }
};
