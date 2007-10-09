#include "reducer.h"

struct counter_base: public np_base {
    Vector *lambda;
    Matrix *den;
    Matrix *num;
    mpq_t count;
    Value tmp;

    counter_base(unsigned dim, unsigned long max_index) : np_base(dim) {
	mpq_init(count);
	num = Matrix_Alloc(max_index, 1);
	den = Matrix_Alloc(dim, 1);
	lambda = Vector_Alloc(dim);
	value_init(tmp);
    }

    virtual void init(Polyhedron *P) {
	vec_ZZ l;
	randomvector(P, l, dim);
	zz2values(l, lambda->p);
    }

    ~counter_base() {
	Matrix_Free(num);
	Matrix_Free(den);
	Vector_Free(lambda);
	mpq_clear(count);
	value_clear(tmp);
    }

    virtual void add_lattice_points(int sign) = 0;

    virtual void handle(const mat_ZZ& rays, Value *vertex, const QQ& c,
			unsigned long det, barvinok_options *options);
    virtual void get_count(Value *result) {
	assert(value_one_p(&count[0]._mp_den));
	value_assign(*result, &count[0]._mp_num);
    }
};

struct counter : public counter_base {
    counter(unsigned dim, unsigned long max_index) :
	counter_base(dim, max_index) {}

    virtual void reset() {
	mpq_set_si(count, 0, 0);
    }

    virtual void add_lattice_points(int sign);
};

struct tcounter : public counter_base {
    mpq_t tcount;
    dpoly todd;
    Vector *todd_denom;
    Value denom;

    tcounter(unsigned dim, unsigned long max_index) :
		counter_base(dim, max_index), todd(dim) {
	mpq_init(tcount);
	setup_todd(dim);
	value_init(denom);
    }

    void setup_todd(unsigned dim);

    void adapt_todd(dpoly& t, const Value c);
    void add_powers(dpoly& n, const Value c);

    ~tcounter() {
	mpq_clear(tcount);
	Vector_Free(todd_denom);
	value_clear(denom);
    }

    virtual void add_lattice_points(int sign);
};
