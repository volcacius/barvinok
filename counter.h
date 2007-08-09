#include "reducer.h"

struct counter : public np_base {
    Vector *lambda;
    Matrix *den;
    Matrix *num;
    mpq_t count;
    Value tmp;

    counter(unsigned dim, unsigned long max_index) : np_base(dim) {
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

    virtual void reset() {
	mpq_set_si(count, 0, 0);
    }

    ~counter() {
	Matrix_Free(num);
	Matrix_Free(den);
	Vector_Free(lambda);
	mpq_clear(count);
	value_clear(tmp);
    }

    void add_falling_powers(dpoly& n, Value c);
    virtual void handle(const mat_ZZ& rays, Value *vertex, const QQ& c,
			unsigned long det, int *closed, barvinok_options *options);
    virtual void get_count(Value *result) {
	assert(value_one_p(&count[0]._mp_den));
	value_assign(*result, &count[0]._mp_num);
    }
};
