#include "reducer.h"

struct counter : public np_base {
    Vector *lambda;
    Matrix *vertex;
    Vector *den;
    ZZ sign;
    Vector *num;
    int j;
    mpq_t count;

    counter(unsigned dim, unsigned long max_index) : np_base(dim) {
	mpq_init(count);
	vertex = Matrix_Alloc(max_index, dim);
	num = Vector_Alloc(max_index);
	den = Vector_Alloc(dim);
	lambda = Vector_Alloc(dim);
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
	Matrix_Free(vertex);
	Vector_Free(num);
	Vector_Free(den);
	Vector_Free(lambda);
	mpq_clear(count);
    }

    virtual void handle(const mat_ZZ& rays, Value *vertex, const QQ& c,
			unsigned long det, int *closed, barvinok_options *options);
    virtual void get_count(Value *result) {
	assert(value_one_p(&count[0]._mp_den));
	value_assign(*result, &count[0]._mp_num);
    }
};
