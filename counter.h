#include "reducer.h"

struct counter : public np_base {
    vec_ZZ lambda;
    mat_ZZ vertex;
    vec_ZZ den;
    ZZ sign;
    vec_ZZ num;
    ZZ offset;
    int j;
    mpq_t count;
    Value tz;

    counter(unsigned dim) : np_base(dim) {
	den.SetLength(dim);
	mpq_init(count);
	value_init(tz);
    }

    virtual void init(Polyhedron *P) {
	randomvector(P, lambda, dim);
    }

    virtual void reset() {
	mpq_set_si(count, 0, 0);
    }

    ~counter() {
	mpq_clear(count);
	value_clear(tz);
    }

    virtual void handle(const mat_ZZ& rays, Value *vertex, const QQ& c,
			unsigned long det, int *closed, barvinok_options *options);
    virtual void get_count(Value *result) {
	assert(value_one_p(&count[0]._mp_den));
	value_assign(*result, &count[0]._mp_num);
    }
};
