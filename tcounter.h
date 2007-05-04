#include "reducer.h"

struct tcounter : public np_base {
    vec_ZZ lambda;
    mat_ZZ vertex;
    vec_ZZ den;
    ZZ sign;
    vec_ZZ num;
    mpq_t count;
    mpq_t tcount;
    Value tz;
    dpoly todd;
    Vector *todd_denom;
    Value denom;

    Value tmp, tmp2;

    tcounter(unsigned dim) : np_base(dim), todd(dim) {
	den.SetLength(dim);
	mpq_init(count);
	mpq_init(tcount);
	value_init(tz);
	setup_todd(dim);
	value_init(denom);
	value_init(tmp);
    }

    void setup_todd(unsigned dim);

    void adapt_todd(dpoly& t, const Value c);
    void add_powers(dpoly& n, const Value c);

    ~tcounter() {
	mpq_clear(count);
	mpq_clear(tcount);
	value_clear(tz);
	Vector_Free(todd_denom);
	value_clear(denom);
	value_clear(tmp);
    }

    virtual void init(Polyhedron *P) {
	randomvector(P, lambda, dim);
    }

    virtual void handle(const mat_ZZ& rays, Value *vertex, const QQ& c,
			unsigned long det, int *closed, barvinok_options *options);
    virtual void get_count(Value *result) {
	assert(value_one_p(&count[0]._mp_den));
	value_assign(*result, &count[0]._mp_num);
    }
};
