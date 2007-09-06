#include "reducer.h"

struct tcounter : public np_base {
    Vector *lambda;
    Matrix *den;
    ZZ sign;
    Matrix *num;
    mpq_t count;
    mpq_t tcount;
    dpoly todd;
    Vector *todd_denom;
    Value denom;

    Value tmp, tmp2;

    tcounter(unsigned dim, unsigned long max_index) : np_base(dim), todd(dim) {
	mpq_init(count);
	mpq_init(tcount);
	setup_todd(dim);
	value_init(denom);
	value_init(tmp);
	num = Matrix_Alloc(max_index, 1);
	den = Matrix_Alloc(dim, 1);
	lambda = Vector_Alloc(dim);
    }

    void setup_todd(unsigned dim);

    void adapt_todd(dpoly& t, const Value c);
    void add_powers(dpoly& n, const Value c);

    ~tcounter() {
	Matrix_Free(num);
	Matrix_Free(den);
	Vector_Free(lambda);
	mpq_clear(count);
	mpq_clear(tcount);
	Vector_Free(todd_denom);
	value_clear(denom);
	value_clear(tmp);
    }

    virtual void init(Polyhedron *P) {
	vec_ZZ l;
	randomvector(P, l, dim);
	zz2values(l, lambda->p);
    }

    virtual void handle(const mat_ZZ& rays, Value *vertex, const QQ& c,
			unsigned long det, barvinok_options *options);
    virtual void get_count(Value *result) {
	assert(value_one_p(&count[0]._mp_den));
	value_assign(*result, &count[0]._mp_num);
    }
};
