#include <bernstein/bernstein.h>
#include "ex_convert.h"

using namespace GiNaC;
using namespace bernstein;

#define ALLOC(type) (type*)malloc(sizeof(type))

static evalue *ex2evalue(const ex& poly, const exvector& params, int pos)
{
    if (pos >= params.size()) {
	evalue *c = ALLOC(evalue);
	value_init(c->d);
	value_init(c->x.n);
	assert(is_a<numeric>(poly));
	numeric2value(ex_to<numeric>(poly).numer(), c->x.n);
	numeric2value(ex_to<numeric>(poly).denom(), c->d);
	return c;
    }

    evalue *v = evalue_var(pos);
    evalue *sum = ex2evalue(poly.coeff(params[pos], poly.degree(params[pos])),
			    params, pos+1);
    for (int i = poly.degree(params[pos])-1; i >= 0; --i) {
	evalue *t = ex2evalue(poly.coeff(params[pos], i), params, pos+1);
	emul(v, sum);
	eadd(t, sum);
	evalue_free(t);
    }
    evalue_free(v);
    return sum;
}

evalue *ex2evalue(const ex& poly, const exvector& params)
{
    return ex2evalue(poly, params, 0);
}
