#include <barvinok/evalue.h>

#if defined(__cplusplus)
extern "C" {
#endif

struct barvinok_options;

struct bernoulli {
    Vector  *b_num;
    Vector  *b_den;
    Vector  **sum;
    int	    size;
    int	    n;
};

struct bernoulli *bernoulli_compute(int n);

evalue *Bernoulli_sum(Polyhedron *P, Polyhedron *C,
			struct barvinok_options *options);
evalue *Bernoulli_sum_evalue(evalue *e, unsigned nvar,
			     struct barvinok_options *options);

#if defined(__cplusplus)
}
#endif
