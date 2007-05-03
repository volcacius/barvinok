#include <sstream>
#include <barvinok/options.h>
#include <barvinok/evalue.h>
#include <barvinok/util.h>
#include "conversion.h"
#include "evalue_read.h"
#include "dpoly.h"
#include "lattice_point.h"

using namespace std;

template <typename T>
void set_from_string(T& v, char *s)
{
    istringstream str(s);
    str >> v;
}

int test_specialization(struct barvinok_options *options)
{
    Value v;
    value_init(v);
    mpq_t count;
    mpq_init(count);
    ZZ sign;

    value_set_si(v, 5);
    dpoly n(2, v);
    assert(value_cmp_si(n.coeff->p[0], 1) == 0);
    assert(value_cmp_si(n.coeff->p[1], 5) == 0);
    assert(value_cmp_si(n.coeff->p[2], 10) == 0);

    value_set_si(v, 1);
    dpoly d(2, v, 1);
    value_set_si(v, 2);
    dpoly d2(2, v, 1);
    d *= d2;
    assert(value_cmp_si(d.coeff->p[0], 2) == 0);
    assert(value_cmp_si(d.coeff->p[1], 1) == 0);
    assert(value_cmp_si(d.coeff->p[2], 0) == 0);

    sign = 1;
    n.div(d, count, sign);
    mpq_canonicalize(count);
    assert(value_cmp_si(mpq_numref(count), 31) == 0);
    assert(value_cmp_si(mpq_denref(count), 8) == 0);

    value_set_si(v, -2);
    dpoly n2(2, v);
    assert(value_cmp_si(n2.coeff->p[0], 1) == 0);
    assert(value_cmp_si(n2.coeff->p[1], -2) == 0);
    assert(value_cmp_si(n2.coeff->p[2], 3) == 0);

    n2.div(d, count, sign);
    mpq_canonicalize(count);
    assert(value_cmp_si(mpq_numref(count), 6) == 0);
    assert(value_cmp_si(mpq_denref(count), 1) == 0);

    value_clear(v);
    mpq_clear(count);
}

int test_lattice_points(struct barvinok_options *options)
{
    Param_Vertices V;
    mat_ZZ tmp;
    set_from_string(tmp, "[[0 0 0 0 4][0 0 0 0 4][-1 0 1 0 4]]");
    V.Vertex = zz2matrix(tmp);
    vec_ZZ lambda;
    set_from_string(lambda, "[3 5 7]");
    mat_ZZ rays;
    set_from_string(rays, "[[0 1 0][4 0 1][0 0 -1]]");
    term_info num;
    evalue *point[4];

    unsigned nvar, nparam;
    char **all_vars;
    point[0] = evalue_read_from_str("( -7/4 * a + ( 7/4 * c + "
		    "( 7 * {( 1/4 * a + ( 3/4 * c + 3/4 ) ) } + -21/4 ) ) )",
		    "a,b,c", &all_vars, &nvar, &nparam, options->MaxRays);
    Free_ParamNames(all_vars, nvar+nparam);
    point[1] = evalue_read_from_str("( -7/4 * a + ( 7/4 * c + "
		    "( 7 * {( 1/4 * a + ( 3/4 * c + 1/2 ) ) } + -1/2 ) ) )",
		    "a,b,c", &all_vars, &nvar, &nparam, options->MaxRays);
    Free_ParamNames(all_vars, nvar+nparam);
    point[2] = evalue_read_from_str("( -7/4 * a + ( 7/4 * c + "
		    "( 7 * {( 1/4 * a + ( 3/4 * c + 1/4 ) ) } + 17/4 ) ) )",
		    "a,b,c", &all_vars, &nvar, &nparam, options->MaxRays);
    Free_ParamNames(all_vars, nvar+nparam);
    point[3] = evalue_read_from_str("( -7/4 * a + ( 7/4 * c + "
		    "( 7 * {( 1/4 * a + ( 3/4 * c + 0 ) ) } + 9 ) ) )",
		    "a,b,c", &all_vars, &nvar, &nparam, options->MaxRays);
    Free_ParamNames(all_vars, nvar+nparam);

    lattice_point(&V, rays, lambda, &num, 4, NULL, options);
    Matrix_Free(V.Vertex);

    for (int i = 0; i < 4; ++i) {
	assert(eequal(num.E[i], point[i]));
	free_evalue_refs(point[i]);
	free(point[i]);
	free_evalue_refs(num.E[i]);
	delete num.E[i];
    }
    delete [] num.E; 
}

int main(int argc, char **argv)
{
    struct barvinok_options *options = barvinok_options_new_with_defaults();
    test_specialization(options);
    test_lattice_points(options);
    barvinok_options_free(options);
}
