#include <sstream>
#include <barvinok/options.h>
#include <barvinok/evalue.h>
#include <barvinok/util.h>
#include "conversion.h"
#include "evalue_read.h"
#include "dpoly.h"
#include "lattice_point.h"
#include "tcounter.h"

using namespace std;

template <typename T>
void set_from_string(T& v, char *s)
{
    istringstream str(s);
    str >> v;
}

int test_evalue(struct barvinok_options *options)
{
    unsigned nvar, nparam;
    char **all_vars;
    evalue *poly1, poly2;

    poly1 = evalue_read_from_str("(1/4 * n^4 + 1/2 * n^3 + 1/4 * n^2)",
				 NULL, &all_vars, &nvar, &nparam,
				 options->MaxRays);
    Free_ParamNames(all_vars, nvar+nparam);

    value_init(poly2.d);
    evalue_copy(&poly2, poly1);
    evalue_negate(poly1);
    eadd(&poly2, poly1);
    reduce_evalue(poly1);
    assert(EVALUE_IS_ZERO(*poly1));
    free_evalue_refs(poly1);
    free_evalue_refs(&poly2);
    free(poly1);
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

int test_todd(struct barvinok_options *options)
{
    tcounter t(2);
    assert(value_cmp_si(t.todd.coeff->p[0], 1) == 0);
    assert(value_cmp_si(t.todd.coeff->p[1], -3) == 0);
    assert(value_cmp_si(t.todd.coeff->p[2], 3) == 0);
    assert(value_cmp_si(t.todd_denom->p[0], 1) == 0);
    assert(value_cmp_si(t.todd_denom->p[1], 6) == 0);
    assert(value_cmp_si(t.todd_denom->p[2], 36) == 0);

    set_from_string(t.lambda, "[1 -1]");

    mat_ZZ rays;
    set_from_string(rays, "[[-1 0][-1 1]]");

    QQ one(1, 1);

    vec_ZZ v;
    set_from_string(v, "[2 0 1]");
    Vector *vertex = Vector_Alloc(3);
    zz2values(v, vertex->p);

    t.handle(rays, vertex->p, one, 1, NULL, options);
    assert(value_cmp_si(mpq_numref(t.count), 71) == 0);
    assert(value_cmp_si(mpq_denref(t.count), 24) == 0);

    set_from_string(rays, "[[0 -1][1 -1]]");
    set_from_string(v, "[0 2 1]");
    zz2values(v, vertex->p);

    t.handle(rays, vertex->p, one, 1, NULL, options);
    assert(value_cmp_si(mpq_numref(t.count), 71) == 0);
    assert(value_cmp_si(mpq_denref(t.count), 12) == 0);

    set_from_string(rays, "[[1 0][0 1]]");
    set_from_string(v, "[0 0 1]");
    zz2values(v, vertex->p);

    t.handle(rays, vertex->p, one, 1, NULL, options);
    assert(value_cmp_si(mpq_numref(t.count), 6) == 0);
    assert(value_cmp_si(mpq_denref(t.count), 1) == 0);

    Vector_Free(vertex);
}

int main(int argc, char **argv)
{
    struct barvinok_options *options = barvinok_options_new_with_defaults();
    test_evalue(options);
    test_specialization(options);
    test_lattice_points(options);
    test_todd(options);
    barvinok_options_free(options);
}
