#include <sstream>
#include <barvinok/options.h>
#include <barvinok/evalue.h>
#include <barvinok/util.h>
#include "conversion.h"
#include "evalue_read.h"
#include "lattice_point.h"

using namespace std;

template <typename T>
void set_from_string(T& v, char *s)
{
    istringstream str(s);
    str >> v;
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
    test_lattice_points(options);
    barvinok_options_free(options);
}
