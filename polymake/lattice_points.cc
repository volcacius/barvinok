#include <gmp.h>
#include <cstdlib>
#include <Rational.h>
#include <Poly.h>
#include <Matrix.h>
extern "C" {
#include <polylib/polylibgmp.h>
}
#include <barvinok/barvinok.h>

namespace polymake { namespace polytope {

void lattice_points(Poly& p)
{
    Value cb;
    Matrix<Rational> F = p.give("FACETS | INEQUALITIES");
    int r = F.rows();
    int c = F.cols();
    ::Matrix *M = Matrix_Alloc(r, c+1);
    for (int i = 0; i < r; ++i) {
	value_set_si(M->p[i][0], 1);
	value_assign(M->p[i][c], mpq_numref(F[i][0].get_rep()));
	for (int j = 1; j < c; ++j)
	    value_assign(M->p[i][j], mpq_numref(F[i][j].get_rep()));
    }

    Polyhedron *A = Constraints2Polyhedron(M, 0);
    Matrix_Free(M);
    value_init(cb);
    barvinok_count(A, &cb, 0);
    Polyhedron_Free(A);
    Integer count(cb);
    value_clear(cb);
    p.take("LATTICE_POINTS") << count;
}

} }

using namespace polymake;

int main(int argc, const char *argv[])
{
    if (argc != 2) {
	return 1;
    }
    try {
	Poly p(argv[1], ios::in | ios::out);
	polytope::lattice_points(p);
    } catch (const std::exception& e) {
	cerr << e.what() << endl;
	return 1;
   }
    return 0;
}
