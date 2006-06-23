#include <gmp.h>
#include <cstdlib>
#include <Rational.h>
#include <Poly.h>
#include <Matrix.h>
extern "C" {
#include <polylib/polylibgmp.h>
}

namespace polymake { namespace polytope {

::Matrix *polymake_constraints2polylib(Matrix<Rational> &F);

} }
