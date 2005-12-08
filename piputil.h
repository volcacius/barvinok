#include <polylib/polylibgmp.h>

struct pip_dep_result {
    Polyhedron *M;
    Polyhedron *D;
};

Polyhedron *pip_lexminmax(Polyhedron *P, int pos, int n, int nparam, int max);
struct pip_dep_result pip_inputdep(Polyhedron *D, int dim, Matrix *M);
Polyhedron *pip_projectout(Polyhedron *P, int pos, int n, int nparam);
