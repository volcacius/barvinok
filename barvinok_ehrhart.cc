#include <unistd.h>
#include <stdlib.h>
#include <strings.h>
#include <polylib/polylibgmp.h>
#include <barvinok/util.h>
#include <barvinok/barvinok.h>
#include "config.h"

/* The input of this example program is a polytope in PolyLib notation,
 * i.e., an n by d+2 matrix of the n constraints A x + b >= 0 defining
 * the polytope * sitting in a d-dimensional space.  The first column
 * is 1 for an inequality and 0 for an equality.  b is placed in the
 * final column.
 * Alternatively, if the matrix is preceded by the word "vertices"
 * on a line by itself, it will be interpreted as a list of vertices
 * in PolyLib notation, i.e., an n by (d+2) matrix, where n is
 * the number of vertices/rays and d the dimension.  The first column is
 * 0 for lines and 1 for vertices/rays.  The final column is the denominator
 * or 0 for rays.  Note that for barvinok_ehrhart, the first column
 * should always be 1.
 */

#ifdef HAVE_GROWING_CHERNIKOVA
#define MAXRAYS    POL_NO_DUAL
#else
#define MAXRAYS  600
#endif

#ifndef HAVE_GETOPT_H
#define getopt_long(a,b,c,d,e) getopt(a,b,c)
#else
#include <getopt.h>
struct option options[] = {
    { "convert",   no_argument,  0,  'c' },
    { "floor",     no_argument,  0,  'f' },
    { "series",    no_argument,  0,  's' },
    { "version",   no_argument,  0,  'V' },
    { 0, 0, 0, 0 }
};
#endif

static Polyhedron *Polyhedron_Read()
{
    int vertices = 0; 
    unsigned NbRows, NbColumns;
    Matrix *M;
    Polyhedron *P;
    char s[128];

    while (fgets(s, sizeof(s), stdin)) {
	if (*s == '#')
	    continue;
	if (strncasecmp(s, "vertices", sizeof("vertices")-1) == 0)
	    vertices = 1;
	if (sscanf(s, "%u %u", &NbRows, &NbColumns) == 2)
	    break;
    }
    if (feof(stdin))
	return NULL;
    M = Matrix_Alloc(NbRows,NbColumns);
    Matrix_Read_Input(M);
    if (vertices)
	P = Rays2Polyhedron(M, MAXRAYS);
    else
	P = Constraints2Polyhedron(M, MAXRAYS);
    Matrix_Free(M);
    return P;
}

/* Construct a cone over P with P placed at x_d = 1, with
 * x_d the coordinate of an extra dimension
 *
 * It's probably a mistake to depend so much on the internal
 * representation.  We should probably simply compute the
 * vertices/facets first.
 */
static Polyhedron *Cone_over_Polyhedron(Polyhedron *P)
{
    unsigned NbConstraints = 0;
    unsigned NbRays = 0;
    Polyhedron *C;
    int i;

    if (POL_HAS(P, POL_INEQUALITIES))
	NbConstraints = P->NbConstraints + 1;
    if (POL_HAS(P, POL_POINTS))
	NbRays = P->NbRays + 1;

    C = Polyhedron_Alloc(P->Dimension+1, NbConstraints, NbRays);
    if (POL_HAS(P, POL_INEQUALITIES)) {
	C->NbEq = P->NbEq;
	for (i = 0; i < P->NbConstraints; ++i)
	    Vector_Copy(P->Constraint[i], C->Constraint[i], P->Dimension+2);
	/* n >= 0 */
	value_set_si(C->Constraint[P->NbConstraints][0], 1);
	value_set_si(C->Constraint[P->NbConstraints][1+P->Dimension], 1);
    }
    if (POL_HAS(P, POL_POINTS)) {
	C->NbBid = P->NbBid;
	for (i = 0; i < P->NbRays; ++i)
	    Vector_Copy(P->Ray[i], C->Ray[i], P->Dimension+2);
	/* vertex 0 */
	value_set_si(C->Ray[P->NbRays][0], 1);
	value_set_si(C->Ray[P->NbRays][1+C->Dimension], 1);
    }
    POL_SET(C, POL_VALID);
    if (POL_HAS(P, POL_INEQUALITIES))
	POL_SET(C, POL_INEQUALITIES);
    if (POL_HAS(P, POL_POINTS))
	POL_SET(C, POL_POINTS);
    if (POL_HAS(P, POL_VERTICES))
	POL_SET(C, POL_VERTICES);
    return C;
}

int main(int argc, char **argv)
{
    Polyhedron *A, *C, *U;
    char **param_name;
    int c, ind = 0;
    int convert = 0;
    int floor = 0;
    int series = 0;

    while ((c = getopt_long(argc, argv, "sfcV", options, &ind)) != -1) {
	switch (c) {
	case 's':
	    series = 1;
	    break;
	case 'c':
	    convert = 1;
	    break;
	case 'f':
	    floor = 1;
	    break;
	case 'V':
	    printf(barvinok_version());
	    exit(0);
	    break;
	}
    }

    A = Polyhedron_Read();
    param_name = Read_ParamNames(stdin, 1);
    Polyhedron_Print(stdout, P_VALUE_FMT, A);
    C = Cone_over_Polyhedron(A);
    U = Universe_Polyhedron(1);
    if (series) {
	gen_fun *gf;
	gf = barvinok_series(C, U, MAXRAYS);
	gf->print(U->Dimension, param_name);
	puts("");
	delete gf;
    } else {
	evalue *EP;
	/* A (conceptually) obvious optimization would be to pass in
	 * the parametric vertices, which are just n times the original
	 * vertices, rather than letting barvinok_enumerate_ev (re)compute
	 * them through Polyhedron2Param_SimplifiedDomain.
	 */
	EP = barvinok_enumerate_ev(C, U, MAXRAYS);
	print_evalue(stdout, EP, param_name);
	if (floor) {
	    fprintf(stderr, "WARNING: floor conversion not supported\n");
	    evalue_frac2floor(EP);
	    print_evalue(stdout, EP, param_name);
	} else if (convert) {
	    evalue_mod2table(EP, C->Dimension);
	    print_evalue(stdout, EP, param_name);
	}
	free_evalue_refs(EP);
	free(EP);
    }
    Free_ParamNames(param_name, 1);
    Polyhedron_Free(A);
    Polyhedron_Free(C);
    Polyhedron_Free(U);
    return 0;
}
