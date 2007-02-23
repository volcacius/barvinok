#ifndef UTIL_H
#define UTIL_H

#include <barvinok/polylib.h>

#if defined(__cplusplus)
extern "C" {
#endif

#include <barvinok/evalue.h>

#ifdef POL_NO_DUAL
#define emptyQ2(P)							\
	((F_ISSET(P, POL_INEQUALITIES) && P->NbEq > P->Dimension) ||	\
	 (F_ISSET(P, POL_POINTS) && P->NbRays == 0))
#define POL_HAS(P,f)		(F_ISSET(P, f))
#define POL_SET(P,f)		(F_SET(P, f))
#else
#define POL_NO_DUAL		0
#define emptyQ2(P)		emptyQ(P)
#define POL_ENSURE_FACETS(P)	/* nothing */
#define POL_ENSURE_VERTICES(P)	/* nothing */
#define	POL_INEQUALITIES	0
#define	POL_POINTS		0
#define	POL_FACETS		0
#define	POL_VERTICES		0
#define	POL_VALID		0
#define POL_HAS(P,f)		(1)
#define POL_SET(P,f)		/* nothing */
#endif

#ifndef POL_INTEGER
#define POL_INTEGER		0
#endif

void value_lcm(const Value i, const Value j, Value* lcm);
int random_int(int max);
Polyhedron *Polyhedron_Read(unsigned MaxRays);
Polyhedron* Polyhedron_Polar(Polyhedron *P, unsigned NbMaxRays);
void Polyhedron_Polarize(Polyhedron *P);
Polyhedron* supporting_cone(Polyhedron *P, int v);
Polyhedron* supporting_cone_p(Polyhedron *P, Param_Vertices *v);
Polyhedron* triangulate_cone(Polyhedron *P, unsigned NbMaxCons);
void check_triangulization(Polyhedron *P, Polyhedron *T);
Polyhedron *remove_equalities(Polyhedron *P);
Polyhedron *remove_equalities_p(Polyhedron *P, unsigned nvar, Matrix **factor);
void manual_count(Polyhedron *P, Value* result);
Polyhedron* Polyhedron_Factor(Polyhedron *P, unsigned nparam, 
			      unsigned NbMaxRays);
void Line_Length(Polyhedron *P, Value *len);
struct barvinok_options;
evalue* ParamLine_Length(Polyhedron *P, Polyhedron *C,
			 struct barvinok_options *options);
void Extended_Euclid(Value a, Value b, Value *x, Value *y, Value *g);
Matrix * unimodular_complete(Vector *row);
Bool isIdentity(Matrix *M);
void Param_Polyhedron_Print(FILE* DST, Param_Polyhedron *PP, char **param_names);
void Enumeration_Print(FILE *Dst, Enumeration *en, char **params);
void Enumeration_Free(Enumeration *en);
void Enumeration_mod2table(Enumeration *en, unsigned nparam);
size_t Enumeration_size(Enumeration *en);
void Free_ParamNames(char **params, int m);
int DomainIncludes(Polyhedron *Pol1, Polyhedron *Pol2);
#ifndef HAVE_DOMAINCONSTRAINTSIMPLIFY
int ConstraintSimplify(Value *old, Value *n, int len, Value* v);
Polyhedron *DomainConstraintSimplify(Polyhedron *P, unsigned MaxRays);
#endif
int line_minmax(Polyhedron *I, Value *min, Value *max);
void count_points_e (int pos, Polyhedron *P, int exist, int nparam,
		     Value *context, Value *res);
int DomainContains(Polyhedron *P, Value *list_args, int len, 
		   unsigned MaxRays, int set);
Polyhedron* Polyhedron_Project(Polyhedron *P, int dim);
Polyhedron *DomainConcat(Polyhedron *head, Polyhedron *tail);

/* only defined if PolyLib has RankingConstraints */
evalue *barvinok_lexsmaller_ev(Polyhedron *P, Polyhedron *D, unsigned dim,
			       Polyhedron *C, unsigned MaxRays);
Enumeration *barvinok_lexsmaller(Polyhedron *P, Polyhedron *D, unsigned dim,
				 Polyhedron *C, unsigned MaxRays);

Matrix *align_matrix(Matrix *M, int nrows);

char **util_generate_names(int n, char *prefix);
void util_free_names(int n, char **names);
void Polyhedron_pprint(FILE *out, Polyhedron *P, int dim, int nparam,
		       char **iter_names, char **param_names);

Polyhedron *Cone_over_Polyhedron(Polyhedron *P);

Matrix *compress_variables(Matrix *Equalities, unsigned nparam);
Matrix *left_inverse(Matrix *M, Matrix **Eq);

int Polyhedron_has_revlex_positive_rays(Polyhedron *P, unsigned nparam);

const char *barvinok_version();

#if defined(__cplusplus)
}
#endif

#endif
