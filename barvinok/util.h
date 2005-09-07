#ifndef UTIL_H
#define UTIL_H

#if defined(__cplusplus)
extern "C" {
#endif

#include <polylib/polylibgmp.h>
#include <barvinok/evalue.h>

#ifdef POL_NO_DUAL
#define emptyQ2(P)							\
	((F_ISSET(P, POL_INEQUALITIES) && P->NbEq > P->Dimension) ||	\
	 (F_ISSET(P, POL_POINTS) && P->NbRays == 0))
#else
#define POL_NO_DUAL		0
#define emptyQ2(P)		emptyQ(P)
#define POL_ENSURE_FACETS(P)	/* nothing */
#define POL_ENSURE_VERTICES(P)	/* nothing */
#endif

void value_lcm(Value i, Value j, Value* lcm);
int random_int(int max);
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
evalue* ParamLine_Length(Polyhedron *P, Polyhedron *C, unsigned MaxRays);
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

/* only defined if PolyLib has RankingConstraints */
evalue *barvinok_ranking_ev(Polyhedron *P, Polyhedron *D, Polyhedron *C, 
			    unsigned MaxRays);
Enumeration *barvinok_ranking(Polyhedron *P, Polyhedron *D, Polyhedron *C, 
			      unsigned MaxRays);

const char *barvinok_version();

#if defined(__cplusplus)
}
#endif

#endif
