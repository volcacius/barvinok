#ifndef UTIL_H
#define UTIL_H

#if defined(__cplusplus)
extern "C" {
#endif

#include <polylib/polylibgmp.h>

void value_lcm(Value i, Value j, Value* lcm);
int random_int(int max);
Polyhedron* Polyhedron_Polar(Polyhedron *P, unsigned NbMaxRays);
void Polyhedron_Polarize(Polyhedron *P);
Polyhedron* supporting_cone(Polyhedron *P, int v);
Polyhedron* supporting_cone_p(Polyhedron *P, Param_Vertices *v);
Polyhedron* triangularize_cone(Polyhedron *P, unsigned NbMaxCons);
void check_triangulization(Polyhedron *P, Polyhedron *T);
Polyhedron *remove_equalities(Polyhedron *P);
Polyhedron *remove_equalities_p(Polyhedron *P, unsigned nvar, Matrix **factor);
void manual_count(Polyhedron *P, Value* result);
Polyhedron* Polyhedron_Reduce(Polyhedron *P, Value* factor);
Matrix * unimodular_complete(Vector *row);
Bool isIdentity(Matrix *M);
void Param_Polyhedron_Print(FILE* DST, Param_Polyhedron *PP, char **param_names);
void Enumeration_Print(FILE *Dst, Enumeration *en, char **params);
void Enumeration_mod2table(Enumeration *en, unsigned nparam);
void Free_ParamNames(char **params, int m);

#include "ev_operations.h"

Polyhedron* ParamPolyhedron_Reduce(Polyhedron *P, unsigned nvar, 
				   evalue* factor);

#if defined(__cplusplus)
}
#endif

#endif
