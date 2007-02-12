#include <barvinok/polylib.h>

/*
 * Scales the parametric polyhedron such that all vertices are integer.
 */
void Param_Polyhedron_Scale_Integer(Param_Polyhedron *PP, Polyhedron **P,
				    Value *det, unsigned MaxRays)
{
  int i;
  int nb_param, nb_vars;
  Vector *denoms;
  Param_Vertices *V;
  Value global_var_lcm;
  Matrix *expansion;

  value_set_si(*det, 1);
  if (!PP->nbV)
    return;

  nb_param = PP->D->Domain->Dimension;
  nb_vars = PP->V->Vertex->NbRows;

  /* Scan the vertices and make an orthogonal expansion of the variable
     space */
  /* a- prepare the array of common denominators */
  denoms = Vector_Alloc(nb_vars);
  value_init(global_var_lcm);

  /* b- scan the vertices and compute the variables' global lcms */
  for (V = PP->V; V; V = V->next)
    for (i = 0; i < nb_vars; i++)
      Lcm3(denoms->p[i], V->Vertex->p[i][nb_param+1], &denoms->p[i]);

  value_set_si(global_var_lcm, 1);
  for (i = 0; i < nb_vars; i++) {
    value_multiply(*det, *det, denoms->p[i]);
    Lcm3(global_var_lcm, denoms->p[i], &global_var_lcm);
  }

  /* scale vertices */
  for (V = PP->V; V; V = V->next)
    for (i = 0; i < nb_vars; i++) {
      Vector_Scale(V->Vertex->p[i], V->Vertex->p[i], denoms->p[i], nb_param+1);
      Vector_AntiScale(V->Vertex->p[i], V->Vertex->p[i],
		       V->Vertex->p[i][nb_param+1], nb_param+2);
    }

  /* the expansion can be actually writen as global_var_lcm.L^{-1} */
  /* this is equivalent to multiply the rows of P by denoms_det */
  for (i = 0; i < nb_vars; i++)
    value_division(denoms->p[i], global_var_lcm, denoms->p[i]);

  /* OPT : we could use a vector instead of a diagonal matrix here (c- and d-).*/
  /* c- make the quick expansion matrix */
  expansion = Matrix_Alloc(nb_vars+nb_param+1, nb_vars+nb_param+1);
  for (i = 0; i < nb_vars; i++)
    value_assign(expansion->p[i][i], denoms->p[i]);
  for (i = nb_vars; i < nb_vars+nb_param+1; i++)
    value_assign(expansion->p[i][i], global_var_lcm);

  /* d- apply the variable expansion to the polyhedron */
  if (P)
    *P = Polyhedron_Preimage(*P, expansion, MaxRays);

  Matrix_Free(expansion);
  value_clear(global_var_lcm);
  Vector_Free(denoms);
}
