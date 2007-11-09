#include <assert.h>
#include <stdlib.h>
#include <barvinok/options.h>
#include <barvinok/util.h>
#include "hilbert.h"
#include "hull.h"
#include "lattice_width.h"
#include "param_util.h"
#include "reduce_domain.h"

#define ALLOC(type) (type*)malloc(sizeof(type))
#define ALLOCN(type,n) (type*)malloc((n) * sizeof(type))

static void clear_width_direction(struct width_direction *wd)
{
    Vector_Free(wd->width);
    Vector_Free(wd->dir);
    if (wd->domain)
	Polyhedron_Free(wd->domain);
}

static struct width_direction_array *new_width_direction_array(void)
{
    struct width_direction_array *dirs = ALLOC(struct width_direction_array);
    assert(dirs);
    dirs->n = 0;
    dirs->alloc = 32;
    dirs->wd = ALLOCN(struct width_direction, dirs->alloc);
    assert(dirs->wd);
    return dirs;
}

static void grow_width_direction_array(struct width_direction_array *dirs,
				       int extra)
{
    if (dirs->n + extra <= dirs->alloc)
	return;
    dirs->alloc = (5*(dirs->n+extra))/4;
    dirs->wd = (struct width_direction*)realloc(dirs->wd,
			    dirs->alloc * sizeof(struct width_direction));
    assert(dirs->wd);
}

void free_width_direction_array(struct width_direction_array *dirs)
{
    int i;

    for (i = 0; i < dirs->n; ++i)
	clear_width_direction(&dirs->wd[i]);
    free(dirs->wd);
    free(dirs);
}

#define INT_BITS (sizeof(unsigned) * 8)

/* For each parametric vertex, compute cone of directions
 * for which this vertex attains the minimal value.
 */
static Matrix **compute_vertex_dirs(Param_Polyhedron *PP)
{
    int i;
    unsigned nvar = PP->V->Vertex->NbRows;
    Param_Vertices *V;
    Matrix **vertex_dirs = ALLOCN(Matrix *, PP->nbV);

    for (i = 0, V = PP->V; V; ++i, V = V->next) {
	int kx;
	unsigned bx;
	int j, k;
	unsigned *facets;
	int n;
	Matrix *M;
	Polyhedron *P;

	if (V->Facets) {
	    int len = (PP->Constraints->NbRows+INT_BITS-1)/INT_BITS;
	    facets = V->Facets;
	    n = bit_vector_count(facets, len);
	} else
	    facets = supporting_constraints(PP->Constraints, V, &n);
	M = Matrix_Alloc(n, 1+nvar+1);
	for (k = 0, j = 0, kx = 0, bx = MSB; j < n; ++k) {
	    if (facets[kx] & bx) {
		value_set_si(M->p[j][0], 1);
		Vector_Copy(PP->Constraints->p[k]+1, M->p[j++]+1, nvar);
	    }
	    NEXT(kx, bx);
	}
	P = Constraints2Polyhedron(M, 0);
	Matrix_Free(M);
	vertex_dirs[i] = Matrix_Alloc(P->NbRays-1, nvar);
	for (k = 0, j = 0; k < P->NbRays; ++k) {
	    if (value_notzero_p(P->Ray[k][1+nvar]))
		continue;
	    Vector_Copy(P->Ray[k]+1, vertex_dirs[i]->p[j++], nvar);
	}
	Polyhedron_Free(P);
	if (!V->Facets)
	    free(facets);
    }
    return vertex_dirs;
}

/* Compute
 * 
 *	 c     a     b
 *	--- = --- - ---
 *	c_d   a_d   b_d
 */
static void Vector_Subtract(Value *a, Value a_d,
			    Value *b, Value b_d,
			    Value *c, Value *c_d, int len)
{
    Value ma, mb;
    value_init(ma);
    value_init(mb);
    value_lcm(*c_d, a_d, b_d);
    value_divexact(ma, *c_d, a_d);
    value_divexact(mb, *c_d, b_d);
    value_oppose(mb, mb);
    Vector_Combine(a, b, c, ma, mb, len);
    value_clear(ma);
    value_clear(mb);
}

/* Compute width for a given direction dir and initialize width_direction
 * structure.
 */
static void compute_width_direction(Matrix *V_min, Matrix *V_max,
				    Value *dir, struct width_direction *wd)
{
    Vector *max = Vector_Alloc(V_min->NbColumns);
    unsigned nvar = V_min->NbRows;
    unsigned nparam = V_min->NbColumns-2;

    wd->width = Vector_Alloc(V_min->NbColumns);
    wd->dir = Vector_Alloc(nvar);
    Vector_Copy(dir, wd->dir->p, nvar);
    wd->domain = NULL;

    V_min->NbColumns--;
    V_max->NbColumns--;

    Vector_Matrix_Product(dir, V_max, max->p);
    Vector_Matrix_Product(dir, V_min, wd->width->p);
    Vector_Subtract(max->p, V_max->p[0][V_max->NbColumns],
		    wd->width->p, V_min->p[0][V_min->NbColumns],
		    wd->width->p, &wd->width->p[nparam+1],
		    nparam+1);

    V_min->NbColumns++;
    V_max->NbColumns++;

    Vector_Normalize(wd->width->p, nparam+2);

    Vector_Free(max);
}

static int Vector_Compare(Value *p1, Value *p2, unsigned len)
{
    int i;

    for (i = 0; i < len; ++i) {
	int sign = mpz_cmp(p1[i], p2[i]);
	if (sign)
	    return sign;
    }
    return 0;
}

static int wd_width_lex_cmp(const void *va, const void *vb)
{
    const struct width_direction *a = (const struct width_direction *)va;
    const struct width_direction *b = (const struct width_direction *)vb;

    return Vector_Compare(a->width->p, b->width->p, a->width->Size);
}

static int wd_dir_lex_cmp(const void *va, const void *vb)
{
    const struct width_direction *a = (const struct width_direction *)va;
    const struct width_direction *b = (const struct width_direction *)vb;

    return Vector_Compare(a->dir->p, b->dir->p, a->dir->Size);
}

static int add_vertex(Matrix *M, int n, Value *v)
{
    if (n >= M->NbRows)
	Matrix_Extend(M, 3*(M->NbRows+10)/2);
    value_set_si(M->p[n][0], 1);
    Vector_Copy(v, M->p[n]+1, M->NbColumns-2);
    value_set_si(M->p[n][M->NbColumns-1], 1);
    return n+1;
}

/* Puts the points in v that lie in P in front of the list
 * and returns their number.
 */
static int valid_vertices(Polyhedron *P, Matrix *v, int n_v)
{
    int i, j, k;
    Value tmp;

    assert(v->NbColumns == P->Dimension+2);
    value_init(tmp);

    for (j = 0, k = 0; j < n_v; ++j) {
	for (i = 0; i < P->NbConstraints; ++i) {
	    Inner_Product(v->p[j]+1, P->Constraint[i]+1, P->Dimension+1, &tmp);
	    if (value_neg_p(tmp))
		break;
	}
	if (i < P->NbConstraints)
	    continue;
	if (j != k)
	    Vector_Exchange(v->p[j]+1, v->p[k]+1, P->Dimension);
	++k;
    }

    value_clear(tmp);
    return k;
}

static struct width_direction_array *
compute_width_directions(Param_Polyhedron *PP, struct barvinok_options *options)
{
    Matrix **vertex_dirs;
    Param_Vertices *V_max, *V_min;
    int i, V_max_i, V_min_i;
    unsigned nvar = PP->V->Vertex->NbRows;
    struct width_direction_array *width_dirs = new_width_direction_array();
    Matrix *all_vertices = Matrix_Alloc(nvar, 1+nvar+1);
    int n_vertices = 0;

    vertex_dirs = compute_vertex_dirs(PP);

    for (V_max = PP->V; V_max; V_max = V_max->next)
	Param_Vertex_Common_Denominator(V_max);

    for (V_max = PP->V, V_max_i = 0; V_max; V_max = V_max->next, V_max_i++) {
	for (V_min = V_max->next, V_min_i = V_max_i+1;
		V_min;
		V_min = V_min->next, V_min_i++) {
	    Matrix *M;
	    Matrix *basis;
	    Polyhedron *C;
	    unsigned V_max_n = vertex_dirs[V_max_i]->NbRows;
	    unsigned V_min_n = vertex_dirs[V_min_i]->NbRows;
	    int sorted_n;
	    int n_valid;

	    if (options->verbose)
		fprintf(stderr, "%d/%d %d/%d %d \r",
				    V_max_i, PP->nbV,
				    V_min_i, PP->nbV,
				    width_dirs->n);

	    M = Matrix_Alloc(V_max_n+V_min_n, 1+nvar+1);
	    for (i = 0; i < V_max_n; ++i) {
		value_set_si(M->p[i][0], 1);
		Vector_Oppose(vertex_dirs[V_max_i]->p[i], M->p[i]+1, nvar);
	    }
	    for (i = 0; i < V_min_n; ++i) {
		value_set_si(M->p[V_max_n+i][0], 1);
		Vector_Copy(vertex_dirs[V_min_i]->p[i], M->p[V_max_n+i]+1, nvar);
	    }
	    C = Constraints2Polyhedron(M, options->MaxRays);
	    Matrix_Free(M);
	    n_valid = valid_vertices(C, all_vertices, n_vertices);
	    basis = Cone_Integer_Hull(C, all_vertices, n_valid, options);
	    grow_width_direction_array(width_dirs, basis->NbRows);
	    qsort(width_dirs->wd, width_dirs->n, sizeof(struct width_direction),
		    wd_dir_lex_cmp);
	    sorted_n = width_dirs->n;
	    for (i = 0; i < basis->NbRows; ++i) {
		Vector v;
		struct width_direction wd;

		v.Size = nvar;
		v.p = basis->p[i];
		wd.dir = &v;
		if (bsearch(&wd, width_dirs->wd, sorted_n,
			    sizeof(struct width_direction),
			    wd_dir_lex_cmp))
		    continue;

		n_vertices = add_vertex(all_vertices, n_vertices, basis->p[i]);
		compute_width_direction(V_min->Vertex, V_max->Vertex,
					basis->p[i],
					&width_dirs->wd[width_dirs->n++]);
	    }
	    Matrix_Free(basis);
	    Polyhedron_Free(C);
	}
    }
    Matrix_Free(all_vertices);

    for (i = 0; i < PP->nbV; ++i)
	Matrix_Free(vertex_dirs[i]);
    free(vertex_dirs);

    return width_dirs;
}

/* Computes the lattice width direction of a parametric polytope.
 * The parameter space is allowed to be unbounded.
 * Currently, the parametric polytope and the parameter space
 * are assumed to be full-dimensional.
 *
 * First, we compute the parametric vertices.
 * Then, for each pair of vertices, we construct a (rational) cone
 * of directions for which one vertex attains the minimal value
 * and the other vertex attians the maximal value.
 * The candidate directions are the elements of the integer hulls
 * of these cones.
 * The minimal direction is then obtained by computing the
 * region in the parameter space where each direction yields
 * a smaller (or equal) width than all the other directions.
 *
 * In principle, we can avoid computing candidate directions
 * for vertices with no overlapping activity domains (possibly
 * after opening some facets of the activity domains in the
 * familiar way).
 *
 * The output is a list of triples, consisting of a direction,
 * the corresponding width and the chamber in the parameter
 * space where this direction leads to the minimal width.
 *
 * The algorithm is described in "Integer points in a parameterised
 * polyhedron" by Friedrich Eisenbrand and Gennady Shmonin.
 */
struct width_direction_array *
Polyhedron_Lattice_Width_Directions(Polyhedron *P, Polyhedron *C,
				    struct barvinok_options *options)
{
    Param_Polyhedron *PP;
    unsigned nparam = C->Dimension;
    int i, j, k;
    struct width_direction_array *width_dirs;
    Polyhedron *TC;
    Vector *inner;

    assert(P->NbEq == 0);
    assert(C->NbEq == 0);

    /* Use true context since the algorithm assumes P is non-empty
     * for every point in the context.
     */
    TC = true_context(P, C, options->MaxRays);
    inner = inner_point(TC);

    /* This is overkill, as we discard the computed chambers. */
    PP = Polyhedron2Param_Polyhedron(P, TC, options);

    width_dirs = compute_width_directions(PP, options);
    Param_Polyhedron_Free(PP);

    qsort(width_dirs->wd, width_dirs->n, sizeof(struct width_direction),
	    wd_width_lex_cmp);

    for (i = 1, j = 1; i < width_dirs->n; ++i) {
	/* We could also weed out width_directions that differ by a
	 * (positive) constant from another width_direction, but then
	 * we'd have to put the two width_directions on a common
	 * denominator first.
	 */
	if (Vector_Equal(width_dirs->wd[j-1].width->p,
			 width_dirs->wd[i].width->p, nparam+2))
	    clear_width_direction(&width_dirs->wd[i]);
	else
	    width_dirs->wd[j++] = width_dirs->wd[i];
    }
    width_dirs->n = j;

    for (i = 0, k = 0; i < width_dirs->n; ++i) {
	Matrix *M = Matrix_Alloc(TC->NbConstraints+width_dirs->n-(i-k)-1, nparam+2);
	for (j = 0; j < TC->NbConstraints; ++j)
	    Vector_Copy(TC->Constraint[j], M->p[j], nparam+2);
	for (j = 0; j < width_dirs->n; ++j) {
	    int pos;
	    if (k <= j && j <= i)
		continue;
	    if (j < k)
		pos = TC->NbConstraints + j;
	    else
		pos = TC->NbConstraints + j - (i-k) - 1;
	    Vector_Subtract(width_dirs->wd[j].width->p,
			    width_dirs->wd[j].width->p[nparam+1],
			    width_dirs->wd[i].width->p,
			    width_dirs->wd[i].width->p[nparam+1],
			    M->p[pos]+1, M->p[pos], nparam+1);
	    value_set_si(M->p[pos][0], 1);
	    Vector_Normalize(M->p[pos]+1, nparam+1);
	    if (!is_internal(inner, M->p[pos]))
		value_decrement(M->p[pos][nparam+1], M->p[pos][nparam+1]);
	}
	width_dirs->wd[i].domain = Constraints2Polyhedron(M, options->MaxRays);
	if (emptyQ(width_dirs->wd[i].domain))
	    clear_width_direction(&width_dirs->wd[i]);
	else
	    width_dirs->wd[k++] = width_dirs->wd[i];
	Matrix_Free(M);
    }
    width_dirs->n = k;
    Vector_Free(inner);
    Polyhedron_Free(TC);

    return width_dirs;
}

/* Construct evalue of chambers with their associated widths */
evalue *Polyhedron_Lattice_Width(Polyhedron *P, Polyhedron *C,
				 struct barvinok_options *options)
{
    evalue *width;
    struct evalue_section *s;
    struct width_direction_array *width_dirs;
    int i;
    unsigned nparam = C->Dimension;

    width_dirs = Polyhedron_Lattice_Width_Directions(P, C, options);
    s = ALLOCN(struct evalue_section, width_dirs->n);
    for (i = 0; i < width_dirs->n; ++i) {
	s[i].D = width_dirs->wd[i].domain;
	width_dirs->wd[i].domain = NULL;
	s[i].E = affine2evalue(width_dirs->wd[i].width->p,
			       width_dirs->wd[i].width->p[nparam+1],
			       nparam);
    }
    free_width_direction_array(width_dirs);

    width = evalue_from_section_array(s, i);
    free(s);

    return width;
}
