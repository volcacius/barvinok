#include <barvinok/util.h>
#include <barvinok/options.h>
#include <unistd.h>
#include "topcom.h"
#include "config.h"

#define ALLOC(type) (type*)malloc(sizeof(type))

void run_points2triangs(pid_t *child, int *in, int *out)
{
    int in_fd[2], out_fd[2];

    if (pipe(in_fd))
	assert(0);
    if (pipe(out_fd))
	assert(0);
    *child = fork();
    assert(*child >= 0);
    if (!*child) {
	int rc;

	dup2(in_fd[0], 0);
	dup2(out_fd[1], 1);
	close(in_fd[0]);
	close(out_fd[1]);
	close(in_fd[1]);
	close(out_fd[0]);

	rc = execl(POINTS2TRIANGS_PATH, "points2triangs", "--regular", NULL);
	assert(0);
    }
    close(in_fd[0]);
    close(out_fd[1]);
    *in = in_fd[1];
    *out = out_fd[0];
}

struct vertex {
    Param_Vertices   vertex;
    unsigned	    *facets;
};

struct domain {
    Param_Domain    domain;
    int		    F_len;
};

static struct vertex *construct_vertex(unsigned *vertex_facets, Polyhedron *P,
				       int d, unsigned nparam, unsigned MaxRays)
{
    unsigned nvar = P->Dimension - nparam;
    Matrix *A = Matrix_Alloc(nvar+1, nvar+1);
    Matrix *inv = Matrix_Alloc(nvar+1, nvar+1);
    Matrix *B = Matrix_Alloc(nvar, nparam+2);
    Matrix *V = Matrix_Alloc(nvar, nparam+2);
    Matrix *Domain = Matrix_Alloc(d-nvar, nparam+2);
    Polyhedron *AD;
    unsigned bx;
    int i, j, ix;
    int ok;
    struct vertex *vertex;

    for (j = 0, i = 0, ix = 0, bx = MSB; i < d; ++i) {
	if ((vertex_facets[ix] & bx) == bx) {
	    Vector_Copy(P->Constraint[i]+1, A->p[j], nvar);
	    Vector_Oppose(P->Constraint[i]+1+nvar, B->p[j++], nparam+1);
	}
	NEXT(ix, bx);
    }
    assert(j == nvar);
    value_set_si(A->p[nvar][nvar], 1);
    ok = Matrix_Inverse(A, inv);
    assert(ok);
    Matrix_Free(A);
    inv->NbRows = nvar;
    inv->NbColumns = nvar;
    Matrix_Product(inv, B, V);
    Matrix_Free(B);
    for (i = 0; i < nvar; ++i) {
	value_assign(V->p[i][nparam+1], inv->p[nvar][nvar]);
	Vector_Normalize(V->p[i], V->NbColumns);
    }
    Matrix_Free(inv);
    for (j = 0, i = 0, ix = 0, bx = MSB; i < d; ++i) {
	if ((vertex_facets[ix] & bx) == bx) {
	    NEXT(ix, bx);
	    continue;
	}
	Param_Inner_Product(P->Constraint[i], V, Domain->p[j]);
	if (First_Non_Zero(Domain->p[j]+1, nparam+1) == -1)
	    vertex_facets[ix] |= bx;
	else
	    value_set_si(Domain->p[j++][0], 1);
	NEXT(ix, bx);
    }
    Domain->NbRows = j;
    A = Matrix_Copy(Domain);
    AD = Constraints2Polyhedron(A, MaxRays);
    Matrix_Free(A);
    POL_ENSURE_VERTICES(AD);
    /* A vertex with a lower-dimensional activity domain
     * saturates more facets than derived above and is actually
     * the superimposition of two or more vertices.
     * We therefore discard the domain and (ultimately)
     * the chamber containing it.
     * We keep the vertex, though, since it may reappear
     * in other chambers, which will then likewise be discarded.
     * The same holds if the activity domain is empty.
     */
    if (AD->NbEq > 0) {
	Matrix_Free(Domain);
	Domain = NULL;
    }
    Polyhedron_Free(AD);
    vertex = calloc(1, sizeof(struct vertex));
    vertex->facets = vertex_facets;
    vertex->vertex.Vertex = V;
    vertex->vertex.Domain = Domain;
    return vertex;
}

static int add_vertex_to_domain(struct vertex **vertices, int words,
				unsigned *vertex_facets,
				Polyhedron *P, int d, unsigned nparam,
				struct domain *domain,
				unsigned MaxRays)
{
    struct vertex *vertex;
    unsigned vbx;
    int vi, vix;

    for (vi = 0, vix = 0, vbx = MSB;
	    *vertices;
	    vertices = (struct vertex **)&(*vertices)->vertex.next, ++vi) {
	int i;
	for (i = 0; i < words; ++i)
	    if (((*vertices)->facets[i] & vertex_facets[i]) != vertex_facets[i])
		break;
	if (i >= words) {
	    if (!(*vertices)->vertex.Domain)
		domain->F_len = 0;
	    else
		domain->domain.F[vix] |= vbx;
	    free(vertex_facets);
	    return 0;
	}
	NEXT(vix, vbx);
    }
    if (domain->F_len <= vix) {
	domain->F_len++;
	domain->domain.F = realloc(domain->domain.F,
				   domain->F_len * sizeof(unsigned));
	domain->domain.F[domain->F_len-1] = 0;
    }
    vertex = construct_vertex(vertex_facets, P, d, nparam, MaxRays);
    if (!vertex->vertex.Domain)
	domain->F_len = 0;
    else
	domain->domain.F[vix] |= vbx;
    vertex->vertex.next = &(*vertices)->vertex;
    *vertices = vertex;
    return 1;
}

static int bit_count(unsigned *F, int F_len)
{
    int i;
    int count = 0;

    for (i = 0; i < F_len; ++i) {
	unsigned v = F[i];
	while (v) {
	    v &= v-1;
	    ++count;
	}
    }
    return count;
}

static void compute_domain(struct domain *domain, struct vertex *vertices,
			   Polyhedron *C, unsigned MaxRays)
{
    unsigned bx;
    int i, ix, j;
    int nbV = bit_count(domain->domain.F, domain->F_len);
    unsigned cols = vertices->vertex.Domain->NbColumns;
    unsigned rows = vertices->vertex.Domain->NbRows;
    Matrix *Constraints = Matrix_Alloc(nbV * rows + C->NbConstraints, cols);

    for (i = 0, j = 0, ix = 0, bx = MSB;
	    vertices;
	    vertices = (struct vertex *)vertices->vertex.next, ++i) {
	if ((domain->domain.F[ix] & bx) == bx)
	    Vector_Copy(vertices->vertex.Domain->p[0],
			Constraints->p[(j++)*rows], rows * cols);
	NEXT(ix, bx);
    }
    Vector_Copy(C->Constraint[0], Constraints->p[j*rows], C->NbConstraints * cols);
    domain->domain.Domain = Constraints2Polyhedron(Constraints, MaxRays);
    Matrix_Free(Constraints);
}

static void add_domain(struct domain **domains, struct domain *domain,
		       struct vertex *vertices, Polyhedron *C,
		       struct barvinok_options *options)
{
    options->stats->topcom_chambers++;

    for (; *domains; domains = (struct domain **)&(*domains)->domain.next) {
	int i;
	for (i = 0; i < (*domains)->F_len; ++i)
	    if (((*domains)->domain.F[i] & domain->domain.F[i])
			!= domain->domain.F[i])
		break;
	if (i < (*domains)->F_len)
	    continue;
	for (; i < domain->F_len; ++i)
	    if (domain->domain.F[i])
		break;
	if (i >= domain->F_len) {
	    Param_Domain_Free(&domain->domain);
	    return;
	}
    }
    options->stats->topcom_distinct_chambers++;
    compute_domain(domain, vertices, C, options->MaxRays);
    *domains = domain;
}

#define INT_BITS (sizeof(unsigned) * 8)

/* Remove any empty or lower-dimensional chamber.  The latter
 * lie on the boundary of the context and are facets of other chambers.
 *
 * While we are examining each chamber, also extend the F vector
 * of each chamber to the maximum.
 */
static void remove_empty_chambers(Param_Domain **PD, unsigned vertex_words)
{
    while (*PD) {
	int remove = 0;
	int i;

	if ((*PD)->Domain->NbEq > 0)
	    remove = 1;
	else {
	    POL_ENSURE_FACETS((*PD)->Domain);
	    if ((*PD)->Domain->NbEq > 0)
		remove = 1;
	}
	if (remove) {
	    Param_Domain *D = *PD;
	    *PD = (*PD)->next;
	    D->next = NULL;
	    Param_Domain_Free(D);
	    continue;
	}
	if ((i = ((struct domain*)(*PD))->F_len) < vertex_words)
	    (*PD)->F = realloc((*PD)->F, vertex_words * sizeof(unsigned));
	for (; i < vertex_words; ++i)
	    (*PD)->F[i] = 0;
	PD = &(*PD)->next;
    }
}

/* Clean up memory in struct vertex not in Param_Vertices.
 * We could also remove the vertices we don't need, but
 * then we'd have to fix up the vertex masks in the domains.
 */
static void clean_up_vertices(Param_Vertices *V)
{
    for (; V; V = V->next)
	free(((struct vertex *)V)->facets);
}

static Param_Polyhedron *points2triangs(Matrix *K, Polyhedron *P, Polyhedron *C,
					struct barvinok_options *options)
{
    int in, out;
    int i, j;
    pid_t child;
    FILE *fin, *fout;
    int d = K->NbRows;
    int words = (d+INT_BITS-1)/INT_BITS;
    struct vertex *vertices = NULL;
    struct domain *domains = NULL;
    int vertex_words = 1;
    Param_Polyhedron *PP = ALLOC(Param_Polyhedron);
    unsigned MaxRays = options->MaxRays;

    PP->nbV = 0;
    /* We need the exact facets, because we may make some of them open later */
    POL_UNSET(options->MaxRays, POL_INTEGER);

    run_points2triangs(&child, &in, &out);

    fin = fdopen(in, "w");
    fprintf(fin, "[\n");
    for (i = 0; i < K->NbRows; ++i) {
	fprintf(fin, "[");
	for (j = 0; j < K->NbColumns; ++j)
	    value_print(fin, P_VALUE_FMT, K->p[i][j]);
	fprintf(fin, "]");
    }
    fprintf(fin, "]\n");
    fclose(fin);

    fout = fdopen(out, "r");
    while (fscanf(fout, "T[%d]:={", &i) == 1) {
	struct domain *domain = ALLOC(struct domain);
	memset(domain, 0, sizeof(*domain));
	domain->domain.F = calloc(vertex_words, sizeof(unsigned));
	domain->F_len = vertex_words;

	while (fgetc(fout) == '{') {	/* '{' or closing '}' */
	    int c;
	    unsigned *F = calloc(words, sizeof(unsigned));

	    for (j = 0; j < K->NbColumns; ++j) {
		unsigned v, shift;
		fscanf(fout, "%d", &v);
		shift = INT_BITS - (v % INT_BITS) - 1;
		F[v / INT_BITS] |= 1u << shift;
		fgetc(fout); /* , or } */
	    }
	    if (!domain->F_len)
		free(F);
	    else if (add_vertex_to_domain(&vertices, words, F, P, d, C->Dimension,
					 domain, options->MaxRays))
		++PP->nbV;
	    if ((c = fgetc(fout)) != ',')	/* , or } */
		ungetc(c, fout);
	}
	if (domain->F_len)
	    vertex_words = domain->F_len;
	fgetc(fout); /* ; */
	fgetc(fout); /* \n */
	if (bit_count(domain->domain.F, domain->F_len) > 0)
	    add_domain(&domains, domain, vertices, C, options);
	else {
	    options->stats->topcom_empty_chambers++;
	    Param_Domain_Free(&domain->domain);
	}
    }
    fclose(fout);

    PP->V = &vertices->vertex;
    PP->D = &domains->domain;

    remove_empty_chambers(&PP->D, vertex_words);
    clean_up_vertices(PP->V);

    options->MaxRays = MaxRays;

    return PP;
}

/* Assumes M is of full row rank */
static Matrix *null_space(Matrix *M)
{
    Matrix *H, *Q, *U;
    Matrix *N;
    int i;

    left_hermite(M, &H, &Q, &U);
    N = Matrix_Alloc(M->NbColumns, M->NbColumns - M->NbRows);
    for (i = 0; i < N->NbRows; ++i)
	Vector_Copy(U->p[i] + M->NbRows, N->p[i], N->NbColumns);
    Matrix_Free(H);
    Matrix_Free(Q);
    Matrix_Free(U);
    return N;
}

/*
 * left_hermite may leave positive entries below the main diagonal in H.
 * This function postprocesses the output of left_hermite to make
 * the non-zero entries below the main diagonal negative.
 */
static void neg_left_hermite(Matrix *A, Matrix **H_p, Matrix **Q_p, Matrix **U_p)
{
    int row, col, i, j;
    Matrix *H, *U, *Q;

    left_hermite(A, &H, &Q, &U);
    *H_p = H;
    *Q_p = Q;
    *U_p = U;

    for (row = 0, col = 0; col < H->NbColumns; ++col, ++row) {
	while (value_zero_p(H->p[row][col]))
	    ++row;
	for (i = 0; i < col; ++i) {
	    if (value_negz_p(H->p[row][i]))
		continue;

	    /* subtract column col from column i in H and U */
	    for (j = 0; j < H->NbRows; ++j)
		value_subtract(H->p[j][i], H->p[j][i], H->p[j][col]);
	    for (j = 0; j < U->NbRows; ++j)
		value_subtract(U->p[j][i], U->p[j][i], U->p[j][col]);

	    /* add row i to row col in Q */
	    for (j = 0; j < Q->NbColumns; ++j)
		value_addto(Q->p[col][j], Q->p[col][j], Q->p[i][j]);
	}
    }
}

static void SwapColumns(Value **V, int n, int i, int j)
{
    int r;

    for (r = 0; r < n; ++r)
	value_swap(V[r][i], V[r][j]);
}

/* C is assumed to be the "true" context, i.e., it has been intersected
 * with the projection of P onto the parameter space.
 * Furthermore, P and C are assumed to be full-dimensional.
 */
Param_Polyhedron *TC_P2PP(Polyhedron *P, Polyhedron *C,
			  struct barvinok_options *options)
{
    unsigned nparam = C->Dimension;
    unsigned nvar = P->Dimension - C->Dimension;
    Matrix *M;
    int i, j, d;
    Matrix *H, *U, *Q;
    Matrix *A;
    Matrix *K;
    int rows;
    Param_Polyhedron *PP;

    assert(P->NbEq == 0);
    assert(C->NbEq == 0);
    /* move constraints only involving parameters down
     * and move unit vectors (if there are any) to the right place.
     */
    for (d = 0, j = P->NbConstraints; d < j; ++d) {
	int index;
	index = First_Non_Zero(P->Constraint[d]+1, nvar);
	if (index != -1) {
	    if (index != d &&
		(value_one_p(P->Constraint[d][1+index]) ||
		 value_mone_p(P->Constraint[d][1+index])) &&
		First_Non_Zero(P->Constraint[d]+1+index+1, nvar-(index+1)) == -1) {
		Vector_Exchange(P->Constraint[d], P->Constraint[index],
				P->Dimension+2);
		--d;
	    }
	    continue;
	}
	while (d < --j && First_Non_Zero(P->Constraint[j]+1, nvar) == -1)
	    ;
	if (d >= j)
	    break;
	Vector_Exchange(P->Constraint[d], P->Constraint[j], P->Dimension+2);
    }
    M = Matrix_Alloc(d+nvar, nvar);
    for (j = 0; j < d; ++j)
	Vector_Copy(P->Constraint[j]+1, M->p[j], nvar);

    neg_left_hermite(M, &H, &Q, &U);
    Matrix_Free(M);
    Matrix_Free(Q);
    Matrix_Free(U);
    H->NbRows -= nvar;

    /* Rearrange rows such that top of H is lower diagonal and
     * add extra unit vector rows that are positive linear
     * combinations of two or more of the top rows.
     */
    for (i = 0; i < H->NbColumns; ++i) {
	for (j = i; j < H->NbRows; ++j)
	    if (value_notzero_p(H->p[j][i]))
		break;
	if (j != i) {
	    Vector_Exchange(P->Constraint[i], P->Constraint[j], P->Dimension+2);
	    Vector_Exchange(H->p[i], H->p[j], H->NbColumns);
	}
	if (First_Non_Zero(H->p[i], i) != -1)
	    value_set_si(H->p[H->NbRows++][i], 1);
    }

    rows = H->NbRows-nvar;
    A = Matrix_Alloc(rows, nvar+rows);
    for (i = nvar; i < d; ++i) {
	Vector_Oppose(H->p[i], A->p[i-nvar], H->NbColumns);
	value_set_si(A->p[i-nvar][i], 1);
    }
    for (; i < H->NbRows; ++i) {
	j = First_Non_Zero(H->p[i], nvar);
	assert(j != -1);
	Vector_Oppose(H->p[j], A->p[i-nvar], H->NbColumns);
	value_set_si(A->p[i-nvar][i], 1);
	SwapColumns(A->p, A->NbRows, i, j);
    }
    Matrix_Free(H);
    K = null_space(A);
    Matrix_Free(A);
    /* Ignore extra constraints */
    K->NbRows = d;
    PP = points2triangs(K, P, C, options);
    Matrix_Free(K);
    return PP;
}
