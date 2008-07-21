#include <assert.h>
#include <iostream>
#include <vector>
#include <barvinok/options.h>
#include "bernoulli.h"
#include "binomial.h"
#include "conversion.h"
#include "decomposer.h"
#include "lattice_point.h"
#include "laurent.h"
#include "param_util.h"
#include "power.h"
#include "reduce_domain.h"
#include "config.h"

using std::cerr;
using std::ostream;
using std::endl;
using std::vector;

#ifdef HAVE_GNUCXX_HASHMAP

#include <ext/hash_map>

#define HASH_MAP __gnu_cxx::hash_map

namespace __gnu_cxx
{
        template<> struct hash< const std::vector<int> >
        {
                size_t operator()( const std::vector<int>& x ) const
                {
			unsigned long __h = 0;
			for (int i = 0; i < x.size(); ++i)
			    __h = 5 * __h + x[i];
                        return size_t(__h);
                }
        };
}

#else

#warning "no hash_map available"
#include <map>
#define HASH_MAP std::map

#endif

#define ALLOC(type) (type*)malloc(sizeof(type))
#define ALLOCN(type,n) (type*)malloc((n) * sizeof(type))

static ostream & operator<< (ostream & os, const vector<int> & v)
{
    os << "[";
    for (int i = 0; i < v.size(); ++i) {
        if (i)
            os << ", ";
        os << v[i];
    }
    os << "]";
    return os;
}

/* Compares the position of the first non-zero element
 * in both vectors and the second non-zero element
 * if the position of the first non-zero element is the same.
 */
static int pos_cmp(const void *a, const void *b)
{
    int pa, pb;
    const Vector *va = *(const Vector **)a;
    const Vector *vb = *(const Vector **)b;

    pa = First_Non_Zero(va->p, va->Size);
    pb = First_Non_Zero(vb->p, vb->Size);

    if (pa == pb) {
	pa = First_Non_Zero(va->p+pa+1, va->Size-pa-1);
	pb = First_Non_Zero(vb->p+pa+1, vb->Size-pa-1);
    }

    return pa - pb;
}

/* Represents the vertex and the rays of a vertex cone */
struct vertex_cone {
    unsigned dim;
    /* The coefficients of the rays, ordered according
     * to the first non-zero coefficient.
     */
    Vector **coeff;
    Matrix Rays;

    /* The powers of the coefficients of the rays */
    struct power ***coeff_power;

    /* The coordinates of the integer point in the fundamental
     * parallelepiped, in the basis formed by the rays of
     * the vertex cone.
     */
    evalue **E_vertex;
    unsigned max_power;
    /* Bernoulli polynomials corresponding to each E_vertex */
    evalue ***bernoulli_t;

    vertex_cone(unsigned dim);
    void init(const mat_ZZ &rays, Param_Vertices *V, unsigned max_power);
    void clear();

    const evalue *get_bernoulli(int n, int i);

    ~vertex_cone() {
	for (int i = 0; i < dim; ++i)
	    Vector_Free(coeff[i]);
	free(coeff);
	delete [] E_vertex;
	free(Rays.p);
	for (int i = 0; i < dim; ++i)
	    delete [] coeff_power[i];
	delete [] coeff_power;
	delete [] bernoulli_t;
    }
};

vertex_cone::vertex_cone(unsigned dim) : dim(dim)
{
    E_vertex = new evalue *[dim];
    bernoulli_t = new evalue **[dim];

    coeff = ALLOCN(Vector *, dim);
    for (int i = 0; i < dim; ++i)
	coeff[i] = Vector_Alloc(dim);

    Rays.NbRows = Rays.NbColumns = dim;
    Rays.p = ALLOCN(Value *, dim);

    coeff_power = new struct power **[dim];
    for (int i = 0; i < dim; ++i)
	coeff_power[i] = new struct power *[dim];
}

void vertex_cone::init(const mat_ZZ &rays, Param_Vertices *V,
			unsigned max_power)
{
    unsigned nparam = V->Vertex->NbColumns - 2;
    this->max_power = max_power;

    for (int i = 0; i < dim; ++i)
	zz2values(rays[i], coeff[i]->p);
    qsort(coeff, dim, sizeof(Vector *), pos_cmp);

    for (int i = 0; i < dim; ++i) {
	for (int j = 0; j < dim; ++j) {
	    if (value_notzero_p(coeff[i]->p[j]))
		coeff_power[i][j] = new struct power(coeff[i]->p[j], max_power);
	    else
		coeff_power[i][j] = NULL;
	}
    }

    for (int i = 0; i < dim; ++i)
	Rays.p[i] = coeff[i]->p;
    Matrix *T = Transpose(&Rays);
    Matrix *L = relative_coordinates(V, T);
    Matrix_Free(T);

    for (int i = 0; i < dim; ++i)
	E_vertex[i] = ceiling(L->p[i], V->Vertex->p[0][nparam+1], nparam, NULL);
    Matrix_Free(L);

    for (int j = 0; j < dim; ++j) {
	bernoulli_t[j] = new evalue *[max_power];
	for (int k = 0; k < max_power; ++k)
	    bernoulli_t[j][k] = NULL;
    }
}

/*
 * Returns b(n, E_vertex[i])
 */
const evalue *vertex_cone::get_bernoulli(int n, int i)
{
    assert(n > 0);
    if (!bernoulli_t[i][n-1]) {
	struct poly_list *bernoulli = bernoulli_compute(n);
	bernoulli_t[i][n-1] = evalue_polynomial(bernoulli->poly[n],
						E_vertex[i]);
    }
    return bernoulli_t[i][n-1];
}

void vertex_cone::clear()
{
    for (int i = 0; i < dim; ++i)
	if (E_vertex[i])
	    evalue_free(E_vertex[i]);

    for (int i = 0; i < dim; ++i) {
	for (int j = 0; j < max_power; ++j) {
	    if (bernoulli_t[i][j])
		evalue_free(bernoulli_t[i][j]);
	}
	delete [] bernoulli_t[i];
    }

    for (int i = 0; i < dim; ++i) {
	for (int j = 0; j < dim; ++j)
	    if (coeff_power[i][j])
		delete coeff_power[i][j];
    }
}

struct todd_product {
    vertex_cone &vc;
    unsigned dim;
    /* The non-zero coefficients in the rays of the vertex cone */
    vector< vector<int> > mask;
    /* For each ray, the power of each variable it contributes */
    vector< vector<int> > selected;
    /* The powers of each variable that still need to be selected */
    vector<int> left;
    /* For each variable, the last ray that has a non-zero coefficient */
    vector<int> last_level;

    HASH_MAP<const vector<int>, const evalue *> cache;

    todd_product(vertex_cone &vc);
    evalue *add(evalue *sum);
    const evalue *get_coefficient(const vector<int> &powers);

    ~todd_product() {
	HASH_MAP<const vector<int>, const evalue *>::iterator j;
	for (j = cache.begin(); j != cache.end(); ++j)
	    if ((*j).second)
		evalue_free(const_cast<evalue *>((*j).second));
    }
};

todd_product::todd_product(vertex_cone &vc) : vc(vc)
{
    dim = vc.dim;
    for (int i = 0; i < dim; ++i) {
	mask.push_back(vector<int>(dim));
	selected.push_back(vector<int>(dim));
    }
    last_level = vector<int>(dim);
    for (int i = 0; i < dim; ++i) {
	for (int j = 0; j < dim; ++j) {
	    if (vc.coeff_power[i][j]) {
		mask[i][j] = 1;
		last_level[j] = i;
	    }
	}
    }
}

void multinomial(const vector<int> &k, Value *m)
{
    int s = 0;
    value_set_si(*m, 1);
    for (int i = 0; i < k.size(); ++i) {
	if (k[i] == 0)
	    continue;
	s += k[i];
	value_multiply(*m, *m, *binomial(s, k[i]));
    }
}

/* Add coefficient of selected powers of variables to sum
 * and return the result.
 * The contribution of each ray j of the vertex cone is
 *
 *                       (      \sum k      )
 *   b(\sum k, \ceil{v}) ( k_1, \ldots, k_n ) c_1^{k_1} \cdots c_n^{k_n}
 *
 *   with k_i the selected powers, c_i the coefficients of the ray
 *   and \ceil{v} the coordinate E_vertex[j] corresponding to the ray.
 */
evalue *todd_product::add(evalue *sum)
{
    evalue *c = NULL;
    for (int i = 0; i < dim; ++i) {
	int s = 0;
	evalue *f = ALLOC(evalue);
	value_init(f->d);
	evalue_set_si(f, 1, 1);
	for (int j = 0; j < dim; ++j) {
	    if (!selected[i][j])
		continue;
	    value_multiply(f->x.n, f->x.n,
				*(*vc.coeff_power[i][j])[selected[i][j]]);
	    value_multiply(f->d, f->d, *factorial(selected[i][j]));
	    s += selected[i][j];
	}
	if (s > 0)
	    emul(vc.get_bernoulli(s, i), f);
	if (!c)
	    c = f;
	else {
	    emul(f, c);
	    evalue_free(f);
	}
    }
    if (!sum)
	sum = c;
    else {
	eadd(c, sum);
	evalue_free(c);
    }
    return sum;
}

static int first_non_zero(const vector<int>& row)
{
    for (int i = 0; i < row.size(); ++i)
	if (row[i] != 0)
	    return i;
    return -1;
}

/* Return coefficient of the term with exponents powers by
 * iterating over all combinations of exponents for each ray
 * that sum up to the given exponents.
 */
const evalue *todd_product::get_coefficient(const vector<int> &powers)
{
    evalue *c = NULL;

    HASH_MAP<const vector<int>, const evalue *>::iterator i;
    i = cache.find(powers);
    if (i != cache.end())
	return (*i).second;

    for (int i = 0; i < vc.dim; ++i)
	for (int j = 0; j < vc.dim; ++j)
	    selected[i][j] = 0;

    left = powers;
    int nz = first_non_zero(left);
    int level = last_level[nz];
    int p = nz;
    while (level >= 0) {
	if (mask[level][p] && left[p] > 0) {
	    selected[level][p]++;
	    left[p]--;
	    /* Fill out remaining powers and make sure we backtrack from
	     * the right position.
	     */
	    for (int i = 0; i < vc.dim; ++i) {
		if (left[i] == 0)
		    continue;
		selected[last_level[i]][i] += left[i];
		left[i] = 0;
		if (last_level[i] >= level) {
		    level = last_level[i];
		    p = i;
		}
	    }

	    c = add(c);
	}
	if (selected[level][p]) {
	    left[p] += selected[level][p];
	    selected[level][p] = 0;
	}
	if (--p < 0) {
	    --level;
	    p = dim-1;
	}
    }
    cache[powers] = c;
    return c;
}

/* Maintains the coefficients of the reciprocals of the
 * (negated) rays of the vertex cone vc.
 */
struct reciprocal {
    vertex_cone &vc;

    /* can_borrow_from[i][j] = 1 if there is a ray
     * with first non-zero coefficient i and a subsequent
     * non-zero coefficient j.
     */
    vector< vector<int> > can_borrow_from;
    /* Total exponent that a variable can borrow from subsequent vars */
    vector<int> can_borrow;
    /* has_borrowed_from[i][j] is the exponent borrowed by i from j */
    vector< vector<int> > has_borrowed_from;
    /* Total exponent borrowed from subsequent vars */
    vector<int> has_borrowed;
    /* The last variable that can borrow from subsequent vars */
    int last;

    /* Position of the first non-zero coefficient in each ray. */
    vector<int> neg;

    /* Power without any "borrowing" */
    vector<int> base_power;
    /* Power after "borrowing" */
    vector<int> power;

    /* The non-zero coefficients in the rays of the vertex cone,
     * except the first.
     */
    vector< vector<int> > mask;
    /* For each ray, the (positive) power of each variable it contributes */
    vector< vector<int> > selected;
    /* The powers of each variable that still need to be selected */
    vector<int> left;

    HASH_MAP<const vector<int>, const evalue *> cache;

    reciprocal(vertex_cone &vc);
    void start(vector<int> &power);
    int next();

    evalue *add(evalue *sum);
    const evalue *get_coefficient();
    ~reciprocal() {
	HASH_MAP<const vector<int>, const evalue *>::iterator j;
	for (j = cache.begin(); j != cache.end(); ++j)
	    if ((*j).second)
		evalue_free(const_cast<evalue *>((*j).second));
    }
};

reciprocal::reciprocal(vertex_cone &vc) : vc(vc)
{
    for (int i = 0; i < vc.dim; ++i) {
	can_borrow_from.push_back(vector<int>(vc.dim));
	has_borrowed_from.push_back(vector<int>(vc.dim));
	mask.push_back(vector<int>(vc.dim));
	selected.push_back(vector<int>(vc.dim));
    }
    can_borrow = vector<int>(vc.dim);
    has_borrowed = vector<int>(vc.dim);
    neg = vector<int>(vc.dim);
    left = vector<int>(vc.dim);

    for (int i = 0; i < vc.dim; ++i) {
	int pos = First_Non_Zero(vc.coeff[i]->p, vc.coeff[i]->Size);
	neg[i] = pos;
	for (int j = pos+1; j < vc.dim; ++j)
	    if (value_notzero_p(vc.coeff[i]->p[j])) {
		mask[i][j] = 1;
		can_borrow_from[neg[i]][j] = 1;
	    }
    }
}

/* Initialize power to the exponent of the todd product
 * required to compute the coefficient in the full product
 * of the term with exponent req_power, without any
 * "borrowing".
 */
void reciprocal::start(vector<int> &req_power)
{
    power = req_power;
    for (int j = 0; j < vc.dim; ++j)
	power[neg[j]]++;

    base_power = power;

    last = -1;
    for (int i = 0; i < vc.dim; ++i) {
	can_borrow[i] = 0;
	has_borrowed[i] = 0;
	for (int j = i+1; j < vc.dim; ++j) {
	    has_borrowed_from[i][j] = 0;
	    if (can_borrow_from[i][j])
		can_borrow[i] += power[j];
	}
	if (can_borrow[i])
	    last = i;
    }
}

/* Set power to the next exponent of the todd product required
 * and return 1 as long as there is any such exponent left.
 */
int reciprocal::next()
{
    int p = last;
    while (p >= 0) {
	if (has_borrowed[p] < can_borrow[p]) {
	    int j;
	    for (j = p+1; j < vc.dim; ++j)
		if (can_borrow_from[p][j]) {
		    if (power[j] > 0)
			break;
		    else if (has_borrowed_from[p][j]) {
			power[j] += has_borrowed_from[p][j];
			power[p] -= has_borrowed_from[p][j];
			has_borrowed[p] -= has_borrowed_from[p][j];
			has_borrowed_from[p][j] = 0;
		    }
		}
	    if (j < vc.dim) {
		has_borrowed_from[p][j]++;
		has_borrowed[p]++;
		power[p]++;
		power[j]--;
		return 1;
	    }
	}
	if (has_borrowed[p]) {
	    for (int j = p+1; j < vc.dim; ++j)
		if (has_borrowed_from[p][j]) {
		    power[j] += has_borrowed_from[p][j];
		    has_borrowed_from[p][j] = 0;
		}
	    power[p] -= has_borrowed[p];
	    has_borrowed[p] = 0;
	}
	--p;
    }
    return 0;
}

/* Add coefficient of selected powers of variables to sum
 * and return the result.
 * The contribution of each ray j of the vertex cone is
 *
 *  (          K           )
 *  ( k_{f+1}, \ldots, k_n ) (-1)^{K+1}
 *  			c_f^{-K-1} c_{f+1}^{k_{f+1}} \cdots c_n^{k_n}
 *
 * K = \sum_{i=f+1}^n k_i
 */
evalue *reciprocal::add(evalue *sum)
{
    evalue *t = NULL;
    for (int i = 0; i < vc.dim; ++i) {
	evalue *c = ALLOC(evalue);
	value_init(c->d);
	value_init(c->x.n);
	value_set_si(c->d, 1);
	multinomial(selected[i], &c->x.n);
	int s = 0;
	for (int j = 0; j < vc.dim; ++j) {
	    if (selected[i][j] == 0)
		continue;
	    value_multiply(c->x.n, c->x.n,
				*(*vc.coeff_power[i][j])[selected[i][j]]);
	    s += selected[i][j];
	}
	value_multiply(c->d, c->d, *(*vc.coeff_power[i][neg[i]])[s+1]);
	if (!(s % 2))
	    value_oppose(c->x.n, c->x.n);
	if (value_neg_p(c->d)) {
	    value_oppose(c->d, c->d);
	    value_oppose(c->x.n, c->x.n);
	}
	if (!t)
	    t = c;
	else {
	    emul(c, t);
	    evalue_free(c);
	}
    }
    if (!sum)
	sum = t;
    else {
	eadd(t, sum);
	evalue_free(t);
    }
    return sum;
}

/* Return coefficient of the term with exponents powers by
 * iterating over all combinations of exponents for each ray
 * that sum up to the given exponents.
 */
const evalue *reciprocal::get_coefficient()
{
    evalue *c = NULL;

    for (int j = 0; j < vc.dim; ++j)
	left[j] = base_power[j] - power[j];

    HASH_MAP<const vector<int>, const evalue *>::iterator i;
    i = cache.find(left);
    if (i != cache.end())
	return (*i).second;

    int nz = first_non_zero(left);
    if (nz == -1)
	return cache[power] = add(c);
    if (left[nz] > 0)
	return NULL;

    for (int i = 0; i < vc.dim; ++i)
	for (int j = 0; j < vc.dim; ++j)
	    selected[i][j] = 0;

    int level = vc.dim-1;
    int p = vc.dim-1;
    while (level >= 0) {
	int nz = first_non_zero(left);
	if (nz < neg[level] || (nz == neg[level] && left[nz] > 0)) {
	    assert(p == vc.dim-1);
	    --level;
	    continue;
	}
	if (nz == neg[level] && mask[level][p]) {
	    selected[level][p]++;
	    left[p]--;
	    left[neg[level]]++;
	    int nz = first_non_zero(left);
	    if (nz == -1)
		c = add(c);
	    else if (left[nz] < 0) {
		level = vc.dim-1;
		p = vc.dim-1;
		continue;
	    }
	}
	if (selected[level][p]) {
	    left[p] += selected[level][p];
	    left[neg[level]] -= selected[level][p];
	    selected[level][p] = 0;
	}
	if (--p < 0) {
	    --level;
	    p = vc.dim-1;
	}
    }
    cache[left] = c;

    return c;
}

/* A term in the input polynomial */
struct term {
    vector<int> powers;
    const evalue *coeff;
};

static void collect_terms_r(vector<struct term> &terms, vector<int> &powers,
			    const evalue *polynomial, unsigned nvar)
{
    if (EVALUE_IS_ZERO(*polynomial))
	return;

    if (value_zero_p(polynomial->d))
	assert(polynomial->x.p->type == ::polynomial);
    if (value_notzero_p(polynomial->d) || polynomial->x.p->pos > nvar) {
	struct term t;
	t.powers = powers;
	t.coeff = polynomial;
	terms.push_back(t);
	return;
    }

    for (int i = polynomial->x.p->size-1; i >= 0; --i) {
	powers[polynomial->x.p->pos-1] = i;
	collect_terms_r(terms, powers, &polynomial->x.p->arr[i], nvar);
    }
}

/* Expand "polynomial" as a sum of powers of the "nvar" variables,
 * collect the terms in "terms" and return the maximal total degree.
 */
static unsigned collect_terms(vector<struct term> &terms,
			  const evalue *polynomial, unsigned nvar)
{
    vector<int> powers(nvar);
    for (int i = 0; i < nvar; ++i)
	powers[i] = 0;
    collect_terms_r(terms, powers, polynomial, nvar);

    unsigned max_degree = 0;
    for (int i = 0; i < terms.size(); ++i) {
	int sum = 0;
	for (int j = 0; j < nvar; ++j)
	    sum += terms[i].powers[j];
	if (sum > max_degree)
	    max_degree = sum;
    }
    return max_degree;
}

/*
static void dump_terms(const vector<struct term> &terms)
{
    for (int i = 0; i < terms.size(); ++i) {
	cerr << terms[i].powers;
	print_evalue(stderr, terms[i].coeff, test);
    }
}
*/

struct laurent_summator : public signed_cone_consumer,
			  public vertex_decomposer {
    const evalue *polynomial;
    unsigned dim;
    vertex_cone vc;
    vector<struct term> terms;
    evalue *result;
    unsigned max_power;

    laurent_summator(const evalue *e, unsigned dim, Param_Polyhedron *PP) :
			polynomial(e), dim(dim), vertex_decomposer(PP, *this),
			vc(dim) {
	max_power = dim + collect_terms(terms, polynomial, dim);
	result = NULL;
    }
    ~laurent_summator() {
	if (result)
	    evalue_free(result);
    }
    virtual void handle(const signed_cone& sc, barvinok_options *options);
};

static int first_non_zero(const vec_ZZ& row)
{
    for (int i = 0; i < row.length(); ++i)
	if (row[i] != 0)
	    return i;
    return -1;
}

void laurent_summator::handle(const signed_cone& sc, barvinok_options *options)
{
    assert(sc.det == 1);

    vc.init(sc.rays, V, max_power);
    reciprocal recip(vc);
    todd_product tp(vc);
    for (int i = 0; i < terms.size(); ++i) {
	recip.start(terms[i].powers);
	do {
	    const evalue *c = recip.get_coefficient();
	    if (!c)
		continue;

	    const evalue *t = tp.get_coefficient(recip.power);

	    evalue *f = evalue_dup(terms[i].coeff);
	    if (sc.sign < 0)
		evalue_negate(f);
	    for (int j = 0; j < dim; ++j)
		evalue_mul(f, *factorial(terms[i].powers[j]));
	    evalue_shift_variables(f, 0, -dim);
	    emul(c, f);
	    emul(t, f);
	    if (!result)
		result = f;
	    else {
		eadd(f, result);
		evalue_free(f);
	    }
	} while (recip.next());
    }
    vc.clear();
}

evalue *laurent_summate(Polyhedron *P, evalue *e, unsigned nvar,
				   struct barvinok_options *options)
{
    Polyhedron *U, *TC;
    Param_Polyhedron *PP;
    struct evalue_section *s;
    int nd = -1;
    Param_Domain *PD;
    evalue *sum;

    U = Universe_Polyhedron(P->Dimension - nvar);
    PP = Polyhedron2Param_Polyhedron(P, U, options);

    for (nd = 0, PD = PP->D; PD; ++nd, PD = PD->next);
    s = ALLOCN(struct evalue_section, nd);

    TC = true_context(P, U, options->MaxRays);
    FORALL_REDUCED_DOMAIN(PP, TC, nd, options, i, PD, rVD)
	Param_Vertices *V;
	laurent_summator ls(e, nvar, PP);

	FORALL_PVertex_in_ParamPolyhedron(V, PD, PP) // _i internal counter
	    ls.decompose_at_vertex(V, _i, options);
	END_FORALL_PVertex_in_ParamPolyhedron;

	s[i].D = rVD;
	s[i].E = ls.result;
	ls.result = NULL;
    END_FORALL_REDUCED_DOMAIN
    Polyhedron_Free(TC);
    Polyhedron_Free(U);
    Param_Polyhedron_Free(PP);

    sum = evalue_from_section_array(s, nd);
    free(s);

    return sum;
}
