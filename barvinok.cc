#include <assert.h>
#include <iostream>
#include <vector>
#include <deque>
#include <string>
#include <sstream>
#include <gmp.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>
#include "ehrhartpolynom.h"
#include <barvinok.h>
#include <util.h>
extern "C" {
#include <polylib/polylibgmp.h>
#include "ev_operations.h"
}

#ifdef NTL_STD_CXX
using namespace NTL;
#endif
using std::cout;
using std::endl;
using std::vector;
using std::deque;
using std::string;
using std::ostringstream;

#define ALLOC(p) (((long *) (p))[0])
#define SIZE(p) (((long *) (p))[1])
#define DATA(p) ((mp_limb_t *) (((long *) (p)) + 2))

static void value2zz(Value v, ZZ& z)
{
    int sa = v[0]._mp_size;
    int abs_sa = sa < 0 ? -sa : sa;

    _ntl_gsetlength(&z.rep, abs_sa);
    mp_limb_t * adata = DATA(z.rep);
    for (int i = 0; i < abs_sa; ++i)
	adata[i] = v[0]._mp_d[i];
    SIZE(z.rep) = sa;
}

static void zz2value(ZZ& z, Value& v)
{
    if (!z.rep) {
	value_set_si(v, 0);
	return;
    }

    int sa = SIZE(z.rep);
    int abs_sa = sa < 0 ? -sa : sa;

    mp_limb_t * adata = DATA(z.rep);
    mpz_realloc2(v, __GMP_BITS_PER_MP_LIMB * abs_sa);
    for (int i = 0; i < abs_sa; ++i)
	v[0]._mp_d[i] = adata[i];
    v[0]._mp_size = sa;
}

/*
 * We just ignore the last column and row
 * If the final element is not equal to one
 * then the result will actually be a multiple of the input
 */
static void matrix2zz(Matrix *M, mat_ZZ& m, unsigned nr, unsigned nc)
{
    m.SetDims(nr, nc);

    for (int i = 0; i < nr; ++i) {
//	assert(value_one_p(M->p[i][M->NbColumns - 1]));
	for (int j = 0; j < nc; ++j) {
	    value2zz(M->p[i][j], m[i][j]);
	}
    }
}

static void values2zz(Value *p, vec_ZZ& v, int len)
{
    v.SetLength(len);

    for (int i = 0; i < len; ++i) {
	value2zz(p[i], v[i]);
    }
}

/*
 * We add a 0 at the end, because we need it afterwards
 */
static Vector * zz2vector(vec_ZZ& v)
{
    Vector *vec = Vector_Alloc(v.length()+1);
    assert(vec);
    for (int i = 0; i < v.length(); ++i)
	zz2value(v[i], vec->p[i]);

    value_set_si(vec->p[v.length()], 0);

    return vec;
}

static void rays(mat_ZZ& r, Polyhedron *C)
{
    unsigned dim = C->NbRays - 1; /* don't count zero vertex */
    assert(C->NbRays - 1 == C->Dimension);
    r.SetDims(dim, dim);
    ZZ tmp;

    int i, c;
    for (i = 0, c = 0; i < dim; ++i)
	if (value_zero_p(C->Ray[i][dim+1])) {
	    for (int j = 0; j < dim; ++j) {
		value2zz(C->Ray[i][j+1], tmp);
		r[j][c] = tmp;
	    }
	    ++c;
	}
}

static Matrix * rays(Polyhedron *C)
{
    unsigned dim = C->NbRays - 1; /* don't count zero vertex */
    assert(C->NbRays - 1 == C->Dimension);

    Matrix *M = Matrix_Alloc(dim+1, dim+1);
    assert(M);

    int i, c;
    for (i = 0, c = 0; i <= dim && c < dim; ++i)
	if (value_zero_p(C->Ray[i][dim+1])) {
	    Vector_Copy(C->Ray[i] + 1, M->p[c], dim);
	    value_set_si(M->p[c++][dim], 0);
	}
    assert(c == dim);
    value_set_si(M->p[dim][dim], 1);

    return M;
}

/*
 * Returns the largest absolute value in the vector
 */
static ZZ max(vec_ZZ& v)
{
    ZZ max = abs(v[0]);
    for (int i = 1; i < v.length(); ++i)
	if (abs(v[i]) > max)
	    max = abs(v[i]);
    return max;
}

class cone {
public:
    cone(Matrix *M) {
	Cone = 0;
	Rays = Matrix_Copy(M);
	set_det();
    }
    cone(Polyhedron *C) {
	Cone = Polyhedron_Copy(C);
	Rays = rays(C);
	set_det();
    }
    void set_det() {
	mat_ZZ A;
	matrix2zz(Rays, A, Rays->NbRows - 1, Rays->NbColumns - 1);
	det = determinant(A);
	Value v;
	value_init(v);
	zz2value(det, v);
	value_clear(v);
    }

    Vector* short_vector(vec_ZZ& lambda) {
	Matrix *M = Matrix_Copy(Rays);
	Matrix *inv = Matrix_Alloc(M->NbRows, M->NbColumns);
	int ok = Matrix_Inverse(M, inv);
	assert(ok);
	Matrix_Free(M);

	ZZ det2;
	mat_ZZ B;
	mat_ZZ U;
	matrix2zz(inv, B, inv->NbRows - 1, inv->NbColumns - 1);
	long r = LLL(det2, B, U);

	ZZ min = max(B[0]);
	int index = 0;
	for (int i = 1; i < B.NumRows(); ++i) {
	    ZZ tmp = max(B[i]);
	    if (tmp < min) {
		min = tmp;
		index = i;
	    }
	}

	Matrix_Free(inv);

	lambda = B[index];

	Vector *z = zz2vector(U[index]);
	Value tmp;
	value_init(tmp);
	Polyhedron *C = poly();
	int i;
	for (i = 0; i < C->NbConstraints; ++i) {
	    Inner_Product(z->p, C->Constraint[i]+1, z->Size-1, &tmp);
	    if (value_pos_p(tmp))
		break;
	}
	if (i == C->NbConstraints) {
	    value_set_si(tmp, -1);
	    Vector_Scale(z->p, z->p, tmp, z->Size-1);
	}
	value_clear(tmp);
	return z;
    }

    ~cone() {
	Polyhedron_Free(Cone);
	Matrix_Free(Rays);
    }

    Polyhedron *poly() {
	if (!Cone) {
	    Matrix *M = Matrix_Alloc(Rays->NbRows+1, Rays->NbColumns+1);
	    for (int i = 0; i < Rays->NbRows; ++i) {
		Vector_Copy(Rays->p[i], M->p[i]+1, Rays->NbColumns);
		value_set_si(M->p[i][0], 1);
	    }
	    Vector_Set(M->p[Rays->NbRows]+1, 0, Rays->NbColumns-1);
	    value_set_si(M->p[Rays->NbRows][0], 1);
	    value_set_si(M->p[Rays->NbRows][Rays->NbColumns], 1);
	    Cone = Rays2Polyhedron(M, M->NbRows+1);
	    assert(Cone->NbConstraints == Cone->NbRays);
	    Matrix_Free(M);
	}
	return Cone;
    }

    ZZ det;
    Polyhedron *Cone;
    Matrix *Rays;
};

class dpoly {
public:
    vec_ZZ coeff;
    dpoly(int d, ZZ& degree, int offset = 0) {
	coeff.SetLength(d+1);

	int min = d + offset;
	if (degree < ZZ(INIT_VAL, min))
	    min = to_int(degree);

	ZZ c = ZZ(INIT_VAL, 1);
	if (!offset)
	    coeff[0] = c;
	for (int i = 1; i <= min; ++i) {
	    c *= (degree -i + 1);
	    c /= i;
	    coeff[i-offset] = c;
	}
    }
    void operator *= (dpoly& f) {
	assert(coeff.length() == f.coeff.length());
	vec_ZZ old = coeff;
	coeff = f.coeff[0] * coeff;
	for (int i = 1; i < coeff.length(); ++i)
	    for (int j = 0; i+j < coeff.length(); ++j)
		coeff[i+j] += f.coeff[i] * old[j];
    }
    void div(dpoly& d, mpq_t count, ZZ& sign) {
	int len = coeff.length();
	Value tmp;
	value_init(tmp);
	mpq_t* c = new mpq_t[coeff.length()];
	mpq_t qtmp;
	mpq_init(qtmp);
	for (int i = 0; i < len; ++i) {
	    mpq_init(c[i]);
	    zz2value(coeff[i], tmp);
	    mpq_set_z(c[i], tmp);

	    for (int j = 1; j <= i; ++j) {
		zz2value(d.coeff[j], tmp);
		mpq_set_z(qtmp, tmp);
		mpq_mul(qtmp, qtmp, c[i-j]);
		mpq_sub(c[i], c[i], qtmp);
	    }

	    zz2value(d.coeff[0], tmp);
	    mpq_set_z(qtmp, tmp);
	    mpq_div(c[i], c[i], qtmp);
	}
	if (sign == -1)
	    mpq_sub(count, count, c[len-1]);
	else
	    mpq_add(count, count, c[len-1]);

	value_clear(tmp);
	mpq_clear(qtmp);
	for (int i = 0; i < len; ++i)
	    mpq_clear(c[i]);
	delete [] c;
    }
};

class dpoly_n {
public:
    Matrix *coeff;
    ~dpoly_n() {
	Matrix_Free(coeff);
    }
    dpoly_n(int d, ZZ& degree_0, ZZ& degree_1, int offset = 0) {
	Value d0, d1;
	value_init(d0);
	value_init(d1);
	zz2value(degree_0, d0);
	zz2value(degree_1, d1);
	coeff = Matrix_Alloc(d+1, d+1+1);
	value_set_si(coeff->p[0][0], 1);
	value_set_si(coeff->p[0][d+1], 1);
	for (int i = 1; i <= d; ++i) {
	    value_multiply(coeff->p[i][0], coeff->p[i-1][0], d0);
	    Vector_Combine(coeff->p[i-1], coeff->p[i-1]+1, coeff->p[i]+1,
			   d1, d0, i);
	    value_set_si(coeff->p[i][d+1], i);
	    value_multiply(coeff->p[i][d+1], coeff->p[i][d+1], coeff->p[i-1][d+1]);
	    value_decrement(d0, d0);
	}
	value_clear(d0);
	value_clear(d1);
    }
    void div(dpoly& d, Vector *count, ZZ& sign) {
	int len = coeff->NbRows;
	Matrix * c = Matrix_Alloc(coeff->NbRows, coeff->NbColumns);
	Value tmp;
	value_init(tmp);
	for (int i = 0; i < len; ++i) {
	    Vector_Copy(coeff->p[i], c->p[i], len+1);
	    for (int j = 1; j <= i; ++j) {
		zz2value(d.coeff[j], tmp);
		value_multiply(tmp, tmp, c->p[i][len]);
		value_oppose(tmp, tmp);
		Vector_Combine(c->p[i], c->p[i-j], c->p[i],
			       c->p[i-j][len], tmp, len);
		value_multiply(c->p[i][len], c->p[i][len], c->p[i-j][len]);
	    }
	    zz2value(d.coeff[0], tmp);
	    value_multiply(c->p[i][len], c->p[i][len], tmp);
	}
	if (sign == -1) {
	    value_set_si(tmp, -1);
	    Vector_Scale(c->p[len-1], count->p, tmp, len);
	    value_assign(count->p[len], c->p[len-1][len]);
	} else
	    Vector_Copy(c->p[len-1], count->p, len+1);
	Vector_Normalize(count->p, len+1);
	value_clear(tmp);
	Matrix_Free(c);
    }
};

/*
 * Barvinok's Decomposition of a simplicial cone
 *
 * Returns two lists of polyhedra
 */
void barvinok_decompose(Polyhedron *C, Polyhedron **ppos, Polyhedron **pneg)
{
    Polyhedron *pos = *ppos, *neg = *pneg;
    vector<cone *> nonuni;
    cone * c = new cone(C);
    ZZ det = c->det;
    int s = sign(det);
    assert(det != 0);
    if (abs(det) > 1) {
	nonuni.push_back(c);
    } else {
	Polyhedron *p = Polyhedron_Copy(c->Cone);
	p->next = pos;
	pos = p;
	delete c;
    }
    vec_ZZ lambda;
    while (!nonuni.empty()) {
	c = nonuni.back();
	nonuni.pop_back();
	Vector* v = c->short_vector(lambda);
	for (int i = 0; i < c->Rays->NbRows - 1; ++i) {
	    if (lambda[i] == 0)
		continue;
	    Matrix* M = Matrix_Copy(c->Rays);
	    Vector_Copy(v->p, M->p[i], v->Size);
	    cone * pc = new cone(M);
	    assert (pc->det != 0);
	    if (abs(pc->det) > 1) {
		assert(abs(pc->det) < abs(c->det));
		nonuni.push_back(pc);
	    } else {
		Polyhedron *p = pc->poly();
		pc->Cone = 0;
		if (sign(pc->det) == s) {
		    p->next = pos;
		    pos = p;
		} else {
		    p->next = neg;
		    neg = p;
		}
		delete pc;
	    }
	    Matrix_Free(M);
	}
	Vector_Free(v);
	delete c;
    }
    *ppos = pos;
    *pneg = neg;
}

/*
 *  Returns a single list of npos "positive" cones followed by nneg
 *  "negative" cones.
 *  The input cone is freed
 */
void decompose(Polyhedron *cone, Polyhedron **parts, int *npos, int *nneg, unsigned MaxRays)
{
    Polyhedron_Polarize(cone);
    if (cone->NbRays - 1 != cone->Dimension) {
	Polyhedron *tmp = cone;
	cone = triangularize_cone(cone, MaxRays);
	Polyhedron_Free(tmp);
    }
    Polyhedron *polpos = NULL, *polneg = NULL;
    *npos = 0; *nneg = 0;
    for (Polyhedron *Polar = cone; Polar; Polar = Polar->next)
	barvinok_decompose(Polar, &polpos, &polneg);

    Polyhedron *last;
    for (Polyhedron *i = polpos; i; i = i->next) {
	Polyhedron_Polarize(i);
	++*npos;
	last = i;
    }
    for (Polyhedron *i = polneg; i; i = i->next) {
	Polyhedron_Polarize(i);
	++*nneg;
    }
    if (last) {
	last->next = polneg;
	*parts = polpos;
    } else
	*parts = polneg;
    Domain_Free(cone);
}

const int MAX_TRY=10;
/*
 * Searches for a vector that is not othogonal to any
 * of the rays in rays.
 */
static void nonorthog(mat_ZZ& rays, vec_ZZ& lambda)
{
    int dim = rays.NumCols();
    bool found = false;
    lambda.SetLength(dim);
    for (int i = 2; !found && i <= 50*dim; i+=4) {
	for (int j = 0; j < MAX_TRY; ++j) {
	    for (int k = 0; k < dim; ++k) {
		int r = random_int(i)+2;
		int v = (2*(r%2)-1) * (r >> 1);
		lambda[k] = v;
	    }
	    int k = 0;
	    for (; k < rays.NumRows(); ++k)
		if (lambda * rays[k] == 0)
		    break;
	    if (k == rays.NumRows()) {
		found = true;
		break;
	    }
	}
    }
    assert(found);
}

static void add_rays(mat_ZZ& rays, Polyhedron *i, int *r)
{
    unsigned dim = i->Dimension;
    for (int k = 0; k < i->NbRays; ++k) {
	if (!value_zero_p(i->Ray[k][dim+1]))
	    continue;
	values2zz(i->Ray[k]+1, rays[(*r)++], dim);
    }
}

void lattice_point(Value* values, Polyhedron *i, vec_ZZ& lambda, ZZ& num)
{
    vec_ZZ vertex;
    unsigned dim = i->Dimension;
    if(!value_one_p(values[dim])) {
	Matrix* Rays = rays(i);
	Matrix *inv = Matrix_Alloc(Rays->NbRows, Rays->NbColumns);
	int ok = Matrix_Inverse(Rays, inv);
	assert(ok);
	Matrix_Free(Rays);
	Rays = rays(i);
	Vector *lambda = Vector_Alloc(dim+1);
	Vector_Matrix_Product(values, inv, lambda->p);
	Matrix_Free(inv);
	for (int j = 0; j < dim; ++j)
	    mpz_cdiv_q(lambda->p[j], lambda->p[j], lambda->p[dim]);
	value_set_si(lambda->p[dim], 1);
	Vector *A = Vector_Alloc(dim+1);
	Vector_Matrix_Product(lambda->p, Rays, A->p);
	Vector_Free(lambda);
	Matrix_Free(Rays);
	values2zz(A->p, vertex, dim);
	Vector_Free(A);
    } else
	values2zz(values, vertex, dim);

    num = vertex * lambda;
}

static EhrhartPolynom *term(string param, ZZ& c, Value *den = NULL)
{
    evalue EP;
    deque<string> params;
    value_init(EP.d);
    value_set_si(EP.d,0);
    EP.x.p = new_enode(polynomial, 2, 1);
    value_init(EP.x.p->arr[0].x.n);
    value_init(EP.x.p->arr[1].x.n);
    value_set_si(EP.x.p->arr[0].d, 1);
    value_set_si(EP.x.p->arr[0].x.n, 0);
    if (den == NULL)
	value_set_si(EP.x.p->arr[1].d, 1);
    else
	value_assign(EP.x.p->arr[1].d, *den);
    zz2value(c, EP.x.p->arr[1].x.n);
    params.push_back(param);
    EhrhartPolynom * ret = new EhrhartPolynom(&EP, params);//--
    free_evalue_refs(&EP);
    return ret;
}

static void vertex_period(deque<string>& params, 
		    Polyhedron *i, vec_ZZ& lambda, Matrix *T, 
		    Value lcm, int p, Vector *val, 
		    EhrhartPolynom *E, evalue* ev,
		    ZZ& offset)
{
    unsigned nparam = T->NbRows - 1;
    unsigned dim = i->Dimension;
    Value tmp;
    ZZ nump;

    if (p == nparam) {
	ZZ num, l;
	Vector * values = Vector_Alloc(dim + 1);
	Vector_Matrix_Product(val->p, T, values->p);
	value_assign(values->p[dim], lcm);
	lattice_point(values->p, i, lambda, num);
	value2zz(lcm, l);
	num *= l;
	num += offset;
	value_init(ev->x.n);
	zz2value(num, ev->x.n);
	value_assign(ev->d, lcm);
	Vector_Free(values);
	return;
    }

    value_init(tmp);
    vec_ZZ vertex;
    values2zz(T->p[p], vertex, dim);
    nump = vertex * lambda;
    if (First_Non_Zero(val->p, p) == -1) {
	value_assign(tmp, lcm);
	EhrhartPolynom * ET = term(params[p], nump, &tmp);
	evalue EV1=(*ET).to_evalue(params);  
	evalue EV2=(*E).to_evalue(params);  
	eadd(&EV1,&EV2);   
	delete ET;
	*E=EhrhartPolynom(&EV2,params);
	free_evalue_refs(&EV1); 
	free_evalue_refs(&EV2);  
    }

    value_assign(tmp, lcm);
    if (First_Non_Zero(T->p[p], dim) != -1)
	Vector_Gcd(T->p[p], dim, &tmp);
    Gcd(tmp, lcm, &tmp);
    if (value_lt(tmp, lcm)) {
	ZZ count;

	value_division(tmp, lcm, tmp);
	value_set_si(ev->d, 0);
	ev->x.p = new_enode(periodic, VALUE_TO_INT(tmp), p+1);
	value2zz(tmp, count);
	do {
	    value_decrement(tmp, tmp);
	    --count;
	    ZZ new_offset = offset - count * nump;
	    value_assign(val->p[p], tmp);
	    vertex_period(params, i, lambda, T, lcm, p+1, val, E, 
			  &ev->x.p->arr[VALUE_TO_INT(tmp)], new_offset);
	} while (value_pos_p(tmp));
    } else
	vertex_period(params, i, lambda, T, lcm, p+1, val, E, ev, offset);
    value_clear(tmp);
}

static void mask_r(Matrix *f, int nr, Vector *lcm, int p, Vector *val, evalue *ev)
{
    unsigned nparam = lcm->Size;

    if (p == nparam) {
	Vector * prod = Vector_Alloc(f->NbRows);
	Matrix_Vector_Product(f, val->p, prod->p);
	int isint = 1;
	for (int i = 0; i < nr; ++i) {
	    value_modulus(prod->p[i], prod->p[i], f->p[i][nparam+1]);
	    isint &= value_zero_p(prod->p[i]);
	}
	value_set_si(ev->d, 1);
	value_init(ev->x.n);
	value_set_si(ev->x.n, isint);
	Vector_Free(prod);
	return;
    }

    Value tmp;
    value_init(tmp);
    if (value_one_p(lcm->p[p]))
	mask_r(f, nr, lcm, p+1, val, ev);
    else { 
	value_assign(tmp, lcm->p[p]);
	value_set_si(ev->d, 0);
	ev->x.p = new_enode(periodic, VALUE_TO_INT(tmp), p+1);
	do {
	    value_decrement(tmp, tmp);
	    value_assign(val->p[p], tmp);
	    mask_r(f, nr, lcm, p+1, val, &ev->x.p->arr[VALUE_TO_INT(tmp)]);
	} while (value_pos_p(tmp));
    }
    value_clear(tmp);
}

/*
 * 
 */
static void mask(Matrix *f, evalue *factor)
{
    int nr = f->NbRows, nc = f->NbColumns;
    int n;
    bool found = false;
    for (n = 0; n < nr && value_notzero_p(f->p[n][nc-1]); ++n)
	if (value_notone_p(f->p[n][nc-1]) &&
	    value_notmone_p(f->p[n][nc-1]))
		found = true;
    if (!found)
	return;

    Value tmp;
    value_init(tmp);
    nr = n;
    unsigned np = nc - 2;
    Vector *lcm = Vector_Alloc(np);
    Vector *val = Vector_Alloc(nc);
    Vector_Set(val->p, 0, nc);
    value_set_si(val->p[np], 1);
    Vector_Set(lcm->p, 1, np);
    for (n = 0; n < nr; ++n) {
	if (value_one_p(f->p[n][nc-1]) ||
	    value_mone_p(f->p[n][nc-1]))
	    continue;
	for (int j = 0; j < np; ++j)
	    if (value_notzero_p(f->p[n][j])) {
		Gcd(f->p[n][j], f->p[n][nc-1], &tmp);
		value_division(tmp, f->p[n][nc-1], tmp);
		value_lcm(tmp, lcm->p[j], &lcm->p[j]);
	    }
    }
    evalue EP;
    value_init(EP.d);
    mask_r(f, nr, lcm, 0, val, &EP);
    value_clear(tmp);
    Vector_Free(val);
    Vector_Free(lcm);
    emul(&EP,factor); 
    free_evalue_refs(&EP);
}

static EhrhartPolynom *multi_mononom(deque<string>& params, vec_ZZ& p)
{
    EhrhartPolynom *X = new EhrhartPolynom();
    evalue EV1=(*X).to_evalue(params);  
    unsigned nparam = p.length()-1;
    for (int i = 0; i < nparam; ++i) {
	EhrhartPolynom *T = term(params[i], p[i]);
	evalue EV2=(*T).to_evalue(params); 
	eadd(&EV2,&EV1); 
	delete T;
	free_evalue_refs(&EV2); 
    }
    *X=EhrhartPolynom(&EV1,params) ;  
    free_evalue_refs(&EV1); 
    return X;
}

struct term_info {
    EhrhartPolynom *E;
    ZZ		    constant;
    ZZ		    coeff;
    int		    pos;
};

void lattice_point(deque<string>& params, 
    Param_Vertices* V, Polyhedron *i, vec_ZZ& lambda, term_info* term)
{
    unsigned nparam = V->Vertex->NbColumns - 2;
    unsigned dim = i->Dimension;
    mat_ZZ vertex;
    vertex.SetDims(V->Vertex->NbRows, nparam+1);
    Value lcm, tmp;
    value_init(lcm);
    value_init(tmp);
    value_set_si(lcm, 1);
    for (int j = 0; j < V->Vertex->NbRows; ++j) {
	value_lcm(lcm, V->Vertex->p[j][nparam+1], &lcm);
    }
    if (value_notone_p(lcm)) {
	Matrix* Rays = rays(i);
	Matrix *inv = Matrix_Alloc(Rays->NbRows, Rays->NbColumns);
	int ok = Matrix_Inverse(Rays, inv);
	assert(ok);
	Matrix_Free(Rays);
	Rays = rays(i);

	Matrix * mv = Matrix_Alloc(dim, nparam+1);
	for (int j = 0 ; j < dim; ++j) {
	    value_division(tmp, lcm, V->Vertex->p[j][nparam+1]);
	    Vector_Scale(V->Vertex->p[j], mv->p[j], tmp, nparam+1);
	}
	Matrix *T = Transpose(mv);

	EhrhartPolynom * EP = new EhrhartPolynom();
	 
	evalue ev;
	Vector *val = Vector_Alloc(nparam+1);
	value_set_si(val->p[nparam], 1);
	ZZ offset(INIT_VAL, 0);
	value_init(ev.d);
	vertex_period(params, i, lambda, T, lcm, 0, val, EP, &ev, offset);
	Vector_Free(val);
        evalue  EV=(*EP).to_evalue(params); 
        eadd(&ev,&EV);
	*EP=EhrhartPolynom(&EV,params);  
	   free_evalue_refs(&ev);   
	   free_evalue_refs(&EV)  ;   

	term->E = EP;
	term->constant = 0;

	Matrix_Free(inv);
	Matrix_Free(Rays);
	Matrix_Free(T);
	Matrix_Free(mv);
	value_clear(lcm);
	value_clear(tmp);
	return;
    }
    for (int i = 0; i < V->Vertex->NbRows; ++i) {
	assert(value_one_p(V->Vertex->p[i][nparam+1]));  // for now
	values2zz(V->Vertex->p[i], vertex[i], nparam+1);
    }

    vec_ZZ num;
    num = lambda * vertex;

    int p = -1;
    int nn = 0;
    for (int j = 0; j < nparam; ++j)
	if (num[j] != 0) {
	    ++nn;
	    p = j;
	}
    if (nn >= 2) {
	term->E = multi_mononom(params, num);
	term->constant = num[nparam];
    } else {
	term->E = NULL;
	term->constant = num[nparam];
	term->pos = p;
	if (p != -1)
	    term->coeff = num[p];
    }

    value_clear(lcm);
    value_clear(tmp);
}

void normalize(Polyhedron *i, vec_ZZ& lambda, ZZ& sign, ZZ& num, vec_ZZ& den)
{
    unsigned dim = i->Dimension;

    int r = 0;
    mat_ZZ rays;
    rays.SetDims(dim, dim);
    add_rays(rays, i, &r);
    den = rays * lambda;
    int change = 0;

    for (int j = 0; j < den.length(); ++j) {
	if (den[j] > 0)
	    change ^= 1;
	else {
	    den[j] = abs(den[j]);
	    num += den[j];
	}
    }
    if (change)
	sign = -sign;
}

void barvinok_count(Polyhedron *P, Value* result, unsigned NbMaxCons)
{
    Polyhedron ** vcone;
    vec_ZZ sign;
    int ncone = 0;
    sign.SetLength(ncone);
    unsigned dim;
    int allocated = 0;
    Value factor;
    Polyhedron *Q;
    int r = 0;

    if (emptyQ(P)) {
	value_set_si(*result, 0);
	return;
    }
    if (P->NbBid == 0)
	for (; r < P->NbRays; ++r)
	    if (value_zero_p(P->Ray[r][P->Dimension+1]))
		break;
    if (P->NbBid !=0 || r < P->NbRays) {
	value_set_si(*result, -1);
	return;
    }
    if (P->NbEq != 0) {
	P = remove_equalities(P);
	if (emptyQ(P)) {
	    Polyhedron_Free(P);
	    value_set_si(*result, 0);
	    return;
	}
	allocated = 1;
    }
    value_init(factor);
    value_set_si(factor, 1);
    Q = Polyhedron_Reduce(P, &factor);
    if (Q) {
	if (allocated)
	    Polyhedron_Free(P);
	P = Q;
	allocated = 1;
    }
    if (P->Dimension == 0) {
	value_assign(*result, factor);
	if (allocated)
	    Polyhedron_Free(P);
	value_clear(factor);
	return;
    }

    dim = P->Dimension;
    vcone = new (Polyhedron *)[P->NbRays];

    for (int j = 0; j < P->NbRays; ++j) {
	int npos, nneg;
	Polyhedron *C = supporting_cone(P, j);
	decompose(C, &vcone[j], &npos, &nneg, NbMaxCons);
	ncone += npos + nneg;
	sign.SetLength(ncone);
	for (int k = 0; k < npos; ++k)
	    sign[ncone-nneg-k-1] = 1;
	for (int k = 0; k < nneg; ++k)
	    sign[ncone-k-1] = -1;
    }

    mat_ZZ rays;
    rays.SetDims(ncone * dim, dim);
    r = 0;
    for (int j = 0; j < P->NbRays; ++j) {
	for (Polyhedron *i = vcone[j]; i; i = i->next) {
	    assert(i->NbRays-1 == dim);
	    add_rays(rays, i, &r);
	}
    }
    vec_ZZ lambda;
    nonorthog(rays, lambda);

    vec_ZZ num;
    mat_ZZ den;
    num.SetLength(ncone);
    den.SetDims(ncone,dim);

    int f = 0;
    for (int j = 0; j < P->NbRays; ++j) {
	for (Polyhedron *i = vcone[j]; i; i = i->next) {
	    lattice_point(P->Ray[j]+1, i, lambda, num[f]);
	    normalize(i, lambda, sign[f], num[f], den[f]);
	    ++f;
	}
    }
    ZZ min = num[0];
    for (int j = 1; j < num.length(); ++j)
	if (num[j] < min)
	    min = num[j];
    for (int j = 0; j < num.length(); ++j)
	num[j] -= min;

    f = 0;
    mpq_t count;
    mpq_init(count);
    for (int j = 0; j < P->NbRays; ++j) {
	for (Polyhedron *i = vcone[j]; i; i = i->next) {
	    dpoly d(dim, num[f]);
	    dpoly n(dim, den[f][0], 1);
	    for (int k = 1; k < dim; ++k) {
		dpoly fact(dim, den[f][k], 1);
		n *= fact;
	    }
	    d.div(n, count, sign[f]);
	    ++f;
	}
    }
    assert(value_one_p(&count[0]._mp_den));
    value_multiply(*result, &count[0]._mp_num, factor);
    mpq_clear(count);

    for (int j = 0; j < P->NbRays; ++j)
	Domain_Free(vcone[j]);

    delete [] vcone;

    if (allocated)
	Polyhedron_Free(P);
    value_clear(factor);
}

static void default_params(deque<string>& params, int n)
{
    for (int i = 1; i <= n; ++i) {
	ostringstream s;
	s << "p" << i;
	params.push_back(s.str());
    }
}

static EhrhartPolynom *uni_polynom(string param, Vector *c)
{ 
    evalue EP;
    deque<string> params;
    unsigned dim = c->Size-2;
    value_init(EP.d);
    value_set_si(EP.d,0);
    EP.x.p = new_enode(polynomial, dim+1, 1);
    for (int j = 0; j <= dim; ++j) {
	value_assign(EP.x.p->arr[j].d, c->p[dim+1]);
	value_init(EP.x.p->arr[j].x.n);
	value_assign(EP.x.p->arr[j].x.n, c->p[j]);
    }
    params.push_back(param);
    EhrhartPolynom * ret = new EhrhartPolynom(&EP, params);
    free_evalue_refs(&EP);
    return ret;
}

static EhrhartPolynom *multi_polynom(deque<string>& params, Vector *c, EhrhartPolynom& X)
{
 	
    unsigned dim = c->Size-2;
    evalue EC;
    value_init(EC.d);
    value_init(EC.x.n);
    value_assign(EC.d, c->p[dim+1]);
        
    EhrhartPolynom *res = new EhrhartPolynom(); 
    value_assign(EC.x.n, c->p[dim]);

    evalue  EV1=(*res).to_evalue(params); 
    eadd(&EC,&EV1);   
    evalue EV2=X.to_evalue(params);  
    for (int i = dim-1; i >= 0; --i) {
	emul(&EV2,&EV1);       
	value_assign(EC.x.n, c->p[i]);
	eadd(&EC,&EV1);
    }
    *res = EhrhartPolynom(&EV1, params);
    free_evalue_refs(&EC);
    free_evalue_refs(&EV2);
    free_evalue_refs(&EV1);
    return res;
}

static EhrhartPolynom *constant(mpq_t c)
{
    evalue EP;
    deque<string> params;
    value_init(EP.d);
    value_init(EP.x.n);
    value_assign(EP.d, &c[0]._mp_den);
    value_assign(EP.x.n, &c[0]._mp_num);
    EhrhartPolynom * ret = new EhrhartPolynom(&EP, params);
    free_evalue_refs(&EP);
    return ret;
}


Enumeration* barvinok_enumerate(Polyhedron *P, Polyhedron* C, unsigned MaxRays)
{
    Polyhedron *CEq = NULL, *rVD, *CA;
    Matrix *CT = NULL;
    Param_Polyhedron *PP;
    Param_Domain *D, *next;
    Param_Vertices *V;
    Enumeration *en, *res;
    int r = 0;
    unsigned nparam = C->Dimension;
    evalue factor;
    value_init(factor.d);
    evalue_set_si(&factor, 1, 1);
    
    res = NULL;

    CA = align_context(C, P->Dimension, MaxRays);
    P = DomainIntersection(P, CA, MaxRays);
    Polyhedron_Free(CA);

    deque<string> params, allparams;
    default_params(params, nparam);
    allparams = params;

    if (C->Dimension == 0 || emptyQ(P)) {
constant:
	res = (Enumeration *)malloc(sizeof(Enumeration));
	res->ValidityDomain = CEq ? CEq : Polyhedron_Copy(C);
	res->next = NULL;
	value_init(res->EP.d);
	value_set_si(res->EP.d, 1);
	value_init(res->EP.x.n);
	if (emptyQ(P))
	    value_set_si(res->EP.x.n, 0);
	else
	    barvinok_count(P, &res->EP.x.n, MaxRays);
	emul(&factor, &res->EP);
out:
	free_evalue_refs(&factor);
	Polyhedron_Free(P);
	if (CT)
	    Matrix_Free(CT);
	   
	return res;
    }

    if (P->NbEq != 0) {
	Matrix *f;
	P = remove_equalities_p(P, P->Dimension-nparam, &f);
	mask(f, &factor);
	Matrix_Free(f);
	if (P->Dimension == nparam) {
	    CEq = P;
	    P = Universe_Polyhedron(0);
	    goto constant;
	}
    }
    Polyhedron *oldP = P;
    PP = Polyhedron2Param_SimplifiedDomain(&P,C,MaxRays,&CEq,&CT);
    if (P != oldP)
	Polyhedron_Free(oldP);

    if (isIdentity(CT)) {
	Matrix_Free(CT);
	CT = NULL;
    } else {
	assert(CT->NbRows != CT->NbColumns);
	if (CT->NbRows == 1) {		// no more parameters
	    assert(PP->D->next == NULL);
	    Param_Polyhedron_Free(PP);
	    goto constant;
	}
	deque<string>::iterator i;
	params.erase(params.begin(), params.end());
	int r = 0, j, p = -1;
	for (i = allparams.begin(), j = 0; i != allparams.end(); ++i, ++j) {
	    if (p < j) {
		if (r >= CT->NbRows - 1)
		    break;
		p = First_Non_Zero(CT->p[r], nparam);
		assert(p != -1);
		assert(First_Non_Zero(CT->p[r]+p+1, nparam-p-1) == -1);
		assert(value_one_p(CT->p[r][p]));
		++r;
	    }
	    if (p == j)
		params.push_back(*i);
	}
	nparam = CT->NbRows - 1;
    }

    unsigned dim = P->Dimension - nparam;
    Polyhedron ** vcone = new (Polyhedron *)[PP->nbV];
    int * npos = new int[PP->nbV];
    int * nneg = new int[PP->nbV];
    vec_ZZ sign;

    int i;
    for (i = 0, V = PP->V; V; ++i, V = V->next) {
	Polyhedron *C = supporting_cone_p(P, V);
	decompose(C, &vcone[i], &npos[i], &nneg[i], MaxRays);
    }

    Vector *c = Vector_Alloc(dim+2);

    for(D=PP->D; D; D=next) {
	next = D->next;
	if (!CEq) {
	    rVD = D->Domain;    
	    D->Domain = NULL;
	} else {
	  Polyhedron *Dt;
	  Dt = CT ? Polyhedron_Preimage(D->Domain,CT,MaxRays) : D->Domain;
	  rVD = DomainIntersection(Dt,CEq,MaxRays);
	  
	  /* if rVD is empty or too small in geometric dimension */
	  if(!rVD || emptyQ(rVD) ||
	     (rVD->Dimension-rVD->NbEq < Dt->Dimension-Dt->NbEq-CEq->NbEq)) {
	    if(rVD)
	      Polyhedron_Free(rVD);
	    if (CT)
		Polyhedron_Free(Dt);
	    continue;		/* empty validity domain */
	  }
	  if (CT)
	      Polyhedron_Free(Dt);
	}
	int ncone = 0;
	sign.SetLength(ncone);
	FORALL_PVertex_in_ParamPolyhedron(V,D,PP) // _i is internal counter
	    ncone += npos[_i] + nneg[_i];
	    sign.SetLength(ncone);
	    for (int k = 0; k < npos[_i]; ++k)
		sign[ncone-nneg[_i]-k-1] = 1;
	    for (int k = 0; k < nneg[_i]; ++k)
		sign[ncone-k-1] = -1;
	END_FORALL_PVertex_in_ParamPolyhedron;

	mat_ZZ rays;
	rays.SetDims(ncone * dim, dim);
	r = 0;
	FORALL_PVertex_in_ParamPolyhedron(V,D,PP) // _i is internal counter
	    for (Polyhedron *i = vcone[_i]; i; i = i->next) {
		assert(i->NbRays-1 == dim);
		add_rays(rays, i, &r);
	    }
	END_FORALL_PVertex_in_ParamPolyhedron;
	vec_ZZ lambda;
	nonorthog(rays, lambda);

	mat_ZZ den;
	den.SetDims(ncone,dim);
	term_info *num = new term_info[ncone];
          
	int f = 0;
	FORALL_PVertex_in_ParamPolyhedron(V,D,PP)
	    for (Polyhedron *i = vcone[_i]; i; i = i->next) {
		lattice_point(params, V, i, lambda, &num[f]);
		normalize(i, lambda, sign[f], num[f].constant, den[f]);
		++f;
	    }
	END_FORALL_PVertex_in_ParamPolyhedron;
	ZZ min = num[0].constant;
	for (int j = 1; j < ncone; ++j)
	    if (num[j].constant < min)
		min = num[j].constant;
	for (int j = 0; j < ncone; ++j)
	    num[j].constant -= min;
	f = 0;
	evalue EP;
	value_init(EP.d);
	value_init(EP.x.n);
	value_set_si(EP.d, 1);
	value_set_si(EP.x.n, 0);
	mpq_t count;
	mpq_init(count);
	FORALL_PVertex_in_ParamPolyhedron(V,D,PP)
	    for (Polyhedron *i = vcone[_i]; i; i = i->next) {
		dpoly n(dim, den[f][0], 1);
		for (int k = 1; k < dim; ++k) {
		    dpoly fact(dim, den[f][k], 1);
		    n *= fact;
		}
		if (num[f].E != NULL) {
		    ZZ one(INIT_VAL, 1);
		    dpoly_n d(dim, num[f].constant, one);
		    d.div(n, c, sign[f]);
		    EhrhartPolynom *ET = multi_polynom(params, c, *num[f].E);
		    evalue EV = ET->to_evalue(params); 
		    eadd(&EV , &EP);
		    delete ET;
		    free_evalue_refs(&EV);
		    delete num[f].E; 
		} else if (num[f].pos != -1) {
		    dpoly_n d(dim, num[f].constant, num[f].coeff);
		    d.div(n, c, sign[f]);
		    EhrhartPolynom *E = uni_polynom(params[num[f].pos], c);
		    evalue EV = E->to_evalue(params);
		    eadd(&EV , &EP);
		    delete E;
		    free_evalue_refs(&EV);
		} else {
		    mpq_set_si(count, 0, 1);
		    dpoly d(dim, num[f].constant);
		    d.div(n, count, sign[f]);
		    EhrhartPolynom *E = constant(count);
		    evalue EV = E->to_evalue(params);
		    eadd(&EV , &EP);
		    delete E;
		    free_evalue_refs(&EV);
		    
		} 
		++f;
	    }
	END_FORALL_PVertex_in_ParamPolyhedron;

	mpq_clear(count);
	delete [] num;

	en = (Enumeration *)malloc(sizeof(Enumeration));
	en->next = res;
	res = en;
	res->ValidityDomain = rVD;
	if (CT)
	    addeliminatedparams_evalue(&EP, CT);
   	emul(&factor, &EP); 
	res->EP = EP;
	reduce_evalue(&res->EP);
    }

    Vector_Free(c);

    for (int j = 0; j < PP->nbV; ++j)
	Domain_Free(vcone[j]);
    delete [] vcone;
    delete [] npos;
    delete [] nneg;

    Param_Polyhedron_Free(PP);
    if (CEq)
	Polyhedron_Free(CEq);

    goto out;
}
