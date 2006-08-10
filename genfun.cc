#include <iostream>
#include <vector>
#include <assert.h>
#include "config.h"
#include <barvinok/genfun.h>
#include <barvinok/barvinok.h>
#include "conversion.h"
#include "genfun_constructor.h"
#include "mat_util.h"

using std::cout;
using std::cerr;
using std::endl;
using std::pair;
using std::vector;

static int lex_cmp(mat_ZZ& a, mat_ZZ& b)
{
    assert(a.NumCols() == b.NumCols());
    int alen = a.NumRows();
    int blen = b.NumRows();
    int len = alen < blen ? alen : blen;

    for (int i = 0; i < len; ++i) {
	int s = lex_cmp(a[i], b[i]);
	if (s)
	    return s;
    }
    return alen-blen;
}

static void lex_order_terms(struct short_rat* rat)
{
    for (int i = 0; i < rat->n.power.NumRows(); ++i) {
	int m = i;
	for (int j = i+1; j < rat->n.power.NumRows(); ++j)
	    if (lex_cmp(rat->n.power[j], rat->n.power[m]) < 0)
		m = j;
	if (m != i) {
	    vec_ZZ tmp = rat->n.power[m];
	    rat->n.power[m] = rat->n.power[i];
	    rat->n.power[i] = tmp;
	    QQ tmp_coeff = rat->n.coeff[m];
	    rat->n.coeff[m] = rat->n.coeff[i];
	    rat->n.coeff[i] = tmp_coeff;
	}
    }
}

void short_rat::add(short_rat *r)
{
    for (int i = 0; i < r->n.power.NumRows(); ++i) {
	int len = n.coeff.length();
	int j;
	for (j = 0; j < len; ++j)
	    if (r->n.power[i] == n.power[j])
		break;
	if (j < len) {
	    n.coeff[j] += r->n.coeff[i];
	    if (n.coeff[j].n == 0) {
		if (j < len-1) {
		    n.power[j] = n.power[len-1];
		    n.coeff[j] = n.coeff[len-1];
		}
		int dim = n.power.NumCols();
		n.coeff.SetLength(len-1);
		n.power.SetDims(len-1, dim);
	    }
	} else {
	    int dim = n.power.NumCols();
	    n.coeff.SetLength(len+1);
	    n.power.SetDims(len+1, dim);
	    n.coeff[len] = r->n.coeff[i];
	    n.power[len] = r->n.power[i];
	}
    }
}

bool short_rat::reduced()
{
    int dim = n.power.NumCols();
    lex_order_terms(this);
    if (n.power.NumRows() % 2 == 0) {
	if (n.coeff[0].n == -n.coeff[1].n &&
	    n.coeff[0].d == n.coeff[1].d) {
	    vec_ZZ step = n.power[1] - n.power[0];
	    int k;
	    for (k = 1; k < n.power.NumRows()/2; ++k) {
		if (n.coeff[2*k].n != -n.coeff[2*k+1].n ||
		    n.coeff[2*k].d != n.coeff[2*k+1].d)
		    break;
		if (step != n.power[2*k+1] - n.power[2*k])
		    break;
	    }
	    if (k == n.power.NumRows()/2) {
		for (k = 0; k < d.power.NumRows(); ++k)
		    if (d.power[k] == step)
			break;
		if (k < d.power.NumRows()) {
		    for (++k; k < d.power.NumRows(); ++k)
			d.power[k-1] = d.power[k];
		    d.power.SetDims(k-1, dim);
		    for (k = 1; k < n.power.NumRows()/2; ++k) {
			n.coeff[k] = n.coeff[2*k];
			n.power[k] = n.power[2*k];
		    }
		    n.coeff.SetLength(k);
		    n.power.SetDims(k, dim);
		    return true;
		}
	    }
	}
    }
    return false;
}

void gen_fun::add(const QQ& c, const vec_ZZ& num, const mat_ZZ& den)
{
    if (c.n == 0)
	return;

    short_rat * r = new short_rat;
    r->n.coeff.SetLength(1);
    ZZ g = GCD(c.n, c.d);
    r->n.coeff[0].n = c.n/g;
    r->n.coeff[0].d = c.d/g;
    r->n.power.SetDims(1, num.length());
    r->n.power[0] = num;
    r->d.power = den;

    /* Make all powers in denominator lexico-positive */
    for (int i = 0; i < r->d.power.NumRows(); ++i) {
	int j;
	for (j = 0; j < r->d.power.NumCols(); ++j)
	    if (r->d.power[i][j] != 0)
		break;
	if (r->d.power[i][j] < 0) {
	    r->d.power[i] = -r->d.power[i];
	    r->n.coeff[0].n = -r->n.coeff[0].n;
	    r->n.power[0] += r->d.power[i];
	}
    }

    /* Order powers in denominator */
    lex_order_rows(r->d.power);

    for (int i = 0; i < term.size(); ++i)
	if (lex_cmp(term[i]->d.power, r->d.power) == 0) {
	    term[i]->add(r);
	    if (term[i]->n.coeff.length() == 0) {
		delete term[i];
		if (i != term.size()-1)
		    term[i] = term[term.size()-1];
		term.pop_back();
	    } else if (term[i]->reduced()) {
		delete r;
		/* we've modified term[i], so removed it
		 * and add it back again
		 */
		r = term[i];
		if (i != term.size()-1)
		    term[i] = term[term.size()-1];
		term.pop_back();
		i = -1;
		continue;
	    }
	    delete r;
	    return;
	}

    term.push_back(r);
}

void gen_fun::add(const QQ& c, const gen_fun *gf)
{
    QQ p;
    for (int i = 0; i < gf->term.size(); ++i) {
	for (int j = 0; j < gf->term[i]->n.power.NumRows(); ++j) {
	    p = c;
	    p *= gf->term[i]->n.coeff[j];
	    add(p, gf->term[i]->n.power[j], gf->term[i]->d.power);
	}
    }
}

static void split_param_compression(Matrix *CP, mat_ZZ& map, vec_ZZ& offset)
{
    Matrix *T = Transpose(CP);
    matrix2zz(T, map, T->NbRows-1, T->NbColumns-1);
    values2zz(T->p[T->NbRows-1], offset, T->NbColumns-1);
    Matrix_Free(T);
}

/*
 * Perform the substitution specified by CP
 *
 * CP is a homogeneous matrix that maps a set of "compressed parameters"
 * to the original set of parameters.
 *
 * This function is applied to a gen_fun computed with the compressed parameters
 * and adapts it to refer to the original parameters.
 *
 * That is, if y are the compressed parameters and x = A y + b are the original
 * parameters, then we want the coefficient of the monomial t^y in the original
 * generating function to be the coefficient of the monomial u^x in the resulting
 * generating function.  
 * The original generating function has the form
 *
 *        a t^m/(1-t^n) = a t^m + a t^{m+n} + a t^{m+2n} + ...
 * 
 * Since each term t^y should correspond to a term u^x, with x = A y + b, we want
 *
 *         a u^{A m + b} + a u^{A (m+n) + b} + a u^{A (m+2n) +b} + ... = 
 *        
 *         = a u^{A m + b}/(1-u^{A n})
 *
 * Therefore, we multiply the powers m and n in both numerator and denominator by A
 * and add b to the power in the numerator.
 * Since the above powers are stored as row vectors m^T and n^T,
 * we compute, say, m'^T = m^T A^T to obtain m' = A m.
 *
 * The pair (map, offset) contains the same information as CP.
 * map is the transpose of the linear part of CP, while offset is the constant part.
 */
void gen_fun::substitute(Matrix *CP)
{
    mat_ZZ map;
    vec_ZZ offset;
    split_param_compression(CP, map, offset);
    Polyhedron *C = Polyhedron_Image(context, CP, 0);
    Polyhedron_Free(context);
    context = C;
    for (int i = 0; i < term.size(); ++i) {
	term[i]->d.power *= map;
	term[i]->n.power *= map;
	for (int j = 0; j < term[i]->n.power.NumRows(); ++j)
	    term[i]->n.power[j] += offset;
    }
}

struct cone {
    int	    *pos;
    vector<pair<Vector *, QQ> >	    vertices;
    cone(int *pos) : pos(pos) {}
};

#ifndef HAVE_COMPRESS_PARMS
static Matrix *compress_parms(Matrix *M, unsigned nparam)
{
    assert(0);
}
#endif

struct parallel_polytopes {
    gf_base *red;
    Matrix *Constraints;
    Matrix *CP, *T;
    int dim;
    int nparam;
    vector<cone>    cones;

    parallel_polytopes(int n, Polyhedron *context, int dim, int nparam) :
			dim(dim), nparam(nparam) {
	red = gf_base::create(Polyhedron_Copy(context), dim, nparam);
	Constraints = NULL;
	CP = NULL;
	T = NULL;
    }
    void add(const QQ& c, Polyhedron *P, unsigned MaxRays) {
	Polyhedron *Q = remove_equalities_p(Polyhedron_Copy(P), P->Dimension-nparam,
					    NULL);
	POL_ENSURE_VERTICES(Q);
	if (emptyQ(Q)) {
	    Polyhedron_Free(Q);
	    return;
	}

	if (Q->NbEq != 0) {
	    Polyhedron *R;
	    if (!CP) {
		Matrix *M;
		M = Matrix_Alloc(Q->NbEq, Q->Dimension+2);
		Vector_Copy(Q->Constraint[0], M->p[0], Q->NbEq * (Q->Dimension+2));
		CP = compress_parms(M, nparam);
		T = align_matrix(CP, Q->Dimension+1);
		Matrix_Free(M);
	    }
	    R = Polyhedron_Preimage(Q, T, MaxRays);
	    Polyhedron_Free(Q);
	    Q = remove_equalities_p(R, R->Dimension-nparam, NULL);
	}
	assert(Q->NbEq == 0);
	assert(Q->Dimension == dim);

	if (First_Non_Zero(Q->Constraint[Q->NbConstraints-1]+1, Q->Dimension) == -1)
	    Q->NbConstraints--;

	if (!Constraints) {
	    red->base->init(Q);
	    Constraints = Matrix_Alloc(Q->NbConstraints, Q->Dimension);
	    for (int i = 0; i < Q->NbConstraints; ++i) {
		Vector_Copy(Q->Constraint[i]+1, Constraints->p[i], Q->Dimension);
	    }
	} else {
	    for (int i = 0; i < Q->NbConstraints; ++i) {
		int j;
		for (j = 0; j < Constraints->NbRows; ++j)
		    if (Vector_Equal(Q->Constraint[i]+1, Constraints->p[j],
					Q->Dimension))
			break;
		assert(j < Constraints->NbRows);
	    }
	}

	for (int i = 0; i < Q->NbRays; ++i) {
	    if (!value_pos_p(Q->Ray[i][dim+1]))
		continue;

	    Polyhedron *C = supporting_cone(Q, i);

	    if (First_Non_Zero(C->Constraint[C->NbConstraints-1]+1, 
			       C->Dimension) == -1)
		C->NbConstraints--;

	    int *pos = new int[1+C->NbConstraints];
	    pos[0] = C->NbConstraints;
	    int l = 0;
	    for (int k = 0; k < Constraints->NbRows; ++k) {
		for (int j = 0; j < C->NbConstraints; ++j) {
		    if (Vector_Equal(C->Constraint[j]+1, Constraints->p[k], 
				     C->Dimension)) {
			pos[1+l++] = k;
			break;
		    }
		}
	    }
	    assert(l == C->NbConstraints);

	    int j;
	    for (j = 0; j < cones.size(); ++j)
		if (!memcmp(pos, cones[j].pos, (1+C->NbConstraints)*sizeof(int)))
		    break;
	    if (j == cones.size())
		cones.push_back(cone(pos));
	    else
		delete [] pos;

	    Polyhedron_Free(C);

	    int k;
	    for (k = 0; k < cones[j].vertices.size(); ++k)
		if (Vector_Equal(Q->Ray[i]+1, cones[j].vertices[k].first->p,
				 Q->Dimension+1))
		    break;

	    if (k == cones[j].vertices.size()) {
		Vector *vertex = Vector_Alloc(Q->Dimension+1);
		Vector_Copy(Q->Ray[i]+1, vertex->p, Q->Dimension+1);
		cones[j].vertices.push_back(pair<Vector*,QQ>(vertex, c));
	    } else {
		cones[j].vertices[k].second += c;
		if (cones[j].vertices[k].second.n == 0) {
		    int size = cones[j].vertices.size();
		    Vector_Free(cones[j].vertices[k].first);
		    if (k < size-1)
			cones[j].vertices[k] = cones[j].vertices[size-1];
		    cones[j].vertices.pop_back();
		}
	    }
	}

	Polyhedron_Free(Q);
    }
    gen_fun *compute(unsigned MaxRays) {
	for (int i = 0; i < cones.size(); ++i) {
	    Matrix *M = Matrix_Alloc(cones[i].pos[0], 1+Constraints->NbColumns+1);
	    Polyhedron *Cone;
	    for (int j = 0; j <cones[i].pos[0]; ++j) {
		value_set_si(M->p[j][0], 1);
		Vector_Copy(Constraints->p[cones[i].pos[1+j]], M->p[j]+1,
			    Constraints->NbColumns);
	    }
	    Cone = Constraints2Polyhedron(M, MaxRays);
	    Matrix_Free(M);
	    for (int j = 0; j < cones[i].vertices.size(); ++j) {
		red->base->do_vertex_cone(cones[i].vertices[j].second,
					  Polyhedron_Copy(Cone),
					  cones[i].vertices[j].first->p,
					  MaxRays);
	    }
	    Polyhedron_Free(Cone);
	}
	if (CP)
	    red->gf->substitute(CP);
	return red->gf;
    }
    void print(std::ostream& os) const {
	for (int i = 0; i < cones.size(); ++i) {
	    os << "[";
	    for (int j = 0; j < cones[i].pos[0]; ++j) {
		if (j)
		    os << ", ";
		os << cones[i].pos[1+j];
	    }
	    os << "]" << endl;
	    for (int j = 0; j < cones[i].vertices.size(); ++j) {
		Vector_Print(stderr, P_VALUE_FMT, cones[i].vertices[j].first);
		os << cones[i].vertices[j].second << endl;
	    }
	}
    }
    ~parallel_polytopes() {
	for (int i = 0; i < cones.size(); ++i) {
	    delete [] cones[i].pos;
	    for (int j = 0; j < cones[i].vertices.size(); ++j)
		Vector_Free(cones[i].vertices[j].first);
	}
	if (Constraints)
	    Matrix_Free(Constraints);
	if (CP)
	    Matrix_Free(CP);
	if (T)
	    Matrix_Free(T);
	delete red;
    }
};

gen_fun *gen_fun::Hadamard_product(const gen_fun *gf, unsigned MaxRays)
{
    QQ one(1, 1);
    Polyhedron *C = DomainIntersection(context, gf->context, MaxRays);
    Polyhedron *U = Universe_Polyhedron(C->Dimension);
    gen_fun *sum = new gen_fun(C);
    for (int i = 0; i < term.size(); ++i) {
	for (int i2 = 0; i2 < gf->term.size(); ++i2) {
	    int d = term[i]->d.power.NumCols();
	    int k1 = term[i]->d.power.NumRows();
	    int k2 = gf->term[i2]->d.power.NumRows();
	    assert(term[i]->d.power.NumCols() == gf->term[i2]->d.power.NumCols());

	    parallel_polytopes pp(term[i]->n.power.NumRows() *
				  gf->term[i2]->n.power.NumRows(),
				  sum->context, k1+k2-d, d);

	    for (int j = 0; j < term[i]->n.power.NumRows(); ++j) {
		for (int j2 = 0; j2 < gf->term[i2]->n.power.NumRows(); ++j2) {
		    Matrix *M = Matrix_Alloc(k1+k2+d+d, 1+k1+k2+d+1);
		    for (int k = 0; k < k1+k2; ++k) {
			value_set_si(M->p[k][0], 1);
			value_set_si(M->p[k][1+k], 1);
		    }
		    for (int k = 0; k < d; ++k) {
			value_set_si(M->p[k1+k2+k][1+k1+k2+k], -1);
			zz2value(term[i]->n.power[j][k], M->p[k1+k2+k][1+k1+k2+d]);
			for (int l = 0; l < k1; ++l)
			    zz2value(term[i]->d.power[l][k], M->p[k1+k2+k][1+l]);
		    }
		    for (int k = 0; k < d; ++k) {
			value_set_si(M->p[k1+k2+d+k][1+k1+k2+k], -1);
			zz2value(gf->term[i2]->n.power[j2][k], 
				 M->p[k1+k2+d+k][1+k1+k2+d]);
			for (int l = 0; l < k2; ++l)
			    zz2value(gf->term[i2]->d.power[l][k], 
				     M->p[k1+k2+d+k][1+k1+l]);
		    }
		    Polyhedron *P = Constraints2Polyhedron(M, MaxRays);
		    Matrix_Free(M);

		    QQ c = term[i]->n.coeff[j];
		    c *= gf->term[i2]->n.coeff[j2];
		    pp.add(c, P, MaxRays);

		    Polyhedron_Free(P);
		}
	    }

	    gen_fun *t = pp.compute(MaxRays);
	    sum->add(one, t);
	    delete t;
	}
    }
    Polyhedron_Free(U);
    return sum;
}

void gen_fun::add_union(gen_fun *gf, unsigned MaxRays)
{
    QQ one(1, 1), mone(-1, 1);

    gen_fun *hp = Hadamard_product(gf, MaxRays);
    add(one, gf);
    add(mone, hp);
    delete hp;
}

static void Polyhedron_Shift(Polyhedron *P, Vector *offset)
{
    Value tmp;
    value_init(tmp);
    for (int i = 0; i < P->NbConstraints; ++i) {
	Inner_Product(P->Constraint[i]+1, offset->p, P->Dimension, &tmp);
	value_subtract(P->Constraint[i][1+P->Dimension],
		       P->Constraint[i][1+P->Dimension], tmp);
    }
    for (int i = 0; i < P->NbRays; ++i) {
	if (value_notone_p(P->Ray[i][0]))
	    continue;
	if (value_zero_p(P->Ray[i][1+P->Dimension]))
	    continue;
	Vector_Combine(P->Ray[i]+1, offset->p, P->Ray[i]+1,
		       P->Ray[i][0], P->Ray[i][1+P->Dimension], P->Dimension);
    }
    value_clear(tmp);
}

void gen_fun::shift(const vec_ZZ& offset)
{
    for (int i = 0; i < term.size(); ++i)
	for (int j = 0; j < term[i]->n.power.NumRows(); ++j)
	    term[i]->n.power[j] += offset;

    Vector *v = Vector_Alloc(offset.length());
    zz2values(offset, v->p);
    Polyhedron_Shift(context, v);
    Vector_Free(v);
}

/* Divide the generating functin by 1/(1-z^power).
 * The effect on the corresponding explicit function f(x) is
 * f'(x) = \sum_{i=0}^\infty f(x - i * power)
 */
void gen_fun::divide(const vec_ZZ& power)
{
    for (int i = 0; i < term.size(); ++i) {
	int r = term[i]->d.power.NumRows();
	int c = term[i]->d.power.NumCols();
	term[i]->d.power.SetDims(r+1, c);
	term[i]->d.power[r] = power;
    }

    Vector *v = Vector_Alloc(1+power.length()+1);
    value_set_si(v->p[0], 1);
    zz2values(power, v->p+1);
    Polyhedron *C = AddRays(v->p, 1, context, context->NbConstraints+1);
    Vector_Free(v);
    Polyhedron_Free(context);
    context = C;
}

static void print_power(std::ostream& os, QQ& c, vec_ZZ& p, 
			unsigned int nparam, char **param_name)
{
    bool first = true;

    for (int i = 0; i < p.length(); ++i) {
	if (p[i] == 0)
	    continue;
	if (first) {
	    if (c.n == -1 && c.d == 1)
		os << "-";
	    else if (c.n != 1 || c.d != 1) {
		os << c.n;
		if (c.d != 1)
		    os << " / " << c.d;
		os << "*";
	    }
	    first = false;
	} else
	    os << "*";
	if (i < nparam)
	    os << param_name[i];
	else
	    os << "x" << i;
	if (p[i] == 1)
	    continue;
	if (p[i] < 0)
	    os << "^(" << p[i] << ")";
	else
	    os << "^" << p[i];
    }
    if (first) {
	os << c.n;
	if (c.d != 1)
	    os << " / " << c.d;
    }
}

void gen_fun::print(std::ostream& os, unsigned int nparam, char **param_name) const
{
    QQ mone(-1, 1);
    for (int i = 0; i < term.size(); ++i) {
	if (i != 0)
	    os << " + ";
	os << "(";
	for (int j = 0; j < term[i]->n.coeff.length(); ++j) {
	    if (j != 0 && term[i]->n.coeff[j].n > 0)
		os << "+";
	    print_power(os, term[i]->n.coeff[j], term[i]->n.power[j],
			nparam, param_name);
	}
	os << ")/(";
	for (int j = 0; j < term[i]->d.power.NumRows(); ++j) {
	    if (j != 0)
		os << " * ";
	    os << "(1";
	    print_power(os, mone, term[i]->d.power[j], nparam, param_name);
	    os << ")";
	}
	os << ")";
    }
}

gen_fun::operator evalue *() const
{
    evalue *EP = NULL;
    evalue factor;
    value_init(factor.d);
    value_init(factor.x.n);
    for (int i = 0; i < term.size(); ++i) {
	unsigned nvar = term[i]->d.power.NumRows();
	unsigned nparam = term[i]->d.power.NumCols();
	Matrix *C = Matrix_Alloc(nparam + nvar, 1 + nvar + nparam + 1); 
	mat_ZZ& d = term[i]->d.power;
	Polyhedron *U = context ? context : Universe_Polyhedron(nparam);

	for (int j = 0; j < term[i]->n.coeff.length(); ++j) {
	    for (int r = 0; r < nparam; ++r) {
		value_set_si(C->p[r][0], 0);
		for (int c = 0; c < nvar; ++c) {
		    zz2value(d[c][r], C->p[r][1+c]);
		}
		Vector_Set(&C->p[r][1+nvar], 0, nparam);
		value_set_si(C->p[r][1+nvar+r], -1);
		zz2value(term[i]->n.power[j][r], C->p[r][1+nvar+nparam]);
	    }
	    for (int r = 0; r < nvar; ++r) {
		value_set_si(C->p[nparam+r][0], 1);
		Vector_Set(&C->p[nparam+r][1], 0, nvar + nparam + 1);
		value_set_si(C->p[nparam+r][1+r], 1);
	    }
	    Polyhedron *P = Constraints2Polyhedron(C, 0);
	    evalue *E = barvinok_enumerate_ev(P, U, 0);
	    Polyhedron_Free(P);
	    if (EVALUE_IS_ZERO(*E)) {
		free_evalue_refs(E);
		free(E);
		continue;
	    }
	    zz2value(term[i]->n.coeff[j].n, factor.x.n);
	    zz2value(term[i]->n.coeff[j].d, factor.d);
	    emul(&factor, E);
	    /*
	    Matrix_Print(stdout, P_VALUE_FMT, C);
	    char *test[] = { "A", "B", "C", "D", "E", "F", "G" };
	    print_evalue(stdout, E, test);
	    */
	    if (!EP)
		EP = E;
	    else {
		eadd(E, EP);
		free_evalue_refs(E);
		free(E);
	    }
	}
	Matrix_Free(C);
	if (!context)
	    Polyhedron_Free(U);
    }
    value_clear(factor.d);
    value_clear(factor.x.n);
    return EP;
}

void gen_fun::coefficient(Value* params, Value* c) const
{
    if (context && !in_domain(context, params)) {
	value_set_si(*c, 0);
	return;
    }

    evalue part;
    value_init(part.d);
    value_init(part.x.n);
    evalue sum;
    value_init(sum.d);
    evalue_set_si(&sum, 0, 1);
    Value tmp;
    value_init(tmp);

    for (int i = 0; i < term.size(); ++i) {
	unsigned nvar = term[i]->d.power.NumRows();
	unsigned nparam = term[i]->d.power.NumCols();
	Matrix *C = Matrix_Alloc(nparam + nvar, 1 + nvar + 1); 
	mat_ZZ& d = term[i]->d.power;

	for (int j = 0; j < term[i]->n.coeff.length(); ++j) {
	    C->NbRows = nparam+nvar;
	    for (int r = 0; r < nparam; ++r) {
		value_set_si(C->p[r][0], 0);
		for (int c = 0; c < nvar; ++c) {
		    zz2value(d[c][r], C->p[r][1+c]);
		}
		zz2value(term[i]->n.power[j][r], C->p[r][1+nvar]);
		value_subtract(C->p[r][1+nvar], C->p[r][1+nvar], params[r]);
	    }
	    for (int r = 0; r < nvar; ++r) {
		value_set_si(C->p[nparam+r][0], 1);
		Vector_Set(&C->p[nparam+r][1], 0, nvar + 1);
		value_set_si(C->p[nparam+r][1+r], 1);
	    }
	    Polyhedron *P = Constraints2Polyhedron(C, 0);
	    if (emptyQ(P)) {
		Polyhedron_Free(P);
		continue;
	    }
	    barvinok_count(P, &tmp, 0);
	    Polyhedron_Free(P);
	    if (value_zero_p(tmp))
		continue;
	    zz2value(term[i]->n.coeff[j].n, part.x.n);
	    zz2value(term[i]->n.coeff[j].d, part.d);
	    value_multiply(part.x.n, part.x.n, tmp);
	    eadd(&part, &sum);
	}
	Matrix_Free(C);
    }

    assert(value_one_p(sum.d));
    value_assign(*c, sum.x.n);

    value_clear(tmp);
    value_clear(part.d);
    value_clear(part.x.n);
    value_clear(sum.d);
    value_clear(sum.x.n);
}

gen_fun *gen_fun::summate(int nvar) const
{
    int dim = context->Dimension;
    int nparam = dim - nvar;

#ifdef USE_INCREMENTAL_DF
    partial_ireducer red(Polyhedron_Project(context, nparam), dim, nparam);
#else
    partial_reducer red(Polyhedron_Project(context, nparam), dim, nparam);
#endif
    red.init(context);
    for (int i = 0; i < term.size(); ++i)
	for (int j = 0; j < term[i]->n.power.NumRows(); ++j)
	    red.reduce(term[i]->n.coeff[j], term[i]->n.power[j], term[i]->d.power);
    return red.gf;
}

/* returns true if the set was finite and false otherwise */
bool gen_fun::summate(Value *sum) const
{
    if (term.size() == 0) {
	value_set_si(*sum, 0);
	return true;
    }

    int maxlen = 0;
    for (int i = 0; i < term.size(); ++i)
	if (term[i]->d.power.NumRows() > maxlen)
	    maxlen = term[i]->d.power.NumRows();

    infinite_icounter cnt(term[0]->d.power.NumCols(), maxlen);
    for (int i = 0; i < term.size(); ++i)
	for (int j = 0; j < term[i]->n.power.NumRows(); ++j)
	    cnt.reduce(term[i]->n.coeff[j], term[i]->n.power[j], term[i]->d.power);

    for (int i = 1; i <= maxlen; ++i)
	if (value_notzero_p(mpq_numref(cnt.count[i]))) {
	    value_set_si(*sum, -1);
	    return false;
	}

    assert(value_one_p(mpq_denref(cnt.count[0])));
    value_assign(*sum, mpq_numref(cnt.count[0]));
    return true;
}
