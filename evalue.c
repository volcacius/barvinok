/***********************************************************************/
/*                copyright 1997, Doran Wilde                          */
/*                copyright 1997-2000, Vincent Loechner                */
/*                copyright 2003-2006, Sven Verdoolaege                */
/*       Permission is granted to copy, use, and distribute            */
/*       for any commercial or noncommercial purpose under the terms   */
/*       of the GNU General Public license, version 2, June 1991       */
/*       (see file : LICENSE).                                         */
/***********************************************************************/

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <barvinok/evalue.h>
#include <barvinok/barvinok.h>
#include <barvinok/util.h>

#ifndef value_pmodulus
#define value_pmodulus(ref,val1,val2)  (mpz_fdiv_r((ref),(val1),(val2)))
#endif

#define ALLOC(type) (type*)malloc(sizeof(type))

#ifdef __GNUC__
#define NALLOC(p,n) p = (typeof(p))malloc((n) * sizeof(*p))
#else
#define NALLOC(p,n) p = (void *)malloc((n) * sizeof(*p))
#endif

void evalue_set_si(evalue *ev, int n, int d) {
    value_set_si(ev->d, d);
    value_init(ev->x.n);
    value_set_si(ev->x.n, n);
}

void evalue_set(evalue *ev, Value n, Value d) {
    value_assign(ev->d, d);
    value_init(ev->x.n);
    value_assign(ev->x.n, n);
}

evalue* evalue_zero()
{
    evalue *EP = ALLOC(evalue);
    value_init(EP->d);
    evalue_set_si(EP, 0, 1);
    return EP;
}

void aep_evalue(evalue *e, int *ref) {
  
    enode *p;
    int i;
  
    if (value_notzero_p(e->d))
        return;	        /* a rational number, its already reduced */
    if(!(p = e->x.p))
        return;	        /* hum... an overflow probably occured */
  
    /* First check the components of p */
    for (i=0;i<p->size;i++)
        aep_evalue(&p->arr[i],ref);
  
    /* Then p itself */
    switch (p->type) {
    case polynomial:
    case periodic:
    case evector:
	p->pos = ref[p->pos-1]+1;
    }
    return;
} /* aep_evalue */

/** Comments */
void addeliminatedparams_evalue(evalue *e,Matrix *CT) {
	
    enode *p;
    int i, j;
    int *ref;

    if (value_notzero_p(e->d))
        return;	         /* a rational number, its already reduced */
    if(!(p = e->x.p))
        return;	         /* hum... an overflow probably occured */
  
    /* Compute ref */
    ref = (int *)malloc(sizeof(int)*(CT->NbRows-1));
    for(i=0;i<CT->NbRows-1;i++)
        for(j=0;j<CT->NbColumns;j++)
            if(value_notzero_p(CT->p[i][j])) {
                ref[i] = j;
                break;
            }
  
    /* Transform the references in e, using ref */
    aep_evalue(e,ref);
    free( ref );
    return;
} /* addeliminatedparams_evalue */

void addeliminatedparams_enum(evalue *e, Matrix *CT, Polyhedron *CEq,
			      unsigned MaxRays, unsigned nparam)
{
    enode *p;
    int i;

    if (CT->NbRows == CT->NbColumns)
	return;

    if (EVALUE_IS_ZERO(*e))
	return;

    if (value_notzero_p(e->d)) {
	evalue res;
	value_init(res.d);
	value_set_si(res.d, 0);
	res.x.p = new_enode(partition, 2, nparam);
	EVALUE_SET_DOMAIN(res.x.p->arr[0], 
	    DomainConstraintSimplify(Polyhedron_Copy(CEq), MaxRays));
	value_clear(res.x.p->arr[1].d);
	res.x.p->arr[1] = *e;
	*e = res;
	return;
    }

    p = e->x.p;
    assert(p);
    assert(p->type == partition);
    p->pos = nparam;

    for (i=0; i<p->size/2; i++) {
	Polyhedron *D = EVALUE_DOMAIN(p->arr[2*i]);
	Polyhedron *T = DomainPreimage(D, CT, MaxRays);
	Domain_Free(D);
	D = T;
	T = DomainIntersection(D, CEq, MaxRays);
	Domain_Free(D);
	EVALUE_SET_DOMAIN(p->arr[2*i], T);
	addeliminatedparams_evalue(&p->arr[2*i+1], CT);
    }
}

static int mod_rational_smaller(evalue *e1, evalue *e2)
{
    int r;
    Value m;
    Value m2;
    value_init(m);
    value_init(m2);

    assert(value_notzero_p(e1->d));
    assert(value_notzero_p(e2->d));
    value_multiply(m, e1->x.n, e2->d);
    value_multiply(m2, e2->x.n, e1->d);
    if (value_lt(m, m2))
	r = 1;
    else if (value_gt(m, m2))
	r = 0;
    else 
	r = -1;
    value_clear(m);
    value_clear(m2);

    return r;
}

static int mod_term_smaller_r(evalue *e1, evalue *e2)
{
    if (value_notzero_p(e1->d)) {
	int r;
	if (value_zero_p(e2->d))
	    return 1;
	r = mod_rational_smaller(e1, e2);
	return r == -1 ? 0 : r;
    }
    if (value_notzero_p(e2->d))
	return 0;
    if (e1->x.p->pos < e2->x.p->pos)
	return 1;
    else if (e1->x.p->pos > e2->x.p->pos)
	return 0;
    else {
	int r = mod_rational_smaller(&e1->x.p->arr[1], &e2->x.p->arr[1]);
	return r == -1 
		 ? mod_term_smaller_r(&e1->x.p->arr[0], &e2->x.p->arr[0])
		 : r;
    }
}

static int mod_term_smaller(const evalue *e1, const evalue *e2)
{
    assert(value_zero_p(e1->d));
    assert(value_zero_p(e2->d));
    assert(e1->x.p->type == fractional || e1->x.p->type == flooring);
    assert(e2->x.p->type == fractional || e2->x.p->type == flooring);
    return mod_term_smaller_r(&e1->x.p->arr[0], &e2->x.p->arr[0]);
}

static void check_order(const evalue *e)
{
    int i;
    evalue *a;

    if (value_notzero_p(e->d))
	return;

    switch (e->x.p->type) {
    case partition:
	for (i = 0; i < e->x.p->size/2; ++i)
	    check_order(&e->x.p->arr[2*i+1]);
	break;
    case relation:
	for (i = 1; i < e->x.p->size; ++i) {
	    a = &e->x.p->arr[i];
	    if (value_notzero_p(a->d))
		continue;
	    switch (a->x.p->type) {
	    case relation:
		assert(mod_term_smaller(&e->x.p->arr[0], &a->x.p->arr[0]));
		break;
	    case partition:
		assert(0);
	    }
	    check_order(a);
	}
	break;
    case polynomial:
	for (i = 0; i < e->x.p->size; ++i) {
	    a = &e->x.p->arr[i];
	    if (value_notzero_p(a->d))
		continue;
	    switch (a->x.p->type) {
	    case polynomial:
		assert(e->x.p->pos < a->x.p->pos);
		break;
	    case relation:
	    case partition:
		assert(0);
	    }
	    check_order(a);
	}
	break;
    case fractional:
    case flooring:
	for (i = 1; i < e->x.p->size; ++i) {
	    a = &e->x.p->arr[i];
	    if (value_notzero_p(a->d))
		continue;
	    switch (a->x.p->type) {
	    case polynomial:
	    case relation:
	    case partition:
		assert(0);
	    }
	}
	break;
    }
}

/* Negative pos means inequality */
/* s is negative of substitution if m is not zero */
struct fixed_param {
    int	    pos;
    evalue  s;
    Value   d;
    Value   m;
};

struct subst {
    struct fixed_param *fixed;
    int			n;
    int			max;
};

static int relations_depth(evalue *e)
{
    int d;

    for (d = 0; 
	 value_zero_p(e->d) && e->x.p->type == relation;
	 e = &e->x.p->arr[1], ++d);
    return d;
}

static void poly_denom_not_constant(evalue **pp, Value *d)
{
    evalue *p = *pp;
    value_set_si(*d, 1);

    while (value_zero_p(p->d)) {
	assert(p->x.p->type == polynomial);
	assert(p->x.p->size == 2);
	assert(value_notzero_p(p->x.p->arr[1].d));
	value_lcm(*d, p->x.p->arr[1].d, d);
	p = &p->x.p->arr[0];
    }
    *pp = p;
}

static void poly_denom(evalue *p, Value *d)
{
    poly_denom_not_constant(&p, d);
    value_lcm(*d, p->d, d);
}

static void realloc_substitution(struct subst *s, int d)
{
    struct fixed_param *n;
    int i;
    NALLOC(n, s->max+d);
    for (i = 0; i < s->n; ++i)
	n[i] = s->fixed[i];
    free(s->fixed);
    s->fixed = n;
    s->max += d;
}

static int add_modulo_substitution(struct subst *s, evalue *r)
{
    evalue *p;
    evalue *f;
    evalue *m;

    assert(value_zero_p(r->d) && r->x.p->type == relation);
    m = &r->x.p->arr[0];

    /* May have been reduced already */
    if (value_notzero_p(m->d))
	return 0;

    assert(value_zero_p(m->d) && m->x.p->type == fractional);
    assert(m->x.p->size == 3);

    /* fractional was inverted during reduction
     * invert it back and move constant in
     */
    if (!EVALUE_IS_ONE(m->x.p->arr[2])) {
	assert(value_pos_p(m->x.p->arr[2].d));
	assert(value_mone_p(m->x.p->arr[2].x.n));
	value_set_si(m->x.p->arr[2].x.n, 1);
	value_increment(m->x.p->arr[1].x.n, m->x.p->arr[1].x.n);
	assert(value_eq(m->x.p->arr[1].x.n, m->x.p->arr[1].d));
	value_set_si(m->x.p->arr[1].x.n, 1);
	eadd(&m->x.p->arr[1], &m->x.p->arr[0]);
	value_set_si(m->x.p->arr[1].x.n, 0);
	value_set_si(m->x.p->arr[1].d, 1);
    }

    /* Oops.  Nested identical relations. */
    if (!EVALUE_IS_ZERO(m->x.p->arr[1]))
	return 0;

    if (s->n >= s->max) {
	int d = relations_depth(r);
	realloc_substitution(s, d);
    }

    p = &m->x.p->arr[0];
    assert(value_zero_p(p->d) && p->x.p->type == polynomial);
    assert(p->x.p->size == 2);
    f = &p->x.p->arr[1];

    assert(value_pos_p(f->x.n));

    value_init(s->fixed[s->n].m);
    value_assign(s->fixed[s->n].m, f->d);
    s->fixed[s->n].pos = p->x.p->pos;
    value_init(s->fixed[s->n].d);
    value_assign(s->fixed[s->n].d, f->x.n);
    value_init(s->fixed[s->n].s.d);
    evalue_copy(&s->fixed[s->n].s, &p->x.p->arr[0]);
    ++s->n;

    return 1;
}

static int type_offset(enode *p)
{
   return p->type == fractional ? 1 : 
	  p->type == flooring ? 1 : 0;
}

static void reorder_terms_about(enode *p, evalue *v)
{
    int i;
    int offset = type_offset(p);

    for (i = p->size-1; i >= offset+1; i--) {
	emul(v, &p->arr[i]);
	eadd(&p->arr[i], &p->arr[i-1]);
	free_evalue_refs(&(p->arr[i]));
    }
    p->size = offset+1;
    free_evalue_refs(v);
}

static void reorder_terms(evalue *e)
{
    enode *p;
    evalue f;

    assert(value_zero_p(e->d));
    p = e->x.p;
    assert(p->type = fractional);  /* for now */

    value_init(f.d);
    value_set_si(f.d, 0);
    f.x.p = new_enode(fractional, 3, -1);
    value_clear(f.x.p->arr[0].d);
    f.x.p->arr[0] = p->arr[0];
    evalue_set_si(&f.x.p->arr[1], 0, 1);
    evalue_set_si(&f.x.p->arr[2], 1, 1);
    reorder_terms_about(p, &f);
    value_clear(e->d);
    *e = p->arr[1];
    free(p);
}

void _reduce_evalue (evalue *e, struct subst *s, int fract) {
  
    enode *p;
    int i, j, k;
    int add;
  
    if (value_notzero_p(e->d)) {
	if (fract)
	    mpz_fdiv_r(e->x.n, e->x.n, e->d);
        return;	/* a rational number, its already reduced */
    }

    if(!(p = e->x.p))
        return;	/* hum... an overflow probably occured */
  
    /* First reduce the components of p */
    add = p->type == relation;
    for (i=0; i<p->size; i++) {
	if (add && i == 1)
	    add = add_modulo_substitution(s, e);

        if (i == 0 && p->type==fractional)
	    _reduce_evalue(&p->arr[i], s, 1);
	else
	    _reduce_evalue(&p->arr[i], s, fract);

	if (add && i == p->size-1) {
	    --s->n;
	    value_clear(s->fixed[s->n].m);
	    value_clear(s->fixed[s->n].d);
	    free_evalue_refs(&s->fixed[s->n].s); 
	} else if (add && i == 1)
	    s->fixed[s->n-1].pos *= -1;
    }

    if (p->type==periodic) {
    
        /* Try to reduce the period */
        for (i=1; i<=(p->size)/2; i++) {
            if ((p->size % i)==0) {
	
                /* Can we reduce the size to i ? */
                for (j=0; j<i; j++)
                    for (k=j+i; k<e->x.p->size; k+=i)
                        if (!eequal(&p->arr[j], &p->arr[k])) goto you_lose;

                /* OK, lets do it */
                for (j=i; j<p->size; j++) free_evalue_refs(&p->arr[j]);
                p->size = i;
                break;

you_lose:   	/* OK, lets not do it */
                continue;
            }
        }

        /* Try to reduce its strength */
        if (p->size == 1) {
	    value_clear(e->d);
            memcpy(e,&p->arr[0],sizeof(evalue));
            free(p);
        }
    }
    else if (p->type==polynomial) {
	for (k = 0; s && k < s->n; ++k) {
	    if (s->fixed[k].pos == p->pos) {
		int divide = value_notone_p(s->fixed[k].d);
		evalue d;

		if (value_notzero_p(s->fixed[k].m)) {
		    if (!fract)
			continue;
		    assert(p->size == 2);
		    if (divide && value_ne(s->fixed[k].d, p->arr[1].x.n))
			continue;
		    if (!mpz_divisible_p(s->fixed[k].m, p->arr[1].d))
			continue;
		    divide = 1;
		}

		if (divide) {
		    value_init(d.d);
		    value_assign(d.d, s->fixed[k].d);
		    value_init(d.x.n);
		    if (value_notzero_p(s->fixed[k].m))
			value_oppose(d.x.n, s->fixed[k].m);
		    else
			value_set_si(d.x.n, 1);
		}

		for (i=p->size-1;i>=1;i--) {
		    emul(&s->fixed[k].s, &p->arr[i]);
		    if (divide)
			emul(&d, &p->arr[i]);
		    eadd(&p->arr[i], &p->arr[i-1]);
		    free_evalue_refs(&(p->arr[i]));
		}
		p->size = 1;
		_reduce_evalue(&p->arr[0], s, fract);

		if (divide)
		    free_evalue_refs(&d);

		break;
	    }
	}

        /* Try to reduce the degree */
        for (i=p->size-1;i>=1;i--) {
            if (!(value_one_p(p->arr[i].d) && value_zero_p(p->arr[i].x.n)))
		break;
	    /* Zero coefficient */
	    free_evalue_refs(&(p->arr[i]));
        }
        if (i+1<p->size)
	    p->size = i+1;

        /* Try to reduce its strength */
        if (p->size == 1) {
	    value_clear(e->d);
            memcpy(e,&p->arr[0],sizeof(evalue));
            free(p);
        }
    }
    else if (p->type==fractional) {
	int reorder = 0;
	evalue v;

	if (value_notzero_p(p->arr[0].d)) {
	    value_init(v.d);
	    value_assign(v.d, p->arr[0].d);
	    value_init(v.x.n);
	    mpz_fdiv_r(v.x.n, p->arr[0].x.n,  p->arr[0].d);

	    reorder = 1;
	} else {
	    evalue *f, *base;
	    evalue *pp = &p->arr[0];
	    assert(value_zero_p(pp->d) && pp->x.p->type == polynomial);
	    assert(pp->x.p->size == 2);

	    /* search for exact duplicate among the modulo inequalities */
	    do {
		f = &pp->x.p->arr[1];
		for (k = 0; s && k < s->n; ++k) {
		    if (-s->fixed[k].pos == pp->x.p->pos &&
			    value_eq(s->fixed[k].d, f->x.n) &&
			    value_eq(s->fixed[k].m, f->d) &&
			    eequal(&s->fixed[k].s, &pp->x.p->arr[0]))
			break;
		}
		if (k < s->n) {
		    Value g;
		    value_init(g);

		    /* replace { E/m } by { (E-1)/m } + 1/m */
		    poly_denom(pp, &g);
		    if (reorder) {
			evalue extra;
			value_init(extra.d);
			evalue_set_si(&extra, 1, 1);
			value_assign(extra.d, g);
			eadd(&extra, &v.x.p->arr[1]);
			free_evalue_refs(&extra); 

			/* We've been going in circles; stop now */
			if (value_ge(v.x.p->arr[1].x.n, v.x.p->arr[1].d)) {
			    free_evalue_refs(&v);
			    value_init(v.d);
			    evalue_set_si(&v, 0, 1);
			    break;
			}
		    } else {
			value_init(v.d);
			value_set_si(v.d, 0);
			v.x.p = new_enode(fractional, 3, -1);
			evalue_set_si(&v.x.p->arr[1], 1, 1);
			value_assign(v.x.p->arr[1].d, g);
			evalue_set_si(&v.x.p->arr[2], 1, 1);
			evalue_copy(&v.x.p->arr[0], &p->arr[0]);
		    }

		    for (f = &v.x.p->arr[0]; value_zero_p(f->d); 
					     f = &f->x.p->arr[0])
			;
		    value_division(f->d, g, f->d);
		    value_multiply(f->x.n, f->x.n, f->d);
		    value_assign(f->d, g);
		    value_decrement(f->x.n, f->x.n);
		    mpz_fdiv_r(f->x.n, f->x.n, f->d);

		    Gcd(f->d, f->x.n, &g);
		    value_division(f->d, f->d, g);
		    value_division(f->x.n, f->x.n, g);

		    value_clear(g);
		    pp = &v.x.p->arr[0];

		    reorder = 1;
		}
	    } while (k < s->n);

	    /* reduction may have made this fractional arg smaller */
	    i = reorder ? p->size : 1;
	    for ( ; i < p->size; ++i)
		if (value_zero_p(p->arr[i].d) && 
			p->arr[i].x.p->type == fractional &&
			!mod_term_smaller(e, &p->arr[i]))
		    break;
	    if (i < p->size) {
		value_init(v.d);
		value_set_si(v.d, 0);
		v.x.p = new_enode(fractional, 3, -1);
		evalue_set_si(&v.x.p->arr[1], 0, 1);
		evalue_set_si(&v.x.p->arr[2], 1, 1);
		evalue_copy(&v.x.p->arr[0], &p->arr[0]);

		reorder = 1;
	    }

	    if (!reorder) {
		Value m;
		Value r;
		evalue *pp = &p->arr[0];
		value_init(m);
		value_init(r);
		poly_denom_not_constant(&pp, &m);
		mpz_fdiv_r(r, m, pp->d);
		if (value_notzero_p(r)) {
		    value_init(v.d);
		    value_set_si(v.d, 0);
		    v.x.p = new_enode(fractional, 3, -1);

		    value_multiply(r, m, pp->x.n);
		    value_multiply(v.x.p->arr[1].d, m, pp->d);
		    value_init(v.x.p->arr[1].x.n);
		    mpz_fdiv_r(v.x.p->arr[1].x.n, r, pp->d);
		    mpz_fdiv_q(r, r, pp->d);

		    evalue_set_si(&v.x.p->arr[2], 1, 1);
		    evalue_copy(&v.x.p->arr[0], &p->arr[0]);
		    pp = &v.x.p->arr[0];
		    while (value_zero_p(pp->d))
			pp = &pp->x.p->arr[0];

		    value_assign(pp->d, m);
		    value_assign(pp->x.n, r);

		    Gcd(pp->d, pp->x.n, &r);
		    value_division(pp->d, pp->d, r);
		    value_division(pp->x.n, pp->x.n, r);

		    reorder = 1;
		}
		value_clear(m);
		value_clear(r);
	    }

	    if (!reorder) {
		int invert = 0;
		Value twice;
		value_init(twice);

		for (pp = &p->arr[0]; value_zero_p(pp->d); 
				      pp = &pp->x.p->arr[0]) {
		    f = &pp->x.p->arr[1];
		    assert(value_pos_p(f->d));
		    mpz_mul_ui(twice, f->x.n, 2);
		    if (value_lt(twice, f->d))
			break;
		    if (value_eq(twice, f->d))
			continue;
		    invert = 1;
		    break;
		}

		if (invert) {
		    value_init(v.d);
		    value_set_si(v.d, 0);
		    v.x.p = new_enode(fractional, 3, -1);
		    evalue_set_si(&v.x.p->arr[1], 0, 1);
		    poly_denom(&p->arr[0], &twice);
		    value_assign(v.x.p->arr[1].d, twice);
		    value_decrement(v.x.p->arr[1].x.n, twice);
		    evalue_set_si(&v.x.p->arr[2], -1, 1);
		    evalue_copy(&v.x.p->arr[0], &p->arr[0]);

		    for (pp = &v.x.p->arr[0]; value_zero_p(pp->d); 
					      pp = &pp->x.p->arr[0]) {
			f = &pp->x.p->arr[1];
			value_oppose(f->x.n, f->x.n);
			mpz_fdiv_r(f->x.n, f->x.n,  f->d);
		    }
		    value_division(pp->d, twice, pp->d);
		    value_multiply(pp->x.n, pp->x.n, pp->d);
		    value_assign(pp->d, twice);
		    value_oppose(pp->x.n, pp->x.n);
		    value_decrement(pp->x.n, pp->x.n);
		    mpz_fdiv_r(pp->x.n, pp->x.n, pp->d);

		    /* Maybe we should do this during reduction of 
		     * the constant.
		     */
		    Gcd(pp->d, pp->x.n, &twice);
		    value_division(pp->d, pp->d, twice);
		    value_division(pp->x.n, pp->x.n, twice);

		    reorder = 1;
		}

		value_clear(twice);
	    }
	}

	if (reorder) {
	    reorder_terms_about(p, &v);
	    _reduce_evalue(&p->arr[1], s, fract);
	}

        /* Try to reduce the degree */
        for (i=p->size-1;i>=2;i--) {
            if (!(value_one_p(p->arr[i].d) && value_zero_p(p->arr[i].x.n)))
		break;
	    /* Zero coefficient */
	    free_evalue_refs(&(p->arr[i]));
        }
        if (i+1<p->size)
	    p->size = i+1;

        /* Try to reduce its strength */
        if (p->size == 2) {
	    value_clear(e->d);
            memcpy(e,&p->arr[1],sizeof(evalue));
	    free_evalue_refs(&(p->arr[0]));
            free(p);
        }
    }
    else if (p->type == flooring) {
        /* Try to reduce the degree */
        for (i=p->size-1;i>=2;i--) {
	    if (!EVALUE_IS_ZERO(p->arr[i]))
		break;
	    /* Zero coefficient */
	    free_evalue_refs(&(p->arr[i]));
        }
        if (i+1<p->size)
	    p->size = i+1;

        /* Try to reduce its strength */
        if (p->size == 2) {
	    value_clear(e->d);
            memcpy(e,&p->arr[1],sizeof(evalue));
	    free_evalue_refs(&(p->arr[0]));
            free(p);
        }
    }
    else if (p->type == relation) {
	if (p->size == 3 && eequal(&p->arr[1], &p->arr[2])) {
	    free_evalue_refs(&(p->arr[2]));
	    free_evalue_refs(&(p->arr[0]));
	    p->size = 2;
	    value_clear(e->d);
	    *e = p->arr[1];
	    free(p);
	    return;
	}
	if (p->size == 3 && EVALUE_IS_ZERO(p->arr[2])) {
	    free_evalue_refs(&(p->arr[2]));
	    p->size = 2;
	}
	if (p->size == 2 && EVALUE_IS_ZERO(p->arr[1])) {
	    free_evalue_refs(&(p->arr[1]));
	    free_evalue_refs(&(p->arr[0]));
	    evalue_set_si(e, 0, 1);
	    free(p);
	} else {
	    int reduced = 0;
	    evalue *m;
	    m = &p->arr[0];

	    /* Relation was reduced by means of an identical 
	     * inequality => remove 
	     */
	    if (value_zero_p(m->d) && !EVALUE_IS_ZERO(m->x.p->arr[1]))
		reduced = 1;

	    if (reduced || value_notzero_p(p->arr[0].d)) {
		if (!reduced && value_zero_p(p->arr[0].x.n)) {
		    value_clear(e->d);
		    memcpy(e,&p->arr[1],sizeof(evalue));
		    if (p->size == 3)
			free_evalue_refs(&(p->arr[2]));
		} else {
		    if (p->size == 3) {
			value_clear(e->d);
			memcpy(e,&p->arr[2],sizeof(evalue));
		    } else
			evalue_set_si(e, 0, 1);
		    free_evalue_refs(&(p->arr[1]));
		}
		free_evalue_refs(&(p->arr[0]));
		free(p);
	    }
	}
    }
    return;
} /* reduce_evalue */

static void add_substitution(struct subst *s, Value *row, unsigned dim)
{
    int k, l;
    evalue *r;

    for (k = 0; k < dim; ++k)
	if (value_notzero_p(row[k+1]))
	    break;

    Vector_Normalize_Positive(row+1, dim+1, k);
    assert(s->n < s->max);
    value_init(s->fixed[s->n].d);
    value_init(s->fixed[s->n].m);
    value_assign(s->fixed[s->n].d, row[k+1]);
    s->fixed[s->n].pos = k+1;
    value_set_si(s->fixed[s->n].m, 0);
    r = &s->fixed[s->n].s;
    value_init(r->d);
    for (l = k+1; l < dim; ++l)
	if (value_notzero_p(row[l+1])) {
	    value_set_si(r->d, 0);
	    r->x.p = new_enode(polynomial, 2, l + 1);
	    value_init(r->x.p->arr[1].x.n);
	    value_oppose(r->x.p->arr[1].x.n, row[l+1]);
	    value_set_si(r->x.p->arr[1].d, 1);
	    r = &r->x.p->arr[0];
	}
    value_init(r->x.n);
    value_oppose(r->x.n, row[dim+1]);
    value_set_si(r->d, 1);
    ++s->n;
}

static void _reduce_evalue_in_domain(evalue *e, Polyhedron *D, struct subst *s)
{
    unsigned dim;
    Polyhedron *orig = D;

    s->n = 0;
    dim = D->Dimension;
    if (D->next)
	D = DomainConvex(D, 0);
    if (!D->next && D->NbEq) {
	int j, k;
	if (s->max < dim) {
	    if (s->max != 0)
		realloc_substitution(s, dim);
	    else {
		int d = relations_depth(e);
		s->max = dim+d;
		NALLOC(s->fixed, s->max);
	    }
	}
	for (j = 0; j < D->NbEq; ++j)
	    add_substitution(s, D->Constraint[j], dim);
    }
    if (D != orig)
	Domain_Free(D);
    _reduce_evalue(e, s, 0);
    if (s->n != 0) {
	int j;
	for (j = 0; j < s->n; ++j) {
	    value_clear(s->fixed[j].d);
	    value_clear(s->fixed[j].m);
	    free_evalue_refs(&s->fixed[j].s); 
	}
    }
}

void reduce_evalue_in_domain(evalue *e, Polyhedron *D)
{
    struct subst s = { NULL, 0, 0 };
    if (emptyQ2(D)) {
	if (EVALUE_IS_ZERO(*e))
	    return;
	free_evalue_refs(e);
	value_init(e->d);
	evalue_set_si(e, 0, 1);
	return;
    }
    _reduce_evalue_in_domain(e, D, &s);
    if (s.max != 0)
	free(s.fixed);
}

void reduce_evalue (evalue *e) {
    struct subst s = { NULL, 0, 0 };

    if (value_notzero_p(e->d))
        return;	/* a rational number, its already reduced */

    if (e->x.p->type == partition) {
	int i;
	unsigned dim = -1;
	for (i = 0; i < e->x.p->size/2; ++i) {
	    Polyhedron *D = EVALUE_DOMAIN(e->x.p->arr[2*i]);

	    /* This shouldn't really happen; 
	     * Empty domains should not be added.
	     */
	    POL_ENSURE_VERTICES(D);
	    if (!emptyQ(D))
		_reduce_evalue_in_domain(&e->x.p->arr[2*i+1], D, &s);

	    if (emptyQ(D) || EVALUE_IS_ZERO(e->x.p->arr[2*i+1])) {
		free_evalue_refs(&e->x.p->arr[2*i+1]);
		Domain_Free(EVALUE_DOMAIN(e->x.p->arr[2*i]));
		value_clear(e->x.p->arr[2*i].d);
		e->x.p->size -= 2;
		e->x.p->arr[2*i] = e->x.p->arr[e->x.p->size];
		e->x.p->arr[2*i+1] = e->x.p->arr[e->x.p->size+1];
		--i;
	    }
	}
	if (e->x.p->size == 0) {
	    free(e->x.p);
	    evalue_set_si(e, 0, 1);
	}
    } else
	_reduce_evalue(e, &s, 0);
    if (s.max != 0)
	free(s.fixed);
}

void print_evalue(FILE *DST, const evalue *e, char **pname)
{
  if(value_notzero_p(e->d)) {    
    if(value_notone_p(e->d)) {
      value_print(DST,VALUE_FMT,e->x.n);
      fprintf(DST,"/");
      value_print(DST,VALUE_FMT,e->d);
    }  
    else {
      value_print(DST,VALUE_FMT,e->x.n);
    }
  }  
  else
    print_enode(DST,e->x.p,pname);
  return;
} /* print_evalue */

void print_enode(FILE *DST,enode *p,char **pname) {
  
  int i;
  
  if (!p) {
    fprintf(DST, "NULL");
    return;
  }
  switch (p->type) {
  case evector:
    fprintf(DST, "{ ");
    for (i=0; i<p->size; i++) {
      print_evalue(DST, &p->arr[i], pname);
      if (i!=(p->size-1))
	fprintf(DST, ", ");
    }
    fprintf(DST, " }\n");
    break;
  case polynomial:
    fprintf(DST, "( ");
    for (i=p->size-1; i>=0; i--) {
      print_evalue(DST, &p->arr[i], pname);
      if (i==1) fprintf(DST, " * %s + ", pname[p->pos-1]);
      else if (i>1) 
	fprintf(DST, " * %s^%d + ", pname[p->pos-1], i);
    }
    fprintf(DST, " )\n");
    break;
  case periodic:
    fprintf(DST, "[ ");
    for (i=0; i<p->size; i++) {
      print_evalue(DST, &p->arr[i], pname);
      if (i!=(p->size-1)) fprintf(DST, ", ");
    }
    fprintf(DST," ]_%s", pname[p->pos-1]);
    break;
  case flooring:
  case fractional:
    fprintf(DST, "( ");
    for (i=p->size-1; i>=1; i--) {
      print_evalue(DST, &p->arr[i], pname);
      if (i >= 2) {
        fprintf(DST, " * ");
	fprintf(DST, p->type == flooring ? "[" : "{");
        print_evalue(DST, &p->arr[0], pname);
	fprintf(DST, p->type == flooring ? "]" : "}");
	if (i>2) 
	  fprintf(DST, "^%d + ", i-1);
	else
	  fprintf(DST, " + ");
      }
    }
    fprintf(DST, " )\n");
    break;
  case relation:
    fprintf(DST, "[ ");
    print_evalue(DST, &p->arr[0], pname);
    fprintf(DST, "= 0 ] * \n");
    print_evalue(DST, &p->arr[1], pname);
    if (p->size > 2) {
	fprintf(DST, " +\n [ ");
	print_evalue(DST, &p->arr[0], pname);
	fprintf(DST, "!= 0 ] * \n");
	print_evalue(DST, &p->arr[2], pname);
    }
    break;
  case partition: {
    char **names = pname;
    int maxdim = EVALUE_DOMAIN(p->arr[0])->Dimension;
    if (!pname || p->pos < maxdim) {
	NALLOC(names, maxdim);
	for (i = 0; i < p->pos; ++i) {
	    if (pname)
		names[i] = pname[i];
	    else {
		NALLOC(names[i], 10);
		snprintf(names[i], 10, "%c", 'P'+i);
	    }
	}
	for ( ; i < maxdim; ++i) {
	    NALLOC(names[i], 10);
	    snprintf(names[i], 10, "_p%d", i);
	}
    }

    for (i=0; i<p->size/2; i++) {
	Print_Domain(DST, EVALUE_DOMAIN(p->arr[2*i]), names);
	print_evalue(DST, &p->arr[2*i+1], names);
    }

    if (!pname || p->pos < maxdim) {
	for (i = pname ? p->pos : 0; i < maxdim; ++i)
	    free(names[i]);
	free(names);
    }

    break;
  }
  default:
    assert(0);
  }
  return;
} /* print_enode */ 

static void eadd_rev(const evalue *e1, evalue *res)
{
    evalue ev;
    value_init(ev.d);
    evalue_copy(&ev, e1);
    eadd(res, &ev);
    free_evalue_refs(res);	  
    *res = ev;
}

static void eadd_rev_cst(const evalue *e1, evalue *res)
{
    evalue ev;
    value_init(ev.d);
    evalue_copy(&ev, e1);
    eadd(res, &ev.x.p->arr[type_offset(ev.x.p)]);
    free_evalue_refs(res);	  
    *res = ev;
}

static int is_zero_on(evalue *e, Polyhedron *D)
{
    int is_zero;
    evalue tmp;
    value_init(tmp.d);
    tmp.x.p = new_enode(partition, 2, D->Dimension);
    EVALUE_SET_DOMAIN(tmp.x.p->arr[0], Domain_Copy(D));
    evalue_copy(&tmp.x.p->arr[1], e);
    reduce_evalue(&tmp);
    is_zero = EVALUE_IS_ZERO(tmp);
    free_evalue_refs(&tmp);
    return is_zero;
}

struct section { Polyhedron * D; evalue E; };

void eadd_partitions(const evalue *e1, evalue *res)
{
    int n, i, j;
    Polyhedron *d, *fd;
    struct section *s;
    s = (struct section *) 
	    malloc((e1->x.p->size/2+1) * (res->x.p->size/2+1) * 
		   sizeof(struct section));
    assert(s);
    assert(e1->x.p->pos == res->x.p->pos);
    assert(e1->x.p->pos == EVALUE_DOMAIN(e1->x.p->arr[0])->Dimension);
    assert(res->x.p->pos == EVALUE_DOMAIN(res->x.p->arr[0])->Dimension);

    n = 0;
    for (j = 0; j < e1->x.p->size/2; ++j) {
	assert(res->x.p->size >= 2);
	fd = DomainDifference(EVALUE_DOMAIN(e1->x.p->arr[2*j]),
			      EVALUE_DOMAIN(res->x.p->arr[0]), 0);
	if (!emptyQ(fd))
	    for (i = 1; i < res->x.p->size/2; ++i) {
		Polyhedron *t = fd;
		fd = DomainDifference(fd, EVALUE_DOMAIN(res->x.p->arr[2*i]), 0);
		Domain_Free(t);
		if (emptyQ(fd))
		    break;
	    }
	if (emptyQ(fd)) {
	    Domain_Free(fd);
	    continue;
	}
	/* See if we can extend one of the domains in res to cover fd */
	for (i = 0; i < res->x.p->size/2; ++i)
	    if (is_zero_on(&res->x.p->arr[2*i+1], fd))
		break;
	if (i < res->x.p->size/2) {
	    EVALUE_SET_DOMAIN(res->x.p->arr[2*i], 
		      DomainConcat(fd, EVALUE_DOMAIN(res->x.p->arr[2*i])));
	    continue;
	}
	value_init(s[n].E.d);
	evalue_copy(&s[n].E, &e1->x.p->arr[2*j+1]);
	s[n].D = fd;
	++n;
    }
    for (i = 0; i < res->x.p->size/2; ++i) {
	fd = EVALUE_DOMAIN(res->x.p->arr[2*i]);
	for (j = 0; j < e1->x.p->size/2; ++j) {
	    Polyhedron *t;
	    d = DomainIntersection(EVALUE_DOMAIN(e1->x.p->arr[2*j]),
				   EVALUE_DOMAIN(res->x.p->arr[2*i]), 0);
	    if (emptyQ(d)) {
		Domain_Free(d);
		continue;
	    }
	    t = fd;
	    fd = DomainDifference(fd, EVALUE_DOMAIN(e1->x.p->arr[2*j]), 0);
	    if (t != EVALUE_DOMAIN(res->x.p->arr[2*i]))
		Domain_Free(t);
	    value_init(s[n].E.d);
	    evalue_copy(&s[n].E, &res->x.p->arr[2*i+1]);
	    eadd(&e1->x.p->arr[2*j+1], &s[n].E);
	    if (!emptyQ(fd) && is_zero_on(&e1->x.p->arr[2*j+1], fd)) {
		d = DomainConcat(fd, d);
		fd = Empty_Polyhedron(fd->Dimension);
	    }
	    s[n].D = d;
	    ++n;
	}
	if (!emptyQ(fd)) {
	    s[n].E = res->x.p->arr[2*i+1];
	    s[n].D = fd;
	    ++n;
	} else {
	    free_evalue_refs(&res->x.p->arr[2*i+1]);
	    Domain_Free(fd);
	}
	if (fd != EVALUE_DOMAIN(res->x.p->arr[2*i]))
	    Domain_Free(EVALUE_DOMAIN(res->x.p->arr[2*i]));
	value_clear(res->x.p->arr[2*i].d);
    }

    free(res->x.p);
    assert(n > 0);
    res->x.p = new_enode(partition, 2*n, e1->x.p->pos);
    for (j = 0; j < n; ++j) {
	s[j].D = DomainConstraintSimplify(s[j].D, 0);
	EVALUE_SET_DOMAIN(res->x.p->arr[2*j], s[j].D);
	value_clear(res->x.p->arr[2*j+1].d);
	res->x.p->arr[2*j+1] = s[j].E;
    }

    free(s);
}

static void explicit_complement(evalue *res)
{
    enode *rel = new_enode(relation, 3, 0);
    assert(rel);
    value_clear(rel->arr[0].d);
    rel->arr[0] = res->x.p->arr[0];
    value_clear(rel->arr[1].d);
    rel->arr[1] = res->x.p->arr[1];
    value_set_si(rel->arr[2].d, 1);
    value_init(rel->arr[2].x.n);
    value_set_si(rel->arr[2].x.n, 0);
    free(res->x.p);
    res->x.p = rel;
}

void eadd(const evalue *e1, evalue *res)
{
 int i; 
    if (value_notzero_p(e1->d) && value_notzero_p(res->d)) {
         /* Add two rational numbers */
	 Value g,m1,m2;
	 value_init(g);
	 value_init(m1);
	 value_init(m2);
	     
         value_multiply(m1,e1->x.n,res->d);
         value_multiply(m2,res->x.n,e1->d);
         value_addto(res->x.n,m1,m2);
         value_multiply(res->d,e1->d,res->d);
         Gcd(res->x.n,res->d,&g);
         if (value_notone_p(g)) {
	      value_division(res->d,res->d,g);
              value_division(res->x.n,res->x.n,g);
         }
         value_clear(g); value_clear(m1); value_clear(m2);
         return ;
     }
     else if (value_notzero_p(e1->d) && value_zero_p(res->d)) {
	  switch (res->x.p->type) {
	  case polynomial:
	      /* Add the constant to the constant term of a polynomial*/
	       eadd(e1, &res->x.p->arr[0]);
	       return ;
	  case periodic:
	      /* Add the constant to all elements of a periodic number */
	      for (i=0; i<res->x.p->size; i++) {
		  eadd(e1, &res->x.p->arr[i]);
	      }
	      return ;
	  case evector:
	      fprintf(stderr, "eadd: cannot add const with vector\n");
	      return;
	  case flooring:
	  case fractional:
	       eadd(e1, &res->x.p->arr[1]);
	       return ;
	  case partition:
		assert(EVALUE_IS_ZERO(*e1));
		break;				/* Do nothing */
	  case relation:
		/* Create (zero) complement if needed */
		if (res->x.p->size < 3 && !EVALUE_IS_ZERO(*e1))
		    explicit_complement(res);
		for (i = 1; i < res->x.p->size; ++i)
		    eadd(e1, &res->x.p->arr[i]);
		break;
	  default:
		assert(0);
	  }
     }
     /* add polynomial or periodic to constant 
      * you have to exchange e1 and res, before doing addition */
     
     else if (value_zero_p(e1->d) && value_notzero_p(res->d)) {
	  eadd_rev(e1, res); 
	  return;
     }
     else {   // ((e1->d==0) && (res->d==0)) 
	assert(!((e1->x.p->type == partition) ^
	         (res->x.p->type == partition)));
	if (e1->x.p->type == partition) {
	    eadd_partitions(e1, res);
	    return;
	}
	if (e1->x.p->type == relation &&
	    (res->x.p->type != relation || 
	     mod_term_smaller(&e1->x.p->arr[0], &res->x.p->arr[0]))) {
		eadd_rev(e1, res);
		return;
	}
	if (res->x.p->type == relation) {
	    if (e1->x.p->type == relation &&
		eequal(&e1->x.p->arr[0], &res->x.p->arr[0])) {
		    if (res->x.p->size < 3 && e1->x.p->size == 3)
			explicit_complement(res);
		    for (i = 1; i < e1->x.p->size; ++i)
			eadd(&e1->x.p->arr[i], &res->x.p->arr[i]);
		    return;
	    }
	    if (res->x.p->size < 3)
		explicit_complement(res);
	    for (i = 1; i < res->x.p->size; ++i)
		eadd(e1, &res->x.p->arr[i]);
	    return;
	}
                 if ((e1->x.p->type != res->x.p->type) ) {
		      /* adding to evalues of different type. two cases are possible  
		       * res is periodic and e1 is polynomial, you have to exchange
		       * e1 and res then to add e1 to the constant term of res */
		     if (e1->x.p->type == polynomial) {
			  eadd_rev_cst(e1, res); 
	             }
                     else if (res->x.p->type == polynomial) {
                          /* res is polynomial and e1 is periodic,
		            add e1 to the constant term of res */
			 
			  eadd(e1,&res->x.p->arr[0]);
		     } else
			assert(0);
	                	 
		     return;
	         }
	         else if (e1->x.p->pos != res->x.p->pos ||
			    ((res->x.p->type == fractional ||
			      res->x.p->type == flooring) &&
			     !eequal(&e1->x.p->arr[0], &res->x.p->arr[0]))) { 
	      	 /* adding evalues of different position (i.e function of different unknowns
		  * to case are possible  */
			   
			switch (res->x.p->type) {
			case flooring:
			case fractional:
			    if (mod_term_smaller(res, e1))
				eadd(e1,&res->x.p->arr[1]);
			    else
				eadd_rev_cst(e1, res);
			    return;
			case polynomial: //  res and e1 are polynomials
			       //  add e1 to the constant term of res
			       
			    if(res->x.p->pos < e1->x.p->pos)
		               eadd(e1,&res->x.p->arr[0]);
			    else
				eadd_rev_cst(e1, res);
		              // value_clear(g); value_clear(m1); value_clear(m2);
		               return;
		        case periodic:  // res and e1 are pointers to periodic numbers
				  //add e1 to all elements of res 
				   
			    if(res->x.p->pos < e1->x.p->pos)
			          for (i=0;i<res->x.p->size;i++) {
			               eadd(e1,&res->x.p->arr[i]);
			          }
			    else
				eadd_rev(e1, res);
			    return;
			default:
			    assert(0);
			}
	         }  
                 
                
		 //same type , same pos  and same size
                 if (e1->x.p->size == res->x.p->size) {
	              // add any element in e1 to the corresponding element in res 
		      i = type_offset(res->x.p);
		      if (i == 1)
			assert(eequal(&e1->x.p->arr[0], &res->x.p->arr[0]));
	              for (; i<res->x.p->size; i++) {
                            eadd(&e1->x.p->arr[i], &res->x.p->arr[i]);
                      }
                      return ;
                }
                
		/* Sizes are different */
		switch(res->x.p->type) {
		case polynomial:
		case flooring:
		case fractional:
                    /* VIN100: if e1-size > res-size you have to copy e1 in a   */
                    /* new enode and add res to that new node. If you do not do */
                    /* that, you lose the the upper weight part of e1 !         */

                     if(e1->x.p->size > res->x.p->size)
			  eadd_rev(e1, res);
                     else {
		        i = type_offset(res->x.p);
		        if (i == 1)
			    assert(eequal(&e1->x.p->arr[0], 
				   &res->x.p->arr[0]));
                        for (; i<e1->x.p->size ; i++) {
                             eadd(&e1->x.p->arr[i], &res->x.p->arr[i]);
                        } 
                      		   
                        return ;
                     } 
		     break;
                
    /* add two periodics of the same pos (unknown) but whith different sizes (periods) */
	        case periodic:
		{
		      /* you have to create a new evalue 'ne' in whitch size equals to the lcm
		       of the sizes of e1 and res, then to copy res periodicaly in ne, after
		       to add periodicaly elements of e1 to elements of ne, and finaly to 
		       return ne. */
		       int x,y,p;
		       Value ex, ey ,ep;
		       evalue *ne;
	               value_init(ex); value_init(ey);value_init(ep);
		       x=e1->x.p->size;
	               y= res->x.p->size;
		       value_set_si(ex,e1->x.p->size);
		       value_set_si(ey,res->x.p->size);
		       value_assign (ep,*Lcm(ex,ey));
		       p=(int)mpz_get_si(ep);
	               ne= (evalue *) malloc (sizeof(evalue)); 
	               value_init(ne->d);
	               value_set_si( ne->d,0);
	    	           	    
	               ne->x.p=new_enode(res->x.p->type,p, res->x.p->pos);
	               for(i=0;i<p;i++)  {
			  evalue_copy(&ne->x.p->arr[i], &res->x.p->arr[i%y]);
	               }
	               for(i=0;i<p;i++)  {
	                    eadd(&e1->x.p->arr[i%x], &ne->x.p->arr[i]);
	               }
      
	              value_assign(res->d, ne->d);
		      res->x.p=ne->x.p;
	    	    
                        return ;
		}
		case evector:
                     fprintf(stderr, "eadd: ?cannot add vectors of different length\n");
                     return ;
		default:
		    assert(0);
                }
     }
     return ;
 }/* eadd  */ 
 
static void emul_rev (evalue *e1, evalue *res)
{
    evalue ev;
    value_init(ev.d);
    evalue_copy(&ev, e1);
    emul(res, &ev);
    free_evalue_refs(res);	  
    *res = ev;
}

static void emul_poly (evalue *e1, evalue *res)
{
    int i, j, o = type_offset(res->x.p);
    evalue tmp;
    int size=(e1->x.p->size + res->x.p->size - o - 1); 
    value_init(tmp.d);
    value_set_si(tmp.d,0);
    tmp.x.p=new_enode(res->x.p->type, size, res->x.p->pos);
    if (o)
	evalue_copy(&tmp.x.p->arr[0], &e1->x.p->arr[0]);
    for (i=o; i < e1->x.p->size; i++) {
	evalue_copy(&tmp.x.p->arr[i], &e1->x.p->arr[i]);
	emul(&res->x.p->arr[o], &tmp.x.p->arr[i]);
    }
    for (; i<size; i++)
	evalue_set_si(&tmp.x.p->arr[i], 0, 1);
    for (i=o+1; i<res->x.p->size; i++)
	for (j=o; j<e1->x.p->size; j++) {
	    evalue ev;
	    value_init(ev.d);
	    evalue_copy(&ev, &e1->x.p->arr[j]);
	    emul(&res->x.p->arr[i], &ev);
	    eadd(&ev, &tmp.x.p->arr[i+j-o]);
	    free_evalue_refs(&ev);
	}
    free_evalue_refs(res);
    *res = tmp;
}

void emul_partitions (evalue *e1,evalue *res)
{
    int n, i, j, k;
    Polyhedron *d;
    struct section *s;
    s = (struct section *) 
	    malloc((e1->x.p->size/2) * (res->x.p->size/2) * 
		   sizeof(struct section));
    assert(s);
    assert(e1->x.p->pos == res->x.p->pos);
    assert(e1->x.p->pos == EVALUE_DOMAIN(e1->x.p->arr[0])->Dimension);
    assert(res->x.p->pos == EVALUE_DOMAIN(res->x.p->arr[0])->Dimension);

    n = 0;
    for (i = 0; i < res->x.p->size/2; ++i) {
	for (j = 0; j < e1->x.p->size/2; ++j) {
	    d = DomainIntersection(EVALUE_DOMAIN(e1->x.p->arr[2*j]),
				   EVALUE_DOMAIN(res->x.p->arr[2*i]), 0);
	    if (emptyQ(d)) {
		Domain_Free(d);
		continue;
	    }

	    /* This code is only needed because the partitions
	       are not true partitions.
	     */
	    for (k = 0; k < n; ++k) {
		if (DomainIncludes(s[k].D, d))
		    break;
		if (DomainIncludes(d, s[k].D)) {
		    Domain_Free(s[k].D);
		    free_evalue_refs(&s[k].E);
		    if (n > k)
			s[k] = s[--n];
		    --k;
		}
	    }
	    if (k < n) {
		Domain_Free(d);
		continue;
	    }

	    value_init(s[n].E.d);
	    evalue_copy(&s[n].E, &res->x.p->arr[2*i+1]);
	    emul(&e1->x.p->arr[2*j+1], &s[n].E);
	    s[n].D = d;
	    ++n;
	}
	Domain_Free(EVALUE_DOMAIN(res->x.p->arr[2*i]));
	value_clear(res->x.p->arr[2*i].d);
	free_evalue_refs(&res->x.p->arr[2*i+1]);
    }

    free(res->x.p);
    if (n == 0)
	evalue_set_si(res, 0, 1);
    else {
	res->x.p = new_enode(partition, 2*n, e1->x.p->pos);
	for (j = 0; j < n; ++j) {
	    s[j].D = DomainConstraintSimplify(s[j].D, 0);
	    EVALUE_SET_DOMAIN(res->x.p->arr[2*j], s[j].D);
	    value_clear(res->x.p->arr[2*j+1].d);
	    res->x.p->arr[2*j+1] = s[j].E;
	}
    }

    free(s);
}

#define value_two_p(val)	(mpz_cmp_si(val,2) == 0)

/* Computes the product of two evalues "e1" and "res" and puts the result in "res". you must
 * do a copy of "res" befor calling this function if you nead it after. The vector type of 
 * evalues is not treated here */

void emul (evalue *e1, evalue *res ){
    int i,j;

if((value_zero_p(e1->d)&&e1->x.p->type==evector)||(value_zero_p(res->d)&&(res->x.p->type==evector))) {    
    fprintf(stderr, "emul: do not proced on evector type !\n");
    return;
}
     
    if (EVALUE_IS_ZERO(*res))
	return;

    if (value_zero_p(e1->d) && e1->x.p->type == partition) {
        if (value_zero_p(res->d) && res->x.p->type == partition)
	    emul_partitions(e1, res);
	else
	    emul_rev(e1, res);
    } else if (value_zero_p(res->d) && res->x.p->type == partition) {
	for (i = 0; i < res->x.p->size/2; ++i)
	    emul(e1, &res->x.p->arr[2*i+1]);
    } else
   if (value_zero_p(res->d) && res->x.p->type == relation) {
	if (value_zero_p(e1->d) && e1->x.p->type == relation &&
	    eequal(&e1->x.p->arr[0], &res->x.p->arr[0])) {
		if (res->x.p->size < 3 && e1->x.p->size == 3)
		    explicit_complement(res);
		if (e1->x.p->size < 3 && res->x.p->size == 3)
		    explicit_complement(e1);
		for (i = 1; i < res->x.p->size; ++i)
		    emul(&e1->x.p->arr[i], &res->x.p->arr[i]);
		return;
	}
	for (i = 1; i < res->x.p->size; ++i)
	    emul(e1, &res->x.p->arr[i]);
   } else
   if(value_zero_p(e1->d)&& value_zero_p(res->d)) {
       switch(e1->x.p->type) {
       case polynomial:
	   switch(res->x.p->type) {
	   case polynomial:
	       if(e1->x.p->pos == res->x.p->pos) {
	       /* Product of two polynomials of the same variable */
		    emul_poly(e1, res);
		    return;
	       }
	       else {
		  /* Product of two polynomials of different variables */     
	          
		if(res->x.p->pos < e1->x.p->pos)
		     for( i=0; i<res->x.p->size ; i++) 
		            emul(e1, &res->x.p->arr[i]);
		else
		    emul_rev(e1, res);
			  
		 return ;
	       }      
	   case periodic:
	   case flooring:
	   case fractional:
	        /* Product of a polynomial and a periodic or fractional */
		emul_rev(e1, res);
		return;
	   default:
		assert(0);
	   }
       case periodic:
	   switch(res->x.p->type) {
	   case periodic:
		if(e1->x.p->pos==res->x.p->pos && e1->x.p->size==res->x.p->size) {
		 /* Product of two periodics of the same parameter and period */	
                   
		     for(i=0; i<res->x.p->size;i++) 
		           emul(&(e1->x.p->arr[i]), &(res->x.p->arr[i]));
	             		     
		     return;
		}
		else{
                  if(e1->x.p->pos==res->x.p->pos && e1->x.p->size!=res->x.p->size) {  
	           /* Product of two periodics of the same parameter and different periods */		  
		    evalue *newp;
	            Value x,y,z;
		    int ix,iy,lcm;
		    value_init(x); value_init(y);value_init(z);
		    ix=e1->x.p->size;
		    iy=res->x.p->size;
		    value_set_si(x,e1->x.p->size);
		    value_set_si(y,res->x.p->size);
		    value_assign (z,*Lcm(x,y));
		    lcm=(int)mpz_get_si(z);
		    newp= (evalue *) malloc (sizeof(evalue));
		    value_init(newp->d);
		    value_set_si( newp->d,0);
		    newp->x.p=new_enode(periodic,lcm, e1->x.p->pos);
		        for(i=0;i<lcm;i++)  {
		           evalue_copy(&newp->x.p->arr[i], 
				       &res->x.p->arr[i%iy]);
			}
			for(i=0;i<lcm;i++)  
		            emul(&e1->x.p->arr[i%ix], &newp->x.p->arr[i]);		
			  
			value_assign(res->d,newp->d);
			res->x.p=newp->x.p;
			
			  value_clear(x); value_clear(y);value_clear(z); 
			  return ;	 
		  }
		  else {
	             /* Product of two periodics of different parameters */
			  
			if(res->x.p->pos < e1->x.p->pos)
			    for(i=0; i<res->x.p->size; i++)
				emul(e1, &(res->x.p->arr[i]));
			else
			    emul_rev(e1, res);
			
			return;
		  }
		}		       
	   case polynomial:
                  /* Product of a periodic and a polynomial */
			  
		       for(i=0; i<res->x.p->size ; i++)
		            emul(e1, &(res->x.p->arr[i]));    
                 
		       return; 
			       
	   }		   
       case flooring:
       case fractional:
	    switch(res->x.p->type) {
	    case polynomial:
	        for(i=0; i<res->x.p->size ; i++)
		    emul(e1, &(res->x.p->arr[i]));    
	        return; 
	    default:
	    case periodic:
		assert(0);
	    case flooring:
	    case fractional:
		assert(e1->x.p->type == res->x.p->type);
	        if (e1->x.p->pos == res->x.p->pos &&
			    eequal(&e1->x.p->arr[0], &res->x.p->arr[0])) {
		    evalue d;
		    value_init(d.d);
		    poly_denom(&e1->x.p->arr[0], &d.d);
		    if (e1->x.p->type != fractional || !value_two_p(d.d))
			emul_poly(e1, res);
		    else {
			evalue tmp;
			value_init(d.x.n);
			value_set_si(d.x.n, 1);
			/* { x }^2 == { x }/2 */
			/* a0 b0 + (a0 b1 + a1 b0 + a1 b1/2) { x } */
			assert(e1->x.p->size == 3);
			assert(res->x.p->size == 3);
			value_init(tmp.d);
			evalue_copy(&tmp, &res->x.p->arr[2]);
			emul(&d, &tmp);
			eadd(&res->x.p->arr[1], &tmp);
			emul(&e1->x.p->arr[2], &tmp);
			emul(&e1->x.p->arr[1], res);
			eadd(&tmp, &res->x.p->arr[2]);
			free_evalue_refs(&tmp);	  
			value_clear(d.x.n);
		    }
		    value_clear(d.d);
		} else {
		    if(mod_term_smaller(res, e1))
			for(i=1; i<res->x.p->size ; i++)
			    emul(e1, &(res->x.p->arr[i]));    
		    else
			emul_rev(e1, res);
		    return; 
		}
	    }
	    break;
       case relation:
	    emul_rev(e1, res);
	    break;
       default:
	    assert(0);
       }		   
   }
   else {
       if (value_notzero_p(e1->d)&& value_notzero_p(res->d)) {
	   /* Product of two rational numbers */
	
	    Value g;
	    value_init(g);
	    value_multiply(res->d,e1->d,res->d);
	    value_multiply(res->x.n,e1->x.n,res->x.n );
	    Gcd(res->x.n, res->d,&g);
	   if (value_notone_p(g)) {
	       value_division(res->d,res->d,g);
	       value_division(res->x.n,res->x.n,g);
	   }
	   value_clear(g);
	   return ;
       }
       else { 
	     if(value_zero_p(e1->d)&& value_notzero_p(res->d)) { 
	      /* Product of an expression (polynomial or peririodic) and a rational number */
		     
		emul_rev(e1, res);
		 return ;
	     }
	     else {
	       /* Product of a rationel number and an expression (polynomial or peririodic) */ 
	         
		   i = type_offset(res->x.p);
		   for (; i<res->x.p->size; i++) 
	              emul(e1, &res->x.p->arr[i]);		 
		  
		 return ;
	     }
       }
   }
   
   return ;
}

/* Frees mask content ! */
void emask(evalue *mask, evalue *res) {
    int n, i, j;
    Polyhedron *d, *fd;
    struct section *s;
    evalue mone;
    int pos;

    if (EVALUE_IS_ZERO(*res)) {
	free_evalue_refs(mask); 
	return;
    }

    assert(value_zero_p(mask->d));
    assert(mask->x.p->type == partition);
    assert(value_zero_p(res->d));
    assert(res->x.p->type == partition);
    assert(mask->x.p->pos == res->x.p->pos);
    assert(res->x.p->pos == EVALUE_DOMAIN(res->x.p->arr[0])->Dimension);
    assert(mask->x.p->pos == EVALUE_DOMAIN(mask->x.p->arr[0])->Dimension);
    pos = res->x.p->pos;

    s = (struct section *) 
	    malloc((mask->x.p->size/2+1) * (res->x.p->size/2) * 
		   sizeof(struct section));
    assert(s);

    value_init(mone.d);
    evalue_set_si(&mone, -1, 1);

    n = 0;
    for (j = 0; j < res->x.p->size/2; ++j) {
	assert(mask->x.p->size >= 2);
	fd = DomainDifference(EVALUE_DOMAIN(res->x.p->arr[2*j]),
			      EVALUE_DOMAIN(mask->x.p->arr[0]), 0);
	if (!emptyQ(fd))
	    for (i = 1; i < mask->x.p->size/2; ++i) {
		Polyhedron *t = fd;
		fd = DomainDifference(fd, EVALUE_DOMAIN(mask->x.p->arr[2*i]), 0);
		Domain_Free(t);
		if (emptyQ(fd))
		    break;
	    }
	if (emptyQ(fd)) {
	    Domain_Free(fd);
	    continue;
	}
	value_init(s[n].E.d);
	evalue_copy(&s[n].E, &res->x.p->arr[2*j+1]);
	s[n].D = fd;
	++n;
    }
    for (i = 0; i < mask->x.p->size/2; ++i) {
	if (EVALUE_IS_ONE(mask->x.p->arr[2*i+1]))
	    continue;

	fd = EVALUE_DOMAIN(mask->x.p->arr[2*i]);
	eadd(&mone, &mask->x.p->arr[2*i+1]);
	emul(&mone, &mask->x.p->arr[2*i+1]);
	for (j = 0; j < res->x.p->size/2; ++j) {
	    Polyhedron *t;
	    d = DomainIntersection(EVALUE_DOMAIN(res->x.p->arr[2*j]),
				   EVALUE_DOMAIN(mask->x.p->arr[2*i]), 0);
	    if (emptyQ(d)) {
		Domain_Free(d);
		continue;
	    }
	    t = fd;
	    fd = DomainDifference(fd, EVALUE_DOMAIN(res->x.p->arr[2*j]), 0);
	    if (t != EVALUE_DOMAIN(mask->x.p->arr[2*i]))
		Domain_Free(t);
	    value_init(s[n].E.d);
	    evalue_copy(&s[n].E, &res->x.p->arr[2*j+1]);
	    emul(&mask->x.p->arr[2*i+1], &s[n].E);
	    s[n].D = d;
	    ++n;
	}

	if (!emptyQ(fd)) {
	    /* Just ignore; this may have been previously masked off */
	}
	if (fd != EVALUE_DOMAIN(mask->x.p->arr[2*i]))
	    Domain_Free(fd);
    }

    free_evalue_refs(&mone);
    free_evalue_refs(mask);
    free_evalue_refs(res);
    value_init(res->d);
    if (n == 0)
	evalue_set_si(res, 0, 1);
    else {
	res->x.p = new_enode(partition, 2*n, pos);
	for (j = 0; j < n; ++j) {
	    EVALUE_SET_DOMAIN(res->x.p->arr[2*j], s[j].D);
	    value_clear(res->x.p->arr[2*j+1].d);
	    res->x.p->arr[2*j+1] = s[j].E;
	}
    }

    free(s);
}

void evalue_copy(evalue *dst, const evalue *src)
{
    value_assign(dst->d, src->d);
    if(value_notzero_p(src->d)) {
	 value_init(dst->x.n);
	 value_assign(dst->x.n, src->x.n);
    } else
	 dst->x.p = ecopy(src->x.p);
}

enode *new_enode(enode_type type,int size,int pos) {
  
  enode *res;
  int i;
  
  if(size == 0) {
    fprintf(stderr, "Allocating enode of size 0 !\n" );
    return NULL;
  }
  res = (enode *) malloc(sizeof(enode) + (size-1)*sizeof(evalue));
  res->type = type;
  res->size = size;
  res->pos = pos;
  for(i=0; i<size; i++) {
    value_init(res->arr[i].d);
    value_set_si(res->arr[i].d,0);
    res->arr[i].x.p = 0;
  }
  return res;
} /* new_enode */

enode *ecopy(enode *e) {
  
  enode *res;
  int i;
  
  res = new_enode(e->type,e->size,e->pos);
  for(i=0;i<e->size;++i) {
    value_assign(res->arr[i].d,e->arr[i].d);
    if(value_zero_p(res->arr[i].d))
      res->arr[i].x.p = ecopy(e->arr[i].x.p);
    else if (EVALUE_IS_DOMAIN(res->arr[i]))
      EVALUE_SET_DOMAIN(res->arr[i], Domain_Copy(EVALUE_DOMAIN(e->arr[i])));
    else {
      value_init(res->arr[i].x.n);
      value_assign(res->arr[i].x.n,e->arr[i].x.n);
    }
  }
  return(res);
} /* ecopy */

int ecmp(const evalue *e1, const evalue *e2)
{
    enode *p1, *p2;
    int i;
    int r;

    if (value_notzero_p(e1->d) && value_notzero_p(e2->d)) {
	Value m, m2;
	value_init(m);
	value_init(m2);
	value_multiply(m, e1->x.n, e2->d);
	value_multiply(m2, e2->x.n, e1->d);

	if (value_lt(m, m2))
	    r = -1;
	else if (value_gt(m, m2))
	    r = 1;
	else 
	    r = 0;

	value_clear(m);
	value_clear(m2);

	return r;
    }
    if (value_notzero_p(e1->d))
	return -1;
    if (value_notzero_p(e2->d))
	return 1;

    p1 = e1->x.p;
    p2 = e2->x.p;

    if (p1->type != p2->type)
	return p1->type - p2->type;
    if (p1->pos != p2->pos)
	return p1->pos - p2->pos;
    if (p1->size != p2->size)
	return p1->size - p2->size;

    for (i = p1->size-1; i >= 0; --i)
	if ((r = ecmp(&p1->arr[i], &p2->arr[i])) != 0)
	    return r;

    return 0;
}

int eequal(const evalue *e1, const evalue *e2)
{ 
    int i;
    enode *p1, *p2;
  
    if (value_ne(e1->d,e2->d))
        return 0;
  
    /* e1->d == e2->d */
    if (value_notzero_p(e1->d)) {    
        if (value_ne(e1->x.n,e2->x.n))
            return 0;
    
        /* e1->d == e2->d != 0  AND e1->n == e2->n */
        return 1;
    }
  
    /* e1->d == e2->d == 0 */
    p1 = e1->x.p;
    p2 = e2->x.p;
    if (p1->type != p2->type) return 0;
    if (p1->size != p2->size) return 0;
    if (p1->pos  != p2->pos) return 0;
    for (i=0; i<p1->size; i++)
        if (!eequal(&p1->arr[i], &p2->arr[i]) ) 
            return 0;
    return 1;
} /* eequal */

void free_evalue_refs(evalue *e) {
  
  enode *p;
  int i;
  
  if (EVALUE_IS_DOMAIN(*e)) {
    Domain_Free(EVALUE_DOMAIN(*e));
    value_clear(e->d);
    return;
  } else if (value_pos_p(e->d)) {
    
    /* 'e' stores a constant */
    value_clear(e->d);
    value_clear(e->x.n);
    return; 
  }  
  assert(value_zero_p(e->d));
  value_clear(e->d);
  p = e->x.p;
  if (!p) return;	/* null pointer */
  for (i=0; i<p->size; i++) {
    free_evalue_refs(&(p->arr[i]));
  }
  free(p);
  return;
} /* free_evalue_refs */

static void mod2table_r(evalue *e, Vector *periods, Value m, int p, 
			Vector * val, evalue *res)
{
    unsigned nparam = periods->Size;

    if (p == nparam) {
	double d = compute_evalue(e, val->p);
	d *= VALUE_TO_DOUBLE(m);
	if (d > 0)
	    d += .25;
	else
	    d -= .25;
	value_assign(res->d, m);
	value_init(res->x.n);
	value_set_double(res->x.n, d);
	mpz_fdiv_r(res->x.n, res->x.n, m);
	return;
    }
    if (value_one_p(periods->p[p]))
	mod2table_r(e, periods, m, p+1, val, res);
    else {
	Value tmp;
	value_init(tmp);

	value_assign(tmp, periods->p[p]);
	value_set_si(res->d, 0);
	res->x.p = new_enode(periodic, VALUE_TO_INT(tmp), p+1);
	do {
	    value_decrement(tmp, tmp);
	    value_assign(val->p[p], tmp);
	    mod2table_r(e, periods, m, p+1, val, 
			&res->x.p->arr[VALUE_TO_INT(tmp)]);
	} while (value_pos_p(tmp));

	value_clear(tmp);
    }
}

static void rel2table(evalue *e, int zero)
{
    if (value_pos_p(e->d)) {
	if (value_zero_p(e->x.n) == zero)
	    value_set_si(e->x.n, 1);
	else
	    value_set_si(e->x.n, 0);
	value_set_si(e->d, 1);
    } else {
	int i;
	for (i = 0; i < e->x.p->size; ++i)
	    rel2table(&e->x.p->arr[i], zero);
    }
}

void evalue_mod2table(evalue *e, int nparam)
{
  enode *p;
  int i;

  if (EVALUE_IS_DOMAIN(*e) || value_pos_p(e->d))
    return;
  p = e->x.p;
  for (i=0; i<p->size; i++) {
    evalue_mod2table(&(p->arr[i]), nparam);
  }
  if (p->type == relation) {
    evalue copy;

    if (p->size > 2) {
      value_init(copy.d);
      evalue_copy(&copy, &p->arr[0]);
    }
    rel2table(&p->arr[0], 1);
    emul(&p->arr[0], &p->arr[1]);
    if (p->size > 2) {
      rel2table(&copy, 0);
      emul(&copy, &p->arr[2]);
      eadd(&p->arr[2], &p->arr[1]);
      free_evalue_refs(&p->arr[2]);	  
      free_evalue_refs(&copy);	  
    }
    free_evalue_refs(&p->arr[0]);	  
    value_clear(e->d);
    *e = p->arr[1];
    free(p);
  } else if (p->type == fractional) {
    Vector *periods = Vector_Alloc(nparam);
    Vector *val = Vector_Alloc(nparam);
    Value tmp;
    evalue *ev;
    evalue EP, res;

    value_init(tmp);
    value_set_si(tmp, 1);
    Vector_Set(periods->p, 1, nparam);
    Vector_Set(val->p, 0, nparam);
    for (ev = &p->arr[0]; value_zero_p(ev->d); ev = &ev->x.p->arr[0]) {
      enode *p = ev->x.p;

      assert(p->type == polynomial);
      assert(p->size == 2);
      value_assign(periods->p[p->pos-1], p->arr[1].d);
      value_lcm(tmp, p->arr[1].d, &tmp);
    }
    value_lcm(tmp, ev->d, &tmp);
    value_init(EP.d);
    mod2table_r(&p->arr[0], periods, tmp, 0, val, &EP);

    value_init(res.d);
    evalue_set_si(&res, 0, 1);
    /* Compute the polynomial using Horner's rule */
    for (i=p->size-1;i>1;i--) {
      eadd(&p->arr[i], &res);
      emul(&EP, &res);
    }
    eadd(&p->arr[1], &res);

    free_evalue_refs(e);	  
    free_evalue_refs(&EP);	  
    *e = res;

    value_clear(tmp);
    Vector_Free(val);
    Vector_Free(periods);
  }
} /* evalue_mod2table */

/********************************************************/
/* function in domain                                   */
/*    check if the parameters in list_args              */
/*    verifies the constraints of Domain P          	*/
/********************************************************/
int in_domain(Polyhedron *P, Value *list_args)
{
  int row, in = 1;
  Value v; /* value of the constraint of a row when
	       parameters are instantiated*/

  value_init(v); 
  
  for (row = 0; row < P->NbConstraints; row++) {
    Inner_Product(P->Constraint[row]+1, list_args, P->Dimension, &v);
    value_addto(v, v, P->Constraint[row][P->Dimension+1]); /*constant part*/
    if (value_neg_p(v) ||
	value_zero_p(P->Constraint[row][0]) && value_notzero_p(v)) {
      in = 0;
      break;
    }
  }
  
  value_clear(v);
  return in || (P->next && in_domain(P->next, list_args));
} /* in_domain */

/****************************************************/
/* function compute enode                           */
/*     compute the value of enode p with parameters */
/*     list "list_args                              */
/*     compute the polynomial or the periodic       */
/****************************************************/

static double compute_enode(enode *p, Value *list_args) {
  
  int i;
  Value m, param;
  double res=0.0;
    
  if (!p)
    return(0.);

  value_init(m);
  value_init(param);

  if (p->type == polynomial) {
    if (p->size > 1)
	value_assign(param,list_args[p->pos-1]);
    
    /* Compute the polynomial using Horner's rule */
    for (i=p->size-1;i>0;i--) {
      res +=compute_evalue(&p->arr[i],list_args);
      res *=VALUE_TO_DOUBLE(param);
    }
    res +=compute_evalue(&p->arr[0],list_args);
  }
  else if (p->type == fractional) {
    double d = compute_evalue(&p->arr[0], list_args);
    d -= floor(d+1e-10);
    
    /* Compute the polynomial using Horner's rule */
    for (i=p->size-1;i>1;i--) {
      res +=compute_evalue(&p->arr[i],list_args);
      res *=d;
    }
    res +=compute_evalue(&p->arr[1],list_args);
  }
  else if (p->type == flooring) {
    double d = compute_evalue(&p->arr[0], list_args);
    d = floor(d+1e-10);
    
    /* Compute the polynomial using Horner's rule */
    for (i=p->size-1;i>1;i--) {
      res +=compute_evalue(&p->arr[i],list_args);
      res *=d;
    }
    res +=compute_evalue(&p->arr[1],list_args);
  }
  else if (p->type == periodic) {
    value_assign(m,list_args[p->pos-1]);
    
    /* Choose the right element of the periodic */
    value_set_si(param,p->size);
    value_pmodulus(m,m,param);
    res = compute_evalue(&p->arr[VALUE_TO_INT(m)],list_args);
  }
  else if (p->type == relation) {
    if (fabs(compute_evalue(&p->arr[0], list_args)) < 1e-10)
      res = compute_evalue(&p->arr[1], list_args);
    else if (p->size > 2)
      res = compute_evalue(&p->arr[2], list_args);
  }
  else if (p->type == partition) {
    int dim = EVALUE_DOMAIN(p->arr[0])->Dimension;
    Value *vals = list_args;
    if (p->pos < dim) {
	NALLOC(vals, dim);
	for (i = 0; i < dim; ++i) {
	    value_init(vals[i]);
	    if (i < p->pos)
		value_assign(vals[i], list_args[i]);
	}
    }
    for (i = 0; i < p->size/2; ++i)
      if (DomainContains(EVALUE_DOMAIN(p->arr[2*i]), vals, p->pos, 0, 1)) {
	res = compute_evalue(&p->arr[2*i+1], vals);
	break;
      }
    if (p->pos < dim) {
	for (i = 0; i < dim; ++i)
	    value_clear(vals[i]);
	free(vals);
    }
  }
  else
    assert(0);
  value_clear(m);
  value_clear(param);
  return res;
} /* compute_enode */

/*************************************************/
/* return the value of Ehrhart Polynomial        */
/* It returns a double, because since it is      */
/* a recursive function, some intermediate value */
/* might not be integral                         */
/*************************************************/

double compute_evalue(const evalue *e, Value *list_args)
{
  double res;
  
  if (value_notzero_p(e->d)) {
    if (value_notone_p(e->d)) 
      res = VALUE_TO_DOUBLE(e->x.n) / VALUE_TO_DOUBLE(e->d);
    else 
      res = VALUE_TO_DOUBLE(e->x.n);
  }
  else 
    res = compute_enode(e->x.p,list_args);
  return res;
} /* compute_evalue */


/****************************************************/
/* function compute_poly :                          */
/* Check for the good validity domain               */
/* return the number of point in the Polyhedron     */
/* in allocated memory                              */
/* Using the Ehrhart pseudo-polynomial              */
/****************************************************/
Value *compute_poly(Enumeration *en,Value *list_args) {

  Value *tmp;
  /*	double d; int i; */

  tmp = (Value *) malloc (sizeof(Value));
  assert(tmp != NULL);
  value_init(*tmp);
  value_set_si(*tmp,0);
  
  if(!en)
    return(tmp);	/* no ehrhart polynomial */
  if(en->ValidityDomain) {
    if(!en->ValidityDomain->Dimension) { /* no parameters */
      value_set_double(*tmp,compute_evalue(&en->EP,list_args)+.25);
      return(tmp);
    }
  }  
  else 
    return(tmp);  /* no Validity Domain */    
  while(en) {
    if(in_domain(en->ValidityDomain,list_args)) {
      
#ifdef EVAL_EHRHART_DEBUG
      Print_Domain(stdout,en->ValidityDomain);
      print_evalue(stdout,&en->EP);
#endif
      
      /*			d = compute_evalue(&en->EP,list_args);
				i = d;
				printf("(double)%lf = %d\n", d, i ); */
      value_set_double(*tmp,compute_evalue(&en->EP,list_args)+.25);
      return(tmp);
    }
    else
      en=en->next;
  }
  value_set_si(*tmp,0);
  return(tmp); /* no compatible domain with the arguments */
} /* compute_poly */ 

static evalue *eval_polynomial(const enode *p, int offset,
			       evalue *base, Value *values)
{
    int i;
    evalue *res, *c;

    res = evalue_zero();
    for (i = p->size-1; i > offset; --i) {
	c = evalue_eval(&p->arr[i], values);
	eadd(c, res);
	free_evalue_refs(c);
	free(c);
	emul(base, res);
    }
    c = evalue_eval(&p->arr[offset], values);
    eadd(c, res);
    free_evalue_refs(c);
    free(c);

    return res;
}

evalue *evalue_eval(const evalue *e, Value *values)
{
    evalue *res = NULL;
    evalue param;
    evalue *param2;
    int i;

    if (value_notzero_p(e->d)) {
	res = ALLOC(evalue);
	value_init(res->d);
	evalue_copy(res, e);
	return res;
    }
    switch (e->x.p->type) {
    case polynomial:
	value_init(param.x.n);
	value_assign(param.x.n, values[e->x.p->pos-1]);
	value_init(param.d);
	value_set_si(param.d, 1);

	res = eval_polynomial(e->x.p, 0, &param, values);
	free_evalue_refs(&param);
	break;
    case fractional:
	param2 = evalue_eval(&e->x.p->arr[0], values);
	mpz_fdiv_r(param2->x.n, param2->x.n, param2->d);

	res = eval_polynomial(e->x.p, 1, param2, values);
	free_evalue_refs(param2);
	free(param2);
	break;
    case flooring:
	param2 = evalue_eval(&e->x.p->arr[0], values);
	mpz_fdiv_q(param2->x.n, param2->x.n, param2->d);
	value_set_si(param2->d, 1);

	res = eval_polynomial(e->x.p, 1, param2, values);
	free_evalue_refs(param2);
	free(param2);
	break;
    case relation:
	param2 = evalue_eval(&e->x.p->arr[0], values);
	if (value_zero_p(param2->x.n))
	    res = evalue_eval(&e->x.p->arr[1], values);
	else if (e->x.p->size > 2)
	    res = evalue_eval(&e->x.p->arr[2], values);
	else
	    res = evalue_zero();
	free_evalue_refs(param2);
	free(param2);
	break;
    case partition:
    	assert(e->x.p->pos == EVALUE_DOMAIN(e->x.p->arr[0])->Dimension);
	for (i = 0; i < e->x.p->size/2; ++i)
	    if (in_domain(EVALUE_DOMAIN(e->x.p->arr[2*i]), values)) {
		res = evalue_eval(&e->x.p->arr[2*i+1], values);
		break;
	    }
	if (!res)
	    res = evalue_zero();
	break;
    default:
	assert(0);
    }
    return res;
}

size_t value_size(Value v) {
    return (v[0]._mp_size > 0 ? v[0]._mp_size : -v[0]._mp_size)
	    * sizeof(v[0]._mp_d[0]);
}

size_t domain_size(Polyhedron *D)
{
    int i, j;
    size_t s = sizeof(*D);

    for (i = 0; i < D->NbConstraints; ++i)
	for (j = 0; j < D->Dimension+2; ++j)
	    s += value_size(D->Constraint[i][j]);

/*
    for (i = 0; i < D->NbRays; ++i)
	for (j = 0; j < D->Dimension+2; ++j)
	    s += value_size(D->Ray[i][j]);
*/

    return D->next ? s+domain_size(D->next) : s;
}

size_t enode_size(enode *p) {
    size_t s = sizeof(*p) - sizeof(p->arr[0]);
    int i;

    if (p->type == partition)
	for (i = 0; i < p->size/2; ++i) {
	    s += domain_size(EVALUE_DOMAIN(p->arr[2*i]));
	    s += evalue_size(&p->arr[2*i+1]);
	}
    else
	for (i = 0; i < p->size; ++i) {
	    s += evalue_size(&p->arr[i]);
	}
    return s;
}

size_t evalue_size(evalue *e)
{
    size_t s = sizeof(*e);
    s += value_size(e->d);
    if (value_notzero_p(e->d))
	s += value_size(e->x.n);
    else
	s += enode_size(e->x.p);
    return s;
}

static evalue *find_second(evalue *base, evalue *cst, evalue *e, Value m)
{
    evalue *found = NULL;
    evalue offset;
    evalue copy;
    int i;

    if (value_pos_p(e->d) || e->x.p->type != fractional)
	return NULL;

    value_init(offset.d);
    value_init(offset.x.n);
    poly_denom(&e->x.p->arr[0], &offset.d);
    value_lcm(m, offset.d, &offset.d);
    value_set_si(offset.x.n, 1);

    value_init(copy.d);
    evalue_copy(&copy, cst);

    eadd(&offset, cst);
    mpz_fdiv_r(cst->x.n, cst->x.n, cst->d);

    if (eequal(base, &e->x.p->arr[0]))
	found = &e->x.p->arr[0];
    else {
	value_set_si(offset.x.n, -2);

	eadd(&offset, cst);
	mpz_fdiv_r(cst->x.n, cst->x.n, cst->d);

	if (eequal(base, &e->x.p->arr[0]))
	    found = base;
    }
    free_evalue_refs(cst);
    free_evalue_refs(&offset);
    *cst = copy;

    for (i = 1; !found && i < e->x.p->size; ++i)
	found = find_second(base, cst, &e->x.p->arr[i], m);

    return found;
}

static evalue *find_relation_pair(evalue *e)
{
    int i;
    evalue *found = NULL;

    if (EVALUE_IS_DOMAIN(*e) || value_pos_p(e->d))
	return NULL;

    if (e->x.p->type == fractional) {
	Value m;
	evalue *cst;

	value_init(m);
	poly_denom(&e->x.p->arr[0], &m);

	for (cst = &e->x.p->arr[0]; value_zero_p(cst->d); 
				    cst = &cst->x.p->arr[0])
	    ;

	for (i = 1; !found && i < e->x.p->size; ++i)
	    found = find_second(&e->x.p->arr[0], cst, &e->x.p->arr[i], m);

	value_clear(m);
    }

    i = e->x.p->type == relation;
    for (; !found && i < e->x.p->size; ++i)
	found = find_relation_pair(&e->x.p->arr[i]);

    return found;
}

void evalue_mod2relation(evalue *e) {
    evalue *d;

    if (value_zero_p(e->d) && e->x.p->type == partition) {
	int i;

	for (i = 0; i < e->x.p->size/2; ++i) {
	    evalue_mod2relation(&e->x.p->arr[2*i+1]);
	    if (EVALUE_IS_ZERO(e->x.p->arr[2*i+1])) {
		value_clear(e->x.p->arr[2*i].d);
		free_evalue_refs(&e->x.p->arr[2*i+1]);
		e->x.p->size -= 2;
		if (2*i < e->x.p->size) {
		    e->x.p->arr[2*i] = e->x.p->arr[e->x.p->size];
		    e->x.p->arr[2*i+1] = e->x.p->arr[e->x.p->size+1];
		}
		--i;
	    }
	}
	if (e->x.p->size == 0) {
	    free(e->x.p);
	    evalue_set_si(e, 0, 1);
	}

	return;
    }

    while ((d = find_relation_pair(e)) != NULL) {
	evalue split;
	evalue *ev;

	value_init(split.d);
	value_set_si(split.d, 0);
	split.x.p = new_enode(relation, 3, 0);
	evalue_set_si(&split.x.p->arr[1], 1, 1);
	evalue_set_si(&split.x.p->arr[2], 1, 1);

	ev = &split.x.p->arr[0];
	value_set_si(ev->d, 0);
	ev->x.p = new_enode(fractional, 3, -1);
	evalue_set_si(&ev->x.p->arr[1], 0, 1);
	evalue_set_si(&ev->x.p->arr[2], 1, 1);
	evalue_copy(&ev->x.p->arr[0], d);

	emul(&split, e);

	reduce_evalue(e);

	free_evalue_refs(&split);	  
    }
}

static int evalue_comp(const void * a, const void * b)
{
    const evalue *e1 = *(const evalue **)a;
    const evalue *e2 = *(const evalue **)b;
    return ecmp(e1, e2);
}

void evalue_combine(evalue *e)
{
    evalue **evs;
    int i, k;
    enode *p;
    evalue tmp;

    if (value_notzero_p(e->d) || e->x.p->type != partition)
	return;

    NALLOC(evs, e->x.p->size/2);
    for (i = 0; i < e->x.p->size/2; ++i)
	evs[i] = &e->x.p->arr[2*i+1];
    qsort(evs, e->x.p->size/2, sizeof(evs[0]), evalue_comp);
    p = new_enode(partition, e->x.p->size, e->x.p->pos);
    for (i = 0, k = 0; i < e->x.p->size/2; ++i) {
	if (k == 0 || ecmp(&p->arr[2*k-1], evs[i]) != 0) {
	    value_clear(p->arr[2*k].d);
	    value_clear(p->arr[2*k+1].d);
	    p->arr[2*k] = *(evs[i]-1);
	    p->arr[2*k+1] = *(evs[i]);
	    ++k;
	} else {
	    Polyhedron *D = EVALUE_DOMAIN(*(evs[i]-1));
	    Polyhedron *L = D;

	    value_clear((evs[i]-1)->d);

	    while (L->next)
		L = L->next;
	    L->next = EVALUE_DOMAIN(p->arr[2*k-2]);
	    EVALUE_SET_DOMAIN(p->arr[2*k-2], D);
	    free_evalue_refs(evs[i]);
	}
    }

    for (i = 2*k ; i < p->size; ++i)
	value_clear(p->arr[i].d);

    free(evs);
    free(e->x.p);
    p->size = 2*k;
    e->x.p = p;

    for (i = 0; i < e->x.p->size/2; ++i) {
	Polyhedron *H;
	if (value_notzero_p(e->x.p->arr[2*i+1].d))
	    continue;
	H = DomainConvex(EVALUE_DOMAIN(e->x.p->arr[2*i]), 0);
	if (H == NULL)
	    continue;
	for (k = 0; k < e->x.p->size/2; ++k) {
	    Polyhedron *D, *N, **P;
	    if (i == k)
		continue;
	    P = &EVALUE_DOMAIN(e->x.p->arr[2*k]);
	    D = *P;
	    if (D == NULL)
		continue;
	    for (; D; D = N) {
		*P = D;
		N = D->next;
		if (D->NbEq <= H->NbEq) {
		    P = &D->next;
		    continue;
		}

		value_init(tmp.d);
		tmp.x.p = new_enode(partition, 2, e->x.p->pos);
		EVALUE_SET_DOMAIN(tmp.x.p->arr[0], Polyhedron_Copy(D));
		evalue_copy(&tmp.x.p->arr[1], &e->x.p->arr[2*i+1]);
		reduce_evalue(&tmp);
		if (value_notzero_p(tmp.d) ||
			ecmp(&tmp.x.p->arr[1], &e->x.p->arr[2*k+1]) != 0)
		    P = &D->next;
		else {
		    D->next = EVALUE_DOMAIN(e->x.p->arr[2*i]);
		    EVALUE_DOMAIN(e->x.p->arr[2*i]) = D;
		    *P = NULL;
		}
		free_evalue_refs(&tmp);
	    }
	}
	Polyhedron_Free(H);
    }

    for (i = 0; i < e->x.p->size/2; ++i) {
	Polyhedron *H, *E;
	Polyhedron *D = EVALUE_DOMAIN(e->x.p->arr[2*i]);
	if (!D) {
	    value_clear(e->x.p->arr[2*i].d);
	    free_evalue_refs(&e->x.p->arr[2*i+1]);
	    e->x.p->size -= 2;
	    if (2*i < e->x.p->size) {
		e->x.p->arr[2*i] = e->x.p->arr[e->x.p->size];
		e->x.p->arr[2*i+1] = e->x.p->arr[e->x.p->size+1];
	    }
	    --i;
	    continue;
	}
	if (!D->next)
	    continue;
	H = DomainConvex(D, 0);
	E = DomainDifference(H, D, 0);
	Domain_Free(D);
	D = DomainDifference(H, E, 0);
	Domain_Free(H);
	Domain_Free(E);
	EVALUE_SET_DOMAIN(p->arr[2*i], D);
    }
}

/* Use smallest representative for coefficients in affine form in
 * argument of fractional.
 * Since any change will make the argument non-standard,
 * the containing evalue will have to be reduced again afterward.
 */
static void fractional_minimal_coefficients(enode *p)
{
    evalue *pp;
    Value twice;
    value_init(twice);

    assert(p->type == fractional);
    pp = &p->arr[0];
    while (value_zero_p(pp->d)) {
	assert(pp->x.p->type == polynomial);
	assert(pp->x.p->size == 2);
	assert(value_notzero_p(pp->x.p->arr[1].d));
	mpz_mul_ui(twice, pp->x.p->arr[1].x.n, 2);
	if (value_gt(twice, pp->x.p->arr[1].d))
	    value_subtract(pp->x.p->arr[1].x.n, 
			   pp->x.p->arr[1].x.n, pp->x.p->arr[1].d);
	pp = &pp->x.p->arr[0];
    }

    value_clear(twice);
}

static Polyhedron *polynomial_projection(enode *p, Polyhedron *D, Value *d,
					 Matrix **R)
{
    Polyhedron *I, *H;
    evalue *pp;
    unsigned dim = D->Dimension;
    Matrix *T = Matrix_Alloc(2, dim+1);
    assert(T);

    assert(p->type == fractional);
    pp = &p->arr[0];
    value_set_si(T->p[1][dim], 1);
    poly_denom(pp, d);
    while (value_zero_p(pp->d)) {
	assert(pp->x.p->type == polynomial);
	assert(pp->x.p->size == 2);
	assert(value_notzero_p(pp->x.p->arr[1].d));
	value_division(T->p[0][pp->x.p->pos-1], *d, pp->x.p->arr[1].d);
	value_multiply(T->p[0][pp->x.p->pos-1], 
		       T->p[0][pp->x.p->pos-1], pp->x.p->arr[1].x.n);
	pp = &pp->x.p->arr[0];
    }
    value_division(T->p[0][dim], *d, pp->d);
    value_multiply(T->p[0][dim], T->p[0][dim], pp->x.n);
    I = DomainImage(D, T, 0);
    H = DomainConvex(I, 0);
    Domain_Free(I);
    if (R)
	*R = T;
    else
	Matrix_Free(T);

    return H;
}

int evalue_range_reduction_in_domain(evalue *e, Polyhedron *D)
{
    int i;
    enode *p;
    Value d, min, max;
    int r = 0;
    Polyhedron *I;
    int bounded;

    if (value_notzero_p(e->d))
	return r;

    p = e->x.p;

    if (p->type == relation) {
	Matrix *T;
	int equal;
	value_init(d);
	value_init(min);
	value_init(max);

	fractional_minimal_coefficients(p->arr[0].x.p);
	I = polynomial_projection(p->arr[0].x.p, D, &d, &T);
	bounded = line_minmax(I, &min, &max); /* frees I */
	equal = value_eq(min, max);
	mpz_cdiv_q(min, min, d);
	mpz_fdiv_q(max, max, d);

	if (bounded && value_gt(min, max)) {
	    /* Never zero */
	    if (p->size == 3) {
		value_clear(e->d);
		*e = p->arr[2];
	    } else {
		evalue_set_si(e, 0, 1);
		r = 1;
	    }
	    free_evalue_refs(&(p->arr[1]));
	    free_evalue_refs(&(p->arr[0]));
	    free(p);
	    value_clear(d);
	    value_clear(min);
	    value_clear(max);
	    Matrix_Free(T);
	    return r ? r : evalue_range_reduction_in_domain(e, D);
	} else if (bounded && equal) {
	    /* Always zero */
	    if (p->size == 3)
		free_evalue_refs(&(p->arr[2]));
	    value_clear(e->d);
	    *e = p->arr[1];
	    free_evalue_refs(&(p->arr[0]));
	    free(p);
	    value_clear(d);
	    value_clear(min);
	    value_clear(max);
	    Matrix_Free(T);
	    return evalue_range_reduction_in_domain(e, D);
	} else if (bounded && value_eq(min, max)) {
	    /* zero for a single value */
	    Polyhedron *E;
	    Matrix *M = Matrix_Alloc(1, D->Dimension+2);
	    Vector_Copy(T->p[0], M->p[0]+1, D->Dimension+1);
	    value_multiply(min, min, d);
	    value_subtract(M->p[0][D->Dimension+1],
			    M->p[0][D->Dimension+1], min);
	    E = DomainAddConstraints(D, M, 0);
	    value_clear(d);
	    value_clear(min);
	    value_clear(max);
	    Matrix_Free(T);
	    Matrix_Free(M);
	    r = evalue_range_reduction_in_domain(&p->arr[1], E);
	    if (p->size == 3)
		r |= evalue_range_reduction_in_domain(&p->arr[2], D);
	    Domain_Free(E);
	    _reduce_evalue(&p->arr[0].x.p->arr[0], 0, 1);
	    return r;
	}

	value_clear(d);
	value_clear(min);
	value_clear(max);
	Matrix_Free(T);
	_reduce_evalue(&p->arr[0].x.p->arr[0], 0, 1);
    }

    i = p->type == relation ? 1 : 
	p->type == fractional ? 1 : 0;
    for (; i<p->size; i++)
	r |= evalue_range_reduction_in_domain(&p->arr[i], D);

    if (p->type != fractional) {
	if (r && p->type == polynomial) {
	    evalue f;
	    value_init(f.d);
	    value_set_si(f.d, 0);
	    f.x.p = new_enode(polynomial, 2, p->pos);
	    evalue_set_si(&f.x.p->arr[0], 0, 1);
	    evalue_set_si(&f.x.p->arr[1], 1, 1);
	    reorder_terms_about(p, &f);
	    value_clear(e->d);
	    *e = p->arr[0];
	    free(p);
	}
	return r;
    }

    value_init(d);
    value_init(min);
    value_init(max);
    fractional_minimal_coefficients(p);
    I = polynomial_projection(p, D, &d, NULL);
    bounded = line_minmax(I, &min, &max); /* frees I */
    mpz_fdiv_q(min, min, d);
    mpz_fdiv_q(max, max, d);
    value_subtract(d, max, min);

    if (bounded && value_eq(min, max)) {
	evalue inc;
	value_init(inc.d);
	value_init(inc.x.n);
	value_set_si(inc.d, 1);
	value_oppose(inc.x.n, min);
	eadd(&inc, &p->arr[0]);
	reorder_terms_about(p, &p->arr[0]); /* frees arr[0] */
	value_clear(e->d);
	*e = p->arr[1];
	free(p);
	free_evalue_refs(&inc);
	r = 1;
    } else if (bounded && value_one_p(d) && p->size > 3) {
	/* replace {g}^2 by -(g-min)^2 + (2{g}+1)*(g-min) - {g}
	 * See pages 199-200 of PhD thesis.
	 */
	evalue rem;
	evalue inc;
	evalue t;
	evalue f;
	evalue factor;
	value_init(rem.d);
	value_set_si(rem.d, 0);
	rem.x.p = new_enode(fractional, 3, -1);
	evalue_copy(&rem.x.p->arr[0], &p->arr[0]);
	value_clear(rem.x.p->arr[1].d);
	value_clear(rem.x.p->arr[2].d);
	rem.x.p->arr[1] = p->arr[1];
	rem.x.p->arr[2] = p->arr[2];
	for (i = 3; i < p->size; ++i)
	    p->arr[i-2] = p->arr[i];
	p->size -= 2;

	value_init(inc.d);
	value_init(inc.x.n);
	value_set_si(inc.d, 1);
	value_oppose(inc.x.n, min);

	value_init(t.d);
	evalue_copy(&t, &p->arr[0]);
	eadd(&inc, &t);

	value_init(f.d);
	value_set_si(f.d, 0);
	f.x.p = new_enode(fractional, 3, -1);
	evalue_copy(&f.x.p->arr[0], &p->arr[0]);
	evalue_set_si(&f.x.p->arr[1], 1, 1);
	evalue_set_si(&f.x.p->arr[2], 2, 1);

	value_init(factor.d);
	evalue_set_si(&factor, -1, 1);
	emul(&t, &factor);

	eadd(&f, &factor);
	emul(&t, &factor);

	value_clear(f.x.p->arr[1].x.n);
	value_clear(f.x.p->arr[2].x.n);
	evalue_set_si(&f.x.p->arr[1], 0, 1);
	evalue_set_si(&f.x.p->arr[2], -1, 1);
	eadd(&f, &factor);

	if (r) {
	    reorder_terms(&rem);
	    reorder_terms(e);
	}

	emul(&factor, e);
	eadd(&rem, e);

	free_evalue_refs(&inc);
	free_evalue_refs(&t);
	free_evalue_refs(&f);
	free_evalue_refs(&factor);
	free_evalue_refs(&rem);

	evalue_range_reduction_in_domain(e, D);

	r = 1;
    } else {
	_reduce_evalue(&p->arr[0], 0, 1);
	if (r)
	    reorder_terms(e);
    }

    value_clear(d);
    value_clear(min);
    value_clear(max);

    return r;
}

void evalue_range_reduction(evalue *e)
{
    int i;
    if (value_notzero_p(e->d) || e->x.p->type != partition)
	return;

    for (i = 0; i < e->x.p->size/2; ++i)
	if (evalue_range_reduction_in_domain(&e->x.p->arr[2*i+1],
			     EVALUE_DOMAIN(e->x.p->arr[2*i]))) {
	    reduce_evalue(&e->x.p->arr[2*i+1]);

	    if (EVALUE_IS_ZERO(e->x.p->arr[2*i+1])) {
		free_evalue_refs(&e->x.p->arr[2*i+1]);
		Domain_Free(EVALUE_DOMAIN(e->x.p->arr[2*i]));
		value_clear(e->x.p->arr[2*i].d);
		e->x.p->size -= 2;
		e->x.p->arr[2*i] = e->x.p->arr[e->x.p->size];
		e->x.p->arr[2*i+1] = e->x.p->arr[e->x.p->size+1];
		--i;
	    }
	}
}

/* Frees EP 
 */
Enumeration* partition2enumeration(evalue *EP)
{
    int i;
    Enumeration *en, *res = NULL;

    if (EVALUE_IS_ZERO(*EP)) {
	free(EP);
	return res;
    }

    for (i = 0; i < EP->x.p->size/2; ++i) {
	assert(EP->x.p->pos == EVALUE_DOMAIN(EP->x.p->arr[2*i])->Dimension);
	en = (Enumeration *)malloc(sizeof(Enumeration));
	en->next = res;
	res = en;
	res->ValidityDomain = EVALUE_DOMAIN(EP->x.p->arr[2*i]);
	value_clear(EP->x.p->arr[2*i].d);
	res->EP = EP->x.p->arr[2*i+1];
    }
    free(EP->x.p);
    value_clear(EP->d);
    free(EP);
    return res;
}

int evalue_frac2floor_in_domain3(evalue *e, Polyhedron *D, int shift)
{
    enode *p;
    int r = 0;
    int i;
    Polyhedron *I;
    Value d, min;
    evalue fl;

    if (value_notzero_p(e->d))
	return r;

    p = e->x.p;

    i = p->type == relation ? 1 : 
	p->type == fractional ? 1 : 0;
    for (; i<p->size; i++)
	r |= evalue_frac2floor_in_domain3(&p->arr[i], D, shift);

    if (p->type != fractional) {
	if (r && p->type == polynomial) {
	    evalue f;
	    value_init(f.d);
	    value_set_si(f.d, 0);
	    f.x.p = new_enode(polynomial, 2, p->pos);
	    evalue_set_si(&f.x.p->arr[0], 0, 1);
	    evalue_set_si(&f.x.p->arr[1], 1, 1);
	    reorder_terms_about(p, &f);
	    value_clear(e->d);
	    *e = p->arr[0];
	    free(p);
	}
	return r;
    }

    if (shift) {
	value_init(d);
	I = polynomial_projection(p, D, &d, NULL);

	/*
	Polyhedron_Print(stderr, P_VALUE_FMT, I);
	*/

	assert(I->NbEq == 0); /* Should have been reduced */

	/* Find minimum */
	for (i = 0; i < I->NbConstraints; ++i)
	    if (value_pos_p(I->Constraint[i][1]))
		break;

	if (i < I->NbConstraints) {
	    value_init(min);
	    value_oppose(I->Constraint[i][2], I->Constraint[i][2]);
	    mpz_cdiv_q(min, I->Constraint[i][2], I->Constraint[i][1]);
	    if (value_neg_p(min)) {
		evalue offset;
		mpz_fdiv_q(min, min, d);
		value_init(offset.d);
		value_set_si(offset.d, 1);
		value_init(offset.x.n);
		value_oppose(offset.x.n, min);
		eadd(&offset, &p->arr[0]);
		free_evalue_refs(&offset);
	    }
	    value_clear(min);
	}

	Polyhedron_Free(I);
	value_clear(d);
    }

    value_init(fl.d);
    value_set_si(fl.d, 0);
    fl.x.p = new_enode(flooring, 3, -1);
    evalue_set_si(&fl.x.p->arr[1], 0, 1);
    evalue_set_si(&fl.x.p->arr[2], -1, 1);
    evalue_copy(&fl.x.p->arr[0], &p->arr[0]);

    eadd(&fl, &p->arr[0]);
    reorder_terms_about(p, &p->arr[0]);
    value_clear(e->d);
    *e = p->arr[1];
    free(p);
    free_evalue_refs(&fl);

    return 1;
}

int evalue_frac2floor_in_domain(evalue *e, Polyhedron *D)
{
    return evalue_frac2floor_in_domain3(e, D, 1);
}

void evalue_frac2floor2(evalue *e, int shift)
{
    int i;
    if (value_notzero_p(e->d) || e->x.p->type != partition) {
	if (!shift) {
	    if (evalue_frac2floor_in_domain3(e, NULL, 0))
		reduce_evalue(e);
	}
	return;
    }

    for (i = 0; i < e->x.p->size/2; ++i)
	if (evalue_frac2floor_in_domain3(&e->x.p->arr[2*i+1],
					EVALUE_DOMAIN(e->x.p->arr[2*i]), shift))
	    reduce_evalue(&e->x.p->arr[2*i+1]);
}

void evalue_frac2floor(evalue *e)
{
    evalue_frac2floor2(e, 1);
}

static Matrix *esum_add_constraint(int nvar, Polyhedron *D, Matrix *C,
				   Vector *row)
{
    int nr, nc;
    int i;
    int nparam = D->Dimension - nvar;

    if (C == 0) {
	nr = D->NbConstraints + 2;
	nc = D->Dimension + 2 + 1;
	C = Matrix_Alloc(nr, nc);
	for (i = 0; i < D->NbConstraints; ++i) {
	    Vector_Copy(D->Constraint[i], C->p[i], 1 + nvar);
	    Vector_Copy(D->Constraint[i] + 1 + nvar, C->p[i] + 1 + nvar + 1,
			D->Dimension + 1 - nvar);
	}
    } else {
	Matrix *oldC = C;
	nr = C->NbRows + 2;
	nc = C->NbColumns + 1;
	C = Matrix_Alloc(nr, nc);
	for (i = 0; i < oldC->NbRows; ++i) {
	    Vector_Copy(oldC->p[i], C->p[i], 1 + nvar);
	    Vector_Copy(oldC->p[i] + 1 + nvar, C->p[i] + 1 + nvar + 1,
			oldC->NbColumns - 1 - nvar);
	}
    }
    value_set_si(C->p[nr-2][0], 1);
    value_set_si(C->p[nr-2][1 + nvar], 1);
    value_set_si(C->p[nr-2][nc - 1], -1);

    Vector_Copy(row->p, C->p[nr-1], 1 + nvar + 1);
    Vector_Copy(row->p + 1 + nvar + 1, C->p[nr-1] + C->NbColumns - 1 - nparam,
		1 + nparam);

    return C;
}

static void floor2frac_r(evalue *e, int nvar)
{
    enode *p;
    int i;
    evalue f;
    evalue *pp;

    if (value_notzero_p(e->d))
	return;

    p = e->x.p;

    assert(p->type == flooring);
    for (i = 1; i < p->size; i++)
	floor2frac_r(&p->arr[i], nvar);

    for (pp = &p->arr[0]; value_zero_p(pp->d); pp = &pp->x.p->arr[0]) {
	assert(pp->x.p->type == polynomial);
	pp->x.p->pos -= nvar;
    }

    value_init(f.d);
    value_set_si(f.d, 0);
    f.x.p = new_enode(fractional, 3, -1);
    evalue_set_si(&f.x.p->arr[1], 0, 1);
    evalue_set_si(&f.x.p->arr[2], -1, 1);
    evalue_copy(&f.x.p->arr[0], &p->arr[0]);

    eadd(&f, &p->arr[0]);
    reorder_terms_about(p, &p->arr[0]);
    value_clear(e->d);
    *e = p->arr[1];
    free(p);
    free_evalue_refs(&f);
}

/* Convert flooring back to fractional and shift position
 * of the parameters by nvar
 */
static void floor2frac(evalue *e, int nvar)
{
    floor2frac_r(e, nvar);
    reduce_evalue(e);
}

evalue *esum_over_domain_cst(int nvar, Polyhedron *D, Matrix *C)
{
    evalue *t;
    int nparam = D->Dimension - nvar;

    if (C != 0) {
	C = Matrix_Copy(C);
	D = Constraints2Polyhedron(C, 0);
	Matrix_Free(C);
    }

    t = barvinok_enumerate_e(D, 0, nparam, 0);

    /* Double check that D was not unbounded. */
    assert(!(value_pos_p(t->d) && value_neg_p(t->x.n)));

    if (C != 0)
	Polyhedron_Free(D);

    return t;
}

evalue *esum_over_domain(evalue *e, int nvar, Polyhedron *D, 
			  Matrix *C)
{
    Vector *row = NULL;
    int i;
    evalue *res;
    Matrix *origC;
    evalue *factor = NULL;
    evalue cum;

    if (EVALUE_IS_ZERO(*e))
	return 0;

    if (D->next) {
	Polyhedron *DD = Disjoint_Domain(D, 0, 0);
	Polyhedron *Q;

	Q = DD;
	DD = Q->next;
	Q->next = 0;

	res = esum_over_domain(e, nvar, Q, C);
	Polyhedron_Free(Q);

	for (Q = DD; Q; Q = DD) {
	    evalue *t;

	    DD = Q->next;
	    Q->next = 0;

	    t = esum_over_domain(e, nvar, Q, C);
	    Polyhedron_Free(Q);

	    if (!res)
		res = t;
	    else if (t) {
		eadd(t, res);
		free_evalue_refs(t);
		free(t);
	    }
	}
	return res;
    }

    if (value_notzero_p(e->d)) {
	evalue *t;

	t = esum_over_domain_cst(nvar, D, C);

	if (!EVALUE_IS_ONE(*e))
	    emul(e, t);

	return t;
    }

    switch (e->x.p->type) {
    case flooring: {
	evalue *pp = &e->x.p->arr[0];

	if (pp->x.p->pos > nvar) {
	    /* remainder is independent of the summated vars */
	    evalue f;
	    evalue *t;

	    value_init(f.d);
	    evalue_copy(&f, e);
	    floor2frac(&f, nvar);

	    t = esum_over_domain_cst(nvar, D, C);

	    emul(&f, t);

	    free_evalue_refs(&f);

	    return t;
	}

	row = Vector_Alloc(1 + D->Dimension + 1 + 1);
	poly_denom(pp, &row->p[1 + nvar]);
	value_set_si(row->p[0], 1);
	for (pp = &e->x.p->arr[0]; value_zero_p(pp->d); 
				   pp = &pp->x.p->arr[0]) {
	    int pos;
	    assert(pp->x.p->type == polynomial);
	    pos = pp->x.p->pos;
	    if (pos >= 1 + nvar)
		++pos;
	    value_assign(row->p[pos], row->p[1+nvar]);
	    value_division(row->p[pos], row->p[pos], pp->x.p->arr[1].d);
	    value_multiply(row->p[pos], row->p[pos], pp->x.p->arr[1].x.n);
	}
	value_assign(row->p[1 + D->Dimension + 1], row->p[1+nvar]);
	value_division(row->p[1 + D->Dimension + 1],
		       row->p[1 + D->Dimension + 1],
		       pp->d);
	value_multiply(row->p[1 + D->Dimension + 1],
		       row->p[1 + D->Dimension + 1],
		       pp->x.n);
	value_oppose(row->p[1 + nvar], row->p[1 + nvar]);
	break;
    }
    case polynomial: {
	int pos = e->x.p->pos;

	if (pos > nvar) {
	    factor = ALLOC(evalue);
	    value_init(factor->d);
	    value_set_si(factor->d, 0);
	    factor->x.p = new_enode(polynomial, 2, pos - nvar);
	    evalue_set_si(&factor->x.p->arr[0], 0, 1);
	    evalue_set_si(&factor->x.p->arr[1], 1, 1);
	    break;
	}

	row = Vector_Alloc(1 + D->Dimension + 1 + 1);
	for (i = 0; i < D->NbRays; ++i)
	    if (value_notzero_p(D->Ray[i][pos]))
		break;
	assert(i < D->NbRays);
	if (value_neg_p(D->Ray[i][pos])) {
	    factor = ALLOC(evalue);
	    value_init(factor->d);
	    evalue_set_si(factor, -1, 1);
	}
	value_set_si(row->p[0], 1);
	value_set_si(row->p[pos], 1);
	value_set_si(row->p[1 + nvar], -1);
	break;
    }
    default:
	assert(0);
    }

    i = type_offset(e->x.p);

    res = esum_over_domain(&e->x.p->arr[i], nvar, D, C);
    ++i;

    if (factor) {
	value_init(cum.d);
	evalue_copy(&cum, factor);
    }

    origC = C;
    for (; i < e->x.p->size; ++i) {
	evalue *t;
	if (row) {
	    Matrix *prevC = C;
	    C = esum_add_constraint(nvar, D, C, row);
	    if (prevC != origC)
		Matrix_Free(prevC);
	}
	/*
	if (row)
	    Vector_Print(stderr, P_VALUE_FMT, row);
	if (C)
	    Matrix_Print(stderr, P_VALUE_FMT, C);
	*/
	t = esum_over_domain(&e->x.p->arr[i], nvar, D, C);

	if (t && factor)
	    emul(&cum, t);

	if (!res)
	    res = t;
	else if (t) {
	    eadd(t, res);
	    free_evalue_refs(t);
	    free(t);
	}
	if (factor && i+1 < e->x.p->size)
	    emul(factor, &cum);
    }
    if (C != origC)
	Matrix_Free(C);

    if (factor) {
	free_evalue_refs(factor);
	free_evalue_refs(&cum);
	free(factor);
    }

    if (row)
	Vector_Free(row);

    reduce_evalue(res);

    return res;
}

evalue *esum(evalue *e, int nvar)
{
    int i;
    evalue *res = ALLOC(evalue);
    value_init(res->d);

    assert(nvar >= 0);
    if (nvar == 0 || EVALUE_IS_ZERO(*e)) {
	evalue_copy(res, e);
	return res;
    }

    evalue_set_si(res, 0, 1);

    assert(value_zero_p(e->d));
    assert(e->x.p->type == partition);

    for (i = 0; i < e->x.p->size/2; ++i) {
	evalue *t;
	t = esum_over_domain(&e->x.p->arr[2*i+1], nvar,
			     EVALUE_DOMAIN(e->x.p->arr[2*i]), 0);
	eadd(t, res);
	free_evalue_refs(t);
	free(t);
    }

    reduce_evalue(res);

    return res;
}

/* Initial silly implementation */
void eor(evalue *e1, evalue *res)
{
    evalue E;
    evalue mone;
    value_init(E.d);
    value_init(mone.d);
    evalue_set_si(&mone, -1, 1);

    evalue_copy(&E, res);
    eadd(e1, &E);
    emul(e1, res);
    emul(&mone, res);
    eadd(&E, res);

    free_evalue_refs(&E); 
    free_evalue_refs(&mone);
}

/* computes denominator of polynomial evalue 
 * d should point to a value initialized to 1
 */
void evalue_denom(const evalue *e, Value *d)
{
    int i, offset;

    if (value_notzero_p(e->d)) {
	value_lcm(*d, e->d, d);
	return;
    }
    assert(e->x.p->type == polynomial ||
	   e->x.p->type == fractional ||
	   e->x.p->type == flooring);
    offset = type_offset(e->x.p);
    for (i = e->x.p->size-1; i >= offset; --i)
	evalue_denom(&e->x.p->arr[i], d);
}

/* Divides the evalue e by the integer n */
void evalue_div(evalue * e, Value n)
{
    int i, offset;

    if (value_notzero_p(e->d)) {
	Value gc;
	value_init(gc);
	value_multiply(e->d, e->d, n);
	Gcd(e->x.n, e->d, &gc);
	if (value_notone_p(gc)) {
	    value_division(e->d, e->d, gc);
	    value_division(e->x.n, e->x.n, gc);
	}
	value_clear(gc);
	return;
    }
    if (e->x.p->type == partition) {
	for (i = 0; i < e->x.p->size/2; ++i)
	    evalue_div(&e->x.p->arr[2*i+1], n);
	return;
    }
    offset = type_offset(e->x.p);
    for (i = e->x.p->size-1; i >= offset; --i)
	evalue_div(&e->x.p->arr[i], n);
}

static void evalue_frac2polynomial_r(evalue *e, int *signs, int sign, int in_frac)
{
    int i, offset;
    Value d;
    enode *p;
    int sign_odd = sign;

    if (value_notzero_p(e->d)) {
	if (in_frac && sign * value_sign(e->x.n) < 0) {
	    value_set_si(e->x.n, 0);
	    value_set_si(e->d, 1);
	}
	return;
    }

    if (e->x.p->type == relation) {
	for (i = e->x.p->size-1; i >= 1; --i)
	    evalue_frac2polynomial_r(&e->x.p->arr[i], signs, sign, in_frac);
	return;
    }

    if (e->x.p->type == polynomial)
	sign_odd *= signs[e->x.p->pos-1];
    offset = type_offset(e->x.p);
    evalue_frac2polynomial_r(&e->x.p->arr[offset], signs, sign, in_frac);
    in_frac |= e->x.p->type == fractional;
    for (i = e->x.p->size-1; i > offset; --i)
	evalue_frac2polynomial_r(&e->x.p->arr[i], signs,
				 (i - offset) % 2 ? sign_odd : sign, in_frac);

    if (e->x.p->type != fractional)
	return;

    /* replace { a/m } by (m-1)/m if sign != 0
     * and by (m-1)/(2m) if sign == 0
     */
    value_init(d);
    value_set_si(d, 1);
    evalue_denom(&e->x.p->arr[0], &d);
    free_evalue_refs(&e->x.p->arr[0]);
    value_init(e->x.p->arr[0].d);
    value_init(e->x.p->arr[0].x.n);
    if (sign == 0)
	value_addto(e->x.p->arr[0].d, d, d);
    else
	value_assign(e->x.p->arr[0].d, d);
    value_decrement(e->x.p->arr[0].x.n, d);
    value_clear(d);

    p = e->x.p;
    reorder_terms_about(p, &p->arr[0]);
    value_clear(e->d);
    *e = p->arr[1];
    free(p);
}

/* Approximate the evalue in fractional representation by a polynomial.
 * If sign > 0, the result is an upper bound;
 * if sign < 0, the result is a lower bound;
 * if sign = 0, the result is an intermediate approximation.
 */
void evalue_frac2polynomial(evalue *e, int sign, unsigned MaxRays)
{
    int i, j, k, dim;
    int *signs;

    if (value_notzero_p(e->d))
	return;
    assert(e->x.p->type == partition);
    /* make sure all variables in the domains have a fixed sign */
    if (sign)
	evalue_split_domains_into_orthants(e, MaxRays);

    assert(e->x.p->size >= 2);
    dim = EVALUE_DOMAIN(e->x.p->arr[0])->Dimension;

    signs = alloca(sizeof(int) * dim);

    for (i = 0; i < e->x.p->size/2; ++i) {
	Polyhedron *D = EVALUE_DOMAIN(e->x.p->arr[2*i]);
	POL_ENSURE_VERTICES(D);
	for (j = 0; j < dim; ++j) {
	    signs[j] = 0;
	    if (!sign)
		continue;
	    for (k = 0; k < D->NbRays; ++k) {
		signs[j] = value_sign(D->Ray[k][1+j]);
		if (signs[j])
		    break;
	    }
	}
	evalue_frac2polynomial_r(&e->x.p->arr[2*i+1], signs, sign, 0);
    }
}

/* Split the domains of e (which is assumed to be a partition)
 * such that each resulting domain lies entirely in one orthant.
 */
void evalue_split_domains_into_orthants(evalue *e, unsigned MaxRays)
{
    int i, dim;
    assert(value_zero_p(e->d));
    assert(e->x.p->type == partition);
    assert(e->x.p->size >= 2);
    dim = EVALUE_DOMAIN(e->x.p->arr[0])->Dimension;

    for (i = 0; i < dim; ++i) {
	evalue split;
	Matrix *C, *C2;
	C = Matrix_Alloc(1, 1 + dim + 1);
	value_set_si(C->p[0][0], 1);
	value_init(split.d);
	value_set_si(split.d, 0);
	split.x.p = new_enode(partition, 4, dim);
	value_set_si(C->p[0][1+i], 1);
	C2 = Matrix_Copy(C);
	EVALUE_SET_DOMAIN(split.x.p->arr[0], Constraints2Polyhedron(C2, MaxRays));
	Matrix_Free(C2);
	evalue_set_si(&split.x.p->arr[1], 1, 1);
	value_set_si(C->p[0][1+i], -1);
	value_set_si(C->p[0][1+dim], -1);
	EVALUE_SET_DOMAIN(split.x.p->arr[2], Constraints2Polyhedron(C, MaxRays));
	evalue_set_si(&split.x.p->arr[3], 1, 1);
	emul(&split, e);
	free_evalue_refs(&split);
	Matrix_Free(C);
    }
}

static Matrix *find_fractional_with_max_periods(evalue *e, Polyhedron *D,
						int max_periods,
						Value *min, Value *max)
{
    Matrix *T;
    Value d;
    int i;

    if (value_notzero_p(e->d))
	return NULL;

    if (e->x.p->type == fractional) {
	Polyhedron *I;
	int bounded;

	value_init(d);
	I = polynomial_projection(e->x.p, D, &d, &T);
	bounded = line_minmax(I, min, max); /* frees I */
	if (bounded) {
	    Value mp;
	    value_init(mp);
	    value_set_si(mp, max_periods);
	    mpz_fdiv_q(*min, *min, d);
	    mpz_fdiv_q(*max, *max, d);
	    value_assign(T->p[1][D->Dimension], d);
	    value_subtract(d, *max, *min);
	    if (value_ge(d, mp)) {
		Matrix_Free(T);
		T = NULL;
	    }
	    value_clear(mp);
	} else {
	    Matrix_Free(T);
	    T = NULL;
	}
	value_clear(d);
	if (T)
	    return T;
    }

    for (i = type_offset(e->x.p); i < e->x.p->size; ++i)
	if ((T = find_fractional_with_max_periods(&e->x.p->arr[i], D, max_periods,
						  min, max)))
	    return T;

    return NULL;
}

/* Look for fractional parts that can be removed by splitting the corresponding
 * domain into at most max_periods parts.
 * We use a very simply strategy that looks for the first fractional part
 * that satisfies the condition, performs the split and then continues
 * looking for other fractional parts in the split domains until no
 * such fractional part can be found anymore.
 */
void evalue_split_periods(evalue *e, int max_periods, unsigned int MaxRays)
{
    int i, j, n;
    Value min;
    Value max;
    Value d;

    if (EVALUE_IS_ZERO(*e))
	return;
    if (value_notzero_p(e->d) || e->x.p->type != partition) {
	fprintf(stderr,
		"WARNING: evalue_split_periods called on incorrect evalue type\n");
	return;
    }

    value_init(min);
    value_init(max);
    value_init(d);

    for (i = 0; i < e->x.p->size/2; ++i) {
	enode *p;
	Matrix *T = NULL;
	Matrix *M;
	Polyhedron *D = EVALUE_DOMAIN(e->x.p->arr[2*i]);
	Polyhedron *E;
	T = find_fractional_with_max_periods(&e->x.p->arr[2*i+1], D, max_periods,
					     &min, &max);
	if (!T)
	    continue;

	M = Matrix_Alloc(2, 2+D->Dimension);

	value_subtract(d, max, min);
	n = VALUE_TO_INT(d)+1;

	value_set_si(M->p[0][0], 1);
	Vector_Copy(T->p[0], M->p[0]+1, D->Dimension+1);
	value_multiply(d, max, T->p[1][D->Dimension]);
	value_subtract(M->p[0][1+D->Dimension], M->p[0][1+D->Dimension], d);
	value_set_si(d, -1);
	value_set_si(M->p[1][0], 1);
	Vector_Scale(T->p[0], M->p[1]+1, d, D->Dimension+1);
	value_addmul(M->p[1][1+D->Dimension], max, T->p[1][D->Dimension]);
	value_addto(M->p[1][1+D->Dimension], M->p[1][1+D->Dimension],
		    T->p[1][D->Dimension]);
	value_decrement(M->p[1][1+D->Dimension], M->p[1][1+D->Dimension]);

	p = new_enode(partition, e->x.p->size + (n-1)*2, e->x.p->pos);
	for (j = 0; j < 2*i; ++j) {
	    value_clear(p->arr[j].d);
	    p->arr[j] = e->x.p->arr[j];
	}
	for (j = 2*i+2; j < e->x.p->size; ++j) {
	    value_clear(p->arr[j+2*(n-1)].d);
	    p->arr[j+2*(n-1)] = e->x.p->arr[j];
	}
	for (j = n-1; j >= 0; --j) {
	    if (j == 0) {
		value_clear(p->arr[2*i+1].d);
		p->arr[2*i+1] = e->x.p->arr[2*i+1];
	    } else
		evalue_copy(&p->arr[2*(i+j)+1], &e->x.p->arr[2*i+1]);
	    if (j != n-1) {
		value_subtract(M->p[1][1+D->Dimension], M->p[1][1+D->Dimension],
			       T->p[1][D->Dimension]);
		value_addto(M->p[0][1+D->Dimension], M->p[0][1+D->Dimension],
			    T->p[1][D->Dimension]);
	    }
	    E = DomainAddConstraints(D, M, MaxRays);
	    EVALUE_SET_DOMAIN(p->arr[2*(i+j)], E);
	    if (evalue_range_reduction_in_domain(&p->arr[2*(i+j)+1], E))
		reduce_evalue(&p->arr[2*(i+j)+1]);
	}
	value_clear(e->x.p->arr[2*i].d);
	Domain_Free(D);
	Matrix_Free(M);
	Matrix_Free(T);
	free(e->x.p);
	e->x.p = p;
	--i;
    }

    value_clear(d);
    value_clear(min);
    value_clear(max);
}

void evalue_extract_affine(const evalue *e, Value *coeff, Value *cst, Value *d)
{
    value_set_si(*d, 1);
    evalue_denom(e, d);
    for ( ; value_zero_p(e->d); e = &e->x.p->arr[0]) {
	assert(e->x.p->type == polynomial);
	assert(e->x.p->size == 2);
	evalue *c = &e->x.p->arr[1];
	value_multiply(coeff[e->x.p->pos-1], *d, c->x.n);
	value_division(coeff[e->x.p->pos-1], coeff[e->x.p->pos-1], c->d);
    }
    value_multiply(*cst, *d, e->x.n);
    value_division(*cst, *cst, e->d);
}

/* returns an evalue that corresponds to
 *
 * c/den x_param
 */
static evalue *term(int param, Value c, Value den)
{
    evalue *EP = ALLOC(evalue);
    value_init(EP->d);
    value_set_si(EP->d,0);
    EP->x.p = new_enode(polynomial, 2, param + 1);
    evalue_set_si(&EP->x.p->arr[0], 0, 1);
    value_init(EP->x.p->arr[1].x.n);
    value_assign(EP->x.p->arr[1].d, den);
    value_assign(EP->x.p->arr[1].x.n, c);
    return EP;
}

evalue *affine2evalue(Value *coeff, Value denom, int nvar)
{
    int i;
    evalue *E = ALLOC(evalue);
    value_init(E->d);
    evalue_set(E, coeff[nvar], denom);
    for (i = 0; i < nvar; ++i) {
	evalue *t = term(i, coeff[i], denom);
	eadd(t, E);
	free_evalue_refs(t);
	free(t);
    }
    return E;
}
