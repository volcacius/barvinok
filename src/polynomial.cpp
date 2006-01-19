/* 
 *	Bernstein Expansion
 *
 *	- polynomial class
 */

#include "polynomial.h"

polynomial::polynomial(void)
{
}

polynomial::polynomial(ex &e)
{
	*this += e;
}

polynomial::polynomial(const ex &e)
{
	*this += e;
}


polynomial& polynomial::operator=(ex &e)
{
	*this = polynomial::polynomial(e);

	return *this;
}


// term i of the polynomial
ex polynomial::term(unsigned int i)
{
	ex e = *this;
	if(!is_a<add>(e)) {
		return e;
	}
	if(this->nops() == 0) {
		return e;
	}
	return this->op(i);
}


// number of terms in the polynomial
unsigned int polynomial::nbTerms(void)
{
	ex e = *this;
	if(!is_a<add>(e)) {
		return 1;
	}

	if(this->nops() == 0) {
		return 1;
	}

	return this->nops();
}
