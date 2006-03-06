/* 
 *	Bernstein Expansion
 *
 *	- polynomial class
 */


#include <ginac/ginac.h>


class polynomial : public GiNaC::ex
{
public:
	polynomial::polynomial(void);
	polynomial::polynomial(GiNaC::ex &);
	polynomial::polynomial(const GiNaC::ex &);

	unsigned int nbTerms(void);
	GiNaC::ex term(unsigned int i);
	polynomial& polynomial::operator=(GiNaC::ex &);
};
