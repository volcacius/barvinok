/* 
 *	Bernstein Expansion
 *
 *	- polynomial class
 */


#include <ginac/ginac.h>


class polynomial : public GiNaC::ex
{
public:
	polynomial(void);
	polynomial(GiNaC::ex &);
	polynomial(const GiNaC::ex &);

	unsigned int nbTerms(void);
	GiNaC::ex term(unsigned int i);
	polynomial& operator=(GiNaC::ex &);
};
