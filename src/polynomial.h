/* 
 *	Bernstein Expansion
 *
 *	- polynomial class
 */


#include <ginac/ginac.h>

using namespace std;
using namespace GiNaC;


class polynomial : public ex
{
public:
	polynomial::polynomial(void);
	polynomial::polynomial(ex &);
	polynomial::polynomial(const ex &);

	unsigned int nbTerms(void);
	ex term(unsigned int i);
	polynomial& polynomial::operator=(ex &);
};
