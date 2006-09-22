#include <iostream>
#include <gmp.h>
extern "C" {
#include <polylib/polylibgmp.h>
}
#include <barvinok/evalue.h>

void evalue_denom(evalue *e, Value *d);
void evalue_print(std::ostream& o, evalue *e, char **p);
