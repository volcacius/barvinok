#include <iostream>
#include <gmp.h>
extern "C" {
#include <polylib/polylibgmp.h>
}
#include <barvinok/evalue.h>

void evalue_print(std::ostream& o, evalue *e, char **p);
