#ifndef EV_OPERATIONS_H
#define EV_OPERATIONS_H

#include <polylib/polylibgmp.h>

#define eadd new_eadd

void evalue_set_si(evalue *ev, int n, int d);
void eadd(evalue *e1,evalue *res);
void emul (evalue *e1, evalue *res );
#endif
