/* Erin Parker (parker@cs.unc.edu), March 2004 */

#include <stdlib.h>
#include "dfa.h"


/* Functions defined in construction.c */
DFA* build_DFA_eq(int, int*, int, int*);
DFA* build_DFA_ineq(int, int*, int, int*);

/* Function defined in count.c */
void count_accepting_paths(DFA*, int, int);
