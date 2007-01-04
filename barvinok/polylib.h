#ifndef BARVINOK_POLYLIB_H
#define BARVINOK_POLYLIB_H

#include <gmp.h>

#if defined(__cplusplus)
extern "C" {
#endif

#include <polylib/polylibgmp.h>

#ifndef value_subtract
#define value_subtract	value_substract
#endif

#ifndef value_addmul
#define value_addmul(ref,val1,val2)						\
	    do {								\
		Value _tmp;							\
		value_init(_tmp);						\
		value_multiply(_tmp,val1,val2);				    	\
		value_addto(ref,ref,_tmp);					\
		value_clear(_tmp);						\
	    } while(0)
#endif

#undef divide
#undef value_compare

#if defined(__cplusplus)
}
#endif

#endif
