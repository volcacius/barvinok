#ifndef LEXMIN_H
#define LEXMIN_H

#include "verify.h"

struct lexmin_options {
		/* Check for integer points in domain
		 */
    #define	BV_LEXMIN_EMPTINESS_CHECK_NONE		0
    #define	BV_LEXMIN_EMPTINESS_CHECK_SAMPLE	1
    #define	BV_LEXMIN_EMPTINESS_CHECK_COUNT		2
    int		emptiness_check;
    int		reduce;

    struct verify_options    verify;
};

#endif
