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
    #define	BV_LEXMIN_POLYSIGN_POLYLIB	0
    #define	BV_LEXMIN_POLYSIGN_CDD		1
    #define	BV_LEXMIN_POLYSIGN_CDDF		2
    int		polysign;

    struct verify_options    verify;
};

#endif
