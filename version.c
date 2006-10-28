#include <barvinok/util.h>
#include "config.h"
#include "version.h"

const char *barvinok_version(void)
{
    return 
	GIT_HEAD_ID"\n"
#ifdef USE_MODULO
	" +MODULO"
#else
	" -MODULO"
#endif
#ifdef USE_INCREMENTAL_BF
	" INCREMENTAL=BF"
#elif defined USE_INCREMENTAL_DF
	" INCREMENTAL=DF"
#else
	" -INCREMENTAL"
#endif
    "\n"
#ifdef HAVE_CORRECT_VERTICES
	" +CORRECT_VERTICES"
#else
	" -CORRECT_VERTICES"
#endif
#ifdef HAVE_PIPLIB
	" +PIPLIB"
#else
	" -PIPLIB"
#endif
#ifdef HAVE_OMEGA
	" +OMEGA"
#else
	" -OMEGA"
#endif
#ifdef HAVE_LIBGLPK
	" +GLPK"
#else
	" -GLPK"
#endif
#ifdef HAVE_GINAC
	" +GINAC"
#else
	" -GINAC"
#endif
    "\n"
    ;
}
