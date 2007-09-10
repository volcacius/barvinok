#include <barvinok/evalue.h>

#if defined(__cplusplus)
extern "C" {
#endif

struct barvinok_options;

evalue *euler_summate(evalue *e, int nvar, struct barvinok_options *options);

#if defined(__cplusplus)
}
#endif
