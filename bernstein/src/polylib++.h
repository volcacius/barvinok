#if defined(__cplusplus)
extern "C" {
#endif
#define matrix polylib_matrix
#define polynomial polylib_polynomial
#include <polylib/polylibgmp.h>
#undef matrix
#undef polynomial
#undef value_compare
#undef divide
#if defined(__cplusplus)
}
#endif
