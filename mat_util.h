#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>

#ifdef NTL_STD_CXX
using namespace NTL;
#endif

int lex_cmp(vec_ZZ& a, vec_ZZ& b);
void lex_order_rows(mat_ZZ& mat);
