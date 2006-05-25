#include <assert.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include "mat_util.h"

int lex_cmp(vec_ZZ& a, vec_ZZ& b)
{
    assert(a.length() == b.length());

    for (int j = 0; j < a.length(); ++j)
	if (a[j] != b[j])
	    return a[j] < b[j] ? -1 : 1;
    return 0;
}

void lex_order_rows(mat_ZZ& mat)
{
    for (int i = 0; i < mat.NumRows(); ++i) {
	int m = i;
	for (int j = i+1; j < mat.NumRows(); ++j)
	    if (lex_cmp(mat[j], mat[m]) < 0)
		m = j;
	if (m != i) {
	    vec_ZZ tmp = mat[m];
	    mat[m] = mat[i];
	    mat[i] = tmp;
	}
    }
}
