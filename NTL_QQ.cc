#include <assert.h>
#include <barvinok/NTL_QQ.h>

NTL_vector_impl(QQ,vec_QQ);

vec_QQ& operator *= (vec_QQ& a, const ZZ& b)
{
    for (int i = 0; i < a.length(); ++i)
	a[i] *= b;
    return a;
}

vec_QQ& operator *= (vec_QQ& a, const QQ& b)
{
    for (int i = 0; i < a.length(); ++i)
	a[i] *= b;
    return a;
}

std::ostream& operator<< (std::ostream& os, const QQ& q)
{
    os << q.n << "/" << q.d;
    return os;
}

std::istream& operator>> (std::istream& is, QQ& q)
{
    char slash;
    is >> q.n >> slash >> q.d;
    assert(slash == '/');
    return is;
}

NTL_io_vector_impl(QQ,vec_QQ);
