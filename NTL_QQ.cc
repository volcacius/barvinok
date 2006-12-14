#include <assert.h>
#include <barvinok/NTL_QQ.h>

NTL_vector_impl(QQ,vec_QQ);

std::ostream& operator<< (std::ostream& os, const QQ& q)
{
    os << q.n << "/" << q.d;
    return os;
}

std::istream& operator>> (std::istream& is, QQ& q)
{
    assert(0);
    return is;
}

NTL_io_vector_impl(QQ,vec_QQ);
