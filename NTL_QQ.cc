#include <assert.h>
#include <barvinok/NTL_QQ.h>

NTL_vector_impl(QQ,vec_QQ);

std::ostream& operator<< (std::ostream& os, const QQ& q)
{
    os << q.n << "/" << q.d;
    return os;
}

std::istream& operator>> (std::istream& os, QQ& q)
{
    assert(0);
}

NTL_io_vector_impl(QQ,vec_QQ);
