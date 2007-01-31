#ifndef NTL_QQ_H
#define NTL_QQ_H

#include <NTL/ZZ.h>
#include <NTL/vector.h>

#ifdef NTL_STD_CXX
using namespace NTL;
#endif

struct QQ {
    ZZ	n;
    ZZ	d;

    QQ() {}
    QQ(int n, int d) {
	this->n = n;
	this->d = d;
    }

    QQ& operator += (const QQ& a) {
	/* This is not thread-safe, but neither is NTL */
	static ZZ tmp;
	n *= a.d;
	mul(tmp, d, a.n);
	n += tmp;
	d *= a.d;
	GCD(tmp, n, d);
	n /= tmp;
	d /= tmp;
	return *this;
    }

    QQ& operator *= (const QQ& a) {
	n *= a.n;
	d *= a.d;
	return *this;
    }

    QQ& operator *= (const ZZ& a) {
	n *= a;
	return *this;
    }
};

NTL_vector_decl(QQ,vec_QQ);

vec_QQ& operator *= (vec_QQ& a, const ZZ& b);

std::ostream& operator<< (std::ostream& os, const QQ& q);
std::istream& operator>> (std::istream& os, QQ& q);

NTL_io_vector_decl(QQ,vec_QQ);

#endif
