#include <ginac/ginac.h>
#include <bernstein/piecewise_lst.h>
#include <barvinok/evalue.h>

struct barvinok_options;

bernstein::piecewise_lst *evalue_range_propagation(
				bernstein::piecewise_lst *pl_all,
				const evalue *e, const GiNaC::exvector& params,
				barvinok_options *options);
