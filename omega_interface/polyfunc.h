#include <ginac/ginac.h>
#include <omega/AST.h>
#include <omega/Relation.h>
#include <barvinok/evalue.h>

struct PolyFunc {
    Relation	domain;
    GiNaC::ex	poly;
};

void maximize(PolyFunc *polyfunc, Map<Variable_Ref *, GiNaC::ex>& variableMap);
evalue *summate(PolyFunc *polyfunc, Map<Variable_Ref *, GiNaC::ex>& variableMap);
