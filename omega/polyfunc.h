#include <ginac/ginac.h>
#include <omega/AST.h>
#include <omega/Relation.h>

struct PolyFunc {
    Relation	domain;
    GiNaC::ex	poly;
};

void maximize(PolyFunc *polyfunc, Map<Variable_Ref *, GiNaC::ex>& variableMap);
