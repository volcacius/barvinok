//-*-c++-*-
/***************************************************************************
                          ehrhartpolynom.h  -  description
                             -------------------
    begin                : Fri 23 Nov  2001
    copyright            : (C) 2001 by Kristof Beyls
    email                : Kristof.Beyls@elis.rug.ac.be
 ***************************************************************************/

#ifndef EHRHARTPOLYNOM_H
#define EHRHARTPOLYNOM_H

#include <set>
#include <map>
#include <string>
#include <deque>
/*
#define swap omega_swap
#define min omega_min
#define max omega_max
#define gcd omega_gcd
#include <omega.h>
#undef gcd
#undef swap
#undef min
#undef max
*/
#include <gmp.h>
#include <assert.h>

extern "C" {
#include <polylib/polylibgmp.h>
#include <polylib/ehrhart.h>
}

using namespace std;
/** class PeriodicNumber represents a multidimensional
 *  periodic number. For more definition and more details
 *  about periodic numbers: see "Counting Solutions to Linear
 *  and Nonlinear Constraints through Ehrhart Polynomials:
 *  Applications to Analyze and Transform Scientific Programs",
 *  By Philippe CLAUSS, in tenth ACM International Conference
 *  on Supercomputing, May 1996
 */
class PeriodicNumber {
  friend class EhrhartPolynom;
  /// parameter_names contains the ordered names of the parameters
  set<string>     parameter_names;
  /// period is the period of the variable
  map<string,int> period;
  /** data contains the multi-dimensional matrix containing
   *  the value's in the periodic number. The data
   *  is layed out column-major, with the dimensions
   *  of the matrix ordered with the alphabetical order
   *  of the corresponding parameter names.
   */
  Value*        data;
  /** stride contains the stride in the data array between
   *  2 consecutive elements of parameter p.
   *  e.g. if the parameters are "a" and "b", and the strides
   *       are 3 and 4;
   *  stride["a"] = 1 and stride["b"] = 3.
   */
  map<string,unsigned long> stride;
  unsigned long datasize;
  
public:
  /** @param parameter_names contains the names of the 
   *  parameters.
   *  @param period contains the periods of the different
   *  parameters.
   */
  PeriodicNumber(const set<string>& _parameter_names,
  		 const map<string, int>& _period) 
    : parameter_names(_parameter_names), period(_period) 
  {
    unsigned long datasize = 1;
    for(map<string,int>::const_iterator i=period.begin();
  	i!=period.end(); i++) {
      stride[(*i).first] = datasize;
      datasize*=(*i).second;
    }
    data = new Value[datasize];
    for(unsigned i=0;i<datasize;i++) {
      value_init(data[i]);
      value_set_si(data[i] , 0);
    }
    this->datasize = datasize;
  }
  PeriodicNumber() {
    parameter_names.insert("");
    period[""]=1;
    datasize=1;
    data = new Value[datasize];
    value_init(data[0]);
    value_set_si(data[0],0);
    stride[""]=1;
  }
  PeriodicNumber(const Value t) {
    parameter_names.insert("");
    period[""]=1;
    datasize=1;
    data = new Value[datasize];
    value_init(data[0]);
    value_assign(data[0],t);
    stride[""]=1;
  }
  PeriodicNumber(const PeriodicNumber& pn) 
    : parameter_names(pn.parameter_names), period(pn.period),
      stride(pn.stride), datasize(pn.datasize)
  {
    data = new Value[datasize];
    for(unsigned i=0;i<datasize; i++) {
      value_init(data[i]);
      value_assign(data[i], pn.data[i]);
    }
    /*
    memcpy(data, pn.data, sizeof(Value)*datasize);
    */
  }
  PeriodicNumber& operator=(const PeriodicNumber& pn) {
    if (*this==pn)
      return *this;
    if (data!=0) {
      for(unsigned i=0; i<datasize; i++)
	value_clear(data[i]);
    }
    delete[] data;
    parameter_names=pn.parameter_names;
    period=pn.period;
    stride=pn.stride;
    datasize=pn.datasize;

    data = new Value[datasize];
    for(unsigned i=0;i<datasize; i++) {
      value_init(data[i]);
      value_assign(data[i], pn.data[i]);
    }
    /* memcpy(data, pn.data, sizeof(Value)*datasize); */
    return *this;
  }
  bool is_zero() const {
    for(unsigned long i=0;i<datasize;i++)
      if (value_notzero_p(data[i]))
	return false;
    return true;
  }
  bool operator==(const PeriodicNumber& pn) const {
    bool t1 = (parameter_names==pn.parameter_names &&
  	       period==pn.period);
    if (!t1)
      return false;
    assert(stride==pn.stride);
    assert(datasize==pn.datasize);
    for(unsigned i=0;i<datasize; i++)
      if (value_ne( data[i], pn.data[i]))
	return false;
    return true;
    /*
    return (0==memcmp(data, pn.data, datasize*sizeof(Value)));	       
    */
  }
  inline bool operator!=(const PeriodicNumber& pn) const { return !((*this)==pn); }
  //  PeriodicNumber() : data(0) {}
  ~PeriodicNumber() {
    for(unsigned i=0;i<datasize; i++)
      value_clear(data[i]);
    delete[] data;
  }
  int dim() const { 
    assert(period.size()==parameter_names.size());
    return period.size(); 
  }
  /** operator[] returns the Value which correspends to the index
   *  @param index The index in the multidimensional matrix,
   *               representing the multi-periodic number.
   */
  Value& operator[](const map<string, int>& index);
  const Value& operator[](const map<string, int>& index) const;
  /** multiply all the values in the periodic number with v */
  PeriodicNumber operator*(const Value v) const;
  /** find the gcd of all the values in the periodic number */
  void get_gcd(Value& v) const;
  /** divide all the values in the periodic number by v.
   *  The gcd of all the values should be a multiple of v.*/
  void divide_by(const Value v);
  /* Calculate the product of 2 periodic numbers.
   *  e.g. product of [1,2]_p and [5,6]_q
   *  is   [ 5  6 ]
   *       [ 10 12]_pq
   */
  PeriodicNumber operator*(const PeriodicNumber& pn) const;
  PeriodicNumber operator+(const PeriodicNumber& pn) const;
  bool smaller_than(const PeriodicNumber& pn) const;
  string to_string() const;
  double to_double() const;
  inline operator string() const { return to_string(); }
  ::evalue to_evalue(const deque<string>& parameter_names) const;
  int get_period(const string& name) const {
    map<string,int>::const_iterator i=period.find(name);
    assert(i!=period.end());
    return (*i).second;
  }
  bool has_param(const string& name) const;
  set<string> get_params() const {
    set<string> result;
    for(map<string,int>::const_iterator i=period.begin();
	i!=period.end(); i++)
      result.insert((*i).first);
    assert(result == parameter_names);
    return result;
  }
};

/** class EhrhartPolynom is PolyAst's own representation of
 *  an Ehrhart polynomial. The first reason to create this
 *  representation is to allow an easy algorithm to add up 2
 *  Ehrhart Polynomials, since the eadd-function in polylib4.19
 *  doesn't handle all cases.
 */
class EhrhartPolynom {
  public:
  struct periodic_rational {
    PeriodicNumber first;
    Value second;
    periodic_rational(const PeriodicNumber& pn, const Value v) {
      first=pn;
      value_init(second);
      value_assign(second, v);
    }
    periodic_rational() {
      value_init(second);
    }
    ~periodic_rational() {
      value_clear(second);
    }
    periodic_rational(const periodic_rational& pr) {
      first = pr.first;
      value_init(second);
      value_assign(second, pr.second);
    }
    periodic_rational& operator=(const periodic_rational& pr) {
      if (&pr == this)
	return (*this);
      first = pr.first;
      value_assign(second, pr.second);
      return (*this);
    }
  };
  /** data represents the Ehrhart polynomial. The first part
   * of the pair in the map represents the coefficients of the
   * parameters. The second part represents the coefficient,
   * which is a periodic number.
   * e.g.
   * \f$2x^2y+ 2/3 y^2-z\f$ is represented as
   * ({(x,2),(y,1)},(2,1)), ({(y,2)},(2,3)), ({(z,1)},(-1,1)).
   */
  map<map<string, int>, periodic_rational> data;
  /** simplify removes the terms with coefficient 0
   *  and removes the exponents with power 0.
   */
  void simplify();
public:
  /** Initialization of EhrhartPolynom from a Polylib representation
      of the Ehrhart Polynomial. */
  EhrhartPolynom(const ::evalue* e, const deque<string>& parameter_names);
  /** Initialization to zero. */
  EhrhartPolynom() {
    Value val, one;
    value_init(val);
    value_set_si(val, 0);
    value_init(one);
    value_set_si(one, 1);
    map<string,int> exponent;
    PeriodicNumber pn(val);
    periodic_rational pr(pn, one);
    data[exponent] = pr;
    value_clear(val);
    value_clear(one);
  }
  /** Initialization to constant value. */
  EhrhartPolynom(long cst) {
    Value val, one;
    value_init(val);
    value_set_si(val, cst);
    value_init(one);
    value_set_si(one, 1);
    map<string,int> exponent;
    PeriodicNumber pn(val);
    periodic_rational pr(pn, one);
    data[exponent] = pr;
    value_clear(val);
    value_clear(one);
  }
  /** single-term polynomial. */
  EhrhartPolynom(const map<string,int>& exponent,
		 const periodic_rational& coef) {
    data[exponent] = coef;
    simplify();
  }
  /** single-term polynomial. */
  EhrhartPolynom(const map<string,int>& exponent,
		 const pair<PeriodicNumber, Value>& coef) {
    data[exponent] = periodic_rational(coef.first, coef.second);
    simplify();
  }
  bool operator==(const EhrhartPolynom& ep) const;
  inline bool operator!=(const EhrhartPolynom& ep) const { return !((*this)==ep); }
  bool smaller_than(const EhrhartPolynom& ep) const;
  ~EhrhartPolynom();
  EhrhartPolynom& operator=(const EhrhartPolynom& e1)
  {
    if (this==&e1)
      return (*this);
    /*
    for(map<map<string, int>, periodic_rational >::iterator
	i=data.begin(); i!=data.end(); i++)
      value_clear((*i).second.second);
      */
    data.clear();
    for(map<map<string, int>, periodic_rational>::const_iterator
	i=e1.data.begin(); i!=e1.data.end(); i++) {
      /*
      Value val;
      value_init(val);
      value_assign(val, (*i).second.second);
      pair<map<string, int>, pair<PeriodicNumber, Value> > new_record=(*i);
      new_record.second.second = val;
      data.insert(new_record);
      */
      data[(*i).first] = (*i).second;
    }
    return (*this);
  }
  /** operator+ adds up the this polynomial and polynomial e1. */
  EhrhartPolynom operator+(const EhrhartPolynom& e1) const;
  /** operator+= adds up the polynomial e1 to this polynomial. */
  EhrhartPolynom operator+=(const EhrhartPolynom& e1) {
    /*
    data = ((*this) + e1).data;
    */
    EhrhartPolynom sum = (*this) + e1;
    (*this) = sum;
    return (*this);
  }
  /** operator* multiplies the this polynomial and polynomial e1. */
  EhrhartPolynom operator*(const EhrhartPolynom& e1) const;
  /** operator/ divides by constant. */
  EhrhartPolynom operator/(const Value v) const;
  /** operator*= multiplies this polynomial with polynomial e1. */
  EhrhartPolynom operator*=(const EhrhartPolynom& e1) {
    EhrhartPolynom prod = (*this) * e1;
    (*this) = prod;
    return (*this);
  }
  ::evalue to_evalue(const deque<string>& parameter_names) const;
  double to_double() const;
  string to_string() const;
  inline operator string() const { return to_string(); }
  bool contains_parameters() const;

  /** subst_var substitutes variable @a var in the polynomial
      with the expression encoded in @a subtitution. 
      @param subtitution is a affine function of variables
      occuring in the polynomial, minus the substituted variable
      @a var.
  */
  EhrhartPolynom subst_var(const string& var,
			   const map<string, int>& substitution,
			   const int divisor=1) const;
#if 0
  /*
   * equal_in_domain tries to find out if this polynomial and ep2 are
   * equal in the domain. e.g. (n-i in domain i=1) and (i-1 in domain
   * i=n) are equal (they both evaluate to n-1).
   * @param ep1 the first Ehrhart polynomial.
   * @param ep2 the second Ehrhart polynomial.
   * @param domain1 the domain corresponding to ep1.
   * @param domain2 the domain corresponding to ep2.
   * @return a null pointer if the polynomials are not equal in the
   * domain, otherwise a pointer to a simplified EhrhartPolynomial 
   * is returned. In the example above, n-1 would be returned.
   */
  static EhrhartPolynom* equal_in_domain
  (const EhrhartPolynom& ep1,
   const EhrhartPolynom& ep2,
   const ::Relation& domain1,
   const ::Relation& domain2);
#endif

  /** struct affine represents an Ehrhart Polynomial as a set of
   *  affine functions, if possible. This is only possible if the
   *  polynomial is univariate, with maximum power 1. 
   */
  struct affine {
    /** period indicates the period of the different parameters in
	the Ehrhart polynomial. */
    map<string,int> period;
    /** @c offset[i] contains the constraints under which @c polynomial[i]
	and @c denumerator[i] are legal. e.g. when @c offset[i] equals
	{("a",2),("b",0)} then @c polynomial[i] and @c denumerator[i] are
	legal iff 
	\f[
	a \textrm{ mod \tt\ period[}a\textrm{\tt ] } = 2 \wedge
	b \textrm{ mod \tt\ period[}b\textrm{\tt ] } = 0
	\f]
     */
    deque<map<string,int> > offset;
    /** @c polynomial[i] indicates the polynomial */
    deque<map<string,int> > polynomial;
    /** @c denumerator[i] indicates the denumerator, the polynomial should
	be divided with. */
    deque<int>              denumerator;
    unsigned size() const { 
      assert(offset.size()==polynomial.size());
      assert(denumerator.size()==offset.size());
      return offset.size(); 
    }
    string to_string() const;
  };
  
  /**
   * get_affine_polynomials returns the Ehrhart polynomials as a set
   * of affine functions when possible. If not possible, it returns
   * the null pointer.
   */
  struct affine* get_affine_polynomials() const;

  /**
   * get_AST_representation returns an PolyAstNode representation of
   * the polynom.
   */
#if 0
  /* KB november 2003: only usefull in PolyAST */
  polyast::Expression* get_AST_representation() const;
#endif
private:
  evalue translate_one_term(const deque<string>& parameter_names,
			    const deque<string>& left_over_var_names,
			    const set<map<string, int> >& terms) const;
};


#endif
