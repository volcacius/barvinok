#include "debug.h"
#include "ehrhartpolynom.h"
#include <iostream>
#include <algorithm>
#if 0
#include "polylib_tools.h"
#include "omega_tools.h"
#endif

// only legal if value==GMP!!
#if defined(GNUMP)
#define value_to_str(str, str2, val) { str = mpz_get_str(str2,10,(val)); }
#define value_gcd(ref, val1, val2)     (mpz_gcd((ref),(val1),(val2)))
#define value_lcm(ref, val1, val2)     (mpz_lcm((ref),(val1),(val2)))
#endif

static long gcd(long a, long b)
{
  if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<"gcd("<<a<<", "<<b<<")=";
  }
  long result;
  while (a!=0) {
    result = b%a;
    b=a;
    a=result;
  }
  result = (b<0)?(-b):b;
  if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<result<<endl;
  }
  return result;
}

static long lcm(long a, long b)
{
  long result=(a*b)/gcd(a,b);
  if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<"lcm("<<a<<", "<<b<<")="<<(long)result<<endl;
  }
  return (long) result;
}

string int2string(int i)
{
	char buf[100];
	buf[99]='\0';
	char* cur=&buf[99];
	bool negative=false;
	if (i==0) return "0";
	if (i<0) {
		negative=true;
		i=-i;
	}
	while(i!=0) {
		*--cur = i%10+'0';
		i /= 10;
	}
	if (negative) *--cur = '-';
	return string(cur);
}


char** get_pname(const deque<string>& parameters)
{
  char** result = new char*[parameters.size()];
  for(unsigned i=0;i<parameters.size();i++) {
    result[i] = new char[parameters[i].size()+1];
    strcpy(result[i], parameters[i].c_str());
  }
  return result;
}

void delete_pname(char** pname, const deque<string>& parameters)
{
  for(unsigned i=0;i<parameters.size();i++)
    delete[] pname[i];
  delete[] pname;
}



string value2string(const Value v) {
  char* str=0;
  value_to_str(str, str, v);
  string result=str;
  free(str);
  return result;
}


inline Value& PeriodicNumber::operator[](const map<string, int>& _index) {
  if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=4) {
    cout<<"in PeriodicNumber::operator[]"<<endl;
    cout<<"_index=";
    for(map<string, int>::const_iterator i=_index.begin();
	i!=_index.end(); i++)
      cout<<"("<<(*i).first<<", "<<(*i).second<<") ";
    cout<<endl;
  }
  unsigned long index=0;
  for(map<string,int>::const_iterator i=_index.begin();
      i!=_index.end(); i++) {
    map<string, unsigned long>::const_iterator j=stride.find((*i).first);
    assert(j!=stride.end());
    index += (*j).second * (*i).second;
  }
  assert(index>=0);
  assert(index<datasize);
  if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=4) {
    cout<<"index="<<index<<endl;
    cout<<"data[index]="<<value2string(data[index])<<endl;
  }
  return data[index];
}

inline const Value& PeriodicNumber::operator[]
(const map<string, int>& _index) const
{
  if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=4) {
    cout<<"in PeriodicNumber::operator[] const"<<endl;
    cout<<"_index=";
    for(map<string, int>::const_iterator i=_index.begin();
	i!=_index.end(); i++)
      cout<<"("<<(*i).first<<", "<<(*i).second<<") ";
    cout<<endl;
  }
  unsigned long index=0;
  for(map<string,int>::const_iterator i=_index.begin();
      i!=_index.end(); i++) {
    map<string, unsigned long>::const_iterator j=stride.find((*i).first);
    assert(j!=stride.end());
    index += (*j).second * (*i).second;
  }
  assert(index>=0);
  assert(index<datasize);
  if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=4) {
    cout<<"index="<<index<<endl;
    cout<<"data[index]="<<value2string(data[index])<<endl;
  }
  return data[index];
}
/*
static long gcd(long a, long b)
{
  if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<"gcd("<<a<<", "<<b<<")=";
  }
  long result;
  while (a!=0) {
    result = b%a;
    b=a;
    a=result;
  }
  result = (b<0)?(-b):b;
  if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<result<<endl;
  }
  return result;
}
*/
/*
static long lcm(long a, long b)
{
  long result=(a*b)/gcd(a,b);
  if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<"lcm("<<a<<", "<<b<<")="<<(long)result<<endl;
  }
  return (long) result;
}
*/

/** operator+ adds up two periodic numbers \f$pn_1\f$ and
 *  \f$pn_2\f$. The following algorithm is used: The parameters in the
 *  result equals the union of the parameters in the two periodic
 *  numbers to be counted together.  The periods of any parameter
 *  \f$p\f$ is the least common multiple of the period of the
 *  parameter in \f$pn_1\f$ and \f$pn_2\f$. If the parameter doesn't
 *  occur in one of \f$pn_1\f$ and \f$pn_2\f$, its period is taken is
 *  being one in \f$pn_1\f$ or \f$pn_2\f$.  After this step, the
 *  values in the matrix are filled in, being the sum of the 
 *  corresponding values in \f$pn_1\f$ and \f$pn_2\f$.*/
PeriodicNumber PeriodicNumber::operator+(const PeriodicNumber& pn2) const
{
  const PeriodicNumber& pn1=(*this);
  if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<"In PeriodicNumber::operator+:"<<endl;
    cout<<"pn1="<<pn1.to_string()<<endl;
    cout<<"pn2="<<pn2.to_string()<<endl;
  }
  set<string> parameters=pn1.parameter_names;
  map<string, int> period;
  for(set<string>::const_iterator i=pn2.parameter_names.begin();
      i!=pn2.parameter_names.end(); i++)
    parameters.insert(*i);

  for(set<string>::const_iterator i=parameters.begin();
      i!=parameters.end(); i++) {
    bool param_in_pn1 = 
      (pn1.parameter_names.find(*i)!=pn1.parameter_names.end());
    bool param_in_pn2 = 
      (pn2.parameter_names.find(*i)!=pn2.parameter_names.end());
    if (param_in_pn1 && param_in_pn2) {
      map<string, int>::const_iterator j=pn1.period.find(*i);
      assert(j!=pn1.period.end());
      map<string, int>::const_iterator k=pn2.period.find(*i);
      assert(k!=pn2.period.end());
      period[*i] = lcm((*j).second, (*k).second);
    } else if (param_in_pn1) {
      map<string, int>::const_iterator j=pn1.period.find(*i);
      assert(j!=pn1.period.end());
      period[*i] = (*j).second;
    } else if (param_in_pn2) {
      map<string, int>::const_iterator j=pn2.period.find(*i);
      assert(j!=pn2.period.end());
      period[*i] = (*j).second;
    } else assert(0);
  }

  PeriodicNumber result(parameters, period);

  deque<unsigned long> strides2;
  for(map<string, unsigned long>::const_iterator i=result.stride.begin();
      i!=result.stride.end(); i++)
    strides2.push_back((*i).second);
  strides2.push_back(result.datasize);
  // now fill in the correct values.
  Value tmp;
  value_init(tmp);
  for(unsigned long i=0;i<result.datasize;i++) {
    // find out which period-values for which parameters
    // map to data[i].
    map<string, int> index;
    int param_nr=0;
    for(set<string>::const_iterator j=parameters.begin();
	j!=parameters.end(); j++, param_nr++) {
      index[*j] = (i%(strides2[param_nr+1])) / strides2[param_nr];
    }
    value_set_si(tmp, 0);
    map<string, int> pn1_index;
    map<string, int> pn2_index;
    for(map<string, int>::const_iterator j=index.begin();j!=index.end(); j++) {
      if (pn1.parameter_names.find((*j).first)==pn1.parameter_names.end())
	continue;
      map<string, int>::const_iterator periodd=pn1.period.find((*j).first);
      assert(periodd!=pn1.period.end());
      pn1_index[(*j).first] = (*j).second % (*periodd).second;
    }
    for(map<string, int>::const_iterator j=index.begin();j!=index.end(); j++) {
      if (pn2.parameter_names.find((*j).first)==pn2.parameter_names.end())
	continue;
      map<string, int>::const_iterator periodd=pn2.period.find((*j).first);
      assert(periodd!=pn2.period.end());
      pn2_index[(*j).first] = (*j).second % (*periodd).second;
    }

    value_addto(tmp, tmp, pn1[pn1_index]);
    value_addto(tmp, tmp, pn2[pn2_index]);
    value_assign(result[index] , tmp);
  }
  value_clear(tmp);

  // if the periodic number has parameters (i.e. not constant coef),
  // then "" should not occur as parameter!
  if (result.parameter_names.find("") != result.parameter_names.end()) {
    assert(result.period[""] == 1);
    result.parameter_names.erase("");
    result.period.erase("");
    result.stride.erase("");
  } 

  if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<"In PeriodicNumber::operator+ result of addition is:"<<endl;
    cout<<result.to_string()<<endl;
  }

  return result;
}

PeriodicNumber PeriodicNumber::operator*(const ::Value v) const
{
  if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<"In PeriodicNumber::operator*(value v):"<<endl;
    cout<<"v="<<value2string(v)<<endl;
    cout<<"pn="<<to_string()<<endl;
  }
  PeriodicNumber result(*this);
  if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<"Before multiplication, result="<<result.to_string()<<endl;
  }
  for(unsigned i=0;i<result.datasize;i++)
    value_multiply(result.data[i],result.data[i],v);
  if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<"After multiplication with v="<<value2string(v)<<":"<<endl;
    cout<<"result="<<result.to_string()<<endl;
  }
  return result;
}

void PeriodicNumber::get_gcd(Value& v) const
{
  if (datasize==0) {
    value_set_si(v, 1);
    return;
  }
  value_assign(v,data[0]);
  for(unsigned i=1;i<datasize;i++)
    value_gcd(v, v, data[i]);
}

void PeriodicNumber::divide_by(const Value v)
{
#ifndef NDEBUG
  Value gcd;
  Value gcd2;
  value_init(gcd); value_init(gcd2);
  get_gcd(gcd); get_gcd(gcd2);
  value_division(gcd, gcd, v);
  value_multiply(gcd, gcd, v);
  if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<"In PeriodicNumber::divide_by("<<value2string(v)<<")"<<endl;
    cout<<" periodicnumber = "<<to_string()<<endl;
    cout<<" gcd = "<<value2string(gcd)<<"; gcd2 = "<<value2string(gcd2)<<endl;
  }
  assert(value_eq(gcd,gcd2));
  value_clear(gcd);
  value_clear(gcd2);
#endif
  for(unsigned i=0; i<datasize; i++)
    value_division(data[i] , data[i], v);
}

/** operator* multiplies two periodic numbers \f$pn_1\f$ and
 *  \f$pn_2\f$. The following algorithm is used: The parameters in the
 *  result equals the union of the parameters in the two periodic
 *  numbers to be multiplied.  The periods of any parameter
 *  \f$p\f$ is the least common multiple of the period of the
 *  parameter in \f$pn_1\f$ and \f$pn_2\f$. If the parameter doesn't
 *  occur in one of \f$pn_1\f$ and \f$pn_2\f$, its period is taken is
 *  being one in \f$pn_1\f$ or \f$pn_2\f$.  After this step, the
 *  values in the matrix are filled in, being the product of the
 *  corresponding values in \f$pn_1\f$ and \f$pn_2\f$.*/
PeriodicNumber PeriodicNumber::operator*(const PeriodicNumber& pn2) const
{
  const PeriodicNumber& pn1=(*this);
  if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<"In PeriodicNumber::operator*:"<<endl;
    cout<<"pn1="<<pn1.to_string()<<endl;
    cout<<"pn2="<<pn2.to_string()<<endl;
  }
  set<string> parameters=pn1.parameter_names;
  map<string, int> period;
  for(set<string>::const_iterator i=pn2.parameter_names.begin();
      i!=pn2.parameter_names.end(); i++)
    parameters.insert(*i);

  for(set<string>::const_iterator i=parameters.begin();
      i!=parameters.end(); i++) {
    bool param_in_pn1 = 
      (pn1.parameter_names.find(*i)!=pn1.parameter_names.end());
    bool param_in_pn2 = 
      (pn2.parameter_names.find(*i)!=pn2.parameter_names.end());
    if (param_in_pn1 && param_in_pn2) {
      map<string, int>::const_iterator j=pn1.period.find(*i);
      assert(j!=pn1.period.end());
      map<string, int>::const_iterator k=pn2.period.find(*i);
      assert(k!=pn2.period.end());
      period[*i] = lcm((*j).second, (*k).second);
    } else if (param_in_pn1) {
      map<string, int>::const_iterator j=pn1.period.find(*i);
      assert(j!=pn1.period.end());
      period[*i] = (*j).second;
    } else if (param_in_pn2) {
      map<string, int>::const_iterator j=pn2.period.find(*i);
      assert(j!=pn2.period.end());
      period[*i] = (*j).second;
    } else assert(0);
  }

  PeriodicNumber result(parameters, period);
  if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<"In PeriodicNumber::operator*:"<<endl;
    cout<<"result.parameters=";
    for(set<string>::const_iterator i=parameters.begin();
        i!=parameters.end(); i++)
      cout<<(*i)<<" ";
    cout<<endl;
    cout<<"result.period=";
    for(map<string, int>::const_iterator i=period.begin();
        i!=period.end(); i++)
      cout<<"("<<(*i).first<<","<<(*i).second<<") ";
    cout<<endl;
  }

  deque<unsigned long> strides2;
  for(map<string, unsigned long>::const_iterator i=result.stride.begin();
      i!=result.stride.end(); i++)
    strides2.push_back((*i).second);
  strides2.push_back(result.datasize);
  // now fill in the correct values.
  Value tmp;
  value_init(tmp);
  for(unsigned long i=0;i<result.datasize;i++) {
    // find out which period-values for which parameters
    // map to data[i].
    map<string, int> index;
    int param_nr=0;
    for(set<string>::const_iterator j=parameters.begin();
	j!=parameters.end(); j++, param_nr++) {
      index[*j] = (i%(strides2[param_nr+1])) / strides2[param_nr];
    }
    value_set_si(tmp, 1);
    map<string, int> pn1_index;
    map<string, int> pn2_index;
    for(map<string, int>::const_iterator j=index.begin();j!=index.end(); j++) {
      if (pn1.parameter_names.find((*j).first)==pn1.parameter_names.end())
	continue;
      map<string, int>::const_iterator periodd=pn1.period.find((*j).first);
      assert(periodd!=pn1.period.end());
      pn1_index[(*j).first] = (*j).second % (*periodd).second;
    }
    for(map<string, int>::const_iterator j=index.begin();j!=index.end(); j++) {
      if (pn2.parameter_names.find((*j).first)==pn2.parameter_names.end())
	continue;
      map<string, int>::const_iterator periodd=pn2.period.find((*j).first);
      assert(periodd!=pn2.period.end());
      pn2_index[(*j).first] = (*j).second % (*periodd).second;
    }

    if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=4) {
      cout<<"In PeriodicNumber::operator*, i="<<i<<", result.datasize="<<datasize<<"\n pn1_index=";
      for(map<string, int>::const_iterator z=pn1_index.begin(); z!=pn1_index.end(); z++)
	cout<<"("<<(*z).first<<","<<(*z).second<<") ";
      cout<<"\npn2_index=";
      for(map<string, int>::const_iterator z=pn2_index.begin(); z!=pn2_index.end(); z++)
	cout<<"("<<(*z).first<<","<<(*z).second<<") ";
      cout<<endl;
      cout<<"pn1[pn1_index]="<<value2string(pn1[pn1_index])<<endl;
      cout<<"pn2[pn2_index]="<<value2string(pn2[pn2_index])<<endl;
      cout<<"tmp="<<value2string(tmp)<<endl;
    }
    value_multiply(tmp, tmp, pn1[pn1_index]);
    if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=4) { cout<<"tmp="<<value2string(tmp)<<endl; }
    value_multiply(tmp, tmp, pn2[pn2_index]);
    if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=4) { cout<<"tmp="<<value2string(tmp)<<endl; }
    value_assign(result[index] , tmp);
  }
  value_clear(tmp);

  // if the periodic number has parameters (i.e. not constant coef),
  // then "" should not occur as parameter!
  if (result.parameter_names.find("") != result.parameter_names.end()) {
    assert(result.period[""] == 1);
    result.parameter_names.erase("");
    result.period.erase("");
    result.stride.erase("");
  } 

  if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<"In PeriodicNumber::operator* result of multiplication is:"<<endl;
    cout<<result.to_string()<<endl;
    cout<<"result.parameter_names=";
    for(set<string>::const_iterator i=result.parameter_names.begin();
        i!=result.parameter_names.end(); i++)
      cout<<(*i)<<" ";
    cout<<endl;
    cout<<"result.period=";
    for(map<string, int>::const_iterator i=result.period.begin();
        i!=result.period.end(); i++)
      cout<<"("<<(*i).first<<","<<(*i).second<<") ";
    cout<<"result.stride =";
    for(map<string, unsigned long>::const_iterator i=result.stride.begin();
        i!=result.stride.end(); i++)
      cout<<"("<<(*i).first<<","<<(*i).second<<") ";
    cout<<"result.data=";
    for(unsigned i=0;i<result.datasize;i++) {
      char* str=0;
      value_to_str(str,0, result.data[i]);
      cout<<str<<" ";
      free(str);
    }
    cout<<endl;
  }
 
  return result;
}

double PeriodicNumber::to_double() const
{
  assert(datasize==1);
  return VALUE_TO_DOUBLE(data[0]);
}

bool PeriodicNumber::smaller_than(const PeriodicNumber& pn) const
{
  if (datasize!=1 || pn.datasize!=1)
    return false;
  return value_lt(data[0], pn.data[0]);
}

string PeriodicNumber::to_string() const
{
  string result;
  deque<unsigned long> strides2;
  for(map<string, unsigned long>::const_iterator i=stride.begin();
      i!=stride.end(); i++)
    strides2.push_back((*i).second);
  strides2.push_back(datasize);
  // now fill in the correct values.
  for(unsigned long i=0;i<datasize;i++) {
    // find out which period-values for which parameters
    // map to data[i].
    map<string, int> index;
    int param_nr=0;
    for(set<string>::const_iterator j=parameter_names.begin();
	j!=parameter_names.end(); j++, param_nr++) {
      index[*j] = (i%(strides2[param_nr+1])) / strides2[param_nr];
    }
    for(map<string, int>::const_iterator j=index.begin();
	j!=index.end(); j++) {
      map<string,int>::const_iterator k=period.find((*j).first);
      assert(k!=period.end());
      if ((*j).second==0 && (*k).second!=1)
	result+="[ ";
      else
	break;
    }
    Value val;
    value_init(val);
    value_assign(val,((*this)[index]));
    char* valstr=0;
    value_to_str(valstr,valstr,val);
    result += string(valstr)+" ";
    free(valstr);
    for(map<string, int>::const_iterator j=index.begin();
	j!=index.end(); j++) {
      map<string,int>::const_iterator k=period.find((*j).first);
      assert(k!=period.end());
      if ((*j).second==(*k).second-1 && (*k).second!=1) {
	result+="]";
	if ((*j).first!=string(""))
	  result+="_";
	result+=(*j).first+", ";
      }
      else
	break;
    }
  }
  if (result.size()>=2 && result.substr(result.size()-2, result.size())
      ==", ")
    result.erase(result.size()-2, result.size());

  if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    string debug_str;
    debug_str+="(\n\tparameter_names=";
    for(set<string>::const_iterator i=parameter_names.begin();
        i!=parameter_names.end(); i++)
      debug_str+=(*i)+" ";
    debug_str+="\n\tperiod=";
    for(map<string, int>::const_iterator i=period.begin();
        i!=period.end(); i++)
      debug_str+="("+(*i).first+","+int2string((*i).second)+") ";
    debug_str+="\n\tstride=";
    for(map<string, unsigned long>::const_iterator i=stride.begin();
        i!=stride.end(); i++)
      debug_str+="("+(*i).first+","+int2string((*i).second)+") ";
    debug_str+="\n\tdata=";
    for(unsigned i=0;i<datasize;i++) {
      char* str=0;
      value_to_str(str,0, data[i]);
      debug_str+=str+string(" ");
      free(str);
    }
    debug_str+=")\n";
    result += debug_str;
  }

  return result;
}

static int find_pos_of(const string& name,
		       const deque<string>& names)
{
  for(unsigned i=0;i<names.size(); i++)
    if (name==names[i])
      return i;
  return -1;
}

bool PeriodicNumber::has_param(const string& name) const 
{
  if (DebugBlackboard.debug("AFFINEPOLYNOM")>=6) {
    cout<<"in PeriodicNumber::has_param("<<name<<") (this="
	<<to_string()<<"): return "
	<<(period.find(name)!=period.end())<<endl;
  }
  return (period.find(name)!=period.end());
}

static void   pn_to_evalue(const PeriodicNumber& pn,
			   const map<string,int>& bound_vars,
			   const set<string>& free_vars,
			   const deque<string>& ep_parameter_names,
			   evalue& ev)
{
  if (free_vars.size()==0) {
    value_set_si(ev.d,1);
    value_assign(ev.x.n , pn[bound_vars]);
    return;
  } else {
    string name;
    // find an appropriate var_name, i.e., different from ""
    // if possible
    for(set<string>::const_iterator i=free_vars.begin();
	i!=free_vars.end(); i++)
      if (name=="")
	name=(*i);
    if (name=="") {
      assert(free_vars.size()==1);
      //      assert(bound_vars.size()==0);
      assert(bound_vars.find("")==bound_vars.end());
      assert(pn.get_period("")==1);
      map<string,int> index=bound_vars;
      index[""]=0;
      value_set_si(ev.d,1);
      value_assign(ev.x.n , pn[index]);
      return;
    }
    set<string> new_free_vars=free_vars;
    new_free_vars.erase(new_free_vars.find(name));
    value_set_si(ev.d,0);
    value_clear(ev.x.n);
    ev.x.p = new_enode(periodic, pn.get_period(name), 
                        find_pos_of(name, ep_parameter_names)+1);
    for(int i=0;i<pn.get_period(name);i++) {
      map<string,int> new_bound_vars=bound_vars;
      assert(new_bound_vars.find(name)==new_bound_vars.end());
      new_bound_vars[name]=i;
      pn_to_evalue(pn, new_bound_vars, new_free_vars,
                   ep_parameter_names, ev.x.p->arr[i]);
    }
  }
}

::evalue PeriodicNumber::to_evalue(const deque<string>& ep_parameter_names)
  const
{
  map<string,int> bound_vars;
  set<string> free_vars=this->parameter_names;

  ::evalue result;
  value_init(result.d);
  value_init(result.x.n);
  pn_to_evalue(*this, bound_vars, free_vars, ep_parameter_names, result);
  return result;
}


static bool check_for_polynomials(::enode* e)
{
  assert(e!=0);
  if (e->type==polynomial)
    return true;
  else
    for(int i=0;i<e->size;i++)
      if (value_zero_p(e->arr[i].d) &&
	  check_for_polynomials(e->arr[i].x.p))
	return true;
  return false;
}

static void find_names(::enode* e, 
		       const deque<string>& parameter_names,
		       set<string>& names)
{
  assert(e!=0);
  assert(e->type==periodic);

  names.insert(parameter_names[e->pos-1]);
  for(int i=0;i<e->size;i++)
    if (value_zero_p(e->arr[i].d))
      find_names(e->arr[i].x.p, parameter_names, names);
}

static void find_periods(::enode* e,
			 const deque<string>& parameter_names,
			 map<string, set<int> >& periods)
{
  assert(e!=0);
  assert(e->type==periodic);

  periods[parameter_names[e->pos-1]].insert(e->size);
  for(int i=0;i<e->size;i++)
    if (value_zero_p(e->arr[i].d))
      find_periods(e->arr[i].x.p, parameter_names, periods);
}

static void find_denominators(::enode* e,
			      set<Value> &denominators)
{
  assert(e!=0);
  assert(e->type==periodic);

  for(int i=0;i<e->size;i++) {
    if (value_zero_p(e->arr[i].d))
      find_denominators(e->arr[i].x.p, denominators);
    else
      denominators.insert(e->arr[i].d);
  }
}

#if 0

static void fill_in_values(PeriodicNumber& pn,
			   ::enode* e,
			   const deque<string>& parameter_names,
			   map<string, int>& name2lcm_period,
			   map<string, set<int> >& name2indices,
			   const int lcm_denom)
{
  assert(e!=0);
  assert(e->type==periodic);
  static int depth=0;
  depth++;
  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<"Entering fill_in_values (depth "<<depth<<"):"<<endl;
    cout<<"e=";
    char **pname=get_pname(parameter_names);
    ::print_enode(stdout, e, pname);
    delete_pname(pname, parameter_names);
    cout<<endl;    
    cout<<"parameter_names= ";
    for(deque<string>::const_iterator i=parameter_names.begin();
	i!=parameter_names.end(); i++)
      cout<<(*i)<<" ";
    cout<<endl;
    cout<<"name2lcm_period= ";
    for(map<string,int>::const_iterator i=name2lcm_period.begin();
	i!=name2lcm_period.end(); i++)
      cout<<"("<<(*i).first<<","<<(*i).second<<")"<<" ";
    cout<<endl;
    cout<<"name2indices= ";
    for(map<string,set<int> >::const_iterator i=name2indices.begin();
	i!=name2indices.end(); i++) {
      cout<<"("<<(*i).first<<",{";
      for(set<int>::const_iterator j=(*i).second.begin();
	  j!=(*i).second.end(); j++) {
	if (j!=(*i).second.begin())
	  cout<<" ";
	cout<<(*j);
      }
      cout<<"}) ";
    }
    cout<<endl;
    cout<<"lcm_denom= "<<lcm_denom<<endl;
  }

  string name=parameter_names[e->pos-1];
  int period =e->size;
  int nr_periods = name2lcm_period[name] / e->size;
  assert(nr_periods*e->size == name2lcm_period[name]);

  assert(name2indices.find(name) == name2indices.end());
  for(int i=0;i<e->size;i++) {
    for(int k=0;k<nr_periods;k++)
      name2indices[name].insert(i+k*period);
    if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
      cout<<" i="<<i<<"; period = "<<period
	  <<"; nr_periods = "<<nr_periods<<"name2indices = ";
      for(map<string,set<int> >::const_iterator i=name2indices.begin();
	  i!=name2indices.end(); i++) {
	cout<<"("<<(*i).first<<",{";
	for(set<int>::const_iterator j=(*i).second.begin();
	    j!=(*i).second.end(); j++) {
	  if (j!=(*i).second.begin())
	    cout<<" ";
	  cout<<(*j);
	}
	cout<<"}) ";
      }
      cout<<endl;
    }

    if (e->arr[i].d==0) {
      fill_in_values(pn, e->arr[i].x.p, parameter_names,
		     name2lcm_period, name2indices, lcm_denom);
    } else {
      // e->arr[i].d!=0
      int multiplier = lcm_denom/e->arr[i].d;
      if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
	cout<<"multiplier="<<multiplier<<endl;
      }
      assert(multiplier * e->arr[i].d == lcm_denom);
      deque<unsigned long> strides;
      map<string, deque<int> > name2indices2;
      unsigned long total=1;
      for(map<string, set<int> >::const_iterator j=name2indices.begin();
	  j!=name2indices.end(); j++) {
	strides.push_back(total);
	total *= (*j).second.size();
	for(set<int>::const_iterator k=(*j).second.begin();
	    k!=(*j).second.end(); k++)
	  name2indices2[(*j).first].push_back(*k);
      }
      strides.push_back(total);
      for(unsigned long j=0;j<total; j++) {
	map<string,int> index;
	unsigned long t=0;
	for(map<string, set<int> >::const_iterator k=name2indices.begin();
	    k!=name2indices.end(); k++,t++) {
	  index[(*k).first] = 
	    name2indices2[(*k).first][(j%strides[t])/strides[t]];
	}
	pn[index] = multiplier * e->arr[i].x.n;
      }
    }
    if (name2indices.find(name) != name2indices.end())
      name2indices.erase(name2indices.find(name));
  }
  depth--;
}


/** handle_periodic_enode processes an enode
 *  of type periodic
 *  @return a pair. The first part is the periodic number
 *          The second part is the denominator of the
 *	    periodic number.
 */
static pair<PeriodicNumber,int> handle_periodic_enode
(const deque<string>& parameter_names,
 ::enode* e)
{
  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<"Entering handle_periodic_enode"<<endl;
    cout<<"e= ";
    char **pname=get_pname(parameter_names);
    ::print_enode(stdout, e, pname);
    delete_pname(pname, parameter_names);
    cout<<endl;
  }
  // check some assumptions first:
  // 1. there are no polynomial enodes further down.
  assert(false==check_for_polynomials(e));
  // 2. find the parameter names occuring in the Ehrhart Polynomial
  set<string> names;
  find_names(e, parameter_names, names);
  // 3. find all the periods for all the names that occur.
  map<string, set<int> > name2periods;
  find_periods(e, parameter_names, name2periods);
  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<"In handle_periodic_enode"<<endl;
    cout<<"name2periods = ";
    for(map<string, set<int> >::const_iterator i=name2periods.begin();
	i!=name2periods.end(); i++) {
      cout<<"("<<(*i).first<<", {";
      for(set<int>::const_iterator j=(*i).second.begin();
	  j!=(*i).second.end(); j++) {
	if (j!=(*i).second.begin())
	  cout<<",";
	cout<<(*j);
      }
      cout<<"}) ";
    }
    cout<<endl;
  }
  // 4. find lcm for all periods of a given name.
  map<string, int> name2lcm_period;
  for(set<string>::const_iterator i=names.begin();
      i!=names.end(); i++) {
    assert(name2periods.find(*i)!=name2periods.end());
    int lcmp=1;
    for(set<int>::const_iterator j=name2periods[*i].begin();
	j!=name2periods[*i].end(); j++)
      lcmp = lcm(lcmp, *j);
    name2lcm_period[*i] = lcmp;
  }
  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<"In handle_periodic_enode"<<endl;
    cout<<"name2lcm_period = ";
    for(map<string,int>::const_iterator i=name2lcm_period.begin();
	i!=name2lcm_period.end(); i++) {
      cout<<"("<<(*i).first<<","<<(*i).second<<") ";
    }
    cout<<endl;
  }
  // 5. find lcm of all denominators
  int lcmd=1;
  set<Value> denominators;
  find_denominators(e, denominators);
  for(set<Value>::const_iterator j=denominators.begin();
      j!=denominators.end(); j++)
    lcmd = lcm(lcmd, *j);
  int common_denom=lcmd;
  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<"common_denom = "<<common_denom<<endl;;
  }
  
  // 5. Construct PeriodicNumber, phase 1
  PeriodicNumber result(names, name2lcm_period);
  // 6. Fill in the values for the PeriodicNumber

  // name2indices indicates for which indices of which
  // parameter names the given value must be filled in.
  map<string, set<int> > name2indices;
  fill_in_values(result, e, parameter_names, name2lcm_period, name2indices,
		 common_denom);
  
  return pair<PeriodicNumber,int>(result,common_denom);
}
					    

static void process_sub_part_of_evalue
(map<map<string, int>, pair<PeriodicNumber, Value> >& data,
 const deque<string>& parameter_names,
 const map<string,int>& current_exponents,
 const PeriodicNumber& current_periodic,
 Value current_denumerator,
 ::evalue* e)
{
  static int depth=0; //only used for debugging
  depth++;

  assert(e!=0);

  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=2) {
    cout<<"Entering process_sub_part_of_evalue (depth "<<depth<<"):"<<endl;
    cout<<"parameter_names= ";
    for(deque<string>::const_iterator i=parameter_names.begin();
	i!=parameter_names.end(); i++)
      cout<<(*i)<<" ";
    cout<<endl;
    cout<<"current_exponent= ";
    for(map<string,int>::const_iterator i=current_exponents.begin();
	i!=current_exponents.end(); i++)
      cout<<"("<<(*i).first<<","<<(*i).second<<")"<<" ";
    cout<<endl;
    cout<<"current_periodic= "<<current_periodic.to_string()<<endl;
    cout<<"current_denumerator= "<<current_denumerator<<endl;
    cout<<"e->d="<<e->d<<endl;
  }

  if (value_notzero_p(e->d)) {
    // the evalue e represents a rational number
    assert(data.find(current_exponents)==data.end());
    value_direct_product(current_denumerator, e->d);
    PeriodicNumber pn;
    if (current_periodic.is_zero())
      pn=PeriodicNumber(e->x.n);
    else
      pn = current_periodic * e->x.n;
    data[current_exponents] = 
      pair<PeriodicNumber, Value>(pn, current_denumerator);
  } else {
    ::enode* en=e->x.p;
    switch(en->type) {
    case polynomial: 
      {
	if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
	  cout<<"The enode is a polynomial enode."<<endl;
	}
	string param_name;
	if (en->pos>=1)
	  param_name = parameter_names[en->pos-1];
	else
	  param_name = "";
	assert(current_exponents.find(param_name)==
	       current_exponents.end());
	const int degree=en->size-1;
	for(int i=0;i<=degree;i++) {
	  map<string,int> exp=current_exponents;
	  exp[param_name] = i;
	  process_sub_part_of_evalue(data, parameter_names,
				     exp, current_periodic,
				     current_denumerator,
				     &(en->arr[i]));
	}
      }
      break;
    case periodic: 
      {
	if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
	  cout<<"The enode is a periodice enode."<<endl;
	}
	pair<PeriodicNumber,int> pn=handle_periodic_enode(parameter_names, en);
	data[current_exponents] =
	  pair<PeriodicNumber, Value>(pn.first, current_denumerator*pn.second);
      }
      break;
    case evector:
      assert(0);
      break;
    }
  }
  depth--;
}

#endif

static EhrhartPolynom convert_evalue(const ::evalue e, 
				     const deque<string>& parameter_names)
{
  if (value_notzero_p(e.d)) {
    // the evalue e represents a rational number
    PeriodicNumber pn(e.x.n);
    map<string,int> exponent;
    exponent[""]=1;
    return EhrhartPolynom(exponent, 
	                  EhrhartPolynom::periodic_rational(pn,e.d));
  } else {
    EhrhartPolynom result;
    ::enode* en=e.x.p;
    switch(en->type) {
    case polynomial:
      {
      	if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
	  cout<<"The enode is a polynomial enode."<<endl;
	}
	string param_name;
	if (en->pos>=1)
	  param_name = parameter_names[en->pos-1];
	else
	  param_name = "";
	const int degree=en->size-1;
	for(int i=0;i<=degree;i++) {
	  Value one;
	  value_init(one);
	  value_set_si(one, 1);
	  PeriodicNumber pn(one);
	  map<string,int> exponent;
	  exponent[param_name]=i;
	  EhrhartPolynom ep1(exponent, EhrhartPolynom::periodic_rational(pn,one));
	  result += ep1 * convert_evalue(en->arr[i],parameter_names); 
	}
      }
      break;
    case periodic:
      {
	if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
	  cout<<"The enode is a periodice enode."<<endl;
	}
	string name=parameter_names[en->pos-1];
	assert(name!=string(""));
	int period=en->size;
	set<string> param_name_set;
	param_name_set.insert(name);
	map<string,int> period_map;
	period_map[name]=period;
	map<string,int> exponent;
	exponent[""]=1;
	for(int i=0;i<period;i++) {
	  // construct PeriodicNumber with parameter name and
	  // period period, and 0 everywhere, except in location i.
	  PeriodicNumber pn(param_name_set, period_map);
	  map<string,int> index;
	  index[name]=i;
	  Value one1, one2;
	  value_init(one1);
	  value_init(one2);
	  value_set_si(one1, 1);
	  value_set_si(one2, 1);
	  value_assign(pn[index],one1);
	  EhrhartPolynom ep(exponent, EhrhartPolynom::periodic_rational(pn, one2));
	  result += ep * convert_evalue(en->arr[i], parameter_names);
	}
      }
      break;
    case evector:
      assert(0);
      break;
    }
    return result;
  }
}

EhrhartPolynom::EhrhartPolynom(const ::evalue* e, 
			       const deque<string>& parameter_names) 
{
  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=2) {
    cout<<"In EhrhartPolynom(e): e="<<endl;
    char** pname=get_pname(parameter_names);
    print_evalue(stdout, const_cast< ::evalue*>(e), pname);
    delete_pname(pname, parameter_names);
    cout<<endl;
  }
#if 0
  map<string,int> current_exponents;
  PeriodicNumber t;
  process_sub_part_of_evalue(data, parameter_names,  
			     current_exponents, t, 1, e);
#endif
  assert(e!=0);
  data = convert_evalue(*e, parameter_names).data;
  simplify();
  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=2) {
    cout<<"In EhrhartPolynom(e): constructed following EhrhartPolynom:"<<endl;
    cout<<to_string()<<endl;
  }
}

EhrhartPolynom::~EhrhartPolynom() {
  data.clear();
  /*
  for(map<map<string, int>, periodic_rational>::iterator i=
      data.begin(); i!=data.end(); i++)
    value_clear((*i).second.second);
    */
}

EhrhartPolynom 
EhrhartPolynom::operator+(const EhrhartPolynom& e2) const
{
  const EhrhartPolynom& e1=(*this);
  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=2) {
    cout<<"In EhrhartPolynom()::operator+: adding following polynoms"<<endl;
    cout<<e1.to_string()<<endl;
    cout<<e2.to_string()<<endl;
  }
  EhrhartPolynom result;
  for(map<map<string, int>, periodic_rational >::const_iterator i=
	e1.data.begin(); i!=e1.data.end(); i++) {
    map<map<string, int>, periodic_rational >::const_iterator j =
      e2.data.find((*i).first);
    if (j==e2.data.end()) {
      // the term was not found in e2.
      result.data[(*i).first] = (*i).second;
    } else {
      Value lcm_val;
      value_init(lcm_val);
      value_lcm(lcm_val, (*i).second.second, (*j).second.second);
      //unsigned long lcm_val = lcm( (*i).second.second, (*j).second.second);
      Value prod1, prod2;
      value_init(prod1); value_init(prod2);
      value_division(prod1, lcm_val, (*i).second.second);
      value_division(prod2, lcm_val, (*j).second.second);
      PeriodicNumber pn1 = (*i).second.first * prod1;
      PeriodicNumber pn2 = (*j).second.first * prod2;
      value_clear(prod1);
      value_clear(prod2);
      result.data[(*i).first] = periodic_rational(pn1+pn2, lcm_val);
    }
  }
  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<"In EhrhartPolynom()::operator+: result of addition (after adding terms to term1)"<<endl;
    cout<<result.to_string()<<endl;
  }
  for(map<map<string, int>, periodic_rational >::const_iterator i=
	e2.data.begin(); i!=e2.data.end(); i++) {
    map<map<string, int>, periodic_rational >::const_iterator j =
      e1.data.find((*i).first);
    if (j==e1.data.end()) {
      result.data[(*i).first] = (*i).second;
    }
  }
  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<"In EhrhartPolynom()::operator+: result of addition (before simplify)"<<endl;
    cout<<result.to_string()<<endl;
  }
  result.simplify();
  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=2) {
    cout<<"In EhrhartPolynom()::operator+: result of addition"<<endl;
    cout<<result.to_string()<<endl;
  }
  return result;
}

EhrhartPolynom EhrhartPolynom::operator*(const EhrhartPolynom& e2) const
{
  const EhrhartPolynom& e1=(*this);
  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=2) {
    cout<<"In EhrhartPolynom()::operator*: multiplying"
      " following polynoms"<<endl;
    cout<<e1.to_string()<<endl;
    cout<<e2.to_string()<<endl;
  }
  EhrhartPolynom result;
  Value div1, div2;
  value_init(div1); value_init(div2);
  Value lcm_val;
  value_init(lcm_val);
  Value gcd_product_denominator;
  value_init(gcd_product_denominator);
  Value gcd_of_product;
  value_init(gcd_of_product);
  for(map<map<string, int>, periodic_rational>::const_iterator i=
	e1.data.begin(); i!=e1.data.end(); i++) {
    for(map<map<string, int>, periodic_rational>::const_iterator j =
	  e2.data.begin(); j!=e2.data.end(); j++) {
      // add up exponents
      map<string,int> exponent=(*i).first;
      for(map<string,int>::const_iterator k=(*j).first.begin();
	  k!=(*j).first.end(); k++)
	exponent[(*k).first]+=(*k).second;
      // multiply coefficients
      value_lcm(lcm_val, (*i).second.second, (*j).second.second);
      Value common_denominator;
      value_init(common_denominator);
      value_multiply(common_denominator, lcm_val, lcm_val);
      /*
      unsigned long lcm_val = lcm( (*i).second.second, (*j).second.second);
      unsigned long common_denominator = lcm_val*lcm_val;
      */
      value_division(div1, lcm_val, (*i).second.second);
      value_division(div2, lcm_val, (*j).second.second);
      PeriodicNumber pn1 = (*i).second.first * div1;
      PeriodicNumber pn2 = (*j).second.first * div2;
      PeriodicNumber pn1_x_pn2 = pn1*pn2;
      pn1_x_pn2.get_gcd(gcd_of_product);
      value_gcd(gcd_product_denominator, gcd_of_product, common_denominator);
      pn1_x_pn2.divide_by(gcd_product_denominator);
      value_division(common_denominator, common_denominator, gcd_product_denominator);
      EhrhartPolynom term(exponent,
			  periodic_rational(pn1_x_pn2, common_denominator));
      result+=term;
    }
  }
  value_clear(div1); value_clear(div2); value_clear(lcm_val); value_clear(gcd_product_denominator);

  result.simplify();
  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=2) {
    cout<<"In EhrhartPolynom()::operator*: result of multiplication:"<<endl;
    cout<<result.to_string()<<endl;
  }
  return result;
}

EhrhartPolynom EhrhartPolynom::operator/(const Value v) const
{
  EhrhartPolynom result(*this);
  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=2) {
    cout<<"In EhrhartPolynom()::operator/: dividing ";
    cout<<this->to_string()<<" by "<<value2string(v)<<endl;
  }
  for(map<map<string, int>, periodic_rational>::iterator i=
	result.data.begin(); i!= result.data.end(); i++) {
    value_multiply ( (*i).second.second, (*i).second.second, v);
  }

  result.simplify();
  return result;
}

 
void EhrhartPolynom::simplify()
{
  Value mult_p1, mult_p2;
  value_init(mult_p1);
  value_init(mult_p2);
  Value lcm_denominator;
  value_init(lcm_denominator);
  map<map<string, int>, periodic_rational> new_data;
  for(map<map<string, int>, periodic_rational>::const_iterator 
	i=data.begin(); i!=data.end();i++) {
    if ((*i).second.first.is_zero())
      continue;
    // remove exponents with power 0, remove exponents with name ""
    map<string,int> exponent;
    for(map<string,int>::const_iterator j=(*i).first.begin();
	j!=(*i).first.end(); j++) {
      if ((*j).second!=0 && (*j).first!=string(""))
	exponent[(*j).first]=(*j).second;
    }
    if (new_data.find(exponent)!=new_data.end()) {
      value_lcm(lcm_denominator, new_data[exponent].second, (*i).second.second);
      value_division(mult_p1, lcm_denominator, new_data[exponent].second);
      value_division(mult_p2, lcm_denominator, (*i).second.second);
      /*
      mult_p1=lcm_denominator / new_data[exponent].second;
      mult_p2=lcm_denominator / (*i).second.second;
      */
      PeriodicNumber pn = (*i).second.first*mult_p2 +
	                  new_data[exponent].first*mult_p1;
      new_data[exponent] = periodic_rational(pn, lcm_denominator);
    } else {
      new_data[exponent] = (*i).second;
    }

  }
  data=new_data;
  value_clear(mult_p1);
  value_clear(mult_p2);
  value_clear(lcm_denominator);
  if (data.size()==0) {
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
    return;
  }
}

double EhrhartPolynom::to_double() const {
  // Should only be called on constant EhrhartPolynomials (i.e. no parameters)!
  EhrhartPolynom ep1=(*this);
  ep1.simplify();

  assert(ep1.data.size()==1);
  const map<string,int> exponent;
  periodic_rational ep1_val = ep1.data[exponent];
  return ep1_val.first.to_double() / VALUE_TO_DOUBLE(ep1_val.second);
}

bool EhrhartPolynom::smaller_than(const EhrhartPolynom& _ep2) const {
  EhrhartPolynom ep1=(*this);
  EhrhartPolynom ep2=_ep2;
  ep1.simplify();
  ep2.simplify();

  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<"In EhrhartPolynom::smaller_than, comparing "
	<<ep1.to_string()<<" and "<<ep2.to_string()<<endl;
  }
  /* only return true if both are comparable, i.e. both
     ep1 and ep2 are constant, and ep1 < ep2 */
  // Value ep1v, ep2v;
  /* check if constant */
  bool result;
  if (ep1.data.size()==1 && ep2.data.size()==1) {
    const map<string,int> exponent;
    periodic_rational ep1_val = ep1.data[exponent];
    periodic_rational ep2_val = ep2.data[exponent];
    PeriodicNumber ep1_scaled=ep1_val.first * ep2_val.second,
      ep2_scaled=ep2_val.first * ep1_val.second;
    result = ep1_scaled.smaller_than(ep2_scaled);
  } else
    result = false;
  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<"result of comparison is "<<result<<endl;
  }
  return result;
}

bool EhrhartPolynom::operator==(const EhrhartPolynom& _ep2) const {
  EhrhartPolynom ep1=(*this);
  EhrhartPolynom ep2=_ep2;
  ep1.simplify();
  ep2.simplify();

  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<"In EhrhartPolynom::operator==, comparing "
	<<ep1.to_string()<<" and "<<ep2.to_string()<<endl;
  }

  map<map<string, int>, periodic_rational>::const_iterator i=data.begin();
  map<map<string, int>, periodic_rational>::const_iterator j=_ep2.data.begin();
  for(; i!=data.end() && j!=_ep2.data.end(); i++, j++) {
    if ((*i).first!=(*j).first ||
	(*i).second.first != (*j).second.first ||
	value_ne((*i).second.second, (*j).second.second))
      {
	if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
	  cout<<"result 1 of comparison is false."<<endl;
	}
	return false;
      }
  }
  if (i!=data.end() || j!=_ep2.data.end()) {
    if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
      cout<<"result 2 of comparison is false."<<endl;
    }
    return false;
  }
  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<"result 3 of comparison is true."<<endl;
  }
  return true;
}

static void fill_in_d(evalue& e, const Value d)
{
  if (value_zero_p(e.d)) {
    for(int i=0;i<e.x.p->size;i++)
      fill_in_d(e.x.p->arr[i], d);
  } else {
    assert(value_one_p(e.d));
    value_gcd(e.d, e.x.n, d);
    value_division(e.x.n, e.x.n, e.d);
    value_division(e.d, d, e.d);
  }
}

static evalue create_coef
(const deque<string>& parameter_names,
 const EhrhartPolynom::periodic_rational& d)
{
  evalue result=d.first.to_evalue(parameter_names);
  fill_in_d(result, d.second);
  return result;
}

/*
static set<map<string,int> > find_terms_with_var
(const set<map<string, int> >& terms,
 const string& varname)
{
  set<map<string, int> > result;
  for(set<map<string,int> >::const_iterator i=terms.begin();
      i!=terms.end(); i++) {
    if ((*i).find(varname)!=(*i).end())
      result.insert(*i);
  }
  return result;
}
*/

static set<map<string,int> > find_terms_with_var_exp
(const set<map<string,int> >&terms,
 const string& varname,
 const int exp)
{
  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<"In find_terms_with_var_exp:"<<endl;
    cout<<" terms("<<terms.size()<<") = ";
    for(set<map<string,int> >::const_iterator i=terms.begin();
	i!=terms.end(); i++) {
      for(map<string,int>::const_iterator j=(*i).begin();
	  j!=(*i).end(); j++)
	cout<<"'"<<(*j).first<<"^"<<(*j).second<<"'";
      cout<<" ";
    }
    cout<<endl;
    cout<<"varname = '"<<varname<<"'"<<endl;
    cout<<"exp = "<<exp<<endl;
  }
  set<map<string, int> > result;
  if (exp!=0) {
    for(set<map<string,int> >::const_iterator i=terms.begin();
	i!=terms.end(); i++) {
      if ((*i).find(varname)!=(*i).end() && 
	  (*( (*i).find(varname) )).second==exp)
	result.insert(*i);
    }
  } else {
    // exp==0
    for(set<map<string,int> >::const_iterator i=terms.begin();
	i!=terms.end(); i++) {
      if ((*i).find(varname)==(*i).end() ||
	  (*( (*i).find(varname) )).second==exp)
	result.insert(*i);
    }
  }
  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<" result("<<result.size()<<") = ";
    for(set<map<string,int> >::const_iterator i=result.begin();
	i!=result.end(); i++) {
      for(map<string,int>::const_iterator j=(*i).begin();
	  j!=(*i).end(); j++)
	cout<<"'"<<(*j).first<<"^"<<(*j).second<<"'";
      cout<<" ";
    }
    cout<<endl;
  }
  return result;
}


static int find_max_exp(const set<map<string,int> >& terms,
			const string& var_name)
{
  int result=0;
  for(set<map<string,int> >::const_iterator i=terms.begin();
      i!=terms.end(); i++) {
    map<string,int>::const_iterator j=(*i).find(var_name);
    if (j!=(*i).end())
      result = std::max(result, (*j).second);
  }
  return result;
}


evalue
EhrhartPolynom::translate_one_term
(const deque<string>& parameter_names,
 const deque<string>& _left_over_var_names,
 const set<map<string, int> >& terms) const
{
  static int depth=0;
  depth++;
  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=2) {
    cout<<"In EhrhartPolynom::translate_one_term(depth = "
	<<depth<<"):"<<endl;
    cout<<to_string()<<endl;
    cout<<"left_over_var_names = ";
    for(unsigned i=0;i<_left_over_var_names.size();i++)
      cout<<_left_over_var_names[i]<<" ";
    cout<<endl;
    cout<<" terms("<<terms.size()<<") = ";
    for(set<map<string,int> >::const_iterator i=terms.begin();
	i!=terms.end(); i++) {
      for(map<string,int>::const_iterator j=(*i).begin();
	  j!=(*i).end(); j++)
	cout<<"'"<<(*j).first<<"^"<<(*j).second<<"'";
      cout<<" ";
    }
    cout<<endl;
  }
  deque<string> left_over_var_names=_left_over_var_names;
  if (left_over_var_names.size()==0) {
    assert(terms.size()<=1);
    if (terms.size()==0) {
      evalue zero_val;
      value_init(zero_val.d);
      value_init(zero_val.x.n);
      value_set_si(zero_val.d,1);
      value_set_si(zero_val.x.n,0);
      depth--;
      return zero_val; // no term with the given coeffs.
    }
    assert(terms.size()==1);
    map<map<string,int>, periodic_rational>::const_iterator t=
      data.find(*(terms.begin()));
    //assert((*t).first==(*(terms.begin())));
    depth--;
    return create_coef(parameter_names, (*t).second);
  } else {
    string var_name = left_over_var_names.front();
    left_over_var_names.pop_front();
    int exp=find_max_exp(terms, var_name);
    assert(exp!=-1);
    int pos=find_pos_of(var_name, parameter_names);
    assert(pos!=-1);
    enode* result = new_enode(polynomial, exp+1, pos+1/*from 1 to m*/);
    for(int i=0;i<exp+1;i++) {
      set<map<string,int> > new_terms = 
	find_terms_with_var_exp(terms, var_name, i);
      result->arr[i] = translate_one_term(parameter_names,
					  left_over_var_names,
					  new_terms);
    }
    evalue t_result;
    value_init(t_result.d);
    value_set_si(t_result.d , 0);
    t_result.x.p=result;
    depth--;
    return t_result;
  }
}

::evalue EhrhartPolynom::to_evalue(const deque<string>& parameter_names) const
{
  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=2) {
    cout<<"In EhrhartPolynom::to_evalue: this="<<endl;
    cout<<to_string()<<endl;
  }
  ::evalue result;
  set<string> var_names;
  for(map<map<string, int>, periodic_rational>::const_iterator 
	i=data.begin(); i!=data.end(); i++) {
    for(map<string,int>::const_iterator j=(*i).first.begin();
	j!=(*i).first.end(); j++)
      var_names.insert((*j).first);
  }
  deque<string> var_names_dq;
  for(set<string>::const_iterator i=var_names.begin();
      i!=var_names.end(); i++)
    var_names_dq.push_back(*i);

  set<map<string,int> > terms;
  for(map<map<string,int>, periodic_rational>::const_iterator i=
	data.begin(); i!=data.end(); i++) {
    terms.insert((*i).first);
  }
  result=translate_one_term(parameter_names, var_names_dq,
			    terms);

  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=2) {
    cout<<"In EhrhartPolynom::to_evalue: result="<<endl;
    char** pname=get_pname(parameter_names);
    print_evalue(stdout, &result, pname);
    delete_pname(pname, parameter_names);
    cout<<endl;
  }

  return result;
}

string EhrhartPolynom::to_string() const
{
  if (data.size()==0)
    return "0_size0";
  string result;
  for(map<map<string, int>, periodic_rational>::const_iterator
	i=data.begin();
      i!=data.end(); i++) {
    if (i!=data.begin())
      result+=" + ";
    // first print coefficient
    if (value_notone_p((*i).second.second)) {
      char *str;
      value_to_str(str, 0, (*i).second.second);
      result+= "(1/"+string(str)+")*";
      free(str);
    }
    result+=(*i).second.first.to_string()+" ";
    for(map<string,int>::const_iterator j=(*i).first.begin();
	j!=(*i).first.end(); j++) 
      if ((*j).second==1)
	result+=string("*")+(*j).first.c_str();
      else
	result+=string("*")+(*j).first+"^"+int2string((*j).second);
  }
  // strip of spaces at the end
  //result+=" ";
  return result.substr(0,result.find_last_not_of(" ")+1);
}

bool EhrhartPolynom::contains_parameters() const
{
  if (data.size()==0)
    return false;
  Value zero;
  value_init(zero);
  value_set_si(zero, 0);
  PeriodicNumber pn_zero(zero);
  value_clear(zero);
  for(map<map<string, int>, periodic_rational>::const_iterator
	i=data.begin();
      i!=data.end(); i++) {
    // first print coefficient
    if ((*i).second.first == PeriodicNumber(zero))
      continue;
    for(map<string,int>::const_iterator j=(*i).first.begin();
	j!=(*i).first.end(); j++) 
      if ((*j).first!="")
	return true;
  }
  return false;
}


EhrhartPolynom EhrhartPolynom::subst_var
(const string& var,
 const map<string,int>& substitution,
 const int divisor) const 
{
  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=2) {
    cout<<"In EhrhartPolynom::subst_var: this="<<(string)(*this)<<endl;
    cout<<"var="<<var<<endl;
    cout<<"substitution=";
    for(map<string,int>::const_iterator i=substitution.begin();
	i!=substitution.end(); i++) {
      cout<<" "<<(((*i).second>=0)?"+":"")<<(*i).second
	  <<(*i).first;
    }
    cout<<" / "<<divisor;
    cout<<endl;
  }
  // split the current polynomial into 2 polynomials:
  // ep_with_var contains all the terms with variable
  // var, ep_without_var contains all the terms without
  // the variable var.
  EhrhartPolynom ep_with_var;
  EhrhartPolynom ep_without_var;

  for(map<map<string, int>, periodic_rational>::const_iterator
	i=data.begin(); i!=data.end(); i++) {
    if ((*i).first.find(var)==(*i).first.end())
      ep_without_var += EhrhartPolynom((*i).first, (*i).second);
    else
      ep_with_var += EhrhartPolynom((*i).first, (*i).second);
  }
  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<"In EhrhartPolynom::subst_var:\n ep_with_var="
	<<(string)(ep_with_var)<<endl;
    cout<<"  ep_without_var="
	<<(string)(ep_without_var)<<endl;
  }

  if (ep_with_var == EhrhartPolynom(0))
    return (*this);

  // find the maximum exponent to which var is powered in the Ehrhartpolynom
  unsigned maxexp=0;
  for(map<map<string, int>, periodic_rational>::const_iterator
	i=ep_with_var.data.begin(); i!=ep_with_var.data.end(); i++) {
    maxexp = std::max(maxexp, (unsigned) (*((*i).first.find(var))).second);
  }
  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<"In EhrhartPolynom::subst_var:\n maxexp="<<maxexp<<endl;
  }
  // create for powers of 1 to n, the powers of substitution.
  deque<EhrhartPolynom> subt_powers;
  
  
  map<string,int> one; one[""]=1;
  // subst^0 = r
  Value val_one;
  value_init(val_one);
  value_set_si(val_one, 1);
  subt_powers.push_back(EhrhartPolynom(one,
				       periodic_rational(PeriodicNumber(val_one),val_one)));
  // subst^1 = subst
  EhrhartPolynom subst;
  Value tmp_val, div_val;
  value_init(tmp_val);
  value_init(div_val);
  value_set_si(div_val, divisor);
  for(map<string,int>::const_iterator i=substitution.begin();
      i!=substitution.end(); i++) {
    map<string,int> exponent;
    exponent[(*i).first]=1;
    value_set_si(tmp_val, (*i).second);
    subst += EhrhartPolynom(exponent,
			    periodic_rational(PeriodicNumber(tmp_val),div_val));
  }
  value_clear(div_val);
  value_clear(tmp_val);
  value_clear(val_one);
  subt_powers.push_back(subst);
  for(unsigned i=2;i<=maxexp;i++) {
    // calculate subt_powers[i]
    subt_powers.push_back(subt_powers[i-1]*subt_powers[1]);
  }
  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<"In EhrhartPolynom::subst_var:\n subt_powers:"<<endl;
    for(unsigned power=0;power<=maxexp;power++) {
      cout<<"subt_powers["<<power<<"]= "
	  <<subt_powers[power].to_string()<<endl;
    }
  }
  
  EhrhartPolynom result=ep_without_var;
  
  for(map<map<string, int>, periodic_rational>::const_iterator
	i=ep_with_var.data.begin(); i!=ep_with_var.data.end(); i++) {
    assert((*i).first.find(var)!=(*i).first.end());
    map<string,int> new_exp=(*i).first;
    unsigned power=new_exp[var];
    assert(power<=maxexp);
    new_exp.erase(new_exp.find(var));
    EhrhartPolynom term(new_exp, (*i).second);
    term *= subt_powers[power];
    result+=term;
  }
  result.simplify();
  
  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<"In EhrhartPolynom::subst_var:\n result=:"<<endl;
    cout<<result.to_string()<<endl;
  }
  return result;
}

static inline bool bit(unsigned long value, int bitnr)
{
  return (value>>bitnr)&1;
}

#if 0
static bool is_constant(const ::Relation& r,
			const int set_var,
			const string& char_name,
			coef_t& value)
{
  // first project away all other variables
  ::Relation t(r);
  for(int i=1;i<=r.n_set(); i++) {
    if (i<set_var) {
      // project away current var 1
      t = Project(t, i, Set_Var);
    } else if (i>set_var) {
      // project away current var 2
      t = Project(t, i, Set_Var);
    }
  }
  // check to see if there's only one EQ-constraint.
  t.simplify();
  
  if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=2) {
    cout<<"in is_constant(,"<<set_var<<","<<char_name<<") after projection, the relation is ";
    t.print();
    cout<<endl;
  }

  // now only the set_var variable remains
  //  assert(t.n_set()==1);

  // there should be only 1 disjunction
  DNF_Iterator di(t.query_DNF());
  if (di==0)
    return false;
  di++;
  if (di!=0)
    return false;
  di = t.query_DNF();

  // There should be no GEQ constraints
  if ((*di)->GEQs() != 0)
    return false;

  // There should be only on EQ constraint
  EQ_Iterator ei=(*di)->EQs();
  if (ei==0)
    return false;
  ei++;
  if (ei!=0)
    return false;
  ei=(*di)->EQs();

  coef_t coeff_of_var=0;
  // now check to see if var equals a constant
  for(Constr_Vars_Iter cvi(*ei); cvi; cvi++) {
    if (string((*cvi).var->char_name()) != char_name) {
      // there are other variables in the EQ constraint
      return false;
    }
    coeff_of_var = (*cvi).coef;
  }

  if (coeff_of_var==0)
    return false;
  
  if (coeff_of_var==1) {
    value = -(*ei).get_const();
    return true;
  }
  if (coeff_of_var==-1) {
    value = (*ei).get_const();
    return true;
  }
  return false;
}


/** equal_in_domain tries to look if the 2 Ehrhart polynomials @a ep1
 * and @a ep2 are equal, taking into account the domain in which the
 * Ehrhart polynomials are defined. e.g. @a ep1 \f$=n-i^2-j\f$ in
 * domain \f$i=n \wedge j=-1\f$ is equal to @a ep2 \f$=-n^2-i-j\f$ in
 * domain \f$i=-n\f$. This is true because, when substituting the
 * value for \f$i\f$ with \f$n/-n\f$ in the 2 polynomials, the same
 * Ehrhart polynomial is generated, i.e. \f$-n^2+n-j\f$.
 *
 * @param ep1 is the first Ehrhart polynomial
 * @param ep2 is the second Ehrhart polynomial
 * @param domain1 is the domain corresponding to @a ep1. The domain must
 * be convex.
 * @param domain2 is the domain corresponding to @a ep2. The domain must
 * be convex.
 *
 * @return If the polynomials are equal inside the domain, a pointer
 * to an Ehrhart polynomials is returned which represents a form of
 * the Ehrhart polynomial which is equal to @a ep1 and @a ep2. In the
 * example above, \f$-n^2+n-j\f$ would be returned. */

EhrhartPolynom* EhrhartPolynom::equal_in_domain
(const EhrhartPolynom& ep1, const EhrhartPolynom& ep2,
 const ::Relation& _domain1, const ::Relation& _domain2)
{
  if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=2) {
    cout<<"In EhrhartPolynom::equal_in_domain"<<endl;
    cout<<"ep1="<<ep1.to_string()<<endl;
    cout<<"ep2="<<ep2.to_string()<<endl;
    cout<<"domain1="; ::copy(_domain1).print();
    cout<<"domain2="; ::copy(_domain2).print();
  }

  if (ep1==ep2) {
    if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=2) {
      cout<<"ep1==ep2; returning ep1";
    }
    return new EhrhartPolynom(ep1);
  }


  // do a very simple check. 
  // If ep2->ep1 in domain1, return ep2.
  // If ep1->ep2 in domain2, return ep1.

  // first for domain1, substitute constant values in both ep1 and ep2
  {
    //    ::Relation tmp(::copy(_domain1));
    EhrhartPolynom tmp_ep1(ep1), tmp_ep2(ep2);
    
    // project to a variable.
    for(int var_to_check=1; var_to_check<= _domain1.n_set(); var_to_check++) {
      ::Relation tmp2(::copy(_domain1));
      ::coef_t lowerbound=0, upperbound=0;
      tmp2.query_variable_bounds(tmp2.set_var(var_to_check),lowerbound, upperbound);
      if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=2) {
	cout <<"Checked variable "<<var_to_check<<" in relation ";
	tmp2.print();
	cout<<" lowerbound="<<lowerbound<<"; upperbound="<<upperbound<<endl;
      }
      if (lowerbound != upperbound) {
	coef_t value;
	bool cst = is_constant(_domain1, var_to_check, tmp2.set_var(var_to_check)->char_name(),
			       value);
	if (cst) 
	  lowerbound = upperbound = value;
      }
      if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=2) {
	cout <<"Checked variable "<<var_to_check<<" in relation ";
	tmp2.print();
	cout<<" lowerbound2="<<lowerbound<<"; upperbound2="<<upperbound<<endl;
      }
      if (lowerbound == upperbound)
	{
	  // the value is constant for this variable.
	  // find the variable's name
	  string var_name = tmp2.set_var(var_to_check)->char_name();
	  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=2)
	    cout <<"variable "<<var_name<<"("<<var_to_check<<") has constant value "
		 <<lowerbound<<endl;
	  map<string, int> substitution;
	  substitution[""] = lowerbound;
	  tmp_ep1=tmp_ep1.subst_var(var_name, substitution);
	  tmp_ep2=tmp_ep2.subst_var(var_name, substitution);
	}
    }
    if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=2) {
      cout<<"In domain1: tmp_ep1="<<tmp_ep1.to_string()<<endl;
      cout<<"            tmp_ep2="<<tmp_ep2.to_string()<<endl;
    }
    
    if (tmp_ep1 == tmp_ep2) {
      // ep1 == ep2 in domain1: return ep2
      if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=2) {
	cout<<"In EhrhartPolynom::equal_in_domain"<<endl;
	cout<<"Found ep1="<<tmp_ep1.to_string()<<endl;
	cout<<"      ep2="<<tmp_ep2.to_string()<<endl;
	cout<<"In domain:";
	::copy(_domain1).print(); cout<<endl;
      }
      return new EhrhartPolynom(ep2);
    }
  }

  // next, check for domain2; substitute constant values in both ep1 and ep2
  {
    //    ::Relation tmp(::copy(_domain2));
    EhrhartPolynom tmp_ep1(ep1), tmp_ep2(ep2);
    
    // project to a variable.
    for(int var_to_check=1; var_to_check<= _domain2.n_set(); var_to_check++) {
      ::Relation tmp2(::copy(_domain2));
      ::coef_t lowerbound, upperbound;
      tmp2.query_variable_bounds(tmp2.set_var(var_to_check),lowerbound, upperbound);
      if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=2) {
	cout <<"Checked variable "<<var_to_check<<" in relation ";
	tmp2.print();
	cout<<" lowerbound="<<lowerbound<<"; upperbound="<<upperbound<<endl;
      }
      if (lowerbound != upperbound) {
	coef_t value;
	bool cst = is_constant(_domain2, var_to_check, tmp2.set_var(var_to_check)->char_name(),
			       value);
	if (cst) 
	  lowerbound = upperbound = value;
      }
      if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=2) {
	cout <<"Checked variable "<<var_to_check<<" in relation ";
	tmp2.print();
	cout<<" lowerbound2="<<lowerbound<<"; upperbound2="<<upperbound<<endl;
      }
      if (lowerbound == upperbound)
	{
	  // the value is constant for this variable.
	  // find the variable's name
	  string var_name = tmp2.set_var(var_to_check)->char_name();
	  if(DebugBlackboard.debug("EHRHARTPOLYNOM")>=2)
	    cout <<"variable "<<var_name<<"("<<var_to_check<<") has constant value "
		 <<lowerbound<<endl;
	  map<string, int> substitution;
	  substitution[""] = lowerbound;
	  tmp_ep1=tmp_ep1.subst_var(var_name, substitution);
	  tmp_ep2=tmp_ep2.subst_var(var_name, substitution);
	}
    }
    if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=2) {
      cout<<"In domain2: tmp_ep1="<<tmp_ep1.to_string()<<endl;
      cout<<"            tmp_ep2="<<tmp_ep2.to_string()<<endl;
    }
    
    if (tmp_ep1 == tmp_ep2) {
      // ep1 == ep2 in domain2: return ep1
      if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=2) {
	cout<<"In EhrhartPolynom::equal_in_domain"<<endl;
	cout<<"Found ep1="<<tmp_ep1.to_string()<<endl;
	cout<<"      ep2="<<tmp_ep2.to_string()<<endl;
	cout<<"In domain:";
	::copy(_domain2).print(); cout<<endl;
      }
      return new EhrhartPolynom(ep1);
    }
  }

  /*
  ::Relation domain1=::copy(_domain1);
  ::Relation domain2=::copy(_domain2);

  map<string, map<string,int> > var2subst_in_domain1;
  map<string, map<string,int> > var2subst_in_domain2;

  set<string> vars_in_domain1;
  set<string> vars_in_domain2;

  // make sure that that GEQ constraints are translated into
  // EQ constraints, where possible.
  domain1.simplify(2,2);
  domain2.simplify(2,2);

  if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=2) {
    cout<<"domain1=";domain1.print(); cout<<endl;
    cout<<"domain2=";domain2.print(); cout<<endl;
  }

  int counter1=0;
  for(DNF_Iterator di(domain1.query_DNF()); di; di++) {
    assert(counter1==0); // domain1 must be convex!
    counter1++;
    for(EQ_Iterator ei=(*di)->EQs(); ei; ei++) {
      map<string,int> equalities;
      int const_coef;
      for(Constr_Vars_Iter cvi(*ei); cvi; cvi++) {
	if (string((*cvi).var->char_name()) ==
	    string(""))
	  goto next_equality;
	equalities[(*cvi).var->char_name()] = (*cvi).coef;
      }
      const_coef = (*ei).get_const();
      // insert equalities into var2subst_in_domain1
      for(map<string,int>::const_iterator i=equalities.begin();
	  i!=equalities.end(); i++) {
	map<string,int> subst_coefs;
	// for now, only handle variables with coeffs 1 or -1
	// TODO: handle other coeffs
	if ((*i).second==1 || (*i).second==-1) {
	  for(map<string,int>::const_iterator j=equalities.begin();
	      j!=equalities.end(); j++) {
	    if (i!=j) {
	      if ((*i).second<0) subst_coefs[(*j).first]=(*j).second;
	      else subst_coefs[(*j).first]= -(*j).second;
	    }
	  }
	  if ((*i).second<0) subst_coefs[""] = const_coef;
	  else subst_coefs[""] = -const_coef;
	  var2subst_in_domain1[(*i).first] = subst_coefs;
	}
      }
      
    next_equality: ;
    }
  }

  int counter2=0;
  for(DNF_Iterator di(domain2.query_DNF()); di; di++) {
    assert(counter2==0); // domain2 must be convex!
    counter2++;
    for(EQ_Iterator ei=(*di)->EQs(); ei; ei++) {
      map<string,int> equalities;
      int const_coef;
      for(Constr_Vars_Iter cvi(*ei); cvi; cvi++) {
	if (string((*cvi).var->char_name()) == string(""))
	  goto next_equality2;
	equalities[(*cvi).var->char_name()] = (*cvi).coef;
      }
      const_coef = (*ei).get_const();
      // insert equalities into var2subst_in_domain1
      for(map<string,int>::const_iterator i=equalities.begin();
	  i!=equalities.end(); i++) {
	map<string,int> subst_coefs;
	// for now, only handle variables with coeffs 1 or -1
	// TODO: handle other coeffs
	if ((*i).second==1 || (*i).second==-1) {
	  for(map<string,int>::const_iterator j=equalities.begin();
	      j!=equalities.end(); j++) {
	    if (i!=j) {
	      if ((*i).second<0) subst_coefs[(*j).first]=(*j).second;
	      else subst_coefs[(*j).first]= -(*j).second;
	    }
	  }
	  if ((*i).second<0) subst_coefs[""] = const_coef;
	  else subst_coefs[""] = -const_coef;
	  var2subst_in_domain2[(*i).first] = subst_coefs;
	}
      }
    next_equality2: ;
    }
  }

  if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=2) {
    cout<<"var2subst_in_domain1("
	<<var2subst_in_domain1.size()<<")"<<endl;
    cout<<"var2subst_in_domain2("
	<<var2subst_in_domain2.size()<<")"<<endl;
    for(map<string, map<string,int> >::const_iterator 
	  i=var2subst_in_domain1.begin();
	i!=var2subst_in_domain1.end() ;i++) {
      cout<<"var2subst_in_domain1["<<(*i).first<<"]= {";
      for(map<string,int>::const_iterator j=(*i).second.begin();
	  j!=(*i).second.end(); j++)
	cout<<"+"<<(*j).second<<(*j).first;
      cout<<"}"<<endl;
    }
    for(map<string, map<string,int> >::const_iterator 
	  i=var2subst_in_domain2.begin();
	i!=var2subst_in_domain2.end() ;i++) {
      cout<<"var2subst_in_domain2["<<(*i).first<<"]= {";
      for(map<string,int>::const_iterator j=(*i).second.begin();
	  j!=(*i).second.end(); j++)
	cout<<"+"<<(*j).second<<(*j).first;
      cout<<"}"<<endl;
    }
  }  

  // now check for every possible combination of substitutions if
  // equality between Ehrhart Polynomials can be reached.
  if (var2subst_in_domain2.size()>=30 ||
      var2subst_in_domain1.size()>=30) {
    // to many combinations, bail out
    if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=3)
      cout<<"too many combinations, bailing out."<<endl;
    return 0;
  }
  unsigned long nr_combinations_ep1 = 1<<var2subst_in_domain1.size();
  unsigned long nr_combinations_ep2 = 1<<var2subst_in_domain2.size();
  // generate combinations
  deque<EhrhartPolynom> ep1_variations;
  deque<EhrhartPolynom> ep2_variations;
  for(unsigned long i=0;i<nr_combinations_ep1;i++) {
    EhrhartPolynom ep(ep1);
    unsigned long counter=0;
    for(map<string, map<string,int> >::const_iterator 
	  j=var2subst_in_domain1.begin();
	j!=var2subst_in_domain1.end(); j++, counter++) {
      if (bit(i,counter))
	ep=ep.subst_var((*j).first, (*j).second);
    }
    ep1_variations.push_back(ep);
  }
  for(unsigned long i=0;i<nr_combinations_ep2;i++) {
    EhrhartPolynom ep(ep2);
    unsigned long counter=0;
    for(map<string, map<string,int> >::const_iterator 
	  j=var2subst_in_domain2.begin();
	j!=var2subst_in_domain2.end(); j++, counter++) {
      if (bit(i,counter))
	ep=ep.subst_var((*j).first, (*j).second);
    }
    ep2_variations.push_back(ep);
  }
  
  if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=3) {
    cout<<"ep1_variations("<<ep1_variations.size()<<"):"<<endl;
    for(unsigned long i=0;i<ep1_variations.size(); i++)
      cout<<ep1_variations[i].to_string()<<endl;
    cout<<"ep2_variations("<<ep2_variations.size()<<"):"<<endl;
    for(unsigned long i=0;i<ep2_variations.size(); i++)
      cout<<ep2_variations[i].to_string()<<endl;
  }

  for(unsigned long i=0;i<ep1_variations.size();i++)
    for(unsigned long j=0;j<ep2_variations.size();j++) {
      if (ep1_variations[i]==ep2_variations[j]) {
	if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=2) {
	  cout<<"Returning "<<ep1_variations[i].to_string()<<endl;
	}
	return new EhrhartPolynom(ep1_variations[i]);
      }
    }
  if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=2) {
    cout<<"Returning null"<<endl;
  }
  */
  if (DebugBlackboard.debug("EHRHARTPOLYNOM")>=2) {
    cout<<"not equal."<<endl;
  }
  return 0;
}
#endif 

/**
 * get_affine_polynomials returns the Ehrhart polynomials as a set
 * of affine functions when possible. If not possible, it returns
 * the null pointer.
 */
struct EhrhartPolynom::affine*
EhrhartPolynom::get_affine_polynomials() const
{
  if (DebugBlackboard.debug("AFFINEPOLYNOM")>=2) {
    cout<<"In EhrhartPolynom::get_affine_polynomials: this="
	<<to_string()<<endl;
  }

  deque<string> var_names;
  // first check if this polynomial is representable as an
  // affine polynomial.
  for(map<map<string,int>, periodic_rational>::const_iterator 
	i=data.begin(); i!=data.end(); i++) {
    unsigned nr_vars=0;
    if ((*i).first.size()==0) {
      //constant term
      var_names.push_back("");
    } else {
      for(map<string,int>::const_iterator j=(*i).first.begin();
	  j!=(*i).first.end(); j++) {
	if ((*j).second > 1) // power >1 => not affine
	  {
	    if (DebugBlackboard.debug("AFFINEPOLYNOM")>=2) {
	      cout<<"Not affine since power of "<<(*j).first<<" is "
		  <<(*j).second<<" >1."<<endl;
	    }
	    return 0;
	  }
	if ((*j).first!=string("")) {
	  nr_vars++;
	  if (std::find(var_names.begin(), var_names.end(), (*j).first) ==
	      var_names.end())
	    var_names.push_back((*j).first);
	}
      }
      if (nr_vars > 1) // more than 1 variable in the term => not affine
	{
	  if (DebugBlackboard.debug("AFFINEPOLYNOM")>=2) {
	    cout<<"Not affine since nr_vars="<<nr_vars<<" >1."<<endl;
	  }
	  return 0;
	}
    }
  }

  if (DebugBlackboard.debug("AFFINEPOLYNOM")>=4) {
    cout<<"var_names("<<var_names.size()<<")=";
    for(deque<string>::const_iterator i=var_names.begin();
	i!=var_names.end(); i++) {
      cout<<"'"<<(*i)<<"' ";
    }
    cout<<endl;
  }

  struct affine* result = new affine;
  deque<string> var_names_as_param;
  set<string> params_set;
  for(map<map<string,int>, periodic_rational>::const_iterator 
	i=data.begin(); i!=data.end(); i++) {
    set<string> p1=(*i).second.first.get_params();
    for(set<string>::const_iterator j=p1.begin(); j!=p1.end(); j++)
      params_set.insert(*j);
  }
  for(set<string>::const_iterator i=params_set.begin();
      i!=params_set.end(); i++)
    var_names_as_param.push_back(*i);

  map<string, set<int> > periods;
  
  for(map<map<string,int>, periodic_rational>::const_iterator 
	i=data.begin(); i!=data.end(); i++) {
    if (DebugBlackboard.debug("AFFINEPOLYNOM")>=6) {
      cout<<"Checking params in "<<(*i).second.first.to_string()<<":";
    }
    for(unsigned  j=0; j!=var_names_as_param.size(); j++) {
      if ((*i).second.first.has_param(var_names_as_param[j])) {
	if (DebugBlackboard.debug("AFFINEPOLYNOM")>=6) {
	  cout<<"'"<<var_names_as_param[j]<<"' ";
	}
	periods[var_names_as_param[j]].insert
	  ((*i).second.first.get_period(var_names_as_param[j]));
      }
    }
    if (DebugBlackboard.debug("AFFINEPOLYNOM")>=6) {
      cout<<endl;
    }
  }


  if (DebugBlackboard.debug("AFFINEPOLYNOM")>=4) {
    cout<<"var_names_as_param("
	<<var_names_as_param.size()<<")=";
    for(unsigned i=0;i<var_names_as_param.size(); i++) {
      cout<<"'"<<var_names_as_param[i]<<"' ";
    }
    cout<<endl;
    cout<<"periods=";
    for(map<string,set<int> >::const_iterator i=periods.begin();
	i!=periods.end(); i++) {
      cout<<"('"<<(*i).first<<"',";
      for(set<int>::const_iterator j=(*i).second.begin();
	  j!=(*i).second.end(); j++)
	cout<<(*j)<<" ";
      cout<<") ";;
    }
    cout<<endl;
  }

  // for every variable, find the lcm of the periods.
  map<string, int> lcm_period;
  for(map<string, set<int> >::const_iterator i=periods.begin();
      i!=periods.end(); i++) {
    int plcm=1;
    for(set<int>::const_iterator j=(*i).second.begin();
	j!=(*i).second.end(); j++)
      plcm=lcm(*j, plcm);
    lcm_period[(*i).first] = plcm;
  }
  
  if (DebugBlackboard.debug("AFFINEPOLYNOM")>=3) {
    cout<<"lcm_period=";
    for(map<string,int>::const_iterator i=lcm_period.begin();
	i!=lcm_period.end(); i++)
      cout<<"('"<<(*i).first<<"',"<<(*i).second<<") ";
    cout<<endl;
  }


  result->period = lcm_period;

  // start creating the affine structure:
  // the number of combinations to be examined are
  // all the compination of periods for the different variables.
  long nr_comb=1;
  for(map<string,int>::const_iterator i=lcm_period.begin();
      i!=lcm_period.end(); i++) {
    nr_comb *= (*i).second;
  }
  if (DebugBlackboard.debug("AFFINEPOLYNOM")>=3) {
    cout<<"nr_comb="<<nr_comb<<endl;
  }
  deque<long> strides;
  for(unsigned j=0;j<var_names_as_param.size();j++) {
    if (j==0)
      strides.push_back(1);
    else
      strides.push_back(strides[j-1] * lcm_period[var_names_as_param[j-1]]);
  }
  if (var_names_as_param.size()==0) {
    assert(nr_comb==1);
  }
  else
    strides.push_back(strides.back() * lcm_period[var_names_as_param.back()]);
  if (DebugBlackboard.debug("AFFINEPOLYNOM")>=3) {
    cout<<"strides=";
    for(deque<long>::const_iterator i=strides.begin(); i!=strides.end(); i++)
      cout<<(*i)<<" ";
    cout<<endl;
  }
  
  for(long i=0;i<nr_comb;i++) {
    // create index
    if (DebugBlackboard.debug("AFFINEPOLYNOM")>=4) {
      cout<<"i="<<i<<endl;
    }
    map<string,int> index;
    for(unsigned j=0;j<var_names_as_param.size();j++) {
      // compute value of
      // var_names_as_param[i] % lcm_period[var_names_as_param[i]]
      assert(strides.size()>j+1);
      int remainder = (i%(strides[j+1]) / strides[j]);
      index[var_names_as_param[j]] = remainder;
    }
    // now find the coefficients for every term and place it
    // in the struct affine result.
    map<string,int> affine_coef;
    Value lcm_denumerator;
    value_init(lcm_denumerator);
    value_set_si(lcm_denumerator, 1);
    for(map<map<string,int>, periodic_rational>::const_iterator
	  j=data.begin(); j!=data.end(); j++) {
      value_lcm(lcm_denumerator ,lcm_denumerator, (*j).second.second);
    }

    if (DebugBlackboard.debug("AFFINEPOLYNOM")>=4) {
      char* str;
      value_to_str(str, 0, lcm_denumerator);
      cout<<"lcm_denumerator="<<str<<endl;
      free(str);
    }
    Value tmp;
    value_init(tmp);
    for(map<map<string,int>, periodic_rational>::const_iterator
	  j=data.begin(); j!=data.end(); j++) {
      string var="";
      for(map<string,int>::const_iterator k=(*j).first.begin();
	  k!=(*j).first.end(); k++) {
	assert(var=="");
	assert((*k).second==1);
	var=(*k).first;
      }
      // filter out the parameters not occuring in the PeriodicNumber
      map<string,int> part_index;
      const PeriodicNumber& pn=(*j).second.first;
      for(map<string,int>::const_iterator k=index.begin();
	  k!=index.end(); k++) {
	if (pn.has_param((*k).first))
	  part_index[(*k).first] = (*k).second;
      }
      if (DebugBlackboard.debug("AFFINEPOLYNOM")>=4) {
	cout<<"part_index= ";
	for(map<string,int>::const_iterator l=part_index.begin();
	    l!=part_index.end(); l++)
	  cout<<"('"<<(*l).first<<"',"<<(*l).second<<") ";
	cout<<endl;
      }
      value_multiply(tmp, pn[part_index], lcm_denumerator);
      value_division(tmp, tmp, (*j).second.second);
      affine_coef[var] = VALUE_TO_INT(tmp); //pn[part_index]*lcm_denumerator/(*j).second.second;
      if (DebugBlackboard.debug("AFFINEPOLYNOM")>=4) {
	cout<<"affine_coef['"<<var<<"']="<<affine_coef[var]<<endl;
      }
    }
    value_clear(tmp);
    // only push back a polynomial if it is different from 0
    bool is_zero=true;
    for(map<string,int>::const_iterator j=affine_coef.begin();
	j!=affine_coef.end(); j++)
      if ((*j).second!=0) {
	is_zero=false;
	break;
      }
    if (!is_zero) {
      result->polynomial.push_back(affine_coef);
      result->offset.push_back(index);
      result->denumerator.push_back(VALUE_TO_INT(lcm_denumerator));
    }
    value_clear(lcm_denumerator);
  }
  
  if (DebugBlackboard.debug("AFFINEPOLYNOM")>=3) {
    cout<<"In EhrhartPolynom::get_affine_polynomials: returning"<<endl;
    cout<<result->to_string()<<endl;
  }
  return result;
}

string EhrhartPolynom::affine::to_string() const {
  string result;
  result+="period=";
  for(map<string,int>::const_iterator i=period.begin();
      i!=period.end(); i++)
    result+="('"+(*i).first+"',"+int2string((*i).second)+") ";
  result+="\n";
  result+="nr_polynomials = "+int2string(size())+"\n";
  for(unsigned i=0;i<size();i++) {
    result+="polynomial["+int2string(i)+"]= (";
    for(map<string,int>::const_iterator j=polynomial[i].begin();
	j!=polynomial[i].end(); j++) {
      result+= (((*j).second>=0)?"+":"")+
	int2string((*j).second)+(*j).first;
    }
    result+=")/"+int2string(denumerator[i]);
    result+="  iff  ";
    for(map<string, int>::const_iterator k=offset[i].begin();
	k!=offset[i].end(); k++) {
      assert(period.find((*k).first)!=period.end());
      int p=(*(period.find((*k).first))).second;
      result+="("+(*k).first+"%"+int2string(p)+"="+int2string((*k).second)+")";
    }
    result+="\n";
  }
  
  return result;
}

#if 0
Expression* EhrhartPolynom::get_AST_representation() const
{
  // Currently, this function is only implemented for Ehrhart
  // Polynomials with period 1.
  // TODO: handle polynomials with larger period.

  if (DebugBlackboard.debug("EHRHARTxTOxAST")>=3) {
    cout<<"In EhrhartPolynom::get_AST_representation: about to convert"<<endl;
    cout<<this->to_string()<<endl;
  }

  for(map<map<string, int>, periodic_rational>::const_iterator i=data.begin();
      i!=data.end(); i++) {
    if ((*i).second.first.datasize > 1) {
      // The period is larger than 1, which cannot handle right now.
      // so return a null pointer.
      return 0;
    }
  }


  Expression* expression = 0;
  // create the different terms in the polytope
  Value coef;
  value_init(coef);
  for(map<map<string, int>, periodic_rational>::const_iterator i=data.begin();
      i!=data.end(); i++) {
    assert((*i).second.first.datasize == 1);
    //    assert((*i).second.second == 1);
    // get the coefficient of the term
    value_assign(coef, (*i).second.first.data[0]);
    Expression* term = 0;
    if (value_notone_p(coef) || (*i).first.size()==0)
      term = new Constant((long long)VALUE_TO_LONG(coef));
    // build the terms
    for(map<string, int>::const_iterator j=(*i).first.begin();
	j!=(*i).first.end(); j++) {
      for(int p=0;p<(*j).second;p++) {
	// (*j).second is the exponent of the variable in the term.
	if (term!=0)
	  term = new BinaryExpression(term,'*', new ScalarRef((*j).first.c_str()));
	else
	  term = new ScalarRef((*j).first.c_str());
      }
    }
    // build the overall polynomial
    if (expression==0)
      expression=term;
    else
      expression = new BinaryExpression(expression, '+', term);
    if (value_notone_p((*i).second.second)) {
      expression = new BinaryExpression(expression, '/',
					new Constant((long long)VALUE_TO_LONG((*i).second.second)));
    }
  }
  value_clear(coef);
  if (data.size()==0)
    expression = new Constant(0LL);

  if (DebugBlackboard.debug("EHRHARTxTOxAST")>=3) {
    cout<<"In EhrhartPolynom::get_AST_representation: converted into"<<endl;
    cout<<expression->to_string()<<endl;
  }

  return expression;
}
#endif
