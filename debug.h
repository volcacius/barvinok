//-*-c++-*-
#ifndef DYNAMIC_DEBUG_H
#define DYNAMIC_DEBUG_H

#include <fstream>
#include <map>
#include <algorithm>
#include <string>

using namespace std;

class DynamicDebug {
  map<string,long> debug_functions;
 public:
  DynamicDebug(const char* DebugConfigFile);
  DynamicDebug(ifstream& DebugConfigStream);
#ifdef DEBUG_RETURNS_FALSE
  inline long debug(const string& s) const { return 0; }
#else
  long debug(const string& s) const { 
    if (debug_functions.find(s) == debug_functions.end())
      return 0;
    else
      return (*(debug_functions.find(s))).second;
  }
#endif
  inline void add(const string& s, const long level) { debug_functions[s] = level; }
  inline void remove(const string& s) { debug_functions.erase(s); }
  map<string,long>::const_iterator begin() const { return debug_functions.begin(); }
  map<string,long>::const_iterator end() const { return debug_functions.end(); }
};

extern DynamicDebug DebugBlackboard;

#endif
