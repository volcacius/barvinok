#include "debug.h"
#include <fstream>
#include <cctype>
#include <iostream>
#include <deque>
#if 0
#include "polyast.h"
#include "constraint.h"
#endif

using namespace std;

#define MAX_LINE_LENGTH 1024

DynamicDebug::DynamicDebug(ifstream& DebugConfigStream)
{
  char line[MAX_LINE_LENGTH];
  while (DebugConfigStream.getline(line, MAX_LINE_LENGTH)) {
    // parse line. Possible formats:
    // DEBUG_PARAMETER
    // DEBUG_PARAMETER level
    string word;
    long level;
    int i;
    for(i=0; !isspace(line[i]) && line[i]!=0; i++)
      word += line[i];
    while (isspace(line[i] && line[i]!=0)) i++;
    if (line[i]!=0) {
      // start parsing level;
      char *p;
      level = strtol(line+i, &p, 10);
      while(isspace(*p) && *p!=0) p++;
      if (*p!=0) {
	cerr<<"Error while parsing Debug Parameters!!\n"
	    <<"for parameter '"<<word<<"', "
	    <<" found level "<<level<<", but level was "
	    <<"not last non-white on line."<<endl;
      }
    } else {
      level=1;
    }
    add(word, level);
  }
}

inline DynamicDebug::DynamicDebug(const char* DebugConfigFile)
{
  ifstream DebugConfigStream(DebugConfigFile);
  char line[MAX_LINE_LENGTH];
  while (DebugConfigStream.getline(line, MAX_LINE_LENGTH)) {
    // parse line. Possible formats:
    // DEBUG_PARAMETER
    // DEBUG_PARAMETER level
    string word;
    long level;
    int i;
    for(i=0; !isspace(line[i]) && line[i]!=0; i++)
      word += line[i];
    while (isspace(line[i] && line[i]!=0)) i++;
    if (line[i]!=0) {
      // start parsing level;
      char *p;
      level = strtol(line+i, &p, 10);
      while(isspace(*p) && *p!=0) p++;
      if (*p!=0) {
	cerr<<"Error while parsing Debug Parameters!!\n"
	    <<"for parameter '"<<word<<"', "
	    <<" found level "<<level<<", but level was "
	    <<"not last non-white on line."<<endl;
      }
    } else {
      level=1;
    }
    add(word, level);
  }  
}

//ifstream DEBUG_CONFIGURATION_FILE("DEBUG_POLYAST");
DynamicDebug DebugBlackboard("DEBUG_POLYAST");

deque<string> split_underscore(const string& s)
{
  deque<string> result;
  const char* word_begin = s.c_str();
  const char* c=s.c_str();
  for(; (*c)!=0; c++) {
    if ((*c)=='_') {
      result.push_back(string(word_begin, c-word_begin));
      word_begin = c+1;
    }
  }
  if (word_begin<c)
    result.push_back(string(word_begin, c-word_begin));
  cout<<"in split_underscore("<<s<<"): returned (";
  for (unsigned i=0;i<result.size();i++)
    cout<<((i!=0)?",":"")<<result[i];
  cout<<")"<<endl;
  return result;
}


