#ifndef MACROSTATE_H
#define MACROSTATE_H

#include "helpinghand.h"

class macrostate{
  
 public:
  macrostate(const int &P, const int &C);
  
  void operator = (const macrostate &equal);
  
  bool operator == (const macrostate &compare) const;
  bool operator < (const macrostate &compare) const;
  
  void add(const int &dP, const int &dC);
  void add_P(const int &dP);
  void add_C(const int &dC);
  
  int at(const unsigned &i) const;
  
 private:
  std::vector<int> state;
};

#endif

