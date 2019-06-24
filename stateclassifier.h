#ifndef STATECLASSIFIER_H
#define STATECLASSIFIER_H

#include "helpinghand.h"

class stateclassifier{
  
 public:
  stateclassifier(const int &A_arg, const int &P_arg, const int &C_arg,
		  const double &cl_arg);
  
  void operator = (const stateclassifier &equal);
  
  bool operator == (const stateclassifier &compare) const;
  bool operator < (const stateclassifier &compare) const;
  
  void add_cl(const double &delta);
  
  int A() const;
  int P() const;
  int C() const;
  double cl() const;
  
 private:
  int A_;
  int P_;
  int C_;
  double cl_; // classifier
};

#endif

