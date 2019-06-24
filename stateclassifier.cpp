#include "stateclassifier.h"

stateclassifier::stateclassifier(const int &A_arg, const int &P_arg, const int &C_arg,
				 const double &cl_arg) :
  A_ (A_arg), P_ (P_arg), C_ (C_arg), cl_ (cl_arg)
{}

void stateclassifier::operator= (const stateclassifier &equal)
{
  A_ = equal.A();
  P_ = equal.P();
  C_ = equal.C();
  cl_ = equal.cl();
}

bool stateclassifier::operator== (const stateclassifier &compare) const {
  return (A_ == compare.A()) && (P_ == compare.P()) && (C_ == compare.C());
}

bool stateclassifier::operator< (const stateclassifier &compare) const {
  return (cl_ < compare.cl());
}

void stateclassifier::add_cl(const double &delta) {
  cl_ += delta;
}

int stateclassifier::A() const {
  return A_;
}

int stateclassifier::P() const {
  return P_;
}

int stateclassifier::C() const {
  return C_;
}

double stateclassifier::cl() const {
  return cl_;
}
