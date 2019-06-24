#include "macrostate.h"

macrostate::macrostate(const int &P, const int &C) :
  state (std::vector<int> (2,0))
{
  state.at(0) = P; state.at(1) = C;
}

void macrostate::operator= (const macrostate &equal)
{
  state.at(0) = equal.at(0);
  state.at(1) = equal.at(1);
}

bool macrostate::operator== (const macrostate &compare) const
{
  return (state.at(0) == compare.at(0)) && (state.at(1) == compare.at(1));
}

bool macrostate::operator< (const macrostate &compare) const
{
  if(state.at(1) < compare.at(1)){
    return true;
  }
  else if(state.at(1) == compare.at(1)){
    return (state.at(0) < compare.at(0));
  }
  else{
    return false;
  }
  
}

void macrostate::add(const int &dP, const int &dC)
{
  state.at(0) += dP;
  state.at(1) += dC;
}

void macrostate::add_P(const int &dP)
{
  state.at(0) += dP;
}

void macrostate::add_C(const int &dC)
{
  state.at(0) += dC;
}

int macrostate::at(const unsigned &i) const
{
  return state.at(i);
}
