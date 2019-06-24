#ifndef HELPING_HAND_H
#define HELPING_HAND_H

#include "initialize.h"

const std::string PercentageBar = "0_____________________0.25_____________________0.50_____________________0.75_______________________1";
bool Cout_One_Percent(const unsigned int &i, unsigned int &p, const unsigned int &N);

template < typename number >
number sum(const std::vector<number> &oughttobesummed)
{
  number result = 0;
  for(unsigned int i = 0; i < oughttobesummed.size(); i++)
    result += oughttobesummed.at(i);
  
  return result;  
};

template < typename number >
number sum_abs(const std::vector<number> &oughttobesummed)
{
  number result = 0;
  for(unsigned int i = 0; i < oughttobesummed.size(); i++)
    result += fabs(oughttobesummed.at(i));
  
  return result;  
};

template < typename number >
void print_elements(const std::vector<number> &oughttobeprinted, const unsigned int &precision = 10)
{
  for(unsigned int i = 0; i < oughttobeprinted.size(); i++)
    std::cout << std::setw(precision+5) << std::setprecision(precision) << oughttobeprinted.at(i); 
  std::cout << std::endl;
};

template < typename number >
void print_rows_and_elements(const std::vector<number> &oughttobeprinted, const unsigned int &precision = 10)
{
  for(unsigned int i = 0; i < oughttobeprinted.size(); i++)
    std::cout << std::setw(5) << i+1 << " " << std::setw(precision+5)
	      << std::setprecision(precision) << oughttobeprinted.at(i) << std::endl;
};

template < typename number >
std::vector<number> stdvector(const unsigned int &quantity, ... )
{
  std::vector<number> vec;
  vec.clear();
  
  va_list variablelist;
  va_start(variablelist,quantity);
  
  for(unsigned int i=0; i < quantity; i++)
    vec.push_back(va_arg(variablelist,number));
  
  va_end(variablelist);
  
  if(vec.size() != quantity){
    std::cerr << "ERROR: stdvector recieved " << vec.size() << " number of arguments "
	      << "instead of " << quantity << " as demanded;"
	      << std::endl;
    exit(-1);
  }
  return vec;
};

std::vector<double> stdvector16(const double &val00, const double &val01, const double &val02, const double &val03, 
				const double &val04, const double &val05, const double &val06, const double &val07, 
				const double &val08, const double &val09, const double &val10, const double &val11, 
				const double &val12, const double &val13, const double &val14, const double &val15);

std::vector<int> stdvector_int16(const int &val00, const int &val01, const int &val02, const int &val03, 
				 const int &val04, const int &val05, const int &val06, const int &val07, 
				 const int &val08, const int &val09, const int &val10, const int &val11, 
				 const int &val12, const int &val13, const int &val14, const int &val15);

std::vector<int> stdvector_int9(const int &val0, const int &val1, const int &val2, const int &val3, 
				const int &val4, const int &val5, const int &val6, const int &val7, 
				const int &val8);

#endif

