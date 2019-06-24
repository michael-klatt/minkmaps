#include "helpinghand.h"

bool Cout_One_Percent(const unsigned int &i, unsigned int &p, const unsigned int &N){
  
  if(N < 100){
    unsigned M = static_cast<unsigned>(100.0/N);
    for(unsigned j = 0; j < M; j++)
      std::cout << "|" << std::flush;
    if(i >= p/(100.0/N-static_cast<unsigned>(100.0/N))){
      p++;
      std::cout << "|" << std::flush;
      if(i == (N-1))
	std::cout << std::endl << std::flush;
      return true;
    }
    else{
      if(i == (N-1))
	std::cout << std::endl << std::flush;
      return false;
    }
    
  }
  else{
    if(i >= p*N/100.0){
      p++;
      std::cout << "|" << std::flush;
      if(i == (N-1))
	std::cout << std::endl << std::flush;
      return true;
    }
    else{
      if(i == (N-1))
	std::cout << std::endl << std::flush;
      return false;
    }
    
  }
  
}

std::vector<double> stdvector16(const double &val00, const double &val01, const double &val02, const double &val03, 
				const double &val04, const double &val05, const double &val06, const double &val07, 
				const double &val08, const double &val09, const double &val10, const double &val11, 
				const double &val12, const double &val13, const double &val14, const double &val15)
{
  std::vector<double> vec(16,0);
  vec.at(0)  = val00; vec.at(1)  = val01; vec.at(2)  = val02; vec.at(3)  = val03;
  vec.at(4)  = val04; vec.at(5)  = val05; vec.at(6)  = val06; vec.at(7)  = val07;
  vec.at(8)  = val08; vec.at(9)  = val09; vec.at(10) = val10; vec.at(11) = val11;
  vec.at(12) = val12; vec.at(13) = val13; vec.at(14) = val14; vec.at(15) = val15;
  return vec;
}

std::vector<int> stdvector_int16(const int &val00, const int &val01, const int &val02, const int &val03, 
			     const int &val04, const int &val05, const int &val06, const int &val07, 
			     const int &val08, const int &val09, const int &val10, const int &val11, 
			     const int &val12, const int &val13, const int &val14, const int &val15)
{
  std::vector<int> vec(16,0);
  vec.at(0)  = val00; vec.at(1)  = val01; vec.at(2)  = val02; vec.at(3)  = val03;
  vec.at(4)  = val04; vec.at(5)  = val05; vec.at(6)  = val06; vec.at(7)  = val07;
  vec.at(8)  = val08; vec.at(9)  = val09; vec.at(10) = val10; vec.at(11) = val11;
  vec.at(12) = val12; vec.at(13) = val13; vec.at(14) = val14; vec.at(15) = val15;
  return vec;
}

std::vector<int> stdvector_int9(const int &val0, const int &val1, const int &val2, const int &val3, 
				const int &val4, const int &val5, const int &val6, const int &val7, 
				const int &val8)
{
  std::vector<int> vec(9,0);
  vec.at(0)  = val0; vec.at(1)  = val1; vec.at(2)  = val2; vec.at(3)  = val3;
  vec.at(4)  = val4; vec.at(5)  = val5; vec.at(6)  = val6; vec.at(7)  = val7;
  vec.at(8)  = val8;
  return vec;
}
