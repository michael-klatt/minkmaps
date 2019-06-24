#include "binfield.h"

unsigned CheckXDimensionOfMatrixFromFile(const std::string &filename, const std::string &prefix_if)
{
  std::ifstream FromFile( (prefix_if + filename).c_str() );
  if(FromFile.fail()){
    std::cerr << "ERROR: ifstream failed to read the matrix from " << prefix_if << filename
	      << " for checking the x-dimension;" << std::endl;
    exit(-1);
  }
  unsigned nmbrofcols = 0;
  std::string line;
  bool comment = true;
  bool nmbr = false;
  
  while(comment && getline(FromFile, line)){
    comment = false;
    
    for (std::string::iterator it=line.begin() ; it < line.end(); it++ ){
      if(*it == '#')
	comment = true;
      if(*it != '\t' && *it != ' '){
	if(!nmbr)
	  nmbrofcols++;
	nmbr = true;}
      if(*it == '\t' || *it == ' ')
	nmbr = false;      
    }
    
    if(comment || nmbrofcols == 0){
      nmbrofcols = 0;
      nmbr = false;
      comment = true;}
  }
  
  if(nmbrofcols == 0){
    std::cerr << "ERROR: ifstream failed to read the matrix from " << prefix_if << filename
	      << " for checking the x-dimension;" << std::endl
	      << "       There are only lines with comments: #" << std::endl;
    exit(-1);}
  
  return nmbrofcols;
}

unsigned CheckYDimensionOfMatrixFromFile(const std::string &filename, const std::string &prefix_if)
{
  std::ifstream FromFile( (prefix_if + filename).c_str() );
  if(FromFile.fail()){
    std::cerr << "ERROR: ifstream failed to read the matrix from " << prefix_if << filename
	      << " for checking the y-dimension;" << std::endl;
    exit(-1);
  }
  unsigned nmbrofrows = 0;
  std::string line;
  bool content = false;
  
  while(getline(FromFile, line)){
    // std::cout << line << std::endl;
    for (std::string::iterator it=line.begin() ; it < line.end(); it++ ){
      if(*it == '#')
	break;
      if(*it != '\t' && *it != ' '){
	content = true;
	break;}
    }
    if(content){
      // std::cout << " content = true :" << line << "!!!" << std::endl; 
      nmbrofrows++;
    }
    content = false;
  }
  
  return nmbrofrows;
}

