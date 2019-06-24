#ifndef BIN_FIELD_H
#define BIN_FIELD_H

#include "helpinghand.h"

unsigned CheckXDimensionOfMatrixFromFile(const std::string &filename, const std::string &prefix_if);
unsigned CheckYDimensionOfMatrixFromFile(const std::string &filename, const std::string &prefix_if);

template < typename valuetype >
class BinField {
  
 public:
  BinField(const unsigned &Nx, const unsigned &Ny, const valuetype &value);
  BinField(const unsigned &N, const valuetype &value);
  BinField(const std::string &filename, const std::string &prefix_if);
  
  BinField<valuetype> operator - (const BinField<valuetype> &Summand);
  void operator += (const BinField<valuetype> &Summand);
  void operator *= (const double &factor);
  
  void out(const double &xbinlength_, const double &ybinlength_,
	   const BinField< double > &xlow_, const BinField< double > &ylow_,
	   const std::string &filename, const std::string &prefix_of) const;
  void MatrixFromFile(const std::string &filename, const std::string &prefix_if);
  void MatrixToFile(const std::string &filename, const std::string &prefix_of) const;
  void fout(const std::string &filename, const std::string &prefix_of) const;
  void xout() const;
  void xout(const unsigned &xlow, const unsigned &xup, const unsigned &ylow, const unsigned &yup) const;
  
  BinField<valuetype> subgrid (const unsigned &start_x, const unsigned &end_x, const unsigned &start_y, const unsigned &end_y);
  
  unsigned call_Nx() const;
  unsigned call_Ny() const;
  valuetype call(const unsigned &xi, const unsigned &yi) const;
  valuetype call(const unsigned &ni) const;
  
  void assign(const unsigned &xi, const unsigned &yi, const valuetype &value);
  void assign(const unsigned &ni, const valuetype &value);
  void add(const unsigned &xi, const unsigned &yi, const valuetype &value);
  void add(const unsigned &ni, const valuetype &value);
  void scale(const unsigned &xi, const unsigned &yi, const valuetype &value);
  void scale(const unsigned &ni, const valuetype &value);
  
  unsigned xi(const unsigned &ni);
  unsigned yi(const unsigned &ni);
  
 private:
  const unsigned Nx_;
  const unsigned Ny_;
  std::vector< std::vector< valuetype > > values_;
  
  void MatrixFromFileError(const std::string &filename, const unsigned &xi, const unsigned &yi,
			   const std::string &prefix_if) const;
  void MatrixFromFileWarning(const std::string &filename, const std::string &prefix_if) const;
};

template < typename valuetype >
BinField<valuetype>::BinField(const unsigned &Nx,
			      const unsigned &Ny, 
			      const valuetype &value) :
Nx_ ( Nx ), Ny_ ( Ny ), values_ ( std::vector< std::vector< valuetype > > (Nx_, std::vector< valuetype > (Ny_, value) ) )
{}

template < typename valuetype >
BinField<valuetype>::BinField(const unsigned &N, 
			      const valuetype &value) :
Nx_ ( N ), Ny_ ( N ), values_ ( std::vector< std::vector< valuetype > > (Nx_, std::vector< valuetype > (Ny_, value) ) )
{}

template < typename valuetype >
BinField<valuetype>::BinField(const std::string &filename, 
			      const std::string &prefix_if) :
Nx_ ( CheckXDimensionOfMatrixFromFile(filename,prefix_if) ), Ny_ ( CheckYDimensionOfMatrixFromFile(filename,prefix_if) ), values_ ( std::vector< std::vector< valuetype > > (Nx_, std::vector< valuetype > (Ny_, 0) ) )
{
  std::ifstream FromFile( (prefix_if + filename).c_str() );
  if(FromFile.fail()){
    std::cerr << "ERROR: ifstream failed to read the matrix from " << prefix_if << filename << ";" << std::endl;
    exit(-1);
  }
  std::string dummy;
  char outlook;
  
  for(unsigned yi = 0; yi < Ny_; yi++){
    for(int xi = 0; xi < int(Nx_); xi++){
      if(FromFile.eof())
	MatrixFromFileError(filename,xi,yi,prefix_if);
      FromFile >> dummy;
      outlook = *(dummy.begin());
      if ( (outlook >= '0') && (outlook <= '9') ){
	std::stringstream dump(dummy);
	dump >> (values_.at(xi)).at(yi);
      }
      else{
	xi--;
      }
    }
  }
  
  FromFile.close();
}

template < typename valuetype >
BinField<valuetype> BinField<valuetype>::operator- (const BinField<valuetype> &Summand)
{
  unsigned Nx = Summand.call_Nx(), Ny = Summand.call_Ny();
  if( !((Nx_ == Nx) && (Ny_ == Ny)) ){
    std::cerr << "ERROR: dimensions mismatch " << Nx_ << "x" << Ny_ << " BinField - " << Nx << "x" << Ny << " BinField;"
	      << std::endl;
    exit(-1);
  }
  BinField<valuetype> tmp(Nx_, Ny_,0);
  for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
    tmp.assign(ni,call(ni) - Summand.call(ni));
  
  return tmp;
}

template < typename valuetype >
void BinField<valuetype>::operator+= (const BinField<valuetype> &Summand)
{
  unsigned Nx = Summand.call_Nx(), Ny = Summand.call_Ny();
  if( !((Nx_ == Nx) && (Ny_ == Ny)) ){
    std::cerr << "ERROR: dimensions mismatch " << Nx_ << "x" << Ny_ << " BinField + " << Nx << "x" << Ny << " BinField;"
	      << std::endl;
    exit(-1);
  }
  for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
    (values_[ni/Ny_])[ni%Ny_] += Summand.call(ni);
}

template < typename valuetype >
void BinField<valuetype>::operator*= (const double &factor)
{
  for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
    (values_[ni/Ny_])[ni%Ny_] *= factor;
}


template < typename valuetype >
BinField<valuetype> BinField<valuetype>::subgrid (const unsigned &start_x, const unsigned &end_x, const unsigned &start_y, const unsigned &end_y)
{
  // includes start and end, i.e., starts from start_x, goes up to end_x (same for y)
  // counting starts at 0
  
  BinField<valuetype> tmp(end_x-start_x+1,end_y-start_y+1,0);
  for(unsigned xi = 0; xi < end_x-start_x+1; xi++)
    for(unsigned yi = 0; yi < end_y-start_y+1; yi++){
      tmp.assign(xi,yi,call(xi+start_x,yi+start_y));
    }
  
  return tmp;
}

template < typename valuetype >
void BinField<valuetype>::out(const double &xbinlength_, const double &ybinlength_,
			      const BinField< double > &xlow_, const BinField< double > &ylow_,
			      const std::string &filename, const std::string &prefix_of) const
{
  double min_D = call(0,0);
  double max_D = min_D;
  for(unsigned xi = 0; xi < Nx_; xi++)
    for(unsigned yi = 0; yi < Ny_; yi++){
      double D = call(xi,yi);
      if ( D < min_D )
	min_D = D;
      else if ( D > max_D )
	max_D = D;
    }
  
  std::stringstream filename_stst;
  filename_stst << prefix_of << filename
		<< "_" << Nx_ << "x" << Ny_ << ".dat";
  std::string filename_st = filename_stst.str();
  std::ofstream Data( filename_st.c_str() );
  
  Data << "# Bin Field" << std::endl
       << "#" << std::setw(13) << "X" << std::setw(14) << "Y" << std::setw(14) << "Data" << std::endl;
  for(unsigned xi = 0; xi < Nx_; xi++)
    {
      for(unsigned yi = 0; yi < Ny_; yi++)
	Data << std::setw(14) << xlow_.call(xi,yi)
	     << std::setw(14) << ylow_.call(xi,yi)
	     << std::setw(14) << call(xi,yi)
	     << std::endl
	     << std::setw(14) << xlow_.call(xi,yi)
	     << std::setw(14) << ylow_.call(xi,yi) + ybinlength_
	     << std::setw(14) << call(xi,yi)
	     << std::endl;
      
      Data << std::endl;
      
      for(unsigned yi = 0; yi < Ny_; yi++)
	Data << std::setw(14) << xlow_.call(xi,yi) + xbinlength_
	     << std::setw(14) << ylow_.call(xi,yi)
	     << std::setw(14) << call(xi,yi)
	     << std::endl
	     << std::setw(14) << xlow_.call(xi,yi) + xbinlength_
	     << std::setw(14) << ylow_.call(xi,yi) + ybinlength_
	     << std::setw(14) << call(xi,yi)
	     << std::endl;
      
      Data << std::endl;
    }
  Data.close();
  
  std::stringstream gpfilename_stst;
  gpfilename_stst << prefix_of << filename
		  << "_" << Nx_ << "x" << Ny_ << ".gp";
  std::string gpfilename_st = gpfilename_stst.str();
  std::ofstream Gnuplot( gpfilename_st.c_str() );
  
  Gnuplot << "set terminal epslatex standalone color colortext 14" << std::endl
	  << "set output '" << filename << "_" << Nx_ << "x" << Ny_ << ".tex'"
	  << std::endl << std::endl
	  << "xmin = " << xlow_.call(0,0) << std::endl
	  << "xmax = " << xlow_.call(Nx_-1,0) + xbinlength_ << std::endl
	  << "ymin = " << ylow_.call(0,0) << std::endl
	  << "ymax = " << ylow_.call(0,Ny_-1) + ybinlength_  << std::endl
	  << std::endl
	  << "unset xtics # border mirror norotate offset 0,0.3" << std::endl
	  << "unset ytics # border mirror norotate offset 0.3,0" << std::endl
	  << "#set cbtics 10" << std::endl
	  << std::endl
	  << "set mxtics 2" << std::endl
	  << "set mytics 2" << std::endl
	  << std::endl;
  if(Ny_ == Nx_)
    Gnuplot << "set size ratio 1" << std::endl;
  else
    Gnuplot << "set size ratio " << Ny_*1.0/Nx_ << std::endl;
  Gnuplot << "set border lw 2" << std::endl
	  << std::endl
	  << "unset xlabel #'$x$' offset 0,0.7" << std::endl
	  << "unset ylabel #'$y$' offset 0.7,0" << std::endl
	  << "set cblabel '\\large $\\langle\\mathcal{D}\\rangle$' offset -0.7, 0" << std::endl
	  << std::endl
	  << "#set cbtics 2 offset -1.1,0" << std::endl
	  << std::endl
	  << "set xrange [xmin:xmax]" << std::endl
	  << "set yrange [ymin:ymax]" << std::endl
	  << "set cbrange [" << min_D << ":" << max_D << "]" << std::endl
	  << std::endl
	  << "set pm3d map" << std::endl
	  << "# set palette defined (" << min_D << " 'black', " << 0.5*min_D << " 'dark-blue', " << 0 << " 'white', " << 0.2*max_D << " 'purple', " << 0.4*max_D << " 'dark-red', " << 0.6*max_D << " 'red', " << 0.8*max_D << " 'orange', " << max_D << " 'yellow')" << std::endl
	  << "splot \"" << filename << "_" << Nx_ << "x" << Ny_ << ".dat" << "\" notitle" << std::endl
	  << "#" << std::endl
	  << std::endl;
  
  Gnuplot.close();   
}

template < typename valuetype >
void BinField<valuetype>::MatrixFromFile(const std::string &filename, const std::string &prefix_if)
{
  std::ifstream FromFile( (prefix_if + filename).c_str() );
  if(FromFile.fail()){
    std::cerr << "ERROR: ifstream failed to read the matrix from " << prefix_if << filename << ";" << std::endl;
    exit(-1);
  }
  std::string dummy;
  char outlook;
  
  for(unsigned yi = 0; yi < Ny_; yi++){
    for(int xi = 0; xi < int(Nx_); xi++){
      if(FromFile.eof())
	MatrixFromFileError(filename,xi,yi,prefix_if);
      FromFile >> dummy;
      outlook = *(dummy.begin());
      if ( (outlook >= '0') && (outlook <= '9') ){
	std::stringstream dump(dummy);
	dump >> (values_.at(xi)).at(yi);
      }
      else{
	xi--;
      }
    }
  }
  
  int endline = 0;
  if(!FromFile.eof()){
    FromFile >> endline;
    if(endline != 0 && !FromFile.eof())
      MatrixFromFileWarning(filename,prefix_if);
  }
  
  FromFile.close();
}

template < typename valuetype >
void BinField<valuetype>::MatrixToFile(const std::string &filename, const std::string &prefix_of) const
{
  std::ofstream Outfile( (prefix_of + filename).c_str() );
  for(unsigned yi = 0; yi < Ny_; yi++)
    {
      for(unsigned xi = 0; xi < Nx_; xi++)
	Outfile << std::setw(15) << call(xi,yi);
      Outfile << std::endl;
    }
  Outfile.close();
}

template < typename valuetype >
void BinField<valuetype>::fout(const std::string &filename, const std::string &prefix_of) const
{
  std::ofstream OutFile( (prefix_of + filename).c_str() );
  for(unsigned yi = Ny_; yi > 0; yi--)
    {
      for(unsigned xi = 0; xi < Nx_; xi++)
	OutFile << std::setw(15) << call(xi,yi-1);
      OutFile << std::endl;
    }
  OutFile.close();
}

template < typename valuetype >
void BinField<valuetype>::xout() const
{
  xout(0,Nx_-1,0,Ny_-1);
}

template < typename valuetype >
void BinField<valuetype>::xout(const unsigned &xlow, const unsigned &xup, const unsigned &ylow, const unsigned &yup) const
{
  std::cout << "Values of the BinField arranged following cartesian coordinates:" << std::endl; 
  std::cout << std::endl;
  for(unsigned yi = yup+1; yi > ylow; yi--)
    {
      for(unsigned xi = xlow; xi <= xup; xi++)
	std::cout << std::setw(15) << call(xi,yi-1);
      std::cout << std::endl;
    }
  std::cout << std::endl;
}

template < typename valuetype >
unsigned BinField<valuetype>::call_Nx() const
{
  return Nx_;
}

template < typename valuetype >
unsigned BinField<valuetype>::call_Ny() const
{
  return Ny_;
}

template < typename valuetype >
valuetype BinField<valuetype>::call(const unsigned &xi, const unsigned &yi) const
{
  /*
  if(xi >= Nx_)
    {
      std::cerr << "ERROR: BinField recieved call at x = " << xi << std::endl
		<< "       but maximum = " << (Nx_-1) << std::endl;
      exit(-1);
    }
  if(yi >= Ny_)
    {
      std::cerr << "ERROR: BinField recieved call at y = " << yi << std::endl
		<< "       but maximum = " << (Ny_-1) << std::endl;
      exit(-1);
    } 
  */
  //return (values_.at(xi)).at(yi);
  return (values_[xi])[yi];
}

template < typename valuetype >
valuetype BinField<valuetype>::call(const unsigned &ni) const
{
  /*
  if(ni >= Nx_*Ny_)
    {
      std::cerr << "ERROR: BinField recieved call at n = " << ni << std::endl
		<< "       but maximum = " << (Nx_*Ny_-1) << std::endl;
      exit(-1);
    }
  // In one loop through the bin-field:
  // one x-column after the other from bottom to top
  // x = 0..Nx_-1 and for each x:
  // y = 0..Ny_-1
  */
  return call(ni/Ny_,ni%Ny_);
}

template < typename valuetype >
void BinField<valuetype>::assign(const unsigned &xi, const unsigned &yi, const valuetype &value)
{
  /*
  if(xi >= Nx_)
    {
      std::cerr << "ERROR: BinField recieved call at x = " << xi << std::endl
		<< "       but maximum = " << (Nx_-1) << std::endl;
      exit(-1);
    }
  if(yi >= Ny_)
    {
      std::cerr << "ERROR: BinField recieved call at y = " << yi << std::endl
		<< "       but maximum = " << (Ny_-1) << std::endl;
      exit(-1);
    } 
  */
  (values_[xi])[yi] = value;
}

template < typename valuetype >
void BinField<valuetype>::assign(const unsigned &ni, const valuetype &value)
{
  /*
  if(ni >= Nx_*Ny_)
    {
      std::cerr << "ERROR: BinField recieved call at n = " << ni << std::endl
		<< "       but maximum = " << (Nx_*Ny_-1) << std::endl;
      exit(-1);
    }
  */
  
  // In one loop through the bin-field:
  // one x-column after the other from bottom to top
  // x = 0..Nx_-1 and for each x:
  // y = 0..Ny_-1
  assign(ni/Ny_,ni%Ny_,value);
}

template < typename valuetype >
void BinField<valuetype>::add(const unsigned &xi, const unsigned &yi, const valuetype &value)
{
  /*
  if(xi >= Nx_)
    {
      std::cerr << "ERROR: BinField recieved call at x = " << xi << std::endl
		<< "       but maximum = " << (Nx_-1) << std::endl;
      exit(-1);
    }
  if(yi >= Ny_)
    {
      std::cerr << "ERROR: BinField recieved call at y = " << yi << std::endl
		<< "       but maximum = " << (Ny_-1) << std::endl;
      exit(-1);
    } 
  */
  (values_[xi])[yi] += value;
}

template < typename valuetype >
void BinField<valuetype>::add(const unsigned &ni, const valuetype &value)
{
  /*
  if(ni >= Nx_*Ny_)
    {
      std::cerr << "ERROR: BinField recieved call at n = " << ni << std::endl
		<< "       but maximum = " << (Nx_*Ny_-1) << std::endl;
      exit(-1);
    }
  */
  (values_[ni/Ny_])[ni%Ny_] += value;
}

template < typename valuetype >
void BinField<valuetype>::scale(const unsigned &xi, const unsigned &yi, const valuetype &value)
{
  /*
  if(xi >= Nx_)
    {
      std::cerr << "ERROR: BinField recieved call at x = " << xi << std::endl
		<< "       but maximum = " << (Nx_-1) << std::endl;
      exit(-1);
    }
  if(yi >= Ny_)
    {
      std::cerr << "ERROR: BinField recieved call at y = " << yi << std::endl
		<< "       but maximum = " << (Ny_-1) << std::endl;
      exit(-1);
    } 
  */
  (values_[xi])[yi] *= value;
}

template < typename valuetype >
void BinField<valuetype>::scale(const unsigned &ni, const valuetype &value)
{
  /*
  if(ni >= Nx_*Ny_)
    {
      std::cerr << "ERROR: BinField recieved call at n = " << ni << std::endl
		<< "       but maximum = " << (Nx_*Ny_-1) << std::endl;
      exit(-1);
    }
  */
  (values_[ni/Ny_])[ni%Ny_] *= value;
}

template < typename valuetype >
unsigned BinField<valuetype>::xi(const unsigned &ni){
  return ni/Ny_;
}

template < typename valuetype >
unsigned BinField<valuetype>::yi(const unsigned &ni){
  return ni - (ni/Ny_)*Ny_;
}

template < typename valuetype >
void BinField<valuetype>::MatrixFromFileError(const std::string &filename, const unsigned &xi, const unsigned &yi,
			   const std::string &prefix_if) const
{
  std::cerr << "ERROR: MatrixFromFile recieved matrix from " << prefix_if << filename << ";" << std::endl
	    << "       Dimension mismatch!" << std::endl
	    << "       Needed " << Nx_ << " x " << Ny_ << " matrix;" << std::endl
	    << "       Failed at position xi = " << xi << ", yi = " << yi << ";" << std::endl;
  exit(-1);
}

template < typename valuetype >
void BinField<valuetype>::MatrixFromFileWarning(const std::string &filename, const std::string &prefix_if) const
{
  std::cerr << "WARNING: MatrixFromFile recieved matrix from " << prefix_if << filename << ";" << std::endl
	    << "         There might be data loss!" << std::endl;
}

#endif

