#include "minkowski.h"

// Minkowski Sky Map Tools
void TurnConfNmbrIntoBinField(const unsigned &conf, BinField<bool> &Configuration,
			      const unsigned &slide, const unsigned &N_conf)
{
  if(conf >= N_conf){
    std::cerr << "ERROR: " << conf << " is greater or equal than " << N_conf
	      << " and thus not a configuration number;" << std::endl;
    exit(-1);
  }
  unsigned x = conf, y = conf, ni = 0;
  
  // Binary code
  while(x > 0){
    // x : 2 = y Rest z
    // Rest z = x - 2*y
    y = x/2;
    // Rest z = true/false Configuration
    Configuration.assign(slide-1-ni+(ni/slide)*slide,ni/slide,x-2*y);
    x = y;
    ni++;
  }
}

bool eq(const double &first, const double &second)
{
  return ( fabs(first - second) < pow(10,-7) );
}

// MINKOWSKI FUNCTIONALS: SUMMING LOOK-UP TABLE (MARCHING SQUARE ALGORITHM)
unsigned convert(const bool &first, const bool &second, const bool &third, const bool &fourth)
{
  return first + 2*second + 4*third + 8*fourth;
}

// Minus sampling boundary condition
double area_mbc(const BinField<bool> &sample,
		const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny)
{
  double total = 0;
  
  if(! ((xi>=0)&&(xi<(sample.call_Nx()-Nx+1))&&(yi>=0)&&(yi<(sample.call_Ny()-Ny+1))&&(Nx>2)&&(Ny>2)) ){
    std::cerr << "ERROR: cannot compute Minkowski functional value; recieved wrong dimensions:" << std::endl
	      << "       sample is a " << sample.call_Nx() << "x" << sample.call_Ny() << "BinField, but" << std::endl
	      << "       SubBinField is " << Nx << "x" << Ny << " at coordinate (" << xi << "," << yi << ");" << std::endl
	      << std::endl;
    exit(-1);
  }
  
  unsigned X;
  unsigned Y;
  
  
  // -----------------------
  // BOTTOM ROW
  Y = 0;
  
  // Range 7
  X = 0;
  total += rg7_area.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
  // Range 8
  X++;
  while(X < (Nx-2))
    {
      total += rg8_area.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
      
      X++;
    }
  
  // Range 9
  total += rg9_area.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
  // -----------------------
  // MID ROWS
  Y++;
  while(Y < (Ny-2))
    {
  
      // Range 4
      X = 0;
      total += rg4_area.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
      // Range 5
      X++;
      while(X < (Nx-2))
	{
	  total += rg5_area.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
      
	  X++;
	}
  
      // Range 6
      total += rg6_area.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
      
      Y++;  
    }
  
  // -----------------------
  // TOP ROW
  
  // Range 1
  X = 0;
  total += rg1_area.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
  // Range 2
  X++;
  while(X < (Nx-2))
    {
      total += rg2_area.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
      
      X++;
    }
  
  // Range 3
  total += rg3_area.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
  return total;
}

double perimeter_mbc(const BinField<bool> &sample,
		     const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny)
{
  double total = 0;
  
  if(! ((xi>=0)&&(xi<(sample.call_Nx()-Nx+1))&&(yi>=0)&&(yi<(sample.call_Ny()-Ny+1))&&(Nx>2)&&(Ny>2)) ){
    std::cerr << "ERROR: cannot compute Minkowski functional value; recieved wrong dimensions:" << std::endl
	      << "       sample is a " << sample.call_Nx() << "x" << sample.call_Ny() << "BinField, but" << std::endl
	      << "       SubBinField is " << Nx << "x" << Ny << " at coordinate (" << xi << "," << yi << ");" << std::endl
	      << std::endl;
    exit(-1);
  }
  
  unsigned X;
  unsigned Y;
  
  
  // -----------------------
  // BOTTOM ROW
  Y = 0;
  
  // Range 7
  X = 0;
  total += rg7_perimeter.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
  // Range 8
  X++;
  while(X < (Nx-2))
    {
      total += rg8_perimeter.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
      
      X++;
    }
  
  // Range 9
  total += rg9_perimeter.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
  // -----------------------
  // MID ROWS
  Y++;
  while(Y < (Ny-2))
    {
  
      // Range 4
      X = 0;
      total += rg4_perimeter.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
      // Range 5
      X++;
      while(X < (Nx-2))
	{
	  total += rg5_perimeter.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
      
	  X++;
	}
  
      // Range 6
      total += rg6_perimeter.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
      
      Y++;  
    }
  
  // -----------------------
  // TOP ROW
  
  // Range 1
  X = 0;
  total += rg1_perimeter.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
  // Range 2
  X++;
  while(X < (Nx-2))
    {
      total += rg2_perimeter.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
      
      X++;
    }
  
  // Range 3
  total += rg3_perimeter.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
  return total;
}

double euler_mbc(const BinField<bool> &sample,
		 const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny)
{
  double total = 0;
  
  if(! ((xi>=0)&&(xi<(sample.call_Nx()-Nx+1))&&(yi>=0)&&(yi<(sample.call_Ny()-Ny+1))&&(Nx>2)&&(Ny>2)) ){
    std::cerr << "ERROR: cannot compute Minkowski functional value; recieved wrong dimensions:" << std::endl
	      << "       sample is a " << sample.call_Nx() << "x" << sample.call_Ny() << "BinField, but" << std::endl
	      << "       SubBinField is " << Nx << "x" << Ny << " at coordinate (" << xi << "," << yi << ");" << std::endl
	      << std::endl;
    exit(-1);
  }
  
  unsigned X;
  unsigned Y;
  
  
  // -----------------------
  // BOTTOM ROW
  Y = 0;
  
  // Range 7
  X = 0;
  total += rg7_euler.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
  // Range 8
  X++;
  while(X < (Nx-2))
    {
      total += rg8_euler.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
      
      X++;
    }
  
  // Range 9
  total += rg9_euler.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
  // -----------------------
  // MID ROWS
  Y++;
  while(Y < (Ny-2))
    {
  
      // Range 4
      X = 0;
      total += rg4_euler.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
      // Range 5
      X++;
      while(X < (Nx-2))
	{
	  total += rg5_euler.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
      
	  X++;
	}
  
      // Range 6
      total += rg6_euler.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
      
      Y++;  
    }
  
  // -----------------------
  // TOP ROW
  
  // Range 1
  X = 0;
  total += rg1_euler.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
  // Range 2
  X++;
  while(X < (Nx-2))
    {
      total += rg2_euler.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
      
      X++;
    }
  
  // Range 3
  total += rg3_euler.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
  return total;
}

double w102_xx_mbc(const BinField<bool> &sample,
		   const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny)
{
  double total = 0;
  
  if(! ((xi>=0)&&(xi<(sample.call_Nx()-Nx+1))&&(yi>=0)&&(yi<(sample.call_Ny()-Ny+1))&&(Nx>2)&&(Ny>2)) ){
    std::cerr << "ERROR: cannot compute Minkowski functional value; recieved wrong dimensions:" << std::endl
	      << "       sample is a " << sample.call_Nx() << "x" << sample.call_Ny() << "BinField, but" << std::endl
	      << "       SubBinField is " << Nx << "x" << Ny << " at coordinate (" << xi << "," << yi << ");" << std::endl
	      << std::endl;
    exit(-1);
  }
  
  unsigned X;
  unsigned Y;
  
  
  // -----------------------
  // BOTTOM ROW
  Y = 0;
  
  // Range 7
  X = 0;
  total += rg7_w102_xx.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
  // Range 8
  X++;
  while(X < (Nx-2))
    {
      total += rg8_w102_xx.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
      
      X++;
    }
  
  // Range 9
  total += rg9_w102_xx.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
  // -----------------------
  // MID ROWS
  Y++;
  while(Y < (Ny-2))
    {
  
      // Range 4
      X = 0;
      total += rg4_w102_xx.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
      // Range 5
      X++;
      while(X < (Nx-2))
	{
	  total += rg5_w102_xx.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
      
	  X++;
	}
  
      // Range 6
      total += rg6_w102_xx.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
      
      Y++;  
    }
  
  // -----------------------
  // TOP ROW
  
  // Range 1
  X = 0;
  total += rg1_w102_xx.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
  // Range 2
  X++;
  while(X < (Nx-2))
    {
      total += rg2_w102_xx.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
      
      X++;
    }
  
  // Range 3
  total += rg3_w102_xx.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
  return total;
}

double w102_xy_mbc(const BinField<bool> &sample,
		   const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny)
{
  double total = 0;
  
  if(! ((xi>=0)&&(xi<(sample.call_Nx()-Nx+1))&&(yi>=0)&&(yi<(sample.call_Ny()-Ny+1))&&(Nx>2)&&(Ny>2)) ){
    std::cerr << "ERROR: cannot compute Minkowski functional value; recieved wrong dimensions:" << std::endl
	      << "       sample is a " << sample.call_Nx() << "x" << sample.call_Ny() << "BinField, but" << std::endl
	      << "       SubBinField is " << Nx << "x" << Ny << " at coordinate (" << xi << "," << yi << ");" << std::endl
	      << std::endl;
    exit(-1);
  }
  
  unsigned X;
  unsigned Y;
  
  
  // -----------------------
  // BOTTOM ROW
  Y = 0;
  
  // Range 7
  X = 0;
  total += rg7_w102_xy.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
  // Range 8
  X++;
  while(X < (Nx-2))
    {
      total += rg8_w102_xy.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
      
      X++;
    }
  
  // Range 9
  total += rg9_w102_xy.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
  // -----------------------
  // MID ROWS
  Y++;
  while(Y < (Ny-2))
    {
  
      // Range 4
      X = 0;
      total += rg4_w102_xy.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
      // Range 5
      X++;
      while(X < (Nx-2))
	{
	  total += rg5_w102_xy.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
      
	  X++;
	}
  
      // Range 6
      total += rg6_w102_xy.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
      
      Y++;  
    }
  
  // -----------------------
  // TOP ROW
  
  // Range 1
  X = 0;
  total += rg1_w102_xy.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
  // Range 2
  X++;
  while(X < (Nx-2))
    {
      total += rg2_w102_xy.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
      
      X++;
    }
  
  // Range 3
  total += rg3_w102_xy.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
  return total;
}

double w102_yy_mbc(const BinField<bool> &sample,
		   const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny)
{
  double total = 0;
  
  if(! ((xi>=0)&&(xi<(sample.call_Nx()-Nx+1))&&(yi>=0)&&(yi<(sample.call_Ny()-Ny+1))&&(Nx>2)&&(Ny>2)) ){
    std::cerr << "ERROR: cannot compute Minkowski functional value; recieved wrong dimensions:" << std::endl
	      << "       sample is a " << sample.call_Nx() << "x" << sample.call_Ny() << "BinField, but" << std::endl
	      << "       SubBinField is " << Nx << "x" << Ny << " at coordinate (" << xi << "," << yi << ");" << std::endl
	      << std::endl;
    exit(-1);
  }
  
  unsigned X;
  unsigned Y;
  
  
  // -----------------------
  // BOTTOM ROW
  Y = 0;
  
  // Range 7
  X = 0;
  total += rg7_w102_yy.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
  // Range 8
  X++;
  while(X < (Nx-2))
    {
      total += rg8_w102_yy.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
      
      X++;
    }
  
  // Range 9
  total += rg9_w102_yy.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
  // -----------------------
  // MID ROWS
  Y++;
  while(Y < (Ny-2))
    {
  
      // Range 4
      X = 0;
      total += rg4_w102_yy.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
      // Range 5
      X++;
      while(X < (Nx-2))
	{
	  total += rg5_w102_yy.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
      
	  X++;
	}
  
      // Range 6
      total += rg6_w102_yy.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
      
      Y++;  
    }
  
  // -----------------------
  // TOP ROW
  
  // Range 1
  X = 0;
  total += rg1_w102_yy.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
  // Range 2
  X++;
  while(X < (Nx-2))
    {
      total += rg2_w102_yy.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
      
      X++;
    }
  
  // Range 3
  total += rg3_w102_yy.at(convert(sample.call(xi+X+1,yi+Y),sample.call(xi+X,yi+Y),sample.call(xi+X+1,yi+Y+1),sample.call(xi+X,yi+Y+1)));
  
  return total;
}

double delta_mbc(const BinField<bool> &sample,
		 const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny)
{
  double w102_xx_ = w102_xx_mbc(sample,xi,yi,Nx,Ny);
  double w102_xy_ = w102_xy_mbc(sample,xi,yi,Nx,Ny);
  double w102_yy_ = w102_yy_mbc(sample,xi,yi,Nx,Ny);
  
  return sqrt((w102_xx_-w102_yy_)*(w102_xx_-w102_yy_)+4*w102_xy_*w102_xy_);
}

// Periodic boundary condition
double area_pbc(const BinField<bool> &sample,
		const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny)
{
  double total = 0;
  
  if(! ((xi>=0)&&(xi<(sample.call_Nx()-Nx+1))&&(yi>=0)&&(yi<(sample.call_Ny()-Ny+1))&&(Nx>2)&&(Ny>2)) ){
    std::cerr << "ERROR: cannot compute Minkowski functional value; recieved wrong dimensions:" << std::endl
	      << "       sample is a " << sample.call_Nx() << "x" << sample.call_Ny() << "BinField, but" << std::endl
	      << "       SubBinField is " << Nx << "x" << Ny << " at coordinate (" << xi << "," << yi << ");" << std::endl
	      << std::endl;
    exit(-1);
  }
  
  unsigned X;
  unsigned Y;
  
  
  // -----------------------
  // BOTTOM ROW
  Y = 0;
  
  // Range 7
  X = 0;
  total += rg7_area.at(convert(sample.call(xi+X,yi+Ny-1),sample.call(xi+Nx-1,yi+Ny-1),sample.call(xi+X,yi+Y),sample.call(xi+Nx-1,yi+Y)));
  
  // Range 8
  X++;
  while(X < Nx)
    {
      total += rg8_area.at(convert(sample.call(xi+X,yi+Ny-1),sample.call(xi+X-1,yi+Ny-1),sample.call(xi+X,yi+Y),sample.call(xi+X-1,yi+Y)));
      X++;
    }
  
  // Range 9
  total += rg9_area.at(convert(sample.call(xi,yi+Ny-1),sample.call(xi+Nx-1,yi+Ny-1),sample.call(xi,yi+Y),sample.call(xi+Nx-1,yi+Y)));
  
  // -----------------------
  // MID ROWS
  Y++;
  while(Y < Ny)
    {
      
      // Range 4
      X = 0;
      total += rg4_area.at(convert(sample.call(xi+X,yi+Y-1),sample.call(xi+Nx-1,yi+Y-1),sample.call(xi+X,yi+Y),sample.call(xi+Nx-1,yi+Y)));
      
      // Range 5
      X++;
      while(X < Nx)
	{
	  total += rg5_area.at(convert(sample.call(xi+X,yi+Y-1),sample.call(xi+X-1,yi+Y-1),sample.call(xi+X,yi+Y),sample.call(xi+X-1,yi+Y)));
      	  X++;
	}
      
      // Range 6
      total += rg6_area.at(convert(sample.call(xi,yi+Y-1),sample.call(xi+Nx-1,yi+Y-1),sample.call(xi,yi+Y),sample.call(xi+Nx-1,yi+Y)));
      
      
      Y++;  
    }
  
  // -----------------------
  // TOP ROW
  
  // Range 1
  X = 0;
  total += rg1_area.at(convert(sample.call(xi+X,yi+Y-1),sample.call(xi+Nx-1,yi+Y-1),sample.call(xi+X,yi),sample.call(xi+Nx-1,yi)));
  
  // Range 2
  X++;
  while(X < Nx)
    {
      total += rg2_area.at(convert(sample.call(xi+X,yi+Y-1),sample.call(xi+X-1,yi+Y-1),sample.call(xi+X,yi),sample.call(xi+X-1,yi)));
      X++;
    }
  
  // Range 3
  total += rg3_area.at(convert(sample.call(xi,yi+Y-1),sample.call(xi+Nx-1,yi+Y-1),sample.call(xi,yi),sample.call(xi+Nx-1,yi)));
  
  return total;
}

double perimeter_pbc(const BinField<bool> &sample,
		     const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny)
{
  double total = 0;
  
  if(! ((xi>=0)&&(xi<(sample.call_Nx()-Nx+1))&&(yi>=0)&&(yi<(sample.call_Ny()-Ny+1))&&(Nx>2)&&(Ny>2)) ){
    std::cerr << "ERROR: cannot compute Minkowski functional value; recieved wrong dimensions:" << std::endl
	      << "       sample is a " << sample.call_Nx() << "x" << sample.call_Ny() << "BinField, but" << std::endl
	      << "       SubBinField is " << Nx << "x" << Ny << " at coordinate (" << xi << "," << yi << ");" << std::endl
	      << std::endl;
    exit(-1);
  }
  
  unsigned X;
  unsigned Y;
  
  
  // -----------------------
  // BOTTOM ROW
  Y = 0;
  
  // Range 7
  X = 0;
  total += rg7_perimeter.at(convert(sample.call(xi+X,yi+Ny-1),sample.call(xi+Nx-1,yi+Ny-1),sample.call(xi+X,yi+Y),sample.call(xi+Nx-1,yi+Y)));
  
  // Range 8
  X++;
  while(X < Nx)
    {
      total += rg8_perimeter.at(convert(sample.call(xi+X,yi+Ny-1),sample.call(xi+X-1,yi+Ny-1),sample.call(xi+X,yi+Y),sample.call(xi+X-1,yi+Y)));
      X++;
    }
  
  // Range 9
  total += rg9_perimeter.at(convert(sample.call(xi,yi+Ny-1),sample.call(xi+Nx-1,yi+Ny-1),sample.call(xi,yi+Y),sample.call(xi+Nx-1,yi+Y)));
  
  // -----------------------
  // MID ROWS
  Y++;
  while(Y < Ny)
    {
  
      // Range 4
      X = 0;
      total += rg4_perimeter.at(convert(sample.call(xi+X,yi+Y-1),sample.call(xi+Nx-1,yi+Y-1),sample.call(xi+X,yi+Y),sample.call(xi+Nx-1,yi+Y)));
  
      // Range 5
      X++;
      while(X < Nx)
	{
	  total += rg5_perimeter.at(convert(sample.call(xi+X,yi+Y-1),sample.call(xi+X-1,yi+Y-1),sample.call(xi+X,yi+Y),sample.call(xi+X-1,yi+Y)));
      	  X++;
	}
  
      // Range 6
      total += rg6_perimeter.at(convert(sample.call(xi,yi+Y-1),sample.call(xi+Nx-1,yi+Y-1),sample.call(xi,yi+Y),sample.call(xi+Nx-1,yi+Y)));
  
      
      Y++;  
    }
  
  // -----------------------
  // TOP ROW
  
  // Range 1
  X = 0;
  total += rg1_perimeter.at(convert(sample.call(xi+X,yi+Y-1),sample.call(xi+Nx-1,yi+Y-1),sample.call(xi+X,yi),sample.call(xi+Nx-1,yi)));
  
  // Range 2
  X++;
  while(X < Nx)
    {
      total += rg2_perimeter.at(convert(sample.call(xi+X,yi+Y-1),sample.call(xi+X-1,yi+Y-1),sample.call(xi+X,yi),sample.call(xi+X-1,yi)));
      X++;
    }
  
  // Range 3
  total += rg3_perimeter.at(convert(sample.call(xi,yi+Y-1),sample.call(xi+Nx-1,yi+Y-1),sample.call(xi,yi),sample.call(xi+Nx-1,yi)));
  
  return total;
}

double euler_pbc(const BinField<bool> &sample,
		 const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny)
{
  double total = 0;
  
  if(! ((xi>=0)&&(xi<(sample.call_Nx()-Nx+1))&&(yi>=0)&&(yi<(sample.call_Ny()-Ny+1))&&(Nx>2)&&(Ny>2)) ){
    std::cerr << "ERROR: cannot compute Minkowski functional value; recieved wrong dimensions:" << std::endl
	      << "       sample is a " << sample.call_Nx() << "x" << sample.call_Ny() << "BinField, but" << std::endl
	      << "       SubBinField is " << Nx << "x" << Ny << " at coordinate (" << xi << "," << yi << ");" << std::endl
	      << std::endl;
    exit(-1);
  }
  
  unsigned X;
  unsigned Y;
  
  
  // -----------------------
  // BOTTOM ROW
  Y = 0;
  
  // Range 7
  X = 0;
  total += rg7_euler.at(convert(sample.call(xi+X,yi+Ny-1),sample.call(xi+Nx-1,yi+Ny-1),sample.call(xi+X,yi+Y),sample.call(xi+Nx-1,yi+Y)));
  
  // Range 8
  X++;
  while(X < Nx)
    {
      total += rg8_euler.at(convert(sample.call(xi+X,yi+Ny-1),sample.call(xi+X-1,yi+Ny-1),sample.call(xi+X,yi+Y),sample.call(xi+X-1,yi+Y)));
      X++;
    }
  
  // Range 9
  total += rg9_euler.at(convert(sample.call(xi,yi+Ny-1),sample.call(xi+Nx-1,yi+Ny-1),sample.call(xi,yi+Y),sample.call(xi+Nx-1,yi+Y)));
  
  // -----------------------
  // MID ROWS
  Y++;
  while(Y < Ny)
    {
  
      // Range 4
      X = 0;
      total += rg4_euler.at(convert(sample.call(xi+X,yi+Y-1),sample.call(xi+Nx-1,yi+Y-1),sample.call(xi+X,yi+Y),sample.call(xi+Nx-1,yi+Y)));
  
      // Range 5
      X++;
      while(X < Nx)
	{
	  total += rg5_euler.at(convert(sample.call(xi+X,yi+Y-1),sample.call(xi+X-1,yi+Y-1),sample.call(xi+X,yi+Y),sample.call(xi+X-1,yi+Y)));
      	  X++;
	}
  
      // Range 6
      total += rg6_euler.at(convert(sample.call(xi,yi+Y-1),sample.call(xi+Nx-1,yi+Y-1),sample.call(xi,yi+Y),sample.call(xi+Nx-1,yi+Y)));
  
      
      Y++;  
    }
  
  // -----------------------
  // TOP ROW
  
  // Range 1
  X = 0;
  total += rg1_euler.at(convert(sample.call(xi+X,yi+Y-1),sample.call(xi+Nx-1,yi+Y-1),sample.call(xi+X,yi),sample.call(xi+Nx-1,yi)));
  
  // Range 2
  X++;
  while(X < Nx)
    {
      total += rg2_euler.at(convert(sample.call(xi+X,yi+Y-1),sample.call(xi+X-1,yi+Y-1),sample.call(xi+X,yi),sample.call(xi+X-1,yi)));
      X++;
    }
  
  // Range 3
  total += rg3_euler.at(convert(sample.call(xi,yi+Y-1),sample.call(xi+Nx-1,yi+Y-1),sample.call(xi,yi),sample.call(xi+Nx-1,yi)));
  
  return total;
}

double w102_xx_pbc(const BinField<bool> &sample,
		   const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny)
{
  double total = 0;
  
  if(! ((xi>=0)&&(xi<(sample.call_Nx()-Nx+1))&&(yi>=0)&&(yi<(sample.call_Ny()-Ny+1))&&(Nx>2)&&(Ny>2)) ){
    std::cerr << "ERROR: cannot compute Minkowski functional value; recieved wrong dimensions:" << std::endl
	      << "       sample is a " << sample.call_Nx() << "x" << sample.call_Ny() << "BinField, but" << std::endl
	      << "       SubBinField is " << Nx << "x" << Ny << " at coordinate (" << xi << "," << yi << ");" << std::endl
	      << std::endl;
    exit(-1);
  }
  
  unsigned X;
  unsigned Y;
  
  
  // -----------------------
  // BOTTOM ROW
  Y = 0;
  
  // Range 7
  X = 0;
  total += rg7_w102_xx.at(convert(sample.call(xi+X,yi+Ny-1),sample.call(xi+Nx-1,yi+Ny-1),sample.call(xi+X,yi+Y),sample.call(xi+Nx-1,yi+Y)));
  
  // Range 8
  X++;
  while(X < Nx)
    {
      total += rg8_w102_xx.at(convert(sample.call(xi+X,yi+Ny-1),sample.call(xi+X-1,yi+Ny-1),sample.call(xi+X,yi+Y),sample.call(xi+X-1,yi+Y)));
      X++;
    }
  
  // Range 9
  total += rg9_w102_xx.at(convert(sample.call(xi,yi+Ny-1),sample.call(xi+Nx-1,yi+Ny-1),sample.call(xi,yi+Y),sample.call(xi+Nx-1,yi+Y)));
  
  // -----------------------
  // MID ROWS
  Y++;
  while(Y < Ny)
    {
  
      // Range 4
      X = 0;
      total += rg4_w102_xx.at(convert(sample.call(xi+X,yi+Y-1),sample.call(xi+Nx-1,yi+Y-1),sample.call(xi+X,yi+Y),sample.call(xi+Nx-1,yi+Y)));
  
      // Range 5
      X++;
      while(X < Nx)
	{
	  total += rg5_w102_xx.at(convert(sample.call(xi+X,yi+Y-1),sample.call(xi+X-1,yi+Y-1),sample.call(xi+X,yi+Y),sample.call(xi+X-1,yi+Y)));
      	  X++;
	}
  
      // Range 6
      total += rg6_w102_xx.at(convert(sample.call(xi,yi+Y-1),sample.call(xi+Nx-1,yi+Y-1),sample.call(xi,yi+Y),sample.call(xi+Nx-1,yi+Y)));
  
      
      Y++;  
    }
  
  // -----------------------
  // TOP ROW
  
  // Range 1
  X = 0;
  total += rg1_w102_xx.at(convert(sample.call(xi+X,yi+Y-1),sample.call(xi+Nx-1,yi+Y-1),sample.call(xi+X,yi),sample.call(xi+Nx-1,yi)));
  
  // Range 2
  X++;
  while(X < Nx)
    {
      total += rg2_w102_xx.at(convert(sample.call(xi+X,yi+Y-1),sample.call(xi+X-1,yi+Y-1),sample.call(xi+X,yi),sample.call(xi+X-1,yi)));
      X++;
    }
  
  // Range 3
  total += rg3_w102_xx.at(convert(sample.call(xi,yi+Y-1),sample.call(xi+Nx-1,yi+Y-1),sample.call(xi,yi),sample.call(xi+Nx-1,yi)));
  
  return total;
}

double w102_xy_pbc(const BinField<bool> &sample,
		   const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny)
{
  double total = 0;
  
  if(! ((xi>=0)&&(xi<(sample.call_Nx()-Nx+1))&&(yi>=0)&&(yi<(sample.call_Ny()-Ny+1))&&(Nx>2)&&(Ny>2)) ){
    std::cerr << "ERROR: cannot compute Minkowski functional value; recieved wrong dimensions:" << std::endl
	      << "       sample is a " << sample.call_Nx() << "x" << sample.call_Ny() << "BinField, but" << std::endl
	      << "       SubBinField is " << Nx << "x" << Ny << " at coordinate (" << xi << "," << yi << ");" << std::endl
	      << std::endl;
    exit(-1);
  }
  
  unsigned X;
  unsigned Y;
  
  
  // -----------------------
  // BOTTOM ROW
  Y = 0;
  
  // Range 7
  X = 0;
  total += rg7_w102_xy.at(convert(sample.call(xi+X,yi+Ny-1),sample.call(xi+Nx-1,yi+Ny-1),sample.call(xi+X,yi+Y),sample.call(xi+Nx-1,yi+Y)));
  
  // Range 8
  X++;
  while(X < Nx)
    {
      total += rg8_w102_xy.at(convert(sample.call(xi+X,yi+Ny-1),sample.call(xi+X-1,yi+Ny-1),sample.call(xi+X,yi+Y),sample.call(xi+X-1,yi+Y)));
      X++;
    }
  
  // Range 9
  total += rg9_w102_xy.at(convert(sample.call(xi,yi+Ny-1),sample.call(xi+Nx-1,yi+Ny-1),sample.call(xi,yi+Y),sample.call(xi+Nx-1,yi+Y)));
  
  // -----------------------
  // MID ROWS
  Y++;
  while(Y < Ny)
    {
  
      // Range 4
      X = 0;
      total += rg4_w102_xy.at(convert(sample.call(xi+X,yi+Y-1),sample.call(xi+Nx-1,yi+Y-1),sample.call(xi+X,yi+Y),sample.call(xi+Nx-1,yi+Y)));
  
      // Range 5
      X++;
      while(X < Nx)
	{
	  total += rg5_w102_xy.at(convert(sample.call(xi+X,yi+Y-1),sample.call(xi+X-1,yi+Y-1),sample.call(xi+X,yi+Y),sample.call(xi+X-1,yi+Y)));
      	  X++;
	}
  
      // Range 6
      total += rg6_w102_xy.at(convert(sample.call(xi,yi+Y-1),sample.call(xi+Nx-1,yi+Y-1),sample.call(xi,yi+Y),sample.call(xi+Nx-1,yi+Y)));
  
      
      Y++;  
    }
  
  // -----------------------
  // TOP ROW
  
  // Range 1
  X = 0;
  total += rg1_w102_xy.at(convert(sample.call(xi+X,yi+Y-1),sample.call(xi+Nx-1,yi+Y-1),sample.call(xi+X,yi),sample.call(xi+Nx-1,yi)));
  
  // Range 2
  X++;
  while(X < Nx)
    {
      total += rg2_w102_xy.at(convert(sample.call(xi+X,yi+Y-1),sample.call(xi+X-1,yi+Y-1),sample.call(xi+X,yi),sample.call(xi+X-1,yi)));
      X++;
    }
  
  // Range 3
  total += rg3_w102_xy.at(convert(sample.call(xi,yi+Y-1),sample.call(xi+Nx-1,yi+Y-1),sample.call(xi,yi),sample.call(xi+Nx-1,yi)));
  
  return total;
}

double w102_yy_pbc(const BinField<bool> &sample,
		   const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny)
{
  double total = 0;
  
  if(! ((xi>=0)&&(xi<(sample.call_Nx()-Nx+1))&&(yi>=0)&&(yi<(sample.call_Ny()-Ny+1))&&(Nx>2)&&(Ny>2)) ){
    std::cerr << "ERROR: cannot compute Minkowski functional value; recieved wrong dimensions:" << std::endl
	      << "       sample is a " << sample.call_Nx() << "x" << sample.call_Ny() << "BinField, but" << std::endl
	      << "       SubBinField is " << Nx << "x" << Ny << " at coordinate (" << xi << "," << yi << ");" << std::endl
	      << std::endl;
    exit(-1);
  }
  
  unsigned X;
  unsigned Y;
  
  
  // -----------------------
  // BOTTOM ROW
  Y = 0;
  
  // Range 7
  X = 0;
  total += rg7_w102_yy.at(convert(sample.call(xi+X,yi+Ny-1),sample.call(xi+Nx-1,yi+Ny-1),sample.call(xi+X,yi+Y),sample.call(xi+Nx-1,yi+Y)));
  
  // Range 8
  X++;
  while(X < Nx)
    {
      total += rg8_w102_yy.at(convert(sample.call(xi+X,yi+Ny-1),sample.call(xi+X-1,yi+Ny-1),sample.call(xi+X,yi+Y),sample.call(xi+X-1,yi+Y)));
      X++;
    }
  
  // Range 9
  total += rg9_w102_yy.at(convert(sample.call(xi,yi+Ny-1),sample.call(xi+Nx-1,yi+Ny-1),sample.call(xi,yi+Y),sample.call(xi+Nx-1,yi+Y)));
  
  // -----------------------
  // MID ROWS
  Y++;
  while(Y < Ny)
    {
  
      // Range 4
      X = 0;
      total += rg4_w102_yy.at(convert(sample.call(xi+X,yi+Y-1),sample.call(xi+Nx-1,yi+Y-1),sample.call(xi+X,yi+Y),sample.call(xi+Nx-1,yi+Y)));
  
      // Range 5
      X++;
      while(X < Nx)
	{
	  total += rg5_w102_yy.at(convert(sample.call(xi+X,yi+Y-1),sample.call(xi+X-1,yi+Y-1),sample.call(xi+X,yi+Y),sample.call(xi+X-1,yi+Y)));
      	  X++;
	}
  
      // Range 6
      total += rg6_w102_yy.at(convert(sample.call(xi,yi+Y-1),sample.call(xi+Nx-1,yi+Y-1),sample.call(xi,yi+Y),sample.call(xi+Nx-1,yi+Y)));
  
      
      Y++;  
    }
  
  // -----------------------
  // TOP ROW
  
  // Range 1
  X = 0;
  total += rg1_w102_yy.at(convert(sample.call(xi+X,yi+Y-1),sample.call(xi+Nx-1,yi+Y-1),sample.call(xi+X,yi),sample.call(xi+Nx-1,yi)));
  
  // Range 2
  X++;
  while(X < Nx)
    {
      total += rg2_w102_yy.at(convert(sample.call(xi+X,yi+Y-1),sample.call(xi+X-1,yi+Y-1),sample.call(xi+X,yi),sample.call(xi+X-1,yi)));
      X++;
    }
  
  // Range 3
  total += rg3_w102_yy.at(convert(sample.call(xi,yi+Y-1),sample.call(xi+Nx-1,yi+Y-1),sample.call(xi,yi),sample.call(xi+Nx-1,yi)));
  
  return total;
}

double delta_pbc(const BinField<bool> &sample,
		 const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny)
{
  double w102_xx_ = w102_xx_pbc(sample,xi,yi,Nx,Ny);
  double w102_xy_ = w102_xy_pbc(sample,xi,yi,Nx,Ny);
  double w102_yy_ = w102_yy_pbc(sample,xi,yi,Nx,Ny);
  
  return sqrt((w102_xx_-w102_yy_)*(w102_xx_-w102_yy_)+4*w102_xy_*w102_xy_);
}

// PIXELIZED DATA
// Minus sampling boundary condition
double area_pixelized_mbc(const BinField<bool> &sample,
			  const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny)
{
  double total = 0;
  
  if(! ((xi>=0)&&(xi<(sample.call_Nx()-Nx+1))&&(yi>=0)&&(yi<(sample.call_Ny()-Ny+1))&&(Nx>2)&&(Ny>2)) ){
    std::cerr << "ERROR: cannot compute Minkowski functional value; recieved wrong dimensions:" << std::endl
	      << "       sample is a " << sample.call_Nx() << "x" << sample.call_Ny() << "BinField, but" << std::endl
	      << "       SubBinField is " << Nx << "x" << Ny << " at coordinate (" << xi << "," << yi << ");" << std::endl
	      << std::endl;
    exit(-1);
  }
  
  for(unsigned ix = xi+1; ix < (xi+Nx-1); ix++)
    for(unsigned iy = yi+1; iy < (yi+Ny-1); iy++)
      total += sample.call(ix,iy);
  
  return total;  
}

double perimeter_pixelized_mbc(const BinField<bool> &sample,
			       const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny)
{
  double total = 0;
  
  if(! ((xi>=0)&&(xi<(sample.call_Nx()-Nx+1))&&(yi>=0)&&(yi<(sample.call_Ny()-Ny+1))&&(Nx>2)&&(Ny>2)) ){
    std::cerr << "ERROR: cannot compute Minkowski functional value; recieved wrong dimensions:" << std::endl
	      << "       sample is a " << sample.call_Nx() << "x" << sample.call_Ny() << "BinField, but" << std::endl
	      << "       SubBinField is " << Nx << "x" << Ny << " at coordinate (" << xi << "," << yi << ");" << std::endl
	      << std::endl;
    exit(-1);
  }
  
  unsigned ix = xi;
  unsigned iy = yi;
  for(iy = yi+1; iy < (yi+Ny-1); iy++)
    total += ( (sample.call(ix,iy)==false) && sample.call(ix+1,iy) );
  
  ix = xi+Nx-2;
  for(iy = yi+1; iy < (yi+Ny-1); iy++)
    total += ( (sample.call(ix+1,iy)==false) && sample.call(ix,iy) );
  
  for(ix = xi+1; ix < (xi+Nx-2); ix++)
    for(iy = yi+1; iy < (yi+Ny-1); iy++)
      total += ( (sample.call(ix+1,iy)==false) && sample.call(ix,iy) ) + ( (sample.call(ix,iy)==false) && sample.call(ix+1,iy) );
  
  iy = yi;
  for(ix = xi+1; ix < (xi+Nx-1); ix++)
    total += ( (sample.call(ix,iy)==false) && sample.call(ix,iy+1) );
  
  iy = yi+Ny-2;
  for(ix = xi+1; ix < (xi+Nx-1); ix++)
    total += ( (sample.call(ix,iy+1)==false) && sample.call(ix,iy) );
  
  for(iy = yi+1; iy < (yi+Ny-2); iy++)
    for(ix = xi+1; ix < (xi+Nx-1); ix++)
    total += ( (sample.call(ix,iy)==false) && sample.call(ix,iy+1) ) + ( (sample.call(ix,iy+1)==false) && sample.call(ix,iy) );
  
  return total;
}

double euler_pixelized_mbc(const BinField<bool> &sample,
			   const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny)
{
  return euler_mbc(sample, xi, yi, Nx, Ny);
}

// Periodic boundary condition
double area_pixelized_pbc(const BinField<bool> &sample,
			  const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny)
{
  double total = 0;
  
  if(! ((xi>=0)&&(xi<(sample.call_Nx()-Nx+1))&&(yi>=0)&&(yi<(sample.call_Ny()-Ny+1))&&(Nx>2)&&(Ny>2)) ){
    std::cerr << "ERROR: cannot compute Minkowski functional value; recieved wrong dimensions:" << std::endl
	      << "       sample is a " << sample.call_Nx() << "x" << sample.call_Ny() << "BinField, but" << std::endl
	      << "       SubBinField is " << Nx << "x" << Ny << " at coordinate (" << xi << "," << yi << ");" << std::endl
	      << std::endl;
    exit(-1);
  }
  
  for(unsigned ix = xi; ix < (xi+Nx); ix++)
    for(unsigned iy = yi; iy < (yi+Ny); iy++)
      total += sample.call(ix,iy);
  
  return total;
}

double perimeter_pixelized_pbc(const BinField<bool> &sample,
			       const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny)
{
  double total = 0;
  
  if(! ((xi>=0)&&(xi<(sample.call_Nx()-Nx+1))&&(yi>=0)&&(yi<(sample.call_Ny()-Ny+1))&&(Nx>2)&&(Ny>2)) ){
    std::cerr << "ERROR: cannot compute Minkowski functional value; recieved wrong dimensions:" << std::endl
	      << "       sample is a " << sample.call_Nx() << "x" << sample.call_Ny() << "BinField, but" << std::endl
	      << "       SubBinField is " << Nx << "x" << Ny << " at coordinate (" << xi << "," << yi << ");" << std::endl
	      << std::endl;
    exit(-1);
  }
  
  unsigned ix = xi;
  unsigned iy = yi;
  for(iy = yi; iy < yi+Ny; iy++){
    total += ( (sample.call(Nx-1,iy)==false) && sample.call(0,iy) );
    total += ( (sample.call(0,iy)==false) && sample.call(Nx-1,iy) );
  }
  
  for(ix = xi+1; ix < xi+Nx; ix++)
    for(iy = yi; iy < yi+Ny; iy++){
      total += ( (sample.call(ix-1,iy)==false) && sample.call(ix,iy) );
      total += ( (sample.call(ix,iy)==false) && sample.call(ix-1,iy) );
    }
  
  for(ix = xi; ix < xi+Nx; ix++){
    total += ( (sample.call(ix,Ny-1)==false) && sample.call(ix,0) );
    total += ( (sample.call(ix,0)==false) && sample.call(ix,Ny-1) );
  }

  for(iy = yi+1; iy < yi+Ny; iy++)
    for(ix = xi; ix < xi+Nx; ix++){
      total += ( (sample.call(ix,iy)==false) && sample.call(ix,iy-1) );
      total += ( (sample.call(ix,iy-1)==false) && sample.call(ix,iy) );
    }
  
  return total;
}

double euler_pixelized_pbc(const BinField<bool> &sample,
			   const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny)
{
  return euler_pbc(sample, xi, yi, Nx, Ny);
}

// PAPAYA
void print_pgm(const BinField<bool> &sample,
	       const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny,
	       const std::string &filename, const std::string &prefix_of, const bool &invert)
{
  
  if(! ((xi>=0)&&(xi<(sample.call_Nx()-Nx+1))&&(yi>=0)&&(yi<(sample.call_Ny()-Ny+1))&&(Nx>2)&&(Ny>2)) ){
    std::cerr << "ERROR: cannot compute Minkowski functional value; recieved wrong dimensions:" << std::endl
	      << "       sample is a " << sample.call_Nx() << "x" << sample.call_Ny() << "BinField, but" << std::endl
	      << "       SubBinField is " << Nx << "x" << Ny << " at coordinate (" << xi << "," << yi << ");" << std::endl
	      << std::endl;
    exit(-1);
  }
  
  std::ofstream OutFile( (prefix_of + filename).c_str() );
  OutFile << "P2" << std::endl;
  OutFile << Nx << " " << Ny << std::endl;
  OutFile << "255" << std::endl;

  unsigned truephase = 255, voidphase = 0;
  if(invert){
    truephase = 0; voidphase = 255;}
  
  for(unsigned j = (Ny-1); j < Ny; j--)
    for(unsigned i = 0; i < Nx; i++)
      {
	if( sample.call(i,j) )
	  OutFile << truephase << std::endl;
	else
	  OutFile << voidphase << std::endl;
      }
  
  OutFile.close();
}

// ---------------------------------------------------------------
// Only BinField as argument:

// Minus sampling boundary condition
double area_mbc(const BinField<bool> &sample)
{
  double total = 0;
  unsigned Nx = sample.call_Nx();
  unsigned Ny = sample.call_Ny();
  
  unsigned X;
  unsigned Y;
  
  
  // -----------------------
  // BOTTOM ROW
  Y = 0;
  
  // Range 7
  X = 0;
  total += rg7_area.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  // Range 8
  X++;
  while(X < (Nx-2))
    {
      total += rg8_area.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
      
      X++;
    }
  
  // Range 9
  total += rg9_area.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  // -----------------------
  // MID ROWS
  Y++;
  while(Y < (Ny-2))
    {
  
      // Range 4
      X = 0;
      total += rg4_area.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
      // Range 5
      X++;
      while(X < (Nx-2))
	{
	  total += rg5_area.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
      
	  X++;
	}
  
      // Range 6
      total += rg6_area.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
      
      Y++;  
    }
  
  // -----------------------
  // TOP ROW
  
  // Range 1
  X = 0;
  total += rg1_area.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  // Range 2
  X++;
  while(X < (Nx-2))
    {
      total += rg2_area.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
      
      X++;
    }
  
  // Range 3
  total += rg3_area.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  return total;
}

double perimeter_mbc(const BinField<bool> &sample)
{
  double total = 0;
  unsigned Nx = sample.call_Nx();
  unsigned Ny = sample.call_Ny();
  
  unsigned X;
  unsigned Y;
  
  
  // -----------------------
  // BOTTOM ROW
  Y = 0;
  
  // Range 7
  X = 0;
  total += rg7_perimeter.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  // Range 8
  X++;
  while(X < (Nx-2))
    {
      total += rg8_perimeter.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
      
      X++;
    }
  
  // Range 9
  total += rg9_perimeter.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  // -----------------------
  // MID ROWS
  Y++;
  while(Y < (Ny-2))
    {
  
      // Range 4
      X = 0;
      total += rg4_perimeter.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
      // Range 5
      X++;
      while(X < (Nx-2))
	{
	  total += rg5_perimeter.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
      
	  X++;
	}
  
      // Range 6
      total += rg6_perimeter.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
      
      Y++;  
    }
  
  // -----------------------
  // TOP ROW
  
  // Range 1
  X = 0;
  total += rg1_perimeter.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  // Range 2
  X++;
  while(X < (Nx-2))
    {
      total += rg2_perimeter.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
      
      X++;
    }
  
  // Range 3
  total += rg3_perimeter.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  return total;
}

double euler_mbc(const BinField<bool> &sample)
{
  double total = 0;
  unsigned Nx = sample.call_Nx();
  unsigned Ny = sample.call_Ny();
  
  unsigned X;
  unsigned Y;
  
  
  // -----------------------
  // BOTTOM ROW
  Y = 0;
  
  // Range 7
  X = 0;
  total += rg7_euler.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  // Range 8
  X++;
  while(X < (Nx-2))
    {
      total += rg8_euler.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
      
      X++;
    }
  
  // Range 9
  total += rg9_euler.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  // -----------------------
  // MID ROWS
  Y++;
  while(Y < (Ny-2))
    {
  
      // Range 4
      X = 0;
      total += rg4_euler.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
      // Range 5
      X++;
      while(X < (Nx-2))
	{
	  total += rg5_euler.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
      
	  X++;
	}
  
      // Range 6
      total += rg6_euler.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
      
      Y++;  
    }
  
  // -----------------------
  // TOP ROW
  
  // Range 1
  X = 0;
  total += rg1_euler.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  // Range 2
  X++;
  while(X < (Nx-2))
    {
      total += rg2_euler.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
      
      X++;
    }
  
  // Range 3
  total += rg3_euler.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  return total;
}

double w102_xx_mbc(const BinField<bool> &sample)
{
  double total = 0;
  unsigned Nx = sample.call_Nx();
  unsigned Ny = sample.call_Ny();
  
  unsigned X;
  unsigned Y;
  
  
  // -----------------------
  // BOTTOM ROW
  Y = 0;
  
  // Range 7
  X = 0;
  total += rg7_w102_xx.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  // Range 8
  X++;
  while(X < (Nx-2))
    {
      total += rg8_w102_xx.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
      
      X++;
    }
  
  // Range 9
  total += rg9_w102_xx.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  // -----------------------
  // MID ROWS
  Y++;
  while(Y < (Ny-2))
    {
  
      // Range 4
      X = 0;
      total += rg4_w102_xx.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
      // Range 5
      X++;
      while(X < (Nx-2))
	{
	  total += rg5_w102_xx.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
      
	  X++;
	}
  
      // Range 6
      total += rg6_w102_xx.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
      
      Y++;  
    }
  
  // -----------------------
  // TOP ROW
  
  // Range 1
  X = 0;
  total += rg1_w102_xx.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  // Range 2
  X++;
  while(X < (Nx-2))
    {
      total += rg2_w102_xx.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
      
      X++;
    }
  
  // Range 3
  total += rg3_w102_xx.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  return total;
}

double w102_xy_mbc(const BinField<bool> &sample)
{
  double total = 0;
  unsigned Nx = sample.call_Nx();
  unsigned Ny = sample.call_Ny();
  
  unsigned X;
  unsigned Y;
  
  
  // -----------------------
  // BOTTOM ROW
  Y = 0;
  
  // Range 7
  X = 0;
  total += rg7_w102_xy.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  // Range 8
  X++;
  while(X < (Nx-2))
    {
      total += rg8_w102_xy.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
      
      X++;
    }
  
  // Range 9
  total += rg9_w102_xy.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  // -----------------------
  // MID ROWS
  Y++;
  while(Y < (Ny-2))
    {
  
      // Range 4
      X = 0;
      total += rg4_w102_xy.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
      // Range 5
      X++;
      while(X < (Nx-2))
	{
	  total += rg5_w102_xy.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
      
	  X++;
	}
  
      // Range 6
      total += rg6_w102_xy.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
      
      Y++;  
    }
  
  // -----------------------
  // TOP ROW
  
  // Range 1
  X = 0;
  total += rg1_w102_xy.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  // Range 2
  X++;
  while(X < (Nx-2))
    {
      total += rg2_w102_xy.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
      
      X++;
    }
  
  // Range 3
  total += rg3_w102_xy.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  return total;
}

double w102_yy_mbc(const BinField<bool> &sample)
{
  double total = 0;
  unsigned Nx = sample.call_Nx();
  unsigned Ny = sample.call_Ny();
  
  unsigned X;
  unsigned Y;
  
  
  // -----------------------
  // BOTTOM ROW
  Y = 0;
  
  // Range 7
  X = 0;
  total += rg7_w102_yy.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  // Range 8
  X++;
  while(X < (Nx-2))
    {
      total += rg8_w102_yy.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
      
      X++;
    }
  
  // Range 9
  total += rg9_w102_yy.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  // -----------------------
  // MID ROWS
  Y++;
  while(Y < (Ny-2))
    {
  
      // Range 4
      X = 0;
      total += rg4_w102_yy.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
      // Range 5
      X++;
      while(X < (Nx-2))
	{
	  total += rg5_w102_yy.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
      
	  X++;
	}
  
      // Range 6
      total += rg6_w102_yy.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
      
      Y++;  
    }
  
  // -----------------------
  // TOP ROW
  
  // Range 1
  X = 0;
  total += rg1_w102_yy.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  // Range 2
  X++;
  while(X < (Nx-2))
    {
      total += rg2_w102_yy.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
      
      X++;
    }
  
  // Range 3
  total += rg3_w102_yy.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  return total;
}

double delta_mbc(const BinField<bool> &sample)
{
  double w102_xx_ = w102_xx_mbc(sample);
  double w102_xy_ = w102_xy_mbc(sample);
  double w102_yy_ = w102_yy_mbc(sample);
  
  return sqrt((w102_xx_-w102_yy_)*(w102_xx_-w102_yy_)+4*w102_xy_*w102_xy_);
}

// Periodic boundary condition
double area_pbc(const BinField<bool> &sample)
{
  double total = 0;
  unsigned Nx = sample.call_Nx();
  unsigned Ny = sample.call_Ny();
  
  unsigned X;
  unsigned Y;
  
  
  // -----------------------
  // BOTTOM ROW
  Y = 0;
  
  // Range 7
  X = 0;
  total += rg7_area.at(convert(sample.call(X,Ny-1),sample.call(Nx-1,Ny-1),sample.call(X,Y),sample.call(Nx-1,Y)));
  
  // Range 8
  X++;
  while(X < Nx)
    {
      total += rg8_area.at(convert(sample.call(X,Ny-1),sample.call(X-1,Ny-1),sample.call(X,Y),sample.call(X-1,Y)));
      X++;
    }
  
  // Range 9
  total += rg9_area.at(convert(sample.call(0,Ny-1),sample.call(Nx-1,Ny-1),sample.call(0,Y),sample.call(Nx-1,Y)));
  
  // -----------------------
  // MID ROWS
  Y++;
  while(Y < Ny)
    {
      
      // Range 4
      X = 0;
      total += rg4_area.at(convert(sample.call(X,Y-1),sample.call(Nx-1,Y-1),sample.call(X,Y),sample.call(Nx-1,Y)));
      
      // Range 5
      X++;
      while(X < Nx)
	{
	  total += rg5_area.at(convert(sample.call(X,Y-1),sample.call(X-1,Y-1),sample.call(X,Y),sample.call(X-1,Y)));
      	  X++;
	}
      
      // Range 6
      total += rg6_area.at(convert(sample.call(0,Y-1),sample.call(Nx-1,Y-1),sample.call(0,Y),sample.call(Nx-1,Y)));
      
      
      Y++;  
    }
  
  // -----------------------
  // TOP ROW
  
  // Range 1
  X = 0;
  total += rg1_area.at(convert(sample.call(X,Y-1),sample.call(Nx-1,Y-1),sample.call(X,0),sample.call(Nx-1,0)));
  
  // Range 2
  X++;
  while(X < Nx)
    {
      total += rg2_area.at(convert(sample.call(X,Y-1),sample.call(X-1,Y-1),sample.call(X,0),sample.call(X-1,0)));
      X++;
    }
  
  // Range 3
  total += rg3_area.at(convert(sample.call(0,Y-1),sample.call(Nx-1,Y-1),sample.call(0,0),sample.call(Nx-1,0)));
  
  return total;
}

double perimeter_pbc(const BinField<bool> &sample)
{
  double total = 0;
  unsigned Nx = sample.call_Nx();
  unsigned Ny = sample.call_Ny();
  
  unsigned X;
  unsigned Y;
  
  
  // -----------------------
  // BOTTOM ROW
  Y = 0;
  
  // Range 7
  X = 0;
  total += rg7_perimeter.at(convert(sample.call(X,Ny-1),sample.call(Nx-1,Ny-1),sample.call(X,Y),sample.call(Nx-1,Y)));
  
  // Range 8
  X++;
  while(X < Nx)
    {
      total += rg8_perimeter.at(convert(sample.call(X,Ny-1),sample.call(X-1,Ny-1),sample.call(X,Y),sample.call(X-1,Y)));
      X++;
    }
  
  // Range 9
  total += rg9_perimeter.at(convert(sample.call(0,Ny-1),sample.call(Nx-1,Ny-1),sample.call(0,Y),sample.call(Nx-1,Y)));
  
  // -----------------------
  // MID ROWS
  Y++;
  while(Y < Ny)
    {
  
      // Range 4
      X = 0;
      total += rg4_perimeter.at(convert(sample.call(X,Y-1),sample.call(Nx-1,Y-1),sample.call(X,Y),sample.call(Nx-1,Y)));
  
      // Range 5
      X++;
      while(X < Nx)
	{
	  total += rg5_perimeter.at(convert(sample.call(X,Y-1),sample.call(X-1,Y-1),sample.call(X,Y),sample.call(X-1,Y)));
      	  X++;
	}
  
      // Range 6
      total += rg6_perimeter.at(convert(sample.call(0,Y-1),sample.call(Nx-1,Y-1),sample.call(0,Y),sample.call(Nx-1,Y)));
  
      
      Y++;  
    }
  
  // -----------------------
  // TOP ROW
  
  // Range 1
  X = 0;
  total += rg1_perimeter.at(convert(sample.call(X,Y-1),sample.call(Nx-1,Y-1),sample.call(X,0),sample.call(Nx-1,0)));
  
  // Range 2
  X++;
  while(X < Nx)
    {
      total += rg2_perimeter.at(convert(sample.call(X,Y-1),sample.call(X-1,Y-1),sample.call(X,0),sample.call(X-1,0)));
      X++;
    }
  
  // Range 3
  total += rg3_perimeter.at(convert(sample.call(0,Y-1),sample.call(Nx-1,Y-1),sample.call(0,0),sample.call(Nx-1,0)));
  
  return total;
}

double euler_pbc(const BinField<bool> &sample)
{
  double total = 0;
  unsigned Nx = sample.call_Nx();
  unsigned Ny = sample.call_Ny();
  
  unsigned X;
  unsigned Y;
  
  
  // -----------------------
  // BOTTOM ROW
  Y = 0;
  
  // Range 7
  X = 0;
  total += rg7_euler.at(convert(sample.call(X,Ny-1),sample.call(Nx-1,Ny-1),sample.call(X,Y),sample.call(Nx-1,Y)));
  
  // Range 8
  X++;
  while(X < Nx)
    {
      total += rg8_euler.at(convert(sample.call(X,Ny-1),sample.call(X-1,Ny-1),sample.call(X,Y),sample.call(X-1,Y)));
      X++;
    }
  
  // Range 9
  total += rg9_euler.at(convert(sample.call(0,Ny-1),sample.call(Nx-1,Ny-1),sample.call(0,Y),sample.call(Nx-1,Y)));
  
  // -----------------------
  // MID ROWS
  Y++;
  while(Y < Ny)
    {
  
      // Range 4
      X = 0;
      total += rg4_euler.at(convert(sample.call(X,Y-1),sample.call(Nx-1,Y-1),sample.call(X,Y),sample.call(Nx-1,Y)));
  
      // Range 5
      X++;
      while(X < Nx)
	{
	  total += rg5_euler.at(convert(sample.call(X,Y-1),sample.call(X-1,Y-1),sample.call(X,Y),sample.call(X-1,Y)));
      	  X++;
	}
  
      // Range 6
      total += rg6_euler.at(convert(sample.call(0,Y-1),sample.call(Nx-1,Y-1),sample.call(0,Y),sample.call(Nx-1,Y)));
  
      
      Y++;  
    }
  
  // -----------------------
  // TOP ROW
  
  // Range 1
  X = 0;
  total += rg1_euler.at(convert(sample.call(X,Y-1),sample.call(Nx-1,Y-1),sample.call(X,0),sample.call(Nx-1,0)));
  
  // Range 2
  X++;
  while(X < Nx)
    {
      total += rg2_euler.at(convert(sample.call(X,Y-1),sample.call(X-1,Y-1),sample.call(X,0),sample.call(X-1,0)));
      X++;
    }
  
  // Range 3
  total += rg3_euler.at(convert(sample.call(0,Y-1),sample.call(Nx-1,Y-1),sample.call(0,0),sample.call(Nx-1,0)));
  
  return total;
}

double w102_xx_pbc(const BinField<bool> &sample)
{
  double total = 0;
  unsigned Nx = sample.call_Nx();
  unsigned Ny = sample.call_Ny();
  
  unsigned X;
  unsigned Y;
  
  
  // -----------------------
  // BOTTOM ROW
  Y = 0;
  
  // Range 7
  X = 0;
  total += rg7_w102_xx.at(convert(sample.call(X,Ny-1),sample.call(Nx-1,Ny-1),sample.call(X,Y),sample.call(Nx-1,Y)));
  
  // Range 8
  X++;
  while(X < Nx)
    {
      total += rg8_w102_xx.at(convert(sample.call(X,Ny-1),sample.call(X-1,Ny-1),sample.call(X,Y),sample.call(X-1,Y)));
      X++;
    }
  
  // Range 9
  total += rg9_w102_xx.at(convert(sample.call(0,Ny-1),sample.call(Nx-1,Ny-1),sample.call(0,Y),sample.call(Nx-1,Y)));
  
  // -----------------------
  // MID ROWS
  Y++;
  while(Y < Ny)
    {
  
      // Range 4
      X = 0;
      total += rg4_w102_xx.at(convert(sample.call(X,Y-1),sample.call(Nx-1,Y-1),sample.call(X,Y),sample.call(Nx-1,Y)));
  
      // Range 5
      X++;
      while(X < Nx)
	{
	  total += rg5_w102_xx.at(convert(sample.call(X,Y-1),sample.call(X-1,Y-1),sample.call(X,Y),sample.call(X-1,Y)));
      	  X++;
	}
  
      // Range 6
      total += rg6_w102_xx.at(convert(sample.call(0,Y-1),sample.call(Nx-1,Y-1),sample.call(0,Y),sample.call(Nx-1,Y)));
  
      
      Y++;  
    }
  
  // -----------------------
  // TOP ROW
  
  // Range 1
  X = 0;
  total += rg1_w102_xx.at(convert(sample.call(X,Y-1),sample.call(Nx-1,Y-1),sample.call(X,0),sample.call(Nx-1,0)));
  
  // Range 2
  X++;
  while(X < Nx)
    {
      total += rg2_w102_xx.at(convert(sample.call(X,Y-1),sample.call(X-1,Y-1),sample.call(X,0),sample.call(X-1,0)));
      X++;
    }
  
  // Range 3
  total += rg3_w102_xx.at(convert(sample.call(0,Y-1),sample.call(Nx-1,Y-1),sample.call(0,0),sample.call(Nx-1,0)));
  
  return total;
}

double w102_xy_pbc(const BinField<bool> &sample)
{
  double total = 0;
  unsigned Nx = sample.call_Nx();
  unsigned Ny = sample.call_Ny();
  
  unsigned X;
  unsigned Y;
  
  
  // -----------------------
  // BOTTOM ROW
  Y = 0;
  
  // Range 7
  X = 0;
  total += rg7_w102_xy.at(convert(sample.call(X,Ny-1),sample.call(Nx-1,Ny-1),sample.call(X,Y),sample.call(Nx-1,Y)));
  
  // Range 8
  X++;
  while(X < Nx)
    {
      total += rg8_w102_xy.at(convert(sample.call(X,Ny-1),sample.call(X-1,Ny-1),sample.call(X,Y),sample.call(X-1,Y)));
      X++;
    }
  
  // Range 9
  total += rg9_w102_xy.at(convert(sample.call(0,Ny-1),sample.call(Nx-1,Ny-1),sample.call(0,Y),sample.call(Nx-1,Y)));
  
  // -----------------------
  // MID ROWS
  Y++;
  while(Y < Ny)
    {
  
      // Range 4
      X = 0;
      total += rg4_w102_xy.at(convert(sample.call(X,Y-1),sample.call(Nx-1,Y-1),sample.call(X,Y),sample.call(Nx-1,Y)));
  
      // Range 5
      X++;
      while(X < Nx)
	{
	  total += rg5_w102_xy.at(convert(sample.call(X,Y-1),sample.call(X-1,Y-1),sample.call(X,Y),sample.call(X-1,Y)));
      	  X++;
	}
  
      // Range 6
      total += rg6_w102_xy.at(convert(sample.call(0,Y-1),sample.call(Nx-1,Y-1),sample.call(0,Y),sample.call(Nx-1,Y)));
  
      
      Y++;  
    }
  
  // -----------------------
  // TOP ROW
  
  // Range 1
  X = 0;
  total += rg1_w102_xy.at(convert(sample.call(X,Y-1),sample.call(Nx-1,Y-1),sample.call(X,0),sample.call(Nx-1,0)));
  
  // Range 2
  X++;
  while(X < Nx)
    {
      total += rg2_w102_xy.at(convert(sample.call(X,Y-1),sample.call(X-1,Y-1),sample.call(X,0),sample.call(X-1,0)));
      X++;
    }
  
  // Range 3
  total += rg3_w102_xy.at(convert(sample.call(0,Y-1),sample.call(Nx-1,Y-1),sample.call(0,0),sample.call(Nx-1,0)));
  
  return total;
}

double w102_yy_pbc(const BinField<bool> &sample)
{
  double total = 0;
  unsigned Nx = sample.call_Nx();
  unsigned Ny = sample.call_Ny();
  
  unsigned X;
  unsigned Y;
  
  
  // -----------------------
  // BOTTOM ROW
  Y = 0;
  
  // Range 7
  X = 0;
  total += rg7_w102_yy.at(convert(sample.call(X,Ny-1),sample.call(Nx-1,Ny-1),sample.call(X,Y),sample.call(Nx-1,Y)));
  
  // Range 8
  X++;
  while(X < Nx)
    {
      total += rg8_w102_yy.at(convert(sample.call(X,Ny-1),sample.call(X-1,Ny-1),sample.call(X,Y),sample.call(X-1,Y)));
      X++;
    }
  
  // Range 9
  total += rg9_w102_yy.at(convert(sample.call(0,Ny-1),sample.call(Nx-1,Ny-1),sample.call(0,Y),sample.call(Nx-1,Y)));
  
  // -----------------------
  // MID ROWS
  Y++;
  while(Y < Ny)
    {
  
      // Range 4
      X = 0;
      total += rg4_w102_yy.at(convert(sample.call(X,Y-1),sample.call(Nx-1,Y-1),sample.call(X,Y),sample.call(Nx-1,Y)));
  
      // Range 5
      X++;
      while(X < Nx)
	{
	  total += rg5_w102_yy.at(convert(sample.call(X,Y-1),sample.call(X-1,Y-1),sample.call(X,Y),sample.call(X-1,Y)));
      	  X++;
	}
  
      // Range 6
      total += rg6_w102_yy.at(convert(sample.call(0,Y-1),sample.call(Nx-1,Y-1),sample.call(0,Y),sample.call(Nx-1,Y)));
  
      
      Y++;  
    }
  
  // -----------------------
  // TOP ROW
  
  // Range 1
  X = 0;
  total += rg1_w102_yy.at(convert(sample.call(X,Y-1),sample.call(Nx-1,Y-1),sample.call(X,0),sample.call(Nx-1,0)));
  
  // Range 2
  X++;
  while(X < Nx)
    {
      total += rg2_w102_yy.at(convert(sample.call(X,Y-1),sample.call(X-1,Y-1),sample.call(X,0),sample.call(X-1,0)));
      X++;
    }
  
  // Range 3
  total += rg3_w102_yy.at(convert(sample.call(0,Y-1),sample.call(Nx-1,Y-1),sample.call(0,0),sample.call(Nx-1,0)));
  
  return total;
}

double delta_pbc(const BinField<bool> &sample)
{
  double w102_xx_ = w102_xx_pbc(sample);
  double w102_xy_ = w102_xy_pbc(sample);
  double w102_yy_ = w102_yy_pbc(sample);
  
  return sqrt((w102_xx_-w102_yy_)*(w102_xx_-w102_yy_)+4*w102_xy_*w102_xy_);
}

// Look-up table but pixelized data: all functional values times 8:
int area_mbc_pix(const BinField<bool> &sample)
{
  int total = 0;
  unsigned Nx = sample.call_Nx();
  unsigned Ny = sample.call_Ny();
  
  unsigned X;
  unsigned Y;
  
  
  // -----------------------
  // BOTTOM ROW
  Y = 0;
  
  // Range 7
  X = 0;
  total += rg7_area_pix.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  // Range 8
  X++;
  while(X < (Nx-2))
    {
      total += rg8_area_pix.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
      
      X++;
    }
  
  // Range 9
  total += rg9_area_pix.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  // -----------------------
  // MID ROWS
  Y++;
  while(Y < (Ny-2))
    {
  
      // Range 4
      X = 0;
      total += rg4_area_pix.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
      // Range 5
      X++;
      while(X < (Nx-2))
	{
	  total += rg5_area_pix.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
      
	  X++;
	}
  
      // Range 6
      total += rg6_area_pix.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
      
      Y++;  
    }
  
  // -----------------------
  // TOP ROW
  
  // Range 1
  X = 0;
  total += rg1_area_pix.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  // Range 2
  X++;
  while(X < (Nx-2))
    {
      total += rg2_area_pix.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
      
      X++;
    }
  
  // Range 3
  total += rg3_area_pix.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  return total;
}

int perimeter_mbc_pix(const BinField<bool> &sample)
{
  int total = 0;
  unsigned Nx = sample.call_Nx();
  unsigned Ny = sample.call_Ny();
  
  unsigned X;
  unsigned Y;
  
  
  // -----------------------
  // BOTTOM ROW
  Y = 0;
  
  // Range 7
  X = 0;
  total += rg7_perimeter_pix.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  // Range 8
  X++;
  while(X < (Nx-2))
    {
      total += rg8_perimeter_pix.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
      
      X++;
    }
  
  // Range 9
  total += rg9_perimeter_pix.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  // -----------------------
  // MID ROWS
  Y++;
  while(Y < (Ny-2))
    {
  
      // Range 4
      X = 0;
      total += rg4_perimeter_pix.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
      // Range 5
      X++;
      while(X < (Nx-2))
	{
	  total += rg5_perimeter_pix.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
      
	  X++;
	}
  
      // Range 6
      total += rg6_perimeter_pix.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
      
      Y++;  
    }
  
  // -----------------------
  // TOP ROW
  
  // Range 1
  X = 0;
  total += rg1_perimeter_pix.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  // Range 2
  X++;
  while(X < (Nx-2))
    {
      total += rg2_perimeter_pix.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
      
      X++;
    }
  
  // Range 3
  total += rg3_perimeter_pix.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  return total;
}

int euler_mbc_pix(const BinField<bool> &sample)
{
  int total = 0;
  unsigned Nx = sample.call_Nx();
  unsigned Ny = sample.call_Ny();
  
  unsigned X;
  unsigned Y;
  
  
  // -----------------------
  // BOTTOM ROW
  Y = 0;
  
  // Range 7
  X = 0;
  total += rg7_euler_pix.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  // Range 8
  X++;
  while(X < (Nx-2))
    {
      total += rg8_euler_pix.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
      
      X++;
    }
  
  // Range 9
  total += rg9_euler_pix.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  // -----------------------
  // MID ROWS
  Y++;
  while(Y < (Ny-2))
    {
  
      // Range 4
      X = 0;
      total += rg4_euler_pix.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
      // Range 5
      X++;
      while(X < (Nx-2))
	{
	  total += rg5_euler_pix.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
      
	  X++;
	}
  
      // Range 6
      total += rg6_euler_pix.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
      
      Y++;  
    }
  
  // -----------------------
  // TOP ROW
  
  // Range 1
  X = 0;
  total += rg1_euler_pix.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  // Range 2
  X++;
  while(X < (Nx-2))
    {
      total += rg2_euler_pix.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
      
      X++;
    }
  
  // Range 3
  total += rg3_euler_pix.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
  
  return total;
}

// White boundary condition
int area_wbc_pix(const BinField<bool> &sample)
{
  int total = 0;
  unsigned Nx = sample.call_Nx();
  unsigned Ny = sample.call_Ny();
  
  unsigned X;
  unsigned Y;
  
  // -----------------------
  // BOTTOM ROW
  
  // Range 7
  total += rg7_area_pix.at(convert(false,false,sample.call(0,0),false));
  
  // Range 8
  X=0;
  while(X < (Nx-1))
    {
      total += rg8_area_pix.at(convert(false,false,sample.call(X+1,0),sample.call(X,0)));
      
      X++;
    }
  
  // Range 9
  total += rg9_area_pix.at(convert(false,false,false,sample.call(Nx-1,0)));
  
  // -----------------------
  // MID ROWS
  Y=0;
  while(Y < (Ny-1))
    {
  
      // Range 4
      total += rg4_area_pix.at(convert(sample.call(0,Y),false,sample.call(0,Y+1),false));
  
      // Range 5
      X=0;
      while(X < (Nx-1))
	{
	  total += rg5_area_pix.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
      
	  X++;
	}
  
      // Range 6
      total += rg6_area_pix.at(convert(false,sample.call(Nx-1,Y),false,sample.call(Nx-1,Y+1)));
      
      Y++;  
    }
  
  // -----------------------
  // TOP ROW
  
  // Range 1
  total += rg1_area_pix.at(convert(sample.call(0,Ny-1),false,false,false));
  
  // Range 2
  X=0;
  while(X < (Nx-1))
    {
      total += rg2_area_pix.at(convert(sample.call(X+1,Ny-1),sample.call(X,Ny-1),false,false));
      
      X++;
    }
  
  // Range 3
  total += rg3_area_pix.at(convert(false,sample.call(Nx-1,Ny-1),false,false));
  
  return total;
}

int perimeter_wbc_pix(const BinField<bool> &sample)
{
  int total = 0;
  unsigned Nx = sample.call_Nx();
  unsigned Ny = sample.call_Ny();
  
  unsigned X;
  unsigned Y;
  
  // -----------------------
  // BOTTOM ROW
  
  // Range 7
  total += rg7_perimeter_pix.at(convert(false,false,sample.call(0,0),false));
  
  // Range 8
  X=0;
  while(X < (Nx-1))
    {
      total += rg8_perimeter_pix.at(convert(false,false,sample.call(X+1,0),sample.call(X,0)));
      
      X++;
    }
  
  // Range 9
  total += rg9_perimeter_pix.at(convert(false,false,false,sample.call(Nx-1,0)));
  
  // -----------------------
  // MID ROWS
  Y=0;
  while(Y < (Ny-1))
    {
  
      // Range 4
      total += rg4_perimeter_pix.at(convert(sample.call(0,Y),false,sample.call(0,Y+1),false));
  
      // Range 5
      X=0;
      while(X < (Nx-1))
	{
	  total += rg5_perimeter_pix.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
      
	  X++;
	}
  
      // Range 6
      total += rg6_perimeter_pix.at(convert(false,sample.call(Nx-1,Y),false,sample.call(Nx-1,Y+1)));
      
      Y++;  
    }
  
  // -----------------------
  // TOP ROW
  
  // Range 1
  total += rg1_perimeter_pix.at(convert(sample.call(0,Ny-1),false,false,false));
  
  // Range 2
  X=0;
  while(X < (Nx-1))
    {
      total += rg2_perimeter_pix.at(convert(sample.call(X+1,Ny-1),sample.call(X,Ny-1),false,false));
      
      X++;
    }
  
  // Range 3
  total += rg3_perimeter_pix.at(convert(false,sample.call(Nx-1,Ny-1),false,false));
  
  return total;
}

int euler_wbc_pix(const BinField<bool> &sample)
{
  int total = 0;
  unsigned Nx = sample.call_Nx();
  unsigned Ny = sample.call_Ny();
  
  unsigned X;
  unsigned Y;
  
  // -----------------------
  // BOTTOM ROW
  
  // Range 7
  total += rg7_euler_pix.at(convert(false,false,sample.call(0,0),false));
  
  // Range 8
  X=0;
  while(X < (Nx-1))
    {
      total += rg8_euler_pix.at(convert(false,false,sample.call(X+1,0),sample.call(X,0)));
      
      X++;
    }
  
  // Range 9
  total += rg9_euler_pix.at(convert(false,false,false,sample.call(Nx-1,0)));
  
  // -----------------------
  // MID ROWS
  Y=0;
  while(Y < (Ny-1))
    {
  
      // Range 4
      total += rg4_euler_pix.at(convert(sample.call(0,Y),false,sample.call(0,Y+1),false));
  
      // Range 5
      X=0;
      while(X < (Nx-1))
	{
	  total += rg5_euler_pix.at(convert(sample.call(X+1,Y),sample.call(X,Y),sample.call(X+1,Y+1),sample.call(X,Y+1)));
      
	  X++;
	}
  
      // Range 6
      total += rg6_euler_pix.at(convert(false,sample.call(Nx-1,Y),false,sample.call(Nx-1,Y+1)));
      
      Y++;  
    }
  
  // -----------------------
  // TOP ROW
  
  // Range 1
  total += rg1_euler_pix.at(convert(sample.call(0,Ny-1),false,false,false));
  
  // Range 2
  X=0;
  while(X < (Nx-1))
    {
      total += rg2_euler_pix.at(convert(sample.call(X+1,Ny-1),sample.call(X,Ny-1),false,false));
      
      X++;
    }
  
  // Range 3
  total += rg3_euler_pix.at(convert(false,sample.call(Nx-1,Ny-1),false,false));
  
  return total;
}

// PIXELIZED DATA
// Minus sampling boundary condition
double area_pixelized_mbc(const BinField<bool> &sample)
{
  double total = 0;
  
  for(unsigned xi = 1; xi < (sample.call_Nx()-1); xi++)
    for(unsigned yi = 1; yi < (sample.call_Ny()-1); yi++)
      total += sample.call(xi,yi);

  return total;
}

double perimeter_pixelized_mbc(const BinField<bool> &sample)
{
  double total = 0;
  unsigned Nx = sample.call_Nx();
  unsigned Ny = sample.call_Ny();

  unsigned xi = 0;
  unsigned yi = 0;
  for(yi = 1; yi < (Ny-1); yi++)
    total += ( (sample.call(xi,yi)==false) && sample.call(xi+1,yi) );
  
  xi = Nx-2;
  for(yi = 1; yi < (Ny-1); yi++)
    total += ( (sample.call(xi+1,yi)==false) && sample.call(xi,yi) );
  
  for(xi = 1; xi < (Nx-2); xi++)
    for(yi = 1; yi < (Ny-1); yi++)
      total += ( (sample.call(xi+1,yi)==false) && sample.call(xi,yi) ) + ( (sample.call(xi,yi)==false) && sample.call(xi+1,yi) );
  
  yi = 0;
  for(xi = 1; xi < (Nx-1); xi++)
    total += ( (sample.call(xi,yi)==false) && sample.call(xi,yi+1) );
  
  yi = Ny-2;
  for(xi = 1; xi < (Nx-1); xi++)
    total += ( (sample.call(xi,yi+1)==false) && sample.call(xi,yi) );
  
  for(yi = 1; yi < (Ny-2); yi++)
    for(xi = 1; xi < (Nx-1); xi++)
    total += ( (sample.call(xi,yi)==false) && sample.call(xi,yi+1) ) + ( (sample.call(xi,yi+1)==false) && sample.call(xi,yi) );
  
  return total;
}

double euler_pixelized_mbc(const BinField<bool> &sample)
{
  return euler_mbc(sample);
}

// Periodic boundary condition
double area_pixelized_pbc(const BinField<bool> &sample)
{
  double total = 0;
  
  for(unsigned xi = 0; xi < sample.call_Nx(); xi++)
    for(unsigned yi = 0; yi < sample.call_Ny(); yi++)
      total += sample.call(xi,yi);
  
  return total;
}

double perimeter_pixelized_pbc(const BinField<bool> &sample)
{
  double total = 0;
  unsigned Nx = sample.call_Nx();
  unsigned Ny = sample.call_Ny();

  unsigned xi = 0;
  unsigned yi = 0;
  for(yi = 0; yi < Ny; yi++){
    total += ( (sample.call(Nx-1,yi)==false) && sample.call(0,yi) );
    total += ( (sample.call(0,yi)==false) && sample.call(Nx-1,yi) );
  }
  
  for(xi = 1; xi < Nx; xi++)
    for(yi = 0; yi < Ny; yi++){
      total += ( (sample.call(xi-1,yi)==false) && sample.call(xi,yi) );
      total += ( (sample.call(xi,yi)==false) && sample.call(xi-1,yi) );
    }
  
  for(xi = 0; xi < Nx; xi++){
    total += ( (sample.call(xi,Ny-1)==false) && sample.call(xi,0) );
    total += ( (sample.call(xi,0)==false) && sample.call(xi,Ny-1) );
  }

  for(yi = 1; yi < Ny; yi++)
    for(xi = 0; xi < Nx; xi++){
      total += ( (sample.call(xi,yi)==false) && sample.call(xi,yi-1) );
      total += ( (sample.call(xi,yi-1)==false) && sample.call(xi,yi) );
    }
  
  return total;
}

double euler_pixelized_pbc(const BinField<bool> &sample)
{
  return euler_pbc(sample);
}

// PAPAYA
void print_pgm(const BinField<bool> &sample,
	       const std::string &filename, const std::string &prefix_of, const bool &invert)
{
  unsigned Nx = sample.call_Nx();
  unsigned Ny = sample.call_Ny();
  
  std::ofstream OutFile( (prefix_of + filename).c_str() );
  OutFile << "P2" << std::endl;
  OutFile << Nx << " " << Ny << std::endl;
  OutFile << "255" << std::endl;

  unsigned truephase = 255, voidphase = 0;
  if(invert){
    truephase = 0; voidphase = 255;}
  
  for( unsigned j = (Ny-1); j < Ny; j--)
    for( unsigned i = 0; i < Nx; i++)
      {
	if( sample.call(i,j) )
	  OutFile << truephase << std::endl;
	else
	  OutFile << voidphase << std::endl;
      }
  
  OutFile.close();
}

// ---------------------------------------------------------------
