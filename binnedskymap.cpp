#include "binnedskymap.h"

// --------------------------------------
// Constructors
BinnedSkyMap::BinnedSkyMap(const unsigned &Nx, const unsigned &Ny,
			   const double &xmin, const double &xmax,
			   const double &ymin, const double &ymax,
			   const std::string &prefix_if, const std::string &prefix_of,
			   const unsigned &seed, const unsigned &seedsub) :
  Nx_ (Nx), Ny_ (Ny), xmin_ (xmin), xmax_ (xmax), ymin_ (ymin), ymax_ (ymax), xbinlength_ ((xmax_ - xmin_)/Nx_), ybinlength_ ((ymax_ - ymin_)/Ny_), xlow_ (BinField< double >(Nx_, Ny_, 0)), ylow_ (BinField< double >(Nx_, Ny_, 0)), intensity_ (BinField< double > (Nx_, Ny_, 0)), countrate_ (BinField< unsigned > (Nx_, Ny_, 0)), srcsub_intensity_ (BinField< double > (Nx_, Ny_, 0)), prefix_if_ (prefix_if), prefix_of_ (prefix_of), seed_ (seed), seedsub_ (seedsub), trial (200)
{
  for(unsigned xi = 0; xi < Nx_; xi++)
    for(unsigned yi = 0; yi < Ny_; yi++)
      {
	xlow_.assign(xi,yi,xmin_ + xi*xbinlength_);
	ylow_.assign(xi,yi,ymin_ + yi*ybinlength_);
      }
}

BinnedSkyMap::BinnedSkyMap(const unsigned &N,
			   const double &xmin, const double &xmax,
			   const double &ymin, const double &ymax,
			   const std::string &prefix_if, const std::string &prefix_of,
			   const unsigned &seed, const unsigned &seedsub) :
  Nx_ (N), Ny_ (N), xmin_ (xmin), xmax_ (xmax), ymin_ (ymin), ymax_ (ymax), xbinlength_ ((xmax_ - xmin_)/Nx_), ybinlength_ ((ymax_ - ymin_)/Ny_), xlow_ (BinField< double >(Nx_, Ny_, 0)), ylow_ (BinField< double >(Nx_, Ny_, 0)), intensity_ (BinField< double > (Nx_, Ny_, 0)), countrate_ (BinField< unsigned > (Nx_, Ny_, 0)), srcsub_intensity_ (BinField< double > (Nx_, Ny_, 0)), prefix_if_ (prefix_if), prefix_of_ (prefix_of), seed_ (seed), seedsub_ (seedsub), trial (200)
{
  for(unsigned xi = 0; xi < Nx_; xi++)
    for(unsigned yi = 0; yi < Ny_; yi++)
      {
	xlow_.assign(xi,yi,xmin_ + xi*xbinlength_);
	ylow_.assign(xi,yi,ymin_ + yi*ybinlength_);
      }
}

BinnedSkyMap::BinnedSkyMap(const std::string &countratefilename,
			   const double &xmin, const double &xmax,
			   const double &ymin, const double &ymax,
			   const std::string &prefix_if, const std::string &prefix_of,
			   const unsigned &seed, const unsigned &seedsub) :
  Nx_ (CheckXDimensionOfMatrixFromFile(countratefilename,prefix_if)), Ny_ (CheckYDimensionOfMatrixFromFile(countratefilename,prefix_if)), xmin_ (xmin), xmax_ (xmax), ymin_ (ymin), ymax_ (ymax), xbinlength_ ((xmax_ - xmin_)/Nx_), ybinlength_ ((ymax_ - ymin_)/Ny_), xlow_ (BinField< double >(Nx_, Ny_, 0)), ylow_ (BinField< double >(Nx_, Ny_, 0)), intensity_ (BinField< double > (Nx_, Ny_, 0)), countrate_ (BinField< unsigned > (countratefilename,prefix_if)), srcsub_intensity_ (BinField< double > (Nx_, Ny_, 0)), prefix_if_ (prefix_if), prefix_of_ (prefix_of), seed_ (seed), seedsub_ (seedsub), trial (200)
{
  for(unsigned xi = 0; xi < Nx_; xi++)
    for(unsigned yi = 0; yi < Ny_; yi++)
      {
	xlow_.assign(xi,yi,xmin_ + xi*xbinlength_);
	ylow_.assign(xi,yi,ymin_ + yi*ybinlength_);
      }
}
/*
BinnedSkyMap::BinnedSkyMap(const std::string &eventlistfilename,
			   const double &xmin, const double &xmax,
			   const double &ymin, const double &ymax,
			   const std::string &prefix_if, const std::string &prefix_of,
			   const unsigned &seed, const unsigned &seedsub) :
  Nx_ (CheckXDimensionOfMatrixFromFile(eventlistfilename,prefix_if)), Ny_ (CheckYDimensionOfMatrixFromFile(eventlistfilename,prefix_if)), xmin_ (xmin), xmax_ (xmax), ymin_ (ymin), ymax_ (ymax), xbinlength_ ((xmax_ - xmin_)/Nx_), ybinlength_ ((ymax_ - ymin_)/Ny_), xlow_ (BinField< double >(Nx_, Ny_, 0)), ylow_ (BinField< double >(Nx_, Ny_, 0)), intensity_ (BinField< double > (Nx_, Ny_, 0)), countrate_ (BinField< unsigned > (eventlistfilename,prefix_if)), srcsub_intensity_ (BinField< double > (Nx_, Ny_, 0)), prefix_if_ (prefix_if), prefix_of_ (prefix_of), seed_ (seed), seedsub_ (seedsub)
{
  for(unsigned xi = 0; xi < Nx_; xi++)
    for(unsigned yi = 0; yi < Ny_; yi++)
      {
	xlow_.assign(xi,yi,xmin_ + xi*xbinlength_);
	ylow_.assign(xi,yi,ymin_ + yi*ybinlength_);
      }
}
*/
// --------------------------------------

// --------------------------------------
// Reset
void BinnedSkyMap::reset()
{
  for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
    { 
      intensity_.assign(ni,0);
      countrate_.assign(ni,0);
      // sources_ = Gaussian();
    }
}
// --------------------------------------

// --------------------------------------
// Manipulate intensity
void BinnedSkyMap::scale_intensity(const double &factor)
{
  for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
    intensity_.assign(ni,intensity_.call(ni)*factor);
}

void BinnedSkyMap::add_bck_intensity(const double &k_bck)
{
  for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
    intensity_.assign(ni,intensity_.call(ni) + k_bck);  
}

void BinnedSkyMap::add_intensity_to_single_bin(const double &lambda, const unsigned &xi, const unsigned &yi)
{
  intensity_.assign(xi, yi,intensity_.call(xi,yi) + lambda);  
}

void BinnedSkyMap::add_src_intensity(const double &k_src, const Gaussian &sources, const bool &verbose)
{
  // Define sources
  sources_ = sources;
  
  // Auxiliary parameters
  double Lambda_i = 0;
  double c1,c2,c3,c4,c5;
  double min_g=1e-14;
  
  if(verbose)
    std::cout << "Intensity computation:" << std::endl
	      << PercentageBar << std::endl;
  unsigned p = 0;
  
  for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
    {
      rx_ = xlow_.call(ni);
      c1 = gauss_pdf(ylow_.call(ni), this);
      
      rx_ = xlow_.call(ni) + xbinlength_;
      c2 = gauss_pdf(ylow_.call(ni), this);
      
      rx_ = xlow_.call(ni);
      c3 = gauss_pdf(ylow_.call(ni) + ybinlength_, this);
      
      rx_ = xlow_.call(ni) + xbinlength_;
      c4 = gauss_pdf(ylow_.call(ni) + ybinlength_, this);
      
      rx_ = xlow_.call(ni) + xbinlength_/2;
      c5 = gauss_pdf(ylow_.call(ni) + ybinlength_/2, this);
      
      if(c1 > min_g && c2 > min_g && c3 > min_g && c4 > min_g && c5 > min_g)
	{
	  Lambda_i = Lambda(k_src, 
			    xlow_.call(ni), xlow_.call(ni) + xbinlength_, 
			    ylow_.call(ni), ylow_.call(ni) + ybinlength_);
	  intensity_.assign(ni, intensity_.call(ni) + Lambda_i);
	}
      if(verbose)
	Cout_One_Percent(ni,p,Nx_*Ny_);
    }
  
}
// --------------------------------------
// Calculate intensity
// Auxiliary function and doubles
double BinnedSkyMap::Lambda(const double &k_src,
			    const double &x_low, const double &x_up,
			    const double &y_low, const double &y_up)
{  
  single_y_low_ = y_low;
  single_y_up_  = y_up;
  
  gsl_function F_margin;
  F_margin.function = &margin_gauss_pdf;
  F_margin.params   = this;
  
  double result;
  double estimated_abs_error;
  size_t nsteps;
  
  gsl_integration_qng (&F_margin, x_low, x_up, 
		       demand_abs_error_, demand_rel_error_, 
		       &result, &estimated_abs_error, &nsteps);
  
  return k_src*result;
}

double margin_gauss_pdf(double x, void * params)
{
  BinnedSkyMap* smap = (BinnedSkyMap*) params;
  
  (*smap).rx_ = x;
  
  gsl_function F;
  F.function = &gauss_pdf;
  F.params   = params;
  
  double result;
  double estimated_abs_error;
  size_t nsteps;
    
  gsl_integration_qng (&F, (*smap).single_y_low_, (*smap).single_y_up_, 
		       demand_abs_error_, demand_rel_error_, 
		       &result, &estimated_abs_error, &nsteps);
  
  return result;
}

double gauss_pdf(double y, void * params)
{
  // Read in parameters
  BinnedSkyMap* smap = (BinnedSkyMap*) params;
  Gaussian sources = (*smap).sources_;
  double x = (*smap).rx_;
  
  // Function value
  double f = 0;
  
  // Construct function
  for(unsigned si = 0; si < sources.number_of_sources(); si++)
    {
      double varx = pow(sources.sigma_xx(si),2), vary = pow(sources.sigma_yy(si),2);
      double cp = cos(sources.phi(si)), sp = sin(sources.phi(si));
      double nx = x - sources.mu_x(si), ny = y - sources.mu_y(si);
      double a = varx*cp*cp + vary*sp*sp, b = (varx-vary)*cp*sp, d = varx*sp*sp+vary*cp*cp;
      f += (sources.ratio(si))/(2*M_PI*sqrt(varx*vary))*exp(-(nx*nx*d+ny*ny*a-2*nx*ny*b)/(2.0*varx*vary));
    }
  
  return f;
}
// --------------------------------------

// --------------------------------------
// Simulation of a countrate-sample with respect to the given intensity
void BinnedSkyMap::add_counts_to_single_bin(const unsigned &k, const unsigned &xi, const unsigned &yi)
{
  countrate_.assign(xi, yi, countrate_.call(xi,yi) + k);
}

void BinnedSkyMap::simulate_countrate()
{
  RandomPoissonBinField(countrate_, intensity_, seed_);
}

void BinnedSkyMap::global_reduction(const double &f)
{
  for(unsigned xi = 0; xi < Nx_; xi++)
    for(unsigned yi = 0; yi < Ny_; yi++)
      countrate_.assign( xi, yi, RandomBinomial(countrate_.call(xi, yi), f, seedsub_) );
}

void BinnedSkyMap::add_global_Poisson(const double &lambda)
{
  for(unsigned xi = 0; xi < Nx_; xi++)
    for(unsigned yi = 0; yi < Ny_; yi++)
      countrate_.add( xi, yi, RandomPoisson(lambda, seedsub_));
}

void BinnedSkyMap::reduction_because_of_detector_acceptance()
{
  for(unsigned xi = 0; xi < Nx_; xi++)
    for(unsigned yi = 0; yi < Ny_; yi++)
      {
	double f_i = exp( (- pow(xi - 0.5*Nx_,2)/pow(Nx_*0.3,2) - pow(yi - 0.5*Ny_,2)/pow(Ny_*0.3,2) )*0.5 );
	
	countrate_.assign( xi, yi, RandomBinomial(countrate_.call(xi, yi), f_i, seedsub_) );
      }
}

void BinnedSkyMap::global_detector_correction(const double &target_f, const double &lambda)
{
  for(unsigned xi = 0; xi < Nx_; xi++)
    for(unsigned yi = 0; yi < Ny_; yi++)
      {
	double f_i = exp( (- pow(xi - 0.5*Nx_,2)/pow(Nx_*0.3,2) - pow(yi - 0.5*Ny_,2)/pow(Ny_*0.3,2) )*0.5 );
	
	if(f_i < target_f){
	  // Poisson fill
	  countrate_.add( xi, yi, RandomPoisson((target_f-f_i)*lambda, seedsub_));
	}
	else if(f_i > target_f){
	  // Post selection
	  countrate_.assign( xi, yi, RandomBinomial(countrate_.call(xi, yi), target_f/f_i, seedsub_) );
	}
	
      }
}

std::vector< BinField< unsigned > > BinnedSkyMap::local_detector_correction(const double &lambda, const unsigned &slide) const
{
  unsigned Nx_mink = Nx_-slide+1;
  unsigned Ny_mink = Ny_-slide+1;
  unsigned center = slide/2;
  
  BinField< unsigned > single_count_map(slide,slide,0);
  std::vector< BinField< unsigned > > countrate_locally_corrected(Nx_mink*Ny_mink,single_count_map);
  
  unsigned ni = 0;
  for(unsigned xi = 0; xi < Nx_mink; xi++)
    for(unsigned yi = 0; yi < Ny_mink; yi++){
      double target_f = exp( (- pow(xi+center - 0.5*Nx_,2)/pow(Nx_*0.3,2) - pow(yi+center - 0.5*Ny_,2)/pow(Ny_*0.3,2) )*0.5 );
      
      // define for each pixel a binfield and use as entry in vector
      // iterate for each binfield over all entries: set to detector corrected bin field (calculate f_i therefore)
      for(unsigned xs = 0; xs < slide; xs++)
	for(unsigned ys = 0; ys < slide; ys++)
	  if ((xs == center) && (ys == center)){
	    countrate_locally_corrected[ni].assign(xs,ys, countrate_.call(xi+xs,yi+ys) );
	  }
	  else{
	    double f_i = exp( (- pow(xi+xs - 0.5*Nx_,2)/pow(Nx_*0.3,2) - pow(yi+ys - 0.5*Ny_,2)/pow(Ny_*0.3,2) )*0.5 );
	    if(f_i < target_f){
	      // Poisson fill
	      countrate_locally_corrected[ni].assign(xs,ys, countrate_.call(xi+xs,yi+ys) + RandomPoisson((target_f-f_i)*lambda, seedsub_) );
	    }
	    else if(f_i > target_f){
	      // Post selection
	      countrate_locally_corrected[ni].assign(xs,ys, RandomBinomial(countrate_.call(xi+xs,yi+ys), target_f/f_i, seedsub_) );
	    }
	    else{
	      countrate_locally_corrected[ni].assign(xs,ys, countrate_.call(xi+xs,yi+ys) );
	    }
	  }
      
      ni++;
    }
  
  return countrate_locally_corrected;
}
// --------------------------------------

// --------------------------------------
// Calling and print out functions
unsigned BinnedSkyMap::call_Nx() const{
  return Nx_;}
unsigned BinnedSkyMap::call_Ny() const{
  return Ny_;}
double BinnedSkyMap::call_xlow(const unsigned &xi, const unsigned &yi) const{
  return xlow_.call(xi,yi);}
double BinnedSkyMap::call_xlow(const unsigned &ni) const{
  return xlow_.call(ni);}
double BinnedSkyMap::call_ylow(const unsigned &xi, const unsigned &yi) const{
  return ylow_.call(xi,yi);}
double BinnedSkyMap::call_ylow(const unsigned &ni) const{
  return ylow_.call(ni);}
double BinnedSkyMap::call_intensity(const unsigned &xi, const unsigned &yi) const{
  return intensity_.call(xi,yi);}
double BinnedSkyMap::call_intensity(const unsigned &ni) const{
  return intensity_.call(ni);}
unsigned BinnedSkyMap::call_countrate(const unsigned &xi, const unsigned &yi) const{
  return countrate_.call(xi,yi);}
unsigned BinnedSkyMap::call_countrate(const unsigned &ni) const{
  return countrate_.call(ni);}
BinField<unsigned> BinnedSkyMap::count_map() const{
  return countrate_;}
  
double BinnedSkyMap::call_xbinlength() const{
  return xbinlength_;}
double BinnedSkyMap::call_ybinlength() const{
  return ybinlength_;}
BinField<double> BinnedSkyMap::call_xlow() const{
  return xlow_;}
BinField<double> BinnedSkyMap::call_ylow() const{
  return ylow_;}

void BinnedSkyMap::xout_xlow(const unsigned &xlow, const unsigned &xup, const unsigned &ylow, const unsigned &yup) const{
  xlow_.xout(xlow, xup, ylow, yup);}
void BinnedSkyMap::xout_ylow(const unsigned &xlow, const unsigned &xup, const unsigned &ylow, const unsigned &yup) const{
  ylow_.xout(xlow, xup, ylow, yup);}
void BinnedSkyMap::xout_intensity(const unsigned &xlow,const unsigned &xup,const unsigned &ylow,const unsigned &yup) const{
  intensity_.xout(xlow, xup, ylow, yup);}
void BinnedSkyMap::xout_countrate(const unsigned &xlow,const unsigned &xup,const unsigned &ylow,const unsigned &yup) const{
  countrate_.xout(xlow, xup, ylow, yup);}

void BinnedSkyMap::fout_xlow(const std::string &filename) const{
  xlow_.fout(filename,prefix_of_);}
void BinnedSkyMap::fout_ylow(const std::string &filename) const{
  ylow_.fout(filename,prefix_of_);}
void BinnedSkyMap::fout_intensity(const std::string &filename) const{
  intensity_.fout(filename,prefix_of_);}
void BinnedSkyMap::fout_countrate(const std::string &filename) const{
  countrate_.fout(filename,prefix_of_);}

void BinnedSkyMap::matrix_out_intensity(const std::string &filename) const{
  intensity_.MatrixToFile(filename,prefix_of_);}
void BinnedSkyMap::matrix_out_countrate(const std::string &filename) const{
  countrate_.MatrixToFile(filename,prefix_of_);}

void BinnedSkyMap::plot_intensity(const std::string &filename) const{
  std::stringstream filename_stst;
  filename_stst << prefix_of_ << filename
		<< "_" << Nx_ << "x" << Ny_ << "_seed_" << seed_ << ".dat";
  std::string filename_st = filename_stst.str();
  std::ofstream Data( filename_st.c_str() );
  
  Data << "# Intensity of a binned sky map" << std::endl
       << "#" << std::setw(13) << "X" << std::setw(14) << "Y" << std::setw(14) << "Lambda" << std::endl;
  for(unsigned xi = 0; xi < Nx_; xi++)
    {
      for(unsigned yi = 0; yi < Ny_; yi++)
	Data << std::setw(14) << xlow_.call(xi,yi)
	     << std::setw(14) << ylow_.call(xi,yi)
	     << std::setw(14) << intensity_.call(xi,yi)
	     << std::endl
	     << std::setw(14) << xlow_.call(xi,yi)
	     << std::setw(14) << ylow_.call(xi,yi) + ybinlength_
	     << std::setw(14) << intensity_.call(xi,yi)
	     << std::endl;
      
      Data << std::endl;
      
      for(unsigned yi = 0; yi < Ny_; yi++)
	Data << std::setw(14) << xlow_.call(xi,yi) + xbinlength_
	     << std::setw(14) << ylow_.call(xi,yi)
	     << std::setw(14) << intensity_.call(xi,yi)
	     << std::endl
	     << std::setw(14) << xlow_.call(xi,yi) + xbinlength_
	     << std::setw(14) << ylow_.call(xi,yi) + ybinlength_
	     << std::setw(14) << intensity_.call(xi,yi)
	     << std::endl;
      
      Data << std::endl;
    }
  Data.close();
  
  std::stringstream gpfilename_stst;
  gpfilename_stst << prefix_of_ << filename
		  << "_" << Nx_ << "x" << Ny_ << "_seed_" << seed_ << ".gp";
  std::string gpfilename_st = gpfilename_stst.str();
  std::ofstream Gnuplot( gpfilename_st.c_str() );
  
  Gnuplot << "set terminal epslatex standalone color colortext 14" << std::endl
	  << "set output '" << filename << "_" << Nx_ << "x" << Ny_ << "_seed_" << seed_ << ".tex'"
	  << std::endl << std::endl
	  << "xmin = " << xmin_ << std::endl
	  << "xmax = " << xmax_ << std::endl
	  << "ymin = " << ymin_ << std::endl
	  << "ymax = " << ymax_ << std::endl
	  << std::endl
	  << "unset xtics #border mirror norotate offset 0,0.3" << std::endl
	  << "unset ytics #border mirror norotate offset 0.3,0" << std::endl
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
	  << "unset xlabel #'$x$' offset 0,0.3" << std::endl
	  << "unset ylabel #'$y$' offset 0.3,0" << std::endl
	  << "set cblabel '\\large $\\lambda$' offset 0,0" << std::endl
	  << std::endl
	  << "set xrange [xmin:xmax]" << std::endl
	  << "set yrange [ymin:ymax]" << std::endl
	  << std::endl
	  << "set pm3d map" << std::endl
	  << "# set palette defined (0 'black', 3 'dark-blue', "
	  << "4 'purple', 6 'dark-red', 9 'red', 12 'orange', 15 'yellow')" << std::endl
	  << "splot '" << filename << "_" << Nx_ << "x" << Ny_ << "_seed_" << seed_ << ".dat" << "' notitle" << std::endl
	  << "#" << std::endl
	  << std::endl;
  
  Gnuplot.close(); 
}
void BinnedSkyMap::plot_countrate(const std::string &filename) const{
  std::stringstream filename_stst;
  filename_stst << prefix_of_ << filename
		<< "_" << Nx_ << "x" << Ny_ << "_seed_" << seed_ << ".dat";
  std::string filename_st = filename_stst.str();
  std::ofstream Data( filename_st.c_str() );
  
  Data << "# Count rate of a binned sky map" << std::endl
       << "#" << std::setw(13) << "X" << std::setw(14) << "Y" << std::setw(14) << "K" << std::endl;
  for(unsigned xi = 0; xi < Nx_; xi++)
    {
      for(unsigned yi = 0; yi < Ny_; yi++)
	Data << std::setw(14) << xlow_.call(xi,yi)
	     << std::setw(14) << ylow_.call(xi,yi)
	     << std::setw(14) << countrate_.call(xi,yi)
	     << std::endl
	     << std::setw(14) << xlow_.call(xi,yi)
	     << std::setw(14) << ylow_.call(xi,yi) + ybinlength_
	     << std::setw(14) << countrate_.call(xi,yi)
	     << std::endl;
      
      Data << std::endl;
      
      for(unsigned yi = 0; yi < Ny_; yi++)
	Data << std::setw(14) << xlow_.call(xi,yi) + xbinlength_
	     << std::setw(14) << ylow_.call(xi,yi)
	     << std::setw(14) << countrate_.call(xi,yi)
	     << std::endl
	     << std::setw(14) << xlow_.call(xi,yi) + xbinlength_
	     << std::setw(14) << ylow_.call(xi,yi) + ybinlength_
	     << std::setw(14) << countrate_.call(xi,yi)
	     << std::endl;
      
      Data << std::endl;
    }
  Data.close();
  
  std::stringstream gpfilename_stst;
  gpfilename_stst << prefix_of_ << filename
		  << "_" << Nx_ << "x" << Ny_ << "_seed_" << seed_ << ".gp";
  std::string gpfilename_st = gpfilename_stst.str();
  std::ofstream Gnuplot( gpfilename_st.c_str() );
  
  Gnuplot << "set terminal epslatex standalone color colortext 14" << std::endl
	  << "set output '" << filename << "_" << Nx_ << "x" << Ny_ << "_seed_" << seed_ << ".tex'"
	  << std::endl << std::endl
	  << "xmin = " << xmin_ << std::endl
	  << "xmax = " << xmax_ << std::endl
	  << "ymin = " << ymin_ << std::endl
	  << "ymax = " << ymax_ << std::endl
	  << std::endl
	  << "unset xtics #border mirror norotate offset 0,0.3" << std::endl
	  << "unset ytics #border mirror norotate offset 0.3,0" << std::endl
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
	  << "unset xlabel #'$x$' offset 0,0.3" << std::endl
	  << "unset ylabel #'$y$' offset 0.3,0" << std::endl
	  << "set cblabel '\\large $k$' offset 0,0" << std::endl
	  << std::endl
	  << "set xrange [xmin:xmax]" << std::endl
	  << "set yrange [ymin:ymax]" << std::endl
	  << std::endl
	  << "set pm3d map" << std::endl
	  << "# set palette defined (0 'black', 3 'dark-blue', "
	  << "4 'purple', 6 'dark-red', 9 'red', 12 'orange', 15 'yellow')" << std::endl
	  << "splot '" << filename << "_" << Nx_ << "x" << Ny_ << "_seed_" << seed_ << ".dat" << "' notitle" << std::endl
	  << "#" << std::endl
	  << std::endl;
  
  Gnuplot.close(); 
}
// --------------------------------------

// --------------------------------------
// Subtract Sources
void BinnedSkyMap::DefineSource(const double &k_src, const Gaussian &sources)
{
  // Define sources
  sources_ = sources;
  
  // Auxiliary parameters
  double Lambda_i = 0;
  double c1,c2,c3,c4, c5;
  
  for(unsigned int ni = 0; ni < Nx_*Ny_; ni++)
    {      
      rx_ = xlow_.call(ni);
      c1 = gauss_pdf(ylow_.call(ni), this);
      
      rx_ = xlow_.call(ni) + xbinlength_;
      c2 = gauss_pdf(ylow_.call(ni), this);
      
      rx_ = xlow_.call(ni);
      c3 = gauss_pdf(ylow_.call(ni) + ybinlength_, this);

      rx_ = xlow_.call(ni) + xbinlength_;
      c4 = gauss_pdf(ylow_.call(ni) + ybinlength_, this);

      rx_ = xlow_.call(ni) + xbinlength_/2;
      c5 = gauss_pdf(ylow_.call(ni) + ybinlength_/2, this);

      if(c1 > (1e-300) && c2 > (1e-300) && c3 > (1e-300) && c4 > (1e-300) && c5 > (1e-300))
	{	  
	  Lambda_i = Lambda(k_src, 
			    xlow_.call(ni), xlow_.call(ni) + xbinlength_, 
			    ylow_.call(ni), ylow_.call(ni) + ybinlength_);      
	  srcsub_intensity_.assign(ni, Lambda_i);
	}
    }
}

void BinnedSkyMap::SubtractSource(const double &lambda){
  SubtractRandomBinomial(countrate_, lambda, srcsub_intensity_, seedsub_);
}
// --------------------------------------

// --------------------------------------
// Global Analysis
double BinnedSkyMap::global_deviation_strength_A_pix_wbc(const double &lambda, const std::string &filename) const
{
  double C = 0, D = 0, maxD = 0, p = 0, mean = 0;
  unsigned A = 0, thresh = 0, M = 0;
  std::vector<unsigned> threshvec; std::vector<double> pvec, Dvec; threshvec.clear(); pvec.clear(); Dvec.clear();
  do{
    A = 0;
    thresh++; threshvec.push_back(thresh);
    p = Bernoulli_p(thresh, lambda); pvec.push_back(p);
    for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
      A += (countrate_.call(ni) >= thresh);
    M = Nx_*Ny_;
    mean = M*p;
    
    if(A > mean){
      C = BinomialCDF_complement(M, A, p);
      
      unsigned bound = M+1;
      for(unsigned m = 0; m <= mean; m++)
	if( BinomialPDF(M, m, p) <= BinomialPDF(M, A, p) )
	  bound = m;
	else
	  break;
      if(bound < M+1)
	C += BinomialCDF(M, bound, p);
    }
    else{      
      C = BinomialCDF(M, A, p);
      
      unsigned bound = M+1;
      for(unsigned m = A+1; m < M+1; m++)
	if( BinomialPDF(M, m, p) <= BinomialPDF(M, A, p) ){
	  bound = m;
	  break;}
      if(bound < M+1)
	C += BinomialCDF_complement(M, bound, p);
    }
    if(C < pow(10,-300)) // Accuracy limit
      C = pow(10,-300);
    
    D = -log10(C); Dvec.push_back(D);
    if(D > maxD)
      maxD = D;
  }while(A > 0);
  
  std::stringstream filename_stst;
  filename_stst << prefix_of_ << filename
		<< "_" << Nx_ << "x" << Ny_ << "_seed_" << seed_ << ".dat";
  std::string filename_st = filename_stst.str();
  std::ofstream Data( filename_st.c_str() );
  Data << "# Deviation strength for a given countrate of a binned sky map" << std::endl
       << "#" << std::setw(13) << "threshold" << std::setw(14) << "p" << std::setw(14) << "D(A_pixelized)" << std::endl;
  for(unsigned i = 0; i < pvec.size(); i++)
    Data << std::setw(14) << threshvec[i] << std::setw(14) << pvec[i] << std::setw(14) <<  Dvec[i] << std::endl;
  Data.close();
  
  std::stringstream gpfilename_stst;
  gpfilename_stst << prefix_of_ << filename
		  << "_" << Nx_ << "x" << Ny_ << "_seed_" << seed_ << ".gp";
  std::string gpfilename_st = gpfilename_stst.str();
  std::ofstream Gnuplot( gpfilename_st.c_str() );
  Gnuplot << "set terminal epslatex standalone color colortext 14" << std::endl
	  << "set output '" << filename << "_" << Nx_ << "x" << Ny_ << "_seed_" << seed_ << ".tex'"
	  << std::endl << std::endl
	  << "xmin = " << 0 << std::endl
	  << "#xmin = " << pvec.back() << std::endl
	  << "xmax = " << threshvec.back() << std::endl
	  << "#xmax = " << 1 << std::endl
	  << "ymin = " << 0 << std::endl
	  << "ymax = " << maxD*1.1 << std::endl
	  << std::endl
	  << "#set xtics border mirror norotate offset 0,0.3" << std::endl
	  << "#set ytics border mirror norotate offset 0.3,0" << std::endl
	  << std::endl
	  << "set mxtics 2" << std::endl
	  << "set mytics 2" << std::endl
	  << std::endl
	  << "set border lw 2" << std::endl
	  << "set pointsize 1.2" << std::endl
	  << std::endl
	  << "set xlabel 'Threshold $\\rho$' offset 0,0.7" << std::endl
	  << "#set xlabel '$p$' offset 0,0.7" << std::endl
	  << "set ylabel '$\\mathcal{D}(A)$' offset 0,0"
	  << std::endl
	  << "set xrange [xmin:xmax]" << std::endl
	  << "set yrange [ymin:ymax]" << std::endl
	  << "#set logscale x" << std::endl
	  << std::endl
	  << "plot '" << filename << "_" << Nx_ << "x" << Ny_ << "_seed_" << seed_ << ".dat'"
	  << " u 1:3 w lp lt 1 pt 13 lw 3 lc rgb 'navy-blue' notitle,\\" << std::endl
	  << "       6.2416 lt 5 lc rgb 'black' lw 3 notitle" << std::endl
	  << "#" << std::endl
	  << std::endl;
  Gnuplot.close(); 
  
  return maxD;
}
// --------------------------------------

// --------------------------------------
// Local Analysis
BinField<double> BinnedSkyMap::significance_map(const double &lambda, const unsigned int &slide,
						const bool &verbose, const bool &matrix_out, const bool &gnuplot_out,
						const std::string &filename) const
{
  // Define significance map
  unsigned Nx_mink = Nx_-slide+1;
  unsigned Ny_mink = Ny_-slide+1;
  BinField< double > SignificanceMap(Nx_mink,Ny_mink,0);
  
  // Significance
  for(unsigned xi = 0; xi < Nx_mink; xi++)
    for(unsigned yi = 0; yi < Ny_mink; yi++){
      unsigned end_x = xi+slide, end_y = yi+slide;
      unsigned count = 0;
      for(unsigned sx = xi; sx < end_x; sx++)
      	for(unsigned sy = yi; sy < end_y; sy++)
	  count += call_countrate(sx, sy);
      
      // sigma = sqrt(2*(count*log(count/(lambda*slide*slide))+(lambda*slide*slide)-count))
      // D(sigma) = -log10(erfc(sigma/sqrt(2)))
      SignificanceMap.assign(xi,yi,-log10(erfc(sqrt(2*(count*log(count/(lambda*slide*slide))+(lambda*slide*slide)-count))/sqrt(2))));
    }
  
  // Output of the significance map:
  // Matrix
  if(matrix_out){
    std::stringstream matfilename_stst;
    matfilename_stst << filename
		     << "_Matrix_" << (Nx_-slide+1) << "x" << (Ny_-slide+1) << "_seed_" << seed_
		     << "_lambda_" << lambda << "_slide_" << slide << ".dat";
    std::string matfilename = matfilename_stst.str();
    
    SignificanceMap.MatrixToFile(matfilename,prefix_of_);
  }
  
  if(gnuplot_out){
    // Data output
    std::stringstream filename_stst;
    filename_stst << prefix_of_ << filename
		  << "_" << (Nx_-slide+1) << "x" << (Ny_-slide+1) << "_seed_" << seed_
		  << "_lambda_" << lambda << "_slide_" << slide << ".dat";
    std::string filename_st = filename_stst.str();
    std::ofstream Data( filename_st.c_str() );
    
    Data << "# Significance map" << std::endl << "#" << std::endl
	 << "# Background intensity: lambda = " << lambda << ";" << std::endl
	 << "# Size of sliding window: slide = " << slide << ";" << std::endl
	 << "#" << std::endl;
    Data << "#" << std::setw(13) << "X" << std::setw(14) << "Y" << std::setw(14) << "sigma" << std::endl;
    for(unsigned xi = 0; xi < (Nx_-slide+1); xi++)
      {
	for(unsigned yi = 0; yi < (Ny_-slide+1); yi++)
	  Data << std::setw(14) << xlow_.call(xi,yi) + (slide-1)*xbinlength_/2.0
	       << std::setw(14) << ylow_.call(xi,yi) + (slide-1)*ybinlength_/2.0
	       << std::setw(14) << SignificanceMap.call(xi,yi)
	       << std::endl
	       << std::setw(14) << xlow_.call(xi,yi) + (slide-1)*xbinlength_/2.0
	       << std::setw(14) << ylow_.call(xi,yi) + ybinlength_ + (slide-1)*ybinlength_/2.0
	       << std::setw(14) << SignificanceMap.call(xi,yi)
	       << std::endl;
      
	Data << std::endl;
      
	for(unsigned yi = 0; yi < (Ny_-slide+1); yi++)
	  Data << std::setw(14) << xlow_.call(xi,yi) + xbinlength_ + (slide-1)*xbinlength_/2.0
	       << std::setw(14) << ylow_.call(xi,yi) + (slide-1)*ybinlength_/2.0
	       << std::setw(14) << SignificanceMap.call(xi,yi)
	       << std::endl
	       << std::setw(14) << xlow_.call(xi,yi) + xbinlength_ + (slide-1)*xbinlength_/2.0
	       << std::setw(14) << ylow_.call(xi,yi) + ybinlength_ + (slide-1)*ybinlength_/2.0
	       << std::setw(14) << SignificanceMap.call(xi,yi)
	       << std::endl;
      
	Data << std::endl;
      }
    Data.close();

    // Gnuplot
    std::stringstream gpfilename_stst;
    gpfilename_stst << prefix_of_ << filename
		    << "_" << (Nx_-slide+1) << "x" << (Ny_-slide+1) << "_seed_" << seed_
		    << "_slide_" << slide << ".gp";
      //<< "_lambda_" << lambda << "_slide_" << slide << ".gp";
    std::string gpfilename_st = gpfilename_stst.str();
    std::ofstream Gnuplot( gpfilename_st.c_str() );
  
    Gnuplot << "set terminal epslatex standalone color colortext 14" << std::endl
	    << "set output '" << filename
	    << "_" << (Nx_-slide+1) << "x" << (Ny_-slide+1) << "_seed_" << seed_
	    << "_slide_" << slide << ".tex'";
    //<< "_lambda_" << lambda << "_slide_" << slide << ".tex'";
    Gnuplot << std::endl << std::endl
	    << "xmin = " << xmin_ + (slide-1)*xbinlength_/2.0 << std::endl
	    << "xmax = " << xmax_ - (slide-1)*xbinlength_/2.0 << std::endl
	    << "ymin = " << ymin_ + (slide-1)*ybinlength_/2.0 << std::endl
	    << "ymax = " << ymax_ - (slide-1)*ybinlength_/2.0 << std::endl
	    << std::endl
	    << "#set title 'Size of sliding observation window: " << slide << "with wbc'" << std::endl
	    << "unset xtics #border mirror norotate offset 0,0.3" << std::endl
	    << "unset ytics #border mirror norotate offset 0.3,0" << std::endl
	    << "#set cbtics 10" << std::endl
	    << std::endl
	    << "set mxtics 2" << std::endl
	    << "set mytics 2" << std::endl
	    << std::endl;
    if(Ny_ == Nx_)
      Gnuplot << "set size ratio 1" << std::endl;
    else
      Gnuplot << "set size ratio " << (Ny_-slide+1)*1.0/(Nx_-slide+1) << std::endl;
    Gnuplot << "set border lw 2" << std::endl
	    << std::endl
	    << "unset xlabel #'$x$' offset 0,0.7" << std::endl
	    << "unset ylabel #'$y$' offset 0.6,0" << std::endl
	    << "set cblabel '\\large $\\mathcal{D}(\\sigma)$' offset 0,0" << std::endl
      //<< "set cblabel '\\large $-\\log_{10}\\left(\\mathrm{erfc}(\\sigma_{\\mathrm{Li\\& Ma}}/\\sqrt{2})\\right)$' offset 0,0" << std::endl
	    << std::endl
	    << "set xrange [xmin:xmax]" << std::endl
	    << "set yrange [ymin:ymax]" << std::endl
	    << "set cbrange [0:9]" << std::endl
	    << std::endl
	    << "set pm3d map" << std::endl
	    << "# set palette defined (0 'black', 3 'dark-blue', "
	    << "4 'purple', 6 'dark-red', 9 'red', 12 'orange', 15 'yellow')" << std::endl
	    << "splot '" << filename << "_" << (Nx_-slide+1) << "x" << (Ny_-slide+1) << "_seed_" << seed_
	    << "_lambda_" << lambda << "_slide_" << slide << ".dat";
    Gnuplot << "' notitle" << std::endl
	    << "#" << std::endl
	    << std::endl;
    
    Gnuplot.close();
  }
  
  if( verbose )
    std::cout << "-------------------" << std::endl;
  
  return SignificanceMap;
}


double BinnedSkyMap::max_of_significance_map(const double &lambda, const unsigned int &slide) const
{
  double max = 0;
  double t = 0;
  
  unsigned Nx_mink = Nx_-slide+1;
  unsigned Ny_mink = Ny_-slide+1;
  BinField< double > tmp_map = significance_map(lambda, slide, false, false, false);
  for(unsigned ni = 0; ni < Nx_mink*Ny_mink; ni++){
    t = tmp_map.call(ni);
    if(t > max)
      max = t;
  }
  
  return max;
}


BinField<double> BinnedSkyMap::significance_map_chessboard(const double &lambda, const unsigned int &slide,
							   const bool &verbose, const bool &matrix_out,
							   const std::string &filename) const
{
  // Define significance map
  unsigned Nx_mink = Nx_/slide;
  unsigned Ny_mink = Ny_/slide;
  BinField< double > SignificanceMap(Nx_mink,Ny_mink,0);
  
  // Significance
  for(unsigned xi = 0; xi < Nx_mink; xi++)
    for(unsigned yi = 0; yi < Ny_mink; yi++){
      unsigned end_x = xi*slide+slide, end_y = yi*slide+slide;
      unsigned count = 0;
      for(unsigned sx = xi*slide; sx < end_x; sx++)
      	for(unsigned sy = yi*slide; sy < end_y; sy++)
	  count += call_countrate(sx, sy);
      
      // sigma = sqrt(2*(count*log(count/(lambda*slide*slide))+(lambda*slide*slide)-count))
      // D(sigma) = -log10(erfc(sigma/sqrt(2)))
      SignificanceMap.assign(xi,yi,-log10(erfc(sqrt((count*log(count/(lambda*slide*slide))+(lambda*slide*slide)-count)))));
    }
  
  // Output of the significance map:
  // Matrix
  if(matrix_out){
    std::stringstream matfilename_stst;
    matfilename_stst << filename
		     << "_Matrix_" << Nx_mink << "x" << Ny_mink << "_seed_" << seed_
		     << "_lambda_" << lambda << "_slide_" << slide << ".dat";
    std::string matfilename = matfilename_stst.str();
    
    SignificanceMap.MatrixToFile(matfilename,prefix_of_);
  }
  
  if( verbose )
    std::cout << "-------------------" << std::endl;
  
  return SignificanceMap;
}

BinField<double> BinnedSkyMap::Minkowski_sky_map_APC_pix_wbc(const double &lambda, const unsigned int &slide,
							     const bool &verbose, const bool &matrix_out, const bool &gnuplot_out,
							     const std::string &filename) const
{  
  // Define MinkowskiSkyMap
  unsigned Nx_mink = Nx_-slide+1;
  unsigned Ny_mink = Ny_-slide+1;
  BinField< double > MinkowskiSkyMap(Nx_mink,Ny_mink,0);
  unsigned slide2 = slide*slide;
  unsigned slide_minus_one = slide - 1;
  unsigned N_A = slide*slide + 1;
  
  // DoS: ready to read in
  std::vector< stateclassifier > DoS_APC; DoS_APC.clear();
  
  // Read in DoS
  if( verbose ){
    std::cout << "Reading DoS ..." << std::endl
	      << PercentageBar << std::endl;}
  unsigned percent = 0;
  for(unsigned Ai = 0; Ai < N_A; Ai++){
    std::stringstream infile_st;
    infile_st << "parameters/structure_" << slide << "x" << slide << "/DoS_averaged_A_"
	      << Ai
	      << "_PC_pix_wbc_" << slide << "x" << slide << "_seed_0_seedwght_0.dat";
    std::string infile = infile_st.str();
    std::ifstream structure_distribution((infile).c_str ());
    if(structure_distribution.fail()){
      std::cerr << "ERROR: ifstream failed to read " << "parameters/structure_" << slide << "x" << slide << "/DoS_averaged_A_"
		<< Ai << "_PC_pix_wbc_" << slide << "x" << slide << "_seed_0_seedwght_0.dat;" << std::endl;
      exit(-1);
    }
    while(!structure_distribution.eof()){
      int P(0), C(0);
      structure_distribution >> P;
      structure_distribution >> C;
      
      double DoS = 0;
      structure_distribution >> DoS;
      if(DoS)
	DoS_APC.push_back(stateclassifier(Ai,P,C,DoS));
    }
    structure_distribution.close();
    if(verbose) 
      Cout_One_Percent(Ai,percent,N_A);
  }
  unsigned nmbr_of_macrostates = DoS_APC.size();
  
  // Determine maximum and minimum thresholds
  unsigned thresh_max = 0;
  for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
    if(countrate_.call(ni) > thresh_max)
      thresh_max = countrate_.call(ni);
  thresh_max++;
  unsigned thresh_min = thresh_max;
  for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
    if(countrate_.call(ni) < thresh_min)
      thresh_min = countrate_.call(ni);
  unsigned nmbr_of_thresh = thresh_max-thresh_min+1;
  
  // Compute Minkowski sky map
  if( verbose ){
    std::cout << "Compute Minkowski sky map ..." << std::endl
	      << PercentageBar << std::endl;
  }
  percent = 0;
  
  // Prepare excursion set and scan-window functional values
  BinField< bool > Excursion_set_old(Nx_,Ny_,true);
  BinField< bool > Excursion_set_older(Nx_,Ny_,true);
  BinField< bool > Excursion_set_new(Nx_,Ny_,true);
  BinField< int > A_Map(Nx_mink,Ny_mink,slide2);
  BinField< int > P_Map(Nx_mink,Ny_mink,4*slide);
  BinField< int > C_Map(Nx_mink,Ny_mink,1);
  
  double p(0), invp(0), Dev(0);
  std::vector< stateclassifier > prob_APC(nmbr_of_macrostates, stateclassifier(0,0,0,1));
  std::vector< stateclassifier > comp_APC(nmbr_of_macrostates, stateclassifier(0,0,0,1));
  
  for(unsigned thresh = thresh_min; thresh <= thresh_max; thresh++){\
    /*
    std::stringstream filename_stst;
    filename_stst << "excursion_set" << thresh << ".dat";
    std::string filename = filename_stst.str();
    Excursion_set_new.fout(filename,prefix_of_);
    */
    p = Bernoulli_p(thresh, lambda);
    invp = 1-p;
    
    // Excursion set
    for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
      Excursion_set_new.assign(ni, countrate_.call(ni) >= thresh);
    
    // determine probability distribution
    for(std::vector< stateclassifier >::size_type sr = 0; sr < nmbr_of_macrostates; sr++){
      int Ai = DoS_APC[sr].A();
      prob_APC[sr] = stateclassifier(Ai,DoS_APC[sr].P(),DoS_APC[sr].C(),DoS_APC[sr].cl()*pow(p,Ai)*pow(invp,slide2-Ai));
    }
    
    // determine compatibility distribution
    boost::sort(prob_APC);
    comp_APC=prob_APC;
    
    double min_comp = comp_APC[0].cl();
    for(std::vector< stateclassifier >::size_type sr = 1; sr < comp_APC.size(); sr++)
      comp_APC[sr].add_cl(comp_APC[sr-1].cl());
    
    // Deviation strength
    for(unsigned xi = 0; xi < Nx_mink; xi++)
      for(unsigned yi = 0; yi < Ny_mink; yi++){
	unsigned end_x = xi+slide_minus_one, end_y = yi+slide_minus_one;
	
	unsigned sx = xi;
	unsigned sy = yi;
	{
	  if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	    A_Map.add(xi,yi,-1);
	    if( !Excursion_set_old.call(sx,sy+1) )
	      P_Map.add(xi,yi,-2);
	    
	    if( !Excursion_set_old.call(sx+1,sy) )
	      P_Map.add(xi,yi,-2);
	    
	    C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx+1,sy+1)*pow_of_2[2]+ Excursion_set_old.call(sx+1,sy)*pow_of_2[4]]);
	    
	    Excursion_set_old.assign(sx,sy,false);
	  } 
	}
	for(sy = yi+1; sy < end_y; sy++)
	  {
	    if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	      A_Map.add(xi,yi,-1);
	      if( Excursion_set_old.call(sx,sy+1) && Excursion_set_old.call(sx,sy-1) )
		P_Map.add(xi,yi,2);    
	      else if( !Excursion_set_old.call(sx,sy+1) && !Excursion_set_old.call(sx,sy-1) )
		P_Map.add(xi,yi,-2);
  	      
	      if( !Excursion_set_old.call(sx+1,sy) )
		P_Map.add(xi,yi,-2);
  
	      C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx+1,sy+1)*pow_of_2[2] + Excursion_set_old.call(sx+1,sy)*pow_of_2[4] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6] + Excursion_set_old.call(sx+1,sy-1)*pow_of_2[7]]);
  
	      Excursion_set_old.assign(sx,sy,false);
	    } 
	  }
	sy = end_y;
	{
	  if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	    A_Map.add(xi,yi,-1);    
	    if( !Excursion_set_old.call(sx,sy-1) )
	      P_Map.add(xi,yi,-2);
	    
	    if( !Excursion_set_old.call(sx+1,sy) )
	      P_Map.add(xi,yi,-2);
  
	    C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx+1,sy)*pow_of_2[4] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6] + Excursion_set_old.call(sx+1,sy-1)*pow_of_2[7]]);
  
	    Excursion_set_old.assign(sx,sy,false);
	  } 
	}
	
	for(sx = xi+1; sx < end_x; sx++){
	  sy = yi;
	  {
	    if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	      A_Map.add(xi,yi,-1);
	      if( !Excursion_set_old.call(sx,sy+1) )
		P_Map.add(xi,yi,-2);
	      
	      if( Excursion_set_old.call(sx+1,sy) && Excursion_set_old.call(sx-1,sy) )
		P_Map.add(xi,yi,2);    
	      else if( !Excursion_set_old.call(sx+1,sy) && !Excursion_set_old.call(sx-1,sy) )
		P_Map.add(xi,yi,-2);
  
	      C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy+1)*pow_of_2[0] + Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx+1,sy+1)*pow_of_2[2] + Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx+1,sy)*pow_of_2[4]]);
  
	      Excursion_set_old.assign(sx,sy,false);
	    } 
	  }
	  for(sy = yi+1; sy < end_y; sy++)
	    {
	      if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
		A_Map.add(xi,yi,-1);
		if( Excursion_set_old.call(sx,sy+1) && Excursion_set_old.call(sx,sy-1) )
		  P_Map.add(xi,yi,2);    
		else if( !Excursion_set_old.call(sx,sy+1) && !Excursion_set_old.call(sx,sy-1) )
		  P_Map.add(xi,yi,-2);
  
		if( Excursion_set_old.call(sx+1,sy) && Excursion_set_old.call(sx-1,sy) )
		  P_Map.add(xi,yi,2);    
		else if( !Excursion_set_old.call(sx+1,sy) && !Excursion_set_old.call(sx-1,sy) )
		  P_Map.add(xi,yi,-2);
  
		C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy+1)*pow_of_2[0] + Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx+1,sy+1)*pow_of_2[2] + Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx+1,sy)*pow_of_2[4] + Excursion_set_old.call(sx-1,sy-1)*pow_of_2[5] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6] + Excursion_set_old.call(sx+1,sy-1)*pow_of_2[7]]);
  
		Excursion_set_old.assign(sx,sy,false);
	      } 
	    }
	  sy = end_y;
	  {
	    if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	      A_Map.add(xi,yi,-1);    
	      if( !Excursion_set_old.call(sx,sy-1) )
		P_Map.add(xi,yi,-2);
  
	      if( Excursion_set_old.call(sx+1,sy) && Excursion_set_old.call(sx-1,sy) )
		P_Map.add(xi,yi,2);    
	      else if( !Excursion_set_old.call(sx+1,sy) && !Excursion_set_old.call(sx-1,sy) )
		P_Map.add(xi,yi,-2);
  
	      C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx+1,sy)*pow_of_2[4] + Excursion_set_old.call(sx-1,sy-1)*pow_of_2[5] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6] + Excursion_set_old.call(sx+1,sy-1)*pow_of_2[7]]);
  
	      Excursion_set_old.assign(sx,sy,false);
	    } 
	  }
	}
	
	sx = end_x;
	sy = yi;
	{
	  if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	    A_Map.add(xi,yi,-1);
	    if( !Excursion_set_old.call(sx,sy+1) )
	      P_Map.add(xi,yi,-2);
	    
	    if( !Excursion_set_old.call(sx-1,sy) )
	      P_Map.add(xi,yi,-2);
  
	    C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy+1)*pow_of_2[0] + Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx-1,sy)*pow_of_2[3]]);
  
	    Excursion_set_old.assign(sx,sy,false);
	  } 
	}
	for(sy = yi+1; sy < end_y; sy++)
	  {
	    if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	      A_Map.add(xi,yi,-1);
	      if( Excursion_set_old.call(sx,sy+1) && Excursion_set_old.call(sx,sy-1) )
		P_Map.add(xi,yi,2);    
	      else if( !Excursion_set_old.call(sx,sy+1) && !Excursion_set_old.call(sx,sy-1) )
		P_Map.add(xi,yi,-2);
	      
	      if( !Excursion_set_old.call(sx-1,sy) )
		P_Map.add(xi,yi,-2);
  
	      C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy+1)*pow_of_2[0] + Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx-1,sy-1)*pow_of_2[5] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6]]);
  
	      Excursion_set_old.assign(sx,sy,false);
	    } 
	  }
	sy = end_y;
	{
	  if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	    A_Map.add(xi,yi,-1);    
	    if( !Excursion_set_old.call(sx,sy-1) )
	      P_Map.add(xi,yi,-2);
	    
	    if( !Excursion_set_old.call(sx-1,sy) )
	      P_Map.add(xi,yi,-2);
  
	    C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx-1,sy-1)*pow_of_2[5] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6]]);
  
	    Excursion_set_old.assign(sx,sy,false);
	  } 
	}
	// std::cout << A_Map.call(xi,yi) << "  " << P_Map.call(xi,yi) << "  " << C_Map.call(xi,yi) << "  " << std::endl;
	
	stateclassifier tmp_state(A_Map.call(xi,yi),P_Map.call(xi,yi),C_Map.call(xi,yi),0);
	int comp_APC_end = nmbr_of_macrostates-1;
	int tmp_sr = 0;
	for(int sr = comp_APC_end; sr >= 0; sr--)
	  if(tmp_state==comp_APC[sr]){
	    tmp_sr = sr;
	    break;
	  }
	
	Dev = comp_APC[tmp_sr].cl();
	if(Dev == 0){
	  Dev = min_comp;
	  std::cerr << "WARNING: unknown structure at (" << xi << "," << yi << "), A = " 
		    << A_Map.call(xi,yi) << ", P = " << P_Map.call(xi,yi) << ", and C = " << C_Map.call(xi,yi) << std::endl;
	}
	Dev = log10(Dev);
	
	if( A_Map.call(xi,yi) >= p*slide2 )
	  Dev *= -1;
	
       	if( fabs(MinkowskiSkyMap.call(xi,yi)) < fabs(Dev) )
	  MinkowskiSkyMap.assign(xi,yi,Dev);
	
	// std::cout << thresh << "    " << Dev << std::endl;
	
	// Excursion set
	for(sx = xi; sx <= end_x; sx++)
	  for(sy = yi; sy <= end_y; sy++)
	    Excursion_set_old.assign(sx, sy, Excursion_set_older.call(sx, sy));
      }// for(xi,yi)
    
    // Excursion set
    for(unsigned ni = 0; ni < Nx_*Ny_; ni++){
      Excursion_set_old.assign(ni, Excursion_set_new.call(ni));
      Excursion_set_older.assign(ni, Excursion_set_new.call(ni));
    }
    
    if(verbose)
      Cout_One_Percent(thresh-thresh_min,percent,nmbr_of_thresh);
  }
  
  // Trial correction
  for(unsigned xi = 0; xi < Nx_mink; xi++)
    for(unsigned yi = 0; yi < Ny_mink; yi++){
      double D = MinkowskiSkyMap.call(xi,yi);
      double Dt = fabs(D) - log10(trial); 
      if(fabs(D) < 9)
	Dt = log10(1-pow(1-pow(0.1, fabs(D) ),trial));
      MinkowskiSkyMap.assign(xi,yi,copysign(Dt,D));
    }
  
  double min_D = MinkowskiSkyMap.call(0,0);
  double max_D = min_D;
  for(unsigned xi = 0; xi < Nx_mink; xi++)
    for(unsigned yi = 0; yi < Ny_mink; yi++){
      double D = MinkowskiSkyMap.call(xi,yi);
      if ( D < min_D )
	min_D = D;
      else if ( D > max_D )
	max_D = D;
    }
  
  // Output of the Minkowski sky map:
  // Matrix
  if(matrix_out){
    std::stringstream matfilename_stst;
    matfilename_stst << filename
		     << "_Matrix_" << (Nx_-slide+1) << "x" << (Ny_-slide+1) << "_seed_" << seed_
		     << "_lambda_" << lambda << "_slide_" << slide << ".dat";
    std::string matfilename = matfilename_stst.str();
    
    MinkowskiSkyMap.MatrixToFile(matfilename,prefix_of_);
  }
  
  if(gnuplot_out){
    // Data output
    std::stringstream filename_stst;
    filename_stst << prefix_of_ << filename
		  << "_" << (Nx_-slide+1) << "x" << (Ny_-slide+1) << "_seed_" << seed_
		  << "_lambda_" << lambda << "_slide_" << slide << ".dat";
    std::string filename_st = filename_stst.str();
    std::ofstream Data( filename_st.c_str() );
    
    Data << "# Minkowski sky map" << std::endl << "#" << std::endl
	 << "# Background intensity: lambda = " << lambda << ";" << std::endl
	 << "# Functionals: APC;" << std::endl
	 << "# Size of sliding window: slide = " << slide << ";" << std::endl
	 << "# White boundary condition;" << std::endl
	 << "# Pixelized data" << std::endl << "#" << std::endl;
    Data << "#" << std::setw(13) << "X" << std::setw(14) << "Y" << std::setw(14) << "D(A, P, chi)" << std::endl;
    for(unsigned xi = 0; xi < (Nx_-slide+1); xi++)
      {
	for(unsigned yi = 0; yi < (Ny_-slide+1); yi++)
	  Data << std::setw(14) << xlow_.call(xi,yi) + (slide-1)*xbinlength_/2.0
	       << std::setw(14) << ylow_.call(xi,yi) + (slide-1)*ybinlength_/2.0
	       << std::setw(14) << MinkowskiSkyMap.call(xi,yi)
	       << std::endl
	       << std::setw(14) << xlow_.call(xi,yi) + (slide-1)*xbinlength_/2.0
	       << std::setw(14) << ylow_.call(xi,yi) + ybinlength_ + (slide-1)*ybinlength_/2.0
	       << std::setw(14) << MinkowskiSkyMap.call(xi,yi)
	       << std::endl;
      
	Data << std::endl;
      
	for(unsigned yi = 0; yi < (Ny_-slide+1); yi++)
	  Data << std::setw(14) << xlow_.call(xi,yi) + xbinlength_ + (slide-1)*xbinlength_/2.0
	       << std::setw(14) << ylow_.call(xi,yi) + (slide-1)*ybinlength_/2.0
	       << std::setw(14) << MinkowskiSkyMap.call(xi,yi)
	       << std::endl
	       << std::setw(14) << xlow_.call(xi,yi) + xbinlength_ + (slide-1)*xbinlength_/2.0
	       << std::setw(14) << ylow_.call(xi,yi) + ybinlength_ + (slide-1)*ybinlength_/2.0
	       << std::setw(14) << MinkowskiSkyMap.call(xi,yi)
	       << std::endl;
      
	Data << std::endl;
      }
    Data.close();

    // Gnuplot
    std::stringstream gpfilename_stst;
    gpfilename_stst << prefix_of_ << filename
		    << "_" << (Nx_-slide+1) << "x" << (Ny_-slide+1) << "_seed_" << seed_
		    << "_slide_" << slide << ".gp";
      //<< "_lambda_" << lambda << "_slide_" << slide << ".gp";
    std::string gpfilename_st = gpfilename_stst.str();
    std::ofstream Gnuplot( gpfilename_st.c_str() );
  
    Gnuplot << "set terminal epslatex standalone color colortext 14" << std::endl
	    << "set output '" << filename
	    << "_" << (Nx_-slide+1) << "x" << (Ny_-slide+1) << "_seed_" << seed_
	    << "_slide_" << slide << ".tex'";
    //<< "_lambda_" << lambda << "_slide_" << slide << ".tex'";
    Gnuplot << std::endl << std::endl
	    << "xmin = " << xmin_ + (slide-1)*xbinlength_/2.0 << std::endl
	    << "xmax = " << xmax_ - (slide-1)*xbinlength_/2.0 << std::endl
	    << "ymin = " << ymin_ + (slide-1)*ybinlength_/2.0 << std::endl
	    << "ymax = " << ymax_ - (slide-1)*ybinlength_/2.0 << std::endl
	    << std::endl
	    << "#set title 'Size of sliding observation window: " << slide << "with wbc'" << std::endl
	    << "unset xtics #border mirror norotate offset 0,0.3" << std::endl
	    << "unset ytics #border mirror norotate offset 0.3,0" << std::endl
	    << "#set cbtics 10" << std::endl
	    << std::endl
	    << "set mxtics 2" << std::endl
	    << "set mytics 2" << std::endl
	    << std::endl;
    if(Ny_ == Nx_)
      Gnuplot << "set size ratio 1" << std::endl;
    else
      Gnuplot << "set size ratio " << (Ny_-slide+1)*1.0/(Nx_-slide+1) << std::endl;
    Gnuplot << "set border lw 2" << std::endl
	    << std::endl
	    << "unset xlabel #'$x$' offset 0,0.7" << std::endl
	    << "unset ylabel #'$y$' offset 0.6,0" << std::endl
	    << "# set cblabel '\\large $\\mathcal{D}(A,P,\\chi)$' offset 0,0" << std::endl
	    << std::endl
	    << "set xrange [xmin:xmax]" << std::endl
	    << "set yrange [ymin:ymax]" << std::endl
	    << "set cbrange [0:9]" << std::endl
	    << "#set cbrange [" << min_D << ":" << max_D << "]" << std::endl
	    << std::endl
	    << "set pm3d map" << std::endl
	    << "# set palette defined (" << min_D << " 'black', " << 0.5*min_D << " 'dark-blue', " << 0 << " 'white', " << 0.2*max_D << " 'purple', " << 0.4*max_D << " 'dark-red', " << 0.6*max_D << " 'red', " << 0.8*max_D << " 'orange', " << max_D << " 'yellow')" << std::endl
	    << "splot '" << filename << "_" << (Nx_-slide+1) << "x" << (Ny_-slide+1) << "_seed_" << seed_
	    << "_lambda_" << lambda << "_slide_" << slide << ".dat";
    Gnuplot << "' notitle" << std::endl
	    << "#" << std::endl
	    << std::endl;
    
    Gnuplot.close();
  }
  
  if( verbose )
    std::cout << "-------------------" << std::endl;
  
  return MinkowskiSkyMap;
}

BinField<double> BinnedSkyMap::Minkowski_sky_map_APC_locally_different_excursion_sets_pix_wbc(const double &lambda, const unsigned int &slide,
											      const std::vector< BinField< unsigned > > &countrate_locally_corrected,
											      const bool &verbose, const bool &matrix_out, const bool &gnuplot_out,
											      const std::string &filename) const
{
  // Define MinkowskiSkyMap
  unsigned Nx_mink = Nx_-slide+1;
  unsigned Ny_mink = Ny_-slide+1;
  unsigned N_mink_prdct = Nx_mink*Ny_mink;
  BinField< double > MinkowskiSkyMap(Nx_mink,Ny_mink,0);
  unsigned slide2 = slide*slide;
  unsigned slide_minus_one = slide - 1;
  unsigned end_x = slide_minus_one, end_y = slide_minus_one;
  unsigned N_A = slide*slide + 1;
  
  // DoS: ready to read in
  std::vector< stateclassifier > DoS_APC; DoS_APC.clear();
  
  // Read in DoS
  if( verbose ){
    std::cout << "Reading DoS ..." << std::endl
	      << PercentageBar << std::endl;}
  unsigned percent = 0;
  for(unsigned Ai = 0; Ai < N_A; Ai++){
    std::stringstream infile_st;
    infile_st << "parameters/structure_" << slide << "x" << slide << "/DoS_averaged_A_"
	      << Ai
	      << "_PC_pix_wbc_" << slide << "x" << slide << "_seed_0_seedwght_0.dat";
    std::string infile = infile_st.str();
    std::ifstream structure_distribution((infile).c_str ());
    if(structure_distribution.fail()){
      std::cerr << "ERROR: ifstream failed to read " << "parameters/structure_" << slide << "x" << slide << "/DoS_averaged_A_"
		<< Ai << "_PC_pix_wbc_" << slide << "x" << slide << "_seed_0_seedwght_0.dat;" << std::endl;
      exit(-1);
    }
    while(!structure_distribution.eof()){
      int P(0), C(0);
      structure_distribution >> P;
      structure_distribution >> C;
      
      double DoS = 0;
      structure_distribution >> DoS;
      if(DoS)
	DoS_APC.push_back(stateclassifier(Ai,P,C,DoS));
    }
    structure_distribution.close();
    if(verbose) 
      Cout_One_Percent(Ai,percent,N_A);
  }
  unsigned nmbr_of_macrostates = DoS_APC.size();
  
  // Compute Minkowski sky map
  if( verbose ){
    std::cout << "Compute Minkowski sky map ..." << std::endl
	      << PercentageBar << std::endl;
  }
  percent = 0;
  
  unsigned center = slide/2;
  for(unsigned xi = 0; xi < Nx_mink; xi++)
    for(unsigned yi = 0; yi < Ny_mink; yi++){
      unsigned Mi = xi*Ny_mink + yi;
      double target_f = exp( (- pow(xi+center - 0.5*Nx_,2)/pow(Nx_*0.3,2) - pow(yi+center - 0.5*Ny_,2)/pow(Ny_*0.3,2) )*0.5 );
      double lambda_i = target_f*lambda;
      
      // Determine maximum and minimum thresholds
      unsigned thresh_max = 0;
      for(unsigned ni = 0; ni < slide2; ni++)
	if(countrate_locally_corrected[Mi].call(ni) > thresh_max)
	  thresh_max = countrate_locally_corrected[Mi].call(ni);
      thresh_max++;
      unsigned thresh_min = thresh_max;
      for(unsigned ni = 0; ni < slide2; ni++)
	if(countrate_locally_corrected[Mi].call(ni) < thresh_min)
	  thresh_min = countrate_locally_corrected[Mi].call(ni);
        
      // Prepare excursion set and scan-window functional values
      BinField< bool > Excursion_set_old(slide,slide,true);
      BinField< bool > Excursion_set_new(slide,slide,true);
    
      int A_Map(slide2);
      int P_Map(4*slide);
      int C_Map(1);
  
      double p(0), invp(0), Dev(0);
      std::vector< stateclassifier > prob_APC(nmbr_of_macrostates, stateclassifier(0,0,0,1));
      std::vector< stateclassifier > comp_APC(nmbr_of_macrostates, stateclassifier(0,0,0,1));
  
      for(unsigned thresh = thresh_min; thresh <= thresh_max; thresh++){
	p = Bernoulli_p(thresh, lambda_i);
	invp = 1-p;
	
	// Excursion set
	for(unsigned ni = 0; ni < slide2; ni++)
	  Excursion_set_new.assign(ni, countrate_locally_corrected[Mi].call(ni) >= thresh);
      
	// determine probability distribution
	for(std::vector< stateclassifier >::size_type sr = 0; sr < nmbr_of_macrostates; sr++){
	  int Ai = DoS_APC[sr].A();
	  prob_APC[sr] = stateclassifier(Ai,DoS_APC[sr].P(),DoS_APC[sr].C(),DoS_APC[sr].cl()*pow(p,Ai)*pow(invp,slide2-Ai));
	}
      
	// determine compatibility distribution
	boost::sort(prob_APC);
	comp_APC=prob_APC;
      
	double min_comp = comp_APC[0].cl();
	for(std::vector< stateclassifier >::size_type sr = 1; sr < comp_APC.size(); sr++)
	  comp_APC[sr].add_cl(comp_APC[sr-1].cl());
	
	unsigned sx = 0;
	unsigned sy = 0;
	{
	  if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	    A_Map -= 1;
	    if( !Excursion_set_old.call(sx,sy+1) )
	      P_Map -= 2;
	    
	    if( !Excursion_set_old.call(sx+1,sy) )
	      P_Map -= 2;
	    
	    C_Map += change_in_C_b2w[Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx+1,sy+1)*pow_of_2[2]+ Excursion_set_old.call(sx+1,sy)*pow_of_2[4]];
	    
	    Excursion_set_old.assign(sx,sy,false);
	  } 
	}
	for(sy = 1; sy < end_y; sy++)
	  {
	    if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	      A_Map -= 1;
	      if( Excursion_set_old.call(sx,sy+1) && Excursion_set_old.call(sx,sy-1) )
		P_Map += 2;    
	      else if( !Excursion_set_old.call(sx,sy+1) && !Excursion_set_old.call(sx,sy-1) )
		P_Map -= 2;
  	      
	      if( !Excursion_set_old.call(sx+1,sy) )
		P_Map -= 2;
  
	      C_Map += change_in_C_b2w[Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx+1,sy+1)*pow_of_2[2] + Excursion_set_old.call(sx+1,sy)*pow_of_2[4] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6] + Excursion_set_old.call(sx+1,sy-1)*pow_of_2[7]];
  
	      Excursion_set_old.assign(sx,sy,false);
	    } 
	  }
	sy = end_y;
	{
	  if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	    A_Map -= 1;    
	    if( !Excursion_set_old.call(sx,sy-1) )
	      P_Map -= 2;
	    
	    if( !Excursion_set_old.call(sx+1,sy) )
	      P_Map -= 2;
  
	    C_Map += change_in_C_b2w[Excursion_set_old.call(sx+1,sy)*pow_of_2[4] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6] + Excursion_set_old.call(sx+1,sy-1)*pow_of_2[7]];
  
	    Excursion_set_old.assign(sx,sy,false);
	  } 
	}
	
	for(sx = 1; sx < end_x; sx++){
	  sy = 0;
	  {
	    if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	      A_Map -= 1;
	      if( !Excursion_set_old.call(sx,sy+1) )
		P_Map -= 2;
	      
	      if( Excursion_set_old.call(sx+1,sy) && Excursion_set_old.call(sx-1,sy) )
		P_Map += 2;    
	      else if( !Excursion_set_old.call(sx+1,sy) && !Excursion_set_old.call(sx-1,sy) )
		P_Map -= 2;
  
	      C_Map += change_in_C_b2w[Excursion_set_old.call(sx-1,sy+1)*pow_of_2[0] + Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx+1,sy+1)*pow_of_2[2] + Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx+1,sy)*pow_of_2[4]];
  
	      Excursion_set_old.assign(sx,sy,false);
	    } 
	  }
	  for(sy = 1; sy < end_y; sy++)
	    {
	      if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
		A_Map -= 1;
		if( Excursion_set_old.call(sx,sy+1) && Excursion_set_old.call(sx,sy-1) )
		  P_Map += 2;    
		else if( !Excursion_set_old.call(sx,sy+1) && !Excursion_set_old.call(sx,sy-1) )
		  P_Map -= 2;
  
		if( Excursion_set_old.call(sx+1,sy) && Excursion_set_old.call(sx-1,sy) )
		  P_Map += 2;    
		else if( !Excursion_set_old.call(sx+1,sy) && !Excursion_set_old.call(sx-1,sy) )
		  P_Map -= 2;
  
		C_Map += change_in_C_b2w[Excursion_set_old.call(sx-1,sy+1)*pow_of_2[0] + Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx+1,sy+1)*pow_of_2[2] + Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx+1,sy)*pow_of_2[4] + Excursion_set_old.call(sx-1,sy-1)*pow_of_2[5] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6] + Excursion_set_old.call(sx+1,sy-1)*pow_of_2[7]];
  
		Excursion_set_old.assign(sx,sy,false);
	      } 
	    }
	  sy = end_y;
	  {
	    if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	      A_Map -= 1;    
	      if( !Excursion_set_old.call(sx,sy-1) )
		P_Map -= 2;
  
	      if( Excursion_set_old.call(sx+1,sy) && Excursion_set_old.call(sx-1,sy) )
		P_Map += 2;    
	      else if( !Excursion_set_old.call(sx+1,sy) && !Excursion_set_old.call(sx-1,sy) )
		P_Map -= 2;
  
	      C_Map += change_in_C_b2w[Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx+1,sy)*pow_of_2[4] + Excursion_set_old.call(sx-1,sy-1)*pow_of_2[5] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6] + Excursion_set_old.call(sx+1,sy-1)*pow_of_2[7]];
  
	      Excursion_set_old.assign(sx,sy,false);
	    } 
	  }
	}
	
	sx = end_x;
	sy = 0;
	{
	  if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	    A_Map -= 1;
	    if( !Excursion_set_old.call(sx,sy+1) )
	      P_Map -= 2;
	    
	    if( !Excursion_set_old.call(sx-1,sy) )
	      P_Map -= 2;
  
	    C_Map += change_in_C_b2w[Excursion_set_old.call(sx-1,sy+1)*pow_of_2[0] + Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx-1,sy)*pow_of_2[3]];
  
	    Excursion_set_old.assign(sx,sy,false);
	  } 
	}
	for(sy = 1; sy < end_y; sy++)
	  {
	    if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	      A_Map -= 1;
	      if( Excursion_set_old.call(sx,sy+1) && Excursion_set_old.call(sx,sy-1) )
		P_Map += 2;    
	      else if( !Excursion_set_old.call(sx,sy+1) && !Excursion_set_old.call(sx,sy-1) )
		P_Map -= 2;
	      
	      if( !Excursion_set_old.call(sx-1,sy) )
		P_Map -= 2;
  
	      C_Map += change_in_C_b2w[Excursion_set_old.call(sx-1,sy+1)*pow_of_2[0] + Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx-1,sy-1)*pow_of_2[5] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6]];
  
	      Excursion_set_old.assign(sx,sy,false);
	    } 
	  }
	sy = end_y;
	{
	  if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	    A_Map -= 1;    
	    if( !Excursion_set_old.call(sx,sy-1) )
	      P_Map -= 2;
	    
	    if( !Excursion_set_old.call(sx-1,sy) )
	      P_Map -= 2;
  
	    C_Map += change_in_C_b2w[Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx-1,sy-1)*pow_of_2[5] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6]];
  
	    Excursion_set_old.assign(sx,sy,false);
	  } 
	}
	// std::cout << A_Map.call(xi,yi) << "  " << P_Map.call(xi,yi) << "  " << C_Map.call(xi,yi) << "  " << std::endl;
	
	stateclassifier tmp_state(A_Map,P_Map,C_Map,0);
	int comp_APC_end = nmbr_of_macrostates-1;
	int tmp_sr = 0;
	for(int sr = comp_APC_end; sr >= 0; sr--)
	  if(tmp_state==comp_APC[sr]){
	    tmp_sr = sr;
	    break;
	  }
	
	Dev = comp_APC[tmp_sr].cl();
	if(Dev == 0){
	  Dev = min_comp;
	  std::cerr << "WARNING: unknown structure at (" << xi << "," << yi << "), A = " 
		    << A_Map << ", P = " << P_Map << ", and C = " << C_Map << std::endl;
	}
	Dev = log10(Dev);
      
	if( A_Map >= p*slide2 )
	  Dev *= -1;
      
	if( fabs(MinkowskiSkyMap.call(xi,yi)) < fabs(Dev) )
	  MinkowskiSkyMap.assign(xi,yi,Dev);
      
	// std::cout << thresh << "    " << Dev << std::endl;
	
	// Excursion set
	for(unsigned ni = 0; ni < slide2; ni++)
	  Excursion_set_old.assign(ni, Excursion_set_new.call(ni));
	
      }// for(thresh)
      
      if(verbose)
	Cout_One_Percent(Mi,percent,N_mink_prdct);
    }// for(xi, yi)
  
  // Trial correction
  for(unsigned xi = 0; xi < Nx_mink; xi++)
    for(unsigned yi = 0; yi < Ny_mink; yi++){
      double D = MinkowskiSkyMap.call(xi,yi);
      double Dt = fabs(D) - log10(trial); 
      if(fabs(D) < 9)
	Dt = log10(1-pow(1-pow(0.1, fabs(D) ),trial));
      MinkowskiSkyMap.assign(xi,yi,copysign(Dt,D));
    }
  
  double min_D = MinkowskiSkyMap.call(0,0);
  double max_D = min_D;
  for(unsigned xi = 0; xi < Nx_mink; xi++)
    for(unsigned yi = 0; yi < Ny_mink; yi++){
      double D = MinkowskiSkyMap.call(xi,yi);
      if ( D < min_D )
	min_D = D;
      else if ( D > max_D )
	max_D = D;
    }
  
  // Output of the Minkowski sky map:
  // Matrix
  if(matrix_out){
    std::stringstream matfilename_stst;
    matfilename_stst << filename
		     << "_Matrix_" << (Nx_-slide+1) << "x" << (Ny_-slide+1) << "_seed_" << seed_
		     << "_lambda_" << lambda << "_slide_" << slide << ".dat";
    std::string matfilename = matfilename_stst.str();
    
    MinkowskiSkyMap.MatrixToFile(matfilename,prefix_of_);
  }
  
  if(gnuplot_out){
    // Data output
    std::stringstream filename_stst;
    filename_stst << prefix_of_ << filename
		  << "_" << (Nx_-slide+1) << "x" << (Ny_-slide+1) << "_seed_" << seed_
		  << "_lambda_" << lambda << "_slide_" << slide << ".dat";
    std::string filename_st = filename_stst.str();
    std::ofstream Data( filename_st.c_str() );
    
    Data << "# Minkowski sky map" << std::endl << "#" << std::endl
	 << "# Background intensity: lambda = " << lambda << ";" << std::endl
	 << "# Functionals: APC;" << std::endl
	 << "# Size of sliding window: slide = " << slide << ";" << std::endl
	 << "# White boundary condition;" << std::endl
	 << "# Pixelized data" << std::endl << "#" << std::endl;
    Data << "#" << std::setw(13) << "X" << std::setw(14) << "Y" << std::setw(14) << "D(A, P, chi)" << std::endl;
    for(unsigned xi = 0; xi < (Nx_-slide+1); xi++)
      {
	for(unsigned yi = 0; yi < (Ny_-slide+1); yi++)
	  Data << std::setw(14) << xlow_.call(xi,yi) + (slide-1)*xbinlength_/2.0
	       << std::setw(14) << ylow_.call(xi,yi) + (slide-1)*ybinlength_/2.0
	       << std::setw(14) << MinkowskiSkyMap.call(xi,yi)
	       << std::endl
	       << std::setw(14) << xlow_.call(xi,yi) + (slide-1)*xbinlength_/2.0
	       << std::setw(14) << ylow_.call(xi,yi) + ybinlength_ + (slide-1)*ybinlength_/2.0
	       << std::setw(14) << MinkowskiSkyMap.call(xi,yi)
	       << std::endl;
      
	Data << std::endl;
      
	for(unsigned yi = 0; yi < (Ny_-slide+1); yi++)
	  Data << std::setw(14) << xlow_.call(xi,yi) + xbinlength_ + (slide-1)*xbinlength_/2.0
	       << std::setw(14) << ylow_.call(xi,yi) + (slide-1)*ybinlength_/2.0
	       << std::setw(14) << MinkowskiSkyMap.call(xi,yi)
	       << std::endl
	       << std::setw(14) << xlow_.call(xi,yi) + xbinlength_ + (slide-1)*xbinlength_/2.0
	       << std::setw(14) << ylow_.call(xi,yi) + ybinlength_ + (slide-1)*ybinlength_/2.0
	       << std::setw(14) << MinkowskiSkyMap.call(xi,yi)
	       << std::endl;
      
	Data << std::endl;
      }
    Data.close();

    // Gnuplot
    std::stringstream gpfilename_stst;
    gpfilename_stst << prefix_of_ << filename
		    << "_" << (Nx_-slide+1) << "x" << (Ny_-slide+1) << "_seed_" << seed_
		    << "_slide_" << slide << ".gp";
      //<< "_lambda_" << lambda << "_slide_" << slide << ".gp";
    std::string gpfilename_st = gpfilename_stst.str();
    std::ofstream Gnuplot( gpfilename_st.c_str() );
  
    Gnuplot << "set terminal epslatex standalone color colortext 14" << std::endl
	    << "set output '" << filename
	    << "_" << (Nx_-slide+1) << "x" << (Ny_-slide+1) << "_seed_" << seed_
	    << "_slide_" << slide << ".tex'";
    //<< "_lambda_" << lambda << "_slide_" << slide << ".tex'";
    Gnuplot << std::endl << std::endl
	    << "xmin = " << xmin_ + (slide-1)*xbinlength_/2.0 << std::endl
	    << "xmax = " << xmax_ - (slide-1)*xbinlength_/2.0 << std::endl
	    << "ymin = " << ymin_ + (slide-1)*ybinlength_/2.0 << std::endl
	    << "ymax = " << ymax_ - (slide-1)*ybinlength_/2.0 << std::endl
	    << std::endl
	    << "#set title 'Size of sliding observation window: " << slide << "with wbc'" << std::endl
	    << "unset xtics #border mirror norotate offset 0,0.3" << std::endl
	    << "unset ytics #border mirror norotate offset 0.3,0" << std::endl
	    << "#set cbtics 10" << std::endl
	    << std::endl
	    << "set mxtics 2" << std::endl
	    << "set mytics 2" << std::endl
	    << std::endl;
    if(Ny_ == Nx_)
      Gnuplot << "set size ratio 1" << std::endl;
    else
      Gnuplot << "set size ratio " << (Ny_-slide+1)*1.0/(Nx_-slide+1) << std::endl;
    Gnuplot << "set border lw 2" << std::endl
	    << std::endl
	    << "unset xlabel #'$x$' offset 0,0.7" << std::endl
	    << "unset ylabel #'$y$' offset 0.6,0" << std::endl
	    << "# set cblabel '\\large $\\mathcal{D}(A,P,\\chi)$' offset 0,0" << std::endl
	    << std::endl
	    << "set xrange [xmin:xmax]" << std::endl
	    << "set yrange [ymin:ymax]" << std::endl
	    << "set cbrange [0:9]" << std::endl
	    << "#set cbrange [" << min_D << ":" << max_D << "]" << std::endl
	    << std::endl
	    << "set pm3d map" << std::endl
	    << "# set palette defined (" << min_D << " 'black', " << 0.5*min_D << " 'dark-blue', " << 0 << " 'white', " << 0.2*max_D << " 'purple', " << 0.4*max_D << " 'dark-red', " << 0.6*max_D << " 'red', " << 0.8*max_D << " 'orange', " << max_D << " 'yellow')" << std::endl
	    << "splot '" << filename << "_" << (Nx_-slide+1) << "x" << (Ny_-slide+1) << "_seed_" << seed_
	    << "_lambda_" << lambda << "_slide_" << slide << ".dat";
    Gnuplot << "' notitle" << std::endl
	    << "#" << std::endl
	    << std::endl;
    
    Gnuplot.close();
  }
  
  if( verbose )
    std::cout << "-------------------" << std::endl;
  
  return MinkowskiSkyMap;
}

BinField<double> BinnedSkyMap::Minkowski_sky_map_APC_pix_wbc_chessboard(const double &lambda, const unsigned int &slide,
									const bool &verbose, const bool &matrix_out,
									const std::string &filename) const
{
  // Define MinkowskiSkyMap
  unsigned Nx_mink = Nx_/slide;
  unsigned Ny_mink = Ny_/slide;
  BinField< double > MinkowskiSkyMap(Nx_mink,Ny_mink,0);
  //BinField< double > MinkowskiSkyMap_sum(Nx_mink,Ny_mink,0);
  //BinField< double > MinkowskiSkyMap_single(Nx_mink,Ny_mink,0);
  unsigned slide2 = slide*slide;
  unsigned slide_minus_one = slide - 1;
  unsigned N_A = slide*slide + 1;
  
  // DoS: ready to read in
  std::vector< stateclassifier > DoS_APC; DoS_APC.clear();
  
  // Read in DoS
  if( verbose ){
    std::cout << "Reading DoS ..." << std::endl
	      << PercentageBar << std::endl;}
  unsigned percent = 0;
  for(unsigned Ai = 0; Ai < N_A; Ai++){
    std::stringstream infile_st;
    infile_st << "parameters/structure_" << slide << "x" << slide << "/DoS_averaged_A_"
	      << Ai
	      << "_PC_pix_wbc_" << slide << "x" << slide << "_seed_0_seedwght_0.dat";
    std::string infile = infile_st.str();
    std::ifstream structure_distribution((infile).c_str ());
    if(structure_distribution.fail()){
      std::cerr << "ERROR: ifstream failed to read " << "parameters/structure_" << slide << "x" << slide << "/DoS_averaged_A_"
		<< Ai << "_PC_pix_wbc_" << slide << "x" << slide << "_seed_0_seedwght_0.dat;" << std::endl;
      exit(-1);
    }
    while(!structure_distribution.eof()){
      int P(0), C(0);
      structure_distribution >> P;
      structure_distribution >> C;
      
      double DoS = 0;
      structure_distribution >> DoS;
      if(DoS)
	DoS_APC.push_back(stateclassifier(Ai,P,C,DoS));
    }
    structure_distribution.close();
    if(verbose) 
      Cout_One_Percent(Ai,percent,N_A);
  }
  unsigned nmbr_of_macrostates = DoS_APC.size();
  
  // Determine maximum and minimum thresholds
  unsigned thresh_max = 0;
  for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
    if(countrate_.call(ni) > thresh_max)
      thresh_max = countrate_.call(ni);
  thresh_max++;
  unsigned thresh_min = thresh_max;
  for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
    if(countrate_.call(ni) < thresh_min)
      thresh_min = countrate_.call(ni);
  unsigned nmbr_of_thresh = thresh_max-thresh_min+1;
  
  // Compute Minkowski sky map
  if( verbose ){
    std::cout << "Compute Minkowski sky map ..." << std::endl
	      << PercentageBar << std::endl;
  }
  percent = 0;
  
  // Prepare excursion set and scan-window functional values
  BinField< bool > Excursion_set_old(Nx_,Ny_,true);
  BinField< bool > Excursion_set_new(Nx_,Ny_,true);
  BinField< int > A_Map(Nx_mink,Ny_mink,slide2);
  BinField< int > P_Map(Nx_mink,Ny_mink,4*slide);
  BinField< int > C_Map(Nx_mink,Ny_mink,1);
  
  double p(0), invp(0), Dev(0);
  std::vector< stateclassifier > prob_APC(nmbr_of_macrostates, stateclassifier(0,0,0,1));
  std::vector< stateclassifier > comp_APC(nmbr_of_macrostates, stateclassifier(0,0,0,1));
  
  for(unsigned thresh = thresh_min; thresh <= thresh_max; thresh++){
    p = Bernoulli_p(thresh, lambda);
    invp = 1-p;
    
    // Excursion set
    for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
      Excursion_set_new.assign(ni, countrate_.call(ni) >= thresh);
    
    // determine probability distribution
    for(std::vector< stateclassifier >::size_type sr = 0; sr < nmbr_of_macrostates; sr++){
      int Ai = DoS_APC[sr].A();
      prob_APC[sr] = stateclassifier(Ai,DoS_APC[sr].P(),DoS_APC[sr].C(),DoS_APC[sr].cl()*pow(p,Ai)*pow(invp,slide2-Ai));
    }
    
    // determine compatibility distribution
    boost::sort(prob_APC);
    comp_APC=prob_APC;
    
    double min_comp = comp_APC[0].cl();
    for(std::vector< stateclassifier >::size_type sr = 1; sr < comp_APC.size(); sr++)
      comp_APC[sr].add_cl(comp_APC[sr-1].cl());
    
    // Deviation strength
    for(unsigned xi = 0; xi < Nx_mink; xi++)
      for(unsigned yi = 0; yi < Ny_mink; yi++){
	unsigned end_x = xi*slide+slide_minus_one, end_y = yi*slide+slide_minus_one;
	
	unsigned sx = xi*slide;
	unsigned sy = yi*slide;
	{
	  if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	    A_Map.add(xi,yi,-1);
	    if( !Excursion_set_old.call(sx,sy+1) )
	      P_Map.add(xi,yi,-2);
	    
	    if( !Excursion_set_old.call(sx+1,sy) )
	      P_Map.add(xi,yi,-2);
	    
	    C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx+1,sy+1)*pow_of_2[2]+ Excursion_set_old.call(sx+1,sy)*pow_of_2[4]]);
	    
	    Excursion_set_old.assign(sx,sy,false);
	  } 
	}
	for(sy = yi*slide+1; sy < end_y; sy++)
	  {
	    if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	      A_Map.add(xi,yi,-1);
	      if( Excursion_set_old.call(sx,sy+1) && Excursion_set_old.call(sx,sy-1) )
		P_Map.add(xi,yi,2);    
	      else if( !Excursion_set_old.call(sx,sy+1) && !Excursion_set_old.call(sx,sy-1) )
		P_Map.add(xi,yi,-2);
  	      
	      if( !Excursion_set_old.call(sx+1,sy) )
		P_Map.add(xi,yi,-2);
  
	      C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx+1,sy+1)*pow_of_2[2] + Excursion_set_old.call(sx+1,sy)*pow_of_2[4] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6] + Excursion_set_old.call(sx+1,sy-1)*pow_of_2[7]]);
  
	      Excursion_set_old.assign(sx,sy,false);
	    } 
	  }
	sy = end_y;
	{
	  if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	    A_Map.add(xi,yi,-1);    
	    if( !Excursion_set_old.call(sx,sy-1) )
	      P_Map.add(xi,yi,-2);
	    
	    if( !Excursion_set_old.call(sx+1,sy) )
	      P_Map.add(xi,yi,-2);
  
	    C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx+1,sy)*pow_of_2[4] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6] + Excursion_set_old.call(sx+1,sy-1)*pow_of_2[7]]);
  
	    Excursion_set_old.assign(sx,sy,false);
	  } 
	}
	
	for(sx = xi*slide+1; sx < end_x; sx++){
	  sy = yi*slide;
	  {
	    if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	      A_Map.add(xi,yi,-1);
	      if( !Excursion_set_old.call(sx,sy+1) )
		P_Map.add(xi,yi,-2);
	      
	      if( Excursion_set_old.call(sx+1,sy) && Excursion_set_old.call(sx-1,sy) )
		P_Map.add(xi,yi,2);    
	      else if( !Excursion_set_old.call(sx+1,sy) && !Excursion_set_old.call(sx-1,sy) )
		P_Map.add(xi,yi,-2);
  
	      C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy+1)*pow_of_2[0] + Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx+1,sy+1)*pow_of_2[2] + Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx+1,sy)*pow_of_2[4]]);
  
	      Excursion_set_old.assign(sx,sy,false);
	    } 
	  }
	  for(sy = yi*slide+1; sy < end_y; sy++)
	    {
	      if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
		A_Map.add(xi,yi,-1);
		if( Excursion_set_old.call(sx,sy+1) && Excursion_set_old.call(sx,sy-1) )
		  P_Map.add(xi,yi,2);    
		else if( !Excursion_set_old.call(sx,sy+1) && !Excursion_set_old.call(sx,sy-1) )
		  P_Map.add(xi,yi,-2);
  
		if( Excursion_set_old.call(sx+1,sy) && Excursion_set_old.call(sx-1,sy) )
		  P_Map.add(xi,yi,2);    
		else if( !Excursion_set_old.call(sx+1,sy) && !Excursion_set_old.call(sx-1,sy) )
		  P_Map.add(xi,yi,-2);
  
		C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy+1)*pow_of_2[0] + Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx+1,sy+1)*pow_of_2[2] + Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx+1,sy)*pow_of_2[4] + Excursion_set_old.call(sx-1,sy-1)*pow_of_2[5] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6] + Excursion_set_old.call(sx+1,sy-1)*pow_of_2[7]]);
  
		Excursion_set_old.assign(sx,sy,false);
	      } 
	    }
	  sy = end_y;
	  {
	    if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	      A_Map.add(xi,yi,-1);    
	      if( !Excursion_set_old.call(sx,sy-1) )
		P_Map.add(xi,yi,-2);
  
	      if( Excursion_set_old.call(sx+1,sy) && Excursion_set_old.call(sx-1,sy) )
		P_Map.add(xi,yi,2);    
	      else if( !Excursion_set_old.call(sx+1,sy) && !Excursion_set_old.call(sx-1,sy) )
		P_Map.add(xi,yi,-2);
  
	      C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx+1,sy)*pow_of_2[4] + Excursion_set_old.call(sx-1,sy-1)*pow_of_2[5] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6] + Excursion_set_old.call(sx+1,sy-1)*pow_of_2[7]]);
  
	      Excursion_set_old.assign(sx,sy,false);
	    } 
	  }
	}
	
	sx = end_x;
	sy = yi*slide;
	{
	  if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	    A_Map.add(xi,yi,-1);
	    if( !Excursion_set_old.call(sx,sy+1) )
	      P_Map.add(xi,yi,-2);
	    
	    if( !Excursion_set_old.call(sx-1,sy) )
	      P_Map.add(xi,yi,-2);
  
	    C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy+1)*pow_of_2[0] + Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx-1,sy)*pow_of_2[3]]);
  
	    Excursion_set_old.assign(sx,sy,false);
	  } 
	}
	for(sy = yi*slide+1; sy < end_y; sy++)
	  {
	    if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	      A_Map.add(xi,yi,-1);
	      if( Excursion_set_old.call(sx,sy+1) && Excursion_set_old.call(sx,sy-1) )
		P_Map.add(xi,yi,2);    
	      else if( !Excursion_set_old.call(sx,sy+1) && !Excursion_set_old.call(sx,sy-1) )
		P_Map.add(xi,yi,-2);
	      
	      if( !Excursion_set_old.call(sx-1,sy) )
		P_Map.add(xi,yi,-2);
  
	      C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy+1)*pow_of_2[0] + Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx-1,sy-1)*pow_of_2[5] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6]]);
  
	      Excursion_set_old.assign(sx,sy,false);
	    } 
	  }
	sy = end_y;
	{
	  if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	    A_Map.add(xi,yi,-1);    
	    if( !Excursion_set_old.call(sx,sy-1) )
	      P_Map.add(xi,yi,-2);
	    
	    if( !Excursion_set_old.call(sx-1,sy) )
	      P_Map.add(xi,yi,-2);
  
	    C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx-1,sy-1)*pow_of_2[5] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6]]);
  
	    Excursion_set_old.assign(sx,sy,false);
	  } 
	}
	// std::cout << A_Map.call(xi,yi) << "  " << P_Map.call(xi,yi) << "  " << C_Map.call(xi,yi) << "  " << std::endl;
	
	stateclassifier tmp_state(A_Map.call(xi,yi),P_Map.call(xi,yi),C_Map.call(xi,yi),0);
	int comp_APC_end = nmbr_of_macrostates-1;
	int tmp_sr = 0;
	for(int sr = comp_APC_end; sr >= 0; sr--)
	  if(tmp_state==comp_APC[sr]){
	    tmp_sr = sr;
	    break;
	  }
	
	Dev = comp_APC[tmp_sr].cl();
	if(Dev == 0){
	  Dev = min_comp;
	  std::cerr << "WARNING: unknown structure A = " << A_Map.call(xi,yi) << ", P = " << P_Map.call(xi,yi) << ", and C = " << C_Map.call(xi,yi) << std::endl;
	}
	Dev = log10(Dev);
	
	if( A_Map.call(xi,yi) >= p*slide2 )
	  Dev *= -1;
	/*
	MinkowskiSkyMap_sum.add(xi,yi,Dev);
	if(thresh==lambda)
	  MinkowskiSkyMap_single.assign(xi,yi,Dev);
	*/
       	if( fabs(MinkowskiSkyMap.call(xi,yi)) < fabs(Dev) )
	  MinkowskiSkyMap.assign(xi,yi,Dev);
	
	// std::cout << thresh << "    " << Dev << std::endl;
      }// for(xi,yi)
    
    // Excursion set
    for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
      if(Excursion_set_old.call(ni) != Excursion_set_new.call(ni)){
	std::cerr << "ERROR: Excursion sets do not match!" << std::endl;
	exit(-1);
      }
    
    if(verbose)
      Cout_One_Percent(thresh-thresh_min,percent,nmbr_of_thresh);
  }
  
  // Trial correction
  for(unsigned xi = 0; xi < Nx_mink; xi++)
    for(unsigned yi = 0; yi < Ny_mink; yi++){
      double D = MinkowskiSkyMap.call(xi,yi);
      double Dt = fabs(D) - log10(trial); 
      if(fabs(D) < 9)
	Dt = log10(1-pow(1-pow(0.1, fabs(D) ),trial));
      MinkowskiSkyMap.assign(xi,yi,copysign(Dt,D));
    }
  
  // Output of the Minkowski sky map:
  // Matrix
  if(matrix_out){
    std::stringstream matfilename_stst;
    matfilename_stst << filename
		     << "_Matrix_" << Nx_mink << "x" << Ny_mink << "_seed_" << seed_
		     << "_lambda_" << lambda << "_slide_" << slide << ".dat";
    std::string matfilename = matfilename_stst.str();
    
    MinkowskiSkyMap.MatrixToFile(matfilename,prefix_of_);
  }
  /*
  if(matrix_out){
    std::stringstream matfilename_stst;
    matfilename_stst << filename
		     << "_single_Matrix_" << Nx_mink << "x" << Ny_mink << "_seed_" << seed_
		     << "_lambda_" << lambda << "_slide_" << slide << ".dat";
    std::string matfilename = matfilename_stst.str();
    
    MinkowskiSkyMap_single.MatrixToFile(matfilename,prefix_of_);
  }
  if(matrix_out){
    std::stringstream matfilename_stst;
    matfilename_stst << filename
		     << "_sum_Matrix_" << Nx_mink << "x" << Ny_mink << "_seed_" << seed_
		     << "_lambda_" << lambda << "_slide_" << slide << ".dat";
    std::string matfilename = matfilename_stst.str();
    
    MinkowskiSkyMap_sum.MatrixToFile(matfilename,prefix_of_);
  }
  */
  if( verbose )
    std::cout << "-------------------" << std::endl;
  
  return MinkowskiSkyMap;
}

BinField<double> BinnedSkyMap::Minkowski_sky_map_sum_D_APC_pix_wbc(const double &lambda, const unsigned int &slide,
								   const bool &verbose, const bool &matrix_out, const bool &gnuplot_out,
								   const std::string &filename) const
{  
  // Define MinkowskiSkyMap
  unsigned Nx_mink = Nx_-slide+1;
  unsigned Ny_mink = Ny_-slide+1;
  BinField< double > MinkowskiSkyMap(Nx_mink,Ny_mink,0);
  unsigned slide2 = slide*slide;
  unsigned slide_minus_one = slide - 1;
  unsigned N_A = slide*slide + 1;
  
  // ECCDF
  double binwidth(0);
  int N_eccdf(0);
  if( slide == 15 && lambda == 100){
    binwidth = 0.020001000050;
    N_eccdf = 18792;
  }
  else if( slide ==  15 && lambda == 50){
    binwidth = 0.009000450023;
    N_eccdf = 19665;
  }
  else if( slide ==  15 && lambda == 25){
    binwidth = 0.007500375019;
    N_eccdf = 18932;
  }
  else if( slide ==  10 && lambda == 100){
    binwidth = 0.013250662533;
    N_eccdf = 19606;
  }
  else if( slide ==  10 && lambda == 50){
    binwidth = 0.009500475024;
    N_eccdf = 19789;
  }
  else if( slide ==  10 && lambda == 25){
    binwidth = 0.007500375019;
    N_eccdf = 19129;
  }
  else if( slide ==  5 && lambda == 100){
    binwidth = 0.015000750038;
    N_eccdf = 16973;
  }
  else if( slide ==  5 && lambda == 50){
    binwidth = 0.011000550028;
    N_eccdf = 18488;
  }
  else if( slide ==  5 && lambda == 25){
    binwidth = 0.006000300015;
    N_eccdf = 19971;
  }
  else{
    std::cerr << "ERROR: unknown parameters w = " << slide << " and lambda = " << lambda << std::endl;
    exit(-1);
  }
  double *neg_log10_eccdf;
  neg_log10_eccdf = load_neg_log10_eccdf(slide, lambda);
  
  // DoS: ready to read in
  std::vector< stateclassifier > DoS_APC; DoS_APC.clear();
  
  // Read in DoS
  if( verbose ){
    std::cout << "Reading DoS ..." << std::endl
	      << PercentageBar << std::endl;}
  unsigned percent = 0;
  for(unsigned Ai = 0; Ai < N_A; Ai++){
    std::stringstream infile_st;
    infile_st << "parameters/structure_" << slide << "x" << slide << "/DoS_averaged_A_"
	      << Ai
	      << "_PC_pix_wbc_" << slide << "x" << slide << "_seed_0_seedwght_0.dat";
    std::string infile = infile_st.str();
    std::ifstream structure_distribution((infile).c_str ());
    if(structure_distribution.fail()){
      std::cerr << "ERROR: ifstream failed to read " << "parameters/structure_" << slide << "x" << slide << "/DoS_averaged_A_"
		<< Ai << "_PC_pix_wbc_" << slide << "x" << slide << "_seed_0_seedwght_0.dat;" << std::endl;
      exit(-1);
    }
    while(!structure_distribution.eof()){
      int P(0), C(0);
      structure_distribution >> P;
      structure_distribution >> C;
      
      double DoS = 0;
      structure_distribution >> DoS;
      if(DoS)
	DoS_APC.push_back(stateclassifier(Ai,P,C,DoS));
    }
    structure_distribution.close();
    if(verbose) 
      Cout_One_Percent(Ai,percent,N_A);
  }
  unsigned nmbr_of_macrostates = DoS_APC.size();
  
  // Determine maximum and minimum thresholds
  unsigned thresh_max = 0;
  for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
    if(countrate_.call(ni) > thresh_max)
      thresh_max = countrate_.call(ni);
  thresh_max++;
  unsigned thresh_min = thresh_max;
  for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
    if(countrate_.call(ni) < thresh_min)
      thresh_min = countrate_.call(ni);
  unsigned nmbr_of_thresh = thresh_max-thresh_min+1;
  
  // Compute Minkowski sky map
  if( verbose ){
    std::cout << "Compute Minkowski sky map ..." << std::endl
	      << PercentageBar << std::endl;
  }
  percent = 0;
  
  // Prepare excursion set and scan-window functional values
  BinField< bool > Excursion_set_old(Nx_,Ny_,true);
  BinField< bool > Excursion_set_older(Nx_,Ny_,true);
  BinField< bool > Excursion_set_new(Nx_,Ny_,true);
  BinField< int > A_Map(Nx_mink,Ny_mink,slide2);
  BinField< int > P_Map(Nx_mink,Ny_mink,4*slide);
  BinField< int > C_Map(Nx_mink,Ny_mink,1);
  
  double p(0), invp(0), Dev(0);
  std::vector< stateclassifier > prob_APC(nmbr_of_macrostates, stateclassifier(0,0,0,1));
  std::vector< stateclassifier > comp_APC(nmbr_of_macrostates, stateclassifier(0,0,0,1));
  
  for(unsigned thresh = thresh_min; thresh <= thresh_max; thresh++){
    p = Bernoulli_p(thresh, lambda);
    invp = 1-p;
    
    // Excursion set
    for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
      Excursion_set_new.assign(ni, countrate_.call(ni) >= thresh);
    
    // determine probability distribution
    for(std::vector< stateclassifier >::size_type sr = 0; sr < nmbr_of_macrostates; sr++){
      int Ai = DoS_APC[sr].A();
      prob_APC[sr] = stateclassifier(Ai,DoS_APC[sr].P(),DoS_APC[sr].C(),DoS_APC[sr].cl()*pow(p,Ai)*pow(invp,slide2-Ai));
    }
    
    // determine compatibility distribution
    boost::sort(prob_APC);
    comp_APC=prob_APC;
    
    double min_comp = comp_APC[0].cl();
    for(std::vector< stateclassifier >::size_type sr = 1; sr < comp_APC.size(); sr++)
      comp_APC[sr].add_cl(comp_APC[sr-1].cl());
    
    // Deviation strength
    for(unsigned xi = 0; xi < Nx_mink; xi++)
      for(unsigned yi = 0; yi < Ny_mink; yi++){
	unsigned end_x = xi+slide_minus_one, end_y = yi+slide_minus_one;
	
	unsigned sx = xi;
	unsigned sy = yi;
	{
	  if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	    A_Map.add(xi,yi,-1);
	    if( !Excursion_set_old.call(sx,sy+1) )
	      P_Map.add(xi,yi,-2);
	    
	    if( !Excursion_set_old.call(sx+1,sy) )
	      P_Map.add(xi,yi,-2);
	    
	    C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx+1,sy+1)*pow_of_2[2]+ Excursion_set_old.call(sx+1,sy)*pow_of_2[4]]);
	    
	    Excursion_set_old.assign(sx,sy,false);
	  } 
	}
	for(sy = yi+1; sy < end_y; sy++)
	  {
	    if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	      A_Map.add(xi,yi,-1);
	      if( Excursion_set_old.call(sx,sy+1) && Excursion_set_old.call(sx,sy-1) )
		P_Map.add(xi,yi,2);    
	      else if( !Excursion_set_old.call(sx,sy+1) && !Excursion_set_old.call(sx,sy-1) )
		P_Map.add(xi,yi,-2);
  	      
	      if( !Excursion_set_old.call(sx+1,sy) )
		P_Map.add(xi,yi,-2);
  
	      C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx+1,sy+1)*pow_of_2[2] + Excursion_set_old.call(sx+1,sy)*pow_of_2[4] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6] + Excursion_set_old.call(sx+1,sy-1)*pow_of_2[7]]);
  
	      Excursion_set_old.assign(sx,sy,false);
	    } 
	  }
	sy = end_y;
	{
	  if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	    A_Map.add(xi,yi,-1);    
	    if( !Excursion_set_old.call(sx,sy-1) )
	      P_Map.add(xi,yi,-2);
	    
	    if( !Excursion_set_old.call(sx+1,sy) )
	      P_Map.add(xi,yi,-2);
  
	    C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx+1,sy)*pow_of_2[4] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6] + Excursion_set_old.call(sx+1,sy-1)*pow_of_2[7]]);
  
	    Excursion_set_old.assign(sx,sy,false);
	  } 
	}
	
	for(sx = xi+1; sx < end_x; sx++){
	  sy = yi;
	  {
	    if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	      A_Map.add(xi,yi,-1);
	      if( !Excursion_set_old.call(sx,sy+1) )
		P_Map.add(xi,yi,-2);
	      
	      if( Excursion_set_old.call(sx+1,sy) && Excursion_set_old.call(sx-1,sy) )
		P_Map.add(xi,yi,2);    
	      else if( !Excursion_set_old.call(sx+1,sy) && !Excursion_set_old.call(sx-1,sy) )
		P_Map.add(xi,yi,-2);
  
	      C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy+1)*pow_of_2[0] + Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx+1,sy+1)*pow_of_2[2] + Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx+1,sy)*pow_of_2[4]]);
  
	      Excursion_set_old.assign(sx,sy,false);
	    } 
	  }
	  for(sy = yi+1; sy < end_y; sy++)
	    {
	      if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
		A_Map.add(xi,yi,-1);
		if( Excursion_set_old.call(sx,sy+1) && Excursion_set_old.call(sx,sy-1) )
		  P_Map.add(xi,yi,2);    
		else if( !Excursion_set_old.call(sx,sy+1) && !Excursion_set_old.call(sx,sy-1) )
		  P_Map.add(xi,yi,-2);
  
		if( Excursion_set_old.call(sx+1,sy) && Excursion_set_old.call(sx-1,sy) )
		  P_Map.add(xi,yi,2);    
		else if( !Excursion_set_old.call(sx+1,sy) && !Excursion_set_old.call(sx-1,sy) )
		  P_Map.add(xi,yi,-2);
  
		C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy+1)*pow_of_2[0] + Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx+1,sy+1)*pow_of_2[2] + Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx+1,sy)*pow_of_2[4] + Excursion_set_old.call(sx-1,sy-1)*pow_of_2[5] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6] + Excursion_set_old.call(sx+1,sy-1)*pow_of_2[7]]);
  
		Excursion_set_old.assign(sx,sy,false);
	      } 
	    }
	  sy = end_y;
	  {
	    if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	      A_Map.add(xi,yi,-1);    
	      if( !Excursion_set_old.call(sx,sy-1) )
		P_Map.add(xi,yi,-2);
  
	      if( Excursion_set_old.call(sx+1,sy) && Excursion_set_old.call(sx-1,sy) )
		P_Map.add(xi,yi,2);    
	      else if( !Excursion_set_old.call(sx+1,sy) && !Excursion_set_old.call(sx-1,sy) )
		P_Map.add(xi,yi,-2);
  
	      C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx+1,sy)*pow_of_2[4] + Excursion_set_old.call(sx-1,sy-1)*pow_of_2[5] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6] + Excursion_set_old.call(sx+1,sy-1)*pow_of_2[7]]);
  
	      Excursion_set_old.assign(sx,sy,false);
	    } 
	  }
	}
	
	sx = end_x;
	sy = yi;
	{
	  if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	    A_Map.add(xi,yi,-1);
	    if( !Excursion_set_old.call(sx,sy+1) )
	      P_Map.add(xi,yi,-2);
	    
	    if( !Excursion_set_old.call(sx-1,sy) )
	      P_Map.add(xi,yi,-2);
  
	    C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy+1)*pow_of_2[0] + Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx-1,sy)*pow_of_2[3]]);
  
	    Excursion_set_old.assign(sx,sy,false);
	  } 
	}
	for(sy = yi+1; sy < end_y; sy++)
	  {
	    if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	      A_Map.add(xi,yi,-1);
	      if( Excursion_set_old.call(sx,sy+1) && Excursion_set_old.call(sx,sy-1) )
		P_Map.add(xi,yi,2);    
	      else if( !Excursion_set_old.call(sx,sy+1) && !Excursion_set_old.call(sx,sy-1) )
		P_Map.add(xi,yi,-2);
	      
	      if( !Excursion_set_old.call(sx-1,sy) )
		P_Map.add(xi,yi,-2);
  
	      C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy+1)*pow_of_2[0] + Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx-1,sy-1)*pow_of_2[5] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6]]);
  
	      Excursion_set_old.assign(sx,sy,false);
	    } 
	  }
	sy = end_y;
	{
	  if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	    A_Map.add(xi,yi,-1);    
	    if( !Excursion_set_old.call(sx,sy-1) )
	      P_Map.add(xi,yi,-2);
	    
	    if( !Excursion_set_old.call(sx-1,sy) )
	      P_Map.add(xi,yi,-2);
  
	    C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx-1,sy-1)*pow_of_2[5] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6]]);
  
	    Excursion_set_old.assign(sx,sy,false);
	  } 
	}
	// std::cout << A_Map.call(xi,yi) << "  " << P_Map.call(xi,yi) << "  " << C_Map.call(xi,yi) << "  " << std::endl;
	
	stateclassifier tmp_state(A_Map.call(xi,yi),P_Map.call(xi,yi),C_Map.call(xi,yi),0);
	int comp_APC_end = nmbr_of_macrostates-1;
	int tmp_sr = 0;
	for(int sr = comp_APC_end; sr >= 0; sr--)
	  if(tmp_state==comp_APC[sr]){
	    tmp_sr = sr;
	    break;
	  }
	
	Dev = comp_APC[tmp_sr].cl();
	if(Dev == 0){
	  Dev = min_comp;
	  std::cerr << "WARNING: unknown structure at (" << xi << "," << yi << "), A = " 
		    << A_Map.call(xi,yi) << ", P = " << P_Map.call(xi,yi) << ", and C = " << C_Map.call(xi,yi) << std::endl;
	}
	Dev = log10(Dev);
	
	//if( A_Map.call(xi,yi) >= p*slide2 )
	Dev *= -1;
	
	MinkowskiSkyMap.add(xi,yi,Dev);
	
	// std::cout << thresh << "    " << Dev << std::endl;
	
	// Excursion set
	for(sx = xi; sx <= end_x; sx++)
	  for(sy = yi; sy <= end_y; sy++)
	    Excursion_set_old.assign(sx, sy, Excursion_set_older.call(sx, sy));
      }// for(xi,yi)
    
    // Excursion set
    for(unsigned ni = 0; ni < Nx_*Ny_; ni++){
      Excursion_set_old.assign(ni, Excursion_set_new.call(ni));
      Excursion_set_older.assign(ni, Excursion_set_new.call(ni));
    }
    
    if(verbose)
      Cout_One_Percent(thresh-thresh_min,percent,nmbr_of_thresh);
  }
  
  double T = 9.0;
  for(unsigned xi = 0; xi < Nx_mink; xi++)
    for(unsigned yi = 0; yi < Ny_mink; yi++){
      double S=MinkowskiSkyMap.call(xi,yi);
      int sum_i = floor(S/binwidth);
      if(sum_i < N_eccdf){
	T = neg_log10_eccdf[sum_i];
      }
      else{
	std::cerr << "WARNING: too significant result" << std::endl;
      }
      /*
      // Determine sign:
      double sum_k = 0;
      for(unsigned txi = 0; txi < slide; txi++)
	for(unsigned tyi = 0; tyi < slide; tyi++)
	  sum_k += countrate_.call(xi+txi,yi+tyi);
      if( sum_k < lambda*slide2)
	MinkowskiSkyMap.assign(xi,yi,(-1)*T);
      else
      */
      MinkowskiSkyMap.assign(xi,yi,T);
    }
  
  double min_D = MinkowskiSkyMap.call(0,0);
  double max_D = min_D;
  for(unsigned xi = 0; xi < Nx_mink; xi++)
    for(unsigned yi = 0; yi < Ny_mink; yi++){
      double D = MinkowskiSkyMap.call(xi,yi);
      if ( D < min_D )
	min_D = D;
      else if ( D > max_D )
	max_D = D;
    }
  
  // Output of the Minkowski sky map:
  // Matrix
  if(matrix_out){
    std::stringstream matfilename_stst;
    matfilename_stst << filename
		     << "_Matrix_" << (Nx_-slide+1) << "x" << (Ny_-slide+1) << "_seed_" << seed_
		     << "_lambda_" << lambda << "_slide_" << slide << ".dat";
    std::string matfilename = matfilename_stst.str();
    
    MinkowskiSkyMap.MatrixToFile(matfilename,prefix_of_);
  }
  
  if(gnuplot_out){
    // Data output
    std::stringstream filename_stst;
    filename_stst << prefix_of_ << filename
		  << "_" << (Nx_-slide+1) << "x" << (Ny_-slide+1) << "_seed_" << seed_
		  << "_lambda_" << lambda << "_slide_" << slide << ".dat";
    std::string filename_st = filename_stst.str();
    std::ofstream Data( filename_st.c_str() );
    
    Data << "# Minkowski sky map" << std::endl << "#" << std::endl
	 << "# Background intensity: lambda = " << lambda << ";" << std::endl
	 << "# Functionals: APC;" << std::endl
	 << "# Size of sliding window: slide = " << slide << ";" << std::endl
	 << "# White boundary condition;" << std::endl
	 << "# Pixelized data" << std::endl << "#" << std::endl;
    Data << "#" << std::setw(13) << "X" << std::setw(14) << "Y" << std::setw(14) << "D(A, P, chi)" << std::endl;
    for(unsigned xi = 0; xi < (Nx_-slide+1); xi++)
      {
	for(unsigned yi = 0; yi < (Ny_-slide+1); yi++)
	  Data << std::setw(14) << xlow_.call(xi,yi) + (slide-1)*xbinlength_/2.0
	       << std::setw(14) << ylow_.call(xi,yi) + (slide-1)*ybinlength_/2.0
	       << std::setw(14) << MinkowskiSkyMap.call(xi,yi)
	       << std::endl
	       << std::setw(14) << xlow_.call(xi,yi) + (slide-1)*xbinlength_/2.0
	       << std::setw(14) << ylow_.call(xi,yi) + ybinlength_ + (slide-1)*ybinlength_/2.0
	       << std::setw(14) << MinkowskiSkyMap.call(xi,yi)
	       << std::endl;
      
	Data << std::endl;
      
	for(unsigned yi = 0; yi < (Ny_-slide+1); yi++)
	  Data << std::setw(14) << xlow_.call(xi,yi) + xbinlength_ + (slide-1)*xbinlength_/2.0
	       << std::setw(14) << ylow_.call(xi,yi) + (slide-1)*ybinlength_/2.0
	       << std::setw(14) << MinkowskiSkyMap.call(xi,yi)
	       << std::endl
	       << std::setw(14) << xlow_.call(xi,yi) + xbinlength_ + (slide-1)*xbinlength_/2.0
	       << std::setw(14) << ylow_.call(xi,yi) + ybinlength_ + (slide-1)*ybinlength_/2.0
	       << std::setw(14) << MinkowskiSkyMap.call(xi,yi)
	       << std::endl;
      
	Data << std::endl;
      }
    Data.close();

    // Gnuplot
    std::stringstream gpfilename_stst;
    gpfilename_stst << prefix_of_ << filename
		    << "_" << (Nx_-slide+1) << "x" << (Ny_-slide+1) << "_seed_" << seed_
		    << "_slide_" << slide << ".gp";
      //<< "_lambda_" << lambda << "_slide_" << slide << ".gp";
    std::string gpfilename_st = gpfilename_stst.str();
    std::ofstream Gnuplot( gpfilename_st.c_str() );
  
    Gnuplot << "set terminal epslatex standalone color colortext 14" << std::endl
	    << "set output '" << filename
	    << "_" << (Nx_-slide+1) << "x" << (Ny_-slide+1) << "_seed_" << seed_
	    << "_slide_" << slide << ".tex'";
    //<< "_lambda_" << lambda << "_slide_" << slide << ".tex'";
    Gnuplot << std::endl << std::endl
	    << "xmin = " << xmin_ + (slide-1)*xbinlength_/2.0 << std::endl
	    << "xmax = " << xmax_ - (slide-1)*xbinlength_/2.0 << std::endl
	    << "ymin = " << ymin_ + (slide-1)*ybinlength_/2.0 << std::endl
	    << "ymax = " << ymax_ - (slide-1)*ybinlength_/2.0 << std::endl
	    << std::endl
	    << "#set title 'Size of sliding observation window: " << slide << "with wbc'" << std::endl
	    << "unset xtics #border mirror norotate offset 0,0.3" << std::endl
	    << "unset ytics #border mirror norotate offset 0.3,0" << std::endl
	    << "#set cbtics 10" << std::endl
	    << std::endl
	    << "set mxtics 2" << std::endl
	    << "set mytics 2" << std::endl
	    << std::endl;
    if(Ny_ == Nx_)
      Gnuplot << "set size ratio 1" << std::endl;
    else
      Gnuplot << "set size ratio " << (Ny_-slide+1)*1.0/(Nx_-slide+1) << std::endl;
    Gnuplot << "set border lw 2" << std::endl
	    << std::endl
	    << "unset xlabel #'$x$' offset 0,0.7" << std::endl
	    << "unset ylabel #'$y$' offset 0.6,0" << std::endl
	    << "# set cblabel '\\large $\\mathcal{D}(A,P,\\chi)$' offset 0,0" << std::endl
	    << std::endl
	    << "set xrange [xmin:xmax]" << std::endl
	    << "set yrange [ymin:ymax]" << std::endl
	    << "set cbrange [0:9]" << std::endl
	    << "#set cbrange [" << min_D << ":" << max_D << "]" << std::endl
	    << std::endl
	    << "set pm3d map" << std::endl
	    << "# set palette defined (" << min_D << " 'black', " << 0.5*min_D << " 'dark-blue', " << 0 << " 'white', " << 0.2*max_D << " 'purple', " << 0.4*max_D << " 'dark-red', " << 0.6*max_D << " 'red', " << 0.8*max_D << " 'orange', " << max_D << " 'yellow')" << std::endl
	    << "splot '" << filename << "_" << (Nx_-slide+1) << "x" << (Ny_-slide+1) << "_seed_" << seed_
	    << "_lambda_" << lambda << "_slide_" << slide << ".dat";
    Gnuplot << "' notitle" << std::endl
	    << "#" << std::endl
	    << std::endl;
    
    Gnuplot.close();
  }
  
  if( verbose )
    std::cout << "-------------------" << std::endl;
  
  return MinkowskiSkyMap;
}

BinField<double> BinnedSkyMap::Minkowski_sky_map_sum_D_APC_pix_wbc_chessboard(const double &lambda, const unsigned int &slide,
									      const bool &verbose, const bool &matrix_out,
									      const std::string &filename) const
{
  // Define MinkowskiSkyMap
  unsigned Nx_mink = Nx_/slide;
  unsigned Ny_mink = Ny_/slide;
  BinField< double > MinkowskiSkyMap(Nx_mink,Ny_mink,0);
  //BinField< double > MinkowskiSkyMap_sum(Nx_mink,Ny_mink,0);
  //BinField< double > MinkowskiSkyMap_single(Nx_mink,Ny_mink,0);
  unsigned slide2 = slide*slide;
  unsigned slide_minus_one = slide - 1;
  unsigned N_A = slide*slide + 1;
  
  // ECCDF
  double binwidth(0);
  int N_eccdf(0);
  if( slide == 15 && lambda == 100){
    binwidth = 0.020001000050;
    N_eccdf = 18792;
  }
  else if( slide ==  15 && lambda == 50){
    binwidth = 0.009000450023;
    N_eccdf = 19665;
  }
  else if( slide ==  15 && lambda == 25){
    binwidth = 0.007500375019;
    N_eccdf = 18932;
  }
  else if( slide ==  10 && lambda == 100){
    binwidth = 0.013250662533;
    N_eccdf = 19606;
  }
  else if( slide ==  10 && lambda == 50){
    binwidth = 0.009500475024;
    N_eccdf = 19789;
  }
  else if( slide ==  10 && lambda == 25){
    binwidth = 0.007500375019;
    N_eccdf = 19129;
  }
  else if( slide ==  5 && lambda == 100){
    binwidth = 0.015000750038;
    N_eccdf = 16973;
  }
  else if( slide ==  5 && lambda == 50){
    binwidth = 0.011000550028;
    N_eccdf = 18488;
  }
  else if( slide ==  5 && lambda == 25){
    binwidth = 0.006000300015;
    N_eccdf = 19971;
  }
  else{
    std::cerr << "ERROR: unknown parameters w = " << slide << " and lambda = " << lambda << std::endl;
    exit(-1);
  }
  double *neg_log10_eccdf;
  neg_log10_eccdf = load_neg_log10_eccdf(slide, lambda);
  
  // DoS: ready to read in
  std::vector< stateclassifier > DoS_APC; DoS_APC.clear();
  
  // Read in DoS
  if( verbose ){
    std::cout << "Reading DoS ..." << std::endl
	      << PercentageBar << std::endl;}
  unsigned percent = 0;
  for(unsigned Ai = 0; Ai < N_A; Ai++){
    std::stringstream infile_st;
    infile_st << "parameters/structure_" << slide << "x" << slide << "/DoS_averaged_A_"
	      << Ai
	      << "_PC_pix_wbc_" << slide << "x" << slide << "_seed_0_seedwght_0.dat";
    std::string infile = infile_st.str();
    std::ifstream structure_distribution((infile).c_str ());
    if(structure_distribution.fail()){
      std::cerr << "ERROR: ifstream failed to read " << "parameters/structure_" << slide << "x" << slide << "/DoS_averaged_A_"
		<< Ai << "_PC_pix_wbc_" << slide << "x" << slide << "_seed_0_seedwght_0.dat;" << std::endl;
      exit(-1);
    }
    while(!structure_distribution.eof()){
      int P(0), C(0);
      structure_distribution >> P;
      structure_distribution >> C;
      
      double DoS = 0;
      structure_distribution >> DoS;
      if(DoS)
	DoS_APC.push_back(stateclassifier(Ai,P,C,DoS));
    }
    structure_distribution.close();
    if(verbose) 
      Cout_One_Percent(Ai,percent,N_A);
  }
  unsigned nmbr_of_macrostates = DoS_APC.size();
  
  // Determine maximum and minimum thresholds
  unsigned thresh_max = 0;
  for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
    if(countrate_.call(ni) > thresh_max)
      thresh_max = countrate_.call(ni);
  thresh_max++;
  unsigned thresh_min = thresh_max;
  for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
    if(countrate_.call(ni) < thresh_min)
      thresh_min = countrate_.call(ni);
  unsigned nmbr_of_thresh = thresh_max-thresh_min+1;
  
  // Compute Minkowski sky map
  if( verbose ){
    std::cout << "Compute Minkowski sky map ..." << std::endl
	      << PercentageBar << std::endl;
  }
  percent = 0;
  
  // Prepare excursion set and scan-window functional values
  BinField< bool > Excursion_set_old(Nx_,Ny_,true);
  BinField< bool > Excursion_set_new(Nx_,Ny_,true);
  BinField< int > A_Map(Nx_mink,Ny_mink,slide2);
  BinField< int > P_Map(Nx_mink,Ny_mink,4*slide);
  BinField< int > C_Map(Nx_mink,Ny_mink,1);
  
  double p(0), invp(0), Dev(0);
  std::vector< stateclassifier > prob_APC(nmbr_of_macrostates, stateclassifier(0,0,0,1));
  std::vector< stateclassifier > comp_APC(nmbr_of_macrostates, stateclassifier(0,0,0,1));
  
  for(unsigned thresh = thresh_min; thresh <= thresh_max; thresh++){
    p = Bernoulli_p(thresh, lambda);
    invp = 1-p;
    
    // Excursion set
    for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
      Excursion_set_new.assign(ni, countrate_.call(ni) >= thresh);
    
    // determine probability distribution
    for(std::vector< stateclassifier >::size_type sr = 0; sr < nmbr_of_macrostates; sr++){
      int Ai = DoS_APC[sr].A();
      prob_APC[sr] = stateclassifier(Ai,DoS_APC[sr].P(),DoS_APC[sr].C(),DoS_APC[sr].cl()*pow(p,Ai)*pow(invp,slide2-Ai));
    }
    
    // determine compatibility distribution
    boost::sort(prob_APC);
    comp_APC=prob_APC;
    
    double min_comp = comp_APC[0].cl();
    for(std::vector< stateclassifier >::size_type sr = 1; sr < comp_APC.size(); sr++)
      comp_APC[sr].add_cl(comp_APC[sr-1].cl());
    
    // Deviation strength
    for(unsigned xi = 0; xi < Nx_mink; xi++)
      for(unsigned yi = 0; yi < Ny_mink; yi++){
	unsigned end_x = xi*slide+slide_minus_one, end_y = yi*slide+slide_minus_one;
	
	unsigned sx = xi*slide;
	unsigned sy = yi*slide;
	{
	  if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	    A_Map.add(xi,yi,-1);
	    if( !Excursion_set_old.call(sx,sy+1) )
	      P_Map.add(xi,yi,-2);
	    
	    if( !Excursion_set_old.call(sx+1,sy) )
	      P_Map.add(xi,yi,-2);
	    
	    C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx+1,sy+1)*pow_of_2[2]+ Excursion_set_old.call(sx+1,sy)*pow_of_2[4]]);
	    
	    Excursion_set_old.assign(sx,sy,false);
	  } 
	}
	for(sy = yi*slide+1; sy < end_y; sy++)
	  {
	    if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	      A_Map.add(xi,yi,-1);
	      if( Excursion_set_old.call(sx,sy+1) && Excursion_set_old.call(sx,sy-1) )
		P_Map.add(xi,yi,2);    
	      else if( !Excursion_set_old.call(sx,sy+1) && !Excursion_set_old.call(sx,sy-1) )
		P_Map.add(xi,yi,-2);
  	      
	      if( !Excursion_set_old.call(sx+1,sy) )
		P_Map.add(xi,yi,-2);
  
	      C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx+1,sy+1)*pow_of_2[2] + Excursion_set_old.call(sx+1,sy)*pow_of_2[4] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6] + Excursion_set_old.call(sx+1,sy-1)*pow_of_2[7]]);
  
	      Excursion_set_old.assign(sx,sy,false);
	    } 
	  }
	sy = end_y;
	{
	  if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	    A_Map.add(xi,yi,-1);    
	    if( !Excursion_set_old.call(sx,sy-1) )
	      P_Map.add(xi,yi,-2);
	    
	    if( !Excursion_set_old.call(sx+1,sy) )
	      P_Map.add(xi,yi,-2);
  
	    C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx+1,sy)*pow_of_2[4] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6] + Excursion_set_old.call(sx+1,sy-1)*pow_of_2[7]]);
  
	    Excursion_set_old.assign(sx,sy,false);
	  } 
	}
	
	for(sx = xi*slide+1; sx < end_x; sx++){
	  sy = yi*slide;
	  {
	    if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	      A_Map.add(xi,yi,-1);
	      if( !Excursion_set_old.call(sx,sy+1) )
		P_Map.add(xi,yi,-2);
	      
	      if( Excursion_set_old.call(sx+1,sy) && Excursion_set_old.call(sx-1,sy) )
		P_Map.add(xi,yi,2);    
	      else if( !Excursion_set_old.call(sx+1,sy) && !Excursion_set_old.call(sx-1,sy) )
		P_Map.add(xi,yi,-2);
  
	      C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy+1)*pow_of_2[0] + Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx+1,sy+1)*pow_of_2[2] + Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx+1,sy)*pow_of_2[4]]);
  
	      Excursion_set_old.assign(sx,sy,false);
	    } 
	  }
	  for(sy = yi*slide+1; sy < end_y; sy++)
	    {
	      if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
		A_Map.add(xi,yi,-1);
		if( Excursion_set_old.call(sx,sy+1) && Excursion_set_old.call(sx,sy-1) )
		  P_Map.add(xi,yi,2);    
		else if( !Excursion_set_old.call(sx,sy+1) && !Excursion_set_old.call(sx,sy-1) )
		  P_Map.add(xi,yi,-2);
  
		if( Excursion_set_old.call(sx+1,sy) && Excursion_set_old.call(sx-1,sy) )
		  P_Map.add(xi,yi,2);    
		else if( !Excursion_set_old.call(sx+1,sy) && !Excursion_set_old.call(sx-1,sy) )
		  P_Map.add(xi,yi,-2);
  
		C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy+1)*pow_of_2[0] + Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx+1,sy+1)*pow_of_2[2] + Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx+1,sy)*pow_of_2[4] + Excursion_set_old.call(sx-1,sy-1)*pow_of_2[5] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6] + Excursion_set_old.call(sx+1,sy-1)*pow_of_2[7]]);
  
		Excursion_set_old.assign(sx,sy,false);
	      } 
	    }
	  sy = end_y;
	  {
	    if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	      A_Map.add(xi,yi,-1);    
	      if( !Excursion_set_old.call(sx,sy-1) )
		P_Map.add(xi,yi,-2);
  
	      if( Excursion_set_old.call(sx+1,sy) && Excursion_set_old.call(sx-1,sy) )
		P_Map.add(xi,yi,2);    
	      else if( !Excursion_set_old.call(sx+1,sy) && !Excursion_set_old.call(sx-1,sy) )
		P_Map.add(xi,yi,-2);
  
	      C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx+1,sy)*pow_of_2[4] + Excursion_set_old.call(sx-1,sy-1)*pow_of_2[5] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6] + Excursion_set_old.call(sx+1,sy-1)*pow_of_2[7]]);
  
	      Excursion_set_old.assign(sx,sy,false);
	    } 
	  }
	}
	
	sx = end_x;
	sy = yi*slide;
	{
	  if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	    A_Map.add(xi,yi,-1);
	    if( !Excursion_set_old.call(sx,sy+1) )
	      P_Map.add(xi,yi,-2);
	    
	    if( !Excursion_set_old.call(sx-1,sy) )
	      P_Map.add(xi,yi,-2);
  
	    C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy+1)*pow_of_2[0] + Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx-1,sy)*pow_of_2[3]]);
  
	    Excursion_set_old.assign(sx,sy,false);
	  } 
	}
	for(sy = yi*slide+1; sy < end_y; sy++)
	  {
	    if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	      A_Map.add(xi,yi,-1);
	      if( Excursion_set_old.call(sx,sy+1) && Excursion_set_old.call(sx,sy-1) )
		P_Map.add(xi,yi,2);    
	      else if( !Excursion_set_old.call(sx,sy+1) && !Excursion_set_old.call(sx,sy-1) )
		P_Map.add(xi,yi,-2);
	      
	      if( !Excursion_set_old.call(sx-1,sy) )
		P_Map.add(xi,yi,-2);
  
	      C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy+1)*pow_of_2[0] + Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx-1,sy-1)*pow_of_2[5] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6]]);
  
	      Excursion_set_old.assign(sx,sy,false);
	    } 
	  }
	sy = end_y;
	{
	  if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	    A_Map.add(xi,yi,-1);    
	    if( !Excursion_set_old.call(sx,sy-1) )
	      P_Map.add(xi,yi,-2);
	    
	    if( !Excursion_set_old.call(sx-1,sy) )
	      P_Map.add(xi,yi,-2);
  
	    C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx-1,sy-1)*pow_of_2[5] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6]]);
  
	    Excursion_set_old.assign(sx,sy,false);
	  } 
	}
	// std::cout << A_Map.call(xi,yi) << "  " << P_Map.call(xi,yi) << "  " << C_Map.call(xi,yi) << "  " << std::endl;
	
	stateclassifier tmp_state(A_Map.call(xi,yi),P_Map.call(xi,yi),C_Map.call(xi,yi),0);
	int comp_APC_end = nmbr_of_macrostates-1;
	int tmp_sr = 0;
	for(int sr = comp_APC_end; sr >= 0; sr--)
	  if(tmp_state==comp_APC[sr]){
	    tmp_sr = sr;
	    break;
	  }
	
	Dev = comp_APC[tmp_sr].cl();
	if(Dev == 0){
	  Dev = min_comp;
	  std::cerr << "WARNING: unknown structure A = " << A_Map.call(xi,yi) << ", P = " << P_Map.call(xi,yi) << ", and C = " << C_Map.call(xi,yi) << std::endl;
	}
	Dev = log10(Dev);
	
	//if( A_Map.call(xi,yi) >= p*slide2 )
	Dev *= -1;
	
	MinkowskiSkyMap.add(xi,yi,Dev);
	
	/*
	MinkowskiSkyMap_sum.add(xi,yi,Dev);
	if(thresh==lambda)
	  MinkowskiSkyMap_single.assign(xi,yi,Dev);
	*/
	/*
       	if( fabs(MinkowskiSkyMap.call(xi,yi)) < fabs(Dev) )
	  MinkowskiSkyMap.assign(xi,yi,Dev);
	*/
	// std::cout << thresh << "    " << Dev << std::endl;
      }// for(xi,yi)
    
    // Excursion set
    for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
      if(Excursion_set_old.call(ni) != Excursion_set_new.call(ni)){
	std::cerr << "ERROR: Excursion sets do not match!" << std::endl;
	exit(-1);
      }
    
    if(verbose)
      Cout_One_Percent(thresh-thresh_min,percent,nmbr_of_thresh);
  }
  
  double T = 9.0;
  for(unsigned xi = 0; xi < Nx_mink; xi++)
    for(unsigned yi = 0; yi < Ny_mink; yi++){
      double S=MinkowskiSkyMap.call(xi,yi);
      int sum_i = floor(S/binwidth);
      if(sum_i < N_eccdf){
	T = neg_log10_eccdf[sum_i];
      }
      else{
	std::cerr << "WARNING: too significant result" << std::endl;
      }
      
      /*
      // Determine sign:
      double sum_k = 0;
      for(unsigned txi = 0; txi < slide; txi++)
	for(unsigned tyi = 0; tyi < slide; tyi++)
	  sum_k += countrate_.call(xi+txi,yi+tyi);
      if( sum_k < lambda*slide2)
	MinkowskiSkyMap.assign(xi,yi,(-1)*T);
      else
      */
      MinkowskiSkyMap.assign(xi,yi,T);
    }
  
  // Output of the Minkowski sky map:
  // Matrix
  if(matrix_out){
    std::stringstream matfilename_stst;
    matfilename_stst << filename
		     << "_Matrix_" << Nx_mink << "x" << Ny_mink << "_seed_" << seed_
		     << "_lambda_" << lambda << "_slide_" << slide << ".dat";
    std::string matfilename = matfilename_stst.str();
    
    MinkowskiSkyMap.MatrixToFile(matfilename,prefix_of_);
  }
  /*
  if(matrix_out){
    std::stringstream matfilename_stst;
    matfilename_stst << filename
		     << "_single_Matrix_" << Nx_mink << "x" << Ny_mink << "_seed_" << seed_
		     << "_lambda_" << lambda << "_slide_" << slide << ".dat";
    std::string matfilename = matfilename_stst.str();
    
    MinkowskiSkyMap_single.MatrixToFile(matfilename,prefix_of_);
  }
  if(matrix_out){
    std::stringstream matfilename_stst;
    matfilename_stst << filename
		     << "_sum_Matrix_" << Nx_mink << "x" << Ny_mink << "_seed_" << seed_
		     << "_lambda_" << lambda << "_slide_" << slide << ".dat";
    std::string matfilename = matfilename_stst.str();
    
    MinkowskiSkyMap_sum.MatrixToFile(matfilename,prefix_of_);
  }
  */
  if( verbose )
    std::cout << "-------------------" << std::endl;
  
  return MinkowskiSkyMap;
}

BinField<double> BinnedSkyMap::Minkowski_sky_map_A_pix_wbc_chessboard(const double &lambda, const unsigned int &slide, const bool &verbose,
								      const std::string &filename) const
{
  // Define MinkowskiSkyMap
  unsigned Nx_mink = Nx_/slide;
  unsigned Ny_mink = Ny_/slide;
  BinField< double > MinkowskiSkyMap(Nx_mink,Ny_mink,0);
  unsigned slide2 = slide*slide, nmbrblckpxl = slide2+1;
  
  // Determine maximum threshold
  unsigned thresh_max = 0;
  for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
    if(countrate_.call(ni) > thresh_max)
      thresh_max = countrate_.call(ni);
  thresh_max++;
  unsigned thresh_min = thresh_max;
  for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
    if(countrate_.call(ni) < thresh_min)
      thresh_min = countrate_.call(ni);
  unsigned nmbr_of_thresh = thresh_max-thresh_min+1;
  
  // Compute Minkowski sky map
  BinField< bool > ExcursionSet(Nx_,Ny_,true); unsigned BlackPixels = 0, percent = 0;
  double p, mean, Dev;
  if(verbose) 
    std::cout << "Compute Minkowski sky map ..." << std::endl
	      << PercentageBar << std::endl;
  for(unsigned thresh = thresh_min; thresh <= thresh_max; thresh++){
    p = Bernoulli_p(thresh, lambda);
    mean = slide2*p;
    
    BlackPixels = 0;
    for(unsigned ni = 0; ni < Nx_*Ny_; ni++){
      ExcursionSet.assign(ni, countrate_.call(ni) >= thresh);
      BlackPixels += ExcursionSet.call(ni);
    }
    
    for(unsigned xi = 0; xi < Nx_mink; xi++)
      for(unsigned yi = 0; yi < Ny_mink; yi++){
	unsigned A = 0;
	double C = 0;
	for(unsigned sx = xi*slide; sx < (xi+1)*slide; sx++)
	  for(unsigned sy = yi*slide; sy < (yi+1)*slide; sy++)
	    A += ExcursionSet.call(sx,sy);
	
	if(A > mean){
	  C = BinomialCDF_complement(slide2, A, p);
	  unsigned bound = nmbrblckpxl;
	  
	  for(unsigned m = 0; m <= mean; m++)
	    if( BinomialPDF(slide2, m, p) <= BinomialPDF(slide2, A, p) )
	      bound = m;
	    else
	      break;
	  if(bound < nmbrblckpxl)
	    C += BinomialCDF(slide2, bound, p);
	}
	else{      
	  C = BinomialCDF(slide2, A, p);
	  unsigned bound = nmbrblckpxl;
	  
	  for(unsigned m = A+1; m < nmbrblckpxl; m++)
	    if( BinomialPDF(slide2, m, p) <= BinomialPDF(slide2, A, p) ){
	      bound = m;
	      break;}
	  if(bound < nmbrblckpxl)
	    C += BinomialCDF_complement(slide2, bound, p);
	}
      
	Dev = log10(C);
	
	if( A >= mean )
	  Dev *= -1;
	
	if( fabs(MinkowskiSkyMap.call(xi,yi)) < fabs(Dev) ){
	  MinkowskiSkyMap.assign(xi,yi,Dev);
	}
      }
    
    if(verbose) 
      Cout_One_Percent(thresh-thresh_min,percent,nmbr_of_thresh);
  }
  
  if(verbose) 
    std::cout << std::endl; 
  
  // Trial correction
  for(unsigned xi = 0; xi < Nx_mink; xi++)
    for(unsigned yi = 0; yi < Ny_mink; yi++){
      double D = MinkowskiSkyMap.call(xi,yi);
      double Dt = fabs(D) - log10(trial); 
      if(fabs(D) < 9)
	Dt = log10(1-pow(1-pow(0.1, fabs(D) ),trial));
      MinkowskiSkyMap.assign(xi,yi,copysign(Dt,D));
    }
  
  // Output of the Minkowski sky map:
  // Matrix
  std::stringstream matfilename_stst;
  matfilename_stst << filename
		   << "_Matrix_" << Nx_mink << "x" << Ny_mink << "_seed_" << seed_
		   << "_lambda_" << lambda << "_slide_" << slide << ".dat";
  std::string matfilename = matfilename_stst.str();
  
  MinkowskiSkyMap.MatrixToFile(matfilename,prefix_of_);
  
  if(verbose) 
    std::cout << "-------------------" << std::endl;
  
  return MinkowskiSkyMap;
}

BinField<double> BinnedSkyMap::Minkowski_sky_map_A_pix_wbc(const double &lambda, const unsigned int &slide, const bool &verbose,
							   const std::string &filename) const
{  
  // Define MinkowskiSkyMap
  BinField< double > MinkowskiSkyMap(Nx_-slide+1,Ny_-slide+1,0);
  unsigned slide2 = slide*slide, nmbrblckpxl = slide2+1;
  
  // Determine maximum threshold
  unsigned thresh_max = 0;
  for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
    if(countrate_.call(ni) > thresh_max)
      thresh_max = countrate_.call(ni);
  thresh_max++;
  unsigned thresh_min = thresh_max;
  for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
    if(countrate_.call(ni) < thresh_min)
      thresh_min = countrate_.call(ni);
  unsigned nmbr_of_thresh = thresh_max-thresh_min+1;
  
  // Compute Minkowski sky map
  BinField< bool > ExcursionSet(Nx_,Ny_,true); unsigned BlackPixels = 0, percent = 0;
  double p, mean, Dev;
  if(verbose) 
    std::cout << "Compute Minkowski sky map ..." << std::endl
	      << PercentageBar << std::endl;
  for(unsigned thresh = thresh_min; thresh <= thresh_max; thresh++){
    BlackPixels = 0;
    p = Bernoulli_p(thresh, lambda);
    mean = slide2*p;
    for(unsigned ni = 0; ni < Nx_*Ny_; ni++){
      ExcursionSet.assign(ni, countrate_.call(ni) >= thresh);
      BlackPixels += ExcursionSet.call(ni);
    }
    
    for(unsigned xi = 0; xi < (Nx_-slide+1); xi++){
      unsigned A = 0;
      for(unsigned i = 0; i < slide2; i++)
	A += ExcursionSet.call(xi+slide-1-i+(i/slide)*slide,i/slide);
      
      for(unsigned yi = 0; yi < (Ny_-slide+1); yi++){
	double C = 0;
	if(yi > 0)
	  for(unsigned i = 0; i < slide; i++){
	    A -= ExcursionSet.call(xi+i,yi-1);
	    A += ExcursionSet.call(xi+i,yi+slide-1);
	  }
	
	if(A > mean){
	  C = BinomialCDF_complement(slide2, A, p);
	  unsigned bound = nmbrblckpxl;
	  
	  for(unsigned m = 0; m <= mean; m++)
	    if( BinomialPDF(slide2, m, p) <= BinomialPDF(slide2, A, p) )
	      bound = m;
	    else
	      break;
	  if(bound < nmbrblckpxl)
	    C += BinomialCDF(slide2, bound, p);
	}
	else{      
	  C = BinomialCDF(slide2, A, p);
	  unsigned bound = nmbrblckpxl;
	  
	  for(unsigned m = A+1; m < nmbrblckpxl; m++)
	    if( BinomialPDF(slide2, m, p) <= BinomialPDF(slide2, A, p) ){
	      bound = m;
	      break;}
	  if(bound < nmbrblckpxl)
	    C += BinomialCDF_complement(slide2, bound, p);
	}
	
	Dev = log10(C);
	
	if( A >= mean )
	  Dev *= -1;
	
       	if( fabs(MinkowskiSkyMap.call(xi,yi)) < fabs(Dev) ){
	  MinkowskiSkyMap.assign(xi,yi,Dev);
	}
      }
    }
    
    if(verbose) 
      Cout_One_Percent(thresh-thresh_min,percent,nmbr_of_thresh);
  }
  
  if(verbose) 
    std::cout << std::endl; 
  
  // Trial correction
  for(unsigned xi = 0; xi < Nx_-slide+1; xi++)
    for(unsigned yi = 0; yi < Ny_-slide+1; yi++){
      double D = MinkowskiSkyMap.call(xi,yi);
      double Dt = fabs(D) - log10(trial); 
      if(fabs(D) < 9)
	Dt = log10(1-pow(1-pow(0.1, fabs(D) ),trial));
      MinkowskiSkyMap.assign(xi,yi,copysign(Dt,D));
    }
  
  double min_D = MinkowskiSkyMap.call(0,0);
  double max_D = min_D;
  for(unsigned xi = 0; xi < Nx_-slide+1; xi++)
    for(unsigned yi = 0; yi < Ny_-slide+1; yi++){
      double D = MinkowskiSkyMap.call(xi,yi);
      if ( D < min_D )
	min_D = D;
      else if ( D > max_D )
	max_D = D;
    }
  
  // Output of the Minkowski sky map:
  // Data output
  std::stringstream filename_stst;
  filename_stst << prefix_of_ << filename
		<< "_" << (Nx_-slide+1) << "x" << (Ny_-slide+1) << "_seed_" << seed_
		<< "_lambda_" << lambda << "_slide_" << slide << ".dat";
  std::string filename_st = filename_stst.str();
  std::ofstream Data( filename_st.c_str() );
  
  Data << "# Minkowski sky map" << std::endl << "#" << std::endl
       << "# Background intensity: lambda = " << lambda << ";" << std::endl
       << "# Size of sliding window: slide = " << slide << ";" << std::endl;
  Data << "#" << std::setw(13) << "X" << std::setw(14) << "Y" << std::setw(14) << "D(A)" << std::endl;
  for(unsigned xi = 0; xi < (Nx_-slide+1); xi++)
    {
      for(unsigned yi = 0; yi < (Ny_-slide+1); yi++)
	Data << std::setw(14) << xlow_.call(xi,yi) + (slide-1)*xbinlength_/2.0
	     << std::setw(14) << ylow_.call(xi,yi) + (slide-1)*ybinlength_/2.0
	     << std::setw(14) << MinkowskiSkyMap.call(xi,yi)
	     << std::endl
	     << std::setw(14) << xlow_.call(xi,yi) + (slide-1)*xbinlength_/2.0
	     << std::setw(14) << ylow_.call(xi,yi) + ybinlength_ + (slide-1)*ybinlength_/2.0
	     << std::setw(14) << MinkowskiSkyMap.call(xi,yi)
	     << std::endl;
      
      Data << std::endl;
      
      for(unsigned yi = 0; yi < (Ny_-slide+1); yi++)
	Data << std::setw(14) << xlow_.call(xi,yi) + xbinlength_ + (slide-1)*xbinlength_/2.0
	     << std::setw(14) << ylow_.call(xi,yi) + (slide-1)*ybinlength_/2.0
	     << std::setw(14) << MinkowskiSkyMap.call(xi,yi)
	     << std::endl
	     << std::setw(14) << xlow_.call(xi,yi) + xbinlength_ + (slide-1)*xbinlength_/2.0
	     << std::setw(14) << ylow_.call(xi,yi) + ybinlength_ + (slide-1)*ybinlength_/2.0
	     << std::setw(14) << MinkowskiSkyMap.call(xi,yi)
	     << std::endl;
      
      Data << std::endl;
    }
  Data.close();
  
  // Matrix
  std::stringstream matfilename_stst;
  matfilename_stst << filename
		   << "_Matrix_" << (Nx_-slide+1) << "x" << (Ny_-slide+1) << "_seed_" << seed_
		   << "_lambda_" << lambda << "_slide_" << slide << ".dat";
  std::string matfilename = matfilename_stst.str();
  
  MinkowskiSkyMap.MatrixToFile(matfilename,prefix_of_);

  // Gnuplot
  std::stringstream gpfilename_stst;
  gpfilename_stst << prefix_of_ << filename
		  << "_" << (Nx_-slide+1) << "x" << (Ny_-slide+1) << "_seed_" << seed_
		  << "_slide_" << slide << ".gp";
      //<< "_lambda_" << lambda << "_slide_" << slide << ".gp";
  std::string gpfilename_st = gpfilename_stst.str();
  std::ofstream Gnuplot( gpfilename_st.c_str() );
  
  Gnuplot << "set terminal epslatex standalone color colortext 14" << std::endl
	  << "set output '" << filename
	  << "_" << (Nx_-slide+1) << "x" << (Ny_-slide+1) << "_seed_" << seed_
	  << "_slide_" << slide << ".tex'";
    //<< "_lambda_" << lambda << "_slide_" << slide << ".tex'";
  Gnuplot << std::endl << std::endl
	  << "xmin = " << xmin_ + (slide-1)*xbinlength_/2.0 << std::endl
	  << "xmax = " << xmax_ - (slide-1)*xbinlength_/2.0 << std::endl
	  << "ymin = " << ymin_ + (slide-1)*ybinlength_/2.0 << std::endl
	  << "ymax = " << ymax_ - (slide-1)*ybinlength_/2.0 << std::endl
	  << std::endl
	  << "#set title 'Size of sliding observation window: " << slide << "with wbc'" << std::endl
	  << std::endl
	  << "unset xtics #border mirror norotate offset 0,0.3" << std::endl
	  << "unset ytics #border mirror norotate offset 0.3,0" << std::endl
	  << "#set cbtics 10" << std::endl
	  << std::endl
	  << "set mxtics 2" << std::endl
	  << "set mytics 2" << std::endl
	  << std::endl;
  if(Ny_ == Nx_)
    Gnuplot << "set size ratio 1" << std::endl;
  else
    Gnuplot << "set size ratio " << (Ny_-slide+1)*1.0/(Nx_-slide+1) << std::endl;
  Gnuplot << "set border lw 2" << std::endl
	  << std::endl
	  << "unset xlabel #'$x$' offset 0,0.7" << std::endl
	  << "unset ylabel #'$y$' offset 0.6,0" << std::endl
	  << "# set cblabel '\\large $\\mathcal{D}(A)$' offset 0,0"
	  << std::endl
	  << "set xrange [xmin:xmax]" << std::endl
	  << "set yrange [ymin:ymax]" << std::endl
	  << "#set cbrange [" << min_D << ":" << max_D << "]" << std::endl
	  << "set cbrange [0:9]" << std::endl
	  << std::endl
	  << "set pm3d map" << std::endl
	  << "# set palette defined (" << min_D << " 'black', " << 0.5*min_D << " 'dark-blue', " << 0 << " 'white', " << 0.2*max_D << " 'purple', " << 0.4*max_D << " 'dark-red', " << 0.6*max_D << " 'red', " << 0.8*max_D << " 'orange', " << max_D << " 'yellow')" << std::endl
	  << "splot '" << filename << "_" << (Nx_-slide+1) << "x" << (Ny_-slide+1) << "_seed_" << seed_
	  << "_lambda_" << lambda << "_slide_" << slide << ".dat";
  Gnuplot << "' notitle" << std::endl
	  << "#" << std::endl
	  << std::endl;
  
  Gnuplot.close();
  
  if(verbose) 
    std::cout << "-------------------" << std::endl;
  
  return MinkowskiSkyMap;
}
// --------------------------------------

BinField<double> BinnedSkyMap::Minkowski_sky_map_sum_D_small_lambda_APC_pix_wbc(const double &lambda, const unsigned int &slide,
								   const bool &verbose, const bool &matrix_out, const bool &gnuplot_out,
								   const std::string &filename) const
{  
  // Define MinkowskiSkyMap
  unsigned Nx_mink = Nx_-slide+1;
  unsigned Ny_mink = Ny_-slide+1;
  BinField< double > MinkowskiSkyMap(Nx_mink,Ny_mink,0);
  unsigned slide2 = slide*slide;
  unsigned slide_minus_one = slide - 1;
  unsigned N_A = slide*slide + 1;
  
  // ECCDF
  double binwidth(0);
  int N_eccdf(0);
  if( slide == 15 && lambda == 100){
    binwidth = 0.020001000050;
    N_eccdf = 18792;
  }
  else if( slide ==  15 && lambda == 50){
    binwidth = 0.009000450023;
    N_eccdf = 19665;
  }
  else if( slide ==  15 && lambda == 25){
    binwidth = 0.007500375019;
    N_eccdf = 18932;
  }
  else if( slide ==  15 && lambda == 5){
    binwidth = 0.003550177509;
    N_eccdf = 19849;
  }
  else if( slide ==  10 && lambda == 100){
    binwidth = 0.013250662533;
    N_eccdf = 19606;
  }
  else if( slide ==  10 && lambda == 50){
    binwidth = 0.009500475024;
    N_eccdf = 19789;
  }
  else if( slide ==  10 && lambda == 25){
    binwidth = 0.007500375019;
    N_eccdf = 19129;
  }
  else if( slide ==  5 && lambda == 100){
    binwidth = 0.015000750038;
    N_eccdf = 16973;
  }
  else if( slide ==  5 && lambda == 50){
    binwidth = 0.011000550028;
    N_eccdf = 18488;
  }
  else if( slide ==  5 && lambda == 25){
    binwidth = 0.006000300015;
    N_eccdf = 19971;
  }
  else{
    std::cerr << "ERROR: unknown parameters w = " << slide << " and lambda = " << lambda << std::endl;
    exit(-1);
  }
  double *neg_log10_eccdf;
  neg_log10_eccdf = load_neg_log10_eccdf(slide, lambda);
  
  // DoS: ready to read in
  std::vector< stateclassifier > DoS_APC; DoS_APC.clear();
  
  // Read in DoS
  if( verbose ){
    std::cout << "Reading DoS ..." << std::endl
	      << PercentageBar << std::endl;}
  unsigned percent = 0;
  for(unsigned Ai = 0; Ai < N_A; Ai++){
    std::stringstream infile_st;
    infile_st << "parameters/structure_" << slide << "x" << slide << "/DoS_averaged_A_"
	      << Ai
	      << "_PC_pix_wbc_" << slide << "x" << slide << "_seed_0_seedwght_0.dat";
    std::string infile = infile_st.str();
    std::ifstream structure_distribution((infile).c_str ());
    if(structure_distribution.fail()){
      std::cerr << "ERROR: ifstream failed to read " << "parameters/structure_" << slide << "x" << slide << "/DoS_averaged_A_"
		<< Ai << "_PC_pix_wbc_" << slide << "x" << slide << "_seed_0_seedwght_0.dat;" << std::endl;
      exit(-1);
    }
    while(!structure_distribution.eof()){
      int P(0), C(0);
      structure_distribution >> P;
      structure_distribution >> C;
      
      double DoS = 0;
      structure_distribution >> DoS;
      if(DoS)
	DoS_APC.push_back(stateclassifier(Ai,P,C,DoS));
    }
    structure_distribution.close();
    if(verbose) 
      Cout_One_Percent(Ai,percent,N_A);
  }
  unsigned nmbr_of_macrostates = DoS_APC.size();
  
  // Determine maximum and minimum thresholds
  unsigned thresh_max = 0;
  for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
    if(countrate_.call(ni) > thresh_max)
      thresh_max = countrate_.call(ni);
  thresh_max++;
  unsigned thresh_min = thresh_max;
  for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
    if(countrate_.call(ni) < thresh_min)
      thresh_min = countrate_.call(ni);
  unsigned nmbr_of_thresh = thresh_max-thresh_min+1;
  
  // Compute Minkowski sky map
  if( verbose ){
    std::cout << "Compute Minkowski sky map ..." << std::endl
	      << PercentageBar << std::endl;
  }
  percent = 0;
  
  // Prepare excursion set and scan-window functional values
  BinField< bool > Excursion_set_old(Nx_,Ny_,true);
  BinField< bool > Excursion_set_older(Nx_,Ny_,true);
  BinField< bool > Excursion_set_new(Nx_,Ny_,true);
  BinField< int > A_Map(Nx_mink,Ny_mink,slide2);
  BinField< int > P_Map(Nx_mink,Ny_mink,4*slide);
  BinField< int > C_Map(Nx_mink,Ny_mink,1);
  
  double p(0), invp(0), Dev(0);
  std::vector< stateclassifier > prob_APC(nmbr_of_macrostates, stateclassifier(0,0,0,1));
  std::vector< stateclassifier > comp_APC(nmbr_of_macrostates, stateclassifier(0,0,0,1));
  
  double relative_difference_in_sum_D = 1;
  unsigned thresh = thresh_min;
  do{
    p = Bernoulli_p(thresh, lambda);
    invp = 1-p;
    
    // Excursion set
    for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
      Excursion_set_new.assign(ni, countrate_.call(ni) >= thresh);
    
    // determine probability distribution
    for(std::vector< stateclassifier >::size_type sr = 0; sr < nmbr_of_macrostates; sr++){
      int Ai = DoS_APC[sr].A();
      prob_APC[sr] = stateclassifier(Ai,DoS_APC[sr].P(),DoS_APC[sr].C(),DoS_APC[sr].cl()*pow(p,Ai)*pow(invp,slide2-Ai));
    }
    
    // determine compatibility distribution
    boost::sort(prob_APC);
    comp_APC=prob_APC;
    
    double min_comp = comp_APC[0].cl();
    for(std::vector< stateclassifier >::size_type sr = 1; sr < comp_APC.size(); sr++)
      comp_APC[sr].add_cl(comp_APC[sr-1].cl());
    
    // Deviation strength
    for(unsigned xi = 0; xi < Nx_mink; xi++)
      for(unsigned yi = 0; yi < Ny_mink; yi++){
	unsigned end_x = xi+slide_minus_one, end_y = yi+slide_minus_one;
	
	unsigned sx = xi;
	unsigned sy = yi;
	{
	  if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	    A_Map.add(xi,yi,-1);
	    if( !Excursion_set_old.call(sx,sy+1) )
	      P_Map.add(xi,yi,-2);
	    
	    if( !Excursion_set_old.call(sx+1,sy) )
	      P_Map.add(xi,yi,-2);
	    
	    C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx+1,sy+1)*pow_of_2[2]+ Excursion_set_old.call(sx+1,sy)*pow_of_2[4]]);
	    
	    Excursion_set_old.assign(sx,sy,false);
	  } 
	}
	for(sy = yi+1; sy < end_y; sy++)
	  {
	    if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	      A_Map.add(xi,yi,-1);
	      if( Excursion_set_old.call(sx,sy+1) && Excursion_set_old.call(sx,sy-1) )
		P_Map.add(xi,yi,2);    
	      else if( !Excursion_set_old.call(sx,sy+1) && !Excursion_set_old.call(sx,sy-1) )
		P_Map.add(xi,yi,-2);
  	      
	      if( !Excursion_set_old.call(sx+1,sy) )
		P_Map.add(xi,yi,-2);
  
	      C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx+1,sy+1)*pow_of_2[2] + Excursion_set_old.call(sx+1,sy)*pow_of_2[4] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6] + Excursion_set_old.call(sx+1,sy-1)*pow_of_2[7]]);
  
	      Excursion_set_old.assign(sx,sy,false);
	    } 
	  }
	sy = end_y;
	{
	  if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	    A_Map.add(xi,yi,-1);    
	    if( !Excursion_set_old.call(sx,sy-1) )
	      P_Map.add(xi,yi,-2);
	    
	    if( !Excursion_set_old.call(sx+1,sy) )
	      P_Map.add(xi,yi,-2);
  
	    C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx+1,sy)*pow_of_2[4] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6] + Excursion_set_old.call(sx+1,sy-1)*pow_of_2[7]]);
  
	    Excursion_set_old.assign(sx,sy,false);
	  } 
	}
	
	for(sx = xi+1; sx < end_x; sx++){
	  sy = yi;
	  {
	    if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	      A_Map.add(xi,yi,-1);
	      if( !Excursion_set_old.call(sx,sy+1) )
		P_Map.add(xi,yi,-2);
	      
	      if( Excursion_set_old.call(sx+1,sy) && Excursion_set_old.call(sx-1,sy) )
		P_Map.add(xi,yi,2);    
	      else if( !Excursion_set_old.call(sx+1,sy) && !Excursion_set_old.call(sx-1,sy) )
		P_Map.add(xi,yi,-2);
  
	      C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy+1)*pow_of_2[0] + Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx+1,sy+1)*pow_of_2[2] + Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx+1,sy)*pow_of_2[4]]);
  
	      Excursion_set_old.assign(sx,sy,false);
	    } 
	  }
	  for(sy = yi+1; sy < end_y; sy++)
	    {
	      if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
		A_Map.add(xi,yi,-1);
		if( Excursion_set_old.call(sx,sy+1) && Excursion_set_old.call(sx,sy-1) )
		  P_Map.add(xi,yi,2);    
		else if( !Excursion_set_old.call(sx,sy+1) && !Excursion_set_old.call(sx,sy-1) )
		  P_Map.add(xi,yi,-2);
  
		if( Excursion_set_old.call(sx+1,sy) && Excursion_set_old.call(sx-1,sy) )
		  P_Map.add(xi,yi,2);    
		else if( !Excursion_set_old.call(sx+1,sy) && !Excursion_set_old.call(sx-1,sy) )
		  P_Map.add(xi,yi,-2);
  
		C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy+1)*pow_of_2[0] + Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx+1,sy+1)*pow_of_2[2] + Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx+1,sy)*pow_of_2[4] + Excursion_set_old.call(sx-1,sy-1)*pow_of_2[5] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6] + Excursion_set_old.call(sx+1,sy-1)*pow_of_2[7]]);
  
		Excursion_set_old.assign(sx,sy,false);
	      } 
	    }
	  sy = end_y;
	  {
	    if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	      A_Map.add(xi,yi,-1);    
	      if( !Excursion_set_old.call(sx,sy-1) )
		P_Map.add(xi,yi,-2);
  
	      if( Excursion_set_old.call(sx+1,sy) && Excursion_set_old.call(sx-1,sy) )
		P_Map.add(xi,yi,2);    
	      else if( !Excursion_set_old.call(sx+1,sy) && !Excursion_set_old.call(sx-1,sy) )
		P_Map.add(xi,yi,-2);
  
	      C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx+1,sy)*pow_of_2[4] + Excursion_set_old.call(sx-1,sy-1)*pow_of_2[5] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6] + Excursion_set_old.call(sx+1,sy-1)*pow_of_2[7]]);
  
	      Excursion_set_old.assign(sx,sy,false);
	    } 
	  }
	}
	
	sx = end_x;
	sy = yi;
	{
	  if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	    A_Map.add(xi,yi,-1);
	    if( !Excursion_set_old.call(sx,sy+1) )
	      P_Map.add(xi,yi,-2);
	    
	    if( !Excursion_set_old.call(sx-1,sy) )
	      P_Map.add(xi,yi,-2);
  
	    C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy+1)*pow_of_2[0] + Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx-1,sy)*pow_of_2[3]]);
  
	    Excursion_set_old.assign(sx,sy,false);
	  } 
	}
	for(sy = yi+1; sy < end_y; sy++)
	  {
	    if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	      A_Map.add(xi,yi,-1);
	      if( Excursion_set_old.call(sx,sy+1) && Excursion_set_old.call(sx,sy-1) )
		P_Map.add(xi,yi,2);    
	      else if( !Excursion_set_old.call(sx,sy+1) && !Excursion_set_old.call(sx,sy-1) )
		P_Map.add(xi,yi,-2);
	      
	      if( !Excursion_set_old.call(sx-1,sy) )
		P_Map.add(xi,yi,-2);
  
	      C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy+1)*pow_of_2[0] + Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx-1,sy-1)*pow_of_2[5] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6]]);
  
	      Excursion_set_old.assign(sx,sy,false);
	    } 
	  }
	sy = end_y;
	{
	  if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	    A_Map.add(xi,yi,-1);    
	    if( !Excursion_set_old.call(sx,sy-1) )
	      P_Map.add(xi,yi,-2);
	    
	    if( !Excursion_set_old.call(sx-1,sy) )
	      P_Map.add(xi,yi,-2);
  
	    C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx-1,sy-1)*pow_of_2[5] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6]]);
  
	    Excursion_set_old.assign(sx,sy,false);
	  } 
	}
	// std::cout << A_Map.call(xi,yi) << "  " << P_Map.call(xi,yi) << "  " << C_Map.call(xi,yi) << "  " << std::endl;
	
	stateclassifier tmp_state(A_Map.call(xi,yi),P_Map.call(xi,yi),C_Map.call(xi,yi),0);
	int comp_APC_end = nmbr_of_macrostates-1;
	int tmp_sr = 0;
	for(int sr = comp_APC_end; sr >= 0; sr--)
	  if(tmp_state==comp_APC[sr]){
	    tmp_sr = sr;
	    break;
	  }
	
	Dev = comp_APC[tmp_sr].cl();
	if(Dev == 0){
	  Dev = min_comp;
	  std::cerr << "WARNING: unknown structure at (" << xi << "," << yi << "), A = " 
		    << A_Map.call(xi,yi) << ", P = " << P_Map.call(xi,yi) << ", and C = " << C_Map.call(xi,yi) << std::endl;
	}
	Dev = log10(Dev);
	
	//if( A_Map.call(xi,yi) >= p*slide2 )
	Dev *= -1;
	
	relative_difference_in_sum_D = Dev/MinkowskiSkyMap.call(xi,yi);
	
	MinkowskiSkyMap.add(xi,yi,Dev);
	
	// std::cout << thresh << "    " << Dev << std::endl;
	
	// Excursion set
	for(sx = xi; sx <= end_x; sx++)
	  for(sy = yi; sy <= end_y; sy++)
	    Excursion_set_old.assign(sx, sy, Excursion_set_older.call(sx, sy));
      }// for(xi,yi)
    
    // Excursion set
    for(unsigned ni = 0; ni < Nx_*Ny_; ni++){
      Excursion_set_old.assign(ni, Excursion_set_new.call(ni));
      Excursion_set_older.assign(ni, Excursion_set_new.call(ni));
    }
    
    if(verbose)
      Cout_One_Percent(thresh-thresh_min,percent,nmbr_of_thresh);
  
    thresh++;
  } while(thresh <= thresh_max || relative_difference_in_sum_D > 0.001);
  
  double T = 9.0;
  for(unsigned xi = 0; xi < Nx_mink; xi++)
    for(unsigned yi = 0; yi < Ny_mink; yi++){
      double S=MinkowskiSkyMap.call(xi,yi);
      int sum_i = floor(S/binwidth);
      if(sum_i < N_eccdf){
	T = neg_log10_eccdf[sum_i];
      }
      else{
	std::cerr << "WARNING: too significant result" << std::endl;
      }
      /*
      // Determine sign:
      double sum_k = 0;
      for(unsigned txi = 0; txi < slide; txi++)
	for(unsigned tyi = 0; tyi < slide; tyi++)
	  sum_k += countrate_.call(xi+txi,yi+tyi);
      if( sum_k < lambda*slide2)
	MinkowskiSkyMap.assign(xi,yi,(-1)*T);
      else
      */
      MinkowskiSkyMap.assign(xi,yi,T);
    }
  
  double min_D = MinkowskiSkyMap.call(0,0);
  double max_D = min_D;
  for(unsigned xi = 0; xi < Nx_mink; xi++)
    for(unsigned yi = 0; yi < Ny_mink; yi++){
      double D = MinkowskiSkyMap.call(xi,yi);
      if ( D < min_D )
	min_D = D;
      else if ( D > max_D )
	max_D = D;
    }
  
  // Output of the Minkowski sky map:
  // Matrix
  if(matrix_out){
    std::stringstream matfilename_stst;
    matfilename_stst << filename
		     << "_Matrix_" << (Nx_-slide+1) << "x" << (Ny_-slide+1) << "_seed_" << seed_
		     << "_lambda_" << lambda << "_slide_" << slide << ".dat";
    std::string matfilename = matfilename_stst.str();
    
    MinkowskiSkyMap.MatrixToFile(matfilename,prefix_of_);
  }
  
  if(gnuplot_out){
    // Data output
    std::stringstream filename_stst;
    filename_stst << prefix_of_ << filename
		  << "_" << (Nx_-slide+1) << "x" << (Ny_-slide+1) << "_seed_" << seed_
		  << "_lambda_" << lambda << "_slide_" << slide << ".dat";
    std::string filename_st = filename_stst.str();
    std::ofstream Data( filename_st.c_str() );
    
    Data << "# Minkowski sky map" << std::endl << "#" << std::endl
	 << "# Background intensity: lambda = " << lambda << ";" << std::endl
	 << "# Functionals: APC;" << std::endl
	 << "# Size of sliding window: slide = " << slide << ";" << std::endl
	 << "# White boundary condition;" << std::endl
	 << "# Pixelized data" << std::endl << "#" << std::endl;
    Data << "#" << std::setw(13) << "X" << std::setw(14) << "Y" << std::setw(14) << "D(A, P, chi)" << std::endl;
    for(unsigned xi = 0; xi < (Nx_-slide+1); xi++)
      {
	for(unsigned yi = 0; yi < (Ny_-slide+1); yi++)
	  Data << std::setw(14) << xlow_.call(xi,yi) + (slide-1)*xbinlength_/2.0
	       << std::setw(14) << ylow_.call(xi,yi) + (slide-1)*ybinlength_/2.0
	       << std::setw(14) << MinkowskiSkyMap.call(xi,yi)
	       << std::endl
	       << std::setw(14) << xlow_.call(xi,yi) + (slide-1)*xbinlength_/2.0
	       << std::setw(14) << ylow_.call(xi,yi) + ybinlength_ + (slide-1)*ybinlength_/2.0
	       << std::setw(14) << MinkowskiSkyMap.call(xi,yi)
	       << std::endl;
      
	Data << std::endl;
      
	for(unsigned yi = 0; yi < (Ny_-slide+1); yi++)
	  Data << std::setw(14) << xlow_.call(xi,yi) + xbinlength_ + (slide-1)*xbinlength_/2.0
	       << std::setw(14) << ylow_.call(xi,yi) + (slide-1)*ybinlength_/2.0
	       << std::setw(14) << MinkowskiSkyMap.call(xi,yi)
	       << std::endl
	       << std::setw(14) << xlow_.call(xi,yi) + xbinlength_ + (slide-1)*xbinlength_/2.0
	       << std::setw(14) << ylow_.call(xi,yi) + ybinlength_ + (slide-1)*ybinlength_/2.0
	       << std::setw(14) << MinkowskiSkyMap.call(xi,yi)
	       << std::endl;
      
	Data << std::endl;
      }
    Data.close();

    // Gnuplot
    std::stringstream gpfilename_stst;
    gpfilename_stst << prefix_of_ << filename
		    << "_" << (Nx_-slide+1) << "x" << (Ny_-slide+1) << "_seed_" << seed_
		    << "_slide_" << slide << ".gp";
      //<< "_lambda_" << lambda << "_slide_" << slide << ".gp";
    std::string gpfilename_st = gpfilename_stst.str();
    std::ofstream Gnuplot( gpfilename_st.c_str() );
  
    Gnuplot << "set terminal epslatex standalone color colortext 14" << std::endl
	    << "set output '" << filename
	    << "_" << (Nx_-slide+1) << "x" << (Ny_-slide+1) << "_seed_" << seed_
	    << "_slide_" << slide << ".tex'";
    //<< "_lambda_" << lambda << "_slide_" << slide << ".tex'";
    Gnuplot << std::endl << std::endl
	    << "xmin = " << xmin_ + (slide-1)*xbinlength_/2.0 << std::endl
	    << "xmax = " << xmax_ - (slide-1)*xbinlength_/2.0 << std::endl
	    << "ymin = " << ymin_ + (slide-1)*ybinlength_/2.0 << std::endl
	    << "ymax = " << ymax_ - (slide-1)*ybinlength_/2.0 << std::endl
	    << std::endl
	    << "#set title 'Size of sliding observation window: " << slide << "with wbc'" << std::endl
	    << "unset xtics #border mirror norotate offset 0,0.3" << std::endl
	    << "unset ytics #border mirror norotate offset 0.3,0" << std::endl
	    << "#set cbtics 10" << std::endl
	    << std::endl
	    << "set mxtics 2" << std::endl
	    << "set mytics 2" << std::endl
	    << std::endl;
    if(Ny_ == Nx_)
      Gnuplot << "set size ratio 1" << std::endl;
    else
      Gnuplot << "set size ratio " << (Ny_-slide+1)*1.0/(Nx_-slide+1) << std::endl;
    Gnuplot << "set border lw 2" << std::endl
	    << std::endl
	    << "unset xlabel #'$x$' offset 0,0.7" << std::endl
	    << "unset ylabel #'$y$' offset 0.6,0" << std::endl
	    << "# set cblabel '\\large $\\mathcal{D}(A,P,\\chi)$' offset 0,0" << std::endl
	    << std::endl
	    << "set xrange [xmin:xmax]" << std::endl
	    << "set yrange [ymin:ymax]" << std::endl
	    << "#set cbrange [" << min_D << ":" << max_D << "]" << std::endl
	    << "set cbrange [0:9]" << std::endl
	    << std::endl
	    << "set pm3d map" << std::endl
	    << "# set palette defined (" << min_D << " 'black', " << 0.5*min_D << " 'dark-blue', " << 0 << " 'white', " << 0.2*max_D << " 'purple', " << 0.4*max_D << " 'dark-red', " << 0.6*max_D << " 'red', " << 0.8*max_D << " 'orange', " << max_D << " 'yellow')" << std::endl
	    << "splot '" << filename << "_" << (Nx_-slide+1) << "x" << (Ny_-slide+1) << "_seed_" << seed_
	    << "_lambda_" << lambda << "_slide_" << slide << ".dat";
    Gnuplot << "' notitle" << std::endl
	    << "#" << std::endl
	    << std::endl;
    
    Gnuplot.close();
  }
  
  if( verbose )
    std::cout << "-------------------" << std::endl;
  
  return MinkowskiSkyMap;
}

BinField<double> BinnedSkyMap::Minkowski_sky_map_sum_D_small_lambda_APC_pix_wbc_chessboard(const double &lambda, const unsigned int &slide,
									      const bool &verbose, const bool &matrix_out,
									      const std::string &filename) const
{
  // Define MinkowskiSkyMap
  unsigned Nx_mink = Nx_/slide;
  unsigned Ny_mink = Ny_/slide;
  BinField< double > MinkowskiSkyMap(Nx_mink,Ny_mink,0);
  //BinField< double > MinkowskiSkyMap_sum(Nx_mink,Ny_mink,0);
  //BinField< double > MinkowskiSkyMap_single(Nx_mink,Ny_mink,0);
  unsigned slide2 = slide*slide;
  unsigned slide_minus_one = slide - 1;
  unsigned N_A = slide*slide + 1;
  
  // ECCDF
  double binwidth(0);
  int N_eccdf(0);
  if( slide == 15 && lambda == 100){
    binwidth = 0.020001000050;
    N_eccdf = 18792;
  }
  else if( slide ==  15 && lambda == 50){
    binwidth = 0.009000450023;
    N_eccdf = 19665;
  }
  else if( slide ==  15 && lambda == 25){
    binwidth = 0.007500375019;
    N_eccdf = 18932;
  }
  else if( slide ==  15 && lambda == 5){
    binwidth = 0.003550177509;
    N_eccdf = 19849;
  }
  else if( slide ==  10 && lambda == 100){
    binwidth = 0.013250662533;
    N_eccdf = 19606;
  }
  else if( slide ==  10 && lambda == 50){
    binwidth = 0.009500475024;
    N_eccdf = 19789;
  }
  else if( slide ==  10 && lambda == 25){
    binwidth = 0.007500375019;
    N_eccdf = 19129;
  }
  else if( slide ==  5 && lambda == 100){
    binwidth = 0.015000750038;
    N_eccdf = 16973;
  }
  else if( slide ==  5 && lambda == 50){
    binwidth = 0.011000550028;
    N_eccdf = 18488;
  }
  else if( slide ==  5 && lambda == 25){
    binwidth = 0.006000300015;
    N_eccdf = 19971;
  }
  else{
    std::cerr << "ERROR: unknown parameters w = " << slide << " and lambda = " << lambda << std::endl;
    exit(-1);
  }
  double *neg_log10_eccdf;
  neg_log10_eccdf = load_neg_log10_eccdf(slide, lambda);
  
  // DoS: ready to read in
  std::vector< stateclassifier > DoS_APC; DoS_APC.clear();
  
  // Read in DoS
  if( verbose ){
    std::cout << "Reading DoS ..." << std::endl
	      << PercentageBar << std::endl;}
  unsigned percent = 0;
  for(unsigned Ai = 0; Ai < N_A; Ai++){
    std::stringstream infile_st;
    infile_st << "parameters/structure_" << slide << "x" << slide << "/DoS_averaged_A_"
	      << Ai
	      << "_PC_pix_wbc_" << slide << "x" << slide << "_seed_0_seedwght_0.dat";
    std::string infile = infile_st.str();
    std::ifstream structure_distribution((infile).c_str ());
    if(structure_distribution.fail()){
      std::cerr << "ERROR: ifstream failed to read " << "parameters/structure_" << slide << "x" << slide << "/DoS_averaged_A_"
		<< Ai << "_PC_pix_wbc_" << slide << "x" << slide << "_seed_0_seedwght_0.dat;" << std::endl;
      exit(-1);
    }
    while(!structure_distribution.eof()){
      int P(0), C(0);
      structure_distribution >> P;
      structure_distribution >> C;
      
      double DoS = 0;
      structure_distribution >> DoS;
      if(DoS)
	DoS_APC.push_back(stateclassifier(Ai,P,C,DoS));
    }
    structure_distribution.close();
    if(verbose) 
      Cout_One_Percent(Ai,percent,N_A);
  }
  unsigned nmbr_of_macrostates = DoS_APC.size();
  
  // Determine maximum and minimum thresholds
  unsigned thresh_max = 0;
  for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
    if(countrate_.call(ni) > thresh_max)
      thresh_max = countrate_.call(ni);
  thresh_max++;
  unsigned thresh_min = thresh_max;
  for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
    if(countrate_.call(ni) < thresh_min)
      thresh_min = countrate_.call(ni);
  unsigned nmbr_of_thresh = thresh_max-thresh_min+1;
  
  // Compute Minkowski sky map
  if( verbose ){
    std::cout << "Compute Minkowski sky map ..." << std::endl
	      << PercentageBar << std::endl;
  }
  percent = 0;
  
  // Prepare excursion set and scan-window functional values
  BinField< bool > Excursion_set_old(Nx_,Ny_,true);
  BinField< bool > Excursion_set_new(Nx_,Ny_,true);
  BinField< int > A_Map(Nx_mink,Ny_mink,slide2);
  BinField< int > P_Map(Nx_mink,Ny_mink,4*slide);
  BinField< int > C_Map(Nx_mink,Ny_mink,1);
  
  double p(0), invp(0), Dev(0);
  std::vector< stateclassifier > prob_APC(nmbr_of_macrostates, stateclassifier(0,0,0,1));
  std::vector< stateclassifier > comp_APC(nmbr_of_macrostates, stateclassifier(0,0,0,1));
  
  double relative_difference_in_sum_D = 1;
  unsigned thresh = thresh_min;
  do{
    p = Bernoulli_p(thresh, lambda);
    invp = 1-p;
    
    // Excursion set
    for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
      Excursion_set_new.assign(ni, countrate_.call(ni) >= thresh);
    
    // determine probability distribution
    for(std::vector< stateclassifier >::size_type sr = 0; sr < nmbr_of_macrostates; sr++){
      int Ai = DoS_APC[sr].A();
      prob_APC[sr] = stateclassifier(Ai,DoS_APC[sr].P(),DoS_APC[sr].C(),DoS_APC[sr].cl()*pow(p,Ai)*pow(invp,slide2-Ai));
    }
    
    // determine compatibility distribution
    boost::sort(prob_APC);
    comp_APC=prob_APC;
    
    double min_comp = comp_APC[0].cl();
    for(std::vector< stateclassifier >::size_type sr = 1; sr < comp_APC.size(); sr++)
      comp_APC[sr].add_cl(comp_APC[sr-1].cl());
    
    // Deviation strength
    for(unsigned xi = 0; xi < Nx_mink; xi++)
      for(unsigned yi = 0; yi < Ny_mink; yi++){
	unsigned end_x = xi*slide+slide_minus_one, end_y = yi*slide+slide_minus_one;
	
	unsigned sx = xi*slide;
	unsigned sy = yi*slide;
	{
	  if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	    A_Map.add(xi,yi,-1);
	    if( !Excursion_set_old.call(sx,sy+1) )
	      P_Map.add(xi,yi,-2);
	    
	    if( !Excursion_set_old.call(sx+1,sy) )
	      P_Map.add(xi,yi,-2);
	    
	    C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx+1,sy+1)*pow_of_2[2]+ Excursion_set_old.call(sx+1,sy)*pow_of_2[4]]);
	    
	    Excursion_set_old.assign(sx,sy,false);
	  } 
	}
	for(sy = yi*slide+1; sy < end_y; sy++)
	  {
	    if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	      A_Map.add(xi,yi,-1);
	      if( Excursion_set_old.call(sx,sy+1) && Excursion_set_old.call(sx,sy-1) )
		P_Map.add(xi,yi,2);    
	      else if( !Excursion_set_old.call(sx,sy+1) && !Excursion_set_old.call(sx,sy-1) )
		P_Map.add(xi,yi,-2);
  	      
	      if( !Excursion_set_old.call(sx+1,sy) )
		P_Map.add(xi,yi,-2);
  
	      C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx+1,sy+1)*pow_of_2[2] + Excursion_set_old.call(sx+1,sy)*pow_of_2[4] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6] + Excursion_set_old.call(sx+1,sy-1)*pow_of_2[7]]);
  
	      Excursion_set_old.assign(sx,sy,false);
	    } 
	  }
	sy = end_y;
	{
	  if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	    A_Map.add(xi,yi,-1);    
	    if( !Excursion_set_old.call(sx,sy-1) )
	      P_Map.add(xi,yi,-2);
	    
	    if( !Excursion_set_old.call(sx+1,sy) )
	      P_Map.add(xi,yi,-2);
  
	    C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx+1,sy)*pow_of_2[4] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6] + Excursion_set_old.call(sx+1,sy-1)*pow_of_2[7]]);
  
	    Excursion_set_old.assign(sx,sy,false);
	  } 
	}
	
	for(sx = xi*slide+1; sx < end_x; sx++){
	  sy = yi*slide;
	  {
	    if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	      A_Map.add(xi,yi,-1);
	      if( !Excursion_set_old.call(sx,sy+1) )
		P_Map.add(xi,yi,-2);
	      
	      if( Excursion_set_old.call(sx+1,sy) && Excursion_set_old.call(sx-1,sy) )
		P_Map.add(xi,yi,2);    
	      else if( !Excursion_set_old.call(sx+1,sy) && !Excursion_set_old.call(sx-1,sy) )
		P_Map.add(xi,yi,-2);
  
	      C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy+1)*pow_of_2[0] + Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx+1,sy+1)*pow_of_2[2] + Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx+1,sy)*pow_of_2[4]]);
  
	      Excursion_set_old.assign(sx,sy,false);
	    } 
	  }
	  for(sy = yi*slide+1; sy < end_y; sy++)
	    {
	      if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
		A_Map.add(xi,yi,-1);
		if( Excursion_set_old.call(sx,sy+1) && Excursion_set_old.call(sx,sy-1) )
		  P_Map.add(xi,yi,2);    
		else if( !Excursion_set_old.call(sx,sy+1) && !Excursion_set_old.call(sx,sy-1) )
		  P_Map.add(xi,yi,-2);
  
		if( Excursion_set_old.call(sx+1,sy) && Excursion_set_old.call(sx-1,sy) )
		  P_Map.add(xi,yi,2);    
		else if( !Excursion_set_old.call(sx+1,sy) && !Excursion_set_old.call(sx-1,sy) )
		  P_Map.add(xi,yi,-2);
  
		C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy+1)*pow_of_2[0] + Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx+1,sy+1)*pow_of_2[2] + Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx+1,sy)*pow_of_2[4] + Excursion_set_old.call(sx-1,sy-1)*pow_of_2[5] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6] + Excursion_set_old.call(sx+1,sy-1)*pow_of_2[7]]);
  
		Excursion_set_old.assign(sx,sy,false);
	      } 
	    }
	  sy = end_y;
	  {
	    if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	      A_Map.add(xi,yi,-1);    
	      if( !Excursion_set_old.call(sx,sy-1) )
		P_Map.add(xi,yi,-2);
  
	      if( Excursion_set_old.call(sx+1,sy) && Excursion_set_old.call(sx-1,sy) )
		P_Map.add(xi,yi,2);    
	      else if( !Excursion_set_old.call(sx+1,sy) && !Excursion_set_old.call(sx-1,sy) )
		P_Map.add(xi,yi,-2);
  
	      C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx+1,sy)*pow_of_2[4] + Excursion_set_old.call(sx-1,sy-1)*pow_of_2[5] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6] + Excursion_set_old.call(sx+1,sy-1)*pow_of_2[7]]);
  
	      Excursion_set_old.assign(sx,sy,false);
	    } 
	  }
	}
	
	sx = end_x;
	sy = yi*slide;
	{
	  if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	    A_Map.add(xi,yi,-1);
	    if( !Excursion_set_old.call(sx,sy+1) )
	      P_Map.add(xi,yi,-2);
	    
	    if( !Excursion_set_old.call(sx-1,sy) )
	      P_Map.add(xi,yi,-2);
  
	    C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy+1)*pow_of_2[0] + Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx-1,sy)*pow_of_2[3]]);
  
	    Excursion_set_old.assign(sx,sy,false);
	  } 
	}
	for(sy = yi*slide+1; sy < end_y; sy++)
	  {
	    if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	      A_Map.add(xi,yi,-1);
	      if( Excursion_set_old.call(sx,sy+1) && Excursion_set_old.call(sx,sy-1) )
		P_Map.add(xi,yi,2);    
	      else if( !Excursion_set_old.call(sx,sy+1) && !Excursion_set_old.call(sx,sy-1) )
		P_Map.add(xi,yi,-2);
	      
	      if( !Excursion_set_old.call(sx-1,sy) )
		P_Map.add(xi,yi,-2);
  
	      C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy+1)*pow_of_2[0] + Excursion_set_old.call(sx,sy+1)*pow_of_2[1] + Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx-1,sy-1)*pow_of_2[5] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6]]);
  
	      Excursion_set_old.assign(sx,sy,false);
	    } 
	  }
	sy = end_y;
	{
	  if(Excursion_set_old.call(sx,sy) && !Excursion_set_new.call(sx,sy)){
	    A_Map.add(xi,yi,-1);    
	    if( !Excursion_set_old.call(sx,sy-1) )
	      P_Map.add(xi,yi,-2);
	    
	    if( !Excursion_set_old.call(sx-1,sy) )
	      P_Map.add(xi,yi,-2);
  
	    C_Map.add(xi,yi,change_in_C_b2w[Excursion_set_old.call(sx-1,sy)*pow_of_2[3] + Excursion_set_old.call(sx-1,sy-1)*pow_of_2[5] + Excursion_set_old.call(sx,sy-1)*pow_of_2[6]]);
  
	    Excursion_set_old.assign(sx,sy,false);
	  } 
	}
	// std::cout << A_Map.call(xi,yi) << "  " << P_Map.call(xi,yi) << "  " << C_Map.call(xi,yi) << "  " << std::endl;
	
	stateclassifier tmp_state(A_Map.call(xi,yi),P_Map.call(xi,yi),C_Map.call(xi,yi),0);
	int comp_APC_end = nmbr_of_macrostates-1;
	int tmp_sr = 0;
	for(int sr = comp_APC_end; sr >= 0; sr--)
	  if(tmp_state==comp_APC[sr]){
	    tmp_sr = sr;
	    break;
	  }
	
	Dev = comp_APC[tmp_sr].cl();
	if(Dev == 0){
	  Dev = min_comp;
	  std::cerr << "WARNING: unknown structure A = " << A_Map.call(xi,yi) << ", P = " << P_Map.call(xi,yi) << ", and C = " << C_Map.call(xi,yi) << std::endl;
	}
	Dev = log10(Dev);
	
	//if( A_Map.call(xi,yi) >= p*slide2 )
	Dev *= -1;
	
	relative_difference_in_sum_D = Dev/MinkowskiSkyMap.call(xi,yi);
	
	MinkowskiSkyMap.add(xi,yi,Dev);
	
	/*
	MinkowskiSkyMap_sum.add(xi,yi,Dev);
	if(thresh==lambda)
	  MinkowskiSkyMap_single.assign(xi,yi,Dev);
	*/
	/*
       	if( fabs(MinkowskiSkyMap.call(xi,yi)) < fabs(Dev) )
	  MinkowskiSkyMap.assign(xi,yi,Dev);
	*/
	// std::cout << thresh << "    " << Dev << std::endl;
      }// for(xi,yi)
    
    // Excursion set
    for(unsigned ni = 0; ni < Nx_*Ny_; ni++)
      if(Excursion_set_old.call(ni) != Excursion_set_new.call(ni)){
	std::cerr << "ERROR: Excursion sets do not match!" << std::endl;
	exit(-1);
      }
    
    if(verbose)
      Cout_One_Percent(thresh-thresh_min,percent,nmbr_of_thresh);
    
    thresh++;
  } while(thresh <= thresh_max || relative_difference_in_sum_D > 0.001);
  
  double T = 9.0;
  for(unsigned xi = 0; xi < Nx_mink; xi++)
    for(unsigned yi = 0; yi < Ny_mink; yi++){
      double S=MinkowskiSkyMap.call(xi,yi);
      int sum_i = floor(S/binwidth);
      if(sum_i < N_eccdf){
	T = neg_log10_eccdf[sum_i];
      }
      else{
	std::cerr << "WARNING: too significant result" << std::endl;
      }
      
      /*
      // Determine sign:
      double sum_k = 0;
      for(unsigned txi = 0; txi < slide; txi++)
	for(unsigned tyi = 0; tyi < slide; tyi++)
	  sum_k += countrate_.call(xi+txi,yi+tyi);
      if( sum_k < lambda*slide2)
	MinkowskiSkyMap.assign(xi,yi,(-1)*T);
      else
      */
      MinkowskiSkyMap.assign(xi,yi,T);
    }
  
  // Output of the Minkowski sky map:
  // Matrix
  if(matrix_out){
    std::stringstream matfilename_stst;
    matfilename_stst << filename
		     << "_Matrix_" << Nx_mink << "x" << Ny_mink << "_seed_" << seed_
		     << "_lambda_" << lambda << "_slide_" << slide << ".dat";
    std::string matfilename = matfilename_stst.str();
    
    MinkowskiSkyMap.MatrixToFile(matfilename,prefix_of_);
  }
  /*
  if(matrix_out){
    std::stringstream matfilename_stst;
    matfilename_stst << filename
		     << "_single_Matrix_" << Nx_mink << "x" << Ny_mink << "_seed_" << seed_
		     << "_lambda_" << lambda << "_slide_" << slide << ".dat";
    std::string matfilename = matfilename_stst.str();
    
    MinkowskiSkyMap_single.MatrixToFile(matfilename,prefix_of_);
  }
  if(matrix_out){
    std::stringstream matfilename_stst;
    matfilename_stst << filename
		     << "_sum_Matrix_" << Nx_mink << "x" << Ny_mink << "_seed_" << seed_
		     << "_lambda_" << lambda << "_slide_" << slide << ".dat";
    std::string matfilename = matfilename_stst.str();
    
    MinkowskiSkyMap_sum.MatrixToFile(matfilename,prefix_of_);
  }
  */
  if( verbose )
    std::cout << "-------------------" << std::endl;
  
  return MinkowskiSkyMap;
}


double * BinnedSkyMap::load_neg_log10_eccdf(const unsigned &slide, const double &lambda) const{
  
  if( slide == 15 && lambda == 100){
    static double neg_log10_eccdf[] = {0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000434, 0.000000000434, 0.000000000434, 0.000000000434, 0.000000000434, 0.000000000869, 0.000000001303, 0.000000001303, 0.000000001303, 0.000000001303, 0.000000001303, 0.000000001737, 0.000000001737, 0.000000001737, 0.000000003040, 0.000000003909, 0.000000003909, 0.000000003909, 0.000000003909, 0.000000003909, 0.000000005212, 0.000000005212, 0.000000005646, 0.000000006514, 0.000000006949, 0.000000007817, 0.000000009120, 0.000000009554, 0.000000009989, 0.000000011726, 0.000000013029, 0.000000014332, 0.000000015635, 0.000000016503, 0.000000018675, 0.000000018675, 0.000000019543, 0.000000022583, 0.000000025189, 0.000000026926, 0.000000031269, 0.000000034744, 0.000000038652, 0.000000045167, 0.000000049075, 0.000000053418, 0.000000057327, 0.000000063841, 0.000000070790, 0.000000079910, 0.000000088596, 0.000000096848, 0.000000105968, 0.000000115088, 0.000000123340, 0.000000136803, 0.000000151134, 0.000000165032, 0.000000175889, 0.000000191958, 0.000000206290, 0.000000223227, 0.000000241902, 0.000000264485, 0.000000279686, 0.000000304875, 0.000000330064, 0.000000353950, 0.000000386522, 0.000000419963, 0.000000449061, 0.000000482067, 0.000000521588, 0.000000560675, 0.000000605841, 0.000000651008, 0.000000699649, 0.000000746553, 0.000000804748, 0.000000865115, 0.000000927220, 0.000001003221, 0.000001072709, 0.000001148276, 0.000001236004, 0.000001315046, 0.000001407116, 0.000001503964, 0.000001600378, 0.000001713729, 0.000001827081, 0.000001932180, 0.000002056389, 0.000002198839, 0.000002336511, 0.000002484606, 0.000002643993, 0.000002798603, 0.000002971453, 0.000003163413, 0.000003343646, 0.000003549938, 0.000003759269, 0.000003989013, 0.000004216586, 0.000004465439, 0.000004726887, 0.000004993981, 0.000005279316, 0.000005583760, 0.000005891245, 0.000006196558, 0.000006557896, 0.000006917063, 0.000007296643, 0.000007696635, 0.000008109222, 0.000008516598, 0.000008968708, 0.000009472935, 0.000009969779, 0.000010476612, 0.000010982578, 0.000011541964, 0.000012125237, 0.000012738478, 0.000013377345, 0.000014028373, 0.000014684613, 0.000015396881, 0.000016129128, 0.000016917837, 0.000017699598, 0.000018541296, 0.000019416872, 0.000020296793, 0.000021250549, 0.000022226023, 0.000023216266, 0.000024285558, 0.000025405233, 0.000026506235, 0.000027697579, 0.000028892836, 0.000030179738, 0.000031505733, 0.000032913820, 0.000034340588, 0.000035786905, 0.000037334861, 0.000038931468, 0.000040504192, 0.000042137729, 0.000043884200, 0.000045671072, 0.000047492698, 0.000049393817, 0.000051351844, 0.000053384587, 0.000055469896, 0.000057678137, 0.000059923309, 0.000062206716, 0.000064602634, 0.000067050689, 0.000069566519, 0.000072157508, 0.000074854933, 0.000077701365, 0.000080536521, 0.000083521991, 0.000086510088, 0.000089612014, 0.000092813002, 0.000096144765, 0.000099543451, 0.000103085514, 0.000106800064, 0.000110552438, 0.000114376957, 0.000118288827, 0.000122311073, 0.000126457602, 0.000130739710, 0.000135107878, 0.000139600338, 0.000144261406, 0.000149061112, 0.000153929080, 0.000158920052, 0.000164022303, 0.000169328812, 0.000174755299, 0.000180281784, 0.000186007330, 0.000191893278, 0.000197909653, 0.000204030391, 0.000210331101, 0.000216775725, 0.000223388603, 0.000230148016, 0.000237089603, 0.000244183390, 0.000251461975, 0.000259054430, 0.000266705248, 0.000274514819, 0.000282554856, 0.000290725417, 0.000299091699, 0.000307741938, 0.000316487963, 0.000325473641, 0.000334641613, 0.000344026664, 0.000353602728, 0.000363388509, 0.000373438353, 0.000383599708, 0.000393967343, 0.000404542576, 0.000415374546, 0.000426439362, 0.000437826162, 0.000449312817, 0.000461084975, 0.000473120487, 0.000485474592, 0.000498053401, 0.000510817801, 0.000523839989, 0.000537341753, 0.000550969606, 0.000564793139, 0.000579045463, 0.000593517873, 0.000608263450, 0.000623340500, 0.000638756886, 0.000654326052, 0.000670429914, 0.000686560471, 0.000703057896, 0.000719773890, 0.000737022137, 0.000754543844, 0.000772205045, 0.000790309446, 0.000808772257, 0.000827409902, 0.000846414301, 0.000865916050, 0.000885753765, 0.000906064579, 0.000926445543, 0.000947152589, 0.000968300675, 0.000989682233, 0.001011514515, 0.001033696159, 0.001056402229, 0.001079468234, 0.001102765791, 0.001126574913, 0.001150886103, 0.001175699012, 0.001200774636, 0.001226286359, 0.001252175459, 0.001278508647, 0.001305192783, 0.001332250150, 0.001359701730, 0.001387696600, 0.001415904160, 0.001444739276, 0.001474062191, 0.001503785844, 0.001533895063, 0.001564336752, 0.001595355518, 0.001626695159, 0.001658647837, 0.001691241153, 0.001724148235, 0.001757660875, 0.001791495755, 0.001825911572, 0.001860672065, 0.001896025947, 0.001931804979, 0.001968154541, 0.002005269369, 0.002042578051, 0.002080343773, 0.002118578001, 0.002157681962, 0.002197321972, 0.002237322250, 0.002277944845, 0.002319142774, 0.002360821447, 0.002403153935, 0.002445844751, 0.002488981797, 0.002533020822, 0.002577254356, 0.002622265314, 0.002667848563, 0.002713999037, 0.002760915781, 0.002808408898, 0.002856602295, 0.002905110147, 0.002954495312, 0.003004350952, 0.003054899846, 0.003105955460, 0.003157513597, 0.003209596313, 0.003262388885, 0.003315961593, 0.003369850348, 0.003424393197, 0.003479631536, 0.003535867306, 0.003592607815, 0.003649706543, 0.003707774229, 0.003766583896, 0.003825974612, 0.003885752929, 0.003946533008, 0.004007995293, 0.004069705197, 0.004132096008, 0.004195566624, 0.004259984191, 0.004324756086, 0.004390178154, 0.004456361697, 0.004523054290, 0.004591129157, 0.004659471938, 0.004728831336, 0.004798809993, 0.004869239558, 0.004940393030, 0.005012630580, 0.005085056790, 0.005158797270, 0.005233210490, 0.005308421651, 0.005384162061, 0.005460839240, 0.005537995778, 0.005616094315, 0.005694983985, 0.005774770843, 0.005855225603, 0.005936853624, 0.006018563484, 0.006101145373, 0.006184594498, 0.006269126816, 0.006354463993, 0.006440511294, 0.006527314987, 0.006614933729, 0.006703539578, 0.006792800075, 0.006882819307, 0.006973887253, 0.007065507973, 0.007158083589, 0.007251266278, 0.007345303782, 0.007440313728, 0.007535836735, 0.007632299221, 0.007729911358, 0.007828337850, 0.007927320503, 0.008027225112, 0.008127779286, 0.008229846545, 0.008332565557, 0.008436546531, 0.008541158322, 0.008646425745, 0.008752580145, 0.008859456835, 0.008967478861, 0.009076679000, 0.009186339535, 0.009296598354, 0.009408184633, 0.009520342503, 0.009633475141, 0.009747432733, 0.009862110184, 0.009978356839, 0.010094776416, 0.010212862994, 0.010331498488, 0.010450670438, 0.010571275496, 0.010692765125, 0.010815564821, 0.010938730607, 0.011062927468, 0.011187852807, 0.011314426393, 0.011441405182, 0.011569211008, 0.011698375063, 0.011828131831, 0.011958859947, 0.012090746046, 0.012224069905, 0.012357647396, 0.012491981134, 0.012627536350, 0.012764063717, 0.012901999898, 0.013040302125, 0.013179729098, 0.013319873474, 0.013460845684, 0.013603376944, 0.013746734970, 0.013890915634, 0.014036481785, 0.014183173201, 0.014330362596, 0.014478891469, 0.014627821510, 0.014777851817, 0.014929309837, 0.015081457387, 0.015234937487, 0.015388904706, 0.015543982060, 0.015700129292, 0.015858140771, 0.016016762209, 0.016176321516, 0.016336464363, 0.016498199604, 0.016660399444, 0.016823410643, 0.016987472094, 0.017152899410, 0.017319301826, 0.017486942226, 0.017655411770, 0.017825394207, 0.017995359458, 0.018166974854, 0.018339808780, 0.018513650108, 0.018688161793, 0.018863911581, 0.019040449004, 0.019217823112, 0.019396637079, 0.019575984923, 0.019756911281, 0.019938696487, 0.020121761907, 0.020305818331, 0.020491042723, 0.020677130923, 0.020864494110, 0.021052117774, 0.021241154690, 0.021430928662, 0.021621879738, 0.021814562000, 0.022007240697, 0.022201653577, 0.022397065338, 0.022593422393, 0.022790752988, 0.022989084053, 0.023188769207, 0.023389120040, 0.023591840898, 0.023794927092, 0.023999322270, 0.024204470784, 0.024410780229, 0.024617371576, 0.024825727714, 0.025034663708, 0.025244925552, 0.025456310195, 0.025669082310, 0.025882764943, 0.026097060537, 0.026312881252, 0.026529578843, 0.026747436780, 0.026966533930, 0.027186887832, 0.027407771309, 0.027629784005, 0.027852740518, 0.028076788646, 0.028301826188, 0.028528321649, 0.028756148715, 0.028984372993, 0.029214348172, 0.029444956897, 0.029676963657, 0.029910045413, 0.030143767646, 0.030378830897, 0.030614932037, 0.030851357434, 0.031089430141, 0.031328414816, 0.031569257714, 0.031810860897, 0.032053364892, 0.032296235009, 0.032540999617, 0.032786554193, 0.033033455779, 0.033281293507, 0.033531192053, 0.033781257949, 0.034033110441, 0.034285252410, 0.034538679311, 0.034793520349, 0.035049187549, 0.035306394547, 0.035563799686, 0.035822409696, 0.036082394172, 0.036343180821, 0.036605243953, 0.036868917202, 0.037132897754, 0.037398850843, 0.037665591630, 0.037932773125, 0.038201315565, 0.038471276418, 0.038742585105, 0.039014324226, 0.039287291750, 0.039561651156, 0.039836521489, 0.040112668099, 0.040389804479, 0.040668099813, 0.040947442159, 0.041227926544, 0.041509732348, 0.041792416380, 0.042076296987, 0.042360871272, 0.042646526637, 0.042933566650, 0.043221554971, 0.043509772793, 0.043800251968, 0.044090199703, 0.044382550002, 0.044675613719, 0.044969983182, 0.045265546197, 0.045562973578, 0.045860083034, 0.046159196009, 0.046458747617, 0.046760292636, 0.047062423725, 0.047365856365, 0.047669610208, 0.047974905117, 0.048281042361, 0.048588144512, 0.048896654799, 0.049207679219, 0.049517947208, 0.049829168633, 0.050142212551, 0.050456282265, 0.050770905845, 0.051086982693, 0.051403484315, 0.051722069444, 0.052041043647, 0.052360954874, 0.052683042557, 0.053005209329, 0.053328304681, 0.053652668592, 0.053977920145, 0.054304743753, 0.054632563097, 0.054961324550, 0.055291077913, 0.055622665425, 0.055954412932, 0.056286849111, 0.056620404813, 0.056955495875, 0.057291706431, 0.057628485931, 0.057966386894, 0.058305954623, 0.058645974966, 0.058987356173, 0.059329725355, 0.059672715925, 0.060017272102, 0.060362857923, 0.060709594497, 0.061057002507, 0.061405530570, 0.061755435492, 0.062105530361, 0.062457730002, 0.062809942982, 0.063163371561, 0.063518691080, 0.063875288636, 0.064232275909, 0.064589631952, 0.064948567368, 0.065309389202, 0.065669878566, 0.066032011822, 0.066395192174, 0.066759486823, 0.067124130899, 0.067490838214, 0.067858168705, 0.068226550557, 0.068595990252, 0.068965996458, 0.069338234298, 0.069711485450, 0.070085154633, 0.070459380053, 0.070835196074, 0.071212023236, 0.071589881814, 0.071968666065, 0.072348365814, 0.072729545423, 0.073111354363, 0.073494161262, 0.073878002959, 0.074263166339, 0.074648769603, 0.075035419125, 0.075422427136, 0.075810535325, 0.076200883995, 0.076591410800, 0.076983498803, 0.077375906856, 0.077769587160, 0.078164165383, 0.078559883559, 0.078956712704, 0.079353913119, 0.079752231719, 0.080152364891, 0.080552885883, 0.080955102595, 0.081357599977, 0.081761018059, 0.082165089342, 0.082570536450, 0.082977546367, 0.083385111735, 0.083793617211, 0.084203330805, 0.084614077737, 0.085025759686, 0.085441207057, 0.085854222993, 0.086268577116, 0.086683473821, 0.087099647772, 0.087516802635, 0.087935672566, 0.088355133178, 0.088775643463, 0.089196487176, 0.089618631657, 0.090042488015, 0.090465698863, 0.090890540547, 0.091317093553, 0.091742953864, 0.092171876781, 0.092600938856, 0.093030578046, 0.093461502793, 0.093892963907, 0.094326944688, 0.094760457456, 0.095195344257, 0.095630223672, 0.096066260144, 0.096503758275, 0.096943343674, 0.097382319611, 0.097822378915, 0.098264134343, 0.098706587057, 0.099149369310, 0.099594279848, 0.100039723737, 0.100485708173, 0.100932961990, 0.101381547878, 0.101829967278, 0.102280353377, 0.102730226590, 0.103181964688, 0.103634511679, 0.104087314378, 0.104542184713, 0.104997368260, 0.105453286279, 0.105910460458, 0.106368580787, 0.106827324288, 0.107287518635, 0.107749348820, 0.108211407650, 0.108675523162, 0.109139273619, 0.109604098902, 0.110070196106, 0.110538302061, 0.111005174139, 0.111473388477, 0.111942527947, 0.112413114844, 0.112883842180, 0.113355606892, 0.113828301181, 0.114300776523, 0.114776235565, 0.115252009574, 0.115728855272, 0.116206655301, 0.116684924756, 0.117163940472, 0.117644244712, 0.118125015173, 0.118605922394, 0.119088590475, 0.119571920839, 0.120057466811, 0.120542226920, 0.121028685070, 0.121515573240, 0.122003333989, 0.122491779054, 0.122980398236, 0.123469431004, 0.123961090939, 0.124453272825, 0.124945492546, 0.125439916036, 0.125934625038, 0.126430406050, 0.126926978806, 0.127424909705, 0.127923224408, 0.128421990710, 0.128921945263, 0.129422739301, 0.129924571878, 0.130426967993, 0.130930138528, 0.131437685112, 0.131942748404, 0.132449251669, 0.132956394120, 0.133463557831, 0.133971660364, 0.134479951397, 0.134989635684, 0.135500680056, 0.136012138193, 0.136525330520, 0.137039899272, 0.137554162776, 0.138069882271, 0.138586639162, 0.139104563257, 0.139623045250, 0.140140782082, 0.140660511788, 0.141180410963, 0.141701184232, 0.142223716156, 0.142747391544, 0.143271331295, 0.143796890256, 0.144322812914, 0.144848264443, 0.145376080497, 0.145904022872, 0.146433419470, 0.146962670853, 0.147493060812, 0.148024154273, 0.148555718850, 0.149088627211, 0.149622434233, 0.150157405696, 0.150692353092, 0.151228770440, 0.151765676846, 0.152304337580, 0.152843357274, 0.153381220477, 0.153920679166, 0.154461109375, 0.155003104678, 0.155545826317, 0.156089264973, 0.156633005148, 0.157176808253, 0.157722390942, 0.158268383521, 0.158815499104, 0.159363370620, 0.159911272021, 0.160460812540, 0.161011133623, 0.161561820292, 0.162114285435, 0.162667811756, 0.163220652513, 0.163774720324, 0.164330343619, 0.164885579114, 0.165442658123, 0.165999879796, 0.166557553489, 0.167116430453, 0.167674984792, 0.168234304486, 0.168795091732, 0.169355529021, 0.169918062303, 0.170480233218, 0.171044135995, 0.171609492748, 0.172174118998, 0.172740499701, 0.173307933291, 0.173876087852, 0.174444588188, 0.175013688767, 0.175584811412, 0.176155428634, 0.176727357635, 0.177300313873, 0.177874395924, 0.178447599721, 0.179021060646, 0.179595560214, 0.180173258448, 0.180749493952, 0.181326785810, 0.181905288410, 0.182483120110, 0.183062400837, 0.183643106865, 0.184223710269, 0.184805982203, 0.185389407226, 0.185974678035, 0.186559285183, 0.187144689691, 0.187730974691, 0.188318497807, 0.188904857639, 0.189492691475, 0.190080822851, 0.190670049534, 0.191259443962, 0.191850816784, 0.192441812206, 0.193033154389, 0.193626692651, 0.194219599194, 0.194813077561, 0.195408096792, 0.196003478151, 0.196599107275, 0.197195251460, 0.197793996986, 0.198391501532, 0.198989764709, 0.199590431342, 0.200190441117, 0.200791194801, 0.201392823463, 0.201995180230, 0.202598768305, 0.203203059821, 0.203808469005, 0.204414436136, 0.205019368510, 0.205625925633, 0.206233979075, 0.206841927084, 0.207451164261, 0.208060300845, 0.208670571767, 0.209281820287, 0.209892890262, 0.210504109741, 0.211116954715, 0.211731087886, 0.212345137548, 0.212960309107, 0.213576422880, 0.214192456837, 0.214810078049, 0.215427474097, 0.216045511255, 0.216664151674, 0.217283574310, 0.217903976356, 0.218525448381, 0.219146957870, 0.219770020801, 0.220393762475, 0.221017738789, 0.221641664415, 0.222268094678, 0.222893939495, 0.223521177233, 0.224149631441, 0.224777557156, 0.225404841298, 0.226034290547, 0.226664152072, 0.227294920355, 0.227925775181, 0.228557735903, 0.229188780095, 0.229821617646, 0.230455696164, 0.231089822667, 0.231724469171, 0.232359869452, 0.232997463964, 0.233635626278, 0.234273176962, 0.234911896939, 0.235549383828, 0.236189149224, 0.236829161677, 0.237469750283, 0.238111534639, 0.238751935135, 0.239395167769, 0.240038569531, 0.240682061900, 0.241327357015, 0.241972751103, 0.242618931087, 0.243264899107, 0.243912210142, 0.244558856048, 0.245206064449, 0.245854481287, 0.246504603142, 0.247153101744, 0.247807373437, 0.248459237738, 0.249111145517, 0.249762383058, 0.250415578825, 0.251069696567, 0.251723955075, 0.252378272738, 0.253033170956, 0.253689606560, 0.254347019486, 0.255004180645, 0.255662613124, 0.256320383307, 0.256978981724, 0.257639396278, 0.258299493346, 0.258960278308, 0.259621683292, 0.260282574800, 0.260945494457, 0.261609115032, 0.262272508705, 0.262937782148, 0.263603486567, 0.264268475524, 0.264934060631, 0.265601340615, 0.266268918617, 0.266937138925, 0.267606502362, 0.268276940946, 0.268946626772, 0.269618164532, 0.270289496061, 0.270959899891, 0.271632147877, 0.272305171392, 0.272978196440, 0.273650631785, 0.274323354376, 0.274998063031, 0.275673921485, 0.276348991926, 0.277025104294, 0.277701496283, 0.278380227854, 0.279057727101, 0.279735488471, 0.280415033164, 0.281095354106, 0.281775945585, 0.282457591159, 0.283139063822, 0.283821731055, 0.284505014010, 0.285188744733, 0.285873328787, 0.286558920572, 0.287243828478, 0.287928432702, 0.288614482468, 0.289300481415, 0.289986467385, 0.290674765845, 0.291363357547, 0.292051874600, 0.292742314101, 0.293433309311, 0.294124945768, 0.294816014903, 0.295506801263, 0.296199250767, 0.296891825290, 0.297586250567, 0.298279114947, 0.298973126241, 0.299669100371, 0.300365289666, 0.301062255388, 0.301758828398, 0.302457831980, 0.303155556849, 0.303854727942, 0.304554352183, 0.305255044780, 0.305954919441, 0.306655567421, 0.307356979186, 0.308058698458, 0.308762131642, 0.309464802827, 0.310166745492, 0.310871398944, 0.311573998244, 0.312278774505, 0.312984240104, 0.313689157071, 0.314395503123, 0.315103007938, 0.315812615312, 0.316520716942, 0.317230148099, 0.317938951936, 0.318647317056, 0.319357588742, 0.320067105472, 0.320776838807, 0.321489562157, 0.322201542832, 0.322912955207, 0.323625827642, 0.324339794245, 0.325053546681, 0.325767208849, 0.326481873457, 0.327197787073, 0.327914796851, 0.328631054177, 0.329348360305, 0.330066540185, 0.330785340391, 0.331504376301, 0.332223450195, 0.332942812684, 0.333662216963, 0.334382593556, 0.335103907747, 0.335825287053, 0.336546587471, 0.337269697796, 0.337991905099, 0.338715048246, 0.339439011301, 0.340164636593, 0.340888637191, 0.341614310251, 0.342339685697, 0.343065558986, 0.343790894489, 0.344517816901, 0.345245887875, 0.345972113384, 0.346700736363, 0.347430037731, 0.348158767999, 0.348889516386, 0.349620883446, 0.350351393140, 0.351082623902, 0.351814070589, 0.352547020215, 0.353280483006, 0.354014346621, 0.354748204055, 0.355482437175, 0.356217223318, 0.356953334783, 0.357689174980, 0.358424581798, 0.359161222119, 0.359900183318, 0.360637892153, 0.361375076664, 0.362113606601, 0.362852475220, 0.363592075355, 0.364330875149, 0.365070268534, 0.365811143716, 0.366550818433, 0.367292869032, 0.368034341089, 0.368778985791, 0.369522931452, 0.370265929830, 0.371010428062, 0.371754295272, 0.372500595855, 0.373246702007, 0.373992421456, 0.374737315657, 0.375483867036, 0.376231739062, 0.376980416990, 0.377728926510, 0.378478238309, 0.379226271298, 0.379975833488, 0.380725871383, 0.381476157324, 0.382226769940, 0.382976561111, 0.383728239543, 0.384479466550, 0.385232278894, 0.385984388384, 0.386737201639, 0.387492567515, 0.388247688673, 0.389006586306, 0.389760161613, 0.390515632742, 0.391271258117, 0.392028141583, 0.392784592205, 0.393542842067, 0.394300826768, 0.395058627809, 0.395817021989, 0.396576406241, 0.397336986192, 0.398096420849, 0.398857259814, 0.399618176237, 0.400379919379, 0.401139888083, 0.401902614469, 0.402665125229, 0.403428597775, 0.404190269017, 0.404954034272, 0.405718401172, 0.406482228896, 0.407246347491, 0.408012500721, 0.408778252460, 0.409544554979, 0.410312546027, 0.411080260340, 0.411846982163, 0.412615891041, 0.413383553571, 0.414151342457, 0.414920589347, 0.415690633409, 0.416461936484, 0.417231464313, 0.418004239945, 0.418775703678, 0.419547253063, 0.420317590927, 0.421091095557, 0.421865223132, 0.422639241291, 0.423413313922, 0.424188195468, 0.424962438877, 0.425736560290, 0.426511289433, 0.427285709401, 0.428060566660, 0.428836777394, 0.429611794582, 0.430388453516, 0.431165587305, 0.431943614009, 0.432723309924, 0.433502072597, 0.434282307421, 0.435061357739, 0.435842251149, 0.436621430881, 0.437401310747, 0.438184352927, 0.438965868002, 0.439747590522, 0.440530613684, 0.441313934182, 0.442098273439, 0.442881937829, 0.443666999561, 0.444451171190, 0.445237047013, 0.446021151634, 0.446805363464, 0.447592003009, 0.448378932238, 0.449166455553, 0.449954043659, 0.450741799758, 0.451530004681, 0.452318785024, 0.453106375956, 0.453896961334, 0.454686395149, 0.455476443505, 0.456266081421, 0.457056254491, 0.457846448290, 0.458637023650, 0.459428766831, 0.460221753029, 0.461012869143, 0.461805369880, 0.462598960307, 0.463392585912, 0.464185678933, 0.464979194148, 0.465774372579, 0.466576038888, 0.467370998652, 0.468165961252, 0.468960713059, 0.469756513285, 0.470554477675, 0.471350327710, 0.472148530158, 0.472947288728, 0.473745560500, 0.474543590046, 0.475342203760, 0.476142059942, 0.476942200317, 0.477743178258, 0.478542844096, 0.479342362521, 0.480145105729, 0.480947658268, 0.481749619961, 0.482551247319, 0.483352661309, 0.484155527770, 0.484957579547, 0.485761389076, 0.486565731717, 0.487370034056, 0.488173719757, 0.488978454968, 0.489783669976, 0.490588598648, 0.491395276432, 0.492201368582, 0.493009059705, 0.493814298175, 0.494622628986, 0.495430360678, 0.496237631415, 0.497046135429, 0.497856107769, 0.498667466225, 0.499477966136, 0.500288319969, 0.501099131261, 0.501908834114, 0.502719551891, 0.503530152109, 0.504341344275, 0.505151967046, 0.505966276375, 0.506777346024, 0.507591975083, 0.508403328171, 0.509217142567, 0.510029962031, 0.510843425557, 0.511657761040, 0.512474024906, 0.513289397338, 0.514107002903, 0.514923094642, 0.515740099071, 0.516557446336, 0.517375660072, 0.518192483832, 0.519010257067, 0.519828569681, 0.520646710344, 0.521465591441, 0.522284164695, 0.523104779091, 0.523925279639, 0.524744637765, 0.525566321035, 0.526385077165, 0.527204889959, 0.528027186379, 0.528851306880, 0.529673738565, 0.530495513396, 0.531318149444, 0.532139993821, 0.532962870390, 0.533786367906, 0.534609796841, 0.535433459307, 0.536258325570, 0.537085382255, 0.537911842439, 0.538736226749, 0.539563042388, 0.540390517237, 0.541215426524, 0.542042190081, 0.542870097013, 0.543696953249, 0.544522930839, 0.545351189716, 0.546178446787, 0.547007576478, 0.547838322788, 0.548667790054, 0.549498935377, 0.550331615782, 0.551160435498, 0.551990024179, 0.552818934600, 0.553651364792, 0.554484523279, 0.555316840350, 0.556150103870, 0.556982166226, 0.557815985817, 0.558647158869, 0.559479900511, 0.560312895962, 0.561147289838, 0.561982464603, 0.562816023585, 0.563651158496, 0.564486935364, 0.565322013934, 0.566157477821, 0.566991414464, 0.567828397267, 0.568663494305, 0.569499073652, 0.570337534602, 0.571175363658, 0.572013006298, 0.572852127978, 0.573690106016, 0.574529787226, 0.575370203089, 0.576208539437, 0.577049031817, 0.577889951182, 0.578730416074, 0.579572368765, 0.580412939030, 0.581254814801, 0.582097625555, 0.582941253877, 0.583783774262, 0.584625809545, 0.585467947304, 0.586311300717, 0.587155497696, 0.587999216371, 0.588844257186, 0.589690830658, 0.590538983152, 0.591384454373, 0.592227950573, 0.593073273693, 0.593918852465, 0.594763750763, 0.595609998161, 0.596457541038, 0.597303909465, 0.598150444795, 0.598996282295, 0.599843692600, 0.600690014870, 0.601538381745, 0.602388668130, 0.603238274437, 0.604087565151, 0.604937287108, 0.605786649367, 0.606639408770, 0.607488955773, 0.608338539349, 0.609188735677, 0.610039077979, 0.610892604307, 0.611743047205, 0.612593580013, 0.613444626202, 0.614296889529, 0.615149555705, 0.616002629196, 0.616855057639, 0.617708000922, 0.618561258491, 0.619414978967, 0.620269034969, 0.621124889981, 0.621980667144, 0.622835859708, 0.623692303211, 0.624547180642, 0.625403610341, 0.626258672272, 0.627114927803, 0.627971860674, 0.628830920101, 0.629688945957, 0.630544815776, 0.631405767563, 0.632268867150, 0.633128887912, 0.633988257427, 0.634846932847, 0.635706592323, 0.636567106678, 0.637426855007, 0.638287455194, 0.639149325181, 0.640010410291, 0.640872431057, 0.641733377821, 0.642593870133, 0.643454245856, 0.644315663128, 0.645177459047, 0.646037875479, 0.646900828172, 0.647763219701, 0.648625553822, 0.649489268388, 0.650352267527, 0.651217802026, 0.652082846702, 0.652947332677, 0.653811437128, 0.654675791639, 0.655541174311, 0.656406298330, 0.657268655385, 0.658133746133, 0.659000569476, 0.659866494949, 0.660735440930, 0.661601436776, 0.662470109201, 0.663337740022, 0.664205379844, 0.665073796486, 0.665941677005, 0.666809327256, 0.667677384902, 0.668545330629, 0.669413679611, 0.670284516315, 0.671153509866, 0.672021962009, 0.672893086748, 0.673764739060, 0.674633489859, 0.675503052088, 0.676374738179, 0.677246033358, 0.678116840682, 0.678991371871, 0.679864648273, 0.680736081989, 0.681608666919, 0.682480943087, 0.683353956676, 0.684228021358, 0.685101078785, 0.685976396421, 0.686850016488, 0.687726595100, 0.688599226148, 0.689474487284, 0.690349440372, 0.691222938784, 0.692097257160, 0.692971848641, 0.693851505529, 0.694728343806, 0.695606451791, 0.696483614075, 0.697362116732, 0.698238532645, 0.699113940348, 0.699990123654, 0.700868065136, 0.701746366677, 0.702627080926, 0.703505957677, 0.704383958357, 0.705262981976, 0.706141395297, 0.707020047081, 0.707899072602, 0.708778419439, 0.709658168280, 0.710536540575, 0.711416657249, 0.712297450524, 0.713178876022, 0.714061086842, 0.714945230827, 0.715825315917, 0.716707099878, 0.717588492843, 0.718469391123, 0.719356340140, 0.720236324436, 0.721121781269, 0.722005443000, 0.722889828008, 0.723772444821, 0.724656951181, 0.725541048815, 0.726425425498, 0.727310603107, 0.728194692245, 0.729077601112, 0.729963283383, 0.730849153784, 0.731737654446, 0.732623736549, 0.733510739087, 0.734395890954, 0.735281592238, 0.736167915954, 0.737054606206, 0.737943129504, 0.738829696960, 0.739716649296, 0.740604444037, 0.741493130541, 0.742379283728, 0.743268544865, 0.744157907734, 0.745047970435, 0.745934136538, 0.746824609264, 0.747717016022, 0.748606914893, 0.749493830503, 0.750384548355, 0.751275357706, 0.752167352492, 0.753057759512, 0.753951265065, 0.754841512937, 0.755730882169, 0.756624300765, 0.757516507170, 0.758407970676, 0.759301043269, 0.760192360822, 0.761085726880, 0.761979478323, 0.762871872882, 0.763766659514, 0.764659726688, 0.765554221519, 0.766448969582, 0.767342171640, 0.768238309703, 0.769132051262, 0.770026976031, 0.770923187482, 0.771818887005, 0.772716200908, 0.773612548969, 0.774510156522, 0.775406181852, 0.776299902959, 0.777198394674, 0.778094760163, 0.778992747176, 0.779889800666, 0.780788218025, 0.781685137760, 0.782582697465, 0.783480293187, 0.784378793624, 0.785277170140, 0.786173448575, 0.787072843964, 0.787972546566, 0.788874090190, 0.789771471139, 0.790668516331, 0.791568950115, 0.792472003359, 0.793370658361, 0.794273045420, 0.795173235636, 0.796075795274, 0.796977483386, 0.797878985979, 0.798778311451, 0.799679749526, 0.800583251881, 0.801485049551, 0.802385185846, 0.803286148007, 0.804190883914, 0.805095052044, 0.805996380581, 0.806898489465, 0.807802462146, 0.808709793691, 0.809614822458, 0.810518473343, 0.811424599189, 0.812325467454, 0.813233044660, 0.814140000130, 0.815045775650, 0.815950317221, 0.816859672275, 0.817767743889, 0.818673227011, 0.819578156776, 0.820481641031, 0.821389069802, 0.822296229355, 0.823203388635, 0.824111351831, 0.825017545321, 0.825924691177, 0.826831634216, 0.827740963002, 0.828647091994, 0.829556274133, 0.830462584249, 0.831368810306, 0.832278156233, 0.833187431400, 0.834097600598, 0.835007162512, 0.835913811490, 0.836824418218, 0.837733261858, 0.838641579306, 0.839550323733, 0.840461987283, 0.841370405459, 0.842279123934, 0.843188476575, 0.844099979917, 0.845009886747, 0.845919194203, 0.846828990247, 0.847736331585, 0.848648591569, 0.849561454180, 0.850477184549, 0.851390818488, 0.852305720140, 0.853216819853, 0.854130167089, 0.855041859031, 0.855956382156, 0.856870142573, 0.857782445888, 0.858693357270, 0.859610354484, 0.860526762543, 0.861440471284, 0.862354939168, 0.863269634383, 0.864184417704, 0.865099139185, 0.866015947655, 0.866931684231, 0.867846011263, 0.868764588363, 0.869683889956, 0.870602281986, 0.871519822390, 0.872436346241, 0.873353156977, 0.874271208200, 0.875191240035, 0.876110074233, 0.877024001711, 0.877945148638, 0.878863166966, 0.879782734789, 0.880698449500, 0.881615474112, 0.882536059567, 0.883455877509, 0.884375461340, 0.885295655825, 0.886216270224, 0.887136094058, 0.888056239102, 0.888974907278, 0.889891471658, 0.890808822785, 0.891730153274, 0.892652906555, 0.893573056256, 0.894492533411, 0.895411111138, 0.896331259593, 0.897252545868, 0.898174365027, 0.899096086322, 0.900023442124, 0.900945210940, 0.901867083305, 0.902789770936, 0.903711646536, 0.904633056396, 0.905553609029, 0.906475462268, 0.907397139461, 0.908321518808, 0.909248155300, 0.910174986204, 0.911097425493, 0.912023378054, 0.912946738432, 0.913871599583, 0.914796828315, 0.915721160129, 0.916647001039, 0.917572617895, 0.918500902914, 0.919426666905, 0.920352015227, 0.921274999108, 0.922200173838, 0.923127272732, 0.924055155362, 0.924981202612, 0.925905757331, 0.926829499298, 0.927753448784, 0.928681144418, 0.929609038471, 0.930534877754, 0.931463273572, 0.932393989086, 0.933318646638, 0.934248442471, 0.935180734475, 0.936109952443, 0.937039404623, 0.937969558867, 0.938901864246, 0.939830363781, 0.940758071406, 0.941687571312, 0.942615308911, 0.943547846990, 0.944478226040, 0.945405424415, 0.946336448855, 0.947266396375, 0.948194773932, 0.949124974201, 0.950053903862, 0.950984863502, 0.951916392345, 0.952845127424, 0.953774404303, 0.954703087490, 0.955633690221, 0.956562043425, 0.957490833767, 0.958424625843, 0.959357289779, 0.960290989907, 0.961223166666, 0.962153467562, 0.963085753564, 0.964019889273, 0.964952821548, 0.965886168319, 0.966821396524, 0.967755175350, 0.968687442599, 0.969624039825, 0.970557235330, 0.971490939640, 0.972424833806, 0.973363538452, 0.974295115428, 0.975231882359, 0.976167015413, 0.977100371937, 0.978038368719, 0.978972271984, 0.979908440778, 0.980840111900, 0.981776013978, 0.982713887194, 0.983649130859, 0.984582909766, 0.985519166981, 0.986456011477, 0.987393771753, 0.988328661260, 0.989267224248, 0.990202588792, 0.991139602026, 0.992077037858, 0.993019578599, 0.993957663533, 0.994895036416, 0.995835275676, 0.996773438120, 0.997709734820, 0.998645651610, 0.999584808884, 1.000522409873, 1.001464758579, 1.002403335167, 1.003344819859, 1.004285749069, 1.005229288320, 1.006167943917, 1.007107643785, 1.008047766502, 1.008986350684, 1.009924652717, 1.010863276425, 1.011806958794, 1.012749163012, 1.013689130965, 1.014626165419, 1.015567481209, 1.016505888392, 1.017445052732, 1.018384992876, 1.019326326975, 1.020269687649, 1.021213451075, 1.022160476413, 1.023104336105, 1.024043320377, 1.024986970493, 1.025929835614, 1.026873453906, 1.027817635986, 1.028759313835, 1.029702159036, 1.030645736675, 1.031585876230, 1.032524315493, 1.033467986116, 1.034412827886, 1.035356290454, 1.036301268809, 1.037244314141, 1.038186789308, 1.039129688889, 1.040072349177, 1.041016554064, 1.041961008068, 1.042904102379, 1.043848706239, 1.044795335417, 1.045741233775, 1.046683853218, 1.047626347035, 1.048569569613, 1.049517142618, 1.050463407111, 1.051406047435, 1.052352918419, 1.053297684806, 1.054242995345, 1.055191674827, 1.056136633908, 1.057082578764, 1.058030136791, 1.058976523282, 1.059923201745, 1.060868995438, 1.061816883393, 1.062762152701, 1.063705903202, 1.064655564669, 1.065602210693, 1.066550975282, 1.067496962020, 1.068444144409, 1.069394619985, 1.070343973358, 1.071290416748, 1.072237793713, 1.073187154438, 1.074140449496, 1.075093492642, 1.076044026901, 1.076995857968, 1.077943554026, 1.078890926912, 1.079840626675, 1.080792454839, 1.081742662066, 1.082695930064, 1.083645050290, 1.084595357496, 1.085548203502, 1.086499074169, 1.087454187867, 1.088405859680, 1.089356356398, 1.090309547471, 1.091259240899, 1.092206450925, 1.093155069382, 1.094109016857, 1.095062458810, 1.096010251521, 1.096962120614, 1.097914105382, 1.098866403766, 1.099820991637, 1.100776323980, 1.101726988673, 1.102679199799, 1.103629517126, 1.104580294013, 1.105532387231, 1.106482115887, 1.107436089552, 1.108388562466, 1.109342497633, 1.110295739033, 1.111249668982, 1.112202594814, 1.113156387584, 1.114113787700, 1.115067721810, 1.116018423335, 1.116975997652, 1.117928525700, 1.118880737713, 1.119840587652, 1.120790559023, 1.121748326591, 1.122698567345, 1.123656474755, 1.124611031640, 1.125567366471, 1.126522952448, 1.127477541163, 1.128431664264, 1.129386870205, 1.130344504170, 1.131298584026, 1.132258680649, 1.133215161777, 1.134171695658, 1.135126754479, 1.136080840716, 1.137035360401, 1.137989846331, 1.138946105715, 1.139902773121, 1.140855804029, 1.141816330860, 1.142768820606, 1.143726046284, 1.144680689996, 1.145637831542, 1.146595918549, 1.147553543488, 1.148511059321, 1.149472817024, 1.150431827431, 1.151385661299, 1.152345560543, 1.153305867659, 1.154262838804, 1.155224984239, 1.156183758977, 1.157133261427, 1.158094731698, 1.159056005107, 1.160019818987, 1.160979786896, 1.161938003502, 1.162897188822, 1.163858028618, 1.164818459896, 1.165776261413, 1.166739807680, 1.167697784006, 1.168654516114, 1.169612904839, 1.170569193703, 1.171531712121, 1.172492763423, 1.173450921580, 1.174406279504, 1.175366800313, 1.176323603600, 1.177282819981, 1.178242300368, 1.179204121737, 1.180163843125, 1.181123633851, 1.182080715613, 1.183044690510, 1.184007525873, 1.184966988603, 1.185930101697, 1.186894240076, 1.187857705716, 1.188821428514, 1.189783159915, 1.190745442340, 1.191703122263, 1.192666816940, 1.193627181015, 1.194594377635, 1.195557914668, 1.196521511580, 1.197491227271, 1.198455561581, 1.199421375151, 1.200380433607, 1.201344901160, 1.202305827273, 1.203269210209, 1.204233073765, 1.205195027065, 1.206158354828, 1.207124069036, 1.208086500734, 1.209054605125, 1.210020547488, 1.210991848072, 1.211958257598, 1.212926014035, 1.213891916399, 1.214859259521, 1.215819860361, 1.216786260801, 1.217752995230, 1.218717028479, 1.219679612518, 1.220644840074, 1.221614069629, 1.222578252868, 1.223542525036, 1.224505607490, 1.225473085757, 1.226437347441, 1.227404766572, 1.228378761567, 1.229345180731, 1.230308013075, 1.231276638912, 1.232245843431, 1.233212311824, 1.234178388997, 1.235145239268, 1.236113072448, 1.237078134510, 1.238043550180, 1.239013277821, 1.239977628745, 1.240942371026, 1.241908230222, 1.242877883276, 1.243846378661, 1.244814291383, 1.245783907211, 1.246751453675, 1.247716720016, 1.248688040264, 1.249658134623, 1.250624522815, 1.251597624020, 1.252564527982, 1.253535201161, 1.254502258514, 1.255469518834, 1.256436508372, 1.257402324947, 1.258368467573, 1.259338760675, 1.260306433904, 1.261280461007, 1.262250934052, 1.263217203266, 1.264185355844, 1.265155743478, 1.266130011437, 1.267097874438, 1.268070845916, 1.269044259143, 1.270011545371, 1.270978178301, 1.271953758631, 1.272931046895, 1.273897409919, 1.274869640879, 1.275844167954, 1.276815629181, 1.277779322610, 1.278751323084, 1.279729837431, 1.280705438975, 1.281676848671, 1.282649428602, 1.283620878648, 1.284589806390, 1.285561235960, 1.286533297555, 1.287505325403, 1.288481693956, 1.289451373706, 1.290415357417, 1.291389148200, 1.292366728052, 1.293336955869, 1.294308286998, 1.295277818207, 1.296250721382, 1.297220195168, 1.298192536925, 1.299165261614, 1.300138098213, 1.301116793746, 1.302090333494, 1.303065947042, 1.304041281959, 1.305019653746, 1.305996166734, 1.306974202334, 1.307951135907, 1.308928272979, 1.309900228803, 1.310878771751, 1.311854867114, 1.312832473996, 1.313804432527, 1.314779171871, 1.315754809957, 1.316727106277, 1.317702062621, 1.318675702700, 1.319652074463, 1.320632481786, 1.321612848704, 1.322589463961, 1.323568335213, 1.324547473886, 1.325526821661, 1.326509543438, 1.327485659374, 1.328460217757, 1.329434343714, 1.330413075960, 1.331390227725, 1.332366203451, 1.333344246401, 1.334323559106, 1.335304126391, 1.336280996654, 1.337255886576, 1.338233565976, 1.339211175099, 1.340190362450, 1.341165941571, 1.342139105297, 1.343116100573, 1.344093581820, 1.345068461623, 1.346049282547, 1.347024163989, 1.348000048325, 1.348981253813, 1.349964010361, 1.350943870547, 1.351921297983, 1.352898718152, 1.353878284200, 1.354858137610, 1.355839980112, 1.356819691804, 1.357795580322, 1.358778131164, 1.359756436690, 1.360741445647, 1.361721112322, 1.362697067288, 1.363679324166, 1.364655048533, 1.365637203185, 1.366612795388, 1.367591171304, 1.368571634508, 1.369546983409, 1.370525516440, 1.371505707562, 1.372481460428, 1.373456485769, 1.374433324488, 1.375411221075, 1.376391283225, 1.377377321010, 1.378354144509, 1.379329654185, 1.380310925472, 1.381296874469, 1.382278259522, 1.383262727975, 1.384242174073, 1.385224825213, 1.386203934673, 1.387188412874, 1.388168259331, 1.389148214817, 1.390132913850, 1.391114282108, 1.392089260377, 1.393066754456, 1.394045130073, 1.395024549993, 1.396009318350, 1.396987040294, 1.397967380944, 1.398949569533, 1.399930908594, 1.400915300959, 1.401897459287, 1.402876627570, 1.403860561968, 1.404836195754, 1.405820162345, 1.406805454674, 1.407786412439, 1.408763869393, 1.409748183535, 1.410726570887, 1.411713174369, 1.412701136874, 1.413683944041, 1.414667784269, 1.415651155396, 1.416635545378, 1.417622512512, 1.418601775894, 1.419588821513, 1.420569262383, 1.421551566272, 1.422534258613, 1.423517705517, 1.424506766333, 1.425490495812, 1.426472168586, 1.427455646995, 1.428441555592, 1.429426380398, 1.430416251598, 1.431394288316, 1.432382630640, 1.433371448480, 1.434356513014, 1.435341414671, 1.436331116919, 1.437319929427, 1.438310057140, 1.439294816317, 1.440279467541, 1.441258342766, 1.442244671497, 1.443234956609, 1.444218881006, 1.445201420030, 1.446187800673, 1.447174639081, 1.448157387028, 1.449141789693, 1.450125845332, 1.451111031348, 1.452089011506, 1.453076755380, 1.454058930071, 1.455038811236, 1.456026431179, 1.457009871028, 1.457993622893, 1.458981270130, 1.459971168537, 1.460956650254, 1.461946474367, 1.462932443679, 1.463921705530, 1.464904510640, 1.465884136328, 1.466879236150, 1.467865372024, 1.468849469841, 1.469837276051, 1.470824496223, 1.471807440278, 1.472796406595, 1.473784126305, 1.474773488462, 1.475768785303, 1.476759572919, 1.477752091100, 1.478735243746, 1.479724061156, 1.480712875505, 1.481699377414, 1.482688336420, 1.483675875475, 1.484658638978, 1.485638502654, 1.486623511857, 1.487612068250, 1.488599789685, 1.489594871226, 1.490581432934, 1.491573258061, 1.492564207848, 1.493553296855, 1.494541754886, 1.495530850220, 1.496516535766, 1.497506771224, 1.498494014203, 1.499478074501, 1.500465428248, 1.501452716845, 1.502448110899, 1.503436392326, 1.504423362135, 1.505410577590, 1.506396000394, 1.507381834281, 1.508365052805, 1.509355792614, 1.510338389913, 1.511330320027, 1.512319180203, 1.513306955100, 1.514302361018, 1.515295956572, 1.516283018815, 1.517280832671, 1.518271862995, 1.519257493636, 1.520248071389, 1.521243134698, 1.522231708919, 1.523226624164, 1.524214733577, 1.525204354006, 1.526189072221, 1.527171437465, 1.528170024124, 1.529168488684, 1.530164587348, 1.531151171860, 1.532143966032, 1.533127621702, 1.534110524344, 1.535104366997, 1.536094325622, 1.537079037173, 1.538067950432, 1.539063913963, 1.540053988872, 1.541050673401, 1.542039952423, 1.543028972841, 1.544024734321, 1.545012029255, 1.546007314607, 1.546997494623, 1.547987023136, 1.548978903653, 1.549966768122, 1.550959263034, 1.551955285082, 1.552947034099, 1.553936434773, 1.554926723131, 1.555915588234, 1.556904079807, 1.557893476893, 1.558887665333, 1.559879847090, 1.560875138164, 1.561869294775, 1.562860557902, 1.563857047786, 1.564851093568, 1.565841922037, 1.566827359487, 1.567815198580, 1.568811211336, 1.569808400739, 1.570804328653, 1.571801346691, 1.572788235477, 1.573788130890, 1.574785406759, 1.575777995734, 1.576769974004, 1.577765833067, 1.578761379552, 1.579752084484, 1.580739480697, 1.581732309786, 1.582717162154, 1.583709615206, 1.584709582496, 1.585708461357, 1.586704595654, 1.587696751097, 1.588689274832, 1.589686874757, 1.590680679355, 1.591674219113, 1.592666398985, 1.593665979542, 1.594660044135, 1.595659042505, 1.596652074490, 1.597651131069, 1.598643442379, 1.599636937861, 1.600630166017, 1.601617774700, 1.602612869911, 1.603605595587, 1.604596594211, 1.605595936570, 1.606590562155, 1.607579008099, 1.608577803001, 1.609576036873, 1.610570635736, 1.611567890490, 1.612567547258, 1.613556952349, 1.614552853883, 1.615550345542, 1.616562185631, 1.617559087492, 1.618555377819, 1.619541407555, 1.620536551084, 1.621530110260, 1.622527167740, 1.623520149637, 1.624520456185, 1.625515462629, 1.626517238197, 1.627520095695, 1.628511556293, 1.629502842843, 1.630494190080, 1.631487861175, 1.632480159086, 1.633478146969, 1.634473248569, 1.635465363692, 1.636464169323, 1.637458661978, 1.638448277775, 1.639442387861, 1.640437260568, 1.641428901673, 1.642423326932, 1.643422728571, 1.644422539326, 1.645426211979, 1.646421915986, 1.647407566297, 1.648405953733, 1.649398408288, 1.650397621308, 1.651391353936, 1.652389530781, 1.653386976677, 1.654388227593, 1.655385133754, 1.656382738936, 1.657385699406, 1.658383307807, 1.659374808234, 1.660376226677, 1.661377350818, 1.662374500638, 1.663369363876, 1.664366511338, 1.665367581469, 1.666364236256, 1.667356460119, 1.668355873584, 1.669355279933, 1.670349042409, 1.671354029424, 1.672343251953, 1.673342102130, 1.674342187980, 1.675334854988, 1.676334619301, 1.677331112382, 1.678327391662, 1.679318656192, 1.680320405137, 1.681317777231, 1.682322000854, 1.683321409215, 1.684315018907, 1.685314210732, 1.686306869628, 1.687298356967, 1.688301880511, 1.689305158603, 1.690298690746, 1.691294351658, 1.692298801062, 1.693289675726, 1.694294910636, 1.695291647180, 1.696284979176, 1.697287164059, 1.698292989483, 1.699295086819, 1.700282534455, 1.701281226462, 1.702270251605, 1.703274429456, 1.704277109720, 1.705281537439, 1.706279526409, 1.707271912153, 1.708275111885, 1.709275675634, 1.710272955743, 1.711269761217, 1.712267695577, 1.713263125843, 1.714264891678, 1.715272896456, 1.716272263847, 1.717281025561, 1.718278945806, 1.719281257751, 1.720284086535, 1.721284961538, 1.722282466223, 1.723286194326, 1.724278897078, 1.725277680792, 1.726271899014, 1.727262442109, 1.728272462999, 1.729273009904, 1.730267279006, 1.731277419965, 1.732275027783, 1.733276342434, 1.734274812774, 1.735281367959, 1.736278855894, 1.737278616419, 1.738275168678, 1.739281780406, 1.740287243896, 1.741293341042, 1.742295200295, 1.743288915179, 1.744291874805, 1.745303389130, 1.746301960659, 1.747294120301, 1.748300082489, 1.749299017120, 1.750299888146, 1.751306304525, 1.752311179240, 1.753314865294, 1.754317521700, 1.755321311580, 1.756329558045, 1.757330538302, 1.758337739577, 1.759333083839, 1.760333090499, 1.761337159854, 1.762337751826, 1.763336373172, 1.764336892183, 1.765343516540, 1.766335666639, 1.767346863667, 1.768360420630, 1.769366568214, 1.770364788621, 1.771363359037, 1.772363870773, 1.773368110246, 1.774363827924, 1.775361108813, 1.776373530720, 1.777370786599, 1.778383008317, 1.779383431144, 1.780386504324, 1.781393973606, 1.782394654547, 1.783391923337, 1.784390879423, 1.785393383891, 1.786392020061, 1.787386862327, 1.788390551786, 1.789394881533, 1.790400011456, 1.791399467258, 1.792401120740, 1.793405386649, 1.794408111824, 1.795419990352, 1.796434721205, 1.797439134031, 1.798438011643, 1.799439109872, 1.800452258861, 1.801456311693, 1.802461010161, 1.803468314510, 1.804467495407, 1.805472809851, 1.806477397281, 1.807479100621, 1.808488289158, 1.809493078345, 1.810504661203, 1.811509207598, 1.812507678439, 1.813518965887, 1.814535928850, 1.815538635874, 1.816540304338, 1.817546542433, 1.818541476197, 1.819542850676, 1.820549297469, 1.821557765288, 1.822563702356, 1.823564944781, 1.824565427183, 1.825566766527, 1.826567972971, 1.827570587751, 1.828572039918, 1.829572052063, 1.830570785284, 1.831587409604, 1.832595549199, 1.833608047723, 1.834616947407, 1.835626292632, 1.836617566938, 1.837623390272, 1.838627145674, 1.839621578576, 1.840632194900, 1.841643630294, 1.842649899468, 1.843652960738, 1.844665572161, 1.845683259426, 1.846691528607, 1.847692999705, 1.848703222805, 1.849705754140, 1.850707371516, 1.851699482339, 1.852708157878, 1.853700265820, 1.854705865308, 1.855716789331, 1.856717332050, 1.857723846956, 1.858734644923, 1.859741292185, 1.860736915482, 1.861746229837, 1.862764861094, 1.863770780455, 1.864775663306, 1.865781218733, 1.866786263593, 1.867791141271, 1.868802715847, 1.869806064667, 1.870817478304, 1.871820389895, 1.872818817898, 1.873817013181, 1.874820275114, 1.875830656606, 1.876840156338, 1.877859220194, 1.878858763952, 1.879877377249, 1.880887392132, 1.881895989295, 1.882898510211, 1.883899295182, 1.884900492527, 1.885912819817, 1.886922424258, 1.887930120110, 1.888929297270, 1.889942542295, 1.890940453100, 1.891948890938, 1.892952887385, 1.893954855606, 1.894968688697, 1.895974093382, 1.896974499790, 1.897968769174, 1.898977985010, 1.899984239117, 1.900977305684, 1.901978816730, 1.902982538458, 1.903995409214, 1.905005413129, 1.906003045678, 1.907008689692, 1.908016773183, 1.909020122556, 1.910029219726, 1.911039640850, 1.912042983897, 1.913056222419, 1.914059679637, 1.915066282203, 1.916072824818, 1.917077256431, 1.918083584900, 1.919088934239, 1.920099940336, 1.921106207244, 1.922121381340, 1.923114738592, 1.924120474501, 1.925138523744, 1.926148779535, 1.927154193090, 1.928158406044, 1.929165758064, 1.930169054697, 1.931176194140, 1.932189649763, 1.933188495557, 1.934204312019, 1.935210089562, 1.936211039816, 1.937224487643, 1.938236877749, 1.939242155624, 1.940245564832, 1.941254598044, 1.942262292842, 1.943268253318, 1.944276625721, 1.945286885412, 1.946297658511, 1.947300211103, 1.948317768439, 1.949332652726, 1.950350107897, 1.951356401288, 1.952354913219, 1.953363254552, 1.954371596599, 1.955374800181, 1.956386532678, 1.957399131487, 1.958409834721, 1.959412136610, 1.960428373528, 1.961431813455, 1.962436581431, 1.963452981655, 1.964464763086, 1.965473171351, 1.966476850706, 1.967487730898, 1.968501090716, 1.969509209350, 1.970525638745, 1.971539693626, 1.972546622515, 1.973554338608, 1.974564971821, 1.975583545802, 1.976594679147, 1.977605367341, 1.978618743836, 1.979626492632, 1.980629690414, 1.981637958569, 1.982649365833, 1.983651966587, 1.984648377208, 1.985655440758, 1.986669941056, 1.987664105446, 1.988673286398, 1.989676717623, 1.990688721296, 1.991692096044, 1.992703047101, 1.993703943494, 1.994717577577, 1.995722100841, 1.996735634060, 1.997745532664, 1.998749990235, 1.999755083555, 2.000766423909, 2.001778642202, 2.002787980076, 2.003795288117, 2.004813413171, 2.005821386044, 2.006825659681, 2.007846765341, 2.008864515176, 2.009877858044, 2.010891166131, 2.011913673483, 2.012912196484, 2.013923693395, 2.014935079444, 2.015950493311, 2.016964764403, 2.017957600903, 2.018966730851, 2.019980393837, 2.020986491821, 2.021999311811, 2.023002547523, 2.023997871011, 2.025002243234, 2.026022684544, 2.027031062570, 2.028052766808, 2.029062206593, 2.030063665590, 2.031086659465, 2.032082889994, 2.033083520050, 2.034089279454, 2.035097703044, 2.036098987604, 2.037107978558, 2.038116284683, 2.039130454028, 2.040132706783, 2.041149118610, 2.042159779062, 2.043163635034, 2.044174768911, 2.045192599704, 2.046199344490, 2.047210752445, 2.048214184588, 2.049222275303, 2.050238123269, 2.051249999925, 2.052245234565, 2.053259005504, 2.054261861736, 2.055280652507, 2.056287797850, 2.057287026606, 2.058298741670, 2.059310728131, 2.060318641033, 2.061326147648, 2.062337099903, 2.063342923785, 2.064370524178, 2.065379912837, 2.066381886726, 2.067405045514, 2.068410131987, 2.069408581626, 2.070415001177, 2.071411881844, 2.072423832201, 2.073432694394, 2.074441843630, 2.075457786874, 2.076480411018, 2.077498803777, 2.078496486846, 2.079518632788, 2.080524160957, 2.081533699389, 2.082531619749, 2.083541998054, 2.084551935977, 2.085566660759, 2.086575068232, 2.087586353958, 2.088614646211, 2.089621569333, 2.090633935908, 2.091646576200, 2.092646799881, 2.093659085021, 2.094669306551, 2.095676578406, 2.096684184315, 2.097689456330, 2.098682724607, 2.099690070874, 2.100697623379, 2.101706915029, 2.102719713092, 2.103733996231, 2.104747503039, 2.105748975568, 2.106756094692, 2.107772456733, 2.108774577070, 2.109778511842, 2.110774740050, 2.111802471718, 2.112812087070, 2.113822587461, 2.114842233509, 2.115855149158, 2.116853324249, 2.117862800000, 2.118875427091, 2.119885212294, 2.120903432055, 2.121914843921, 2.122926772296, 2.123941988394, 2.124961841689, 2.125982354183, 2.126991248678, 2.127998468582, 2.129009315796, 2.130034004661, 2.131045143200, 2.132061761154, 2.133085425607, 2.134105771881, 2.135125616297, 2.136133184755, 2.137144049121, 2.138163241458, 2.139167598551, 2.140174343698, 2.141181624636, 2.142186366813, 2.143209080607, 2.144224159490, 2.145232272285, 2.146252157405, 2.147261763336, 2.148281359634, 2.149294841975, 2.150303451025, 2.151312623516, 2.152322172898, 2.153330612698, 2.154338053628, 2.155341999066, 2.156349391157, 2.157357128829, 2.158364521124, 2.159359399689, 2.160374594045, 2.161385302700, 2.162387827715, 2.163400137878, 2.164427179484, 2.165442353005, 2.166456846751, 2.167472885627, 2.168491435146, 2.169514817351, 2.170531097983, 2.171541252871, 2.172566298452, 2.173578030150, 2.174596214242, 2.175620955725, 2.176643424489, 2.177673536680, 2.178675427520, 2.179680226169, 2.180687354946, 2.181701839989, 2.182701965762, 2.183728532508, 2.184742313134, 2.185759531552, 2.186774798277, 2.187797797882, 2.188808857537, 2.189824764250, 2.190840828929, 2.191839146065, 2.192856960874, 2.193886736057, 2.194902699822, 2.195917568295, 2.196940144984, 2.197975617133, 2.198988773099, 2.199992046289, 2.201000678071, 2.202021961377, 2.203036294503, 2.204059949936, 2.205079825801, 2.206098751705, 2.207096285061, 2.208116381883, 2.209119057731, 2.210124194769, 2.211144445571, 2.212158108993, 2.213172866773, 2.214205860568, 2.215217578498, 2.216250665514, 2.217264656986, 2.218269678760, 2.219288472705, 2.220290693401, 2.221298919089, 2.222323330781, 2.223339197108, 2.224352276347, 2.225367943237, 2.226375092888, 2.227380111650, 2.228394515744, 2.229388315381, 2.230401963655, 2.231417095392, 2.232423628701, 2.233432574539, 2.234451097561, 2.235471118608, 2.236477669951, 2.237490986727, 2.238501709346, 2.239503330939, 2.240532354556, 2.241553293862, 2.242564795360, 2.243584289370, 2.244590698402, 2.245598221830, 2.246624793471, 2.247641814684, 2.248657065579, 2.249663659953, 2.250663929756, 2.251694419184, 2.252703423628, 2.253706131720, 2.254724588841, 2.255729710335, 2.256731123843, 2.257736817230, 2.258757216762, 2.259775912041, 2.260789797671, 2.261801770508, 2.262823345353, 2.263825164000, 2.264839607767, 2.265844812100, 2.266845925661, 2.267858284525, 2.268852360329, 2.269872808099, 2.270872078288, 2.271888516348, 2.272905873473, 2.273926680256, 2.274929198097, 2.275948710509, 2.276961500224, 2.277968914755, 2.278990312514, 2.279996905252, 2.281009651942, 2.282032830451, 2.283032841225, 2.284048190133, 2.285057713538, 2.286053225407, 2.287046819122, 2.288054324768, 2.289058934261, 2.290080186445, 2.291089498568, 2.292102608381, 2.293111007992, 2.294110555675, 2.295135629993, 2.296161325835, 2.297161814500, 2.298175830872, 2.299188241752, 2.300201197678, 2.301223647312, 2.302237883329, 2.303264271766, 2.304272613875, 2.305278478367, 2.306287205439, 2.307304978214, 2.308311627434, 2.309339473311, 2.310345086578, 2.311349031106, 2.312362969088, 2.313360603114, 2.314384538563, 2.315372468328, 2.316385686096, 2.317404159520, 2.318427558450, 2.319432169452, 2.320441925412, 2.321443655995, 2.322452447499, 2.323463862099, 2.324485338614, 2.325484412309, 2.326489013215, 2.327506375244, 2.328510950476, 2.329532694651, 2.330543274949, 2.331556305484, 2.332577121877, 2.333583116996, 2.334592855365, 2.335600996498, 2.336623644930, 2.337651069627, 2.338666911606, 2.339670988787, 2.340683197645, 2.341699011179, 2.342714337991, 2.343717669321, 2.344724092391, 2.345719952371, 2.346726110275, 2.347737989857, 2.348739435684, 2.349756702792, 2.350775481682, 2.351789529176, 2.352810549230, 2.353820341251, 2.354844087890, 2.355858624784, 2.356861609783, 2.357854145213, 2.358863044750, 2.359880360742, 2.360872450012, 2.361894089706, 2.362906119126, 2.363924528124, 2.364932953502, 2.365942414675, 2.366964337501, 2.367990900095, 2.369008213538, 2.370021399467, 2.371033688983, 2.372055913058, 2.373082189415, 2.374104113481, 2.375115158768, 2.376132073979, 2.377155515896, 2.378174839452, 2.379186162452, 2.380227050318, 2.381244318397, 2.382254445041, 2.383251181911, 2.384270832208, 2.385285078327, 2.386287746120, 2.387281292339, 2.388293788029, 2.389324083376, 2.390359069337, 2.391371716016, 2.392386836603, 2.393403905712, 2.394424654653, 2.395453529889, 2.396467210451, 2.397466451671, 2.398475715459, 2.399483625492, 2.400498794781, 2.401508132325, 2.402528928362, 2.403529582909, 2.404540375089, 2.405567226757, 2.406577683382, 2.407597157328, 2.408611463311, 2.409627697638, 2.410654141366, 2.411668112518, 2.412673673709, 2.413690687745, 2.414722163894, 2.415739580452, 2.416745213593, 2.417749771503, 2.418769188102, 2.419779014875, 2.420795086139, 2.421830288020, 2.422843356626, 2.423830558447, 2.424836180538, 2.425847609957, 2.426850375771, 2.427883146394, 2.428871624317, 2.429881521817, 2.430886393886, 2.431887139327, 2.432906435206, 2.433907133723, 2.434926103492, 2.435950432258, 2.436955564397, 2.437961718612, 2.438978204528, 2.439985113939, 2.440991126305, 2.441992745161, 2.443007037612, 2.444034932334, 2.445035254265, 2.446046497205, 2.447066908548, 2.448070712775, 2.449072323324, 2.450103795116, 2.451117964764, 2.452124176055, 2.453136176087, 2.454159684542, 2.455169755605, 2.456180567323, 2.457184901657, 2.458201417859, 2.459215942986, 2.460222567582, 2.461243840511, 2.462263743484, 2.463280115740, 2.464315949724, 2.465338536563, 2.466360613557, 2.467357201899, 2.468366425023, 2.469381582487, 2.470385263629, 2.471399113203, 2.472415464008, 2.473444017207, 2.474448077629, 2.475446028532, 2.476450960998, 2.477456137987, 2.478453583719, 2.479466555739, 2.480489511487, 2.481520673456, 2.482513393672, 2.483528485969, 2.484542245477, 2.485536059652, 2.486551193015, 2.487558428398, 2.488544596466, 2.489570280791, 2.490604709841, 2.491619246015, 2.492640883720, 2.493662629507, 2.494672133096, 2.495692419424, 2.496669722451, 2.497681877024, 2.498699134504, 2.499707251722, 2.500700106782, 2.501721433522, 2.502724851842, 2.503752480720, 2.504773799061, 2.505787503911, 2.506810276691, 2.507825395716, 2.508849480655, 2.509852102678, 2.510858452861, 2.511868834071, 2.512905480577, 2.513939785687, 2.514953389376, 2.515975348456, 2.517004716701, 2.518014198219, 2.519014983359, 2.520033468443, 2.521037337333, 2.522064338620, 2.523081608176, 2.524090958974, 2.525120413430, 2.526136561300, 2.527135940871, 2.528144805429, 2.529160718978, 2.530167826066, 2.531166649998, 2.532183453500, 2.533185002131, 2.534179950853, 2.535209352730, 2.536228212502, 2.537245128753, 2.538256182625, 2.539259824238, 2.540272571176, 2.541264577854, 2.542268392131, 2.543298505274, 2.544312512364, 2.545324928985, 2.546338641643, 2.547356717060, 2.548365977938, 2.549394362172, 2.550396807292, 2.551419504622, 2.552433148799, 2.553426487913, 2.554430043520, 2.555456832594, 2.556474793502, 2.557475863044, 2.558516958654, 2.559513774082, 2.560540668571, 2.561545943001, 2.562564970031, 2.563574787588, 2.564572138685, 2.565573702087, 2.566590387940, 2.567605126780, 2.568608892234, 2.569614015748, 2.570643770444, 2.571663014421, 2.572658192825, 2.573673068253, 2.574663082693, 2.575672196918, 2.576688249406, 2.577704877886, 2.578720763797, 2.579728470774, 2.580723636094, 2.581751921302, 2.582773008904, 2.583803498626, 2.584818574605, 2.585815612712, 2.586816119134, 2.587817927114, 2.588799978091, 2.589806547006, 2.590820194079, 2.591813137221, 2.592841518007, 2.593842337429, 2.594857599472, 2.595850237035, 2.596854932438, 2.597873656822, 2.598872357705, 2.599874397562, 2.600909940145, 2.601909751427, 2.602917961842, 2.603913165445, 2.604942131104, 2.605973891044, 2.607010567913, 2.608027709563, 2.609034881713, 2.610041740918, 2.611058388663, 2.612088977173, 2.613109186936, 2.614099828164, 2.615106697243, 2.616113035171, 2.617129084217, 2.618156169139, 2.619184785302, 2.620191028843, 2.621210138955, 2.622238925257, 2.623249177481, 2.624248804531, 2.625261182918, 2.626280702452, 2.627285130489, 2.628310525259, 2.629327248174, 2.630373798306, 2.631391279880, 2.632391031719, 2.633416990390, 2.634398402061, 2.635408484891, 2.636416409982, 2.637426302829, 2.638415693700, 2.639407722342, 2.640409802095, 2.641429034813, 2.642431410024, 2.643479481969, 2.644518212499, 2.645531975567, 2.646549072859, 2.647590742858, 2.648598562786, 2.649594772255, 2.650622994541, 2.651616659435, 2.652620214460, 2.653648003691, 2.654657445649, 2.655648406481, 2.656657786435, 2.657655104337, 2.658660258911, 2.659676473153, 2.660699047429, 2.661735594686, 2.662724079626, 2.663717222162, 2.664721470882, 2.665724426557, 2.666755712397, 2.667767828449, 2.668779067688, 2.669792464020, 2.670804363842, 2.671794149107, 2.672806844462, 2.673823136341, 2.674832363167, 2.675840852440, 2.676838069135, 2.677841510750, 2.678859093336, 2.679855376263, 2.680854574773, 2.681845639772, 2.682862615164, 2.683867925701, 2.684871995136, 2.685870174127, 2.686886069102, 2.687902440861, 2.688930957212, 2.689910659960, 2.690927748258, 2.691948933416, 2.692950680715, 2.693962686267, 2.694957906722, 2.695970939773, 2.696977478972, 2.697987006357, 2.699009309765, 2.700029671718, 2.701050254724, 2.702048528833, 2.703059843630, 2.704068465517, 2.705079435320, 2.706070249305, 2.707079257920, 2.708077754609, 2.709081219395, 2.710115970160, 2.711119471435, 2.712131787808, 2.713157462756, 2.714171847715, 2.715177337182, 2.716195552785, 2.717189892126, 2.718232363893, 2.719246626921, 2.720270334125, 2.721279086651, 2.722292249842, 2.723273104351, 2.724280117022, 2.725301467402, 2.726321525182, 2.727291367271, 2.728282196647, 2.729281345650, 2.730284899065, 2.731300133897, 2.732306961716, 2.733336106054, 2.734347199324, 2.735360415916, 2.736406534595, 2.737416506236, 2.738434063734, 2.739462353759, 2.740484960601, 2.741468311075, 2.742484136082, 2.743508598215, 2.744524630251, 2.745525882506, 2.746507159458, 2.747527082333, 2.748552083585, 2.749560723251, 2.750574889993, 2.751579665046, 2.752586770120, 2.753573807650, 2.754569263646, 2.755570717193, 2.756589610391, 2.757588283129, 2.758591499544, 2.759595790355, 2.760594651031, 2.761607101859, 2.762636501304, 2.763647428505, 2.764676123667, 2.765694347406, 2.766732223076, 2.767771058465, 2.768795042803, 2.769828349179, 2.770840033965, 2.771855108269, 2.772867926969, 2.773887241661, 2.774903781633, 2.775947853871, 2.776964816464, 2.777991719895, 2.779010613390, 2.780038183948, 2.781051927251, 2.782076456479, 2.783103671921, 2.784125397096, 2.785153239011, 2.786180864885, 2.787204808533, 2.788218904244, 2.789207573817, 2.790219664336, 2.791222303110, 2.792218379784, 2.793226728038, 2.794260949185, 2.795254812091, 2.796254757965, 2.797231416415, 2.798226104327, 2.799213775684, 2.800222889236, 2.801210445946, 2.802212371758, 2.803232902128, 2.804252515625, 2.805294777002, 2.806283661139, 2.807304061335, 2.808329378608, 2.809342563063, 2.810359800563, 2.811380551495, 2.812402579379, 2.813423909212, 2.814427250135, 2.815426667921, 2.816439205300, 2.817465519366, 2.818434792589, 2.819454662458, 2.820484976062, 2.821483476266, 2.822500149722, 2.823527019136, 2.824544723977, 2.825541569020, 2.826563137286, 2.827598502051, 2.828614974569, 2.829618870041, 2.830648909508, 2.831683166335, 2.832674097622, 2.833698091520, 2.834713819759, 2.835696823206, 2.836723802720, 2.837755009949, 2.838772791899, 2.839794766551, 2.840817947669, 2.841818803368, 2.842839209037, 2.843845345371, 2.844885114205, 2.845874992541, 2.846878121272, 2.847882960461, 2.848887063255, 2.849883349475, 2.850913351961, 2.851950435667, 2.852953472808, 2.853955418922, 2.854977097631, 2.855991833561, 2.856997697609, 2.858000259580, 2.858985680026, 2.859991903242, 2.860995733070, 2.862056256340, 2.863053785145, 2.864091402784, 2.865085664303, 2.866081888225, 2.867087118270, 2.868082500189, 2.869085951539, 2.870083354459, 2.871109519667, 2.872127115196, 2.873153262177, 2.874151934801, 2.875171153797, 2.876172850014, 2.877176534664, 2.878173358404, 2.879194506815, 2.880204877632, 2.881236105694, 2.882267801197, 2.883272077991, 2.884319607480, 2.885334315453, 2.886344379563, 2.887349761813, 2.888363522194, 2.889363159524, 2.890418078637, 2.891416035488, 2.892399680090, 2.893415117596, 2.894430891427, 2.895447339750, 2.896462066510, 2.897474368092, 2.898490753646, 2.899502631839, 2.900542432257, 2.901563262226, 2.902562205728, 2.903601365693, 2.904609197411, 2.905634400617, 2.906622096730, 2.907621874367, 2.908640849978, 2.909649523797, 2.910649232477, 2.911609436124, 2.912657000862, 2.913658324638, 2.914641624308, 2.915641817670, 2.916635000224, 2.917655249266, 2.918671418398, 2.919650987461, 2.920606725246, 2.921618229297, 2.922631731290, 2.923636675823, 2.924642490710, 2.925644785109, 2.926667738595, 2.927667002333, 2.928692524580, 2.929719735366, 2.930723462482, 2.931735081583, 2.932714840692, 2.933717692669, 2.934747903448, 2.935771573605, 2.936779639878, 2.937806986873, 2.938811495338, 2.939800184332, 2.940776350674, 2.941757374607, 2.942775640722, 2.943782181853, 2.944818981004, 2.945804206401, 2.946803967029, 2.947851480390, 2.948850565880, 2.949875171673, 2.950877765336, 2.951880735171, 2.952894988602, 2.953902242757, 2.954910272588, 2.955893181965, 2.956905848927, 2.957889348671, 2.958907871638, 2.959913344535, 2.960931073863, 2.961937666558, 2.962933438253, 2.963943089315, 2.964956294973, 2.965950989067, 2.966942734543, 2.967938363481, 2.968938301746, 2.969952300086, 2.970925208935, 2.971919028488, 2.972895542139, 2.973911471482, 2.974915845389, 2.975911043180, 2.976913056625, 2.977903766289, 2.978911220714, 2.979952119173, 2.980990114345, 2.981981012857, 2.983017609773, 2.984047476826, 2.985046222952, 2.986063252895, 2.987069180162, 2.988060119250, 2.989035538567, 2.990006787317, 2.991005310441, 2.992046215626, 2.993074661519, 2.994071703212, 2.995080485912, 2.996045566529, 2.997063694958, 2.998026280377, 2.999006169341, 3.000009988888, 3.001017004644, 3.002005417713, 3.003024073610, 3.003996031047, 3.005025523091, 3.006042488187, 3.007032266514, 3.008068104901, 3.009093558556, 3.010027652526, 3.011040820562, 3.012024656007, 3.013044736747, 3.014068564897, 3.015084920085, 3.016116729883, 3.017112144839, 3.018129769651, 3.019152961650, 3.020155368182, 3.021199765693, 3.022195943396, 3.023213194673, 3.024194720994, 3.025209766529, 3.026217963519, 3.027220183420, 3.028256701171, 3.029257603068, 3.030260817035, 3.031270554223, 3.032266271583, 3.033299912360, 3.034323798659, 3.035313358515, 3.036315565930, 3.037357481003, 3.038357306170, 3.039369899392, 3.040390578514, 3.041429905790, 3.042413302913, 3.043379254924, 3.044350246104, 3.045358605424, 3.046381875590, 3.047396906078, 3.048400720367, 3.049421945740, 3.050414848673, 3.051428604325, 3.052427580919, 3.053447033576, 3.054440330907, 3.055453175449, 3.056430799779, 3.057429465738, 3.058438879764, 3.059414790415, 3.060413861198, 3.061407731689, 3.062419927081, 3.063471178113, 3.064468047604, 3.065509123726, 3.066508662804, 3.067532830967, 3.068517211490, 3.069514531478, 3.070522831637, 3.071512484408, 3.072514148466, 3.073486236384, 3.074515153492, 3.075540829198, 3.076566860876, 3.077571956531, 3.078612692598, 3.079581850267, 3.080614347692, 3.081617335286, 3.082619492821, 3.083671883473, 3.084676159778, 3.085711860526, 3.086682162024, 3.087684928295, 3.088714519687, 3.089744421812, 3.090756968092, 3.091760614613, 3.092759057993, 3.093713994211, 3.094678595930, 3.095643179722, 3.096587126546, 3.097581516853, 3.098628328179, 3.099645436020, 3.100634267982, 3.101675848544, 3.102718282629, 3.103757709893, 3.104800736502, 3.105822446161, 3.106787137981, 3.107759543401, 3.108728552491, 3.109739425911, 3.110738646726, 3.111792977291, 3.112768785755, 3.113812258590, 3.114798276758, 3.115750253203, 3.116723639642, 3.117761860684, 3.118735775347, 3.119726182479, 3.120718853394, 3.121704027395, 3.122708147663, 3.123731918139, 3.124749425481, 3.125764681020, 3.126775337281, 3.127715498149, 3.128721947760, 3.129701461958, 3.130688471319, 3.131710076454, 3.132725836832, 3.133717387956, 3.134712392578, 3.135692468011, 3.136697962717, 3.137701616462, 3.138682491848, 3.139717706501, 3.140709155781, 3.141696855174, 3.142689218836, 3.143695342797, 3.144707439190, 3.145743160780, 3.146788665497, 3.147779928861, 3.148814447584, 3.149790729231, 3.150782730673, 3.151806570886, 3.152817393966, 3.153816340948, 3.154803324298, 3.155833589550, 3.156846984884, 3.157843386258, 3.158818288733, 3.159858134205, 3.160850783326, 3.161875967721, 3.162857447417, 3.163836083667, 3.164888027730, 3.165926616404, 3.166877122523, 3.167847610536, 3.168836288035, 3.169836853213, 3.170807549503, 3.171797835817, 3.172807841462, 3.173818257348, 3.174816739564, 3.175791482202, 3.176786686717, 3.177779599147, 3.178791174138, 3.179841906425, 3.180831957613, 3.181837473122, 3.182820180040, 3.183840925373, 3.184825522566, 3.185827679620, 3.186835493336, 3.187836950206, 3.188801814192, 3.189782944777, 3.190833011225, 3.191845765986, 3.192856825457, 3.193839705431, 3.194818693501, 3.195842159968, 3.196857793709, 3.197832658332, 3.198826190714, 3.199842642546, 3.200868375757, 3.201846073692, 3.202850919273, 3.203836567715, 3.204832810285, 3.205797161647, 3.206780438177, 3.207749128798, 3.208795836723, 3.209725394308, 3.210762770078, 3.211766560293, 3.212770549399, 3.213751287800, 3.214753472084, 3.215720861287, 3.216740482487, 3.217768956610, 3.218758904682, 3.219818831348, 3.220834409823, 3.221808217356, 3.222809602555, 3.223837297215, 3.224836089134, 3.225809425181, 3.226836197446, 3.227877879839, 3.228944874376, 3.229955495821, 3.230957386887, 3.231948259026, 3.232962188975, 3.234007521305, 3.234997924245, 3.236007042461, 3.237077729696, 3.238070658833, 3.239045530203, 3.240059579461, 3.241046496013, 3.242091775502, 3.243140337082, 3.244093906725, 3.245051101796, 3.246000463269, 3.247048544170, 3.248008436953, 3.249022080468, 3.250012608077, 3.251088235218, 3.252064868042, 3.253059256898, 3.254046573379, 3.255000200990, 3.256009954877, 3.257056598299, 3.258049907695, 3.259122779055, 3.260162729720, 3.261198044322, 3.262254898090, 3.263181361526, 3.264171240534, 3.265173777626, 3.266193063571, 3.267221175645, 3.268205819406, 3.269213688158, 3.270219855780, 3.271294865730, 3.272248153081, 3.273223905040, 3.274254933617, 3.275266314025, 3.276302208564, 3.277311794049, 3.278310542651, 3.279285981299, 3.280290113625, 3.281318983710, 3.282326169083, 3.283301506652, 3.284254804217, 3.285299822667, 3.286287746113, 3.287300643509, 3.288298195036, 3.289268453799, 3.290260372874, 3.291271548119, 3.292293596021, 3.293279657978, 3.294307303462, 3.295337386372, 3.296351872981, 3.297406633652, 3.298387120022, 3.299434725471, 3.300367767436, 3.301394087015, 3.302407153186, 3.303395513870, 3.304379126419, 3.305326372176, 3.306324923353, 3.307331063406, 3.308396076734, 3.309435369243, 3.310441648649, 3.311426244327, 3.312444286356, 3.313462038968, 3.314397078725, 3.315366453426, 3.316361390832, 3.317351397921, 3.318332820108, 3.319348174301, 3.320375897489, 3.321380570479, 3.322420419643, 3.323396001428, 3.324305959951, 3.325297731754, 3.326311106849, 3.327285327710, 3.328308906052, 3.329289477788, 3.330279701282, 3.331236802458, 3.332210949210, 3.333215347181, 3.334167697430, 3.335170058857, 3.336124826536, 3.337110011847, 3.338055816000, 3.339082371860, 3.340059092413, 3.341056110645, 3.342070697531, 3.343091487883, 3.344080154194, 3.345079728045, 3.346115332738, 3.347120573027, 3.348132017844, 3.349138060991, 3.350218417831, 3.351190312337, 3.352199564504, 3.353215085463, 3.354179973340, 3.355159800027, 3.356131981329, 3.357119191917, 3.358107661096, 3.359072571921, 3.360102320998, 3.361035774881, 3.361970239975, 3.362981846011, 3.364002842700, 3.365008128852, 3.365961276242, 3.366944825752, 3.367979246304, 3.368994805670, 3.369994419706, 3.371010626031, 3.371984212524, 3.372955886302, 3.373961587637, 3.375009787449, 3.375982068544, 3.376950324124, 3.377950812184, 3.378999341691, 3.380002483106, 3.380925461368, 3.381888071923, 3.382887428190, 3.383845993507, 3.384791931960, 3.385748380628, 3.386810638050, 3.387812913579, 3.388763289363, 3.389766890159, 3.390788835299, 3.391769301491, 3.392738039461, 3.393702491919, 3.394703575175, 3.395749099932, 3.396710525786, 3.397714234238, 3.398764865282, 3.399790783962, 3.400778694314, 3.401785287607, 3.402809590799, 3.403784591507, 3.404782742125, 3.405717971188, 3.406646356698, 3.407663327279, 3.408678233162, 3.409679899256, 3.410682762888, 3.411691309654, 3.412643797556, 3.413713205042, 3.414682554903, 3.415699313037, 3.416719590878, 3.417767272531, 3.418832297614, 3.419801735892, 3.420806530330, 3.421802184117, 3.422812772645, 3.423846462283, 3.424772885538, 3.425771896720, 3.426750006537, 3.427731487035, 3.428703536146, 3.429759537124, 3.430699844767, 3.431683261815, 3.432713605791, 3.433714566954, 3.434623315855, 3.435654745774, 3.436669638569, 3.437651216293, 3.438688678560, 3.439655711874, 3.440623705398, 3.441643086395, 3.442638389344, 3.443638390951, 3.444661253891, 3.445657446698, 3.446626779956, 3.447685932991, 3.448652488098, 3.449634651704, 3.450616590022, 3.451606896215, 3.452609316693, 3.453617758996, 3.454614940237, 3.455636735765, 3.456699473429, 3.457697534434, 3.458669173645, 3.459671776780, 3.460659137606, 3.461652520161, 3.462662042336, 3.463716866175, 3.464714743200, 3.465713649218, 3.466735211124, 3.467709455683, 3.468650113131, 3.469590251436, 3.470572212393, 3.471549967214, 3.472526060756, 3.473513397688, 3.474448595584, 3.475496129105, 3.476451221172, 3.477473610444, 3.478465737934, 3.479458827171, 3.480431873405, 3.481438685126, 3.482387168338, 3.483329797607, 3.484322162657, 3.485296885363, 3.486336344871, 3.487332943896, 3.488351889891, 3.489410758349, 3.490308331590, 3.491269676582, 3.492253390038, 3.493229871827, 3.494186871984, 3.495160924927, 3.496149418386, 3.497166090225, 3.498108566975, 3.499120246728, 3.500124672286, 3.501124541584, 3.502099117026, 3.503026094872, 3.503991091049, 3.505005469161, 3.505965131859, 3.506942269940, 3.507889443957, 3.508844294650, 3.509871490329, 3.510847610395, 3.511786415557, 3.512814943911, 3.513763686383, 3.514750024639, 3.515735760021, 3.516705183660, 3.517718259993, 3.518742308519, 3.519768777417, 3.520796237431, 3.521743825593, 3.522739793277, 3.523661179750, 3.524657198516, 3.525687560547, 3.526656115215, 3.527626834763, 3.528624665960, 3.529592449416, 3.530591865069, 3.531580293205, 3.532581338433, 3.533565407705, 3.534584428967, 3.535552184995, 3.536487745184, 3.537443288984, 3.538370934372, 3.539381758556, 3.540399462643, 3.541405960279, 3.542410253199, 3.543436605140, 3.544364986120, 3.545365480880, 3.546363701542, 3.547388728361, 3.548360911979, 3.549289121853, 3.550324168987, 3.551278225050, 3.552279301533, 3.553296664757, 3.554316416818, 3.555318289514, 3.556280265459, 3.557260045770, 3.558290728422, 3.559287653830, 3.560323165476, 3.561261511483, 3.562304916296, 3.563328587949, 3.564322824993, 3.565306573490, 3.566314953753, 3.567344925049, 3.568308226132, 3.569328441363, 3.570284854571, 3.571248233239, 3.572210510092, 3.573166796916, 3.574185469641, 3.575180409870, 3.576177634676, 3.577209964700, 3.578142808544, 3.579099080237, 3.580054159482, 3.581008033751, 3.581925860991, 3.582848956057, 3.583838978752, 3.584832933087, 3.585852596441, 3.586889757097, 3.587855418013, 3.588877155173, 3.589872586028, 3.590878768755, 3.591873712574, 3.592891351160, 3.593861940960, 3.594841538873, 3.595838763332, 3.596889778420, 3.597840108114, 3.598785625321, 3.599791957964, 3.600781574657, 3.601783867835, 3.602744979344, 3.603741355698, 3.604682345338, 3.605716465713, 3.606686329514, 3.607725238428, 3.608704895500, 3.609713287551, 3.610660232577, 3.611612798466, 3.612498054724, 3.613504601564, 3.614479521440, 3.615467384506, 3.616502395377, 3.617529077630, 3.618590673875, 3.619511983769, 3.620438876807, 3.621344142971, 3.622304078921, 3.623293503827, 3.624223027117, 3.625280974438, 3.626181709499, 3.627128483996, 3.628127129403, 3.629129925354, 3.630120216897, 3.631075625999, 3.632116908335, 3.633050604791, 3.633986312940, 3.634965265615, 3.635940795536, 3.636961818118, 3.637973924983, 3.638994070132, 3.640016617216, 3.640963671592, 3.641952784317, 3.642999509559, 3.644044936664, 3.645037264278, 3.645985732255, 3.647030674385, 3.648022131055, 3.649033276444, 3.649961424143, 3.650953769206, 3.651926952367, 3.652904274105, 3.653919033836, 3.654942056509, 3.655902596406, 3.656853440632, 3.657871554636, 3.658818805301, 3.659770110608, 3.660751343041, 3.661702909767, 3.662670547263, 3.663718423467, 3.664702607298, 3.665757415552, 3.666695834500, 3.667727209575, 3.668710403744, 3.669667410777, 3.670669253312, 3.671679529883, 3.672700338082, 3.673705112182, 3.674734805478, 3.675657859041, 3.676752026161, 3.677584283140, 3.678476132928, 3.679515137764, 3.680406804362, 3.681352430996, 3.682327288099, 3.683268732659, 3.684226915187, 3.685237706732, 3.686248747523, 3.687143795597, 3.688210125773, 3.689228118256, 3.690188914229, 3.691128381086, 3.692106219283, 3.693154821584, 3.694203814083, 3.695173540977, 3.696115235515, 3.697106541643, 3.698156464638, 3.699165480838, 3.700113703448, 3.701050909294, 3.702062306646, 3.703027841343, 3.703951594729, 3.704930149419, 3.705897675440, 3.706935918084, 3.707945612190, 3.708982101783, 3.709967615794, 3.711024576013, 3.711992364805, 3.712978013346, 3.713916456009, 3.714962805555, 3.715995873834, 3.717049515208, 3.718024032275, 3.719039400227, 3.720038911949, 3.721043014000, 3.722012804637, 3.723051323703, 3.724170565956, 3.725115085412, 3.726068597311, 3.726933878136, 3.727926212500, 3.728941758582, 3.729913045810, 3.730926242089, 3.731911352960, 3.732901051365, 3.733904776972, 3.734906110109, 3.735926307402, 3.736813842586, 3.737826662754, 3.738877555959, 3.739849873450, 3.740793288723, 3.741726776998, 3.742592643926, 3.743647927092, 3.744696130317, 3.745717849549, 3.746659571162, 3.747649488453, 3.748658709757, 3.749611716986, 3.750574156357, 3.751604912282, 3.752611093192, 3.753486647534, 3.754504605605, 3.755490327903, 3.756406409848, 3.757376590207, 3.758289195932, 3.759293535623, 3.760322711384, 3.761316729831, 3.762310516177, 3.763336801297, 3.764284745649, 3.765219587438, 3.766197004285, 3.767156300377, 3.768125358947, 3.769106792777, 3.770103238940, 3.771112231981, 3.772036208533, 3.772998204012, 3.773952012222, 3.774962237184, 3.775977410288, 3.776865039460, 3.777934157411, 3.778891049339, 3.779857902278, 3.780819046947, 3.781771813010, 3.782710873971, 3.783702113370, 3.784655942219, 3.785715268630, 3.786793133791, 3.787697842222, 3.788655157902, 3.789646695194, 3.790616366004, 3.791663467379, 3.792648433412, 3.793673444932, 3.794641331887, 3.795698194596, 3.796749476057, 3.797778771697, 3.798802314530, 3.799833753637, 3.800864902536, 3.801870984220, 3.802832511029, 3.803771293660, 3.804745352533, 3.805827122808, 3.806844779396, 3.807739281358, 3.808618858350, 3.809637482025, 3.810588298751, 3.811667850583, 3.812631590947, 3.813654025669, 3.814571178779, 3.815643644482, 3.816570706962, 3.817708063079, 3.818728228441, 3.819759399044, 3.820735535661, 3.821736917442, 3.822726176155, 3.823656923165, 3.824572270477, 3.825544765854, 3.826595181982, 3.827560543713, 3.828607072167, 3.829556388527, 3.830522482907, 3.831531981146, 3.832582227803, 3.833578767525, 3.834497489391, 3.835477623980, 3.836454015239, 3.837486373655, 3.838431370111, 3.839318425717, 3.840282453820, 3.841245613368, 3.842286418244, 3.843202578310, 3.844169209487, 3.845128876017, 3.846072386111, 3.846972149914, 3.847965578017, 3.848906077215, 3.849977723662, 3.850891785200, 3.851881878408, 3.852855664451, 3.853788215604, 3.854772507348, 3.855709188045, 3.856829006800, 3.857839005330, 3.858826255921, 3.859784308501, 3.860769692614, 3.861713094647, 3.862674378269, 3.863602896609, 3.864523866876, 3.865520076023, 3.866566484696, 3.867538583342, 3.868560992384, 3.869502195683, 3.870438997691, 3.871481187500, 3.872558267995, 3.873521166225, 3.874495962901, 3.875521859702, 3.876439074425, 3.877443375239, 3.878400765809, 3.879376719952, 3.880315308176, 3.881255929251, 3.882291314964, 3.883156579082, 3.884189858590, 3.885222268028, 3.886170247260, 3.887204031310, 3.888126139776, 3.889208339229, 3.890259508661, 3.891211798140, 3.892233946712, 3.893238127626, 3.894186766408, 3.895151127332, 3.896192859416, 3.897233668454, 3.898156732234, 3.899123073319, 3.900115723218, 3.901172905668, 3.902104389145, 3.903110833627, 3.904140512006, 3.905193584322, 3.906189725128, 3.907132042918, 3.908090467901, 3.909100328687, 3.910080765392, 3.911112965562, 3.912105055989, 3.913159863257, 3.914085371293, 3.915109285186, 3.916132038525, 3.917150028698, 3.918170410656, 3.919182379001, 3.919907691821, 3.920902001706, 3.921938504542, 3.922882931318, 3.923858571695, 3.924891203583, 3.925955593104, 3.926905135484, 3.927816299647, 3.928692524580, 3.929655463463, 3.930509504433, 3.931550671962, 3.932549720738, 3.933584615461, 3.934569681380, 3.935489599085, 3.936430226995, 3.937421775913, 3.938411825089, 3.939445689364, 3.940459300800, 3.941444921185, 3.942489844171, 3.943480088883, 3.944373243728, 3.945421414627, 3.946379987054, 3.947286824910, 3.948218688516, 3.949164145397, 3.950092307210, 3.951014697429, 3.951981820688, 3.953017357242, 3.953871780231, 3.954763106128, 3.955789613000, 3.956759577343, 3.957755354601, 3.958702080340, 3.959714201067, 3.960716783818, 3.961618308532, 3.962629286585, 3.963558754001, 3.964538241283, 3.965568083035, 3.966532013068, 3.967437644693, 3.968446121943, 3.969485283409, 3.970462014272, 3.971554824866, 3.972519919587, 3.973523937012, 3.974534376530, 3.975481494640, 3.976348419889, 3.977270655283, 3.978170072281, 3.979174833678, 3.980144583675, 3.981025034099, 3.982003095560, 3.983016774494, 3.983923998122, 3.984778603673, 3.985802992505, 3.986753971985, 3.987639495839, 3.988577588512, 3.989606747927, 3.990549104201, 3.991502027199, 3.992491192360, 3.993375672354, 3.994261957349, 3.995188708882, 3.996203536641, 3.997074047915, 3.998210057482, 3.999162619462, 4.000212856450, 4.001191600119, 4.002133273761, 4.003146980479, 4.004070988347, 4.005010146141, 4.006030598111, 4.006863690551, 4.007901780068, 4.008853700008, 4.009749967072, 4.010746014975, 4.011677427847, 4.012668968191, 4.013559709617, 4.014560052667, 4.015427684378, 4.016342145369, 4.017403168183, 4.018353514423, 4.019337728503, 4.020337831379, 4.021299189476, 4.022235252972, 4.023200825089, 4.024131818459, 4.025069412870, 4.026004424869, 4.026913729778, 4.027774012212, 4.028821600059, 4.029890330064, 4.030845115548, 4.031890798872, 4.032798474342, 4.033839487365, 4.034798298974, 4.035721508346, 4.036755380332, 4.037654345155, 4.038531444649, 4.039495926246, 4.040467321578, 4.041531681633, 4.042440598103, 4.043389813164, 4.044418070757, 4.045453590776, 4.046370761128, 4.047328614652, 4.048283731758, 4.049172858393, 4.050219791638, 4.051166637247, 4.052012734006, 4.052934064570, 4.054083564254, 4.055127624746, 4.056105009876, 4.056955840502, 4.057917502876, 4.058791784901, 4.059662847975, 4.060545646808, 4.061565340829, 4.062572386426, 4.063601883711, 4.064487695985, 4.065562152215, 4.066593707019, 4.067500868395, 4.068526873333, 4.069596088171, 4.070432928498, 4.071317444156, 4.072342255812, 4.073266658483, 4.074265167195, 4.075142040326, 4.076067253711, 4.076942591878, 4.077871653204, 4.079005826509, 4.080007190459, 4.081220260509, 4.082200509064, 4.083072533887, 4.084041166070, 4.084885214596, 4.085857902695, 4.086885819108, 4.087836420346, 4.088783778782, 4.089775928046, 4.090706124121, 4.091493535010, 4.092352206600, 4.093357935411, 4.094371395680, 4.095295271655, 4.096242798032, 4.097127214554, 4.098095080860, 4.099070564688, 4.100004502312, 4.100945932002, 4.101889406883, 4.102840439147, 4.103799074030, 4.104748774375, 4.105589776527, 4.106532304368, 4.107493569663, 4.108445817573, 4.109344290849, 4.110457400927, 4.111505995270, 4.112399576721, 4.113232993510, 4.114282668818, 4.115380201245, 4.116412372794, 4.117321807041, 4.118352905507, 4.119403609265, 4.120491249408, 4.121472459901, 4.122565300140, 4.123585847236, 4.124510443976, 4.125431216146, 4.126382993083, 4.127226244081, 4.128298663026, 4.129157339524, 4.130205236350, 4.131208666802, 4.132184979935, 4.133222509833, 4.134238859956, 4.135180513266, 4.136225232878, 4.137272471682, 4.138262527785, 4.139200987505, 4.140129485913, 4.141084012649, 4.142004504439, 4.143011456199, 4.144020748097, 4.144911111909, 4.145894448076, 4.146928743793, 4.147947191302, 4.148925194771, 4.149887006293, 4.150801782571, 4.151700013464, 4.152630964124, 4.153551544585, 4.154430696045, 4.155404791724, 4.156381077151, 4.157365801654, 4.158252716325, 4.159235432140, 4.160214096070, 4.161025045045, 4.161894251561, 4.162885251485, 4.163884851730, 4.164835973004, 4.165700127111, 4.166591497253, 4.167650783740, 4.168501357447, 4.169283050563, 4.170297502735, 4.171230577126, 4.172197940721, 4.173044537823, 4.174028921800, 4.174983051787, 4.175958819435, 4.176995532349, 4.177910414680, 4.178892789177, 4.179857676921, 4.180752275325, 4.181721307835, 4.182798347337, 4.183831641084, 4.184860751809, 4.185912296697, 4.186766051944, 4.187855691030, 4.188726692288, 4.189639767114, 4.190507623063, 4.191390712520, 4.192390566478, 4.193453744550, 4.194451569152, 4.195465314614, 4.196351727599, 4.197274152329, 4.198321940964, 4.199303537451, 4.200280469974, 4.201218188432, 4.202240952582, 4.203252261125, 4.204238126823, 4.205052109623, 4.206021136101, 4.207027305088, 4.208021787501, 4.209018552387, 4.209848593923, 4.210750775362, 4.211739686991, 4.212900997543, 4.213809548887, 4.214663045466, 4.215753691933, 4.216754074544, 4.217613383257, 4.218697101510, 4.219711497736, 4.220692172416, 4.221624422080, 4.222522428299, 4.223494946429, 4.224731880349, 4.225753183752, 4.226732971705, 4.227729648701, 4.228794807681, 4.229678298489, 4.230497130741, 4.231406292280, 4.232324775972, 4.233267498932, 4.234115469108, 4.235017322617, 4.235973392308, 4.237043996484, 4.237936935218, 4.238794080818, 4.239894305146, 4.240846059675, 4.241868116559, 4.242793826721, 4.243630178398, 4.244658816188, 4.245491127342, 4.246363327814, 4.247367761003, 4.248412994918, 4.249391335583, 4.250217320947, 4.251068105586, 4.251982622206, 4.252906843400, 4.253833035623, 4.254644101548, 4.255472326527, 4.256521593406, 4.257424109437, 4.258359996064, 4.259132243454, 4.259961178031, 4.260847122311, 4.261600012139, 4.262560849815, 4.263523817972, 4.264488926080, 4.265416171650, 4.266337382081, 4.267140026063, 4.268161526857, 4.269072445649, 4.270049977240, 4.271013502097, 4.272003541706, 4.273158732355, 4.273974094719, 4.274872764935, 4.276002823808, 4.277004316111, 4.278139942310, 4.279171158619, 4.280113767528, 4.281099906259, 4.282030086795, 4.282978928103, 4.283896446124, 4.284799171805, 4.285645088793, 4.286534658478, 4.287451308835, 4.288386771041, 4.289256618619, 4.290246816052, 4.291077953222, 4.292021268386, 4.292915483301, 4.293965339286, 4.294906389452, 4.295772243975, 4.296493677813, 4.297414386073, 4.298466552658, 4.299382804318, 4.300353024173, 4.301473202109, 4.302369686853, 4.303425208204, 4.304422018880, 4.305263216614, 4.306132411670, 4.306976932076, 4.308008416669, 4.308936198440, 4.309919156244, 4.310833264573, 4.311767107239, 4.312613746248, 4.313479915583, 4.314500100906, 4.315387999730, 4.316142794600, 4.317025049815, 4.317809781016, 4.318867349819, 4.319818647684, 4.320753854586, 4.321672861889, 4.322648589789, 4.323589917198, 4.324698361060, 4.325625762283, 4.326472250810, 4.327449602823, 4.328540191795, 4.329401657739, 4.330394925378, 4.331614083310, 4.332612429695, 4.333491379341, 4.334278331613, 4.335132475500, 4.335903585148, 4.336808090839, 4.337799555729, 4.338850142755, 4.339979277732, 4.340892320800, 4.341969513895, 4.342934581291, 4.343892211018, 4.344900000029, 4.345890870480, 4.346777835343, 4.347734314462, 4.348867421450, 4.350110455631, 4.350976854549, 4.351971936350, 4.352949725131, 4.353674708757, 4.354793949140, 4.355551562696, 4.356320365257, 4.357199254307, 4.358238438071, 4.359240385651, 4.359995853469, 4.360912128916, 4.362020214600, 4.362970826728, 4.363803068068, 4.364707304349, 4.365704143610, 4.366602249277, 4.367573079121, 4.368606969287, 4.369806117112, 4.370620986093, 4.371621288764, 4.372541968040, 4.373454340973, 4.374533206391, 4.375449776582, 4.376430277153, 4.377340509568, 4.378429084411, 4.379385131226, 4.380134812608, 4.381146850513, 4.382046097848, 4.383094080987, 4.384060464076, 4.384955233240, 4.386010268454, 4.386961996497, 4.387926424529, 4.388637869334, 4.389765824666, 4.390693793285, 4.391591648358, 4.392544976785, 4.393618634889, 4.394565654652, 4.395568731847, 4.396444271375, 4.397419167736, 4.398352784529, 4.399419125837, 4.400346138757, 4.401307959067, 4.402447408902, 4.403292970319, 4.404393565134, 4.405231886552, 4.406226732959, 4.407112958105, 4.408123242642, 4.408879717763, 4.409882965215, 4.410720778764, 4.411649847046, 4.412580907107, 4.413525221429, 4.414426481377, 4.415205323167, 4.416132723194, 4.417141537775, 4.418152701170, 4.418972447192, 4.419907935299, 4.420719562037, 4.421612960289, 4.422565650323, 4.423451343824, 4.424500406591, 4.425390058660, 4.426304716602, 4.427093562741, 4.427779165635, 4.428943682024, 4.429713976918, 4.430672917112, 4.431610513508, 4.432644212860, 4.433621439141, 4.434470938249, 4.435487800283, 4.436471453321, 4.437397885222, 4.438183335681, 4.439292523449, 4.440177131287, 4.441063544648, 4.442035886446, 4.442950190956, 4.443938840990, 4.444808783489, 4.445983551414, 4.446979138865, 4.447843019358, 4.448696418104, 4.449563724720, 4.450641105658, 4.451610581867, 4.452545289427, 4.453704164879, 4.454767128253, 4.455795484543, 4.456614955190, 4.457759833710, 4.458670422333, 4.459645494548, 4.460522427195, 4.461388568656, 4.462294213919, 4.463428932772, 4.464440066740, 4.465250672325, 4.466113601772, 4.467067354879, 4.467831868306, 4.468930023178, 4.469992510024, 4.470864922499, 4.472112434047, 4.472976214490, 4.473945175227, 4.474747818251, 4.475564929465, 4.476565713254, 4.477725343989, 4.478665743222, 4.479831002476, 4.480684014409, 4.481657179535, 4.482513789420, 4.483715888886, 4.484748914794, 4.485717952140, 4.486555986379, 4.487475692676, 4.488624726586, 4.489656098610, 4.490797477669, 4.491739694488, 4.492751486081, 4.493494967595, 4.494537984646, 4.495556324240, 4.496549806558, 4.497613853126, 4.498557303139, 4.499543962742, 4.500519118770, 4.501634301206, 4.502697074459, 4.503596249519, 4.504733255612, 4.505706231335, 4.506583777133, 4.507505016657, 4.508316207153, 4.509199047989, 4.510041520575, 4.510772987777, 4.511787826701, 4.512706041994, 4.513782115332, 4.514505714989, 4.515358555118, 4.516369917505, 4.517040708833, 4.517869924608, 4.518901503435, 4.519518755938, 4.520309545798, 4.521188291048, 4.522227791651, 4.523197344734, 4.524125510681, 4.524910196611, 4.525900340699, 4.526921970156, 4.527902071631, 4.529075246700, 4.529883646849, 4.530678815485, 4.531726535475, 4.532687937019, 4.533770208294, 4.534706396313, 4.535465743699, 4.536286138458, 4.537257698166, 4.538186446377, 4.539267491108, 4.540155357612, 4.541029948371, 4.541891179383, 4.542738967535, 4.543618782015, 4.544393887418, 4.545246580034, 4.546284248299, 4.547232523297, 4.548198218652, 4.549089173971, 4.550059011227, 4.551031019106, 4.552005207346, 4.552966070525, 4.553758054675, 4.554707230574, 4.555783394930, 4.556815270850, 4.557786847277, 4.558713434659, 4.559799584639, 4.560777867928, 4.561679205814, 4.562693468676, 4.563408043205, 4.564474148501, 4.565622738750, 4.566774374763, 4.567415494575, 4.568346802653, 4.569328441363, 4.570183145693, 4.570893991667, 4.571784188439, 4.572741159730, 4.573618883821, 4.574563604591, 4.575363335350, 4.576262750018, 4.577246058699, 4.578100064157, 4.578922809983, 4.579730615922, 4.580589525914, 4.581433571626, 4.582445275637, 4.583176495233, 4.583925609078, 4.584943141889, 4.585896110742, 4.587018947242, 4.587875593827, 4.588717086983, 4.589644616566, 4.590489547731, 4.591319178190, 4.592184357616, 4.593153366902, 4.594175717850, 4.595234681598, 4.596176251098, 4.597119866398, 4.598099963239, 4.598961522167, 4.599842076610, 4.600862989661, 4.601834216894, 4.602738019719, 4.603608838270, 4.604638655059, 4.605530807653, 4.606635414815, 4.607355798262, 4.608394549114, 4.609153173105, 4.610072352619, 4.611082154843, 4.612058756293, 4.613019742672, 4.613965005259, 4.615001806738, 4.616023149569, 4.617082864912, 4.618181228374, 4.618812693366, 4.619824951003, 4.620893995264, 4.621965677543, 4.623058242853, 4.623916009638, 4.624812084639, 4.625930196274, 4.626959199970, 4.627990647547, 4.628710426935, 4.629746047039, 4.630821263363, 4.631880541385, 4.633128970016, 4.634100379276, 4.635036480173, 4.636049740603, 4.637479827711, 4.638517704175, 4.639293010751, 4.640240383584, 4.641208837608, 4.642065153000, 4.643133153318, 4.644127225107, 4.645066033441, 4.646372241014, 4.647393045863, 4.648126423379, 4.649248182184, 4.650275777366, 4.651052934074, 4.652006850213, 4.652943334633, 4.654018878174, 4.655057837235, 4.655823944804, 4.656847528338, 4.657597060285, 4.658506122601, 4.659278342687, 4.660051938306, 4.660886585229, 4.661742769754, 4.662800394626, 4.663800520533, 4.664802954922, 4.665928435259, 4.666814791947, 4.667702961311, 4.668552457122, 4.669586226651, 4.670622462777, 4.671498077315, 4.672559323757, 4.673377435484, 4.674484336637, 4.675388346671, 4.676314852898, 4.677222685297, 4.678277325654, 4.679293025690, 4.680373515844, 4.681248186143, 4.682041075299, 4.682940041260, 4.683945741011, 4.684764590382, 4.685690277755, 4.686617942422, 4.687484144337, 4.688542831619, 4.689391643541, 4.690391122041, 4.691222226335, 4.691990813762, 4.693039137117, 4.693617867850, 4.694476666789, 4.695681845510, 4.696695922561, 4.697690721630, 4.698535926857, 4.699730685470, 4.700558079198, 4.701561774689, 4.702545900381, 4.703685896647, 4.704894927190, 4.705798399305, 4.707101824982, 4.708098599535, 4.709075440617, 4.709965388637, 4.710790217355, 4.711773138262, 4.712691048430, 4.713408816563, 4.714330194040, 4.715208444050, 4.716065886543, 4.716766635111, 4.717468516188, 4.718216928494, 4.719285080509, 4.720036632541, 4.720926510426, 4.722047152961, 4.723193654371, 4.724251115521, 4.725126616630, 4.726096333579, 4.726975566182, 4.727949425011, 4.728972005738, 4.729880401381, 4.730860803208, 4.731960550066, 4.732875224890, 4.733791830192, 4.734969802277, 4.736103669395, 4.736979716896, 4.737786294524, 4.738641953806, 4.739642358364, 4.740430010356, 4.741386659197, 4.742585441238, 4.743498733789, 4.744172917297, 4.745258623908, 4.746153353801, 4.747074188334, 4.747851145181, 4.748629494503, 4.749482414390, 4.750385897655, 4.751046384504, 4.752051270743, 4.753058487517, 4.753747687701, 4.754734160543, 4.755500221166, 4.756416324002, 4.757458571702, 4.758378819209, 4.759350925379, 4.760475296843, 4.761276624605, 4.762255246125, 4.762959208621, 4.763739928999, 4.764623077598, 4.765356192238, 4.766369419836, 4.767232525824, 4.768097350543, 4.769091480679, 4.770318157682, 4.771317390287, 4.772061902631, 4.773297010415, 4.774096555074, 4.775001039461, 4.775959261372, 4.776893619071, 4.777803953698, 4.778898880038, 4.779970406008, 4.780834778170, 4.781963668799, 4.783148216716, 4.784018943208, 4.785156151952, 4.786322898118, 4.787067008292, 4.787865686532, 4.788745932382, 4.789601211348, 4.790753151247, 4.791800582866, 4.792715794088, 4.793497942411, 4.794416741359, 4.795554418567, 4.796341700544, 4.797239312607, 4.798193357043, 4.798821443617, 4.799669434875, 4.800656281311, 4.801562864674, 4.802609158867, 4.803740889495, 4.804487789326, 4.805652192720, 4.806541457391, 4.807488300496, 4.808465152993, 4.809416204023, 4.810313217407, 4.811352704000, 4.812056470937, 4.812874271898, 4.813891620187, 4.814911357242, 4.815762970984, 4.816616257947, 4.817556814043, 4.818327877932, 4.819243507696, 4.820074967926, 4.821051814883, 4.821973211730, 4.822925454090, 4.823966650736, 4.824720181558, 4.825649402521, 4.826464104990, 4.827367883314, 4.828185818949, 4.829210409554, 4.830178671138, 4.831207979686, 4.832328283770, 4.832919080327, 4.834162375310, 4.834874430912, 4.836093566542, 4.836957953086, 4.838003474578, 4.838991562694, 4.839831707041, 4.840823968207, 4.841848674607, 4.842875804495, 4.843693200506, 4.844542498548, 4.845667338758, 4.846581514962, 4.847405922073, 4.848323769153, 4.849028767937, 4.850257585287, 4.851243148678, 4.852292766220, 4.853251986369, 4.854461764288, 4.855612625898, 4.856329566530, 4.857078942484, 4.857860919868, 4.858895290672, 4.859774874734, 4.860593229558, 4.861413129346, 4.862392729495, 4.863342838126, 4.864072664995, 4.865090120868, 4.866141874797, 4.866940357246, 4.867932518696, 4.869216085411, 4.870213465480, 4.870987206874, 4.871762329231, 4.872765580837, 4.873608808383, 4.874323589618, 4.875658829225, 4.876540380557, 4.877456475931, 4.878276054363, 4.879228707063, 4.880216413838, 4.881239409558, 4.881966871171, 4.882861333596, 4.883624738820, 4.884389488326, 4.885589197602, 4.886591493182, 4.887394998465, 4.888502251192, 4.889477627399, 4.890623968596, 4.891807194865, 4.892756102421, 4.893809102737, 4.894762399254, 4.895649480757, 4.896298880410, 4.897463393126, 4.898252926054, 4.899388247166, 4.900561084446, 4.901563954660, 4.902083590763, 4.903124731940, 4.904133546521, 4.905109802993, 4.906648367984, 4.907279355316, 4.908122100371, 4.909213072051, 4.910094888561, 4.911013881052, 4.911934822310, 4.913035426123, 4.913640169325, 4.914602377983, 4.915495216754, 4.916354072304, 4.917142851404, 4.917789283399, 4.918724720471, 4.919662176752, 4.920529320765, 4.921543181947, 4.922559415529, 4.923250359376, 4.923869505457, 4.924599044486, 4.925402955410, 4.926208357192, 4.927162125337, 4.927860436760, 4.928817844010, 4.929592678260, 4.930590929328, 4.931480193999, 4.932371283272, 4.933338697699, 4.934233612377, 4.935317433771, 4.936403966709, 4.937305226658, 4.938283706840, 4.939604382680, 4.940436582099, 4.941232444630, 4.942105780275, 4.942828558453, 4.943743264150, 4.944851110111, 4.945770090137, 4.946729433319, 4.947690900353, 4.948654500663, 4.949542905214, 4.950045851978, 4.950898321791, 4.951869072971, 4.952413642926, 4.953231780339, 4.953739019808, 4.954755279522, 4.955695480824, 4.956677026425, 4.957621401860, 4.958765268829, 4.959635044140, 4.960903289559, 4.962095642688, 4.962852496367, 4.963770455914, 4.964410182757, 4.965251256853, 4.966415424485, 4.967623028790, 4.968429967859, 4.969117046890, 4.970210529168, 4.971144190610, 4.971957704909, 4.972813538163, 4.973834591161, 4.974899038953, 4.975678557858, 4.976541762356, 4.977489149566, 4.978273335586, 4.979058940138, 4.979970366457, 4.981049978499, 4.982257335838, 4.983008421794, 4.983551681741, 4.984472595686, 4.985353475316, 4.986572872929, 4.987669044420, 4.988598740075, 4.989233766482, 4.990379159186, 4.991399828238, 4.992294885635, 4.993106292052, 4.993833574515, 4.994562076949, 4.995463682149, 4.996324097451, 4.997272527310, 4.998266287191, 4.999175623394, 5.000347574634, 5.001479117165, 5.002264223402, 5.003313244400, 5.004277078005, 5.005155150447, 5.006123085059, 5.007137328611, 5.007755855953, 5.008640997362, 5.009394788576, 5.010327752376, 5.011619434778, 5.013049212141, 5.013765869865, 5.014798141635, 5.015787833239, 5.016238439714, 5.017050711425, 5.017909760704, 5.018861217359, 5.019632973815, 5.020633757604, 5.021682503253, 5.022230681908, 5.022825323980, 5.023650020997, 5.024798037742, 5.025764722557, 5.026918350399, 5.027843464141, 5.028724151262, 5.029374223312, 5.030397735151, 5.031610458153, 5.032498825277, 5.033295223342, 5.034422042089, 5.035457533921, 5.036070577973, 5.037346595107, 5.038057116859, 5.038768802955, 5.039576792922, 5.040433953362, 5.041483896577, 5.042440598103, 5.043639446327, 5.044793462458, 5.045902279521, 5.046868774817, 5.048128442872, 5.048662481204, 5.049245818406, 5.050122295963, 5.051342567859, 5.052223291535, 5.053105804898, 5.054186873413, 5.054925511521, 5.055862926842, 5.056851870641, 5.057743849581, 5.058737090681, 5.059284343093, 5.060181337179, 5.061630402548, 5.062331685601, 5.063385738025, 5.064089863569, 5.064895978849, 5.065754118977, 5.066917837630, 5.067729224101, 5.068440432909, 5.069356558958, 5.070121476343, 5.071041159119, 5.071911540335, 5.072886388066, 5.073400346093, 5.074223946163, 5.075152367024, 5.075927570090, 5.077170780333, 5.077430231014, 5.078053545771, 5.078886026163, 5.079928875702, 5.080974235413, 5.081864773834, 5.082757142093, 5.083335535459, 5.084441884588, 5.085233863074, 5.085815566577, 5.086875209610, 5.087459117209, 5.088309841246, 5.089642442727, 5.090230085267, 5.091407761152, 5.091997798827, 5.092481153893, 5.093718844228, 5.094258072608, 5.094960071924, 5.095933948085, 5.096692920036, 5.097888276552, 5.098759697927, 5.100015437451, 5.100726812682, 5.101603954070, 5.102702877941, 5.103418672494, 5.104025267641, 5.104964402548, 5.105628546144, 5.106515653782, 5.107460195341, 5.108518296160, 5.109690783100, 5.110194248132, 5.110754339107, 5.111427403142, 5.112213965162, 5.112889296875, 5.113791375833, 5.114865033937, 5.115657852353, 5.116622510252, 5.117589315626, 5.118672415899, 5.119643800581, 5.120617362826, 5.121880515303, 5.122859110215, 5.123551213122, 5.124360062996, 5.125402231297, 5.126214539182, 5.126912014410, 5.127668878770, 5.128543825977, 5.129654628919, 5.130533589919, 5.131531900979, 5.132414673453, 5.133063182268, 5.134126471792, 5.135014539340, 5.136796140971, 5.137391636035, 5.138465589141, 5.139362583226, 5.140261433803, 5.140981856111, 5.141703475466, 5.142546882965, 5.143210711247, 5.143754599210, 5.144723196170, 5.145997766873, 5.146423456380, 5.147642016332, 5.148680487351, 5.149537567238, 5.150457747995, 5.151318345976, 5.152797636019, 5.153662887870, 5.154343940016, 5.154901959986, 5.155895769302, 5.156767221902, 5.157390760389, 5.158390287832, 5.159141445958, 5.160333443118, 5.160899217293, 5.161591721506, 5.162222230446, 5.163422725159, 5.164436248331, 5.165770097148, 5.166534139829, 5.167810538931, 5.168450148004, 5.169668006548, 5.170374646642, 5.171533452647, 5.172114017210, 5.173601217812, 5.174508970121, 5.175418623767, 5.176525770830, 5.177243667049, 5.178290002702, 5.179338865356, 5.180258702727, 5.181048688360, 5.181906130853, 5.182699121607, 5.183957659078, 5.184887341810, 5.185885643871, 5.186552455735, 5.187420844591, 5.188156993824, 5.189297139053, 5.190305641283, 5.191316490871, 5.192194467729, 5.193074223116, 5.193413065672, 5.194431182451, 5.195383583013, 5.196338076764, 5.197363081917, 5.198184831419, 5.198802165541, 5.199970640756, 5.200797343699, 5.201902067938, 5.202732459169, 5.204050501097, 5.205024255949, 5.206419132632, 5.207118254615, 5.207818503850, 5.208519883980, 5.209363038068, 5.209785229756, 5.210701388841, 5.211831628859, 5.213319546803, 5.213532523260, 5.214456630043, 5.215525356237, 5.216954427885, 5.217957583379, 5.218747405752, 5.219754716135, 5.220259248823, 5.221197795987, 5.222065951162, 5.223080997158, 5.224316750974, 5.225264117448, 5.226652458019, 5.227825039175, 5.228486010020, 5.229589868486, 5.230548820598, 5.231361898752, 5.232250632654, 5.232992636050, 5.233884716779, 5.235226283089, 5.236048173967, 5.237021509132, 5.238147305534, 5.238974748289, 5.240332155310, 5.240861183719, 5.242148656314, 5.243135775939, 5.244048958996, 5.244811414392, 5.245422343955, 5.246187216435, 5.247337056879, 5.247874692702, 5.249029015563, 5.249800272171, 5.250572900878, 5.251269443902, 5.251889532505, 5.253443638590, 5.254378778693, 5.255081457559, 5.255550542553, 5.256098449515, 5.256725476488, 5.258139605935, 5.259006068415, 5.259795264493, 5.260585897301, 5.261377972082, 5.262012673667, 5.262727823464, 5.263523817972, 5.264002115908, 5.264720551940, 5.265680319141, 5.267043630424, 5.268169579712, 5.268411234813, 5.269056306572, 5.269944847625, 5.271240524832, 5.272214825817, 5.273109859258, 5.273843533827, 5.274415027729, 5.275642195774, 5.276462238468, 5.277530614116, 5.278354233710, 5.278931698203, 5.279675282583, 5.280751601552, 5.281249265260, 5.282495925236, 5.283162276700, 5.283829652140, 5.284498054707, 5.285335007137, 5.286089645871, 5.286761538454, 5.287855585785, 5.288867927693, 5.289628735739, 5.290560425868, 5.292089334287, 5.292344676469, 5.293709042741, 5.294820755126, 5.296107046367, 5.297052753818, 5.297914278564, 5.298432014944, 5.299209778626, 5.300248968310, 5.301638433945, 5.302596276800, 5.303905840005, 5.305219363988, 5.305833704067, 5.306712842994, 5.307593765166, 5.308476477832, 5.309006967900, 5.309892560544, 5.310868802765, 5.312114475151, 5.313453100937, 5.314079207805, 5.315154638356, 5.316142794600, 5.317043073699, 5.317584138323, 5.318216233532, 5.319483190619, 5.320208829020, 5.320935681879, 5.322119418488, 5.322575562299, 5.323214969581, 5.324771746407, 5.325874017257, 5.326979092871, 5.327532686932, 5.328734567053, 5.329939782527, 5.330962199115, 5.331893762067, 5.332547047110, 5.333482019445, 5.334137699797, 5.335264031481, 5.335922409815, 5.336581787747, 5.337147766735, 5.337525496250, 5.338850142755, 5.339419087573, 5.340749531227, 5.341320971418, 5.342179543984, 5.342657263185, 5.343422708604, 5.344669441991, 5.345342245350, 5.346112441929, 5.347366931917, 5.348140730753, 5.349012905617, 5.349789645340, 5.350665141288, 5.351542405717, 5.352519226826, 5.353596273777, 5.354283060630, 5.355069292086, 5.355955507185, 5.356646038023, 5.357238796735, 5.358823453387, 5.359717370303, 5.361311113310, 5.361910278015, 5.363111093016, 5.363913484897, 5.364918563989, 5.366329593949, 5.367137959900, 5.368657713516, 5.369572124975, 5.370488465800, 5.371202514433, 5.372020010170, 5.372736583432, 5.373659632625, 5.374481771028, 5.375408540873, 5.376233999866, 5.377371573871, 5.378408324141, 5.379135524735, 5.380697924124, 5.381847266621, 5.383524486111, 5.385102783967, 5.386263858738, 5.387004343968, 5.388170520502, 5.389126999620, 5.389659288548, 5.390405590775, 5.391153177674, 5.391902053675, 5.392866795608, 5.393618634889, 5.394371777992, 5.394587201847, 5.395234115296, 5.396314450385, 5.397180657567, 5.398157210218, 5.398591939465, 5.399027104313, 5.400226060854, 5.400772137226, 5.401647290131, 5.402633949734, 5.403622856002, 5.404063093731, 5.405386490840, 5.406049704736, 5.406603157700, 5.407268233606, 5.408712734941, 5.409715596282, 5.410273743746, 5.411056357260, 5.411952503014, 5.413300198376, 5.414426481377, 5.415895029601, 5.416801226031, 5.418050341627, 5.418619311290, 5.419759491735, 5.421017157297, 5.421590029669, 5.422163658707, 5.423083044035, 5.423658649794, 5.424581208786, 5.425505731715, 5.426084559578, 5.427012291802, 5.427593132442, 5.428407616639, 5.429223631205, 5.430743166671, 5.431563585583, 5.432738307646, 5.434034182553, 5.434388275098, 5.435807539374, 5.436993812938, 5.438183335681, 5.438660058541, 5.439017944414, 5.439854160151, 5.440452444420, 5.441171475183, 5.442011851775, 5.443938840990, 5.444784594874, 5.445026541667, 5.445510839996, 5.445995678988, 5.446723953863, 5.447331783888, 5.448184177649, 5.449527042893, 5.450751443146, 5.451733454829, 5.452471423540, 5.453210648369, 5.453951133598, 5.454692883534, 5.455559862682, 5.456180194857, 5.456925764966, 5.458795309317, 5.459545386329, 5.460547508451, 5.461300620458, 5.462432742847, 5.463062977295, 5.464073258604, 5.464452720823, 5.464959186749, 5.465466243995, 5.466100899203, 5.467372998771, 5.467882883751, 5.469160221383, 5.469672210222, 5.470184803355, 5.470826396738, 5.471597562046, 5.472498989019, 5.473531487531, 5.474307475495, 5.475733731233, 5.477034404508, 5.477686204843, 5.478208350361, 5.479516467259, 5.479909671887, 5.480828536178, 5.481486060122, 5.482408269288, 5.483861423283, 5.485053994692, 5.485850865525, 5.486649201194, 5.487582451399, 5.489053051327, 5.489723155583, 5.490394295388, 5.491201034596, 5.491874463917, 5.492819022740, 5.493224463393, 5.494307492588, 5.495529137506, 5.496481687276, 5.497436330893, 5.498256270372, 5.498803757973, 5.499900808084, 5.500312917382, 5.501000636420, 5.502103257087, 5.503208684300, 5.503623945988, 5.504733255612, 5.505706231335, 5.506681391768, 5.507239610973, 5.507518989871, 5.508498233763, 5.509620079997, 5.510463370518, 5.511308301683, 5.512861624523, 5.513286224018, 5.513853003193, 5.514136670403, 5.514846650096, 5.515700160653, 5.516270101000, 5.516840790283, 5.517698232777, 5.518270803040, 5.519131076313, 5.519849274727, 5.520856752021, 5.521577812260, 5.521722168080, 5.522589312093, 5.524038410808, 5.524619406857, 5.525201181199, 5.525492360883, 5.526075306584, 5.526805090795, 5.527828853308, 5.528414945815, 5.528708288941, 5.529295570277, 5.530325227448, 5.531652669588, 5.532539890493, 5.532984181562, 5.534171184643, 5.534319788402, 5.534617148552, 5.535361440905, 5.536107011014, 5.536853863274, 5.537602002101, 5.538801711377, 5.539703673243, 5.541060138110, 5.541362150974, 5.542118103266, 5.542875373697, 5.544241796896, 5.545765104252, 5.546528766277, 5.548213564476, 5.548674191510, 5.550059011227, 5.550830267835, 5.551448260798, 5.551602896542, 5.552531869050, 5.553929064299, 5.555018887912, 5.555799011136, 5.556737012541, 5.558461961298, 5.559720786764, 5.561457651214, 5.562090964460, 5.563360368307, 5.564474148501, 5.565431095966, 5.565750547604, 5.567030709126, 5.567190994967, 5.568475415813, 5.568958054664, 5.570732333567, 5.571379327328, 5.572351628813, 5.573488738635, 5.573976984310, 5.574302786637, 5.575445023393, 5.576099081472, 5.576590272267, 5.577574323629, 5.578396073130, 5.578560609780, 5.579549140894, 5.581201709410, 5.582528306797, 5.583525920900, 5.584692707774, 5.586030028252, 5.586365002801, 5.586700235919, 5.587203571283, 5.587875593827, 5.588380294037, 5.588885581449, 5.589560213790, 5.590404980603, 5.591251393816, 5.592609095529, 5.594312213327, 5.594995334950, 5.596707854842, 5.597051170656, 5.597738617545, 5.598599459218, 5.598944274228, 5.599462010608, 5.600672467841, 5.601886308269, 5.602233743874, 5.603277721496, 5.603626272463, 5.604324214730, 5.605023280445, 5.605723473232, 5.607478910068, 5.608359296508, 5.609948503541, 5.611011214875, 5.612076533027, 5.613144470815, 5.613857891069, 5.614572485195, 5.615288257062, 5.616364131638, 5.617262734239, 5.618343517414, 5.618884919290, 5.619969752032, 5.620331965966, 5.621057301387, 5.621420423884, 5.622693748932, 5.623240604595, 5.623788149717, 5.624519285381, 5.624885315308, 5.625801742071, 5.626536278368, 5.627087997030, 5.628193541493, 5.629116983222, 5.629671992220, 5.630598586303, 5.631899148291, 5.633017024022, 5.633763876282, 5.634886568372, 5.635824367229, 5.636576067083, 5.637706062036, 5.638649975648, 5.639028116274, 5.639974910811, 5.640543979879, 5.641113795594, 5.642446280257, 5.643209539648, 5.643974142807, 5.645315446045, 5.646468440922, 5.647046088290, 5.647624504999, 5.648010544564, 5.648783654661, 5.649558143465, 5.649751981666, 5.650916831220, 5.652474840001, 5.653451441452, 5.653451441452, 5.654430243944, 5.655018586073, 5.656394491896, 5.657379957447, 5.657774770639, 5.658565475422, 5.660151216962, 5.661145253748, 5.661942124580, 5.662740460250, 5.663740447986, 5.664542099311, 5.665546248849, 5.665948559653, 5.667157733006, 5.668370282370, 5.669586226651, 5.670398751643, 5.671009144551, 5.672028376377, 5.672641065614, 5.672641065614, 5.673869043289, 5.674689628289, 5.675923420261, 5.676954264518, 5.677367288308, 5.677987561418, 5.679438319805, 5.680269505670, 5.680477550935, 5.681519274825, 5.682354456779, 5.683400697906, 5.684029654543, 5.685290307045, 5.686132779631, 5.686765708305, 5.687611050629, 5.688246138944, 5.689944262249, 5.690582774222, 5.691222226335, 5.691862621362, 5.692290076595, 5.693146251307, 5.694004117229, 5.694648630553, 5.695294101787, 5.695509472227, 5.695509472227, 5.696372023616, 5.697452627513, 5.697885623044, 5.698319050706, 5.699187205882, 5.700057099977, 5.700492701300, 5.700492701300, 5.701146923590, 5.701365216876, 5.702239488901, 5.702896349851, 5.703773712739, 5.704872914748, 5.706195640081, 5.706858516549, 5.708409174342, 5.708853238268, 5.709520186669, 5.710634048480, 5.711750774428, 5.712870379281, 5.714442690992, 5.714667772356, 5.716020715762, 5.717377887122, 5.718058066559, 5.719421629632, 5.720105019988, 5.721246399047, 5.722620025333, 5.722849386036, 5.724458311599, 5.725149679983, 5.725842150736, 5.726535727379, 5.727230413448, 5.727694155598, 5.727926212500, 5.729787145104, 5.730487055782, 5.731422028117, 5.732828271597, 5.733298033116, 5.734474664781, 5.735418270762, 5.735890843694, 5.736600668666, 5.737311655698, 5.738499226802, 5.739213331345, 5.739690054205, 5.741123370628, 5.741602195904, 5.741841806659, 5.742801573861, 5.742801573861, 5.743763466794, 5.744968836654, 5.745451922891, 5.746419710438, 5.747875447494, 5.748118545447, 5.749092300299, 5.749579997691, 5.749824051916, 5.750801642609, 5.751781438810, 5.752763450493, 5.753747687701, 5.754734160543, 5.755475488430, 5.755970410970, 5.756713853917, 5.757707095017, 5.758702612890, 5.759700417997, 5.760450279160, 5.760950906860, 5.762707662433, 5.763462738511, 5.764471553092, 5.765735875621, 5.765989182413, 5.766242637034, 5.767003889608, 5.767257937279, 5.768785352037, 5.769551078622, 5.770062314092, 5.771343041891, 5.772627557710, 5.773657912836, 5.774690718274, 5.774949303862, 5.775985188627, 5.777023550107, 5.778585762158, 5.779107750780, 5.779891911960, 5.781202001888, 5.782252926737, 5.783042792639, 5.783570169124, 5.784362436565, 5.784891418947, 5.784891418947, 5.786482243004, 5.788078915691, 5.788345599447, 5.789413975095, 5.791021482724, 5.791827473333, 5.791827473333, 5.792634962531, 5.793713955588, 5.794254459057, 5.795608668081, 5.795880017344, 5.796423225022, 5.797239312607, 5.797784224199, 5.798876102793, 5.799696817018, 5.800793520838, 5.801617869992, 5.803823814960, 5.804377056413, 5.805485658118, 5.806318970459, 5.807153884811, 5.807711387432, 5.808269606637, 5.808548985535, 5.809668301830, 5.810509686301, 5.810790510418, 5.811352704000, 5.811634073937, 5.811634073937, 5.812479279164, 5.814174640387, 5.815592514588, 5.816730156317, 5.817585347565, 5.818442226137, 5.820448208835, 5.821598658466, 5.822752163744, 5.823330067332, 5.824778199657, 5.825649402521, 5.825940192275, 5.826813731588, 5.827689031478, 5.828858848972, 5.828858848972, 5.829151796357, 5.830619504688, 5.830913642513, 5.831502516477, 5.831797253157, 5.832387327272, 5.833273944420, 5.833866029695, 5.834458923278, 5.835647144216, 5.836540448230, 5.838631997765, 5.839831707041, 5.840432806766, 5.842241113953, 5.843148098930, 5.844967771209, 5.846185135655, 5.847100403606, 5.847711655617, 5.848936746646, 5.849550590539, 5.850472986246, 5.851089006891, 5.852942328972, 5.854182285508, 5.854492828590, 5.855114581713, 5.856360764725, 5.857297754262, 5.857923538927, 5.858550226600, 5.859177819891, 5.860435733824, 5.862962545210, 5.863596551866, 5.865185629680, 5.865822892423, 5.866780543268, 5.868381335651, 5.869666231505, 5.869988050328, 5.870632404277, 5.870632404277, 5.870954940112, 5.871923987331, 5.871923987331, 5.873219422988, 5.873868592738, 5.874844170419, 5.877129077136, 5.877456475931, 5.877456475931, 5.877784121727, 5.878440155813, 5.879097182385, 5.879426068794, 5.880084589742, 5.881404634776, 5.882066164960, 5.883724412419, 5.884056823061, 5.885389015768, 5.886056647693, 5.887060023916, 5.887730231583, 5.889410289701, 5.890084136976, 5.891096872333, 5.891773343625, 5.892790030352, 5.894149325615, 5.895854449446, 5.895854449446, 5.897223385117, 5.898940645092, 5.900319358891, 5.900319358891, 5.901702463505, 5.902048929006, 5.902048929006, 5.902395671126, 5.903089986992, 5.903785414654, 5.905528871358, 5.906228218501, 5.906928693624, 5.907279355316, 5.909389292172, 5.909741947069, 5.911863911299, 5.912928794093, 5.913996294382, 5.915066425063, 5.915781313261, 5.916139199133, 5.917573699139, 5.917933065715, 5.919373513078, 5.919734372660, 5.920456992597, 5.921543181947, 5.921905849594, 5.921905849594, 5.921905849594, 5.922995673207, 5.924818145381, 5.926648297613, 5.927015255372, 5.927382523455, 5.928117992694, 5.930331903088, 5.931072388318, 5.932185488838, 5.932929143955, 5.933674074638, 5.934047019686, 5.935167780261, 5.936666641048, 5.936666641048, 5.937418015772, 5.938547520913, 5.940436582099, 5.941573975543, 5.941953769605, 5.942714355582, 5.944239535312, 5.944621668625, 5.945386945443, 5.946153573148, 5.947306058075, 5.948461609485, 5.949233688767, 5.950007143080, 5.950781977330, 5.950781977330, 5.951558196450, 5.952335805398, 5.953895212754, 5.956637721979, 5.957818405484, 5.958607314842, 5.960189445852, 5.960982678003, 5.961777361631, 5.963770455914, 5.964170174747, 5.964170174747, 5.964170174747, 5.964170174747, 5.964570261815, 5.965772739229, 5.966576244513, 5.967381239149, 5.968187728670, 5.968591535748, 5.969400278034, 5.970210529168, 5.971022294791, 5.972650392225, 5.974284616099, 5.975104039893, 5.975925012693, 5.976336081802, 5.976747540366, 5.977984260182, 5.977984260182, 5.979224511806, 5.979224511806, 5.980053318321, 5.981299501334, 5.982549270489, 5.982549270489, 5.982966660701, 5.983384452443, 5.983384452443, 5.984221243611, 5.984640244591, 5.986741334716, 5.987584625238, 5.989276134608, 5.990124366288, 5.990549104201, 5.991399828238, 5.992252221999, 5.993106292052, 5.993106292052, 5.993962045003, 5.994819487496, 5.996108833763, 5.996108833763, 5.996539467890, 5.997402019280, 5.998699066980, 5.999565922521, 6.001304841688, 6.001304841688, 6.001740661576, 6.001740661576, 6.002176919254, 6.002613615603, 6.004364805402, 6.004803708403, 6.006563769502, 6.007446482168, 6.007446482168, 6.007888512213, 6.008773924308, 6.010105436281, 6.010105436281, 6.010550182333, 6.010550182333, 6.010995384301, 6.010995384301, 6.010995384301, 6.012333735074, 6.013676222949, 6.014124642692, 6.015022873585, 6.015922966097, 6.018181392829, 6.018634490921, 6.019542107724, 6.020451625296, 6.021363051616, 6.022276394711, 6.022276394711, 6.023191662662, 6.024108863598, 6.025028005702, 6.026410376573, 6.026872146400, 6.026872146400, 6.027797161621, 6.029188389127, 6.031050319019, 6.031517051446, 6.031517051446, 6.032452023781, 6.033389013318, 6.034798298974, 6.034798298974, 6.036684488614, 6.037157318799, 6.038578905934, 6.039053804266, 6.040005161672, 6.040481623027, 6.040481623027, 6.040958607679, 6.042392712940, 6.043831569525, 6.044312249686, 6.044312249686, 6.045275209021, 6.045275209021, 6.045275209021, 6.046240308267, 6.046240308267, 6.046723663333, 6.046723663333, 6.048176964684, 6.049635145624, 6.050122295963, 6.051098239030, 6.051098239030, 6.052076380168, 6.052566278113, 6.053547734987, 6.054531414868, 6.056011124926, 6.056505484094, 6.057991946978, 6.058488567366, 6.058985756294, 6.060980223551, 6.062482107983, 6.063486257521, 6.063486257521, 6.063486257521, 6.063989204285, 6.064492734175, 6.064492734175, 6.066512712151, 6.068542129311, 6.068542129311, 6.069050968832, 6.069050968832, 6.071092309756, 6.071604147743, 6.072116589669, 6.074687908500, 6.074687908500, 6.074687908500, 6.076238039171, 6.076755981370, 6.077274542007, 6.077274542007, 6.077274542007, 6.078313524516, 6.078313524516, 6.078313524516, 6.079876673709, 6.079876673709, 6.080398976216, 6.080921907624, 6.080921907624, 6.081445469450, 6.083546051450, 6.084600164788, 6.086186147616, 6.086186147616, 6.086186147616, 6.087246696329, 6.087246696329, 6.088309841246, 6.088842391260, 6.089909454406, 6.091514981121, 6.091514981121, 6.093126465278, 6.094743951252, 6.095284454721, 6.096367483916, 6.096367483916, 6.096367483916, 6.097997108649, 6.098541678604, 6.099086932262, 6.100179497573, 6.100179497573, 6.101274818411, 6.101274818411, 6.102922996791, 6.104025267641, 6.105130343255, 6.105683937316, 6.106238237942, 6.106238237942, 6.107348966123, 6.108462542327, 6.109020403010, 6.111259039317, 6.113509274828, 6.113509274828, 6.114073660199, 6.114638779968, 6.116338564846, 6.118045028660, 6.118615343229, 6.119758224105, 6.120330794368, 6.120904120500, 6.121478204499, 6.122628654130, 6.123782159408, 6.124938736608, 6.125518182301, 6.125518182301, 6.126098402136, 6.126679398185, 6.127261172527, 6.127261172527, 6.128427064454, 6.130181792021, 6.130768280269, 6.130768280269, 6.131355561605, 6.131355561605, 6.133712660916, 6.135488918942, 6.135488918942, 6.136082623042, 6.136677139880, 6.136677139880, 6.137272471682, 6.137272471682, 6.137868620687, 6.139063379300, 6.139063379300, 6.140861702705, 6.142064735281, 6.142064735281, 6.142667503569, 6.143875555758, 6.143875555758, 6.145693958199, 6.148130399270, 6.148130399270, 6.148741651281, 6.149353764817, 6.149966742310, 6.151195298948, 6.151810883009, 6.153662887870, 6.154901959986, 6.155522824254, 6.155522824254, 6.156767221902, 6.158015195410, 6.159893905543, 6.162411561764, 6.163675884293, 6.164943898280, 6.165579296318, 6.167491087294, 6.167491087294, 6.169411331315, 6.170053304058, 6.170053304058, 6.170053304058, 6.170696227169, 6.171984935776, 6.173277479831, 6.173925197299, 6.173925197299, 6.174573882232, 6.175874166083, 6.175874166083, 6.176525770830, 6.177178354697, 6.177831920632, 6.177831920632, 6.179142010560, 6.181114585406, 6.183096160624, 6.184422251676, 6.185752404268, 6.186419011432, 6.187086643357, 6.189767482005, 6.191114132640, 6.191114132640, 6.191114132640, 6.191114132640, 6.191789027076, 6.193141970481, 6.194499141842, 6.194499141842, 6.194499141842, 6.195179321279, 6.196542884352, 6.197910742118, 6.198596289983, 6.201349354555, 6.201349354555, 6.202732459169, 6.202732459169, 6.204119982656, 6.206209615309, 6.206908399823, 6.207608310502, 6.210419287836, 6.210419287836, 6.211831628859, 6.211831628859, 6.213248577854, 6.214670164989, 6.216811308925, 6.218244625348, 6.218244625348, 6.218244625348, 6.218963061379, 6.220403508742, 6.221848749616, 6.221848749616, 6.222573177611, 6.223298816012, 6.223298816012, 6.225483034271, 6.225483034271, 6.225483034271, 6.225483034271, 6.225483034271, 6.226213555019, 6.226213555019, 6.228412519119, 6.229884705213, 6.229884705213, 6.230622673924, 6.230622673924, 6.232844133918, 6.232844133918, 6.234331445241, 6.235823867610, 6.237321436273, 6.238072161579, 6.239577516577, 6.241088107602, 6.242603971207, 6.243363891754, 6.244125144328, 6.244125144328, 6.244125144328, 6.244887733605, 6.244887733605, 6.246416941107, 6.247183568812, 6.247951552181, 6.248720896017, 6.251037138744, 6.251811972994, 6.253365801062, 6.254925208418, 6.255707016877, 6.255707016877, 6.258848401148, 6.258848401148, 6.259637310506, 6.260427655550, 6.260427655550, 6.262012673667, 6.263603497723, 6.264401100302, 6.265200170411, 6.266000713462, 6.266000713462, 6.266802734893, 6.268411234813, 6.269217724334, 6.270835210307, 6.270835210307, 6.272458742971, 6.273272790973, 6.274088367705, 6.274905478919, 6.276544327965, 6.277366077466, 6.279014255846, 6.279840696594, 6.280668713016, 6.280668713016, 6.281498311133, 6.283162276700, 6.284832642152, 6.285670240255, 6.286509456906, 6.287350298373, 6.288192770959, 6.288192770959, 6.288192770959, 6.289036881005, 6.289036881005, 6.289882634888, 6.289882634888, 6.290730039024, 6.294136287716, 6.295849483160, 6.295849483160, 6.297569463554, 6.298432014944, 6.298432014944, 6.299296282855, 6.299296282855, 6.301899454377, 6.301899454377, 6.302770657240, 6.304518323510, 6.306273051076, 6.306273051076, 6.308034897233, 6.308918507877, 6.309803919971, 6.309803919971, 6.310691140876, 6.311580177997, 6.312471038785, 6.312471038785, 6.312471038785, 6.313363730738, 6.313363730738, 6.313363730738, 6.314258261398, 6.314258261398, 6.316052869248, 6.316952961761, 6.318758762624, 6.318758762624, 6.319664486585, 6.319664486585, 6.319664486585, 6.319664486585, 6.319664486585, 6.319664486585, 6.320572103388, 6.321481620960, 6.321481620960, 6.321481620960, 6.323306390375, 6.324221658326, 6.325138859262, 6.325138859262, 6.325138859262, 6.325138859262, 6.325138859262, 6.326058001366, 6.326058001366, 6.327902142064, 6.328827157285, 6.328827157285, 6.328827157285, 6.328827157285, 6.329754146926, 6.330683119434, 6.330683119434, 6.331614083310, 6.331614083310, 6.331614083310, 6.331614083310, 6.334419008982, 6.335358024444, 6.335358024444, 6.337242168318, 6.339134521996, 6.339134521996, 6.341035157336, 6.341035157336, 6.341988603343, 6.343901797987, 6.343901797987, 6.343901797987, 6.344861565189, 6.347753658997, 6.348721986002, 6.348721986002, 6.349692476868, 6.350665141288, 6.350665141288, 6.350665141288, 6.353596273777, 6.355561410532, 6.355561410532, 6.356547323514, 6.356547323514, 6.357535479758, 6.358525889496, 6.359518563030, 6.362510270487, 6.363512103647, 6.364516253185, 6.365522729839, 6.365522729839, 6.365522729839, 6.366531544420, 6.366531544420, 6.368556230987, 6.370590400897, 6.372634143407, 6.372634143407, 6.372634143407, 6.373659632625, 6.374687549038, 6.374687549038, 6.374687549038, 6.374687549038, 6.374687549038, 6.374687549038, 6.374687549038, 6.374687549038, 6.376750709602, 6.376750709602, 6.376750709602, 6.377785977034, 6.377785977034, 6.377785977034, 6.377785977034, 6.377785977034, 6.378823718225, 6.378823718225, 6.378823718225, 6.379863945026, 6.381951903288, 6.382999658879, 6.385102783967, 6.387216143280, 6.387216143280, 6.387216143280, 6.387216143280, 6.388276691993, 6.390405590775, 6.391473966423, 6.392544976785, 6.393618634889, 6.394694953859, 6.396855627380, 6.401209493237, 6.402304814074, 6.402304814074, 6.402304814074, 6.403402904374, 6.403402904374, 6.403402904374, 6.404503778174, 6.404503778174, 6.404503778174, 6.404503778174, 6.404503778174, 6.406713932980, 6.408935392973, 6.408935392973, 6.410050398674, 6.410050398674, 6.411168274406, 6.412289034981, 6.412289034981, 6.412289034981, 6.412289034981, 6.412289034981, 6.412289034981, 6.414539270491, 6.416801226031, 6.417936637088, 6.419075024324, 6.419075024324, 6.419075024324, 6.419075024324, 6.419075024324, 6.419075024324, 6.419075024324, 6.419075024324, 6.420216403383, 6.420216403383, 6.420216403383, 6.420216403383, 6.420216403383, 6.421360790032, 6.421360790032, 6.421360790032, 6.423658649794, 6.423658649794, 6.425968732272, 6.425968732272, 6.428291168191, 6.429457060118, 6.429457060118, 6.429457060118, 6.429457060118, 6.429457060118, 6.429457060118, 6.429457060118, 6.429457060118, 6.430626090385, 6.430626090385, 6.431798275933, 6.431798275933, 6.432973633841, 6.432973633841, 6.432973633841, 6.435333935748, 6.435333935748, 6.435333935748, 6.435333935748, 6.435333935748, 6.435333935748, 6.436518914606, 6.436518914606, 6.436518914606, 6.436518914606, 6.436518914606, 6.436518914606, 6.437707135544, 6.437707135544, 6.437707135544, 6.437707135544, 6.437707135544, 6.438898616351, 6.438898616351, 6.441291429467, 6.442492798094, 6.443697499233, 6.444905551422, 6.447331783888, 6.447331783888, 6.448550002027, 6.450996737974, 6.450996737974, 6.452225294612, 6.452225294612, 6.452225294612, 6.452225294612, 6.452225294612, 6.453457336522, 6.453457336522, 6.453457336522, 6.454692883534, 6.454692883534, 6.454692883534, 6.454692883534, 6.454692883534, 6.454692883534, 6.457174573041, 6.458420756053, 6.458420756053, 6.458420756053, 6.458420756053, 6.459670525209, 6.462180904927, 6.462180904927, 6.462180904927, 6.462180904927, 6.463441557428, 6.465973893944, 6.467245621007, 6.468521082958, 6.471083299722, 6.471083299722, 6.471083299722, 6.471083299722, 6.471083299722, 6.471083299722, 6.471083299722, 6.471083299722, 6.471083299722, 6.471083299722, 6.472370099129, 6.472370099129, 6.472370099129, 6.472370099129, 6.473660722610, 6.474955192963, 6.477555766494, 6.477555766494, 6.477555766494, 6.477555766494, 6.478861916296, 6.478861916296, 6.480172006224, 6.480172006224, 6.484126156288, 6.485452247340, 6.486782399932, 6.486782399932, 6.488116639021, 6.489454989793, 6.489454989793, 6.490797477669, 6.492144128304, 6.492144128304, 6.492144128304, 6.494850021680, 6.496209316943, 6.497572880016, 6.498940737782, 6.501689446210, 6.501689446210, 6.501689446210, 6.503070351927, 6.505845405982, 6.507239610973, 6.508638306166, 6.510041520575, 6.511449283500, 6.512861624523, 6.514278573518, 6.514278573518, 6.514278573518, 6.515700160653, 6.515700160653, 6.515700160653, 6.515700160653, 6.515700160653, 6.515700160653, 6.517126416391, 6.517126416391, 6.517126416391, 6.517126416391, 6.519993057043, 6.519993057043, 6.519993057043, 6.521433504406, 6.521433504406, 6.521433504406, 6.522878745280, 6.522878745280, 6.522878745280, 6.522878745280, 6.524328811676, 6.524328811676, 6.524328811676, 6.524328811676, 6.524328811676, 6.524328811676, 6.525783735924, 6.525783735924, 6.525783735924, 6.525783735924, 6.525783735924, 6.525783735924, 6.525783735924, 6.525783735924, 6.525783735924, 6.527243550683, 6.527243550683, 6.528708288941, 6.530177984022, 6.533132379646, 6.533132379646, 6.533132379646, 6.533132379646, 6.534617148552, 6.534617148552, 6.534617148552, 6.534617148552, 6.536107011014, 6.537602002101, 6.537602002101, 6.537602002101, 6.537602002101, 6.537602002101, 6.540607512241, 6.540607512241, 6.543633966871, 6.545155139991, 6.548213564476, 6.548213564476, 6.548213564476, 6.549750891681, 6.549750891681, 6.551293680095, 6.554395796726, 6.554395796726, 6.555955204082, 6.555955204082, 6.555955204082, 6.555955204082, 6.555955204082, 6.555955204082, 6.557520230936, 6.559090917935, 6.562249437180, 6.562249437180, 6.562249437180, 6.562249437180, 6.562249437180, 6.562249437180, 6.563837352959, 6.565431095966, 6.565431095966, 6.565431095966, 6.567030709126, 6.567030709126, 6.567030709126, 6.567030709126, 6.567030709126, 6.568636235841, 6.570247719998, 6.570247719998, 6.570247719998, 6.570247719998, 6.570247719998, 6.570247719998, 6.570247719998, 6.573488738635, 6.576754126063, 6.576754126063, 6.578396073130, 6.578396073130, 6.578396073130, 6.580044251510, 6.581698708680, 6.583359492662, 6.583359492662, 6.583359492662, 6.583359492662, 6.583359492662, 6.585026652029, 6.586700235919, 6.586700235919, 6.586700235919, 6.588380294037, 6.590066876669, 6.590066876669, 6.590066876669, 6.590066876669, 6.590066876669, 6.591760034688, 6.591760034688, 6.593459819566, 6.593459819566, 6.596879478824, 6.596879478824, 6.596879478824, 6.598599459218, 6.600326278519, 6.600326278519, 6.605548319174, 6.605548319174, 6.607303046740, 6.607303046740, 6.607303046740, 6.607303046740, 6.610833915635, 6.610833915635, 6.610833915635, 6.610833915635, 6.610833915635, 6.612610173661, 6.612610173661, 6.614393726402, 6.614393726402, 6.614393726402, 6.614393726402, 6.616184634020, 6.616184634020, 6.616184634020, 6.616184634020, 6.616184634020, 6.616184634020, 6.616184634020, 6.617982957425, 6.619788758288, 6.621602099052, 6.621602099052, 6.621602099052, 6.621602099052, 6.621602099052, 6.621602099052, 6.621602099052, 6.623423042943, 6.623423042943, 6.625251653990, 6.625251653990, 6.625251653990, 6.625251653990, 6.627087997030, 6.627087997030, 6.627087997030, 6.627087997030, 6.628932137728, 6.628932137728, 6.630784142590, 6.632644078974, 6.634512015109, 6.634512015109, 6.634512015109, 6.634512015109, 6.636388020108, 6.638272163982, 6.640164517660, 6.640164517660, 6.642065153000, 6.642065153000, 6.642065153000, 6.642065153000, 6.643974142807, 6.645891560853, 6.647817481889, 6.647817481889, 6.647817481889, 6.649751981666, 6.649751981666, 6.649751981666, 6.649751981666, 6.651695136952, 6.653647025549, 6.655607726315, 6.657577319178, 6.657577319178, 6.659555885160, 6.663540266151, 6.663540266151, 6.663540266151, 6.663540266151, 6.663540266151, 6.665546248849, 6.665546248849, 6.665546248849, 6.665546248849, 6.665546248849, 6.665546248849, 6.665546248849, 6.665546248849, 6.667561540084, 6.667561540084, 6.669586226651, 6.669586226651, 6.669586226651, 6.669586226651, 6.669586226651, 6.669586226651, 6.669586226651, 6.669586226651, 6.673664139071, 6.675717544702, 6.675717544702, 6.675717544702, 6.675717544702, 6.675717544702, 6.675717544702, 6.677780705266, 6.677780705266, 6.677780705266, 6.677780705266, 6.677780705266, 6.677780705266, 6.677780705266, 6.677780705266, 6.679853713889, 6.679853713889, 6.679853713889, 6.679853713889, 6.679853713889, 6.679853713889, 6.679853713889, 6.679853713889, 6.679853713889, 6.681936665037, 6.681936665037, 6.681936665037, 6.681936665037, 6.681936665037, 6.684029654543, 6.684029654543, 6.684029654543, 6.684029654543, 6.686132779631, 6.686132779631, 6.686132779631, 6.688246138944, 6.688246138944, 6.688246138944, 6.692503962087, 6.692503962087, 6.692503962087, 6.692503962087, 6.692503962087, 6.692503962087, 6.692503962087, 6.692503962087, 6.692503962087, 6.692503962087, 6.692503962087, 6.692503962087, 6.692503962087, 6.692503962087, 6.692503962087, 6.692503962087, 6.692503962087, 6.692503962087, 6.692503962087, 6.694648630553, 6.694648630553, 6.694648630553, 6.694648630553, 6.694648630553, 6.694648630553, 6.694648630553, 6.694648630553, 6.694648630553, 6.694648630553, 6.694648630553, 6.694648630553, 6.696803942580, 6.698970004336, 6.701146923590, 6.703334809738, 6.703334809738, 6.703334809738, 6.705533773838, 6.705533773838, 6.705533773838, 6.705533773838, 6.705533773838, 6.707743928644, 6.707743928644, 6.709965388637, 6.709965388637, 6.709965388637, 6.709965388637, 6.709965388637, 6.709965388637, 6.709965388637, 6.709965388637, 6.709965388637, 6.709965388637, 6.709965388637, 6.709965388637, 6.709965388637, 6.709965388637, 6.709965388637, 6.709965388637, 6.709965388637, 6.709965388637, 6.712198270070, 6.714442690992, 6.714442690992, 6.714442690992, 6.714442690992, 6.714442690992, 6.714442690992, 6.714442690992, 6.714442690992, 6.714442690992, 6.714442690992, 6.714442690992, 6.714442690992, 6.714442690992, 6.714442690992, 6.714442690992, 6.714442690992, 6.714442690992, 6.716698771296, 6.716698771296, 6.716698771296, 6.718966632752, 6.723538195827, 6.725842150736, 6.730487055782, 6.730487055782, 6.730487055782, 6.730487055782, 6.730487055782, 6.730487055782, 6.732828271597, 6.732828271597, 6.735182176990, 6.735182176990, 6.735182176990, 6.737548910270, 6.742321425131, 6.744727494897, 6.744727494897, 6.744727494897, 6.747146969020, 6.747146969020, 6.749579997691, 6.754487332186, 6.754487332186, 6.754487332186, 6.756961951314, 6.756961951314, 6.756961951314, 6.756961951314, 6.759450751717, 6.759450751717, 6.759450751717, 6.759450751717, 6.761953896871, 6.761953896871, 6.761953896871, 6.767003889608, 6.767003889608, 6.769551078622, 6.769551078622, 6.769551078622, 6.769551078622, 6.769551078622, 6.774690718274, 6.774690718274, 6.774690718274, 6.774690718274, 6.774690718274, 6.774690718274, 6.777283528852, 6.777283528852, 6.777283528852, 6.777283528852, 6.777283528852, 6.777283528852, 6.777283528852, 6.777283528852, 6.777283528852, 6.779891911960, 6.779891911960, 6.779891911960, 6.782516055786, 6.785156151952, 6.785156151952, 6.785156151952, 6.787812395596, 6.787812395596, 6.787812395596, 6.787812395596, 6.787812395596, 6.787812395596, 6.787812395596, 6.787812395596, 6.787812395596, 6.787812395596, 6.787812395596, 6.787812395596, 6.787812395596, 6.787812395596, 6.787812395596, 6.787812395596, 6.787812395596, 6.790484985457, 6.793174123968, 6.793174123968, 6.793174123968, 6.793174123968, 6.793174123968, 6.793174123968, 6.798602875680, 6.801342913046, 6.804100347591, 6.804100347591, 6.804100347591, 6.804100347591, 6.804100347591, 6.806875401646, 6.806875401646, 6.806875401646, 6.806875401646, 6.806875401646, 6.809668301830, 6.809668301830, 6.812479279164, 6.815308569182, 6.815308569182, 6.815308569182, 6.815308569182, 6.815308569182, 6.815308569182, 6.815308569182, 6.815308569182, 6.815308569182, 6.815308569182, 6.815308569182, 6.818156412055, 6.818156412055, 6.818156412055, 6.818156412055, 6.818156412055, 6.821023052707, 6.821023052707, 6.821023052707, 6.821023052707, 6.823908740944, 6.823908740944, 6.823908740944, 6.823908740944, 6.823908740944, 6.826813731588, 6.826813731588, 6.826813731588, 6.826813731588, 6.826813731588, 6.832682665252, 6.832682665252, 6.832682665252, 6.835647144216, 6.835647144216, 6.838631997765, 6.838631997765, 6.838631997765, 6.838631997765, 6.838631997765, 6.838631997765, 6.838631997765, 6.838631997765, 6.838631997765, 6.838631997765, 6.838631997765, 6.844663962535, 6.844663962535, 6.844663962535, 6.844663962535, 6.844663962535, 6.844663962535, 6.844663962535, 6.844663962535, 6.844663962535, 6.844663962535, 6.844663962535, 6.847711655617, 6.847711655617, 6.847711655617, 6.847711655617, 6.847711655617, 6.847711655617, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.853871964322, 6.853871964322, 6.853871964322, 6.856985199746, 6.860120913599, 6.860120913599, 6.860120913599, 6.860120913599, 6.860120913599, 6.860120913599, 6.860120913599, 6.860120913599, 6.863279432844, 6.863279432844, 6.863279432844, 6.863279432844, 6.872895201635, 6.872895201635, 6.872895201635, 6.872895201635, 6.872895201635, 6.872895201635, 6.872895201635, 6.872895201635, 6.872895201635, 6.876148359033, 6.876148359033, 6.879426068794, 6.879426068794, 6.879426068794, 6.879426068794, 6.879426068794, 6.886056647693, 6.889410289701, 6.889410289701, 6.889410289701, 6.889410289701, 6.889410289701, 6.889410289701, 6.889410289701, 6.892790030352, 6.899629454882, 6.899629454882, 6.899629454882, 6.899629454882, 6.899629454882, 6.899629454882, 6.899629454882, 6.899629454882, 6.899629454882, 6.899629454882, 6.899629454882, 6.899629454882, 6.899629454882, 6.903089986992, 6.903089986992, 6.906578314838, 6.906578314838, 6.906578314838, 6.906578314838, 6.906578314838, 6.906578314838, 6.906578314838, 6.906578314838, 6.906578314838, 6.906578314838, 6.906578314838, 6.906578314838, 6.910094888561, 6.913640169325, 6.913640169325, 6.913640169325, 6.917214629684, 6.917214629684, 6.917214629684, 6.920818753952, 6.920818753952, 6.920818753952, 6.924453038607, 6.924453038607, 6.924453038607, 6.928117992694, 6.931814138254, 6.931814138254, 6.935542010773, 6.939302159646, 6.939302159646, 6.939302159646, 6.943095148664, 6.943095148664, 6.943095148664, 6.943095148664, 6.943095148664, 6.943095148664, 6.946921556517, 6.946921556517, 6.946921556517, 6.946921556517, 6.946921556517, 6.946921556517, 6.946921556517, 6.946921556517, 6.946921556517, 6.950781977330, 6.954677021213, 6.954677021213, 6.958607314842, 6.958607314842, 6.958607314842, 6.962573502059, 6.962573502059, 6.962573502059, 6.962573502059, 6.962573502059, 6.966576244513, 6.966576244513, 6.966576244513, 6.966576244513, 6.970616222315, 6.970616222315, 6.970616222315, 6.970616222315, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.978810700930, 6.978810700930, 6.978810700930, 6.978810700930, 6.978810700930, 6.978810700930, 6.978810700930, 6.978810700930, 6.978810700930, 6.978810700930, 6.982966660701, 6.982966660701, 6.982966660701, 6.982966660701, 6.982966660701, 6.982966660701, 6.982966660701, 6.982966660701, 6.982966660701, 6.982966660701, 6.982966660701, 6.982966660701, 6.982966660701, 6.982966660701, 6.982966660701, 6.987162775295, 6.987162775295, 6.991399828238, 6.991399828238, 6.995678626217, 6.995678626217, 6.995678626217, 6.995678626217, 6.995678626217, 6.995678626217, 6.995678626217, 6.995678626217, 7.000000000000, 7.004364805402, 7.004364805402, 7.004364805402, 7.008773924308, 7.008773924308, 7.008773924308, 7.008773924308, 7.008773924308, 7.008773924308, 7.008773924308, 7.008773924308, 7.008773924308, 7.008773924308, 7.013228265734, 7.017728766960, 7.017728766960, 7.017728766960, 7.017728766960, 7.017728766960, 7.017728766960, 7.017728766960, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.031517051446, 7.031517051446, 7.031517051446, 7.036212172654, 7.036212172654, 7.036212172654, 7.036212172654, 7.036212172654, 7.036212172654, 7.036212172654, 7.036212172654, 7.036212172654, 7.036212172654, 7.036212172654, 7.036212172654, 7.036212172654, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.045757490561, 7.050609993355, 7.050609993355, 7.050609993355, 7.050609993355, 7.050609993355, 7.050609993355, 7.055517327850, 7.055517327850, 7.055517327850, 7.055517327850, 7.055517327850, 7.055517327850, 7.055517327850, 7.055517327850, 7.065501548756, 7.065501548756, 7.065501548756, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.075720713938, 7.075720713938, 7.075720713938, 7.075720713938, 7.075720713938, 7.075720713938, 7.075720713938, 7.080921907624, 7.080921907624, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.096910013008, 7.096910013008, 7.096910013008, 7.096910013008, 7.096910013008, 7.096910013008, 7.096910013008, 7.096910013008, 7.096910013008, 7.096910013008, 7.096910013008, 7.096910013008, 7.096910013008, 7.102372908710, 7.102372908710, 7.102372908710, 7.107905397310, 7.107905397310, 7.107905397310, 7.107905397310, 7.107905397310, 7.107905397310, 7.113509274828, 7.113509274828, 7.119186407719, 7.119186407719, 7.119186407719, 7.119186407719, 7.119186407719, 7.119186407719, 7.119186407719, 7.119186407719, 7.119186407719, 7.119186407719, 7.119186407719, 7.119186407719, 7.124938736608, 7.124938736608, 7.124938736608, 7.124938736608, 7.124938736608, 7.124938736608, 7.124938736608, 7.124938736608, 7.124938736608, 7.124938736608, 7.124938736608, 7.124938736608, 7.124938736608, 7.124938736608, 7.130768280269, 7.130768280269, 7.136677139880, 7.136677139880, 7.136677139880, 7.136677139880, 7.136677139880, 7.136677139880, 7.136677139880, 7.136677139880, 7.136677139880, 7.136677139880, 7.142667503569, 7.142667503569, 7.142667503569, 7.142667503569, 7.142667503569, 7.154901959986, 7.154901959986, 7.154901959986, 7.154901959986, 7.154901959986, 7.154901959986, 7.167491087294, 7.167491087294, 7.167491087294, 7.173925197299, 7.173925197299, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.187086643357, 7.187086643357, 7.187086643357, 7.187086643357, 7.187086643357, 7.193820026016, 7.193820026016, 7.207608310502, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.229147988358, 7.229147988358, 7.229147988358, 7.229147988358, 7.229147988358, 7.229147988358, 7.229147988358, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.244125144328, 7.244125144328, 7.244125144328, 7.244125144328, 7.244125144328, 7.244125144328, 7.244125144328, 7.244125144328, 7.244125144328, 7.244125144328, 7.244125144328, 7.244125144328, 7.244125144328, 7.244125144328, 7.244125144328, 7.244125144328, 7.244125144328, 7.244125144328, 7.251811972994, 7.251811972994, 7.251811972994, 7.251811972994, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.267606240177, 7.267606240177, 7.267606240177, 7.275724130399, 7.275724130399, 7.275724130399, 7.275724130399, 7.275724130399, 7.275724130399, 7.275724130399, 7.275724130399, 7.275724130399, 7.275724130399, 7.275724130399, 7.275724130399, 7.275724130399, 7.275724130399, 7.283996656365, 7.283996656365, 7.283996656365, 7.283996656365, 7.283996656365, 7.283996656365, 7.283996656365, 7.283996656365, 7.283996656365, 7.283996656365, 7.283996656365, 7.283996656365, 7.283996656365, 7.283996656365, 7.283996656365, 7.292429823902, 7.292429823902, 7.292429823902, 7.292429823902, 7.292429823902, 7.292429823902, 7.292429823902, 7.301029995664, 7.301029995664, 7.301029995664, 7.301029995664, 7.301029995664, 7.301029995664, 7.301029995664, 7.301029995664, 7.309803919971, 7.309803919971, 7.309803919971, 7.309803919971, 7.309803919971, 7.309803919971, 7.309803919971, 7.309803919971, 7.309803919971, 7.309803919971, 7.309803919971, 7.309803919971, 7.309803919971, 7.309803919971, 7.309803919971, 7.309803919971, 7.309803919971, 7.309803919971, 7.309803919971, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.327902142064, 7.327902142064, 7.327902142064, 7.327902142064, 7.327902142064, 7.327902142064, 7.327902142064, 7.327902142064, 7.327902142064, 7.327902142064, 7.327902142064, 7.327902142064, 7.327902142064, 7.327902142064, 7.327902142064, 7.327902142064, 7.327902142064, 7.327902142064, 7.327902142064, 7.327902142064, 7.327902142064, 7.327902142064, 7.327902142064, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.356547323514, 7.356547323514, 7.366531544420, 7.366531544420, 7.366531544420, 7.366531544420, 7.376750709602, 7.376750709602, 7.376750709602, 7.376750709602, 7.376750709602, 7.376750709602, 7.376750709602, 7.376750709602, 7.387216143280, 7.387216143280, 7.387216143280, 7.387216143280, 7.387216143280, 7.387216143280, 7.387216143280, 7.387216143280, 7.387216143280, 7.387216143280, 7.387216143280, 7.387216143280, 7.387216143280, 7.387216143280, 7.387216143280, 7.387216143280, 7.387216143280, 7.408935392973, 7.408935392973, 7.408935392973, 7.408935392973, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.431798275933, 7.443697499233, 7.443697499233, 7.443697499233, 7.443697499233, 7.443697499233, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.468521082958, 7.468521082958, 7.468521082958, 7.468521082958, 7.468521082958, 7.468521082958, 7.468521082958, 7.468521082958, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.522878745280, 7.522878745280, 7.522878745280, 7.522878745280, 7.522878745280, 7.522878745280, 7.522878745280, 7.522878745280, 7.522878745280, 7.522878745280, 7.522878745280, 7.522878745280, 7.522878745280, 7.522878745280, 7.522878745280, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.552841968658, 7.552841968658, 7.552841968658, 7.552841968658, 7.552841968658, 7.552841968658, 7.552841968658, 7.552841968658, 7.552841968658, 7.552841968658, 7.552841968658, 7.552841968658, 7.552841968658, 7.552841968658, 7.552841968658, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.602059991328, 7.602059991328, 7.602059991328, 7.602059991328, 7.602059991328, 7.602059991328, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.795880017344, 7.795880017344, 7.795880017344, 7.795880017344, 7.795880017344, 7.795880017344, 7.795880017344, 7.795880017344, 7.795880017344, 7.795880017344, 7.795880017344, 7.795880017344, 7.795880017344, 7.795880017344, 7.795880017344, 7.795880017344, 7.795880017344, 7.795880017344, 7.795880017344, 7.795880017344, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.853871964322, 7.853871964322, 7.853871964322, 7.853871964322, 7.853871964322, 7.853871964322, 7.853871964322, 7.853871964322, 7.853871964322, 7.853871964322, 7.853871964322, 7.853871964322, 7.853871964322, 7.853871964322, 7.853871964322, 7.853871964322, 7.853871964322, 7.853871964322, 7.853871964322, 7.853871964322, 7.853871964322, 7.853871964322, 7.853871964322, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.958607314842, 7.958607314842, 7.958607314842, 7.958607314842, 7.958607314842, 7.958607314842, 7.958607314842, 7.958607314842, 7.958607314842, 7.958607314842, 7.958607314842, 7.958607314842, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000};
    return neg_log10_eccdf;
  }
  else if( slide ==  15 && lambda == 50){
    static double neg_log10_eccdf[] = {0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000434, 0.000000000434, 0.000000000434, 0.000000000869, 0.000000000869, 0.000000000869, 0.000000000869, 0.000000000869, 0.000000000869, 0.000000001303, 0.000000001303, 0.000000001303, 0.000000001303, 0.000000001303, 0.000000001303, 0.000000001303, 0.000000001303, 0.000000002172, 0.000000002172, 0.000000002606, 0.000000003040, 0.000000003040, 0.000000003909, 0.000000004343, 0.000000004777, 0.000000005212, 0.000000005646, 0.000000005646, 0.000000006514, 0.000000006514, 0.000000006514, 0.000000006514, 0.000000007383, 0.000000008252, 0.000000008252, 0.000000009555, 0.000000010857, 0.000000011726, 0.000000012595, 0.000000013029, 0.000000013897, 0.000000014766, 0.000000015635, 0.000000016069, 0.000000016938, 0.000000017372, 0.000000017806, 0.000000017806, 0.000000019109, 0.000000020846, 0.000000021715, 0.000000023886, 0.000000026058, 0.000000027361, 0.000000028229, 0.000000030401, 0.000000032138, 0.000000033875, 0.000000035612, 0.000000037349, 0.000000039955, 0.000000043429, 0.000000045601, 0.000000048641, 0.000000052550, 0.000000053418, 0.000000056893, 0.000000059064, 0.000000062973, 0.000000066447, 0.000000069921, 0.000000072527, 0.000000077304, 0.000000080779, 0.000000085556, 0.000000089899, 0.000000098151, 0.000000102928, 0.000000106836, 0.000000111614, 0.000000120300, 0.000000127683, 0.000000133763, 0.000000144620, 0.000000153740, 0.000000162426, 0.000000172849, 0.000000181101, 0.000000188050, 0.000000198907, 0.000000207159, 0.000000221056, 0.000000231045, 0.000000243205, 0.000000257971, 0.000000270566, 0.000000288372, 0.000000300966, 0.000000316167, 0.000000334407, 0.000000350910, 0.000000367848, 0.000000386522, 0.000000410409, 0.000000424740, 0.000000441678, 0.000000469473, 0.000000489450, 0.000000509428, 0.000000531143, 0.000000553726, 0.000000576743, 0.000000605407, 0.000000625819, 0.000000648402, 0.000000674894, 0.000000704861, 0.000000733090, 0.000000765662, 0.000000812131, 0.000000844704, 0.000000882053, 0.000000916797, 0.000000959792, 0.000000994970, 0.000001042742, 0.000001081829, 0.000001125258, 0.000001173031, 0.000001212118, 0.000001271182, 0.000001319823, 0.000001368464, 0.000001425357, 0.000001490501, 0.000001550868, 0.000001619053, 0.000001693318, 0.000001760633, 0.000001834898, 0.000001900043, 0.000001972136, 0.000002052481, 0.000002127614, 0.000002209262, 0.000002289607, 0.000002370820, 0.000002463760, 0.000002555396, 0.000002660062, 0.000002762990, 0.000002860273, 0.000002975362, 0.000003089582, 0.000003206842, 0.000003335829, 0.000003462209, 0.000003601619, 0.000003728000, 0.000003873055, 0.000004025494, 0.000004180539, 0.000004336452, 0.000004492365, 0.000004646107, 0.000004820696, 0.000004977043, 0.000005148157, 0.000005327523, 0.000005511232, 0.000005714050, 0.000005914697, 0.000006124898, 0.000006325980, 0.000006564845, 0.000006788944, 0.000006995672, 0.000007234538, 0.000007483393, 0.000007724865, 0.000007985446, 0.000008255582, 0.000008523981, 0.000008801935, 0.000009094222, 0.000009373913, 0.000009673149, 0.000010000180, 0.000010312879, 0.000010642951, 0.000010979538, 0.000011327416, 0.000011694405, 0.000012050102, 0.000012420566, 0.000012780606, 0.000013163231, 0.000013551068, 0.000013961055, 0.000014375385, 0.000014797099, 0.000015246175, 0.000015708281, 0.000016159530, 0.000016606436, 0.000017106328, 0.000017618816, 0.000018136516, 0.000018667247, 0.000019214047, 0.000019807320, 0.000020377575, 0.000020997777, 0.000021599739, 0.000022256425, 0.000022861429, 0.000023489452, 0.000024123991, 0.000024813256, 0.000025513379, 0.000026210463, 0.000026931436, 0.000027672823, 0.000028445049, 0.000029227699, 0.000030027724, 0.000030844255, 0.000031650364, 0.000032490351, 0.000033355966, 0.000034237218, 0.000035150612, 0.000036052281, 0.000036986527, 0.000037966380, 0.000038929296, 0.000039919578, 0.000040966325, 0.000041982237, 0.000043131058, 0.000044255994, 0.000045391791, 0.000046529329, 0.000047707698, 0.000048925161, 0.000050157395, 0.000051400925, 0.000052647065, 0.000053973998, 0.000055320480, 0.000056682169, 0.000058088166, 0.000059555411, 0.000060993560, 0.000062430410, 0.000063903751, 0.000065379704, 0.000066907350, 0.000068545764, 0.000070195478, 0.000071866048, 0.000073513168, 0.000075209813, 0.000076889090, 0.000078667410, 0.000080486134, 0.000082356991, 0.000084260434, 0.000086131308, 0.000088109047, 0.000090048135, 0.000092071069, 0.000094045360, 0.000096173001, 0.000098380145, 0.000100524314, 0.000102737128, 0.000104972541, 0.000107275734, 0.000109647574, 0.000112078940, 0.000114529435, 0.000117039458, 0.000119559052, 0.000122084308, 0.000124709930, 0.000127372493, 0.000130050278, 0.000132836252, 0.000135651784, 0.000138454302, 0.000141308536, 0.000144168871, 0.000147197788, 0.000150209783, 0.000153307819, 0.000156446716, 0.000159634294, 0.000162912263, 0.000166261943, 0.000169662915, 0.000173161235, 0.000176644376, 0.000180127111, 0.000183678956, 0.000187282967, 0.000190950444, 0.000194642282, 0.000198538364, 0.000202389728, 0.000206354096, 0.000210346308, 0.000214404603, 0.000218516380, 0.000222654702, 0.000226875188, 0.000231166543, 0.000235586996, 0.000240035303, 0.000244530587, 0.000249174531, 0.000253805923, 0.000258529924, 0.000263232683, 0.000268081506, 0.000273057276, 0.000278074823, 0.000283126325, 0.000288239162, 0.000293465486, 0.000298862669, 0.000304146488, 0.000309556841, 0.000315069395, 0.000320741523, 0.000326492826, 0.000332228993, 0.000338103884, 0.000344068824, 0.000350072530, 0.000356075015, 0.000362320990, 0.000368610522, 0.000374870152, 0.000381315914, 0.000387926951, 0.000394519398, 0.000401264957, 0.000407970628, 0.000414805077, 0.000421766135, 0.000428855114, 0.000436004637, 0.000443326433, 0.000450754865, 0.000458101691, 0.000465692107, 0.000473460911, 0.000481215071, 0.000489089804, 0.000497032073, 0.000505045795, 0.000513286634, 0.000521475887, 0.000529830962, 0.000538389699, 0.000546951215, 0.000555657269, 0.000564441773, 0.000573238196, 0.000582094810, 0.000591160784, 0.000600314797, 0.000609709505, 0.000619190529, 0.000628707424, 0.000638380234, 0.000648159385, 0.000658013568, 0.000667921911, 0.000678061450, 0.000688303011, 0.000698686622, 0.000709156613, 0.000719851758, 0.000730639827, 0.000741511692, 0.000752556977, 0.000763757424, 0.000774870711, 0.000786388038, 0.000797843887, 0.000809515851, 0.000821279502, 0.000833038251, 0.000845083640, 0.000857322136, 0.000869638436, 0.000882130897, 0.000894759497, 0.000907484643, 0.000920395565, 0.000933460071, 0.000946562401, 0.000959700817, 0.000973026361, 0.000986724044, 0.001000491807, 0.001014418024, 0.001028428270, 0.001042530824, 0.001056834969, 0.001071269324, 0.001085997171, 0.001100791698, 0.001115649428, 0.001130848606, 0.001146213783, 0.001161643517, 0.001177332469, 0.001193065538, 0.001209019977, 0.001225268109, 0.001241405350, 0.001257866373, 0.001274355460, 0.001290883070, 0.001307756306, 0.001324807061, 0.001342145571, 0.001359477803, 0.001377026384, 0.001394786986, 0.001412782725, 0.001430880735, 0.001449152925, 0.001467613264, 0.001486394251, 0.001505169079, 0.001524118610, 0.001543345291, 0.001562553646, 0.001582150347, 0.001601723958, 0.001621552592, 0.001641577868, 0.001661933216, 0.001682465219, 0.001703061412, 0.001723959002, 0.001744899020, 0.001766191953, 0.001787473283, 0.001809231449, 0.001830830209, 0.001852855727, 0.001875014523, 0.001897325809, 0.001919789173, 0.001942490142, 0.001965493756, 0.001988602862, 0.002012037376, 0.002035652488, 0.002059568227, 0.002083614452, 0.002107861011, 0.002132606455, 0.002157171744, 0.002182221608, 0.002207292559, 0.002232665729, 0.002258195796, 0.002284008525, 0.002310086937, 0.002336243327, 0.002362579014, 0.002389396657, 0.002416300243, 0.002443296775, 0.002470569279, 0.002498039181, 0.002525752824, 0.002553681420, 0.002581841177, 0.002610343114, 0.002638871391, 0.002667632700, 0.002696674718, 0.002726077480, 0.002755840621, 0.002785540675, 0.002815596828, 0.002846008277, 0.002876505357, 0.002907404628, 0.002938454322, 0.002969817576, 0.003001460348, 0.003033176710, 0.003065320076, 0.003097641220, 0.003130137120, 0.003162867309, 0.003195876029, 0.003229089399, 0.003262461960, 0.003296149533, 0.003330251537, 0.003364561486, 0.003398952934, 0.003433610188, 0.003468563081, 0.003503805992, 0.003539210264, 0.003575073260, 0.003610976878, 0.003647413393, 0.003683927863, 0.003720829993, 0.003757568794, 0.003794581890, 0.003832049867, 0.003870101667, 0.003908070910, 0.003946458504, 0.003985333372, 0.004024368646, 0.004063712100, 0.004103390124, 0.004143135349, 0.004183142049, 0.004223234002, 0.004263931388, 0.004304872074, 0.004345907862, 0.004387694641, 0.004429526244, 0.004471430328, 0.004513601754, 0.004556045429, 0.004598634144, 0.004641983365, 0.004685485298, 0.004729218138, 0.004773097210, 0.004817475185, 0.004862187084, 0.004907041943, 0.004952222546, 0.004997709663, 0.005043669470, 0.005089505833, 0.005135886253, 0.005182606531, 0.005229723919, 0.005277042483, 0.005324753097, 0.005372694079, 0.005421016340, 0.005469718247, 0.005518611671, 0.005567839204, 0.005617431756, 0.005667102131, 0.005717098543, 0.005767623991, 0.005818590628, 0.005869681221, 0.005920792354, 0.005972487187, 0.006024444024, 0.006076779227, 0.006129133070, 0.006182038149, 0.006235199920, 0.006288573092, 0.006342568458, 0.006396459473, 0.006451192479, 0.006505870224, 0.006560941527, 0.006616221334, 0.006671868931, 0.006728117786, 0.006784429952, 0.006841177374, 0.006898209315, 0.006955780080, 0.007013371718, 0.007071502909, 0.007129569038, 0.007188172783, 0.007247517054, 0.007306846470, 0.007366771484, 0.007426966747, 0.007487325884, 0.007548107479, 0.007609199890, 0.007670447651, 0.007732252692, 0.007794016123, 0.007856475111, 0.007919465864, 0.007982438330, 0.008046168494, 0.008109745179, 0.008173785224, 0.008238446403, 0.008303220797, 0.008368204530, 0.008433695707, 0.008499520949, 0.008565807087, 0.008632414781, 0.008699118208, 0.008766516991, 0.008834359703, 0.008901989155, 0.008970078260, 0.009038747404, 0.009107676209, 0.009176975691, 0.009246611420, 0.009316981123, 0.009387640039, 0.009458377032, 0.009529640506, 0.009600824655, 0.009672830467, 0.009745234629, 0.009818197261, 0.009891340098, 0.009964696112, 0.010038317855, 0.010112513072, 0.010187202470, 0.010262357843, 0.010337667657, 0.010413582161, 0.010489581123, 0.010566058846, 0.010642877164, 0.010719845288, 0.010797331942, 0.010874863595, 0.010953100777, 0.011031669228, 0.011110513878, 0.011189816247, 0.011269396068, 0.011349166116, 0.011429736452, 0.011510300332, 0.011591746265, 0.011673016533, 0.011755069493, 0.011836930432, 0.011919428162, 0.012002363029, 0.012085635252, 0.012169816271, 0.012254287451, 0.012339192831, 0.012424252459, 0.012509760091, 0.012595689611, 0.012682186595, 0.012768802335, 0.012855725232, 0.012943196632, 0.013030722349, 0.013118927730, 0.013207558887, 0.013296383673, 0.013386147488, 0.013476281526, 0.013566569593, 0.013657364939, 0.013748311059, 0.013839588285, 0.013930875282, 0.014022803604, 0.014115037173, 0.014207808177, 0.014300965250, 0.014394686414, 0.014488655912, 0.014583133947, 0.014677796056, 0.014773268710, 0.014868967742, 0.014964860033, 0.015061129150, 0.015157702486, 0.015254658507, 0.015352206678, 0.015449917621, 0.015548419365, 0.015646962365, 0.015746468468, 0.015846014041, 0.015946276713, 0.016047103758, 0.016147993432, 0.016249259566, 0.016350722484, 0.016452856839, 0.016555754686, 0.016658510848, 0.016761992327, 0.016865635271, 0.016970055321, 0.017074659937, 0.017179517922, 0.017284599175, 0.017390117214, 0.017496330074, 0.017603019366, 0.017710241529, 0.017817602386, 0.017925637454, 0.018033696686, 0.018141940847, 0.018251173119, 0.018361067594, 0.018471374620, 0.018581765423, 0.018692519805, 0.018803785428, 0.018915373062, 0.019026863234, 0.019139173114, 0.019251946504, 0.019365242335, 0.019478756684, 0.019592526516, 0.019707377314, 0.019822319406, 0.019937784818, 0.020053237275, 0.020169605289, 0.020286125993, 0.020403147718, 0.020520462329, 0.020638748194, 0.020757094528, 0.020875518688, 0.020994517557, 0.021114523820, 0.021234508981, 0.021354612151, 0.021475537960, 0.021596705128, 0.021718034405, 0.021840133777, 0.021962760853, 0.022085732376, 0.022209546849, 0.022333467956, 0.022457605998, 0.022582575080, 0.022707700484, 0.022833192892, 0.022958925784, 0.023085398123, 0.023212051159, 0.023339359351, 0.023466852479, 0.023594843880, 0.023723171153, 0.023851806138, 0.023981282397, 0.024110752737, 0.024240774640, 0.024371029332, 0.024502022461, 0.024633402601, 0.024765018363, 0.024897250762, 0.025029604866, 0.025163039370, 0.025296491406, 0.025430391991, 0.025565056135, 0.025700076294, 0.025835337096, 0.025971045308, 0.026107323555, 0.026243980900, 0.026380989999, 0.026518215911, 0.026656460026, 0.026794946330, 0.026933373302, 0.027072623727, 0.027211886448, 0.027351404752, 0.027492229606, 0.027632697017, 0.027774154809, 0.027916240842, 0.028058267748, 0.028201165096, 0.028344137757, 0.028487443163, 0.028631109460, 0.028774968081, 0.028919522855, 0.029064969027, 0.029211020408, 0.029357384852, 0.029503834897, 0.029650604962, 0.029797748846, 0.029945867611, 0.030093943352, 0.030242542198, 0.030391359584, 0.030540988387, 0.030690935367, 0.030841352369, 0.030991941822, 0.031142954908, 0.031295058561, 0.031447371451, 0.031600124515, 0.031753351371, 0.031906439270, 0.032060174494, 0.032214484278, 0.032369275635, 0.032524426425, 0.032680234346, 0.032836221843, 0.032993180978, 0.033150091862, 0.033307659197, 0.033465757964, 0.033624214158, 0.033782902815, 0.033941996509, 0.034101682164, 0.034261490938, 0.034421930202, 0.034582919271, 0.034744362761, 0.034906176911, 0.035068287724, 0.035231116108, 0.035394042786, 0.035557431727, 0.035721761556, 0.035886219155, 0.036051027375, 0.036215955281, 0.036382178400, 0.036548445321, 0.036714933771, 0.036881542805, 0.037049027867, 0.037217040762, 0.037385645502, 0.037554797296, 0.037723639823, 0.037893259325, 0.038063632510, 0.038234194918, 0.038405402691, 0.038577175879, 0.038748916845, 0.038921225941, 0.039093325912, 0.039266683060, 0.039440371004, 0.039613940509, 0.039788368080, 0.039963394742, 0.040139463241, 0.040315734199, 0.040492460979, 0.040668910591, 0.040846137586, 0.041023618793, 0.041201798102, 0.041379818376, 0.041558736040, 0.041737428629, 0.041917067671, 0.042097750010, 0.042278455852, 0.042459342272, 0.042641111824, 0.042823205764, 0.043005658991, 0.043188633107, 0.043371920458, 0.043555259270, 0.043739192807, 0.043923877961, 0.044108131642, 0.044293668725, 0.044479348135, 0.044664872072, 0.044852238227, 0.045038737164, 0.045226552930, 0.045414284089, 0.045603034175, 0.045791523692, 0.045980570125, 0.046170175613, 0.046359796748, 0.046550203136, 0.046741635178, 0.046932836163, 0.047124724031, 0.047316778563, 0.047509339579, 0.047701754803, 0.047894950219, 0.048089151454, 0.048283625948, 0.048478522617, 0.048674159219, 0.048869611818, 0.049065688740, 0.049262851470, 0.049460115917, 0.049657875595, 0.049855741928, 0.050053633633, 0.050252666264, 0.050452355016, 0.050652138556, 0.050851601965, 0.051052267313, 0.051253899208, 0.051455085478, 0.051656454993, 0.051858586409, 0.052061335822, 0.052264641845, 0.052469136420, 0.052672869313, 0.052877415954, 0.053081917150, 0.053287155426, 0.053493231011, 0.053699507140, 0.053906522453, 0.054113941483, 0.054321857782, 0.054529054322, 0.054736971459, 0.054945321361, 0.055154033597, 0.055363175242, 0.055573464040, 0.055784074460, 0.055994864633, 0.056205880739, 0.056417697716, 0.056629340974, 0.056841801740, 0.057054426422, 0.057268265195, 0.057481730907, 0.057696046583, 0.057910842736, 0.058126101160, 0.058342172192, 0.058558255357, 0.058774804592, 0.058991558368, 0.059208893222, 0.059426168574, 0.059644433052, 0.059863149734, 0.060081796579, 0.060300435798, 0.060519615482, 0.060739732487, 0.060959810691, 0.061180286472, 0.061401630614, 0.061623971013, 0.061846279078, 0.062068879359, 0.062291807337, 0.062515046908, 0.062738766090, 0.062962574471, 0.063188329168, 0.063413968586, 0.063640149689, 0.063865675894, 0.064091758704, 0.064318601888, 0.064546297777, 0.064773186989, 0.065000899963, 0.065229124530, 0.065458364862, 0.065687273089, 0.065916699336, 0.066146733373, 0.066377809250, 0.066608772719, 0.066840003956, 0.067072001025, 0.067304442068, 0.067537242494, 0.067770882533, 0.068004360854, 0.068238078595, 0.068472154546, 0.068707188000, 0.068942742183, 0.069178222010, 0.069414292780, 0.069650841183, 0.069887854155, 0.070124859760, 0.070362312927, 0.070599961908, 0.070838585056, 0.071077551156, 0.071316886805, 0.071556836292, 0.071796897932, 0.072037068760, 0.072279110304, 0.072520449753, 0.072761928525, 0.073003963978, 0.073246609404, 0.073489072596, 0.073732665543, 0.073976578001, 0.074220997979, 0.074465659214, 0.074711027818, 0.074956545452, 0.075202369275, 0.075448143211, 0.075695192132, 0.075942539955, 0.076189779261, 0.076437392958, 0.076685091950, 0.076933903379, 0.077183337290, 0.077432743774, 0.077683200368, 0.077933873716, 0.078184663757, 0.078436234983, 0.078687954101, 0.078939774408, 0.079191491175, 0.079444620049, 0.079698007160, 0.079951690988, 0.080206301439, 0.080460297595, 0.080714344064, 0.080969271328, 0.081224868790, 0.081480503078, 0.081737155001, 0.081993849580, 0.082250889642, 0.082508663278, 0.082766288909, 0.083024822484, 0.083283878324, 0.083542993501, 0.083802015806, 0.084061739233, 0.084321870154, 0.084582366466, 0.084842987971, 0.085105294370, 0.085367058859, 0.085629470495, 0.085892889731, 0.086156333794, 0.086419952067, 0.086684302578, 0.086949325553, 0.087213995929, 0.087479213887, 0.087745016476, 0.088011182888, 0.088278282596, 0.088544412959, 0.088811338983, 0.089079182530, 0.089347086268, 0.089614806248, 0.089883435978, 0.090152461264, 0.090421964032, 0.090692092233, 0.090962693248, 0.091232957697, 0.091504022026, 0.091775679470, 0.092047991162, 0.092320388297, 0.092592443619, 0.092865094361, 0.093139044394, 0.093412510892, 0.093686717109, 0.093961089016, 0.094235684007, 0.094511022840, 0.094787087370, 0.095063304252, 0.095339772109, 0.095615840724, 0.095893052774, 0.096170824495, 0.096448446990, 0.096727513588, 0.097005935908, 0.097284702012, 0.097564145182, 0.097843271568, 0.098123152340, 0.098404312334, 0.098685985884, 0.098966938981, 0.099248737641, 0.099530582730, 0.099813014172, 0.100096100802, 0.100379195325, 0.100662736793, 0.100947004863, 0.101231637319, 0.101516586475, 0.101801439421, 0.102087472600, 0.102373059888, 0.102660138843, 0.102946888598, 0.103233560106, 0.103521211380, 0.103808640737, 0.104095937573, 0.104383567079, 0.104673187214, 0.104961751883, 0.105250861480, 0.105540792188, 0.105831300041, 0.106122679957, 0.106414125103, 0.106706048585, 0.106997498340, 0.107289519670, 0.107582114709, 0.107878141175, 0.108171708143, 0.108465694458, 0.108759782844, 0.109054066603, 0.109348798503, 0.109644536666, 0.109939702140, 0.110236063093, 0.110532246627, 0.110829630629, 0.111126816777, 0.111424656043, 0.111722116099, 0.112020708602, 0.112319510476, 0.112618834402, 0.112918290984, 0.113217729862, 0.113517393954, 0.113818400005, 0.114119227365, 0.114421116705, 0.114723364793, 0.115025090420, 0.115327096614, 0.115630286128, 0.115933407834, 0.116236824690, 0.116540649055, 0.116844781036, 0.117150114094, 0.117454601612, 0.117759862072, 0.118066052556, 0.118371590963, 0.118677927234, 0.118984700209, 0.119290938980, 0.119598494909, 0.119906379839, 0.120215238710, 0.120523508027, 0.120832439141, 0.121141632644, 0.121451402578, 0.121761011927, 0.122071642325, 0.122381886588, 0.122692141219, 0.123003126686, 0.123314914782, 0.123627161844, 0.123939325636, 0.124251781610, 0.124565116319, 0.124878541195, 0.125192660583, 0.125506206607, 0.125820289010, 0.126135164281, 0.126450993197, 0.126766480924, 0.127082790398, 0.127398838916, 0.127715178207, 0.128031468151, 0.128349254493, 0.128666505511, 0.128984332123, 0.129302989273, 0.129621221900, 0.129939963782, 0.130258773291, 0.130578809579, 0.130898298740, 0.131218251054, 0.131538746742, 0.131859069032, 0.132180751768, 0.132502260489, 0.132824421346, 0.133146955906, 0.133469603802, 0.133793077836, 0.134116248870, 0.134440270779, 0.134764360478, 0.135088664941, 0.135413185062, 0.135739101299, 0.136065419143, 0.136391458574, 0.136717701904, 0.137045206293, 0.137372563370, 0.137700616419, 0.138028177453, 0.138356244928, 0.138684806659, 0.139012758828, 0.139341324568, 0.139671234715, 0.140001182247, 0.140330837687, 0.140661124202, 0.140991445785, 0.141322469791, 0.141653858832, 0.141985029379, 0.142316125990, 0.142648177894, 0.142981568582, 0.143314475996, 0.143648227626, 0.143981257630, 0.144314282851, 0.144648045364, 0.144982907487, 0.145318365428, 0.145652964557, 0.145988840982, 0.146323914692, 0.146660817109, 0.146997402310, 0.147334043712, 0.147670905991, 0.148008603162, 0.148346651125, 0.148684797291, 0.149023719870, 0.149362581268, 0.149701438763, 0.150041275602, 0.150381605137, 0.150721404599, 0.151062022982, 0.151402951804, 0.151745201705, 0.152087496542, 0.152430627690, 0.152772557688, 0.153114417895, 0.153457475572, 0.153800900992, 0.154145312914, 0.154488903562, 0.154833852441, 0.155178930872, 0.155523395237, 0.155868517308, 0.156214208204, 0.156560101615, 0.156906553084, 0.157253202552, 0.157600122066, 0.157947313926, 0.158293557383, 0.158640665340, 0.158988964696, 0.159337462126, 0.159686401272, 0.160035428269, 0.160384349586, 0.160733398030, 0.161083622119, 0.161433312615, 0.161784293434, 0.162134956941, 0.162486031978, 0.162836920310, 0.163187759113, 0.163540742876, 0.163893984000, 0.164246337393, 0.164599417186, 0.164952936017, 0.165306386989, 0.165660380881, 0.166015744329, 0.166371104487, 0.166727299483, 0.167081676767, 0.167436946920, 0.167793789357, 0.168150652112, 0.168507519630, 0.168864122592, 0.169221128253, 0.169578774169, 0.169936980755, 0.170295162274, 0.170653982987, 0.171012956825, 0.171371897046, 0.171731486309, 0.172091000471, 0.172450934584, 0.172811102597, 0.173172036742, 0.173532305507, 0.173894088668, 0.174256185136, 0.174617663140, 0.174980397438, 0.175342715699, 0.175704729229, 0.176068529980, 0.176432346269, 0.176796755993, 0.177160814126, 0.177525453512, 0.177890282874, 0.178254704696, 0.178619525611, 0.178985496599, 0.179351384402, 0.179717946616, 0.180084030850, 0.180450635811, 0.180818448790, 0.181185209803, 0.181556836616, 0.181924443836, 0.182292626804, 0.182660788179, 0.183030033070, 0.183399820065, 0.183768962740, 0.184138659679, 0.184508494255, 0.184878386116, 0.185249838107, 0.185620777037, 0.185992136365, 0.186363326575, 0.186734798260, 0.187106439619, 0.187478360490, 0.187850937477, 0.188223848437, 0.188597454000, 0.188971470486, 0.189345872489, 0.189719821137, 0.190094250829, 0.190468935606, 0.190844257980, 0.191218514105, 0.191593689152, 0.191969108850, 0.192345111019, 0.192721314461, 0.193097482305, 0.193474707864, 0.193851953941, 0.194229252230, 0.194605737177, 0.194983399283, 0.195362173238, 0.195740950668, 0.196119564839, 0.196498622030, 0.196877620829, 0.197257258450, 0.197637739599, 0.198018329639, 0.198400039400, 0.198781114403, 0.199162390803, 0.199543831018, 0.199924802050, 0.200306000118, 0.200688555456, 0.201071882088, 0.201454601209, 0.201838516404, 0.202221731441, 0.202605285613, 0.202988873878, 0.203374367570, 0.203758942903, 0.204143716629, 0.204528007361, 0.204913869216, 0.205299656846, 0.205685677294, 0.206072468336, 0.206459390390, 0.206846566566, 0.207233938447, 0.207622218655, 0.208010504867, 0.208398900648, 0.208787306946, 0.209176071467, 0.209565037951, 0.209953965768, 0.210342628935, 0.210732247703, 0.211122254456, 0.211512922752, 0.211903574213, 0.212294818131, 0.212686053384, 0.213077668355, 0.213470229619, 0.213862932140, 0.214254869102, 0.214646581303, 0.215039287010, 0.215432336013, 0.215827046702, 0.216220200864, 0.216613804943, 0.217007753904, 0.217402030450, 0.217796653077, 0.218191671927, 0.218586959881, 0.218982797060, 0.219378458444, 0.219773928814, 0.220170290309, 0.220566857269, 0.220964018551, 0.221360838432, 0.221758698583, 0.222156246294, 0.222553798649, 0.222951600601, 0.223349921245, 0.223749059016, 0.224147832647, 0.224546782677, 0.224947475161, 0.225346412851, 0.225746613476, 0.226147632073, 0.226547192789, 0.226947738838, 0.227348678106, 0.227749309165, 0.228150642828, 0.228553566477, 0.228954981466, 0.229358349712, 0.229761304218, 0.230164174763, 0.230566910547, 0.230969247699, 0.231372416649, 0.231776300878, 0.232180923532, 0.232586879880, 0.232991190880, 0.233396123926, 0.233801303934, 0.234207073067, 0.234613143403, 0.235019580345, 0.235426253875, 0.235832657490, 0.236239897425, 0.236646763188, 0.237054294581, 0.237461382650, 0.237869593923, 0.238278526027, 0.238686697579, 0.239095608598, 0.239504097610, 0.239912166834, 0.240322365290, 0.240732776161, 0.241143404225, 0.241553156011, 0.241964896705, 0.242375690260, 0.242786902372, 0.243196958529, 0.243608721820, 0.244020358664, 0.244432592665, 0.244844123164, 0.245256243373, 0.245669111350, 0.246081344320, 0.246495098168, 0.246909304855, 0.247323431084, 0.247738407410, 0.248151568166, 0.248565685808, 0.248979829643, 0.249395096817, 0.249810289004, 0.250226250942, 0.250641452279, 0.251057445002, 0.251473538223, 0.251889419258, 0.252305993182, 0.252722212449, 0.253139223057, 0.253556996536, 0.253974741275, 0.254392540281, 0.254810237924, 0.255228881719, 0.255648624237, 0.256066633225, 0.256485657173, 0.256904618922, 0.257324432927, 0.257743581562, 0.258164146352, 0.258584802931, 0.259006267900, 0.259426962360, 0.259847552809, 0.260269003296, 0.260690875039, 0.261113290886, 0.261536466161, 0.261959683445, 0.262382665928, 0.262805651937, 0.263228494584, 0.263651220942, 0.264075362718, 0.264497949988, 0.264922187757, 0.265346769139, 0.265771749194, 0.266196452101, 0.266620992188, 0.267046259332, 0.267472689424, 0.267898148762, 0.268323722415, 0.268748663614, 0.269175931555, 0.269602468929, 0.270028833630, 0.270461348114, 0.270888289036, 0.271315711729, 0.271744693037, 0.272173495434, 0.272602140752, 0.273029910658, 0.273458425125, 0.273887210236, 0.274316392975, 0.274745710024, 0.275175280027, 0.275605184417, 0.276035685350, 0.276465329677, 0.276895943393, 0.277326892382, 0.277759397308, 0.278190564883, 0.278620552377, 0.279052314426, 0.279482878673, 0.279915434765, 0.280346858480, 0.280779549330, 0.281211669252, 0.281644876618, 0.282078242975, 0.282512148743, 0.282945826930, 0.283379345641, 0.283813706596, 0.284248373720, 0.284682831376, 0.285117547410, 0.285552890961, 0.285988110879, 0.286423433117, 0.286858365734, 0.287293610651, 0.287729809432, 0.288167263021, 0.288603631575, 0.289041539979, 0.289478182727, 0.289915496054, 0.290352523027, 0.290790708772, 0.291227928243, 0.291665714981, 0.292104779074, 0.292544424640, 0.292983462702, 0.293423034666, 0.293861594436, 0.294300967838, 0.294741268187, 0.295182190200, 0.295622029072, 0.296062193630, 0.296502962917, 0.296943663727, 0.297384624423, 0.297825845343, 0.298267430394, 0.298708631180, 0.299151500061, 0.299594337077, 0.300037362653, 0.300479926287, 0.300922404717, 0.301365563908, 0.301808691963, 0.302252570516, 0.302696666915, 0.303141585356, 0.303585988494, 0.304030360560, 0.304474558288, 0.304918970673, 0.305363683885, 0.305809739928, 0.306255553933, 0.306701150202, 0.307146249227, 0.307592834851, 0.308039157223, 0.308486129615, 0.308932094165, 0.309378442757, 0.309826057146, 0.310273705687, 0.310720584198, 0.311167794980, 0.311616276667, 0.312064419226, 0.312512716112, 0.312960716469, 0.313409595908, 0.313857795552, 0.314307551693, 0.314756579973, 0.315205911461, 0.315656673136, 0.316107004777, 0.316557493292, 0.317009392174, 0.317461045510, 0.317911075844, 0.318361061362, 0.318812614860, 0.319264033247, 0.319715636613, 0.320167259840, 0.320620219158, 0.321071198169, 0.321522130597, 0.321974704825, 0.322426995654, 0.322880123371, 0.323333500340, 0.323787884732, 0.324241368779, 0.324695670805, 0.325151564168, 0.325606012798, 0.326061403046, 0.326516592464, 0.326972112865, 0.327427367596, 0.327883064019, 0.328338887612, 0.328795240128, 0.329251632399, 0.329709215547, 0.330165496053, 0.330622045371, 0.331078599492, 0.331536277943, 0.331993054982, 0.332451730399, 0.332910054153, 0.333367871158, 0.333825813458, 0.334284113492, 0.334742028434, 0.335201196313, 0.335661484195, 0.336120510780, 0.336579052087, 0.337039435054, 0.337499873915, 0.337958814577, 0.338419264099, 0.338879949286, 0.339342243131, 0.339804696251, 0.340264913781, 0.340726131542, 0.341187888224, 0.341650109678, 0.342112020643, 0.342575640127, 0.343038334210, 0.343501046688, 0.343964357196, 0.344426908909, 0.344889884610, 0.345353984446, 0.345817626482, 0.346280283356, 0.346743478982, 0.347207540112, 0.347672686584, 0.348137429508, 0.348602237085, 0.349067772591, 0.349532816042, 0.349998744955, 0.350464231141, 0.350929705247, 0.351395666101, 0.351863596975, 0.352330633767, 0.352798674364, 0.353264880622, 0.353731885018, 0.354199356791, 0.354666933286, 0.355134701859, 0.355601974152, 0.356069801011, 0.356538061296, 0.357006883329, 0.357475352436, 0.357943131247, 0.358411698942, 0.358880217014, 0.359348552629, 0.359817959701, 0.360288116605, 0.360758726232, 0.361228191139, 0.361697758562, 0.362167162314, 0.362637228116, 0.363108631954, 0.363578454469, 0.364048952493, 0.364519910523, 0.364990427760, 0.365461910702, 0.365934048269, 0.366405975721, 0.366877848486, 0.367349923860, 0.367821714697, 0.368294661553, 0.368767448901, 0.369241003545, 0.369714413797, 0.370188017787, 0.370666373395, 0.371139242112, 0.371613632806, 0.372087691125, 0.372561585416, 0.373037094372, 0.373511955517, 0.373987221370, 0.374462157243, 0.374938202069, 0.375414728042, 0.375890737207, 0.376366593014, 0.376842297489, 0.377318092877, 0.377794280533, 0.378270918249, 0.378746962939, 0.379224121736, 0.379700618528, 0.380178202539, 0.380655737418, 0.381133921208, 0.381611291534, 0.382088756912, 0.382566365284, 0.383046006028, 0.383525404065, 0.384005621022, 0.384485123264, 0.384963914158, 0.385443812643, 0.385922035801, 0.386401427897, 0.386882248365, 0.387361394583, 0.387841896363, 0.388320972112, 0.388801018117, 0.389281913553, 0.389762623940, 0.390243275004, 0.390724850514, 0.391205173200, 0.391686832505, 0.392167874835, 0.392649337956, 0.393131087380, 0.393614204893, 0.394097753893, 0.394580892718, 0.395064974082, 0.395548424068, 0.396032653882, 0.396517649299, 0.397001972452, 0.397484979484, 0.397970439667, 0.398454807137, 0.398939828600, 0.399424483737, 0.399909066287, 0.400394771020, 0.400880574664, 0.401365093705, 0.401851261511, 0.402336565891, 0.402822847996, 0.403308407778, 0.403795372714, 0.404280413159, 0.404766190042, 0.405252444646, 0.405739567071, 0.406226631147, 0.406713500901, 0.407200854966, 0.407688324602, 0.408174572296, 0.408662080594, 0.409150517781, 0.409639932108, 0.410128661349, 0.410617726560, 0.411106829441, 0.411595829448, 0.412085998718, 0.412574867858, 0.413063380680, 0.413553154924, 0.414042899602, 0.414533367511, 0.415022605665, 0.415513513707, 0.416003128969, 0.416492980696, 0.416984671979, 0.417476069932, 0.417966353170, 0.418457517204, 0.418948561595, 0.419439557180, 0.419929948090, 0.420421463926, 0.420913339778, 0.421406015162, 0.421898730350, 0.422390112247, 0.422883236313, 0.423377021095, 0.423870139335, 0.424363031198, 0.424855959807, 0.425349278530, 0.425842101200, 0.426335003877, 0.426829218559, 0.427324932646, 0.427820160624, 0.428314504287, 0.428807456382, 0.429302496302, 0.429797091672, 0.430292210014, 0.430787492965, 0.431282505492, 0.431777783586, 0.432272775236, 0.432768481103, 0.433264182216, 0.433760558241, 0.434256758574, 0.434754295828, 0.435251988443, 0.435748908838, 0.436246366444, 0.436743643005, 0.437241206731, 0.437739080565, 0.438237602058, 0.438736711961, 0.439233414603, 0.439731999725, 0.440231169867, 0.440729600025, 0.441228328184, 0.441726941600, 0.442226998593, 0.442726578806, 0.443225235261, 0.443725064537, 0.444224926201, 0.444723303328, 0.445222783279, 0.445723009234, 0.446223138578, 0.446723009957, 0.447224300181, 0.447725222391, 0.448225613733, 0.448726587169, 0.449227710298, 0.449727802543, 0.450227772034, 0.450728806964, 0.451230954547, 0.451733346686, 0.452235618148, 0.452736862509, 0.453239150958, 0.453741645690, 0.454244146510, 0.454747423705, 0.455250218103, 0.455752305331, 0.456254894221, 0.456758558960, 0.457262778620, 0.457768337009, 0.458272333824, 0.458775519829, 0.459279524599, 0.459783580394, 0.460287662897, 0.460792523151, 0.461296468487, 0.461801621841, 0.462305931756, 0.462811283060, 0.463318376665, 0.463824171449, 0.464328743144, 0.464834286770, 0.465340562849, 0.465847040070, 0.466352409341, 0.466859544398, 0.467367320749, 0.467875805000, 0.468381818107, 0.468889139949, 0.469396265405, 0.469904152870, 0.470411145523, 0.470918624122, 0.471426315781, 0.471934797314, 0.472443944489, 0.472951321314, 0.473460226958, 0.473968667704, 0.474477146260, 0.474985391092, 0.475494073051, 0.476003172184, 0.476513502422, 0.477022968932, 0.477532034840, 0.478042067657, 0.478551113238, 0.479063160299, 0.479572592475, 0.480086993856, 0.480596596829, 0.481106756399, 0.481619149658, 0.482130309557, 0.482641116441, 0.483152900025, 0.483665229222, 0.484176476456, 0.484687949713, 0.485199403522, 0.485711097570, 0.486224681792, 0.486738358577, 0.487250049738, 0.487763603587, 0.488275299070, 0.488787537928, 0.489300263758, 0.489813862588, 0.490326514201, 0.490840161636, 0.491353224480, 0.491867476442, 0.492381612036, 0.492896767735, 0.493411535599, 0.493925424624, 0.494440125825, 0.494954554012, 0.495470398222, 0.495984633805, 0.496500241910, 0.497013372172, 0.497528775417, 0.498044601009, 0.498559906614, 0.499078916098, 0.499592765993, 0.500109532407, 0.500625603753, 0.501141699729, 0.501658761287, 0.502174824281, 0.502690380514, 0.503208260936, 0.503724410364, 0.504241307090, 0.504759701453, 0.505274998078, 0.505793797603, 0.506311169222, 0.506827750237, 0.507345920002, 0.507864763289, 0.508382983855, 0.508899987171, 0.509417550546, 0.509935255112, 0.510454271110, 0.510971361547, 0.511492395007, 0.512011134552, 0.512529736775, 0.513050411063, 0.513569705349, 0.514089632658, 0.514608402046, 0.515128737509, 0.515648704789, 0.516168718048, 0.516689113513, 0.517209157371, 0.517728532883, 0.518249628835, 0.518770467429, 0.519292395146, 0.519813027529, 0.520333376667, 0.520854495540, 0.521374155844, 0.521896076619, 0.522418845185, 0.522941519846, 0.523465267885, 0.523988027152, 0.524510886039, 0.525033150187, 0.525555476526, 0.526079906221, 0.526600826147, 0.527123878830, 0.527646277157, 0.528169314874, 0.528693222412, 0.529216091082, 0.529740847468, 0.530264464343, 0.530786075854, 0.531311316841, 0.531833200778, 0.532356436126, 0.532879486390, 0.533404625305, 0.533929672345, 0.534454662090, 0.534979440241, 0.535504410639, 0.536030760970, 0.536557513947, 0.537081870161, 0.537607171734, 0.538132958015, 0.538660109713, 0.539184738150, 0.539710317126, 0.540236444029, 0.540762314539, 0.541288384580, 0.541815068439, 0.542341042827, 0.542868447778, 0.543395906643, 0.543923119493, 0.544450831645, 0.544977459993, 0.545503481731, 0.546031530859, 0.546558730726, 0.547086946355, 0.547615558484, 0.548143541287, 0.548672515478, 0.549200590494, 0.549728496831, 0.550256814369, 0.550784711760, 0.551312738468, 0.551842560821, 0.552371869845, 0.552898743964, 0.553426404044, 0.553956674675, 0.554485354680, 0.555015081227, 0.555544711782, 0.556075220291, 0.556604331198, 0.557134121981, 0.557664161474, 0.558194325744, 0.558724317312, 0.559255090246, 0.559786219512, 0.560316211278, 0.560846788995, 0.561377115645, 0.561907960812, 0.562440830447, 0.562973424333, 0.563504367234, 0.564037413117, 0.564568903859, 0.565101506926, 0.565632081900, 0.566165145142, 0.566696781746, 0.567228835865, 0.567762979314, 0.568296073630, 0.568829060353, 0.569360944181, 0.569893659271, 0.570426117664, 0.570960410170, 0.571494740679, 0.572027915382, 0.572562204795, 0.573096877684, 0.573630540226, 0.574165908486, 0.574700420585, 0.575235586452, 0.575769898459, 0.576304808064, 0.576839329887, 0.577375079894, 0.577908177261, 0.578442680035, 0.578978955002, 0.579512930854, 0.580047412129, 0.580583530783, 0.581120365045, 0.581657414543, 0.582193150849, 0.582729738262, 0.583265703501, 0.583800245593, 0.584337374299, 0.584874602085, 0.585411413567, 0.585947120074, 0.586485259708, 0.587020633680, 0.587556697007, 0.588094936816, 0.588633135428, 0.589169661136, 0.589708410627, 0.590247088808, 0.590784558810, 0.591321938932, 0.591860982564, 0.592398443261, 0.592937047915, 0.593475222817, 0.594012856400, 0.594551052219, 0.595089678196, 0.595628271252, 0.596167524507, 0.596712817255, 0.597251790973, 0.597792249775, 0.598331771611, 0.598869731339, 0.599410404294, 0.599949778679, 0.600489636839, 0.601029543012, 0.601570077841, 0.602111288010, 0.602653102146, 0.603194337296, 0.603734319052, 0.604275828603, 0.604816500213, 0.605357420415, 0.605899900246, 0.606441926695, 0.606984026070, 0.607523191366, 0.608063993839, 0.608605666339, 0.609148382560, 0.609689260247, 0.610233534922, 0.610775908403, 0.611319300804, 0.611860923703, 0.612403824239, 0.612947719577, 0.613490001548, 0.614031882903, 0.614576059333, 0.615118648493, 0.615660489603, 0.616204066440, 0.616746038792, 0.617289802027, 0.617832951710, 0.618376855484, 0.618920953684, 0.619463781534, 0.620007862638, 0.620551810483, 0.621095845144, 0.621637941637, 0.622182866352, 0.622727096523, 0.623272750143, 0.623818744997, 0.624362300021, 0.624908387357, 0.625453325227, 0.625999085397, 0.626544850392, 0.627090796030, 0.627636826316, 0.628181501703, 0.628729092468, 0.629275129223, 0.629819057086, 0.630364694247, 0.630910876702, 0.631456901148, 0.632002684229, 0.632549511862, 0.633097461764, 0.633644816680, 0.634191212464, 0.634738245978, 0.635285538081, 0.635830903288, 0.636379522308, 0.636925733049, 0.637472098235, 0.638021006703, 0.638568843108, 0.639116826577, 0.639664197179, 0.640212010074, 0.640760729442, 0.641309030608, 0.641857288079, 0.642405710506, 0.642956248340, 0.643504308657, 0.644054534903, 0.644604514565, 0.645154257329, 0.645703617351, 0.646253348134, 0.646803169055, 0.647352314106, 0.647900567475, 0.648450648511, 0.648999694929, 0.649550633976, 0.650099630034, 0.650650816939, 0.651201733571, 0.651753248662, 0.652303616188, 0.652855121423, 0.653404434209, 0.653953675263, 0.654504601681, 0.655056918717, 0.655607148566, 0.656159032681, 0.656709349411, 0.657261528211, 0.657813463873, 0.658365680632, 0.658916388596, 0.659468513503, 0.660020942202, 0.660572909681, 0.661123985383, 0.661674964119, 0.662227081726, 0.662781080869, 0.663334407329, 0.663887524318, 0.664441894160, 0.664994548804, 0.665547443161, 0.666099530392, 0.666651870833, 0.667205958112, 0.667759334547, 0.668313133705, 0.668864746751, 0.669418690300, 0.669972203824, 0.670526069824, 0.671079624882, 0.671633509190, 0.672187043254, 0.672740721573, 0.673293798757, 0.673847909095, 0.674403273163, 0.674957162153, 0.675511540390, 0.676066340868, 0.676620450522, 0.677175061540, 0.677732294470, 0.678286506061, 0.678840101069, 0.679396314433, 0.679951763299, 0.680507942207, 0.681063857013, 0.681619257502, 0.682177533495, 0.682734745850, 0.683290795295, 0.683847545008, 0.684404112734, 0.684960181473, 0.685516085230, 0.686071592692, 0.686628531341, 0.687184663519, 0.687742788972, 0.688298229865, 0.688854040491, 0.689411345061, 0.689970299430, 0.690527418516, 0.691085044223, 0.691642763370, 0.692198781711, 0.692755666948, 0.693313022837, 0.693869740730, 0.694427854385, 0.694985829843, 0.695542771568, 0.696101043229, 0.696658283897, 0.697218424829, 0.697776030043, 0.698333333007, 0.698893251740, 0.699451923744, 0.700009963624, 0.700567581605, 0.701126342011, 0.701683691780, 0.702241576145, 0.702801538508, 0.703363256959, 0.703921918511, 0.704481776859, 0.705042703587, 0.705601002096, 0.706160067808, 0.706720318350, 0.707279836130, 0.707839432967, 0.708399694088, 0.708962027666, 0.709519904111, 0.710080935248, 0.710641690524, 0.711203054634, 0.711764080771, 0.712323069311, 0.712887231220, 0.713448236199, 0.714008294300, 0.714569507743, 0.715130730640, 0.715693478602, 0.716255185143, 0.716819257226, 0.717379830917, 0.717943148022, 0.718505002590, 0.719067123295, 0.719635854737, 0.720197031772, 0.720758813865, 0.721322759301, 0.721885027537, 0.722447850469, 0.723008729989, 0.723572345584, 0.724134261399, 0.724697568746, 0.725260583385, 0.725821815604, 0.726383674577, 0.726947919612, 0.727511154811, 0.728075576526, 0.728641413963, 0.729208583238, 0.729771633864, 0.730335487775, 0.730898667747, 0.731464299008, 0.732030218016, 0.732593944919, 0.733157483574, 0.733721281597, 0.734286439036, 0.734850582879, 0.735415989482, 0.735978919602, 0.736542270167, 0.737107058931, 0.737671799768, 0.738236130238, 0.738802044641, 0.739366371473, 0.739933162586, 0.740499960837, 0.741067138517, 0.741634207438, 0.742200343434, 0.742766447391, 0.743332342655, 0.743897350758, 0.744462335331, 0.745028577378, 0.745595169437, 0.746162588800, 0.746729417340, 0.747296163880, 0.747861626602, 0.748428427580, 0.748996639319, 0.749560413403, 0.750128645764, 0.750697196948, 0.751264842560, 0.751831315650, 0.752399852233, 0.752966215276, 0.753535863766, 0.754101751199, 0.754668828704, 0.755236257084, 0.755801314172, 0.756370839725, 0.756938862351, 0.757505892005, 0.758073244969, 0.758643841350, 0.759212027758, 0.759780773656, 0.760345840893, 0.760913307201, 0.761481879562, 0.762048734045, 0.762614871134, 0.763185206223, 0.763753546131, 0.764321174367, 0.764889141087, 0.765459607878, 0.766028823095, 0.766599158336, 0.767168348189, 0.767739350957, 0.768308616664, 0.768879522288, 0.769448415877, 0.770020326649, 0.770591367997, 0.771162248335, 0.771730377965, 0.772301719839, 0.772870346779, 0.773441595740, 0.774013021506, 0.774582455434, 0.775151863204, 0.775723754640, 0.776294664327, 0.776864366562, 0.777434619396, 0.778007252654, 0.778579251702, 0.779151563719, 0.779723508994, 0.780294296862, 0.780864328190, 0.781435528780, 0.782007024178, 0.782579746452, 0.783152402580, 0.783722320212, 0.784293576083, 0.784864292930, 0.785436518612, 0.786007668447, 0.786581188417, 0.787155235315, 0.787725918324, 0.788297907035, 0.788871256369, 0.789442592942, 0.790016792340, 0.790588485875, 0.791162707757, 0.791734845375, 0.792308246527, 0.792880626606, 0.793454180431, 0.794027225117, 0.794599538430, 0.795175352094, 0.795751789116, 0.796326343045, 0.796901791409, 0.797476921701, 0.798054841514, 0.798630086857, 0.799205963862, 0.799780370608, 0.800353714352, 0.800927667722, 0.801504448405, 0.802079091315, 0.802652116283, 0.803227331097, 0.803800494719, 0.804374883562, 0.804948031978, 0.805525726060, 0.806103158600, 0.806678098604, 0.807253037286, 0.807830543277, 0.808404214099, 0.808979348658, 0.809552234629, 0.810128867286, 0.810704912838, 0.811281366302, 0.811857560883, 0.812435426095, 0.813011539791, 0.813586705402, 0.814162704531, 0.814737002337, 0.815314439294, 0.815891272181, 0.816464512006, 0.817043305804, 0.817619142241, 0.818192382842, 0.818772921874, 0.819351441654, 0.819930342933, 0.820505837102, 0.821081677790, 0.821659334335, 0.822238363040, 0.822818335167, 0.823396656547, 0.823975219145, 0.824552454219, 0.825128445447, 0.825707411161, 0.826282832453, 0.826859582681, 0.827436880942, 0.828016257033, 0.828592994541, 0.829172600205, 0.829751586544, 0.830330047014, 0.830907846074, 0.831486968777, 0.832065047442, 0.832644788743, 0.833222586513, 0.833803052674, 0.834381905133, 0.834961019322, 0.835540389311, 0.836119738152, 0.836697561884, 0.837276024053, 0.837855511870, 0.838435654190, 0.839014447091, 0.839593039815, 0.840173149837, 0.840751041033, 0.841332209708, 0.841912129145, 0.842490207052, 0.843066492526, 0.843647249296, 0.844228528864, 0.844809044159, 0.845387939286, 0.845967713698, 0.846554250776, 0.847134064260, 0.847717377943, 0.848297056882, 0.848878114715, 0.849460282648, 0.850040965842, 0.850622343360, 0.851205332633, 0.851787259375, 0.852369113684, 0.852951829074, 0.853533377666, 0.854115200108, 0.854695496925, 0.855277245500, 0.855856648715, 0.856437521826, 0.857020450891, 0.857602045197, 0.858182010084, 0.858763933242, 0.859343904123, 0.859924257348, 0.860507891280, 0.861087756086, 0.861673541048, 0.862256413806, 0.862840117400, 0.863423744026, 0.864008070281, 0.864590716391, 0.865172301756, 0.865755591559, 0.866339566853, 0.866920965479, 0.867502916208, 0.868086955559, 0.868670985396, 0.869253005546, 0.869837837461, 0.870422529876, 0.871006429324, 0.871591056705, 0.872173550369, 0.872757717340, 0.873342038482, 0.873927202086, 0.874509329567, 0.875092837748, 0.875676123068, 0.876259441603, 0.876841248713, 0.877425529458, 0.878009626626, 0.878593653348, 0.879179584512, 0.879764614818, 0.880348146076, 0.880930643294, 0.881514385670, 0.882097503483, 0.882683788683, 0.883268522954, 0.883854115374, 0.884436853967, 0.885023931606, 0.885609157453, 0.886193919805, 0.886778888326, 0.887364521847, 0.887950738115, 0.888535414993, 0.889122996265, 0.889709621731, 0.890293974150, 0.890877701874, 0.891463240152, 0.892050280249, 0.892637660436, 0.893222433015, 0.893807589340, 0.894393306986, 0.894984020279, 0.895570747258, 0.896156513853, 0.896742633317, 0.897328896872, 0.897914153908, 0.898502067449, 0.899088822600, 0.899675309839, 0.900262490060, 0.900848947849, 0.901436925485, 0.902025852742, 0.902615093823, 0.903201104157, 0.903788727554, 0.904375160856, 0.904961671779, 0.905548787169, 0.906136515408, 0.906726595676, 0.907317969951, 0.907907287400, 0.908495847191, 0.909081866227, 0.909670842901, 0.910259951828, 0.910852510166, 0.911439694079, 0.912029414346, 0.912617219631, 0.913205469505, 0.913794396166, 0.914383320202, 0.914971929865, 0.915563748310, 0.916152503857, 0.916743915834, 0.917332784622, 0.917923056924, 0.918512239011, 0.919097643364, 0.919687602782, 0.920276615237, 0.920866384192, 0.921456444093, 0.922044958578, 0.922634493336, 0.923223719486, 0.923813789880, 0.924403867577, 0.924992256025, 0.925581764707, 0.926171034095, 0.926759676884, 0.927349860741, 0.927940417281, 0.928529276549, 0.929121333162, 0.929711601120, 0.930299295234, 0.930893115805, 0.931483209480, 0.932075992818, 0.932666302299, 0.933256905018, 0.933848233682, 0.934439408840, 0.935031187887, 0.935622557405, 0.936213289635, 0.936802314514, 0.937394496824, 0.937985895091, 0.938578065858, 0.939169954095, 0.939760544381, 0.940350477650, 0.940941630341, 0.941531421324, 0.942125070433, 0.942716277748, 0.943305859262, 0.943898612394, 0.944491384532, 0.945079930428, 0.945669983956, 0.946260767338, 0.946852889611, 0.947444662059, 0.948039550034, 0.948633046930, 0.949223785884, 0.949817275586, 0.950409500806, 0.951002003211, 0.951592956937, 0.952184999866, 0.952776144764, 0.953366382961, 0.953960154789, 0.954553210160, 0.955144095818, 0.955738343734, 0.956329427324, 0.956922099123, 0.957516277897, 0.958112481405, 0.958706226803, 0.959300385652, 0.959894384361, 0.960487364660, 0.961080111426, 0.961673362147, 0.962264446177, 0.962857308599, 0.963450446470, 0.964044691370, 0.964639858579, 0.965232771769, 0.965829638745, 0.966424002729, 0.967018891435, 0.967613914909, 0.968207502435, 0.968800220948, 0.969394073293, 0.969988896848, 0.970587316311, 0.971180388066, 0.971774482448, 0.972370038602, 0.972966718638, 0.973561677807, 0.974156794333, 0.974751776813, 0.975349984172, 0.975945195323, 0.976541935169, 0.977142681057, 0.977739217153, 0.978337631461, 0.978934434576, 0.979528620253, 0.980124453868, 0.980720021774, 0.981316374252, 0.981911430520, 0.982511061765, 0.983107862767, 0.983704803181, 0.984296403581, 0.984892539830, 0.985486731745, 0.986084614627, 0.986682483588, 0.987281088166, 0.987874314795, 0.988469752703, 0.989066435838, 0.989662837283, 0.990258691460, 0.990855934104, 0.991451239826, 0.992049179207, 0.992646401487, 0.993244095547, 0.993840922002, 0.994439410158, 0.995038135986, 0.995632455858, 0.996229635218, 0.996824804325, 0.997423946093, 0.998020413763, 0.998615441951, 0.999212855805, 0.999811261892, 1.000408585298, 1.001005686662, 1.001603736509, 1.002201157491, 1.002797919649, 1.003395857483, 1.003991428828, 1.004587760975, 1.005182174934, 1.005780114665, 1.006377583062, 1.006977220591, 1.007573409249, 1.008173488719, 1.008770494281, 1.009368858606, 1.009968284006, 1.010562032149, 1.011157007554, 1.011754856205, 1.012352675552, 1.012951005741, 1.013550394381, 1.014149808274, 1.014749489000, 1.015348397211, 1.015949952661, 1.016548485047, 1.017146325455, 1.017745740968, 1.018342573613, 1.018942772593, 1.019539809058, 1.020139714499, 1.020741415497, 1.021339936998, 1.021940750810, 1.022538847295, 1.023137919755, 1.023739058199, 1.024340193895, 1.024940456343, 1.025539168289, 1.026141681732, 1.026739752967, 1.027339028208, 1.027938890682, 1.028539786922, 1.029139105306, 1.029739084455, 1.030338272924, 1.030935812853, 1.031533260739, 1.032132953294, 1.032734669224, 1.033334645491, 1.033931530703, 1.034529975658, 1.035129769071, 1.035730335406, 1.036329844634, 1.036930527756, 1.037532540018, 1.038132533540, 1.038733689486, 1.039334570879, 1.039936784722, 1.040540521357, 1.041143317557, 1.041745239951, 1.042346135320, 1.042947153679, 1.043545721055, 1.044146003956, 1.044744633555, 1.045343467546, 1.045945015969, 1.046546286879, 1.047146329139, 1.047744995918, 1.048343722104, 1.048945914388, 1.049548134775, 1.050154491489, 1.050757507739, 1.051357647451, 1.051957702283, 1.052559753941, 1.053158351378, 1.053761328724, 1.054364120619, 1.054968627689, 1.055571864842, 1.056175787838, 1.056778081996, 1.057379007007, 1.057981330469, 1.058586120663, 1.059187603373, 1.059789681041, 1.060393073673, 1.060997735617, 1.061599281811, 1.062202850147, 1.062806330010, 1.063409267553, 1.064012791695, 1.064616243523, 1.065217466650, 1.065820650206, 1.066423933819, 1.067027027881, 1.067631666002, 1.068230855813, 1.068834298234, 1.069437372559, 1.070043586955, 1.070647076734, 1.071251349992, 1.071855860525, 1.072459499750, 1.073061584560, 1.073665760794, 1.074270320099, 1.074874592087, 1.075475882236, 1.076079113367, 1.076682644656, 1.077285127091, 1.077886534351, 1.078489863066, 1.079094859570, 1.079702067101, 1.080305558076, 1.080909689993, 1.081514275738, 1.082116545674, 1.082721275526, 1.083324154714, 1.083927313478, 1.084528134869, 1.085132583544, 1.085739102131, 1.086342818477, 1.086948309003, 1.087551595234, 1.088155619580, 1.088763548657, 1.089365404782, 1.089969079172, 1.090574412402, 1.091180017287, 1.091785646995, 1.092391966639, 1.092993291294, 1.093598407396, 1.094205522325, 1.094808576270, 1.095415422636, 1.096020940273, 1.096625643243, 1.097234085001, 1.097839000907, 1.098446019024, 1.099050487973, 1.099654073080, 1.100261539939, 1.100867195183, 1.101476455718, 1.102082880368, 1.102691676874, 1.103299273081, 1.103898551605, 1.104505707366, 1.105110377232, 1.105718964800, 1.106324534117, 1.106930771232, 1.107536259158, 1.108145745578, 1.108754018822, 1.109361530686, 1.109966352310, 1.110582168667, 1.111187796177, 1.111794364925, 1.112400516230, 1.113011390717, 1.113618640563, 1.114226480805, 1.114833554951, 1.115441614846, 1.116046754565, 1.116654783877, 1.117260798320, 1.117866844897, 1.118472443338, 1.119079715814, 1.119688244818, 1.120294705881, 1.120899771840, 1.121503918251, 1.122115786884, 1.122724641514, 1.123333814380, 1.123938146494, 1.124544992685, 1.125152803887, 1.125759204329, 1.126369037971, 1.126977378728, 1.127586753425, 1.128194528168, 1.128802418520, 1.129409084006, 1.130018736558, 1.130625702496, 1.131237842031, 1.131842861773, 1.132450693347, 1.133058969733, 1.133666332733, 1.134274694235, 1.134881088967, 1.135488580762, 1.136094844059, 1.136702859203, 1.137310404106, 1.137919903871, 1.138528605242, 1.139133859169, 1.139741725234, 1.140350875283, 1.140960844890, 1.141571401658, 1.142178672715, 1.142787035445, 1.143395599035, 1.144006982991, 1.144616508345, 1.145226908564, 1.145836363319, 1.146447453365, 1.147054219217, 1.147666635999, 1.148272426403, 1.148882783216, 1.149491799289, 1.150100553840, 1.150708657458, 1.151316469247, 1.151924017595, 1.152533348792, 1.153143913031, 1.153748511714, 1.154356648691, 1.154965005508, 1.155572339133, 1.156185706638, 1.156794003182, 1.157403770699, 1.158013589461, 1.158625442142, 1.159234774166, 1.159849147904, 1.160461444724, 1.161074014477, 1.161686567210, 1.162296987523, 1.162906080425, 1.163516250270, 1.164125237951, 1.164735943876, 1.165346391267, 1.165958283419, 1.166567157661, 1.167177575998, 1.167785900849, 1.168397811792, 1.169006433064, 1.169619810623, 1.170230475806, 1.170839091728, 1.171450424408, 1.172061334443, 1.172672284195, 1.173287893962, 1.173897831127, 1.174510774546, 1.175124330771, 1.175737218881, 1.176348893801, 1.176958017579, 1.177571820821, 1.178184692675, 1.178795880684, 1.179401949948, 1.180014749455, 1.180627861901, 1.181236653148, 1.181849962740, 1.182461753209, 1.183072189025, 1.183685327015, 1.184296961847, 1.184908634927, 1.185521570166, 1.186132551623, 1.186743345682, 1.187355828988, 1.187967202262, 1.188579128998, 1.189195168715, 1.189807410697, 1.190420348676, 1.191033768622, 1.191647732124, 1.192267575162, 1.192881532875, 1.193496224146, 1.194110123023, 1.194723333711, 1.195333802384, 1.195949508282, 1.196563821128, 1.197178621197, 1.197794648938, 1.198409233739, 1.199025410598, 1.199640929265, 1.200252982631, 1.200869562022, 1.201482770425, 1.202097751933, 1.202712850558, 1.203324798695, 1.203937367219, 1.204548018752, 1.205161682420, 1.205778599976, 1.206394208742, 1.207007445644, 1.207612583420, 1.208224864716, 1.208843026045, 1.209460028438, 1.210074083401, 1.210686665689, 1.211303377092, 1.211922387498, 1.212539426295, 1.213153398315, 1.213768637419, 1.214382586407, 1.214998843787, 1.215613429645, 1.216231237229, 1.216846762252, 1.217458460344, 1.218068961842, 1.218684792231, 1.219298438858, 1.219914668821, 1.220526838471, 1.221141100733, 1.221755075162, 1.222369455009, 1.222981671730, 1.223599164033, 1.224215454101, 1.224830047257, 1.225448489162, 1.226063588468, 1.226678733118, 1.227297015154, 1.227912493959, 1.228529007964, 1.229143954591, 1.229761925621, 1.230378931777, 1.230996283519, 1.231612907664, 1.232227420983, 1.232840629881, 1.233457761333, 1.234073187596, 1.234684380247, 1.235300219471, 1.235916050977, 1.236533924973, 1.237151839518, 1.237768321863, 1.238384891038, 1.239002163617, 1.239617404134, 1.240230957605, 1.240845250553, 1.241461072443, 1.242075251127, 1.242691719705, 1.243309026562, 1.243923244332, 1.244536799099, 1.245152940278, 1.245773153753, 1.246397088099, 1.247012762005, 1.247628910547, 1.248245749879, 1.248858159383, 1.249474341708, 1.250093029414, 1.250708507609, 1.251323728295, 1.251939713132, 1.252553069288, 1.253173609913, 1.253792607758, 1.254410390428, 1.255027974963, 1.255646298146, 1.256261773380, 1.256877219776, 1.257492431680, 1.258105612843, 1.258721378737, 1.259339155218, 1.259954398035, 1.260569571997, 1.261185404623, 1.261802072065, 1.262421213730, 1.263037005473, 1.263654117894, 1.264274399000, 1.264889660911, 1.265503978861, 1.266120585665, 1.266737210338, 1.267354141074, 1.267975273787, 1.268595224487, 1.269214035268, 1.269832233520, 1.270449353933, 1.271064344877, 1.271683812461, 1.272302652732, 1.272918606548, 1.273534293780, 1.274152414537, 1.274771849646, 1.275389148167, 1.276007546743, 1.276628600764, 1.277244952519, 1.277862196715, 1.278477713566, 1.279091791763, 1.279708079173, 1.280330260996, 1.280951850919, 1.281569473166, 1.282186452973, 1.282803952396, 1.283418135666, 1.284038617921, 1.284656918144, 1.285274583754, 1.285892038590, 1.286509801318, 1.287125516577, 1.287742729387, 1.288357986076, 1.288980198290, 1.289596298123, 1.290215789469, 1.290834350084, 1.291453096256, 1.292075082106, 1.292691526842, 1.293309103798, 1.293930337347, 1.294550364275, 1.295168252635, 1.295784549832, 1.296399144707, 1.297014472869, 1.297633862826, 1.298250590277, 1.298874616387, 1.299495948484, 1.300114400356, 1.300733517183, 1.301352700784, 1.301969538973, 1.302588788728, 1.303207613282, 1.303826752624, 1.304446084281, 1.305064485642, 1.305685234987, 1.306306459629, 1.306930431947, 1.307550073371, 1.308172383789, 1.308792280090, 1.309412805673, 1.310036835165, 1.310651292411, 1.311271235768, 1.311896268933, 1.312516735746, 1.313132785021, 1.313750550148, 1.314372258512, 1.314995871730, 1.315614722661, 1.316235895947, 1.316856238393, 1.317473111006, 1.318089379515, 1.318713688815, 1.319335100788, 1.319950888958, 1.320573139154, 1.321197146560, 1.321816575663, 1.322432746894, 1.323050488044, 1.323671113126, 1.324288649293, 1.324908083430, 1.325530221940, 1.326149718955, 1.326769308313, 1.327390936554, 1.328016163953, 1.328640691454, 1.329258174555, 1.329882431220, 1.330500270666, 1.331119111338, 1.331741902123, 1.332365017809, 1.332983307153, 1.333603114633, 1.334225992506, 1.334846582045, 1.335466696262, 1.336087603024, 1.336710332295, 1.337332888693, 1.337952357538, 1.338574567397, 1.339197148362, 1.339819132016, 1.340448066286, 1.341065197076, 1.341684398334, 1.342306269892, 1.342929482777, 1.343552192670, 1.344172554248, 1.344793044306, 1.345414162391, 1.346032720989, 1.346649749859, 1.347272661398, 1.347890788089, 1.348512285974, 1.349137032435, 1.349756382125, 1.350377054232, 1.350999004444, 1.351621534336, 1.352235898355, 1.352857826842, 1.353479725954, 1.354100818910, 1.354723823586, 1.355345606998, 1.355964329190, 1.356593568468, 1.357218422049, 1.357841681198, 1.358462356266, 1.359086431474, 1.359708014381, 1.360331314607, 1.360949258893, 1.361573146161, 1.362192821205, 1.362816576103, 1.363443003163, 1.364063003908, 1.364685389530, 1.365311871110, 1.365938904699, 1.366562077406, 1.367186105133, 1.367814880094, 1.368439484376, 1.369060670724, 1.369684303385, 1.370303759207, 1.370927405488, 1.371550814493, 1.372172909558, 1.372793909169, 1.373417316689, 1.374039452125, 1.374662562387, 1.375289103122, 1.375916383922, 1.376536810211, 1.377161291227, 1.377787054997, 1.378412829035, 1.379039339731, 1.379661446462, 1.380281453696, 1.380909300204, 1.381536518808, 1.382160069178, 1.382781192379, 1.383403562175, 1.384027214201, 1.384656533287, 1.385279172857, 1.385905114190, 1.386533270393, 1.387158608369, 1.387782896485, 1.388403877135, 1.389027065936, 1.389653483080, 1.390278810232, 1.390898330259, 1.391520928390, 1.392148009365, 1.392768690692, 1.393392387707, 1.394011870871, 1.394632594511, 1.395261393083, 1.395885440946, 1.396512561972, 1.397137035219, 1.397761778222, 1.398389997082, 1.399017166811, 1.399643009044, 1.400274000704, 1.400896093861, 1.401520546324, 1.402146172029, 1.402769615235, 1.403394658354, 1.404021934607, 1.404644042668, 1.405270676013, 1.405896169067, 1.406517957502, 1.407148888355, 1.407771830068, 1.408396233816, 1.409025045062, 1.409651722988, 1.410276938413, 1.410904028416, 1.411529459827, 1.412154895737, 1.412777155433, 1.413401973112, 1.414028209296, 1.414655033858, 1.415282753135, 1.415908382268, 1.416531457613, 1.417156619756, 1.417782899044, 1.418410242110, 1.419039336131, 1.419667756144, 1.420293771836, 1.420918207136, 1.421542544204, 1.422168468921, 1.422795205077, 1.423422225262, 1.424047880597, 1.424671090030, 1.425297831476, 1.425925698717, 1.426553419681, 1.427185846915, 1.427807135640, 1.428433938314, 1.429061833563, 1.429691525672, 1.430321406844, 1.430944998325, 1.431575281894, 1.432206411019, 1.432832305166, 1.433457076025, 1.434084457911, 1.434707418114, 1.435333604406, 1.435956689339, 1.436581073002, 1.437203883606, 1.437828052820, 1.438454383791, 1.439077370085, 1.439709367647, 1.440334038383, 1.440959620901, 1.441586994232, 1.442213423684, 1.442835075183, 1.443461427845, 1.444091570997, 1.444720429007, 1.445344567934, 1.445968477312, 1.446595093817, 1.447224330587, 1.447848987325, 1.448475739141, 1.449103579983, 1.449728121809, 1.450356209076, 1.450977931222, 1.451603407160, 1.452225823640, 1.452853026768, 1.453482531268, 1.454111664472, 1.454736822501, 1.455362224944, 1.455990452852, 1.456618658733, 1.457242522421, 1.457868816633, 1.458485318477, 1.459116308712, 1.459737451766, 1.460365890192, 1.460994912905, 1.461626381806, 1.462260155165, 1.462884994659, 1.463512666346, 1.464141739663, 1.464762265930, 1.465389982796, 1.466019255959, 1.466648272138, 1.467276710473, 1.467904873256, 1.468529564556, 1.469154503317, 1.469783609324, 1.470405058168, 1.471029465962, 1.471658066467, 1.472287114284, 1.472914003580, 1.473547110008, 1.474173531515, 1.474799043557, 1.475426054811, 1.476053322778, 1.476679103169, 1.477304691671, 1.477935477124, 1.478564683287, 1.479193087424, 1.479818521573, 1.480446013065, 1.481073794547, 1.481701220784, 1.482334446820, 1.482964107656, 1.483594748853, 1.484224638431, 1.484851357946, 1.485483339688, 1.486110948460, 1.486733897542, 1.487360436065, 1.487983952262, 1.488619789929, 1.489243617571, 1.489871991994, 1.490507659296, 1.491136413371, 1.491768235162, 1.492402974804, 1.493030993958, 1.493661086560, 1.494288285919, 1.494913243272, 1.495543546263, 1.496167849926, 1.496794088401, 1.497417695267, 1.498046478165, 1.498674242213, 1.499296525449, 1.499929340920, 1.500559284445, 1.501188380308, 1.501815368502, 1.502446025408, 1.503071347762, 1.503697294770, 1.504323243639, 1.504948151252, 1.505574821824, 1.506200726231, 1.506829068629, 1.507457916264, 1.508085870844, 1.508709368160, 1.509337929401, 1.509967879470, 1.510603796812, 1.511232824732, 1.511858502674, 1.512491952563, 1.513121529208, 1.513751495355, 1.514381155771, 1.515011630907, 1.515645998253, 1.516275134055, 1.516907767026, 1.517536761326, 1.518161269222, 1.518786676443, 1.519413947800, 1.520043334573, 1.520676889915, 1.521303250025, 1.521938055115, 1.522569854631, 1.523201560503, 1.523839873882, 1.524465483429, 1.525098397722, 1.525726319576, 1.526358448614, 1.526991046035, 1.527623951555, 1.528255846061, 1.528890319947, 1.529518313500, 1.530150410878, 1.530783650724, 1.531411186282, 1.532041877252, 1.532669013860, 1.533293706280, 1.533926158950, 1.534553288161, 1.535184377227, 1.535814654640, 1.536452255614, 1.537077930309, 1.537710050037, 1.538341576039, 1.538971542910, 1.539604140020, 1.540227143399, 1.540858736785, 1.541485190811, 1.542114713664, 1.542740482923, 1.543370933936, 1.544002225495, 1.544630783218, 1.545262065788, 1.545887154182, 1.546516109093, 1.547148548014, 1.547782307866, 1.548413708316, 1.549046719966, 1.549678900247, 1.550307884868, 1.550937287598, 1.551561433297, 1.552189094609, 1.552821495258, 1.553457179058, 1.554082439470, 1.554709349089, 1.555337367767, 1.555969561033, 1.556600485558, 1.557231795282, 1.557865640209, 1.558495744563, 1.559128841681, 1.559759695368, 1.560388231362, 1.561015149492, 1.561647611549, 1.562275590737, 1.562909527272, 1.563539653086, 1.564172891515, 1.564797743538, 1.565434017883, 1.566063854339, 1.566700226223, 1.567328711717, 1.567955281506, 1.568588353493, 1.569219450460, 1.569850869032, 1.570484353860, 1.571118861135, 1.571753081722, 1.572381448038, 1.573015128016, 1.573647325737, 1.574277414151, 1.574903375514, 1.575535077759, 1.576163068530, 1.576797787145, 1.577431694861, 1.578066890836, 1.578697255480, 1.579330844442, 1.579955337579, 1.580587145047, 1.581218995436, 1.581849096748, 1.582479731641, 1.583112231523, 1.583743205717, 1.584370694844, 1.585000611845, 1.585626726865, 1.586258268613, 1.586894957516, 1.587528834746, 1.588153123182, 1.588786768367, 1.589424680374, 1.590060590428, 1.590693388352, 1.591326177529, 1.591956716223, 1.592591180189, 1.593225976634, 1.593859315740, 1.594494996761, 1.595126070257, 1.595756538288, 1.596384596520, 1.597019814661, 1.597649978918, 1.598283280426, 1.598921197557, 1.599551261682, 1.600177000309, 1.600808180094, 1.601443002026, 1.602074688102, 1.602702492575, 1.603334951796, 1.603970496998, 1.604599949143, 1.605232117327, 1.605863910227, 1.606495640764, 1.607126481270, 1.607759559538, 1.608392380996, 1.609024007585, 1.609658127506, 1.610292289457, 1.610923620004, 1.611551040098, 1.612177732053, 1.612811277202, 1.613443590096, 1.614079932456, 1.614710805607, 1.615349259573, 1.615979863999, 1.616609068097, 1.617234867604, 1.617866326275, 1.618496737850, 1.619130234023, 1.619768527875, 1.620399570423, 1.621037066240, 1.621673281935, 1.622305589582, 1.622940496094, 1.623574524998, 1.624206464565, 1.624841650035, 1.625480002647, 1.626111950468, 1.626747227991, 1.627370932656, 1.628010121254, 1.628643806484, 1.629282875089, 1.629912160893, 1.630544826905, 1.631172768705, 1.631806679833, 1.632442821734, 1.633078628064, 1.633712657221, 1.634347257816, 1.634983030680, 1.635614405722, 1.636249970247, 1.636882570066, 1.637522199972, 1.638146217964, 1.638777693634, 1.639408479615, 1.640035557037, 1.640667756258, 1.641303082797, 1.641937873844, 1.642573822942, 1.643210628251, 1.643847603622, 1.644482890459, 1.645113967098, 1.645750860705, 1.646384341791, 1.647016128014, 1.647647850644, 1.648284203461, 1.648921702930, 1.649554093169, 1.650189288137, 1.650817114348, 1.651453343145, 1.652085183857, 1.652722356924, 1.653354679454, 1.653990547457, 1.654622740200, 1.655247372135, 1.655875479657, 1.656510109152, 1.657141979560, 1.657779234080, 1.658411629864, 1.659043482100, 1.659681372748, 1.660317082798, 1.660952928962, 1.661587595500, 1.662221514834, 1.662852004752, 1.663485292353, 1.664126098063, 1.664756149905, 1.665392041370, 1.666026631002, 1.666658702255, 1.667289776985, 1.667923670422, 1.668559725456, 1.669193124577, 1.669826859976, 1.670462111209, 1.671095727027, 1.671730635714, 1.672364737877, 1.673003612651, 1.673633063545, 1.674265069180, 1.674902433379, 1.675536784092, 1.676166932299, 1.676795850237, 1.677427684647, 1.678062922911, 1.678704376464, 1.679334886167, 1.679968224802, 1.680602009639, 1.681240994171, 1.681876932774, 1.682515581075, 1.683148302241, 1.683782135276, 1.684411077995, 1.685037442012, 1.685668227671, 1.686304887106, 1.686943959743, 1.687573905549, 1.688214044709, 1.688851373803, 1.689486537673, 1.690125993595, 1.690766115432, 1.691397835404, 1.692031073983, 1.692661042500, 1.693293640767, 1.693923105143, 1.694557416792, 1.695189664047, 1.695824277455, 1.696462971640, 1.697100898403, 1.697735108149, 1.698366038577, 1.699006139140, 1.699647141034, 1.700282055286, 1.700914540019, 1.701544321125, 1.702177773148, 1.702804504746, 1.703441137165, 1.704079846722, 1.704709331222, 1.705344423047, 1.705977068555, 1.706607786088, 1.707241877608, 1.707874635497, 1.708512955826, 1.709148570052, 1.709789523849, 1.710423220358, 1.711058624338, 1.711692879811, 1.712328421210, 1.712965073424, 1.713599965255, 1.714235426724, 1.714873013168, 1.715508830105, 1.716146031072, 1.716776136892, 1.717417127186, 1.718047561026, 1.718682319795, 1.719311635950, 1.719946149704, 1.720590424305, 1.721230010454, 1.721870517021, 1.722507338930, 1.723139723925, 1.723782342534, 1.724423979283, 1.725057987334, 1.725696386239, 1.726335586190, 1.726978924559, 1.727611016579, 1.728242287752, 1.728877013482, 1.729513110869, 1.730146081709, 1.730783948485, 1.731416810063, 1.732056617624, 1.732695984520, 1.733330559419, 1.733969711136, 1.734602756947, 1.735240905128, 1.735869945105, 1.736505721528, 1.737135625406, 1.737774327252, 1.738415491673, 1.739052483974, 1.739687240005, 1.740325576274, 1.740954853997, 1.741594243094, 1.742230040951, 1.742868260478, 1.743501260023, 1.744138484505, 1.744782050536, 1.745410767162, 1.746041314904, 1.746673700425, 1.747307420730, 1.747942042800, 1.748572311333, 1.749206859968, 1.749842752126, 1.750478109774, 1.751117556856, 1.751749732301, 1.752386414577, 1.753022383892, 1.753660295653, 1.754294015470, 1.754929303582, 1.755571434580, 1.756206316246, 1.756847585615, 1.757482572849, 1.758114682783, 1.758747938357, 1.759388008185, 1.760026523487, 1.760655341578, 1.761291713761, 1.761934592534, 1.762573396419, 1.763213040598, 1.763849820993, 1.764487283927, 1.765120929845, 1.765764137336, 1.766400664802, 1.767035205655, 1.767673549300, 1.768307355715, 1.768938721025, 1.769576293980, 1.770209124300, 1.770841417592, 1.771479174655, 1.772117252881, 1.772756115553, 1.773395455470, 1.774037751302, 1.774665566114, 1.775302987863, 1.775934891164, 1.776569870151, 1.777206870940, 1.777847801949, 1.778484803687, 1.779121043581, 1.779744930769, 1.780374770408, 1.781014390463, 1.781658736642, 1.782292070245, 1.782920295654, 1.783567926418, 1.784201092965, 1.784833807955, 1.785473912215, 1.786111404982, 1.786740824795, 1.787382656745, 1.788020693741, 1.788660657227, 1.789288277871, 1.789921251215, 1.790563996708, 1.791198994591, 1.791829732208, 1.792472051835, 1.793107123702, 1.793749581080, 1.794398129910, 1.795036785164, 1.795663087151, 1.796294015818, 1.796926841931, 1.797564896810, 1.798199087639, 1.798833522702, 1.799479860893, 1.800113375340, 1.800750615102, 1.801381616456, 1.802013591006, 1.802646927644, 1.803284668100, 1.803923069916, 1.804560251616, 1.805207521228, 1.805836565484, 1.806479900539, 1.807115805624, 1.807754819052, 1.808391868665, 1.809023671056, 1.809655805536, 1.810296044637, 1.810932030323, 1.811569905546, 1.812207760800, 1.812843054400, 1.813483546859, 1.814125268384, 1.814765813292, 1.815400206470, 1.816037318961, 1.816676079342, 1.817306284620, 1.817945144898, 1.818581828932, 1.819223142564, 1.819864544073, 1.820503188481, 1.821137364854, 1.821774427779, 1.822414273200, 1.823054744838, 1.823689130330, 1.824320850251, 1.824959178889, 1.825601615315, 1.826242412777, 1.826882087340, 1.827518121897, 1.828151960913, 1.828796974469, 1.829429370497, 1.830059898954, 1.830706166404, 1.831348743351, 1.831983246495, 1.832612829224, 1.833251520257, 1.833888278229, 1.834523508618, 1.835168137741, 1.835807624851, 1.836448233784, 1.837084297571, 1.837720397640, 1.838361651246, 1.839000495679, 1.839646645723, 1.840286933585, 1.840924764399, 1.841564950596, 1.842207622019, 1.842845983501, 1.843496491358, 1.844127590949, 1.844766656469, 1.845405020964, 1.846042253420, 1.846686280513, 1.847322005293, 1.847971453642, 1.848612772191, 1.849239539481, 1.849875941816, 1.850510138327, 1.851155034534, 1.851796289671, 1.852442728947, 1.853084403340, 1.853715832503, 1.854356659167, 1.854993705387, 1.855624181076, 1.856261375078, 1.856902160582, 1.857535946767, 1.858171285811, 1.858816403674, 1.859450792070, 1.860095580493, 1.860735465781, 1.861382071094, 1.862019968233, 1.862656207994, 1.863295314998, 1.863934252710, 1.864574290792, 1.865206962850, 1.865842152312, 1.866478367952, 1.867115932739, 1.867751584336, 1.868388360116, 1.869023051415, 1.869661566886, 1.870301505762, 1.870946325292, 1.871588484609, 1.872225122784, 1.872861528913, 1.873504516215, 1.874143061312, 1.874779226438, 1.875417335364, 1.876053771506, 1.876693037970, 1.877336390132, 1.877968727918, 1.878611149838, 1.879242486578, 1.879878200632, 1.880523257521, 1.881161907546, 1.881797826082, 1.882434147055, 1.883073326076, 1.883699392196, 1.884338540751, 1.884974065898, 1.885604348450, 1.886233976262, 1.886876098806, 1.887511998949, 1.888148093044, 1.888794264308, 1.889431028599, 1.890066131662, 1.890703279144, 1.891340652608, 1.891978217811, 1.892613057484, 1.893265673533, 1.893907127635, 1.894545851349, 1.895188176955, 1.895824688822, 1.896459294913, 1.897098667618, 1.897751612248, 1.898390004701, 1.899036806195, 1.899667509018, 1.900304721549, 1.900937995740, 1.901573129504, 1.902209262827, 1.902848968330, 1.903485514000, 1.904115889519, 1.904751017773, 1.905389416392, 1.906030328896, 1.906665848006, 1.907299702384, 1.907941931548, 1.908586061860, 1.909220894707, 1.909860362259, 1.910504059542, 1.911149490983, 1.911793047312, 1.912426979993, 1.913066780920, 1.913707739411, 1.914348040745, 1.914995322306, 1.915630588435, 1.916270760520, 1.916919409799, 1.917566407076, 1.918205807743, 1.918840422884, 1.919479466441, 1.920114429251, 1.920755062099, 1.921393778445, 1.922021350276, 1.922657208202, 1.923301569328, 1.923942221975, 1.924593677901, 1.925234485928, 1.925879279939, 1.926518432170, 1.927150080130, 1.927776874286, 1.928422399926, 1.929060218534, 1.929699898712, 1.930333012871, 1.930967940457, 1.931603649245, 1.932240178464, 1.932880618816, 1.933523048434, 1.934157771437, 1.934803552111, 1.935434274179, 1.936073335540, 1.936714352250, 1.937351015506, 1.937986617980, 1.938618438849, 1.939248423749, 1.939889987779, 1.940526820222, 1.941166446291, 1.941810548236, 1.942451993088, 1.943090996219, 1.943729681929, 1.944367627146, 1.945010145934, 1.945653501744, 1.946291402948, 1.946938928643, 1.947579876842, 1.948225010451, 1.948860371465, 1.949497629900, 1.950141400306, 1.950788220358, 1.951421442413, 1.952061539322, 1.952700867333, 1.953343868282, 1.953988603919, 1.954632733260, 1.955275782023, 1.955911622969, 1.956547610310, 1.957183546563, 1.957820927710, 1.958459245668, 1.959092692912, 1.959727302931, 1.960368906777, 1.961000383605, 1.961649796421, 1.962296956494, 1.962936548636, 1.963572970746, 1.964207927216, 1.964846056316, 1.965481835294, 1.966118546376, 1.966766571535, 1.967401544301, 1.968044427063, 1.968690202502, 1.969325446597, 1.969966646622, 1.970606603001, 1.971250227162, 1.971894603049, 1.972534106587, 1.973168224885, 1.973811570606, 1.974450547467, 1.975080131799, 1.975719295409, 1.976367750263, 1.977002099172, 1.977637252245, 1.978273913930, 1.978906587192, 1.979547143915, 1.980182962336, 1.980813729363, 1.981453902994, 1.982103315144, 1.982744183533, 1.983380272535, 1.984026252733, 1.984682376553, 1.985319381409, 1.985955850305, 1.986585463188, 1.987223538608, 1.987867916444, 1.988508853043, 1.989147517748, 1.989791534754, 1.990429328533, 1.991072229804, 1.991709352041, 1.992354792454, 1.993000381569, 1.993628186342, 1.994266886210, 1.994907986159, 1.995542640238, 1.996183905620, 1.996829094098, 1.997460864557, 1.998093468916, 1.998732928547, 1.999371553078, 2.000002996642, 2.000645493627, 2.001287025961, 2.001929507358, 2.002565208010, 2.003209715755, 2.003849660525, 2.004489277121, 2.005132254124, 2.005771695380, 2.006412696574, 2.007061310552, 2.007695157118, 2.008335154075, 2.008966607669, 2.009612522939, 2.010254508884, 2.010887202772, 2.011520329166, 2.012162910433, 2.012805280915, 2.013446721542, 2.014080991112, 2.014717760784, 2.015359154867, 2.016000461238, 2.016645152814, 2.017281583473, 2.017920305915, 2.018557974793, 2.019193131539, 2.019835264563, 2.020481171072, 2.021120790176, 2.021758339069, 2.022398791570, 2.023042662913, 2.023688636859, 2.024330153064, 2.024963372403, 2.025612395968, 2.026255147502, 2.026897696604, 2.027538652786, 2.028171195557, 2.028807723871, 2.029455922596, 2.030107789946, 2.030748983096, 2.031388416595, 2.032025146204, 2.032670630001, 2.033298740101, 2.033933911514, 2.034569919194, 2.035205682232, 2.035843792301, 2.036490068396, 2.037135273473, 2.037767983715, 2.038411580538, 2.039044538782, 2.039686748186, 2.040327239951, 2.040973784389, 2.041615509322, 2.042265699182, 2.042912261663, 2.043555563179, 2.044205733151, 2.044847487815, 2.045480981274, 2.046129260951, 2.046773576213, 2.047414343908, 2.048054894141, 2.048696341948, 2.049335089496, 2.049975460040, 2.050616190635, 2.051257281482, 2.051898928546, 2.052537064704, 2.053175649032, 2.053814485494, 2.054455395368, 2.055101689838, 2.055740849612, 2.056377193262, 2.057013628800, 2.057649510530, 2.058279470630, 2.058920442918, 2.059564554408, 2.060204434318, 2.060844059345, 2.061488230291, 2.062134360313, 2.062766598834, 2.063404734427, 2.064042953440, 2.064683724667, 2.065324382649, 2.065959162453, 2.066600086007, 2.067242818759, 2.067881781855, 2.068527178447, 2.069167424243, 2.069797649861, 2.070431447295, 2.071066017576, 2.071709302757, 2.072354722155, 2.072995193618, 2.073646593237, 2.074289283312, 2.074922810826, 2.075565584774, 2.076199735984, 2.076836006723, 2.077473782029, 2.078116810269, 2.078757043348, 2.079388836300, 2.080033716258, 2.080670508199, 2.081309963472, 2.081945221726, 2.082584562065, 2.083230631338, 2.083887146012, 2.084519962277, 2.085158721516, 2.085797574976, 2.086429632619, 2.087089094083, 2.087735579230, 2.088378450291, 2.089017689769, 2.089666253501, 2.090307820458, 2.090943214896, 2.091582543296, 2.092227808750, 2.092862685940, 2.093506141067, 2.094151360104, 2.094781494499, 2.095422282610, 2.096068839730, 2.096714624450, 2.097353817218, 2.097994224247, 2.098638738191, 2.099283445865, 2.099921841575, 2.100565885190, 2.101198988002, 2.101847071165, 2.102482210745, 2.103125714853, 2.103763941838, 2.104409017941, 2.105056602440, 2.105694961002, 2.106334647596, 2.106968166218, 2.107616242330, 2.108258544315, 2.108895826370, 2.109541254553, 2.110174378173, 2.110804614588, 2.111443794041, 2.112087176205, 2.112731963326, 2.113370717512, 2.114010412552, 2.114642003260, 2.115295128214, 2.115938630421, 2.116582349100, 2.117222301750, 2.117864850929, 2.118507096979, 2.119153151178, 2.119784774671, 2.120430498741, 2.121071731992, 2.121702705443, 2.122350944228, 2.123001823833, 2.123647618263, 2.124288592448, 2.124922465347, 2.125558250569, 2.126196361684, 2.126832910913, 2.127477966355, 2.128117273169, 2.128753374856, 2.129398366013, 2.130040390942, 2.130689996941, 2.131334346926, 2.131978948107, 2.132632346526, 2.133269847283, 2.133909112805, 2.134546419819, 2.135185374887, 2.135829546481, 2.136468906935, 2.137116475388, 2.137750697112, 2.138390385645, 2.139025215704, 2.139673794419, 2.140310864607, 2.140962929744, 2.141614049862, 2.142260784274, 2.142899973643, 2.143531764511, 2.144173857595, 2.144811991378, 2.145453128214, 2.146092416306, 2.146724305534, 2.147364066527, 2.148001412809, 2.148648502046, 2.149279714427, 2.149916936251, 2.150556261479, 2.151196160156, 2.151846244027, 2.152489775177, 2.153126352263, 2.153767700353, 2.154405349099, 2.155039963471, 2.155669415613, 2.156314470766, 2.156941784791, 2.157583115567, 2.158222894069, 2.158853285821, 2.159500770888, 2.160149097132, 2.160796946191, 2.161434804035, 2.162079970724, 2.162725276059, 2.163363697010, 2.164000840382, 2.164645899286, 2.165287914454, 2.165935016593, 2.166584486495, 2.167221142225, 2.167874649236, 2.168518193870, 2.169166283886, 2.169796978935, 2.170433220876, 2.171068722004, 2.171723598828, 2.172367386518, 2.173018857143, 2.173661006239, 2.174296061595, 2.174948549758, 2.175575862017, 2.176219199339, 2.176855464156, 2.177501812261, 2.178141727549, 2.178793402858, 2.179439426967, 2.180087005202, 2.180735879709, 2.181382823651, 2.182026770292, 2.182670813355, 2.183312302387, 2.183950163301, 2.184601651017, 2.185239613506, 2.185881979235, 2.186522560647, 2.187159477308, 2.187787491703, 2.188421576084, 2.189057527261, 2.189699721217, 2.190337884739, 2.190979818845, 2.191623513398, 2.192262077986, 2.192906390673, 2.193540674523, 2.194191982045, 2.194838894347, 2.195482480204, 2.196128317395, 2.196775936300, 2.197416653680, 2.198070172523, 2.198708685941, 2.199362228932, 2.199998653998, 2.200642699511, 2.201290048754, 2.201928754155, 2.202571100609, 2.203209336574, 2.203840525915, 2.204482717874, 2.205115135154, 2.205767376863, 2.206407466844, 2.207054937738, 2.207698400927, 2.208342678572, 2.208981589139, 2.209620737293, 2.210266677063, 2.210914496565, 2.211554730637, 2.212192795123, 2.212838179035, 2.213478417344, 2.214120027522, 2.214755038696, 2.215390551838, 2.216030138733, 2.216668094064, 2.217310282862, 2.217951701003, 2.218585734278, 2.219224939360, 2.219851540957, 2.220496076716, 2.221138968933, 2.221796420284, 2.222435008393, 2.223083973279, 2.223734272850, 2.224379723143, 2.225027009104, 2.225675407285, 2.226311830747, 2.226968376854, 2.227626210407, 2.228266749972, 2.228897936180, 2.229544923470, 2.230189924591, 2.230834037681, 2.231467414795, 2.232114242014, 2.232747560533, 2.233387899081, 2.234026279962, 2.234666271552, 2.235309746115, 2.235958287976, 2.236601434629, 2.237245910104, 2.237880227804, 2.238516526325, 2.239148862389, 2.239785590455, 2.240426199689, 2.241062686230, 2.241705562265, 2.242342941499, 2.242986349006, 2.243630330608, 2.244271152768, 2.244916967416, 2.245560838955, 2.246199006010, 2.246838802215, 2.247475933574, 2.248120537136, 2.248761863315, 2.249408456901, 2.250053620172, 2.250694792644, 2.251335208826, 2.251987044463, 2.252621911456, 2.253257863501, 2.253900047117, 2.254532568165, 2.255176406367, 2.255817834796, 2.256475184171, 2.257126699817, 2.257766221299, 2.258406607160, 2.259047859837, 2.259691245333, 2.260335822656, 2.260973041384, 2.261615955811, 2.262243220568, 2.262894940893, 2.263528837365, 2.264174033238, 2.264812837338, 2.265457624173, 2.266093672044, 2.266736030457, 2.267378215022, 2.268026985744, 2.268673824651, 2.269312262288, 2.269946303092, 2.270585238833, 2.271204598689, 2.271847744898, 2.272483548889, 2.273120692320, 2.273763421513, 2.274408001954, 2.275065649429, 2.275708643570, 2.276349636920, 2.276987222122, 2.277627802289, 2.278268916581, 2.278902311469, 2.279535722167, 2.280176846498, 2.280819664868, 2.281463519144, 2.282097186635, 2.282751683291, 2.283388902500, 2.284040923251, 2.284680958452, 2.285324200175, 2.285965962989, 2.286612960690, 2.287244346811, 2.287895275744, 2.288534184689, 2.289179190682, 2.289822955244, 2.290466403914, 2.291087291346, 2.291733892664, 2.292371240676, 2.293009184325, 2.293637477591, 2.294283699284, 2.294919921220, 2.295561708230, 2.296204616855, 2.296849166835, 2.297495967156, 2.298142524234, 2.298782354706, 2.299426244592, 2.300065370527, 2.300703268569, 2.301349667165, 2.301990935654, 2.302627747291, 2.303259819137, 2.303889664651, 2.304543891599, 2.305180601953, 2.305811836200, 2.306462370881, 2.307107803699, 2.307748551491, 2.308376287977, 2.309016168777, 2.309643172135, 2.310294595293, 2.310939621285, 2.311586140688, 2.312232108791, 2.312865382367, 2.313513436220, 2.314144285247, 2.314782506965, 2.315414844329, 2.316063119798, 2.316706150811, 2.317359965608, 2.318013050047, 2.318656715411, 2.319290918310, 2.319923145614, 2.320558384260, 2.321201923427, 2.321839947882, 2.322490592015, 2.323124207771, 2.323770922383, 2.324410718357, 2.325060546409, 2.325705096499, 2.326334123754, 2.326974574757, 2.327601751155, 2.328242873449, 2.328883647021, 2.329525831156, 2.330170916909, 2.330815101771, 2.331460243567, 2.332101120266, 2.332753876564, 2.333406305305, 2.334048843969, 2.334699186449, 2.335343172203, 2.335992444741, 2.336636184016, 2.337277763042, 2.337925492066, 2.338563771510, 2.339210198239, 2.339874210769, 2.340507752764, 2.341152508891, 2.341784389783, 2.342432765257, 2.343074550723, 2.343721598373, 2.344389094661, 2.345019458010, 2.345665272780, 2.346307036341, 2.346945212014, 2.347576881992, 2.348215668990, 2.348867712369, 2.349497136867, 2.350142841473, 2.350792721858, 2.351422018138, 2.352076454637, 2.352711822066, 2.353361249759, 2.354010276237, 2.354652510415, 2.355282901421, 2.355915194385, 2.356566571082, 2.357194608308, 2.357832266599, 2.358467690003, 2.359106625913, 2.359743221935, 2.360373383722, 2.361014333743, 2.361657229818, 2.362291977303, 2.362912529356, 2.363555234340, 2.364202206977, 2.364864632713, 2.365506204795, 2.366137020687, 2.366792300373, 2.367425494494, 2.368074208473, 2.368715975202, 2.369347001136, 2.370000222759, 2.370651981288, 2.371294916483, 2.371930726112, 2.372568389645, 2.373234683751, 2.373876217843, 2.374522301585, 2.375172129669, 2.375811994301, 2.376446499824, 2.377100871400, 2.377747524499, 2.378389329363, 2.379024392007, 2.379675479208, 2.380317848033, 2.380958349315, 2.381594150166, 2.382237794531, 2.382880716245, 2.383515033614, 2.384153223862, 2.384788350555, 2.385430209394, 2.386070377094, 2.386723974691, 2.387373259189, 2.388029883244, 2.388641482624, 2.389281721978, 2.389921094582, 2.390546786108, 2.391177549301, 2.391820470511, 2.392467346883, 2.393092103392, 2.393728834970, 2.394343673060, 2.394991086297, 2.395649942217, 2.396278648256, 2.396917905991, 2.397562987612, 2.398212179431, 2.398846784370, 2.399478831269, 2.400116164092, 2.400767000910, 2.401420018457, 2.402070621418, 2.402722420042, 2.403357718109, 2.404001323478, 2.404643789064, 2.405281906060, 2.405914769188, 2.406562953262, 2.407208002303, 2.407840903532, 2.408461379904, 2.409102350043, 2.409744379146, 2.410386464940, 2.411024689943, 2.411665423075, 2.412305419526, 2.412947259625, 2.413582058191, 2.414240894260, 2.414888425645, 2.415526522251, 2.416151970750, 2.416791360969, 2.417422268744, 2.418057392288, 2.418695609891, 2.419339557156, 2.419985945601, 2.420629064844, 2.421269242613, 2.421897515344, 2.422540371316, 2.423172558188, 2.423811543754, 2.424444546434, 2.425089568402, 2.425742147577, 2.426384116382, 2.427028660868, 2.427676256177, 2.428320393817, 2.428955109749, 2.429598111982, 2.430244172867, 2.430878429072, 2.431509038243, 2.432156775495, 2.432785127397, 2.433435007763, 2.434083030156, 2.434712168405, 2.435361864029, 2.435997007013, 2.436644118629, 2.437290650577, 2.437933027411, 2.438576713479, 2.439215982015, 2.439855475472, 2.440491600940, 2.441133336663, 2.441778544259, 2.442411119916, 2.443043172829, 2.443693880127, 2.444339523015, 2.444977536556, 2.445616367598, 2.446267789495, 2.446904876713, 2.447531824385, 2.448171866922, 2.448823351978, 2.449456133480, 2.450098897940, 2.450736362035, 2.451368009849, 2.452011398823, 2.452643057935, 2.453294383287, 2.453926059868, 2.454552843051, 2.455201962104, 2.455835056749, 2.456485848051, 2.457134380884, 2.457784755923, 2.458432487155, 2.459077936238, 2.459716835735, 2.460366328943, 2.461000599866, 2.461635295575, 2.462278348881, 2.462925508208, 2.463564667062, 2.464202238562, 2.464854046496, 2.465500745394, 2.466154760648, 2.466820577084, 2.467455686935, 2.468091854536, 2.468731766876, 2.469378766782, 2.470020707718, 2.470650762579, 2.471306029115, 2.471960869571, 2.472608962334, 2.473245497817, 2.473896417162, 2.474538340740, 2.475190294167, 2.475846215602, 2.476481792809, 2.477108920932, 2.477759006475, 2.478404186216, 2.479053204963, 2.479684059967, 2.480321082723, 2.480972185914, 2.481625714830, 2.482263880157, 2.482917112520, 2.483554138331, 2.484210110348, 2.484853015801, 2.485482795121, 2.486102981033, 2.486753492335, 2.487385368604, 2.488012287254, 2.488664061304, 2.489295509686, 2.489941699045, 2.490572052640, 2.491207898307, 2.491854111052, 2.492490353120, 2.493125771317, 2.493767806451, 2.494417164436, 2.495055545880, 2.495702346116, 2.496321104297, 2.496952746205, 2.497608526403, 2.498248063310, 2.498890735145, 2.499532438781, 2.500168772030, 2.500818147196, 2.501476900740, 2.502127546240, 2.502766590356, 2.503399932193, 2.504055962512, 2.504697574955, 2.505331099510, 2.505965131859, 2.506602043346, 2.507236399161, 2.507875179121, 2.508503415637, 2.509134245314, 2.509773717956, 2.510420463724, 2.511065215562, 2.511706270085, 2.512351239502, 2.513004951203, 2.513631304461, 2.514258988131, 2.514892554721, 2.515539999602, 2.516192260049, 2.516823944210, 2.517478137303, 2.518118568829, 2.518752919698, 2.519403420985, 2.520037351587, 2.520688052595, 2.521326603208, 2.521977506639, 2.522616365923, 2.523271959743, 2.523918676998, 2.524567812016, 2.525213552365, 2.525859671228, 2.526489381602, 2.527127169508, 2.527777754232, 2.528428141999, 2.529066289239, 2.529715375629, 2.530376774444, 2.531023103690, 2.531670839431, 2.532318655015, 2.532953510875, 2.533576535690, 2.534244588545, 2.534888667461, 2.535528486445, 2.536161487794, 2.536790480340, 2.537433857861, 2.538080738581, 2.538709215903, 2.539352136344, 2.539986673823, 2.540640236356, 2.541264577854, 2.541921580304, 2.542563218321, 2.543212632338, 2.543843875973, 2.544477712061, 2.545116893284, 2.545771818625, 2.546400835638, 2.547029693716, 2.547694409590, 2.548323456416, 2.548970786699, 2.549592141404, 2.550237204067, 2.550879983735, 2.551510263119, 2.552155705451, 2.552800712642, 2.553421983829, 2.554048500041, 2.554680905774, 2.555322813309, 2.555970670240, 2.556608855892, 2.557260985918, 2.557907976526, 2.558546031500, 2.559189432131, 2.559830792549, 2.560471522977, 2.561101817944, 2.561753135244, 2.562394173288, 2.563032824375, 2.563663669901, 2.564295592360, 2.564933539109, 2.565584403427, 2.566245202682, 2.566888905483, 2.567518482504, 2.568163273336, 2.568805482813, 2.569439618998, 2.570085979461, 2.570733141708, 2.571374633153, 2.572012048256, 2.572663550287, 2.573313917434, 2.573960700608, 2.574583498389, 2.575222541231, 2.575871193942, 2.576515739248, 2.577146479353, 2.577771237482, 2.578401172830, 2.579030046260, 2.579666266098, 2.580313828922, 2.580952760988, 2.581591308670, 2.582229966775, 2.582860922683, 2.583488136094, 2.584126091262, 2.584785852398, 2.585437254355, 2.586093486091, 2.586732934997, 2.587384411077, 2.588050658357, 2.588672783793, 2.589314863683, 2.589935594252, 2.590572101129, 2.591232585527, 2.591883554745, 2.592523434929, 2.593168514203, 2.593814382621, 2.594448581018, 2.595110888719, 2.595748183094, 2.596399959823, 2.597038634855, 2.597713334603, 2.598331210132, 2.598976873409, 2.599620215860, 2.600262609732, 2.600912192509, 2.601547478932, 2.602165624592, 2.602795091932, 2.603439588630, 2.604082250723, 2.604716077240, 2.605365183587, 2.605999834802, 2.606616103904, 2.607254696807, 2.607885426817, 2.608520953230, 2.609174892786, 2.609804882843, 2.610449247952, 2.611094747898, 2.611738011296, 2.612386142621, 2.613022415043, 2.613664438387, 2.614307412249, 2.614956171168, 2.615591025533, 2.616232373144, 2.616867838981, 2.617521336971, 2.618151662236, 2.618773876894, 2.619413979142, 2.620059190784, 2.620700647643, 2.621327254009, 2.621955856924, 2.622587556512, 2.623224188902, 2.623868697080, 2.624503553010, 2.625141903135, 2.625770735249, 2.626410584571, 2.627058922356, 2.627718365313, 2.628342083060, 2.628980374815, 2.629614053130, 2.630252364756, 2.630888460031, 2.631527533387, 2.632159728710, 2.632796201878, 2.633455830823, 2.634090093378, 2.634738581233, 2.635382036535, 2.636023816713, 2.636666734814, 2.637307026385, 2.637944867069, 2.638591015386, 2.639226204266, 2.639859291268, 2.640498806491, 2.641152190064, 2.641797230249, 2.642445708316, 2.643103175354, 2.643731809060, 2.644366907262, 2.645024031104, 2.645669666165, 2.646303567317, 2.646954961520, 2.647597302154, 2.648239628208, 2.648877487881, 2.649508922337, 2.650151949367, 2.650814198995, 2.651460136221, 2.652107035598, 2.652767004219, 2.653398648843, 2.654056864204, 2.654690388980, 2.655321893487, 2.655972018278, 2.656618784353, 2.657262372498, 2.657892296870, 2.658528082199, 2.659159450188, 2.659801459860, 2.660459522545, 2.661087736427, 2.661734797464, 2.662367055651, 2.663024622538, 2.663676379759, 2.664306861679, 2.664949704891, 2.665576208159, 2.666217913793, 2.666862585723, 2.667499127362, 2.668131951272, 2.668763267943, 2.669409500700, 2.670062791370, 2.670704451253, 2.671349302006, 2.671976133567, 2.672605506015, 2.673251345842, 2.673886464607, 2.674543861226, 2.675195470287, 2.675837558289, 2.676476267058, 2.677110547826, 2.677738311922, 2.678363050130, 2.679004449716, 2.679639320928, 2.680279073565, 2.680900814470, 2.681542429471, 2.682182904636, 2.682823488847, 2.683483458478, 2.684115472897, 2.684767742785, 2.685382057237, 2.686003775383, 2.686641370515, 2.687265952221, 2.687916200499, 2.688537319765, 2.689162088484, 2.689791796780, 2.690446479540, 2.691075495426, 2.691697309037, 2.692362785844, 2.692999083197, 2.693642963623, 2.694280711118, 2.694925205088, 2.695552559134, 2.696182546787, 2.696815610349, 2.697442240923, 2.698067609979, 2.698690408710, 2.699349309009, 2.699990491496, 2.700610170635, 2.701236410567, 2.701876452346, 2.702495108675, 2.703141394363, 2.703783812728, 2.704417727483, 2.705033851439, 2.705661656015, 2.706299644755, 2.706942995484, 2.707568695156, 2.708201286482, 2.708841020437, 2.709480363335, 2.710114410617, 2.710771697254, 2.711401824188, 2.712057706771, 2.712706961049, 2.713339006406, 2.713976693335, 2.714621621170, 2.715265253646, 2.715919550256, 2.716569859257, 2.717231787578, 2.717852540068, 2.718476679437, 2.719098305331, 2.719721277813, 2.720380502478, 2.721032961479, 2.721644532374, 2.722293395479, 2.722938640085, 2.723579788983, 2.724221195367, 2.724849721396, 2.725478927553, 2.726106272839, 2.726736377514, 2.727376670633, 2.727998871797, 2.728663351031, 2.729298809326, 2.729939163119, 2.730575791676, 2.731204467389, 2.731832883396, 2.732472296052, 2.733117584591, 2.733771832261, 2.734412458849, 2.735028548890, 2.735666780019, 2.736316836090, 2.736961230622, 2.737585696333, 2.738217003325, 2.738858512654, 2.739504785210, 2.740157512111, 2.740810743294, 2.741448914370, 2.742077472948, 2.742714867744, 2.743338768886, 2.743966457823, 2.744595537678, 2.745225046981, 2.745865630093, 2.746519515171, 2.747143814938, 2.747783591896, 2.748416282671, 2.749047703324, 2.749693465846, 2.750337012639, 2.750974171302, 2.751624522797, 2.752262838083, 2.752890291851, 2.753549430602, 2.754203159699, 2.754835894644, 2.755490822560, 2.756140051692, 2.756774624672, 2.757421553898, 2.758063477029, 2.758694141365, 2.759334456239, 2.759982463829, 2.760625182947, 2.761242036821, 2.761874074355, 2.762506530331, 2.763123295236, 2.763756818137, 2.764396315566, 2.765028918518, 2.765671052828, 2.766323773763, 2.766961478124, 2.767593253599, 2.768221874463, 2.768865434856, 2.769504586054, 2.770150307282, 2.770789047850, 2.771421545212, 2.772051881752, 2.772692655305, 2.773315821489, 2.773958722075, 2.774594563713, 2.775238843791, 2.775905337640, 2.776499126851, 2.777135581863, 2.777754486521, 2.778379488668, 2.779013485420, 2.779665405362, 2.780299450410, 2.780914491789, 2.781557980561, 2.782210051809, 2.782831494428, 2.783465433935, 2.784076261290, 2.784720487181, 2.785374148337, 2.785991647910, 2.786625969117, 2.787271862301, 2.787906192017, 2.788543317795, 2.789184854281, 2.789797093353, 2.790442630407, 2.791102282971, 2.791723685321, 2.792355939848, 2.793005564105, 2.793616736362, 2.794250673194, 2.794896911037, 2.795549807671, 2.796207218603, 2.796858280888, 2.797484983096, 2.798134963847, 2.798756954636, 2.799394340156, 2.800046090761, 2.800692233874, 2.801332193251, 2.801981355189, 2.802609434539, 2.803238423538, 2.803872471937, 2.804502186679, 2.805110356789, 2.805744646988, 2.806391822378, 2.807046647651, 2.807700509286, 2.808316250414, 2.808960558536, 2.809575008585, 2.810181352433, 2.810810736764, 2.811456789632, 2.812084645350, 2.812714820878, 2.813344216707, 2.813989241318, 2.814614536715, 2.815276211191, 2.815938895307, 2.816588926783, 2.817217979142, 2.817860506909, 2.818503414792, 2.819136393789, 2.819760546152, 2.820375832452, 2.821007809526, 2.821649060286, 2.822299624503, 2.822945964690, 2.823577355952, 2.824194021331, 2.824844639900, 2.825484611185, 2.826126690962, 2.826758646721, 2.827415456647, 2.828078523193, 2.828703959683, 2.829322969071, 2.829954311796, 2.830600098598, 2.831235068729, 2.831871263479, 2.832488014959, 2.833142610908, 2.833792863459, 2.834440828051, 2.835078174864, 2.835694740739, 2.836333037944, 2.836937366655, 2.837593330156, 2.838244899459, 2.838896848384, 2.839563584450, 2.840201284867, 2.840840826156, 2.841478597265, 2.842102512173, 2.842743047243, 2.843375444248, 2.844006337510, 2.844626305066, 2.845260844097, 2.845886870172, 2.846517154953, 2.847150188457, 2.847811065378, 2.848441700097, 2.849084603691, 2.849742287548, 2.850373580814, 2.851001478712, 2.851642013194, 2.852278239496, 2.852917256386, 2.853555045011, 2.854213950735, 2.854851471507, 2.855502384407, 2.856133381246, 2.856785282839, 2.857410641725, 2.858042226815, 2.858660930926, 2.859295594666, 2.859917660997, 2.860549755029, 2.861191603816, 2.861846724838, 2.862488279730, 2.863132368229, 2.863753611503, 2.864378922775, 2.864992723514, 2.865632571572, 2.866271767735, 2.866890808623, 2.867513614120, 2.868142443177, 2.868777962780, 2.869408948037, 2.870046647013, 2.870666904170, 2.871303548072, 2.871936599193, 2.872557296449, 2.873189260192, 2.873834811983, 2.874470590464, 2.875081892870, 2.875720479607, 2.876352819444, 2.876995896462, 2.877631735685, 2.878265553981, 2.878910485756, 2.879519517963, 2.880126110016, 2.880782393686, 2.881414550289, 2.882067489031, 2.882696216354, 2.883331830861, 2.883974361494, 2.884612849660, 2.885246942658, 2.885871944784, 2.886490156247, 2.887152785997, 2.887798986344, 2.888440775006, 2.889070057332, 2.889693851985, 2.890310784532, 2.890931635287, 2.891572659970, 2.892211582634, 2.892844999226, 2.893491235269, 2.894118694639, 2.894756605064, 2.895391357867, 2.896036269298, 2.896683851561, 2.897333429467, 2.897972992555, 2.898611779524, 2.899236010695, 2.899899421773, 2.900531379485, 2.901165641717, 2.901808799015, 2.902449787867, 2.903125774332, 2.903763839662, 2.904389949973, 2.905027782767, 2.905700107281, 2.906344769242, 2.906988286140, 2.907620821128, 2.908264475004, 2.908904154553, 2.909533846046, 2.910150679565, 2.910779353186, 2.911414604694, 2.912067456977, 2.912699623848, 2.913337692463, 2.913982756290, 2.914618432869, 2.915268261745, 2.915912980969, 2.916518522082, 2.917186634743, 2.917807253585, 2.918447117186, 2.919081075989, 2.919750618393, 2.920363345516, 2.920993592670, 2.921621129938, 2.922233961989, 2.922838206355, 2.923490264781, 2.924125798774, 2.924760438090, 2.925426364199, 2.926067304016, 2.926718362594, 2.927391341706, 2.928031510455, 2.928676309414, 2.929327972260, 2.929947718751, 2.930589448739, 2.931231016293, 2.931857941103, 2.932509194129, 2.933155840077, 2.933802704497, 2.934429247589, 2.935072777046, 2.935748724231, 2.936378833351, 2.937029392736, 2.937667759892, 2.938313095333, 2.938952221551, 2.939594179223, 2.940229896513, 2.940859344594, 2.941505647870, 2.942157095888, 2.942804573359, 2.943448061041, 2.944096703801, 2.944724519901, 2.945360137969, 2.946003974492, 2.946647230555, 2.947277977940, 2.947906561079, 2.948539912973, 2.949177666847, 2.949815584826, 2.950451341442, 2.951087254024, 2.951725264919, 2.952343199364, 2.952963183432, 2.953594591328, 2.954226918547, 2.954860950612, 2.955495517668, 2.956119628941, 2.956756432259, 2.957382359495, 2.958030875531, 2.958676807518, 2.959329238186, 2.959985422926, 2.960619197843, 2.961244762879, 2.961887142242, 2.962506171712, 2.963134861612, 2.963752876732, 2.964396178646, 2.965033622938, 2.965655551725, 2.966296455991, 2.966940722228, 2.967570629610, 2.968192975754, 2.968823086524, 2.969435897309, 2.970072679784, 2.970718922778, 2.971347832934, 2.971956890633, 2.972581073756, 2.973224121297, 2.973848494041, 2.974466804208, 2.975103219696, 2.975738104962, 2.976375153613, 2.977027143104, 2.977678050187, 2.978321671249, 2.978954662388, 2.979588577468, 2.980251222497, 2.980884125146, 2.981530021451, 2.982143534494, 2.982791725373, 2.983424163571, 2.984072595436, 2.984692648148, 2.985311907979, 2.985940040734, 2.986597716935, 2.987231930602, 2.987863270836, 2.988489186126, 2.989146501166, 2.989774694140, 2.990404646699, 2.991030834728, 2.991669855992, 2.992290192430, 2.992922952476, 2.993548933681, 2.994180960602, 2.994834507822, 2.995485601572, 2.996128633566, 2.996758393693, 2.997381729654, 2.998003800046, 2.998621134783, 2.999256690016, 2.999883624672, 3.000547556080, 3.001172874551, 3.001791245004, 3.002399579525, 3.003027572257, 3.003669175597, 3.004291113197, 3.004922727546, 3.005557021363, 3.006174181533, 3.006806336687, 3.007461945624, 3.008103058909, 3.008734041952, 3.009390794470, 3.010021874860, 3.010643191828, 3.011251582241, 3.011891176863, 3.012522774735, 3.013189762981, 3.013859570669, 3.014482814103, 3.015123591229, 3.015763063441, 3.016411145893, 3.017015481477, 3.017659556792, 3.018267445099, 3.018901133465, 3.019502132666, 3.020108058572, 3.020743073663, 3.021363507807, 3.022000362591, 3.022634492389, 3.023275964287, 3.023897276929, 3.024551647254, 3.025178470636, 3.025788225972, 3.026400684668, 3.027049134360, 3.027664763657, 3.028322058275, 3.028969671269, 3.029585241501, 3.030245916987, 3.030876356183, 3.031546939407, 3.032157290276, 3.032828450623, 3.033430756370, 3.034081340627, 3.034708908443, 3.035327020145, 3.035972904553, 3.036578172729, 3.037224029284, 3.037838626084, 3.038465007564, 3.039088016995, 3.039721439116, 3.040349114302, 3.040996311804, 3.041632525356, 3.042277810087, 3.042930287416, 3.043572703409, 3.044205011909, 3.044863764150, 3.045503261439, 3.046164953359, 3.046818947741, 3.047448252804, 3.048059551211, 3.048706204326, 3.049333873145, 3.049945884432, 3.050594378553, 3.051195465437, 3.051819404798, 3.052494718552, 3.053176483499, 3.053797870232, 3.054427531448, 3.055078813182, 3.055718729256, 3.056347227927, 3.056975151931, 3.057597538567, 3.058237703091, 3.058846983325, 3.059494472177, 3.060126966474, 3.060792353116, 3.061467766522, 3.062095126473, 3.062754505058, 3.063376189883, 3.064029968484, 3.064673139882, 3.065282436824, 3.065928983601, 3.066533971842, 3.067156532122, 3.067775924943, 3.068370276274, 3.069006167023, 3.069625146443, 3.070255730785, 3.070908195651, 3.071529381780, 3.072168892744, 3.072783669409, 3.073420402784, 3.074058586109, 3.074693582218, 3.075325892289, 3.075974126305, 3.076638872396, 3.077270391064, 3.077908025623, 3.078565852683, 3.079177250557, 3.079812995774, 3.080471103310, 3.081104035542, 3.081735270399, 3.082402599886, 3.083038881181, 3.083698739378, 3.084320048133, 3.084916370267, 3.085578031667, 3.086205214622, 3.086846034828, 3.087461241943, 3.088091152404, 3.088714519687, 3.089343584412, 3.089971958847, 3.090615690083, 3.091243229985, 3.091904948853, 3.092524682248, 3.093172211220, 3.093815855939, 3.094443182257, 3.095101148668, 3.095733042638, 3.096361519863, 3.096981677510, 3.097603265723, 3.098236635307, 3.098904741391, 3.099567868607, 3.100192625060, 3.100805135143, 3.101419607371, 3.102021217643, 3.102656117054, 3.103294700521, 3.103914914577, 3.104558115578, 3.105223851114, 3.105842395074, 3.106468480585, 3.107082687903, 3.107735053655, 3.108366662474, 3.109009796976, 3.109613077728, 3.110258062036, 3.110923628690, 3.111529575682, 3.112170667210, 3.112804259854, 3.113438214272, 3.114055023774, 3.114691373500, 3.115338285685, 3.115960637665, 3.116629895309, 3.117250119973, 3.117882626054, 3.118500079308, 3.119127553376, 3.119750213476, 3.120402991800, 3.121045274800, 3.121648853711, 3.122292983245, 3.122926541743, 3.123589888063, 3.124244420742, 3.124862307505, 3.125476436373, 3.126101305183, 3.126718934154, 3.127364810117, 3.127980157952, 3.128617985843, 3.129287749771, 3.129944489435, 3.130591663972, 3.131212779180, 3.131814781290, 3.132434114634, 3.133086782664, 3.133709115669, 3.134324057733, 3.134959426316, 3.135582670219, 3.136199084235, 3.136823515835, 3.137455998384, 3.138086419090, 3.138720147154, 3.139306317152, 3.139908263285, 3.140529651463, 3.141183788697, 3.141813025861, 3.142444381785, 3.143043451606, 3.143646371660, 3.144278584106, 3.144905655144, 3.145522703795, 3.146150358988, 3.146830684242, 3.147455354975, 3.148083979275, 3.148715961381, 3.149323138641, 3.149950793936, 3.150587957601, 3.151240822285, 3.151860783539, 3.152536545103, 3.153159595978, 3.153768071635, 3.154404667417, 3.155008685457, 3.155616030496, 3.156254720459, 3.156865058221, 3.157506837357, 3.158180196165, 3.158830183670, 3.159443527248, 3.160068411503, 3.160724377251, 3.161363702741, 3.162025413394, 3.162681185768, 3.163301261173, 3.163923490057, 3.164528850192, 3.165148394990, 3.165766916447, 3.166409891392, 3.167043611384, 3.167655894997, 3.168301674758, 3.168900357727, 3.169503716924, 3.170113055782, 3.170726467986, 3.171357501367, 3.172001068857, 3.172600353245, 3.173217291175, 3.173848716349, 3.174468729483, 3.175094828594, 3.175714021037, 3.176390808382, 3.177017728100, 3.177671052173, 3.178318157787, 3.178989837345, 3.179590313057, 3.180221871802, 3.180828006081, 3.181463347228, 3.182091693065, 3.182749393143, 3.183363042421, 3.183987510801, 3.184604906262, 3.185226507101, 3.185868986685, 3.186523761750, 3.187164155135, 3.187792778608, 3.188427004757, 3.189011150897, 3.189626997691, 3.190243719020, 3.190892995487, 3.191529743555, 3.192154582614, 3.192792507450, 3.193421201143, 3.194051485212, 3.194706483687, 3.195329784770, 3.195941707180, 3.196566101973, 3.197196866786, 3.197815536572, 3.198461837308, 3.199112535792, 3.199782097397, 3.200452003847, 3.201121565052, 3.201754839275, 3.202398035753, 3.203044958832, 3.203667856950, 3.204280527567, 3.204905897761, 3.205544021212, 3.206160044331, 3.206791624577, 3.207427625551, 3.208066662919, 3.208694001684, 3.209325764190, 3.209929572704, 3.210544800155, 3.211136889279, 3.211811824350, 3.212439642153, 3.213064112844, 3.213711504232, 3.214353459111, 3.214977127415, 3.215632373647, 3.216251451524, 3.216869266421, 3.217491544434, 3.218121891689, 3.218786215862, 3.219401171920, 3.220003305423, 3.220634424330, 3.221255619038, 3.221876979675, 3.222500680410, 3.223126004099, 3.223766769782, 3.224368438442, 3.225000833742, 3.225618818141, 3.226253769333, 3.226861824658, 3.227502263626, 3.228148789196, 3.228784510887, 3.229414535565, 3.230037361849, 3.230668469294, 3.231308632509, 3.231959371977, 3.232584377310, 3.233196909430, 3.233860904163, 3.234484182849, 3.235091193580, 3.235724462178, 3.236376617879, 3.237007267745, 3.237674115543, 3.238320185986, 3.238950654950, 3.239595612557, 3.240244549956, 3.240879333267, 3.241515803118, 3.242175203243, 3.242795345900, 3.243403444016, 3.244036770583, 3.244685517004, 3.245305437193, 3.245951492123, 3.246549470930, 3.247158248493, 3.247787089485, 3.248446083816, 3.249094519138, 3.249747011450, 3.250397393527, 3.251032493908, 3.251655345017, 3.252244155868, 3.252868748382, 3.253480227047, 3.254125312693, 3.254776042948, 3.255440262028, 3.256096883081, 3.256713711401, 3.257366763506, 3.258009785533, 3.258653761038, 3.259335777949, 3.259980143614, 3.260554244338, 3.261198836792, 3.261852323957, 3.262494873164, 3.263094595684, 3.263704712084, 3.264359588843, 3.265007460006, 3.265629878574, 3.266255595704, 3.266865354878, 3.267493659904, 3.268128511714, 3.268765099306, 3.269396160763, 3.270028140533, 3.270613258286, 3.271262422585, 3.271868698671, 3.272503475980, 3.273139182456, 3.273772557813, 3.274390519703, 3.275037177626, 3.275672509783, 3.276309593262, 3.276987879571, 3.277593977350, 3.278239657077, 3.278843378901, 3.279434717664, 3.280084797319, 3.280698551067, 3.281337245154, 3.281968566932, 3.282594147820, 3.283224799600, 3.283894776149, 3.284540700271, 3.285178374375, 3.285826211838, 3.286437220505, 3.287038156207, 3.287671087744, 3.288312534687, 3.288902557766, 3.289477313256, 3.290109577060, 3.290734280202, 3.291353087585, 3.292017014621, 3.292615503191, 3.293224201843, 3.293864511265, 3.294486086902, 3.295117121783, 3.295721616851, 3.296382808268, 3.297011443585, 3.297621168356, 3.298252463095, 3.298831092979, 3.299434725471, 3.300085130867, 3.300707868653, 3.301313248017, 3.301922953912, 3.302533516977, 3.303183345807, 3.303826280556, 3.304512194397, 3.305115012545, 3.305788915641, 3.306423409807, 3.307064116727, 3.307700478037, 3.308340423153, 3.308978658732, 3.309635552711, 3.310260612374, 3.310942553794, 3.311557929874, 3.312198241295, 3.312812723031, 3.313448631990, 3.314061306506, 3.314696358582, 3.315354782957, 3.315979144352, 3.316631414581, 3.317280157384, 3.317883817360, 3.318516349206, 3.319168821569, 3.319770579824, 3.320424032397, 3.321067554285, 3.321731161585, 3.322353816567, 3.322961834418, 3.323569789837, 3.324175849088, 3.324786424965, 3.325390510779, 3.326007397966, 3.326660174135, 3.327266873714, 3.327916002748, 3.328550371198, 3.329221812986, 3.329850667536, 3.330483222589, 3.331105529707, 3.331762290619, 3.332364028220, 3.332989038161, 3.333620566506, 3.334232386878, 3.334876994890, 3.335479305420, 3.336072094954, 3.336743959656, 3.337380975909, 3.337965965850, 3.338604778367, 3.339222717385, 3.339844386263, 3.340446021397, 3.341061825359, 3.341703303557, 3.342339043759, 3.342970932651, 3.343605658370, 3.344256664693, 3.344868291040, 3.345526005697, 3.346103769616, 3.346719934559, 3.347347604600, 3.347956828653, 3.348560125145, 3.349195314619, 3.349820743843, 3.350449021585, 3.351069437509, 3.351712214781, 3.352343236804, 3.352970282315, 3.353596273777, 3.354201570555, 3.354840155230, 3.355491495733, 3.356086623585, 3.356760574971, 3.357400958340, 3.358047239583, 3.358661751713, 3.359286073991, 3.359950090403, 3.360544390565, 3.361154468232, 3.361794375029, 3.362420218757, 3.363050973345, 3.363676625094, 3.364293131401, 3.364918563989, 3.365517691649, 3.366155991619, 3.366775019564, 3.367384811256, 3.368014715336, 3.368632340661, 3.369268123553, 3.369922142570, 3.370572050814, 3.371218849277, 3.371846163524, 3.372460050801, 3.373086085712, 3.373708917065, 3.374317218793, 3.374939760050, 3.375550820618, 3.376140024784, 3.376748641538, 3.377389177962, 3.377996437996, 3.378602471422, 3.379243669935, 3.379920188241, 3.380517440380, 3.381101936986, 3.381769851716, 3.382374892272, 3.382928331344, 3.383549693778, 3.384179309548, 3.384843549598, 3.385471880252, 3.386116969654, 3.386698476311, 3.387298773017, 3.387915814772, 3.388522047412, 3.389123807855, 3.389694442578, 3.390341571737, 3.390996077642, 3.391629099299, 3.392281263183, 3.392883965979, 3.393514373730, 3.394146774184, 3.394785486214, 3.395440251777, 3.396035466596, 3.396643405628, 3.397220761480, 3.397814081531, 3.398471258469, 3.399073910501, 3.399669769419, 3.400307928550, 3.400963426289, 3.401626484792, 3.402210507836, 3.402838139182, 3.403493071068, 3.404128066579, 3.404777227268, 3.405418522896, 3.406028687427, 3.406650787607, 3.407310389395, 3.407946550869, 3.408563617898, 3.409150366260, 3.409775841668, 3.410398867175, 3.411021668800, 3.411630797986, 3.412256492169, 3.412870728092, 3.413494836508, 3.414142382658, 3.414755907737, 3.415413250421, 3.416036498211, 3.416726393221, 3.417385592001, 3.417988937479, 3.418683082335, 3.419314465229, 3.419958189232, 3.420607444412, 3.421250797972, 3.421883633082, 3.422484073362, 3.423115257800, 3.423746208701, 3.424353846692, 3.425017800609, 3.425690872709, 3.426308193661, 3.426956571946, 3.427574533505, 3.428208508744, 3.428855069057, 3.429452390310, 3.430038843851, 3.430626090385, 3.431252816092, 3.431852272658, 3.432451382260, 3.433067800034, 3.433706314316, 3.434304447114, 3.434906951416, 3.435586065391, 3.436247268878, 3.436851288904, 3.437501340624, 3.438160704359, 3.438786477859, 3.439398820199, 3.439974947163, 3.440596155368, 3.441221851947, 3.441855658829, 3.442530093682, 3.443162199405, 3.443783160203, 3.444384468829, 3.444973301806, 3.445596859416, 3.446212819865, 3.446819934438, 3.447446150994, 3.448083020635, 3.448696418104, 3.449324126105, 3.449934384942, 3.450530796195, 3.451148890375, 3.451796133428, 3.452407416599, 3.453025724167, 3.453672069258, 3.454305779261, 3.454916893722, 3.455547465072, 3.456151640127, 3.456797684394, 3.457430994466, 3.458061488572, 3.458692899340, 3.459298968169, 3.459889605052, 3.460528697373, 3.461195106124, 3.461788329142, 3.462389920097, 3.463016308598, 3.463641075717, 3.464271802478, 3.464893312769, 3.465482733292, 3.466097088505, 3.466744115880, 3.467319495560, 3.467969623763, 3.468550462696, 3.469187085791, 3.469814394177, 3.470404121660, 3.471044754585, 3.471659314467, 3.472277322088, 3.472866535205, 3.473503060973, 3.474095234384, 3.474715424896, 3.475337799874, 3.475963666862, 3.476589135807, 3.477255906881, 3.477901513807, 3.478550696248, 3.479147179230, 3.479760212188, 3.480341295977, 3.480966533654, 3.481611102291, 3.482205215185, 3.482793541826, 3.483436863874, 3.484054662472, 3.484683947194, 3.485355305569, 3.485966533268, 3.486587943490, 3.487246250251, 3.487872166643, 3.488432118187, 3.489069121725, 3.489692308089, 3.490333850317, 3.490892952275, 3.491520109921, 3.492111759735, 3.492750135455, 3.493331292422, 3.493967401330, 3.494603086909, 3.495220687295, 3.495847328872, 3.496485774132, 3.497080139811, 3.497742260575, 3.498372553752, 3.498990061092, 3.499596098770, 3.500241457081, 3.500796958092, 3.501433081340, 3.502050819229, 3.502633511013, 3.503288937270, 3.503946740071, 3.504551412047, 3.505136081119, 3.505775813083, 3.506413701194, 3.507058111018, 3.507586067558, 3.508246217228, 3.508866721063, 3.509468461448, 3.510065414456, 3.510671634192, 3.511313940077, 3.511926141092, 3.512519416500, 3.513147475857, 3.513746675740, 3.514386450905, 3.514990207010, 3.515566333139, 3.516247289020, 3.516850784025, 3.517508020267, 3.518134749235, 3.518755214041, 3.519416777014, 3.520077910943, 3.520678114195, 3.521293571588, 3.521908457990, 3.522476485412, 3.523075670085, 3.523662630011, 3.524264906822, 3.524888380624, 3.525505468560, 3.526142395534, 3.526730596617, 3.527353234856, 3.527976767039, 3.528636401548, 3.529232399423, 3.529860108491, 3.530451896734, 3.531043017352, 3.531646760860, 3.532223240385, 3.532831590386, 3.533452660201, 3.534108786224, 3.534700445892, 3.535286954243, 3.535911547869, 3.536563932089, 3.537196350903, 3.537723321986, 3.538307929660, 3.538993972786, 3.539575781365, 3.540179460645, 3.540749284272, 3.541378767777, 3.541994036851, 3.542623813464, 3.543243884480, 3.543860284143, 3.544494296348, 3.545123140524, 3.545711698165, 3.546359117306, 3.547009033377, 3.547627738947, 3.548259605200, 3.548878559242, 3.549521484629, 3.550156114176, 3.550783953794, 3.551425070185, 3.552085714248, 3.552700845882, 3.553285796158, 3.553919734765, 3.554492317336, 3.555162325424, 3.555724058470, 3.556331860067, 3.556926421883, 3.557550021109, 3.558176087260, 3.558733873267, 3.559378969565, 3.559947763358, 3.560531511673, 3.561136597673, 3.561780525573, 3.562411138948, 3.563039493514, 3.563663987938, 3.564227279411, 3.564878998914, 3.565518921744, 3.566166186334, 3.566758358887, 3.567362566356, 3.567966009572, 3.568550993809, 3.569164143932, 3.569733006766, 3.570370437634, 3.570947352048, 3.571555800369, 3.572129428124, 3.572710309127, 3.573370015175, 3.574068184328, 3.574664713671, 3.575333931415, 3.575979642369, 3.576582081200, 3.577173873700, 3.577774687145, 3.578415814237, 3.578963988217, 3.579616772722, 3.580212717782, 3.580903782622, 3.581594291645, 3.582179682550, 3.582737615754, 3.583339525558, 3.583972263335, 3.584574219296, 3.585210430881, 3.585814106655, 3.586425325669, 3.587008879097, 3.587628504669, 3.588279307062, 3.588979966912, 3.589533208365, 3.590136166530, 3.590792433120, 3.591444607179, 3.592055319581, 3.592743399127, 3.593388294600, 3.594049497716, 3.594646808663, 3.595231261431, 3.595756572530, 3.596318519937, 3.596927545696, 3.597535707533, 3.598113734736, 3.598725284973, 3.599354961147, 3.599978636059, 3.600589357233, 3.601147195663, 3.601776923596, 3.602379749752, 3.602967746232, 3.603558283066, 3.604238656493, 3.604893868429, 3.605557075199, 3.606193234751, 3.606858430265, 3.607436696374, 3.608089708706, 3.608687256239, 3.609289159731, 3.609957348733, 3.610557473272, 3.611172620111, 3.611820618244, 3.612439337306, 3.613087447791, 3.613700838030, 3.614350835238, 3.614974961704, 3.615616116415, 3.616202580443, 3.616874309779, 3.617484075726, 3.618026208822, 3.618578042033, 3.619184785302, 3.619772472551, 3.620440689082, 3.621086340016, 3.621685695156, 3.622302258853, 3.622877780127, 3.623519766350, 3.624148080209, 3.624788287330, 3.625429439593, 3.626001780100, 3.626626344539, 3.627231558762, 3.627883705257, 3.628522061713, 3.629222380368, 3.629855300951, 3.630511402229, 3.631090483896, 3.631688936635, 3.632316152043, 3.632938678992, 3.633586382541, 3.634214475579, 3.634819124941, 3.635392728230, 3.636040347711, 3.636698341294, 3.637219803848, 3.637785271900, 3.638466695942, 3.639098108394, 3.639830864948, 3.640388360266, 3.640950371968, 3.641631075520, 3.642244241153, 3.642835374936, 3.643406296312, 3.644004754951, 3.644569553342, 3.645196477552, 3.645782040214, 3.646406870527, 3.647017187667, 3.647643798830, 3.648255857940, 3.648841693436, 3.649434132404, 3.650017680533, 3.650650581857, 3.651290243903, 3.651870448131, 3.652468987015, 3.653033187183, 3.653615726174, 3.654281323718, 3.654936170178, 3.655517339585, 3.656150442095, 3.656756883235, 3.657346414870, 3.657954529518, 3.658516014175, 3.659131713287, 3.659797888312, 3.660453163519, 3.661043757732, 3.661665046168, 3.662277246834, 3.662888313343, 3.663492236238, 3.664137080443, 3.664831056413, 3.665449749686, 3.666073351713, 3.666754301038, 3.667327286250, 3.667939439651, 3.668501845359, 3.669109576151, 3.669791245795, 3.670443484652, 3.671074303791, 3.671701961869, 3.672340739885, 3.672953870264, 3.673547386910, 3.674182733973, 3.674767664568, 3.675376006598, 3.675985201964, 3.676599377946, 3.677206162304, 3.677861367446, 3.678471990181, 3.679087621251, 3.679745673199, 3.680385998739, 3.680943934299, 3.681571426670, 3.682210273212, 3.682816584568, 3.683419553358, 3.684044341078, 3.684661624464, 3.685288202910, 3.685940973459, 3.686582064008, 3.687190287148, 3.687873441207, 3.688491955620, 3.689104984011, 3.689706124896, 3.690346415366, 3.690985519305, 3.691629835268, 3.692240897828, 3.692839974964, 3.693472040771, 3.694111469716, 3.694738938795, 3.695373776402, 3.695990133303, 3.696540420920, 3.697162761834, 3.697809817656, 3.698442656839, 3.699074247521, 3.699654557372, 3.700259600247, 3.700920014907, 3.701474404678, 3.702112610308, 3.702732041450, 3.703372099249, 3.703964774072, 3.704542861927, 3.705148150502, 3.705805017008, 3.706465087792, 3.707046515640, 3.707599926289, 3.708202839171, 3.708839909742, 3.709471242303, 3.710050028630, 3.710631817903, 3.711279165829, 3.711896159941, 3.712455789009, 3.713065485671, 3.713729953103, 3.714251463508, 3.714906486032, 3.715521891160, 3.716172056378, 3.716764372813, 3.717377887122, 3.717944629316, 3.718532555273, 3.719219097153, 3.719870268057, 3.720508725273, 3.721152692980, 3.721774730942, 3.722432039902, 3.722982470767, 3.723604838717, 3.724235004413, 3.724829209123, 3.725408073257, 3.726052418640, 3.726690777591, 3.727367165897, 3.728014426645, 3.728634752846, 3.729241996742, 3.729866411663, 3.730468376851, 3.731040783862, 3.731646721309, 3.732281639462, 3.732875224890, 3.733429657448, 3.734034231278, 3.734625504032, 3.735264795212, 3.735857746846, 3.736491753501, 3.737119574462, 3.737703195274, 3.738294733193, 3.738894219536, 3.739504070039, 3.740090906208, 3.740738299455, 3.741319625143, 3.741892141718, 3.742563833890, 3.743195692321, 3.743847733876, 3.744466996349, 3.745135440698, 3.745754124895, 3.746354316862, 3.746928663543, 3.747588766820, 3.748184205232, 3.748812122090, 3.749433630960, 3.750029163506, 3.750640188478, 3.751127179606, 3.751705425536, 3.752338457910, 3.752967496393, 3.753555586293, 3.754188857319, 3.754744036595, 3.755391407661, 3.755948127332, 3.756612174936, 3.757222606207, 3.757781679752, 3.758413677701, 3.759084002657, 3.759720397501, 3.760302704153, 3.760878280061, 3.761539882774, 3.762159796236, 3.762725265924, 3.763346874930, 3.763954242889, 3.764557410488, 3.765196824149, 3.765842246479, 3.766432825118, 3.767059767346, 3.767679984790, 3.768306183928, 3.768925634213, 3.769510205770, 3.770190217047, 3.770820048745, 3.771471322643, 3.772105586085, 3.772679017451, 3.773204255390, 3.773830736591, 3.774512383946, 3.775052781223, 3.775640482680, 3.776239357670, 3.776802690470, 3.777486420367, 3.778085242639, 3.778737074279, 3.779410791901, 3.780004424501, 3.780609351745, 3.781236116994, 3.781816484185, 3.782431837144, 3.782992725166, 3.783615025670, 3.784227649093, 3.784886125932, 3.785479344904, 3.786115837132, 3.786692136310, 3.787375656403, 3.787974953398, 3.788628463135, 3.789159997798, 3.789769793021, 3.790388486294, 3.790978538545, 3.791590894349, 3.792220263874, 3.792837069639, 3.793578934762, 3.794205786189, 3.794809174121, 3.795399844620, 3.795869160118, 3.796458557049, 3.797035146898, 3.797620678920, 3.798239749758, 3.798862437354, 3.799409390585, 3.800047461015, 3.800708427222, 3.801389643447, 3.802052657085, 3.802675325311, 3.803243668942, 3.803884637033, 3.804449029595, 3.804983665277, 3.805557815272, 3.806168856187, 3.806791891506, 3.807426971674, 3.808043441402, 3.808649605951, 3.809279012272, 3.809872888702, 3.810473192499, 3.811034975199, 3.811625630186, 3.812214271195, 3.812846046215, 3.813419387261, 3.814018956453, 3.814715722069, 3.815422125068, 3.815996036307, 3.816593481857, 3.817188899576, 3.817799408762, 3.818390765719, 3.818939991468, 3.819570167721, 3.820169683321, 3.820729787261, 3.821359689463, 3.821984741570, 3.822558737873, 3.823095923716, 3.823737952034, 3.824334558221, 3.824940691581, 3.825617427531, 3.826216623003, 3.826813731588, 3.827397068790, 3.827969500437, 3.828569025580, 3.829125422999, 3.829729481428, 3.830366705391, 3.830963665761, 3.831526088052, 3.832095140379, 3.832732892688, 3.833336075459, 3.834031958326, 3.834639917021, 3.835227925231, 3.835897084181, 3.836570256699, 3.837184766773, 3.837821073913, 3.838479272770, 3.839123477793, 3.839738611011, 3.840366645057, 3.840995588619, 3.841649571808, 3.842244134092, 3.842887902845, 3.843581095971, 3.844129774913, 3.844703445646, 3.845244421921, 3.845901793811, 3.846499160258, 3.847091241372, 3.847684130783, 3.848290080479, 3.848949015025, 3.849630454010, 3.850297580245, 3.850891785200, 3.851492974535, 3.852057925586, 3.852654545605, 3.853261279503, 3.853847148203, 3.854458657980, 3.855123914789, 3.855712303310, 3.856264057613, 3.857000822131, 3.857635565985, 3.858249303656, 3.858885876196, 3.859520240447, 3.860199597263, 3.860832732515, 3.861460478517, 3.862003784968, 3.862690206654, 3.863336497181, 3.863885331550, 3.864453931665, 3.865086937603, 3.865711303920, 3.866336569154, 3.867017088850, 3.867583403495, 3.868211371812, 3.868837038583, 3.869511843134, 3.870213465480, 3.870845251050, 3.871397203442, 3.871975730635, 3.872603609223, 3.873206449486, 3.873865344474, 3.874479698348, 3.875052577589, 3.875603383976, 3.876151624417, 3.876739793907, 3.877302567756, 3.877964432334, 3.878581332670, 3.879235284343, 3.879824354645, 3.880503269481, 3.881107274626, 3.881791507606, 3.882453627798, 3.883083578374, 3.883741026910, 3.884436081804, 3.885078916160, 3.885702674946, 3.886334016570, 3.886906022340, 3.887508948694, 3.888092574233, 3.888724040899, 3.889383357545, 3.889939171556, 3.890566579616, 3.891269272450, 3.891966331045, 3.892630592117, 3.893234731212, 3.893812503646, 3.894445536481, 3.895059028706, 3.895710961207, 3.896326244869, 3.896973232256, 3.897607461623, 3.898149862062, 3.898692940764, 3.899257361915, 3.899843207936, 3.900495459920, 3.901110647177, 3.901650517511, 3.902336705431, 3.902985768822, 3.903621889028, 3.904164892082, 3.904743449596, 3.905343732970, 3.905944847206, 3.906490764292, 3.907061911716, 3.907732127505, 3.908385793926, 3.909001702247, 3.909660811081, 3.910282063911, 3.911024496168, 3.911633530446, 3.912307295189, 3.912974998092, 3.913590335113, 3.914178030667, 3.914805783742, 3.915323649271, 3.915920853672, 3.916558300335, 3.917236165481, 3.917911495333, 3.918580675865, 3.919211224038, 3.919828245286, 3.920464224872, 3.921072166319, 3.921622942848, 3.922210724637, 3.922777489589, 3.923432369555, 3.924011669414, 3.924628251553, 3.925209149851, 3.925823774565, 3.926398943326, 3.926949180089, 3.927536868714, 3.928125353680, 3.928681468719, 3.929293649826, 3.929858654270, 3.930476198648, 3.931150212561, 3.931754751616, 3.932453058010, 3.933029635023, 3.933677802504, 3.934278406555, 3.934954613072, 3.935571963150, 3.936126455049, 3.936685409582, 3.937248843083, 3.937809244976, 3.938423131455, 3.939083179385, 3.939574150911, 3.940205676071, 3.940678976207, 3.941224858103, 3.941854991193, 3.942444195034, 3.943003727869, 3.943567796253, 3.944170787010, 3.944793738374, 3.945459716929, 3.946057670645, 3.946556595078, 3.947167598295, 3.947679350121, 3.948288080643, 3.948905387344, 3.949527439162, 3.950150383240, 3.950805243729, 3.951387309169, 3.952051816754, 3.952592722199, 3.953251278599, 3.953824918977, 3.954528369267, 3.955103700099, 3.955762155675, 3.956405898122, 3.956995523608, 3.957648976988, 3.958256073199, 3.958891673394, 3.959524248674, 3.960300411172, 3.960859632099, 3.961439443396, 3.962016050553, 3.962633271468, 3.963147628221, 3.963770455914, 3.964378174024, 3.964958698744, 3.965564071132, 3.966182341172, 3.966753215354, 3.967332897371, 3.967929489041, 3.968555177732, 3.969262685551, 3.969902455637, 3.970563460770, 3.971209215705, 3.971843721129, 3.972552534072, 3.973144096517, 3.973744642105, 3.974370582987, 3.974956429463, 3.975633389639, 3.976262060642, 3.976833996214, 3.977398440481, 3.977938851536, 3.978463404245, 3.979112744007, 3.979713316323, 3.980269067347, 3.980775668863, 3.981320301376, 3.981944766325, 3.982545098612, 3.983183862247, 3.983764992529, 3.984380416367, 3.985017691415, 3.985655902961, 3.986244560277, 3.986775034917, 3.987360993733, 3.988023784287, 3.988662200948, 3.989428694097, 3.989963074734, 3.990578851417, 3.991135925515, 3.991706497298, 3.992329019562, 3.992922525202, 3.993589585237, 3.994266243312, 3.994905324927, 3.995536751197, 3.996151877966, 3.996742013534, 3.997276843087, 3.997903085158, 3.998486950996, 3.999110607628, 3.999704780116, 4.000364960670, 4.000917329480, 4.001492190645, 4.002076541176, 4.002652939821, 4.003348255443, 4.003970171589, 4.004553479337, 4.005106810298, 4.005836879778, 4.006435924971, 4.007018142456, 4.007609980853, 4.008202626883, 4.008937923593, 4.009541264510, 4.010127662770, 4.010719305324, 4.011320670164, 4.011945188565, 4.012467797270, 4.012986560820, 4.013537306790, 4.014120156202, 4.014762195313, 4.015355690432, 4.015940986989, 4.016504516995, 4.017231422749, 4.017859979905, 4.018426007086, 4.019019996322, 4.019605711992, 4.020137626474, 4.020825078769, 4.021372175545, 4.021947369876, 4.022564496064, 4.023168757439, 4.023819736264, 4.024384402033, 4.024963602398, 4.025483701831, 4.026041313893, 4.026682761383, 4.027325157684, 4.028005561846, 4.028682394185, 4.029206968968, 4.029815898455, 4.030355816858, 4.031026995559, 4.031638484083, 4.032180674518, 4.032718859892, 4.033271779013, 4.033877047969, 4.034389135116, 4.035024209672, 4.035603642794, 4.036282987310, 4.036859375790, 4.037474403264, 4.038057116859, 4.038635866322, 4.039115579209, 4.039719535563, 4.040233797887, 4.040891798282, 4.041455228067, 4.042052880021, 4.042660937880, 4.043255454717, 4.043740300408, 4.044268966679, 4.044865690384, 4.045415015539, 4.046037459449, 4.046781702108, 4.047435173345, 4.048021713963, 4.048681913159, 4.049275005844, 4.049888395604, 4.050492896024, 4.051098239030, 4.051797385411, 4.052482956465, 4.053110712762, 4.053857353839, 4.054398487923, 4.055073390014, 4.055690092969, 4.056381841536, 4.056990502654, 4.057585141594, 4.058200458376, 4.058866379031, 4.059478534655, 4.060056643824, 4.060710435076, 4.061265182213, 4.061930813777, 4.062657666636, 4.063270169359, 4.063863412978, 4.064447392571, 4.064996848546, 4.065617712815, 4.066168652381, 4.066902638603, 4.067516088374, 4.068135486508, 4.068689631406, 4.069315801182, 4.069871455134, 4.070535092598, 4.071199745707, 4.071819299710, 4.072506450325, 4.073014820562, 4.073626680456, 4.074177577179, 4.074770442610, 4.075467449352, 4.075984472808, 4.076548730360, 4.077098161910, 4.077882045013, 4.078479992625, 4.079042294045, 4.079641842358, 4.080273566316, 4.080817270957, 4.081408799570, 4.081938193741, 4.082447230016, 4.083030468412, 4.083730336687, 4.084262570071, 4.084795456513, 4.085371286330, 4.086032582819, 4.086509340568, 4.087124601352, 4.087703529736, 4.088261943761, 4.088842391260, 4.089391601345, 4.089962876447, 4.090465365103, 4.091059478955, 4.091557876522, 4.092207186643, 4.092835956490, 4.093465637988, 4.094187935046, 4.094846595161, 4.095392636185, 4.095933948085, 4.096595263592, 4.097224990701, 4.097746835586, 4.098307429836, 4.098923284274, 4.099545474959, 4.100223257408, 4.100880184597, 4.101450326387, 4.102197027576, 4.102834935862, 4.103402140853, 4.103986640879, 4.104638237340, 4.105241005628, 4.105900030415, 4.106454607091, 4.107148825289, 4.107833020927, 4.108390073043, 4.108925516125, 4.109595749650, 4.110177456557, 4.110659072667, 4.111152442644, 4.111826124417, 4.112495225765, 4.113114641014, 4.113718011941, 4.114418295779, 4.114995183362, 4.115606841901, 4.116185311484, 4.116764552602, 4.117281979431, 4.117868383108, 4.118524042555, 4.119214980664, 4.119752502212, 4.120388092934, 4.121001661257, 4.121616097652, 4.122214141098, 4.122795722600, 4.123343466415, 4.123932340379, 4.124510443976, 4.125106696282, 4.125668964902, 4.126313280095, 4.126882930567, 4.127558177796, 4.128286992051, 4.128859238981, 4.129367886700, 4.129824422721, 4.130592250618, 4.131173419441, 4.131714193583, 4.132343985584, 4.133098583343, 4.133665393340, 4.134156044805, 4.134706470233, 4.135251664370, 4.135963817274, 4.136605754873, 4.137266514323, 4.137922314258, 4.138644839807, 4.139380542019, 4.140021558358, 4.140615492098, 4.141204227584, 4.141703475466, 4.142257530126, 4.142872635512, 4.143494656652, 4.144141778804, 4.144705015090, 4.145256846740, 4.145797229304, 4.146393036146, 4.147026216218, 4.147648117730, 4.148270911071, 4.148839531489, 4.149415023648, 4.150009683153, 4.150506879109, 4.151139939158, 4.151780083074, 4.152408834385, 4.153044674980, 4.153761883650, 4.154412104009, 4.154957801438, 4.155640888873, 4.156194356101, 4.156698687291, 4.157309649742, 4.157877742652, 4.158552928543, 4.159141445958, 4.159561408470, 4.160082224321, 4.160584807316, 4.161201265165, 4.161749262945, 4.162241159864, 4.162771518824, 4.163353132261, 4.163967199850, 4.164518701080, 4.165185240038, 4.165795543594, 4.166336633332, 4.167025109470, 4.167625228363, 4.168123829677, 4.168699851092, 4.169385672133, 4.169969793918, 4.170496819193, 4.171204810236, 4.171681745338, 4.172243137021, 4.172811720690, 4.173342208124, 4.173996505232, 4.174502479443, 4.175139027360, 4.175769999990, 4.176284563113, 4.176780161909, 4.177505014721, 4.178113256658, 4.178656816560, 4.179187935366, 4.179851105811, 4.180350793707, 4.181035510150, 4.181642122042, 4.182216546425, 4.182937301457, 4.183400803212, 4.184043902989, 4.184575012072, 4.185313003378, 4.185812357031, 4.186412340293, 4.187113369994, 4.187735228418, 4.188357978547, 4.188847430260, 4.189438188461, 4.190029751348, 4.190561500739, 4.191087158665, 4.191809290122, 4.192431149690, 4.192992940226, 4.193521551149, 4.194084754884, 4.194560314387, 4.195152093638, 4.195874203498, 4.196433641619, 4.197034816987, 4.197650517410, 4.198376796905, 4.199021873760, 4.199619767974, 4.200246033588, 4.200749076116, 4.201238896239, 4.201819116074, 4.202448563563, 4.202877941016, 4.203446480123, 4.204085240487, 4.204648406093, 4.205240168661, 4.205895528724, 4.206523929210, 4.206971345614, 4.207713394460, 4.208456713363, 4.209004497549, 4.209496687702, 4.210024654022, 4.210659062006, 4.211259077234, 4.211711401101, 4.212348279652, 4.212943543455, 4.213532523260, 4.214164968641, 4.214748487450, 4.215339921841, 4.215989288922, 4.216611021582, 4.217204999611, 4.217821304876, 4.218330774946, 4.218833655152, 4.219437881043, 4.220057365407, 4.220684953486, 4.221248390736, 4.221928377571, 4.222544177278, 4.223117292715, 4.223771131971, 4.224302192234, 4.224943238921, 4.225672851487, 4.226403691862, 4.226967278245, 4.227568266384, 4.228126023370, 4.228618324990, 4.229177433050, 4.229921573882, 4.230474979829, 4.231110421109, 4.231709769565, 4.232258046417, 4.232858981846, 4.233505358548, 4.234160145014, 4.234718981956, 4.235390536879, 4.235936006307, 4.236609447230, 4.237193941928, 4.237734174590, 4.238350259178, 4.238952160744, 4.239502124773, 4.240150922390, 4.240634378344, 4.241262162913, 4.241951501342, 4.242725469174, 4.243280235376, 4.243965170564, 4.244536776219, 4.244941165035, 4.245544633064, 4.246110668540, 4.246754090586, 4.247383113879, 4.247851637511, 4.248389911133, 4.248967374163, 4.249522461981, 4.249931522115, 4.250565167783, 4.251222972928, 4.251796462754, 4.252487205139, 4.253062367935, 4.253708390582, 4.254433390738, 4.254964265433, 4.255574010108, 4.256168944872, 4.256654890743, 4.257361264979, 4.257880003925, 4.258430859506, 4.258974530382, 4.259653103320, 4.260214120604, 4.260783780606, 4.261512768578, 4.262274758543, 4.263022171421, 4.263635373718, 4.264409083730, 4.264944307980, 4.265560232183, 4.266104892529, 4.266722466067, 4.267180197228, 4.267815395284, 4.268419292300, 4.269032099082, 4.269637691508, 4.270171312273, 4.270770395072, 4.271386530711, 4.271987293339, 4.272597023527, 4.273183170966, 4.273745638910, 4.274382350668, 4.274930016023, 4.275453804777, 4.276134035557, 4.276757833302, 4.277300280256, 4.277859874543, 4.278345989778, 4.278931698203, 4.279427279636, 4.280039276608, 4.280743311987, 4.281257564489, 4.281797354818, 4.282337816895, 4.282887282824, 4.283570922440, 4.284163724835, 4.284707141650, 4.285435549905, 4.285930224021, 4.286593467832, 4.287148347901, 4.287864012225, 4.288530218126, 4.289188995203, 4.289730277587, 4.290297656316, 4.290831838724, 4.291460131372, 4.292097843275, 4.292804671059, 4.293341947927, 4.293862802517, 4.294418499821, 4.295017739325, 4.295557762732, 4.296107046367, 4.296734422320, 4.297414386073, 4.297897031310, 4.298423380945, 4.298941724960, 4.299555898972, 4.300101598504, 4.300561210722, 4.301168992139, 4.301847237738, 4.302369686853, 4.302988731301, 4.303731003254, 4.304334487744, 4.304851176641, 4.305579086304, 4.306149989106, 4.306721643377, 4.307320488139, 4.307796630267, 4.308432299562, 4.309104294741, 4.309582397660, 4.310078764655, 4.310708903795, 4.311268807945, 4.311882865675, 4.312417535574, 4.312926082404, 4.313426287946, 4.313980759787, 4.314661402094, 4.315262328024, 4.315837124345, 4.316241733997, 4.316835844188, 4.317403708539, 4.318035541005, 4.318740667398, 4.319320088932, 4.319891213023, 4.320490340157, 4.321199466746, 4.321773069327, 4.322420419643, 4.323096151311, 4.323763783237, 4.324294963159, 4.324927731454, 4.325598187155, 4.326232858493, 4.326960651883, 4.327523454579, 4.328133211227, 4.328808637659, 4.329299674768, 4.329865518764, 4.330404218976, 4.331064573282, 4.331660683930, 4.332313618197, 4.332939490353, 4.333388431581, 4.333987742704, 4.334653572575, 4.335273429862, 4.335884761297, 4.336468680407, 4.337025075346, 4.337525496250, 4.338102167030, 4.338745917006, 4.339352672259, 4.339979277732, 4.340568731804, 4.341092305106, 4.341730966568, 4.342408783621, 4.342972845964, 4.343643024082, 4.344266259710, 4.344900000029, 4.345342245350, 4.345919764693, 4.346584862742, 4.347357268151, 4.347908446304, 4.348508768738, 4.348954706101, 4.349488493815, 4.350120180855, 4.350830711100, 4.351483866259, 4.351991470567, 4.352617029885, 4.353165134416, 4.353821812536, 4.354351799088, 4.354921765227, 4.355531867695, 4.356113246035, 4.356784276042, 4.357416782226, 4.357901668169, 4.358535805008, 4.359131150412, 4.359747199245, 4.360493589423, 4.361031785752, 4.361570649864, 4.362160174151, 4.362630368515, 4.363201288130, 4.363762923453, 4.364275045469, 4.364848132690, 4.365492501577, 4.366188284473, 4.366763903617, 4.367360524572, 4.367978232945, 4.368647564196, 4.369206129323, 4.369857001654, 4.370468081651, 4.371080022687, 4.371662166424, 4.372163230268, 4.372798059047, 4.373495391541, 4.374029402321, 4.374574359067, 4.375099397796, 4.375790120652, 4.376492277889, 4.377029988690, 4.377599446675, 4.378159278268, 4.378823718225, 4.379468365280, 4.380041031634, 4.380551862126, 4.381251318541, 4.381889118263, 4.382433556932, 4.383073096522, 4.383661045143, 4.384197191770, 4.384786665199, 4.385408584321, 4.385968017799, 4.386528172837, 4.387173775154, 4.387777911515, 4.388372268532, 4.388914266585, 4.389595379419, 4.390224228005, 4.390928766450, 4.391591648358, 4.392105543035, 4.392738039461, 4.393178113983, 4.393908978233, 4.394641074511, 4.395212535982, 4.395719933548, 4.396238739392, 4.396801479324, 4.397343263928, 4.397853158460, 4.398483216350, 4.399037989023, 4.399746080969, 4.400313386945, 4.400859573200, 4.401384559039, 4.402074567440, 4.402611999595, 4.403128121396, 4.403688863219, 4.404404585180, 4.404922842992, 4.405419627492, 4.406038642868, 4.406647464421, 4.407312608234, 4.407834350022, 4.408523565327, 4.409135882851, 4.409760221691, 4.410351941664, 4.411011598476, 4.411661052766, 4.412244148982, 4.412850501746, 4.413480207490, 4.414133366694, 4.414708509105, 4.415329615536, 4.416098759001, 4.416687848017, 4.417164234709, 4.417743407668, 4.418391633968, 4.418949655563, 4.419462756734, 4.419999310411, 4.420547965179, 4.421062959289, 4.421647358490, 4.422152178698, 4.422588632515, 4.423301684364, 4.424027438252, 4.424523491398, 4.425297542393, 4.426015059488, 4.426629362091, 4.427232919744, 4.427767536225, 4.428314455383, 4.429176960475, 4.429643893626, 4.430403732354, 4.431071148279, 4.431739591455, 4.432420819416, 4.433008943771, 4.433586079380, 4.434270212169, 4.434825387218, 4.435215615496, 4.435831233118, 4.436483318156, 4.437065092499, 4.437754732067, 4.438230984426, 4.438922479341, 4.439495584805, 4.440057484341, 4.440775860953, 4.441423417458, 4.442071940946, 4.442781621746, 4.443287525790, 4.443806086427, 4.444409843520, 4.445026541667, 4.445535069095, 4.446080581489, 4.446663217623, 4.447234473898, 4.447916102149, 4.448550002027, 4.449074888264, 4.449612638641, 4.450249021460, 4.450910869024, 4.451561442414, 4.452003897777, 4.452693058902, 4.453259974791, 4.453889378251, 4.454618651506, 4.455237639783, 4.455795484543, 4.456440998826, 4.456963077092, 4.457460878656, 4.458034056743, 4.458620477594, 4.459320227169, 4.459996055150, 4.460522427195, 4.461036882675, 4.461639947404, 4.462193493355, 4.462760350619, 4.463391061006, 4.464060615572, 4.464579281944, 4.465225317985, 4.465846925815, 4.466418574941, 4.467092816990, 4.467691606725, 4.468380598644, 4.468827752017, 4.469467342251, 4.470005326931, 4.470556812061, 4.471096148862, 4.471790567523, 4.472460318035, 4.473118193518, 4.473789996158, 4.474372203788, 4.475007052117, 4.475668799711, 4.476227528309, 4.476656808059, 4.477125597686, 4.477699250833, 4.478208350361, 4.478770357834, 4.479228343290, 4.479844113055, 4.480460757137, 4.481025688856, 4.481657179535, 4.482052328005, 4.482513789420, 4.483028567531, 4.483610076312, 4.484086436048, 4.484629594016, 4.485107073949, 4.485571794383, 4.486130116108, 4.486729115534, 4.487248917490, 4.487876173277, 4.488504336326, 4.489093228438, 4.489843884231, 4.490461466448, 4.491053020178, 4.491780120925, 4.492319499989, 4.492751486081, 4.493467909592, 4.494022935890, 4.494524422903, 4.494958608875, 4.495542730660, 4.496073195806, 4.496767860088, 4.497490945390, 4.498064812651, 4.498721590821, 4.499406792569, 4.500024399818, 4.500629132888, 4.501220936713, 4.501703233556, 4.502227477276, 4.502738533893, 4.503291704879, 4.503901007867, 4.504622197056, 4.505288974795, 4.505831486510, 4.506541949005, 4.507058111018, 4.507421186812, 4.508008335818, 4.508470224703, 4.509142941215, 4.509774560727, 4.510491508422, 4.511111003995, 4.511632628308, 4.512112511762, 4.512578788641, 4.513187113665, 4.513909721815, 4.514349542467, 4.515102527184, 4.515614734076, 4.516284359096, 4.516840790283, 4.517469415868, 4.518127589689, 4.518586038746, 4.519260264406, 4.519964296771, 4.520525465897, 4.521173871044, 4.521866571899, 4.522459130063, 4.522994572587, 4.523588672866, 4.524183586969, 4.524633941721, 4.525361305866, 4.525973234585, 4.526498432047, 4.527141203775, 4.527726368389, 4.528312322513, 4.528884389970, 4.529471909651, 4.530045507607, 4.530472520813, 4.531106452141, 4.531756085348, 4.532377097582, 4.532939731996, 4.533488261514, 4.534007779216, 4.534602275708, 4.535227374110, 4.535942875735, 4.536674501563, 4.537377425099, 4.537871647885, 4.538366433732, 4.538936886296, 4.539643484138, 4.540110168029, 4.540667835108, 4.541211118289, 4.541891179383, 4.542557159272, 4.543133146733, 4.543633966871, 4.544196180123, 4.544713450989, 4.545338039334, 4.545856672825, 4.546406490079, 4.547048824536, 4.547676783022, 4.548321000427, 4.549058420928, 4.549550731137, 4.550151489734, 4.550660473825, 4.551324591834, 4.551927810935, 4.552547368773, 4.553214380699, 4.553851324304, 4.554551485734, 4.555205989724, 4.555814627903, 4.556408481551, 4.557050130422, 4.557645677061, 4.558085024070, 4.558650552666, 4.559311268128, 4.559862633233, 4.560525196683, 4.560999075920, 4.561647548236, 4.562011749780, 4.562709335343, 4.563408043205, 4.563964637392, 4.564633493387, 4.565239537711, 4.565670662662, 4.566294155243, 4.567014683795, 4.567656158719, 4.568330728686, 4.568829298517, 4.569586288107, 4.570150862142, 4.570732333567, 4.571346954738, 4.571962446965, 4.572481433636, 4.573163546196, 4.573732792860, 4.574253900704, 4.574922485226, 4.575608445590, 4.576279120265, 4.576639421915, 4.577262466054, 4.578165826495, 4.578741672122, 4.579367741666, 4.579994715041, 4.580374639112, 4.580969970915, 4.581383876443, 4.581997181476, 4.582594743156, 4.583060082459, 4.583775682901, 4.584325665434, 4.584909742329, 4.585310706680, 4.585879373956, 4.586314740145, 4.586717004369, 4.587388273844, 4.588060582476, 4.588582338461, 4.589071000897, 4.589509579997, 4.590134476401, 4.590844909597, 4.591302231104, 4.591810931555, 4.592405169489, 4.592983208815, 4.593562018527, 4.594141600682, 4.594721957346, 4.595251782835, 4.595885002494, 4.596450545968, 4.596982485778, 4.597583847093, 4.598203260112, 4.598633928402, 4.599306624877, 4.599928501652, 4.600360885036, 4.601122926913, 4.601764771444, 4.602216365490, 4.602790219621, 4.603382257417, 4.604027452205, 4.604708561616, 4.605303222108, 4.605951279167, 4.606652970934, 4.607461320530, 4.608147840941, 4.608729591669, 4.609524146062, 4.610249341162, 4.610763016201, 4.611277299525, 4.611969883540, 4.612716980519, 4.613233584283, 4.613857891069, 4.614608245785, 4.615091302170, 4.615664510066, 4.616346178537, 4.616921046188, 4.617622696532, 4.617982957425, 4.618541953078, 4.619029407184, 4.619481241930, 4.619987855556, 4.620640085512, 4.621220668939, 4.621656616425, 4.622256759646, 4.622803064991, 4.623441290994, 4.623861207920, 4.624647360767, 4.625324958823, 4.626113768211, 4.626775270406, 4.627364119418, 4.627916891426, 4.628581147890, 4.629116983222, 4.629597950027, 4.630172107475, 4.630654244865, 4.631229801410, 4.631787518813, 4.632290077575, 4.632961061844, 4.633521008391, 4.634119081534, 4.634624347174, 4.635242693429, 4.635843143877, 4.636369219888, 4.637084201000, 4.637743779217, 4.638272163982, 4.638801192391, 4.639330866012, 4.639918044894, 4.640524998892, 4.640961771621, 4.641570187070, 4.642179456062, 4.642789580994, 4.643496108106, 4.643955011313, 4.644452704268, 4.645200314737, 4.645641817265, 4.646256829246, 4.646795680093, 4.647412329414, 4.648087793672, 4.648590248075, 4.649132003713, 4.649868326127, 4.650372845934, 4.650877952524, 4.651539364180, 4.652182287291, 4.652709024144, 4.653334133258, 4.653823126578, 4.654489042303, 4.655254146275, 4.655882932170, 4.656433867614, 4.657143241698, 4.657972311913, 4.658605048476, 4.659238708237, 4.659714560530, 4.660230656231, 4.660687714229, 4.661264691200, 4.661902246293, 4.662740460250, 4.663200168652, 4.663880630191, 4.664341547711, 4.664883249719, 4.665224668322, 4.665968684980, 4.666552725503, 4.667056839874, 4.667460553110, 4.668127501511, 4.668653698347, 4.669342765289, 4.669870437512, 4.670622462777, 4.671049867930, 4.671640786452, 4.672252928748, 4.672825041062, 4.673356963897, 4.673889539029, 4.674176580997, 4.674710162792, 4.675511766692, 4.675944013185, 4.676438536537, 4.677016193052, 4.677718667628, 4.678235918924, 4.678878167291, 4.679459080075, 4.680144726332, 4.680852283745, 4.681435844899, 4.681894907960, 4.682500778893, 4.683065626440, 4.683652171267, 4.684260500392, 4.684764590382, 4.685269266155, 4.685816660863, 4.686449128666, 4.687061390987, 4.687420705094, 4.688161406824, 4.688882157337, 4.689349163522, 4.689901728144, 4.690369832574, 4.691008971000, 4.691520959838, 4.692119044012, 4.692824987960, 4.693339123449, 4.693939718472, 4.694584136162, 4.695121882302, 4.695789613565, 4.696156225111, 4.696501554208, 4.697084920431, 4.697798989395, 4.698123954727, 4.698818027862, 4.699295847404, 4.700100640454, 4.700776076633, 4.701452564921, 4.701780276647, 4.702458332182, 4.703005923330, 4.703663945398, 4.704103181321, 4.704740863518, 4.705401521546, 4.706306049254, 4.707035454617, 4.707677460085, 4.708276043656, 4.708831024285, 4.709475691563, 4.710232729047, 4.710812531779, 4.711303739410, 4.711929717357, 4.712444589358, 4.713184386563, 4.713902971186, 4.714645258970, 4.715208444050, 4.715727220142, 4.716472635138, 4.717106113345, 4.717604495257, 4.718216928494, 4.718784767389, 4.719444391994, 4.720127818197, 4.720652507589, 4.721200686244, 4.721635152262, 4.722253299973, 4.722964111828, 4.723699075429, 4.724435284930, 4.725011318229, 4.725518860311, 4.726235041495, 4.726697721432, 4.727346302570, 4.728111947340, 4.728599879559, 4.729204746623, 4.729903718580, 4.730370325642, 4.730860803208, 4.731305046438, 4.732077708469, 4.732640509172, 4.733227536487, 4.734050710049, 4.734686797446, 4.735182176990, 4.735725384668, 4.736174634479, 4.736648031606, 4.737335375324, 4.737952540739, 4.738594372926, 4.739332463010, 4.739952475005, 4.740501685969, 4.741099443216, 4.741698024343, 4.742081549686, 4.742753535093, 4.743306302002, 4.743883853346, 4.744630995734, 4.745379425690, 4.745790403270, 4.746201770132, 4.746758946227, 4.747438213794, 4.748191501376, 4.748824299993, 4.749262932096, 4.749824051916, 4.750459235498, 4.750948471195, 4.751634328387, 4.752149433026, 4.752566870437, 4.753033893441, 4.753550661733, 4.754191323213, 4.754684783647, 4.755252957164, 4.755648647170, 4.756168538018, 4.756639452326, 4.757235022248, 4.757582815582, 4.758179681482, 4.758727530098, 4.759301020814, 4.759950227887, 4.760525336535, 4.761076154008, 4.761452112319, 4.762129658523, 4.762883729465, 4.763588712256, 4.764092972959, 4.764724123313, 4.765179120193, 4.765761199653, 4.766293345714, 4.766876921479, 4.767461282469, 4.768122812654, 4.768785352037, 4.769346752821, 4.770011163455, 4.770522941127, 4.771086594005, 4.771522645555, 4.772216099146, 4.772936407657, 4.773400094793, 4.773915884024, 4.774561483195, 4.775208043507, 4.775777813915, 4.776659841122, 4.777413576610, 4.777908113546, 4.778585762158, 4.779212223834, 4.779813432097, 4.780336899263, 4.780834778170, 4.781569536381, 4.782095125797, 4.782832021229, 4.783174576734, 4.783886902685, 4.784547507117, 4.785182634131, 4.786004383632, 4.786880861189, 4.787306457579, 4.787839040247, 4.788372276831, 4.788906168941, 4.789494209319, 4.790324165429, 4.790726327215, 4.791263122886, 4.791773694064, 4.792338708803, 4.792823592905, 4.793336003009, 4.793822002480, 4.794416741359, 4.794876870796, 4.795418824422, 4.796233025039, 4.796722279108, 4.797375476421, 4.797920558993, 4.798439021883, 4.798903434961, 4.799505178261, 4.800025537470, 4.800519085138, 4.801205499824, 4.801865480168, 4.802333573806, 4.802967681082, 4.803464585578, 4.804211009964, 4.804792450497, 4.805125052070, 4.805652192720, 4.806207769720, 4.806958927847, 4.807571944669, 4.808046232847, 4.808493095376, 4.808968391151, 4.809388202186, 4.810032701649, 4.810481613874, 4.811212087363, 4.811915626285, 4.812592097577, 4.813128385635, 4.813750179221, 4.814344541080, 4.814939717479, 4.815422125068, 4.816018781085, 4.816730156317, 4.817471224721, 4.817985013849, 4.818671012910, 4.819444059296, 4.819988880942, 4.820390766274, 4.821023052707, 4.821569860052, 4.821973211730, 4.822463500070, 4.823272242356, 4.823879788944, 4.824372235633, 4.824836225508, 4.825300711629, 4.825882018746, 4.826231176863, 4.826755440938, 4.827046972454, 4.827630623236, 4.828068876644, 4.828449054323, 4.828888134820, 4.829532923614, 4.829826326194, 4.830443117568, 4.830913642513, 4.831296319751, 4.831649859814, 4.832357804518, 4.832889521403, 4.833599491425, 4.834162375310, 4.834844738608, 4.835260615709, 4.836093566542, 4.836600067213, 4.837137006678, 4.837525209562, 4.838153040472, 4.838901643549, 4.839351425564, 4.839831707042, 4.840492962523, 4.841215482277, 4.841878849663, 4.842361935899, 4.843027058104, 4.843753809603, 4.844208647642, 4.844694333840, 4.844998163769, 4.845606462043, 4.846032778355, 4.846673038509, 4.847436485743, 4.848231897105, 4.848844744849, 4.849550590539, 4.850073040886, 4.850626909509, 4.851273983602, 4.851922023238, 4.852601960652, 4.853159065965, 4.853934010823, 4.854492828590, 4.855301281357, 4.855799539816, 4.856298370575, 4.856891477238, 4.857454115914, 4.857892228269, 4.858612945133, 4.859177819891, 4.859995049381, 4.860435733824, 4.860876865890, 4.861381566101, 4.862202956942, 4.862899195421, 4.863437963410, 4.864072664995, 4.864612891618, 4.865408565289, 4.865950457260, 4.866844461802, 4.867484165087, 4.867932518696, 4.868381335651, 4.868862726221, 4.869312506018, 4.869955857712, 4.870890413779, 4.871600731282, 4.872150409812, 4.872700784942, 4.873089705400, 4.873641272931, 4.874096031592, 4.874909286917, 4.875821944525, 4.876344325562, 4.877096350827, 4.877948037367, 4.878571481632, 4.879097182386, 4.879722279661, 4.880183454063, 4.880909147578, 4.881437687364, 4.881735273911, 4.882364182611, 4.882927663414, 4.883392256012, 4.884189858590, 4.884622505721, 4.885022255215, 4.885689323132, 4.886357417228, 4.887026540667, 4.887361489382, 4.887897945229, 4.888300722426, 4.889141043268, 4.889780776311, 4.890421453096, 4.890894136245, 4.891536458796, 4.892044227452, 4.892620417196, 4.893163412633, 4.893639091193, 4.894319537054, 4.895035154472, 4.895512888688, 4.896025330614, 4.896743766645, 4.897257663855, 4.897634908662, 4.898149862062, 4.898699819532, 4.899457149988, 4.899870796640, 4.900699273767, 4.901425489974, 4.901841016539, 4.902083590763, 4.902881575655, 4.903472334394, 4.903994260285, 4.904621400440, 4.905214532180, 4.905493933546, 4.906088258951, 4.906788508186, 4.907209200323, 4.907735639212, 4.908016667763, 4.908860846174, 4.909424544778, 4.909918381961, 4.910660192229, 4.911155437273, 4.911863911299, 4.912538033828, 4.912964335194, 4.913497801030, 4.914031922954, 4.914816491947, 4.915316502097, 4.915960219332, 4.916461548770, 4.917035206222, 4.917358221843, 4.918040933700, 4.918688704840, 4.919337443606, 4.920059402848, 4.920529320765, 4.920963542913, 4.921470684746, 4.921942133021, 4.922305134113, 4.923068425444, 4.923578032639, 4.924015315256, 4.924599044486, 4.925220117668, 4.925842080303, 4.926281649654, 4.926868435059, 4.927529517696, 4.928338876558, 4.929002203007, 4.929592678260, 4.930110003493, 4.930553916120, 4.931220636990, 4.931777020654, 4.932445623307, 4.932929143955, 4.933152489026, 4.933711354734, 4.934308271907, 4.934868627860, 4.935317433772, 4.935691793317, 4.936066475837, 4.936779264418, 4.937342819774, 4.937718930027, 4.938170692705, 4.938886946082, 4.939302159646, 4.939906816176, 4.940550187492, 4.941194513324, 4.941498056570, 4.942029768289, 4.942676294631, 4.943323784879, 4.943972241914, 4.944774616146, 4.945348649558, 4.946000139307, 4.946767851160, 4.947306058075, 4.947767909748, 4.948384476995, 4.948963304859, 4.949620243739, 4.950200722081, 4.950820755354, 4.951441675102, 4.952180172183, 4.952764084540, 4.953036845877, 4.953621911952, 4.954207767270, 4.954794413965, 4.955303482890, 4.955773922887, 4.956480539754, 4.957030926607, 4.957542625426, 4.957857816935, 4.958212681028, 4.958765268829, 4.959239475771, 4.959951758453, 4.960347976418, 4.961022377667, 4.961658066346, 4.962135444226, 4.962812629063, 4.963171566623, 4.963610671334, 4.964290162172, 4.964770443650, 4.965331444166, 4.965973476225, 4.966495827055, 4.966898063336, 4.967381239149, 4.967703655261, 4.968147368604, 4.968470354197, 4.968753163767, 4.969197951228, 4.969764703988, 4.970494474573, 4.971144190610, 4.971876284871, 4.972405788174, 4.972976745411, 4.973711937976, 4.974080001498, 4.974530280939, 4.975104039893, 4.975391203873, 4.975760693931, 4.976212720210, 4.976829878879, 4.977530387232, 4.977984260182, 4.978603943258, 4.979224511806, 4.979804510247, 4.980219269596, 4.980925270822, 4.981424316533, 4.981798977504, 4.982299028776, 4.982883150561, 4.983342655178, 4.983970036924, 4.984346901465, 4.984766023753, 4.985647499349, 4.986109939672, 4.986614982250, 4.987162775295, 4.987837932029, 4.988387270781, 4.988937305270, 4.989572827283, 4.990081915334, 4.990761629032, 4.991442408216, 4.992038966664, 4.992550955502, 4.992978074421, 4.993448390913, 4.993790759462, 4.994433428867, 4.995077050704, 4.995506662453, 4.995764633641, 4.996237979172, 4.996884282900, 4.997402019280, 4.998223032923, 4.998828985059, 4.999522538626, 5.000304112589, 5.000695427662, 5.001130636118, 5.001697059901, 5.002395212540, 5.002963289117, 5.003488327846, 5.004145520125, 5.004452550525, 5.004759798137, 5.005594853311, 5.006299305180, 5.006651960077, 5.007049039430, 5.007358130122, 5.007711646203, 5.008153946355, 5.008640997362, 5.009039900318, 5.009661145212, 5.010194349063, 5.010550182333, 5.011307297450, 5.012065734768, 5.012512493867, 5.012914970376, 5.013362604390, 5.013855535290, 5.014304140310, 5.014843078723, 5.015337693810, 5.016058141016, 5.016418813294, 5.016960383954, 5.017367005605, 5.017864505196, 5.018589159834, 5.019405844224, 5.019905686215, 5.020406104151, 5.020861530453, 5.021545566635, 5.022093572363, 5.022505030927, 5.023100049017, 5.023971159909, 5.024476287040, 5.025074013910, 5.025902996206, 5.026548855975, 5.027334407734, 5.027889771597, 5.028399479770, 5.029188389128, 5.029885677715, 5.030537492025, 5.031190286087, 5.031657169041, 5.031844062850, 5.032498825277, 5.032873418924, 5.033154576345, 5.033576654056, 5.034140063193, 5.034798298974, 5.035504660445, 5.035881856849, 5.036306594762, 5.036779013477, 5.037393927076, 5.037772768650, 5.038531444649, 5.039386542352, 5.039767126872, 5.040100412134, 5.040290975754, 5.040576978057, 5.041436116778, 5.041961983902, 5.042680103145, 5.043495409683, 5.044071843103, 5.044889769029, 5.045178816948, 5.046095406587, 5.046820398934, 5.047498152137, 5.048322562657, 5.049099910634, 5.049683836075, 5.050366076201, 5.051000545973, 5.051587034221, 5.052468254304, 5.052762392129, 5.053253064966, 5.053990115234, 5.054728418492, 5.055122694036, 5.055764156207, 5.056208801071, 5.056802369931, 5.057545473658, 5.057991946978, 5.058637664288, 5.059383917663, 5.060031709449, 5.060880282352, 5.061230177217, 5.061830653610, 5.062482107983, 5.063084319055, 5.063737658097, 5.064291252157, 5.064744718216, 5.065400561782, 5.065956280964, 5.066563332174, 5.067221929939, 5.067830754079, 5.068542129311, 5.069101885590, 5.069203737017, 5.069713352754, 5.070121476343, 5.070683273247, 5.071399340155, 5.071962793593, 5.072578304950, 5.073246094810, 5.073811950893, 5.074481641823, 5.075100735986, 5.075617322798, 5.075824129698, 5.076341578207, 5.076807809579, 5.077378328477, 5.077793722561, 5.078469586499, 5.078990198503, 5.079354998593, 5.079824477990, 5.080398976216, 5.080764962077, 5.081235968972, 5.081864773834, 5.082599532045, 5.083125121461, 5.084020085860, 5.084600164788, 5.085339569411, 5.085656842881, 5.086663074067, 5.087459117209, 5.088150203501, 5.088416299019, 5.088948979672, 5.089268901957, 5.089962876447, 5.090497458595, 5.090979145789, 5.091622227568, 5.092212556889, 5.092749917119, 5.093395628275, 5.093988374729, 5.094258072608, 5.094635931233, 5.095068172604, 5.095663207798, 5.096042291477, 5.096638663745, 5.097290187030, 5.097779471721, 5.098105968026, 5.098541678604, 5.098868748645, 5.099523628611, 5.100124803979, 5.100891141807, 5.101549080802, 5.102043189993, 5.102427886174, 5.102812923420, 5.103363569470, 5.103639154531, 5.104522203724, 5.104964402548, 5.105628546144, 5.106349183014, 5.106793246940, 5.107182175690, 5.107571453055, 5.107738393084, 5.108351056129, 5.108852969551, 5.109690783100, 5.110306208556, 5.110754339107, 5.111371274615, 5.111708154644, 5.112157734893, 5.112551499750, 5.112945621949, 5.113565680371, 5.113847818029, 5.114582234889, 5.114921615851, 5.115657852353, 5.116338564846, 5.116679321617, 5.117418546446, 5.117703199062, 5.118444170207, 5.119186407719, 5.119701008574, 5.120330794368, 5.121076266032, 5.121650577782, 5.121995529732, 5.122340755888, 5.122628654130, 5.122916743349, 5.123493495734, 5.124013228572, 5.124591439923, 5.124996646400, 5.125518182301, 5.126330707293, 5.127261172527, 5.127610611582, 5.128427064454, 5.128835867197, 5.129420539447, 5.130299026326, 5.130768280269, 5.131238041788, 5.131767131578, 5.132473588800, 5.132709330145, 5.133122185663, 5.133594501622, 5.134185619832, 5.134718315004, 5.135370275455, 5.136260892655, 5.136736636350, 5.137451230475, 5.138226703281, 5.138525331143, 5.139183036135, 5.140081514799, 5.140681534903, 5.140801638466, 5.141342515909, 5.141884067810, 5.142366014850, 5.142908845327, 5.143512787131, 5.144238627660, 5.144541419614, 5.145026326274, 5.145511774956, 5.146180154143, 5.146666894998, 5.147032308971, 5.147764060588, 5.148252580867, 5.148864004988, 5.149231273071, 5.149844077578, 5.150212175762, 5.150519162756, 5.151379882566, 5.152057361155, 5.152427340858, 5.152921137934, 5.153477331584, 5.153848523471, 5.154405906840, 5.155026061853, 5.155957957959, 5.156580334795, 5.156891858000, 5.157640426669, 5.158265221025, 5.158640529545, 5.159329438667, 5.159705668389, 5.160145015398, 5.160899217293, 5.161591721506, 5.162348442154, 5.162853656091, 5.163486001109, 5.163929191085, 5.164372833790, 5.164689999131, 5.164943898280, 5.165261481096, 5.165960981841, 5.166343003107, 5.166725360709, 5.167171870461, 5.167682730078, 5.168386144691, 5.168834366091, 5.169411331315, 5.170246081075, 5.170631892011, 5.171146840323, 5.171404543628, 5.172114017210, 5.172824651701, 5.173342208124, 5.173990022209, 5.174768676800, 5.175158528246, 5.175483671993, 5.176134690676, 5.176917203467, 5.177308989224, 5.177635747558, 5.178028182358, 5.178355482458, 5.178945244953, 5.179470151476, 5.179864248130, 5.180521871638, 5.181246409502, 5.181444220702, 5.182302445237, 5.182765269575, 5.183294816334, 5.183957659078, 5.184488662637, 5.185020316239, 5.185486047632, 5.185752404268, 5.186285608119, 5.187153463033, 5.187688390869, 5.187956102070, 5.188424994129, 5.188961491396, 5.189700258960, 5.190372958106, 5.190979279516, 5.191383964573, 5.191789027076, 5.192059278785, 5.192397330084, 5.193006486318, 5.193955764252, 5.194567111868, 5.194975155570, 5.195519810894, 5.196065150136, 5.196611175016, 5.196952789509, 5.197705288603, 5.198321940964, 5.199351644636, 5.199695422444, 5.200246033588, 5.200728391650, 5.201487466969, 5.202178688636, 5.202801730161, 5.203078924670, 5.203633845022, 5.203981030653, 5.204606665069, 5.205372555336, 5.205790883654, 5.206488994207, 5.207538268653, 5.208098919990, 5.208519883980, 5.209011524911, 5.209503723033, 5.210137369954, 5.210630846409, 5.211054272976, 5.211336786879, 5.211690187893, 5.212043876716, 5.212681243375, 5.213319546803, 5.213958789757, 5.214385475053, 5.214670164989, 5.215026290046, 5.215382707367, 5.216310763653, 5.216954427885, 5.217455715990, 5.218172847068, 5.218388217507, 5.219322725567, 5.220042948753, 5.220908796155, 5.221342368053, 5.221704008911, 5.222065951162, 5.222428195309, 5.223008415144, 5.223589411193, 5.223807485253, 5.224243962156, 5.224608028303, 5.225191169689, 5.225410049735, 5.225994269742, 5.226506107729, 5.227091805029, 5.227458267359, 5.227825039175, 5.228339040651, 5.228927216779, 5.229442525149, 5.230032198671, 5.230474979829, 5.231435890486, 5.231657941355, 5.232250632654, 5.232695682547, 5.233438447362, 5.233959139619, 5.234182484690, 5.234629519708, 5.235077015350, 5.235674394374, 5.236572006437, 5.237546517636, 5.238673677579, 5.238899461042, 5.239652921470, 5.239954672034, 5.240256632402, 5.240785568766, 5.241466577763, 5.242452146531, 5.242679903453, 5.243363891754, 5.244048958996, 5.244887733605, 5.245345930745, 5.246110668540, 5.246646787359, 5.247260306065, 5.247874692702, 5.248412994918, 5.249029015563, 5.249491605149, 5.250186414707, 5.250882337644, 5.251191995041, 5.251579377532, 5.252432837262, 5.252899068635, 5.253833035623, 5.254378778693, 5.254847104923, 5.255159603215, 5.255707016877, 5.256333478554, 5.257039329786, 5.257510535418, 5.258218303857, 5.258454483224, 5.258927227627, 5.259716280318, 5.260269471435, 5.261060968797, 5.261774551858, 5.262250926108, 5.263125638352, 5.263364502313, 5.264401100302, 5.264800451578, 5.265040238728, 5.265440178421, 5.265840486756, 5.266561972909, 5.267043630424, 5.267847581935, 5.268330668171, 5.268814292366, 5.269298455718, 5.269863996003, 5.270268404713, 5.270916242956, 5.271240524832, 5.271889815900, 5.272377422031, 5.272865576240, 5.273272790973, 5.273843533827, 5.274333339686, 5.274741933640, 5.275150912371, 5.275642195774, 5.276216063035, 5.276872841204, 5.277448337999, 5.277695213126, 5.278271801427, 5.278684119394, 5.279096829187, 5.279509931550, 5.280006173632, 5.280834505911, 5.281332264684, 5.281830594609, 5.282246306789, 5.282828973168, 5.283579266153, 5.284247283177, 5.284748971121, 5.285418791160, 5.285921835018, 5.286257521591, 5.286845598163, 5.287266140930, 5.287434472127, 5.287855585785, 5.288445831750, 5.288952396133, 5.289713352297, 5.290560425868, 5.290984583028, 5.291239276310, 5.291749111409, 5.292004253577, 5.292514988033, 5.293196902963, 5.293879890297, 5.294649537114, 5.294992040667, 5.295849483160, 5.296021174992, 5.296622631488, 5.297138829427, 5.297914278564, 5.298518364379, 5.298950369271, 5.299469343021, 5.300335679798, 5.301029995664, 5.301464507438, 5.302073555194, 5.302596276800, 5.303468880030, 5.304168227173, 5.304781081095, 5.305394801066, 5.305745887975, 5.306009389539, 5.306624848975, 5.306976932076, 5.307241181845, 5.308034897233, 5.308564847856, 5.309095445945, 5.309272456130, 5.310158590863, 5.310691140876, 5.311224344727, 5.311580177997, 5.312025379965, 5.312560225454, 5.313185045493, 5.313721321933, 5.314079207805, 5.315064917359, 5.315603521581, 5.315962962513, 5.316592700868, 5.317043073699, 5.317584138323, 5.318035541005, 5.318849250068, 5.319392571008, 5.319845858266, 5.320390428220, 5.320662969479, 5.321208565634, 5.321572677566, 5.321937095026, 5.322119418488, 5.322301818526, 5.323123568027, 5.323946875348, 5.324496615272, 5.325230685985, 5.325598187155, 5.325874017257, 5.326426203577, 5.327348077160, 5.327717375211, 5.328364403398, 5.328734567053, 5.329475842218, 5.330125497510, 5.330218384791, 5.330590132712, 5.331148351917, 5.331614083310, 5.332453660488, 5.333014281670, 5.333575627481, 5.334419008982, 5.334794371565, 5.335452037753, 5.335640125449, 5.336676066372, 5.336959025106, 5.337242168318, 5.337619979984, 5.337998120610, 5.338281942305, 5.339419087573, 5.340273904762, 5.340939927759, 5.342466112442, 5.343039817257, 5.343326954115, 5.343518484210, 5.343997679317, 5.344381416459, 5.345342245350, 5.345919764693, 5.346884006834, 5.347560252411, 5.348334396077, 5.348721986002, 5.349498205122, 5.349498205122, 5.350081281265, 5.350275814070, 5.350762527650, 5.351347304869, 5.351932870551, 5.352421445788, 5.353400248280, 5.354086724966, 5.354577730651, 5.354774288465, 5.354872600742, 5.355364496232, 5.356349961783, 5.356843534380, 5.357436562896, 5.358030402298, 5.358228529346, 5.358525889496, 5.358724242768, 5.359319846722, 5.359717370303, 5.360015751958, 5.360413913327, 5.360812440064, 5.361410916707, 5.361910278015, 5.362710452322, 5.363311552047, 5.364114314719, 5.364516253185, 5.365019199949, 5.365825128237, 5.366127737342, 5.366834646316, 5.368049173741, 5.368759219764, 5.369165482172, 5.369470428573, 5.369572124975, 5.369877357141, 5.370488465800, 5.371100435579, 5.371815491927, 5.372429335819, 5.373146585333, 5.373865021365, 5.374584647846, 5.374996398985, 5.375408540873, 5.375821074252, 5.376647318462, 5.377164520478, 5.377578726024, 5.377682339166, 5.377889639639, 5.378304537671, 5.379031564356, 5.379343520180, 5.380072289709, 5.380593589113, 5.381324461115, 5.381742655160, 5.382370702242, 5.382475465114, 5.382580253263, 5.383314479104, 5.384365531123, 5.384470776363, 5.385102783967, 5.386158178124, 5.386898483033, 5.387216143280, 5.387534036047, 5.388064374960, 5.388807939132, 5.389446294683, 5.389659288548, 5.389765824666, 5.390619055749, 5.391153177674, 5.391580948683, 5.391687957303, 5.392652223232, 5.393081474052, 5.393833685392, 5.394156460942, 5.394479476563, 5.394802732612, 5.395234115296, 5.395342027952, 5.395881993808, 5.396098168268, 5.396855627380, 5.397397479580, 5.398157210218, 5.398374520446, 5.399027104313, 5.399244850360, 5.399353764338, 5.399789693591, 5.400444409014, 5.401318901093, 5.401428336518, 5.401756808346, 5.401975927666, 5.402633949734, 5.402963335022, 5.403292970319, 5.403952992455, 5.404283380057, 5.405165644417, 5.405496956180, 5.406271001292, 5.406713932980, 5.406824736522, 5.407268233606, 5.407601153884, 5.407712184048, 5.408045444953, 5.408601448719, 5.409046764812, 5.409604052816, 5.410162056853, 5.410608976863, 5.410944468948, 5.411392195257, 5.411728293158, 5.412289034981, 5.412738150307, 5.413187730557, 5.413750361134, 5.414200990987, 5.415329615536, 5.415781887883, 5.416008200802, 5.416461180746, 5.416574499593, 5.416914633652, 5.417595701981, 5.418277840051, 5.418619311290, 5.419417123186, 5.419873674588, 5.420445039599, 5.420788219768, 5.421246215574, 5.421819390372, 5.422393322637, 5.423083044035, 5.423658649794, 5.424696666578, 5.425390058660, 5.426084559578, 5.426780172886, 5.427128397800, 5.427476902150, 5.427942010074, 5.428757149440, 5.428990327691, 5.429457060118, 5.430041181903, 5.430509045651, 5.431211787685, 5.431680914905, 5.432503108896, 5.433326862394, 5.433916215832, 5.434270212169, 5.434742656580, 5.435452288244, 5.436281660034, 5.436756298860, 5.437588167050, 5.438064236686, 5.438898616351, 5.439495584805, 5.439854160151, 5.440093374964, 5.440572200241, 5.441051554022, 5.441651491238, 5.442252258359, 5.442733471130, 5.443456291516, 5.443576878629, 5.444180316939, 5.444784594874, 5.445631999010, 5.446116973356, 5.446723953863, 5.447088549783, 5.447575154296, 5.448672011996, 5.448916134814, 5.449527042893, 5.449894000652, 5.450383760481, 5.450506286785, 5.451119437362, 5.451733454829, 5.452594540333, 5.453087356819, 5.454074670644, 5.454569170535, 5.454816631785, 5.455435902504, 5.456304367691, 5.456304367691, 5.456801414362, 5.458046525542, 5.458420756053, 5.458920232223, 5.459920911196, 5.460422116655, 5.460547508451, 5.461300620458, 5.461929212957, 5.462432742847, 5.462810773756, 5.463189134008, 5.463946844841, 5.464452720823, 5.464705879957, 5.464959186749, 5.465593100862, 5.466227941615, 5.466863711721, 5.467755356373, 5.468010448587, 5.468776625467, 5.469416140355, 5.470184803355, 5.470954829234, 5.471340354765, 5.471726222833, 5.472241247479, 5.472370099129, 5.473144012874, 5.473660722610, 5.473919308198, 5.474307475495, 5.474695990042, 5.475084852460, 5.475603877896, 5.475733731233, 5.476513667657, 5.476773958034, 5.477164686339, 5.477425367309, 5.477686204843, 5.478077755116, 5.478208350361, 5.478469658721, 5.478992747591, 5.479516467259, 5.479778564118, 5.480565805086, 5.481091426309, 5.481222931073, 5.482144581070, 5.482408269288, 5.482804102050, 5.483200295918, 5.483464626104, 5.483596851553, 5.484126156288, 5.484523558618, 5.485053994692, 5.485585079420, 5.485717952140, 5.486649201194, 5.487449007096, 5.488250288655, 5.488383979431, 5.488651484510, 5.489320968968, 5.490259984430, 5.490528647897, 5.490931954983, 5.491604966867, 5.492009275180, 5.492548939098, 5.493089274448, 5.493494967595, 5.494307492588, 5.494850021680, 5.495257363728, 5.495665088198, 5.496617936526, 5.496890563329, 5.497436330893, 5.498119506245, 5.498393077581, 5.499351936628, 5.499900808084, 5.500725418108, 5.501138311007, 5.501276029252, 5.501827339363, 5.502241281713, 5.502379350219, 5.502932063601, 5.503070351927, 5.503208684300, 5.503901007867, 5.504594436854, 5.505149978320, 5.505567101274, 5.506123889147, 5.506960411682, 5.507239610973, 5.507938395487, 5.508498233763, 5.508918586577, 5.509339346644, 5.510041520575, 5.510885630621, 5.511590311097, 5.512013668871, 5.512296136836, 5.512437439744, 5.513144644723, 5.513569521146, 5.514420523015, 5.514988785421, 5.516127545777, 5.516412703031, 5.516698047642, 5.517126416391, 5.517698232777, 5.518127589689, 5.518700726667, 5.519705539997, 5.520136886977, 5.520856752021, 5.521289244487, 5.522011023749, 5.522734004575, 5.522878745280, 5.523748203993, 5.523893283160, 5.525055664535, 5.525201181199, 5.525783735924, 5.525929496785, 5.526367073126, 5.527389802404, 5.527828853308, 5.528561592611, 5.528855034839, 5.528855034839, 5.529442514783, 5.529589509024, 5.530030790500, 5.530177984022, 5.530767257493, 5.531652669588, 5.532096053477, 5.532539890493, 5.532836034031, 5.533428927614, 5.534468443026, 5.535063570878, 5.535510452566, 5.536107011014, 5.536554968230, 5.536554968230, 5.536853863274, 5.537602002101, 5.537901618865, 5.538351431937, 5.538951908329, 5.539252458156, 5.540155357612, 5.540607512241, 5.540909210399, 5.541362150974, 5.541664374008, 5.542572307054, 5.542572307054, 5.543785844642, 5.544089759617, 5.544546031221, 5.545155139991, 5.545307550761, 5.545612532853, 5.546375926409, 5.547140664204, 5.547600154089, 5.548060130635, 5.549289121853, 5.550059011227, 5.550675906901, 5.551139154393, 5.551757587366, 5.552376902240, 5.552686891176, 5.553618187778, 5.554551485734, 5.555330769061, 5.555955204082, 5.556424120250, 5.557206774060, 5.558304864359, 5.558461961298, 5.558461961298, 5.559720786764, 5.559878396812, 5.560667306170, 5.560983271612, 5.561615892965, 5.562090964461, 5.562883906952, 5.563678299860, 5.564155634016, 5.564792896759, 5.565271458220, 5.566550206238, 5.567351339987, 5.567832730557, 5.568153954301, 5.568475415813, 5.568958054664, 5.569602408613, 5.570247719998, 5.570732333567, 5.571055709964, 5.571217488503, 5.571703186017, 5.572838597074, 5.572838597074, 5.573814174755, 5.574139854922, 5.575281662668, 5.575608445590, 5.576262750018, 5.577082019232, 5.577574323629, 5.578231598793, 5.578560609780, 5.579219380451, 5.579714115058, 5.580374639112, 5.581036169296, 5.581201709410, 5.582030357785, 5.582694416755, 5.583193128177, 5.583525920900, 5.584192272364, 5.584525831891, 5.585026652029, 5.585862637816, 5.586365002801, 5.587203571283, 5.587539452570, 5.588043762070, 5.588043762070, 5.588380294037, 5.588885581449, 5.589222766623, 5.590066876669, 5.590235895734, 5.590235895734, 5.590574131329, 5.591590421532, 5.592099459857, 5.592439150514, 5.592609095529, 5.593119329951, 5.593289541390, 5.593800576337, 5.594312213327, 5.594824453782, 5.595508382241, 5.596364810209, 5.596536298655, 5.597394758080, 5.597910649428, 5.598599459218, 5.598599459218, 5.598771832502, 5.599289363227, 5.599634726650, 5.599807511407, 5.600326278519, 5.601018933342, 5.601192269797, 5.601886308270, 5.602407565962, 5.602581457649, 5.602755418990, 5.603975103391, 5.604498875694, 5.604848408495, 5.605023280445, 5.605373235728, 5.606073993414, 5.606600304707, 5.607303046740, 5.607303046740, 5.607830850510, 5.608183076387, 5.608359296508, 5.608888386297, 5.609594843520, 5.610125441609, 5.610302451794, 5.610479534154, 5.610833915635, 5.611543547300, 5.611898798429, 5.612254340391, 5.612610173661, 5.613322716039, 5.614393726402, 5.615825861193, 5.616903070051, 5.617262734239, 5.617802789623, 5.618523909725, 5.618884919290, 5.619426996933, 5.620150821237, 5.620513186283, 5.620875853930, 5.621420423884, 5.621965677543, 5.622329560666, 5.623240604595, 5.623423042944, 5.624336386039, 5.624702261783, 5.624702261783, 5.624702261783, 5.625434939277, 5.625985259708, 5.626168854926, 5.627087997030, 5.628193541493, 5.628378072824, 5.628747370875, 5.629301907424, 5.629857152949, 5.630042392654, 5.630413109264, 5.630969778191, 5.631341287608, 5.631527161560, 5.631899148291, 5.632644078974, 5.633390289608, 5.633763876282, 5.634324859544, 5.634512015109, 5.634512015109, 5.634886568372, 5.635261444945, 5.636200054521, 5.636764195516, 5.636764195516, 5.637517525249, 5.637706062036, 5.638272163982, 5.639028116274, 5.639785386705, 5.639974910811, 5.640354207325, 5.640923773941, 5.641303900426, 5.641494088510, 5.641874714723, 5.641874714723, 5.641874714723, 5.642065153000, 5.642255674820, 5.643400564275, 5.643591672961, 5.644356949779, 5.644740094473, 5.645699437655, 5.645891560853, 5.646083769080, 5.646083769080, 5.646276062411, 5.647238808276, 5.647238808276, 5.647431613821, 5.647817481889, 5.648203693103, 5.648783654661, 5.648783654661, 5.649364391741, 5.649751981666, 5.649945906421, 5.649945906421, 5.650139917808, 5.650722472532, 5.651305809734, 5.651305809734, 5.651500429716, 5.651889931520, 5.652865217090, 5.653647025549, 5.654038458187, 5.654626269441, 5.655018586073, 5.655214877367, 5.655607726315, 5.655607726315, 5.655607726315, 5.655804284129, 5.656197666838, 5.656985502849, 5.657379957447, 5.657577319178, 5.658169943080, 5.658565475422, 5.659952682339, 5.660747365967, 5.661344334421, 5.661742769754, 5.661942124580, 5.662341108974, 5.662540738709, 5.662940273679, 5.663140179083, 5.663740447986, 5.664341547711, 5.664943480561, 5.665345233117, 5.665345233117, 5.665948559653, 5.666351243485, 5.666955970177, 5.666955970177, 5.667359589613, 5.667763584509, 5.668370282370, 5.668775218979, 5.669586226651, 5.670398751644, 5.671009144551, 5.671620396561, 5.672028376377, 5.672436739813, 5.673254620435, 5.673664139071, 5.673664139071, 5.674689628289, 5.675306086138, 5.676129393460, 5.676129393460, 5.676747899828, 5.677780705266, 5.678815972698, 5.679230771661, 5.679438319805, 5.679853713889, 5.680685695910, 5.681727919788, 5.683191247947, 5.684029654543, 5.684659523373, 5.685079944008, 5.685711339053, 5.686132779631, 5.686765708305, 5.686976889677, 5.687611050629, 5.688670047696, 5.688882157338, 5.689519108537, 5.689944262249, 5.690156995284, 5.691008971000, 5.691649051413, 5.692290076595, 5.692717952967, 5.693360558976, 5.693574972449, 5.694004117229, 5.694218848745, 5.694648630553, 5.695294101787, 5.695724949523, 5.696156225111, 5.696372023616, 5.696587929403, 5.696803942580, 5.697020063252, 5.697885623044, 5.698752911364, 5.698970004336, 5.700057099977, 5.700492701300, 5.700710665912, 5.700928739973, 5.701146923590, 5.702020755841, 5.702020755841, 5.702239488901, 5.703115524461, 5.703115524461, 5.703554205794, 5.703993330686, 5.704432900038, 5.705093089395, 5.705974905905, 5.706195640081, 5.707743928644, 5.708187312533, 5.708853238268, 5.708853238268, 5.708853238268, 5.709520186669, 5.709742730606, 5.710634048480, 5.711303739410, 5.711527199400, 5.712198270070, 5.712422190921, 5.713094647028, 5.713992877921, 5.715118285345, 5.715343717212, 5.715569266156, 5.716020715762, 5.716020715762, 5.716472635138, 5.717151397165, 5.717831221695, 5.717831221695, 5.718512112060, 5.718739312945, 5.719421629632, 5.719877103698, 5.719877103698, 5.720105019988, 5.720561211713, 5.721703791909, 5.722161666998, 5.722620025333, 5.723768042078, 5.723768042078, 5.723768042078, 5.723998010038, 5.724228099835, 5.724688645458, 5.725380380909, 5.726073219899, 5.726304412070, 5.726767165957, 5.727694155598, 5.727926212500, 5.728623128106, 5.730020323355, 5.730720610228, 5.731188096260, 5.731656086049, 5.731890270192, 5.732124580681, 5.733063088841, 5.733533104560, 5.734003629505, 5.735182176990, 5.735890843694, 5.735890843694, 5.736363931412, 5.736363931412, 5.736600668666, 5.737074530668, 5.737548910270, 5.737548910270, 5.738261452647, 5.738975166008, 5.738975166008, 5.739451627363, 5.739451627363, 5.740645072692, 5.741362717276, 5.741602195905, 5.742321425131, 5.742321425131, 5.743041847439, 5.743522793758, 5.744245213357, 5.745210312603, 5.745451922891, 5.745693667669, 5.746177561292, 5.746177561292, 5.746661994674, 5.747146969020, 5.747632485540, 5.748361779552, 5.748848656825, 5.749579997691, 5.750312572195, 5.750801642609, 5.751046384504, 5.751291264399, 5.751291264399, 5.752272167090, 5.752763450493, 5.753009300758, 5.753994095924, 5.754240644033, 5.754487332186, 5.754734160543, 5.755228238505, 5.755722879198, 5.756218083906, 5.756713853917, 5.756961951314, 5.757458571702, 5.757707095017, 5.757955760630, 5.758453519403, 5.758702612890, 5.758951849328, 5.759950227887, 5.760450279160, 5.761201437286, 5.762205006726, 5.762707662433, 5.762959208621, 5.763210900591, 5.763966852882, 5.763966852882, 5.764724123313, 5.765229704839, 5.765735875621, 5.767003889608, 5.767512133647, 5.768020973169, 5.768785352037, 5.769551078622, 5.770062314092, 5.770574152079, 5.771086594005, 5.772113295386, 5.772370350429, 5.772884917411, 5.773142429711, 5.773657912836, 5.773915884024, 5.773915884024, 5.775466937394, 5.775466937394, 5.775466937394, 5.775725985706, 5.775985188627, 5.777023550107, 5.777283528852, 5.777803953698, 5.778585762158, 5.779368980552, 5.780153613976, 5.780415473786, 5.780677491581, 5.781464494783, 5.782516055786, 5.782779344356, 5.783306400830, 5.783834097714, 5.783834097714, 5.784098186796, 5.784098186796, 5.784626847217, 5.785421046430, 5.785951320588, 5.786216700665, 5.786482243004, 5.787013815263, 5.787013815263, 5.787546038960, 5.788078915692, 5.788612447063, 5.788879458742, 5.789146634685, 5.789681480174, 5.790484985457, 5.791021482724, 5.791558643561, 5.792096469614, 5.792365632611, 5.792634962531, 5.792634962531, 5.792634962531, 5.793174123968, 5.793713955588, 5.794254459057, 5.795337488252, 5.796151536254, 5.796695083862, 5.797239312607, 5.797511682940, 5.797784224199, 5.798876102793, 5.798876102793, 5.799149501909, 5.799696817018, 5.799970733446, 5.800244822747, 5.801342913046, 5.801342913046, 5.801617869992, 5.801893001127, 5.802443786846, 5.802995271977, 5.803547458297, 5.804931003531, 5.804931003531, 5.805485658118, 5.806041021981, 5.806875401646, 5.807432546663, 5.807711387432, 5.808828544271, 5.809388202186, 5.809668301830, 5.810790510418, 5.811071516239, 5.811352704000, 5.812197361282, 5.813326132500, 5.814174640387, 5.814741234703, 5.815024809302, 5.815024809302, 5.815876645760, 5.815876645760, 5.816160962944, 5.816445466381, 5.816730156317, 5.817585347565, 5.819014419213, 5.819014419213, 5.820448208835, 5.821310760224, 5.822175028135, 5.822752163744, 5.822752163744, 5.823908740944, 5.824778199657, 5.826231176863, 5.826522356547, 5.827105302248, 5.827397068790, 5.828273546347, 5.828566099057, 5.828858848972, 5.829444941479, 5.829738284605, 5.829738284605, 5.830619504688, 5.830619504688, 5.830913642513, 5.831502516477, 5.831502516477, 5.832387327272, 5.832682665252, 5.833569886157, 5.833866029695, 5.834162375310, 5.834755673875, 5.834755673875, 5.835647144216, 5.836242476018, 5.837435593477, 5.837734385702, 5.838033383636, 5.838332587562, 5.839231438139, 5.839231438139, 5.839531468881, 5.840132152907, 5.841034739617, 5.841939206063, 5.842845560094, 5.843148098930, 5.843148098930, 5.843450848668, 5.843753809603, 5.843753809603, 5.844360366240, 5.845576026885, 5.845880474484, 5.846185135655, 5.846490010699, 5.847711655617, 5.848017604543, 5.848630149753, 5.848936746646, 5.849550590539, 5.850165303284, 5.850165303284, 5.850780887345, 5.851089006891, 5.852323675759, 5.852942328972, 5.853871964322, 5.854492828590, 5.855114581713, 5.855425792390, 5.855737226238, 5.856048883576, 5.856360764725, 5.856360764725, 5.857297754262, 5.858236769724, 5.859177819891, 5.859491956962, 5.859806321421, 5.859806321421, 5.860120913599, 5.860435733824, 5.860435733824, 5.861066059743, 5.861066059743, 5.862013267276, 5.862329462763, 5.862962545211, 5.863913902616, 5.865185629680, 5.865185629680, 5.865504144165, 5.866141874797, 5.867100230056, 5.867100230056, 5.867740310469, 5.868381335651, 5.868702203402, 5.869023308394, 5.869344650978, 5.869344650978, 5.870310107801, 5.870954940112, 5.870954940112, 5.871923987331, 5.871923987331, 5.872247484167, 5.872571222148, 5.873219422988, 5.873543886568, 5.874518734299, 5.875495775166, 5.876475019057, 5.877129077136, 5.877456475931, 5.877784121727, 5.877784121727, 5.878768544850, 5.879097182385, 5.880744110722, 5.881074247174, 5.882066164960, 5.882397308310, 5.883724412419, 5.885055584287, 5.885055584287, 5.885389015768, 5.885389015768, 5.885722703438, 5.886390848927, 5.887394998465, 5.887730231583, 5.888401475120, 5.889410289701, 5.889410289701, 5.889410289701, 5.889747082647, 5.890759031412, 5.891096872333, 5.891773343625, 5.892111974817, 5.892111974817, 5.893129455521, 5.893129455521, 5.894489815230, 5.894830572001, 5.895512888688, 5.895854449446, 5.896538377905, 5.896538377905, 5.896880746454, 5.897223385117, 5.897566294319, 5.897909474488, 5.897909474488, 5.898252926054, 5.898940645092, 5.899629454882, 5.899974269892, 5.900319358891, 5.900664722314, 5.900664722314, 5.901010360599, 5.902048929006, 5.902395671126, 5.903089986992, 5.903437561626, 5.903437561626, 5.903437561626, 5.903437561626, 5.904133546521, 5.905179619645, 5.906928693624, 5.907279355316, 5.907981529247, 5.908684840303, 5.909036923404, 5.909036923404, 5.910094888561, 5.910094888561, 5.910094888561, 5.910448117114, 5.910801633195, 5.910801633195, 5.910801633195, 5.912573542964, 5.913284336055, 5.913640169325, 5.915066425063, 5.915066425063, 5.915066425063, 5.915423722066, 5.915423722066, 5.915781313261, 5.916139199133, 5.916855856857, 5.917214629684, 5.917933065715, 5.919373513078, 5.919373513078, 5.919373513078, 5.919734372660, 5.920456992597, 5.920456992597, 5.920456992597, 5.920818753952, 5.921180816901, 5.921180816901, 5.921905849594, 5.922632094716, 5.923359556330, 5.925183559355, 5.925549281045, 5.927015255372, 5.927382523455, 5.927382523455, 5.928117992694, 5.928486194905, 5.929223537157, 5.929223537157, 5.930331903088, 5.930701987884, 5.930701987884, 5.931072388318, 5.931443104928, 5.931814138254, 5.932557157224, 5.932557157224, 5.932929143955, 5.934047019686, 5.934047019686, 5.934793871946, 5.935167780261, 5.935542010773, 5.935542010773, 5.935916564036, 5.937042165916, 5.937042165916, 5.937794191180, 5.938924676370, 5.939679971312, 5.939679971312, 5.940058111938, 5.940436582099, 5.941573975543, 5.942333896090, 5.942714355582, 5.943095148664, 5.944239535312, 5.945004138471, 5.945004138471, 5.945386945443, 5.945386945443, 5.945386945443, 5.946153573148, 5.946921556517, 5.946921556517, 5.947690900353, 5.948847477553, 5.949620243739, 5.950394387405, 5.951558196450, 5.952725132616, 5.954285941059, 5.955460239608, 5.956244873031, 5.956637721979, 5.956637721979, 5.957030926607, 5.957030926607, 5.958212681028, 5.958212681028, 5.958212681028, 5.958212681028, 5.958607314842, 5.959397659886, 5.959793372425, 5.959793372425, 5.960585880824, 5.960982678003, 5.962175249412, 5.962175249412, 5.963371104638, 5.963770455914, 5.964570261815, 5.964570261815, 5.964970717798, 5.965772739229, 5.965772739229, 5.966174306047, 5.967784296702, 5.968187728670, 5.969805214643, 5.970210529168, 5.970210529168, 5.970210529168, 5.971428747307, 5.971428747307, 5.971835580576, 5.971835580576, 5.972650392225, 5.973058372041, 5.973058372041, 5.973466735477, 5.973466735477, 5.973875483255, 5.974284616099, 5.974284616099, 5.974284616099, 5.974694134735, 5.974694134735, 5.975104039893, 5.975514332301, 5.975514332301, 5.976336081802, 5.977984260182, 5.978397283972, 5.979224511806, 5.979638717352, 5.980468315469, 5.981299501334, 5.981299501334, 5.981299501334, 5.981299501334, 5.982549270489, 5.982549270489, 5.982549270489, 5.982966660701, 5.983384452443, 5.983802646488, 5.986320302709, 5.986320302709, 5.986741334717, 5.988006885341, 5.988006885341, 5.988852639224, 5.989276134608, 5.989276134608, 5.989276134608, 5.989700043360, 5.990124366288, 5.990974257913, 5.991399828238, 5.991825815994, 5.992679047077, 5.992679047077, 5.993106292052, 5.993533957751, 5.993533957751, 5.993533957751, 5.994390554640, 5.994819487496, 5.994819487496, 5.995678626217, 5.996539467891, 5.996970529446, 5.996970529446, 5.999565922521, 6.000000000000, 6.000000000000, 6.001304841688, 6.001740661576, 6.002176919254, 6.002176919254, 6.002176919254, 6.002613615603, 6.003050751505, 6.003050751505, 6.003926345515, 6.004364805402, 6.005682847330, 6.005682847330, 6.006563769502, 6.007446482168, 6.008773924308, 6.008773924308, 6.009217308197, 6.009661145212, 6.010995384301, 6.011887159732, 6.011887159732, 6.011887159732, 6.012780770092, 6.013676222949, 6.014573525917, 6.015022873585, 6.015472686656, 6.015922966097, 6.015922966097, 6.015922966097, 6.016373712875, 6.016373712875, 6.016373712875, 6.016824927962, 6.017276612331, 6.017728766960, 6.018181392829, 6.019088062223, 6.020451625296, 6.020451625296, 6.021363051616, 6.022276394711, 6.023191662662, 6.024108863598, 6.024108863598, 6.024568191491, 6.025028005702, 6.025488307263, 6.026410376573, 6.027334407734, 6.028260409112, 6.029188389128, 6.031050319019, 6.031050319019, 6.032452023781, 6.032452023781, 6.032452023781, 6.033389013318, 6.033858267261, 6.034328028780, 6.035740369803, 6.036212172654, 6.036684488614, 6.036684488614, 6.038104526332, 6.039053804266, 6.039053804266, 6.040005161672, 6.040005161672, 6.040005161672, 6.040958607679, 6.040958607679, 6.040958607679, 6.040958607679, 6.040958607679, 6.041436116778, 6.043351420795, 6.043831569525, 6.044793462458, 6.045757490561, 6.047207556956, 6.047691990338, 6.047691990338, 6.049635145624, 6.050609993355, 6.050609993355, 6.051098239030, 6.051587034221, 6.051587034221, 6.051587034221, 6.052566278113, 6.053547734987, 6.053547734987, 6.054039296422, 6.055024091588, 6.055517327850, 6.055517327850, 6.056505484094, 6.056505484094, 6.056505484094, 6.056505484094, 6.057991946978, 6.057991946978, 6.058488567366, 6.058985756294, 6.059981844992, 6.059981844992, 6.060480747381, 6.061480274824, 6.062482107983, 6.063989204285, 6.063989204285, 6.064996848546, 6.064996848546, 6.065501548756, 6.065501548756, 6.066006836169, 6.066006836169, 6.066006836169, 6.066006836169, 6.066512712151, 6.068033885272, 6.068033885272, 6.068033885272, 6.069050968832, 6.069560405233, 6.070070439915, 6.071092309756, 6.071092309756, 6.071092309756, 6.071604147743, 6.072116589669, 6.072116589669, 6.073143291050, 6.073657553374, 6.074687908500, 6.075204004202, 6.075720713938, 6.075720713938, 6.076755981370, 6.076755981370, 6.076755981370, 6.077274542007, 6.077793722561, 6.077793722561, 6.077793722561, 6.077793722561, 6.077793722561, 6.077793722561, 6.077793722561, 6.077793722561, 6.078313524516, 6.078833949362, 6.079876673709, 6.081445469450, 6.082494490447, 6.083019952680, 6.083546051450, 6.083546051450, 6.084072788303, 6.085128182460, 6.085128182460, 6.086186147616, 6.086186147616, 6.086186147616, 6.086186147616, 6.086186147616, 6.086716098240, 6.087246696329, 6.088309841246, 6.088309841246, 6.089375595111, 6.089375595111, 6.089909454406, 6.090443970759, 6.090443970759, 6.091514981121, 6.091514981121, 6.092588639225, 6.093126465278, 6.093664958195, 6.093664958195, 6.094743951252, 6.094743951252, 6.095284454721, 6.095825631716, 6.095825631716, 6.095825631716, 6.096367483916, 6.096910013008, 6.099086932262, 6.099086932262, 6.099086932262, 6.099086932262, 6.099632871344, 6.100179497573, 6.100179497573, 6.100179497573, 6.100179497573, 6.100179497573, 6.100726812682, 6.101274818411, 6.102372908710, 6.102922996791, 6.103473782510, 6.103473782510, 6.103473782510, 6.104025267641, 6.105130343255, 6.105130343255, 6.105683937316, 6.105683937316, 6.106793246940, 6.107905397310, 6.107905397310, 6.107905397310, 6.108462542327, 6.108462542327, 6.109020403010, 6.109020403010, 6.109578981199, 6.110138278742, 6.111259039317, 6.111820506082, 6.112382699664, 6.112945621949, 6.112945621949, 6.114073660199, 6.114638779969, 6.115204636051, 6.115204636051, 6.115204636051, 6.115771230367, 6.115771230367, 6.116338564846, 6.116906641424, 6.116906641424, 6.118045028660, 6.119186407719, 6.119758224105, 6.120904120500, 6.121478204499, 6.122053048371, 6.122053048371, 6.122053048371, 6.122628654130, 6.124360062996, 6.126098402136, 6.126098402136, 6.126679398185, 6.126679398185, 6.127261172527, 6.127261172527, 6.128427064454, 6.128427064454, 6.129011186239, 6.130181792021, 6.130181792021, 6.130181792021, 6.131355561605, 6.131943638177, 6.131943638177, 6.133122185663, 6.133122185663, 6.134303940084, 6.134303940084, 6.134896025359, 6.135488918942, 6.136082623042, 6.136677139880, 6.137272471682, 6.137272471682, 6.137272471682, 6.137868620687, 6.138465589141, 6.140261433803, 6.140261433803, 6.140861702705, 6.142667503569, 6.143271109617, 6.143875555758, 6.144480844332, 6.145086977692, 6.145693958199, 6.146301788224, 6.146910470148, 6.147520006363, 6.148741651281, 6.149353764817, 6.151195298948, 6.151810883009, 6.153044674980, 6.153044674980, 6.153662887870, 6.153662887870, 6.154281982033, 6.156144577377, 6.158015195410, 6.158640529545, 6.158640529545, 6.159893905543, 6.159893905543, 6.160521952626, 6.161150909263, 6.161150909263, 6.161780778092, 6.161780778092, 6.162411561765, 6.162411561765, 6.163043262940, 6.164309428508, 6.164309428508, 6.164309428508, 6.164943898280, 6.165579296318, 6.165579296318, 6.165579296318, 6.166215625344, 6.166852888087, 6.166852888087, 6.166852888087, 6.166852888087, 6.166852888087, 6.167491087294, 6.167491087294, 6.167491087294, 6.170053304058, 6.170053304058, 6.170696227169, 6.170696227169, 6.170696227169, 6.171984935776, 6.172630726946, 6.172630726946, 6.173277479831, 6.173277479831, 6.173277479831, 6.173277479831, 6.173925197299, 6.175874166083, 6.175874166083, 6.175874166083, 6.175874166083, 6.176525770830, 6.177831920632, 6.177831920632, 6.177831920632, 6.177831920632, 6.177831920632, 6.178486471595, 6.179142010560, 6.179142010560, 6.179798540514, 6.180456064458, 6.180456064458, 6.181774106386, 6.181774106386, 6.181774106386, 6.181774106386, 6.182434630440, 6.184422251676, 6.185752404268, 6.185752404268, 6.187086643357, 6.187755303200, 6.187755303200, 6.188424994129, 6.188424994129, 6.188424994129, 6.189095719331, 6.189767482005, 6.189767482005, 6.189767482005, 6.189767482005, 6.191114132640, 6.191114132640, 6.191114132640, 6.191789027076, 6.191789027076, 6.192464971931, 6.192464971931, 6.192464971931, 6.193820026016, 6.193820026016, 6.195179321279, 6.195860567665, 6.195860567665, 6.195860567665, 6.197226274708, 6.197226274708, 6.197910742118, 6.198596289983, 6.198596289983, 6.199282921718, 6.199970640756, 6.199970640756, 6.200659450546, 6.201349354555, 6.201349354555, 6.202040356263, 6.202732459169, 6.202732459169, 6.204119982656, 6.205511953341, 6.206209615309, 6.206908399823, 6.206908399823, 6.208309350980, 6.208309350980, 6.209011524911, 6.209714835967, 6.211124884225, 6.211831628859, 6.211831628859, 6.211831628859, 6.211831628859, 6.211831628859, 6.211831628859, 6.214670164989, 6.215382707367, 6.215382707367, 6.215382707367, 6.217527375834, 6.217527375834, 6.218244625348, 6.218963061379, 6.219682687860, 6.219682687860, 6.220403508742, 6.221125527997, 6.221848749616, 6.221848749616, 6.221848749616, 6.223298816012, 6.223298816012, 6.224025668871, 6.224025668871, 6.224753740260, 6.224753740260, 6.225483034271, 6.225483034271, 6.225483034271, 6.226945306636, 6.226945306636, 6.227678293277, 6.227678293277, 6.227678293277, 6.227678293277, 6.227678293277, 6.228412519119, 6.228412519119, 6.228412519119, 6.229884705213, 6.229884705213, 6.229884705213, 6.230622673924, 6.232102383982, 6.233587152888, 6.235077015350, 6.235077015350, 6.237321436273, 6.237321436273, 6.237321436273, 6.237321436273, 6.237321436273, 6.238072161579, 6.239577516577, 6.239577516577, 6.240332155310, 6.240332155310, 6.240332155310, 6.240332155310, 6.241845378033, 6.243363891754, 6.243363891754, 6.243363891754, 6.246416941107, 6.247183568812, 6.248720896017, 6.248720896017, 6.249491605149, 6.249491605149, 6.250263684431, 6.251037138744, 6.251037138744, 6.251037138744, 6.251037138744, 6.252588192114, 6.253365801062, 6.253365801062, 6.254144804826, 6.254144804826, 6.254144804826, 6.254144804826, 6.254144804826, 6.254144804826, 6.254144804826, 6.254925208418, 6.255707016877, 6.255707016877, 6.256490235272, 6.258060922271, 6.258060922271, 6.258848401148, 6.259637310506, 6.259637310506, 6.260427655550, 6.261219441516, 6.262012673667, 6.262807357295, 6.262807357295, 6.263603497723, 6.263603497723, 6.263603497723, 6.265200170411, 6.265200170411, 6.265200170411, 6.266802734893, 6.266802734893, 6.267606240177, 6.267606240177, 6.268411234813, 6.268411234813, 6.269217724334, 6.269217724334, 6.269217724334, 6.269217724334, 6.269217724334, 6.270025714300, 6.270835210307, 6.271646217979, 6.272458742971, 6.273272790973, 6.273272790973, 6.274905478919, 6.274905478919, 6.275724130399, 6.276544327965, 6.276544327965, 6.276544327965, 6.277366077466, 6.278189384787, 6.279014255846, 6.279014255846, 6.279840696594, 6.279840696594, 6.279840696594, 6.279840696594, 6.279840696594, 6.279840696594, 6.280668713016, 6.282329496998, 6.283162276700, 6.283162276700, 6.283162276700, 6.284832642152, 6.285670240255, 6.288192770959, 6.288192770959, 6.289036881005, 6.289036881005, 6.289036881005, 6.291579099865, 6.291579099865, 6.292429823902, 6.292429823902, 6.293282217663, 6.294136287716, 6.294992040667, 6.294992040667, 6.297569463554, 6.299296282855, 6.300162274133, 6.301029995664, 6.301029995664, 6.301029995664, 6.301029995664, 6.302770657240, 6.302770657240, 6.304518323510, 6.304518323510, 6.304518323510, 6.304518323510, 6.305394801066, 6.306273051076, 6.306273051076, 6.307153080723, 6.307153080723, 6.307153080723, 6.308918507877, 6.308918507877, 6.309803919971, 6.310691140876, 6.310691140876, 6.310691140876, 6.310691140876, 6.310691140876, 6.310691140876, 6.310691140876, 6.311580177997, 6.311580177997, 6.311580177997, 6.312471038785, 6.313363730738, 6.313363730738, 6.313363730738, 6.313363730738, 6.313363730738, 6.313363730738, 6.314258261398, 6.314258261398, 6.315154638356, 6.315154638356, 6.315154638356, 6.315154638356, 6.316052869248, 6.317854923626, 6.317854923626, 6.319664486585, 6.320572103388, 6.320572103388, 6.320572103388, 6.321481620960, 6.323306390375, 6.323306390375, 6.323306390375, 6.323306390375, 6.324221658326, 6.326058001366, 6.326058001366, 6.326058001366, 6.326058001366, 6.326979092871, 6.326979092871, 6.326979092871, 6.327902142064, 6.327902142064, 6.328827157285, 6.329754146926, 6.330683119434, 6.331614083310, 6.331614083310, 6.332547047110, 6.333482019445, 6.334419008982, 6.334419008982, 6.336299074610, 6.336299074610, 6.336299074610, 6.336299074610, 6.336299074610, 6.336299074610, 6.336299074610, 6.337242168318, 6.337242168318, 6.337242168318, 6.337242168318, 6.338187314463, 6.340083799930, 6.340083799930, 6.340083799930, 6.340083799930, 6.341035157336, 6.341035157336, 6.341035157336, 6.341988603343, 6.341988603343, 6.342944147143, 6.342944147143, 6.342944147143, 6.342944147143, 6.342944147143, 6.342944147143, 6.342944147143, 6.342944147143, 6.343901797987, 6.344861565189, 6.344861565189, 6.344861565189, 6.344861565189, 6.345823458122, 6.345823458122, 6.346787486225, 6.347753658997, 6.347753658997, 6.347753658997, 6.347753658997, 6.347753658997, 6.347753658997, 6.348721986002, 6.348721986002, 6.348721986002, 6.349692476868, 6.349692476868, 6.349692476868, 6.349692476868, 6.349692476868, 6.350665141288, 6.351639989019, 6.351639989019, 6.351639989019, 6.351639989019, 6.352617029885, 6.352617029885, 6.352617029885, 6.353596273777, 6.353596273777, 6.355561410532, 6.355561410532, 6.356547323514, 6.356547323514, 6.356547323514, 6.356547323514, 6.356547323514, 6.356547323514, 6.356547323514, 6.356547323514, 6.357535479758, 6.357535479758, 6.357535479758, 6.357535479758, 6.357535479758, 6.358525889496, 6.360513510731, 6.361510743045, 6.362510270488, 6.362510270488, 6.362510270488, 6.364516253185, 6.366531544420, 6.366531544420, 6.367542707815, 6.367542707815, 6.367542707815, 6.367542707815, 6.369572124975, 6.370590400897, 6.370590400897, 6.370590400897, 6.371611069950, 6.371611069950, 6.371611069950, 6.372634143407, 6.372634143407, 6.372634143407, 6.372634143407, 6.372634143407, 6.373659632625, 6.373659632625, 6.374687549038, 6.374687549038, 6.375717904164, 6.375717904164, 6.375717904164, 6.376750709602, 6.377785977034, 6.378823718225, 6.378823718225, 6.380906669373, 6.380906669373, 6.380906669373, 6.380906669373, 6.380906669373, 6.381951903288, 6.381951903288, 6.382999658879, 6.384049948344, 6.384049948344, 6.384049948344, 6.384049948344, 6.386158178124, 6.386158178124, 6.387216143280, 6.389339836910, 6.389339836910, 6.389339836910, 6.390405590775, 6.391473966423, 6.392544976785, 6.392544976785, 6.392544976785, 6.393618634889, 6.393618634889, 6.393618634889, 6.394694953859, 6.394694953859, 6.394694953859, 6.394694953859, 6.395773946916, 6.395773946916, 6.395773946916, 6.395773946916, 6.396855627380, 6.397940008672, 6.397940008672, 6.397940008672, 6.399027104313, 6.399027104313, 6.399027104313, 6.401209493237, 6.402304814074, 6.402304814074, 6.402304814074, 6.402304814074, 6.402304814074, 6.404503778174, 6.404503778174, 6.404503778174, 6.404503778174, 6.404503778174, 6.404503778174, 6.405607449625, 6.405607449625, 6.406713932980, 6.406713932980, 6.408935392974, 6.408935392974, 6.411168274406, 6.413412695328, 6.414539270492, 6.414539270492, 6.415668775632, 6.416801226031, 6.416801226031, 6.416801226031, 6.416801226031, 6.416801226031, 6.416801226031, 6.416801226031, 6.416801226031, 6.416801226031, 6.416801226031, 6.416801226031, 6.417936637088, 6.417936637088, 6.417936637088, 6.419075024324, 6.420216403383, 6.420216403383, 6.420216403383, 6.420216403383, 6.420216403383, 6.421360790032, 6.421360790032, 6.421360790032, 6.422508200163, 6.424812155072, 6.425968732272, 6.428291168191, 6.429457060118, 6.430626090385, 6.431798275933, 6.431798275933, 6.432973633841, 6.432973633841, 6.432973633841, 6.434152181326, 6.434152181326, 6.434152181326, 6.435333935748, 6.436518914606, 6.437707135544, 6.437707135544, 6.437707135544, 6.437707135544, 6.440093374964, 6.441291429467, 6.441291429467, 6.441291429467, 6.441291429467, 6.443697499233, 6.443697499233, 6.443697499233, 6.444905551422, 6.444905551422, 6.446116973356, 6.447331783888, 6.447331783888, 6.448550002027, 6.449771646945, 6.450996737974, 6.450996737974, 6.452225294612, 6.452225294612, 6.453457336522, 6.453457336522, 6.454692883534, 6.454692883534, 6.454692883534, 6.454692883534, 6.457174573041, 6.459670525209, 6.460923901207, 6.460923901207, 6.460923901207, 6.460923901207, 6.462180904927, 6.462180904927, 6.463441557428, 6.463441557428, 6.464705879957, 6.464705879957, 6.464705879957, 6.464705879957, 6.465973893944, 6.465973893944, 6.465973893944, 6.465973893944, 6.465973893944, 6.467245621008, 6.467245621008, 6.467245621008, 6.467245621008, 6.467245621008, 6.467245621008, 6.467245621008, 6.467245621008, 6.467245621008, 6.467245621008, 6.467245621008, 6.471083299722, 6.471083299722, 6.472370099129, 6.473660722610, 6.473660722610, 6.473660722610, 6.473660722610, 6.476253533188, 6.476253533188, 6.476253533188, 6.476253533188, 6.476253533188, 6.476253533188, 6.476253533188, 6.476253533188, 6.476253533188, 6.476253533188, 6.476253533188, 6.476253533188, 6.476253533188, 6.476253533188, 6.477555766494, 6.478861916296, 6.478861916296, 6.478861916296, 6.480172006224, 6.480172006224, 6.480172006224, 6.480172006224, 6.480172006224, 6.481486060122, 6.481486060122, 6.481486060122, 6.481486060122, 6.481486060122, 6.481486060122, 6.481486060122, 6.484126156288, 6.484126156288, 6.485452247340, 6.485452247340, 6.485452247340, 6.485452247340, 6.485452247340, 6.485452247340, 6.485452247340, 6.485452247340, 6.485452247340, 6.486782399932, 6.486782399932, 6.486782399932, 6.486782399932, 6.486782399932, 6.486782399932, 6.488116639021, 6.488116639021, 6.488116639021, 6.488116639021, 6.488116639021, 6.489454989793, 6.489454989793, 6.490797477669, 6.492144128304, 6.492144128304, 6.493494967595, 6.493494967595, 6.493494967595, 6.494850021680, 6.494850021680, 6.494850021680, 6.494850021680, 6.494850021680, 6.494850021680, 6.494850021680, 6.494850021680, 6.494850021680, 6.494850021680, 6.494850021680, 6.496209316943, 6.496209316943, 6.496209316943, 6.496209316943, 6.496209316943, 6.497572880016, 6.497572880016, 6.497572880016, 6.497572880016, 6.497572880016, 6.497572880016, 6.503070351927, 6.503070351927, 6.504455662454, 6.504455662454, 6.504455662454, 6.504455662454, 6.504455662454, 6.504455662454, 6.507239610973, 6.507239610973, 6.507239610973, 6.507239610973, 6.507239610973, 6.507239610973, 6.507239610973, 6.508638306166, 6.508638306166, 6.508638306166, 6.508638306166, 6.508638306166, 6.508638306166, 6.508638306166, 6.508638306166, 6.510041520575, 6.510041520575, 6.510041520575, 6.510041520575, 6.510041520575, 6.510041520575, 6.510041520575, 6.510041520575, 6.511449283500, 6.511449283500, 6.511449283500, 6.511449283500, 6.512861624523, 6.514278573518, 6.514278573518, 6.514278573518, 6.514278573518, 6.514278573518, 6.514278573518, 6.514278573518, 6.514278573518, 6.514278573518, 6.515700160653, 6.515700160653, 6.517126416391, 6.517126416391, 6.517126416391, 6.518557371498, 6.518557371498, 6.518557371498, 6.518557371498, 6.518557371498, 6.519993057043, 6.519993057043, 6.519993057043, 6.519993057043, 6.519993057043, 6.519993057043, 6.521433504406, 6.522878745280, 6.522878745280, 6.522878745280, 6.522878745280, 6.524328811676, 6.524328811676, 6.525783735924, 6.525783735924, 6.527243550683, 6.527243550683, 6.527243550683, 6.527243550683, 6.527243550683, 6.528708288941, 6.528708288941, 6.528708288941, 6.528708288941, 6.530177984022, 6.531652669588, 6.531652669588, 6.531652669588, 6.531652669588, 6.531652669588, 6.531652669588, 6.531652669588, 6.531652669588, 6.533132379646, 6.533132379646, 6.536107011014, 6.536107011014, 6.536107011014, 6.536107011014, 6.536107011014, 6.537602002101, 6.537602002101, 6.537602002101, 6.537602002101, 6.537602002101, 6.539102157243, 6.539102157243, 6.539102157243, 6.539102157243, 6.540607512241, 6.540607512241, 6.540607512241, 6.540607512241, 6.540607512241, 6.543633966871, 6.543633966871, 6.543633966871, 6.543633966871, 6.543633966871, 6.543633966871, 6.543633966871, 6.545155139992, 6.546681659953, 6.546681659953, 6.546681659953, 6.546681659953, 6.548213564476, 6.548213564476, 6.548213564476, 6.549750891681, 6.549750891681, 6.549750891681, 6.549750891681, 6.549750891681, 6.549750891681, 6.549750891681, 6.551293680095, 6.551293680095, 6.551293680095, 6.551293680095, 6.551293680095, 6.552841968658, 6.552841968658, 6.552841968658, 6.552841968658, 6.554395796726, 6.554395796726, 6.555955204082, 6.555955204082, 6.555955204082, 6.555955204082, 6.555955204082, 6.557520230936, 6.557520230936, 6.557520230936, 6.557520230936, 6.559090917935, 6.559090917935, 6.560667306170, 6.560667306170, 6.560667306170, 6.560667306170, 6.562249437180, 6.562249437180, 6.563837352959, 6.563837352959, 6.563837352959, 6.563837352959, 6.563837352959, 6.563837352959, 6.565431095966, 6.565431095966, 6.565431095966, 6.565431095966, 6.567030709126, 6.567030709126, 6.567030709126, 6.567030709126, 6.568636235841, 6.568636235841, 6.568636235841, 6.568636235841, 6.568636235841, 6.568636235841, 6.568636235841, 6.568636235841, 6.570247719998, 6.570247719998, 6.570247719998, 6.570247719998, 6.570247719998, 6.570247719998, 6.570247719998, 6.570247719998, 6.570247719998, 6.571865205971, 6.571865205971, 6.571865205971, 6.571865205971, 6.575118363369, 6.575118363369, 6.576754126063, 6.576754126063, 6.580044251510, 6.580044251510, 6.581698708680, 6.581698708680, 6.581698708680, 6.583359492662, 6.583359492662, 6.583359492662, 6.585026652029, 6.585026652029, 6.585026652029, 6.586700235919, 6.586700235919, 6.586700235919, 6.586700235919, 6.586700235919, 6.586700235919, 6.586700235919, 6.586700235919, 6.588380294037, 6.588380294037, 6.588380294037, 6.588380294037, 6.590066876669, 6.591760034688, 6.591760034688, 6.591760034688, 6.591760034688, 6.591760034688, 6.591760034688, 6.591760034688, 6.591760034688, 6.591760034688, 6.591760034688, 6.591760034688, 6.593459819566, 6.593459819566, 6.593459819566, 6.593459819566, 6.593459819566, 6.593459819566, 6.593459819566, 6.593459819566, 6.593459819566, 6.593459819566, 6.593459819566, 6.593459819566, 6.596879478824, 6.598599459218, 6.598599459218, 6.598599459218, 6.598599459218, 6.600326278519, 6.600326278519, 6.602059991328, 6.602059991328, 6.602059991328, 6.602059991328, 6.603800652904, 6.603800652904, 6.603800652904, 6.603800652904, 6.603800652904, 6.603800652904, 6.603800652904, 6.603800652904, 6.603800652904, 6.605548319174, 6.605548319174, 6.605548319174, 6.605548319174, 6.605548319174, 6.605548319174, 6.605548319174, 6.605548319174, 6.607303046740, 6.610833915635, 6.610833915635, 6.610833915635, 6.614393726402, 6.614393726402, 6.614393726402, 6.614393726402, 6.614393726402, 6.614393726402, 6.614393726402, 6.614393726402, 6.614393726402, 6.614393726402, 6.614393726402, 6.614393726402, 6.614393726402, 6.614393726402, 6.614393726402, 6.616184634020, 6.616184634020, 6.616184634020, 6.616184634020, 6.616184634020, 6.616184634020, 6.616184634020, 6.616184634020, 6.616184634020, 6.616184634020, 6.616184634020, 6.616184634020, 6.616184634020, 6.616184634020, 6.616184634020, 6.616184634020, 6.616184634020, 6.616184634020, 6.616184634020, 6.617982957425, 6.617982957425, 6.617982957425, 6.617982957425, 6.617982957425, 6.619788758288, 6.619788758288, 6.619788758288, 6.619788758288, 6.619788758288, 6.619788758288, 6.619788758288, 6.619788758288, 6.621602099052, 6.621602099052, 6.621602099052, 6.621602099052, 6.621602099052, 6.623423042944, 6.623423042944, 6.623423042944, 6.623423042944, 6.623423042944, 6.623423042944, 6.623423042944, 6.623423042944, 6.625251653990, 6.625251653990, 6.625251653990, 6.627087997030, 6.627087997030, 6.627087997030, 6.627087997030, 6.628932137728, 6.628932137728, 6.628932137728, 6.628932137728, 6.628932137728, 6.628932137728, 6.628932137728, 6.628932137728, 6.628932137728, 6.628932137728, 6.630784142590, 6.630784142590, 6.630784142590, 6.632644078974, 6.632644078974, 6.632644078974, 6.634512015109, 6.634512015109, 6.634512015109, 6.634512015109, 6.636388020108, 6.636388020108, 6.636388020108, 6.636388020108, 6.636388020108, 6.636388020108, 6.636388020108, 6.636388020108, 6.636388020108, 6.636388020108, 6.638272163982, 6.638272163982, 6.638272163982, 6.638272163982, 6.638272163982, 6.638272163982, 6.638272163982, 6.638272163982, 6.638272163982, 6.640164517660, 6.640164517660, 6.640164517660, 6.640164517660, 6.640164517660, 6.642065153000, 6.642065153000, 6.642065153000, 6.642065153000, 6.642065153000, 6.642065153000, 6.643974142807, 6.643974142807, 6.643974142807, 6.647817481889, 6.647817481889, 6.649751981666, 6.649751981666, 6.651695136952, 6.651695136952, 6.651695136952, 6.651695136952, 6.653647025549, 6.653647025549, 6.653647025549, 6.653647025549, 6.653647025549, 6.653647025549, 6.653647025549, 6.655607726315, 6.655607726315, 6.657577319178, 6.657577319178, 6.657577319178, 6.657577319178, 6.657577319178, 6.659555885160, 6.659555885160, 6.659555885160, 6.659555885160, 6.659555885160, 6.659555885160, 6.659555885160, 6.659555885160, 6.661543506395, 6.661543506395, 6.661543506395, 6.661543506395, 6.661543506395, 6.661543506395, 6.663540266151, 6.663540266151, 6.663540266151, 6.663540266151, 6.663540266151, 6.663540266151, 6.663540266151, 6.663540266151, 6.663540266151, 6.663540266151, 6.663540266151, 6.663540266151, 6.663540266151, 6.663540266151, 6.663540266151, 6.663540266151, 6.663540266151, 6.663540266151, 6.663540266151, 6.663540266151, 6.663540266151, 6.663540266151, 6.667561540084, 6.667561540084, 6.667561540084, 6.667561540084, 6.669586226651, 6.669586226651, 6.673664139071, 6.673664139071, 6.673664139071, 6.673664139071, 6.673664139071, 6.673664139071, 6.677780705266, 6.679853713889, 6.681936665037, 6.681936665037, 6.681936665037, 6.681936665037, 6.681936665037, 6.681936665037, 6.681936665037, 6.681936665037, 6.684029654543, 6.686132779631, 6.686132779631, 6.686132779631, 6.686132779631, 6.686132779631, 6.686132779631, 6.686132779631, 6.686132779631, 6.688246138944, 6.688246138944, 6.692503962087, 6.692503962087, 6.692503962087, 6.694648630553, 6.694648630553, 6.694648630553, 6.694648630553, 6.696803942580, 6.696803942580, 6.696803942580, 6.696803942580, 6.696803942580, 6.696803942580, 6.696803942580, 6.696803942580, 6.698970004336, 6.698970004336, 6.698970004336, 6.698970004336, 6.698970004336, 6.698970004336, 6.698970004336, 6.698970004336, 6.698970004336, 6.698970004336, 6.698970004336, 6.701146923590, 6.701146923590, 6.701146923590, 6.701146923590, 6.701146923590, 6.701146923590, 6.701146923590, 6.701146923590, 6.701146923590, 6.701146923590, 6.701146923590, 6.701146923590, 6.701146923590, 6.701146923590, 6.701146923590, 6.701146923590, 6.701146923590, 6.701146923590, 6.701146923590, 6.701146923590, 6.701146923590, 6.703334809738, 6.703334809738, 6.703334809738, 6.707743928644, 6.707743928644, 6.707743928644, 6.707743928644, 6.707743928644, 6.709965388637, 6.709965388637, 6.712198270070, 6.712198270070, 6.712198270070, 6.712198270070, 6.712198270070, 6.712198270070, 6.712198270070, 6.712198270070, 6.712198270070, 6.712198270070, 6.712198270070, 6.712198270070, 6.712198270070, 6.712198270070, 6.712198270070, 6.712198270070, 6.712198270070, 6.712198270070, 6.712198270070, 6.714442690992, 6.714442690992, 6.714442690992, 6.714442690992, 6.714442690992, 6.714442690992, 6.716698771296, 6.716698771296, 6.716698771296, 6.716698771296, 6.718966632752, 6.718966632752, 6.718966632752, 6.718966632752, 6.718966632752, 6.718966632752, 6.718966632752, 6.718966632752, 6.718966632752, 6.718966632752, 6.718966632752, 6.721246399047, 6.721246399047, 6.721246399047, 6.721246399047, 6.723538195827, 6.723538195827, 6.723538195827, 6.723538195827, 6.723538195827, 6.725842150736, 6.725842150736, 6.725842150736, 6.728158393464, 6.728158393464, 6.730487055782, 6.730487055782, 6.730487055782, 6.730487055782, 6.730487055782, 6.730487055782, 6.730487055782, 6.730487055782, 6.730487055782, 6.732828271597, 6.732828271597, 6.732828271597, 6.732828271597, 6.732828271597, 6.732828271597, 6.735182176990, 6.735182176990, 6.735182176990, 6.735182176990, 6.735182176990, 6.737548910270, 6.737548910270, 6.739928612015, 6.739928612015, 6.739928612015, 6.744727494897, 6.744727494897, 6.744727494897, 6.744727494897, 6.744727494897, 6.744727494897, 6.747146969020, 6.747146969020, 6.747146969020, 6.747146969020, 6.747146969020, 6.749579997691, 6.749579997691, 6.749579997691, 6.752026733638, 6.752026733638, 6.752026733638, 6.752026733638, 6.752026733638, 6.752026733638, 6.752026733638, 6.752026733638, 6.752026733638, 6.754487332186, 6.756961951314, 6.756961951314, 6.756961951314, 6.756961951314, 6.759450751717, 6.759450751717, 6.759450751717, 6.759450751717, 6.759450751717, 6.759450751717, 6.759450751717, 6.761953896871, 6.761953896871, 6.761953896871, 6.761953896871, 6.761953896871, 6.761953896871, 6.761953896871, 6.764471553092, 6.767003889608, 6.767003889608, 6.767003889608, 6.767003889608, 6.767003889608, 6.769551078622, 6.769551078622, 6.769551078622, 6.769551078622, 6.769551078622, 6.769551078622, 6.769551078622, 6.769551078622, 6.769551078622, 6.769551078622, 6.769551078622, 6.769551078622, 6.769551078622, 6.769551078622, 6.769551078622, 6.772113295386, 6.772113295386, 6.772113295386, 6.774690718274, 6.774690718274, 6.774690718274, 6.774690718274, 6.777283528852, 6.777283528852, 6.777283528852, 6.777283528852, 6.777283528852, 6.777283528852, 6.777283528852, 6.777283528852, 6.777283528852, 6.777283528852, 6.777283528852, 6.777283528852, 6.777283528852, 6.777283528852, 6.777283528852, 6.779891911960, 6.782516055786, 6.785156151952, 6.785156151952, 6.785156151952, 6.785156151952, 6.785156151952, 6.785156151952, 6.785156151952, 6.785156151952, 6.785156151952, 6.785156151952, 6.785156151952, 6.787812395596, 6.787812395596, 6.787812395596, 6.787812395596, 6.787812395596, 6.787812395596, 6.787812395596, 6.787812395596, 6.790484985457, 6.790484985457, 6.790484985457, 6.790484985457, 6.790484985457, 6.790484985457, 6.790484985457, 6.790484985457, 6.790484985457, 6.793174123968, 6.793174123968, 6.793174123968, 6.793174123968, 6.793174123968, 6.793174123968, 6.793174123968, 6.793174123968, 6.795880017344, 6.795880017344, 6.795880017344, 6.795880017344, 6.795880017344, 6.795880017344, 6.795880017344, 6.795880017344, 6.795880017344, 6.795880017344, 6.795880017344, 6.795880017344, 6.795880017344, 6.795880017344, 6.795880017344, 6.795880017344, 6.795880017344, 6.798602875680, 6.798602875680, 6.801342913046, 6.801342913046, 6.801342913046, 6.801342913046, 6.801342913046, 6.801342913046, 6.801342913046, 6.801342913046, 6.801342913046, 6.801342913046, 6.801342913046, 6.801342913046, 6.804100347591, 6.804100347591, 6.806875401646, 6.806875401646, 6.806875401646, 6.806875401646, 6.806875401646, 6.806875401646, 6.806875401646, 6.806875401646, 6.806875401646, 6.806875401646, 6.809668301830, 6.809668301830, 6.809668301830, 6.809668301830, 6.809668301830, 6.809668301830, 6.812479279164, 6.812479279164, 6.812479279164, 6.812479279164, 6.812479279164, 6.812479279164, 6.812479279164, 6.812479279164, 6.812479279164, 6.812479279164, 6.815308569182, 6.815308569182, 6.815308569182, 6.815308569182, 6.818156412055, 6.818156412055, 6.821023052707, 6.821023052707, 6.821023052707, 6.821023052707, 6.826813731588, 6.826813731588, 6.826813731588, 6.826813731588, 6.826813731588, 6.826813731588, 6.826813731588, 6.826813731588, 6.826813731588, 6.829738284605, 6.829738284605, 6.829738284605, 6.829738284605, 6.829738284605, 6.829738284605, 6.829738284605, 6.832682665252, 6.835647144216, 6.835647144216, 6.838631997765, 6.838631997765, 6.838631997765, 6.838631997765, 6.838631997765, 6.838631997765, 6.838631997765, 6.838631997765, 6.838631997765, 6.838631997765, 6.838631997765, 6.838631997765, 6.838631997765, 6.838631997765, 6.838631997765, 6.838631997765, 6.838631997765, 6.838631997765, 6.841637507905, 6.841637507905, 6.841637507905, 6.841637507905, 6.841637507905, 6.841637507905, 6.841637507905, 6.841637507905, 6.841637507905, 6.841637507905, 6.841637507905, 6.844663962535, 6.847711655617, 6.847711655617, 6.847711655617, 6.847711655617, 6.847711655617, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.850780887345, 6.853871964322, 6.853871964322, 6.853871964322, 6.853871964322, 6.853871964322, 6.853871964322, 6.853871964322, 6.853871964322, 6.853871964322, 6.856985199746, 6.856985199746, 6.856985199746, 6.856985199746, 6.856985199746, 6.856985199746, 6.856985199746, 6.856985199746, 6.856985199746, 6.856985199746, 6.856985199746, 6.856985199746, 6.856985199746, 6.856985199746, 6.856985199746, 6.856985199746, 6.856985199746, 6.856985199746, 6.856985199746, 6.856985199746, 6.856985199746, 6.856985199746, 6.856985199746, 6.860120913599, 6.860120913599, 6.860120913599, 6.860120913599, 6.860120913599, 6.860120913599, 6.860120913599, 6.860120913599, 6.860120913599, 6.863279432844, 6.863279432844, 6.863279432844, 6.863279432844, 6.863279432844, 6.866461091630, 6.866461091630, 6.866461091630, 6.866461091630, 6.866461091630, 6.869666231505, 6.869666231505, 6.869666231505, 6.869666231505, 6.869666231505, 6.869666231505, 6.869666231505, 6.869666231505, 6.869666231505, 6.869666231505, 6.869666231505, 6.869666231505, 6.869666231505, 6.869666231505, 6.869666231505, 6.869666231505, 6.869666231505, 6.869666231505, 6.869666231505, 6.869666231505, 6.869666231505, 6.869666231505, 6.869666231505, 6.869666231505, 6.869666231505, 6.869666231505, 6.869666231505, 6.869666231505, 6.869666231505, 6.869666231505, 6.869666231505, 6.869666231505, 6.869666231505, 6.872895201635, 6.872895201635, 6.872895201635, 6.872895201635, 6.872895201635, 6.872895201635, 6.872895201635, 6.872895201635, 6.872895201635, 6.876148359033, 6.876148359033, 6.879426068794, 6.879426068794, 6.882728704344, 6.882728704344, 6.886056647693, 6.889410289701, 6.889410289701, 6.889410289701, 6.889410289701, 6.892790030352, 6.892790030352, 6.892790030352, 6.892790030352, 6.892790030352, 6.892790030352, 6.892790030352, 6.892790030352, 6.892790030352, 6.892790030352, 6.892790030352, 6.892790030352, 6.892790030352, 6.892790030352, 6.892790030352, 6.892790030352, 6.892790030352, 6.892790030352, 6.892790030352, 6.899629454882, 6.899629454882, 6.899629454882, 6.899629454882, 6.899629454882, 6.899629454882, 6.899629454882, 6.899629454882, 6.899629454882, 6.899629454882, 6.899629454882, 6.899629454882, 6.899629454882, 6.899629454882, 6.899629454882, 6.899629454882, 6.899629454882, 6.903089986992, 6.903089986992, 6.906578314838, 6.906578314838, 6.906578314838, 6.906578314838, 6.906578314838, 6.906578314838, 6.906578314838, 6.910094888561, 6.910094888561, 6.910094888561, 6.910094888561, 6.910094888561, 6.910094888561, 6.910094888561, 6.910094888561, 6.913640169325, 6.913640169325, 6.913640169325, 6.913640169325, 6.917214629684, 6.917214629684, 6.917214629684, 6.917214629684, 6.920818753952, 6.920818753952, 6.920818753952, 6.920818753952, 6.924453038607, 6.924453038607, 6.924453038607, 6.928117992694, 6.928117992694, 6.928117992694, 6.928117992694, 6.928117992694, 6.928117992694, 6.928117992694, 6.928117992694, 6.928117992694, 6.928117992694, 6.928117992694, 6.928117992694, 6.928117992694, 6.928117992694, 6.928117992694, 6.928117992694, 6.928117992694, 6.928117992694, 6.928117992694, 6.928117992694, 6.931814138254, 6.931814138254, 6.931814138254, 6.935542010773, 6.935542010773, 6.943095148664, 6.943095148664, 6.943095148664, 6.943095148664, 6.943095148664, 6.943095148664, 6.943095148664, 6.943095148664, 6.943095148664, 6.943095148664, 6.943095148664, 6.943095148664, 6.943095148664, 6.943095148664, 6.943095148664, 6.943095148664, 6.943095148664, 6.943095148664, 6.943095148664, 6.943095148664, 6.943095148664, 6.943095148664, 6.943095148664, 6.943095148664, 6.943095148664, 6.943095148664, 6.946921556517, 6.946921556517, 6.946921556517, 6.946921556517, 6.950781977330, 6.950781977330, 6.950781977330, 6.954677021213, 6.954677021213, 6.954677021213, 6.954677021213, 6.954677021213, 6.954677021213, 6.954677021213, 6.954677021213, 6.954677021213, 6.954677021213, 6.954677021213, 6.954677021213, 6.954677021213, 6.954677021213, 6.954677021213, 6.954677021213, 6.954677021213, 6.954677021213, 6.954677021213, 6.954677021213, 6.954677021213, 6.954677021213, 6.954677021213, 6.954677021213, 6.954677021213, 6.954677021213, 6.954677021213, 6.954677021213, 6.954677021213, 6.954677021213, 6.954677021213, 6.958607314842, 6.958607314842, 6.958607314842, 6.958607314842, 6.958607314842, 6.958607314842, 6.962573502059, 6.962573502059, 6.962573502059, 6.962573502059, 6.962573502059, 6.962573502059, 6.966576244513, 6.966576244513, 6.966576244513, 6.966576244513, 6.966576244513, 6.966576244513, 6.966576244513, 6.966576244513, 6.966576244513, 6.966576244513, 6.966576244513, 6.966576244513, 6.966576244513, 6.966576244513, 6.966576244513, 6.966576244513, 6.970616222315, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.974694134735, 6.978810700930, 6.978810700930, 6.978810700930, 6.978810700930, 6.978810700930, 6.978810700930, 6.978810700930, 6.978810700930, 6.978810700930, 6.978810700930, 6.978810700930, 6.978810700930, 6.978810700930, 6.978810700930, 6.978810700930, 6.978810700930, 6.978810700930, 6.978810700930, 6.978810700930, 6.978810700930, 6.982966660701, 6.982966660701, 6.982966660701, 6.982966660701, 6.982966660701, 6.982966660701, 6.982966660701, 6.982966660701, 6.982966660701, 6.982966660701, 6.982966660701, 6.982966660701, 6.982966660701, 6.982966660701, 6.982966660701, 6.982966660701, 6.982966660701, 6.982966660701, 6.982966660701, 6.982966660701, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.987162775295, 6.991399828238, 6.991399828238, 6.991399828238, 6.991399828238, 6.991399828238, 6.991399828238, 6.991399828238, 6.991399828238, 6.991399828238, 6.991399828238, 6.991399828238, 6.991399828238, 6.991399828238, 6.991399828238, 6.995678626217, 6.995678626217, 7.000000000000, 7.000000000000, 7.004364805402, 7.004364805402, 7.004364805402, 7.008773924308, 7.017728766960, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.022276394711, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.026872146400, 7.036212172654, 7.036212172654, 7.036212172654, 7.036212172654, 7.036212172654, 7.036212172654, 7.036212172654, 7.036212172654, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.040958607679, 7.045757490561, 7.045757490561, 7.045757490561, 7.045757490561, 7.045757490561, 7.045757490561, 7.045757490561, 7.045757490561, 7.045757490561, 7.045757490561, 7.045757490561, 7.045757490561, 7.050609993355, 7.050609993355, 7.055517327850, 7.055517327850, 7.060480747381, 7.060480747381, 7.060480747381, 7.060480747381, 7.060480747381, 7.060480747381, 7.060480747381, 7.060480747381, 7.060480747381, 7.060480747381, 7.060480747381, 7.060480747381, 7.060480747381, 7.060480747381, 7.060480747381, 7.060480747381, 7.060480747381, 7.060480747381, 7.060480747381, 7.060480747381, 7.065501548756, 7.065501548756, 7.065501548756, 7.065501548756, 7.065501548756, 7.065501548756, 7.065501548756, 7.065501548756, 7.065501548756, 7.065501548756, 7.065501548756, 7.065501548756, 7.065501548756, 7.065501548756, 7.065501548756, 7.065501548756, 7.065501548756, 7.065501548756, 7.065501548756, 7.065501548756, 7.065501548756, 7.065501548756, 7.065501548756, 7.065501548756, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.070581074286, 7.075720713938, 7.075720713938, 7.075720713938, 7.075720713938, 7.075720713938, 7.075720713938, 7.075720713938, 7.075720713938, 7.075720713938, 7.080921907624, 7.080921907624, 7.086186147616, 7.086186147616, 7.086186147616, 7.086186147616, 7.086186147616, 7.086186147616, 7.086186147616, 7.086186147616, 7.086186147616, 7.086186147616, 7.086186147616, 7.086186147616, 7.086186147616, 7.086186147616, 7.086186147616, 7.086186147616, 7.086186147616, 7.086186147616, 7.086186147616, 7.086186147616, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.091514981121, 7.096910013008, 7.096910013008, 7.096910013008, 7.096910013008, 7.096910013008, 7.096910013008, 7.096910013008, 7.096910013008, 7.096910013008, 7.096910013008, 7.096910013008, 7.102372908710, 7.102372908710, 7.102372908710, 7.102372908710, 7.102372908710, 7.102372908710, 7.102372908710, 7.102372908710, 7.102372908710, 7.102372908710, 7.102372908710, 7.102372908710, 7.107905397310, 7.107905397310, 7.107905397310, 7.107905397310, 7.113509274828, 7.113509274828, 7.113509274828, 7.113509274828, 7.113509274828, 7.113509274828, 7.113509274828, 7.119186407719, 7.119186407719, 7.119186407719, 7.119186407719, 7.119186407719, 7.119186407719, 7.119186407719, 7.119186407719, 7.119186407719, 7.119186407719, 7.119186407719, 7.119186407719, 7.124938736608, 7.124938736608, 7.124938736608, 7.124938736608, 7.124938736608, 7.130768280269, 7.130768280269, 7.130768280269, 7.130768280269, 7.130768280269, 7.130768280269, 7.130768280269, 7.130768280269, 7.136677139880, 7.142667503569, 7.142667503569, 7.142667503569, 7.142667503569, 7.142667503569, 7.142667503569, 7.142667503569, 7.148741651281, 7.148741651281, 7.154901959986, 7.154901959986, 7.154901959986, 7.154901959986, 7.154901959986, 7.154901959986, 7.154901959986, 7.154901959986, 7.154901959986, 7.154901959986, 7.154901959986, 7.154901959986, 7.154901959986, 7.161150909263, 7.161150909263, 7.161150909263, 7.161150909263, 7.161150909263, 7.161150909263, 7.161150909263, 7.161150909263, 7.161150909263, 7.161150909263, 7.161150909263, 7.161150909263, 7.161150909263, 7.161150909263, 7.161150909263, 7.161150909263, 7.161150909263, 7.161150909263, 7.161150909263, 7.161150909263, 7.161150909263, 7.161150909263, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.167491087294, 7.173925197299, 7.173925197299, 7.173925197299, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.180456064458, 7.187086643357, 7.187086643357, 7.187086643357, 7.193820026016, 7.193820026016, 7.193820026016, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.200659450546, 7.207608310502, 7.207608310502, 7.207608310502, 7.207608310502, 7.207608310502, 7.207608310502, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.214670164989, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.221848749616, 7.229147988358, 7.229147988358, 7.229147988358, 7.229147988358, 7.229147988358, 7.229147988358, 7.229147988358, 7.229147988358, 7.229147988358, 7.229147988358, 7.229147988358, 7.229147988358, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.236572006437, 7.244125144328, 7.244125144328, 7.244125144328, 7.244125144328, 7.244125144328, 7.244125144328, 7.244125144328, 7.244125144328, 7.244125144328, 7.244125144328, 7.244125144328, 7.244125144328, 7.251811972994, 7.251811972994, 7.251811972994, 7.251811972994, 7.251811972994, 7.251811972994, 7.251811972994, 7.251811972994, 7.251811972994, 7.251811972994, 7.251811972994, 7.251811972994, 7.251811972994, 7.251811972994, 7.251811972994, 7.251811972994, 7.251811972994, 7.251811972994, 7.251811972994, 7.251811972994, 7.251811972994, 7.251811972994, 7.251811972994, 7.251811972994, 7.251811972994, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.259637310506, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.267606240177, 7.275724130399, 7.275724130399, 7.275724130399, 7.275724130399, 7.275724130399, 7.275724130399, 7.275724130399, 7.275724130399, 7.275724130399, 7.275724130399, 7.275724130399, 7.275724130399, 7.275724130399, 7.275724130399, 7.283996656365, 7.283996656365, 7.283996656365, 7.283996656365, 7.283996656365, 7.283996656365, 7.283996656365, 7.283996656365, 7.283996656365, 7.283996656365, 7.283996656365, 7.283996656365, 7.283996656365, 7.292429823902, 7.292429823902, 7.292429823902, 7.292429823902, 7.292429823902, 7.292429823902, 7.292429823902, 7.292429823902, 7.292429823902, 7.292429823902, 7.292429823902, 7.292429823902, 7.292429823902, 7.292429823902, 7.292429823902, 7.292429823902, 7.292429823902, 7.292429823902, 7.292429823902, 7.292429823902, 7.292429823902, 7.292429823902, 7.292429823902, 7.301029995664, 7.301029995664, 7.301029995664, 7.301029995664, 7.301029995664, 7.301029995664, 7.301029995664, 7.301029995664, 7.301029995664, 7.301029995664, 7.301029995664, 7.301029995664, 7.301029995664, 7.301029995664, 7.301029995664, 7.309803919972, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.318758762624, 7.327902142064, 7.327902142064, 7.327902142064, 7.327902142064, 7.327902142064, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.337242168318, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.346787486225, 7.356547323514, 7.356547323514, 7.356547323514, 7.356547323514, 7.356547323514, 7.356547323514, 7.356547323514, 7.356547323514, 7.356547323514, 7.366531544420, 7.366531544420, 7.366531544420, 7.366531544420, 7.366531544420, 7.366531544420, 7.366531544420, 7.366531544420, 7.366531544420, 7.366531544420, 7.376750709602, 7.376750709602, 7.376750709602, 7.376750709602, 7.376750709602, 7.376750709602, 7.376750709602, 7.376750709602, 7.376750709602, 7.376750709602, 7.376750709602, 7.376750709602, 7.376750709602, 7.376750709602, 7.376750709602, 7.376750709602, 7.376750709602, 7.387216143280, 7.397940008672, 7.397940008672, 7.397940008672, 7.397940008672, 7.397940008672, 7.397940008672, 7.397940008672, 7.397940008672, 7.397940008672, 7.397940008672, 7.397940008672, 7.397940008672, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.408935392974, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.420216403383, 7.431798275933, 7.443697499233, 7.443697499233, 7.443697499233, 7.443697499233, 7.443697499233, 7.443697499233, 7.443697499233, 7.443697499233, 7.443697499233, 7.443697499233, 7.443697499233, 7.443697499233, 7.443697499233, 7.443697499233, 7.443697499233, 7.443697499233, 7.443697499233, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.455931955650, 7.468521082958, 7.468521082958, 7.468521082958, 7.468521082958, 7.468521082958, 7.468521082958, 7.468521082958, 7.468521082958, 7.468521082958, 7.468521082958, 7.468521082958, 7.468521082958, 7.468521082958, 7.468521082958, 7.468521082958, 7.468521082958, 7.468521082958, 7.468521082958, 7.468521082958, 7.468521082958, 7.468521082958, 7.468521082958, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.481486060122, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.494850021680, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.508638306166, 7.522878745280, 7.522878745280, 7.522878745280, 7.522878745280, 7.522878745280, 7.522878745280, 7.522878745280, 7.522878745280, 7.522878745280, 7.522878745280, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.537602002101, 7.552841968658, 7.552841968658, 7.552841968658, 7.552841968658, 7.552841968658, 7.552841968658, 7.552841968658, 7.552841968658, 7.552841968658, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.568636235841, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.585026652029, 7.602059991328, 7.602059991328, 7.602059991328, 7.602059991328, 7.602059991328, 7.602059991328, 7.602059991328, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.619788758288, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.638272163982, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.657577319178, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.677780705266, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.698970004336, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.721246399047, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.744727494897, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.769551078622, 7.795880017344, 7.795880017344, 7.795880017344, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.823908740944, 7.853871964322, 7.853871964322, 7.853871964322, 7.853871964322, 7.853871964322, 7.853871964322, 7.853871964322, 7.853871964322, 7.853871964322, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.886056647693, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.920818753952, 7.958607314842, 7.958607314842, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.000000000000, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.045757490561, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.096910013008, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.154901959986, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.221848749616, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.301029995664, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.397940008672, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.522878745280, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 8.698970004336, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000, 9.000000000000};
    return neg_log10_eccdf;
  }
  else if( slide ==  15 && lambda == 25){
    return neg_log10_eccdf;
  }
  else if( slide ==  15 && lambda == 5){
    return neg_log10_eccdf;
  }
  else if( slide ==  10 && lambda == 100){
    return neg_log10_eccdf;
  }
  else if( slide ==  10 && lambda == 50){
    return neg_log10_eccdf;
  }
  else if( slide ==  10 && lambda == 25){
    return neg_log10_eccdf;
  }
  else if( slide ==  5 && lambda == 100){
    return neg_log10_eccdf;
  }
  else if( slide ==  5 && lambda == 50){
    return neg_log10_eccdf;
  }
  else if( slide ==  5 && lambda == 25){
    return neg_log10_eccdf;
  }
  else{
    std::cerr << "ERROR: unknown parameters w = " << slide << " and lambda = " << lambda << std::endl;
    exit(-1);
  }
}