#include "gaussian.h"

Gaussian::Gaussian() :
  mu_x_(std::vector<double>(0,0)), mu_y_(std::vector<double>(0,0)), sigma_xx_(std::vector<double>(0,0)), sigma_yy_(std::vector<double>(0,0)), phi_(std::vector<double>(0,0)), ratio_(std::vector<double>(0,0)), Nsources_(0)
{
  
}

Gaussian::Gaussian( const double &mu_x_par,
		    const double &mu_y_par,
		    const double &sigma_xx_par,
		    const double &sigma_yy_par,
		    const double &phi_par) :
  mu_x_(std::vector<double>(1,mu_x_par)), mu_y_(std::vector<double>(1,mu_y_par)), sigma_xx_(std::vector<double>(1,sigma_xx_par)), sigma_yy_(std::vector<double>(1,sigma_yy_par)), phi_(std::vector<double>(1,phi_par)), ratio_(std::vector<double>(1,1)), Nsources_(1)
{
  
}

Gaussian::Gaussian( const std::vector< double > &mu_x_par,
		    const std::vector< double > &mu_y_par,
		    const std::vector< double > &sigma_xx_par,
		    const std::vector< double > &sigma_yy_par,
		    const std::vector< double > &phi_par,
		    const std::vector< double > &ratio_par,
		    const unsigned &Nsources) :
  mu_x_(mu_x_par), mu_y_(mu_y_par), sigma_xx_(sigma_xx_par), sigma_yy_(sigma_yy_par), phi_(phi_par), ratio_ (ratio_par), Nsources_(Nsources)
{
  // Normalization constant of the vector ratio_par
  double summation = sum_abs<double>(ratio_par);
  
  // Define ratio
  std::vector< double > ratio_(ratio_par.size(),1.0/ratio_par.size());
  if(fabs(summation-1) < 1e-6)
    {
      for(unsigned i = 0; i < ratio_par.size(); i++)
	ratio_.at(i) = fabs(ratio_par.at(i));
    }
  else if(summation == 0)
    {

    }
  else
    {
      std::cerr << "WARNING: Normalization of ratio vector is not 1 but: " << summation << ";" << std::endl;
      std::cerr << "Renormalization: ratio = [";
      for(unsigned i = 0; i < ratio_par.size(); i++)
	ratio_.at(i) = fabs(ratio_par.at(i))/summation;      
      for(unsigned i = 0; i < ratio_.size(); i++)
	{
	  std::cerr << ratio_.at(i);
	  if(i<(ratio_.size()-1))
	    std::cerr <<  "; ";
	}
      std::cerr << "]" << std::endl;
      // THERE IS A PROBLEM
      // Most likely summation == 1 is not true!
      
    }
  
  // Check the number of sources;
  CheckNumberOfSources();
}

void Gaussian::AddGaussian( const double &mu_x_par,
			    const double &mu_y_par,
			    const double &sigma_xx_par,
			    const double &sigma_yy_par,
			    const double &phi_par,
			    const double &ratio_par)
{
  mu_x_.push_back(mu_x_par);
  mu_y_.push_back(mu_y_par);
  sigma_xx_.push_back(sigma_xx_par);
  sigma_yy_.push_back(sigma_yy_par);
  phi_.push_back(phi_par);
  for(unsigned i = 0; i < Nsources_; i++)
    ratio_.at(i) *= (1-ratio_par);
  ratio_.push_back(ratio_par);
  Nsources_++;  
}

double Gaussian::setratio(const std::vector< double > &ratio_par)
{
  ratio_ = ratio_par;
  return sum(ratio_);  
}

unsigned Gaussian::number_of_sources() const 
{
  return Nsources_;
}

std::vector< double > Gaussian::all_mu_x() const
{
  return mu_x_;
}
std::vector< double > Gaussian::all_mu_y() const
{
  return mu_y_;
}
std::vector< double > Gaussian::all_sigma_xx() const
{
  return sigma_xx_;
}
std::vector< double > Gaussian::all_sigma_yy() const
{
  return sigma_yy_;
}
std::vector< double > Gaussian::all_phi() const
{
  return phi_;
}
std::vector< double > Gaussian::all_ratio() const
{
  return ratio_;
}

double Gaussian::mu_x(const unsigned &numbrofsource) const
{
  return mu_x_.at(numbrofsource);
}
double Gaussian::mu_y(const unsigned &numbrofsource) const
{
  return mu_y_.at(numbrofsource);
}
double Gaussian::sigma_xx(const unsigned &numbrofsource) const
{
  return sigma_xx_.at(numbrofsource);
}
double Gaussian::sigma_yy(const unsigned &numbrofsource) const
{
  return sigma_yy_.at(numbrofsource);
}
double Gaussian::phi(const unsigned &numbrofsource) const
{
  return phi_.at(numbrofsource);
}
double Gaussian::ratio(const unsigned &numbrofsource) const
{
  return ratio_.at(numbrofsource);
}

void Gaussian::CheckNumberOfSources() const
{
  if(!(Nsources_==mu_x_.size()))
    {
      std::cerr << "Error in a Gaussian: Number of x-coordinates for " 
		<< " the centers of the distributions is " <<  mu_x_.size()
		<< " not " << Nsources_
		<< " as demanded!" << std::endl;
      exit(-1);
    }
  
  if(!(Nsources_==mu_y_.size()))
    {
      std::cerr << "Error in a Gaussian: Number of y-coordinates for " 
		<< " the centers of the distributions is " <<  mu_y_.size()
		<< " not " << Nsources_
		<< " as demanded!" << std::endl;
      exit(-1);
    }

  if(!(Nsources_==sigma_xx_.size()))
    {
      std::cerr << "Error in a Gaussian: Number of values sigma_xx for" 
		<< " the standard deviations of the distributions is " <<  sigma_xx_.size()
		<< " not " << Nsources_
		<< " as demanded!" << std::endl;
      exit(-1);
    }

  if(!(Nsources_==sigma_yy_.size()))
    {
      std::cerr << "Error in a Gaussian: Number of values sigma_yy for" 
		<< " the standard deviations of the distributions is " <<  sigma_yy_.size()
		<< " not " << Nsources_
		<< " as demanded!" << std::endl;
      exit(-1);
    }
  
  if(!(Nsources_==phi_.size()))
    {
      std::cerr << "Error in a Gaussian: Number of angles for " 
		<< " the orientation of the distributions is " <<  phi_.size()
		<< " not " << Nsources_
		<< " as demanded!" << std::endl;
      exit(-1);
    }
  
  if(!(Nsources_==ratio_.size()))
    {
      std::cerr << "Error in a Gaussian: Number of ratios of " 
		<< " the distributions is " <<  ratio_.size()
		<< " not " << Nsources_
		<< " as demanded!" << std::endl;
      exit(-1);
    }
}

