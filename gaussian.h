#ifndef GAUSSIAN_H
#define GAUSSIAN_H

#include "helpinghand.h"

class Gaussian {
  
 public:
  Gaussian();
  Gaussian( const double &mu_x_par,
  	    const double &mu_y_par,
  	    const double &sigma_xx_par,
	    const double &sigma_yy_par,
	    const double &phi_par);
  Gaussian( const std::vector< double > &mu_x_par,
	    const std::vector< double > &mu_y_par,
	    const std::vector< double > &sigma_xx_par,
	    const std::vector< double > &sigma_yy_par,
	    const std::vector< double > &phi_par,
	    const std::vector< double > &ratio_par,
	    const unsigned &Nsources);
  
  void AddGaussian( const double &mu_x_par,
		    const double &mu_y_par,
		    const double &sigma_xx_par,
		    const double &sigma_yy_par,
		    const double &phi_par,
		    const double &ratio);
  
  double setratio(const std::vector< double > &ratio_par);
  
  unsigned number_of_sources() const;
  
  std::vector< double > all_mu_x() const;
  std::vector< double > all_mu_y() const;
  std::vector< double > all_sigma_xx() const;
  std::vector< double > all_sigma_yy() const;
  std::vector< double > all_phi() const;
  std::vector< double > all_ratio() const;
  
  double mu_x(const unsigned &numbrofsource) const;
  double mu_y(const unsigned &numbrofsource) const;
  double sigma_xx(const unsigned &numbrofsource) const;
  double sigma_yy(const unsigned &numbrofsource) const;
  double phi(const unsigned &numbrofsource) const;
  double ratio(const unsigned &numbrofsource) const;
  
 private:
  std::vector< double > mu_x_;
  std::vector< double > mu_y_;
  std::vector< double > sigma_xx_;
  std::vector< double > sigma_yy_;
  std::vector< double > phi_;
  std::vector< double > ratio_;
  unsigned Nsources_;
  
  void CheckNumberOfSources() const;
};

#endif

