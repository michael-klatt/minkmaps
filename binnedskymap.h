#ifndef BINNED_SKY_MAP_H
#define BINNED_SKY_MAP_H

#include "gaussian.h"
#include "minkowski.h"
#include "randomnumbers.h"
#include "stateclassifier.h"
#include "lookup_phase_space.h"

#include <gsl/gsl_integration.h>
#include <boost/range/algorithm/sort.hpp>

// --------------------
// Global variables
// Accuracy of integration
const double demand_abs_error_ = 0;
const double demand_rel_error_ = 1e-7;

// External functions: Calculate intensity
// params points to std::vector<double>; content = {double y_low; double y_up}
double margin_gauss_pdf(double x, void * params);

// params points to double; content = x 
double gauss_pdf(double y, void * params);
// --------------------

// --------------------
// Class of a binned sky map
class BinnedSkyMap {
  
 public:
  BinnedSkyMap(const unsigned &Nx, const unsigned &Ny,
	       const double &xmin, const double &xmax,
	       const double &ymin, const double &ymax,
	       const std::string &prefix_if, const std::string &prefix_of,
	       const unsigned &seed, const unsigned &seedsub);
  BinnedSkyMap(const unsigned &N,
	       const double &xmin, const double &xmax,
	       const double &ymin, const double &ymax,
	       const std::string &prefix_if, const std::string &prefix_of,
	       const unsigned &seed, const unsigned &seedsub);
  BinnedSkyMap(const std::string &countratefilename,
	       const double &xmin, const double &xmax,
	       const double &ymin, const double &ymax,
	       const std::string &prefix_if, const std::string &prefix_of,
	       const unsigned &seed, const unsigned &seedsub);
  /*
  BinnedSkyMap(const std::string &eventlistfilename,
	       const double &xmin, const double &xmax,
	       const double &ymin, const double &ymax,
	       const std::string &prefix_if, const std::string &prefix_of,
	       const unsigned &seed, const unsigned &seedsub);
  */
  // --------------------
  // Reset
  void reset();
  // --------------------
  
  // --------------------
  // Manipulate intensity
  void scale_intensity(const double &factor);
  
  void add_bck_intensity(const double &k_bck);
  
  void add_intensity_to_single_bin(const double &lambda, const unsigned &xi, const unsigned &yi);
  
  void add_src_intensity(const double &k_src, const Gaussian &sources, const bool &verbose=true);
  
  // Calculate source intensity
  // Auxiliary function and doubles
  Gaussian sources_;                   // Cache source parameters
  double single_y_low_, single_y_up_;  // Cache integration boundaries
  double rx_;                          // Cache position of integration 
  
  double Lambda(const double &k_src, const double &x_low, const double &x_up, const double &y_low, const double &y_up);
  // --------------------
  
  // --------------------
  // Simulation of a countrate-sample with respect to the given intensity
  void add_counts_to_single_bin(const unsigned &k, const unsigned &xi, const unsigned &yi);
  
  void simulate_countrate();
  
  void reduction_because_of_detector_acceptance();
  
  void global_reduction(const double &f);
  
  void add_global_Poisson(const double &lambda);
  
  void global_detector_correction(const double &target_f, const double &lambda);
  
  std::vector< BinField< unsigned > > local_detector_correction(const double &lambda, const unsigned &slide) const;
  // --------------------
  
  // -------------------------------
  // Calling and print out functions
  unsigned call_Nx() const;
  unsigned call_Ny() const;
  double call_xlow(const unsigned &xi, const unsigned &yi) const;
  double call_xlow(const unsigned &ni) const;
  double call_ylow(const unsigned &xi, const unsigned &yi) const;
  double call_ylow(const unsigned &ni) const;
  double call_intensity(const unsigned &xi, const unsigned &yi) const;
  double call_intensity(const unsigned &ni) const;
  unsigned call_countrate(const unsigned &xi, const unsigned &yi) const;
  unsigned call_countrate(const unsigned &ni) const;
  BinField<unsigned> count_map() const;
  
  double call_xbinlength() const;
  double call_ybinlength() const;
  BinField<double> call_xlow() const;
  BinField<double> call_ylow() const;
  
  void xout_xlow(const unsigned &xlow, const unsigned &xup, const unsigned &ylow, const unsigned &yup) const;
  void xout_ylow(const unsigned &xlow, const unsigned &xup, const unsigned &ylow, const unsigned &yup) const;
  void xout_intensity(const unsigned &xlow, const unsigned &xup, const unsigned &ylow, const unsigned &yup) const;
  void xout_countrate(const unsigned &xlow, const unsigned &xup, const unsigned &ylow, const unsigned &yup) const;
  
  void fout_xlow(const std::string &filename) const;
  void fout_ylow(const std::string &filename) const;
  void fout_intensity(const std::string &filename) const;
  void fout_countrate(const std::string &filename) const;
  
  void matrix_out_intensity(const std::string &filename = "a_Intensity_Matrix") const;
  void matrix_out_countrate(const std::string &filename = "b_Countrate_Matrix") const;
  
  void plot_intensity(const std::string &filename = "a_Intensity") const;
  void plot_countrate(const std::string &filename = "b_Countrate") const;
  // -------------------------------
  
  // -------------------------------
  // Subtract Sources
  void DefineSource(const double &k_src, const Gaussian &sources);
  
  void SubtractSource(const double &lambda);
  // -------------------------------
  
  // -------------------------------
  // Global Analysis
  double global_deviation_strength_A_pix_wbc(const double &lambda, const std::string &filename="global_deviation_strength") const;
  // -------------------------------
  
  // -------------------------------
  // Local Analysis
  BinField<double> significance_map(const double &lambda, const unsigned int &slide=15, const bool &verbose=true,
				    const bool &matrix_out=true, const bool &gnuplot_out=true,
				    const std::string &filename="significance_map") const;
  
  double max_of_significance_map(const double &lambda, const unsigned int &slide=15) const;
  
  BinField<double> significance_map_chessboard(const double &lambda, const unsigned int &slide=15, const bool &verbose=true,
					       const bool &matrix_out=true,
					       const std::string &filename="significance_map") const;

  BinField<double> Minkowski_sky_map_APC_pix_wbc(const double &lambda, const unsigned int &slide=15,
						 const bool &verbose=true, const bool &matrix_out=true, const bool &gnuplot_out=true,
						 const std::string &filename="Minkowski_sky_map_APC_pix_wbc") const;
  
  BinField<double> Minkowski_sky_map_APC_pix_wbc_chessboard(const double &lambda, const unsigned int &slide=15,
							    const bool &verbose=true, const bool &matrix_out=true,
							    const std::string &filename="Minkowski_sky_map_APC_pix_wbc") const;
  
  BinField<double> Minkowski_sky_map_APC_locally_different_excursion_sets_pix_wbc(const double &lambda, const unsigned int &slide,
										  const std::vector< BinField< unsigned > > &countrate_locally_corrected,
										  const bool &verbose=true, const bool &matrix_out=true, const bool &gnuplot_out=true,
										  const std::string &filename="Minkowski_sky_map_locally_corrected_APC_pix_wbc") const;
  
  BinField<double> Minkowski_sky_map_sum_D_APC_pix_wbc(const double &lambda, const unsigned int &slide,
						       const bool &verbose=true, const bool &matrix_out=true, const bool &gnuplot_out=true,
						       const std::string &filename="Minkowski_sky_map_sum_D_APC_pix_wbc") const;
  
  BinField<double> Minkowski_sky_map_sum_D_APC_pix_wbc_chessboard(const double &lambda, const unsigned int &slide,
								  const bool &verbose=true, const bool &matrix_out=true,
								  const std::string &filename="Minkowski_sky_map_sum_D_APC_pix_wbc") const;

  BinField<double> Minkowski_sky_map_sum_D_small_lambda_APC_pix_wbc(const double &lambda, const unsigned int &slide,
								    const bool &verbose=true, const bool &matrix_out=true, const bool &gnuplot_out=true,
								    const std::string &filename="Minkowski_sky_map_sum_D_APC_pix_wbc") const;
  
  BinField<double> Minkowski_sky_map_sum_D_small_lambda_APC_pix_wbc_chessboard(const double &lambda, const unsigned int &slide,
									       const bool &verbose=true, const bool &matrix_out=true,
									       const std::string &filename="Minkowski_sky_map_sum_D_APC_pix_wbc") const;
  
  BinField<double> Minkowski_sky_map_A_pix_wbc(const double &lambda, const unsigned int &slide=15, const bool &verbose=true,
					       const std::string &filename="Minkowski_sky_map_A_pix_wbc") const;
  
  BinField<double> Minkowski_sky_map_A_pix_wbc_chessboard(const double &lambda, const unsigned int &slide=15, const bool &verbose=true,
							  const std::string &filename="Minkowski_sky_map_A_pix_wbc") const;
  // -------------------------------
  
 private:
  double * load_neg_log10_eccdf(const unsigned &slide, const double &lambda) const;
  
  const unsigned Nx_;
  const unsigned Ny_;
  
  const double xmin_;
  const double xmax_;
  const double ymin_;
  const double ymax_;
  const double xbinlength_;
  const double ybinlength_;
  
  BinField< double > xlow_;
  BinField< double > ylow_;
  
  BinField< double > intensity_;
  BinField< unsigned > countrate_;
  BinField< double > srcsub_intensity_;
  
  const std::string prefix_if_;
  const std::string prefix_of_;
  const unsigned seed_;
  const unsigned seedsub_;

  const unsigned trial;
};
// --------------------

#endif

