#ifndef PHASE_SPACE_H
#define PHASE_SPACE_H

#include "lookup_phase_space.h"
#include "gaussian.h"
#include "randomnumbers.h"
#include "stateclassifier.h"
#include <gsl/gsl_histogram.h>
#include <ctime>
#include <sys/resource.h>

// --------------------
// Class of determining phase space
class PhaseSpace {
  
 public:
  PhaseSpace(const unsigned &N, const unsigned &A,
	     const double &nmbr_of_microstates, const unsigned &NmbrMacrostates,
	     const double &f_initial, const double &f_mod, const double &f_epsilon, 
	     const double &hist_dev,
	     const std::string &prefix_if, const std::string &prefix_of,
	     const bool &verbose,
	     const unsigned &seed, const unsigned &seedwght);
  void reset();
  
  void CHANGE_STATE();
  void CHANGE_STATE_wo_hist();
  
  void EVOLVE(const unsigned &slowdown);
  void EVOLVE_classic_20times();
  void EVOLVE_classic();
  void EVOLVE_analytic();
  
  unsigned evaluate_NmbrMacrostates() const;
  
  void NORMALIZE();
  void NORMALIZE_with_total_sum();
  
  void display_sample() const;
  
  void output_DoS(const bool &numerical);
  void output_ammendment_DoS();
  void output_ammendment_log_DoS();
  void output_ammendment_log_DoS_all();
  
  void evaluate_single_deviations(const unsigned &start_seed, const unsigned &start_seedwght,
				  const unsigned &end_seed, const unsigned &end_seedwght,
				  const unsigned &bins=61, const double &max_dev=0.25);
  void evaluate_all_maximum_deviations(const unsigned &start_seed, const unsigned &start_seedwght,
				       const unsigned &end_seed, const unsigned &end_seedwght,
				       const unsigned &end_actual_seed, const unsigned &end_actual_seedwght);
  
 private:
  unsigned n_bin(const int &P, const int &C) const;
  void change_macrostate(const unsigned &xb_, const unsigned &yb_, const unsigned &xw_, const unsigned &yw_);
  
  const unsigned N, N2;
  const unsigned A, W;
  
  const double nmbr_of_microstates;
  const unsigned NmbrMacrostates;
  const double timescale;
  
  const double f_initial;
  const double f_mod;
  double f;
  const double f_epsilon, hist_dev;
  
  BinField<bool> sample;
  std::vector<unsigned> black_pixels_x, black_pixels_y, white_pixels_x, white_pixels_y;
  unsigned xb, yb, xw, yw;
  
  int P0, C0, bin0, P1, C1, bin1;
  
  const unsigned n_bins_P, n_bins_C, abs_min_C;
  const int max_P, max_C, min_C;
  std::vector<double> log_DoS;
  std::vector< std::vector<double> > log_DoS_all;
  std::vector<double> DoS;
  std::vector<unsigned> hist;
  
  unsigned non_zero_hist, hist_min;
  double hist_mean;
  double total_count;
  
  const std::string prefix_if;
  const std::string prefix_of;
  const bool verbose;
  const unsigned seed;
  const unsigned seedwght;
  const unsigned original_seed;
  const unsigned original_seedwght;
};
// --------------------

#endif

