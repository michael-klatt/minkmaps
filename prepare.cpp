// Preparations for Morphometric Minowski Numerical Analysis 
#include "binnedskymap.h"
#include "phasespace.h"
#include <boost/math/special_functions/binomial.hpp>

static std::string config_file = "minkmap.conf";
static std::string prefix_if = "input/";
static std::string prefix_of = "output/";
static std::string countratefilename = "no_filename";
static unsigned Nx = 100;
static unsigned Ny = Nx;
static double xmin = -1;
static double xmax = 1;
static double ymin = xmin;
static double ymax = xmax;
static double lambda = 5;
static double k_src = 0;
static std::string functionals = "APC";
static unsigned slide = 5;
static unsigned A_blackpix = 0;
static bool mbc = true;
static bool pix = true;
static bool verbose = false;
static unsigned seed = 17;
static unsigned seedsub = 23;

int main(int clc, char* clv[]){
  // Read in parameters
  intialize(clc, clv, config_file, prefix_if, prefix_of, countratefilename,
	    Nx, Ny, xmin, xmax, ymin, ymax, lambda, k_src, functionals, slide, A_blackpix, mbc, pix, verbose,
	    seed, seedsub);
  
  // Morphometric Minkowski Numerical Analysis
  double nmbr_of_microstates = boost::math::binomial_coefficient<double>(slide*slide, A_blackpix);
  std::cout << nmbr_of_microstates << " microstates;" << std::endl;
  
  double f_mod = 0.7, f_epsilon = 1e-7, Hdev = 0.8, f_initial = 1.;
  const unsigned NmbrMacrostates = 1;
  
  if(nmbr_of_microstates < 0.5*pow(10,12)){
    /*
    PhaseSpace united_states(slide, A_blackpix, nmbr_of_microstates,
			     NmbrMacrostates, f_initial, f_mod, f_epsilon, Hdev,
			     prefix_if, prefix_of, verbose, seed, seedsub);
      united_states.EVOLVE_analytic();
      united_states.output_ammendment_DoS();
    */
  }
  else{
    PhaseSpace united_states(slide, A_blackpix, nmbr_of_microstates,
			     NmbrMacrostates, f_initial, f_mod, f_epsilon, Hdev,
			     prefix_if, prefix_of, verbose, seed, seedsub);
    united_states.EVOLVE_classic();
    united_states.output_ammendment_log_DoS();

    // united_states.EVOLVE_classic_20times();
    // united_states.output_ammendment_log_DoS_all();
  }
  
  return 0;
}
