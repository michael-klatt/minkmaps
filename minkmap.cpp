// Minkowski Maps for a Morphometric Statistical Analysis
#include "binnedskymap.h"
#include "phasespace.h"

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
static unsigned slide = 15;
static unsigned A_blackpix = 0;
static bool mbc = false;
static bool pix = false;
static bool verbose = true;
static unsigned seed = 17;
static unsigned seedsub = 23;

void add_ellipse(BinnedSkyMap &skylab, double x0, double y0, double phi, double a, double b);

int main(int clc, char* clv[]){
  // Read in parameters
  intialize(clc, clv, config_file, prefix_if, prefix_of, countratefilename,
	    Nx, Ny, xmin, xmax, ymin, ymax, lambda, k_src, functionals, slide, A_blackpix, mbc, pix, verbose,
	    seed, seedsub);
    
  BinnedSkyMap skylab(Nx,Ny,xmin,xmax,ymin,ymax,prefix_if,prefix_of,seed,seedsub);
  skylab.add_bck_intensity(lambda);
    
  add_ellipse(skylab, 0, 0, 0, 0.07, 0.07);
  add_ellipse(skylab, -7, -7, 0, 0.07, 0.07);
  add_ellipse(skylab, -7, 6, 0, 0.07, 0.07);
  add_ellipse(skylab, 6, -7, 0, 0.07, 0.07);
  add_ellipse(skylab, 6, 6, 0, 0.07, 0.07);
  
  int N_runs = 1;
  
  skylab.plot_intensity("intensity");
  skylab.simulate_countrate();
  skylab.plot_countrate("counts_map");
  
  for(int run=0;run<N_runs;run++){
    std::cout << "run " << run << " of 0 ... " << N_runs-1 << std::endl;
    skylab.simulate_countrate();
    
    std::stringstream name;
    name << "singleMSM-joint-" << run;
    skylab.Minkowski_sky_map_APC_pix_wbc(lambda, slide, verbose, true, true, name.str());

    std::stringstream other;
    other << "singleMSM-simple-" << run;
    skylab.Minkowski_sky_map_A_pix_wbc(lambda, slide, verbose, other.str());
  }

  return 0;
}

void add_ellipse(BinnedSkyMap &skylab, double x0, double y0, double phi, double a, double b){
  
  for(unsigned xi = 0; xi < Nx; xi++)
    for(unsigned yi = 0; yi < Nx; yi++){
      double x = xi  - 0.5*Nx - x0;
      double y = yi  - 0.5*Ny - y0;
      
      double xprime = cos(phi*M_PI)*x - sin(phi*M_PI)*y;
      double yprime = sin(phi*M_PI)*x + cos(phi*M_PI)*y;
      
      double semi_a = a*Nx;
      double semi_b = b*Ny;
      if(xprime*xprime/semi_a/semi_a + yprime*yprime/semi_b/semi_b < 1)
	skylab.add_intensity_to_single_bin(0.7*k_src, xi, yi);
    }
  
}

