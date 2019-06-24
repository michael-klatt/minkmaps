#include "initialize.h"

void intialize(int clc, char* clv[], std::string &config_file, std::string &prefix_if, std::string &prefix_of,
	       std::string &countratefilename, unsigned &Nx, unsigned &Ny,
	       double &xmin, double &xmax, double &ymin, double &ymax,
	       double &lambda, double &k_src, std::string &functionals, unsigned &slide, unsigned &A_blackpix, bool &mbc, bool &pix, bool &verbose,
	       unsigned &seed, unsigned &seedsub)
{
  try{
    progopt::options_description generic("Generic options");
    generic.add_options()
      ("help,h", "Print options")
      ("version", "Print version number")
      ("config", progopt::value<std::string>(&config_file)->default_value(config_file),
       "Configuration file")
      ;
    progopt::options_description config("Configuration");
    config.add_options()
      ("prefix_if,i", progopt::value<std::string>(&prefix_if)->default_value(prefix_if), "Set prefix for input files")
      ("prefix_of,o", progopt::value<std::string>(&prefix_of)->default_value(prefix_of), "Set prefix for output files")
      ("countrate,K", progopt::value<std::string>(&countratefilename)->default_value(countratefilename),
       "Set filename for count rate")
      ("NmbrBinsX,U", progopt::value<unsigned>(&Nx)->default_value(Nx),
       "Set observation window number of bins in x-direction")
      ("NmbrBinsY,V", progopt::value<unsigned>(&Ny)->default_value(Ny),
       "Set observation window number of bins in y-direction")
      ("xmin,x", progopt::value<double>(&xmin)->default_value(xmin), "Set observation window minimum x-coordinate")
      ("xmax,X", progopt::value<double>(&xmax)->default_value(xmax), "Set observation window maximum x-coordinate")
      ("ymin,y", progopt::value<double>(&ymin)->default_value(ymin), "Set observation window minimum y-coordinate")
      ("ymax,Y", progopt::value<double>(&ymax)->default_value(ymax), "Set observation window maximum y-coordinate")
      ("lambda,l", progopt::value<double>(&lambda)->default_value(lambda), "Set background intensity")
      ("k_src,k", progopt::value<double>(&k_src)->default_value(k_src), "Set total number of source signals")
      ("Mink,M", progopt::value<std::string>(&functionals)->default_value(functionals),
       "Set applied Minkowski functionals")
      ("slide,w", progopt::value<unsigned>(&slide)->default_value(slide), "Set linear size of sliding observation window")
      ("A_blackpix,A", progopt::value<unsigned>(&A_blackpix)->default_value(A_blackpix), "Set number of black pixels")
      ("mbc,m", progopt::value<bool>(&mbc)->default_value(mbc), "Set boundary conditions to minus sampling")
      ("pix,p", progopt::value<bool>(&pix)->default_value(pix), "Set data format to pixelized")
      ("verbose,v", progopt::value<bool>(&verbose)->default_value(verbose), "Set output to verbose or quite")
      ("seed,s", progopt::value<unsigned>(&seed)->default_value(seed), "Set seed of random number generators")
      ("seedsub,S", progopt::value<unsigned>(&seedsub)->default_value(seedsub), "Set seed of source subtraction")
      ;
    
    progopt::options_description cmdlineopt("\nMorphometric Minkowski Analyses - Options");
    cmdlineopt.add(generic).add(config);
    progopt::options_description confopt;
    confopt.add(config);
    
    progopt::variables_map varmap;
    progopt::store(progopt::parse_command_line(clc, clv, cmdlineopt), varmap);
    progopt::notify(varmap);
    // Reading config file
    std::ifstream readconf( config_file.c_str() );
    if(readconf.fail()){
      std::cerr << "ERROR: ifstream failed to read the configuration file " << config_file << ";" << std::endl;
      exit(-1);
    }
    else{
      progopt::store(parse_config_file(readconf, confopt), varmap);
      progopt::notify(varmap);
    }
    readconf.close();
    
    if(varmap.count("help")){
      std::cout << cmdlineopt << std::endl;
      exit(0);
    }
    if(varmap.count("version")) {
      std::cout << std::endl
		<< "# ---------------------------------------------------------- " << std::endl
		<< "# - Minkowski Maps for a Morphometric Statistical Analysis - " << std::endl
		<< "# ---------------------------------------------------------- " << std::endl << std::endl;
      
      std::cout << "Version: " << __MINKMAP_VERSION << std::endl << std::endl;
      exit(0);
    }
    
    std::cout << std::endl
	      << "# ---------------------------------------------------------- " << std::endl
	      << "# - Minkowski Maps for a Morphometric Statistical Analysis - " << std::endl
	      << "# ---------------------------------------------------------- " << std::endl << std::endl;
    
    std::cout << "# Version: " << __MINKMAP_VERSION << std::endl
	      << std::endl
	      << "# Observation window number of bins in x-direction: Nx = " << Nx << std::endl
	      << "# Observation window number of bins in y-direction: Ny = " << Ny << std::endl
	      << "# Observation window minimum x-coordinate:          xmin = " << xmin << std::endl
	      << "# Observation window maximum x-coordinate:          xmax = " << xmax << std::endl
	      << "# Observation window minimum y-coordinate:          ymin = " << ymin << std::endl
	      << "# Observation window maximum y-coordinate:          ymax = " << ymax << std::endl
	      << "# Background intensity:                             lambda = " << lambda << std::endl
	      << "# Total number of source signals:                   k_src = " << k_src << std::endl
	      << "# Applied Minkowski functionals:                    Mink = " << functionals << std::endl
	      << "# Linear size of sliding observation window:        slide = " << slide << std::endl
	      << "# Number of black pixels:                           A_blackpix = " << A_blackpix << std::endl;
    if(mbc) std::cout << "# Boundary condition:                               mbc" << std::endl;
    else    std::cout << "# Boundary condition:                               pbc" << std::endl;
    if(pix) std::cout << "# Data format:                                      pix" << std::endl;
    else    std::cout << "# Data format:                                      msa" << std::endl;
    if(verbose) std::cout << "# Output format:                                    verbose" << std::endl;
    else    std::cout << "# Output format:                                    quite" << std::endl;
    std::cout << "# Seed of random number generators:                 seed = " << seed << std::endl
	      << "# Seed of source subtraction:                       seedsub = " << seedsub << std::endl
	      << std::endl
	      << "# Configuration file:                               config = " << config_file << std::endl
	      << "# Prefix for input files:                           prefix_if = " << prefix_if << std::endl
	      << "# Prefix for output files:                          prefix_of = " << prefix_of << std::endl
	      << "# Filename for count rate:                          countrate = " << countratefilename << std::endl
	      << std::endl;
    
  }
  catch(std::exception& error){
    std::cerr << "# ERROR: " << error.what() << std::endl;
    exit(-1);
  }
  catch(...) {
    std::cerr << "# Exception of unknown type!\n";
    exit(-1);
  }
  
}
