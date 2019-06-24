#ifndef INITIALIZE_H
#define INITIALIZE_H

#ifndef __MINKMAP_VERSION
#define __MINKMAP_VERSION "_unspecified version_"
#endif

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <map>
#include <cstdlib>
#include <stdarg.h>

#include <boost/program_options.hpp>
namespace progopt = boost::program_options;

void intialize(int clc, char* clv[], std::string &config_file, std::string &prefix_if, std::string &prefix_of,
	       std::string &countratefilename, unsigned &Nx, unsigned &Ny,
	       double &xmin, double &xmax, double &ymin, double &ymax,
	       double &lambda, double &k_src, std::string &functionals, unsigned &slide, unsigned &A_blackpix, bool &mbc, bool &pix, bool &verbose,
	       unsigned &seed, unsigned &seedsub);


#endif
