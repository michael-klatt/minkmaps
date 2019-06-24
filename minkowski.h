#ifndef MINKOWSKI_H
#define MINKOWSKI_H

#include "lookuptable.h"
#include "randomnumbers.h"

typedef double (*Minkowski)(const BinField<bool>&);

// Minkowski Sky Map Tools
void TurnConfNmbrIntoBinField(const unsigned &conf, BinField<bool> &Configuration,
			      const unsigned &slide, const unsigned &N_conf);

bool eq(const double &first, const double &second);

// MINKOWSKI FUNCTIONALS: SUMMING LOOK-UP TABLE (MARCHING SQUARE ALGORITHM)
unsigned convert(const bool &right_low, const bool &left_low, const bool &right_up, const bool &left_up);

// Minus sampling boundary condition
double area_mbc(const BinField<bool> &sample);
double area_mbc(const BinField<bool> &sample,
		const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny);

double perimeter_mbc(const BinField<bool> &sample);
double perimeter_mbc(const BinField<bool> &sample,
		     const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny);

double euler_mbc(const BinField<bool> &sample);
double euler_mbc(const BinField<bool> &sample,
		 const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny);

double w102_xx_mbc(const BinField<bool> &sample);
double w102_xx_mbc(const BinField<bool> &sample,
		   const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny);

double w102_xy_mbc(const BinField<bool> &sample);
double w102_xy_mbc(const BinField<bool> &sample,
		   const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny);

double w102_yy_mbc(const BinField<bool> &sample);
double w102_yy_mbc(const BinField<bool> &sample,
		   const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny);

double delta_mbc(const BinField<bool> &sample);
double delta_mbc(const BinField<bool> &sample,
		 const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny);

// Periodic boundary condition
double area_pbc(const BinField<bool> &sample);
double area_pbc(const BinField<bool> &sample,
		const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny);

double perimeter_pbc(const BinField<bool> &sample);
double perimeter_pbc(const BinField<bool> &sample,
		     const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny);

double euler_pbc(const BinField<bool> &sample);
double euler_pbc(const BinField<bool> &sample,
		 const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny);

double w102_xx_pbc(const BinField<bool> &sample);
double w102_xx_pbc(const BinField<bool> &sample,
		   const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny);

double w102_xy_pbc(const BinField<bool> &sample);
double w102_xy_pbc(const BinField<bool> &sample,
		   const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny);

double w102_yy_pbc(const BinField<bool> &sample);
double w102_yy_pbc(const BinField<bool> &sample,
		   const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny);

double delta_pbc(const BinField<bool> &sample);
double delta_pbc(const BinField<bool> &sample,
		 const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny);

// Look-up table but pixelized data: all functional values times 8:
int area_mbc_pix(const BinField<bool> &sample);

int perimeter_mbc_pix(const BinField<bool> &sample);

int euler_mbc_pix(const BinField<bool> &sample);

int area_wbc_pix(const BinField<bool> &sample);

int perimeter_wbc_pix(const BinField<bool> &sample);

int euler_wbc_pix(const BinField<bool> &sample);

// PIXELIZED DATA
// Minus sampling boundary condition
double area_pixelized_mbc(const BinField<bool> &sample);
double area_pixelized_mbc(const BinField<bool> &sample,
			  const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny);

double perimeter_pixelized_mbc(const BinField<bool> &sample);
double perimeter_pixelized_mbc(const BinField<bool> &sample,
			       const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny);

double euler_pixelized_mbc(const BinField<bool> &sample);
double euler_pixelized_mbc(const BinField<bool> &sample,
			   const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny);

// Periodic boundary condition
double area_pixelized_pbc(const BinField<bool> &sample);
double area_pixelized_pbc(const BinField<bool> &sample,
			  const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny);

double perimeter_pixelized_pbc(const BinField<bool> &sample);
double perimeter_pixelized_pbc(const BinField<bool> &sample,
			       const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny);

double euler_pixelized_pbc(const BinField<bool> &sample);
double euler_pixelized_pbc(const BinField<bool> &sample,
			   const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny);

// PAPAYA
void print_pgm(const BinField<bool> &sample,
	       const std::string &filename, const std::string &prefix_of, const bool &invert = false);
void print_pgm(const BinField<bool> &sample,
	       const unsigned &xi, const unsigned &yi, const unsigned &Nx, const unsigned &Ny,
	       const std::string &filename, const std::string &prefix_of, const bool &invert = false);

#endif


