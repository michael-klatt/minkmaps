#ifndef RANDOM_NUMBERS_H
#define RANDOM_NUMBERS_H

// RANDOMNUMBER GENERATORS
// DEFAULT: MERSENNE TWISTER

#include "binfield.h"

// Boost
#include <boost/random.hpp>
#include <boost/math/distributions/poisson.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
// GSL
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

// -------------------------
// Random Number Generators:
// -------------------------


// -------------------------
// Random Number Generator of Poisson Distributed numbers
unsigned RandomPoisson (const double &mean, const unsigned &seed);


// Generator of Random Poisson Bin Fields
void RandomPoissonBinField (BinField<unsigned> &countrate, const BinField<double> &mean, const unsigned &seed);
// -------------------------


// -------------------------
// Uniform Random Number Generator: 0 <= sample < 1; sample in [0;1)
double RandomRational (const unsigned &seed);


// Random numbers r with r >=min and r < max
double RandomUniform (const double &min, const double &max, const unsigned &seed);


// Random numbers n with n >=0 and n <= N
unsigned RandomNumber (const unsigned &N, const unsigned &seed);


// Random integer numbers n with n >=0 and n <= N (different seed to RandomNumber)
int RandomInteger (const unsigned &N, const unsigned &seed);


// Random integer numbers i with i >=n and i <= N (different seed to RandomNumber)
int RandomIntegerBetweenUpperLowerBound (const unsigned &n, const unsigned &N, const unsigned &seed);


// Random numbers n with n >=0 and n < N
int RandomGSLInteger (const unsigned &N, const unsigned &seed);

// Bernoulli Experiment
bool RandomBernoulli (const double &p, const unsigned &seed);


// Bernoulli Experiment on a BinField
void RandomBernoulliBinField (BinField<bool> &sample, const BinField<double> &p, const unsigned &seed);
// -------------------------


// -------------------------
// Binomial Random Number Generator
unsigned RandomBinomial (const unsigned &n, const double &p, const unsigned &seed);


// Subtracting countrate BinField by Binomial Random Number Generator
void SubtractRandomBinomial (BinField<unsigned> &countrate, const double &lambda, const BinField<double> &inty, const unsigned &seed);
// -------------------------


// Random number Gaussian distributed
double RandomGaussian(const double &mean, const double &sigma, const unsigned int &seed);


// -----------------------------------------
// PDFs for Random Number Distributions:
// -----------------------------------------


// Function value of the Poisson-pdf
double PoissonPDF(const unsigned &k, const double &mean);


// Function value of the Poisson-CDF, which is the sum of PoissonPDF(k) for k <= K
double PoissonCDF(const unsigned &K, const double &mean);


// Function value of the Poisson-CDF, which is the sum of PoissonPDF(k) for k > K
double PoissonComplementaryCDF(const unsigned &K, const double &mean);


// Given an intensity lambda: what is the probability p for recieving a filled bin at threshold thresh?
double Bernoulli_p(const unsigned &thresh, const double &lambda);


// Functional value of the Binomial-pdf
double BinomialPDF(const unsigned &M, const unsigned &m, const double &p);


// Functional value of the Binomial-CDF, which is the sum of BinomialPDF(mm) for mm <= m
double BinomialCDF(const unsigned &M, const unsigned &m, const double &p);


// Functional value of the Binomial-CDF, which is the sum of BinomialPDF(mm) for mm >= m
double BinomialCDF_complement(const unsigned &M, const unsigned &m, const double &p);


#endif
