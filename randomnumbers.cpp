#include "randomnumbers.h"

// -------------------------
// Random Number Generators:
// -------------------------

// -------------------------
// Random Number Generator of Poisson Distributed numbers
unsigned RandomPoisson (const double &mean, const unsigned &seed)
{
  
  const gsl_rng_type * GSLrngType;
  static gsl_rng * rngPoisson = 0;
  
  /* create a generator chosen by the 
     environment variable GSL_RNG_TYPE */
  
  if (! rngPoisson) {
    gsl_rng_env_setup();
    
    GSLrngType = gsl_rng_mt19937;//gsl_rng_default;
    rngPoisson    = gsl_rng_alloc (GSLrngType);
    // seed = 0 is replaced by default seed
    gsl_rng_set (rngPoisson, seed);
  }
  
  return gsl_ran_poisson (rngPoisson, mean);
}


// Generator of Random Poisson Bin Fields
void RandomPoissonBinField (BinField<unsigned> &countrate, const BinField<double> &mean, const unsigned &seed)
{
  unsigned Nx = mean.call_Nx();
  unsigned Ny = mean.call_Ny();
  
  const gsl_rng_type * GSLrngType;
  static gsl_rng * rngPoissonBinField = 0;
  
  if (! rngPoissonBinField) {
    gsl_rng_env_setup();
    
    GSLrngType = gsl_rng_mt19937;//gsl_rng_default;
    rngPoissonBinField    = gsl_rng_alloc (GSLrngType);
    // seed = 0 is replaced by default seed
    gsl_rng_set (rngPoissonBinField, seed);
  }
  
  for(unsigned ni = 0; ni < Nx*Ny; ni++)
    countrate.assign(ni, gsl_ran_poisson (rngPoissonBinField, mean.call(ni)) );
}
// -------------------------


// -------------------------
// Uniform Random Number Generator: 0 <= sample < 1; sample in [0;1)
double RandomRational (const unsigned &seed)
{
  
  const gsl_rng_type * GSLrngType;
  static gsl_rng * rngUniform = 0;
  
  if (! rngUniform) {
    gsl_rng_env_setup();
    
    GSLrngType = gsl_rng_mt19937;//gsl_rng_default;
    rngUniform    = gsl_rng_alloc (GSLrngType);
    // seed = 0 is replaced by default seed
    gsl_rng_set (rngUniform, seed);
  }
  
  return gsl_rng_uniform (rngUniform);
}


// Random numbers r with r >=min and r < max
double RandomUniform (const double &min, const double &max, const unsigned &seed)
{
  
  // Mersenne twister random number generator
  static boost::mt19937 rngenerator(seed+1);
  
  // Uniform distribution
  boost::uniform_real<double> uniform_pdf(min,max);
  
  // Binding random number generator to distribution, forming a function
  boost::variate_generator< boost::mt19937&, boost::uniform_real<double> >  uniform_rng (rngenerator, uniform_pdf);
  
  return uniform_rng();
}


// Random numbers n with n >=0 and n <= N
unsigned RandomNumber (const unsigned &N, const unsigned &seed)
{
  
  // Mersenne twister random number generator
  static boost::mt19937 rngenerator(seed+1);
  
  // Uniform distribution
  boost::uniform_int<> uniform_pdf(0,N);
  
  // Binding random number generator to distribution, forming a function
  boost::variate_generator< boost::mt19937&, boost::uniform_int<> >  uniform_rng (rngenerator, uniform_pdf);
  
  return uniform_rng();
}


// Random numbers n with n >=0 and n <= N
int RandomInteger (const unsigned &N, const unsigned &seed)
{
  
  // Mersenne twister random number generator
  static boost::mt19937 rngenerator(seed);
  
  // Uniform distribution
  boost::uniform_int<> uniform_integer_pdf(0,N);
  
  // Binding random number generator to distribution, forming a function
  boost::variate_generator< boost::mt19937&, boost::uniform_int<> >  uniform_rng (rngenerator, uniform_integer_pdf);
  
  return uniform_rng();
}


// Random integer numbers i with i >=n and i <= N (different seed to RandomNumber)
int RandomIntegerBetweenUpperLowerBound (const unsigned &n, const unsigned &N, const unsigned &seed)
{
  
  // Mersenne twister random number generator
  static boost::mt19937 rngenerator(seed);
  
  // Uniform distribution
  boost::uniform_int<> uniform_integer_pdf(n,N);
  
  // Binding random number generator to distribution, forming a function
  boost::variate_generator< boost::mt19937&, boost::uniform_int<> >  uniform_rng (rngenerator, uniform_integer_pdf);
  
  return uniform_rng();
}


// Random numbers n with n >=0 and n < N
int RandomGSLInteger (const unsigned &N, const unsigned &seed)
{
  
  const gsl_rng_type * GSLrngType;
  static gsl_rng * rngInteger = 0;
  
  /* create a generator chosen by the 
     environment variable GSL_RNG_TYPE */
  
  if (! rngInteger) {
    gsl_rng_env_setup();
    
    GSLrngType = gsl_rng_mt19937;//gsl_rng_default;
    rngInteger    = gsl_rng_alloc (GSLrngType);
    // seed = 0 is replaced by default seed
    gsl_rng_set (rngInteger, seed);
  }
  
  return gsl_rng_uniform_int (rngInteger, N);
}


// Bernoulli Experiment
bool RandomBernoulli (const double &p, const unsigned &seed)
{
  if( p < 0 || p > 1 )
    {
      std::cerr << "ERROR: RandomBernoulli recieved as probability p: " << p << std::endl;
      exit(-1);
    }
  
  const gsl_rng_type * GSLrngType;
  static gsl_rng * rngThresh = 0;
  
  if (! rngThresh) {
    gsl_rng_env_setup();
    
    GSLrngType = gsl_rng_mt19937;//gsl_rng_default;
    rngThresh    = gsl_rng_alloc (GSLrngType);
    // seed = 0 is replaced by default seed
    gsl_rng_set (rngThresh, seed);
  }
  
  return ( gsl_rng_uniform(rngThresh) < p );
}


// Bernoulli Experiment on a BinField
void RandomBernoulliBinField (BinField<bool> &sample, const BinField<double> &p, const unsigned &seed)
{
  unsigned Nx = p.call_Nx();
  unsigned Ny = p.call_Ny();
  
  if(Nx != sample.call_Nx() || Ny != sample.call_Ny() )
    {
      std::cerr << "ERROR: RandomBernoulliBinField recieved per reference sample and p BinFields;" << std::endl
		<< "       Dimensions do not match." << std::endl;
      exit(-1);
    }
  
  const gsl_rng_type * GSLrngType;
  static gsl_rng * rngThreshBinField = 0;
  
  if (! rngThreshBinField) {
    gsl_rng_env_setup();
    
    GSLrngType = gsl_rng_mt19937;//gsl_rng_default;
    rngThreshBinField    = gsl_rng_alloc (GSLrngType);
    // seed = 0 is replaced by default seed
    gsl_rng_set (rngThreshBinField, seed);
  }
  
  for(unsigned ni = 0; ni < Nx*Ny; ni++)
    sample.assign(ni, gsl_rng_uniform (rngThreshBinField) < p.call(ni));
}
// -------------------------


// -------------------------
// Binomial Random Number Generator
unsigned RandomBinomial (const unsigned &n, const double &p, const unsigned &seed)
{
  const gsl_rng_type * GSLrngType;
  static gsl_rng * rngBino = 0;
  
  if (! rngBino) {
    gsl_rng_env_setup();
    
    GSLrngType = gsl_rng_mt19937;//gsl_rng_default;
    rngBino    = gsl_rng_alloc (GSLrngType);
    // seed = 0 is replaced by default seed
    gsl_rng_set (rngBino, seed);
  }
  
  return gsl_ran_binomial(rngBino, p, n);
}


// Subtracting countrate BinField by Binomial Random Number Generator
void SubtractRandomBinomial (BinField<unsigned> &countrate, const double &lambda, const BinField<double> &inty, const unsigned &seed)
{
  unsigned Nx = inty.call_Nx();
  unsigned Ny = inty.call_Ny();
  
  if(Nx != countrate.call_Nx() || Ny != countrate.call_Ny() )
    {
      std::cerr << "ERROR: RandomBernoulliBinField recieved per reference countrate and p BinFields;" << std::endl
		<< "       Dimensions do not match." << std::endl;
      exit(-1);
    }
  const gsl_rng_type * GSLrngType;
  static gsl_rng * rngBino = 0;
  
  if (! rngBino) {
    gsl_rng_env_setup();
    
    GSLrngType = gsl_rng_mt19937;//gsl_rng_default;
    rngBino    = gsl_rng_alloc (GSLrngType);
    // seed = 0 is replaced by default seed
    gsl_rng_set (rngBino, seed);
  }
  
  for(unsigned ni = 0; ni < Nx*Ny; ni++)
    countrate.assign( ni, gsl_ran_binomial( rngBino, lambda/(lambda + inty.call(ni)), countrate.call(ni)) );
  
}
// -------------------------


// Random number Gaussian distributed
double RandomGaussian(const double &mean, const double &sigma, const unsigned int &seed)
{
  const gsl_rng_type * GSLrngTypeGauss;
  static gsl_rng * rngGauss = 0;
  
  if (! rngGauss) {
    gsl_rng_env_setup();
    
    GSLrngTypeGauss = gsl_rng_mt19937;//gsl_rng_default;
    rngGauss    = gsl_rng_alloc (GSLrngTypeGauss);
    // seed = 0 is replaced by default seed
    gsl_rng_set (rngGauss, seed);
  }
  
  return (mean + gsl_ran_gaussian (rngGauss, sigma)) ;
}

// -------------------------



// -----------------------------------------
// PDFs for Random Number Distributions:
// -----------------------------------------


// Function value of the Poisson-pdf
double PoissonPDF(const unsigned &k, const double &mean)
{
  return gsl_ran_poisson_pdf(k, mean);
}


// Function value of the Poisson-CDF, which is the sum of PoissonPDF(k) for k <= K
double PoissonCDF(const unsigned &K, const double &mean)
{
  // P  
  return gsl_cdf_poisson_P(K,mean);
}


// Function value of the Poisson-CDF, which is the sum of PoissonPDF(k) for k > K
double PoissonComplementaryCDF(const unsigned &K, const double &mean)
{
  // Q  
  return gsl_cdf_poisson_Q(K,mean);
}


// Given an intensity lambda: what is the probability p for recieving a filled bin at threshold thresh?
double Bernoulli_p(const unsigned &thresh, const double &lambda)
{
  if(thresh)
    return gsl_cdf_poisson_Q(thresh-1,lambda);
  else
    return 1;
}


// Functional value of the Binomial-pdf
double BinomialPDF(const unsigned &M, const unsigned &m, const double &p)
{
  return gsl_ran_binomial_pdf(m, p, M);
}


// Functional value of the Binomial-CDF, which is the sum of BinomialPDF(mm) for mm <= m
double BinomialCDF(const unsigned &M, const unsigned &m, const double &p)
{
  // P  
  return gsl_cdf_binomial_P(m, p, M);
}


// Functional value of the Binomial-CDF, which is the sum of BinomialPDF(mm) for mm >= m
double BinomialCDF_complement(const unsigned &M, const unsigned &m, const double &p)
{
  // Q
  return (gsl_cdf_binomial_Q(m, p, M) +  gsl_ran_binomial_pdf(m, p, M));
}
