#include <complex>
#include <cmath>
#include "complex_def.h"
#include "Inverse_Mellin.h"

long double InverseMellin(double z, std::function<  std::complex<long double>(std::complex<long double>)> f)
{ // Talbot algorithm for a function f : C -> C
  if ( z == 1. ) return 0.;
  long double log1ox = -std::log(z);
  // prec ~ 10^{-0.6*M}
  unsigned int M = 50;
  long double r = 1.+2.*M/5./log1ox;
  std::complex<long double> integrand = f(r) ;

  long double t, cott, sigma;
  std::complex<long double> N;
  std::complex<long double> I(0,1);
  long double res = std::real(0.5 * integrand * exp( r*log1ox )) ;

  for(unsigned k = 1; k < M; k++) {
    t = k*M_PI/M;
    cott = 1./std::tan(t);
    sigma = t + cott * ( t * cott - 1. );
    N = r * t * ( cott + I );
  	integrand =  f(N);
    res += std::real( exp( N*log1ox ) * integrand * (1.+I*sigma) );
  }
  return res*r/M;
}


std::vector<long double> InverseMellin(double z, std::function< std::vector< std::complex<long double> >(std::complex<long double>)> f)
{ // Talbot algorithm for a function f : C -> C^size
  long double log1ox = -std::log(z);
  // prec ~ 10^{-0.6*M}
  unsigned int M = 50;
  long double r = 1.+2.*M/5./log1ox;
  std::vector< std::complex<long double> > integrand = f(r) ;
  std::size_t size = integrand.size(); // = res.size() = dimension of the codomain of f
  if(z==1.) return std::vector<long double>(size,0.);
  
  long double t, cott, sigma;
  std::complex<long double> N;
  std::complex<long double> I(0,1);
  std::vector<long double> res; 


  for (auto f_r : integrand )
  	res.push_back( std::real(0.5 * f_r * std::exp( r*log1ox ))*r/M );

  for(unsigned k = 1; k < M; k++) {
    t = k*M_PI/M;
    cott = 1./std::tan(t);
    sigma = t + cott * ( t * cott - 1. );
    N = r * t * ( cott + I );
  	integrand =  f(N);
  	for (unsigned i = 0; i < size ; ++i)
    	res[i] += std::real( std::exp( N*log1ox ) * integrand[i] * (1.+I*sigma) )*r/M;
  }
  return res;
}

