#ifndef __supporting_functions__
#define __supporting_functions__

#include <complex>
#include <functional>
#include <vector>

long double InverseMellin(double z, std::function<  std::complex<long double>(std::complex<long double>)>& f);
std::vector<long double> InverseMellin(double z, std::function< std::vector< std::complex<long double> >(std::complex<long double>)>& f);

#endif
