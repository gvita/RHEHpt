#ifndef __COMPLEX_DEF_H__
#define __COMPLEX_DEF_H__

#include <complex>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_result.h>

std::complex<long double> operator*(const std::complex<long double>& a, const double& b);

std::complex<long double> operator*(const double& a,const std::complex<long double>& b);

std::complex<long double> operator/(const std::complex<long double>& a, const double& b);

std::complex<long double> operator/(const double& a,const std::complex<long double>& b);

std::complex<long double> operator-(const std::complex<long double>& a, const double& b);

std::complex<long double> operator-(const double& a,const std::complex<long double>& b);

std::complex<long double> operator+(const std::complex<long double>& a, const double& b);

std::complex<long double> operator+(const double& a,const std::complex<long double>& b);

//Useful Math functions
long double log_r(long double x);
long double dilog_r(long double x);

const std::complex<long double> II(0.0,1.0);
std::complex<long double> log_c(std::complex<long double> z);
std::complex<long double> dilog_c(std::complex<long double> z);
#endif /* __COMPLEX_DEF_H__ */
