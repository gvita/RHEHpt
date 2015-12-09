#include "complex_def.h"

std::complex<long double> operator*(const std::complex<long double>& a, const double& b){
	return a*static_cast<long double>(b);
}

std::complex<long double> operator*(const double& a,const std::complex<long double>& b){
	return static_cast<long double>(a)*b;
}

std::complex<long double> operator/(const std::complex<long double>& a, const double& b){
	return a/static_cast<long double>(b);
}

std::complex<long double> operator/(const double& a,const std::complex<long double>& b){
	return static_cast<long double>(a)/b;
}

std::complex<long double> operator-(const std::complex<long double>& a, const double& b){
	return a-static_cast<long double>(b);
}

std::complex<long double> operator-(const double& a,const std::complex<long double>& b){
	return static_cast<long double>(a)-b;
}

std::complex<long double> operator+(const std::complex<long double>& a, const double& b){
	return a+static_cast<long double>(b);
}

std::complex<long double> operator+(const double& a,const std::complex<long double>& b){
	return static_cast<long double>(a)+b;
}

//Useful Math functions
long double log_r(long double x){
//	std::cout << "log_r(" << x << ")" << std::endl;
	return gsl_sf_log(x);
}
long double dilog_r(long double x){
//	std::cout << "dilog_r(" << x << ")" << std::endl;
	return gsl_sf_dilog(x);
}

std::complex<long double> log_c(std::complex<long double> z){
	if ((imag(z)>-1E-10)&&(imag(z)<1E-10)){
	if(real(z) < 0){
		std::complex<long double> ris(log_r(abs(z)),(long double) (M_PI));
		return ris;
	}
	else {
		std::complex<long double> ris(log_r(abs(z)),0.0);
		return ris;
	} 
	}
	gsl_sf_result lnr;
	gsl_sf_result theta;
	gsl_sf_complex_log_e(std::real(z),std::imag(z),&lnr, &theta);
	std::complex<long double> ris(lnr.val,theta.val);
	return ris;
}
std::complex<long double> dilog_c(std::complex<long double> z){
	gsl_sf_result ris_r;
	gsl_sf_result ris_i;
	gsl_sf_complex_dilog_e(std::abs(z),std::arg(z),&ris_r, &ris_i);
	std::complex<long double> ris(ris_r.val,ris_i.val);
	return ris;
}
