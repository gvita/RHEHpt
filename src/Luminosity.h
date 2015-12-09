#ifndef LUM_H
#define LUM_H

//#include "costant.h"
//#include "complex.h"
#include <LHAPDF/LHAPDF.h>
#include <cmath>
#include <iostream>
#include <gsl/gsl_chebyshev.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_dilog.h>
#include <algorithm>
#include <gsl/gsl_integration.h>


class Luminosity {
public:
	typedef std::vector<double> Rvector;
	typedef std::complex<double> dcomplex;

	Luminosity(std::shared_ptr<LHAPDF::PDF> thePDF, std::size_t order=15);
	virtual ~Luminosity();
	
/*	double* Cogg(){
	return gsl_cheb_coeffs(Chebgg);
	}
	double* Coqq(){
	return gsl_cheb_coeffs(Chebqq);
	}
	double* Cogq(){
	return gsl_cheb_coeffs(Chebgq);
	}
	double* Coqq2(){
	return gsl_cheb_coeffs(Chebqq2);
	}
	double* Coqqprime(){
	return gsl_cheb_coeffs(Chebqqprime);
	}
	double* Coqqbarprime(){
	return gsl_cheb_coeffs(Chebqqbarprime);
	}	
	double* Coeffgg(){
	return C_largegg;
	}
	double* Coeffqq(){
	return C_largeqq;
	}
	double* Coeffgq(){
	return C_largegq;
	}
	double GetT(int i, int j){
	return T[i][j];
	}*/
	void Cheb_Lum(double );
	void c_large();
/*	double Lum_gg(double );
	double Lum_gg_x(double );
	double Lum_gg_N(double );
	double Lum_qbarq(double );
	double Lum_qq(double );
	double Lum_qbarq_x(double );
	double Lum_gq_x(double );
	double Lum_qq_x(double );
	double Lum_qqbarprime_x(double );
	double Lum_qqprime_x(double );	
	double Lum_qqprimo(double );
	double Lum_qbarqprimo(double );
	double Lum_gq(double );*/
	dcomplex CLum_gg_N(dcomplex );	
	dcomplex CLum_gq_N(dcomplex );
	dcomplex CLum_qbarq_N(dcomplex );
	dcomplex CLum_qq_N(dcomplex );
	dcomplex CLum_qqprime_N(dcomplex );
	dcomplex CLum_qqbarprime_N(dcomplex);
	//void cheb_approx(double , double);

private:
	size_t _n;	//Order of Chebishev approximation
	gsl_cheb_series *Chebgg;
	gsl_cheb_series *Chebqq;
	gsl_cheb_series *Chebgq;
	gsl_cheb_series *Chebqq2;
	gsl_cheb_series *Chebqqprime;
	gsl_cheb_series *Chebqqbarprime;
	std::vector<Rvector> _T;
	Rvector C_largegg;
	Rvector C_largeqq;
	Rvector C_largegq;
	Rvector C_largeqq2;
	Rvector C_largeqqprime;
	Rvector C_largeqqbarprime;
	std::shared_ptr< LHAPDF::PDF > _PDF;
};

//double Lum_gg(double );

#endif
