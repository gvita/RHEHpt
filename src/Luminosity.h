#ifndef LUM_H
#define LUM_H

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
		typedef std::complex<long double> dcomplex;
		typedef LHAPDF::PDF* PDF_ptr;
	
		Luminosity(PDF_ptr thePDF, double MUF=125.09, double MUR=125.09, std::size_t order=20);
		Luminosity(const Luminosity& Lumi);
		virtual ~Luminosity();
	
		
		
		
		void Cheb_Lum(double );
		void c_large();

		dcomplex CLum_N(dcomplex N, unsigned channel) const;
		long double CLum_x(long double z, unsigned channel) const;

		dcomplex CLum_gg_N(dcomplex N) const {	return CLum_N(N,0);	}
		long double CLum_gg_x(long double z) const { return CLum_x(z,0);  }
		//dcomplex CLum_gq_N(dcomplex N) const {	return CLum_N(N,0);	}
		//dcomplex CLum_qbarq_N(dcomplex N) const {	return CLum_N(N,0);	}
		//dcomplex CLum_qq_N(dcomplex N) const {	return CLum_N(N,0);	}
		//dcomplex CLum_qqprime_N(dcomplex N) const {	return CLum_N(N,0);	}
		//dcomplex CLum_qqbarprime_N(dcomplex N) const {	return CLum_N(N,0);	}

		double xfxQ(int channel, double x){
		  return (_PDF->xfxQ(channel,x,_Muf)/x);
		}
		double get_alphas(){
		  return(_PDF->alphasQ(_Mur));
		}
		
		inline PDF_ptr get_PDF(){
		  return _PDF;
		}
		inline double get_muf(){
		  return _Muf;
		}
		inline double get_mur(){
		  return _Mur;
		}
		inline size_t get_n(){
		  return _n;
		}
		
		struct par_struct{
			double Muf;
			PDF_ptr PDF;
			double t;
		};
		//void cheb_approx(double , double);

	private:
		size_t _n;	//Order of Chebishev approximation
		std::array< gsl_cheb_series*,1>  Chebseries;
		std::array< std::vector<double> ,1>  C_larges;

		std::vector< std::vector< double > > _T;
		LHAPDF::PDF* _PDF;
	
//		static double _Lgg(double x,void* p);
//		static double _Lum_gg(double tau);
		
		double _Muf;
		double _Mur;
	//	gsl_cheb_series *Chebgg;
	//	gsl_cheb_series *Chebqq;
	//	gsl_cheb_series *Chebgq;
	//	gsl_cheb_series *Chebqq2;
	//	gsl_cheb_series *Chebqqprime;
	//	gsl_cheb_series *Chebqqbarprime;

	//	std::vector < double > C_largegg;
	//	std::vector < double > C_largeqq;
	//	std::vector < double > C_largegq;
	//	std::vector < double > C_largeqq2;
	//	std::vector < double > C_largeqqprime;
	//	std::vector < double > C_largeqqbarprime;
};

//double Lum_gg(double );

#endif
