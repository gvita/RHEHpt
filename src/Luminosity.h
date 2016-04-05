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
		
		dcomplex CLum_gq_N(dcomplex N) const {	return CLum_N(N,1);	}
		long double CLum_gq_x(long double z) const { return CLum_x(z,1);  }
		
		dcomplex CLum_qqbar_N(dcomplex N) const {	return CLum_N(N,2);	}
		long double CLum_qqbar_x(long double z) const { return CLum_x(z,2);  }
		
		dcomplex CLum_qq_N(dcomplex N) const {	return CLum_N(N,3);	}
		long double CLum_qq_x(long double z) const { return CLum_x(z,3);  }
		
		dcomplex CLum_qqbarprime_N(dcomplex N) const {	return CLum_N(N,4);	}
		long double CLum_qqbarprime_x(long double z) const {	return CLum_x(z,4);	}
		
		dcomplex CLum_qqprime_N(dcomplex N) const {	return CLum_N(N,5);	}
		long double CLum_qqprime_x(long double z) const { return CLum_x(z,5);  }
		

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
			unsigned int Nf;
		};
		
		void setNf(unsigned int NF){
		  if ((NF < 1)|| (NF>6)){
		    std::cout << "ERROR: INCONGRUENT NUMBER OF FLAVOUR (1-6 ACCEPTED)" << std::endl;
		    NF=5.;
		  }
		  _Nf=NF;
		}
		//void cheb_approx(double , double);

	private:
		size_t _n;	//Order of Chebishev approximation
		std::array< gsl_cheb_series*,6>  Chebseries;
		std::array< std::vector<double> ,6>  C_larges;

		std::vector< std::vector< double > > _T;
		LHAPDF::PDF* _PDF;
		
		double _Muf;
		double _Mur;
		unsigned int _Nf;

};

#endif
