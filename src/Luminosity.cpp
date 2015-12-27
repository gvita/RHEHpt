#include "Luminosity.h"
#include "complex_def.h"
#include <functional>

Luminosity::Luminosity(PDF_ptr thePDF, std::size_t n):_n(n),_PDF(thePDF)
{	
	for (auto& it : Chebseries)
		it = gsl_cheb_alloc(_n);
	for (auto& it : C_larges)	
		it = std::vector < double >( _n+1 , 0 );

	for(std::size_t i = 0; i < (_n + 1) ; ++i)
		_T.push_back( std::vector < double >(_n+1,0) );

	_T[0][0]=1.0;
	_T[1][0]=0.0;
	_T[1][1]=1.0;
	for(std::size_t i = 2 ; i < (_n+1) ; i++){
		_T[i][0] = - _T[i-2][0];
		for(std::size_t j = 1 ; j < (i-1) ; j++){
			_T[i][j] = 2*_T[i-1][j-1] - _T[i-2][j];
		}
		_T[i][i-1] = 2*_T[i-1][i-2];
		_T[i][i] = 2*_T[i-1][i-1];
	}

}

Luminosity::~Luminosity(){
//	std::cout << "distruggendo " << std::endl;
 	for (auto& it : Chebseries)
		gsl_cheb_free(it);
//	std::cout << "distrutto " << std::endl;
}

double f_gg(double x1,void *p){
//	std::cout.precision(15);
	Luminosity::par_struct* p_struct = (Luminosity::par_struct *) p;
	double x2 = p_struct -> t /x1;
	double Muf = p_struct -> Muf; 
	Luminosity::PDF_ptr PDF = p_struct -> PDF;
//	std::cout << " pars : " << std::endl;
//	std::cout << " muf : " <<  Muf << std::endl;
//	std::cout << " tau : " <<  p_struct -> t << std::endl;
//	std::cout << " PDF : " <<  PDF << std::endl;
//	std::cout << " x1 = " <<  x1 << std::endl;
//	std::cout << " x2 = " <<  x2 << std::endl;			
//	std::cout << "PDF: " << x1 << " " << ( PDF-> xfxQ(0,x1,Muf)/x1 * PDF -> xfxQ(0,x2,Muf)/x2 ) / x1 << std::endl;
	return ( PDF-> xfxQ(0,x1,Muf)/x1 * PDF -> xfxQ(0,x2,Muf)/x2 ) / x1;
}


double _Lum_gg(double tau,void* p){
	Luminosity::par_struct* par = (Luminosity::par_struct*) p;
	par -> t = tau;
	double precision = 1e-4;
	double ris = 0.0, error = 0.0;
	gsl_function Integrand;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	Integrand.function = f_gg;
	Integrand.params = par;
	gsl_integration_qags (&Integrand, tau, 1.0, 0, precision, 1000,	w, &ris, &error);
	gsl_integration_workspace_free (w);
	return ris;
}

double _Lgg(double x, void *p){
	double z = std::exp(x);
	return (z*_Lum_gg(z,p));
}

void Luminosity::Cheb_Lum(double muf){
	double umin=-12.0;
	for (auto& it : C_larges)	
		it = std::vector < double >( _n+1 , 0 );
	
	for (auto& it : Chebseries){
		gsl_cheb_free(it);
		it = gsl_cheb_alloc(_n);
	}
	
	gsl_function Largegg;

	Largegg.function = &_Lgg;
	par_struct par;
	par.PDF = _PDF;
	par.Muf = muf;

	Largegg.params = &par;
	gsl_cheb_init(Chebseries[0],&Largegg,umin,0.0);
	
	c_large();

 }

//void Luminosity::Resize(size_t n){
//	_n = n;
//	for (auto it : Chebseries)
//		gsl_cheb_free(it);
//		it = gsl_cheb_alloc(_n);
//}

void Luminosity::c_large(){
	double umin = -12.0;

	std::vector<double*> co_channels;
	for (auto it : Chebseries)
		co_channels.push_back(gsl_cheb_coeffs(it));

	for(unsigned k = 0; k < Chebseries.size() ; ++k){
		C_larges[k][0] = - co_channels[k][0]/2.;
		for(unsigned l = 0 ; l < (_n+1) ; l++ )
		 	for (unsigned i = l ; i < (_n+1) ; i++ )
				for(unsigned j = i ; j < (_n+1) ; j++ )
					C_larges[k][l] += pow(2./umin,l)*gsl_sf_fact(i)/gsl_sf_fact(i-l)*co_channels[k][j]*_T[j][i];
	}

}

Luminosity::dcomplex Luminosity::CLum_N(dcomplex N, unsigned channel) const {
	dcomplex ris(0.,0.);
	for (std::size_t i = 0 ; i < _n + 1 ; ++i)
		ris += C_larges[channel][i]/pow(N-1.,i+1.);
	return ris;
}
