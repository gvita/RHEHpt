#include "RHEHpt.h"
#include "LHAPDF/LHAPDF.h"
#include "Mfinite_HardPart_2.h"
#include "cuba.h"
#include <cmath>

using namespace std::placeholders;

namespace RHEHpt{

RHEHpt::RHEHpt(double CME, const std::string& PDFfile, double MH, double MT, double MB, unsigned choice,bool verbose):Exact_FO_ME(CME,MT,MH),_CME(CME),_mH(MH),_mT(MT),_mB(MB),_choice(choice),_verbose(verbose)
{
		_s = std::pow(CME,2);	 // GeV^2
	 _PDF = LHAPDF::mkPDF(PDFfile,0);	// create pdf choosing the reference replica from the set
	 _as = _PDF -> alphasQ(_mH);
	std::cout << _s << " " << _as << " " << _mH << " " << _mT << " "<< _mB << std::endl;		 
}


RHEHpt::~RHEHpt(){
	delete _PDF;
}


/******************************************/
/*		RESUMMED DISTRIBUTIONS				*/
/******************************************/

double RHEHpt::core_hpt_infinite(double N) const{
		const double sigma0 = _as * _as * Gf * std::sqrt(2.) / (576. * M_PI);
		const double as_N = _as / N;
//		const double M1 = _an_dim(as_N);
//		const double M2 = _an_dim(as_N);
		const double M1 = N;
		const double M2 = M1;	
		std::cout << "M=" << M1 << std::endl;
		const double gamma_factor = tgamma( 1. + M1 ) * tgamma( 1. + M2 ) * tgamma( 2. - M1 - M2 ) / ( tgamma(2.-M1) * tgamma(2.-M2) * tgamma(M1 + M2 ) );
//		return R(as_N) * R(as_N) * k * xp_factor * gamma_factor * ( 1 + (2 * M1 * M2)/(1 - M1 - M2));
	return	sigma0 * gamma_factor * ( 1. + (2. * M1 * M2)/(1. - M1 - M2));
}

double RHEHpt::hpt_infinite(double xp,double N) const{
//		const double xp = (pt*pt) / (_mH*_mH);
//		const double xp_factor = std::pow(xp, 2.*_an_dim(_as/N)	- 1. );
		const double xp_factor = std::pow(xp, 2.*N	- 1. );
		return core_hpt_infinite(N)*xp_factor;
}

int _core_hptfinite(const int *ndim, const double x[], const int *ncomp, double res[], void *pars){
	res[0] = 0.;
	double * p = (double *) pars;
	const long double xp = p[0];
	const long double yt = p[1];
	const long double M = p[2];
	
	const long double x1 = 0.25*x[0];
	const long double sqrtx1 = std::sqrt(x1);
	const long double r = x[1];
	const long double x2 = 1. + x1 + 2.* sqrtx1 *(1.-2.*r);
	const long double x2_red = 2.25 - 2.*r; // x2(1/4,r)
	long double DF0,F0,F0red;
	DF0 = D_F_0( x1, x2, xp, yt, r );		
	F0 = F_0( x1, x2, xp, yt, r);
	F0red = F_0(0.25,x2_red,xp,yt,r);
	const long double derivative_term = std::pow(x1*x2,M)/x2	* ( DF0 - (1.-M) * F0* (1.+sqrtx1 - 2.*r)/(sqrtx1*x2) );

	const long double surface_term = std::pow(x2_red/4.,M-1)*F0red;

	/** P1 Total **/	
	res[0] = 0.25*(surface_term - derivative_term)/std::sqrt(r*(1.-r));

	// and P_2, i.e. the term not integrated by parts
	const long double x1p = 0.25/x[0];
	const long double x2p = x[1] + x[1]/std::sqrt(x[0]) + x1p;
	const long double F_inf = F_Infinitive(x1p,x2p,xp,yt);
	res[0] += M * (std::pow( x1p*x2p,M)/(x[0]*x2p)) * (1. + 2.*sqrt(x1p) )*F_inf/std::sqrt(2.*x1p*x2p + 2.*x1p + 2.*x2p - x1p*x1p - x2p*x2p -1.);
		return 0;
}

//FIXME aggiungere massa bottom a private classe e cambiare F_costant

double RHEHpt::hpt_finite(double xp, double N) const
{
	double M = _an_dim(_as/N);
	double p[3] = { xp, get_yt(), M};
	double the_integral[1], error[1], prob[1];
	double epsrel=1.e-3, epsabs=1.e-15;
	int last = 4;
	int verbose = 0;
	int nregions, neval, fail;
		Cuhre(2, 1, _core_hptfinite, &p, 1,
			epsrel, epsabs, verbose | last,
			0, 500000, 9,
			NULL, NULL,
			&nregions, &neval, &fail, the_integral, error, prob);
	const double F_constant = 4608.*std::pow(M_PI,3.);
	std::cout << "hpt_finite integral (" << xp << ") = " << the_integral[0] << " ± " << error[0] << "\tp = " << prob[0] << std::endl;
	std::cout << std::endl;

		const double sigma0 = _as * _as * Gf * std::sqrt(2.) / (576. * M_PI);

		return sigma0 * M * std::pow(xp,2.*M-1.) * F_constant * 2. * the_integral[0];

}

/******************************************/
/*		RESUMMED Hpt EXPANSION				*/
/******************************************/


unsigned int factorial(unsigned int n)
{
	return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

double I_1_2( RHEHpt::_parameters* par);
double I_3( RHEHpt::_parameters* par);
double half_C( RHEHpt::_parameters* par){
	return (I_1_2(par) + I_3(par))/(factorial(par->j -1) * factorial(par -> k));
	
}

double RHEHpt::C(unsigned j, unsigned k,double xp) const
{
	const double yt = (_mT*_mT) / (_mH*_mH);
	const double yb = (_mB*_mB) / (_mH*_mH);
	const double F_constant = 4608.*std::pow(M_PI,3.);
	if ( k*j != 0 ) {
		_parameters par;
		switch (_choice){
			case(0): { // total
				par.set(j,k,xp,
						_parameters::HardF1( std::bind( &F_0_tot, _1, _2, _3, yt, yb, _4)) ,
						_parameters::HardF1( std::bind( &D_F_0_tot, _1, _2, _3, yt, yb, _4)) ,
						_parameters::HardF2( std::bind( &F_Infinitive_tot, _1, _2, _3, yt, yb)) 
				);
				break;
			}
			case(1): { // only top
				par.set(j,k,xp,
						_parameters::HardF1( std::bind( &F_0, _1, _2, _3, yt, _4)) ,
						_parameters::HardF1( std::bind( &D_F_0, _1, _2, _3, yt, _4)) ,
						_parameters::HardF2( std::bind( &F_Infinitive, _1, _2, _3, yt))
				);
				break; 
			}
			case(2): { // only bottom
				par.set(j,k,xp,
						_parameters::HardF1( std::bind( &F_0, _1, _2, _3, yb, _4)) ,
						_parameters::HardF1( std::bind( &D_F_0, _1, _2, _3, yb, _4)) ,
						_parameters::HardF2( std::bind( &F_Infinitive, _1, _2, _3, yb))
				);
				break;
			}
			case(3): { // interference
				par.set(j,k,xp,
						_parameters::HardF1( std::bind( &F_0_inter, _1, _2, _3, yt, yb, _4)) ,
						_parameters::HardF1( std::bind( &D_F_0_inter, _1, _2, _3, yt, yb, _4)) ,
						_parameters::HardF2( std::bind( &F_Infinitive_inter, _1, _2, _3, yt, yb))
				);
				break;
			}
		}
		if( j != k ){
			double C_j_k = half_C(&par);
			par.switch_jk();
			double C_k_j = half_C(&par);
			return F_constant * (C_j_k + C_k_j);
		}
		else
			return 2. * F_constant * half_C(&par);
		
	}
	else if ( k + j > 1.) return 0.;
	else switch (_choice){
		case(0): {
			return M_PI * F_constant * std::norm( A1x_0(xp,yt) + A1x_0(xp,yb) );
			break;
		}
		case(1): {
			return M_PI * F_constant * std::norm(A1x_0(xp,yt)); // checked against c0 infinite
			break; 
		}
		case(2): {
			return M_PI* F_constant * std::norm(A1x_0(xp,yb));
			break;
		}
		case(3): {
			return M_PI * F_constant*real(A1x_0(xp,yt)*std::conj(A1x_0(xp,yb))+std::conj(A1x_0(xp,yt))*A1x_0(xp,yb));
			break;
		}
	}
}


std::vector< std::array<unsigned,2> > RHEHpt::partition_of(unsigned n) const
{
	/* returns all the couples of int (x,y) s.t. x + y = n and x >= y	*/
	std::vector < std::array<unsigned,2> > result;
	const std::size_t max_i = n/2; // deliberate integer division
	for(unsigned i = 1; i <= max_i ; ++i) 
		result.push_back({n-i,i});
	return result;
}


std::vector< std::array<unsigned,2> > RHEHpt::symmetric_partition_of(unsigned n) const
{
	/* returns all the couples of int (x,y) s.t. x + y = n */
	std::vector < std::array<unsigned,2> > result;
	const std::size_t max_i = n/2; // deliberate integer division
	for(unsigned i = 0; i <= max_i ; ++i){
		result.push_back({n-i,i});
		if(n-i != i)
			result.push_back({i,n-i});
	}
	return result;
}

std::vector<double> RHEHpt::Integral_coeffs(unsigned order, double xp, bool HQ) const
{ 
	if(order < 1 ) return std::vector<double>(1,0.); 
	std::vector<double> coeff_list( order + 1 );

	coeff_list[0] = 0;
	if ( HQ ){
		coeff_list[1] = 2.;
		if (order == 1) return coeff_list;

		coeff_list[2] = 0.;
		if (order == 2) return coeff_list;
		
		coeff_list[3] = 2.;
		if (order == 3) return coeff_list;
		std::cerr << "N3LO for heavy_quark not implemented." << std::endl;
		return std::vector<double>(1,0.);
	}
	coeff_list[1] = C(1,0,xp);

	if(order > 1){
		for(unsigned i = 2 ; i < order + 1 ; ++i){
			// get couples of int (x,y) s.t. x + y = i and x >= y
			std::vector< std::array<unsigned,2> > indices = partition_of(i); 
			std::vector<double> tmp(indices.size());
			std::cout << "number of couples at order " << i << " = " << indices.size() << std::endl;
			// compute all coefficients of the i-th order (for i=4 they are c_{1,3} and c_{2,2}), store result in tmp
			std::transform(indices.cbegin(),indices.cend(),tmp.begin(),[&](const std::array<unsigned,2>& index_pair){
				std::cout << "computing C(" << index_pair[0] << "," << index_pair[1] << ",xp)" << std::endl;
				return C(index_pair[0],index_pair[1],xp); 
			});
			coeff_list[i] = 2.*std::accumulate(tmp.cbegin(),tmp.cend(),0.); // sum of the partial coefficients. 2 comes from M1=M2 symm
		}
	}
	return coeff_list; // coeff_list[i] should go with M^i
}

std::vector < double > RHEHpt::Xp_prefactor_coeffs(unsigned order, double xp) const
{
	// compute the analytic expansion of xp^(2 M-1) up to the n-1-th order (cause the O(M^0) term of Integral_coeffs is 0)
	std::vector<double> xp_coeff_list( order +1 ); // this starts at M^0
	for(unsigned i = 0 ; i < order+1.	; ++i)
		xp_coeff_list[i] = std::pow(2.*std::log(xp),i)/factorial(i)/xp;
	return xp_coeff_list;
}


std::vector<double> RHEHpt::pt_distr_series_terms(unsigned order, double xp,double N,bool HQ) const
{
/*
	IMPORTANT NOTE: This is formally correct up to order 4 where you need to expand also R(as_N).
	Also at order 4, O(_as) != O(M) since M ~ _as + O(_as^4)
*/

	//FIXME need to add R(M)^2
	if(order < 1 ) return std::vector<double>(1,0.); 
	
	// compute the analytic expansion of xp^(2 M-1) up to the n-1-th order (cause the O(M^0) term of Integral_coeffs is 0)
	std::vector<double> xp_coeff_list( order +1 ); // this starts at M^0
	for(unsigned i = 0 ; i < order+1.	; ++i)
		xp_coeff_list[i] = std::pow(2.*std::log(xp),i)/factorial(i)/xp;
		
	// compute the coeff_list for the double integral
	// this starts at M^0. This is the longest part of the computation
	std::vector<double> integral_coeff_list = Integral_coeffs(order,xp,HQ);

	if( _verbose ){
		/* print stuff */
		std::cout << "integral_coeff_list " << std::endl;
		for (auto it : integral_coeff_list){
			std::cout << it << " " ;
		}
		std::cout << std::endl;
	
		std::cout << "xp_coeff_list " << std::endl;
		for (auto it : xp_coeff_list){
			std::cout << it << " " ;
		}
		std::cout << std::endl;
	}	
	return pt_distr_series_terms(N,integral_coeff_list,xp_coeff_list,order);	

}

double RHEHpt::pt_distr_series(unsigned order, double xp,double N, bool HQ) const
{
		const double sigma0 = _as * _as * Gf * std::sqrt(2.) / (576. * M_PI);
		std::vector<double> terms = pt_distr_series_terms(order, xp, N, HQ);
	return sigma0 * std::accumulate(terms.cbegin(),terms.cend(),0.);
}


/*************************** I_i integrals definitions ************************/


/************************************* I_1_2 *************************************/

int _core_I_1_2(const int *ndim, const double x[], const int *ncomp, double res[], void *pars){
	res[0] = 0.;
	if(x[0] < 0.25){
		// I_2 calculation
		RHEHpt::_parameters * p = (RHEHpt::_parameters *) pars;
		const long double xp = p -> xp;
		const unsigned j = p -> j;
		const unsigned k = p -> k;
		RHEHpt::_parameters::HardF1 F_0f = p -> F_0f;
		RHEHpt::_parameters::HardF1 D_0f = p -> D_0f;
		const long double x1 = x[0];
		const long double sqrtx1 = std::sqrt(x1);
		const long double r = x[1];
		const long double x2 = 1. + x1 + 2.* sqrtx1 *(1.-2.*r);
		const long double logjx2 = std::pow(std::log(x2),j-1);

		const long double x2_red = 2.25 - 2.*r; // x2(1/4,r)
		long double DF0,F0,F0red;
		DF0 = D_0f( x1, x2, xp, r );
		F0 = F_0f( x1, x2, xp, r);
		F0red = F_0f(0.25, x2_red, xp, r);
		
		const long double C1 = (logjx2 / x2) * DF0 ;
		//const long double C1=0.0;
		long double C2 = logjx2 * F0* (1.+sqrtx1 - 2.*r)/(sqrtx1*x2*x2);
		//long double C2=0.0;
		if( j > 1 ) 
			C2 -= (j-1.) * std::pow(std::log(x2),j-2) * F0* (1.+sqrtx1 - 2.*r)/(sqrtx1*x2*x2);
		
		// I_1 calculation

		const long double I1_numerator = std::pow(std::log(0.25),k)*std::pow(std::log(x2_red),j-1)*F0red/x2_red;
		
		/** Total **/
		
		// 4*I1_numerator because I_1 is an integral on dr, here we are integrating in dr dx	
		res[0] = (4. * I1_numerator -	(C1 - C2 )*std::pow(std::log(x1),k) )/std::sqrt(r*(1.-r));
		/*if (res[0]!=res[0]){
		  std::cout << res[0] << " YES " << x[0]<< " " << x[1] << " " << x2 << " " << x2_red<< " " << DF0<< std::endl;
		}*/
		
	}
	return 0;
}

double I_1_2(RHEHpt::_parameters* par){
	double the_integral[1], error[1], prob[1];
	double epsrel=1.e-3, epsabs=1.e-15;
	int last = 4;
	int verbose = 0;
	int nregions, neval, fail;
		Cuhre(2, 1, _core_I_1_2, par, 1,
			epsrel, epsabs, verbose | last,
			0, 500000, 9,
			NULL, NULL,
			&nregions, &neval, &fail, the_integral, error, prob);
	std::cout << "I_2 integral (" << par -> j << "," << par -> k << "," << par -> xp << ") = " << the_integral[0] << " ± " << error[0] << "\tp = " << prob[0] << std::endl;
	std::cout << std::endl;
		return the_integral[0];
}


/************************************* I_3 *************************************/

int _core_I_3(const int *ndim, const double x[], const int *ncomp, double res[], void *pars){ // TO DO rescale x1 to 1/4 -> 1 and put this inside _core_I_1_2
	res[0] = 0.;
	RHEHpt::_parameters * p = (RHEHpt::_parameters *) pars;
	const long double xp = p -> xp;
	const unsigned j = p -> j;
	const unsigned k = p -> k;
	RHEHpt::_parameters::HardF2 F_inf = p -> F_inf;
			
	const long double x1 = 0.25/x[0];
	const long double x2 = x[1] + x[1]/std::sqrt(x[0]) + x1;
	res[0] = k * ( std::pow(std::log(x1),k-1) * std::pow(std::log(x2),j-1)/(x[0] * x2))
	*(1. + 2.*sqrt(x1) )*F_inf(x1,x2,xp)/std::sqrt(2.*x1*x2 + 2.*x1 + 2.*x2 - x1*x1 - x2*x2 -1.);
	return 0;
}

double I_3(RHEHpt::_parameters* par){
	double the_integral[1], error[1], prob[1];
	double epsrel=1.e-3, epsabs=1.e-15;
	int last = 4;
	int verbose = 0;
	int nregions, neval, fail;
		Cuhre(2, 1, _core_I_3, par, 1,
			epsrel, epsabs, verbose | last,
			0, 500000, 9,
			NULL, NULL,
			&nregions, &neval, &fail, the_integral, error, prob);
	std::cout << "I_3 integral (" << par -> j << "," << par -> k << "," << par -> xp << ") = " << the_integral[0] << " ± " << error[0] << "\tp = " << prob[0] << std::endl;
	std::cout << std::endl;
		return the_integral[0];
}


/******************************************/
/*		FIXED ORDER DISTRIBUTION			*/
/******************************************/

double RHEHpt::ptd_FO(double xp)
{
	const double MatrixElementFO = Exact_FO_ME(xp);
	const double K = 1./(16.*M_PIl*256.) * Gf * sqrt(2.)*_as*_as*_as;
	return K *0.5 * 3./M_PIl * MatrixElementFO;
}

}
