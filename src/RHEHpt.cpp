#include "RHEHpt.h"
#include "LHAPDF/LHAPDF.h"
#include "Mfinite_HardPart_2.h"
#include "cuba.h"
#include <cmath>

using namespace std::placeholders;

namespace RHEHpt{

RHEHpt::RHEHpt(double CME, const std::string& PDFfile, double MH, double MT, double MB, double MUR, double MUF, 
				unsigned choice,unsigned channel,bool verbose):
					_PDF(LHAPDF::mkPDF(PDFfile,0)),
		 			Exact_FO_fullmass(CME,0.1,MT,MB,MH,choice),Exact_FO_PL(CME,0.1,MH,MUF,MUR),
		 			_CME(CME), _mH(MH), _mT(MT), _mB(MB), _muR(MUR), _muF(MUF), _choice(choice), _channel(channel), _verbose(verbose),
		 			_Lum(_PDF), _Integral_coeff_interpolator(choice)
{
	_s = std::pow(CME,2);	 // GeV^2
	_as=_PDF -> alphasQ(_muR);
	_Lum.setNf(5);
	_Lum.Cheb_Lum( _muF );
	Exact_FO_fullmass.SetAS(_as);
	Exact_FO_PL.SetAS(_as);
	

	std::cout << _s << " " << _as << " " << _mH << " " << _mT << " "<< _mB << " " << _muR << " " << _muF << std::endl;		 
}


RHEHpt::~RHEHpt(){
	std::cout << "deleting RHEHpt" << std::endl;
	delete _PDF;
}


/******************************************/
/*		RESUMMED DISTRIBUTIONS				*/
/******************************************/

double RHEHpt::core_hpt_infinite(double N) const{
//		const double sigma0 = _as * _as * Gf * std::sqrt(2.) / (576. * M_PI);
		const double as_N = _as / N;
//		const double M1 = _an_dim(as_N);
//		const double M2 = _an_dim(as_N);
		const double M1 = N;
		const double M2 = M1;	
		std::cout << "M=" << M1 << std::endl;
		const double gamma_factor = tgamma( 1. + M1 ) * tgamma( 1. + M2 ) * tgamma( 2. - M1 - M2 ) / ( tgamma(2.-M1) * tgamma(2.-M2) * tgamma(M1 + M2 ) );
//		return R(as_N) * R(as_N) * k * xp_factor * gamma_factor * ( 1 + (2 * M1 * M2)/(1 - M1 - M2));
	return	SIGMA0() * gamma_factor * ( 1. + (2. * M1 * M2)/(1. - M1 - M2));
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

//	const double sigma0 = _as * _as * Gf * std::sqrt(2.) / (576. * M_PI);

	return SIGMA0() * M * std::pow(xp,2.*M-1.) * F_constant * 2. * the_integral[0];

}

/******************************************/
/*		RESUMMED Hpt EXPANSION			  */
/******************************************/


unsigned int factorial(unsigned int n)
{
	return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

long double RHEHpt::M(long double x,unsigned i) const{ //inverse mellin of _an_dim^i at x
	return std::pow(3.*_as/M_PIl,i)*std::pow(std::log(1./x),i-1)/factorial(i-1);
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

double RHEHpt::_Integral_coeff_from_integrals(unsigned i, double xp) const
{
	// get couples of int (x,y) s.t. x + y = i and x >= y
	std::vector< std::array<unsigned,2> > indices = partition_of(i); 
	std::vector<double> tmp(indices.size());
	std::cout << "number of couples at order " << i << " = " << indices.size() << std::endl;
	// compute all coefficients of the i-th order (for i=4 they are c_{1,3} and c_{2,2}), store result in tmp
	std::transform(indices.cbegin(),indices.cend(),tmp.begin(),[&](const std::array<unsigned,2>& index_pair){
		std::cout << "computing C(" << index_pair[0] << "," << index_pair[1] << ",xp)" << std::endl;
		return C(index_pair[0],index_pair[1],xp); 
	});
	
	return 2.*std::accumulate(tmp.cbegin(),tmp.cend(),0.); // sum of the partial coefficients. 2 comes from M1=M2 symm
}

double RHEHpt::_Integral_coeff_from_file(unsigned order,double xp) const
{
	if (order > 3) return _Integral_coeff_from_integrals(order,xp); // N3LO and beyond not precomputed
	else if(order > 1) return _Integral_coeff_interpolator(order,xp);
	else if (order == 1) return C(1,0,xp);
	else return 0;
}


std::vector<double> RHEHpt::Integral_coeffs(unsigned order, double xp, bool HQ) const
{ 
	if(order < 1 ) return std::vector<double>(1,0.); 
	std::vector<double> coeff_list( order + 1 );

	coeff_list[0] = 0; // as^2

	if ( HQ ){
		coeff_list[1] = 2.; // as^3
		if (order == 1) return coeff_list;

		coeff_list[2] = 0.; // as^4
		if (order == 2) return coeff_list;
		
		coeff_list[3] = 2.; // as^5
		if (order == 3) return coeff_list;
		std::cerr << "N3LO not implemented." << std::endl;
		return std::vector<double>(1,0.);
	}
	
	/**** finite mass case ****/
	
	coeff_list[1] = C(1,0,xp);

	if(order > 1){
		const double yt = get_yt();
		const double yb = get_yb();
		
		/* If yb ~ yb_def and yt ~ yt_def use precomputed coefficients.
		   We give a 1% tolerance in the check, since, experimentally there is about a 1% uncertainty on the quark/higgs mass.*/
		if (std::abs(yb - yb_def)/yb_def < 0.01 && std::abs(yt - yt_def)/yt_def < 0.01 ){
			for(unsigned i = 2 ; i < order + 1 ; ++i)
				coeff_list[i] = _Integral_coeff_from_file(i,xp);
		}
		else{
			for(unsigned i = 2 ; i < order + 1 ; ++i)
				coeff_list[i] = _Integral_coeff_from_integrals(i,xp);
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


std::vector<double> RHEHpt::pt_distr_series_terms(unsigned order, double xp,double N,bool HQ,bool Nspace) const
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
	return pt_distr_series_terms< double >(N,integral_coeff_list,xp_coeff_list,order,Nspace);	

}

double RHEHpt::pt_distr_series(unsigned order, double xp,double N, bool HQ, bool Nspace) const
{
	std::vector<double> terms = pt_distr_series_terms(order, xp, N, HQ, Nspace);
	return SIGMA0() * std::accumulate(terms.cbegin(),terms.cend(),0.);
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

//CORE FUNCTION

double NLO_PL_sing_function(double x1,void *p){
  RHEHpt::_par_part *pars=(RHEHpt::_par_part *)p;
  double ris=0.;
  switch (pars->_channel_int){
    case(1):{
      ris=(pars->_Pointlike_int).NLO_PL_sing_doublediff(pars->_xp_int,x1);
      break;
    }
    case(2):{
      ris=(pars->_Pointlike_int).NLO_PL_sing_doublediff_gq(pars->_xp_int,x1);
      break;
    }
    case(3):{
      ris=(pars->_Pointlike_int).NLO_PL_sing_doublediff_qqbar(pars->_xp_int,x1);
      break;
    }
    case(4):{
      ris=0.;
      break;
    }
    case(5):{
      ris=0.;
      break;
    }
    case(6):{
      ris=0.;
      break;
    }
  }
  return(ris);
}
double NLO_PL_notsing_function(double x1,void *p){
  RHEHpt::_par_part *pars=(RHEHpt::_par_part *)p;
  double ris=0.;
  switch (pars->_channel_int){
    case(1):{
      ris=(pars->_Pointlike_int).NLO_PL_notsing_doublediff(pars->_xp_int,x1);
      break;
    }
    case(2):{
      ris=(pars->_Pointlike_int).NLO_PL_notsing_doublediff_gq(pars->_xp_int,x1);
      break;
    }
    case(3):{
      ris=(pars->_Pointlike_int).NLO_PL_notsing_doublediff_qqbar(pars->_xp_int,x1);
      break;
    }
    case(4):{
      ris=(pars->_Pointlike_int).NLO_PL_notsing_doublediff_qq(pars->_xp_int,x1);
      break;
    }
    case(5):{
      ris=(pars->_Pointlike_int).NLO_PL_notsing_doublediff_qqbarprime(pars->_xp_int,x1);
      break;
    }
    case(6):{
      ris=(pars->_Pointlike_int).NLO_PL_notsing_doublediff_qqprime(pars->_xp_int,x1);
      break;
    }
  }
  return(ris);
}


//--------------------------------------------------------------------------------------------------------------------------------------------------------------------

long double RHEHpt::sigma_part(long double CME_part,long double pt, unsigned int order, bool heavyquark){
	long double sigma = 0.0;
	long double sigma0 = SIGMA0();
	if (_channel==0){
	  std::cout<< "TOT channel is not a partonic channel: switch to gg channel major contribution" << std::endl;
	  _channel=1;
	}
	if ( heavyquark == false ){
		Exact_FO_fullmass.setchannel(_channel);
		Exact_FO_fullmass.SetCME(CME_part);
		sigma = sigma0 * Exact_FO_fullmass(get_xp(pt));
	}
	else {
		Exact_FO_PL.setchannel(_channel);
		if ( order > 1 ){
			std::cout << "Error order must be or 0 (LO) or 1(NLO)" << std::endl;
			return 0;
		}
		if ( order == 0 ){
			Exact_FO_PL.SetCME(CME_part);
			sigma = sigma0 * Exact_FO_PL.LO_PL(get_xp(pt));
		}
		if( order == 1 ){
			Exact_FO_PL.SetCME(CME_part);
			_par_part par(Exact_FO_PL,get_xp(pt),_channel);
			double precision = 1e-5;
			double NLO_PL_sing_ris = 0.0, NLO_PL_sing_error = 0.0;
			gsl_function NLO_PL_sing;
			gsl_integration_workspace * w = gsl_integration_workspace_alloc (100000);
			NLO_PL_sing.function = NLO_PL_sing_function;
			NLO_PL_sing.params = &par;
			gsl_integration_qags (&NLO_PL_sing, 0.0, 1.0, 0, precision, 100000,	w, &NLO_PL_sing_ris, &NLO_PL_sing_error);
			gsl_integration_workspace_free (w);
			double NLO_PL_notsing_ris = 0.0, NLO_PL_notsing_error = 0.0;
			gsl_function NLO_PL_notsing;
			w = gsl_integration_workspace_alloc (100000);
			NLO_PL_notsing.function = NLO_PL_notsing_function;
			NLO_PL_notsing.params = &par;
			gsl_integration_qags (&NLO_PL_notsing, 0.0, 1.0, 0, precision, 100000, w, &NLO_PL_notsing_ris, &NLO_PL_notsing_error);
			gsl_integration_workspace_free (w);
			long double NLO_PL_delta_ris= 0.;
			switch(_channel){
			  case(1):{
			    NLO_PL_delta_ris=Exact_FO_PL.NLO_PL_delta(get_xp(pt));
			    break;
			  }
			  case(2):{
			    NLO_PL_delta_ris=Exact_FO_PL.NLO_PL_delta_gq(get_xp(pt));
			    break;
			  }
			  case(3):{
			    NLO_PL_delta_ris=Exact_FO_PL.NLO_PL_delta_qqbar(get_xp(pt));
			    break;
			  }
			  case(4):{
			    NLO_PL_delta_ris=0.;
			    break;
			  }
			  case(5):{
			    NLO_PL_delta_ris=0.;
			    break;
			  }
			  case(6):{
			    NLO_PL_delta_ris=0.;
			    break;
			  }
			}
			sigma = sigma0*(_as*_as/(4.*M_PIl*M_PIl)*(NLO_PL_delta_ris + NLO_PL_sing_ris + NLO_PL_notsing_ris));
			//sigma=(NLO_PL_delta_ris + NLO_PL_sing_ris + NLO_PL_notsing_ris);
			std::cout << "Sigma_part_NLO( "<< CME_part << " , " << pt << ")= " << sigma << " +- " << sigma0*(_as*_as/(4.*M_PIl*M_PIl)*( NLO_PL_notsing_error + NLO_PL_sing_error ))<< std::endl;
		}
	}
	return sigma;
}

double rho(double xp){
	return std::pow( std::sqrt(1.+xp) + std::sqrt(xp) ,2);
}

double LO_fullmass_function(double z,void *p){
	RHEHpt::_par_hadro *pars=(RHEHpt::_par_hadro *)p;
	long double tauprime=(pars->_tau_int)* rho(pars -> _xp_int);
	unsigned int channel=(pars->_channel_int);
	pars->_Finite_int.SetCME( pars -> _mH_int / std::sqrt(z) * std::sqrt( rho(pars -> _xp_int) ) );
	long double ris=0.0;
	switch (channel){
	  case(0):{
	    (pars->_Finite_int).setchannel(1);
	    ris+= 1./z*pars->_Lumi_int.CLum_gg_x(tauprime/z)*pars->_Finite_int(pars->_xp_int)/z;
	    (pars->_Finite_int).setchannel(2);
	    ris+= 1./z*pars->_Lumi_int.CLum_gq_x(tauprime/z)*pars->_Finite_int(pars->_xp_int)/z;
	    (pars->_Finite_int).setchannel(3);
	    ris+= 1./z*pars->_Lumi_int.CLum_qqbar_x(tauprime/z)*pars->_Finite_int(pars->_xp_int)/z;
	    break;
	  }
	  case(1):{
	    (pars->_Finite_int).setchannel(1);
	    ris+= 1./z*pars->_Lumi_int.CLum_gg_x(tauprime/z)*pars->_Finite_int(pars->_xp_int)/z;
	    break;
	  }
	  case(2):{
	    (pars->_Finite_int).setchannel(2);
	    ris+= 1./z*pars->_Lumi_int.CLum_gq_x(tauprime/z)*pars->_Finite_int(pars->_xp_int)/z;
	    break;
	  }
	  case(3):{
	    (pars->_Finite_int).setchannel(3);
	    ris+= 1./z*pars->_Lumi_int.CLum_qqbar_x(tauprime/z)*pars->_Finite_int(pars->_xp_int)/z;
	    break;
	  }
	}
	return (ris);
}

double LO_PL_function(double z,void *p){
	RHEHpt::_par_hadro *pars=(RHEHpt::_par_hadro *)p;
	long double tauprime = (pars->_tau_int) * rho(pars -> _xp_int);
	unsigned int channel=(pars->_channel_int);
	pars->_Pointlike_int.SetCME(pars->_mH_int/std::sqrt(z)* std::sqrt( rho(pars -> _xp_int) ) );
	long double ris=0.;
	switch(channel){
	  case(0):{
	    (pars->_Pointlike_int).setchannel(1);
	    ris+=(1./z) * pars->_Lumi_int.CLum_gg_x(tauprime/z) * pars->_Pointlike_int.LO_PL(pars->_xp_int)/z;
	    (pars->_Pointlike_int).setchannel(2);
	    ris+=(1./z) * pars->_Lumi_int.CLum_gq_x(tauprime/z) * pars->_Pointlike_int.LO_PL(pars->_xp_int)/z;
	    (pars->_Pointlike_int).setchannel(3);
	    ris+=(1./z) * pars->_Lumi_int.CLum_qqbar_x(tauprime/z) * pars->_Pointlike_int.LO_PL(pars->_xp_int)/z;
	    break;
	  }
	  case(1):{
	    (pars->_Pointlike_int).setchannel(1);
	    ris+=(1./z) * pars->_Lumi_int.CLum_gg_x(tauprime/z) * pars->_Pointlike_int.LO_PL(pars->_xp_int)/z;
	    break;
	  }
	  case(2):{
	    (pars->_Pointlike_int).setchannel(2);
	    ris+=(1./z) * pars->_Lumi_int.CLum_gq_x(tauprime/z) * pars->_Pointlike_int.LO_PL(pars->_xp_int)/z;
	    break;
	  }
	  case(3):{
	    (pars->_Pointlike_int).setchannel(3);
	    ris+=(1./z) * pars->_Lumi_int.CLum_qqbar_x(tauprime/z) * pars->_Pointlike_int.LO_PL(pars->_xp_int)/z;
	    break;
	  }
	}
	return (ris);
}

int _core_NLO_PL_notsing(const int *ndim, const double x[], const int *ncomp, double res[], void *pars){
	res[0] = 0.;
	RHEHpt::_par_hadro *p = (RHEHpt::_par_hadro *) pars;
	long double tauprime = ( p -> _tau_int ) * rho(p -> _xp_int);

	//Setting integration limit
	long double z = tauprime+(1.-tauprime)*x[0];
	p->_Pointlike_int.SetCME(p->_mH_int/std::sqrt(z)* std::sqrt( rho(p -> _xp_int) ) );
	res[0]+=( (1./z) * p->_Lumi_int.CLum_gg_x(tauprime/z) * p->_Pointlike_int.NLO_PL_notsing_doublediff(p->_xp_int,x[1]) / z );
	return 0;
}

int _core_NLO_PL_sing(const int *ndim, const double x[], const int *ncomp, double res[], void *pars){
	res[0] = 0.;
	RHEHpt::_par_hadro *p = (RHEHpt::_par_hadro *) pars;
	long double tauprime = (p->_tau_int) * rho(p -> _xp_int) ;

	//Setting integration limit
	long double z = tauprime+(1.-tauprime)*x[0];
	p->_Pointlike_int.SetCME(p->_mH_int/std::sqrt(z)* std::sqrt( rho(p -> _xp_int) ) );
	res[0] += ( (1./z) * p->_Lumi_int.CLum_gg_x(tauprime/z) * p->_Pointlike_int.NLO_PL_sing_doublediff(p->_xp_int,x[1]) / z );
	return 0;
}

double NLO_delta_function(double z,void *p){
  RHEHpt::_par_hadro *pars=(RHEHpt::_par_hadro *)p;
  long double tauprime = (pars->_tau_int) * rho(pars -> _xp_int) ;
  pars->_Pointlike_int.SetCME( pars -> _mH_int / std::sqrt(z) * std::sqrt( rho(pars -> _xp_int) ) );
  return (1./z*pars->_Lumi_int.CLum_gg_x(tauprime/z)*pars->_Pointlike_int.NLO_PL_delta(pars->_xp_int)/z);
}

double pt_hadro(double z,void *p){
  RHEHpt::_par_expansion *pars=(RHEHpt::_par_expansion *)p;
  long double tauprime = (pars->_tau_int) * rho(pars -> _xp_int) ;
  
  return (1./z) * pars->_Lumi_int.CLum_gg_x(tauprime/z) * pars -> ds_xp_part( z )/z;
  
}



long double RHEHpt::sigma_hadro_FO_fullmass(long double pt){
	long double sigma = 0.0;
	long double sigma_error = 0.0;
	std::string name;	// A string to identify what has been computed, i.e. a human readable equivalent of order and heavyquark

	_par_hadro par(_Lum,Exact_FO_fullmass,Exact_FO_PL, get_xp(pt),get_tau(),_mH,_channel);
	long double tauprime=(get_tau()) * rho( get_xp(pt) ) ;
	double precision = 1e-6;
	double LO_fullmass_ris = 0.0, LO_fullmass_error = 0.0;
	gsl_function LO_fullmass;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (100000);
	LO_fullmass.function = LO_fullmass_function;
	LO_fullmass.params = &par;
	gsl_integration_qags (&LO_fullmass, tauprime, 1.0, 0, precision, 100000,	w, &LO_fullmass_ris, &LO_fullmass_error);
	gsl_integration_workspace_free (w);
	sigma = LO_fullmass_ris * tauprime * gev2_to_pb;
	sigma_error = LO_fullmass_error*tauprime*gev2_to_pb;
	name = "Sigma_hadro_LO_fullmass";
	sigma *= SIGMA0();
	sigma_error *= SIGMA0();

	std::cout << name + "( " << _CME << " , " << pt << ")= " << sigma << " ± " << sigma_error << std::endl;
	return sigma;

}

std::vector<double> RHEHpt::sigma_hadro_FO_pointlike(std::vector<double>& ptgrid, unsigned int order){
//	int gridsize = ptgrid.size();
//	std::vector< double > result( gridsize, 0. );
	std::vector< double > result;
	if (order == 1 ){
		for (auto pt : ptgrid){
			_par_hadro par(_Lum,Exact_FO_fullmass,Exact_FO_PL, get_xp(pt),get_tau(),_mH,_channel);
			long double tauprime=(get_tau()) * rho( get_xp(pt) ) ;
			double precision = 1e-6;
			double LO_PL_ris = 0.0, LO_PL_error = 0.0;
			gsl_function LO_PL;
			gsl_integration_workspace * w = gsl_integration_workspace_alloc (100000);
			LO_PL.function = LO_PL_function;
			LO_PL.params = &par;
			gsl_integration_qags (&LO_PL, tauprime, 1.0, 0, precision, 100000,      w, &LO_PL_ris, &LO_PL_error);
			gsl_integration_workspace_free (w);
			result.push_back( SIGMA0() * LO_PL_ris * tauprime * gev2_to_pb );
		}
		return result;
	}
	// in this way it always returns a vector of size ptgrid. In caso of an error or a choice of order == 0, the vector is all zero.
	int gridsize = ptgrid.size();
	result.resize( gridsize, 0. );
	if ( order == 2){
	   	if (_channel > 4){
		  std::cout << "Error channel must be or 0 (all channel) or 1 (GG) or 2 (GQ) or 3 (QQ)" << std::endl;
		  return result;
		}
		_channel += 1; // channel in hqt starts from 1
		std::cout << _channel << std::endl;
		double CME = _CME;
		double MH = _mH;
		int len = ((_PDF -> set()).name()).length();
		int hqt_order = order;
		const char* pdfsetname = (_PDF -> set()).name().c_str();
		int pdfmem = 0;
		int channel=_channel;
		hqt_( &CME, &MH, &MH, &_muR, &_muR, &hqt_order, pdfsetname, &len, &pdfmem, &ptgrid[0],&result[0],&gridsize, &channel);
	   	return result;
	
	}
	else if ( order > 2 ){
			std::cout << "Error order must be or 1 (LO,as^3) or 2 (NLO,as^4)" << std::endl;
	 		return result;
	}
	return result; // 
}

long double RHEHpt::pt_distr_hadro(long double pt, unsigned int order, bool heavyquark){
	double xp = get_xp(pt);
	std::vector<double> int_coeff = Integral_coeffs( order, xp, heavyquark);
	std::vector<double> xp_coeff = Xp_prefactor_coeffs( order, xp);
	
	// the partonic cross section (at fixed xp, as a function of x)
	_par_expansion::Partonic_Distr ds_xp_part ( [&](double tau){ 
		return pt_distr_series(pt_distr_series_terms(tau, int_coeff, xp_coeff, order, false), order); // false --> not N space but x space
	}); 

	_par_expansion par( _Lum, ds_xp_part, xp , get_tau(), _mH);

	long double tauprime=(get_tau()) * rho( xp ) ;

	double precision = 1e-6;
	double ris = 0.0, error = 0.0;
	gsl_function f;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (100000);
	f.function = pt_hadro;
	f.params = &par;
	gsl_integration_qags (&f, tauprime, 1.0, 0, precision, 100000,	w, &ris, &error);
	gsl_integration_workspace_free (w);
	long double sigma = ris * tauprime * gev2_to_pb;
	long double sigma_error = error*tauprime*gev2_to_pb;
	return sigma;
}

std::vector<double> RHEHpt::pt_distr_hadro(const std::vector< double >& ptgrid, unsigned int order, bool heavyquark){
	std::vector< double > result;
	for (auto pt : ptgrid)
		result.push_back( pt_distr_hadro(pt, order, heavyquark) );
	return result;
}

}
