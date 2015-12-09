#include <complex>
#include <cmath>
#include <iostream>
#include "supporting_functions.h"
#include "gsl/gsl_sf_dilog.h"

typedef std::complex<long double> ldcomplex;
typedef std::complex<double> dcomplex;
typedef std::vector<double> Rvector;

//namespace RHEHpt{
//    
//long double B0(const long double x, const long double yt){ // CORRECT ONLY IF ARGUMENT IS < 0
//    const long double a = x/(x - 4.*yt);
//    const long double one = 1.;
//    return - std::sqrt( 1./a ) * std::log( (one + std::sqrt(a)) / (one - std::sqrt(a)) ) / (16*M_PIl*M_PIl);
//}


//long double B1(const long double x, const long double yt){ // CORRECT ONLY IF ARGUMENT IS > 0 && < 4. yt
//    const long double a = x/(4.*yt - x);
//    const long double one = 1.;
//    return - std::sqrt( 1./a ) * std::atan( std::sqrt(a) ) / ( 8. * M_PIl * M_PIl);
//}

//ldcomplex Li2(std::complex <long double> z){
//    gsl_sf_result result_re;
//    gsl_sf_result result_im;
//    gsl_sf_complex_dilog_e (std::abs(z), std::arg(z), &result_re , &result_im);
//    return ldcomplex (result_re.val,result_im.val);
//}

//long double Delta3(long double x1, long double x2){
//    return 1. + x1*x1 + x2*x2 - 2.*x1*x2 + 2.*x1 + 2.*x2;
//}

//long double C0(long double x1,long double x2, long double yt){
//    const long double sdelta = std::sqrt(Delta3(x1,x2));
//    const long double d1p = ( 1. + ( - x1 + x2 -1. ) / sdelta) / 2. ;
//    const long double d1m = ( 1. - ( - x1 + x2 -1. ) / sdelta) / 2. ;
//    const long double d2p = ( 1. + ( + x1 - x2 -1. ) / sdelta) / 2. ;
//    const long double d2m = ( 1. - ( + x1 - x2 -1. ) / sdelta) / 2. ;
//    const long double d3p = ( 1. + ( + x1 + x2 +1. ) / sdelta) / 2. ;
//    const long double d3m = ( 1. - ( + x1 + x2 +1. ) / sdelta) / 2. ;
//    const long double xp = - x2 * ( 1. + std::sqrt(1. + 4. * yt / x2)) / ( 2 * yt);
//    const long double xm = - x2 * ( 1. - std::sqrt(1. + 4. * yt / x2)) / ( 2 * yt);
//    const long double yp = - x1 * ( 1. + std::sqrt(1. + 4. * yt / x1)) / ( 2 * yt);
//    const long double ym = - x1 * ( 1. - std::sqrt(1. + 4. * yt / x1)) / ( 2 * yt);
//    ldcomplex I(0.,1.);
//    long double one = 1.; // stupid compiler cannot do 1. + I because 1. is double (and not long double)
//    const ldcomplex zp = ( one + I * std::sqrt<long double>( 4. * yt - 1.))/(2. * yt);
//    const ldcomplex zm = ( one - I * std::sqrt<long double>( 4. * yt - 1.))/(2. * yt);
//    const ldcomplex T1 = 
//        std::log(1. - ym) * std::log( (1.-ym*d1p)/(1.-ym*d1m) ) + 
//        std::log(1. - xm) * std::log( (1.-xm*d2p)/(1.-xm*d2m) ) +
//        std::log(one - zm) * std::log( (one-zm*d3p)/(one-zm*d3m) );
//    const ldcomplex T2 = 
//        Li2(yp*d1p) + Li2(ym*d1p) - Li2(yp*d1m) - Li2(ym*d1m) +
//        Li2(xp*d2p) + Li2(xm*d2p) - Li2(xp*d2m) - Li2(xm*d2m) +
//        Li2(zp*d3p) + Li2(zm*d3p) - Li2(zp*d3m) - Li2(zm*d3m);
////    if(std::imag(T1+T2)/std::real( T1 + T2 ) > 1E-3) std::cerr << "Error in RHEHpt::C0(long double x1,long double x2, long double yt): Im(C0) > Re(C0)/1000 and the Im part is being ignored as negligible." << std::endl;
//    return std::real( T1 + T2 ) / ( 16. * M_PIl * M_PIl * sdelta);
//}

//long double A1x_0(long double x1,long double yt){
//        const long double one = 1;
//        const long double L1sq = -std::pow(atan(std::sqrt(4.*yt -1)/(2*yt -1)),2);
////		if (x1 < 1E-15){
////			std::cout << "A1(0,0), x = " << x1 << std::endl;
////         return (4. + L1sq * (4.*yt - 1.))/(32. * M_PIl * M_PIl );
////        }
//        
//        const long double L2 = std::log((std::sqrt(1.+4.*yt/x1) - one)/(std::sqrt(1.+4.*yt/x1) + one)); 
//        const long double C0term = L1sq - std::pow(L2,2);
//        const long double L3 = (4.*yt -1. - x1) /((32.*M_PIl*M_PIl)*std::pow(1. + x1,2));		
//        const long double B0term = 2*x1*(std::real(B0(-x1,yt)) - std::real(B1(1,yt)))/std::pow(1+x1,2);
//        return C0term * L3 + B0term + 1./(8.*M_PIl*M_PIl*(1.+x1));
//}

//long double A1(long double x1, long double x2, long double yt){
////    if(x1 < 1E-16 )	return A1x_0(x2,yt);
////    if(x2 < 1E-16 )	return A1x_0(x1,yt); 
//    
//    const long double Delta = Delta3(x1,x2);
//    const long double Delta2 = Delta * Delta;
//    const long double L1 = 4. * yt * ( 1. + x1 + x2)/ Delta - 1. - 4.*x1*x2/Delta + 12.*x1*x2*(1. + x1 + x2)/Delta2 ;

//    const long double L2 = (B0(-x2,yt) - B1(1.,yt)) * ( - 2.*x2/Delta + 12* x1 * x2 * ( 1. + x1 - x2)/Delta2);

//    const long double L3 = (B0(-x1,yt) - B1(1.,yt)) * ( - 2.*x1/Delta + 12* x1 * x2 * ( 1. - x1 + x2)/Delta2);

//    const long double L4 = 2. * (1. + x1 + x2) / ( Delta * (16. * M_PIl * M_PIl) );
////    std::cout << "A1 = " << C0(x1,x2,yt) << "*" << L1 << "\t"<< L2 << "\t"<< L3 << "\t"<< L4 << " yt=" << yt << std::endl;
////	std::cout << "A1(" << x1 << "," << x2 << ") = " << std::real(C0(x1,x2,yt) * L1 - L2 - L3 + L4) << std::endl;
//    return C0(x1,x2,yt) * L1 - L2 - L3 + L4; // C0 is complex
//}


//long double A2x_0(long double x1,long double yt){
//        const long double one = 1.;
//        const long double L1sq = -std::pow(atan(std::sqrt(4.*yt -1)/(2*yt -1)),2);
////		if (x1 < 1E-15){
////			std::cout << "A2(0,0), x = " << x1 << std::endl;
////         return (4. + L1sq * (4.*yt - 1.))/(64. * M_PIl * M_PIl );
////        }
//        
//        const long double L2 = std::log((std::sqrt(1.+4.*yt/x1) - one)/(std::sqrt(1.+4.*yt/x1) + one)); 
//        const long double C0term = (L1sq - std::pow(L2,2))/((32.*M_PIl*M_PIl)*(1.+x1));
//        const long double L3 = (2.*yt -0.5*(1.+x1));		
//        const long double B0term = (std::real(B0(-x1,yt)) - std::real(B1(1,yt)))*x1/(1.+x1);
////		std::cout << " fooo = " << C0term * L3 + B0term + 1./(4.*M_PIl*M_PIl) << std::endl;
//        return C0term * L3 + B0term + 1./(16.*M_PIl*M_PIl);
//}


//long double A2(long double x1, long double x2, long double yt){
////    if(x1 < 1E-16 )	return A2x_0(x2,yt);
////    if(x2 < 1E-16 )	return A2x_0(x1,yt); 
//    const long double Delta = Delta3(x1,x2);
//    const long double L1 = 2. * yt  - ( 1. + x1 + x2)/ 2. + 2.*x1*x2/Delta ;
//    const long double L2 = (B0(-x2,yt) - B1(1.,yt)) * ( x2 * ( 1. - x1 + x2)/Delta);
//    const long double L3 = (B0(-x1,yt) - B1(1.,yt)) * ( x1 * ( 1. + x1 - x2)/Delta);
//    const long double L4 = 1. / (16. * M_PIl * M_PIl);
////    if(x1 > 0.99 || x2 > 0.99)
////	    std::cout << "A2(" << x1 << "," << x2 << ") = " << C0(x1,x2,yt) << "*" << L1 << "\t"<< std::real(L2) << "\t"<< std::real(L3) << "\t"<< L4 << std::endl;
////	std::cout << "A2(" << x1 << "," << x2 << ") = " << std::real(C0(x1,x2,yt) * L1 + L2 + L3 + L4) << std::endl;
//    return C0(x1,x2,yt) * L1 + L2 + L3 + L4; // C0 is complex
//}

//long double A3(long double x1, long double x2, long double yt){
////	if(x2 < 1E-15 || x1 < 1E-15 ){	
////			const long double tau = std::sqrt( 4.*yt  - 1. ) * std::atan( 1./std::sqrt( 4.*yt - 1. ) ); // this is the leading order of A3(xp*x1,xp*x2) when xp -> 0. It is not to be used unless A3 is multiplied by some xp^2n * x1^n * x2^n factor, as in MA2
////			return (10. - 8.*tau + tau*tau )/(8.*M_PIl*M_PIl);
////	}
//	return (A1(x1,x2,yt)*((1.+(x1+x2))/2.) - A2(x1,x2,yt))/(x1*x2);
//}

//long double MA2(long double x1, long double x2, long double xp, long double yt, long double costh){
////	if(costh < -5.){	// no costh2 passed
////		costh = (1. - x1 - x2)/sqrt(4.*x1*x2);
////	}
//	const long double rA1 = std::real(RHEHpt::A1(xp*x1,xp*x2,yt));
//	const long double rA3 = std::real(RHEHpt::A3(xp*x1,xp*x2,yt));
//	const long double L1 = xp*xp * x1*x2 * rA3*rA3 ;
//	const long double L2 = sqrt(4.*x1*x2) * rA1*rA3 * xp*costh;
//	const long double L3 = rA1*costh * rA1*costh;
//	return L1 + L2 + L3;
//}
//};

dcomplex InverseMellinIntegrand(dcomplex N, double log1ox, std::function< dcomplex(dcomplex,const Rvector&)> f, const Rvector& par) {
  return exp( N*log1ox ) * f(N,par) / M_PI;
}


double InverseMellin(double z, std::function< dcomplex(dcomplex,const Rvector&)> f, const Rvector& par)
{ // Talbot algorithm
  if(z==1.) return 0;
  dcomplex I(0,1);
  double log1ox = -std::log(z);
  double res = 0.; 
  // prec ~ 10^{-0.6*M}
  int M = 50;
  double r = 1.+2.*M/5./log1ox;
  double t, cott, sigma;
  dcomplex N;
  res = std::real(0.5 * InverseMellinIntegrand(r, log1ox, f, par));
  for(int k = 1; k < M; k++) {
    t = k*M_PI/M;
    cott = 1./tan(t);
    sigma = t + cott*(t*cott - 1.);
    N = r*t*(cott+I);
    res += std::real(InverseMellinIntegrand(N, log1ox, f, par) * (1.+I*sigma));
  }
  res *= M_PI*r/M;
  return res;
}

