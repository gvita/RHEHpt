#ifndef __LOBAUR_H__
#define __LOBAUR_H__

#include <iostream>
#include <cmath>
#include <vector>
#include <complex>

class LOBaur{
	public:
		LOBaur(double CME,double MT = 173.3, double MH = 125.09):_CME(CME),_mq(MT),_mH(MH){}
		
		void SetMandelstam(double xp){
			x = pow(_mH/_CME,2);
			y = pow(_mq/_mH,2);
			t = 0.5*(-1.+x-sqrt((1.-x)*(1.-x)-4.*x*xp));
			u = 0.5*(-1.+x+sqrt((1.-x)*(1.-x)-4.*x*xp));
			s1 = 1.-x;
			t1 = t-x;
			u1 = u-x;
		}
		
		void SetCME(double CME){
			_CME = CME;
		}
		
		long double MatrixElementFO(double pt);
		
		long double operator ()(double pt){
			return MatrixElementFO(pt);
		}
				
	private:
		std::complex<long double> B1(double);
		std::complex<long double> C(double);
		std::complex<long double> C1(double);
		std::complex<long double> D(double, double );
		std::complex<long double> E(double, double);
		std::complex<long double> Appp();
		std::complex<long double> Ampm();
		std::complex<long double> Ampp();
		std::complex<long double> Appm();

		long double t,u,x,y,s1,t1,u1;
		long double _CME;	//	sqrt(s)
		long double _mq;		//	Mass of quark
		long double _mH;		//	Mass of Higgs
};

#endif /* __LOBAUR_H__ */
