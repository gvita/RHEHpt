#ifndef __LOBAUR_H__
#define __LOBAUR_H__

#include <iostream>
#include <cmath>
#include <vector>
#include <complex>


class LOBaur{
	public:
		LOBaur(double CME, double AS,double MT = 173.3, double MB=4.18, double MH = 125.09, unsigned int Choice=1, unsigned int channel=1):_CME(CME),_as(AS),_mT(MT),_mB(MB),_mH(MH),_choice(Choice),_channel(channel){
			x = std::pow(_mH/_CME,2);
			yt = std::pow(_mT/_mH,2);
			yb = std::pow(_mB/_mH,2);
		}
		
		
		void SetMandelstam(double xp){
			t = 0.5*(-1.+x-sqrt((1.-x)*(1.-x)-4.*x*xp));
			u = 0.5*(-1.+x+sqrt((1.-x)*(1.-x)-4.*x*xp));
			s1 = 1.-x;
			t1 = t-x;
			u1 = u-x;
		}
		
		void SetCME(double CME){
			_CME = CME;
			x = std::pow(_mH/_CME,2);
		}
		void SetAS(double AS){
		  _as=AS;
		}
		
		
		long double MatrixElementFO(long double xp);
		long double sigmapartLO(long double xp);
		
		long double operator ()(long double xp){
			return sigmapartLO(xp);
		}
		
		void setchoice(unsigned int CHOICE){
		if(CHOICE > 3 )
		    std::cout << "Error invalid choice; set a number from 0 to 3" << std::endl;  
		  else{
		    _choice=CHOICE;
		  }
		}
		void setchannel(unsigned int CHANNEL){
		if((CHANNEL > 6 )||(CHANNEL<1))
		    std::cout << "Error invalid choice; set a number from 1 to 6" << std::endl;  
		  else{
		    _channel=CHANNEL;
		  }
		}
		
		inline double get_alphas() const{ return _as; }
		inline double get_scale() const{ return _CME; }
		inline double get_mH() const{ return _mH; }
		inline double get_x() const{ return std::pow(_mH/_CME,2.); }
		inline double get_mt()  const{ return _mT; }
		inline double get_mb() const{return _mB; }
		inline unsigned int get_choice() const{return _choice; }
				
	private:
		std::complex<long double> B1(long double,long double);
		std::complex<long double> C(long double, long double );
		std::complex<long double> C1(long double, long double );
		std::complex<long double> D(long double,long double , long double );
		std::complex<long double> E(long double, long double, long double );
		std::complex<long double> Appp(long double );
		std::complex<long double> Ampm(long double );
		std::complex<long double> Ampp(long double );
		std::complex<long double> Appm(long double );
		std::complex<long double> Aqqbar(long double );
		std::complex<long double> Aqg(long double );
		std::complex<long double> Agqbar(long double );
		long double t,u,x,yt,yb,s1,t1,u1;
		
		long double _CME;	//	sqrt(s)
		long double _mT;		//	Mass of Top
		long double _mB;		// 	Mass of quark
		long double _mH;		//	Mass of Higgs
		unsigned int _choice;
		unsigned int _channel;
		long double _as;
		long double epsilon=1e-20; // Epsilon

		
		
		
};

#endif /* __LOBAUR_H__ */
