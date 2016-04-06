#ifndef __NLOPL_H__
#define __NLOPL_H__

#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include "complex_def.h"
#include "LHAPDF/LHAPDF.h"



class NLOPL{
public:
  NLOPL(double CME, double AS, double MH=125.09, double MUF=125.09, double MUR=125.09, unsigned int channel=1):_CME(CME),_mH(MH),_as(AS),_channel(channel){
      x=std::pow(_mH/_CME,2);
    _muF=std::pow(MUF/_mH,2);
    _muR=std::pow(MUR/_mH,2);
  }
  NLOPL( const NLOPL& NLO){
    _mH=NLO._mH;
    _CME=NLO._CME;
    _muF=NLO._muF;
    _muR=NLO._muR;
    _as=NLO._as;
    x=NLO.x;
    _channel=NLO._channel;
  }
  void SetCME(double CME){
    _CME = CME;
    x = std::pow(_mH/_CME,2);
  }
  void SetAS(double AS){
    _as=AS;
  }
  
  void setchannel(unsigned int CHANNEL){
    if ((CHANNEL <1) || (CHANNEL>6)){
      std::cout << "ERROR: invalid Channel; options are (1)=GG, (2)=GQ, (3)=QQbar, (4)=QQ, (5)=QQbar', (6)=QQ'" << std::endl;
      CHANNEL=1;
    }
    _channel=CHANNEL;
  }
  
  inline double get_alphas() const{ return _as; }
  inline double get_scale() const{ return _CME; }
  inline double get_mH() const{ return _mH; }
  inline double get_x() const{ return std::pow(_mH/_CME,2.); }
  inline double get_muf()  const{ return _muF; }
  inline double get_mur() const{return _muR; }
  
  
		
  long double LO_PL(double xp);
  
  long double NLO_PL_sing_doublediff(double xp, double z);
  long double NLO_PL_notsing_doublediff(double xp, double z);
  long double NLO_PL_delta(double xp);

  long double NLO_PL_sing_doublediff_gq(double xp, double z);
  long double NLO_PL_notsing_doublediff_gq(double xp, double z);
  long double NLO_PL_delta_gq(double xp);
  
  long double NLO_PL_sing_doublediff_qqbar(double xp, double z);
  long double NLO_PL_notsing_doublediff_qqbar(double xp, double z);
  long double NLO_PL_delta_qqbar(double xp);
  
  long double NLO_PL_notsing_doublediff_qq(double xp, double z);
  long double NLO_PL_notsing_doublediff_qqprime(double xp, double z);
  long double NLO_PL_notsing_doublediff_qqbarprime(double xp, double z);
  
 
  
    
  
private:
  long double _mH;
  long double _CME;
  long double _muF;
  long double _muR;
  unsigned int _channel;
  long double _as; //VALUE OF as(mH)
  long double _Nc=3.0;
  long double _Nf=5.0;
  long double _Cf=4./3.;
  long double x; 
  
};

#endif /* __NLOPL_H__ */