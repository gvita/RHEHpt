#ifndef __NLOPL_H__
#define __NLOPL_H__

#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include "LHAPDF/LHAPDF.h"



class NLOPL{
public:
  NLOPL(double CME, double AS, double MH=125.09, double MUF=125.09, double MUR=125.09):_CME(CME),_mH(MH),_as(AS){
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
  }
  void SetCME(double CME){
    _CME = CME;
    x = std::pow(_mH/_CME,2);
  }
  void SetAS(double AS){
    _as=AS;
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
		  
 
  
    
  
private:
  long double _mH;
  long double _CME;
  long double _muF;
  long double _muR;
  long double _as; //VALUE OF as(mH)
  long double _Nc=3.0;
  long double _Nf=5.0;
  long double x; 
  
};

#endif /* __NLOPL_H__ */