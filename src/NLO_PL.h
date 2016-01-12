#ifndef __NLOPL_H__
#define __NLOPL_H__

#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include "RHEHpt.h"
#include "LHAPDF/LHAPDF.h"


using namespace RHEHpt;

class NLOPL{
public:
  NLOPL(double CME, double MH=125.09, double MUF=125.09, double MUR=125.09,const std::string& PDFname = "NNPDF30_nnlo_as_0118"):_CME(CME),_mH(MH),_PDF(LHAPDF::mkPDF(PDFname,0)){
    x=std::pow(_mH/_CME,2);
    _muF=std::pow(MUF/_mH,2);
    _muR=std::pow(MUR/_mH,2);
    _as = _PDF -> alphasQ(_mH);
  }
  void SetCME(double CME){
    _CME = CME;
    x = std::pow(_mH/_CME,2);
  }
  ~NLOPL(){
    delete _PDF;
  }
  
  inline double get_alphas() const{ return _as; }
  inline LHAPDF::PDF* get_PDF() const { return _PDF; }
  inline double get_scale() const{ return _CME; }
  inline double get_mH() const{ return _mH; }
  inline double get_x() const{ return std::pow(_mH/_CME,2.); }
		
  long double LO_PL(double xp);
  long double NLO_PL(double xp);
  
  long double G2_sing(double xp);
  long double G2_notsing(double xp);
  long double G2_delta(double xp);
		
  long double operator ()(double xp, unsigned int order){
    if (order>1){ //Only available possibility are 0 (LO) or 1 (NLO)
      std::cout << "ERROR: ONLY AVAILABLE OPTIONS ARE 0 (LO) or 1 (NLO)" << std::endl;
      return 0;
    }
    if(order==0){	
      return LO_PL(xp);
    }
    if(order==1){
      return NLO_PL(xp);
    }
  }
    
  
private:
  long double _mH;
  long double _CME;
  long double _muF;
  long double _muR;
  long double _as; //VALUE OF as(mH)
  long double _Nc=3.0;
  long double _Nf=5.0;
  long double x; 
  LHAPDF::PDF* _PDF;
  
};

#endif /* __NLOPL_H__ */