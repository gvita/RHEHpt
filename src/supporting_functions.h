#ifndef __supporting_functions__
#define __supporting_functions__

#include <complex>
#include <functional>
#include <vector>
//namespace RHEHpt{

//    long double B0(const long double x, const long double yt); // CORRECT ONLY IF ARGUMENT IS < 0
//    long double B1(const long double x, const long double yt); // CORRECT ONLY IF ARGUMENT IS > 0 && < 4 yt
//    long double Delta(long double x1, long double x2);
//    std::complex<long double> Li2(std::complex <long double> z);
//    long double C0(long double x1,long double x2, long double yt);

//    long double A1(long double x1, long double x2, long double yt);
//    long double A2(long double x1, long double x2, long double yt);
//    long double A3(long double x1, long double x2, long double yt); 
//	long double MA2(long double x1, long double x2, long double xp, long double yt, long double costh=-10.);    
//};
double InverseMellin(double z, std::function< std::complex<double>(std::complex<double>,const std::vector<double>&)> f,const std::vector<double>& par);
    

#endif
