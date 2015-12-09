#include <cmath>
#include <complex>

//Hard Part Functions
std::complex<long double> C0(long double x1,long double x2,long double y);
std::complex<long double> B0(long double x, long double y);
std::complex<long double> A1(long double x1, long double x2, long double y);
std::complex<long double> A1x_0(long double x1,long double y);
std::complex<long double> A2(long double x1, long double x2, long double y);
std::complex<long double> A3(long double x1, long double x2, long double y);
long double F_0(long double x1, long double x2, long double xp, long double y,long double r);
long double F_Infinitive(long double x1, long double x2, long double xp, long double y);



//First derivative wrt x1 of the Hard Part Functions
std::complex<long double> D_C0(long double x1, long double x2, long double y);
std::complex<long double> D_B0(long double x1, long double y);
std::complex<long double> D_A1(long double x1, long double x2, long double y);
std::complex<long double> D_A2(long double x1, long double x2, long double y);
std::complex<long double> D_A3(long double x1, long double x2, long double y);
long double D_F_0(long double x1, long double x2, long double xp, long double y,long double r);

//Function for Mtop infinitive case

long double F_MtopInf_0(long double x1, long double r, long double xp, long double y);
long double F_MtopInf_Inf(long double x1, long double x2, long double xp, long double y);
long double D_F_MtopInf(long double x1, long double r, long double xp, long double y);


