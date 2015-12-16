#include <iostream> 			// just for debugging purpose 
#include "complex_def.h"
#include "Mfinite_HardPart_2.h"

using namespace std;


//NOTe IN C0,A1,A2 x1 and x2 are xi and xi_bar
std::complex<long double> C0(long double x1, long double x2,long double y)
{
  //Useful quantities
  long double epsilon=1e-20;
  long double D3=(1.+x1+x2)*(1.+x1+x2)-4.*x1*x2;
  long double delta1= (-x1+x2-1.)/sqrt(D3);
  long double delta2=(x1-x2-1.)/sqrt(D3);
  long double delta3=(x1+x2+1.)/sqrt(D3);
  long double delta1p=(1.+delta1)/2.;
  long double delta1m=(1.-delta1)/2.;
  long double delta2p=(1.+delta2)/2.;
  long double delta2m=(1.-delta2)/2.;
  long double delta3p=(1.+delta3)/2.;
  long double delta3m=(1.-delta3)/2.;
  std::complex<long double> rad_x1 = static_cast<std::complex<long double> > (1.+4.*y/x1); 
  std::complex<long double> rad_x2 = static_cast<std::complex<long double> > (1.+4.*y/x2); 
  std::complex<long double> rad_y = static_cast<std::complex<long double> > (1.-4.*(y-std::complex<long double>(0.,1.)*epsilon)); 
  std::complex<long double> xmi=-x2/(2.*y)*(1.-sqrt(rad_x2));
  std::complex<long double> xpi=-x2/(2.*y)*(1.+sqrt(rad_x2));
  std::complex<long double> ymi=-x1/(2.*y)*(1.-sqrt(rad_x1));
  std::complex<long double> ypi=-x1/(2.*y)*(1.+sqrt(rad_x1));
  std::complex<long double> zmi=1./(2.*y)*(1.- sqrt(rad_y));
  std::complex<long double> zpi=1./(2.*y)*(1.+ sqrt(rad_y));
  
  std::complex<long double> C0_ris;
  //C0 expression
  C0_ris=1./(16.*M_PI*M_PI*sqrt(D3))*(log_c(1.-ymi)*(log_c((1.-ymi*delta1p))-log_c((1.-ymi*delta1m)))
  +log_c(1.-xmi)*(log_c((1.-xmi*delta2p))-log_c((1.-xmi*delta2m)))+dilog_c(ypi*delta1p)+dilog_c(ymi*delta1p)
  -dilog_c(ypi*delta1m)-dilog_c(ymi*delta1m)+dilog_c(xpi*delta2p)+dilog_c(xmi*delta2p)
  -dilog_c(xpi*delta2m)-dilog_c(xmi*delta2m)+log_c(1.-zmi)*(log_c((1.-zmi*delta3p))-log_c((1.-zmi*delta3m)))
  +dilog_c(zpi*delta3p)+dilog_c(zmi*delta3p)-dilog_c(zpi*delta3m)-dilog_c(zmi*delta3m));
  return C0_ris;
}

std::complex<long double> B0(long double x,long double y)
{
    long double epsilon=1e-20;
    std::complex<long double> rad_x = static_cast<std::complex<long double> > (1.-4.*(y-std::complex<long double>(0.,1.)*epsilon)/x);
   // cout <<log_c((1.+1./sqrt(rad_x))/(1.-1./sqrt(rad_x)),-1.)<< endl;
    std::complex<long double> ris= (-1./(16.*M_PI*M_PI)*sqrt(rad_x)*(log_c((1.+1./sqrt(rad_x))/(1.-1./sqrt(rad_x)))));
    return ris;
}
std::complex<long double> A1(long double x1, long double x2, long double y)
{
  if(x1 < 1E-15) return A1x_0(x2,y);
  if(x2 < 1E-15) return A1x_0(x1,y);

  //Useful quantities
  long double D3=(1.+x1+x2)*(1.+x1+x2)-4.*x1*x2;
  //A1 Expression
  std::complex<long double> A1;
  A1=C0(x1,x2,y)*(4.*y/D3*(1.+x1+x2)-1.-4.*x1*x2/D3+12.*x1*x2/D3/D3*(1.+x1+x2))
  -(B0(-x2,y)-B0(1.,y))*(-2.*x2/D3+12.*x1*x2/D3/D3*(1.+x1-x2))
  -(B0(-x1,y)-B0(1.,y))*(-2.*x1/D3+12.*x1*x2/D3/D3*(1.-x1+x2))
  +2./D3*1./(16.*M_PI*M_PI)*(1.+x1+x2);
  return A1;
}
std::complex<long double> A2(long double x1, long double x2, long double y)
{
  //Useful quantities
  long double D3=(1.+x1+x2)*(1.+x1+x2)-4.*x1*x2;
  //A1 Expression
  std::complex<long double> A2;
  A2=C0(x1,x2,y)*(2.*y-1./2.*(1.+x1+x2)+2.*x1*x2/D3)
  +(B0(-x2,y)-B0(1.,y))*(x2/D3*(1.-x1+x2))
  +(B0(-x1,y)-B0(1.,y))*(x1/D3*(1.+x1-x2))
  +1./(16.*M_PI*M_PI);
  return A2;
}
std::complex<long double> A3(long double x1, long double x2, long double y)
{
  //Definition of A3=1/(x1 x2)((1+x1+x2)/2 A1-A2)
  //Useful quantities
  long double D3=(1.+x1+x2)*(1.+x1+x2)-4.*x1*x2;
  //A3 Expression
  std::complex<long double> A3;
  A3=C0(x1,x2,y)*(8.*y-2.*(x1+x2-1.)+24.*x1*x2/D3)+(B0(-x2,y)-B0(1.,y))*(2.-6./D3*(1.+x1-x2)*(1.+x1+x2))
  +(B0(-x1,y)-B0(1.,y))*(2.-6./D3*(1.-x1+x2)*(1.+x1+x2))+1./(4.*M_PI*M_PI);
  return A3/D3;
}

long double F_0(long double x1, long double x2, long double xp, long double y,long double r)
{
  //definition of x2
  const long double sqrtx1 = std::sqrt(x1);
  //Expression of F
  const std::complex<long double> A1_v=A1(xp*x1,xp*x2,y);
  const std::complex<long double> A3_v=A3(xp*x1,xp*x2,y);
  long double F_ris=(pow(1.+sqrtx1-2.*r,2)/(x2)*real(A1_v*std::conj(A1_v))+xp*xp*x1*x2*real(A3_v*std::conj(A3_v))
  -xp*sqrtx1*(1.+sqrtx1-2.*r)*real(std::conj(A1_v)*A3_v+std::conj(A3_v)*A1_v));
  //std::cout << F_ris << endl;
  return y*y*F_ris;
}

long double F_Infinitive(long double x1, long double x2, long double xp, long double y)
{ // i.e. F_0 without the definition of A3
  //definition of x2

//  long double x2=(1.+sqrt(x1))*(1.+sqrt(x1))-(1.+2.*sqrt(x1))*r;

  //Expression of F
  const std::complex<long double> A1_v=A1(xp*x1,xp*x2,y);
  const std::complex<long double> A2_v=A2(xp*x1,xp*x2,y);
  long double F_ris=1./(xp*xp*x1*x2)*pow(abs((1.+xp)/2.*A1_v-A2_v),2);
  return y*y*F_ris;
}

long double F_0_inter(long double x1, long double x2, long double xp, long double y, long double y2, long double r){
  const long double sqrtx1 = std::sqrt(x1);
  const std::complex<long double> A1_v_yt=A1(xp*x1,xp*x2,y);
  const std::complex<long double> A3_v_yt=A3(xp*x1,xp*x2,y);
  const std::complex<long double> A1_v_yb=A1(xp*x1,xp*x2,y2);
  const std::complex<long double> A3_v_yb=A3(xp*x1,xp*x2,y2);
  long double F_ris=(pow(1.+sqrtx1-2.*r,2)/x2*real(A1_v_yt*std::conj(A1_v_yb)+std::conj(A1_v_yt)*A1_v_yb)
  +xp*xp*x1*x2*real(A3_v_yt*std::conj(A3_v_yb)+std::conj(A3_v_yt)*A3_v_yb)
  -xp*sqrtx1*(1.+sqrtx1-2.*r)*real(std::conj(A1_v_yt)*A3_v_yb
  +std::conj(A3_v_yb)*A1_v_yt+std::conj(A1_v_yb)*A3_v_yt+std::conj(A3_v_yt)*A1_v_yb));
  return y*y2*F_ris;  
}

long double F_Infinitive_inter(long double x1, long double x2, long double xp, long double y, long double y2)
{
  const std::complex<long double> A1_v_yt=A1(xp*x1,xp*x2,y);
  const std::complex<long double> A2_v_yt=A2(xp*x1,xp*x2,y);
  const std::complex<long double> A1_v_yb=A1(xp*x1,xp*x2,y2);
  const std::complex<long double> A2_v_yb=A2(xp*x1,xp*x2,y2);
  long double F_ris=1./(xp*xp*x1*x2)*((1.+xp)*(1.+xp)/4.*real(A1_v_yt*std::conj(A1_v_yb)+A1_v_yb*std::conj(A1_v_yt))
  +real(A2_v_yt*std::conj(A2_v_yb)+A2_v_yb*std::conj(A2_v_yt))
  -(1.+xp)/2.*real(std::conj(A1_v_yt)*A2_v_yb+A1_v_yt*std::conj(A2_v_yb)+std::conj(A2_v_yt)*A1_v_yb+A2_v_yt*std::conj(A1_v_yb)));
  return y*y2*F_ris;
}

long double F_0_tot(long double x1, long double x2, long double xp, long double y, long double y2,long double r){
  const long double sqrtx1 = std::sqrt(x1);
  const std::complex<long double> A1_v_yt=A1(xp*x1,xp*x2,y);
  const std::complex<long double> A3_v_yt=A3(xp*x1,xp*x2,y);
  const std::complex<long double> A1_v_yb=A1(xp*x1,xp*x2,y2);
  const std::complex<long double> A3_v_yb=A3(xp*x1,xp*x2,y2);
  
  std::complex<long double> Amp=y*(-(1.+sqrtx1-2.*r)/sqrt(x2)*A1_v_yt+xp*sqrtx1*sqrt(x2)*A3_v_yt)
  +y2*(-(1.+sqrtx1-2.*r)/sqrt(x2)*A1_v_yb+xp*sqrtx1*sqrt(x2)*A3_v_yb); 
  
  long double F_ris=real(Amp*std::conj(Amp));
  return F_ris;
}

long double F_Infinitive_tot(long double x1, long double x2, long double xp, long double y, long double y2){
  const std::complex<long double> A1_v_yt=A1(xp*x1,xp*x2,y);
  const std::complex<long double> A2_v_yt=A2(xp*x1,xp*x2,y);
  const std::complex<long double> A1_v_yb=A1(xp*x1,xp*x2,y2);
  const std::complex<long double> A2_v_yb=A2(xp*x1,xp*x2,y2);
  
  std::complex<long double> Amp=y*((1.+xp)/2.*A1_v_yt-A2_v_yt)+y2*((1.+xp)/2.*A1_v_yb-A2_v_yb);
  
  long double F_ris=1./(xp*xp*x1*x2)*real(Amp*std::conj(Amp));
  return F_ris;
}





std::complex<long double> D_C0(long double x1, long double x2, long double y){
  long double epsilon=1e-20;
   const long double tiny=1e-6;
  //Useful quantities
  long double D3=(1.+x1+x2)*(1.+x1+x2)-4.*x1*x2;
  long double delta1= (-x1+x2-1.)/sqrt(D3);
  long double delta2=(x1-x2-1.)/sqrt(D3);
  long double delta3=(x1+x2+1.)/sqrt(D3);
  long double delta1p=(1.+delta1)/2.;
  long double delta1m=(1.-delta1)/2.;
  long double delta2p=(1.+delta2)/2.;
  long double delta2m=(1.-delta2)/2.;
  long double delta3p=(1.+delta3)/2.;
  long double delta3m=(1.-delta3)/2.;
  
  std::complex<long double> rad_x1 = static_cast<std::complex<long double> > (1.+4.*y/x1); 
  std::complex<long double> rad_x2 = static_cast<std::complex<long double> > (1.+4.*y/x2); 
  std::complex<long double> rad_y = static_cast<std::complex<long double> > (1.-4.*(y-std::complex<long double>(0.,1.)*epsilon)); 
  std::complex<long double> xmi=-x2/(2.*y)*(1.-sqrt(rad_x2));
  std::complex<long double> xpi=-x2/(2.*y)*(1.+sqrt(rad_x2));
  std::complex<long double> ymi=-x1/(2.*y)*(1.-sqrt(rad_x1));
  std::complex<long double> ypi=-x1/(2.*y)*(1.+sqrt(rad_x1));
  std::complex<long double> zmi=1./(2.*y)*(1.- sqrt(rad_y));
  std::complex<long double> zpi=1./(2.*y)*(1.+ sqrt(rad_y));
  
  //Useful Derivatives
  long double D_D3=2.*(1.+x1-x2);
  long double D_delta1=-4.*x2/pow(sqrt(D3),3);
  long double D_delta2=2.*(1.+x1+x2)/pow(sqrt(D3),3);
  long double D_delta3=2.*x2*(1.-x1+x2)/pow(sqrt(D3),3);
  
  std::complex<long double> D_ymi=(1.-ymi)/(x1*sqrt(rad_x1));
  std::complex<long double> D_ypi=-(1.-ypi)/(x1*sqrt(rad_x1));
  
  //Derivate of C0
  std::complex<long double> D_C0_ris=-1./D3*(1.+x1-x2)*C0(x1,x2,y)+1./(16.*M_PI*M_PI)*1./sqrt(D3)*
  (-1./(1.-ymi)*D_ymi*(log_c((1.-ymi*delta1p))-log_c((1.-ymi*delta1m)))+
  log_c(1.-ymi)/((1.-ymi*delta1p)*(1.-ymi*delta1m))*(-(1.-ymi*delta1m)*(D_ymi*delta1p+ymi/2.*D_delta1)
  +(1.-ymi*delta1p)*(D_ymi*delta1m-ymi/2.*D_delta1))-log_c(1.-xmi)/((1.-xmi*delta2p)*(1.-xmi*delta2m))*
  ((1.-xmi*delta2p)+(1.-xmi*delta2m))*xmi/2.*D_delta2-log_c(1.-zmi)/((1.-zmi*delta3p)*(1.-zmi*delta3m))*
  ((1.-zmi*delta3p)+(1.-zmi*delta3m))*zmi/2.*D_delta3-log_c(1.-xmi*delta2m)/(xmi*delta2m)*(xmi/2.*D_delta2)
  -log_c(1.-xpi*delta2m)/(xpi*delta2m)*(xpi/2.*D_delta2)
  -log_c(1.-zmi*delta3p)/(zmi*delta3p)*(zmi/2.*D_delta3)
  -log_c(1.-zpi*delta3p)/(zpi*delta3p)*(zpi/2.*D_delta3)
  );
  //0/0 part in D_C0 when x1->0
  std::complex<long double> D_C0_diver=0.0;
  if ((x1 > tiny)&&(x2>tiny)) {
    D_C0_diver= +1./(16.*M_PI*M_PI)*1./sqrt(D3)*(-log_c(1.-ypi*delta1p)/(ypi*delta1p)*(D_ypi*delta1p+ypi/2.*D_delta1)
  -log_c(1.-ymi*delta1p)/(ymi*delta1p)*(D_ymi*delta1p+ymi/2.*D_delta1)
  +log_c(1.-ypi*delta1m)/(ypi*delta1m)*(D_ypi*delta1m-ypi/2.*D_delta1)
  +log_c(1.-ymi*delta1m)/(ymi*delta1m)*(D_ymi*delta1m-ymi/2.*D_delta1)
  -log_c(1.-xmi*delta2p)/(xmi*delta2p)*(xmi/2.*D_delta2)
  -log_c(1.-xpi*delta2p)/(xpi*delta2p)*(xpi/2.*D_delta2)
  -log_c(1.-zmi*delta3m)/(zmi*delta3m)*(zmi/2.*D_delta3)
  -log_c(1.-zpi*delta3m)/(zpi*delta3m)*(zpi/2.*D_delta3));
  }
  if ((x1< tiny)&&(x2 >tiny)){
    D_C0_diver= +1./(16.*M_PI*M_PI)*1./sqrt(D3)*((1.-x2)/(2.*(1.+x2)*y));
  } 
  if ((x1> tiny)&&(x2 < tiny)){
    D_C0_diver=+1./(16.*M_PI*M_PI)*1./sqrt(D3)*(log_c((1.+2.*y/x1+sqrt(rad_x1))/(1.+2.*y/x1-sqrt(rad_x1)))/(2.*x1*sqrt(rad_x1)));
  }
  if ((x1 < tiny)&&(x2 <tiny)){
    D_C0_diver=+1./(16.*M_PI*M_PI)*1./sqrt(D3)*(1./(2.*y));
  }
  std::complex<long double> ris=D_C0_ris+D_C0_diver;
  /*if(real(ris) != real(ris)){
    cout << "DC0nan " << x1 << " " << x2 << " " << y << endl;
  }*/
  return D_C0_ris+D_C0_diver;
}
std::complex<long double> D_B0(long double x, long double y){ 
  std::complex<long double> rad_x = static_cast<std::complex<long double> > (1.+4.*y/x);
  return (-1./(16.*M_PI*M_PI)*(1./x-2.*y/x/x/sqrt(rad_x)*(log_c((1.+sqrt(rad_x)))-log_c((sqrt(rad_x)-1.)))));
}
std::complex<long double> D_A1(long double x1, long double x2, long double y){
  //Useful quantities
  long double D3=(1.+x1+x2)*(1.+x1+x2)-4.*x1*x2;
  //Derivate of A1
  std::complex<long double> D_A1_ris=D_C0(x1,x2,y)*(4.*y/D3*(1.+x1+x2)-1.-4.*x1*x2/D3+12.*x1*x2*(1.+x1+x2)/D3/D3)
  +C0(x1,x2,y)*(4.*y/D3/D3*(3.*x2*x2+2.*x2*(1.-x1)-(1.+x1)*(1.+x1))
		-4.*x2/D3/D3*((1.+x2)*(1.+x2)-x1*x1)+12.*x2/D3/D3/D3*(pow(1.+x2,3)+4.*x1*x2*(1.+x2)-3*x1*x1*(1.+x2)-2.*x1*x1*x1))
  -(B0(-x2,y)-B0(1.,y))*(-4.*x2/pow(D3,3)*(5.*x1*x1*x1-6.*x1*x1*(-1.+x2)+4.*(-1.+x2)*(1.+x2)*(1.+x2)-x1*(3.+22.*x2+3.*x2*x2)))
  -(B0(-x1,y)-B0(1.,y))*(2./pow(D3,3)*(x1*x1*x1*x1-18.*x1*x1*x2*(1.+x2)+pow(1.+x2,3)*(-1.+5.*x2)+2.*x1*x1*x1*(1.+5.*x2)
  +2.*x1*(-1.-13.*x2-11.*x2*x2+pow(x2,3))))
  -D_B0(x1,y)*(-2.*x1/D3+12.*x1*x2/D3/D3*(1.-x1+x2))+2./D3/D3*(3.*x2*x2+2.*x2*(1.-x1)-(1.+x1)*(1.+x1))/(16.*M_PI*M_PI);
  return D_A1_ris;
}
std::complex<long double> D_A2(long double x1, long double x2, long double y){
  //Useful quantities
  long double D3=(1.+x1+x2)*(1.+x1+x2)-4.*x1*x2;
  //Derivate of A2
  std::complex<long double> D_A2_ris=D_C0(x1,x2,y)*(2.*y-1./2.*(1.+x1+x2)+2.*x1*x2/D3)-C0(x1,x2,y)*
  (1./2.+4.*x1*x2/D3/D3*(1.+x1-x2)-2.*x2/D3)
  +(B0(-x2,y)-B0(1.,y))*(-2.*x2/D3/D3*(1.+x1-x2)*(1.-x1+x2)-x2/D3)
  +(B0(-x1,y)-B0(1.,y))*(-2.*x1/D3/D3*(1.+x1-x2)*(1.+x1-x2)+1./D3*(1.+x1-x2)+x1/D3)+D_B0(x1,y)*(x1/D3*(1.+x1-x2));
  return D_A2_ris;
}
std::complex<long double> D_A3(long double x1, long double x2, long double y){
//Useful quantities
  long double D3=(1.+x1+x2)*(1.+x1+x2)-4.*x1*x2;
  long double D_D3=2.*(1.+x1-x2);
  //Derivate of A3
  std::complex<long double> D_A3_ris;
  D_A3_ris=-1./D3*D_D3*A3(x1,x2,y)+1./D3*(D_C0(x1,x2,y)*(8.*y-2.*(x1+x2-1.)+24.*x1*x2/D3)
  +C0(x1,x2,y)*(-2.+24.*x2/D3-24.*D_D3*x1*x2/D3/D3)
  +(B0(-x2,y)-B0(1.,y))*(6./D3/D3*D_D3*(1.+x1-x2)*(1.+x1+x2)-6./D3*(1.+x1+x2)-6./D3*(1.+x1-x2))
  +(B0(-x1,y)-B0(1.,y))*(6./D3/D3*D_D3*(1.-x1+x2)*(1.+x1+x2)+6./D3*(1.+x1+x2)-6./D3*(1.-x1+x2))
  +D_B0(x1,y)*(2.-6./D3*(1.-x1+x2)*(1.+x1+x2)));
  return D_A3_ris;
}

long double D_F_0(long double x1, long double x2, long double xp, long double y,long double r){
  const long double sqrtx1 = std::sqrt(x1);
 const long double tiny=1e-8;
 
  //Compute one time function and derivatives
  const std::complex<long double> A1_v=A1(xp*x1,xp*x2,y);
  const std::complex<long double> A3_v=A3(xp*x1,xp*x2,y);
  //times xp to pass from derivates wrt xi to derivates wrt x1
   const std::complex<long double> D_A1_x1=xp*D_A1(xp*x1,xp*x2,y);
   const std::complex<long double> D_A1_x2=xp*D_A1(xp*x2,xp*x1,y);
   const std::complex<long double> D_A3_x1=xp*D_A3(xp*x1,xp*x2,y);
   const std::complex<long double> D_A3_x2=xp*D_A3(xp*x2,xp*x1,y);
   
   //Total derivates of A1 and A3
   const std::complex<long double> D_A1_v = D_A1_x1+D_A1_x2*(1.+sqrtx1-2.*r)/sqrtx1;
   const std::complex<long double> D_A3_v= D_A3_x1+D_A3_x2*(1.+sqrtx1-2.*r)/sqrtx1;
   
  //Derivate of F
  long double D_F_ris=0.0;
  
  D_F_ris=(xp > tiny)?((4.*(1.-r)*r*(1+sqrtx1-2.*r))/(sqrtx1*x2*x2)*real(A1_v*std::conj(A1_v))
  +pow(1.+sqrtx1-2.*r,2)/x2*real(A1_v*std::conj(D_A1_v)+std::conj(A1_v)*D_A1_v))
  +(xp*xp*x2*real(A3_v*std::conj(A3_v))+xp*xp*sqrtx1*(1.+sqrtx1-2.*r)*real(A3_v*std::conj(A3_v))
  +xp*xp*x1*x2*real(A3_v*std::conj(D_A3_v)+std::conj(A3_v)*D_A3_v))
  +(0.5*(-xp*(1.+2.*sqrtx1-2.*r)/sqrtx1*real(std::conj(A1_v)*A3_v+std::conj(A3_v)*A1_v))-xp*sqrtx1*(1.+sqrtx1-2.*r)
  *real(std::conj(D_A1_v)*A3_v+std::conj(A1_v)*D_A3_v+std::conj(D_A3_v)*A1_v+std::conj(A3_v)*D_A1_v))
  :((4.*(1.-r)*r*(1+sqrtx1-2.*r))/(sqrtx1*x2*x2)*real(A1_v*std::conj(A1_v)));
  return y*y*D_F_ris;
}

long double D_F_0_inter(long double x1, long double x2, long double xp, long double y, long double y2, long double r){
  const long double sqrtx1 = std::sqrt(x1);
  const long double tiny=1e-8;
  
  const std::complex<long double> A1_v_yt=A1(xp*x1,xp*x2,y);
  const std::complex<long double> A3_v_yt=A3(xp*x1,xp*x2,y);
  const std::complex<long double> A1_v_yb=A1(xp*x1,xp*x2,y2);
  const std::complex<long double> A3_v_yb=A3(xp*x1,xp*x2,y2);
  
  const std::complex<long double> D_A1_x1_yt=xp*D_A1(xp*x1,xp*x2,y);
  const std::complex<long double> D_A1_x2_yt=xp*D_A1(xp*x2,xp*x1,y);
  const std::complex<long double> D_A3_x1_yt=xp*D_A3(xp*x1,xp*x2,y);
  const std::complex<long double> D_A3_x2_yt=xp*D_A3(xp*x2,xp*x1,y);
  
  const std::complex<long double> D_A1_x1_yb=xp*D_A1(xp*x1,xp*x2,y2);
  const std::complex<long double> D_A1_x2_yb=xp*D_A1(xp*x2,xp*x1,y2);
  const std::complex<long double> D_A3_x1_yb=xp*D_A3(xp*x1,xp*x2,y2);
  const std::complex<long double> D_A3_x2_yb=xp*D_A3(xp*x2,xp*x1,y2);
  
  const std::complex<long double> D_A1_v_yt = D_A1_x1_yt+D_A1_x2_yt*(1.+sqrtx1-2.*r)/sqrtx1;
  const std::complex<long double> D_A3_v_yt= D_A3_x1_yt+D_A3_x2_yt*(1.+sqrtx1-2.*r)/sqrtx1;
  const std::complex<long double> D_A1_v_yb = D_A1_x1_yb+D_A1_x2_yb*(1.+sqrtx1-2.*r)/sqrtx1;
  const std::complex<long double> D_A3_v_yb= D_A3_x1_yb+D_A3_x2_yb*(1.+sqrtx1-2.*r)/sqrtx1;
  
  long double D_F_ris=0.0;
  
  D_F_ris=(xp > tiny) ? ((4.*(1.-r)*r*(1+sqrtx1-2.*r))/(sqrtx1*x2*x2)*real(A1_v_yt*std::conj(A1_v_yb)+A1_v_yb*std::conj(A1_v_yt))
  +pow(1.+sqrtx1-2.*r,2)/x2*real(D_A1_v_yt*std::conj(A1_v_yb)+A1_v_yt*std::conj(D_A1_v_yb)
  +std::conj(D_A1_v_yt)*A1_v_yb+D_A1_v_yb*std::conj(A1_v_yt))
  +xp*xp*x2*real(A3_v_yt*std::conj(A3_v_yb)+std::conj(A3_v_yt)*A3_v_yb)
  +xp*xp*sqrtx1*(1.+sqrtx1-2.*r)*real(A3_v_yt*std::conj(A3_v_yb)+std::conj(A3_v_yt)*A3_v_yb)
  +xp*xp*x1*x2*real(D_A3_v_yt*std::conj(A3_v_yb)+A3_v_yt*std::conj(D_A3_v_yb)+A3_v_yt*std::conj(D_A3_v_yb)+D_A3_v_yt*std::conj(A3_v_yb))
  +0.5*(-xp*(1.+2.*sqrtx1-2.*r)/sqrtx1)*real(std::conj(A1_v_yt)*A3_v_yb
  +std::conj(A3_v_yb)*A1_v_yt+std::conj(A1_v_yb)*A3_v_yt+std::conj(A3_v_yt)*A1_v_yb)
  -xp*sqrtx1*(1.+sqrtx1-2.*r)*real(A1_v_yt*std::conj(D_A3_v_yb)
  +D_A1_v_yt*std::conj(A3_v_yb)+std::conj(A1_v_yt)*D_A3_v_yb+std::conj(D_A1_v_yt)*A3_v_yb
  +std::conj(A3_v_yt)*D_A1_v_yb+std::conj(D_A3_v_yt)*A1_v_yb+A3_v_yt*std::conj(D_A1_v_yb)
  +D_A3_v_yt*std::conj(A1_v_yb)))
  :((4.*(1.-r)*r*(1+sqrtx1-2.*r))/(sqrtx1*x2*x2)*real(A1_v_yt*std::conj(A1_v_yb)+A1_v_yb*std::conj(A1_v_yt)));

  return y*y2*D_F_ris; 
}

long double D_F_0_tot(long double x1, long double x2, long double xp, long double y, long double y2, long double r){
  const long double sqrtx1 = std::sqrt(x1);
  const long double tiny = 1e-8;
  
  const std::complex<long double> A1_v_yt=A1(xp*x1,xp*x2,y);
  const std::complex<long double> A3_v_yt=A3(xp*x1,xp*x2,y);
  const std::complex<long double> A1_v_yb=A1(xp*x1,xp*x2,y2);
  const std::complex<long double> A3_v_yb=A3(xp*x1,xp*x2,y2);
  
  const std::complex<long double> D_A1_x1_yt=xp*D_A1(xp*x1,xp*x2,y);
  const std::complex<long double> D_A1_x2_yt=xp*D_A1(xp*x2,xp*x1,y);
  const std::complex<long double> D_A3_x1_yt=xp*D_A3(xp*x1,xp*x2,y);
  const std::complex<long double> D_A3_x2_yt=xp*D_A3(xp*x2,xp*x1,y);
  
  const std::complex<long double> D_A1_x1_yb=xp*D_A1(xp*x1,xp*x2,y2);
  const std::complex<long double> D_A1_x2_yb=xp*D_A1(xp*x2,xp*x1,y2);
  const std::complex<long double> D_A3_x1_yb=xp*D_A3(xp*x1,xp*x2,y2);
  const std::complex<long double> D_A3_x2_yb=xp*D_A3(xp*x2,xp*x1,y2);
  
  const std::complex<long double> D_A1_v_yt = D_A1_x1_yt+D_A1_x2_yt*(1.+sqrtx1-2.*r)/sqrtx1;
  const std::complex<long double> D_A3_v_yt= D_A3_x1_yt+D_A3_x2_yt*(1.+sqrtx1-2.*r)/sqrtx1;
  const std::complex<long double> D_A1_v_yb = D_A1_x1_yb+D_A1_x2_yb*(1.+sqrtx1-2.*r)/sqrtx1;
  const std::complex<long double> D_A3_v_yb= D_A3_x1_yb+D_A3_x2_yb*(1.+sqrtx1-2.*r)/sqrtx1;
  
  std::complex<long double> Amp= y*(-(1.+sqrtx1-2.*r)/sqrt(x2)*A1_v_yt+xp*sqrtx1*sqrt(x2)*A3_v_yt)
  +y2*(-(1.+sqrtx1-2.*r)/sqrt(x2)*A1_v_yb+xp*sqrtx1*sqrt(x2)*A3_v_yb);
  std::complex<long double> DAmp=y*(-(1.+sqrtx1-2.*r)/sqrt(x2)*D_A1_v_yt-2.*(1.-r)*r/sqrtx1/pow(sqrt(x2),3.)*A1_v_yt
  +xp*sqrtx1*sqrt(x2)*D_A3_v_yt+xp*(1.+(3.-6.*r)*sqrtx1+2.*x1)/(2.*sqrtx1*sqrt(x2))*A3_v_yt)
  +y2*(-(1.+sqrtx1-2.*r)/sqrt(x2)*D_A1_v_yb-2.*(1.-r)*r/sqrtx1/pow(sqrt(x2),3.)*A1_v_yb
  +xp*sqrtx1*sqrt(x2)*D_A3_v_yb+xp*(1.+(3.-6.*r)*sqrtx1+2.*x1)/(2.*sqrtx1*sqrt(x2))*A3_v_yb);
  std::complex<long double> Amp0=y*(-(1.+sqrtx1-2.*r)/sqrt(x2)*A1_v_yt)+y2*(-(1.+sqrtx1-2.*r)/sqrt(x2)*A1_v_yb);
  std::complex<long double> DAmp0=y*(-2.*(1.-r)*r/sqrtx1/pow(sqrt(x2),3.)*A1_v_yt)+y2*(-2.*(1.-r)*r/sqrtx1/pow(sqrt(x2),3.)*A1_v_yb);
  
  long double D_F_ris=0.0;
  D_F_ris=(xp > tiny) ? real(Amp*std::conj(DAmp)+DAmp*std::conj(Amp)):real(Amp0*std::conj(DAmp0)+DAmp0*std::conj(Amp0));
  return D_F_ris;
  
}


std::complex<long double> A1x_0(long double x1,long double yt){
	std::complex<long double> rad_y=static_cast<std::complex<long double> > (1.-4.*yt);
	std::complex<long double> rad_x1=static_cast<std::complex<long double> > (1.+4.*yt/x1);
        const std::complex<long double> L1sq = std::pow(log_c((sqrt(rad_y)-1.)/(sqrt(rad_y)+1.)),2);
		if (x1 < 1E-15){
			std::cout << "A1(0,0), x = " << x1 << std::endl;
         return (4. + L1sq * (4.*yt - 1.))/(32. * M_PIl * M_PIl );
        }
        
        const std::complex<long double> L2 = log_c((sqrt(rad_x1) - 1.))-log_c((sqrt(rad_x1) + 1.)); 
        const std::complex<long double> C0term = L1sq - std::pow(L2,2);
        const long double L3 = (4.*yt -1. - x1) /((32.*M_PIl*M_PIl)*std::pow(1. + x1,2));		
        const std::complex<long double> B0term = 2*x1*(B0(-x1,yt) - B0(1.,yt))/std::pow(1+x1,2);
        return yt*(C0term * L3 + B0term + 1./(8.*M_PIl*M_PIl*(1.+x1)));
	//return B0term;

  
}




//Functions Mtop Infinitive case (Cross Checking)
long double F_MtopInf_0(long double x1, long double r, long double xp, long double y)
{
  //definition of x2
  long double x2=(1.+sqrt(x1))*(1.+sqrt(x1))-(4.*sqrt(x1))*r;
  //Overall Constant terms
  const long double K=2./M_PI;
  //Expression of F
  
  long double F_ris=K*pow(1.+sqrt(x1)-2.*r,2)/x2;
  return F_ris;
}
long double F_MtopInf_Inf(long double x1, long double x2, long double xp, long double y)
{
  //Overall Constant terms
  const long double K=2./M_PI;
  //Expression of F
  
  long double F_ris=K*pow(1.-x1-x2,2)/(4.*x1*x2);
  return F_ris;
}
long double D_F_MtopInf(long double x1, long double r, long double xp, long double y)
{
  //definition of x2
  long double x2=(1.+sqrt(x1))*(1.+sqrt(x1))-(4.*sqrt(x1))*r;
  //Overall Constant terms
  const long double K=2./M_PI;
  //Expression of F
  
  long double F_ris=K*(4.*r*(1.-r)*(1.+sqrt(x1)-2.*r))/(sqrt(x1)*x2*x2);
  return F_ris;
}   
  
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  



