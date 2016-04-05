#include "NLO_PL.h"
#include "cuba.h"
#include <gsl/gsl_sf_dilog.h>

/*double B2pp_ANALITIC(long double x, long double xp)
{
  return(1./3.*3.*5.*(2.+3.*xp)/(xp*(1.+xp))*std::log(x));
}*/ //Used to cross-check the complete form above 

long double B2pp_ANALITICS(long double x, long double xp){
  long double ymax=0.5*std::log((1.+std::sqrt(1.-4.*x*(1.+xp)/std::pow(1.+x,2)))/(1.-std::sqrt(1.-4.*x*(1.+xp)/std::pow(1.+x,2))));
  long double wmax=std::exp(ymax);
  long double R1=1./12.*x*xp*(1./std::sqrt(x)*(1./wmax*(2.+3.*x)*std::sqrt(1.+xp))
		-2.*std::pow(wmax,2)*x*(1.+xp)-8.*std::pow(x,5)*xp*xp*std::pow(1.+xp,3)/(3.*std::pow(x+x*xp-wmax*std::sqrt(x*(1.+xp)),3))
		+4.*x/(wmax*std::sqrt(x*(1.+xp))-x)+wmax*std::sqrt(1.+xp)*(2.+3.*x-8.*x*x*(1.+2.*xp))/std::sqrt(x)
		+2.*x*x*x*xp*std::pow(1.+xp,2)*(1.+x*x*(1.+7.*xp)-x*(4.+7.*xp))/((x-1.)*std::pow(x+x*xp-wmax*std::sqrt(x*(1.+xp)),2))
		-(4.*x*(2.*(1.+xp)+2.*x*xp*(1.+xp)+x*x*(-3.+4.*xp+18.*xp*xp+11.*xp*xp*xp)+x*x*x*x*(2.+9.*xp+18.*xp*xp+11.*xp*xp*xp)-x*x*x*(1.+15.*xp+36.*xp*xp+22.*xp*xp*xp)))
		/(std::pow(x-1.,2)*(x+x*xp-wmax*std::sqrt(x*(1.+xp))))
		-(-2.+3.*x+x*x+2.*xp-3.*x*xp)*std::log(wmax*std::sqrt(x*(1.+xp)))/x
		-((2.+2.*xp-2.*xp*xp*xp+std::pow(x,6)*(2.+xp)-3.*std::pow(x,5)*(4.+3.*xp)+std::pow(x,4)*(30.+26.*xp+4.*xp*xp-3.*xp*xp*xp)
		-x*(12.+23.*xp+4.*xp*xp+3.*xp*xp*xp)-std::pow(x,3)*(40.+48.*xp+9.*xp*xp*xp)+x*x*(30.+51.*xp+9.*xp*xp*xp))*std::log(1.+xp-wmax*std::sqrt(x*(1.+xp))))
		/(std::pow(x-1.,3)*x*(x-xp-1.)*xp)
		-1./(x*xp*(x-xp-1.))*(4.*x*x*x*x+x*x*x*(10.+xp)-2.*x*x*(3.+5.*xp)+2.*(1.-xp+xp*xp*xp)-x*(2.+7.*xp+3.*xp*xp*xp))*std::log(-x+wmax*std::sqrt(x*(1.+xp)))
		-1./(std::pow(x-1.,3)*x*xp)*4.*(-1.+2.*x*(2.+xp)-x*x*(5.+6.*xp+3.*xp*xp)-6.*std::pow(x,5)*xp*(2.+6.*xp+5.*xp*xp)
		-x*x*x*xp*(1.+8.*xp+10.*xp*xp)+std::pow(x,6)*(-1.+3.*xp+12.*xp*xp+10.*xp*xp*xp)
		+std::pow(x,4)*(3.+14.*xp+33.*xp*xp+30.*xp*xp*xp))*std::log(wmax*std::sqrt(x*(1.+xp))-x*(1.+xp)));
  
  ymax=std::log(((1.+x)-std::sqrt(std::pow(1.-x,2)-4.*x*xp))/(2.*std::sqrt(x*(1.+xp))));
  wmax=std::exp(ymax);
  
  long double R2=1./12.*x*xp*(1./std::sqrt(x)*(1./wmax*(2.+3.*x)*std::sqrt(1.+xp))
		-2.*std::pow(wmax,2)*x*(1.+xp)-8.*std::pow(x,5)*xp*xp*std::pow(1.+xp,3)/(3.*std::pow(x+x*xp-wmax*std::sqrt(x*(1.+xp)),3))
		+4.*x/(wmax*std::sqrt(x*(1.+xp))-x)+wmax*std::sqrt(1.+xp)*(2.+3.*x-8.*x*x*(1.+2.*xp))/std::sqrt(x)
		+2.*x*x*x*xp*std::pow(1.+xp,2)*(1.+x*x*(1.+7.*xp)-x*(4.+7.*xp))/((x-1.)*std::pow(x+x*xp-wmax*std::sqrt(x*(1.+xp)),2))
		-(4.*x*(2.*(1.+xp)+2.*x*xp*(1.+xp)+x*x*(-3.+4.*xp+18.*xp*xp+11.*xp*xp*xp)+x*x*x*x*(2.+9.*xp+18.*xp*xp+11.*xp*xp*xp)-x*x*x*(1.+15.*xp+36.*xp*xp+22.*xp*xp*xp)))
		/(std::pow(x-1.,2)*(x+x*xp-wmax*std::sqrt(x*(1.+xp))))
		-(-2.+3.*x+x*x+2.*xp-3.*x*xp)*std::log(wmax*std::sqrt(x*(1.+xp)))/x
		-((2.+2.*xp-2.*xp*xp*xp+std::pow(x,6)*(2.+xp)-3.*std::pow(x,5)*(4.+3.*xp)+std::pow(x,4)*(30.+26.*xp+4.*xp*xp-3.*xp*xp*xp)
		-x*(12.+23.*xp+4.*xp*xp+3.*xp*xp*xp)-std::pow(x,3)*(40.+48.*xp+9.*xp*xp*xp)+x*x*(30.+51.*xp+9.*xp*xp*xp))*std::log(1.+xp-wmax*std::sqrt(x*(1.+xp))))
		/(std::pow(x-1.,3)*x*(x-xp-1.)*xp)
		-1./(x*xp*(x-xp-1.))*(4.*x*x*x*x+x*x*x*(10.+xp)-2.*x*x*(3.+5.*xp)+2.*(1.-xp+xp*xp*xp)-x*(2.+7.*xp+3.*xp*xp*xp))*std::log(-x+wmax*std::sqrt(x*(1.+xp)))
		-1./(std::pow(x-1.,3)*x*xp)*4.*(-1.+2.*x*(2.+xp)-x*x*(5.+6.*xp+3.*xp*xp)-6.*std::pow(x,5)*xp*(2.+6.*xp+5.*xp*xp)
		-x*x*x*xp*(1.+8.*xp+10.*xp*xp)+std::pow(x,6)*(-1.+3.*xp+12.*xp*xp+10.*xp*xp*xp)
		+std::pow(x,4)*(3.+14.*xp+33.*xp*xp+30.*xp*xp*xp))*std::log(wmax*std::sqrt(x*(1.+xp))-x*(1.+xp)));
  
  
  long double R3=-1./6.*(1.+x*x-x*(2.+xp))*(std::log(x)+std::log(1.+xp)+2.*std::log(2)-2.*std::log(1.+x+std::sqrt((1.-x)*(1.-x)-4.*x*xp)));
  
  return(2./xp*(R1-R2+R3));
  
  
}

//Delta contribution

//gg channel

long double NLOPL::NLO_PL_delta(double xp){
  const long double rad=std::sqrt((1.-x)*(1.-x)-4.*x*xp);
  const long double t=0.5*(x-1.+rad);
  const long double u=0.5*(x-1.-rad);
  const long double MUR=pow(_muR/_mH,2.);
  const long double b0=11./6.*_Nc-1./3.*_Nf;
  const long double Delta=(5.*_Nc-3.*(_Nc*_Nc-1.)/(2.*_Nc));
  const long double delta=3./2.*b0*(std::log(MUR*x/(-t))+std::log(MUR*x/(-u)))+(67./18.*_Nc-5./9.*_Nf);
  const long double U=1./2.*std::pow(std::log(u/t),2.)+M_PIl*M_PIl/3.-std::log(x)*std::log(x/(-t))-std::log(x)*std::log(x/(-u))
		      -std::log(x/(-t))*std::log(x/(-u))+std::pow(std::log(x),2.)+std::pow(std::log(x/(x-t)),2.)+std::pow(std::log(x/(x-u)),2.)
		      +2.*gsl_sf_dilog(1.-x)+2.*gsl_sf_dilog(x/(x-t))+2.*gsl_sf_dilog(x/(x-u));
  const long double ris=x*(Delta+delta+_Nc*U)*_Nc*(pow(x,4)+1.+pow(t,4)+pow(u,4))/(u*t)
			+(_Nc-_Nf)*_Nc/3.*(x*x+(x*x/t)+(x*x/u)+x);			  
  const long double Jac1=-(1.-x-rad)/(2.*rad);
  const long double Jac2=(1.-x+rad)/(2.*rad);
  const long double Si5z1=1./t*(std::pow(x,4)+1.+std::pow(t,4)+std::pow(u,4))/(u*t)*std::log(x/(-t));
  const long double Si5z2=1./u*(std::pow(x,4)+1.+std::pow(u,4)+std::pow(t,4))/(u*t)*std::log(x/(-u));

  return (ris/rad*2.+2.*x*_Nc*b0*(Jac2*Si5z2-Jac1*Si5z1)+_Nc*_Nf*B2pp_ANALITICS(x,xp));
}


long double NLOPL::NLO_PL_sing_doublediff(double xp, double zz)
{
  long double ris=0.0;
  //Set Energy, Transverse Momentum, Integration variable za
  const long double xx=x;
  const long double Nc=_Nc;
  const long double Nf=_Nf;
  const long double muF=_muF;
  const long double muR=_muR;
  const long double b0=11./6.*Nc-1./3.*Nf;
  const long double zmin=xx*std::pow(std::sqrt(1.+xp)+std::sqrt(xp),2);
  const long double z=zmin+(1.-zmin)*zz;
  
  //Set Mandelstam variable
  const long double  rad=std::sqrt((z-xx)*(z-xx)-4.*xx*xp*z);
  const long double  t1=0.5*(xx-z+rad);
  const long double  u1=xx-0.5/z*(xx+z+rad);
  const long double  Q1=1.+t1+u1-xx;
  const long double  Qt1=Q1+xx*xp;
  const long double  zb1=(2.*xx*z-z-xx-rad)/(z*(-2.+xx+z-rad));
  const long double  Jac1=(xx-z+rad)/(2.*z*rad);
  const long double  Jac1z1=(xx-1.+std::sqrt((1.-xx)*(1.-xx)-4.*xx*xp))/(2.*std::sqrt((1.-xx)*(1.-xx)-4.*xx*xp));
  const long double  t2=0.5*(xx-z-rad);
  const long double  u2=xx-0.5/z*(xx+z-rad);
  const long double  Q2=1.+t2+u2-xx;
  const long double  Qt2=Q2+xx*xp;
  const long double  zb2=(2.*xx*z-xx-z+rad)/(z*(-2.+xx+z+rad));
  const long double  Jac2=-(xx-z-rad)/(2.*z*rad);
  const long double  Jac2z1=-(xx-1.-std::sqrt((1.-xx)*(1.-xx)-4.*xx*xp))/(2.*std::sqrt((1.-xx)*(1.-xx)-4.*xx*xp));
  
  //Define and add the different singular parts 
  const long double Si12=2.*xx/(-t2)*(1.+std::pow(z,4)+std::pow(1.-z,4))/(z)*Nc*Nc*(std::pow(xx,4)+std::pow(z,4)+std::pow(t2,4)+std::pow(xx*xp*z/t2,4))/(z*z*xx*xp);
  const long double Si11=2.*xx/(-t1)*(1.+std::pow(z,4)+std::pow(1.-z,4))/(z)*Nc*Nc*(std::pow(xx,4)+std::pow(z,4)+std::pow(t1,4)+std::pow(xx*xp*z/t1,4))/(z*z*xx*xp);
  const long double Si12z1=16.*Nc*Nc*(1.+std::pow(xx,4)-2.*xx*(1.+xp)-2*xx*xx*xx*(1.+xp)+xx*xx*(3.+4.*xp+xp*xp))/(xp*(1.-xx+sqrt(1.-2.*xx+xx*xx-4.*xx*xp)));
  const long double Si11z1=16.*Nc*Nc*(1.+std::pow(xx,4)-2.*xx*(1.+xp)-2*xx*xx*xx*(1.+xp)+xx*xx*(3.+4.*xp+xp*xp))/(xp*(1.-xx-sqrt(1.-2.*xx+xx*xx-4.*xx*xp)));
  ris+=std::log(1.-z)/(1.-z)*(Jac2*Si12-Jac2z1*Si12z1-Jac1*Si11+Jac1z1*Si11z1);
  
  const long double Si21=2.*xx*(z/(-t1))*Nc*Nc/2.*((std::pow(xx,4)+1.+std::pow(Q1,4)+std::pow(u1,4)+std::pow(t1,4))+z*zb1*(std::pow(xx,4)
		  +1.+std::pow(Q1,4)+std::pow(u1/zb1,4)+std::pow(t1/z,4)))/(u1*t1);
  const long double Si22=2.*xx*(z/(-t2))*Nc*Nc/2.*((std::pow(xx,4)+1.+std::pow(Q2,4)+std::pow(u2,4)+std::pow(t2,4))+z*zb2*(std::pow(xx,4)
		  +1.+std::pow(Q2,4)+std::pow(u2/zb2,4)+std::pow(t2/z,4)))/(u2*t2);
  const long double Si21z1=(8.*Nc*Nc*(-std::pow(1.+(xx-1.)*xx,2)+2.*std::pow(1.-xx,2)*xx*xp-xx*xx*xp*xp)/(xp*(-1.+xx+sqrt(std::pow(1.-xx,2)-4.*xx*xp))));
  const long double Si22z1=(8.*Nc*Nc*(-std::pow(1.+(xx-1.)*xx,2)+2.*std::pow(1.-xx,2)*xx*xp-xx*xx*xp*xp)/(xp*(-1.+xx-sqrt(std::pow(1.-xx,2)-4.*xx*xp))));
 ris+=std::log(1.-z)/(1.-z)*(Jac2*Si22-Jac2z1*Si22z1-Jac1*Si21+Jac1z1*Si21z1);
  
  const long double Si31=2.*xx*(z/(-t1))*Nc*Nc/2.*((std::pow(xx,4)+1.+std::pow(Q1,4)+std::pow(u1,4)+std::pow(t1,4))+z*zb1*(std::pow(xx,4)
		  +1.+std::pow(Q1,4)+std::pow(u1/zb1,4)+std::pow(t1/z,4)))/(u1*t1)*(std::log(Qt1*z/(-t1)));
  const long double Si32=2.*xx*(z/(-t2))*Nc*Nc/2.*((std::pow(xx,4)+1.+std::pow(Q2,4)+std::pow(u2,4)+std::pow(t2,4))+z*zb2*(std::pow(xx,4)
		  +1.+std::pow(Q2,4)+std::pow(u2/zb2,4)+std::pow(t2/z,4)))/(u2*t2)*(std::log(Qt2*z/(-t2)));
  const long double Si31z1=-(8.*Nc*Nc*(std::pow(1.+(xx-1.)*xx,2)-2.*std::pow(1.-xx,2)*xx*xp+xx*xx*xp*xp)/(xp*(-1.+xx+sqrt(std::pow(1.-xx,2)-4.*xx*xp)))
		    *std::log(2.*xx*xp/(1.-xx-sqrt(std::pow(1.-xx,2)-4.*xx*xp))));
  const long double Si32z1=-(8.*Nc*Nc*(std::pow(1.+(xx-1.)*xx,2)-2.*std::pow(1.-xx,2)*xx*xp+xx*xx*xp*xp)/(xp*(-1.+xx-sqrt(std::pow(1.-xx,2)-4.*xx*xp)))
		    *std::log(2.*xx*xp/(1.-xx+sqrt(std::pow(1.-xx,2)-4.*xx*xp))));
  ris+=-1./(1.-z)*(Jac2*Si32-Jac2z1*Si32z1-Jac1*Si31+Jac1z1*Si31z1);
  
  const long double Si41=2.*xx*(z/t1)*b0*Nc/2.*(std::pow(xx,4)+1.+z*zb1*(std::pow(u1/zb1,4)+std::pow(t1/z,4)))/(u1*t1);
  const long double Si42=2.*xx*(z/t2)*b0*Nc/2.*(std::pow(xx,4)+1.+z*zb2*(std::pow(u2/zb2,4)+std::pow(t2/z,4)))/(u2*t2);
  const long double Si41z1=(4.*Nc*b0*(std::pow(1.+(xx-1.)*xx,2)-2.*std::pow(1.-xx,2)*xx*xp+xx*xx*xp*xp)/(xp*(-1.+xx+sqrt(std::pow(1.-xx,2)-4.*xx*xp))));
  const long double Si42z1=(4.*Nc*b0*(std::pow(1.+(xx-1.)*xx,2)-2.*std::pow(1.-xx,2)*xx*xp+xx*xx*xp*xp)/(xp*(-1.+xx-sqrt(std::pow(1.-xx,2)-4.*xx*xp)))); 
  ris+=1./(1.-z)*(Jac2*Si42-Jac2z1*Si42z1-Jac1*Si41+Jac1z1*Si41z1);
  
  const long double Si51=2.*xx/t1*(1.+std::pow(z,4)+std::pow(1.-z,4))/z*std::log(muF*xx*z/(-t1))*Nc*Nc*(std::pow(xx,4)+std::pow(z,4)+std::pow(t1,4)+std::pow(xx*xp*z/t1,4))/(z*z*xx*xp);
  const long double Si52=2.*xx/t2*(1.+std::pow(z,4)+std::pow(1.-z,4))/z*std::log(muF*xx*z/(-t2))*Nc*Nc*(std::pow(xx,4)+std::pow(z,4)+std::pow(t2,4)+std::pow(xx*xp*z/t2,4))/(z*z*xx*xp);
  const long double Si51z1=(16.*Nc*Nc*(std::pow(1.+(xx-1.)*xx,2)-2.*std::pow(1.-xx,2)*xx*xp+xx*xx*xp*xp)/(xp*(-1.+xx+sqrt(std::pow(1.-xx,2)-4.*xx*xp)))*std::log(muF*2.*xx/(1.-xx-sqrt(std::pow(1.-xx,2)-4.*xx*xp))));
  const long double Si52z1=(16.*Nc*Nc*(std::pow(1.+(xx-1.)*xx,2)-2.*std::pow(1.-xx,2)*xx*xp+xx*xx*xp*xp)/(xp*(-1.+xx-sqrt(std::pow(1.-xx,2)-4.*xx*xp)))*std::log(muF*2.*xx/(1.-xx+sqrt(std::pow(1.-xx,2)-4.*xx*xp))));
  ris+=1./(1.-z)*(Jac2*Si52-Jac2z1*Si52z1-Jac1*Si51+Jac1z1*Si51z1);
  ris*=(1.-zmin);
  return ris;
}

long double NLOPL::NLO_PL_notsing_doublediff(double xp, double zz){
  
  long double ris=0.0;
  
  //Set Energy, Transverse Momentum, Integration variable z
  const long double xx=x;
  const long double Nc=_Nc;
  const long double Nf=_Nf;
  const long double muF=_muF;
  const long double muR=_muR;
  const long double b0=11./6.*Nc-1./3.*Nf;
  const long double Cf=(Nc*Nc-1.)/(2.*Nc);
  const long double zmin=xx*std::pow(std::sqrt(1.+xp)+std::sqrt(xp),2);
  const long double z=zmin+(1.-zmin)*zz;
  
  //Set Mandelstam variable
  const long double  rad=std::sqrt((z-xx)*(z-xx)-4.*xx*xp*z);
  const long double  t1=0.5*(xx-z+rad);
  const long double  u1=xx-0.5/z*(xx+z+rad);
  long double  Q1=1.+t1+u1-xx;
  const long double  Qt1=Q1+xx*xp;
  const long double  zb1=(2.*xx*z-z-xx-rad)/(z*(-2.+xx+z-rad));
  const long double  Jac1=(xx-z+rad)/(2.*z*rad);
  const long double  t2=0.5*(xx-z-rad);
  const long double  u2=xx-0.5/z*(xx+z-rad);
  long double  Q2=1.+t2+u2-xx;
  const long double  Qt2=Q2+xx*xp;
  const long double  zb2=(2.*xx*z-xx-z+rad)/(z*(-2.+xx+z+rad));
  const long double  Jac2=-(xx-z-rad)/(2.*z*rad);
  const long double  A1=1.+xx-Q1;
  const long double  B1=std::sqrt(A1*A1-4.*xx);
  const long double  L1a1=std::log(xx/(z*z));
  const long double  L1b1=std::log(xx/(zb1*zb1));
  const long double  L2a1=std::log(xx/std::pow(A1-z,2));
  const long double  L2b1=std::log(xx/std::pow(A1-zb1,2));
  const long double  L31=std::log((A1+B1)/(A1-B1));
  const long double  A2=1.+xx-Q2;
  const long double  B2=std::sqrt(A2*A2-4.*xx);
  const long double  L1a2=std::log(xx/(z*z));
  const long double  L1b2=std::log(xx/(zb2*zb2));
  const long double  L2a2=std::log(xx/std::pow(A2-z,2));
  const long double  L2b2=std::log(xx/std::pow(A2-zb2,2));
  const long double  L32=std::log((A2+B2)/(A2-B2));
  
  //Define and add regular parts in G2s (paper Glover for definition of G2s)
  const long double Fi1_1=xx/(-t1)*(-2.*Nf*std::log(muF*xx/Q1)*0.5*(z*z+(1.-z)*(1.-z))+2.*Nf*z*(1.-z))*Cf*(z*z*xx*xx*xp*xp/(t1*t1)+z*z)/(-t1);
  const long double Fi1_2=xx/(-t2)*(-2.*Nf*std::log(muF*xx/Q2)*0.5*(z*z+(1.-z)*(1.-z))+2.*Nf*z*(1.-z))*Cf*(z*z*xx*xx*xp*xp/(t2*t2)+z*z)/(-t2);			
  ris+=2.*(Jac2*Fi1_2-Jac1*Fi1_1);
  
  
  //Separate divergent part Q->0 (ref Grazzini)
  const long double tiny=1e-4;
  const long double Fi2_1 = (Q1 > tiny*(xx*xp) ) ? Nc*Nc*((std::pow(xx,4)+1.+std::pow(Q1,4)+std::pow(u1/zb1,4)+std::pow(t1/z,4))*(Q1+Qt1)/(Q1*Qt1)
			+(2.*xx*xx*(std::pow(xx-t1,4)+std::pow(xx-u1,4)+std::pow(u1,4)+std::pow(t1,4)))/(u1*t1*(xx-u1)*(xx-t1)))*1./(xp)
			*std::log(xx*xp/Qt1) : -Nc*Nc*((std::pow(xx,4)+1.+std::pow(u1,4)+std::pow(t1,4))/(xx*xp*xp));
  const long double Fi2_2 = ( Q2 > tiny*(xx*xp)) ? Nc*Nc*((std::pow(xx,4)+1.+std::pow(Q2,4)+std::pow(u2/zb2,4)+std::pow(t2/z,4))*(Q2+Qt2)/(Q2*Qt2)
			+(2.*xx*xx*(std::pow(xx-t2,4)+std::pow(xx-u2,4)+std::pow(u2,4)+std::pow(t2,4)))/(u2*t2*(xx-u2)*(xx-t2)))*1./(xp)
			*std::log(xx*xp/Qt2) : -Nc*Nc*((std::pow(xx,4)+1.+std::pow(u2,4)+std::pow(t2,4))/(xx*xp*xp));
  ris+=(Jac2*Fi2_2-Jac1*Fi2_1);
  
  //Define and add G2ns (paper Glover for definition of G2ns)
  //A1234 
  const long double A1234_1_1=-0.5*((pow(xx*xp/t1,4)+pow(xx*Q1/t1,4))*L1a1+pow(u1,4)*L1b1);
  const long double A1234_1_2=(pow(xx,2)*Q1*pow(u1,3))/(2.*t1*(u1-xx)*(t1-xx))*(L1a1+L1b1)
			     +(xx*xp*Q1*pow(u1,3))/(2.*A1*(t1-xx))*(L2b1-L1b1)
			     +(xx*xp*Q1*(pow(xx,4)+pow(u1-xx,4)))/(2.*A1*t1*(u1-xx))*(L2a1-L1a1);
  const long double A1234_1_3=pow(xx,2)*xp*Q1*(pow(u1,4)/(2.*pow(B1*t1,2))+pow(u1,2)/(2.*pow(B1,2)))
			     +pow(xx,4)*pow(xp*Q1,2)*(-6./pow(B1,4)-4./pow(t1,4)+8./pow(B1*t1,2));
  const long double A1234_1_4=L31*(xx*xp*pow(u1,3)*(u1+t1)/(B1*t1)+pow(xx,4)*pow(Q1*xp,2)*(3.*A1/pow(B1,5)-1./(A1*pow(B1,3)))
			     -pow(xx,2)*xp*Q1*(1./(B1*t1)*(t1*t1+t1*u1+4.*pow(u1,2)-2.*xx*Q1)
			     +A1/(2.*pow(B1,3))*(t1*t1+3.*t1*u1+3.*u1*u1-6*xx*Q1)+1./(2.*A1*B1)*(t1*t1+t1*u1+7.*u1*u1-2.*xx*Q1)));
  const long double A1234_2_1=-0.5*((pow(xx*xp/t2,4)+pow(xx*Q2/t2,4))*L1a2+pow(u2,4)*L1b2);
  const long double A1234_2_2=(pow(xx,2)*Q2*pow(u2,3))/(2.*t2*(u2-xx)*(t2-xx))*(L1a2+L1b2)
			     +(xx*xp*Q2*pow(u2,3))/(2.*A2*(t2-xx))*(L2b2-L1b2)
			     +(xx*xp*Q2*(pow(xx,4)+pow(u2-xx,4)))/(2.*A2*t2*(u2-xx))*(L2a2-L1a2);
  const long double A1234_2_3=pow(xx,2)*xp*Q2*(pow(u2,4)/(2.*pow(B2*t2,2))+pow(u2,2)/(2.*pow(B2,2)))
			     +pow(xx,4)*pow(xp*Q2,2)*(-6./pow(B2,4)-4./pow(t2,4)+8./pow(B2*t2,2));
  const long double A1234_2_4=L32*(xx*xp*pow(u2,3)*(u2+t2)/(B2*t2)+pow(xx,4)*pow(Q2*xp,2)*(3.*A2/pow(B2,5)-1./(A2*pow(B2,3)))
			     -pow(xx,2)*xp*Q2*(1./(B2*t2)*(t2*t2+t2*u2+4.*pow(u2,2)-2.*xx*Q2)
			     +A2/(2.*pow(B2,3))*(t2*t2+3.*t2*u2+3.*u2*u2-6*xx*Q2)+1./(2.*A2*B2)*(t2*t2+t2*u2+7.*u2*u2-2.*xx*Q2)));
  ris+=2.*Nc*Nc*(Jac2/(xp*Q2)*(A1234_2_1+A1234_2_2+A1234_2_3+A1234_2_4)-Jac1/(xp*Q1)*(A1234_1_1+A1234_1_2+A1234_1_3+A1234_1_4));
  
  //A3412
  const long double A3412_1_1=xx*xp*Q1*pow(A1,3)/(2.*t1*(u1-xx))*(L2a1+L1b1)
			      +xx*xp*(u1+t1)/(16.*u1*t1*B1)*(pow(A1,4)+6.*pow(A1*B1,2)+pow(B1,4))*L31;
  const long double A3412_1_2=(-xx*xp/(2.*u1*t1)*(pow(1.-Q1,4)+pow(xx,4)+2.*Q1*A1*pow(1.-Q1,2)-2.*Q1*pow(xx,3))
			      -pow(xx,4)*pow(Q1*xp,2)/pow(u1,4)+(2.*pow(xx,2)*xp*Q1*(A1*A1-xx))/(u1*u1))*L1b1;
  const long double A3412_1_3=xx*xp*pow(Q1-u1,3)/(8.*u1*t1*(Q1-t1))*((Q1-t1)+Q1*u1)*
			      (4./3.+2.*xx*xp/Qt1+4.*pow(xx*xp/Qt1,2)-44./3.*pow(xx*xp/Qt1,3))
			      +xx*xp*pow(Q1-t1,3)/(8.*u1*t1*(Q1-u1))*((Q1-u1)+Q1*t1)*
			      (4./3.+2.*xx*xp/Qt1+4.*pow(xx*xp/Qt1,2)-44./3.*pow(xx*xp/Qt1,3));
  const long double A3412_1_4=xx*xp*pow(Q1-u1,2)/(4.*u1*t1*(Q1-t1))*
			      (-3.*(t1-xx)*((Q1-t1)+Q1*u1)-Q1*(xx*(t1-xx)+Q1*(u1-xx)))*(1.+2.*xx*xp/Qt1-6.*pow(xx*xp/Qt1,2))
			      +xx*xp*pow(Q1-t1,2)/(4.*u1*t1*(Q1-u1))*
			      (-3.*(u1-xx)*((Q1-u1)+Q1*t1)+3.*Q1*(xx*(u1-xx)+Q1*(t1-xx))+4.*u1)*(1.+2.*xx*xp/Qt1-6.*pow(xx*xp/Qt1,2));
  const long double A3412_1_5=xx*xp*(Q1-t1)/(2.*u1*t1*(Q1-u1))*(3.*pow(u1-xx,2)*((Q1-u1)+Q1*t1)+8.*u1*t1+2.*u1-2.*Q1*u1*pow(u1-Q1,2)
			      -3.*xx*Q1*pow(t1-xx,2)-3.*Q1*(xx-Q1)*pow(t1,2)-Q1*u1*(4.*u1*t1-u1*xx-Q1*t1+2.*pow(t1,2)-4.*xx*xx)+3.*xx*pow(Q1,2)*(t1-xx)
			      +xx*Q1*u1*(t1-Q1))*(1.-2.*xx*xp/Qt1)+xx*xp*(Q1-u1)/(2.*u1*t1*(Q1-t1))*
			      (3.*pow(t1-xx,2)*((Q1-t1)+Q1*u1)+3.*(t1-xx)*Q1*(xx*(t1-xx)+Q1*(u1-xx))+Q1*u1*(xx*(t1-Q1)+Q1*(u1-xx)))*(1.-2.*xp*xx/Qt1);
  const long double A3412_1_6=-4.*pow(xx,4)*pow(Q1*xp,2)/pow(u1,4)+xp*Q1*pow(xx*B1,2)/(2.*u1*u1)
			      +pow(xx,3)*xp/6.*((1.+Q1)/(u1*t1)+Q1/(u1*u1)+Q1/(t1*t1))+(2.*xp*pow(xx,3)*Q1)/(u1*u1)+pow(xx,3)*xp/u1
			      -xx*xp/(12.*t1*u1)*(30.*pow(xx,3)+54.*pow(Q1,2)*xx+8.*pow(Q1,3))
			      +xx*xp/(12*u1*t1)*(11.+17.*pow(xx,4)+Q1*(61.*u1*u1*t1+17.*pow(u1,3)+73.*u1*t1*t1+29.*pow(t1,3))
			      +xx*(24.*u1*u1*t1+6.*pow(u1,3)+36.*u1*t1*t1+18.*pow(t1,3))+Q1*Q1*(-21.*u1*u1-33.*t1*t1-52.*u1*t1)
			      +xx*Q1*(-73.*u1*u1-109.*t1*t1-170.*u1*t1)+xx*xx*(-23.*u1*u1-35.*t1*t1-52.*u1*t1)
			      +xx*xx*Q1*(134.*t1+110.*u1)+4.*pow(Q1,4)+52.*xx*pow(Q1,3)+20.*xx*xx*Q1*Q1-22.*pow(xx,3)*Q1);
  const long double A3412_2_1=xx*xp*Q2*pow(A2,3)/(2.*t2*(u2-xx))*(L2a2+L1b2)
			      +xx*xp*(u2+t2)/(16.*u2*t2*B2)*(pow(A2,4)+6.*pow(A2*B2,2)+pow(B2,4))*L32;
  const long double A3412_2_2=(-xx*xp/(2.*u2*t2)*(pow(1.-Q2,4)+pow(xx,4)+2.*Q2*A2*pow(1.-Q2,2)-2.*Q2*pow(xx,3))
			      -pow(xx,4)*pow(Q2*xp,2)/pow(u2,4)+(2.*pow(xx,2)*xp*Q2*(A2*A2-xx))/(u2*u2))*L1b2;
  const long double A3412_2_3=xx*xp*pow(Q2-u2,3)/(8.*u2*t2*(Q2-t2))*((Q2-t2)+Q2*u2)*
			      (4./3.+2.*xx*xp/Qt2+4.*pow(xx*xp/Qt2,2)-44./3.*pow(xx*xp/Qt2,3))
			      +xx*xp*pow(Q2-t2,3)/(8.*u2*t2*(Q2-u2))*((Q2-u2)+Q2*t2)*
			      (4./3.+2.*xx*xp/Qt2+4.*pow(xx*xp/Qt2,2)-44./3.*pow(xx*xp/Qt2,3));
  const long double A3412_2_4=xx*xp*pow(Q2-u2,2)/(4.*u2*t2*(Q2-t2))*
			      (-3.*(t2-xx)*((Q2-t2)+Q2*u2)-Q2*(xx*(t2-xx)+Q2*(u2-xx)))*(1.+2.*xx*xp/Qt2-6.*pow(xx*xp/Qt2,2))
			      +xx*xp*pow(Q2-t2,2)/(4.*u2*t2*(Q2-u2))*
			      (-3.*(u2-xx)*((Q2-u2)+Q2*t2)+3.*Q2*(xx*(u2-xx)+Q2*(t2-xx))+4.*u2)*(1.+2.*xx*xp/Qt2-6.*pow(xx*xp/Qt2,2));
  const long double A3412_2_5=xx*xp*(Q2-t2)/(2.*u2*t2*(Q2-u2))*(3.*pow(u2-xx,2)*((Q2-u2)+Q2*t2)+8.*u2*t2+2.*u2-2.*Q2*u2*pow(u2-Q2,2)
			      -3*xx*Q2*pow(t2-xx,2)-3.*Q2*(xx-Q2)*pow(t2,2)-Q2*u2*(4.*u2*t2-u2*xx-Q2*t2+2.*pow(t2,2)-4.*xx*xx)+3.*xx*pow(Q2,2)*(t2-xx)
			      +xx*Q2*u2*(t2-Q2))*(1-2.*xx*xp/Qt2)+xx*xp*(Q2-u2)/(2.*u2*t2*(Q2-t2))*
			      (3.*pow(t2-xx,2)*((Q2-t2)+Q2*u2)+3.*(t2-xx)*Q2*(xx*(t2-xx)+Q2*(u2-xx))+Q2*u2*(xx*(t2-Q2)+Q2*(u2-xx)))*(1-2.*xp*xx/Qt2);
  const long double A3412_2_6=-4.*pow(xx,4)*pow(Q2*xp,2)/pow(u2,4)+xp*Q2*pow(xx*B2,2)/(2.*u2*u2)
			      +pow(xx,3)*xp/6.*((1.+Q2)/(u2*t2)+Q2/(u2*u2)+Q2/(t2*t2))+(2.*xp*pow(xx,3)*Q2)/(u2*u2)+pow(xx,3)*xp/u2
			      -xx*xp/(12.*t2*u2)*(30.*pow(xx,3)+54.*pow(Q2,2)*xx+8.*pow(Q2,3))
			      +xx*xp/(12*u2*t2)*(11.+17.*pow(xx,4)+Q2*(61.*u2*u2*t2+17.*pow(u2,3)+73.*u2*t2*t2+29.*pow(t2,3))
			      +xx*(24.*u2*u2*t2+6.*pow(u2,3)+36.*u2*t2*t2+18.*pow(t2,3))+Q2*Q2*(-21.*u2*u2-33.*t2*t2-52.*u2*t2)
			      +xx*Q2*(-73.*u2*u2-109.*t2*t2-170.*u2*t2)+xx*xx*(-23.*u2*u2-35.*t2*t2-52.*u2*t2)
			      +xx*xx*Q2*(134.*t2+110.*u2)+4.*pow(Q2,4)+52.*xx*pow(Q2,3)+20.*xx*xx*Q2*Q2-22.*pow(xx,3)*Q2);
  ris+=2.*Nc*Nc*(Jac2/(xp*Q2)*(A3412_2_1+A3412_2_2+A3412_2_3+A3412_2_4+A3412_2_5+A3412_2_6)-Jac1/(xp*Q1)*(A3412_1_1+A3412_1_2+A3412_1_3+A3412_1_4+A3412_1_5+A3412_1_6));
  
  //A1324
  const long double A1324_1_1=-0.5*((pow(xx*xp/t1,4)+pow(xx*Q1/t1,4))*L1a1+pow(u1,4)*L1b1)
			      +pow(xx,3)*xp*Q1/(t1*t1)*L1a1+(pow(xx,2)*Q1*pow(u1,3)/(2.*t1*(u1-xx)*(t1-xx))+xx*xp*pow(u1,3)/(2.*t1))
			      *(L1a1+L1b1);
  const long double A1324_1_2=xx*xp*(1.-zb1)*pow(u1,3)/(2.*A1*(t1-xx))*(L2b1-L1b1)
			      +(xx*xp*(1.-z)*(pow(xx,4)+pow(u1-xx,4)))/(2.*A1*t1*(u1-xx))*(L2a1-L1a1)
			      +pow(xx,3)*xp*Q1/(A1*B1)*L31+xx*xx*xp*Q1/(2.*pow(t1,4))*(pow(xx*xp,2)-6.*xx*xx*xp*Q1+pow(xx*Q1,2));
  const long double A1324_2_1=-0.5*((pow(xx*xp/t2,4)+pow(xx*Q2/t2,4))*L1a2+pow(u2,4)*L1b2)
			      +pow(xx,3)*xp*Q2/(t2*t2)*L1a2+(pow(xx,2)*Q2*pow(u2,3)/(2.*t2*(u2-xx)*(t2-xx))+xx*xp*pow(u2,3)/(2.*t2))
			      *(L1a2+L1b2);
  const long double A1324_2_2=xx*xp*(1.-zb2)*pow(u2,3)/(2.*A2*(t2-xx))*(L2b2-L1b2)
			      +(xx*xp*(1.-z)*(pow(xx,4)+pow(u2-xx,4)))/(2.*A2*t2*(u2-xx))*(L2a2-L1a2)
			      +pow(xx,3)*xp*Q2/(A2*B2)*L32+xx*xx*xp*Q2/(2.*pow(t2,4))*(pow(xx*xp,2)-6.*xx*xx*xp*Q2+pow(xx*Q2,2));
  ris+=2.*Nc*Nc*(Jac2/(xp*Q2)*(A1324_2_1+A1324_2_2)-Jac1/(xp*Q1)*(A1324_1_1+A1324_1_2));
  
  //A3241
  const long double A3241_1_1=xx*xp*pow(A1,3)*(1.-z)/(2.*t1*(u1-xx))*(L2a1-L1a1)
			      +(-pow(xx,4)*pow(xp*Q1,2)/pow(t1,4)+pow(xx,3)*xp*pow(Q1,2)/(u1*t1)
			      -xx*xx*xp*Q1*pow(A1,4)/(2.*u1*t1*(u1-xx)*(t1-xx))+xx*xx*xp*Q1*(u1+t1)*(2.*A1*A1-xx)/(u1*t1*t1))*L1a1;
  const long double A3241_1_2=xx*xp*Q1*pow(Q1-u1,2)/(2.*u1*t1*pow(Q1-t1,2))*
			      (-u1*t1-pow(Q1-t1,2))*(-3.+10.*Q1/Qt1-6.*pow(Q1/Qt1,2))
			      +xx*xp*Q1*(Q1-u1)/(u1*t1*pow(Q1-t1,2))*(u1*t1*(Q1-u1)-pow(Q1-t1,3)-xx*pow(Q1-t1,2)-xx*(Q1-t1)*(Q1-u1))
			      *(-1.+2.*Q1/Qt1);
  const long double A3241_1_3=xx*xx*xp*Q1*(B1*B1/(2.*t1*t1)-2.*xx*Q1/(t1*t1)+(u1+t1)*(u1+t1)/(2.*u1*t1))
			      -4.*pow(xx,4)*pow(xp*Q1,2)/pow(t1,4)+xx*xp*Q1/(4.*u1*t1)
			      *(pow(t1+u1,2)-(t1+u1)*(6.*Q1+4.*xx)+6.*Q1*Q1+8.*xx*Q1)+pow(xx,3)*xp*Q1*pow(t1+u1,2)/(4.*u1*u1*t1*t1);
  const long double A3241_2_1=xx*xp*pow(A2,3)*(1.-z)/(2.*t2*(u2-xx))*(L2a2-L1a2)
			      +(-pow(xx,4)*pow(xp*Q2,2)/pow(t2,4)+pow(xx,3)*xp*pow(Q2,2)/(u2*t2)
			      -xx*xx*xp*Q2*pow(A2,4)/(2.*u2*t2*(u2-xx)*(t2-xx))+xx*xx*xp*Q2*(u2+t2)*(2.*A2*A2-xx)/(u2*t2*t2))*L1a2;
  const long double A3241_2_2=xx*xp*Q2*pow(Q2-u2,2)/(2.*u2*t2*pow(Q2-t2,2))*
			      (-u2*t2-pow(Q2-t2,2))*(-3.+10.*Q2/Qt2-6.*pow(Q2/Qt2,2))
			      +xx*xp*Q2*(Q2-u2)/(u2*t2*pow(Q2-t2,2))*(u2*t2*(Q2-u2)-pow(Q2-t2,3)-xx*pow(Q2-t2,2)-xx*(Q2-t2)*(Q2-u2))
			      *(-1.+2.*Q2/Qt2);
  const long double A3241_2_3=xx*xx*xp*Q2*(B2*B2/(2.*t2*t2)-2.*xx*Q2/(t2*t2)+(u2+t2)*(u2+t2)/(2.*u2*t2))
			      -4.*pow(xx,4)*pow(xp*Q2,2)/pow(t2,4)+xx*xp*Q2/(4.*u2*t2)
			      *(pow(t2+u2,2)-(t2+u2)*(6.*Q2+4.*xx)+6.*Q2*Q2+8.*xx*Q2)+pow(xx,3)*xp*Q2*pow(t2+u2,2)/(4.*u2*u2*t2*t2);
  ris+=2.*Nc*Nc*(Jac2/(xp*Q2)*(A3241_2_1+A3241_2_2+A3241_2_3)-Jac1/(xp*Q1)*(A3241_1_1+A3241_1_2+A3241_1_3));
  
  //Aepsilon
  const long double Aepsilon_1=4.*pow(xx,4)*pow(xp*Q1,2)*(1./pow(t1,4)+1./pow(u1,4));
  const long double Aepsilon_2=4.*pow(xx,4)*pow(xp*Q2,2)*(1./pow(t2,4)+1./pow(u2,4));
  ris+=Nc*Nc*(Jac2/(xp*Q2)*(Aepsilon_2)-Jac1/(xp*Q1)*(Aepsilon_1));
  
  //A0
  const long double A0_1= (pow(t1/z,4)+pow(u1/zb1,4))*xx*xp/(Qt1*Qt1)*(5.-7.*Q1/Qt1+20./3.*pow(Q1/Qt1,2))
			 +xx*xp*(17./3.+4.*std::log(xx*xp/Qt1));
  const long double A0_2=(pow(t2/z,4)+pow(u2/zb2,4))*xx*xp/(Qt2*Qt2)*(5.-7.*Q2/Qt2+20./3.*pow(Q2/Qt2,2))
			 +xx*xp*(17./3.+4.*std::log(xx*xp/Qt2));
  ris+=Nc*Nc*(Jac2/(xp)*(A0_2)-Jac1/(xp)*(A0_1));
  
  //B1pm
  const long double B1pm_1=xx*xp*z*pow(1.-z,3)/t1+pow(xx*xp*z/t1,3)*(1.-z)
			   +4.*pow(xx*xp*z/t1*(1.-z),2)-xx*xp*Q1*(1.+std::log(xx*xp/Qt1));
  const long double B1pm_2=xx*xp*z*pow(1.-z,3)/t2+pow(xx*xp*z/t2,3)*(1.-z)
			   +4.*pow(xx*xp*z/t2*(1.-z),2)-xx*xp*Q2*(1.+std::log(xx*xp/Qt2));
  ris+=2.*Nf*Cf*(Jac2/(xp*Q2)*B1pm_2-Jac1/(xp*Q1)*B1pm_1);
  
  //B2pm //FIXME? C'è una differenza con risultato mathematica ma non capisco da dove dipenda. Contributo completamente negiglible comunque
  const long double B2pm_1=1./3.*pow(t1/z,4.)*xx*xp/Qt1*(-3./Qt1+3.*Q1/(Qt1*Qt1)-2.*Q1*Q1/(Qt1*Qt1*Qt1))
			   -1./3.*xx*xp;
  const long double B2pm_2=1./3.*pow(t2/z,4.)*xx*xp/Qt2*(-3./Qt2+3.*Q1/(Qt2*Qt2)-2.*Q2*Q2/(Qt2*Qt2*Qt2))
			   -1./3.*xx*xp;
  ris+=2.*Nf*Nc*(Jac2/xp*B2pm_2-Jac1/xp*B2pm_1);
  
  //B1pp
  const long double B1pp_1=xx*xp*pow(z,3)*(1.-z)/t1+pow(xx*xp*(1.-z)/t1,3)*z
			   +4.*pow(xx*xp*z*(1.-z)/t1,2)-xx*xp*Q1/(Qt1*Qt1*u1*t1)*(pow(u1*t1+xx*xp*Q1,2)+2.*xx*xp*Q1*Qt1)
			   +xx*xp*Q1/(u1*t1)*(1.+Q1*Q1)*std::log(1.+(xx*xp/Q1));
  const long double B1pp_2=xx*xp*pow(z,3)*(1.-z)/t2+pow(xx*xp*(1.-z)/t2,3)*z
			   +4.*pow(xx*xp*z*(1.-z)/t2,2)-xx*xp*Q2/(Qt2*Qt2*u2*t2)*(pow(u2*t2+xx*xp*Q2,2)+2.*xx*xp*Q2*Qt2)
			   +xx*xp*Q2/(u2*t2)*(1.+Q2*Q2)*std::log(1.+(xx*xp/Q2));
  ris+=2.*Nf*Cf*(Jac2/(xp*Q2)*B1pp_2-Jac1/(xp*Q1)*B1pp_1);
  
  //B2pp 
  const long double B2pp_1_1=-xx*xp*Q1/(2.*u1*t1)*(1.+Q1*Q1)*std::log(Qt1/Q1);
  const long double B2pp_1_2=+xx*xp*std::pow(Q1-u1,3.)/(2.*u1*t1*(Q1-t1))*((Q1-t1)+Q1*u1)
			     *(2./3.+Q1/Qt1-19./3.*std::pow(Q1/Qt1,3.))
			     -xx*xp*std::pow(Q1-u1,2.)/(2.*u1*t1*std::pow(Q1-t1,2.))
			     *(3.*std::pow(Q1-t1,3.)*Q1+(Q1-t1)*Q1*(2.*u1*t1+xx*xx)
			     +std::pow(Q1-t1,2.)*(1+4.*xx*Q1-u1*(Q1+xx))-u1*u1*Q1*Q1+u1*t1*t1*xx)*(1.-2.*Q1*Q1/(Qt1*Qt1))
			     +xx*xp*(Q1-u1)/(2.*u1*t1*(Q1-t1))*(3.*Q1*(Q1+1.)*(Q1-t1)-t1+xx*Q1+Q1*u1*(xx-Q1)*(xx-Q1))*(1.-2.*Q1/Qt1);
  const long double B2pp_1_3=xx*xp/(12.*u1*t1)*(-2.+6.*xx*t1*(t1-xx)+2.*xx*xx*xx+8.*Q1*(1.-Q1)*(1.-Q1)
			     -2.*u1*t1*Q1+7.*xx*xp*Q1-2.*Q1*Q1*xx-xx*xx*xx*Q1+3.*xx*Q1*Q1*Q1-4.*u1*t1*xx*Q1)
			     +11./6.*xx*xp*Q1*Q1/(u1*t1)-xx*xp*xx*xx*Q1/(3.*t1*t1);
  const long double B2pp_2_1=-xx*xp*Q2/(2.*u2*t2)*(1.+Q2*Q2)*std::log(Qt2/Q2);
  const long double B2pp_2_2=+xx*xp*std::pow(Q2-u2,3.)/(2.*u2*t2*(Q2-t2))*((Q2-t2)+Q2*u2)
			     *(2./3.+Q2/Qt2-19./3.*std::pow(Q2/Qt2,3.))
			     -xx*xp*std::pow(Q2-u2,2.)/(2.*u2*t2*std::pow(Q2-t2,2.))
			     *(3.*std::pow(Q2-t2,3.)*Q2+(Q2-t2)*Q2*(2.*u2*t2+xx*xx)
			     +std::pow(Q2-t2,2.)*(1+4.*xx*Q2-u2*(Q2+xx))-u2*u2*Q2*Q2+u2*t2*t2*xx)*(1.-2.*Q2*Q2/(Qt2*Qt2))
			     +xx*xp*(Q2-u2)/(2.*u2*t2*(Q2-t2))*(3.*Q2*(Q2+1.)*(Q2-t2)-t2+xx*Q2+Q2*u2*(xx-Q2)*(xx-Q2))*(1.-2.*Q2/Qt2);
  const long double B2pp_2_3=xx*xp/(12.*u2*t2)*(-2.+6.*xx*t2*(t2-xx)+2.*xx*xx*xx+8.*Q2*(1.-Q2)*(1.-Q2)
			     -2.*u2*t2*Q2+7.*xx*xp*Q2-2.*Q2*Q2*xx-xx*xx*xx*Q2+3.*xx*Q2*Q2*Q2-4.*u2*t2*xx*Q2)
			     +11./6.*xx*xp*Q2*Q2/(u2*t2)-xx*xp*xx*xx*Q2/(3.*t2*t2);
  //res[0]+=2.*Nc*Nf*(Jac2/(xp*Q2)*(B2pp_2_1+B2pp_2_2+B2pp_2_3)-Jac1/(xp*Q1)*(B2pp_1_1+B2pp_1_2+B2pp_1_3));
  ris+=2.*Nc*Nf*(Jac2/(xp*Q2)*(B2pp_2_1)-Jac1/(xp*Q1)*(B2pp_1_1)); // Il resto è integrato analiticamente in B2pp_ANALITICS
  /*if (res[0]!=res[0]) {
    std::cout << "nan Yes \t z= " << z << " x= " << xx << std::endl;
  }*/
  
  ris*=(1.-zmin);
  return ris;
}

long double NLOPL::LO_PL(double xp){ //we have to multiply later for sigma0 normalization
  const long double rad=std::sqrt((1.-x)*(1.-x)-4.*x*xp);
  const long double t=0.5*(x-1.+rad);
  const long double u=0.5*(x-1.-rad);
  switch(_channel){
    case(1):{
      return (4.*_Nc*(std::pow(1.+(x-1.)*x,2)-2.*(x-1.)*(x-1.)*x*xp+x*x*xp*xp)/(xp*rad)*_as/(2.*M_PIl));
      break;
    }
    case(2):{
      return(_as/(2.*M_PIl)*_Cf*x/rad*((1.+u*u)/(-t)+(1.+t*t)/(-u)));
      break;
    }
    case(3):{
      return(_as/(2.*M_PIl)*2.*_Cf*_Cf*x/rad*(u*u+t*t));
      break;
    }
  }
}




















//VECCHIA VERSIONE (Backup)

/*int _core_G2_sing(const int *ndim, const double x[], const int *ncomp, double res[], void *pars){
  res[0] = 0.;
  long double *p=(long double *) pars;
  
  //Set Energy, Transverse Momentum, Integration variable za
  const long double xx=p[0];
  const long double xp=p[1];
  const long double Nc=p[2];
  const long double Nf=p[3];
  const long double muF=p[4];
  const long double muR=p[5];
  const long double b0=11./6.*Nc-1./3.*Nf;
  const long double zmin=xx*std::pow(std::sqrt(1.+xp)+std::sqrt(xp),2);
  const long double z=zmin+(1.-zmin)*x[0];
  
  //Set Mandelstam variable
  const long double  rad=std::sqrt((z-xx)*(z-xx)-4.*xx*xp*z);
  const long double  t1=0.5*(xx-z+rad);
  const long double  u1=xx-0.5/z*(xx+z+rad);
  const long double  Q1=1.+t1+u1-xx;
  const long double  Qt1=Q1+xx*xp;
  const long double  zb1=(2.*xx*z-z-xx-rad)/(z*(-2.+xx+z-rad));
  const long double  Jac1=(xx-z+rad)/(2.*z*rad);
  const long double  Jac1z1=(xx-1.+std::sqrt((1.-xx)*(1.-xx)-4.*xx*xp))/(2.*std::sqrt((1.-xx)*(1.-xx)-4.*xx*xp));
  const long double  t2=0.5*(xx-z-rad);
  const long double  u2=xx-0.5/z*(xx+z-rad);
  const long double  Q2=1.+t2+u2-xx;
  const long double  Qt2=Q2+xx*xp;
  const long double  zb2=(2.*xx*z-xx-z+rad)/(z*(-2.+xx+z+rad));
  const long double  Jac2=-(xx-z-rad)/(2.*z*rad);
  const long double  Jac2z1=-(xx-1.-std::sqrt((1.-xx)*(1.-xx)-4.*xx*xp))/(2.*std::sqrt((1.-xx)*(1.-xx)-4.*xx*xp));
  
  //Define and add the different singular parts 
  const long double Si12=2.*xx/(-t2)*(1.+std::pow(z,4)+std::pow(1.-z,4))/(z)*Nc*Nc*(std::pow(xx,4)+std::pow(z,4)+std::pow(t2,4)+std::pow(xx*xp*z/t2,4))/(z*z*xx*xp);
  const long double Si11=2.*xx/(-t1)*(1.+std::pow(z,4)+std::pow(1.-z,4))/(z)*Nc*Nc*(std::pow(xx,4)+std::pow(z,4)+std::pow(t1,4)+std::pow(xx*xp*z/t1,4))/(z*z*xx*xp);
  const long double Si12z1=16.*Nc*Nc*(1.+std::pow(xx,4)-2.*xx*(1.+xp)-2*xx*xx*xx*(1.+xp)+xx*xx*(3.+4.*xp+xp*xp))/(xp*(1.-xx+sqrt(1.-2.*xx+xx*xx-4.*xx*xp)));
  const long double Si11z1=16.*Nc*Nc*(1.+std::pow(xx,4)-2.*xx*(1.+xp)-2*xx*xx*xx*(1.+xp)+xx*xx*(3.+4.*xp+xp*xp))/(xp*(1.-xx-sqrt(1.-2.*xx+xx*xx-4.*xx*xp)));
  res[0]+=std::log(1.-z)/(1.-z)*(Jac2*Si12-Jac2z1*Si12z1-Jac1*Si11+Jac1z1*Si11z1);
  
  const long double Si21=2.*xx*(z/(-t1))*Nc*Nc/2.*((std::pow(xx,4)+1.+std::pow(Q1,4)+std::pow(u1,4)+std::pow(t1,4))+z*zb1*(std::pow(xx,4)
		  +1.+std::pow(Q1,4)+std::pow(u1/zb1,4)+std::pow(t1/z,4)))/(u1*t1);
  const long double Si22=2.*xx*(z/(-t2))*Nc*Nc/2.*((std::pow(xx,4)+1.+std::pow(Q2,4)+std::pow(u2,4)+std::pow(t2,4))+z*zb2*(std::pow(xx,4)
		  +1.+std::pow(Q2,4)+std::pow(u2/zb2,4)+std::pow(t2/z,4)))/(u2*t2);
  const long double Si21z1=(8.*Nc*Nc*(-std::pow(1.+(xx-1.)*xx,2)+2.*std::pow(1.-xx,2)*xx*xp-xx*xx*xp*xp)/(xp*(-1.+xx+sqrt(std::pow(1.-xx,2)-4.*xx*xp))));
  const long double Si22z1=(8.*Nc*Nc*(-std::pow(1.+(xx-1.)*xx,2)+2.*std::pow(1.-xx,2)*xx*xp-xx*xx*xp*xp)/(xp*(-1.+xx-sqrt(std::pow(1.-xx,2)-4.*xx*xp))));
  res[0]+=std::log(1.-z)/(1.-z)*(Jac2*Si22-Jac2z1*Si22z1-Jac1*Si21+Jac1z1*Si21z1);
  
  const long double Si31=2.*xx*(z/(-t1))*Nc*Nc/2.*((std::pow(xx,4)+1.+std::pow(Q1,4)+std::pow(u1,4)+std::pow(t1,4))+z*zb1*(std::pow(xx,4)
		  +1.+std::pow(Q1,4)+std::pow(u1/zb1,4)+std::pow(t1/z,4)))/(u1*t1)*(std::log(Qt1*z/(-t1)));
  const long double Si32=2.*xx*(z/(-t2))*Nc*Nc/2.*((std::pow(xx,4)+1.+std::pow(Q2,4)+std::pow(u2,4)+std::pow(t2,4))+z*zb2*(std::pow(xx,4)
		  +1.+std::pow(Q2,4)+std::pow(u2/zb2,4)+std::pow(t2/z,4)))/(u2*t2)*(std::log(Qt2*z/(-t2)));
  const long double Si31z1=-(8.*Nc*Nc*(std::pow(1.+(xx-1.)*xx,2)-2.*std::pow(1.-xx,2)*xx*xp+xx*xx*xp*xp)/(xp*(-1.+xx+sqrt(std::pow(1.-xx,2)-4.*xx*xp)))
		    *std::log(2.*xx*xp/(1.-xx-sqrt(std::pow(1.-xx,2)-4.*xx*xp))));
  const long double Si32z1=-(8.*Nc*Nc*(std::pow(1.+(xx-1.)*xx,2)-2.*std::pow(1.-xx,2)*xx*xp+xx*xx*xp*xp)/(xp*(-1.+xx-sqrt(std::pow(1.-xx,2)-4.*xx*xp)))
		    *std::log(2.*xx*xp/(1.-xx+sqrt(std::pow(1.-xx,2)-4.*xx*xp))));
  res[0]+=-1./(1.-z)*(Jac2*Si32-Jac2z1*Si32z1-Jac1*Si31+Jac1z1*Si31z1);
  
  const long double Si41=2.*xx*(z/t1)*b0*Nc/2.*(std::pow(xx,4)+1.+z*zb1*(std::pow(u1/zb1,4)+std::pow(t1/z,4)))/(u1*t1);
  const long double Si42=2.*xx*(z/t2)*b0*Nc/2.*(std::pow(xx,4)+1.+z*zb2*(std::pow(u2/zb2,4)+std::pow(t2/z,4)))/(u2*t2);
  const long double Si41z1=(4.*Nc*b0*(std::pow(1.+(xx-1.)*xx,2)-2.*std::pow(1.-xx,2)*xx*xp+xx*xx*xp*xp)/(xp*(-1.+xx+sqrt(std::pow(1.-xx,2)-4.*xx*xp))));
  const long double Si42z1=(4.*Nc*b0*(std::pow(1.+(xx-1.)*xx,2)-2.*std::pow(1.-xx,2)*xx*xp+xx*xx*xp*xp)/(xp*(-1.+xx-sqrt(std::pow(1.-xx,2)-4.*xx*xp)))); 
  res[0]+=1./(1.-z)*(Jac2*Si42-Jac2z1*Si42z1-Jac1*Si41+Jac1z1*Si41z1);
  
  const long double Si51=2.*xx/t1*(1.+std::pow(z,4)+std::pow(1.-z,4))/z*std::log(muF*xx*z/(-t1))*Nc*Nc*(std::pow(xx,4)+std::pow(z,4)+std::pow(t1,4)+std::pow(xx*xp*z/t1,4))/(z*z*xx*xp);
  const long double Si52=2.*xx/t2*(1.+std::pow(z,4)+std::pow(1.-z,4))/z*std::log(muF*xx*z/(-t2))*Nc*Nc*(std::pow(xx,4)+std::pow(z,4)+std::pow(t2,4)+std::pow(xx*xp*z/t2,4))/(z*z*xx*xp);
  const long double Si51z1=(16.*Nc*Nc*(std::pow(1.+(xx-1.)*xx,2)-2.*std::pow(1.-xx,2)*xx*xp+xx*xx*xp*xp)/(xp*(-1.+xx+sqrt(std::pow(1.-xx,2)-4.*xx*xp)))*std::log(muF*2.*xx/(1.-xx-sqrt(std::pow(1.-xx,2)-4.*xx*xp))));
  const long double Si52z1=(16.*Nc*Nc*(std::pow(1.+(xx-1.)*xx,2)-2.*std::pow(1.-xx,2)*xx*xp+xx*xx*xp*xp)/(xp*(-1.+xx-sqrt(std::pow(1.-xx,2)-4.*xx*xp)))*std::log(muF*2.*xx/(1.-xx+sqrt(std::pow(1.-xx,2)-4.*xx*xp))));
  res[0]+=1./(1.-z)*(Jac2*Si52-Jac2z1*Si52z1-Jac1*Si51+Jac1z1*Si51z1);
  res[0]*=(1.-zmin);
  return 0;
}

long double NLOPL::G2_sing(double xp){
  double the_integral[1], error[1], prob[1];
  double epsrel=1.e-6, epsabs=1.e-15;
  int last = 4;
  int verbose = 0;
  int nregions, neval, fail;
  long double par[6]={x,xp,_Nc,_Nf,_muF,_muR};
  Cuhre(2, 1, _core_G2_sing, &par, 1,
	epsrel, epsabs, verbose | last,
	0, 500000, 9,
	NULL, NULL,
	&nregions, &neval, &fail, the_integral, error, prob);
  std::cout << "G2_sing("<< par[0] <<"," << par[1] << ") = " << the_integral[0] << " ± " << error[0] << "\tp = " << prob[0] << std::endl;
  std::cout << std::endl;
  return static_cast<long double>(the_integral[0]);
}

int _core_G2_notsing(const int *ndim, const double x[], const int *ncomp, double res[], void *pars){
  res[0] = 0.;
  long double *p=(long double *) pars;
  
  //Set Energy, Transverse Momentum, Integration variable y
  const long double xx=p[0];
  const long double xp=p[1];
  const long double Nc=p[2];
  const long double Nf=p[3];
  const long double muF=p[4];
  const long double muR=p[5];
  const long double b0=11./6.*Nc-1./3.*Nf;
  const long double Cf=(Nc*Nc-1.)/(2.*Nc);
  const long double zmin=xx*std::pow(std::sqrt(1.+xp)+std::sqrt(xp),2);
  const long double z=zmin+(1.-zmin)*x[0];
  
  //Set Mandelstam variable
  const long double  rad=std::sqrt((z-xx)*(z-xx)-4.*xx*xp*z);
  const long double  t1=0.5*(xx-z+rad);
  const long double  u1=xx-0.5/z*(xx+z+rad);
  long double  Q1=1.+t1+u1-xx;
  const long double  Qt1=Q1+xx*xp;
  const long double  zb1=(2.*xx*z-z-xx-rad)/(z*(-2.+xx+z-rad));
  const long double  Jac1=(xx-z+rad)/(2.*z*rad);
  const long double  t2=0.5*(xx-z-rad);
  const long double  u2=xx-0.5/z*(xx+z-rad);
  long double  Q2=1.+t2+u2-xx;
  const long double  Qt2=Q2+xx*xp;
  const long double  zb2=(2.*xx*z-xx-z+rad)/(z*(-2.+xx+z+rad));
  const long double  Jac2=-(xx-z-rad)/(2.*z*rad);
  const long double  A1=1.+xx-Q1;
  const long double  B1=std::sqrt(A1*A1-4.*xx);
  const long double  L1a1=std::log(xx/(z*z));
  const long double  L1b1=std::log(xx/(zb1*zb1));
  const long double  L2a1=std::log(xx/std::pow(A1-z,2));
  const long double  L2b1=std::log(xx/std::pow(A1-zb1,2));
  const long double  L31=std::log((A1+B1)/(A1-B1));
  const long double  A2=1.+xx-Q2;
  const long double  B2=std::sqrt(A2*A2-4.*xx);
  const long double  L1a2=std::log(xx/(z*z));
  const long double  L1b2=std::log(xx/(zb2*zb2));
  const long double  L2a2=std::log(xx/std::pow(A2-z,2));
  const long double  L2b2=std::log(xx/std::pow(A2-zb2,2));
  const long double  L32=std::log((A2+B2)/(A2-B2));
  
  //Define and add regular parts in G2s (paper Glover for definition of G2s)
  const long double Fi1_1=xx/(-t1)*(-2.*Nf*std::log(muF*xx/Q1)*0.5*(z*z+(1.-z)*(1.-z))+2.*Nf*z*(1.-z))*Cf*(z*z*xx*xx*xp*xp/(t1*t1)+z*z)/(-t1);
  const long double Fi1_2=xx/(-t2)*(-2.*Nf*std::log(muF*xx/Q2)*0.5*(z*z+(1.-z)*(1.-z))+2.*Nf*z*(1.-z))*Cf*(z*z*xx*xx*xp*xp/(t2*t2)+z*z)/(-t2);			
  res[0]+=2.*(Jac2*Fi1_2-Jac1*Fi1_1);
  
  
  //Separate divergent part Q->0 (ref Grazzini)
  const long double tiny=1e-4;
  const long double Fi2_1 = (Q1 > tiny*(xx*xp) ) ? Nc*Nc*((std::pow(xx,4)+1.+std::pow(Q1,4)+std::pow(u1/zb1,4)+std::pow(t1/z,4))*(Q1+Qt1)/(Q1*Qt1)
			+(2.*xx*xx*(std::pow(xx-t1,4)+std::pow(xx-u1,4)+std::pow(u1,4)+std::pow(t1,4)))/(u1*t1*(xx-u1)*(xx-t1)))*1./(xp)
			*std::log(xx*xp/Qt1) : -Nc*Nc*((std::pow(xx,4)+1.+std::pow(u1,4)+std::pow(t1,4))/(xx*xp*xp));
  const long double Fi2_2 = ( Q2 > tiny*(xx*xp)) ? Nc*Nc*((std::pow(xx,4)+1.+std::pow(Q2,4)+std::pow(u2/zb2,4)+std::pow(t2/z,4))*(Q2+Qt2)/(Q2*Qt2)
			+(2.*xx*xx*(std::pow(xx-t2,4)+std::pow(xx-u2,4)+std::pow(u2,4)+std::pow(t2,4)))/(u2*t2*(xx-u2)*(xx-t2)))*1./(xp)
			*std::log(xx*xp/Qt2) : -Nc*Nc*((std::pow(xx,4)+1.+std::pow(u2,4)+std::pow(t2,4))/(xx*xp*xp));
  res[0]+=(Jac2*Fi2_2-Jac1*Fi2_1);
  
  //Define and add G2ns (paper Glover for definition of G2ns)
  //A1234 
  const long double A1234_1_1=-0.5*((pow(xx*xp/t1,4)+pow(xx*Q1/t1,4))*L1a1+pow(u1,4)*L1b1);
  const long double A1234_1_2=(pow(xx,2)*Q1*pow(u1,3))/(2.*t1*(u1-xx)*(t1-xx))*(L1a1+L1b1)
			     +(xx*xp*Q1*pow(u1,3))/(2.*A1*(t1-xx))*(L2b1-L1b1)
			     +(xx*xp*Q1*(pow(xx,4)+pow(u1-xx,4)))/(2.*A1*t1*(u1-xx))*(L2a1-L1a1);
  const long double A1234_1_3=pow(xx,2)*xp*Q1*(pow(u1,4)/(2.*pow(B1*t1,2))+pow(u1,2)/(2.*pow(B1,2)))
			     +pow(xx,4)*pow(xp*Q1,2)*(-6./pow(B1,4)-4./pow(t1,4)+8./pow(B1*t1,2));
  const long double A1234_1_4=L31*(xx*xp*pow(u1,3)*(u1+t1)/(B1*t1)+pow(xx,4)*pow(Q1*xp,2)*(3.*A1/pow(B1,5)-1./(A1*pow(B1,3)))
			     -pow(xx,2)*xp*Q1*(1./(B1*t1)*(t1*t1+t1*u1+4.*pow(u1,2)-2.*xx*Q1)
			     +A1/(2.*pow(B1,3))*(t1*t1+3.*t1*u1+3.*u1*u1-6*xx*Q1)+1./(2.*A1*B1)*(t1*t1+t1*u1+7.*u1*u1-2.*xx*Q1)));
  const long double A1234_2_1=-0.5*((pow(xx*xp/t2,4)+pow(xx*Q2/t2,4))*L1a2+pow(u2,4)*L1b2);
  const long double A1234_2_2=(pow(xx,2)*Q2*pow(u2,3))/(2.*t2*(u2-xx)*(t2-xx))*(L1a2+L1b2)
			     +(xx*xp*Q2*pow(u2,3))/(2.*A2*(t2-xx))*(L2b2-L1b2)
			     +(xx*xp*Q2*(pow(xx,4)+pow(u2-xx,4)))/(2.*A2*t2*(u2-xx))*(L2a2-L1a2);
  const long double A1234_2_3=pow(xx,2)*xp*Q2*(pow(u2,4)/(2.*pow(B2*t2,2))+pow(u2,2)/(2.*pow(B2,2)))
			     +pow(xx,4)*pow(xp*Q2,2)*(-6./pow(B2,4)-4./pow(t2,4)+8./pow(B2*t2,2));
  const long double A1234_2_4=L32*(xx*xp*pow(u2,3)*(u2+t2)/(B2*t2)+pow(xx,4)*pow(Q2*xp,2)*(3.*A2/pow(B2,5)-1./(A2*pow(B2,3)))
			     -pow(xx,2)*xp*Q2*(1./(B2*t2)*(t2*t2+t2*u2+4.*pow(u2,2)-2.*xx*Q2)
			     +A2/(2.*pow(B2,3))*(t2*t2+3.*t2*u2+3.*u2*u2-6*xx*Q2)+1./(2.*A2*B2)*(t2*t2+t2*u2+7.*u2*u2-2.*xx*Q2)));
  res[0]+=2.*Nc*Nc*(Jac2/(xp*Q2)*(A1234_2_1+A1234_2_2+A1234_2_3+A1234_2_4)-Jac1/(xp*Q1)*(A1234_1_1+A1234_1_2+A1234_1_3+A1234_1_4));
  
  //A3412
  const long double A3412_1_1=xx*xp*Q1*pow(A1,3)/(2.*t1*(u1-xx))*(L2a1+L1b1)
			      +xx*xp*(u1+t1)/(16.*u1*t1*B1)*(pow(A1,4)+6.*pow(A1*B1,2)+pow(B1,4))*L31;
  const long double A3412_1_2=(-xx*xp/(2.*u1*t1)*(pow(1.-Q1,4)+pow(xx,4)+2.*Q1*A1*pow(1.-Q1,2)-2.*Q1*pow(xx,3))
			      -pow(xx,4)*pow(Q1*xp,2)/pow(u1,4)+(2.*pow(xx,2)*xp*Q1*(A1*A1-xx))/(u1*u1))*L1b1;
  const long double A3412_1_3=xx*xp*pow(Q1-u1,3)/(8.*u1*t1*(Q1-t1))*((Q1-t1)+Q1*u1)*
			      (4./3.+2.*xx*xp/Qt1+4.*pow(xx*xp/Qt1,2)-44./3.*pow(xx*xp/Qt1,3))
			      +xx*xp*pow(Q1-t1,3)/(8.*u1*t1*(Q1-u1))*((Q1-u1)+Q1*t1)*
			      (4./3.+2.*xx*xp/Qt1+4.*pow(xx*xp/Qt1,2)-44./3.*pow(xx*xp/Qt1,3));
  const long double A3412_1_4=xx*xp*pow(Q1-u1,2)/(4.*u1*t1*(Q1-t1))*
			      (-3.*(t1-xx)*((Q1-t1)+Q1*u1)-Q1*(xx*(t1-xx)+Q1*(u1-xx)))*(1.+2.*xx*xp/Qt1-6.*pow(xx*xp/Qt1,2))
			      +xx*xp*pow(Q1-t1,2)/(4.*u1*t1*(Q1-u1))*
			      (-3.*(u1-xx)*((Q1-u1)+Q1*t1)+3.*Q1*(xx*(u1-xx)+Q1*(t1-xx))+4.*u1)*(1.+2.*xx*xp/Qt1-6.*pow(xx*xp/Qt1,2));
  const long double A3412_1_5=xx*xp*(Q1-t1)/(2.*u1*t1*(Q1-u1))*(3.*pow(u1-xx,2)*((Q1-u1)+Q1*t1)+8.*u1*t1+2.*u1-2.*Q1*u1*pow(u1-Q1,2)
			      -3.*xx*Q1*pow(t1-xx,2)-3.*Q1*(xx-Q1)*pow(t1,2)-Q1*u1*(4.*u1*t1-u1*xx-Q1*t1+2.*pow(t1,2)-4.*xx*xx)+3.*xx*pow(Q1,2)*(t1-xx)
			      +xx*Q1*u1*(t1-Q1))*(1.-2.*xx*xp/Qt1)+xx*xp*(Q1-u1)/(2.*u1*t1*(Q1-t1))*
			      (3.*pow(t1-xx,2)*((Q1-t1)+Q1*u1)+3.*(t1-xx)*Q1*(xx*(t1-xx)+Q1*(u1-xx))+Q1*u1*(xx*(t1-Q1)+Q1*(u1-xx)))*(1.-2.*xp*xx/Qt1);
  const long double A3412_1_6=-4.*pow(xx,4)*pow(Q1*xp,2)/pow(u1,4)+xp*Q1*pow(xx*B1,2)/(2.*u1*u1)
			      +pow(xx,3)*xp/6.*((1.+Q1)/(u1*t1)+Q1/(u1*u1)+Q1/(t1*t1))+(2.*xp*pow(xx,3)*Q1)/(u1*u1)+pow(xx,3)*xp/u1
			      -xx*xp/(12.*t1*u1)*(30.*pow(xx,3)+54.*pow(Q1,2)*xx+8.*pow(Q1,3))
			      +xx*xp/(12*u1*t1)*(11.+17.*pow(xx,4)+Q1*(61.*u1*u1*t1+17.*pow(u1,3)+73.*u1*t1*t1+29.*pow(t1,3))
			      +xx*(24.*u1*u1*t1+6.*pow(u1,3)+36.*u1*t1*t1+18.*pow(t1,3))+Q1*Q1*(-21.*u1*u1-33.*t1*t1-52.*u1*t1)
			      +xx*Q1*(-73.*u1*u1-109.*t1*t1-170.*u1*t1)+xx*xx*(-23.*u1*u1-35.*t1*t1-52.*u1*t1)
			      +xx*xx*Q1*(134.*t1+110.*u1)+4.*pow(Q1,4)+52.*xx*pow(Q1,3)+20.*xx*xx*Q1*Q1-22.*pow(xx,3)*Q1);
  const long double A3412_2_1=xx*xp*Q2*pow(A2,3)/(2.*t2*(u2-xx))*(L2a2+L1b2)
			      +xx*xp*(u2+t2)/(16.*u2*t2*B2)*(pow(A2,4)+6.*pow(A2*B2,2)+pow(B2,4))*L32;
  const long double A3412_2_2=(-xx*xp/(2.*u2*t2)*(pow(1.-Q2,4)+pow(xx,4)+2.*Q2*A2*pow(1.-Q2,2)-2.*Q2*pow(xx,3))
			      -pow(xx,4)*pow(Q2*xp,2)/pow(u2,4)+(2.*pow(xx,2)*xp*Q2*(A2*A2-xx))/(u2*u2))*L1b2;
  const long double A3412_2_3=xx*xp*pow(Q2-u2,3)/(8.*u2*t2*(Q2-t2))*((Q2-t2)+Q2*u2)*
			      (4./3.+2.*xx*xp/Qt2+4.*pow(xx*xp/Qt2,2)-44./3.*pow(xx*xp/Qt2,3))
			      +xx*xp*pow(Q2-t2,3)/(8.*u2*t2*(Q2-u2))*((Q2-u2)+Q2*t2)*
			      (4./3.+2.*xx*xp/Qt2+4.*pow(xx*xp/Qt2,2)-44./3.*pow(xx*xp/Qt2,3));
  const long double A3412_2_4=xx*xp*pow(Q2-u2,2)/(4.*u2*t2*(Q2-t2))*
			      (-3.*(t2-xx)*((Q2-t2)+Q2*u2)-Q2*(xx*(t2-xx)+Q2*(u2-xx)))*(1.+2.*xx*xp/Qt2-6.*pow(xx*xp/Qt2,2))
			      +xx*xp*pow(Q2-t2,2)/(4.*u2*t2*(Q2-u2))*
			      (-3.*(u2-xx)*((Q2-u2)+Q2*t2)+3.*Q2*(xx*(u2-xx)+Q2*(t2-xx))+4.*u2)*(1.+2.*xx*xp/Qt2-6.*pow(xx*xp/Qt2,2));
  const long double A3412_2_5=xx*xp*(Q2-t2)/(2.*u2*t2*(Q2-u2))*(3.*pow(u2-xx,2)*((Q2-u2)+Q2*t2)+8.*u2*t2+2.*u2-2.*Q2*u2*pow(u2-Q2,2)
			      -3*xx*Q2*pow(t2-xx,2)-3.*Q2*(xx-Q2)*pow(t2,2)-Q2*u2*(4.*u2*t2-u2*xx-Q2*t2+2.*pow(t2,2)-4.*xx*xx)+3.*xx*pow(Q2,2)*(t2-xx)
			      +xx*Q2*u2*(t2-Q2))*(1-2.*xx*xp/Qt2)+xx*xp*(Q2-u2)/(2.*u2*t2*(Q2-t2))*
			      (3.*pow(t2-xx,2)*((Q2-t2)+Q2*u2)+3.*(t2-xx)*Q2*(xx*(t2-xx)+Q2*(u2-xx))+Q2*u2*(xx*(t2-Q2)+Q2*(u2-xx)))*(1-2.*xp*xx/Qt2);
  const long double A3412_2_6=-4.*pow(xx,4)*pow(Q2*xp,2)/pow(u2,4)+xp*Q2*pow(xx*B2,2)/(2.*u2*u2)
			      +pow(xx,3)*xp/6.*((1.+Q2)/(u2*t2)+Q2/(u2*u2)+Q2/(t2*t2))+(2.*xp*pow(xx,3)*Q2)/(u2*u2)+pow(xx,3)*xp/u2
			      -xx*xp/(12.*t2*u2)*(30.*pow(xx,3)+54.*pow(Q2,2)*xx+8.*pow(Q2,3))
			      +xx*xp/(12*u2*t2)*(11.+17.*pow(xx,4)+Q2*(61.*u2*u2*t2+17.*pow(u2,3)+73.*u2*t2*t2+29.*pow(t2,3))
			      +xx*(24.*u2*u2*t2+6.*pow(u2,3)+36.*u2*t2*t2+18.*pow(t2,3))+Q2*Q2*(-21.*u2*u2-33.*t2*t2-52.*u2*t2)
			      +xx*Q2*(-73.*u2*u2-109.*t2*t2-170.*u2*t2)+xx*xx*(-23.*u2*u2-35.*t2*t2-52.*u2*t2)
			      +xx*xx*Q2*(134.*t2+110.*u2)+4.*pow(Q2,4)+52.*xx*pow(Q2,3)+20.*xx*xx*Q2*Q2-22.*pow(xx,3)*Q2);
  res[0]+=2.*Nc*Nc*(Jac2/(xp*Q2)*(A3412_2_1+A3412_2_2+A3412_2_3+A3412_2_4+A3412_2_5+A3412_2_6)-Jac1/(xp*Q1)*(A3412_1_1+A3412_1_2+A3412_1_3+A3412_1_4+A3412_1_5+A3412_1_6));
  
  //A1324
  const long double A1324_1_1=-0.5*((pow(xx*xp/t1,4)+pow(xx*Q1/t1,4))*L1a1+pow(u1,4)*L1b1)
			      +pow(xx,3)*xp*Q1/(t1*t1)*L1a1+(pow(xx,2)*Q1*pow(u1,3)/(2.*t1*(u1-xx)*(t1-xx))+xx*xp*pow(u1,3)/(2.*t1))
			      *(L1a1+L1b1);
  const long double A1324_1_2=xx*xp*(1.-zb1)*pow(u1,3)/(2.*A1*(t1-xx))*(L2b1-L1b1)
			      +(xx*xp*(1.-z)*(pow(xx,4)+pow(u1-xx,4)))/(2.*A1*t1*(u1-xx))*(L2a1-L1a1)
			      +pow(xx,3)*xp*Q1/(A1*B1)*L31+xx*xx*xp*Q1/(2.*pow(t1,4))*(pow(xx*xp,2)-6.*xx*xx*xp*Q1+pow(xx*Q1,2));
  const long double A1324_2_1=-0.5*((pow(xx*xp/t2,4)+pow(xx*Q2/t2,4))*L1a2+pow(u2,4)*L1b2)
			      +pow(xx,3)*xp*Q2/(t2*t2)*L1a2+(pow(xx,2)*Q2*pow(u2,3)/(2.*t2*(u2-xx)*(t2-xx))+xx*xp*pow(u2,3)/(2.*t2))
			      *(L1a2+L1b2);
  const long double A1324_2_2=xx*xp*(1.-zb2)*pow(u2,3)/(2.*A2*(t2-xx))*(L2b2-L1b2)
			      +(xx*xp*(1.-z)*(pow(xx,4)+pow(u2-xx,4)))/(2.*A2*t2*(u2-xx))*(L2a2-L1a2)
			      +pow(xx,3)*xp*Q2/(A2*B2)*L32+xx*xx*xp*Q2/(2.*pow(t2,4))*(pow(xx*xp,2)-6.*xx*xx*xp*Q2+pow(xx*Q2,2));
  res[0]+=2.*Nc*Nc*(Jac2/(xp*Q2)*(A1324_2_1+A1324_2_2)-Jac1/(xp*Q1)*(A1324_1_1+A1324_1_2));
  
  //A3241
  const long double A3241_1_1=xx*xp*pow(A1,3)*(1.-z)/(2.*t1*(u1-xx))*(L2a1-L1a1)
			      +(-pow(xx,4)*pow(xp*Q1,2)/pow(t1,4)+pow(xx,3)*xp*pow(Q1,2)/(u1*t1)
			      -xx*xx*xp*Q1*pow(A1,4)/(2.*u1*t1*(u1-xx)*(t1-xx))+xx*xx*xp*Q1*(u1+t1)*(2.*A1*A1-xx)/(u1*t1*t1))*L1a1;
  const long double A3241_1_2=xx*xp*Q1*pow(Q1-u1,2)/(2.*u1*t1*pow(Q1-t1,2))*
			      (-u1*t1-pow(Q1-t1,2))*(-3.+10.*Q1/Qt1-6.*pow(Q1/Qt1,2))
			      +xx*xp*Q1*(Q1-u1)/(u1*t1*pow(Q1-t1,2))*(u1*t1*(Q1-u1)-pow(Q1-t1,3)-xx*pow(Q1-t1,2)-xx*(Q1-t1)*(Q1-u1))
			      *(-1.+2.*Q1/Qt1);
  const long double A3241_1_3=xx*xx*xp*Q1*(B1*B1/(2.*t1*t1)-2.*xx*Q1/(t1*t1)+(u1+t1)*(u1+t1)/(2.*u1*t1))
			      -4.*pow(xx,4)*pow(xp*Q1,2)/pow(t1,4)+xx*xp*Q1/(4.*u1*t1)
			      *(pow(t1+u1,2)-(t1+u1)*(6.*Q1+4.*xx)+6.*Q1*Q1+8.*xx*Q1)+pow(xx,3)*xp*Q1*pow(t1+u1,2)/(4.*u1*u1*t1*t1);
  const long double A3241_2_1=xx*xp*pow(A2,3)*(1.-z)/(2.*t2*(u2-xx))*(L2a2-L1a2)
			      +(-pow(xx,4)*pow(xp*Q2,2)/pow(t2,4)+pow(xx,3)*xp*pow(Q2,2)/(u2*t2)
			      -xx*xx*xp*Q2*pow(A2,4)/(2.*u2*t2*(u2-xx)*(t2-xx))+xx*xx*xp*Q2*(u2+t2)*(2.*A2*A2-xx)/(u2*t2*t2))*L1a2;
  const long double A3241_2_2=xx*xp*Q2*pow(Q2-u2,2)/(2.*u2*t2*pow(Q2-t2,2))*
			      (-u2*t2-pow(Q2-t2,2))*(-3.+10.*Q2/Qt2-6.*pow(Q2/Qt2,2))
			      +xx*xp*Q2*(Q2-u2)/(u2*t2*pow(Q2-t2,2))*(u2*t2*(Q2-u2)-pow(Q2-t2,3)-xx*pow(Q2-t2,2)-xx*(Q2-t2)*(Q2-u2))
			      *(-1.+2.*Q2/Qt2);
  const long double A3241_2_3=xx*xx*xp*Q2*(B2*B2/(2.*t2*t2)-2.*xx*Q2/(t2*t2)+(u2+t2)*(u2+t2)/(2.*u2*t2))
			      -4.*pow(xx,4)*pow(xp*Q2,2)/pow(t2,4)+xx*xp*Q2/(4.*u2*t2)
			      *(pow(t2+u2,2)-(t2+u2)*(6.*Q2+4.*xx)+6.*Q2*Q2+8.*xx*Q2)+pow(xx,3)*xp*Q2*pow(t2+u2,2)/(4.*u2*u2*t2*t2);
  res[0]+=2.*Nc*Nc*(Jac2/(xp*Q2)*(A3241_2_1+A3241_2_2+A3241_2_3)-Jac1/(xp*Q1)*(A3241_1_1+A3241_1_2+A3241_1_3));
  
  //Aepsilon
  const long double Aepsilon_1=4.*pow(xx,4)*pow(xp*Q1,2)*(1./pow(t1,4)+1./pow(u1,4));
  const long double Aepsilon_2=4.*pow(xx,4)*pow(xp*Q2,2)*(1./pow(t2,4)+1./pow(u2,4));
  res[0]+=Nc*Nc*(Jac2/(xp*Q2)*(Aepsilon_2)-Jac1/(xp*Q1)*(Aepsilon_1));
  
  //A0
  const long double A0_1= (pow(t1/z,4)+pow(u1/zb1,4))*xx*xp/(Qt1*Qt1)*(5.-7.*Q1/Qt1+20./3.*pow(Q1/Qt1,2))
			 +xx*xp*(17./3.+4.*std::log(xx*xp/Qt1));
  const long double A0_2=(pow(t2/z,4)+pow(u2/zb2,4))*xx*xp/(Qt2*Qt2)*(5.-7.*Q2/Qt2+20./3.*pow(Q2/Qt2,2))
			 +xx*xp*(17./3.+4.*std::log(xx*xp/Qt2));
  res[0]+=Nc*Nc*(Jac2/(xp)*(A0_2)-Jac1/(xp)*(A0_1));
  
  //B1pm
  const long double B1pm_1=xx*xp*z*pow(1.-z,3)/t1+pow(xx*xp*z/t1,3)*(1.-z)
			   +4.*pow(xx*xp*z/t1*(1.-z),2)-xx*xp*Q1*(1.+std::log(xx*xp/Qt1));
  const long double B1pm_2=xx*xp*z*pow(1.-z,3)/t2+pow(xx*xp*z/t2,3)*(1.-z)
			   +4.*pow(xx*xp*z/t2*(1.-z),2)-xx*xp*Q2*(1.+std::log(xx*xp/Qt2));
  res[0]+=2.*Nf*Cf*(Jac2/(xp*Q2)*B1pm_2-Jac1/(xp*Q1)*B1pm_1);
  
  //B2pm //FIXME? C'è una differenza con risultato mathematica ma non capisco da dove dipenda. Contributo completamente negiglible comunque
  const long double B2pm_1=1./3.*pow(t1/z,4.)*xx*xp/Qt1*(-3./Qt1+3.*Q1/(Qt1*Qt1)-2.*Q1*Q1/(Qt1*Qt1*Qt1))
			   -1./3.*xx*xp;
  const long double B2pm_2=1./3.*pow(t2/z,4.)*xx*xp/Qt2*(-3./Qt2+3.*Q1/(Qt2*Qt2)-2.*Q2*Q2/(Qt2*Qt2*Qt2))
			   -1./3.*xx*xp;
  res[0]+=2.*Nf*Nc*(Jac2/xp*B2pm_2-Jac1/xp*B2pm_1);
  
  //B1pp
  const long double B1pp_1=xx*xp*pow(z,3)*(1.-z)/t1+pow(xx*xp*(1.-z)/t1,3)*z
			   +4.*pow(xx*xp*z*(1.-z)/t1,2)-xx*xp*Q1/(Qt1*Qt1*u1*t1)*(pow(u1*t1+xx*xp*Q1,2)+2.*xx*xp*Q1*Qt1)
			   +xx*xp*Q1/(u1*t1)*(1.+Q1*Q1)*std::log(1.+(xx*xp/Q1));
  const long double B1pp_2=xx*xp*pow(z,3)*(1.-z)/t2+pow(xx*xp*(1.-z)/t2,3)*z
			   +4.*pow(xx*xp*z*(1.-z)/t2,2)-xx*xp*Q2/(Qt2*Qt2*u2*t2)*(pow(u2*t2+xx*xp*Q2,2)+2.*xx*xp*Q2*Qt2)
			   +xx*xp*Q2/(u2*t2)*(1.+Q2*Q2)*std::log(1.+(xx*xp/Q2));
  res[0]+=2.*Nf*Cf*(Jac2/(xp*Q2)*B1pp_2-Jac1/(xp*Q1)*B1pp_1);
  
  //B2pp //FIXME
  const long double B2pp_1_1=-xx*xp*Q1/(2.*u1*t1)*(1.+Q1*Q1)*std::log(Qt1/Q1);
  const long double B2pp_1_2=+xx*xp*std::pow(Q1-u1,3.)/(2.*u1*t1*(Q1-t1))*((Q1-t1)+Q1*u1)
			     *(2./3.+Q1/Qt1-19./3.*std::pow(Q1/Qt1,3.))
			     -xx*xp*std::pow(Q1-u1,2.)/(2.*u1*t1*std::pow(Q1-t1,2.))
			     *(3.*std::pow(Q1-t1,3.)*Q1+(Q1-t1)*Q1*(2.*u1*t1+xx*xx)
			     +std::pow(Q1-t1,2.)*(1+4.*xx*Q1-u1*(Q1+xx))-u1*u1*Q1*Q1+u1*t1*t1*xx)*(1.-2.*Q1*Q1/(Qt1*Qt1))
			     +xx*xp*(Q1-u1)/(2.*u1*t1*(Q1-t1))*(3.*Q1*(Q1+1.)*(Q1-t1)-t1+xx*Q1+Q1*u1*(xx-Q1)*(xx-Q1))*(1.-2.*Q1/Qt1);
  const long double B2pp_1_3=xx*xp/(12.*u1*t1)*(-2.+6.*xx*t1*(t1-xx)+2.*xx*xx*xx+8.*Q1*(1.-Q1)*(1.-Q1)
			     -2.*u1*t1*Q1+7.*xx*xp*Q1-2.*Q1*Q1*xx-xx*xx*xx*Q1+3.*xx*Q1*Q1*Q1-4.*u1*t1*xx*Q1)
			     +11./6.*xx*xp*Q1*Q1/(u1*t1)-xx*xp*xx*xx*Q1/(3.*t1*t1);
  const long double B2pp_2_1=-xx*xp*Q2/(2.*u2*t2)*(1.+Q2*Q2)*std::log(Qt2/Q2);
  const long double B2pp_2_2=+xx*xp*std::pow(Q2-u2,3.)/(2.*u2*t2*(Q2-t2))*((Q2-t2)+Q2*u2)
			     *(2./3.+Q2/Qt2-19./3.*std::pow(Q2/Qt2,3.))
			     -xx*xp*std::pow(Q2-u2,2.)/(2.*u2*t2*std::pow(Q2-t2,2.))
			     *(3.*std::pow(Q2-t2,3.)*Q2+(Q2-t2)*Q2*(2.*u2*t2+xx*xx)
			     +std::pow(Q2-t2,2.)*(1+4.*xx*Q2-u2*(Q2+xx))-u2*u2*Q2*Q2+u2*t2*t2*xx)*(1.-2.*Q2*Q2/(Qt2*Qt2))
			     +xx*xp*(Q2-u2)/(2.*u2*t2*(Q2-t2))*(3.*Q2*(Q2+1.)*(Q2-t2)-t2+xx*Q2+Q2*u2*(xx-Q2)*(xx-Q2))*(1.-2.*Q2/Qt2);
  const long double B2pp_2_3=xx*xp/(12.*u2*t2)*(-2.+6.*xx*t2*(t2-xx)+2.*xx*xx*xx+8.*Q2*(1.-Q2)*(1.-Q2)
			     -2.*u2*t2*Q2+7.*xx*xp*Q2-2.*Q2*Q2*xx-xx*xx*xx*Q2+3.*xx*Q2*Q2*Q2-4.*u2*t2*xx*Q2)
			     +11./6.*xx*xp*Q2*Q2/(u2*t2)-xx*xp*xx*xx*Q2/(3.*t2*t2);
  //res[0]+=2.*Nc*Nf*(Jac2/(xp*Q2)*(B2pp_2_1+B2pp_2_2+B2pp_2_3)-Jac1/(xp*Q1)*(B2pp_1_1+B2pp_1_2+B2pp_1_3));
  res[0]+=2.*Nc*Nf*(Jac2/(xp*Q2)*(B2pp_2_1)-Jac1/(xp*Q1)*(B2pp_1_1));
  /*if (res[0]!=res[0]) {
    std::cout << "nan Yes \t z= " << z << " x= " << xx << std::endl;
  }*//*
  
  res[0]*=(1.-zmin);
  return 0;
}

long double NLOPL::G2_notsing(double xp){
  double the_integral[1], error[1], prob[1];
  double epsrel=1.e-6, epsabs=1.e-15;
  int last = 4;
  int verbose = 0;
  int nregions, neval, fail;
  long double par[6]={x,xp,_Nc,_Nf,_muF,_muR};
  Cuhre(2, 1, _core_G2_notsing, &par, 1,
	epsrel, epsabs, verbose | last,
	0, 500000, 9,
	NULL, NULL,
	&nregions, &neval, &fail, the_integral, error, prob);
  std::cout << "G2_notsing("<< par[0] <<"," << par[1] << ") = " << the_integral[0] << " ± " << error[0] << "\tp = " << prob[0] << std::endl;
  std::cout << std::endl;
  return (static_cast<long double>(the_integral[0]));
}
		
long double NLOPL::NLO_PL(double xp){
  long double sigma0=_as*_as*std::sqrt(2)*RHEHpt::Gf/(576.*M_PIl);
  return (sigma0*(_as*_as/(4.*M_PIl*M_PIl)*(G2_sing(xp)+G2_notsing(xp)+G2_delta(xp)))); 
}

long double NLOPL::LO_PL(double xp){
  long double sigma0=_as*_as*std::sqrt(2)*RHEHpt::Gf/(576.*M_PIl);
  const long double rad=std::sqrt((1.-x)*(1.-x)-4.*x*xp);
  const long double t=0.5*(x-1.+rad);
  const long double u=0.5*(x-1.-rad);
  const long double ris=4.*_Nc*(std::pow(1.+(x-1.)*x,2)-2.*(x-1.)*(x-1.)*x*xp+x*x*xp*xp)/(xp*rad);
  return (sigma0*ris*_as/(2.*M_PIl));
  
}*/
  
		
		
//I check all the terms with other programs all return the same value for the integral, in many times with better convergence.		
//There are still problem with B2pp, I'm working on it


		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		