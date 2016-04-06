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
  const long double MUF=pow(_muF/_mH,2.);
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
  const long double Si5z1=1./t*(std::pow(x,4)+1.+std::pow(t,4)+std::pow(u,4))/(u*t)*std::log(x*MUF/(-t));
  const long double Si5z2=1./u*(std::pow(x,4)+1.+std::pow(u,4)+std::pow(t,4))/(u*t)*std::log(x*MUF/(-u));

  return (ris/rad*2.+2.*x*_Nc*b0*(Jac2*Si5z2-Jac1*Si5z1)+_Nc*_Nf*B2pp_ANALITICS(x,xp));
}

//gq channel

long double NLOPL::NLO_PL_delta_gq(double xp){
  const long double rad=std::sqrt((1.-x)*(1.-x)-4.*x*xp);
  const long double t=0.5*(x-1.+rad);
  const long double u=0.5*(x-1.-rad);
  const long double MUR=pow(_muR/_mH,2.);
  const long double MUF=pow(_muF/_mH,2.);
  const long double b0=11./6.*_Nc-1./3.*_Nf;
  const long double Delta=(5.*_Nc-3.*(_Nc*_Nc-1.)/(2.*_Nc));
  const long double V11=0.5*std::pow(std::log(u/t),2.)+0.5*std::pow(std::log(1./(-u)),2.)-0.5*std::pow(std::log(1./(-t)),2.)
		      -std::log(x)*std::log((-t)/x)+std::log(x)*std::log((-u)/x)-std::log((-t)/x)*std::log((-u)/x)
		      +2.*dilog_r(x/(x-u))+std::pow(std::log(x/(x-u)),2.)+M_PIl*M_PIl;
  const long double V21=std::pow(std::log(x),2)+std::pow(std::log(x/(x-t)),2)-2.*std::log(1./x)*std::log((-t)/x)+2.*dilog_r(1.-x)
		      +2.*dilog_r(x/(x-t))-7./2.-2.*M_PIl*M_PIl/3.;
  const long double V31=b0*(2.*std::log(MUR*x/(-u))+std::log(MUR*x/(-t)))+(67./9.*_Nc-10./9.*_Nf);
  const long double V12=0.5*std::pow(std::log(t/u),2.)+0.5*std::pow(std::log(1./(-t)),2.)-0.5*std::pow(std::log(1./(-u)),2.)
		      -std::log(x)*std::log((-u)/x)+std::log(x)*std::log((-t)/x)-std::log((-u)/x)*std::log((-t)/x)
		      +2.*dilog_r(x/(x-t))+std::pow(std::log(x/(x-t)),2.)+M_PIl*M_PIl;
  const long double V22=std::pow(std::log(x),2)+std::pow(std::log(x/(x-u)),2)-2.*std::log(1./x)*std::log((-u)/x)+2.*dilog_r(1.-x)
		      +2.*dilog_r(x/(x-u))-7./2.-2.*M_PIl*M_PIl/3.;
  const long double V32=b0*(2.*std::log(MUR*x/(-t))+std::log(MUR*x/(-u)))+(67./9.*_Nc-10./9.*_Nf);
  long double ris=0.;
  ris+=x*((Delta+_Nc*V11+_Cf*V21+V31)*_Cf*(1.+t*t)/(-u)+(_Nc-_Cf)*_Cf*((1.+t*t+u*u-u*x)/(-u)));
  ris+=x*((Delta+_Nc*V12+_Cf*V22+V32)*_Cf*(1.+u*u)/(-t)+(_Nc-_Cf)*_Cf*((1.+u*u+t*t-t*x)/(-t)));
  const long double Jac1=-(1.-x-rad)/(2.*rad);
  const long double Jac2=(1.-x+rad)/(2.*rad);
  const long double Sideltaz1=x/t*b0*std::log(MUF*x/(-t))*(1.+t*t)/(-x*xp/(t));
  const long double Sideltaz2=x/u*b0*std::log(MUF*x/(-u))*(1.+u*u)/(-x*xp/(u));
  const long double Sideltazb1=x/(t)*3./2.*_Cf*std::log(MUF*x/(-t))*_Cf*(1.+u*u)/(-t);
  const long double Sideltazb2=x/(u)*3./2.*_Cf*std::log(MUF*x/(-u))*_Cf*(1.+t*t)/(-u);
  return(ris/rad+Jac2*(Sideltaz2+Sideltazb2)-Jac1*(Sideltaz1+Sideltazb1)); // OK
}

//qqbar channel
long double NLOPL::NLO_PL_delta_qqbar(double xp){
  const long double rad=std::sqrt((1.-x)*(1.-x)-4.*x*xp);
  const long double t=0.5*(x-1.+rad);
  const long double u=0.5*(x-1.-rad);
  const long double MUR=pow(_muR/_mH,2.);
  const long double MUF=pow(_muF/_mH,2.);
  const long double b0=11./6.*_Nc-1./3.*_Nf;
  const long double Delta=(5.*_Nc-3.*(_Nc*_Nc-1.)/(2.*_Nc));
  const long double W1=std::log(-u/x)*std::log(-t/x)-std::log(1./x)*std::log(-u/x)-std::log(1./x)*std::log(-t/x)
		      +2.*dilog_r(1.-x)+std::pow(std::log(x),2)-0.5*std::pow(std::log(u/t),2)-5./3.*M_PIl*M_PIl;
  const long double W2=3./2.*(std::log(1./(-t))+std::log(1./(-u)))+std::pow(std::log(u/t),2)-2.*std::log(-u/x)*std::log(-t/x)
		      +std::pow(std::log(x/(x-u)),2)+std::pow(std::log(x/(x-t)),2)+2.*dilog_r(x/(x-u))+2.*dilog_r(x/(x-t))-7.+2.*M_PIl*M_PIl;
  const long double W3=b0/2.*(4.*std::log(MUR*x)+std::log(MUR*x/(-u))+std::log(MUR*x/(-t)))+(67./6.*_Nc-5./3.*_Nf);
  const long double ris=x*((Delta+_Nc*W1+_Cf*W2+W3)*2.*_Cf*_Cf*(t*t+u*u)+(_Nc-_Cf)*2.*_Cf*_Cf*((t*t+u*u+1.-x)));
  const long double Jac1=-(1.-x-rad)/(2.*rad);
  const long double Jac2=(1.-x+rad)/(2.*rad);
  const long double Sidelta1=2.*x/(t)*_Cf*3./2.*std::log(MUF*x/(-t))*2.*_Cf*_Cf*(x*x*xp*xp/(t*t)+t*t);
  const long double Sidelta2=2.*x/(u)*_Cf*3./2.*std::log(MUF*x/(-u))*2.*_Cf*_Cf*(x*x*xp*xp/(u*u)+u*u);
  
  //return (2.*ris/rad+Jac2*Sidelta2-Jac1*Sidelta1);
  return 0;

  
}


//Singular Part
//gg channel
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

//gq channel
long double NLOPL::NLO_PL_sing_doublediff_gq(double xp, double zz)
{
  long double ris=0.0;
  //Set Energy, Transverse Momentum, Integration variable za
  const long double xx=x;
  const long double Nc=_Nc;
  const long double Nf=_Nf;
  const long double MUF=std::pow(_muF/_mH,2);
  const long double MUR=std::pow(_muR/_mH,2);
  const long double b0=11./6.*Nc-1./3.*Nf;
  const long double zmin=xx*std::pow(std::sqrt(1.+xp)+std::sqrt(xp),2);
  const long double z=zmin+(1.-zmin)*zz;
  
  //Set Mandelstam variable
  const long double  rad=std::sqrt((z-xx)*(z-xx)-4.*xx*xp*z);
  const long double  t1=0.5*(xx-z+rad);
  const long double  u1=xx-0.5/z*(xx+z+rad);
  const long double rad1=std::sqrt((1.-xx)*(1.-xx)-4.*xx*xp);
  const long double t=0.5*(xx-1.+rad1);
  const long double u=0.5*(xx-1.-rad1);
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
  const long double Si11=xx/(t1)*_Nc*_Cf*((1.+std::pow(z,4)+std::pow(1.-z,4))/z*std::log(MUF*xx*z/(-t1)))*(z*z+t1*t1)/(-xx*xp*z/(t1));
  const long double Si12=xx/(t2)*_Nc*_Cf*((1.+std::pow(z,4)+std::pow(1.-z,4))/z*std::log(MUF*xx*z/(-t2)))*(z*z+t2*t2)/(-xx*xp*z/(t2));
  const long double Si11z1=xx/(t)*_Nc*_Cf*(2.*std::log(MUF*xx/(-t)))*(1.+t*t)/(-xx*xp/(t));
  const long double Si12z1=xx/(u)*_Nc*_Cf*(2.*std::log(MUF*xx/(-u)))*(1.+u*u)/(-xx*xp/(u));
  ris+=1./(1.-z)*(Jac2*Si12-Jac2z1*Si12z1-Jac1*Si11+Jac1z1*Si11z1); // OK 4 Nc Cf Log(x)^2/xp
  
  const long double Si21=xx/(-t1)*_Nc*_Cf*(1.+std::pow(z,4)+std::pow(1.-z,4))/z*(z*z+t1*t1)/(-xx*xp*z/(t1));
  const long double Si22=xx/(-t2)*_Nc*_Cf*(1.+std::pow(z,4)+std::pow(1.-z,4))/z*(z*z+t2*t2)/(-xx*xp*z/(t2));
  const long double Si21z1=xx/(-t)*_Nc*_Cf*(2.)*(1.+t*t)/(-xx*xp/(t));
  const long double Si22z1=xx/(-u)*_Nc*_Cf*(2.)*(1.+u*u)/(-xx*xp*z/(u));
  ris+=std::log(1.-z)/(1.-z)*(Jac2*Si22-Jac2z1*Si22z1-Jac1*Si21+Jac1z1*Si21z1); //OK no logs
  
  const long double Si31=xx/(t1)*(_Cf*(1.+z*z)*std::log(MUF*xx*z/(-t1)))*_Cf*(z*z+xx*xx*xp*xp*z*z/(t1*t1))/(-t1);
  const long double Si32=xx/(t2)*(_Cf*(1.+z*z)*std::log(MUF*xx*z/(-t2)))*_Cf*(z*z+xx*xx*xp*xp*z*z/(t2*t2))/(-t2);
  const long double Si31z1=xx/(t)*(_Cf*(2.)*std::log(MUF*xx/(-t)))*_Cf*(1.+xx*xx*xp*xp/(t*t))/(-t);
  const long double Si32z1=xx/(u)*(_Cf*(2.)*std::log(MUF*xx/(-u)))*_Cf*(1.+xx*xx*xp*xp/(u*u))/(-u);
  ris+=1./(1.-z)*(Jac2*Si32-Jac2z1*Si32z1-Jac1*Si31+Jac1z1*Si31z1); //OK no logs
  
  const long double Si41=xx/(-t1)*(_Cf*(1.+z*z))*_Cf*(z*z+xx*xx*xp*xp*z*z/(t1*t1))/(-t1);
  const long double Si42=xx/(-t2)*(_Cf*(1.+z*z))*_Cf*(z*z+xx*xx*xp*xp*z*z/(t2*t2))/(-t2);
  const long double Si41z1=xx/(-t)*(_Cf*(2.))*_Cf*(1.+xx*xx*xp*xp/(t*t))/(-t);
  const long double Si42z1=xx/(-u)*(_Cf*(2.))*_Cf*(1.+xx*xx*xp*xp/(u*u))/(-u);
  ris+=std::log(1.-z)/(1.-z)*(Jac2*Si42-Jac2z1*Si42z1-Jac1*Si41+Jac1z1*Si41z1); //OK no logs
  
  const long double Si51=xx*z/(-t1)*_Nc*_Cf*((-t1-t1*t1*t1+Q1*Q1*Q1*t1+Q1*t1*t1*t1)/(u1*t1)
			 +(z*zb1*(-(t1/z)-std::pow(t1/z,3)-Q1*Q1*Q1*(u1/zb1)-Q1*std::pow(u1/zb1,3)))/(u1*t1));
  const long double Si52=xx*z/(-t2)*_Nc*_Cf*((-t2-t2*t2*t2+Q2*Q2*Q2*t2+Q2*t2*t2*t2)/(u2*t2)
			 +(z*zb2*(-(t2/z)-std::pow(t2/z,3)-Q2*Q2*Q2*(u2/zb2)-Q2*std::pow(u2/zb2,3)))/(u2*t2));
  const long double Si51z1=xx/(-t)*_Nc*_Cf*((-t-t*t*t)/(u*t)+(-(t)-std::pow(t,3))/(u*t));
  const long double Si52z1=xx*z/(-u)*_Nc*_Cf*((-u-u*u*u)/(u*t)+(-(u)-std::pow(u,3))/(u*t));
  ris+=std::log(1.-z)/(1.-z)*(Jac2*Si52-Jac2z1*Si52z1-Jac1*Si51+Jac1z1*Si51z1); //OK +Nc*Cf*Log(x)^2/xp
  
  const long double Si61=xx*z/(t1)*std::log(Qt1*z/(-t1))*_Nc*_Cf*((-t1-t1*t1*t1+Q1*Q1*Q1*t1+Q1*t1*t1*t1)/(u1*t1)
			 +(z*zb1*(-(t1/z)-std::pow(t1/z,3)-Q1*Q1*Q1*(u1/zb1)-Q1*std::pow(u1/zb1,3)))/(u1*t1));
  const long double Si62=xx*z/(t2)*std::log(Qt2*z/(-t2))*_Nc*_Cf*((-t2-t2*t2*t2+Q2*Q2*Q2*t2+Q2*t2*t2*t2)/(u2*t2)
			 +(z*zb2*(-(t2/z)-std::pow(t2/z,3)-Q2*Q2*Q2*(u2/zb2)-Q2*std::pow(u2/zb2,3)))/(u2*t2));
  const long double Si61z1=xx/(t)*std::log(xx*xp/(-t))*_Nc*_Cf*((-t-t*t*t)/(u*t)+(-(t)-std::pow(t,3))/(u*t));
  const long double Si62z1=xx/(u)*std::log(xx*xp/(-u))*_Nc*_Cf*((-u-u*u*u)/(u*t)+(-(u)-std::pow(u,3))/(u*t));
  ris+=1./(1.-z)*(Jac2*Si62-Jac2z1*Si62z1-Jac1*Si61+Jac1z1*Si61z1); //OK -3Nc*Cf*Log(x)^2/xp
  
  const long double Si71=xx*z/(t1)*3./2.*_Cf*_Cf*(u1*u1+1.)/(-t1);
  const long double Si72=xx*z/(t2)*3./2.*_Cf*_Cf*(u2*u2+1.)/(-t2);
  const long double Si71z1=xx*z/(t)*3./2.*_Cf*_Cf*(u*u+1.)/(-t);
  const long double Si72z1=xx*z/(u)*3./2.*_Cf*_Cf*(t*t+1.)/(-u);
  ris+=1./(1.-z)*(Jac2*Si72-Jac2z1*Si72z1-Jac1*Si71+Jac1z1*Si71z1);//OK (Log semplice non Log squared)
  
  ris*=(1.-zmin);
  return ris;
}

//qqbar channel

long double NLOPL::NLO_PL_sing_doublediff_qqbar(double xp, double zz)
{
  long double ris=0.0;
  //Set Energy, Transverse Momentum, Integration variable za
  const long double xx=x;
  const long double Nc=_Nc;
  const long double Nf=_Nf;
  const long double MUF=std::pow(_muF/_mH,2);
  const long double MUR=std::pow(_muR/_mH,2);
  const long double b0=11./6.*Nc-1./3.*Nf;
  const long double zmin=xx*std::pow(std::sqrt(1.+xp)+std::sqrt(xp),2);
  const long double z=zmin+(1.-zmin)*zz;
  
  //Set Mandelstam variable
  const long double  rad=std::sqrt((z-xx)*(z-xx)-4.*xx*xp*z);
  const long double  t1=0.5*(xx-z+rad);
  const long double  u1=xx-0.5/z*(xx+z+rad);
  const long double rad1=std::sqrt((1.-xx)*(1.-xx)-4.*xx*xp);
  const long double t=0.5*(xx-1.+rad1);
  const long double u=0.5*(xx-1.-rad1);
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
  const long double Si11=xx/(-t1)*(-_Cf*(1.+z*z)*std::log(MUF*z*xx/(-t1)))*2.*_Cf*_Cf*(xx*xx*xp*xp/(t1*t1)*z*z+t1*t1)/z;
  const long double Si12=xx/(-t2)*(-_Cf*(1.+z*z)*std::log(MUF*z*xx/(-t2)))*2.*_Cf*_Cf*(xx*xx*xp*xp/(t2*t2)*z*z+t2*t2)/z;
  const long double Si11z1=xx/(-t)*(-_Cf*(2.)*std::log(MUF*xx/(-t)))*2.*_Cf*_Cf*(xx*xx*xp*xp/(t*t)+t*t);
  const long double Si12z1=xx/(-t2)*(-_Cf*(2.)*std::log(MUF*xx/(-u)))*2.*_Cf*_Cf*(xx*xx*xp*xp/(u*u)+u*u);
  ris+=2./(1.-z)*(Jac2*Si12-Jac2z1*Si12z1-Jac1*Si11+Jac1z1*Si11z1);
  
  const long double Si21= xx/(-t1)*(_Cf*(1.+z*z))*2.*_Cf*_Cf*(xx*xx*xp*xp/(t1*t1)*z*z+t1*t1)/z;
  const long double Si22=xx/(-t2)*(_Cf*(1.+z*z))*2.*_Cf*_Cf*(xx*xx*xp*xp/(t2*t2)*z*z+t2*t2)/z;
  const long double Si21z1=xx/(-t)*(_Cf*(2.))*2.*_Cf*_Cf*(xx*xx*xp*xp/(t*t)+t*t);
  const long double Si22z1=xx/(-u)*(_Cf*(2.))*2.*_Cf*_Cf*(xx*xx*xp*xp/(u*u)+u*u);
  ris+=2.*std::log(1.-z)/(1.-z)*(Jac2*Si22-Jac2z1*Si22z1-Jac1*Si21+Jac1z1*Si21z1);
  
  const long double Si31=xx*z/(-t1)*(2.*_Cf-_Nc)*_Cf*_Cf*((t1*t1+u1*u1+std::pow(t1/z,2)+std::pow(u1/zb1,2)));
  const long double Si32=xx*z/(-t2)*(2.*_Cf-_Nc)*_Cf*_Cf*((t2*t2+u2*u2+std::pow(t2/z,2)+std::pow(u2/zb1,2)));
  const long double Si31z1=xx/(-t)*(2.*_Cf-_Nc)*_Cf*_Cf*((t*t+u*u+std::pow(t,2)+std::pow(u,2)));
  const long double Si32z1=xx/(-u)*(2.*_Cf-_Nc)*_Cf*_Cf*((u*u+t*t+std::pow(u,2)+std::pow(t,2)));
  ris+=2.*std::log(1.-z)/(1.-z)*(Jac2*Si32-Jac2z1*Si32z1-Jac1*Si31+Jac1z1*Si31z1);
  
  const long double Si41=xx*z/(t1)*std::log(Qt1*z/(-t1))*(2.*_Cf-_Nc)*_Cf*_Cf*((t1*t1+u1*u1+std::pow(t1/z,2)+std::pow(u1/zb1,2)));
  const long double Si42=xx*z/(t2)*std::log(Qt2*z/(-t2))*(2.*_Cf-_Nc)*_Cf*_Cf*((t2*t2+u2*u2+std::pow(t2/z,2)+std::pow(u2/zb1,2)));
  const long double Si41z1=xx/t*std::log(xx*xp/(-t))*(2.*_Cf-_Nc)*_Cf*_Cf*((t*t+u*u+std::pow(t,2)+std::pow(u,2)));
  const long double Si42z1=xx/u*std::log(xx*xp/(-u))*(2.*_Cf-_Nc)*_Cf*_Cf*((u*u+t*t+std::pow(u,2)+std::pow(t,2)));
  ris+=2./(1.-z)*(Jac2*Si42-Jac2z1*Si42z1-Jac1*Si41+Jac1z1*Si41z1);
  
  const long double Si51=-xx*z/(-t1)*b0*_Cf*_Cf*((t1*t1+u1*u1));
  const long double Si52=-xx*z/(-t2)*b0*_Cf*_Cf*((t2*t2+u2*u2));
  const long double Si51z1=-xx/(-t)*b0*_Cf*_Cf*((t*t+u*u));
  const long double Si52z1=-xx/(-u)*b0*_Cf*_Cf*((t*t+u*u));
  ris+=2./(1.-z)*(Jac2*Si52-Jac2z1*Si52z1-Jac1*Si51+Jac1z1*Si51z1);
  
  ris*=(1.-zmin);
  return ris;
}


//Not-Singular Part
//gg channel

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

//gq channel

long double NLOPL::NLO_PL_notsing_doublediff_gq(double xp, double zz){
  
  long double ris=0.0;
  
  //Set Energy, Transverse Momentum, Integration variable z
  const long double xx=x;
  const long double Nc=_Nc;
  const long double Nf=_Nf;
  const long double MUF=std::pow(_muF/_mH,2);
  const long double MUR=std::pow(_muR/_mH,2);
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
  
  //Define and add regular parts in G2s (refer to Glover paper for definition of G2s)
  const long double Fi11=xx/(-t1)*(-0.5*(z*z+(1.-z)*(1.-z))*std::log(MUF*xx/Q1)+z*(1.-z))
			*2.*_Cf*_Cf*(xx*xx*xp*xp*z*z/(t1*t1)+t1*t1)/(z)+xx/(-t1)*_Cf*(1.-z)*_Cf*(xx*xx*xp*xp*z*z/(t1*t1)+z*z)/(-t1)
			+xx/(-t1)*(-_Cf*(1.+(1.-z)*(1.-z))/z*std::log(MUF*xx/Q1)+_Cf*z)
			*_Nc*(std::pow(xx,4)+std::pow(z,4)+std::pow(z*xx*xp/t1,4)+std::pow(t1,4))/(z*z*xx*xp);
  const long double Fi12=xx/(-t2)*(-0.5*(z*z+(1.-z)*(1.-z))*std::log(MUF*xx/Q2)+z*(1.-z))
			*2.*_Cf*_Cf*(xx*xx*xp*xp*z*z/(t2*t2)+t2*t2)/(z)+xx/(-t2)*_Cf*(1.-z)*_Cf*(xx*xx*xp*xp*z*z/(t2*t2)+z*z)/(-t2)
			+xx/(-t2)*(-_Cf*(1.+(1.-z)*(1.-z))/z*std::log(MUF*xx/Q2)+_Cf*z)
			*_Nc*(std::pow(xx,4)+std::pow(z,4)+std::pow(z*xx*xp/t2,4)+std::pow(t2,4))/(z*z*xx*xp);
  ris+=(Jac2*Fi12-Jac1*Fi11); // OK 6 Nc Cf Log(x)^2/xp
  
  //Separate Divergent part (Q^2->0) see Grazzini
  const long double tiny=1e-4;
  const long double Fi21=(Q1>tiny*xx*xp) ? xx* _Nc*_Cf*(((-(t1/z)-std::pow(t1/z,3)-Q1*Q1*Q1*u1/zb1-Q1*std::pow(u1/zb1,3))*(Q1+Qt1))/(Q1*Qt1)
			-(2.*xx*xx*((xx-t1)*(xx-t1)+t1*t1))/(u1*(xx-u1)))/(xx*xp)*std::log(xx*xp/Qt1) : (t1*t1*t1+t1*z*z)/(z*z*z*xx*xp*xp) ;
  const long double Fi22=(Q2>tiny*xx*xp) ? xx* _Nc*_Cf*(((-(t2/z)-std::pow(t2/z,3)-Q2*Q2*Q2*u2/zb2-Q2*std::pow(u2/zb2,3))*(Q2+Qt2))/(Q2*Qt2)
			-(2.*xx*xx*((xx-t2)*(xx-t2)+t2*t2))/(u2*(xx-u2)))/(xx*xp)*std::log(xx*xp/Qt2) : (t2*t2*t2+t2*z*z)/(z*z*z*xx*xp*xp) ;
  ris+=(Jac2*Fi22-Jac1*Fi21); // OK -10 Nc Cf Log(x)^2/xp
  
  const long double C1pm1=-2.*xx*gsl_sf_log(xx*xp/Qt1);
  const long double C1pm2=-2.*xx*gsl_sf_log(xx*xp/Qt2);
  ris+=_Cf*_Cf*(Jac2*(C1pm2)-Jac1*(C1pm1));
  
  const long double C1mp1=xx*xp*Q1-3.*xx*xp*t1*t1/2./u1+xx*xp/(2.*Qt1)*(pow(Q1-t1,3)+Q1*pow(Q1-u1,3))*(-3.+10.*Q1/Qt1-6.*Q1*Q1/Qt1/Qt1);
  const long double C1mp2=xx*xp*Q2-3.*xx*xp*t2*t2/2./u2+xx*xp/(2.*Qt2)*(pow(Q2-t2,3)+Q1*pow(Q2-u2,3))*(-3.+10.*Q2/Qt2-6.*Q2*Q2/Qt2/Qt2);
  ris+=_Cf*_Cf*(1./xp/Q2*Jac2*(C1mp2)-1./xp/Q1*Jac1*(C1mp1));
  
  const long double C1pp1=-3.*xx*xp/(2.*u1)-xx*xp*Q1*A1*A1/u1*L2b1+xx*xp/(u1*t1*t1)*L1a1*((A1-xx)*(A1-xx)*t1*t1-2.*Q1*xx*u1*t1*A1-Q1*xx*xx*(Q1-t1)*u1)
			  +xx*xp/(u1*B1)*L31*((1.+Q1-xx)*((A1-xx)*(A1-xx)+Q1*A1*A1)-4.*Q1*A1*(A1-xx))+1./2.*xx*xp*(Q1-t1)*(Q1-t1)*(1./(Q1-u1)-Q1/u1)*(-3.+10.*Q1/Qt1-6.*Q1*Q1/Qt1/Qt1)
			  +1./2.*xx*xp*(Q1-u1)*(Q1-u1)*(Q1/(Q1-t1)+1./u1)*(-3.+10.*Q1/Qt1-6.*Q1*Q1/Qt1/Qt1)+xx*xp*(Q1-t1)/(u1*(Q1-u1))*(2.*u1*(1.+t1)
			  -Q1*(4.*xx*Q1-Q1*t1-xx*u1-2.*u1*t1))*(-1.+2.*Q1/Qt1)+xx*xp*(Q1-u1)/(u1*(Q1-t1))*(t1-2.*u1*t1+2.*Q1*u1*(Q1-t1))*(-1.+2.*Q1/Qt1)
			  +xx*xp*Q1*xx*(t1+u1)/t1-2.*xx*xp*Q1*Q1*xx*xx/(t1*t1)+xx*xp*Q1*xx*xx/(2.*u1*u1)+xx*xp/(2.*u1)*(-2.*(Q1+xx)*u1+2.*xx+xx*xx+xx*Q1*(2.*(1.-Q1)+3.*xx-u1)+5.*Q1*(1.-Q1));
  const long double C1pp2=-3.*xx*xp/(2.*u2)-xx*xp*Q2*A2*A2/u2*L2b2+xx*xp/(u2*t2*t2)*L1a2*((A2-xx)*(A2-xx)*t2*t2-2.*Q2*xx*u2*t2*A2-Q2*xx*xx*(Q2-t2)*u2)
			  +xx*xp/(u2*B2)*L32*((1.+Q2-xx)*((A2-xx)*(A2-xx)+Q2*A2*A2)-4.*Q2*A2*(A2-xx))+1./2.*xx*xp*(Q2-t2)*(Q2-t2)*(1./(Q2-u2)-Q2/u2)*(-3.+10.*Q2/Qt2-6.*Q2*Q2/Qt2/Qt2)
			  +1./2.*xx*xp*(Q2-u2)*(Q2-u2)*(Q2/(Q2-t2)+1./u2)*(-3.+10.*Q2/Qt2-6.*Q2*Q2/Qt2/Qt2)+xx*xp*(Q2-t2)/(u2*(Q2-u2))*(2.*u2*(1.+t2)
			  -Q2*(4.*xx*Q2-Q2*t2-xx*u2-2.*u2*t2))*(-1.+2.*Q2/Qt2)+xx*xp*(Q2-u2)/(u2*(Q2-t2))*(t2-2.*u2*t2+2.*Q2*u2*(Q2-t2))*(-1.+2.*Q2/Qt2)
			  +xx*xp*Q2*xx*(t2+u2)/t2-2.*xx*xp*Q2*Q2*xx*xx/(t2*t2)+xx*xp*Q2*xx*xx/(2.*u2*u2)+xx*xp/(2.*u2)*(-2.*(Q2+xx)*u2+2.*xx+xx*xx+xx*Q2*(2.*(1.-Q2)+3.*xx-u2)+5.*Q2*(1.-Q2));
  ris+=_Cf*_Cf*(1./xp/Q2*Jac2*(C1pp2)-1./xp/Q1*Jac1*(C1pp1));
  
  const long double C1mm1=xx*xp*t1*t1/u1*L1a1-xx*xp*Q1*(xx-t1)*(xx-t1)/u1*L2b1+xx*xp*xx*Q1/B1/B1*(t1*(u1+t1)-2.*xx*Q1)
			  +xx*xp/(u1*B1)*(t1*t1*B1*B1-xx*t1*t1*(u1+t1)+2.*Q1*Q1*xx*xx+Q1*xx*xx*(3.*t1-u1)+Q1*xx*xx*u1/(B1*B1)*(-t1*(u1+t1)+2.*xx*Q1+Q1*(t1-u1)))*L31;
  const long double C1mm2=xx*xp*t2*t2/u2*L1a2-xx*xp*Q2*(xx-t2)*(xx-t2)/u2*L2b2+xx*xp*xx*Q2/B2/B2*(t2*(u2+t2)-2.*xx*Q2)
			  +xx*xp/(u2*B2)*(t2*t2*B2*B2-xx*t2*t2*(u2+t2)+2.*Q2*Q2*xx*xx+Q2*xx*xx*(3.*t2-u2)+Q2*xx*xx*u2/(B2*B2)*(-t2*(u2+t2)+2.*xx*Q2+Q2*(t2-u2)))*L32;
  ris+=_Cf*_Cf*(1./xp/Q2*Jac2*(C1mm2)-1./xp/Q1*Jac1*(C1mm1));
  
  const long double C2mp1=xx*xp*Q1/Qt1/Qt1*(pow(Q1-t1,3)+Q1*pow(Q1-u1,3))*(-2.+3.*Q1/Qt1)+2.*xx*xp*Q1+4.*xx*xp*Q1*gsl_sf_log(xx*xp/Qt1);
  const long double C2mp2=xx*xp*Q2/Qt2/Qt2*(pow(Q2-t2,3)+Q2*pow(Q2-u2,3))*(-2.+3.*Q2/Qt2)+2.*xx*xp*Q2+4.*xx*xp*Q2*gsl_sf_log(xx*xp/Qt2);
  ris+=_Nc*_Cf*(1./xp/Q2*Jac2*(C2mp2)-1./xp/Q1*Jac1*(C2mp1));
  
  const long double C2pp1=1./2.*xx*xp*A1*A1*(1.-z)*(L2a1-L1a1)+(xx*xp*(xx-t1)*A1*A1*(1.-zb1))/(2.*u1)*(L1b1-L2b1)
			  +xx*xp*pow(1.-Q1,3)/(2.*u1)*(L1b1-L1a1)+xx*xp*Q1*A1*A1/(xx-u1)*(L1b1+L2a1)
			  +xx*xp*Q1/(u1*u1)*L1b1*(2.*u1*(1.-Q1)*(1.-Q1)+4.*xx*(xx-t1)*A1-2.*xx*xx*(Q1-t1)-xx*xx*xx)
			  -1./2.*xx*xp*pow(Q1-t1,2)*(1./(Q1-u1)-Q1/u1)*(-3.+10.*Q1/Qt1-6.*Q1*Q1/Qt1/Qt1)
			  -1./2.*xx*xp*pow(Q1-u1,2)*(Q1/(Q1-t1)+1./u1)*(-3.+10.*Q1/Qt1-6.*Q1*Q1/Qt1/Qt1)
			  +xx*xp*(Q1-t1)/(2.*u1*(Q1-u1))*((-3.*xx*xp+u1*u1-Q1*Q1)*(1.+Q1)-xx*u1+2.*Q1*Q1*(Q1-u1))*(-1.+2.*Q1/Qt1)
			  +xx*xp*(Q1-u1)/(2.*u1*(Q1-t1))*(3.*xx*xp*(1.+Q1)-u1*Q1*(Q1-t1)+3.*Q1*xx*(Q1-u1)+t1*(u1+1.))*(-1.+2.*Q1/Qt1)
			  +xx*xp*xx*Q1/(2.*u1*u1)*(2.*(1.-Q1)*(1.-Q1)-2.*xx*(1.-xx)-u1*(Q1-u1)-4.*xx*Q1)+xx*xp*(u1-t1)*(xx+Q1)/(2.*u1)
			  -8.*pow(xx*xx*xp*Q1,2)/pow(u1,4)-2.*pow(xx*xx*xp*Q1,2)/pow(u1,4)*L1b1;
  const long double C2pp2=1./2.*xx*xp*A2*A2*(1.-z)*(L2a2-L1a2)+(xx*xp*(xx-t2)*A2*A2*(1.-zb2))/(2.*u2)*(L1b2-L2b2)
			  +xx*xp*pow(1.-Q2,3)/(2.*u2)*(L1b2-L1a2)+xx*xp*Q2*A2*A2/(xx-u2)*(L1b2+L2a2)
			  +xx*xp*Q2/(u2*u2)*L1b2*(2.*u2*(1.-Q2)*(1.-Q2)+4.*xx*(xx-t2)*A2-2.*xx*xx*(Q2-t2)-xx*xx*xx)
			  -1./2.*xx*xp*pow(Q2-t2,2)*(1./(Q2-u2)-Q2/u2)*(-3.+10.*Q2/Qt2-6.*Q2*Q2/Qt2/Qt2)
			  -1./2.*xx*xp*pow(Q2-u2,2)*(Q2/(Q2-t2)+1./u2)*(-3.+10.*Q2/Qt2-6.*Q2*Q2/Qt2/Qt2)
			  +xx*xp*(Q2-t2)/(2.*u2*(Q2-u2))*((-3.*xx*xp+u2*u2-Q2*Q2)*(1.+Q2)-xx*u2+2.*Q2*Q2*(Q2-u2))*(-1.+2.*Q2/Qt2)
			  +xx*xp*(Q2-u2)/(2.*u2*(Q2-t2))*(3.*xx*xp*(1.+Q2)-u2*Q2*(Q2-t2)+3.*Q2*xx*(Q2-u2)+t2*(u2+1.))*(-1.+2.*Q2/Qt2)
			  +xx*xp*xx*Q2/(2.*u2*u2)*(2.*(1.-Q2)*(1.-Q2)-2.*xx*(1.-xx)-u2*(Q2-u2)-4.*xx*Q2)+xx*xp*(u2-t2)*(xx+Q2)/(2.*u2)
			  -8.*pow(xx*xx*xp*Q2,2)/pow(u2,4)-2.*pow(xx*xx*xp*Q2,2)/pow(u2,4)*L1b2;
  ris+=_Nc*_Cf*(1./xp/Q2*Jac2*(C2pp2)-1./xp/Q1*Jac1*(C2pp1));
  
  const long double C2mm1=xx*xp*t1*t1*Q1/(2.*u1)*(L1a1+3.*L1b1)+1./2.*xx*xp*t1*t1*(1.-z)*(L2a1-L1a1)+xx*xp*Q1*pow(xx-t1,3)*zb1/(2.*u1*u1)*(L2b1-L1b1)+xx*xp*t1*t1*Q1/(xx-u1)*(L1b1+L2a1)
			  +xx*xp*t1*t1/(2.*u1)*(L1b1-L1a1)+xx*xp*xx*Q1/(u1*u1)*(4.*t1*(t1-xx)+xx*xx)*L1b1-2.*pow(xx*xx*xp*Q1,2)/pow(u1,4)*L1b1+xx*xp*t1*t1*xx*Q1/(u1*u1)-8.*pow(xx*xx*xp*Q1,2)/(pow(u1,4));
  const long double C2mm2=xx*xp*t2*t2*Q2/(2.*u2)*(L1a2+3.*L1b2)+1./2.*xx*xp*t2*t2*(1.-z)*(L2a2-L1a2)+xx*xp*Q2*pow(xx-t2,3)*zb2/(2.*u2*u2)*(L2b2-L1b2)+xx*xp*t2*t2*Q2/(xx-u2)*(L1b2+L2a2)
			  +xx*xp*t2*t2/(2.*u2)*(L1b2-L1a2)+xx*xp*xx*Q2/(u2*u2)*(4.*t2*(t2-xx)+xx*xx)*L1b2-2.*pow(xx*xx*xp*Q2,2)/pow(u2,4)*L1b2+xx*xp*t2*t2*xx*Q2/(u2*u2)-8.*pow(xx*xx*xp*Q2,2)/(pow(u2,4));
  ris+=_Nc*_Cf*(1./xp/Q2*Jac2*(C2mm2)-1./xp/Q1*Jac1*(C2mm1));
  
  const long double C2epsilon1=4.*pow(xx*xx*xp*Q1,2)/(pow(u1,4));
  const long double C2epsilon2=4.*pow(xx*xx*xp*Q2,2)/(pow(u2,4));
  ris+=_Nc*_Cf*(1./xp/Q2*Jac2*(C2epsilon2)-1./xp/Q1*Jac1*(C2epsilon1));// All C2's OK -2 Nc Cf Log(x)^2/xp
  
  ris*=(1.-zmin);
  return ris;
}

// qqbar channel

long double NLOPL::NLO_PL_notsing_doublediff_qqbar(double xp, double zz){
  
  long double ris=0.0;
  
  //Set Energy, Transverse Momentum, Integration variable z
  const long double xx=x;
  const long double Nc=_Nc;
  const long double Nf=_Nf;
  const long double MUF=std::pow(_muF/_mH,2);
  const long double MUR=std::pow(_muR/_mH,2);
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
  
  //Define and add regular parts in G2s (refer to Glover paper for definition of G2s)
  const long double Fi11= xx/(-t1)*_Cf*(1.-z)*2.*_Cf*_Cf*(t1*t1+z*z*xx*xx*xp*xp/(t1*t1))/z
			+ xx/(-t1)*(-_Cf*(1.+(1.-z)*(1.-z))/z*std::log(MUF*xx/Q1)+_Cf*z)*_Cf*(z*z+t1*t1)/(-xx*xp*z/(t1));
  const long double Fi12= xx/(-t2)*_Cf*(1.-z)*2.*_Cf*_Cf*(t2*t2+z*z*xx*xx*xp*xp/(t2*t2))/z
			+ xx/(-t2)*(-_Cf*(1.+(1.-z)*(1.-z))/z*std::log(MUF*xx/Q2)+_Cf*z)*_Cf*(z*z+t2*t2)/(-xx*xp*z/(t2));
  ris+=2.*(Jac2*Fi12-Jac1*Fi11);
  
  const long double Fi21=1./xp*2.*_Cf*_Cf*((1.-Q1)*(1.-Q1)+(u1+t1-2.*Q1)*(u1+t1-2.*Q1))*std::log(xx*xp/Qt1);
  const long double Fi22=1./xp*2.*_Cf*_Cf*((1.-Q2)*(1.-Q2)+(u2+t2-2.*Q2)*(u2+t2-2.*Q2))*std::log(xx*xp/Qt2);
  ris+=2.*(Jac2*Fi22-Jac1*Fi21);
  
  const long double D1pm1=-xx*xp*Q1*(1.+std::log(xx*xp/Qt1))-(xx*xp*xx*xp*z*(1.-z))/(t1)-xx*xp*(1.-z)*(1.-z);
  const long double D1pm2=-xx*xp*Q2*(1.+std::log(xx*xp/Qt2))-(xx*xp*xx*xp*z*(1.-z))/(t2)-xx*xp*(1.-z)*(1.-z);
  ris+=2.*2.*_Cf*_Cf*_Cf*(Jac2/xp/Q2*D1pm2-Jac1/xp/Q1*D1pm1);
  
  const long double D2pm1=-1./3.*xx*xp*Q1-xx*xp*t1*t1/(6.*z*z)*(11.-12.*Q1/Qt1+3.*Q1*Q1/Qt1/Qt1)+11.*xx*xp*t1*t1/6.;
  const long double D2pm2=-1./3.*xx*xp*Q2-xx*xp*t2*t2/(6.*z*z)*(11.-12.*Q2/Qt2+3.*Q2*Q2/Qt2/Qt2)+11.*xx*xp*t2*t2/6.;
  ris+=2.*2.*_Nc*_Cf*_Cf*(Jac2/xp/Q2*D2pm2-Jac1/xp/Q1*D2pm1);
  
  const long double D1pp1=xx*xp*u1*u1*(1.-zb1)/A1*(L2b1-L1b1)+xx*xp*(xx-u1)*(xx-u1)*(1.-z)/A1*(L2a1-L1a1)
			 +xx*xp*xx*Q1*(xx*xp+u1*t1)/(t1*t1)*L1a1-2.*xx*xp*xx*xx*Q1/(A1*B1)*L31+xx*xp*xx*Q1*(2.*xx*xp-u1*t1)/(t1*t1);
  const long double D1pp2=xx*xp*u2*u2*(1.-zb2)/A2*(L2b2-L1b2)+xx*xp*(xx-u2)*(xx-u2)*(1.-z)/A2*(L2a2-L1a2)
			 +xx*xp*xx*Q2*(xx*xp+u2*t2)/(t2*t2)*L1a2-2.*xx*xp*xx*xx*Q2/(A2*B2)*L32+xx*xp*xx*Q2*(2.*xx*xp-u2*t2)/(t2*t2);
  ris+=2.*2.*_Cf*_Cf*_Cf*(Jac2/xp/Q2*D1pp2-Jac1/xp/Q1*D1pp1);
  
  const long double D2pp1=xx*xp*u1*u1*(xx-t1)*(1.-zb1)/(2.*A1)*(L1b1-L2b1)+xx*xp*std::pow(xx-u1,3)*(1.-z)/(2.*A1)*(L1a1-L2a1)
			 -0.5*xx*xp*u1*u1*(L1a1+L1b1)+6.*std::pow(xx*xp*xx*Q1,2)/std::pow(B1,4)-xx*xp*xx*Q1*u1*u1/B1/B1
			 +L31*(xx*xp*u1*u1*(u1+t1)/B1+std::pow(xx*xp*xx*Q1,2)*(1./(A1*B1*B1*B1)-3.*A1/std::pow(B1,5))
			 +xx*xp*xx*Q1*((t1-3.*u1)/(2.*B1)+A1*(B1*B1+2.*u1*u1)/(4.*B1*B1*B1)+(t1*t1-6.*t1*u1+7.*u1+u1)/(4.*A1*B1)));
  const long double D2pp2=xx*xp*u2*u2*(xx-t2)*(1.-zb2)/(2.*A2)*(L1b2-L2b2)+xx*xp*std::pow(xx-u2,3)*(1.-z)/(2.*A2)*(L1a2-L2a2)
			 -0.5*xx*xp*u2*u2*(L1a2+L1b2)+6.*std::pow(xx*xp*xx*Q2,2)/std::pow(B2,4)-xx*xp*xx*Q2*u2*u2/B2/B2
			 +L32*(xx*xp*u2*u2*(u2+t2)/B2+std::pow(xx*xp*xx*Q2,2)*(1./(A2*B2*B2*B2)-3.*A2/std::pow(B2,5))
			 +xx*xp*xx*Q2*((t2-3.*u2)/(2.*B2)+A2*(B2*B2+2.*u2*u2)/(4.*B2*B2*B2)+(t2*t2-6.*t2*u2+7.*u2+u2)/(4.*A2*B2)));
  ris+=2.*2.*_Nc*_Cf*_Cf*(Jac2/xp/Q2*D2pp2-Jac1/xp/Q1*D2pp1);
  
  const long double E11=8./3.*xx-4./3.*xx*xx;
  ris+=_Nf*_Cf*_Cf*E11*(Jac2-Jac1);
  
  const long double E21=xx*(Q1-xx*xp)/(Qt1*Qt1)*(std::pow(t1/z,2)+std::pow(u1/zb1,2))+2.*xx;
  const long double E22=xx*(Q2-xx*xp)/(Qt2*Qt2)*(std::pow(t2/z,2)+std::pow(u2/zb2,2))+2.*xx;
  ris+=_Cf*_Cf*(Jac2*E22-Jac1*E21);
  
  const long double E31=-2.*xx*xp*(std::pow(u1+t1-2.*Q1,2)-2.*xx*xp)*std::log(xx*xp/Qt1)
			-xx*xp*Q1*(2.*Qt1+Q1)/(Qt1*Qt1)*(std::pow(t1/z,2)+std::pow(u1/zb1,2))-6.*xp*xx*Q1;
  const long double E32=-2.*xx*xp*(std::pow(u2+t2-2.*Q2,2)-2.*xx*xp)*std::log(xx*xp/Qt2)
			-xx*xp*Q2*(2.*Qt2+Q2)/(Qt2*Qt2)*(std::pow(t2/z,2)+std::pow(u2/zb2,2))-6.*xp*xx*Q2;
  ris+=_Cf*_Cf/_Nc*(Jac2/(xp*Q2)*E32-Jac1/(xp*Q1)*E31);
  			 
  ris*=(1.-zmin);
  return ris;  
}

//qq channel

long double NLOPL::NLO_PL_notsing_doublediff_qq(double xp, double zz){
  
  long double ris=0.0;
  
  //Set Energy, Transverse Momentum, Integration variable z
  const long double xx=x;
  const long double Nc=_Nc;
  const long double Nf=_Nf;
  const long double MUF=std::pow(_muF/_mH,2);
  const long double MUR=std::pow(_muR/_mH,2);
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
  
  //Define and add regular parts in G2s (refer to Glover paper for definition of G2s)
  const long double Fi11=xx/(-t1)*(-_Cf*(1.+(1.-z)*(1.-z))/z*std::log(MUF*xx/Q1)+_Cf*z)*_Cf*(t1*t1+z*z)/(-xx*xp*z/(t1));
  const long double Fi12=xx/(-t2)*(-_Cf*(1.+(1.-z)*(1.-z))/z*std::log(MUF*xx/Q2)+_Cf*z)*_Cf*(t2*t2+z*z)/(-xx*xp*z/(t2));
  ris+=2.*(Jac2*Fi12-Jac1*Fi11);
  
  const long double Fi21=xx*2.*_Cf*_Cf*((1.-Q1)*(1.-Q1)+std::pow(u1+t1-2.*Q1,2))/(xx*xp)*std::log(xx*xp/Qt1);
  const long double Fi22=xx*2.*_Cf*_Cf*((1.-Q2)*(1.-Q2)+std::pow(u2+t2-2.*Q2,2))/(xx*xp)*std::log(xx*xp/Qt2);
  ris+=(Jac2*Fi22-Jac1*Fi21);
  
  const long double E41=2.*xx*(1.+Q1*Q1)/Qt1*std::log(x*xp/Q1)+4.*xx*std::log(xx*xp/Qt1);
  const long double E42=2.*xx*(1.+Q2*Q2)/Qt2*std::log(x*xp/Q2)+4.*xx*std::log(xx*xp/Qt2);
  ris+=_Cf*_Cf/_Nc*(Jac2*E42-Jac1*E41);
  
  const long double E21=xx*(Q1-xx*xp)/(Qt1*Qt1)*(std::pow(t1/z,2)+std::pow(u1/zb1,2))+2.*xx;
  const long double E22=xx*(Q2-xx*xp)/(Qt2*Qt2)*(std::pow(t2/z,2)+std::pow(u2/zb2,2))+2.*xx;
  ris+=_Cf*_Cf*(Jac2*E22-Jac1*E21);
  
  ris*=(1.-zmin);
  return ris;
  
  
}

//qqprime channel

long double NLOPL::NLO_PL_notsing_doublediff_qqprime(double xp, double zz){
  
  long double ris=0.0;
  
  //Set Energy, Transverse Momentum, Integration variable z
  const long double xx=x;
  const long double Nc=_Nc;
  const long double Nf=_Nf;
  const long double MUF=std::pow(_muF/_mH,2);
  const long double MUR=std::pow(_muR/_mH,2);
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
  
  //Define and add regular parts in G2s (refer to Glover paper for definition of G2s)
  const long double Fi11=xx/(-t1)*(-_Cf*(1.+(1.-z)*(1.-z))/z*std::log(MUF*xx/Q1)+_Cf*z)*_Cf*(t1*t1+z*z)/(-xx*xp*z/(t1));
  const long double Fi12=xx/(-t2)*(-_Cf*(1.+(1.-z)*(1.-z))/z*std::log(MUF*xx/Q2)+_Cf*z)*_Cf*(t2*t2+z*z)/(-xx*xp*z/(t2));
  ris+=2.*(Jac2*Fi12-Jac1*Fi11);
  
  const long double Fi21=xx*2.*_Cf*_Cf*((1.-Q1)*(1.-Q1)+std::pow(u1+t1-2.*Q1,2))/(xx*xp)*std::log(xx*xp/Qt1);
  const long double Fi22=xx*2.*_Cf*_Cf*((1.-Q2)*(1.-Q2)+std::pow(u2+t2-2.*Q2,2))/(xx*xp)*std::log(xx*xp/Qt2);
  ris+=(Jac2*Fi22-Jac1*Fi21);
  
  const long double E21=xx*(Q1-xx*xp)/(Qt1*Qt1)*(std::pow(t1/z,2)+std::pow(u1/zb1,2))+2.*xx;
  const long double E22=xx*(Q2-xx*xp)/(Qt2*Qt2)*(std::pow(t2/z,2)+std::pow(u2/zb2,2))+2.*xx;
  ris+=_Cf*_Cf*(Jac2*E22-Jac1*E21);
  
  ris*=(1.-zmin);
  return ris;
}

long double NLOPL::NLO_PL_notsing_doublediff_qqbarprime(double xp, double zz)
{
  return(NLO_PL_notsing_doublediff_qqprime(xp,zz));
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
    case(4):{
      return 0;
      break;
    }
    case(5):{
      return 0;
      break;
    }
    case(6):{
      return 0;
      break;
    }
  }
}



















		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		