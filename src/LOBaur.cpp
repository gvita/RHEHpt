#include "complex_def.h"
#include "LOBaur.h"

//Auxiliary functions

std::complex<long double> B_Aux(std::complex<long double> v){
  std::complex<long double> res;
  if (std::real(v)<0.){
    res=-std::sqrt(1.-4./v)*log_c((1.+std::sqrt(1.-4./v))/(std::sqrt(1.-4./v)-1.));
  }
  else {
    if ((std::real(v)<4.)&&(std::real(v)>0.)){
      res=std::sqrt(4./v-1.)*(2.*std::acos(std::real(std::sqrt(v)/2.))-M_PIl);
    }
    else {
      res=-std::sqrt(1.-4./v)*(log_c((1.+std::sqrt(1.-4./v))/(1.-std::sqrt(1.-4./v)))-II*M_PIl);    }
  }
  return res;
}

std::complex<long double> C_Aux(std::complex<long double> v){
  std::complex<long double> res;
  if (std::real(v)<0.){
    res=std::pow(log_c((1.+std::sqrt(1.-4./v))/(std::sqrt(1.-4./v)-1.)),2);
  }
  else {
    if ((std::real(v)<4.)&&(std::real(v)>0.)){
      res=-std::pow(2.*std::acos(std::real(std::sqrt(v)*0.5))-M_PIl,2);
    }
    else {
      res=std::pow(log_c((1.+std::sqrt(1.-4./v))/(1.-std::sqrt(1.-4./v))),2)-M_PIl*M_PIl
      -2.*II*M_PIl*log_c((1.+std::sqrt(1.-4./v))/(1.-std::sqrt(1.-4./v)));
    }
  }
  return res;
}

std::complex<long double> D_Aux(std::complex<long double> s, std::complex<long double> t, std::complex<long double> u, std::complex<long double> v){
  std::complex<long double> beta=0.5*(1.+std::sqrt(1.+4.*t/(u*s)));
  std::complex<long double> res;
  if (std::real(v)<0.){
    std::complex<long double> gamm=0.5*(1.+std::sqrt(1.-4./v));
    res =-dilog_c(gamm/(gamm+beta-1.))+dilog_c((gamm-1.)/(gamm+beta-1))+dilog_c((beta-gamm)/beta)-dilog_c((beta-gamm)/(beta-1.))
    +(std::pow(log_c(beta),2)-std::pow(log_c(beta-1.),2))/2.+log_c(gamm)*log_c((gamm+beta-1.)/beta)+log_c(gamm-1.)*log_c((beta-1.)/(gamm+beta-1));
  }
  else {
    if ((std::real(v)<4.)&&(std::real(v)>0.)){
      std::complex<long double> a=std::sqrt(4./v-1.);
      std::complex<long double> r=std::sqrt((a*a+1.)/(a*a+std::pow(2.*beta-1.,2)));
      std::complex<long double> phi=std::acos(std::real(r*(a*a+2.*beta-1.)/(1.+a*a)));
      if(std::abs(beta)-1. < 1e-3) phi = acos(1.-real(2*pow(a*(beta-1.)/(1+a*a),2))); // expansion of the argument of acos
      std::complex<long double> theta=std::acos(std::real(r*(a*a-2.*beta+1.)/(1.+a*a)));
      res=dilog_c(r*std::exp(II*theta))+dilog_c(r*std::exp(-II*theta))-dilog_c(r*std::exp(II*phi))-dilog_c(r*std::exp(-II*phi))
      +(phi-theta)*(phi+theta-M_PIl);
    }
    else{
      std::complex<long double> gamm=0.5*(1.+std::sqrt(1.-4./v));
      res=-dilog_c(gamm/(gamm+beta-1.))+dilog_c((gamm-1.)/(gamm+beta-1.))+dilog_c(gamm/(gamm-beta))
      -dilog_c((gamm-1.)/(gamm-beta))+log_c(gamm/(1.-gamm))*log_c((gamm+beta-1.)/(beta-gamm))
      -II*M_PIl*log_c((gamm+beta-1.)/(beta-gamm));
      
    }
  }
  return (res*2./(2.*beta-1.));
}



//Functions B1, C, C1, D, E

std::complex<long double> LOBaur::B1(long double v, long double y) {
  std::complex<long double> v_c=static_cast<std::complex<long double> > (v/(y*x));
  std::complex<long double> y_c=static_cast<std::complex<long double> > (1./y);
  return (B_Aux(v_c)-B_Aux(y_c));
}
std::complex<long double> LOBaur::C(long double v, long double y) {
  std::complex<long double> v_c=static_cast<std::complex<long double> > (v/(y*x));
  return (1./(2.*v)*C_Aux(v_c));
}
std::complex<long double> LOBaur::C1(long double v, long double y){
  return(v/(v-x)*C(v,y)-x/(v-x)*C(x,y));
}  
std::complex<long double> LOBaur::D(long double v, long double w, long double y){
  std::complex<long double> v_c=static_cast<std::complex<long double> > (v/(y*x));
  std::complex<long double> w_c=static_cast<std::complex<long double> > (w/(y*x));
  std::complex<long double> u_c=static_cast<std::complex<long double> > (1./y*(1.-w/x-v/x));
  std::complex<long double> y_c=static_cast<std::complex<long double> > (1./y);
  return(1./(v*w)*(D_Aux(w_c,u_c,v_c,w_c)+D_Aux(w_c,u_c,v_c,v_c)-D_Aux(w_c,u_c,v_c,y_c)));
}

std::complex<long double> LOBaur::E(long double v, long double w, long double y){
  return(v*C(v,y)+w*C(w,y)+(v-x)*C1(v,y)+(w-x)*C1(w,y)-v*w*D(v,w,y));
}

//Helicity Amplitudes

std::complex<long double> LOBaur::Appp(long double y){
  std::complex<long double> tc= static_cast<std::complex<long double> > (t);
  std::complex<long double> uc = static_cast<std::complex<long double> > (u);
  std::complex<long double> t1c = static_cast<std::complex<long double> > (t1);
  std::complex<long double> u1c = static_cast<std::complex<long double> > (u1);
  std::complex<long double> s1c = static_cast<std::complex<long double> > (s1);
  std::complex<long double> sc = static_cast<std::complex<long double> > (1.);
  std::complex<long double> Mppp=std::sqrt(1./8.*sc*tc*uc)*(-64.*(1./(uc*tc)+1./(tc*t1c)+1./(uc*u1c))
  -64./sc*((2.*sc+tc)/(u1c*u1c)*B1(u,y)+(2.*sc+uc)/(t1c*t1c)*B1(t,y))
  -16.*(sc-4.*y*x)/(tc*uc*sc)*(s1c*C1(1.,y)+(uc-sc)*C1(t,y)+(tc-sc)*C1(u,y))
  -128.*y*x*(1./(tc*t1c)*C1(t,y)+1./(uc*u1c)*C1(u,y))+64.*y*x/sc*D(u,t,y)
  +8.*(sc-4.*y*x)/(sc*tc*uc)*(sc*tc*D(1.,t,y)+uc*sc*D(u,1.,y)-uc*tc*D(u,t,y))-32./sc/sc*E(u,t,y));
  return (Mppp*y*x);
  
  
  
  
}
std::complex<long double> LOBaur::Ampm(long double y){
  std::complex<long double> tc= static_cast<std::complex<long double> > (t);
  std::complex<long double> uc = static_cast<std::complex<long double> > (u);
  std::complex<long double> t1c = static_cast<std::complex<long double> > (t1);
  std::complex<long double> u1c = static_cast<std::complex<long double> > (u1);
  std::complex<long double> s1c = static_cast<std::complex<long double> > (s1);
  std::complex<long double> sc = static_cast<std::complex<long double> > (1.);
  std::complex<long double> Mmpm=std::sqrt(1./8.*sc*uc*tc)*(-64.*(1./(uc*sc)+1./(sc*s1c)+1./(uc*u1c))
  -64./tc*((2.*tc+sc)/(u1c*u1c)*B1(u,y)+(2.*tc+uc)/(s1c*s1c)*B1(1.,y))
  -16.*(tc-4.*y*x)/(sc*uc*tc)*(t1c*C1(t,y)+(uc-tc)*C1(1.,y)+(sc-tc)*C1(u,y))
  -128.*y*x*(1./(sc*s1c)*C1(1.,y)+1./(uc*u1c)*C1(u,y))+64.*y*x/tc*D(u,1.,y)
  +8.*(tc-4.*y*x)/(tc*sc*uc)*(tc*sc*D(t,1.,y)+uc*tc*D(u,t,y)-uc*sc*D(u,1.,y))-32./tc/tc*E(u,1.,y));
  return (Mmpm*y*x);
  
  
  
  
}
std::complex<long double> LOBaur::Ampp(long double y){
  std::complex<long double> tc= static_cast<std::complex<long double> > (t);
  std::complex<long double> uc = static_cast<std::complex<long double> > (u);
  std::complex<long double> t1c = static_cast<std::complex<long double> > (t1);
  std::complex<long double> u1c = static_cast<std::complex<long double> > (u1);
  std::complex<long double> s1c = static_cast<std::complex<long double> > (s1);
  std::complex<long double> sc = static_cast<std::complex<long double> > (1.);
  std::complex<long double> Mmpp=std::sqrt(1./8.*sc*tc*uc)*(-64.*(1./(sc*tc)+1./(tc*t1c)+1./(sc*s1c))
  -64./uc*((2.*uc+tc)/(s1c*s1c)*B1(1.,y)+(2.*uc+sc)/(t1c*t1c)*B1(t,y))
  -16.*(uc-4.*y*x)/(tc*sc*uc)*(u1c*C1(u,y)+(sc-uc)*C1(t,y)+(tc-uc)*C1(1.,y))
  -128.*y*x*(1./(tc*t1c)*C1(t,y)+1./(sc*s1c)*C1(1.,y))+64.*y*x/uc*D(1.,t,y)
  +8.*(uc-4.*y*x)/(uc*tc*sc)*(uc*tc*D(u,t,y)+sc*uc*D(1.,u,y)-sc*tc*D(1.,t,y))-32./uc/uc*E(1.,t,y));
  return (Mmpp*y*x);
  
  
  
  
}
std::complex<long double> LOBaur::Appm(long double y){
  std::complex<long double> tc = static_cast<std::complex<long double> > (t);
  std::complex<long double> uc = static_cast<std::complex<long double> > (u);
  std::complex<long double> t1c = static_cast<std::complex<long double> > (t1);
  std::complex<long double> u1c = static_cast<std::complex<long double> > (u1);
  std::complex<long double> s1c = static_cast<std::complex<long double> > (s1);
  std::complex<long double> sc = static_cast<std::complex<long double> > (1.);
  std::complex<long double> Mppm=std::sqrt(1./8.*sc*tc*uc)*(64.*x/(sc*tc*uc)
  +16.*(x-4.*y*x)/(sc*tc*uc)*(s1c*C1(1.,y)+u1c*C1(u,y)+t1c*C1(t,y))
  -8.*(x-4.*y*x)/(sc*tc*uc)*(sc*tc*D(1.,t,y)+uc*sc*D(u,1.,y)+uc*tc*D(u,t,y)));
  return (Mppm*y*x);
}

std::complex<long double> LOBaur::Aqqbar(long double y){
  std::complex<long double> tc = static_cast<std::complex<long double> > (t);
  std::complex<long double> uc = static_cast<std::complex<long double> > (u);
  std::complex<long double> t1c = static_cast<std::complex<long double> > (t1);
  std::complex<long double> u1c = static_cast<std::complex<long double> > (u1);
  std::complex<long double> s1c = static_cast<std::complex<long double> > (s1);
  std::complex<long double> sc = static_cast<std::complex<long double> > (1.);
  std::complex<long double> Aqqbar=y*x*(2.+2.*sc/s1c*B1(1.,y)+(4.*y*x-uc-tc)*C1(1.,y));
  return Aqqbar;
}

std::complex<long double> LOBaur::Agqbar(long double y){
  std::complex<long double> tc = static_cast<std::complex<long double> > (t);
  std::complex<long double> uc = static_cast<std::complex<long double> > (u);
  std::complex<long double> t1c = static_cast<std::complex<long double> > (t1);
  std::complex<long double> u1c = static_cast<std::complex<long double> > (u1);
  std::complex<long double> s1c = static_cast<std::complex<long double> > (s1);
  std::complex<long double> sc = static_cast<std::complex<long double> > (1.);
  std::complex<long double> Agqbar=y*x*(2.+2.*tc/t1c*B1(t,y)+(4.*y*x-sc-uc)*C1(t,y));
  return Agqbar;
}

std::complex<long double> LOBaur::Aqg(long double y){
  std::complex<long double> tc = static_cast<std::complex<long double> > (t);
  std::complex<long double> uc = static_cast<std::complex<long double> > (u);
  std::complex<long double> t1c = static_cast<std::complex<long double> > (t1);
  std::complex<long double> u1c = static_cast<std::complex<long double> > (u1);
  std::complex<long double> s1c = static_cast<std::complex<long double> > (s1);
  std::complex<long double> sc = static_cast<std::complex<long double> > (1.);
  std::complex<long double> Aqg=y*x*(2.+2.*uc/u1c*B1(u,y)+(4.*y*x-sc-tc)*C1(u,y));
  return Aqg;
}



long double LOBaur::MatrixElementFO(long double xp){
  if( xp > std::pow(1.-x,2.)/(4.*x) )  // process is kinematically forbidden
    return 0.;
  SetMandelstam(xp);
  long double rad=std::sqrt((1.-x)*(1.-x)-4.*x*xp);
  switch (_channel){
    case (1):{
      switch (_choice){
	case(0): {
	  return x/rad*( std::norm( Appp(yt)+Appp(yb) ) + std::norm( Appm(yt)+Appm(yb) ) +
	  std::norm( Ampp(yt)+Ampp(yb) ) + std::norm( Ampm(yt)+Ampm(yb) ) );
	  break;
	}
	case(1): {
	  return x/rad*( std::norm( Appp(yt) ) + std::norm( Appm(yt) ) +
	  std::norm( Ampp(yt) ) + std::norm( Ampm(yt) ) );
	  break;
	}
	case(2): {
	  return x/rad*( std::norm( Appp(yb) ) + std::norm( Appm(yb) ) +
	  std::norm( Ampp(yb) ) + std::norm( Ampm(yb) ) );
	  break;
	}
	case(3):{
	  return x/rad*std::real( std::conj( Appp(yb) )*Appp(yt)+std::conj( Appp(yt) )*Appp(yb) 
	  + std::conj( Appm(yb) )*Appm(yt) +std::conj( Appm(yt) )*Appm(yb)+
	  std::conj( Ampp(yb) )*Ampp(yt)+std::conj( Ampp(yt) )*Ampp(yb) 
	  + std::conj( Ampm(yb) )*Ampm(yt)+std::conj( Ampm(yt) )*Ampm(yb));
	  break;    
	}
      }
      break;
    }
    case(2):{
      switch (_choice){
	case(0): {
	  return -256./9.*x/rad*std::real((1.+t*t)/(u*(u-x)*(u-x))*(Aqg(yt)+Aqg(yb))*std::conj(Aqg(yt)+Aqg(yb))+(1.+u*u)/(t*(t-x)*(t-x))*(Agqbar(yt)+Agqbar(yb))*std::conj(Agqbar(yt)+Agqbar(yb)));
	  break;
	}
	case(1): {
	  return -256./9.*x/rad*std::real((1.+t*t)/(u*(u-x)*(u-x))*(Aqg(yt))*std::conj(Aqg(yt))+(1.+u*u)/(t*(t-x)*(t-x))*(Agqbar(yt))*std::conj(Agqbar(yt)));
	  break;
	}
	case(2): {
	  return -256./9.*x/rad*std::real((1.+t*t)/(u*(u-x)*(u-x))*(Aqg(yb))*std::conj(Aqg(yb))+(1.+u*u)/(t*(t-x)*(t-x))*(Agqbar(yb))*std::conj(Agqbar(yb)));
	  break;
	}
	case(3): {
	  return -256./9.*x/rad*std::real((1.+t*t)/(u*(u-x)*(u-x))*((Aqg(yt))*std::conj(Aqg(yb))+(Aqg(yb))*std::conj(Aqg(yt)))+(1.+u*u)/(t*(t-x)*(t-x))*(Agqbar(yt)*std::conj(Agqbar(yb))+Agqbar(yb)*std::conj(Agqbar(yt))));
	  break;
	}
      }
      break;
    }
    case(3):{
      switch (_choice){
	case(0): {
	  return 2048./27.*x/rad*std::real((u*u+t*t)/std::pow(1.-x,2)*(Aqqbar(yt)+Aqqbar(yb))*std::conj(Aqqbar(yt)+Aqqbar(yb)));
	  break;
	}
	case(1): {
	  return 2048./27.*x/rad*std::real((u*u+t*t)/std::pow(1.-x,2)*(Aqqbar(yt))*std::conj(Aqqbar(yt)));
	  break;
	}
	case(2): {
	  return 2048./27.*x/rad*std::real((u*u+t*t)/std::pow(1.-x,2)*(Aqqbar(yb))*std::conj(Aqqbar(yb)));
	  break;
	}
	case(3): {
	  return 2048./27.*x/rad*std::real(((u*u+t*t)/std::pow(1.-x,2)*(Aqqbar(yt)*std::conj(Aqqbar(yb))+Aqqbar(yb)*std::conj(Aqqbar(yt)))));
	  break;
	}
      }
      break;
    }
	  
  }
}
long double LOBaur::sigmapartLO(long double xp)
{
  const double K = 9./(16.*4.)*_as;
  return (K *0.5 * 3./M_PIl * MatrixElementFO(xp));
}

