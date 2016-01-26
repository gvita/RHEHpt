#include "complex_def.h"
#include "LOBaur.h"

//Functions B1, C, C1, D, E

std::complex<long double> LOBaur::B1(double v, double y) {
  std::complex<long double> rad_v =  static_cast<std::complex<long double> > (1.-4.*(y-std::complex<long double>(0.,1.)*epsilon)*x/v);
  std::complex<long double> rad_y =  static_cast<std::complex<long double> > (1.-4.*(y-std::complex<long double>(0.,1.)*epsilon));
  std::complex<long double> z=0.5*(1.+sqrt(rad_v));
  std::complex<long double> zy=0.5*(1.+sqrt(rad_y));
  return (-sqrt(rad_v)*log_c(-z/(1.-z))+sqrt(rad_y)*log_c(-zy/(1.-zy)));
}
std::complex<long double> LOBaur::C(double v, double y) {
  std::complex<long double> rad_v = static_cast<std::complex<long double> >  (1.-4.*(y-std::complex<long double>(0.,1.)*epsilon)*x/v);
  std::complex<long double> z=0.5*(1.+sqrt(rad_v));
  return (1./(2.*v)*pow(log_c(-z/(1.-z)),2));
}
std::complex<long double> LOBaur::C1(double v, double y){
  return(v/(v-x)*C(v,y)-x/(v-x)*C(x,y));
}  
std::complex<long double> LOBaur::D(double v, double w, double y){
  std::complex<long double> rad_vw = static_cast<std::complex<long double> > (1+4.*y*x*(x-v-w)/(v*w));
  std::complex<long double> rad_v = static_cast<std::complex<long double> > (1.-4.*(y-std::complex<long double>(0.,1.)*epsilon)*x/v);
  std::complex<long double> rad_w = static_cast<std::complex<long double> > (1.-4.*(y-std::complex<long double>(0.,1.)*epsilon)*x/w);
  std::complex<long double> rad_y = static_cast<std::complex<long double> > (1.-4.*(y-std::complex<long double>(0.,1.)*epsilon));
  std::complex<long double> xpi = 0.5*(1.+sqrt(rad_vw));
  std::complex<long double> xmi = 0.5*(1.-sqrt(rad_vw));
  std::complex<long double> yv = 0.5*(1.+sqrt(rad_v));
  std::complex<long double> yw = 0.5*(1.+sqrt(rad_w));
  std::complex<long double> yy = 0.5*(1.+sqrt(rad_y));
  std::complex<long double> log1 = (pow(real(xmi),2)  >1.e-34) ? log_c(-xmi/xpi)*log_c(1.-std::complex<long double>(0.,1.)*epsilon-v/(y*x)*xmi*xpi) :0.;
  std::complex<long double> log2 = (pow(real(xmi),2)  >1.e-34) ? log_c(-xmi/xpi)*log_c(1.-std::complex<long double>(0.,1.)*epsilon-w/(y*x)*xmi*xpi) :0.;
  std::complex<long double> log3 = (pow(real(xmi),2)  >1.e-34) ? log_c(-xmi/xpi)*log_c(1.-std::complex<long double>(0.,1.)*epsilon-1/(y)*xmi*xpi) :0.;
  
  std::complex<long double> Int1=2./(sqrt(rad_vw))*
  (dilog_c(xmi/(xmi-yv))-dilog_c(xpi/(xpi-yv))+dilog_c(xmi/(yv-xpi))
  -dilog_c(xpi/(yv-xmi))+log1);
  std::complex<long double> Int2=2./(sqrt(rad_vw))*
  (dilog_c(xmi/(xmi-yw))-dilog_c(xpi/(xpi-yw))+dilog_c(xmi/(yw-xpi))
  -dilog_c(xpi/(yw-xmi))+log2);
  std::complex<long double> Int3=2./(sqrt(rad_vw))*
  (dilog_c(xmi/(xmi-yy))-dilog_c(xpi/(xpi-yy))+dilog_c(xmi/(yy-xpi))
  -dilog_c(xpi/(yy-xmi))+log3); 
  
  return(1./(v*w)*(-Int3+Int1+Int2));
}

std::complex<long double> LOBaur::E(double v, double w, double y){
  return(v*C(v,y)+w*C(w,y)+(v-x)*C1(v,y)+(w-x)*C1(w,y)-v*w*D(v,w,y));
}

//Helicity Amplitudes

std::complex<long double> LOBaur::Appp(double y){
  std::complex<long double> tc= static_cast<std::complex<long double> > (t);
  std::complex<long double> uc = static_cast<std::complex<long double> > (u);
  std::complex<long double> t1c = static_cast<std::complex<long double> > (t1);
  std::complex<long double> u1c = static_cast<std::complex<long double> > (u1);
  std::complex<long double> s1c = static_cast<std::complex<long double> > (s1);
  std::complex<long double> sc = static_cast<std::complex<long double> > (1.);
  std::complex<long double> Mppp=sqrt(1./8.*sc*tc*uc)*(-64.*(1./(uc*tc)+1./(tc*t1c)+1./(uc*u1c))
  -64./sc*((2.*sc+tc)/(u1c*u1c)*B1(u,y)+(2.*sc+uc)/(t1c*t1c)*B1(t,y))
  -16.*(sc-4.*y*x)/(tc*uc*sc)*(s1c*C1(1.,y)+(uc-sc)*C1(t,y)+(tc-sc)*C1(u,y))
  -128.*y*x*(1./(tc*t1c)*C1(t,y)+1./(uc*u1c)*C1(u,y))+64.*y*x/sc*D(u,t,y)
  +8.*(sc-4.*y*x)/(sc*tc*uc)*(sc*tc*D(1.,t,y)+uc*sc*D(u,1.,y)-uc*tc*D(u,t,y))-32./sc/sc*E(u,t,y));
  return (Mppp*y*x);
  
  
  
  
}
std::complex<long double> LOBaur::Ampm(double y){
  std::complex<long double> tc= static_cast<std::complex<long double> > (t);
  std::complex<long double> uc = static_cast<std::complex<long double> > (u);
  std::complex<long double> t1c = static_cast<std::complex<long double> > (t1);
  std::complex<long double> u1c = static_cast<std::complex<long double> > (u1);
  std::complex<long double> s1c = static_cast<std::complex<long double> > (s1);
  std::complex<long double> sc = static_cast<std::complex<long double> > (1.);
  std::complex<long double> Mmpm=sqrt(1./8.*sc*uc*tc)*(-64.*(1./(uc*sc)+1./(sc*s1c)+1./(uc*u1c))
  -64./tc*((2.*tc+sc)/(u1c*u1c)*B1(u,y)+(2.*tc+uc)/(s1c*s1c)*B1(1.,y))
  -16.*(tc-4.*y*x)/(sc*uc*tc)*(t1c*C1(t,y)+(uc-tc)*C1(1.,y)+(sc-tc)*C1(u,y))
  -128.*y*x*(1./(sc*s1c)*C1(1.,y)+1./(uc*u1c)*C1(u,y))+64.*y*x/tc*D(u,1.,y)
  +8.*(tc-4.*y*x)/(tc*sc*uc)*(tc*sc*D(t,1.,y)+uc*tc*D(u,t,y)-uc*sc*D(u,1.,y))-32./tc/tc*E(u,1.,y));
  return (Mmpm*y*x);
  
  
  
  
}
std::complex<long double> LOBaur::Ampp(double y){
  std::complex<long double> tc= static_cast<std::complex<long double> > (t);
  std::complex<long double> uc = static_cast<std::complex<long double> > (u);
  std::complex<long double> t1c = static_cast<std::complex<long double> > (t1);
  std::complex<long double> u1c = static_cast<std::complex<long double> > (u1);
  std::complex<long double> s1c = static_cast<std::complex<long double> > (s1);
  std::complex<long double> sc = static_cast<std::complex<long double> > (1.);
  std::complex<long double> Mmpp=sqrt(1./8.*sc*tc*uc)*(-64.*(1./(sc*tc)+1./(tc*t1c)+1./(sc*s1c))
  -64./uc*((2.*uc+tc)/(s1c*s1c)*B1(1.,y)+(2.*uc+sc)/(t1c*t1c)*B1(t,y))
  -16.*(uc-4.*y*x)/(tc*sc*uc)*(u1c*C1(u,y)+(sc-uc)*C1(t,y)+(tc-uc)*C1(1.,y))
  -128.*y*x*(1./(tc*t1c)*C1(t,y)+1./(sc*s1c)*C1(1.,y))+64.*y*x/uc*D(1.,t,y)
  +8.*(uc-4.*y*x)/(uc*tc*sc)*(uc*tc*D(u,t,y)+sc*uc*D(1.,u,y)-sc*tc*D(1.,t,y))-32./uc/uc*E(1.,t,y));
  return (Mmpp*y*x);
  
  
  
  
}
std::complex<long double> LOBaur::Appm(double y){
  std::complex<long double> tc = static_cast<std::complex<long double> > (t);
  std::complex<long double> uc = static_cast<std::complex<long double> > (u);
  std::complex<long double> t1c = static_cast<std::complex<long double> > (t1);
  std::complex<long double> u1c = static_cast<std::complex<long double> > (u1);
  std::complex<long double> s1c = static_cast<std::complex<long double> > (s1);
  std::complex<long double> sc = static_cast<std::complex<long double> > (1.);
  std::complex<long double> Mppm=sqrt(1./8.*sc*tc*uc)*(64.*x/(sc*tc*uc)
  +16.*(x-4.*y*x)/(sc*tc*uc)*(s1c*C1(1.,y)+u1c*C1(u,y)+t1c*C1(t,y))
  -8.*(x-4.*y*x)/(sc*tc*uc)*(sc*tc*D(1.,t,y)+uc*sc*D(u,1.,y)+uc*tc*D(u,t,y)));
  return (Mppm*y*x);
}

long double LOBaur::MatrixElementFO(double xp){
  if( xp > std::pow(1.-x,2.)/(4.*x) )  // process is kinematically forbidden
    return 0.;
  SetMandelstam(xp);
  switch (_choice){
    case(0): {
        return x*( std::norm( Appp(yt)+Appp(yb) ) + std::norm( Appm(yt)+Appm(yb) ) +
    std::norm( Ampp(yt)+Ampp(yb) ) + std::norm( Ampm(yt)+Ampm(yb) ) );
      break;
    }
    case(1): {
        return x*( std::norm( Appp(yt) ) + std::norm( Appm(yt) ) +
    std::norm( Ampp(yt) ) + std::norm( Ampm(yt) ) );
      break;
    }
    case(2): {
      return x*( std::norm( Appp(yb) ) + std::norm( Appm(yb) ) +
    std::norm( Ampp(yb) ) + std::norm( Ampm(yb) ) );
      break;
    }
    case(3):{
      return x*real( std::conj( Appp(yb) )*Appp(yt)+std::conj( Appp(yt) )*Appp(yb) 
      + std::conj( Appm(yb) )*Appm(yt) +std::conj( Appm(yt) )*Appm(yb)+
    std::conj( Ampp(yb) )*Ampp(yt)+std::conj( Ampp(yt) )*Ampp(yb) 
    + std::conj( Ampm(yb) )*Ampm(yt)+std::conj( Ampm(yt) )*Ampm(yb));
   
    break;
      
    }
  }
}
long double LOBaur::sigmapartLO(double xp)
{
  const long double Gf = 0.00001166364;
  const double K = 9./(16.*4.)*_as;
  return (K *0.5 * 3./M_PIl * MatrixElementFO(xp));
}

