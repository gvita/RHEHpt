#include "complex_def.h"
#include "LOBaur.h"

//Functions B1, C, C1, D, E

std::complex<long double> LOBaur::B1(double v) {
  std::complex<long double> rad_v =  static_cast<std::complex<long double> > (1.-4.*(y-std::complex<long double>(0.,1.)*epsilon)*x/v);
  std::complex<long double> rad_y =  static_cast<std::complex<long double> > (1.-4.*(y-std::complex<long double>(0.,1.)*epsilon));
  std::complex<long double> z=0.5*(1.+sqrt(rad_v));
  std::complex<long double> zy=0.5*(1.+sqrt(rad_y));
  return (-sqrt(rad_v)*log_c(-z/(1.-z))+sqrt(rad_y)*log_c(-zy/(1.-zy)));
}
std::complex<long double> LOBaur::C(double v) {
  std::complex<long double> rad_v = static_cast<std::complex<long double> >  (1.-4.*(y-std::complex<long double>(0.,1.)*epsilon)*x/v);
  std::complex<long double> z=0.5*(1.+sqrt(rad_v));
  return (1./(2.*v)*pow(log_c(-z/(1.-z)),2));
}
std::complex<long double> LOBaur::C1(double v){
  return(v/(v-x)*C(v)-x/(v-x)*C(x));
}  
std::complex<long double> LOBaur::D(double v, double w){
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

std::complex<long double> LOBaur::E(double v, double w){
  return(v*C(v)+w*C(w)+(v-x)*C1(v)+(w-x)*C1(w)-v*w*D(v,w));
}

//Helicity Amplitudes

std::complex<long double> LOBaur::Appp(){
  std::complex<long double> tc= static_cast<std::complex<long double> > (t);
  std::complex<long double> uc = static_cast<std::complex<long double> > (u);
  std::complex<long double> t1c = static_cast<std::complex<long double> > (t1);
  std::complex<long double> u1c = static_cast<std::complex<long double> > (u1);
  std::complex<long double> s1c = static_cast<std::complex<long double> > (s1);
  std::complex<long double> sc = static_cast<std::complex<long double> > (1.);
  std::complex<long double> Mppp=sqrt(1./8.*sc*tc*uc)*(-64.*(1./(uc*tc)+1./(tc*t1c)+1./(uc*u1c))
  -64./sc*((2.*sc+tc)/(u1c*u1c)*B1(u)+(2.*sc+uc)/(t1c*t1c)*B1(t))
  -16.*(sc-4.*y*x)/(tc*uc*sc)*(s1c*C1(1.)+(uc-sc)*C1(t)+(tc-sc)*C1(u))
  -128.*y*x*(1./(tc*t1c)*C1(t)+1./(uc*u1c)*C1(u))+64.*y*x/sc*D(u,t)
  +8.*(sc-4.*y*x)/(sc*tc*uc)*(sc*tc*D(1.,t)+uc*sc*D(u,1.)-uc*tc*D(u,t))-32./sc/sc*E(u,t));
  return (Mppp*y*x);
  
  
  
  
}
std::complex<long double> LOBaur::Ampm(){
  std::complex<long double> tc= static_cast<std::complex<long double> > (t);
  std::complex<long double> uc = static_cast<std::complex<long double> > (u);
  std::complex<long double> t1c = static_cast<std::complex<long double> > (t1);
  std::complex<long double> u1c = static_cast<std::complex<long double> > (u1);
  std::complex<long double> s1c = static_cast<std::complex<long double> > (s1);
  std::complex<long double> sc = static_cast<std::complex<long double> > (1.);
  std::complex<long double> Mmpm=sqrt(1./8.*sc*uc*tc)*(-64.*(1./(uc*sc)+1./(sc*s1c)+1./(uc*u1c))
  -64./tc*((2.*tc+sc)/(u1c*u1c)*B1(u)+(2.*tc+uc)/(s1c*s1c)*B1(1.))
  -16.*(tc-4.*y*x)/(sc*uc*tc)*(t1c*C1(t)+(uc-tc)*C1(1.)+(sc-tc)*C1(u))
  -128.*y*x*(1./(sc*s1c)*C1(1.)+1./(uc*u1c)*C1(u))+64.*y*x/tc*D(u,1.)
  +8.*(tc-4.*y*x)/(tc*sc*uc)*(tc*sc*D(t,1.)+uc*tc*D(u,t)-uc*sc*D(u,1.))-32./tc/tc*E(u,1.));
  return (Mmpm*y*x);
  
  
  
  
}
std::complex<long double> LOBaur::Ampp(){
  std::complex<long double> tc= static_cast<std::complex<long double> > (t);
  std::complex<long double> uc = static_cast<std::complex<long double> > (u);
  std::complex<long double> t1c = static_cast<std::complex<long double> > (t1);
  std::complex<long double> u1c = static_cast<std::complex<long double> > (u1);
  std::complex<long double> s1c = static_cast<std::complex<long double> > (s1);
  std::complex<long double> sc = static_cast<std::complex<long double> > (1.);
  std::complex<long double> Mmpp=sqrt(1./8.*sc*tc*uc)*(-64.*(1./(sc*tc)+1./(tc*t1c)+1./(sc*s1c))
  -64./uc*((2.*uc+tc)/(s1c*s1c)*B1(1.)+(2.*uc+sc)/(t1c*t1c)*B1(t))
  -16.*(uc-4.*y*x)/(tc*sc*uc)*(u1c*C1(u)+(sc-uc)*C1(t)+(tc-uc)*C1(1.))
  -128.*y*x*(1./(tc*t1c)*C1(t)+1./(sc*s1c)*C1(1.))+64.*y*x/uc*D(1.,t)
  +8.*(uc-4.*y*x)/(uc*tc*sc)*(uc*tc*D(u,t)+sc*uc*D(1.,u)-sc*tc*D(1.,t))-32./uc/uc*E(1.,t));
  return (Mmpp*y*x);
  
  
  
  
}
std::complex<long double> LOBaur::Appm(){
  std::complex<long double> tc = static_cast<std::complex<long double> > (t);
  std::complex<long double> uc = static_cast<std::complex<long double> > (u);
  std::complex<long double> t1c = static_cast<std::complex<long double> > (t1);
  std::complex<long double> u1c = static_cast<std::complex<long double> > (u1);
  std::complex<long double> s1c = static_cast<std::complex<long double> > (s1);
  std::complex<long double> sc = static_cast<std::complex<long double> > (1.);
  std::complex<long double> Mppm=sqrt(1./8.*sc*tc*uc)*(64.*x/(sc*tc*uc)
  +16.*(x-4.*y*x)/(sc*tc*uc)*(s1c*C1(1.)+u1c*C1(u)+t1c*C1(t))
  -8.*(x-4.*y*x)/(sc*tc*uc)*(sc*tc*D(1.,t)+uc*sc*D(u,1.)+uc*tc*D(u,t)));
  return (Mppm*y*x);
}


long double LOBaur::MatrixElementFO( double xp ){
  SetMandelstam(xp);
  return x*( std::norm( Appp() ) + std::norm( Appm() ) +
    std::norm( Ampp() ) + std::norm( Ampm() ) );   
}










  
  
