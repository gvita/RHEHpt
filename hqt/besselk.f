c.....Function BK(n,z)
c.....BK(n,z) is the n-derivative of BesselK[nu,z]
c.....with respect to nu in nu=1

c     NOW b0->b0p INCLUDED

      function IK(m)
      implicit none
      external BK
      real *8 BK,IK,argum,dloqt,a_param,b0p
      integer m
      include 'scales.F'
      include 'const.F'
      common/a_param/a_param,b0p

      argum=b0p*qt/q
      dloqt=DLOG(a_param*qt/q)
      
      if (m.eq.1) then
         IK=-(2*b0p/q/qt)*BK(0,argum)
      elseif (m.eq.2) then
         IK=-(4*b0p/q/qt)*(BK(1,argum)-BK(0,argum)*DLOqt)
      elseif (m.eq.3) then
         IK=(b0p/q/qt)*(-6*BK(2,argum)+12*BK(1,argum)*DLOqt+
     /                 BK(0,argum)*(pi**2-6*(DLOqt**2)))
      elseif (m.eq.4) then
         IK=-(4*b0p/q/qt)*(2*BK(3,argum)-6*DLOqt*BK(2,argum)+
     /                    BK(1,argum)*(6*(DLOqt)**2-pi**2)+
     /                    BK(0,argum)*(pi**2*DLOqt-
     /                                 2*(DLOqt)**3-4*Z3))
      endif


      return
      end



      function BK(n,z)
      implicit none
      external fb
      real *8 bk,fb,errest,z,zz,max,adpint
      integer n,nn,ifail
      common/nuorder/nn
      common/zz/zz
      nn=n
      zz=z
      max=10d0
      bk=adpint(fb,0d0,max,1d-10,1d-5,errest,ifail)
      return
      end
      
      
      function fb(t)
      implicit none
      integer nn,nu
      real *8 fb,t,zz
      common/nuorder/nn
      common/zz/zz
      nu=1
      if(nn.eq.0) then
         fb=dexp(-zz*dcosh(t))*dcosh(nu*t)
      elseif(nn.eq.1) then
         fb=dexp(-zz*dcosh(t))*t*dsinh(nu*t)
      elseif(nn.eq.2) then
         fb=dexp(-zz*dcosh(t))*t*t*dcosh(nu*t)
      elseif(nn.eq.3) then
         fb=dexp(-zz*dcosh(t))*t*t*t*dsinh(nu*t)
      endif      
      return
      end
      
