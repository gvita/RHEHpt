
c***************************************************************************

      subroutine higgsqty(ris,y1,y2,ord)
c       
c     Written by M.Grazzini  
c     Returns dsigma/dqt/dy if y1=y2
c     and dsigma/dqt integrated between y1 and y2 if y1<y2  
c     Based on C.Glosser, C. Schmidt, JHEP 0212 (2002) 016
c
c     q <-> xmh, q2 <-> xmh^2
c     QQ2=Q2 in Glosser-Schmidt
c     qt Higgs transverse momentum (GeV)
c     y1,y2 minimum and maximum rapidities   
c
c     ord order of the calculation (1) -> LO (2) -> NLO
c
c     the normalization is sigma0 and is given through a common
c

      implicit none
      real *8 xmuf2
      real *8 pi,alphas
      real *8 y1,y2
      real *8 ris,av3,d3
      real *8 ALLGG,ALLqqb,ALLqq,ALLqg,ALLqqpb,ALLqqp,abisq
      real *8 delta,distr,distrcross
      real *8 tm,ymax
      real *8 xmb,xmt,sigma0,xlam,aass,alfa,yy1,yy2,amt
      real *8 xr,tmpx,tmp,tmp1,tmp2,beta0
      real *8 rdelta,rdistr,rdistrcross
      integer nloop,nf,isetproton,ih1,ih2
      integer i,k,ord,iord,inorm,ndim,iwseed,n3,ncl3
      EXTERNAL ALLGG,ALLqqb,ALLqq,ALLqg,ALLqqpb,ALLqqp,abisq,
     #delta,distr,distrcross


      common/nf/nf
      common/fixvar/xmuf2
      common/isetproton/isetproton
      common/pdf/ih1,ih2
      common/nloop/nloop
      common/internal/yy1,yy2
      common/tm/tm
      common/beta0/beta0
      common/pi/pi
      common/aass/aass
      common/iord/iord
      common/sigma0/sigma0
      common/amt/amt
      common/xmb/xmb

      include 'scales.F'


      PI=3.14159265358979312D0
      beta0=(33d0-2*nf)/6d0

      xmuf2=muf2
      iord=ord-1



      xmt=amt

c.....Transverse mass

      tm=dsqrt(q2+qt**2)

c.....Check y limits: y2 cannot be smaller than y1

      if(y1.gt.y2) then
      write(*,*)'Wrong y1,y2 !'
      stop
      endif


c.....Kinematical limits: qt

      xr=(1-xtau)**2-4*xtau*(qt/q)**2


      if(xr.lt.0) then
       ris=0
       return
      endif


c.....Kinematical limits: y


      tmpx=(q2+shad)/sroot/tm

      ymax=dlog((tmpx+dsqrt(tmpx**2-4))/2)


c     Always assume y1 not larger than y2

      if(y1.gt.ymax.or.y2.lt.(-ymax)) then
       ris=0
       return
      elseif(y1.lt.-ymax.and.y2.gt.ymax) then
       yy1=-ymax
       yy2=ymax     
      elseif(y1.lt.-ymax.and.y2.gt.-ymax.and.y2.lt.ymax) then
       yy1=-ymax
       yy2=y2
      elseif(y1.gt.-ymax.and.y1.lt.ymax.and.y2.gt.ymax) then
       yy1=y1
       yy2=ymax
      elseif(y1.gt.-ymax.and.y1.lt.ymax.and.y2.lt.ymax) then
       yy1=y1
       yy2=y2
      endif

c      write(*,*)yy1,yy2


      ndim=3
      

c     Number of iterations and calls to vegas

      if(qt/q.le.0.04d0) then
      ncl3=80000
      n3=15
      else
      ncl3=40000
      n3=10
      endif

c     If Rapidity integral is symmetric -> compute only first contribution

      if(iord.eq.1) then
       call xinteg(distr,ndim,n3,ncl3,av3,d3)      
       rdistr=aass/2*av3
       if(yy2.eq.-yy1) then
        rdistrcross=rdistr
        else
        call xinteg(distrcross,ndim,n3,ncl3,av3,d3)
        rdistrcross=aass/2*av3
       endif
      elseif(iord.eq.0) then
       rdistr=0d0
       rdistrcross=0d0      
      endif


c     Number of iterations
      n3=10

c     Number of calls to vegas
      ncl3=10000

      iwseed=1
      ndim=2

      call xinteg(delta,ndim,n3,ncl3,av3,d3)
      write(*,*) 'int=',av3,' +- ',d3

      rdelta=av3


      ris=sigma0*(aass/2d0)*(rdelta+rdistr+rdistrcross)

        
      return
      end 
      
   

c**********************************************************
      function delta(xx)
      implicit none
      real *8 delta,gg0,qg0,gq0,qqb0
      real *8 lo
      real *8 xx(1:2)
      real *8 sh,th,uh,s,x1,x2,x,t,jj,tiny
      real *8 y,yh,y1,y2,tm
      real *8 x2min,jac,tmp,esp
      real *8 allgg,allqg,allgq,allqqb
      real *8 aass
      real *8 beta0,de,uu,ddgg,Li2
      real *8 pi
      real *8 V1,V2,V3,V1c,V2c,V3c,W1,W2,W3
      real *8 ddqg,ddgq,ddqqb
      integer j,nf,iord
      external gg0,qg0,gq0,qqb0,Li2    

      common/internal/y1,y2
      common/tm/tm
      common/beta0/beta0
      common/nf/nf
      common/rap/yh
      common/pi/pi 
      common/aass/aass
      common/iord/iord

      include 'scales.F'


      tiny=1d-8
     
      

      jac=1d0

c     x and y have to vary between 0 and 1

      t=xx(1)

      if(y1.ne.y2) then
       yh=y1+(y2-y1)*xx(2)
       jac=jac*(y2-y1)
      else
       yh=y1
      endif
 
      x2min=(q2-sroot*tm*dexp(-yh))/(sroot*tm*dexp(yh)-shad)
      if(x2min.gt.1d0.or.x2min.lt.0d0) then
       write(*,*)'ERROR IN x2min !'
       stop
      endif

c     First possibility
c     esp=0.3d0
c      x2=x2min**(t**esp)

c     Second possibility: better small qt behaviour
c     The function is now peaked at central values of
c     the integration variable x

      esp=8d0

      x2=dexp((1-t**esp)*dlog(x2min))
      jac=jac*x2*dlog(x2min)*esp*t**(esp-1d0)
 
c      x2=x2min**t
c      jac=x2*dlog(x2min)

      x1=(sroot*tm*dexp(yh)*x2-q2)/(x2*shad-sroot*tm*dexp(-yh))

c     Impose x1,x2 < 1

      if(x1.gt.1-tiny.or.x2.gt.1-tiny) then
      delta=0d0
      return
      endif

       sh=x1*x2*shad
       th=q2-sroot*x2*tm*dexp(yh)
       uh=q2-sroot*x1*tm*dexp(-yh)


       jj=abs(x2*shad-sroot*tm*dexp(-yh))

c      Leading order contribution
     

       lo=0d0
     #   +gg0(x1,x2)*allgg(x1,x2)
     #   +qg0(x1,x2)*allqg(x1,x2)
     #   +gq0(x1,x2)*allgq(x1,x2)
     #   +qqb0(x1,x2)*allqqb(x1,x2)    



c      NLO delta terms

 
c     gq coefficients

      V1=0.5d0*((dlog(uh/th))**2+(dlog(-sh/uh))**2-(dlog(-sh/th))**2)
     #+dlog(sh/q2)*dlog(-th/q2)-dlog(sh/q2)*dlog(-uh/q2)
     #-dlog(-th/q2)*dlog(-uh/q2)+2*Li2(q2/(q2-uh))
     #+(dlog(q2/(q2-uh)))**2+pi**2


      V2=(dlog(q2/sh))**2+(dlog(q2/(q2-th)))**2
     #-2*dlog(sh/q2)*dlog(-th/q2)+2*Li2(1-q2/sh)
     #+2*Li2(q2/(q2-th))-3.5d0-2*pi**2/3d0


      V3=beta0*(2*dlog(-mur2/uh)+dlog(-mur2/th))+67d0/3-10d0*nf/9d0

c     Coefficient of delta(QQ2) in gq


      ddgq=((11d0+3*V1+4d0*V2/3+V3)*gq0(x1,x2)
     #+20d0/9*(sh**2+th**2+uh**2-uh*q2)/(-uh))*allgq(x1,x2)


c     qg coefficients

      V1c=0.5d0*((dlog(th/uh))**2+(dlog(-sh/th))**2-(dlog(-sh/uh))**2)
     #+dlog(sh/q2)*dlog(-uh/q2)-dlog(sh/q2)*dlog(-th/q2)
     #-dlog(-uh/q2)*dlog(-th/q2)+2*Li2(q2/(q2-th))
     #+(dlog(q2/(q2-th)))**2+pi**2


      V2c=(dlog(q2/sh))**2+(dlog(q2/(q2-uh)))**2
     #-2*dlog(sh/q2)*dlog(-uh/q2)+2*Li2(1-q2/sh)
     #+2*Li2(q2/(q2-uh))-3.5d0-2*pi**2/3d0


      V3c=beta0*(2*dlog(-mur2/th)+dlog(-mur2/uh))+67d0/3-10d0*nf/9d0


c     Coefficient of delta(QQ2) in gq


      ddqg=((11d0+3*V1c+4d0*V2c/3+V3c)*qg0(x1,x2)
     #+20d0/9*(sh**2+th**2+uh**2-th*q2)/(-th))*allqg(x1,x2)


c     qqb coefficients

      W1=dlog(-uh/q2)*dlog(-th/q2)-dlog(sh/q2)*dlog(-uh/q2)
     #-dlog(sh/q2)*dlog(-th/q2)+2*Li2(1-q2/sh)
     #+(dlog(q2/sh))**2-0.5d0*(dlog(uh/th))**2-5*Pi**2/3d0


      W2=1.5d0*dlog(sh**2/th/uh)+(dlog(uh/th))**2
     #-2*dlog(-uh/q2)*dlog(-th/q2)+(dlog(q2/(q2-uh)))**2
     #+(dlog(q2/(q2-th)))**2
     #+2*Li2(q2/(q2-uh))+2*Li2(q2/(q2-th))-7+2*Pi**2


      W3=beta0/2*(4*dlog(mur2/sh)+dlog(-mur2/uh)+dlog(-mur2/th))
     #+(67d0/2-5d0*nf/3d0)


c     Coefficient of delta(QQ2) in qqb

      ddqqb=((11d0+3*W1+4d0*W2/3+W3)*qqb0(x1,x2)
     #+160d0/27*(th**2+uh**2+sh**2-sh*q2)/sh)


c     Coefficient of delta(QQ2) in gg


       uu=0.5d0*(dlog(uh/th))**2+pi**2/3
     #    -dlog(sh/q2)*dlog(-th/q2)-dlog(sh/q2)*dlog(-uh/q2)
     #    -dlog(-uh/q2)*dlog(-th/q2)
     #    +(dlog(q2/sh))**2+(dlog(q2/(q2-th)))**2
     #    +(dlog(q2/(q2-uh)))**2+2*Li2(1d0-q2/sh)
     #    +2*Li2(q2/(q2-th))+2*Li2(q2/(q2-uh))



       de=1.5d0*beta0*(dlog(-mur2/th)+dlog(-mur2/uh))
     #   +67d0/6-5d0/9*nf


       ddgg=((11d0+de+3*uu)*gg0(x1,x2)
     #  +(3-nf)*(q2**2/sh+q2**2/th+q2**2/uh+q2))*allgg(x1,x2)



c     Total result

       delta=lo+iord*aass/2*(ddgg+ddgq+ddqg+ddqqb)


       delta=-2*qt*delta*jac/sh/jj



      return
      end


c*********************************************************
c     Singular terms
c*********************************************************

      function distr(zz)
c
c     Here use parametrization of Eq. (B.1)
c
      implicit none
      real *8 zz(1:3)
      real *8 y1,y2,yh
      real *8 distr,distr1,distr2,z1,z2,tm,z2min
      real *8 sh,uh,th,x1,x2,jac,jacy,tiny
      real *8 x10,x20,dcut,la,lb,QQ2
      real *8 tmp,a1,a2,a10,a20
      real *8 b1,b2,b10,b20,c1,c2,d10,d20,m10,m20,r1,r2
      real *8 Pgg,Pqq,Pqg,Pgq
      real *8 gg0,qg0,gq0,qqb0
      real *8 allgg,allgq,allqg,allqqb,allqq,allqqpb,allqqp
      real *8 beta0,u,t,s,temp0
      real *8 pre1,pre2,pre10,pre20
      real *8 big1,big2,big3,big4,big5
      real *8 REGgg,REGgq,REGqg,REGqqb,REGqqpb,REGqq
      integer nf


      common/internal/y1,y2
      common/tm/tm
      common/nf/nf
      common/rap/yh
      common/beta0/beta0
      common/coef/u,t,s
      common/big/big1,big2,big3,big4,big5
      common/regu/uh,th,sh,REGgg,REGgq,REGqg,REGqqb,REGqqpb,REGqq


      include 'scales.F'

      tiny=1d-7
      

c     zz(1) <-> z1      
c     zz(2) <-> z2
c     zz(3) <-> yh

      jacy=1d0
      if(y1.ne.y2) then
       yh=y1+(y2-y1)*zz(3)
       jacy=jacy*(y2-y1)
      elseif(y1.eq.y2) then
       yh=y1
      endif

 
      x10=tm/sroot*dexp(yh)
      x20=tm/sroot*dexp(-yh)
c     d=qt/(tm+qt)
      dcut=x10*(qt/tm)**2/(1-x10*(1-(qt/tm)**2))


c*********************************************************
c     First double integral zbp->z1,za->z2
c     Set both integrals from 0 to 1

      jac=1d0
      z1=x20+(1-dcut-x20)*zz(1)
      jac=jac*(1-dcut-x20)

      lb=qt**2/tm**2*z1/(1-z1)

      z2=x10*(1+lb)+(1-x10*(1+lb))*zz(2)
      jac=jac*(1-x10*(1+lb))     

      x1=x10/z2*(1+lb)
      x2=x20/z1


      sh=tm**2/z1/z2*(1+lb)
      th=-tm**2/z1*(1-z1)*(1+lb)
      uh=q2-tm**2/z2*(1+lb)      
      QQ2=tm**2/z1/z2*(1-z2)*(1-z1)*(1+lb)

      if(QQ2.lt.0) then
       write(*,*)'QQ2 < 0 !'
       stop
      endif



c     Phase space prefactor

      pre1=tm**2*(1+lb)/(z1*z2)**2

c     Impose x1,x2 < 1

      if(x1.gt.1-tiny.or.x2.gt.1-tiny) then
       a1=0d0
       b1=0d0
       c1=0d0
       r1=0d0
      else

      u=uh
      t=th
      s=sh

      a1=0
      b1=0
      c1=0
      r1=0

c     gg contribution


      call coeff   
 

c     a):   (Log(1-z2)/(1-z2)+ terms

      a1=a1+(1/(-th)*Pgg(z2)*gg0(z2*x1,x2)+(-z2/th)*big1)
     #   *allgg(x1,x2)/sh

c     b):    1/(1-z2)+ term

      b1=b1+(1d0/th*Pgg(z2)*dlog(-muf2*z2/th)*gg0(z2*x1,x2)
     #+z2/th*big1*dlog((QQ2+qt**2)*z2/(-th))+z2/th*big2)
     #*allgg(x1,x2)/sh


c     c):  Regular (non delta or plus distributions) terms: 
c          term in big3 symmetrized: half here and half in distrcross

      
      c1=c1+(1/(-th)*(-2*nf*Pqg(z2)*dlog(muf2/QQ2)+2*nf*z2*(1-z2))
     #*qg0(z2*x1,x2)+0.5d0*big3)*allgg(x1,x2)/sh

      


c     add gq contribution: term in big5 symmetrized
       
      a1=a1+(1/(-th)*Pgg(z2)*gq0(z2*x1,x2)+(-z2/th)*big4)
     #   *allgq(x1,x2)/sh


      b1=b1+(1d0/th*Pgg(z2)*dlog(-muf2*z2/th)*gq0(z2*x1,x2)
     #+z2/th*big4*dlog((QQ2+qt**2)*z2/(-th)))*allgq(x1,x2)/sh


      c1=c1+(1/(-th)*(-Pqg(z2)*dlog(muf2/QQ2)+z2*(1-z2))
     #*qqb0(z2*x1,x2)+0.5d0*big5)*allgq(x1,x2)/sh


c     add qg contribution: term in big5 symmetrized


      u=th
      t=uh
      s=sh

      call coeff

   
      a1=a1-1d0/th*Pqq(z2)*qg0(z2*x1,x2)*allqg(x1,x2)/sh


      b1=b1+(1d0/th*Pqq(z2)*dlog(-muf2*z2/th)*qg0(z2*x1,x2)
     #+z2/th*8d0/3*(uh**2+sh**2)/(-th))*allqg(x1,x2)/sh


      c1=c1+(-1/th*(4d0/3*(1-z2)*qg0(z2*x1,x2)
     #+(-Pgq(z2)*dlog(muf2/QQ2)+4d0/3*z2)*gg0(z2*x1,x2))
     #+0.5d0*big5)*allqg(x1,x2)/sh


c     SWITCH OFF gg,qg,gq      

c       a1=0
c       b1=0
c       c1=0
c       r1=0


c     add qqb contribution: last regular term symmetrized


      a1=a1+1/(-th)*(Pqq(z2)*qqb0(z2*x1,x2)-z2*16d0/27
     #*(th**2+uh**2+(QQ2-th)**2+(QQ2-uh)**2)/sh)*allqqb(x1,x2)/sh



      b1=b1+(1/th*Pqq(z2)*dlog(-muf2*z2/th)*qqb0(z2*x1,x2)
     #-z2/th*dlog((QQ2+qt**2)*z2/(-th))*16d0/27
     #*(th**2+uh**2+(QQ2-th)**2+(QQ2-uh)**2)/sh
     #+z2/th*16d0/9*beta0*(uh**2+th**2)/sh
     #)*allqqb(x1,x2)/sh



      c1=c1+(
     #-1/th*(4d0/3*(1-z2)*qqb0(z2*x1,x2)+(-Pgq(z2)*dlog(muf2/QQ2)
     #+4d0/3*z2)*gq0(z2*x1,x2))+16d0/9*((sh-QQ2)**2+(uh+th-2*QQ2)**2)/sh
     #*dlog(qt**2/(qt**2+QQ2))/qt**2
     #)*allqqb(x1,x2)/sh



c    Add qQb, qq and qQ contributions


      c1=c1+(-1/th*(-Pgq(z2)*dlog(muf2/QQ2)+4d0/3*z2)*gq0(z2*x1,x2)
     #+16d0/9*((sh-QQ2)**2+(uh+th-2*QQ2)**2)/sh
     #*dlog(qt**2/(qt**2+QQ2))/qt**2)
     #*(allqq(x1,x2)+allqqpb(x1,x2)+allqqp(x1,x2))/sh


c     d): Regular terms in gg,qg,gq
c         half here and half in distrcross

      call REG

      r1=0.5d0*(0d0
     #+REGgg*allgg(x1,x2)
     #+REGgq*allgq(x1,x2)
     #+REGqg*allqg(x1,x2)
     #+REGqqb*allqqb(x1,x2)
     #+REGqqpb*(allqqpb(x1,x2)+allqqp(x1,x2))
     #+REGqq*allqq(x1,x2)
     #)/sh



      endif

c*********************************************************
c     First double integral with zbp->z1,za->z2=1

      x1=x10*(1+lb)
      x2=x20/z1    

      sh=tm**2/z1*(1+lb)
      th=-tm**2/z1*(1-z1)*(1+lb)
      uh=q2-tm**2*(1+lb) 
      QQ2=0d0



      pre10=tm**2*(1+lb)/(z1)**2

      if(x1.gt.1-tiny.or.x2.gt.1-tiny) then
       a10=0d0
       b10=0d0
       d10=0d0
       m10=0d0
      else


 
      u=uh
      t=th
      s=sh

      a10=0
      b10=0
      d10=0



c     gg contribution


      call coeff   

      
      a10=a10+(1/(-th)*pgg(1d0)*gg0(x1,x2)+(-1/th)*big1)
     #   *allgg(x1,x2)/sh


      b10=b10+(6d0/th*dlog(-muf2/th)*gg0(x1,x2)
     #+big1/th*dlog((QQ2+qt**2)/(-th))+big2/th)
     #*allgg(x1,x2)/sh


c     Contribution from delta(1-z) term in Pgg

      d10=d10+1/th*beta0*dlog(muf2/(-th))*gg0(x1,x2)*allgg(x1,x2)/sh




c     add gq contribution


      a10=a10+(1/(-th)*Pgg(1d0)*gq0(x1,x2)+(-1d0/th)*big4)
     #   *allgq(x1,x2)/sh



      
      b10=b10+(1d0/th*Pgg(1d0)*dlog(-muf2/th)*gq0(x1,x2)
     #+1/th*big4*dlog((QQ2+qt**2)/(-th)))*allgq(x1,x2)/sh

c    add contribution from delta(1-z) term in Pgg


      d10=d10+1/th*beta0*dlog(muf2/(-th))*gq0(x1,x2)
     #*allgq(x1,x2)/sh
   


c      add qg contribution


      u=th
      t=uh
      s=sh

      call coeff

    
     
      a10=a10-1d0/th*Pqq(1d0)*qg0(x1,x2)*allqg(x1,x2)/sh


      b10=b10+(1d0/th*Pqq(1d0)*dlog(-muf2/th)*qg0(x1,x2)
     #+1/th*8d0/3*(uh**2+sh**2)/(-th))*allqg(x1,x2)/sh


c     add contribution from delta(1-z) term in Pqq


      d10=d10+1/th*2*dlog(muf2/(-th))*qg0(x1,x2)
     #*allqg(x1,x2)/sh
  

c     SWITCH OFF gg,qg,gq


c      a10=0
c      b10=0
c      d10=0



c     add qqb contribution

      a10=a10+1/(-th)*(Pqq(1d0)*qqb0(x1,x2)-16d0/27
     #*(th**2+uh**2+(QQ2-th)**2+(QQ2-uh)**2)/sh)*allqqb(x1,x2)/sh



      b10=b10+(1/th*Pqq(1d0)*dlog(-muf2/th)*qqb0(x1,x2)
     #-1d0/th*dlog((QQ2+qt**2)/(-th))*16d0/27
     #*(th**2+uh**2+(QQ2-th)**2+(QQ2-uh)**2)/sh
     #+1d0/th*16d0/9*beta0*(uh**2+th**2)/sh
     #)*allqqb(x1,x2)/sh


c     add contribution from delta(1-z) term in Pqq


      d10=d10+1/th*2*dlog(muf2/(-th))*qqb0(x1,x2)
     #*allqqb(x1,x2)/sh
  


c     Mismatch from plus distributions in z2: integral from 0 to z2min

c     \int_0^z2min 1/(1-z2) =  -log(1-z2min)      
c     \int_0^z2min Log(1-z2)/(1-z2) ->-1/2 Log[1-z2min]^2

      z2min=x10*(1+lb)


      m10=-0.5d0*a10*(dlog(1-z2min))**2-dlog(1-z2min)*b10      

c     Multiply by jacobian from z1 integration

      m10=m10*pre10*(1-dcut-x20)

      endif

c*********************************************************
c     Result

      a1=dlog(1-z2)/(1-z2)*(pre1*a1-pre10*a10)*jac

      b1=(pre1*b1-pre10*b10)/(1-z2)*jac

      c1=pre1*c1*jac

c     For d10 only jacobian from z1 integration

      d10=pre10*d10*(1-dcut-x20)

      r1=r1*pre1*jac

      distr1=(a1+b1+c1+d10-m10+r1)/shad


      distr=2*qt*jacy*distr1

      return
      end


c*********************************************************

      function distrcross(zz)
c
c     Here use parametrization of Eq. (B.1)
c     but obtained with (a<->b),(u<->t)
c
      implicit none
      real *8 zz(1:3)
      real *8 y1,y2,yh
      real *8 distrcross,dcross1,dcross2,z1,z2,tm,z2min
      real *8 sh,uh,th,x1,x2,jac,jacy,tiny
      real *8 x10,x20,d,la,lb,QQ2
      real *8 tmp,a1,a2,a10,a20
      real *8 b1,b2,b10,b20,c1,c2,d10,d20,m10,m20,r1,r2
      real *8 Pgg,Pqg,Pqq,Pgq
      real *8 gg0,qg0,gq0,qqb0
      real *8 allgg,allgq,allqg,allqqb,allqq,allqqpb,allqqp
      real *8 beta0,u,t,s
      real *8 pre1,pre2,pre10,pre20
      real *8 big1,big2,big3,big4,big5,xxa,xxb
      real *8 REGgg,REGgq,REGqg,REGqqb,REGqqpb,REGqq
      integer nf

      common/internal/y1,y2
      common/tm/tm
      common/nf/nf
      common/rap/yh
      common/beta0/beta0
      common/coef/u,t,s
      common/big/big1,big2,big3,big4,big5
      common/regu/uh,th,sh,REGgg,REGgq,REGqg,REGqqb,REGqqpb,REGqq

      include 'scales.F'

      tiny=1d-6

      muf2=muf**2

c     zz(1) <-> z1      
c     zz(2) <-> z2
c     zz(3) <-> yh

      jacy=1d0
      if(y1.ne.y2) then
       yh=y1+(y2-y1)*zz(3)
       jacy=jacy*(y2-y1)
      elseif(y1.eq.y2) then
       yh=y1
      endif

 
      x10=tm/sroot*dexp(yh)
      x20=tm/sroot*dexp(-yh)
c      d=qt/(tm+qt)
      d=x20*(qt/tm)**2/(1-x20*(1-(qt/tm)**2))


c*********************************************************
c     First double integral zap->z1,zb->z2
c     Set both integrals from 0 to 1

      jac=1d0
      z1=x10+(1-d-x10)*zz(1)
      jac=jac*(1-d-x10)

      la=qt**2/tm**2*z1/(1-z1)

      z2=x20*(1+la)+(1-x20*(1+la))*zz(2)
      jac=jac*(1-x20*(1+la))     

            

      x1=x10/z1
      x2=x20/z2*(1+la)

      sh=tm**2/z1/z2*(1+la)
      uh=-tm**2/z1*(1-z1)*(1+la)
      th=q2-tm**2/z2*(1+la)      
      QQ2=tm**2/z1/z2*(1-z2)*(1-z1)*(1+la)
     
      if(QQ2.lt.0) then
       write(*,*)'QQ2 < 0 !'
       stop
      endif



c     Phase space prefactor

      pre1=tm**2*(1+la)/(z1*z2)**2

c     Impose x1,x2 < 1

      if(x1.gt.1-tiny.or.x2.gt.1-tiny) then
       a1=0d0
       b1=0d0
       c1=0d0
       r1=0d0
      else


      u=uh
      t=th
      s=sh

      a1=0
      b1=0
      c1=0
      r1=0

c     gg contribution


      call coeff   
      

c     a):   (Log(1-z2)/(1-z2)+ terms

      a1=a1+(1/(-uh)*pgg(z2)*gg0(x1,x2*z2)+(-z2/uh)*big1)
     #   *allgg(x1,x2)/sh

c     b):    1/(1-z2)+ term

      b1=b1+(1/uh*Pgg(z2)*dlog(-muf2*z2/uh)*gg0(x1,x2*z2)
     #+z2/uh*big1*dlog((QQ2+qt**2)*z2/(-uh))+z2/uh*big2)
     #*allgg(x1,x2)/sh


c     c): Regular term

      c1=c1+(1/(-uh)*(-2*nf*Pqg(z2)*dlog(muf2/QQ2)+2*nf*z2*(1-z2))
     #*gq0(x1,x2*z2)+0.5d0*big3)*allgg(x1,x2)/sh



c     add gq contribution: term in big5 symmetrized
       
      a1=a1-1d0/uh*Pqq(z2)*gq0(x1,z2*x2)*allgq(x1,x2)/sh

      b1=b1+(1d0/uh*Pqq(z2)*dlog(-muf2*z2/uh)*gq0(x1,z2*x2)
     #+z2/uh*8d0/3*(th**2+sh**2)/(-uh))*allgq(x1,x2)/sh


      c1=c1+(-1/uh*(4d0/3*(1-z2)*gq0(x1,z2*x2)
     #+(-Pgq(z2)*dlog(muf2/QQ2)+4d0/3*z2)*gg0(x1,z2*x2))
     #+0.5d0*big5)*allgq(x1,x2)/sh


c     add qg contribution: term in big5 symmetrized


      u=th
      t=uh
      s=sh

      call coeff

      a1=a1+(1/(-uh)*Pgg(z2)*qg0(x1,z2*x2)+(-z2/uh)*big4)
     #   *allqg(x1,x2)/sh


      b1=b1+(1d0/uh*Pgg(z2)*dlog(-muf2*z2/uh)*qg0(x1,z2*x2)
     #+z2/uh*big4*dlog((QQ2+qt**2)*z2/(-uh)))*allqg(x1,x2)/sh


      c1=c1+(1/(-uh)*(-Pqg(z2)*dlog(muf2/QQ2)+z2*(1-z2))
     #*qqb0(x1,z2*x2)+0.5d0*big5)*allqg(x1,x2)/sh


c     SWITCH OFF gg,qg,gq      

c      a1=0
c      b1=0
c      c1=0
c      r1=0


c     add qqb contribution: last regular term symmetrized


      a1=a1+1/(-uh)*(Pqq(z2)*qqb0(x1,z2*x2)-z2*16d0/27
     #*(th**2+uh**2+(QQ2-th)**2+(QQ2-uh)**2)/sh)*allqqb(x1,x2)/sh



      b1=b1+(1/uh*Pqq(z2)*dlog(-muf2*z2/uh)*qqb0(x1,z2*x2)
     #-z2/uh*dlog((QQ2+qt**2)*z2/(-uh))*16d0/27
     #*(th**2+uh**2+(QQ2-th)**2+(QQ2-uh)**2)/sh
     #+z2/uh*16d0/9*beta0*(uh**2+th**2)/sh
     #)*allqqb(x1,x2)/sh



      c1=c1+(
     #-1/uh*(4d0/3*(1-z2)*qqb0(x1,z2*x2)+(-Pgq(z2)*dlog(muf2/QQ2)
     #+4d0/3*z2)*qg0(x1,z2*x2))+16d0/9*((sh-QQ2)**2+(uh+th-2*QQ2)**2)/sh
     #*dlog(qt**2/(qt**2+QQ2))/qt**2
     #)*allqqb(x1,x2)/sh


c    Add qQb, qq and qQ contributions


      c1=c1+(-1/uh*(-Pgq(z2)*dlog(muf2/QQ2)+4d0/3*z2)*qg0(x1,z2*x2)
     #+16d0/9*((sh-QQ2)**2+(uh+th-2*QQ2)**2)/sh
     #*dlog(qt**2/(qt**2+QQ2))/qt**2)
     #*(allqq(x1,x2)+allqqpb(x1,x2)+allqqp(x1,x2))/sh



c     d): Regular terms in gg,qg,gq
c         half here and half in distr

      call REG

      r1=0.5d0*(0d0
     #+REGgg*allgg(x1,x2)
     #+REGgq*allgq(x1,x2)
     #+REGqg*allqg(x1,x2)
     #+REGqqb*allqqb(x1,x2)
     #+REGqqpb*(allqqpb(x1,x2)+allqqp(x1,x2))
     #+REGqq*allqq(x1,x2)
     #)/sh

      endif


c*********************************************************
c     First double integral with zap->z1,zb->z2=1

      x1=x10/z1
      x2=x20*(1+la)
      sh=tm**2/z1*(1+la)
      uh=-tm**2/z1*(1-z1)*(1+la)
      th=q2-tm**2*(1+la)      
     
      QQ2=0d0


      pre10=tm**2*(1+la)/(z1)**2

      if(x1.gt.1-tiny.or.x2.gt.1-tiny) then
       a10=0d0
       b10=0d0
       d10=0d0
       m10=0d0
      else


      u=uh
      t=th
      s=sh

      a10=0
      b10=0
      d10=0


c     gg contribution


      call coeff   
      

      a10=a10+(1/(-uh)*pgg(1d0)*gg0(x1,x2)+(-1/uh)*big1)
     #   *allgg(x1,x2)/sh

      b10=b10+(1/uh*Pgg(1d0)*dlog(-muf2/uh)*gg0(x1,x2)
     #+1/uh*big1*dlog((QQ2+qt**2)/(-uh))+1/uh*big2)
     #*allgg(x1,x2)/sh

c     Contribution from delta(1-z) term in Pgg

      d10=d10+1/uh*beta0*dlog(muf2/(-uh))*gg0(x1,x2)*allgg(x1,x2)/sh



c     add gq contribution


      a10=a10-1d0/uh*Pqq(1d0)*gq0(x1,x2)*allgq(x1,x2)/sh

      b10=b10+(1d0/uh*Pqq(1d0)*dlog(-muf2/uh)*gq0(x1,x2)
     #+1d0/uh*8d0/3*(th**2+sh**2)/(-uh))*allgq(x1,x2)/sh


c    add contribution from delta(1-z) term in Pqq


      d10=d10+1/uh*2*dlog(muf2/(-uh))*gq0(x1,x2)
     #*allgq(x1,x2)/sh


c      add qg contribution


      u=th
      t=uh
      s=sh

      call coeff



      a10=a10+(1/(-uh)*Pgg(1d0)*qg0(x1,x2)+(-1d0/uh)*big4)
     #   *allqg(x1,x2)/sh


      b10=b10+(1d0/uh*Pgg(1d0)*dlog(-muf2/uh)*qg0(x1,x2)
     #+1d0/uh*big4*dlog((QQ2+qt**2)/(-uh)))*allqg(x1,x2)/sh

c    add contribution from delta(1-z) term in Pgg


      d10=d10+1/uh*beta0*dlog(muf2/(-uh))*qg0(x1,x2)
     #*allqg(x1,x2)/sh
   

c     SWITCH OFF gg,gq,qg


c      a10=0
c      b10=0
c      d10=0


c     add qqb contribution


      a10=a10+1/(-uh)*(Pqq(1d0)*qqb0(x1,x2)-16d0/27
     #*(th**2+uh**2+(QQ2-th)**2+(QQ2-uh)**2)/sh)*allqqb(x1,x2)/sh



      b10=b10+(1/uh*Pqq(1d0)*dlog(-muf2/uh)*qqb0(x1,x2)
     #-1d0/uh*dlog((QQ2+qt**2)/(-uh))*16d0/27
     #*(th**2+uh**2+(QQ2-th)**2+(QQ2-uh)**2)/sh
     #+1d0/uh*16d0/9*beta0*(uh**2+th**2)/sh
     #)*allqqb(x1,x2)/sh


c     add contribution from delta(1-z) term in Pqq


      d10=d10+1/uh*2*dlog(muf2/(-uh))*qqb0(x1,x2)
     #*allqqb(x1,x2)/sh
  


c     Mismatch from plus distributions in z2: integral from 0 to z2min

c     \int_0^z2min 1/(1-z2) =  -log(1-z2min)      
c     \int_0^z2min Log(1-z2)/(1-z2) ->-1/2 Log[1-z2min]^2

      z2min=x20*(1+la)


      m10=-0.5d0*(dlog(1-z2min))**2*a10-dlog(1-z2min)*b10      

c     Multiply by jacobian from z1 integration

      m10=m10*pre10*(1-d-x10)

      endif

c*********************************************************
c     Result

      a1=dlog(1-z2)/(1-z2)*(pre1*a1-pre10*a10)*jac

      b1=(pre1*b1-pre10*b10)/(1-z2)*jac

      c1=pre1*c1*jac
      

c     For d10 only jacobian from z1 integration

      d10=pre10*d10*(1-d-x10)


      r1=r1*pre1*jac

      dcross1=(a1+b1+c1+d10-m10+r1)/shad

      

      distrcross=2*qt*jacy*dcross1

      return
      end

c*********************************************************

c*********************************************************
c     Pgg (with no delta term and to be divided by (1-z)+

      function Pgg(z)
      implicit none
      real *8 z,Pgg
      Pgg=3*(1+z**4+(1-z)**4)/z
      return 
      end

c*********************************************************
c     Pqq (with no delta term and to be divided by (1-z)+

      function Pqq(z)
      implicit none
      real *8 z,Pqq
      Pqq=4d0/3*(1+z**2)
      return 
      end


c*********************************************************
c     Pqg 

      function Pqg(z)
      implicit none
      real *8 z,Pqg
      Pqg=0.5d0*(z**2+(1-z)**2)
      return
      end

c*********************************************************
c     Pgq 

      function Pgq(z)
      implicit none
      real *8 z,Pgq
      Pgq=4d0/3*(1+(1-z)**2)/z
      return
      end



c*********************************************************

      function gg0(x1,x2)
      implicit none
      real *8 gg0,x1,x2
      real *8 yh
      real *8 uh,th,sh,tm
      common/tm/tm
      common/rap/yh

      include 'scales.F'
      
      sh=x1*x2*shad
      th=q2-sroot*x2*tm*dexp(yh)
      uh=q2-sroot*x1*tm*dexp(-yh)

      gg0=3*(q2**4+sh**4+th**4+uh**4)/(uh*th*sh)
      return
      end

c**********************************************************

      function qg0(x1,x2)
      implicit none
      real *8 qg0,x1,x2
      real *8 yh
      real *8 uh,th,sh,tm
      common/tm/tm
      common/rap/yh
      
      include 'scales.F'      

      sh=x1*x2*shad
      th=q2-sroot*x2*tm*dexp(yh)
      uh=q2-sroot*x1*tm*dexp(-yh)

      qg0=4d0/3*(uh**2+sh**2)/(-th)
      return
      end

      function gq0(x1,x2)
      implicit none
      real *8 gq0,x1,x2
      real *8 y1,y2,yh
      real *8 uh,th,sh,tm
      common/tm/tm
      common/internal/y1,y2
      common/rap/yh

      include 'scales.F'
      
      
      sh=x1*x2*shad
      th=q2-sroot*x2*tm*dexp(yh)
      uh=q2-sroot*x1*tm*dexp(-yh)

      gq0=4d0/3*(th**2+sh**2)/(-uh)
      return
      end

      function qqb0(x1,x2)
      implicit none
      real *8 qqb0,x1,x2
      real *8 yh
      real *8 uh,th,sh,tm
      common/tm/tm
      common/rap/yh
      
      include 'scales.F'

      sh=x1*x2*shad
      th=q2-sroot*x2*tm*dexp(yh)
      uh=q2-sroot*x1*tm*dexp(-yh)

      qqb0=2*16d0/9*(th**2+uh**2)/sh
      return
      end

C**********************************************************


      subroutine coeff
c     Returns the values of some lenghty coefficients

      implicit none
      real *8 beta0,tmp,tiny
      real *8 uh,th,sh,QQ2,xa,xb
      real *8 big1,big2,big3,big4,big5

      common/coef/uh,th,sh
      common/beta0/beta0
      common/big/big1,big2,big3,big4,big5

      include 'scales.F'

      QQ2=uh+sh+th-q2
      xa=th/(th-QQ2)
      xb=uh/(uh-QQ2)


c     gg coefficients


      big1=0.5d0*9*((q2**4+sh**4+QQ2**4+uh**4+th**4)
     #+xa*xb*(q2**4+sh**4+QQ2**4+(uh/xb)**4+(th/xa)**4))/sh/uh/th


      big2=beta0*3d0/2*
     #(q2**4+sh**4+xa*xb*((uh/xb)**4+(th/xa)**4))/sh/uh/th

c     be careful when QQ2=0 to avoid numerical problems

      tiny=1d-4
      tmp=tiny*qt2


      if(QQ2.gt.tmp) then
      big3=9*((q2**4+sh**4+QQ2**4+(uh/xb)**4+(th/xa)**4)*(2*QQ2+qt**2)
     #/sh**2/QQ2/(QQ2+qt**2)+2*q2**2*((q2-th)**4+(q2-uh)**4
     #+uh**4+th**4)/sh/uh/th/(q2-uh)/(q2-th))
     #/qt**2*dlog(qt**2/(qt**2+QQ2))
      else
      big3=-9*((q2**4+sh**4+uh**4+th**4)/sh**2/qt**4)
      endif      
   

c     qg coefficients: be careful when QQ2=0 to avoid numerical problems     


      big4=4*((-sh**3*th-sh*th**3+QQ2**3*th+QQ2*th**3)/sh/uh/th
     #+xa*xb*(-sh**3*th/xa-sh*(th/xa)**3-QQ2**3*uh/xb-QQ2*(uh/xb)**3)
     #/(sh*uh*th))


      if(QQ2.gt.tmp) then
      big5=4*((-sh**3*th/xa-sh*(th/xa)**3-QQ2**3*uh/xb-QQ2*(uh/xb)**3)
     #*(2*QQ2+qt**2)/sh**2/(QQ2+qt**2)
     #-2*QQ2*q2**2*((q2-th)**2+th**2)/sh/uh/(q2-uh))
     #*dlog(qt**2/(qt**2+QQ2))/(QQ2*qt**2)
      else
      big5=-4*(-sh**3*th/xa-sh*(th/xa)**3)/sh**2/qt**4
      endif

      
      return 
      end



C**********************************************************
C     Returns regular part of coefficient functions

      subroutine REG

      implicit none      
      real *8 A,B,L1a,L1b,L2a,L2b
      real *8 Rgq,Rqqb,A0,Aep,REGgg,REGgq,REGqg,REGqqb,REGqqpb,REGqq
      real *8 A1234,A3412,A1324,A3241,B1pmB1pp,B2pmB2pp
      real *8 uh,th,sh,QQ2,xa,xb
      real *8 XQT2,E2,E4
      integer nf

      common/nf/nf
      common/regu/uh,th,sh,REGgg,REGgq,REGqg,REGqqb,REGqqpb,REGqq


      include 'scales.F'
      

      QQ2=sh+th+uh-q2

      xQT2=QQ2+qt2

      xa=th/(th-QQ2)
      xb=uh/(uh-QQ2)


c     Eq. (A.4)

      A0=((th/xa)**4+(uh/xb)**4)*qt2*QQ2/xQT2**2
     #*(5-7*QQ2/xQT2+20d0/3*(QQ2/xQT2)**2)
     #+sh**2*QQ2*qt2*(17d0/3+4*dlog(qt2/xQT2))


c     Eq. (A.12)

      Aep=4*(sh*qt2*q2*QQ2)**2*(1d0/th**4+1d0/uh**4)


c     Regular part of gg coefficient

      REGgg=1/sh**2/qt2/QQ2*(9*(A0+Aep
     #+A1234(uh,th,sh)+A1234(th,uh,sh)
     #+A3412(uh,th,sh)+A3412(th,uh,sh)
     #+A1324(uh,th,sh)+A1324(th,uh,sh)
     #+A3241(uh,th,sh)+A3241(th,uh,sh))
     #+4d0/3*nf
     #*(B1pmB1pp(uh,th,sh)+B1pmB1pp(th,uh,sh))
     #+3*nf
     #*(B2pmB2pp(uh,th,sh)+B2pmB2pp(th,uh,sh)))

c     Regular part of gq coefficient

      REGgq=Rgq(uh,th,sh)

c     Regular part of qg coefficient

      REGqg=Rgq(th,uh,sh)

c     Regular part of qqb coefficient

      REGqqb=Rqqb(uh,th,sh)+Rqqb(th,uh,sh)

c     Regular part of qqpb coefficient (also qqp)


      E2=sh*qt2*QQ2*(QQ2-qt2)/(qt2+QQ2)**2*((QQ2-th)**2+(QQ2-uh)**2)
     #+2*sh**2*qt2*QQ2

      REGqqpb=16d0/9/sh**2/qt2/QQ2*E2
    

c     Regular part of qq coefficient

   
      E4=2*sh*qt2*QQ2*(sh**2+QQ2**2)/(qt2+QQ2)*dlog(qt2/QQ2)
     #+4*sh**2*qt2*QQ2*dlog(qt2/(qt2+QQ2))



      REGqq=16d0/9/sh**2/qt2/QQ2*(E2+E4/3)



      return
      end

c********************************************************
c     Eq. (A.5)

      function A1234(u,t,s)
      
      implicit none
      real *8 x1,x2,u,t,s,A1234
      real *8 u2,t2,s2,u3,u4,t4,QQ2,q4
      real *8 xQt2
      real *8 A,B,L1a,L1b,L2a,L2b,L3


      include 'scales.F'
      
      q4=q2**2

      QQ2=u+t+s-q2

      x1=t/(t-QQ2)
      x2=u/(u-QQ2)

  
      xQT2=QQ2+qt2


      t2=t**2
      t4=t**4

      u2=u**2
      u3=u**3
      u4=u2**2

      s2=s**2
     

c     Definitions in Eqs. (A.2)-(A.3)

      A=s+q2-QQ2
      B=dsqrt(A**2-4*q2*s)

      L1a=dlog(q2/s/x1**2)
      L1b=dlog(q2/s/x2**2)
      L2a=dlog(q2*s/(A-s*x1)**2)
      L2b=dlog(q2*s/(A-s*x2)**2)

      L3=dlog((A+B)/(A-B))

   

      A1234=-0.5d0*(((s*qt2/t)**4+(q2*QQ2/t)**4)*L1a+u**4*L1b)
     #+0.5d0*s*q4*QQ2*u3/t/(u-q2)/(t-q2)*(L1a+L1b)
     #+0.5d0*s*qt2*QQ2*u3/A/(t-q2)*(L2b-L1b)
     #+0.5d0*s*qt2*QQ2*(q4**2+(u-q2)**4)/A/t/(u-q2)*(L2a-L1a)
     #+s*qt2*q2*QQ2*(u4/(B*t)**2+u2/B**2)/2
     #+(s*qt2*q2*QQ2)**2*(-6/B**4-4/t4+8/(B*t)**2)
     #+L3*(s*qt2*u3*(u+t)/B/t
     #+(s*qt2*q2*QQ2)**2*(3*A/B**5-1d0/A/B**3)
     #-s*qt2*q2*QQ2*((t2+t*u+4*u**2-2*q2*QQ2)/B/t
     #+A*(t2+3*t*u+3*u**2-6*q2*QQ2)/2/B**3
     #+(t2+t*u+7*u2-2*q2*QQ2)/(2*A*B)))


      return
      end

c****************************************************
c     Eq. (A.7)

      function A3412(u,t,s)

      implicit none
      real *8 x1,x2,u,t,s,A3412
      real *8 s2,s3,s4,u2,u3,u4,t2,t3,t4,QQ2,q4,q6,q8
      real *8 xQt2
      real *8 A,B,L1a,L1b,L2a,L2b,L3
      

      include 'scales.F'
      
      q4=q2**2   
      q6=q2**3
      q8=q2**4
      
      QQ2=u+t+s-q2

      x1=t/(t-QQ2)
      x2=u/(u-QQ2)

      xQT2=QQ2+qt2
      

      t2=t**2
      t3=t**3
      t4=t**4

      u2=u**2
      u3=u**3
      u4=u**4

      s2=s**2
      s3=s**3
      s4=s2**2

    

c     Definitions in Eqs. (A.2)-(A.3)

      A=s+q2-QQ2
      B=dsqrt(A**2-4*q2*s)

      L1a=dlog(q2/s/x1**2)
      L1b=dlog(q2/s/x2**2)
      L2a=dlog(q2*s/(A-s*x1)**2)
      L2b=dlog(q2*s/(A-s*x2)**2)

      L3=dlog((A+B)/(A-B))


      A3412=s*qt2*QQ2*A**3/2/t/(u-q2)*(L2a+L1b)
     #+s*qt2*(u+t)/16/u/t/B*(A**4+6*A**2*B**2+B**4)*L3
     #+(-s*qt2/2/u/t*((s-QQ2)**4+q8+2*QQ2*A*(s-QQ2)**2-2*QQ2*q6)
     #-(s*qt2*q2*QQ2)**2/u4+2*s*qt2*q2*QQ2*(A**2-s*q2)/u2)*L1b
     #+s*qt2/(8*u*t)*((QQ2-u)**3/(QQ2-t)*((QQ2-t)*s+QQ2*u)
     #+(QQ2-t)**3/(QQ2-u)*((QQ2-u)*s+QQ2*t))
     #*(4d0/3+2*qt2/xQt2+4*(qt2/xQt2)**2-44d0/3*(qt2/xQt2)**3)
     #+s*qt2*(QQ2-u)**2/4/u/t/(QQ2-t)*
     #(-3*(t-q2)*((QQ2-t)*s+QQ2*u)-QQ2*(q2*(t-q2)+QQ2*(u-q2)))
     #*(1+2*qt2/xQt2-6*(qt2/xQt2)**2)
     #+s*qt2*(QQ2-t)**2/4/u/t/(QQ2-u)
     #*(-3*(u-q2)*((QQ2-u)*s+QQ2*t)+3*QQ2*(q2*(u-q2)
     #+QQ2*(t-q2))+4*u*s2)*(1+2*qt2/xQt2-6*(qt2/xQt2)**2)
     #+s*qt2*(QQ2-u)/2/u/t/(QQ2-t)*(3*(t-q2)**2*((QQ2-t)*s+QQ2*u)
     #+3*(t-q2)*QQ2*(q2*(t-q2)+QQ2*(u-q2))
     #+QQ2*u*(q2*(t-QQ2)+QQ2*(u-q2)))*(1-2*qt2/xQt2)
     #+s*qt2*(QQ2-t)/2/u/t/(QQ2-u)*(3*(u-q2)**2*((QQ2-u)*s+QQ2*t)
     #+8*u*t*s2+2*u*s3-2*QQ2*u*(u-QQ2)**2
     #-3*q2*QQ2*(t-q2)**2-3*QQ2*(q2-QQ2)*t2
     #-QQ2*u*(4*u*t-u*q2-QQ2*t+2*t2-4*q4)
     #+3*q2*QQ2**2*(t-q2)+q2*QQ2*u*(t-QQ2))*(1-2*qt2/xQt2)
     #-4*(s*qt2*q2*QQ2)**2/u4+s*qt2*q2*QQ2*(B**2)/(2*u2)
     #+s2*qt2*q4/6*((s+QQ2)/u/t+QQ2/u2+QQ2/t2)
     #+2*s2*qt2*QQ2*q4/u2+s2*qt2*q4/u-s2*qt2/12/u/t
     #*(30*q6+54*(QQ2**2)*q2+8*QQ2**3)
     #+s*qt2/12/u/t*(11*s4+17*q8+QQ2*(61*u2*t+17*u3+73*u*t2
     #+29*t3)+q2*(24*u2*t+6*u3+36*u*t2+18*t3)
     #+QQ2**2*(-21*u2-33*t2-52*u*t)
     #+q2*QQ2*(-73*u2-109*t2-170*u*t)
     #+q4*(-23*u2-35*t2-52*u*t)
     #+q4*QQ2*(134*t+110*u)+4*QQ2**4+52*q2*QQ2**3+20*q4*QQ2**2
     #-22*q6*QQ2)



      return
      end

c********************************************************
c     Eq. (A.9)

      function A1324(u,t,s)
      
      implicit none
      real *8 x1,x2,u,t,s,A1324
      real *8 u2,t2,s2,u3,u4,t4,QQ2,q4,q8
      real *8 xQt2
      real *8 A,B,L1a,L1b,L2a,L2b,L3

      include 'scales.F'
      
      q4=q2**2
      q8=q2**4

      QQ2=u+t+s-q2
      x1=t/(t-QQ2)
      x2=u/(u-QQ2)


      xQT2=QQ2+qt2


      t2=t**2
      t4=t**4

      u2=u**2
      u3=u**3
      u4=u2**2

      s2=s**2
     


c     Definitions in Eqs. (A.2)-(A.3)

      A=s+q2-QQ2
      B=dsqrt(A**2-4*q2*s)

      L1a=dlog(q2/s/x1**2)
      L1b=dlog(q2/s/x2**2)
      L2a=dlog(q2*s/(A-s*x1)**2)
      L2b=dlog(q2*s/(A-s*x2)**2)

      L3=dlog((A+B)/(A-B))


      A1324=-0.5d0*(((s*qt2/t)**4+(q2*QQ2/t)**4)*L1a+u4*L1b)
     #+s2*qt2*q4*QQ2/t2*L1a
     #+(s*q4*QQ2*u3/t/(u-q2)/(t-q2)+s*qt2*u3/t)/2*(L1a+L1b)
     #+s2*qt2*(1-x2)*u3/A/2/(t-q2)*(L2b-L1b)
     #+s2*qt2*(1-x1)*(q8+(u-q2)**4)/2/A/t/(u-q2)*(L2a-L1a)
     #+s2*qt2*q4*QQ2/A/B*L3+s*qt2*q2*QQ2/2/t4
     #*((s*qt2)**2-6*s*qt2*q2*QQ2+q4*QQ2**2)

      return
      end

c********************************************************
c     Eq. (A.9)
c     TO BE SYMMETRIZED !

      function A3241(u,t,s)
      
      implicit none
      real *8 x1,x2,u,t,s,A3241
      real *8 u2,t2,s2,u3,t4,QQ2,q4
      real *8 xQt2
      real *8 A,B,L1a,L1b,L2a,L2b,L3


      include 'scales.F'
      
      q4=q2**2

      QQ2=u+t+s-q2
      x1=t/(t-QQ2)
      x2=u/(u-QQ2)

      xQT2=QQ2+qt2


      t2=t**2
      t4=t**4

      u2=u**2
      u3=u**3

      s2=s**2
     

c     Definitions in Eqs. (A.2)-(A.3)

      A=s+q2-QQ2
      B=dsqrt(A**2-4*q2*s)

      L1a=dlog(q2/s/x1**2)
      L1b=dlog(q2/s/x2**2)
      L2a=dlog(q2*s/(A-s*x1)**2)
      L2b=dlog(q2*s/(A-s*x2)**2)

      L3=dlog((A+B)/(A-B))


      A3241=s2*qt2*A**3*(1-x1)/2/t/(u-q2)*(L2a-L1a)
     #+(-(s*qt2*q2*QQ2)**2/t4+s*qt2*q4*QQ2**2/u/t
     #-s*qt2*q2*QQ2*A**4/2/u/t/(u-q2)/(t-q2)
     #+s*qt2*QQ2*q2*(u+t)*(2*A**2-s*q2)/u/t2)*L1a
     #+s2*qt2*QQ2*(QQ2-u)**2/2/u/t/(QQ2-t)**2
     #*(-u*t-(QQ2-t)**2)*(-3+10*QQ2/xQt2-6*(QQ2/xQt2)**2)
     #+s2*qt2*QQ2*(QQ2-u)/u/t/(QQ2-t)**2*(u*t*(QQ2-u)-(QQ2-t)**3
     #-q2*(QQ2-t)**2-q2*(QQ2-t)*(QQ2-u))*(-1+2*QQ2/xQt2)
     #+s*qt2*q2*QQ2*(B**2/2/t2-2*q2*QQ2/t2+(u+t)**2/2/u/t)
     #-4*(s*qt2*q2*QQ2)**2/t4+s2*qt2*QQ2/4/u/t
     #*((t+u)**2-(t+u)*(6*QQ2+4*q2)+6*QQ2**2+8*q2*QQ2)
     #+s2*qt2*QQ2*q4*(t+u)**2/4/u2/t2


      return
      end



c********************************************************
c     Eq. (A.13)+(A.15)
c     TO BE SYMMETRIZED !

      function B1pmB1pp(u,t,s)
      
      implicit none
      real *8 x1,x2,u,t,s,B1pm,B1pp,B1pmB1pp
      real *8 u2,t2,s2,u3,t4,QQ2,q4
      real *8 xQt2

      include 'scales.F'
      
      q4=q2**2

      QQ2=u+t+s-q2
      x1=t/(t-QQ2)
      x2=u/(u-QQ2)

      xQT2=QQ2+qt2


      t2=t**2
      t4=t**4

      u2=u**2
      u3=u**3

      s2=s**2
     


      B1pm=s2**2*qt2*x1*(1-x1)**3/t+s2**2*qt2**3*x1**3*(1-x1)/t**3
     #+4*s2**2*qt2**2*x1**2*(1-x1)**2/t2-s2*qt2*QQ2*(1+dlog(qt2/xQt2))

      B1pp=s2**2*qt2*x1**3*(1-x1)/t+s2**2*qt2**3*x1*(1-x1)**3/t**3
     #+4*s2**2*qt2**2*x1**2*(1-x1)**2/t2
     #-s2*qt2*QQ2/xQt2**2/u/t*((u*t+qt2*QQ2)**2+2*s*qt2*QQ2*xQt2)
     #+s2*qt2*QQ2/u/t*(s2+QQ2**2)*dlog(xQt2/QQ2)

      B1pmB1pp=B1pm+B1pp

      return
      end



c*******************************************************
c     Eq. (A.14)+(A.16)
c     TO BE SYMMETRIZED !

      function B2pmB2pp(u,t,s)
      
      implicit none
      real *8 x1,x2,u,t,s,B2pm,B2pp,B2pmB2pp
      real *8 u2,t2,s2,s3,s4,u3,t4,QQ2,q4,q6
      real *8 xQt2

      include 'scales.F'
      
      q4=q2**2
      q6=q2**3

      QQ2=u+t+s-q2
      x1=t/(t-QQ2)
      x2=u/(u-QQ2)

      xQT2=QQ2+qt2


      t2=t**2
      t4=t**4

      u2=u**2
      u3=u**3

      s2=s**2
      s3=s**3
      s4=s2**2
     

      B2pm=1/3d0*(t/x1)**4*qt2/xQt2*(qt2**3-QQ2**3-xQt2**3)/xQt2**3
     #-s2*qt2*QQ2/3

      B2pp=-s2*qt2*QQ2/(2*u*t)*(s2+QQ2**2)*dlog(xQt2/QQ2)
     #+s*qt2*((QQ2-u)**3)/(2*u*t)/(QQ2-t)*((QQ2-t)*s+QQ2*u)
     #*(2d0/3+QQ2/xQt2-10d0/3*(QQ2/xQt2)**3)
     #-s*qt2*(QQ2-u)**2/2/u/t/(QQ2-t)**2*(3*(QQ2-t)**3*QQ2
     #+(QQ2-t)*QQ2*(2*u*t+q4)+(QQ2-t)**2
     #*(s2+4*q2*QQ2-u*(QQ2+q2))-u2*QQ2**2+u*t2*q2)
     #*(1-2*(QQ2/xQt2)**2)
     #+s*qt2*(QQ2-u)/2/u/t/(QQ2-t)
     #*(3*QQ2*s*(QQ2+s)*(QQ2-t)-t*s3+q2*QQ2*s2+QQ2*u*(q2-QQ2)**2)
     #*(1-2*QQ2/xQt2)
     #+s*qt2/12/u/t*(-2*s4+6*s*q2*t*(t-q2)
     #+2*s*q6+8*QQ2*s*(s-QQ2)**2-2*u*t*s*QQ2+7*s2*q2*QQ2
     #-2*s*QQ2**2*q2-q6*QQ2+3*q2*QQ2**3-4*u*t*q2*QQ2)
     #+11d0/6*s3*qt2*QQ2**2/u/t-s2*qt2*q4*QQ2/3/t2


      B2pmB2pp=B2pm+B2pp

      return
      end

C************************************************************************
c     Regular part of qg coefficient function

      function Rgq(u,t,s)
      implicit none
      real *8 Rgq,x1,x2,u,t,s
      real *8 C1pm,C1mp,C1pp,C1mm,C2pm,C2mp,C2pp,C2mm,C2ep
      real *8 u2,u4,t2,t4,s2,s4,QQ2,q4,q6
      real *8 xQt2
      real *8 A,B,L1a,L1b,L2a,L2b,L3

      include 'scales.F'

      
      q4=q2**2
      q6=q2**3

      QQ2=u+t+s-q2
      x1=t/(t-QQ2)
      x2=u/(u-QQ2)

      xQT2=QQ2+qt2

      s2=s**2
      s4=s2**2
      u2=u**2
      u4=u2**2
      t2=t**2
      t4=t2**2

c     Definitions in Eqs. (A.2)-(A.3)

      A=s+q2-QQ2
      B=dsqrt(A**2-4*q2*s)

      L1a=dlog(q2/s/x1**2)
      L1b=dlog(q2/s/x2**2)
      L2a=dlog(q2*s/(A-s*x1)**2)
      L2b=dlog(q2*s/(A-s*x2)**2)

      L3=dlog((A+B)/(A-B))


c     Eqs (A.18)-(A.26)
      

      C1pm=-2*s2*qt2*QQ2*dlog(qt2/XQT2)

      C2pm=0d0
      
      C1mp=s2*qt2*QQ2-1.5d0*s2*qt2*t2/u+qt2/2/xQT2*(s*(QQ2-t)**3
     #+QQ2*(QQ2-u)**3)*(-3+10*QQ2/xQT2-6*(QQ2/xQT2)**2)
      
      C2mp=qt2*QQ2/xQT2**2*(s*(QQ2-t)**3+QQ2*(QQ2-u)**3)
     #*(-2+3*QQ2/xQT2)+2*s2*qt2*QQ2+4*s2*qt2*QQ2*dlog(qt2/xQT2)


      C1pp=-1.5d0*s4*qt2/u-s*qt2*QQ2*A**2/u*L2b+s*qt2/u/t2*L1a
     #*((A-q2)**2*s*t2-2*QQ2*q2*u*t*A-QQ2*q4*(QQ2-t)*u)
     #+s*qt2/u/B*L3*((s+QQ2-q2)*(s*(A-q2)**2+QQ2*A**2)
     #-4*s*QQ2*A*(A-q2))+0.5d0*s*qt2*((QQ2-t)**2*(s/(QQ2-u)-QQ2/u)
     #+(QQ2-u)**2*(QQ2/(QQ2-t)+s/u))*(-3+10*QQ2/xQT2-6*(QQ2/xQT2)**2)
     #+s*qt2*(QQ2-t)/u/(QQ2-u)*(2*s*u*(s+t)
     #-QQ2*(4*q2*QQ2-QQ2*t-q2*u-2*u*t))*(-1+2*QQ2/xQT2)
     #+s*qt2*(QQ2-u)/u/(QQ2-t)*(t*s2-2*u*t*s+2*QQ2*u*(QQ2-t))
     #*(-1+2*QQ2/xQT2)+s*qt2*QQ2*q2*(u+t)/t-2*s*qt2*QQ2**2*q4/t2
     #+s2*qt2*QQ2*q4/2/u2+s*qt2/2/u*(-2*(QQ2+q2)*u*s+2*q2*s2
     #+q4*s+q2*QQ2*(2*(s-QQ2)+3*q2-u)+5*QQ2*s*(s-QQ2))

 

      C2pp=0.5d0*s*qt2*A**2*(1-x1)*(L2a-L1a)+s*qt2*(q2-t)*A**2
     #*(1-x2)/2/u*(L1b-L2b)+0.5d0*s*qt2/u*(L1b-L1a)*(s-QQ2)**3
     #+s*qt2*QQ2*A**2/(q2-u)*(L1b-L2a)+s*qt2*QQ2/u2*L1b
     #*(2*u*(s-QQ2)**2+4*q2*(q2-t)*A-2*q4*(QQ2-t)-q6)
     #-0.5d0*s*qt2*((QQ2-t)**2*(s/(QQ2-u)-QQ2/u)
     #+(QQ2-u)**2*(QQ2/(QQ2-t)+s/u))*(-3+10*QQ2/xQT2-6*(QQ2/xQT2)**2)
     #+0.5d0*s*qt2*(QQ2-t)/u/(QQ2-u)*((-3*s*qt2+u2-QQ2**2)*(s+QQ2)
     #-q2*s*u+2*QQ2**2*(QQ2-u))*(-1+2*QQ2/xQT2)
     #+0.5d0*s*qt2*(QQ2-u)/u/(QQ2-t)*(3*s*qt2*(s+QQ2)-u*QQ2*(QQ2-t)
     #+3*QQ2*q2*(QQ2-u)+s*t*(u+s))*(-1+2*QQ2/xQT2)
     #+0.5d0*s*qt2*q2*QQ2/u2*(2*(s-QQ2)**2-2*q2*(s-q2)
     #-u*(QQ2-u)-4*q2*QQ2)+s2*qt2*(u-t)*(q2+QQ2)/(2*u)
     #-2*(s*qt2*q2*QQ2)**2/u4*(4+L1b)




      C1mm=s2*qt2*t2/u*L1a-s*qt2*QQ2*(q2-t)**2/u*L2b
     #+s*qt2*q2*QQ2/B**2*(t*(u+t)-2*q2*QQ2)+s*qt2/u/B
     #*(t2*B**2-q2*t2*(u+t)+2*QQ2**2*q4+QQ2*q4*(3*t-u)
     #+QQ2*q4*u/B**2*(-t*(u+t)+2*q2*QQ2+QQ2*(t-u)))*L3



      C2mm=s*qt2*t2*QQ2/(2*u)*(L1a+3*L1b)
     #+0.5d0*s*qt2*t2*(1-x1)*(L2a-L1a)+s*qt2*QQ2*(q2-t)**3*x2/(2*u2)
     #*(L2b-L1b)+s*qt2*t2*QQ2/(q2-u)*(L1b+L2a)
     #+s2*qt2*t2/2/u*(L1b-L1a)+s*qt2*q2*QQ2/u2
     #*(4*t*(t-q2)+q4)*L1b-2*(s*qt2*q2*QQ2)**2/u4*(L1b+4)
     #+s*qt2*t2*q2*QQ2/u2

 
      C2ep=4/u4*(s*qt2*q2*QQ2)**2


      Rgq=1/s2/qt2/QQ2*(16d0/9*(C1pm+C1mp+C1pp+C1mm)
     #+4*(C2pm+C2mp+C2pp+C2mm+C2ep))


      return
      end

C*************************************************************************
c     Regular part of qqb coefficient function
c     Remember: ADD (u<->t) PERMUTATION !

      function Rqqb(u,t,s)
      implicit none
      real *8 Rqqb,x1,x2,u,t,s
      real *8 D1pm,D1pp,D2pm,D2pp,E1,E2,E3
      real *8 u2,t2,s2,s3,QQ2,q4
      real *8 xQt2
      real *8 A,B,L1a,L1b,L2a,L2b,L3
      integer nf
      common/nf/nf

      include 'scales.F'
      
      q4=q2**2


      QQ2=u+t+s-q2
      x1=t/(t-QQ2)
      x2=u/(u-QQ2)

      xQT2=QQ2+qt2

      s2=s**2
      s3=s**3
      u2=u**2
      t2=t**2


c     Definitions in Eqs. (A.2)-(A.3)

      A=s+q2-QQ2
      B=dsqrt(A**2-4*q2*s)

      L1a=dlog(q2/s/x1**2)
      L1b=dlog(q2/s/x2**2)
      L2a=dlog(q2*s/(A-s*x1)**2)
      L2b=dlog(q2*s/(A-s*x2)**2)

      L3=dlog((A+B)/(A-B))


c     Definitions in Eqs (A.28)-(A.34)


      D1pm=-s2*qt2*QQ2*(1+dlog(qt2/xQT2))-s3*qt2**2*x1*(1-x1)/t
     #-s3*qt2*(1-x1)**2


      D2pm=-s2*qt2*QQ2/3-s*qt2*t2/6/x1**2
     #*(11-12*QQ2/xQT2+3*(QQ2/xQT2)**2)+11d0*s*qt2*t2/6


      D1pp=s2*qt2*u2*(1-x2)*(L2b-L1b)/A+s2*qt2*(q2-u)**2*(1-x1)/A
     #*(L2a-L1a)+s*qt2*q2*QQ2*(s*qt2+u*t)/t2*L1a
     #-2*s2*qt2*q4*QQ2/A/B*L3+s*qt2*q2*QQ2*(2*s*qt2-u*t)/t2


      D2pp=s*qt2*u2*(q2-t)*(1-x2)/(2*A)*(L1b-L2b)+s*qt2*(q2-u)**3
     #*(1-x1)/(2*A)*(L1a-L2a)-0.5d0*s*qt2*u2*(L1a+L1b)
     #+6*(s*qt2*q2*QQ2)**2/B**4-s*qt2*q2*QQ2*u2/B**2
     #+L3*(s*qt2*u2*(u+t)/B+(s*qt2*q2*QQ2)**2*(1/A/B**3-3*A/B**5)
     #+s*qt2*q2*QQ2*((t-3*u)/2/B+A*(B**2+2*u2)/4/B**3
     #+(t2-6*u*t+7*u2)/4/A/B))


      E1=4d0/3*(2*s2*qt2*QQ2-s*qt2*q2*QQ2)

      E2=s*qt2*QQ2*(QQ2-qt2)/xQT2**2*
     #((QQ2-t)**2+(QQ2-u)**2)+2*s2*qt2*QQ2
      
      E3=-2*s*qt2*((u+t-2*QQ2)**2-2*s*qt2)*dlog(qt2/xQT2)
     #-s*qt2*QQ2*(2*xQT2+QQ2)/xQT2**2*
     #((QQ2-t)**2+(QQ2-u)**2)-6*s2*qt2*QQ2
      

c     Eq. (A.27) (add (u<->t) permutation)


      Rqqb=(128d0/27*(D1pm+D1pp)+32d0/3*(D2pm+D2pp)
     #+0.5d0*(16d0/9*nf*E1+16d0/9*E2+16d0/27*E3))/s2/qt2/QQ2


      return
      end


C*************************************************************************

C.....Luminosities

! Channels
!      if(channel.eq.1) then
!         write(*,*) 'higgsqtchannel='
!      endif

C.....COSTRUCT GG Luminosity (jinst=1)

      FUNCTION ALLGG(X1,X2)
      IMPLICIT REAL*8(A-H, O-Z)
      integer channel
      common/channel/channel
      common/isetproton/isetproton
      common/fixvar/xmu2
      common/nf/nf

      dimension sf(1:6,1:2)


      CALL STRFUN(X1,X2,SF)
      if(channel.eq.1 .or. channel.eq.2)then
      ALLGG=SF(1,1)+SF(1,2)                      !TOT,GG
!         write(*,*) 'ALLGG=',ALLGG
      else
      ALLGG=0.0d0                                !GQ,QQ
!         write(*,*) 'ALLGG=0'
      endif
      RETURN
      END

C.....COSTRUCT qg Luminosity (jinst=2)

      FUNCTION ALLqg(X1,X2)
      IMPLICIT REAL*8(A-H, O-Z)
      integer channel
      common/channel/channel
      common/isetproton/isetproton
      common/fixvar/xmu2
      common/nf/nf
      dimension sf(1:6,1:2)

      CALL STRFUN(X1,X2,SF)
      if(channel.eq.1 .or. channel.eq.3)then
      ALLqg=SF(2,1)                              !TOT,GQ
!         write(*,*) 'ALLqG='
      else 
      ALLqg=0.0d0                                !GG,QQ
!         write(*,*) 'ALLqG=0'
      endif 
      RETURN
      END

C.....COSTRUCT gq Luminosity (jinst=2)

      FUNCTION ALLgq(X1,X2)
      IMPLICIT REAL*8(A-H, O-Z)
      integer channel
      common/channel/channel
      common/isetproton/isetproton
      common/fixvar/xmu2
      common/nf/nf
      dimension sf(1:6,1:2)

      CALL STRFUN(X1,X2,SF)
      if(channel.eq.1 .or. channel.eq.3)then
      ALLgq=SF(2,2)                              !TOT,GQ
!         write(*,*) 'ALLgq='
      else
      ALLgq=0.0d0                                !GG,QQ
!         write(*,*) 'ALLgq=0'
      endif
      RETURN
      END



C.....COSTRUCT qbarq Luminosity (jinst=3)

      FUNCTION ALLqqb(X1,X2)
      IMPLICIT REAL*8(A-H, O-Z)
      integer channel
      common/channel/channel
      common/isetproton/isetproton
      common/fixvar/xmu2
      common/nf/nf
      dimension sf(1:6,1:2)

      CALL STRFUN(X1,X2,SF)
      if(channel.eq.1 .or. channel.eq.4)then
      ALLqqb=SF(3,1)+SF(3,2)                     !TOT,QQ
!         write(*,*) 'ALLqqb='
      else
      ALLqqb=0.0d0                               !GG,GQ
!         write(*,*) 'ALLqqb=0'
      endif
      RETURN
      END

C.....qq luminosity (jinst=4)

      FUNCTION ALLqq(X1,X2)
      IMPLICIT REAL*8(A-H, O-Z)
      integer channel
      common/channel/channel
      common/isetproton/isetproton
      common/fixvar/xmu2
      common/nf/nf
      dimension sf(1:6,1:2)

      CALL STRFUN(X1,X2,SF)
      if(channel.eq.1 .or. channel.eq.4)then
      ALLqq=SF(4,1)+SF(4,2)                      !TOT,QQ
!         write(*,*) 'ALLqq='
      else
      ALLqq=0.0d0                                !GG,GQ
!         write(*,*) 'ALLqq=0'
      endif
      RETURN
      END

C.....qQb luminosity (jinst=5)

      FUNCTION ALLqqpb(X1,X2)
      IMPLICIT REAL*8(A-H, O-Z)
      integer channel
      common/channel/channel
      common/isetproton/isetproton
      common/fixvar/xmu2
      common/nf/nf
      dimension sf(1:6,1:2)

      CALL STRFUN(X1,X2,SF)
      if(channel.eq.1 .or. channel.eq.4)then
      ALLqqpb=SF(5,1)+SF(5,2)                    !TOT,QQ
!         write(*,*) 'ALLqqpb='
      else
      ALLqqpb=0.0d0                              !GG,GQ
!         write(*,*) 'ALLqqpb=0'
      endif
      RETURN
      END

C.....qQ luminosity (jinst=6)

      FUNCTION ALLqqp(X1,X2)
      IMPLICIT REAL*8(A-H, O-Z)
      integer channel
      common/channel/channel
      common/isetproton/isetproton
      common/fixvar/xmu2
      common/nf/nf
      dimension sf(1:6,1:2)

      CALL STRFUN(X1,X2,SF)
      if(channel.eq.1 .or. channel.eq.4)then
      ALLqqp=SF(6,1)+SF(6,2)                     !TOT,QQ
!         write(*,*) 'ALLqqp='
      else
      ALLqqp=0.0d0                               !GG,GQ
!         write(*,*) 'ALLqqp=0'
      endif
      RETURN
      END



c**********************************************************

C.....Inclusion of two heavy quarks: amplitude level

      FUNCTION ABI(X)
      IMPLICIT NONE
      REAL*8 X,FR,FI,root,dasin,PI,etap,etam,rlog
      REAL*8 AR,AI
      COMPLEX*16 ABI
      PI=4.D0*DATAN(1.D0)
      IF(X.GE.0.25D0) THEN
                FR=-2.D0*(DASIN(0.5D0/DSQRT(X)))**2
                FI=0.D0
      ELSEIF(X.LT.0.25D0) THEN
                ROOT=DSQRT(0.25D0-X)
                ETAP=0.5D0+ROOT
                ETAM=0.5D0-ROOT
                RLOG=DLOG(ETAP/ETAM)
                FR=0.5D0*(RLOG**2-PI**2)
                FI=PI*RLOG
                ENDIF
      AR=2.D0*X+X*(4.D0*X-1.D0)*FR
      AI=       X*(4.D0*X-1.D0)*FI
      ABI=3D0*(AR+(0,1)*AI)
      RETURN
      END


      FUNCTION ABISQ(x,y)
      implicit none
      complex *16 abi,t1,t2
      real*8 abisq,x,y
      external ABI
      t1=abi(x)+abi(y)      
      t2=t1*dconjg(t1)
      abisq=dreal(t2)
      return
      end


c
c
c Integration routines
c
c
      subroutine xinteg(sig,indim,initn,incalls,avg,stddev)
      implicit real * 8 (a-h,o-z)
      common/bveg1/xl(10),xu(10),acc,ndim,ncall,itmx,nprn
      common/bveg2/xi(50,10),si,si2,swgt,schi,ndo,it
      common/cwww/w1max,w1vgs,w1evt
      common/ciweight/iweight
      external sig
c
      do j=1,10
        xl(j)=0.d0
        xu(j)=1.d0
      enddo
      nprn=1
      acc=-1.d0
      ndim=indim
      ncall=incalls
      itmx=1
      do j=1,initn
        if(j.eq.1)then
          call vegas(sig,avgi,sd,chi2a)
        else
          call vegas3(sig,avgi,sd,chi2a)
        endif
      enddo
      avg=avgi
      stddev=sd
      return
      end


  


      subroutine strfun(x1,x2,sf)
c This subroutines returns the luminosity for a given initial state
c Assuming one quark flavour to simplify the writing, the values
c returned are 
c
c   jinst=1     -->  sf(1,1) = ( g_h1(x1)*g_h2(x2) )/2
c    (g,g)           sf(1,2) = ( g_h2(x2)*g_h1(x1) )/2
c
c   jinst=2     -->  sf(2,1) = ( q_h1(x1)+qb_h1(x1) )*g_h2(x2)
c    (q,g)           sf(2,2) = ( q_h2(x2)+qb_h2(x2) )*g_h1(x1)
c
c   jinst=3     -->  sf(3,1) = ( q_h1(x1)*qb_h2(x2) + qb_h1(x1)*q_h2(x2) )/2
c    (q,qb)          sf(3,2) = ( q_h2(x2)*qb_h1(x1) + qb_h2(x2)*q_h1(x1) )/2
c
c   jinst=4     -->  sf(4,1) = ( q_h1(x1)*q_h2(x2) + qb_h1(x1)*qb_h2(x2) )/2
c    (q,q)           sf(4,2) = ( q_h2(x2)*q_h1(x1) + qb_h2(x2)*qb_h1(x1) )/2
c
c   jinst=5     -->  sf(5,1) = q_h1(x1)*Qb_h2(x2) + qb_h1(x1)*Q_h2(x2)
c    (q,Qb)          sf(5,2) = q_h2(x2)*Qb_h1(x1) + qb_h2(x2)*Q_h1(x1)
c
c   jinst=6     -->  sf(6,1) = q_h1(x1)*Q_h2(x2) + qb_h1(x1)*Qb_h2(x2)
c    (q,Q)           sf(6,2) = q_h2(x2)*Q_h1(x1) + qb_h2(x2)*Qb_h1(x1)
c
c sf(*,1) is for the event
c sf(*,2) is for the event after permutation of the initial state
c
c x1 and x2 are the Bjorken variables
c
      implicit none
      real * 8 pi,sf,x1,x2
      real * 8 fh1(-5:5),fh2(-5:5)
      parameter(pi=3.14159265358979312D0)
      dimension sf(1:6,1:2)
      integer nl
      common/nf/nl
      integer i,j
      real *8 muf2
      common/fixvar/muf2
      integer isetproton,ih1,ih2
      common/isetproton/isetproton
      common/pdf/ih1,ih2


c      call partonsgrid(muf2,x1,fh1,5,isetproton,ih1)
c      call partonsgrid(muf2,x2,fh2,5,isetproton,ih2)

       call partons(muf2,x1,fh1,5,isetproton,ih1)
       call partons(muf2,x2,fh2,5,isetproton,ih2)
      
c
      do i=1,6
        do j=1,2
          sf(i,j) = 0
        enddo
      enddo

c jinst=1
      sf(1,1) = (fh1(0)*fh2(0))/2
      sf(1,2) = (fh2(0)*fh1(0))/2
c jinst=2
      do i=1,nl
        sf(2,1) = sf(2,1) + (fh1(i)+fh1(-i))*fh2(0)
        sf(2,2) = sf(2,2) + (fh2(i)+fh2(-i))*fh1(0)
      enddo
c jinst=3
      do i=1,nl
        sf(3,1) = sf(3,1) + ( fh1(i)*fh2(-i)
     #                      + fh1(-i)*fh2(i) )/2
        sf(3,2) = sf(3,2) + ( fh2(i)*fh1(-i)
     #                      + fh2(-i)*fh1(i) )/2
      enddo
c jinst=4
      do i=1,nl
        sf(4,1) = sf(4,1) + ( fh1( i)*fh2( i)
     #                      + fh1(-i)*fh2(-i) )/2
        sf(4,2) = sf(4,2) + ( fh2( i)*fh1( i)
     #                      + fh2(-i)*fh1(-i) )/2
      enddo
c jinst=5
      do i=1,nl-1
        do j=i+1,nl
          sf(5,1) = sf(5,1) + fh1( i)*fh2(-j)
     #                      + fh1(-i)*fh2( j)
          sf(5,2) = sf(5,2) + fh2( i)*fh1(-j)
     #                      + fh2(-i)*fh1( j)
        enddo
      enddo
c jinst=6
      do i=1,nl-1
        do j=i+1,nl
          sf(6,1) = sf(6,1) + fh1( i)*fh2( j)
     #                      + fh1(-i)*fh2(-j)
          sf(6,2) = sf(6,2) + fh2( i)*fh1( j)
     #                      + fh2(-i)*fh1(-j)
        enddo
      enddo
      return
      end
c


        FUNCTION LI2(x)                                                     
        implicit none                                                           
*      !! Dilogarithm for arguments x < = 1.0                             
        real*8 X,Y,T,S,A,PI3,PI6,ZERO,ONE,HALF,MALF,MONE,MTWO  
        real*8 C(0:18),H,ALFA,B0,B1,B2,LI2OLD                 
        real*8 Li2                                             
        integer  i                                                       
                                                                           
        DATA ZERO /0.0d0/, ONE /1.0d0/                               
        DATA HALF /0.5d0/, MALF /-0.5d0/                             
        DATA MONE /-1.0d0/, MTWO /-2.0d0/                            
        DATA PI3 /3.289868133696453d0/, PI6 /1.644934066848226d0/    
                                                                           
        DATA C( 0) / 0.4299669356081370d0/                              
        DATA C( 1) / 0.4097598753307711d0/                              
        DATA C( 2) /-0.0185884366501460d0/                              
        DATA C( 3) / 0.0014575108406227d0/                              
        DATA C( 4) /-0.0001430418444234d0/                              
        DATA C( 5) / 0.0000158841554188d0/                              
        DATA C( 6) /-0.0000019078495939d0/                              
        DATA C( 7) / 0.0000002419518085d0/                              
        DATA C( 8) /-0.0000000319334127d0/                              
        DATA C( 9) / 0.0000000043454506d0/                              
        DATA C(10) /-0.0000000006057848d0/                              
        DATA C(11) / 0.0000000000861210d0/                              
        DATA C(12) /-0.0000000000124433d0/                              
        DATA C(13) / 0.0000000000018226d0/                              
        DATA C(14) /-0.0000000000002701d0/                              
        DATA C(15) / 0.0000000000000404d0/                              
        DATA C(16) /-0.0000000000000061d0/                              
        DATA C(17) / 0.0000000000000009d0/                              
        DATA C(18) /-0.0000000000000001d0/                              
                                                                           
        if(x .gt. 1.00000000001d0) then                                    
          write(6,*)'problems in LI2'
          write(6,*)'x=',x 
          stop                                               
        elseif(x .gt. 1.0d0) then                                          
          x = 1.d0                                                      
        endif                                                              
        IF(X .EQ. ONE) THEN                                                
         LI2OLD=PI6
         LI2=LI2OLD                                                       
         RETURN                                                            
        ELSE IF(X .EQ. MONE) THEN                                          
         LI2OLD=MALF*PI6
         LI2=LI2OLD                                                  
         RETURN                                                            
        END IF                                                             
        T=-X                                                               
        IF(T .LE. MTWO) THEN                                               
         Y=MONE/(ONE+T)                                                    
         S=ONE                                                             
         A=-PI3+HALF*(LOG(-T)**2-LOG(ONE+ONE/T)**2)                        
        ELSE IF(T .LT. MONE) THEN                                          
         Y=MONE-T                                                          
         S=MONE                                                            
         A=LOG(-T)                                                         
         A=-PI6+A*(A+LOG(ONE+ONE/T))                                       
        ELSE IF(T .LE. MALF) THEN                                          
         Y=(MONE-T)/T                                                      
         S=ONE                                                             
         A=LOG(-T)                                                         
         A=-PI6+A*(MALF*A+LOG(ONE+T))                                      
        ELSE IF(T .LT. ZERO) THEN                                          
         Y=-T/(ONE+T)                                                      
         S=MONE                                                            
         A=HALF*LOG(ONE+T)**2                                              
        ELSE IF(T .LE. ONE) THEN                                           
         Y=T                                                               
         S=ONE                                                             
         A=ZERO                                                            
        ELSE                                                               
         Y=ONE/T                                                           
         S=MONE                                                            
         A=PI6+HALF*LOG(T)**2                                              
        END IF                                                             
                                                                           
        H=Y+Y-ONE                                                          
        ALFA=H+H                                                           
        B1=ZERO                                                            
        B2=ZERO                                                            
        DO  I = 18,0,-1                                                    
          B0=C(I)+ALFA*B1-B2                                               
          B2=B1                                                            
          B1=B0                                                            
        ENDDO                                                              
        LI2OLD=-(S*(B0-H*B2)+A) 
         LI2=LI2OLD
        end     
c
c                                   
C......Function Li3(x) for -1 < x < 1 (From Daniel)

       function LI3(x)
       implicit none
       double precision LI3,xlog,x,PI,Z3
       PI=3.14159265358979312D0
       Z3=1.20205690315959429D0
  
       if (x.lt.0.5d0) then
       xlog=dlog(1d0-x)
       LI3=  -xlog -(3*xlog**2)/8.-(17*xlog**3)/216.-(5*xlog**4)/576
     #   - (7*xlog**5)/54000. + (7*xlog**6)/86400. + 19*xlog**7/5556600
     #   - xlog**8/752640-11*xlog**9/127008000+11*xlog**10/435456000
       elseif (x.lt.1d0) then
       xlog=dlog(x)
       LI3=Z3+(Pi**2*xlog)/6+(3d0/4-dlog(-xlog)/2)*xlog**2 
     # -xlog**3/12-xlog**4/288+xlog**6/86400-xlog**8/10160640  
       elseif (x.eq.1d0) then
         LI3=1.20205690315959429D0
       else
        write(6,*)'wrong argument of Li3!!' 
       endif
       return
       end   
          

      FUNCTION S12(z)                                                     
      implicit none
      external LI2,LI3
      real *8 S12,z,Z3,LI2,LI3,t
      Z3=1.20205690315959429D0
c
      if(z.gt.0d0) then
      S12=(dlog(1-z)**2*dlog(z)+2*dlog(1-z)*LI2(1-z)
     #-2*LI3(1-z)+2*Z3)/2d0
      elseif(z.eq.0d0) then
      S12=0d0
      elseif(z.lt.0d0) then
      t=-z
      S12=dlog(t)*dlog(1+t)**2/2-dlog(1 + t)**3/3- 
     #dlog(1+t)*LI2(1/(1+t))-LI3(1/(1 + t))+Z3
      endif
      RETURN
      END
