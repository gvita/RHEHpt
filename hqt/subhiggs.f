      subroutine subtraction(factscale,subqqb,subqg,subgg)

c.....computes the expansion of the resummed part to 1st or 2nd order
c.....minor bug corrected in sigma2gg

c.....January 2006: Bug corrected in xlog3

c.....May 2010: a dependence for general muf and mur added


      implicit double precision (A-H,O-Z)
      integer f,flag,flag2,flag4
      integer iord,imod
      real *8 q2s,lambda
      real*8 IK
      complex *16 alphas,all,factscale
      common/flag/flag
      common/flag2/flag2
      common/flag4/flag4
      common/isetproton/isetproton
      common/order/iord 
      common/scales/q2s,lambda
      common/nloop/nloop
!      common/xlam/xlam
      common/nflavors/f
      common/alfa/nll
      common/idistr/idist
      common/ippbar/ippbar
      common/logs/xlog0,xlog1,xlog2,xlog3
      common/modified/imod
      external xqq,xqg,xgg
 
      include 'scales.F'
 
      Z3=1.20205690316d0

      if (imod.eq.1) then
            xlog0=(-qt2/2)* IK(1)
            xlog1=(-qt2/4) * Ik(2)
            xlog2=(-qt2/6) * IK(3)
            xlog3= -(qt2/8)   * IK(4) + 4*Z3*xlog0 
      elseif (imod.eq.0) then
            xlog0=1d0
            xlog1=DLOG(q2/qt2)
            xlog2=(DLOG(q2/qt2))**2
            xlog3=(DLOG(q2/qt2))**3
      else
           write(6,*)'wrong imod= ',imod
           stop
      endif
      
c.....program disadvantage: single precision!
      q2s=sngl(muf2)
c.....compute the expansion of the resummed part

      call sub (xtau,factscale,
     /     uv,dv,us,ds,ss,gl,subqqb,subqg,subgg)
      return
      end
      
      subroutine sub (x,q2,uv,dv,us,ds,ss,gl,subqqb,subqg,subgg)
      implicit double precision (a-h, o-z) 
      complex *16 q2
      integer iord,flag2,flag4
      dimension pa(9)  
      common/order/iord 
      common/flag2/flag2
      common/flag4/flag4
c.....call polarized parton distributions

      call subinv(pa,x,q2)
c.....*** u,d - valence 
      uv=pa(1)
      dv=pa(2)
c.....*** u,d - sea
      us=pa(3)
      ds=pa(4)
c.....*** s - sea
      ss=pa(5)
c.....*** gluon
      gl=pa(6)
c.....*** subtractions
      subqqb=pa(7) 
      subqg=pa(8) 
      subgg=pa(9) 
      end
      
*     Mellin inversion of the SIGMAs contributions.
*     The integration parameters (for the contour in the complex n-plane
*     and for the Gauss quadrature) are taken from the common-blocks
*     CONT, MOMS and WEIGHTS. 
                  
      subroutine subinv(pa,xb,q2)    
      implicit real *8 (a-z)
      integer flag2,flag4
      double precision c,ex,pa(9),xb,wn(136) 
      complex *16 cc(2),xnm(2),cdex(2),n(2,136),
     .     fn(2,137,9),q2,alps,alpq,alphas,II
      integer nmax, i1, i2, m, iord, kp
      common/scales/q2s,lambda
      common/coupl/alps,alpq
      common/contin/c,cc                              
      common/weights/wn
      common/moms/n
      common/order/iord
      common/flag2/flag2
      common/flag4/flag4
      II=cmplx(0.0,1.0)
*...  coupling constants at static point, input scale :
!      alps = alphas (cmplx(q2s,0d0))
!      alpq = alphas (q2)
*...  q2 and x dependent quantities : 
      x  = xb                                    
      ax = dlog (x)                                        
      ex = dble (exp (-c * ax))                                  
*...  integration length parameter for mellin inversion
      if ( X .le. 0.001 ) then                          
         nmax = 24              ! zmax = 2 
      else if ( X .le. 0.05 ) then                   
         nmax = 40              ! zmax = 4  
      else if ( X .le. 0.2 ) then                   
         nmax = 56              ! zmax = 8  
      else if ( X .le. 0.4 ) then                   
         nmax = 72              ! zmax = 12 
      else if ( X .le. 0.7 ) then                   
         nmax = 88              ! zmax = 18 
      else                                                          
         nmax = 136             ! zmax = 36
      end if              
c      nmax=136
      nmax=136
      xxx=1                               
*...  calculation of the parton densities and output:
*     index of func :   1     2    3     4    5   6   
*     distribution  :  UV    DV    UB    DB   S   G  
      call subreno (fn, n, nmax,iord)
      do 2 i2 = 1, 9
         fun = 0.d0
         if(i2.ge.7) xxx=x
         do 1 i1 = 1, nmax
            do 3 kp=1,2
               xnm(kp) = (c - n(kp,i1)+1.) * ax
               cdex(kp) = cdexp (xnm(kp)) / 3.1415927 * cc(kp)/2./II
 3          continue
c     fz = dble (aimag (fn(kp,i1+1,i2)*cex(kp)))
           fz = dble (((fn(1,i1+1,i2)*cdex(1) - fn(2,i1+1,i2)*cdex(2))))
           fun = fun + wn(i1) * fz
 1       continue
         pa(i2) = fun * ex/xxx
 2    continue
      return                         
      end
      
*...  Mellin-n space Q**2 - evolution of parton densities and structure 
*     functions in next to leading order :
*     (for valence input scale not equal sea/gluon input scale)
*     Output for FN: Array of parton moments 
*     2. index of FN :   1     2      3        4      5   6   7   8 
*     Distribution   :  uv    dv    db-ub    db+ub    s   g  G1P G1N
      
      subroutine subreno (fn, n, nmax, iord)
      implicit complex*16 (a - z)
      dimension FX1(-5:5), FX2(-5:5)
      real *8 beta0N,beta1N,beta2N,A1gN,A2gN,A3gN,B1gN,B2gN,
     /     A1qN,A2qN,A3qN,B1qN,B2qN
      real*8 aass,xlog0,xlog1,xlog2,xlog3
      dimension ans(2,136), am(2,136), ap(2,136), al(2,136),
     1     be(2,136), ab(2,136), ac(2,136), rmin(2,136),
     2     rplus(2,136), C1Q(2,136), C1G(2,136),C3Q(2,136), 
     3     rmmqq(2,136),rmmqg(2,136),rmmgq(2,136) ,rmmgg(2,136),
     4     rmpqq(2,136),rmpqg(2,136),rmpgq(2,136),rmpgg(2,136),
     5     rpmqq(2,136),rpmqg(2,136),rpmgq(2,136),rpmgg(2,136),
     6     rppqq(2,136),rppqg(2,136),rppgq(2,136),rppgg(2,136),
     7     uvi(2,136), dvi(2,136), usi(2,136), dsi(2,136),
     8     ssi(2,136), gli(2,136), chi(2,136), boi(2,136),
     9     fn(2,137,9), n(2,136)
      integer  F, K1, KK1,NMAX, iord, iiord, KP, flag2, flag4
      integer i,j,ippbar
      integer channel
      common/channel/channel
      common / ANOMS  / ANS, AM, AP, AL, BE, AB, AC, RMIN, RPLUS, C1Q, 
     1     C1G, C3Q, RMMQQ, RMMQG, RMMGQ, RMMGG, RMPQQ, RMPQG, RMPGQ, 
     2     RMPGG, RPMQQ, RPMQG, RPMGQ, RPMGG, RPPQQ, RPPQG, RPPGQ, 
     3     RPPGG
      common / INPUT  / UVI, DVI, USI, DSI, SSI, CHI, BOI, GLI
      common / order/ iiord
      common / COUPL  / ALPS,  ALPQ
      common/flag2/flag2
      common/flag4/flag4
      common/ippbar/ippbar
      common/aass/aass
      common/logs/xlog0,xlog1,xlog2,xlog3
      include 'const.F'
      include 'scales.F'

      real*8 a_param,b0p
      common/a_param/a_param,b0p

c.....this is because SIGMA IS NORMALIZED TO AS/2*PI
c.....while our coefficients are normalized to AS/PI
      beta0N=beta0*2d0
      beta1N=beta1*4d0
      beta2N=beta2*8d0
      A1gN=A1g*2d0
      A2gN=A2g*4d0
      A3gN=A3g*8d0
      B1gN=B1g*2d0
      B2gN=B2g*4d0
      A1qN=A1q*2d0
      A2qN=A2q*4d0
      A3qN=A3q*8d0
      B1qN=B1q*2d0
      B2qN=B2q*4d0

      iord=iiord
      oswi=cmplx(Real(iord),0.)
      silly=0
      DO 1 KP=1,2
      DO 1 KK1 = 0, NMAX
         IF(KK1.EQ.0) THEN
            K1=1
         ELSE
            K1=KK1
         ENDIF 
         
      UVN = UVI(KP,K1)
      DVN = DVI(KP,K1)
      NS3N = UVI(KP,K1) + 2*USI(KP,K1) - DVI(KP,K1) - 2.*DSI(KP,K1)
      NS8N = UVI(KP,K1) + 2*USI(KP,K1) + DVI(KP,K1) + 2.*DSI(KP,K1)
     .     - 4.*SSI(KP,K1)
      GLN = GLI(KP,K1)
      SIN = UVI(KP,K1) + DVI(KP,K1) + 2*USI(KP,K1) + 2*DSI(KP,K1)
     .     + 2*SSI(KP,K1) + 2*CHI(KP,K1) + 2*BOI(KP,K1)
      NS15N = UVI(KP,K1) + DVI(KP,K1) + 2*USI(KP,K1) + 2*DSI(KP,K1) +
     .     2*SSI(KP,K1) - 6*CHI(KP,K1)
      NS24N = UVI(KP,K1) + DVI(KP,K1) + 2*USI(KP,K1) + 2*DSI(KP,K1) +
     .     2*SSI(KP,K1) + 2*CHI(KP,K1) - 8*BOI(KP,K1)
      NS35N = SIN
      F = 5

C...  FLAVOUR DECOMPOSITION OF THE QUARK SEA :
      SSN = (10.* SIN + 2.* NS35N + 3.* NS24N + 5.* NS15N - 20.* NS8N)
     1     / 120.
      DSN = (10.* SIN + 2.* NS35N + 3.* NS24N + 5.* NS15N + 10.* NS8N
     1     - 30.* NS3N - 60.* DVN) / 120.
      USN = (10.* SIN + 2.* NS35N + 3.* NS24N + 5.* NS15N + 10.* NS8N
     1     + 30.* NS3N - 60.* UVN) / 120.
      CHN = (UVN + DVN + 2*USN + 2*DSN + 2*SSN - NS15N)/6.
      BON = (UVN + DVN + 2*USN + 2*DSN + 2*SSN + 2*CHN - NS24N)/8.
      
      if (nf.eq.3) then        !GRV
         SSN= (20.* SIN - 20.* NS8N)/120.
         DSN = (20.* SIN + 10.* NS8N - 30.* NS3N - 60.* DVN) / 120.
         USN = (20.* SIN + 10.* NS8N + 30.* NS3N - 60.* UVN) / 120.
         CHN=CMPLX(0D0,0D0)
         BON=CMPLX(0D0,0D0)
      endif



c     NEW: Correct luminosities



c h1=proton

      FX1(0)=GLN
      FX1(1)=UVN+USN
      FX1(2)=DVN+DSN
      FX1(3)=SSN
      IF(NF.GE.4) FX1(4)=CHN
      IF(NF.GE.5) FX1(5)=BON
      FX1(-1)=USN
      FX1(-2)=DSN
      FX1(-3)=SSN
      IF(NF.GE.4) FX1(-4)=CHN
      IF(NF.GE.5) FX1(-5)=BON

      if(ippbar.eq.1) then
c h2=proton
       FX2(0)=GLN
       FX2(1)=UVN+USN
       FX2(2)=DVN+DSN
       FX2(3)=SSN
       IF(NF.GE.4) FX2(4)=CHN
       IF(NF.GE.5) FX2(5)=BON
       FX2(-1)=USN
       FX2(-2)=DSN
       FX2(-3)=SSN
       IF(NF.GE.4) FX2(-4)=CHN
       IF(NF.GE.5) FX2(-5)=BON
      elseif(ippbar.eq.-1) then
c h2=antiproton
       FX2(0)=GLN
       FX2(-1)=UVN+USN
       FX2(-2)=DVN+DSN
       FX2(-3)=SSN
       IF(NF.GE.4) FX2(-4)=CHN
       IF(NF.GE.5) FX2(-5)=BON
       FX2(1)=USN
       FX2(2)=DSN
       FX2(3)=SSN
       IF(NF.GE.4) FX2(4)=CHN
       IF(NF.GE.5) FX2(5)=BON
      endif 
    

       GGN=0
       QGN=0
       GQN=0
       QQBN=0
       QQN=0
       QQPBN1=0
       QQPBN2=0
       QQPN1=0
       QQPN2=0

c      gg

       GGN=FX1(0)*FX2(0)


c      qg

       do i=1,nf
       QGN=QGN+(FX1(i)+FX1(-i))*FX2(0)
       enddo

c      gq

       do i=1,nf
       GQN=GQN+(FX2(i)+FX2(-i))*FX1(0)
       enddo
  
c      qqb

       do i=1,nf
       QQBN=QQBN+(FX1(i)*FX2(-i)+FX1(-i)*FX2(i))
       enddo


c      qq


       do i=1,nf
        QQN=QQN + ( FX1(i)*FX2(i)+ FX1(-i)*FX2(-i) )
       enddo
       
c      qqpb


       do i=1,nf-1
        do j=i+1,nf
           QQPBN1 = QQPBN1 + FX1( i)*FX2(-j)
     #                      + FX1(-i)*FX2( j)
           QQPBN2 = QQPBN2  + FX2( i)*FX1(-j)
     #                      + FX2(-i)*FX1( j)
        enddo
       enddo

       QQPBN=QQPBN1+QQPBN2

c      qqp

       
       do i=1,nf-1
        do j=i+1,nf
          QQPN1 = QQPN1 + FX1(i)*FX2(j)
     #                      + FX1(-i)*FX2(-j)
          QQPN2 = QQPN2 + FX2(i)*FX1(j)
     #                      + FX2(-i)*FX1(-j)
        enddo
       enddo

       QQPN=QQPN1+QQPN2


C.....Coefficient functions and anomalous dimensions
c.....normalized to alphas/2*pi      
      XN=N(KP,K1)

      C1QQ=2*pi**2/9d0+2/3d0*(2*pi**2/3d0-8d0) + 
     /     4/3d0/(XN*(XN+1))
      C1QG=1d0/((XN+1)*(XN+2))
      C1GQ=4/3d0/(XN+1)
      if (flag4.eq.1) then
         C1qq=2*pi**2/3d0-16/3d0+4/3d0/(N(KP,K1)*(N(KP,K1)+1))
         C1gg=pi**2/2d0+11/2d0+pi**2
         H1g=0d0
         H1q=0d0
      elseif (flag4.eq.2) then
         C1qq=-2/9d0*pi**2+4/3d0/(N(KP,K1)*(N(KP,K1)+1))
         C1gg=-pi**2/2d0
         H1g=11d0+4*pi**2
         H1q=16/3d0*(pi**2/3d0-2d0)
      elseif (flag4.eq.3) then
         C1qq=-2/3d0+4/3d0/(N(KP,K1)*(N(KP,K1)+1))
         C1gg=0d0
         H1g=11d0+3*pi**2
         H1q=4/3d0*(pi**2-7d0)
      endif

      call ancalc(QQIN, QGFN, GQIN, GGIN, GGFN, NS1MIN, NS1PIN, NS1FN,
     /     QQ1FN, QG1FN, GQ1IN, GQ1FN, GG1IN, GG1FN, XN)

      gamma1qq=-1d0*(QQIN/4d0)
      gamma1qg=-1d0*(QGFN/8d0)
      gamma1gq=-1d0*(GQIN/4d0)
      gamma1gg=-1d0*((GGIN+nf*GGFN)/4d0)
      gamma2qq=-1d0*(((NS1PIN+nf*NS1FN)+nf*QQ1FN)/8d0)
      gamma2qg=-1d0*(QG1FN/16d0)
      gamma2gq=-1d0*((GQ1IN+nf*GQ1FN)/8d0)
      gamma2gg=-1d0*((GG1IN+nf*GG1FN)/8d0)
      
C.....SIGMA1 and SIGMA2
c.....The gq and qq contributions are taken from Arnold-Kaufmann
c.....but have to be DIVIDED BY 2 in the end because of 
c.....luminosities normalization!


c.....NOTE: the factor 3 in fron of beta0*DLOG(mur2/q2) is due to the
c.....fact that the Born is alphas^2

      SIGMA1gg=A1gN*xlog1+ (B1gN+2*gamma1gg)*xlog0
      SIGMA1gq=(2*gamma1gq)*xlog0
      SIGMA2gg=(-1/2d0*A1gN**2)*xlog3 +
     /       (A1gN*(beta0N-3/2d0*B1gN)-3*A1gN*gamma1gg)*xlog2 +
     /       (A2gN+beta0N*B1gN-B1gN**2+3*beta0N*A1gN*DLOG(mur2/q2)+
     /       2*A1gN*C1GG-4*gamma1gg**2-4*nf*gamma1gq*gamma1qg+
     /      (beta0N-2*B1gN-A1gN*DLOG(muf2/q2))*2*gamma1gg)* xlog1 +
     /      ( (B2gN+2*gamma2gg+2*(B1gN+2*gamma1gg)*C1GG+2*Z3*A1gN**2-
     /       2*beta0N*C1GG+4*nf*C1GQ*gamma1qg+
     /       3*beta0N*B1gN*DLOG(mur2/q2)+(3*beta0N*DLOG(mur2/q2)-
     /       B1gN*DLOG(muf2/q2))*2*gamma1gg-(2*gamma1gg**2+
c    /       2*nf*gamma1gq*gamma1qg)*2*DLOG(muf2/q2))  +
c    /       DLOG(q2/qt2)*(A1gN*H1g) + H1g*(B1gN+2*gamma1gg))*xlog0
     /       2*nf*gamma1gq*gamma1qg)*2*DLOG(muf2/q2)))*xlog0 +
     /       H1g*(A1gN*xlog1+(B1gN+2*gamma1gg)*xlog0)



      SIGMA2gq=-3*A1gN*gamma1gq*xlog2 +
     /       ((beta0N-2*B1gN-A1gN*DLOG(muf2/q2))*2*gamma1gq+
     /       2*A1gN*C1GQ-6*gamma1gg*gamma1gq-2*gamma1gq*gamma1qq)*
     /      xlog1+
     /      ( (2*gamma1gg*C1GQ+4*gamma1gq*C1GG+2*gamma1qq*C1GQ+
     /       2*gamma2gq+(B1gN-beta0N)*2*C1GQ+
     /       2*(3*beta0N*DLOG(mur2/q2)-B1gN*DLOG(muf2/q2))*gamma1gq-
     /       2*(3*gamma1gg*gamma1gq+gamma1gq*gamma1qq)*DLOG(muf2/q2)) 
     /       + H1g*2*gamma1gq )*xlog0


      SIGMA2qq=-2*gamma1gq**2* xlog1 +
     /       ((2*gamma1gq*C1GQ-2*gamma1gq**2*DLOG(muf2/q2)))*xlog0



c.....NOW a DEPENDENCE (muf=mur=mH)!


      SIGMA1gg=Sigma1gg+2*A1gN*log(a_param)*xlog0
      SIGMA1gq=SIGMA1gq


      SIGMA2gg=SIGMA2gg-3*A1gN**2*log(a_param)*xlog2
     /     -6*A1gN**2*log(a_param)**2*xlog1
     /     +(-6*A1gN*(B1gN+2*gamma1gg)+4*A1gN*beta0N)*xlog1*log(a_param)
     /     +(-4*A1gN**2*log(a_param)**3+log(a_param)**2
     /     *(-6*A1gN*(B1gN+2*gamma1gg)+4*A1gN*beta0N)
     /     +2*log(a_param)*(A2gN+beta0N*(B1gN+2*gamma1gg)
     /     -(B1gN+2*gamma1gg)**2+2*A1gN*C1gg
     /     -2*2*nf*gamma1gq*gamma1qg))*xlog0


      SIGMA2gq=SIGMA2gq+2d0*(-6*A1gN*gamma1gq*log(a_param)*xlog1
     /    +(-6*A1gN*gamma1gq*log(a_param)**2
     /    +log(a_param)*(2*A1g*C1gq+2*gamma1gq*(beta0N-2*B1gN)
     /                 -6*gamma1gg*gamma1gq-2*gamma1gq*gamma1qq))*xlog0)



      SIGMA2qq=SIGMA2qq-4*gamma1gq**2*log(a_param)*xlog0


c.....Finally add a dependence for muf and/or mur different from mH

      sigma2gg=sigma2gg-4*xlog0*A1gN*log(a_param)*
     &     (0.5d0*beta0N*dlog(q2/mur2)
     &      -(gamma1gg*log(q2/muf2)-beta0N*dlog(q2/mur2)))


      sigma2gq=sigma2gq-8*xlog0*
     &     (-0.5d0*A1gN*log(a_param)*gamma1gq*dlog(q2/muf2))



! Channels
!         write(*,*) 'subhiggschannel=',channel
      if(channel.eq.1)then
       goto 444
      elseif(channel.eq.2)then ! GG channel
       SIGMA1gq=dcmplx(0d0,0d0)
       SIGMA2gq=dcmplx(0d0,0d0)
       SIGMA2qq=dcmplx(0d0,0d0)
      elseif(channel.eq.3)then ! GQ channel
       SIGMA1gg=dcmplx(0d0,0d0)
       SIGMA2gg=dcmplx(0d0,0d0)
       SIGMA2qq=dcmplx(0d0,0d0)
      elseif(channel.eq.4)then ! QQ channel
       SIGMA1gg=dcmplx(0d0,0d0)
       SIGMA1gq=dcmplx(0d0,0d0)
       SIGMA2gg=dcmplx(0d0,0d0)
       SIGMA2gq=dcmplx(0d0,0d0)
      endif
 444  continue


*.....OUTPUT
c.....See above for the 1/2 factor in gq and qq contributions

      FN(KP,KK1+1,1) = UVN
      FN(KP,KK1+1,2) = DVN   
      FN(KP,KK1+1,3) = USN
      FN(KP,KK1+1,4) = DSN
      FN(KP,KK1+1,5) = SSN
      FN(KP,KK1+1,6) = GLN
      
      alpr=aass/2d0

      if (flag2.eq.1) then 
         FN(KP,KK1+1,7) = 0d0
         FN(KP,KK1+1,8) = (QGN+GQN)*(ALPr*SIGMA1gq/2d0)
         FN(KP,KK1+1,9) = GGN*(ALPr*SIGMA1gg)
      elseif (flag2.eq.2) then  
         FN(KP,KK1+1,7) = (QQBN+QQN+QQPBN+QQPN)
     /                         *((ALPr)**2*SIGMA2qq)
         FN(KP,KK1+1,8) = (QGN+GQN)*((ALPr*SIGMA1gq/2d0)+
     /                        ((ALPr)**2*SIGMA2gq/2d0))
         FN(KP,KK1+1,9) = GGN*((ALPr*SIGMA1gg)+
     /                        ((ALPr)**2*SIGMA2gg))
      endif
      
 1    CONTINUE
      RETURN
      END
      




