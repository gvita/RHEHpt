      subroutine hqt(CME,higgsmass,rscale,Murin,Mufact,order,pdfsetname,
     .pdfsetnamelen,pdfmember,ptgrid,result,gridsize)

c.....This subroutine computes the Higgs boson qt distribution
c.....in pp collisions at LO and NLO accuracy.
c.....
c.....This subroutine has been derived from the HqT2.0 program by
c.....G. Bozzi, S. Catani, D. de Florian, M. Grazzini
c.....Phys. Lett. B564 (2003) 65 [hep-ph/0302104]
c.....Nucl. Phys. B737 (2006) 73 [hep-ph/0508068]
c.....and included as part of RHEHpt by G.Vita (vita@mit.edu).
c.....In RHEHpt this subroutine is called by a C++ method and
c.....it returns by filling the vector result

      implicit none
c.....Argoument definitions
      real *8 CME,higgsmass,rscale,Murin,Mufact
      integer order,pdfsetnamelen,gridsize,i
      character*64 pdfsetname
      integer pdfmember
	  real* 8 ptgrid (gridsize)
	  real* 8 result (gridsize) ! this is going to be the array where
	  
c.....other defs
      logical NEWA3
      character*64 pdfname
      integer j,ifail,flag,flag1,flag2,flag4,flagalphas,imod
      integer ih1,ih2,inorm,icoll,nset,mem,channel,cicle
      integer isetproton,ippbar,nloop,iord,flagrealcomplex
      real *8 v,fu,fd,max,min,tmp1,tmp2,xlam,phi,abisq
      real *8 errest,adpint,sigma0,gevpb,gf,amt,xmb,amh,factor
      real *8 arg1,arg2,aass
      real *8 integral
      real *8 subqqb,subqg,subgg,a_param,b0p,resummscale
      complex *16 S,ru,rd,rr,jac
      complex *16 alphas
      real *8 allres,expres,ris,matched,switched,g,LT
      real *8 dg2,Gamma2g,dq2,Gamma2q

c.....this is for LL calculation
      complex *16 qqb1,qg1,gg1,LLgg
      common/LLgg/LLgg

      common/flagrealcomplex/flagrealcomplex
      common/aass/aass
      common/v/v
      common/nloop/nloop
      common/order/iord
      common/jac/jac
      common/j/j
      common/isetproton/isetproton
      common/flag/flag
      common/flag1/flag1
      common/flag2/flag2
      common/flag4/flag4
      common/flagalphas/flagalphas
      common/xlam/xlam
      common/ippbar/ippbar
      common/modified/imod
      common/a_param/a_param,b0p
      common/sigma0/sigma0
      common/amt/amt
      common/xmb/xmb
      common/NP/g
      common/nset/nset
      common/LT/LT
      common/channel/channel
      common/pdf/ih1,ih2
      common/pdfname/pdfname
      common/mem/mem
      include 'const.F'
      include 'scales.F'

      external fu,fd,adpint,alphas,S,abisq
      data gevpb/3.8937966d8/


      write (6,*)
      write (6,*)'****************************************************'
      write (6,*)'*                                                  *'
      write (6,*)'*               HqT version 2.0                    *'
      write (6,*)'*                                                  *'
      write (6,*)'****************************************************'

c   Choose pp or ppbar collider

      ih1=1
      ih2=1

c  Choose CM energy
      sroot=CME

c  imod=0 original sudakov
c  imod=1 modified sudakov
      imod=1


c      read (*,*) amh
       amh=higgsmass
       
       
! Change of the evolution scale
      resummscale=rscale
      a_param=amh/resummscale

c.....Renormalization and factorization scales

      mur=Murin
      muf=Mufact
      
      mur2=mur**2
      muf2=muf**2

c......Parameters(amt=mtop, amh=mhiggs, gf=gfermi)
c................(nf=number of active flavours)
c................(flag=0,1 needs to call just once the pdfs fitting routine and evolution code initialization)
c................(flag1=1,2 set the order of resummation: NLL,NNLL)
c................(ippbar=1,-1 sets p-p or p-pbar events)
c................(iord sets the order for the evolution: LO(0),NLO(1))
c................(flag2 sets the order of matching: LO(1), NLO(2))
c................(nloop sets the order of the calculation)
c................(Flag channel: channel that contribuite to the cross-section 1=Tot, 2=GG, 3=GQ, 4=QQ)
c................(NEWA3: Flag for the new/old version of A3)

      amt=172d0
      xmb=4.67d0
      gf=1.16637d-5
      arg1=(amt/amh)**2
      arg2=(xmb/amh)**2
      ippbar=1
      nf=5
      flag=0
      flagalphas=0
      channel=1
      NEWA3=.true.


c     Order of the calculation: (1)-> LO (2)-> NLO

      flag1=order
      
      if (flag1.eq.1) then
         write (6,*) '  Computing LO pt distribution'
         iord=0 ! LO evolution
         flag2=1 ! LO matching
      elseif(flag1.eq.2) then
         write(6,*) '  Computing NLO pt distribution'
         iord=1 ! NLO evolution
         flag2=2 ! NLO matching
      else
  41  format(A)
      write(*,*)'Wrong order !'
      stop
      endif
      write (6,42)'   CM energy=',sroot/1d3,'TeV, collider: pp'
  42  format(A,F8.3,A)

       pdfname=pdfsetname
       mem=pdfmember

! Initializzation for native PDFs, dummy for LHAPFD
      call lambda5qcd(xlam,isetproton)
! Initializzation for LHAPFD, dummy for native PDFs
      call inizializzazPDF(pdfsetnamelen)

c.....Constants
c..............(CA,Cf=color factors)
c..............(beta_i=beta function coefficients defined to aass/pi)
c..............(Euler=Euler's Gamma function)
      pi=dacos(-1d0)
      II=(0d0,1d0)
      ifail=0
      CA=3d0
      Cf=4/3d0
      TR=1/2d0
      Euler=0.57721566d0
      Z2=1.644934d0
      Z3=1.202057d0
      b0=2*dexp(-Euler)

      b0p=b0*a_param
      beta0=(33-2*nf)/12d0
      beta1=(153-19*nf)/24d0
      beta2=2857/128d0-5033*nf/1152d0+325*nf**2/3456d0

c.....Sudakov coefficients defined 
c.....according to Catani-De Florian-Grazzini
c.....WARNING!!!
c.....they are defined according to alphas/pi normalization!!!
c.....Now I added B2-dependence from resummation scheme 
c.....both for Higgs prodution and DY process

c.....choice of resummation scheme

      flag4=1 


c.....Choosing processes and scales

      shad=sroot**2
      q=amh
      q2=q**2
      xtau=q2/shad

      LT=2.0D0*dlog(q/amt)
      
      

c.....gluon coefficients
      A1g=CA
      A2g=CA/2d0*(67/6d0-(pi**2)/2d0-5d0/9*nf)

c.....exact result for A3 for threshold resummation

      if(NEWA3)then
!Thomas Becher & Matthias Neubert arXiv: 1007.4005
!their normalization: aass/4/Pi, our normalization: aass/Pi ->  A3g = A3g/64d0 & Beta0=4d0*Beta0

      Gamma2g = 4d0*CA*(CA**2 *(245d0/6d0 - 134d0*pi**2/27d0
     . + 11d0*Pi**4 /45d0 + 22d0*Z3/3d0) + CA*TR*nf*(-418d0/27d0
     . + 40d0*Pi**2/27d0 - 56d0*Z3/3d0) + CF*TR*nf*(-55d0/3d0 + 16d0
     . * Z3) - 16d0/27d0 * TR**2 * nf**2)

      dg2 = CA*(CA*(808d0/27d0 - 28d0*Z3) - 224d0/27d0*TR*nf)

      A3g = Gamma2g + 2d0 *(4d0*Beta0) * dg2
      A3g = A3g/64d0
      else
      A3g=CA*(13.7683d0-2.14673d0*nf-nf**2/108d0)
      endif

      B1g=-(11*CA-2*nf)/6d0

      if (flag4.eq.1) then
         B2g=CA**2*(23/24d0+(11*pi**2)/18d0-3*Z3/2d0)+
     /        Cf*nf/2d0-CA*nf*(1/12d0+pi**2/9d0)-11/8d0*Cf*CA 
         C1ggn=(pi**2/2d0+11/2d0+pi**2)/2d0
      elseif (flag4.eq.2) then
         B2g=(Cf*nf/2d0+2*CA*nf/3d0-CA**2*(8/3d0+3*Z3))/2d0
         C1ggn=-pi**2/2d0/2d0
      elseif (flag4.eq.3) then
         B2g=(Cf*nf/2d0+2*CA*nf/3d0-CA**2*(8/3d0+3*Z3))/2d0
     /        +pi**2/2d0*beta0
         C1ggn=0d0
      endif


c.....quark coefficients
      A1q=Cf
      A2q=Cf/2d0*(67/6d0-(pi**2)/2d0-5/9d0*nf)

      if(NEWA3)then

      Gamma2q = 4d0*CF*(CA**2 *(245d0/6d0 - 134d0*pi**2/27d0
     . + 11d0*Pi**4 /45d0 + 22d0*Z3/3d0) + CA*TR*nf*(-418d0/27d0
     . + 40d0*Pi**2/27d0 - 56d0*Z3/3d0) + CF*TR*nf*(-55d0/3d0 + 16d0
     . * Z3) - 16d0/27d0 * TR**2 * nf**2)

      dq2 = CF*(CA*(808d0/27d0 - 28d0*Z3) - 224d0/27d0*TR*nf)

      A3q = Gamma2q + 2d0 *(4d0*Beta0) * dq2
      A3q = A3q/64d0
      else
      A3q=Cf*(13.81-2.15*nf-nf**2/108d0)
      endif
      B1q=-(3*Cf)/2d0
      if (flag4.eq.1) then
         B2q=Cf**2*(pi**2/4d0-3/16d0-3*Z3)+
     /        CA*Cf*(11*pi**2/36d0-193/48d0+3*Z3/2d0)+
     /        Cf*nf*(17/24d0-pi**2/18d0)
      elseif (flag4.eq.2) then
         B2q=1/2d0*Cf*nf*(1/6d0+2*pi**2/9d0)-
     /        CA*Cf*(17/24d0+11*pi**2/18d0-3*Z3)-
     /        Cf**2*(3/8d0-pi**2/2d0+6*Z3)
      elseif (flag4.eq.3) then
         B2q=Cf**2*(pi**2/4d0-3/16d0-3*Z3)+
     /        CA*Cf*(-11*pi**2/72d0-13/16d0+3*Z3/2d0)+
     /        Cf*nf*(1/8d0+pi**2/36d0)
      endif


c.....here alphas gives alphas/4pi at the renormalization scale "mur2".
c.....if you need alphas,    you have to write aass=alphas*4*pi
c.....if you need alphas/pi, you have to write aass=alphas*4
c     aass=alphas(cmplx(q2))*4*pih
      write(6,*)
      aass=dreal(alphas(dcmplx(mur2)))*4d0
      write(6,*)
      write(6,44) '   alpha_S(',sqrt(mur2),' GeV) = ',aass*pi
  44  format(A,F7.2,A,F8.6)


c.....Here you choose NP smearing

      g=0d0

c.....mtop->infinity limit:
         factor=gf/288d0/pi/dsqrt(2d0)

c.....Please note here you need aass, not aass/pi or aass/4pi!!!!!!
      sigma0=factor*gevpb*(aass*pi)**2


c.....Here starts the integration: you can choose 
c.....the parameters v,phi,min,max
c.....To integrate on complex plane set phi=pi/8,max=10
c........and set the Sudakov without the b-cutoff
c.....To integrate on real axis set phi=0,max=1
c........and set the Sudakov with the b-cutoff
c     flagrealcomplex=0(1) real(complex)

      flagrealcomplex=0
      flagrealcomplex=1

      if (flagrealcomplex.eq.0) then
C     real axis
         v=1d0
         phi=0d0
         min=1d-5
         max=3d0
      write (6,45) '   Mh= ',amh, 'GeV, real axis'
      write (12,45) '( Mh= ',amh, 'GeV, real axis'
      else
c  complex plane
         v=1d0
         phi=pi/9
         min=1d-5
         max=10d0
         phi=pi/4
         min=1d-5
         max=3d0
      write (6,45) '   Mh= ',amh, 'GeV, complex plane'
      write (12,45) '( Mh= ',amh, 'GeV, complex plane'
      endif
  45  format(A,F7.2,A)
      close(11)


      write (6,45) '   Q=  ',resummscale, 'GeV'
      write (12,45) '( Q=  ',resummscale, 'GeV'
      write (6,46) '   mur=',mur, 'GeV, muf=',muf,'GeV'
      write (12,46) '( mur=',mur, 'GeV, muf=',muf,'GeV'
  46  format(A,F7.2,A,F7.2,A)
      write (6,47) '   g=',g, ', norm=',inorm
      write (12,47) '( g=',g, ', norm=',inorm
  47  format(A,F5.2,A,I3)

c.....this is one-time calculation of luminosities for LL result
      if (flag1.eq.0) then
         call evolfast (dcmplx(muf2),qqb1,qg1,gg1)
         LLgg=gg1
      endif

      do i=1,gridsize 
		 qt = ptgrid(i)
         qt2 = qt**2

         if(flag1.ne.0) then
            call higgsqty(ris,-10d0,10d0,flag2)
            write(6,90)ptgrid(i),allres,ris
         else
            write(6,90) ptgrid(i),allres
         endif
		result(i) = ris
 1    enddo !Steps cicle
  5   format(A,1(F10.5))
  
  90  format(6(1PE15.7))

      end
      
c.....Function to be integrated in the upper complex b plane
      function fu(t)
      implicit none
      integer j,flag1,iord,nloop
      real *8 t,fu
      complex *16 qqb,qg,gg,LLgg
      complex *16 tmp,b,jac,h1,S
      complex *16 final
      common/order/iord
      common/nloop/nloop
      common/jac/jac
      common/j/j
      common/flag1/flag1
      common/LLgg/LLgg

      real *8 a_param,b0p
      common/a_param/a_param,b0p

      external h1,S
      include 'scales.F'
      include 'const.F'
      if (flag1.eq.0) then
         b=jac*t
         tmp=jac*b*h1(qt*b)*S(b)*LLgg
         if(j.eq.1) then
            fu=dreal(tmp)
         elseif(j.eq.2) then
            fu=dimag(tmp)
         endif
      else
         b=jac*t
         final=(b0**2)/(b**2)
         call evolfast(final,qqb,qg,gg)
         tmp=jac*b*h1(qt*b)*S(b)*gg
         if(j.eq.1) then
            fu=dreal(tmp)
         elseif(j.eq.2) then
            fu=dimag(tmp)
         endif
      endif
      return
      end
      
c.....Function to be integrated in the lower complex b plane
      function fd(t)
      implicit none
      integer j,flag1,iord,nloop
      real *8 t,fd
      complex *16 qqb,qg,gg,LLgg
      complex *16 tmp,b,jac,h2,S
      complex *16 final
      common/order/iord
      common/nloop/nloop
      common/jac/jac
      common/j/j
      common/flag1/flag1
      common/LLgg/LLgg

      real *8 a_param,b0p
      common/a_param/a_param,b0p

      external h2,S
      include 'scales.F'
      include 'const.F'
      if (flag1.eq.0) then
         b=jac*t
         tmp=jac*b*h2(qt*b)*S(b)*LLgg
         if(j.eq.1) then
            fd=dreal(tmp)
         elseif(j.eq.2) then
            fd=dimag(tmp)
         endif
      else
         b=jac*t
c         final=(b0p**2)/(b**2)      ! Change of the evolution scale
         final=(b0**2)/(b**2)
         call evolfast(final,qqb,qg,gg)
         tmp=jac*b*h2(qt*b)*S(b)*gg
         if(j.eq.1) then
            fd=dreal(tmp)
         elseif(j.eq.2) then
            fd=dimag(tmp)
         endif
      endif
      return
      end

c.....Sudakov form factor
      function S(b)
      implicit none
      integer flag1,flagrealcomplex,imod
      complex *16 S,f0,f1,f2,b,bstar,blim,blog
      real*8 a_param,b0p,aass,g
      common/aass/aass
      common/flag1/flag1
      common/flagrealcomplex/flagrealcomplex
      common/modified/imod
      common/a_param/a_param,b0p
      common/NP/g
      include 'scales.F'
      include 'const.F'
      blim=(1/q)*exp(1/(2*aass*beta0))

      bstar=b
c.....choose bstar (b) for real axis (complex plane) integration
CC      if (flagrealcomplex.eq.0) bstar=b/sqrt(1+(b**2)/(blim**2))

      if (imod.eq.1) blog=log( (q*bstar/b0p)**2+1d0)    !modified sudakov
      if (imod.eq.0) blog= log( (q*bstar/b0p)**2 )    !normal sudakov

      if (flag1.eq.0) then
         S=exp(blog*f0(beta0*aass*blog))

      elseif (flag1.eq.1) then
         S=exp(blog*f0(beta0*aass*blog) +
     /                           f1(beta0*aass*blog))

      elseif (flag1.eq.2) then
         S=exp(blog*f0(beta0*aass*blog) +
     /                           f1(beta0*aass*blog) +
     /                        aass*f2(beta0*aass*blog))

      endif
      S=S*exp(-g*b**2)

      return
      end
      
      
c.....Soft-gluon-Resummation of LL
      function f0(y)
      implicit none
      complex *16 f0,y
      include 'const.F'
      include 'scales.F'
      f0=(A1g/beta0)*(y+log(1-y))/(y)
      return
      end
      
c.....Soft-gluon-Resummation of NLL
c.....Now we have mu_r dependence!
      function f1(y)
      implicit none
      complex *16 f1,y
      include 'const.F'
      include 'scales.F'
      real *8 a_param,b0p
      common/a_param/a_param,b0p

      f1=((A1g*beta1)/(beta0**3))*((1d0/2)*log(1-y)*log(1-y) +
     \     (y)/(1-y)+log(1-y)/(1-y)) -
     \     (A2g/(beta0**2))*(log(1-y)+(y)/(1-y)) + 
     \     (B1g/beta0)*log(1-y) +
     \     (A1g/beta0)*(y/(1-y)+log(1-y))*log(q2/mur2)


c    a dependence
      
      f1=f1-2*log(a_param)*A1g/beta0*y/(1-y)
      return
      end

c.....Soft-gluon-Resummation of NNLL
c.....Now we have mu_r dependence!
      function f2(y)
      implicit none
      complex *16 f2,y
      include 'const.F'
      include 'scales.F'
      real *8 a_param,b0p
      common/a_param/a_param,b0p

      f2=((A2g*beta1)/(beta0**3))*((y/2)*((3*y-2d0)/(1-y)**2)-
     \     ((1-2*y)*log(1-y)/(1-y)/(1-y))) - 
     \     (B2g/beta0)*((y)/(1-y))+
     \     (B1g*beta1/beta0**2)*((y)/(1-y)+log(1-y)/(1-y))-
     \     (A3g/2/beta0**2)*(y)*(y)/(1-y)/(1-y) +
     \     A1g*((beta1**2/2/beta0**4)*(1-2*y)/(1-y)
     \     /(1-y)*log(1-y)*log(1-y) +
     \     log(1-y)*((beta0*beta2-beta1**2)/(beta0**4)+
     \     beta1**2/beta0**4/(1-y)) +
     \     (y)/(2*beta0**4*(1-y)*(1-y))*
     \     (beta0*beta2*(2d0-3*y)+beta1**2*y)) -
     \     (A1g/2)*(y)*(y)/(1-y)/(1-y)*log(q2/mur2)*log(q2/mur2) +
     \     log(q2/mur2)*(B1g*y/(1-y)+A2g/beta0*y*y/(1-y)/(1-y)+
     \     A1g*beta1/beta0**2*(y/(1-y)+(1-2*y)/(1-y)/(1-y)*log(1-y))) 
     \     +2d0*C1ggn*((y)/(1-y))

      f2=f2+2*A1g*y*(y-2)/(1-y)**2*log(a_param)**2-log(a_param)
     \ *(2*B1g*y/(1-y)+2*y/beta0*A2g/(1-y)**2
     \ -2*A1g*beta1/beta0**2*y*log(1-y)/(1-y)**2)
     \ +2d0*A1g*y*log(q2/mur2)*log(a_param)/(1-y)**2    !Double log-dependent term
      return
      end
      

