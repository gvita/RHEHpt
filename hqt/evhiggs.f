CC    New: May 2010

CC    Now evolution from resummation scale to b0^2/b^2

CC    Now the variable final must be b0^2/b^2 (to be set in auto.f)

CC    Full exponentiation

      subroutine evolfast(final,qqb,qg,gg)
      implicit double precision (A-H,O-Z)
      integer f,flag,flag1,flag4
      integer iord,xi
      real*8 q2s,lambda
      complex *16 alphas,all,final,ff
      complex *16 qqb,qg,gg,uv,dv,us,ds,ss,gl
      common/flag/flag
      common/flag1/flag1
      common/flag4/flag4
      common/isetproton/isetproton
      common/order/iord 
      common/CUV/A1UV,A2UV,A3UV,A4UV,A5UV,A6UV,A7UV,A8UV
      common/CDV/A1DV,A2DV,A3DV,A4DV,A5DV,A6DV,A7DV,A8DV
      common/CUS/A1US,A2US,A3US,A4US,A5US,A6US,A7US,A8US
      common/CDS/A1DS,A2DS,A3DS,A4DS,A5DS,A6DS,A7DS,A8DS
      common/CSS/A1SS,A2SS,A3SS,A4SS,A5SS,A6SS,A7SS,A8SS
      common/CGL/A1GL,A2GL,A3GL,A4GL,A5GL,A6GL,A7GL,A8GL
      common/CCH/A1CH,A2CH,A3CH,A4CH,A5CH,A6CH,A7CH,A8CH
      common/CBO/A1BO,A2BO,A3BO,A4BO,A5BO,A6BO,A7BO,A8BO
      save/CUV/
      save/CDV/
      save/CUS/
      save/CDS/
      save/CSS/
      save/CGL/
      save/CCH/
      save/CBO/
      common/scales/q2s,lambda
      common/scheme/N
      common/nloop/nloop
      common/xlam/xlam
      common/nf/f
      common/idistr/idist
c     common/x/x,q2,q20
      common/ippbar/ippbar
!
      common/a_param/a_param,b0p
!      external xqq,xqg,xgg
      include 'scales.F'
      
      data N/1/              !MSBAR SCHEME


c.....NEW: Start evolution from resummation scale


      q2s=q2/a_param**2


C.....Here call the program to fit the pdfs at q2 around 'x'

      xtauf=muf2/shad

      if (flag.eq.0) then
         call fiter(1,xtauf,muf2,
     /        A1UV,A2UV,A3UV,A4UV,A5UV,A6UV,A7UV,A8UV)
         call fiter(2,xtauf,muf2,
     /        A1DV,A2DV,A3DV,A4DV,A5DV,A6DV,A7DV,A8DV)
         call fiter(3,xtauf,muf2,
     /        A1US,A2US,A3US,A4US,A5US,A6US,A7US,A8US)
         call fiter(4,xtauf,muf2,
     /        A1DS,A2DS,A3DS,A4DS,A5DS,A6DS,A7DS,A8DS)
         call fiter(5,xtauf,muf2,
     /        A1SS,A2SS,A3SS,A4SS,A5SS,A6SS,A7SS,A8SS)
         call fiter(6,xtauf,muf2,
     /        A1GL,A2GL,A3GL,A4GL,A5GL,A6GL,A7GL,A8GL)
         call fiter(7,xtauf,muf2,
     /        A1CH,A2CH,A3CH,A4CH,A5CH,A6CH,A7CH,A8CH)
         call fiter(8,xtauf,muf2,
     /        A1BO,A2BO,A3BO,A4BO,A5BO,A6BO,A7BO,A8BO)

c.....NEW: Make grid

c      call gridmaker

C.....Initialize evolution code
      call init 
      call inp200 
C   CHECK FOR FIT TO PDF:
!      do xi=1,100,1
!      xx=sqrt(xtau)*real(xi)/10d0
!      cola=muf2**(0.5)
!      ipp=1
!      ff=muf2
!      call  sub (xx,ff,xuv,xdv,xus,xds,xss,xgl,subqqb,subqg,subgg)
!      call distri(xx,cola,upv,dnv,usea,dsea,str,chm,bot,glu,ipp)
!      write(6,*)xi,xx,xgl/glu,xus/usea
      write(*,*)'Fit completed'
      endif
      flag=1 

c.....Compute distributions/luminosities at Q2 from evolution from Q20
      call evol (xtau,final,uv,dv,us,ds,ss,gl,qqb,qg,gg)



      return
      end
      
      subroutine fiter(idist,x,q2,a1,a2,a3,a4,a5,a6,a7,a8)
      implicit double precision (A-H,O-Z)
      external fcng1,chi2
      dimension nprm(8), vstrt(8),stp(8),bl(8),bu(8),arglis(16) 
      dimension parsal(8)
      character*10 pnam(8) 
      data nprm /1,2,3,4,5,6,7,8/  
      data pnam  / 'A1','A2','A3','A4','A5', 'A6', 'A7', 'A8'/
      data vstrt /  1D0, -.1D0, 5D0, 3D0, -1D0, 0D0, 0D0, 0D0/
      data stp   /  1D0, 1D0, 1D0, 10D0, 15D0, 25D0, 25D0, 25D0/
c      data bl    /0D0, -0.8D0, 4D0, -20D0, -40D0, -80D0,-80D0,-80D0/
c      data bu    /0D0, 1.5D0, 19D0, 600D0, 40D0, 80D0, 80D0, 80D0/

c      data bl    /0D0, -0.8D0, 5D0, -5D0, -30D0, -30D0,-30D0,-30D0/
c      data bu    /0D0, 1.2D0, 12D0, 50D0, 30D0, 30D0, 30D0, 30D0/

      data bl    /0D0, -0.8D0, 4D0, -20D0, -60D0, -90D0,-90D0,-90D0/
      data bu    /0D0, 1.5D0, 25D0, 600D0, 60D0, 90D0, 90D0, 90D0/

      common/sal/parsal
      common/distribucion/id
      common/xq2/rx,rq2
      id=idist
      rx=x
      rq2=q2
c.....initialization :
      call mninit(5,1,1)
c     call mninit(5,6,6)
c.....definitions of the parameters :
      do 11 i = 1, 8
         call mnparm (nprm(i),pnam(i),VSTRT(i),STP(i),BL(i),BU(i),
     .        ierflg,chi2)
         if (ierflg .ne. 0) then
            write (6,*) ' unable to define parameter no.', i
            stop
         end if
 11   continue
c.....output  
      arglis(1) = 0.
      call mnexcm(fcng1,'set print',arglis,1,ierflg,chi2)
c.....first call :
      arglis(1) = 1.            !   IFLAG = 1 
      call mnexcm(fcng1,'CALL FCN',arglis,1,ierflg,chi2)
c.....simplex fit :
   
c      goto 100

      arglis(1)=5000.
      call mnexcm(fcng1,'SIMPLEX',arglis,1,ierflg,chi2)
      call mnexcm(fcng1,'MIGRAD',arglis,1,ierflg,chi2)
      
      arglis(1)=20000.
      call mnexcm(fcng1,'SIMPLEX',arglis,1,ierflg,chi2)
      call mnexcm(fcng1,'MIGRAD',arglis,1,ierflg,chi2)
      arglis(1)=28000.
      call mnexcm(fcng1,'SIMPLEX',arglis,1,ierflg,chi2)
      call mnexcm(fcng1,'MIGRAD',arglis,1,ierflg,chi2)
 100  continue     

c      arglis(1)=500.
c      call mnexcm(fcng1,'SIMPLEX',arglis,1,ierflg,chi2)
c      call mnexcm(fcng1,'MIGRAD',arglis,1,ierflg,chi2)

cc.....last call :
      arglis(1) = 3             !   iflag = 3
      call mnexcm (fcng1, 'call fcn', arglis, 1, ierflg,chi2)
c.....stop :
      call mnexcm (fcng1,'stop',arglis,1,ierflg,chi2)
      a1=parsal(1)
      a2=parsal(2)
      a3=parsal(3)
      a4=parsal(4)
      a5=parsal(5)
      a6=parsal(6)
      a7=parsal(7)
      a8=parsal(8)
 1200 format(4F8.4)

      
c      call distri(rx,q,upv,dnv,usea,dsea,str,chm,bot,glu,ippbar)
c      call evol (rx,cmplx(q*q),uv,dv,us,ds,ss,gl,qqb,qg,gg)
c       write(6,*)rx,glu/gl
      

      return
      end
      
      subroutine fcng1 (npar, g, f, x, iflag ,chi2)
      implicit double precision (a-h, o-z)
      dimension x(*), g(*)
      external chi2 
      f = chi2 (x)
      return
      end
      
c     NEW version of chi2 from Daniel

      
      double precision function chi2(param)
      implicit double precision (a-h,o-z)
      dimension param(8),xx1(222)!,xx1(63),xx1(74)
      dimension parsal(8)
      common/sal/parsal
      common/distribucion/id
      common/xq2/rx,rq2
      common/expp/aa
c      data npoints/63/
!      data npoints/74/
      data npoints/222/
c      data xx1/1.d-5,2.d-5,4.d-5,6.d-5,8.d-5,9D-5,
c     .     1.D-4,2.D-4,3.D-4,4.D-4,5.D-4,6.D-4,8.D-4,
c     .     1.D-3,2.D-3,3D-3,4.D-3,5.D-3,6.D-3,7D-3,8d-3,9d-3,
c     .     1.D-2,2.D-2,3D-2,4.D-2,5D-2,
c     .     6.D-2,6.5d-2,7D-2,7.5d-2,8.D-2,8.5d-2,9d-2,9.5d-2,
c     .     .1D0,.11d0,.125D0,.15D0,.175D0,.2D0,.225D0,.25D0,.275D0,
c     .     .3D0,.325D0,.35D0,.375D0,.4D0,.425D0,.45D0,.475D0,
c     .     .5D0,.525D0,.55D0,.575D0,.6D0,.65D0,.7D0,.75D0,
c     .     .8D0,.9D0,1.D0/


!      data xx1/1.d-5,2.d-5,4.d-5,6.d-5,8.d-5,9D-5,
!     .     1.D-4,2.D-4,3.D-4,4.D-4,5.D-4,6.D-4,8.D-4,
!     .     1.D-3,2.D-3,3D-3,4.D-3,5.D-3,6.D-3,7D-3,8d-3,9d-3,
!     .     1.D-2,2.D-2,3D-2,4.D-2,5D-2,
!     .     6.D-2,6.5d-2,7D-2,7.5d-2,8.D-2,8.5d-2,9d-2,9.5d-2,
!     .     .1D0,.11d0,.125D0,.15D0,.175D0,.2D0,.225D0,.25D0,.275D0,
!     .     .3D0,.325D0,.35D0,.375D0,.4D0,.425D0,.45D0,.475D0,
!     .     .5D0,.525D0,.55D0,.575D0,.6D0,
!     .     .625d0,.65d0,.675d0,.7d0,.725d0,.75d0,.775d0,.8d0,
!     .     .825d0,.85d0,.875d0,.9d0,.92d0,.94d0,.96d0,.98d0,1d0/

        data xx1/1.d-5,2.d-5,4.d-5,6.d-5,8.d-5,9D-5,
     .     1.D-4,2.D-4,3.D-4,4.D-4,5.D-4,6.D-4,8.D-4,
     .     1.D-3,2.D-3,3D-3,4.D-3,5.D-3,6.D-3,7D-3,8d-3,9d-3,
     . 0.01d0,0.02d0,0.03d0,
     . 0.04d0,0.06d0,0.08d0,0.1d0,0.105d0,0.11d0,0.115d0,0.12d0,
     . 0.125d0,0.13d0,0.135d0,
     . 0.14d0,0.145d0,0.15d0,0.155d0,0.16d0,0.165d0,0.17d0,0.1725d0,
     . 0.175d0,0.1775d0,0.18d0,
     . 0.1825d0,0.185d0,0.1875d0,0.19d0,0.1925d0,0.195d0,0.1975d0,
     . 0.20d0,0.2025d0,0.205d0,
     . 0.2075d0,0.21d0,0.2125d0,0.215d0,0.2175d0,0.22d0,0.2225d0,
     . 0.225d0,0.23d0,0.235d0,
     . 0.24d0,0.245d0,0.25d0,0.255d0,0.26d0,0.265d0,0.27d0,0.275d0,
     . 0.28d0,0.285d0,0.29d0,0.295d0,
     . 0.30d0,0.305d0,0.31d0,0.315d0,0.32d0,0.325d0,0.33d0,0.335d0,
     . 0.34d0,0.345d0,0.35d0,0.355d0,
     . 0.36d0,0.365d0,0.37d0,0.375d0,0.38d0,0.385d0,0.39d0,
     . 0.395d0,0.40d0,0.405d0,0.41d0,0.415d0,0.42d0,0.425d0,0.43d0,
     . 0.435d0,0.44d0,.445d0,0.45d0,
     . 0.455d0,0.46d0,0.465d0,0.47d0,0.475d0,0.48d0,0.485d0,0.49d0,
     . 0.495d0,0.50d0,0.505d0,
     . 0.51d0,0.515d0,0.52d0,0.525d0,0.53d0,0.535d0,0.54d0,0.545d0,
     . 0.55d0,0.555d0,0.56d0,
     . 0.565d0,0.57d0,0.575d0,0.58d0,0.585d0,0.59d0,0.595d0,0.60d0,
     . 0.605d0,0.61d0,0.615d0,0.62d0,
     . 0.625d0,0.63d0,0.635d0,0.64d0,0.645d0,0.65d0,0.655d0,0.66d0,
     . 0.665d0,0.67d0,0.675d0,0.68d0,
     . 0.685d0,0.69d0,0.70d0,0.705d0,0.71d0,0.715d0,0.72d0,0.725d0,
     . 0.73d0,0.735d0,0.74d0,
     . 0.745d0,0.75d0,0.755d0,0.76d0,0.765d0,0.77d0,0.775d0,0.78d0,
     . 0.7825d0,0.785d0,0.7875d0,
     . 0.79d0,0.7925d0,0.795d0,0.7975d0,0.80d0,0.8025d0,0.805d0,
     . 0.8075d0,0.81d0,0.8125d0,
     . 0.815d0,0.8175d0,0.82d0,0.8225d0,0.825d0,0.8275d0,0.83d0,
     . 0.835d0,0.84d0,0.845d0,0.85d0,
     . 0.855d0,0.86d0,0.865d0,0.87d0,0.8725d0,0.875d0,
     . 0.8775d0,0.88d0,0.8825d0,0.885d0,0.8875d0,0.89d0,0.8925d0,
     . 0.895d0,0.8975d0,0.90d0,
     . 0.9d0,0.92d0,0.94d0,0.96d0,0.98d0,0.99d0,1d0/! 222


      a1=(param(1))
      a2=(param(2))
      a3=(param(3))
      a4=(param(4))
      a5=(param(5))
      a6=(param(6))
      a7=(param(7))
      a8=(param(8))
      chi2=0.d0
      
c.....aa can be changed to try to improve the fit
c.....the common changes it at the same time in the subroutine f0moments
      aa=2.5
c      aa=3.5
      
c.....take less/more points to improve speed/accuracy
      do i=1,npoints,1
         x=rx**( (1-xx1(i))/2.)   
         f=a1*x**a2*(1-x)**a3*(1+a4*x+a5*x**(0.5)+a6*x**(1.5)+A7*X**2
     .        +A8*X**(aa))   

         call distri(x,sqrt(rq2),unv,dnv,us,ds,str,chm,bot,glu,1)
         if(id.eq.1) then
            f2=unv
         elseif(id.eq.2) then
            f2=dnv
         elseif(id.eq.3) then
            f2=us
         elseif(id.eq.4) then
            f2=ds
         elseif(id.eq.5) then
            f2=str
         elseif(id.eq.6) then
            f2=glu
         elseif(id.eq.7) then
            f2=chm
         elseif(id.eq.8) then
            f2=bot
         endif
        aaa=1d0
         if (x.gt.rx) aaa=30d0
         chi2=chi2+(f2-f)**2*100*aaa*x**1.2
      enddo
      
      do i=1,npoints,1
         x=rx**( (1+xx1(i))/2.)     
         f=a1*x**a2*(1-x)**a3*(1+a4*x+a5*x**(0.5)+a6*x**(1.5)+A7*X**2
     .        +A8*X**(aa))         
         call distri(x,sqrt(rq2),unv,dnv,us,ds,str,chm,bot,glu,1)
         if(id.eq.1) then
            f2=unv
         elseif(id.eq.2) then
            f2=dnv
         elseif(id.eq.3) then
            f2=us
         elseif(id.eq.4) then
            f2=ds
         elseif(id.eq.5) then
            f2=str
         elseif(id.eq.6) then
            f2=glu
         elseif(id.eq.7) then
            f2=chm
         elseif(id.eq.8) then
            f2=bot
         endif
         aaaa=1d0
         if (x.gt.rx) aaaa=30d0
c        Modifica mia per x< Sqrt(tau)
         if(id.gt.2) then         
         chi2=chi2+(f2-f)**2*100*aaaa*x**1.2
         elseif(id.eq.1) then
         chi2=chi2+(f2-f)**2*100*aaaa
         elseif(id.eq.2) then
         chi2=chi2+(f2-f)**2*100*aaaa
         endif
        enddo
      
C.....return the parameters to save the last set
      do i=1,8
         parsal(i)=param(i)
      enddo
      return
      end  
      
       
     
     
      subroutine distri(x,q,upv,dnv,usea,dsea,str,chm,bot,glu,ippbar)
C.....here I call the new  sets !!!!! (in program prog_pdf.f)
C.....gives always x*distribution!!!!!
C.....MODIFICA MIA sett. 2004
      implicit none
      real*8 upv,dnv,usea,dsea,str,chm,bot,glu,x,q,q2
      integer isetproton,ippbar
      common/isetproton/isetproton
      REAL *8 FX(-5:5)
      q2=q*q
c      sx=sngl(x)
      call partons(q2,x,fx,5,isetproton,ippbar)
      usea=fx(-1)*x
      dsea=fx(-2)*x
      str=fx(-3)*x
      chm=fx(-4)*x
      bot=fx(-5)*x
      glu=fx(0)*x
      upv=fx(1)*x-usea 
      dnv=fx(2)*x-dsea
!      write(*,*)x,q,upv,dnv,usea,dsea,str,chm,bot,glu,ippbar
!      stop
      return
      end
      
      
C.....PSI - FUNCTION FOR COMPLEX ARGUMENT
      COMPLEX*16 FUNCTION PSIFN (Z)
      COMPLEX*16 Z, ZZ, RZ, DZ, SUB, CDLOG
      SUB = DCMPLX (0.,0.)
      ZZ = Z
 1    CONTINUE
      IF (DREAL (ZZ) .LT. 10.) THEN
         SUB = SUB - 1./ ZZ
         ZZ = ZZ + 1.
         GOTO 1
      END IF
      RZ = 1./ ZZ
      DZ = RZ * RZ
      PSIFN = SUB + CDLOG(ZZ) - RZ/2.- DZ/2520. * ( 210.+ DZ * (-21.+
     1     10.*DZ ))
      RETURN
      END
C     
C.....FIRST DERIVATIVE OF THE PSI - FUNCTION FOR COMPLEX ARGUMENT :
      COMPLEX*16 FUNCTION PSIFN1 (Z)
      COMPLEX*16 Z, ZZ, RZ, DZ, SUB
      SUB = DCMPLX (0.,0.)
      ZZ = Z
 1    CONTINUE
      IF (DREAL (ZZ) .LT. 10.) THEN
         SUB = SUB + 1./ (ZZ * ZZ)
         ZZ = ZZ + 1.
         GOTO 1
      END IF
      RZ = 1./ ZZ
      DZ = RZ * RZ
      PSIFN1 = SUB + RZ + DZ/2. * ( 1 + RZ/630. * ( 210.- DZ * ( 42.-
     1     DZ * ( 30.- 42.*DZ ))))
      RETURN
      END
C     
C.....SECOND DERIVATIVE OF THE PSI - FUNCTION FOR COMPLEX ARGUMENT :
      COMPLEX*16 FUNCTION PSIFN2 (Z)
      COMPLEX*16 Z, ZZ, RZ, DZ, SUB
      SUB = DCMPLX (0.,0.)
      ZZ = Z
 1    CONTINUE
      IF (DREAL (ZZ) .LT. 10.) THEN
         SUB = SUB - 2./ (ZZ * ZZ * ZZ)
         ZZ = ZZ + 1.
         GOTO 1
      END IF
      RZ = 1./ ZZ
      DZ = RZ * RZ
      PSIFN2 = SUB - DZ/60. * ( 60.+ RZ * ( 60.+ RZ * ( 30.- DZ *
     1     ( 10.- DZ * ( 10.- DZ * ( 18.- 50.* DZ ))))))
      RETURN
      END
C     
C.....BETA FUNCTION FOR COMPLEX ARGUMENT :
      COMPLEX*16 FUNCTION CBETA (Z1, Z2)
      implicit COMPLEX*16 (A - Z)
      LNGAM (X) = (X - 0.5) * CDLOG (X) - X + 0.91893853 + 1./(12.* X)
     1     * (1.- 1./(30.* X*X) * (1.- 1./(3.5 * X*X)
     2     * (1.- 4./(3.* X*X))))
      SUB = DCMPLX (0., 0.)
      ZZ1 = Z1
 1    CONTINUE
      IF ( DREAL (ZZ1) .LT. 15.) THEN
         SUB = SUB + CDLOG ((ZZ1+Z2) / ZZ1)
         ZZ1 = ZZ1 + 1.
         GOTO 1
      END IF
      ZZ2 = Z2
 2    CONTINUE
      IF ( DREAL (ZZ2) .LT. 15.) THEN
         SUB = SUB + CDLOG ((ZZ1+ZZ2) / ZZ2)
         ZZ2 = ZZ2 + 1.
         GOTO 2
      END IF
      LG1 = LNGAM (ZZ1)
      LG2 = LNGAM (ZZ2)
      LG12 = LNGAM (ZZ1 + ZZ2)
      CBETA = CDEXP (LG1 + LG2 - LG12 + SUB)
      RETURN
      END

C.....DOUBLE PRECISE GAUSS INTEGRATION :
      FUNCTION DINTEG (F, ALFA, BETA, EPS, K)
      implicit double precision (A-H, O-Z)
      dimension W(12), X(12)
      common / DIST / L
      data CONST / 1.0 D-12 /
      data W
     1 /0.101228536290376d0, 0.222381034453374d0, 0.313706645877887d0,
     2  0.362683783378362d0, 0.027152459411754d0, 0.062253523938647d0,
     3  0.095158511682492d0, 0.124628971255533d0, 0.149595988816576d0,
     4  0.169156519395002d0, 0.182603415044923d0, 0.189450610455068d0/
      data X
     1 /0.960289856497536d0, 0.796666477413627d0, 0.525532409916329d0,
     2  0.183434642495650d0, 0.989400934991649d0, 0.944575023073232d0,
     3  0.865631202387831d0, 0.755404408355003d0, 0.617876244402643d0,
     4  0.458016777657227d0, 0.281603550779258d0, 0.095012509837637d0/
      DINTEG = 0.0 D0
      M = 0
      L = K
      IF ( ALFA . EQ. BETA ) RETURN
      A = ALFA
      B = BETA
      DELTA = CONST * (DABS(A-B))
      AA = A
 1    Y = B - AA
      IF( DABS(Y) .LE. DELTA ) RETURN
C     IF( DABS(Y) .LE. DELTA ) THEN
C     WRITE (*,*) M
C     RETURN
C     END IF
 2    BB = AA + Y
      C1 = 0.5 D0 * (AA + BB)
      C2 = C1 - AA
      S8 = 0.0 D0
      S16 = 0.0 D0
      DO 15 I = 1, 4
         C3 = X(I) * C2
         S8 = S8 + W(I) * (F(C1+C3) + F(C1-C3))
 15   CONTINUE
      DO 16 I = 5, 12
         C3 = X(I) * C2
         S16 = S16 + W(I) * (F(C1+C3) + F(C1-C3))
  16   CONTINUE
       S8 = S8 * C2
       S16= S16 * C2
       IF( DABS(S16-S8) .GT. EPS * DABS(S8)) THEN
          Y = 0.5 * Y
          M = M + 1
          IF ( DABS(Y) .LE. DELTA ) THEN
             DINTEG = 0.0
             WRITE (*,10)
 10          FORMAT (1X,' DINTEG : TOO HIGH ACCURACY ')
          ELSE
             GOTO 2
          END IF
       ELSE
          DINTEG = DINTEG + S16
          AA = BB
          GOTO 1
       END IF
       RETURN
       END

      function alphasl(nq2)
      implicit none
      real *8 a_param,b0p,aass
      complex*16 xlambda,aa1,all,alphasl,nq2,qq,t,xlt,bstar,b,blog
      integer flagrealcomplex,imod,iord
      common/order/iord
      common/aass/aass
      common/flagrealcomplex/flagrealcomplex
      common/modified/imod
      common/a_param/a_param,b0p
      include 'scales.F'
      include 'const.F'

c.....Here computes NLL expression for alphas
c     here nq2=b0^2/b^2 and the result is now  alpha(nq2)/alpha(muf2)  


      b=b0/nq2**(0.5)

      blog=cdlog((q*b/(b0p))**2+1)

      xlambda=beta0*aass*blog


c     HERE now a dependence (without constant term)!


      aa1=cdlog(1-xlambda)+aass*iord*
     .      (beta1/beta0*cdlog(1-xlambda)/(1-xlambda) 
     .       + beta0*xlambda/(1-xlambda)*dlog(q2/mur2) !+
!     .          beta0*dlog(q2/muf2)                          ! TOGLIERE  (FATTO)
     .       -2*beta0*xlambda*dlog(a_param)/(1-xlambda)   )  
      alphasl=CDExp(-aa1)
      return
      end
      

      function alphaslB(nq2)
      implicit none
      real *8 a_param,b0p,aass
      complex*16 xlambda,aa2,all,alphaslB,nq2,qq,t,xlt,bstar,b,blog
      integer flagrealcomplex,imod,iord
      common/order/iord
      common/aass/aass
      common/flagrealcomplex/flagrealcomplex
      common/modified/imod
      common/a_param/a_param,b0p
      include 'scales.F'
      include 'const.F'

c.....This function computes the factor Exp(lambda/(1-lambda))
c.....where lambda=beta0*aass*blog

c.....This function is needed to resum the logarithms which multiply the 
c.....N-dependent part of the C coefficients (the N-independent part 
c.....is in the f2 function in auto.f).

      b=b0/nq2**(0.5)

      blog=cdlog((q*b/(b0p))**2+1)

      xlambda=beta0*aass*blog

      aa2=xlambda/(1-xlambda)

      alphaslB=CDExp(aa2)
      return
      end





      subroutine evol (x,q2,uv,dv,us,ds,ss,gl,qqb,qg,gg)
      implicit double precision (a-h, o-z) 
      complex*16 q2
      complex *16 uv,dv,us,ds,ss,gl,qqb,qg,gg,pa
      integer iord,flag4
      dimension pa(9)  
      common/order/iord 
      common/flag4/flag4
c.....call polarized parton distributions
      call inv(pa,x,q2)
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
c.....*** luminosities 
      qqb=pa(7) 
      qg=pa(8) 
      gg=pa(9) 
      end
      
*.....Parton density  calculation for the proton at externally fixed 
*     Bjorken-x, XB, and scale Q**2 in GeV**2, QS, calculated via a 
*     Mellin inversion of the analytic RGE solution in RENO.    
*
*.....The integration parameters (for the contour in the complex n-plane
*     and for the Gauss quadrature) are taken from the common-blocks
*     CONT, MOMS and WEIGHTS. 
                  
      subroutine inv(pa,xb,q2)    
      implicit none
      double precision c,xb,wn(136),xxx,q2s,lambda,x,ax,ex
      complex*16 cc(2),xnm(2),cdex(2),n(2,136),alps,
     .     fn(2,137,9),q2,II,fz,alphasl,alphaslB,alphas,aexp,aexpB
      complex *16 pa(9),fun
      integer nmax, i1, i2, m, iord, kp, flag4
      common/scales/q2s,lambda
      common/coupl/alps,aexp,aexpB
      common/contin/c,cc                              
      common/weights/wn
      common/moms/n
      common/order/iord
      common/flag4/flag4
      II=dcmplx(0.0,1.0)

*...  coupling constants at static point, input scale :
      alps=alphas(dcmplx(q2s,0d0))
*... alphasl gives the LL/NLL evolution of alpha from Q2 to b0^2/b^2
      aexp=alphasl(q2)
*... needed to resum xlam/(1-xlam)
      aexpB=alphaslB(q2)

*...  q2 and x dependent quantities : 
      x=xb                                    
      ax=dlog(x)                                        
      ex=dexp(-c*ax)                                   
*...  integration length parameter for mellin inversion
c      if ( X .le. 0.001 ) then                          
c         nmax = 24              ! zmax = 2 
c      else if ( X .le. 0.05 ) then                   
c         nmax = 40              ! zmax = 4  
c      else if ( X .le. 0.2 ) then                   
c         nmax = 56              ! zmax = 8  
c      else if ( X .le. 0.4 ) then                   
c         nmax = 72              ! zmax = 12 
c      else if ( X .le. 0.7 ) then                   
c         nmax = 88              ! zmax = 18 
c      else                                                          
c         nmax = 136             ! zmax = 36
c      end if              
      nmax=136
      xxx=1d0                               
*...  calculation of the parton densities and output:
*     index of func :   1     2    3     4    5   6   
*     distribution  :  UV    DV    UB    DB   S   G  
      call reno (fn, n, nmax,iord)
CC      do 2 i2 = 1, 9
         i2=9
         fun=(0d0,0d0)
         if(i2.ge.7) xxx=xb
         do 1 i1 = 1, nmax
         do 3 kp=1,2
               xnm(kp) = (c - n(kp,i1)+1.) * ax
               cdex(kp) = cdexp (xnm(kp)) / 3.1415927 * cc(kp)/2./II
 3          continue
c     fz = dble (aimag (fn(kp,i1+1,i2)*cdex(kp)))
           fz=dcmplx ((( fn(1,i1+1,i2)*cdex(1) -fn(2,i1+1,i2)*cdex(2))))
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
      
      subroutine reno (fn, n, nmax, iord)
      implicit none
      integer channel
      real *8 pi,beta0,beta1,aass,GE,zeta4,zeta2,zeta3,Cf,CA,log2
      real *8 a_param,b0p,loga,A1g,B1g,H1g,A2g,B2g
      complex*16 ans(2,136), am(2,136), ap(2,136), al(2,136),
     1     be(2,136), ab(2,136), ac(2,136), rmin(2,136),
     2     rplus(2,136), C1Q(2,136), C1G(2,136),C3Q(2,136), 
     3     rmmqq(2,136),rmmqg(2,136),rmmgq(2,136) ,rmmgg(2,136),
     4     rmpqq(2,136),rmpqg(2,136),rmpgq(2,136),rmpgg(2,136),
     5     rpmqq(2,136),rpmqg(2,136),rpmgq(2,136),rpmgg(2,136),
     6     rppqq(2,136),rppqg(2,136),rppgq(2,136),rppgg(2,136),
     7     uvi(2,136), dvi(2,136), usi(2,136), dsi(2,136),
     8     ssi(2,136), gli(2,136), chi(2,136), boi(2,136),
     9     fn(2,137,9), n(2,136)
C
C ADDED: Anomalous dimensions used for scale dependence
C
      complex*16 QQIM(2,136),QGFM(2,136),GQIM(2,136),GGIM(2,136),
     #      GGFM(2,136),NS1MIM(2,136),NS1PIM(2,136),NS1FM(2,136),
     #      QQ1FM(2,136),QG1FM(2,136),GQ1IM(2,136),GQ1FM(2,136),
     #      GG1IM(2,136),GG1FM(2,136)

      common/anomdim/QQIM,QGFM,GQIM,GGIM,GGFM,NS1MIM,NS1PIM,NS1FM,
     #      QQ1FM,QG1FM,GQ1IM,GQ1FM,GG1IM,GG1FM

      common/channel/channel
      complex*16 H2ggM(2,136),H2gqM(2,136),H2qqM(2,136)
      common / H2mat /H2ggM,H2gqM,H2qqM
C


C
      complex *16 S,XL,XL1,ENS,EM,EP,EMP,EPM,SG,GL,SSN,DSN,USN,CHN
      complex *16 BON,alp,alpn,alps,aexp,aexpB,oswi
      complex *16 QQ1FN,QG1FN,GQ1IN,GQ1FN,GG1IN,GG1FN
      complex *16 gamma1qq,gamma1gq,gamma1qg,gamma1gg
      complex *16 gamma2qqV,gamma2qqS,gamma2qqbV,gamma2qqbS,
     / gamma2gq,gamma2qg,gamma2gg
      complex *16 UVN,DVN,NS3N,NS8N,GLN,SIN,NS15N,NS24N,NS35N
      complex *16 C1gq,C1gg,H1stgg,H1stgq,C2gq
      complex *16 Hgg,Hgq,Hqq
      complex *16 H2stgg,H2stgq,H2stqq
      complex *16 QGN,GGN,QQBARN,qqn1,gqn1,qgn1,ggn1,G1N
C
      INTEGER  F,NF,K1,KK1,NMAX,iord,iiord,KP,flag1,flag4
   
      common / ANOMS  / ANS, AM, AP, AL, BE, AB, AC, RMIN, RPLUS, C1Q, 
     1     C1G, C3Q, RMMQQ, RMMQG, RMMGQ, RMMGG, RMPQQ, RMPQG, RMPGQ, 
     2     RMPGG, RPMQQ, RPMQG, RPMGQ, RPMGG, RPPQQ, RPPQG, RPPGQ, 
     3     RPPGG
      common / INPUT  / UVI, DVI, USI, DSI, SSI, CHI, BOI, GLI
      common/ order/ iiord
      common/nf/nf
      common/coupl/alps,aexp,aexpB
      common/flag1/flag1
      common/flag4/flag4
      common/aass/aass
      common/a_param/a_param,b0p
!

      integer si,sj,sk,ippbar,flagch,pcf
      common/ippbar/ippbar
      complex *16 FX1(-5:5), FX2(-5:5)
      common/flagch/flagch
!

      include 'scales.F'

      data pi/3.1415926536d0/
      data zeta3/1.20205690315959429d0/

      iord=iiord
      oswi=dcmplx(Real(iord),0.)

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
C
      XL = 1d0/aexp 
      ALP = alps*aexp
C
      S   = CDLOG (XL)
      XL1 = 1.- XL
      ENS = CDEXP (-ANS(KP,K1)*S)
      EM  = CDEXP (-AM(KP,K1)*S)
      EP  = CDEXP (-AP(KP,K1)*S)
      EMP = EM / EP
      EPM = EP / EM
      
*...  EVOLUTION OF LIGHT PARTON DESITIES BETWEEN THRESHOLDS :
      UVN  = UVN  * ENS * (1.+ oswi * ALP * XL1 * RMIN(KP,K1))
      DVN  = DVN  * ENS * (1.+ oswi * ALP * XL1 * RMIN(KP,K1))
      NS3N = NS3N * ENS * (1.+ oswi * ALP * XL1 * RPLUS(KP,K1))
      NS8N = NS8N * ENS * (1.+ oswi * ALP * XL1 * RPLUS(KP,K1))
      SG = SIN
      GL = GLN
      SIN = EM * ((AL(KP,K1) + oswi * ALP * (RMMQQ(KP,K1) * XL1 
     1     + RMPQQ(KP,K1) * (EPM-XL)))* SG  
     2     + (BE(KP,K1) + oswi * ALP * (RMMQG(KP,K1) * XL1 
     3     + RMPQG(KP,K1) * (EPM-XL))) * GL)
     4     + EP * ((AC(KP,K1) + oswi * ALP * (RPPQQ(KP,K1) * XL1 
     5     + RPMQQ(KP,K1) * (EMP-XL))) * SG
     6     +(-BE(KP,K1) + oswi * ALP * (RPPQG(KP,K1) * XL1 
     7     + RPMQG(KP,K1) * (EMP-XL))) * GL)
      GLN = EM * ((AB(KP,K1) + oswi * ALP * (RMMGQ(KP,K1) * XL1 
     1     + RMPGQ(KP,K1) * (EPM-XL))) * SG 
     2     + (AC(KP,K1) + oswi * ALP * (RMMGG(KP,K1) * XL1 
     3     + RMPGG(KP,K1) * (EPM-XL))) * GL)
     4     + EP *((-AB(KP,K1) + oswi * ALP * (RPPGQ(KP,K1) * XL1 
     5     + RPMGQ(KP,K1) * (EMP-XL))) * SG
     6     + (AL(KP,K1) + oswi * ALP * (RPPGG(KP,K1) * XL1 
     7     + RPMGG(KP,K1) * (EMP-XL))) * GL)
      
      NS15N = NS15N * ENS * (1.+  oswi * ALP * XL1 * RPLUS(KP,K1))
      NS24N = NS24N * ENS * (1.+  oswi * ALP * XL1 * RPLUS(KP,K1))
      NS35N = SIN
      
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


c     INCLUDE HERE THE CORRESPONDING COEFFICIENTS C1 AND H1

c     This is alpha_S(b0^2/b^2)

       alpn=alp*4*pi

       C1gq=4/3d0/(N(KP,K1)+1)
       C1gg=pi**2/2d0+11/2d0+pi**2
       H1g=0d0

       pcf=2

       if(flag1.eq.1) then   


c Full hard coefficient: expand C*H*C at order alpha_s.
c Add renorm. scale dep terms due to the fact that
c the LO x-section begins at order alphas^2

      beta0=(33-2*nf)/12d0
      beta1=(153-19*nf)/24d0 


      Hgg=1d0+aass/2d0*(H1g+2*C1gg-2*pcf*beta0*dlog(q2/mur2))
!      Hgq=alpn/2d0/pi*C1gq
      Hgq=aass/2d0*aexp*C1gq
      Hqq=0d0


      if(muf2.ne.q2.or.a_param.ne.1) then 

      A1g=3d0
      B1g=-(33d0-2*nf)/6d0


c      call ancalc(QQIN, QGFN, GQIN, GGIN, GGFN, NS1MIN, NS1PIN, NS1FN,
c     /     QQ1FN, QG1FN, GQ1IN, GQ1FN, GG1IN, GG1FN, N(KP,K1))

      gamma1qq=-1d0*(QQIM(KP,K1)/4d0)
      gamma1qg=-1d0*(QGFM(KP,K1)/8d0)
      gamma1gq=-1d0*(GQIM(KP,K1)/4d0)
      gamma1gg=-1d0*((GGIM(KP,K1)+nf*GGFM(KP,K1))/4d0)
            
      
      Hgg=Hgg-aass*gamma1gg*(dlog(muf2/q2)+2d0*log(a_param))
      Hgq=Hgq-aass/2d0*gamma1gq*(dlog(muf2/q2)+2d0*log(a_param))

      Hgg=Hgg+aass/2d0*(-4*log(a_param))*(B1g+A1g*log(a_param))

      endif

CC    INCLUSION OF H2

       elseif(flag1.eq.2) then

       H1stgg=(H1g+2d0*C1gg)
       H1stgq=C1gq

       A1g=3d0
       B1g=-(33d0-2*nf)/6d0
       
       CA=3d0
       Cf=4d0/3
       A2g=CA/2d0*(67/6d0-(pi**2)/2d0-5d0/9*nf)
       B2g=CA**2*(23/24d0+(11*pi**2)/18d0-3*zeta3/2d0)+
     /        Cf*nf/2d0-CA*nf*(1/12d0+pi**2/9d0)-11/8d0*Cf*CA
       beta0=(33-2*nf)/12d0
       beta1=(153-19*nf)/24d0


C     Mellin transform of (1-z)/z

      G1N=1/N(KP,K1)/(N(KP,K1)-1)

C     H2 coefficients without effect of spin correlations

      H2stgg=4.0D0*(H2ggM(KP,K1)-CA**2*G1N*G1N)
      H2stgq=4.0D0*(H2gqM(KP,K1)-CA*CF*G1N*G1N)
      H2stqq=4.0D0*(H2qqM(KP,K1)-CF**2*G1N*G1N)


C     C2gq coefficient (alpha/(2pi) normalization)

      C2gq=H2stgq-C1gg*C1gq


CC H2 scale dependence: here is from DY

      loga=log(a_param)

c      call ancalc(QQIN, QGFN, GQIN, GGIN, GGFN, NS1MIN, NS1PIN, NS1FN,
c     /     QQ1FN, QG1FN, GQ1IN, GQ1FN, GG1IN, GG1FN, N(KP,K1))

      gamma1qq=-1d0*(QQIM(KP,K1)/4d0)
      gamma1qg=-1d0*(QGFM(KP,K1)/8d0)
      gamma1gq=-1d0*(GQIM(KP,K1)/4d0)
      gamma1gg=-1d0*((GGIM(KP,K1)+nf*GGFM(KP,K1))/4d0)
!
      gamma2qg=-1d0*(QG1FM(KP,K1)/16d0)
      gamma2gq=-1d0*((GQ1IM(KP,K1)+nf*GQ1FM(KP,K1))/8d0)
      gamma2gg=-1d0*((GG1IM(KP,K1)+nf*GG1FM(KP,K1))/8d0)  !OK
!
      gamma2qqV=-1d0*((NS1PIM(KP,K1)+2*nf*NS1FM(KP,K1)
     #          +NS1MIM(KP,K1)))/16d0
      gamma2qqbV=-1d0*(NS1PIM(KP,K1)-NS1MIM(KP,K1))/16d0
      gamma2qqS=-1d0*QQ1FM(KP,K1)/16d0
      gamma2qqbS=gamma2qqS


C    Scale dependence: first order

      H1stgg=H1stgg+(-4*log(a_param))*(B1g+A1g*log(a_param))
     / -2d0*gamma1gg*(log(muf2/q2)+2d0*log(a_param))
      H1stgq=H1stgq-gamma1gq*(log(muf2/q2)+2d0*log(a_param))
      H1stgg=H1stgg-2*pcf*beta0*dlog(q2/mur2)

      

c    Second order: muF,muR and Q
c    pcF=2 

C da qui:check

      H2stgg=H2stgg+4d0*( (A1g/6d0*beta0*8d0*log(a_param)**3
     / +2d0*log(a_param)**2*(A2g+beta0*(-(B1g+2d0*A1g*log(a_param)+
     / gamma1gg)))
     / -2d0*log(a_param)*(B2g+2d0*A2g*log(a_param)-beta0*(C1gg)
     / +2/4d0*gamma2gg) 
     / +beta0/2d0*(gamma1gg)*log(q2/muf2)**2
     / +(2/4d0*gamma2gg)*log(q2/muf2)
     /+1/4d0*(H1stgg+H1g+2d0*C1gg)*(gamma1gg*(log(q2/muf2)
     / -2d0*log(a_param))-
     / ((B1g+A1g*log(a_param))*2d0*log(a_param)+pcf*beta0*log(q2/mur2)))
     /+1/4d0*(H1stgq+C1gq)*(1/2d0*gamma1qg*(log(q2/muf2)
     / -2d0*log(a_param)))
     /+1/4d0*(H1stgq+C1gq)*(1/2d0*gamma1qg*(log(q2/muf2)
     / -2d0*log(a_param))))
     / -pcf*(0.5d0*(beta0*log(q2/mur2))**2+beta1*log(q2/mur2))
     /) 

      

      H2stgq=H2stgq+4d0*(
     / +2d0*log(a_param)**2*(beta0*(-1/2d0*gamma1gq)) 
     / -2d0*log(a_param)*(-beta0*(1/2d0*C1gq)
     / +1/4d0*(gamma2gq))  
     / +beta0/2d0*(1/2d0*gamma1gq)*log(q2/muf2)**2
     / +(1/4d0*gamma2gq)*log(q2/muf2)
     /+1/4d0*(H1stgg+H1g+2d0*C1gg)
     /*(1/2d0*gamma1gq*(log(q2/muf2)
     / -2d0*log(a_param)))
     /+1/4d0*(H1stgq+C1gq)*(1/2d0*(gamma1gg+gamma1qq)
     /*(log(q2/muf2)
     / -2d0*log(a_param))
     /-((B1g+A1g*log(a_param))*2d0*log(a_param)+pcf*beta0*log(q2/mur2))) 
     /)


      H2stqq=H2stqq+4d0*(
     /1/2d0*(H1stgq+C1gq)*(gamma1gq/2d0*(log(q2/muf2)
     /-2d0*log(a_param)))
     /)

  


c    Second order: additional muR dependent terms (fourth line of eq. (70))



       H2stgg=H2stgg-4d0*(beta0*log(q2/mur2)*H1stgg/2d0)
       H2stgq=H2stgq-4d0*(beta0*log(q2/mur2)*H1stgq/2d0)



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



CC H coefficient up to NNLO  with exponentiated logarithms up to NNLL
    

       Hgg=(1+(aass/2d0)*H1stgg+(aass/2d0)**2*H2stgg)  ! no exponential: here C1gg-C1ggdelta=0 !
       Hgg=Hgg+aass**2*(CA*G1N)**2*aexp*aexp           ! spin correlations

       Hgq=(aass/2d0)*H1stgq+(aass/2d0)**2*H2stgq
     /     *aexp*aexpB**(aass/2d0*(C2gq/C1gq-C1gg))    ! g<-q factor
       Hgq=Hgq+aass**2*CA*CF*(G1N)**2*aexp*aexp        ! spin correlations

       Hqq=(aass/2d0)**2*H2stqq
     /     *aexp*aexpB**(aass/2d0*(C2gq/C1gq-C1gg))    ! g<-q factor
     /     *aexp*aexpB**(aass/2d0*(C2gq/C1gq-C1gg))    ! g<-q factor
       Hqq=Hqq+aass**2*(CF*G1N)**2*aexp*aexp           ! spin correlations

    
      endif


 
         QQN1=(UVN+2*USN+DVN+2*DSN+2*SSN+2*CHN+2*BON)**2
         QGN1=2*(UVN+2*USN+DVN+2*DSN+2*SSN+2*CHN+2*BON)*GLN
         GGN1=GLN*GLN

!         GGN=Hgg*GGN1 + QGN1*Hgq + QQN1*Hqq
!  Channels
      if(channel.eq.1)then
       GGN=Hgg*GGN1 + QGN1*Hgq + QQN1*Hqq           ! Tot
      elseif(channel.eq.2)then                      ! GG channel
       GGN=Hgg*GGN1                                 ! Only channel gluon-gluon
      elseif(channel.eq.3)then                      ! GQ channel
       GGN=QGN1*Hgq                                 ! Only channel gluon-quark
      elseif(channel.eq.4)then                      ! QQ channel
       GGN=QQN1*Hqq                                 ! Only channel quark-quark
      endif


*...  OUTPUT
      FN(KP,KK1+1,1) = UVN
      FN(KP,KK1+1,2) = DVN
      FN(KP,KK1+1,3) = USN
      FN(KP,KK1+1,4) = DSN
      FN(KP,KK1+1,5) = SSN
      FN(KP,KK1+1,6) = GLN
      FN(KP,KK1+1,7) = QQBARN
      FN(KP,KK1+1,8) = QGN
      FN(KP,KK1+1,9) = GGN
      
 1    CONTINUE
      RETURN
      END
      
*     
*...Moments of NLO-input distributions
*     
      
      SUBROUTINE INP200
      implicit COMPLEX *16(A - Z)            
      dimension UV(2,136), DV(2,136),US(2,136),DS(2,136),GL(2,136),
     1     SS(2,136), CH(2,136), BO(2,136), N(2,136)
      INTEGER K1,nnf, KP
      common/nf/nnF
      common / MOMS / N
      common / INPUT  / UV, DV, US, DS, SS, CH, BO, GL
      save / INPUT /
      
      DO 1 KP=1,2
         DO 1 K1 = 1, 136
            XN = N(KP,K1)
C     Notacion
C     la transformada de  x^a (1-x)^b  (*x^(n-1)) es
C     Gamma[1+b] Gamma[a+n]/ Gamma[1+a+b+n] =  BETA(a+n ,1+b)
C     aqui va la transformada de f(x) (y no x*f(x))
c     el programa devuelve x*f(x) 
            call F0MOMENTS(XN,UVi,DVi,USi,DSi,SSi,GLi,CHi,BOi) 
            UV(KP,K1)=UVI
            DV(KP,K1)=DVI
            US(KP,K1)=USI
            DS(KP,K1)=DSI
            SS(KP,K1)=SSI
            GL(KP,K1)=GLI
            CH(KP,K1)=CHI
            BO(KP,K1)=BOI
            if (nnf.eq.3) then  !GRV
               CH(KP,K1)=CMPLX(0D0,0D0)
               BO(KP,K1)=CMPLX(0D0,0D0)
            endif
            
 1       CONTINUE     
         RETURN
         END
      
*     
*...  Anomalous dimensions for leading and next to leading order
*     evolution of parton densities and Wilson coefficients for NLO 
*     structure functions.
*     The moments are calculated on an externally given array of mellin 
*     moments, N, suitable for a (non-adaptive) quadrature.
*     
*     

      SUBROUTINE ANOM (ANS, AM, AP, AL, BE, AB, AC, RMIN, RPLUS, 
     1     RMMQQ, RMMQG, RMMGQ, RMMGG, RMPQQ, RMPQG, RMPGQ, RMPGG, 
     2     RPMQQ, RPMQG, RPMGQ, RPMGG, RPPQQ, RPPQG, RPPGQ, RPPGG, 
     1     C1Q, C1G, C3Q,
C
C    ADDED
C
     #      QQIM,QGFM,GQIM,GGIM,GGFM,NS1MIM,NS1PIM,NS1FM,
     #      QQ1FM,QG1FM,GQ1IM,GQ1FM,GG1IM,GG1FM,N)
      implicit COMPLEX*16 (A - Z)
      dimension ANS(2,136), AM(2,136), AP(2,136), AL(2,136),
     1     BE(2,136), AB(2,136), AC(2,136), RMIN(2,136), 
     2     RPLUS(2,136), N(2,136), C1Q(2,136), C1G(2,136), C3Q(2,136), 
     3     RMMQQ(2,136),RMMQG(2,136),RMMGQ(2,136),RMMGG(2,136),
     4     RMPQQ(2,136),RMPQG(2,136),RMPGQ(2,136),RMPGG(2,136),
     5     RPMQQ(2,136),RPMQG(2,136),RPMGQ(2,136),RPMGG(2,136),
     6     RPPQQ(2,136),RPPQG(2,136),RPPGQ(2,136),RPPGG(2,136)

C    ADDED

           dimension QQIM(2,136),QGFM(2,136),GQIM(2,136),GGIM(2,136),
     #      GGFM(2,136),NS1MIM(2,136),NS1PIM(2,136),NS1FM(2,136),
     #      QQ1FM(2,136),QG1FM(2,136),GQ1IM(2,136),GQ1FM(2,136),
     #      GG1IM(2,136),GG1FM(2,136)

C

      REAL*8 B0, B1, B0F
      INTEGER F, FR, K1, K2, KP
      common/nf/F

      DO 1 KP=1,2
         DO 1 K1 = 1, 136
            XN = N(KP,K1)
            CALL ANCALC (QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F,
     1           QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, XN)
C
C
            QQIM(KP,K1)=QQI
            QGFM(KP,K1)=QGF
            GQIM(KP,K1)=GQI
            GGIM(KP,K1)=GGI
            GGFM(KP,K1)=GGF
            NS1MIM(KP,K1)=NS1MI
            NS1PIM(KP,K1)=NS1PI
            NS1FM(KP,K1)=NS1F
            QQ1FM(KP,K1)=QQ1F
            QG1FM(KP,K1)=QG1F
            GQ1IM(KP,K1)=GQ1I
            GQ1FM(KP,K1)=GQ1F
            GG1IM(KP,K1)=GG1I
            GG1FM(KP,K1)=GG1F
C
            FR = 5
*...  ANOMALOUS DIMENSIONS AND RELATED QUANTITIES IN LEADING ORDER :
            B0 = 11.- 2./3.* FR
            B0F = 11.- 2./3.* F
            B02 = 2.* B0
            B02F = 2.* B0F
            QQ = QQI
            QG = F * QGF
            GQ = GQI
            GG = GGI + F * GGF
            SQ = CDSQRT ((GG - QQ) * (GG - QQ) + 4.* QG * GQ)
            GP = 0.5 * (QQ + GG + SQ)
            GM = 0.5 * (QQ + GG - SQ)
            XANS = QQ / B02
            XAM = GM / B02
            XAP = GP / B02
            XAL = (QQ - GP) / (GM - GP)
            XBE = QG / (GM - GP)
            XAB = GQ / (GM - GP)
            XAC = 1.- XAL
            ANS(KP,K1) = XANS
            AM(KP,K1) = XAM
            AP(KP,K1) = XAP
            AL(KP,K1) = XAL
            BE(KP,K1) = XBE
            AB(KP,K1) = XAB
            AC(KP,K1) = XAC
*...  NEXT TO LEADING ORDER : ANOMALOUS DIMENSIONS IN THE MS-BAR 
*     FACTORIZATION SCHEME OF BARDEEN ET AL. (1981) :
            NS1M = NS1MI + F * NS1F
            NS1P = NS1PI + F * NS1F
            QQ1 = NS1P + F * QQ1F
            QG1 = F * QG1F
            GQ1 = GQ1I + F * GQ1F
            GG1 = GG1I + F * GG1F
*...  COMBINATIONS OF ANOMALOUS DIMENSIONS FOR NLO - SINGLET EVOLUTION :
            B1 = 102 - 38./3.* FR
            B10 = B1 / B0
            RMIN(KP,K1) = (NS1M - QQ * B10) / B02
            RPLUS(KP,K1) = (NS1P - QQ * B10) / B02
            RQQ = (QQ1 - QQ * B10) / B02
            RQG = (QG1 - QG * B10) / B02
            RGQ = (GQ1 - GQ * B10) / B02
            RGG = (GG1 - GG * B10) / B02
            NMP = 1.- XAM + XAP
            NPM = 1.- XAP + XAM
            DMQQ =  XAL * RQQ + XBE * RGQ
            DMQG =  XAL * RQG + XBE * RGG
            DMGQ =  XAB * RQQ + XAC * RGQ
            DMGG =  XAB * RQG + XAC * RGG
            DPQQ =  XAC * RQQ - XBE * RGQ
            DPQG =  XAC * RQG - XBE * RGG
            DPGQ = -XAB * RQQ + XAL * RGQ
            DPGG = -XAB * RQG + XAL * RGG
            RMMQQ(KP,K1) =   XAL * DMQQ + XAB * DMQG
            RMMQG(KP,K1) =   XBE * DMQQ + XAC * DMQG
            RMMGQ(KP,K1) =   XAL * DMGQ + XAB * DMGG
            RMMGG(KP,K1) =   XBE * DMGQ + XAC * DMGG
            RMPQQ(KP,K1) =  (XAC * DMQQ - XAB * DMQG) / NMP
            RMPQG(KP,K1) = (-XBE * DMQQ + XAL * DMQG) / NMP
            RMPGQ(KP,K1) =  (XAC * DMGQ - XAB * DMGG) / NMP
            RMPGG(KP,K1) = (-XBE * DMGQ + XAL * DMGG) / NMP
            RPMQQ(KP,K1) =  (XAL * DPQQ + XAB * DPQG) / NPM
            RPMQG(KP,K1) =  (XBE * DPQQ + XAC * DPQG) / NPM
            RPMGQ(KP,K1) =  (XAL * DPGQ + XAB * DPGG) / NPM
            RPMGG(KP,K1) =  (XBE * DPGQ + XAC * DPGG) / NPM
            RPPQQ(KP,K1) =   XAC * DPQQ - XAB * DPQG
            RPPQG(KP,K1) =  -XBE * DPQQ + XAL * DPQG
            RPPGQ(KP,K1) =   XAC * DPGQ - XAB * DPGG
            RPPGG(KP,K1) =  -XBE * DPGQ + XAL * DPGG
 2          CONTINUE
 1          CONTINUE
            RETURN
            END


*     
*...  CALCULATION OF ANOMALOUS DIMENSIONS AND WILSON COEFFICIENTS
*     UP TO THEIR DEPENDENCE OF THE NUMBER OF ACTIVE FLAVOURS F :
*     
      SUBROUTINE ANCALC (QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F,
     1     QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F,  XN)
      implicit COMPLEX*16 (A - Z)
      REAL*8 ZETA2, ZETA3
      XNS = XN * XN
      XN1 = XN + 1.
      XN2 = XN + 2.
      XNM = XN - 1.
C...  LEADING ORDER :
      CPSI = PSIFN (XN1) + 0.577216
      QQI = (8./3.) * (-3.- 2./(XN * XN1) + 4.* CPSI)
      QGF = -4.* (XNS + XN +2.) / (XN * XN1 * XN2)
      GQI = -(16./3.) * (XNS + XN + 2.) / (XN * XN1 * XNM)
      GGI = -22.- 24./(XN * XNM) - 24./(XN1 * XN2) + 24.* CPSI
      GGF = 4./3.
C...  NEXT TO LEADING ORDER :
      XNT = XNS * XN
      XNFO = XNT * XN
      XN1S = XN1 * XN1
      XN1T = XN1S * XN1
C...  ANALYTIC CONTINUATIONS OF N-SUMS AS GIVEN IN GLUECK ET AL. (1990) :
      ZETA2 = 1.644934
      ZETA3 = 1.202057
      CPSI1 = ZETA2 - PSIFN1 (XN1)
      CPSI2 = 0.5 * PSIFN2 (XN1) - ZETA2
      SPMOM = 1.01/XN1 - 0.846/XN2 + 1.155/(XN+3.) - 1.074/(XN+4.) +
     1     0.55/(XN+5.)
      SLC = -5./8.* ZETA3
      SLV = - ZETA2/2.* (PSIFN (XN1/2.) - PSIFN (XN/2.))
     1     + CPSI/XNS + SPMOM
      SSCHLM = SLC - SLV
      SSTR2M = ZETA2 - PSIFN1 (XN1/2.)
      SSTR3M = 0.5 * PSIFN2 (XN1/2.) + ZETA3
      SSCHLP = SLC + SLV
      SSTR2P = ZETA2 - PSIFN1 (XN2/2.)
      SSTR3P = 0.5 * PSIFN2 (XN2/2.) + ZETA3
C...  NON-SINGLET PIECES AS GIVEN IN CURCI ET AL. (1980) :
      NS1MA = 16.* CPSI * (2.* XN + 1.) / (XNS * XN1S) +
     1     16.* (2.* CPSI - 1./(XN * XN1)) * ( CPSI1 - SSTR2M ) +
     2     64.* SSCHLM + 24.* CPSI1 - 3. - 8.* SSTR3M -
     3     8.* (3.* XNT + XNS -1.) / (XNT * XN1T) +
     4     16.* (2.* XNS + 2.* XN +1.) / (XNT * XN1T)
      NS1PA = 16.* CPSI * (2.* XN + 1.) / (XNS * XN1S) +
     1     16.* (2.* CPSI - 1./(XN * XN1)) * ( CPSI1 - SSTR2P ) +
     2     64.* SSCHLP + 24.* CPSI1 - 3. - 8.* SSTR3P -
     3     8.* (3.* XNT + XNS -1.) / (XNT * XN1T) -
     4     16.* (2.* XNS + 2.* XN +1.) / (XNT * XN1T)
      NS1B = CPSI * (536./9. + 8.* (2.* XN + 1.) / (XNS * XN1S)) -
     1     (16.* CPSI + 52./3.- 8./(XN * XN1)) * CPSI1 - 43./6. -
     2     (151.* XNFO + 263.* XNT + 97.* XNS + 3.* XN + 9.) *
     3     4./ (9.* XNT * XN1T)
      NS1C = -160./9.* CPSI + 32./3.* CPSI1 + 4./3. +
     1     16.* (11.* XNS + 5.* XN - 3.) / (9.* XNS * XN1S)
      NS1MI = -2./9.* NS1MA + 4.* NS1B
      NS1PI = -2./9.* NS1PA + 4.* NS1B
      NS1F = 2./3. * NS1C
C...  SINGLET PIECES AS GIVEN IN FLORATOS ET AL. (1981) :
      XNFI = XNFO * XN
      XNSI = XNFI * XN
      XNSE = XNSI * XN
      XNE = XNSE * XN
      XNN = XNE * XN
      XNMS = XNM * XNM
      XN2S = XN2 * XN2
      XN2T = XN2S * XN2
      QQ1F = (5.* XNFI + 32.* XNFO + 49.* XNT + 38.* XNS + 28.* XN
     1     + 8.) / (XNM * XNT * XN1T * XN2S) * (-32./3.)
      QG1A = (-2.* CPSI * CPSI + 2.* CPSI1 - 2.* SSTR2P)
     1     * (XNS + XN + 2.) / (XN * XN1 * XN2)
     2     + (8.* CPSI * (2.* XN + 3.)) / (XN1S * XN2S)
     3     + 2.* (XNN + 6.* XNE + 15. * XNSE + 25.* XNSI + 36.* XNFI
     4     + 85.* XNFO + 128.* XNT + 104.* XNS + 64.* XN + 16.)
     5     / (XNM * XNT * XN1T * XN2T)
      QG1B = (2.* CPSI * CPSI - 2.* CPSI1 + 5.) * (XNS + XN + 2.)
     1     / (XN * XN1 * XN2)   -   4.* CPSI / XNS
     2     + (11.* XNFO + 26.* XNT + 15.* XNS + 8.* XN + 4.)
     3     / (XNT * XN1T * XN2)
      QG1F = - 12.* QG1A - 16./3.* QG1B
      GQ1A = (-2.* CPSI * CPSI + 10.* CPSI - 2.* CPSI1)
     1     * (XNS + XN + 2.) / (XNM * XN * XN1)  -  4.* CPSI / XN1S
     2     - (12.* XNSI + 30.* XNFI + 43.* XNFO + 28.* XNT - XNS
     3     - 12.* XN - 4.) / (XNM * XNT * XN1T)
      GQ1B = (CPSI * CPSI + CPSI1 - SSTR2P) * (XNS + XN + 2.)
     1     / (XNM * XN * XN1)
     2     - CPSI * (17.* XNFO + 41.* XNS - 22.* XN - 12.)
     3     / (3.* XNMS * XNS * XN1)
     4     + (109.* XNN + 621.* XNE + 1400.* XNSE + 1678.* XNSI
     5     + 695.* XNFI - 1031.* XNFO - 1304.* XNT - 152.* XNS
     6     + 432.* XN + 144.) / (9.* XNMS * XNT * XN1T * XN2S)
      GQ1C = (CPSI - 8./3.) * (XNS + XN + 2.) / (XNM * XN * XN1)
     1     + 1./ XN1S
      GQ1I = - 64./9.* GQ1A - 32.* GQ1B
      GQ1F = - 64./9.* GQ1C
      GG1A = 16./9.* (38.* XNFO + 76.* XNT + 94.* XNS + 56.* XN + 12.)
     1     / (XNM * XNS * XN1S * XN2)   -   160./9.* CPSI + 32./3.
      GG1B = (2.* XNSI + 4.* XNFI + XNFO - 10.* XNT - 5.* XNS - 4.* XN
     1     - 4.) * 16. / (XNM * XNT * XN1T * XN2)   +   8.
      GG1C = (2.* XNFI + 5.* XNFO + 8.* XNT + 7.* XNS - 2.* XN - 2.)
     1     * 64.* CPSI / (XNMS * XNS * XN1S * XN2S)
     2     + 536./9.* CPSI - 64./3.
     3     + 32.* SSTR2P * (XNS + XN + 1.) / (XNM * XN * XN1 * XN2)
     4     - 16.* CPSI * SSTR2P + 32.* SSCHLP - 4.* SSTR3P
     5     - 4.* (457.* XNN + 2742.* XNE + 6040.* XNSE + 6098.* XNSI
     6     + 1567.* XNFI - 2344.* XNFO - 1632.* XNT + 560.* XNS
     7     + 1488.* XN + 576.) / (9.* XNMS * XNT * XN1T * XN2T)
      GG1I = 9.* GG1C
      GG1F = 3./2.* GG1A + 2./3.* GG1B
      
      RETURN
      END
*     
*     
*...  Initialization of support points in n-space and weights for the 
*     Gauss quadrature and of the anomalous dimensions for the RG 
*     evolution at these n-values. The outputs are written into the 
*     common-blocks CONT, MOMS, WEIGHTS and ANOMS, respectively. 
*     


      SUBROUTINE INIT
      implicit double precision (A - Z)
      INTEGER I1, I2, I3, K
      dimension ZS(8), WZ(8), DOWN(17), UP(17), WN(136)
      COMPLEX*16 ANS(2,136), AM(2,136), AP(2,136), AL(2,136),
     1     BE(2,136), AB(2,136), AC(2,136), RMIN(2,136),
     2     RPLUS(2,136), N(2,136), C1Q(2,136), C1G(2,136), C3Q(2,136),
     3     RMMQQ(2,136),RMMQG(2,136),RMMGQ(2,136),RMMGG(2,136),
     4     RMPQQ(2,136),RMPQG(2,136),RMPGQ(2,136),RMPGG(2,136),
     5     RPMQQ(2,136),RPMQG(2,136),RPMGQ(2,136),RPMGG(2,136),
     6     RPPQQ(2,136),RPPQG(2,136),RPPGQ(2,136),RPPGG(2,136),
     7     CC(2)
C ADDED
      complex*16 QQIM(2,136),QGFM(2,136),GQIM(2,136),GGIM(2,136),
     #      GGFM(2,136),NS1MIM(2,136),NS1PIM(2,136),NS1FM(2,136),
     #      QQ1FM(2,136),QG1FM(2,136),GQ1IM(2,136),GQ1FM(2,136),
     #      GG1IM(2,136),GG1FM(2,136)
C
      common/anomdim/QQIM,QGFM,GQIM,GGIM,GGFM,NS1MIM,NS1PIM,NS1FM,
     #      QQ1FM,QG1FM,GQ1IM,GQ1FM,GG1IM,GG1FM
C      
      common / WEIGHTS / WN
      common / MOMS / N
      common / CONTin / C, CC
      common / ANOMS / ANS, AM, AP, AL, BE, AB, AC, RMIN, RPLUS, C1Q,
     1     C1G, C3Q, RMMQQ, RMMQG, RMMGQ, RMMGG, RMPQQ, RMPQG, RMPGQ,
     2     RMPGG, RPMQQ, RPMQG, RPMGQ, RPMGG, RPPQQ, RPPQG, RPPGQ,
     3     RPPGG
C ADDED
      complex*16 H2ggM(2,136),H2gqM(2,136),H2qqM(2,136)
      common / H2mat /H2ggM,H2gqM,H2qqM
C
      save / WEIGHTS /
      save / MOMS /
      save / CONTin /
      save / ANOMS /
      save / H2mat /

*...  WEIGHTS AND SUPPORT POINTS FOR NOMALIZED 8 POINT GAUSS QUADRATURE :
      data WZ
     1     / 0.10122 85362 90376,  0.22238 10344 53374, 
     2     0.31370 66458 77887,  0.36268 37833 78362, 
     3     0.36268 37833 78362,  0.31370 66458 77887,
     4     0.22238 10344 53374,  0.10122 85362 90376/
      data ZS
     1     /-0.96028 98564 97536, -0.79666 64774 13627, 
     2     -0.52553 24099 16329, -0.18343 46424 95650, 
     3     0.18343 46424 95650,  0.52553 24099 16329,
     4     0.79666 64774 13627,  0.96028 98564 97536/
*...  INTEGRATION CONTOUR PARAMETERS :
      data DOWN / 0.D0, 0.5D0, 1.D0, 2.D0, 3.D0, 4.D0, 6.D0, 8.D0,
     1     1.D1, 1.2D1, 1.5D1, 1.8D1, 2.1D1, 2.4D1, 2.7D1, 3.D1, 3.3D1/
      
! Modifica Massimiliano  19/12/07
!      C = 1.8
      C=1.0
      
      PHI = 3.141592654 * 3./4.
      CO = DCOS (PHI)
      SI = DSIN (PHI)
      CC(1) = DCMPLX (CO, SI)
      CC(2) = DCMPLX (CO,-SI)
      DO 1 I1 = 1, 16
         UP(I1) = DOWN(I1+1)
 1    CONTINUE
      UP(17) = 36.D0 
*...  SUPPORT POINTS AND WEIGHTS FOR THE GAUSS INTEGRATION : 
*     (THE FACTOR (UP-DOWN)/2 IS INCLUDED IN THE WEIGHTS)
      K = 0
      DO 2 I2 = 1, 17
         SUM  = UP(I2) + DOWN(I2) 
         DIFF = UP(I2) - DOWN(I2) 
         DO 3 I3 = 1, 8
            K = K + 1
            Z = 0.5 * (SUM + DIFF * ZS(I3))
            WN(K) = DIFF / 2.* WZ(I3) 
            N(1,K)  = DCMPLX (C+CO*Z+1.,SI*Z)
            N(2,K)  = DCMPLX (C+CO*Z+1.,-SI*Z)
            
 3       CONTINUE
 2    CONTINUE 
      CALL ANOM (ANS, AM, AP, AL, BE, AB, AC, RMIN, RPLUS,
     1     RMMQQ, RMMQG, RMMGQ, RMMGG, RMPQQ, RMPQG, RMPGQ, RMPGG,
     2     RPMQQ, RPMQG, RPMGQ, RPMGG, RPPQQ, RPPQG, RPPGQ, RPPGG,
     3     C1Q, C1G, C3Q,
     #     QQIM,QGFM,GQIM,GGIM,GGFM,NS1MIM,NS1PIM,NS1FM,
     #     QQ1FM,QG1FM,GQ1IM,GQ1FM,GG1IM,GG1FM,N)


C    ADDED
C
      call H2matrix(H2ggM,H2gqM,H2qqM,N)

      RETURN
      END 


C ADDED: Computes and stores matrix of C2 values

      subroutine H2matrix(H2ggM,H2gqM,H2qqM,N)
      implicit none
      integer KP,K1
      complex*16 H2gg,H2gq,H2qq,XN
      complex*16 H2ggM(2,136),H2gqM(2,136),H2qqM(2,136),N(2,136)
      complex*16 MellinH2gg,MellinH2gq,MellinH2qq

      do KP=1,2
      do K1=1,136
      XN=N(KP,K1)
      H2ggM(KP,K1)=MellinH2gg(XN)
      H2gqM(KP,K1)=MellinH2gq(XN)
      H2qqM(KP,K1)=MellinH2qq(XN)
      enddo
      enddo

      return
      end


C     Computes the moments of L^(0)
      SUBROUTINE F0MOMENTS(N,UV,DV,US,DS,SS,GL,CH,BO)
      implicit COMPLEX*16 (A-Z)
      double precision A1UV,A2UV,A3UV,A4UV,A5UV,A6UV,A7UV,A8UV
      double precision A1DV,A2DV,A3DV,A4DV,A5DV,A6DV,A7DV,A8DV
      double precision A1US,A2US,A3US,A4US,A5US,A6US,A7US,A8US 
      double precision A1DS,A2DS,A3DS,A4DS,A5DS,A6DS,A7DS,A8DS
      double precision A1SS,A2SS,A3SS,A4SS,A5SS,A6SS,A7SS,A8SS
      double precision A1GL,A2GL,A3GL,A4GL,A5GL,A6GL,A7GL,A8GL
      double precision A1CH,A2CH,A3CH,A4CH,A5CH,A6CH,A7CH,A8CH
      double precision A1BO,A2BO,A3BO,A4BO,A5BO,A6BO,A7BO,A8BO
      double precision GE,ZETA2,ZETA3,PI,aa 
      common/ CUV/ A1UV,A2UV,A3UV,A4UV,A5UV,A6UV,A7UV,A8UV
      common/ CDV/ A1DV,A2DV,A3DV,A4DV,A5DV,A6DV,A7DV,A8DV
      common/ CUS/ A1US,A2US,A3US,A4US,A5US,A6US,A7US,A8US
      common/ CDS/ A1DS,A2DS,A3DS,A4DS,A5DS,A6DS,A7DS,A8DS
      common/ CSS/ A1SS,A2SS,A3SS,A4SS,A5SS,A6SS,A7SS,A8SS
      common/ CGL/ A1GL,A2GL,A3GL,A4GL,A5GL,A6GL,A7GL,A8GL
      common/ CCH/ A1CH,A2CH,A3CH,A4CH,A5CH,A6CH,A7CH,A8CH
      common/ CBO/ A1BO,A2BO,A3BO,A4BO,A5BO,A6BO,A7BO,A8BO
      common/expp/aa
      data GE/ 0.577216d0 /
      data ZETA2/1.64493d0 / 
      data ZETA3/1.20206d0 /
      data PI/3.1415926536d0 /
      AAA=CMPLX(aa,0.d0)
C     U-VALENCE
      CA1=CMPLX(A1UV,0.d0)
      CA2=CMPLX(A2UV,0.d0)
      CA3=CMPLX(A3UV,0.d0)
      CA4=CMPLX(A4UV,0.d0)
      CA5=CMPLX(A5UV,0.d0)
      CA6=CMPLX(A6UV,0.d0)
      CA7=CMPLX(A7UV,0.d0)
      CA8=CMPLX(A8UV,0.d0)
      UV=CA1*CBETA(CA2-1+N,CA3+1) + CA1*CA4*CBETA(CA2+N,CA3+1)+
     .     CA1*CA5*CBETA(CA2-1+0.5+N,CA3+1)
     .     + CA1*CA6*CBETA(CA2-1+1.5+N,CA3+1) 
     .     + CA1*CA7*CBETA(CA2-1+2+N,CA3+1)
     .     + CA1*CA8*CBETA(CA2-1+AAA+N,CA3+1)
C     D-VALENCE
      CA1=CMPLX(A1DV,0.d0)
      CA2=CMPLX(A2DV,0.d0)
      CA3=CMPLX(A3DV,0.d0)
      CA4=CMPLX(A4DV,0.d0)
      CA5=CMPLX(A5DV,0.d0)
      CA6=CMPLX(A6DV,0.d0)
      CA7=CMPLX(A7DV,0.d0)
      CA8=CMPLX(A8DV,0.d0)
      DV=CA1*CBETA(CA2-1+N,CA3+1) + CA1*CA4*CBETA(CA2+N,CA3+1)+
     .     CA1*CA5*CBETA(CA2-1+0.5+N,CA3+1) 
     .     + CA1*CA6*CBETA(CA2-1+1.5+N,CA3+1) 
     .     + CA1*CA7*CBETA(CA2-1+2+N,CA3+1)
     .     + CA1*CA8*CBETA(CA2-1+AAA+N,CA3+1)
C     U-BAR
      CA1=CMPLX(A1US,0.d0)
      CA2=CMPLX(A2US,0.d0)
      CA3=CMPLX(A3US,0.d0)
      CA4=CMPLX(A4US,0.d0)
      CA5=CMPLX(A5US,0.d0)
      CA6=CMPLX(A6US,0.d0)
      CA7=CMPLX(A7US,0.d0)
      CA8=CMPLX(A8US,0.d0)
      US=CA1*CBETA(CA2-1+N,CA3+1) + CA1*CA4*CBETA(CA2+N,CA3+1)+
     .     CA1*CA5*CBETA(CA2-1+0.5+N,CA3+1) 
     .     + CA1*CA6*CBETA(CA2-1+1.5+N,CA3+1) 
     .     + CA1*CA7*CBETA(CA2-1+2+N,CA3+1)
     .     + CA1*CA8*CBETA(CA2-1+AAA+N,CA3+1)
C     D-BAR
      CA1=CMPLX(A1DS,0.d0)
      CA2=CMPLX(A2DS,0.d0)
      CA3=CMPLX(A3DS,0.d0)
      CA4=CMPLX(A4DS,0.d0)
      CA5=CMPLX(A5DS,0.d0)
      CA6=CMPLX(A6DS,0.d0)
      CA7=CMPLX(A7DS,0.d0)
      CA8=CMPLX(A8DS,0.d0)
      DS=CA1*CBETA(CA2-1+N,CA3+1) + CA1*CA4*CBETA(CA2+N,CA3+1)+
     .     CA1*CA5*CBETA(CA2-1+0.5+N,CA3+1) 
     .     + CA1*CA6*CBETA(CA2-1+1.5+N,CA3+1) 
     .     + CA1*CA7*CBETA(CA2-1+2+N,CA3+1)
     .     + CA1*CA8*CBETA(CA2-1+AAA+N,CA3+1)
C     S-BAR
      CA1=CMPLX(A1SS,0.d0)
      CA2=CMPLX(A2SS,0.d0)
      CA3=CMPLX(A3SS,0.d0)
      CA4=CMPLX(A4SS,0.d0)
      CA5=CMPLX(A5SS,0.d0)
      CA6=CMPLX(A6SS,0.d0)
      CA7=CMPLX(A7SS,0.d0)
      CA8=CMPLX(A8SS,0.d0)
      SS=CA1*CBETA(CA2-1+N,CA3+1) + CA1*CA4*CBETA(CA2+N,CA3+1)+
     .     CA1*CA5*CBETA(CA2-1+0.5+N,CA3+1)
     .     + CA1*CA6*CBETA(CA2-1+1.5+N,CA3+1) 
     .     + CA1*CA7*CBETA(CA2-1+2+N,CA3+1)
     .     + CA1*CA8*CBETA(CA2-1+AAA+N,CA3+1)
C     GLUON
      CA1=CMPLX(A1GL,0.d0)
      CA2=CMPLX(A2GL,0.d0)
      CA3=CMPLX(A3GL,0.d0)
      CA4=CMPLX(A4GL,0.d0)
      CA5=CMPLX(A5GL,0.d0)
      CA6=CMPLX(A6GL,0.d0)
      CA7=CMPLX(A7GL,0.d0)
      CA8=CMPLX(A8GL,0.d0)
      GL=CA1*CBETA(CA2-1+N,CA3+1) + CA1*CA4*CBETA(CA2+N,CA3+1)+
     .     CA1*CA5*CBETA(CA2-1+0.5+N,CA3+1)
     .     + CA1*CA6*CBETA(CA2-1+1.5+N,CA3+1) 
     .     + CA1*CA7*CBETA(CA2-1+2+N,CA3+1)
     .     + CA1*CA8*CBETA(CA2-1+AAA+N,CA3+1)
C     CHARM
      CA1=CMPLX(A1CH,0.d0)
      CA2=CMPLX(A2CH,0.d0)
      CA3=CMPLX(A3CH,0.d0)
      CA4=CMPLX(A4CH,0.d0)
      CA5=CMPLX(A5CH,0.d0)
      CA6=CMPLX(A6CH,0.d0)
      CA7=CMPLX(A7CH,0.d0)
      CA8=CMPLX(A8CH,0.d0)
      CH=CA1*CBETA(CA2-1+N,CA3+1) + CA1*CA4*CBETA(CA2+N,CA3+1)+
     .     CA1*CA5*CBETA(CA2-1+0.5+N,CA3+1)
     .     + CA1*CA6*CBETA(CA2-1+1.5+N,CA3+1) 
     .     + CA1*CA7*CBETA(CA2-1+2+N,CA3+1)
     .     + CA1*CA8*CBETA(CA2-1+AAA+N,CA3+1)
C     BOTTOM
      CA1=CMPLX(A1BO,0.d0)
      CA2=CMPLX(A2BO,0.d0)
      CA3=CMPLX(A3BO,0.d0)
      CA4=CMPLX(A4BO,0.d0)
      CA5=CMPLX(A5BO,0.d0)
      CA6=CMPLX(A6BO,0.d0)
      CA7=CMPLX(A7BO,0.d0)
      CA8=CMPLX(A8BO,0.d0)
      BO=CA1*CBETA(CA2-1+N,CA3+1) + CA1*CA4*CBETA(CA2+N,CA3+1)+
     .     CA1*CA5*CBETA(CA2-1+0.5+N,CA3+1)
     .     + CA1*CA6*CBETA(CA2-1+1.5+N,CA3+1) 
     .     + CA1*CA7*CBETA(CA2-1+2+N,CA3+1)
     .     + CA1*CA8*CBETA(CA2-1+AAA+N,CA3+1)
      RETURN
      END               
      
      
      
      function xqq(xx)
      implicit double precision (A-H,J-Z)
      real aa,uv0
c     common/x/x,Q2,Q20
      common/ippbar/ippbar
      include 'scales.F'
      z=xtau**xx
      y=xtau/z
      call distri(z,sqrt(q2),unv1,dnv1,us1,ds1,str1,chm,bot,glu1,ippbar)
      call distri(y,sqrt(q2),unv2,dnv2,us2,ds2,str2,chm,bot,glu2,ippbar)
      xqq=((unv1+us1)*us2 + (dnv1+ds1)*ds2 + str1*str2) +
     .     ((unv2+us2)*us1 + (dnv2+ds2)*ds1 + str2*str1)  
      return
      end
      
      function xqg(xx)
      implicit double precision (A-H,J-Z)
      real aa,uv0
c     common/x/x,Q2,Q20
      common/ippbar/ippbar
      include 'scales.F'
      z=xtau**xx
      y=xtau/z
      q=sqrt(q2)
      call distri(z,q,unv1,dnv1,us1,ds1,str1,chm1,bot1,glu1,ippbar)
      call distri(y,q,unv2,dnv2,us2,ds2,str2,chm2,bot2,glu2,ippbar)
      xqg=((unv1+2*us1)+(dnv1+2*ds1)+ 2*str1+0*chm1+0*bot1)*glu2 +
     .     ((unv2+2*us2)+(dnv2+2*ds2)+ 2*str2+0*chm2+0*bot2)*glu1
      return
      end
      
      function xgg(xx)
      implicit double precision (A-H,J-Z)
      real aa,uv0
c     common/x/x,Q2,Q20
      common/ippbar/ippbar
      include 'scales.F'
      z=xtau**xx
      y=xtau/z
      call distri(z,sqrt(q2),unv1,dnv1,us1,ds1,str1,chm,bot,glu1,ippbar)
      call distri(y,sqrt(q2),unv2,dnv2,us2,ds2,str2,chm,bot,glu2,ippbar)
      xgg=glu1*glu2
      return
      end

!-----------------------------------------------------------------------
! Written by D.Tommasini: Mellin of H2gg and H2gq
!-----------------------------------------------------------------------

      COMPLEX*16 FUNCTION MellinH2gq(ZN)
      IMPLICIT none
      COMPLEX*16 ZN,ZN2,ZN3,ZN4,ZN5,ZN6,ZNp1,ZNm1,ZNp2,ZNp3
      COMPLEX*16 ZNb,plinomZN,plinomZNm, temp
      COMPLEX*16 PG0ZN,PG1ZN,PG1ZNmez,PG1ZNp1mez,PG2ZN,PG2ZNmezp1
      COMPLEX*16 PG2ZNp1mez,PG0ZNmezp1,PG0ZNp1mez
      COMPLEX*16 DACG1ZN,DACG1ZNp1,ACG4ZN,ACG4ZNp1,ACG4ZNp2,DACG1ZNp2
      COMPLEX*16 DACG1,ACG4
      Real*8 GE,GE2,GE3,ZET3,Pi,Pi2,DL
      Real*8 CA,CF,nf
      CA=3.0D0
      CF=4.0D0/3.0D0
      nf=5.0D0
      GE   = 0.57721566490153D0
      GE2  =  GE*GE
      GE3  =  GE2*GE
      ZET3 = 1.20205690315959428540D0
      Pi = 3.141592653589793238462643D0
      Pi2 = Pi*Pi
      DL = DLOG(2.0D0)

      ZN2=ZN*ZN
      ZN3=ZN2*ZN
      ZN4=ZN3*ZN
      ZN5=ZN4*ZN
      ZN6=ZN5*ZN
      ZNp1 = ZN+(1.0D0,0.0D0)
      ZNm1 = ZN+(-1.0D0,0.0D0)
      ZNp2 = ZN+(2.0D0,0.0D0)
      ZNp3 = ZN+(3.0D0,0.0D0)
      ZNb = (-1.0D0 + ZN2)
      plinomZN = 2.0D0 + ZN + ZN2
      plinomZNm=-2.0D0 + ZN + ZN2
!      write(*,*) ZN,ZN2,zn3,zn4,zn5,zn6,znp1,znm1,znp2,znp3,znb
!----------------------------------------------------------------------
      CALL PSI0(ZN,PG0ZN)                        !PolyGamma[0, ZN]
      temp=ZN/2.0D0+1.0D0
      CALL PSI0(temp,PG0ZNmezp1)                 !PolyGamma[0, 1 + ZN/2]
      temp=ZNp1/2.0D0
      CALL PSI0(temp,PG0ZNp1mez)                 !PolyGamma[0, ZNp1/2]
      CALL PSI1(ZN,PG1ZN)                        !PolyGamma[1, ZN]
      temp=ZN/2.0D0
      CALL PSI1(temp,PG1ZNmez)                   !PolyGamma[1, ZN/2]
      temp=ZNp1/2.0D0
      CALL PSI1(temp,PG1ZNp1mez)                 !PolyGamma[1, ZNp1/2]
      CALL PSI2(ZN,PG2ZN)                        !PolyGamma[2, ZN]
      temp=ZN/2.0D0+1.0D0
      CALL PSI2(temp,PG2ZNmezp1)                 !PolyGamma[2, 1 + ZN/2]
      temp=ZNp1/2.0D0
      CALL PSI2(temp,PG2ZNp1mez)                 !PolyGamma[2, ZNp1/2]
!--------------------------------------------------------------------
! NB:The Blumlein's convention to calculate the Mellin transform is z^N (not z^(N-1) )
!--------------------------------------------------------------------
      DACG1ZN=DACG1(ZNm1)
      DACG1ZNp1=DACG1(ZN)
      DACG1ZNp2=DACG1(ZNp1)
      ACG4ZN=ACG4(ZNm1)
      ACG4ZNp1=ACG4(ZN)
      ACG4ZNp2=ACG4(ZNp1)

!      write (*,*) 'H2gqZN',ZN

      MellinH2gq=(CF*((4.0D0*ZN*nf*(72.0D0 + ZN*ZNp1*(-12.0D0 +
     # 18.0D0*GE2*ZN*ZNp1*plinomZN +
     # 4.0D0*ZN*(-13.0D0 + ZN*(144.0D0 + ZN*(74.0D0 + 13.0D0*ZN))) -
     # 12.0D0*GE*
     # (-6.0D0+ZN*(7.0D0 + ZN*(27.0D0 +2.0D0*ZN*(5.0D0+ZN))))+3.0D0*
     # ZN*ZNp1*plinomZN*
     #    Pi2) + 6.0D0*ZN*ZNp1*
     #  (2.0D0*(6.0D0 + ZN*(-7.0D0 + 3.0D0*GE*ZNp1*plinomZN -
     #       ZN*(27.0D0 + 2.0D0*ZN*(5.0D0 + ZN))))*PG0ZN +
     #   3.0D0*ZN*ZNp1*plinomZN*PG0ZN**2 -
     #   3.0D0*ZN*ZNp1*plinomZN*PG1ZN)))/(ZNm1*ZNp1**3) +
     #    (CA*(8.0D0*ZNm1**2*(5184.0D0 +
     #  ZN*(24048.0D0 + ZN*(36744.0D0 + ZN*(19532.0D0+ZN*(20342.0D0+
     # ZN*(57126.0D0 +
     # ZN*(48911.0D0 + ZN*(-18621.0D0 +ZN*(-54864.0D0+ZN*(-34946.0D0+
     # ZN*(-9377.0D0 +
     #                     ZN*(-807.0D0 + 40.0D0*ZN)))))))))))) -
     # 3.0D0*ZN*ZNp1*plinomZNm**2*(576.0D0*ZN3*ZNp1**3*plinomZNm*
     #    DACG1ZN + 576.0D0*ZNm1**2*ZN2*ZNp1**3*ZNp2*DACG1ZNp1 -
     #   576.0D0*PG0ZN + ZN*(288.0D0*ZN2*ZNp2*ZNb**2*DACG1ZNp2 +
     #     ZNp1*(576.0D0*ZN2*ZNp1**2*plinomZNm*ACG4ZN +
     #       576.0D0*ZN*ZNp2*ZNb**2*ACG4ZNp1 + 576.0D0*ZN2*ACG4ZNp2 -
     #    288.0D0*ZN3*ACG4ZNp2 - 864.0D0*ZN4*ACG4ZNp2 +288.0D0*ZN5*
     # ACG4ZNp2+
     #       288.0D0*ZN6*ACG4ZNp2 - 8.0D0*(240.0D0 + ZN*(-484.0D0 +
     # ZN*(-721.0D0 +
     #        ZN*(211.0D0 + ZN*(462.0D0 + ZN*(213.0D0 + 43.0D0*
     # ZN))))))*PG0ZN +
     #       12.0D0*ZN*ZNp1*ZNp2*(-28.0D0 + ZN*(29.0D0 + ZN*(18.0D0 +
     # 5.0D0*ZN)))*
     #        PG0ZN**2 - 24.0D0*ZNm1*ZN*ZNp1*ZNp2*plinomZN*
     #     PG0ZN**3 - 36.0D0*ZNp2*(-4.0D0 + ZN2*(15.0D0 +ZN*(3.0D0 +
     # ZN + ZN2)))*
     #        PG1ZNmez - 12.0D0*(-48.0D0 + ZN*ZNp1*(-248.0D0 +
     # ZN*(-30.0D0 + ZN*(17.0D0 + ZN*(52.0D0 + 17.0D0*ZN)))))*PG1ZN -
     #       72.0D0*ZNm1*ZN*ZNp1*ZNp2*plinomZN*PG0ZN*
     # PG1ZN + 36.0D0*ZNp2*(-4.0D0 + ZN2*(15.0D0 + ZN*(3.0D0 +ZN+
     # ZN2)))*
     #        PG1ZNp1mez + 9.0D0*ZNm1*ZN*ZNp1*ZNp2*
     #        plinomZN*PG2ZNmezp1 + 120.0D0*ZNm1*ZN*ZNp1*
     #        ZNp2*plinomZN*PG2ZN - 9.0D0*ZNm1*ZN*ZNp1*
     #        ZNp2*plinomZN*PG2ZNp1mez)))))/
     # (ZNp2**3*ZNb**4) +
     # (6.0D0*CA*ZN*(2.0D0*GE*(72.0D0 + ZN*(204.0D0 + ZN*(-346.0D0 +
     # 6.0D0*GE2*ZNp1**2*
     #        (-2.0D0 + ZN + ZN3) - 3.0D0*GE*ZNp1**2*
     #        (-28.0D0 + ZN*(29.0D0 + ZN*(18.0D0 + 5.0D0*ZN))) +
     #    2.0D0*ZN*(-516.0D0 + ZN*(3.0D0 + ZN*(335.0D0+ZN*(170.0D0+
     # 43.0D0*ZN))))))) +
     # ZN*ZNp1**2*(-24.0D0 + ZN*(76.0D0 + 6.0D0*GE*(-2.0D0 + ZN + ZN3) +
     #     ZN*(85.0D0 + ZN*(-174.0D0 + 61.0D0*ZN))))*Pi2 + 6.0D0*ZN2*
     # ZNp1**2*
     #  ((-2.0D0 + ZN + ZN3)*Pi2*PG0ZNmezp1 +
     #   (2.0D0*GE*(28.0D0 + 3.0D0*GE*(-2.0D0 + ZN + ZN3) -
     # ZN*(29.0D0 + ZN*(18.0D0 + 5.0D0*ZN))) + (-2.0D0+ ZN+ZN3)*Pi2)*
     # PG0ZN +
     #   6.0D0*GE*(-2.0D0 + ZN + ZN3)*PG0ZN**2 -
     #   (-2.0D0 + ZN + ZN3)*(Pi2*PG0ZNp1mez - 6.0D0*GE*
     #      PG1ZN - 46.0D0*ZET3))))/(ZNm1**2*ZNp1**3) -
     # (18.0D0*CF*(24.0D0 - 2.0D0*ZN*(-2.0D0*GE3*ZN2*ZNp1**3*plinomZN -
     # 6.0D0*GE*ZNp1*(4.0D0 + ZN2*ZNp1*(1.0D0 + ZN*(36.0D0 +5.0D0*ZN*
     # (4.0D0+ZN)))) +
     # 3.0D0*GE2*ZN*ZNp1**2*(-4.0D0 + ZN*(4.0D0 + ZN*(17.0D0 + ZN*
     # (6.0D0+ZN)))) +
     # 3.0D0*(8.0D0 + ZN*(3.0D0 + ZN*(-99.0D0 + ZN*(-107.0D0 + ZN*
     # (63.0D0 + ZN*(113.0D0 +
     # (26.0D0 - 3.0D0*ZN)*
     #                ZN))))))) + ZN2*ZNp1**2*
     # (4.0D0 + ZN*(-4.0D0 + 2.0D0*GE*ZNp1*plinomZN - ZN*(17.0D0+ ZN*
     # (6.0D0+ ZN))))*
     #  Pi2 + 2.0D0*ZN*ZNp1*(3.0D0*ZN*ZNp1*
     #    (4.0D0 + ZN*(-4.0D0 + 2.0D0*GE*ZNp1*plinomZN -
     # ZN*(17.0D0 + ZN*(6.0D0 + ZN))))*PG0ZN**2 + 2.0D0*ZN2*ZNp1**2*
     #    plinomZN*PG0ZN**3 + PG0ZN*
     #    (24.0D0 + ZN*ZNp1*(6.0D0*GE2*ZN*ZNp1*plinomZN +
     #       6.0D0*ZN*(1.0D0 + ZN*(36.0D0 + 5.0D0*ZN*(4.0D0 + ZN))) -
     # 6.0D0*GE*
     #  (-4.0D0 + ZN*(4.0D0 + ZN*(17.0D0 + ZN*(6.0D0 + ZN)))) + ZN*
     # ZNp1*plinomZN*
     #        Pi2) - 6.0D0*ZN2*ZNp1**2*plinomZN*PG1ZN) +
     # ZN*ZNp1*(-3.0D0*(4.0D0 + ZN*(-4.0D0 + 2.0D0*GE*ZNp1*plinomZN -
     #         ZN*(17.0D0 + ZN*(6.0D0 + ZN))))*PG1ZN + 2.0D0*ZN*ZNp1*
     #      plinomZN*(PG2ZN + 2.0D0*ZET3)))))/
     # (ZNm1*ZNp1**4)))/(1728.0D0*ZN4)

      RETURN
      END

      COMPLEX*16 FUNCTION MellinH2gg(ZN)
      IMPLICIT none
      COMPLEX*16 ZN,ZN2,ZN3,ZN4,ZN5,ZN6,ZNp1,ZNm1,ZNp2,ZNp3,ZNm2
      COMPLEX*16 ZNb,plinomZN,plinomZNm, temp
      COMPLEX*16 PG0ZNm1,PG0ZNmez,PG1ZNm1,PG0ZNp1,PG0ZNp2mez
      COMPLEX*16 PG1ZNmezp1,PG1ZNp1,PG0ZNp2,PG1ZNp2,PG0ZNp3,PG1ZNp3
      COMPLEX*16 PG2ZNm1,PG2ZNmez,PG2ZNp3,PG3ZNmezp1,PG3ZNm1,PG3ZN
      COMPLEX*16 PG3ZNp1mez,PG0ZNm1mez
      COMPLEX*16 PG0ZN,PG1ZN,PG1ZNmez,PG1ZNp1mez,PG2ZN,PG2ZNmezp1
      COMPLEX*16 PG2ZNp1mez,PG0ZNmezp1,PG0ZNp1mez
      COMPLEX*16 DDACG1ZNm1,DACG2ZNm1,DACG2ZN,DACG2ZNp1,DACG2ZNp2
      COMPLEX*16 DACG2ZNp3,DACG4ZNm1,ACG1ZNp1,ACG1ZNp2,ACG1ZNp3
      COMPLEX*16 ACG13ZNm1,ACG13ZN,ACG13ZNp1,ACG13ZNp2,ACG13ZNp3
      COMPLEX*16 ACG20ZNm1,ACG22ZNm1,ACG6ZNm1,ACG7ZNm1,ACG9ZNm1
      COMPLEX*16 DDACG1,DACG2,DACG4,ACG1,ACG13,ACG20,ACG22,ACG5ZNm1
      COMPLEX*16 ACG1ZNm1
      COMPLEX*16 ACG5,ACG6,ACG7,ACG9
      Real*8 GE,GE2,GE3,ZET3,Pi,Pi2,DL
      Real*8 CA,CF,nf,LT
      common/LT/LT
      CA=3.0D0
      CF=4.0D0/3.0D0
      nf=5.0D0
      GE   = 0.57721566490153
      GE2  =  GE*GE
      GE3  =  GE2*GE
      ZET3 = 1.20205690315959428540
      Pi = 3.141592653589793238462643
      Pi2 = Pi*Pi
      DL = DLOG(2.0D0)
!      write (*,*) 'H2ggZN',ZN
      ZN2=ZN*ZN
      ZN3=ZN2*ZN
      ZN4=ZN3*ZN
      ZN5=ZN4*ZN
      ZN6=ZN5*ZN
      ZNp1 = ZN+(1.0D0,0.0D0)
      ZNm1 = ZN+(-1.0D0,0.0D0)
      ZNp2 = ZN+(2.0D0,0.0D0)
      ZNm2 = ZN+(-2.0D0,0.0D0)
      ZNp3 = ZN+(3.0D0,0.0D0)
      ZNb = (-1.0D0 + ZN2)
      plinomZN = 2.0D0 + ZN + ZN2
      plinomZNm=-2.0D0 + ZN + ZN2
!      write(*,*) zn,zn2,zn3,zn4,zn5,zn6,znp1,znm1,znp2,znp3,znb
!------------------------------------------------------------------
      CALL PSI0(ZN,PG0ZN)                        !PolyGamma[0, ZN]
      temp=ZN/2.0D0+1.0D0
      CALL PSI0(temp,PG0ZNmezp1)                 !PolyGamma[0, 1 + ZN/2]
      temp=ZNp1/2.0D0
      CALL PSI0(temp,PG0ZNp1mez)                 !PolyGamma[0, ZNp1/2]
      CALL PSI1(ZN,PG1ZN)                        !PolyGamma[1, ZN]
      temp=ZN/2.0D0
      CALL PSI1(temp,PG1ZNmez)                   !PolyGamma[1, ZN/2]
      temp=ZNp1/2.0D0
      CALL PSI1(temp,PG1ZNp1mez)                 !PolyGamma[1, ZNp1/2]
      CALL PSI2(ZN,PG2ZN)                        !PolyGamma[2, ZN]
      temp=ZN/2.0D0+1.0D0
      CALL PSI2(temp,PG2ZNmezp1)                 !PolyGamma[2, 1 + ZN/2]
      temp=ZNp1/2.0D0
      CALL PSI2(temp,PG2ZNp1mez)                 !PolyGamma[2, ZNp1/2]
      temp=ZNm1
      CALL PSI0(temp,PG0ZNm1)                    !PolyGamma[0, -1 + ZN]
      temp=ZN/2.0D0
      CALL PSI0(temp,PG0ZNmez)                   !PolyGamma[0, ZN/2]
      temp=ZNm1
      CALL PSI1(temp,PG1ZNm1)                    !PolyGamma[1, -1 + ZN]
      temp=ZNp1
      CALL PSI0(temp,PG0ZNp1)                    !PolyGamma[0, 1 + ZN]
      temp=ZNp2/2.0D0
      CALL PSI0(temp,PG0ZNp2mez)                 !PolyGamma[0, ZNp2/2]
      temp=ZN/2.0D0+1.0D0
      CALL PSI1(temp,PG1ZNmezp1)                 !PolyGamma[1, 1 + ZN/2]
      CALL PSI1(ZNp1,PG1ZNp1)                    !PolyGamma[1, 1 + ZN]
      CALL PSI0(ZNp2,PG0ZNp2)                    !PolyGamma[0, 2 + ZN]
      CALL PSI1(ZNp2,PG1ZNp2)                    !PolyGamma[1, 2 + ZN]
      CALL PSI0(ZNp3,PG0ZNp3)                    !PolyGamma[0, 3 + ZN]
      CALL PSI1(ZNp3,PG1ZNp3)                    !PolyGamma[1, 3 + ZN]
      CALL PSI2(ZNm1,PG2ZNm1)                    !PolyGamma[2, -1 + ZN]
      temp=ZN/2.0D0
      CALL PSI2(temp,PG2ZNmez)                   !PolyGamma[2, ZN/2]
      CALL PSI2(ZNp3,PG2ZNp3)                    !PolyGamma[2, 3 + ZN]
      temp=ZN/2.0D0+1.0D0
      CALL PSI3(temp,PG3ZNmezp1)                 !PolyGamma[3, 1 + ZN/2]
      CALL PSI3(ZNm1,PG3ZNm1)                    !PolyGamma[3, -1 + ZN]
      CALL PSI3(ZN,PG3ZN)                        !PolyGamma[3, ZN]
      temp=ZNp1/2.0D0
      CALL PSI3(temp,PG3ZNp1mez)                 !PolyGamma[3, ZNp1/2]
      temp=ZNm1/2.0D0
      CALL PSI0(temp,PG0ZNm1mez)                 !PolyGamma[0, ZNm1/2]
!--------------------------------------------------------------------
! NB:The Blumlein's convention to calculate the Mellin transform is z^N (not z^(N-1) )
!--------------------------------------------------------------------
      DDACG1ZNm1=DDACG1(ZNm2)                    !DDACG1[-1 + ZN]
      DACG2ZNm1=DACG2(ZNm2)                      !DACG2[-1 + ZN]
      DACG2ZN=DACG2(ZNm1)                        !DACG2[ZN]
      DACG2ZNp1=DACG2(ZN)                        !DACG2[1 + ZN]
      DACG2ZNp2=DACG2(ZNp1)                      !DACG2[2 + ZN]
      DACG2ZNp3=DACG2(ZNp2)                      !DACG2[3 + ZN]
      DACG4ZNm1=DACG4(ZNm2)                      !DACG4[-1 + ZN]
      DACG4ZNm1=DACG4(ZNm2)                      !DACG4[-1 + ZN]
      ACG1ZNp1=ACG1(ZN)                          !ACG1[1 + ZN]
      ACG1ZNm1=ACG1(ZNm2)                        !ACG1[-1 + ZN]
      ACG1ZNp2=ACG1(ZNp1)                        !ACG1[2 + ZN]
      ACG1ZNp3=ACG1(ZNp2)                        !ACG1[3 + ZN]
      ACG5ZNm1=ACG5(ZNm2)                        !ACG5[-1 + ZN]
      ACG6ZNm1=ACG6(ZNm2)                        !ACG6[-1 + ZN]
      ACG7ZNm1=ACG7(ZNm2)                        !ACG7[-1 + ZN]
      ACG9ZNm1=ACG9(ZNm2)                        !ACG9[-1 + ZN]
      ACG13ZNm1=ACG13(ZNm2)                      !ACG13[-1 + ZN]
      ACG13ZN=ACG13(ZNm1)                        !ACG13[ZN]
      ACG13ZNp1=ACG13(ZN)                        !ACG13[1 + ZN]
      ACG13ZNp2=ACG13(ZNp1)                      !ACG13[2 + ZN]
      ACG13ZNp3=ACG13(ZNp2)                      !ACG13[3 + ZN]
      ACG20ZNm1=ACG20(ZNm2)                      !ACG20[-1 + ZN]
      ACG22ZNm1=ACG22(ZNm2)                      !ACG22[-1 + ZN]


!      write(*,*) ZN,DDACG1ZNm1,DACG2ZNm1,DACG2ZN


      MellinH2gg=(9.0D0*CF**2)/4.0D0-(CA**2*(8.0D0+ ZN*(24.0D0 +
     # ZN*(31.0D0 + 2.0D0*ZN*
     # (9.0D0 + 2.0D0*ZN))))*DL**2)/
     # (ZN2*ZNp1**2*ZNp2**2) + CA**2*(-DDACG1ZNm1/2.0D0 + DACG2ZNm1 +
     # 2.0D0*DACG2ZN +3.0D0*DACG2ZNp1 +2.0D0*DACG2ZNp2 + DACG2ZNp3 +
     # DACG4ZNm1+
     # (4.0D0*ACG1ZNp1)/ZN2 + (2.0D0*ACG1ZNp2)/ZNp1**2 +
     # (2.0D0*ACG1ZNp3)/ZNp2**2 +
     # 2.0D0*ACG13ZNm1 + 4.0D0*ACG13ZN + 6.0D0*ACG13ZNp1 +
     # 4.0D0*ACG13ZNp2 +
     # 2.0D0*ACG13ZNp3 + 3.0D0*ACG20ZNm1 - 2.0D0*ACG22ZNm1 +
     # ACG5ZNm1 -
     # 2.0D0*ACG6ZNm1 - 3.0D0*ACG7ZNm1 + 2.0D0*ACG9ZNm1) +
     # (CA**2*(32.0D0*ZN*ZNb*(-54.0D0 +
     # ZN2*(162.0D0+ZN*(-47.0D0 + ZN*(-648.0D0 + ZN*(465.0D0 + ZN*
     # (-108.0D0+101.0D0*ZN*
     # (-3.0D0 + ZN2)))))))*PG0ZNm1 +
     # (186624.0D0+ZN*(312192.0D0 +ZN*(-463296.0D0 +ZN*(-948160.0D0 +
     # ZN*(707600.0D0 + 756.0D0*LT*ZNm1**4*ZNp1**4*ZNp2**4 +
     # ZN*(1232712.0D0 + ZN*(-1566444.0D0 + ZN*(-2736488.0D0 + ZN*
     # (204381.0D0 +
     # ZN*(2124608.0D0 + ZN*(390184.0D0 + ZN*(-1291248.0D0 + ZN*
     # (-807822.0D0 +
     # ZN*(86264.0D0 + ZN*(236748.0D0 + ZN*(83096.0D0 + 9561.0D0*
     # ZN))))))))))))))) + 432.0D0*ZNm1**4*ZN*ZNp1*ZNp2*
     # (16.0D0+ZN*(72.0D0+ZN*(132.0D0+ZN*(119.0D0 + ZN*(57.0D0 +
     # ZN*(15.0D0+2.0D0*ZN))))))*
     # PG0ZNmez)/ZNp2**4))/(864.0D0*ZN4*ZNb**4) +
     # (CA**2*Pi2*(-288.0D0- 816.0D0*ZN + 356.0D0*ZN2+3312.0D0*ZN3 +
     # 3986.0D0*ZN4 +
     # 1884.0D0*ZN5 +
     # 314.0D0*ZN6+52.0D0*ZN2*Pi2+156.0D0*ZN3*Pi2+169.0D0*ZN4*Pi2 +
     # 78.0D0*ZN5*Pi2+
     # 13.0D0*ZN6*Pi2 + 192.0D0*ZN*DL+720.0D0*ZN2*DL+936.0D0*ZN3*DL +
     # 504.0D0*ZN4*DL+96.0D0*ZN5*DL+24.0D0*ZN2*(2.0D0+3.0D0*ZN+
     # ZN2)**2*ACG1ZNm1-
     # 12.0D0*ZN*(8.0D0 + 22.0D0*ZN + 23.0D0*ZN2 + 11.0D0*ZN3 +
     # 2.0D0*ZN4)*PG0ZNmez +
     # 96.0D0*ZN*PG0ZNp1mez + 264.0D0*ZN2*PG0ZNp1mez +
     # 276.0D0*ZN3*PG0ZNp1mez + 132.0D0*ZN4*PG0ZNp1mez +
     # 24.0D0*ZN5*PG0ZNp1mez - 48.0D0*ZN2*PG1ZNm1 -
     # 144.0D0*ZN3*PG1ZNm1 - 156.0D0*ZN4*PG1ZNm1 -
     # 72.0D0*ZN5*PG1ZNm1 - 12.0D0*ZN6*PG1ZNm1))/
     # (144.0D0*ZN2*ZNp1**2*ZNp2**2) +
     # CA**2*(-1.0D0/(12.0D0*ZNp1**2)+(-2.0D0/ZN3 + ZNp1**(-3)-1.0D0/
     # (12.0D0*ZNp1) -
     # ZNp2**(-3))*PG0ZNp1 +
     # ((-2.0D0/ZN3 + ZNp1**(-3) - ZNp2**(-3))*PG0ZNp2mez)/2.0D0 -
     # (PG0ZNm1**2*PG1ZNm1)/2.0D0 + PG1ZNm1**2 -
     # ((8.0D0 + ZN*(24.0D0 + ZN*(23.0D0+2.0D0*ZN*(5.0D0 + ZN))))*
     # PG1ZNmez)/
     # (4.0D0*ZN2*ZNp1**2*ZNp2**2) + PG0ZN**2*PG1ZN) +
     # (CA**2*((-2.0D0/ZN2 + ZNp1**(-2) - ZNp2**(-2))*PG1ZNmezp1 +
     # 2.0D0*(-4.0D0*PG1ZN**2 + ((8.0D0+ZN*(24.0D0+ZN*(23.0D0 +
     # 2.0D0*ZN*(5.0D0 + ZN))))*
     # PG1ZNp1mez)/(ZN2*ZNp1**2*ZNp2**2) +
     # (2.0D0*(-8.0D0 + ZN2*(47.0D0 + ZN*(56.0D0+ZN*(30.0D0+7.0D0*
     # ZN))))*PG1ZNp1)/
     # (ZNm1*ZN2*ZNp1**2*ZNp2**2) - 3.0D0*PG0ZNp1**2*
     # PG1ZNp1 + 6.0D0*PG1ZNp1**2 +
     # 2.0D0*PG0ZNp2**2*PG1ZNp2 - 4.0D0*PG1ZNp2**2 -
     # PG0ZNp3**2*PG1ZNp3 + 2.0D0*PG1ZNp3**2 +
     # PG2ZNm1/ZNb)))/4.0D0 +
     # (CA**2*(-96.0D0*(4.0D0 + ZN2*(7.0D0 + 4.0D0*ZN))*PG2ZNm1 +
     # ZNm1*(12.0D0*(4.0D0 + ZN*(5.0D0 + 2.0D0*ZN))*PG2ZNmez -
     # 12.0D0*(4.0D0 + ZN*(5.0D0 + 2.0D0*ZN))*PG2ZNp1mez +
     # ZN*ZNp1*ZNp2*(192.0D0*PG0ZNp3*PG2ZNp3 +
     # PG3ZNmezp1 - 80.0D0*PG3ZNm1 +
     # 16.0D0*PG3ZN - PG3ZNp1mez))))/
     # (192.0D0*ZN*ZNp2*ZNb) +
     # (CA**2*(504.0D0+2.0D0*ZN*(299.0D0+ZN*(55.0D0 - ZN*(173.0D0 +
     # 55.0D0*ZN))) +
     # 9.0D0*ZNm1*ZN*ZNp1*ZNp2*(PG0ZNm1mez -
     # PG0ZNmez - 14.0D0*PG0ZN))*ZET3)/
     # (36.0D0*ZNm1*ZN*ZNp1*ZNp2) +
     # (CA*(-45.0D0 -
     # (2.0D0*(144.0D0 + ZN*(48.0D0 + ZN*(-1964.0D0+ZN*(-6382.0D0+
     # 4.0D0*GE*ZNm1*ZNp1**2*
     # ZNp2*(47.0D0 + 56.0D0*ZN) + ZN*(-8669.0D0 + ZN*(-3142.0D0 +
     # ZN*(3748.0D0+
     # ZN*(3740.0D0 + 861.0D0*ZN))))))))*nf)/(ZN3*ZNp1**3*
     # plinomZNm) - 448.0D0*nf*PG0ZN + (72.0D0*nf*PG0ZNp2)/
     # ZNp1-12.0D0*(CF*(435.0D0+99.0D0*LT+54.0D0*Pi2) + 2.0D0*nf*
     # (5.0D0*Pi2 +16.0D0
     # *ZET3))))/
     # 864.0D0 + (CA**2*GE*(-2592.0D0 + ZN*ZNp1*
     # (-1296.0D0 + 54.0D0*GE*ZNm1**3*ZNp2*
     # (8.0D0 + ZN*(24.0D0 + ZN*(23.0D0 + 2.0D0*ZN*(5.0D0 + ZN)))) +
     # ZN*(11016.0D0 + ZN*(-9640.0D0+ ZN*(-14936.0D0+ ZN*(10346.0D0+
     # ZN*(6269.0D0+
     # ZN*(-9818.0D0+ ZN*(-6304.0D0+ ZN*(4.0D0+ ZN)*(395.0D0+
     # 404.0D0*ZN))))))))) -
     # 54.0D0*ZN*ZNp2**3*ZNb*(2.0D0*(-1.0D0+ ZN2*(5.0D0 + 2.0D0*ZNm1*
     # ZN))*
     # PG0ZNm1 + ZNm1*ZN*ZNp1*
     # ((-2.0D0 + ZN*(2.0D0 - 4.0D0*ZN + GE*ZNb))*PG1ZNm1 +
     # ZN*ZNb*(2.0D0*PG0ZNp3*PG1ZNp3 -
     # 2.0D0*PG2ZNm1 + 7.0D0*ZET3)))))/
     # (108.0D0*ZN3*ZNp2**3*ZNb**3) +
     # (CF*(-2.0D0 +
     # (nf*(24.0D0+ ZN*ZNp1*(48.0D0+ ZN*(-18.0D0+ ZN*(-72.0D0+ ZN*
     # (-62.0D0+211.0D0*ZN +
     # 142.0D0*ZN2-152.0D0*ZN3-164.0D0*ZN4-41.0D0*ZN5+12.0D0*LT*
     # ZNp1**3*plinomZNm +
     # 24.0D0*ZNp1**3*plinomZNm*ZET3))))))/
     # (ZN4*ZNp1**4*plinomZNm)))/24.0D0

      RETURN
      END

      complex*16 function MellinH2qq(ZN)
      implicit none
      complex*16 ZN
      real *8 CF

      CF=4d0/3

      MellinH2qq=-(CF**2*(-(1+ZN)**(-2)+(4d0*(-3+2*ZN))/(-1+ZN)**2
     &           - (4 + 8*ZN)/ZN**2))/4d0

      return
      end




      COMPLEX*16 FUNCTION ACG1(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(1+X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
!      IMPLICIT COMPLEX*16 (A-H,O-Z)
      IMPLICIT none
      COMPLEX*16 ZN,T
      REAL*8 AK2,DL,ZERO,ONE
      INTEGER L,K
C
      DIMENSION AK2(10)
      DATA AK2/0.999999980543793D+0,
     &        -0.999995797779624D+0,
     &         0.916516447393493D+0,
     &        -0.831229921350708D+0,
     &         0.745873737923571D+0,
     &        -0.634523908078600D+0,
     &         0.467104011423750D+0,
     &        -0.261348046799178D+0,
     &         0.936814286867420D-1,
     &        -0.156249375012462D-1/

!      COMMON/ACLOG2/ AK2(10)
!      COMMON/ACCON1/ DL
!      COMMON/VAL   / ZERO,ONE
      DL = DLOG(2.0D0)
      ZERO=0.0D0
      ONE =1.0D0
C
      T=DCMPLX(ZERO,ZERO)
      DO 1 L=2,11
      K=L-1
      T=T+AK2(K)/(ZN+DBLE(K+1))
1     CONTINUE
C
      ACG1=(DL*DL- ZN*T)/2.0D0
C
      RETURN
      END
      COMPLEX*16   FUNCTION DACG1(ZN)

C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(1+X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT none
      COMPLEX*16 ZN,T1,T2
      REAL*8 AK2,DL,ZERO
      INTEGER L,K

      DIMENSION AK2(10)
      DATA AK2/0.999999980543793D+0,
     &        -0.999995797779624D+0,
     &         0.916516447393493D+0,
     &        -0.831229921350708D+0,
     &         0.745873737923571D+0,
     &        -0.634523908078600D+0,
     &         0.467104011423750D+0,
     &        -0.261348046799178D+0,
     &         0.936814286867420D-1,
     &        -0.156249375012462D-1/
      DL = DLOG(2.0D0)
      ZERO=0.0D0

      T1=DCMPLX(ZERO,ZERO)
      T2=DCMPLX(ZERO,ZERO)
      DO 1 L=2,11
      K=L-1
      T1=T1+AK2(K)/(ZN+DBLE(L))
      T2=T2+AK2(K)/(ZN+DBLE(L))/(ZN+DBLE(L))
1     CONTINUE

      DACG1=(ZN*T2 - T1)/(2.0D0)

      RETURN
      END





C --------------------------------------------------------------------------------------------------------------------------------

      COMPLEX*16   FUNCTION DDACG1(ZN)

C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(1+X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT none
      COMPLEX*16 ZN,T1,T2,ZNPKP1
      REAL*8 AK2,DL,ZERO
      INTEGER L,K

      DIMENSION AK2(10)
      DATA AK2/0.999999980543793D+0,
     &        -0.999995797779624D+0,
     &         0.916516447393493D+0,
     &        -0.831229921350708D+0,
     &         0.745873737923571D+0,
     &        -0.634523908078600D+0,
     &         0.467104011423750D+0,
     &        -0.261348046799178D+0,
     &         0.936814286867420D-1,
     &        -0.156249375012462D-1/
      DL = DLOG(2.0D0)
      ZERO=0.0D0

      T1=DCMPLX(ZERO,ZERO)
      T2=DCMPLX(ZERO,ZERO)
      DO 1 L=2,11
      K=L-1
      ZNPKP1=ZN+DBLE(L)
      T1=T1+AK2(K)/ZNPKP1/ZNPKP1
      T2=T2+AK2(K)/ZNPKP1/ZNPKP1/ZNPKP1
1     CONTINUE

      DDACG1=(T1-ZN*T2)
c      write(*,*) ZN,T1,T2,DACG1
      RETURN
      END




      COMPLEX*16 FUNCTION DACG2(ZN)
C     ----------------------------
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(X)*LOG(1+X)**2/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
!      IMPLICIT COMPLEX*16 (A-H,O-Z)
      IMPLICIT none
      COMPLEX*16 ZN,T1,T2
      REAL*8 AK3,DL,ZERO,ONE
      INTEGER K,L
C
      DIMENSION AK3(11)
      DATA AK3/9.99999989322696D-1,
     &        -1.49999722020708D+0,
     &         1.74988008499745D+0,
     &        -1.87296689068405D+0,
     &         1.91539974617231D+0,
     &        -1.85963744001295D+0,
     &         1.62987195424434D+0,
     &        -1.17982353224299D+0,
     &         6.28710122994999D-1,
     &        -2.11307487211713D-1,
     &         3.28953352932140D-2/
!      COMMON/ACLOG3/ AK3(11)
!      COMMON/ACCON1/ DL
!      COMMON/VAL   / ZERO,ONE
      DL = DLOG(2.0D0)
      ZERO=0.0D0
      ONE =1.0D0

C
      T1=DCMPLX(ZERO,ZERO)
      T2=DCMPLX(ZERO,ZERO)
      DO 1 L=3,13
      K=L-2
      T1=T1+AK3(K)/(ZN+DBLE(L))
      T2=T2+AK3(K)/(ZN+DBLE(L))/(ZN+DBLE(L))
1     CONTINUE
C
      DACG2=(ZN*T2 - T1)/3.0D0
C
      RETURN
      END


      COMPLEX*16 FUNCTION ACG4(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LI2(-X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
!      IMPLICIT COMPLEX*16 (A-H,O-Z)
      IMPLICIT none
      COMPLEX*16 ZN,ZL,ZNL,ZNL1,T,V1
      REAL*8 AK1,DL,ZERO,ONE,ZET2,PI
      INTEGER L


C
!      COMMON/ACLOG1/ AK1(9)
!      COMMON/ACCON1/ DL
!      COMMON/ACCON2/ ZET2
!      COMMON/VAL   / ZERO,ONE
      DIMENSION AK1(9)
      DATA AK1/0.999999974532240D+0,
     &        -0.499995525890027D+0,
     &         0.333203435554182D+0,
     &        -0.248529457735332D+0,
     &         0.191451164493502D+0,
     &        -0.137466222203386D+0,
     &         0.792107405737825D-1,
     &        -0.301109652783781D-1,
     &         0.538406198111749D-2/

      DL = DLOG(2.0D0)
      ZERO=0.0D0
      ONE =1.0D0
      PI = 3.141592653589793238462643D0
      ZET2 = PI**2/6.0D0


C
      T=DCMPLX(-ZET2/2.0D0*DL,ZERO)
C
      DO 1 L=1,9
      ZL =DCMPLX(DBLE(L),ZERO)
      ZNL =ZN+ZL
      ZNL1=ZNL+DCMPLX(ONE,ZERO)
      CALL BET(ZNL1,V1)
C
      T=T+AK1(L)*(ZN/ZNL*ZET2/2.0D0+ZL/ZNL**2*(DL-V1))
1     CONTINUE
C
      ACG4=T
C
      RETURN
      END

      COMPLEX*16 FUNCTION DACG4(ZN)
C     ----------------------------
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: Log(X)*LI2(-X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
!      IMPLICIT COMPLEX*16 (A-H,O-Z)
      IMPLICIT none
      COMPLEX*16 ZN,ZL,ZNL,ZNL1,T,V1,V2
      REAL*8 AK1,DL,ZERO,ONE,ZET2,PI
      INTEGER L

C
!      COMMON/ACLOG1/ AK1(9)
!      COMMON/ACCON1/ DL
!      COMMON/ACCON2/ ZET2
!      COMMON/VAL   / ZERO,ONE
C
      DIMENSION AK1(9)
      DATA AK1/0.999999974532240D+0,
     &        -0.499995525890027D+0,
     &         0.333203435554182D+0,
     &        -0.248529457735332D+0,
     &         0.191451164493502D+0,
     &        -0.137466222203386D+0,
     &         0.792107405737825D-1,
     &        -0.301109652783781D-1,
     &         0.538406198111749D-2/

      DL = DLOG(2.0D0)
      ZERO=0.0D0
      ONE =1.0D0
      PI = 3.141592653589793238462643D0
      ZET2 = PI**2/6.0D0

      T=DCMPLX(ZERO,ZERO)
C
      DO 1 L=1,9
      ZL =DCMPLX(DBLE(L),ZERO)
      ZNL =ZN+ZL
      ZNL1=ZNL+DCMPLX(ONE,ZERO)
      CALL BET(ZNL1,V1)
      CALL BET1(ZNL1,V2)
C
      T=T+AK1(L)*L*(ZET2/2.0D0-V2-2*(DL-V1)/ZNL)/ZNL/ZNL
1     CONTINUE
C
      DACG4=T
C
      RETURN
      END


      COMPLEX*16 FUNCTION ACG5(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(X) LI2(X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
!      IMPLICIT COMPLEX*16 (A-H,O-Z)
      IMPLICIT none
      COMPLEX*16 ZN,ZL,ZNL,ZL1,T,Z,PS,S1P,S1
      REAL*8 AK1,DL,ZERO,ONE,ZET2,PI,GE
      INTEGER L

C
!      COMMON/ACLOG1/ AK1(9)
!      COMMON/VAL   / ZERO,ONE
!      COMMON/ACCON1/ DL
!      COMMON/ACCON2/ ZET2
!      COMMON/ACCON4/ GE
C
      DIMENSION AK1(9)
      DATA AK1/0.999999974532240D+0,
     &        -0.499995525890027D+0,
     &         0.333203435554182D+0,
     &        -0.248529457735332D+0,
     &         0.191451164493502D+0,
     &        -0.137466222203386D+0,
     &         0.792107405737825D-1,
     &        -0.301109652783781D-1,
     &         0.538406198111749D-2/

      GE   = 0.57721566490153D+0
      DL = DLOG(2.0D0)
      ZERO=0.0D0
      ONE =1.0D0
      PI = 3.141592653589793238462643D0
      ZET2 = PI**2/6.0D0

      T=DCMPLX(ZERO,ZERO)
      Z=ZN
      DO 1 L=1,9
      ZL =DCMPLX(DBLE(L),ZERO)
      ZNL=Z+ZL
      ZL1=ZNL+ONE
      CALL PSI0(ZL1,PS)
      CALL PSI1(ZL1,S1P)
      S1=PS+GE
      T=T-AK1(L)*ZL/ZNL**2*(ZET2+S1P-2.0D0*S1/ZNL)
1     CONTINUE
C
      ACG5=T
C
      RETURN
      END
      COMPLEX*16 FUNCTION ACG6(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LI3(X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
!      IMPLICIT COMPLEX*16 (A-H,O-Z)
      IMPLICIT none
      COMPLEX*16 ZN,ZL,ZNL,ZNL1,T,V1,S1
      REAL*8 AK1,DL,ZERO,ONE,ZET3,GE,ZET2,PI
      INTEGER L

C
!      COMMON/ACLOG1/ AK1(9)
!      COMMON/ACCON1/ DL
!      COMMON/ACCON2/ ZET2
!      COMMON/ACCON3/ ZET3
!      COMMON/ACCON4/ GE
!      COMMON/VAL   / ZERO,ONE
C
      DIMENSION AK1(9)
      DATA AK1/0.999999974532240D+0,
     &        -0.499995525890027D+0,
     &         0.333203435554182D+0,
     &        -0.248529457735332D+0,
     &         0.191451164493502D+0,
     &        -0.137466222203386D+0,
     &         0.792107405737825D-1,
     &        -0.301109652783781D-1,
     &         0.538406198111749D-2/

      ZET3 = 1.20205690315959428540D0
      GE   = 0.57721566490153D+0
      DL = DLOG(2.0D0)
      ZERO=0.0D0
      ONE =1.0D0
      PI = 3.141592653589793238462643D0
      ZET2 = PI**2/6.0D0


      T=DCMPLX(DL*ZET3,ZERO)
      DO 1 L=1,9
      ZL=DCMPLX(DBLE(L),ZERO)
      ZNL=ZN+ZL
      ZNL1=ZNL+ONE
      CALL PSI0(ZNL1,V1)
      S1=V1+GE
C
      T=T-AK1(L)*(ZN/ZNL*ZET3+ZL/ZNL**2*(ZET2-S1/ZNL))
1     CONTINUE
C
      ACG6=T
C
      RETURN
      END
      COMPLEX*16 FUNCTION ACG7(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LI3(-X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
!      IMPLICIT COMPLEX*16 (A-H,O-Z)
      IMPLICIT none
      COMPLEX*16 ZN,ZL,ZNL,ZNL1,T,V1
      REAL*8 AK1,DL,ZERO,ONE,ZET2,ZET3,PI
      INTEGER L

C
!      COMMON/ACLOG1/ AK1(9)
!      COMMON/ACCON1/ DL
!      COMMON/ACCON2/ ZET2
!      COMMON/ACCON3/ ZET3
!      COMMON/VAL   / ZERO,ONE
C
      DIMENSION AK1(9)
      DATA AK1/0.999999974532240D+0,
     &        -0.499995525890027D+0,
     &         0.333203435554182D+0,
     &        -0.248529457735332D+0,
     &         0.191451164493502D+0,
     &        -0.137466222203386D+0,
     &         0.792107405737825D-1,
     &        -0.301109652783781D-1,
     &         0.538406198111749D-2/

      ZET3 = 1.20205690315959428540D0
      DL = DLOG(2.0D0)
      ZERO=0.0D0
      ONE =1.0D0
      PI = 3.141592653589793238462643D0
      ZET2 = PI**2/6.0D0


      T=DCMPLX(-3.0D0*ZET3/4.0D0*DL,ZERO)
C
      DO 1 L=1,9
      ZL =DCMPLX(DBLE(L),ZERO)
      ZNL =ZN+ZL
      ZNL1=ZNL+DCMPLX(ONE,ZERO)
      CALL BET(ZNL1,V1)
C
      T=T+AK1(L)*(ZN/ZNL*3.0D0*ZET3/4.0D0+ZL/ZNL**2/2.0D0*ZET2
     & -ZL/ZNL**3*(DL-V1))
1     CONTINUE
C
      ACG7=T
C
      RETURN
      END


      COMPLEX*16 FUNCTION ACG9(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: S12(-X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
!      IMPLICIT COMPLEX*16 (A-H,O-Z)
      IMPLICIT none
      COMPLEX*16 ZN,T,T1,T2,ZNKL,ZK,ZNK
      REAL*8 AK1,AK2,AK3,DL,ZERO,ONE,ZET3,PI
      INTEGER K,L,L1
C
!      COMMON/ACLOG1/ AK1(9)
!      COMMON/ACLOG2/ AK2(10)
!      COMMON/ACLOG3/ AK3(11)
!      COMMON/ACCON1/ DL
!      COMMON/ACCON3/ ZET3
!      COMMON/VAL   / ZERO,ONE
C
      DIMENSION AK1(9),AK2(10),AK3(11)
      DATA AK1/0.999999974532240D+0,
     &        -0.499995525890027D+0,
     &         0.333203435554182D+0,
     &        -0.248529457735332D+0,
     &         0.191451164493502D+0,
     &        -0.137466222203386D+0,
     &         0.792107405737825D-1,
     &        -0.301109652783781D-1,
     &         0.538406198111749D-2/
c
      DATA AK2/0.999999980543793D+0,
     &        -0.999995797779624D+0,
     &         0.916516447393493D+0,
     &        -0.831229921350708D+0,
     &         0.745873737923571D+0,
     &        -0.634523908078600D+0,
     &         0.467104011423750D+0,
     &        -0.261348046799178D+0,
     &         0.936814286867420D-1,
     &        -0.156249375012462D-1/
C
      DATA AK3/9.99999989322696D-1,
     &        -1.49999722020708D+0,
     &         1.74988008499745D+0,
     &        -1.87296689068405D+0,
     &         1.91539974617231D+0,
     &        -1.85963744001295D+0,
     &         1.62987195424434D+0,
     &        -1.17982353224299D+0,
     &         6.28710122994999D-1,
     &        -2.11307487211713D-1,
     &         3.28953352932140D-2/

      ZET3 = 1.20205690315959428540D0
      DL = DLOG(2.0D0)
      ZERO=0.0D0
      ONE =1.0D0
      PI = 3.141592653589793238462643D0


      T=DCMPLX(ZET3*DL/8.0D0,ZERO)
      DO 1 K=1,9
      T1=DCMPLX(ZERO,ZERO)
      DO 2 L=2,11
      L1=L-1
      ZNKL=ZN+DCMPLX(DBLE(K+L),ZERO)
      T1=T1+AK2(L1)/ZNKL
2     CONTINUE
      ZK=DCMPLX(DBLE(K),ZERO)
      ZNK=ZN+ZK
      T=T-AK1(K)*ZN/ZNK*(ZET3/8.0D0-T1/2.0D0)
1     CONTINUE
      T2=DCMPLX(ZERO,ZERO)
      DO 3 K=3,13
      L=K-2
      ZK=DCMPLX(DBLE(K),ZERO)
      ZNK=ZN+ZK
      T2=T2+AK3(L)/ZNK
3     CONTINUE
      T=T-T2/2.0D0
C
      ACG9=T
C
      RETURN
      END


      COMPLEX*16 FUNCTION ACG13(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(1+X)/(1+X)*LI2(-X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
!      IMPLICIT COMPLEX*16 (A-H,O-Z)
      IMPLICIT none
      COMPLEX*16 ZN,T0,T1,T2,V1,ZNK1
      REAL*8 AK2,AK3,DL,ZERO,ONE,ZET2,PI
      INTEGER L,K

C
!      COMMON/ACLOG2/ AK2(10)
!      COMMON/ACLOG3/ AK3(11)
!      COMMON/ACCON1/ DL
!      COMMON/ACCON2/ ZET2
!      COMMON/VAL   / ZERO,ONE
C
      DIMENSION AK2(10),AK3(11)
      DATA AK2/0.999999980543793D+0,
     &        -0.999995797779624D+0,
     &         0.916516447393493D+0,
     &        -0.831229921350708D+0,
     &         0.745873737923571D+0,
     &        -0.634523908078600D+0,
     &         0.467104011423750D+0,
     &        -0.261348046799178D+0,
     &         0.936814286867420D-1,
     &        -0.156249375012462D-1/
C
      DATA AK3/9.99999989322696D-1,
     &        -1.49999722020708D+0,
     &         1.74988008499745D+0,
     &        -1.87296689068405D+0,
     &         1.91539974617231D+0,
     &        -1.85963744001295D+0,
     &         1.62987195424434D+0,
     &        -1.17982353224299D+0,
     &         6.28710122994999D-1,
     &        -2.11307487211713D-1,
     &         3.28953352932140D-2/


      DL = DLOG(2.0D0)
      ZERO=0.0D0
      ONE =1.0D0
      PI = 3.141592653589793238462643D0
      ZET2 = PI**2/6.0D0

      T0=DCMPLX(-1.0D0/4.0D0*ZET2*DL**2,ZERO)
C
      T1=DCMPLX(ZERO,ZERO)
      DO 1 L=3,13
      K=L-2
      T1=T1+AK3(K)/(ZN+DBLE(L))
1     CONTINUE
C
      T2=DCMPLX(ZERO,ZERO)
      DO 2 L=2,11
      K=L-1
      ZNK1=ZN+DBLE(L+1)
      CALL BET(ZNK1,V1)
      T2=T2+AK2(K)*ZN/(ZN+DBLE(L))*(ZET2/2.0D0-(DL-V1)/(ZN+DBLE(L)))
2     CONTINUE
C
      ACG13=T0+(T1+T2)/2.0D0
C
      RETURN
      END

      COMPLEX*16 FUNCTION ACG20(ZZ)
C     -----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: (LI3(X)-ZETA3)/(X-1)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
!      IMPLICIT COMPLEX*16 (A-H,O-Z)
      IMPLICIT none
      COMPLEX*16 ZN,S1,S2,ZN1,ZZ,ZNK,ZNK1,R,T,T1,T2,PS,PS1
      REAL*8 DL,ZERO,ONE,ZETA2,ZET2,ZETA3,ZET3,GE,P24,P34,P22,CK2,CK4,PI
      INTEGER L,K

C
!      COMMON/ACLOG7/ CK2(13)
!      COMMON/ACLOG9/ CK4(13)
!      COMMON/POLY2 / P22(4)
!      COMMON/POLY4 / P24(5),P34(5)
!      COMMON/ACCON1/ DL
!      COMMON/ACCON2/ ZET2
!      COMMON/ACCON3/ ZET3
!      COMMON/ACCON4/ GE
!      COMMON/VAL   / ZERO,ONE
C
      DIMENSION CK2(13),CK4(13),P22(4),P24(5),P34(5)
      DL = DLOG(2.0D0)
      ZERO=0.0D0
      ONE =1.0D0
      ZET3 = 1.20205690315959428540D0
      GE   = 0.57721566490153D+0
      PI = 3.141592653589793238462643D0
      ZET2 = PI**2/6.0D0
      ZETA2=ZET2
      ZETA3=ZET3

      DATA CK2/ 1.9891961855478257D-8,!0.480322239287449D+0,
     &         -6.050451121009104D-6,!-0.168480825099580D+1,
     &          2.1253046159349207,!0.209270571620726D+1,
     &         -1.0523034804772378D0,!-0.101728150275998D+1,
     &          0.160179976133047D+0,
     &         -0.351982236677917D+0,
     &          0.141033316846244D+1,
     &         -0.353343997843391D+1,
     &          0.593934696819832D+1,
     &         -0.660019784804042D+1,
     &          0.466330349413063D+1,
     &         -0.189825467489058D+1,
     &          0.339772909487512D+0/
C
      DATA CK4/-1.844613928270178D-8,!0.192962504274437D+0,
     &          2.2150086978693064D0,!0.000005641557253D+0,
     &         -0.9133677154535804D0,!-0.196891075399448D+1,
     &          3.4783104357500143D0,!0.392919138747074D+1,
     &         -2.823955592989266D0,!-0.290306105685546D+1,
     &          0.992890266001707D+0,
     &         -0.130026190226546D+1,
     &          0.341870577921103D+1,
     &         -0.576763902370864D+1,
     &          0.645554138192407D+1,
     &         -0.459405622046138D+1,
     &          0.188510809558304D+1,
     &         -0.340476080290674D+0/



C----------------------------------------------------------------------

!      P12(1)=ZETA3-11.0D0/6.0*ZETA2+4.0D0/3.0D0!=-0.48032221939548725
!      P12(2)=3.0D0*ZETA2-13.0D0/4.0D0!=1.684802200544679
!      P12(3)=-3.0D0/2.0D0*ZETA2+5.0D0/2.0D0!=0.0325988997276605
!      P12(4)=1.0D0/3.0D0*ZETA2-7.0D0/12.0D0!=-0.03502197771725779
C-------------------------------
      P22(1)=-1.0D0
      P22(2)=5.0D0/2.0D0
      P22(3)=-2.0D0
      P22(4)=1.0D0/2.0D0

C-------------------------------
!      P14(1)=257.D0/144.0D0-205.0D0/72.0D0*ZET2+ZET2**2
!      P14(2)=-167.0D0/36.0D0+25.0D0/6.0D0*ZET2
!      P14(3)=101.0D0/24.0D0-23.0D0/12.0D0*ZET2
!      P14(4)=-59.0D0/36.0D0+13.0D0/18.0D0*ZET2
!      P14(5)=41.0D0/144.0D0-ZET2/8.0D0
C-------------------------------
      P24(1)=-167.0D0/36.0D0+25.0D0/6.0D0*ZET2
      P24(2)=235.0D0/18.0D0-8.0D0*ZET2
      P24(3)=-40.0D0/3.0D0+6.0D0*ZET2
      P24(4)=109.0D0/18.0D0-8.0D0/3.0D0*ZET2
      P24(5)=-41.0D0/36.0D0+ZET2/2.0D0
C-------------------------------
      P34(1)=35.0D0/12.0D0
      P34(2)=-26.0D0/3.0D0
      P34(3)=19.0D0/2.0D0
      P34(4)=-14.0D0/3.0D0
      P34(5)=11.0D0/12.0D0
C----------------------------------------------------------------------
C
C >>> ACCOUNT FOR POLYNOM PARTS
C
!      CK2(1)=CK2(1)+P12(1)
!      CK2(2)=CK2(2)+P12(2)
!      CK2(3)=CK2(3)+P12(3)
!      CK2(4)=CK2(4)+P12(4)
C
!      write(*,*) CK2(1),CK2(2),CK2(3),CK2(4),CK2(5)
!      CK4(1)=CK4(1)+P14(1)
!      CK4(2)=CK4(2)+P14(2)
!      CK4(3)=CK4(3)+P14(3)
!      CK4(4)=CK4(4)+P14(4)
!      CK4(5)=CK4(5)+P14(5)
!      write(*,*) CK4(1),CK4(2),CK4(3),CK4(4),CK4(5),CK4(6)
C
      T=DCMPLX(ZET2**2/2.0D0,ZERO)
C
      ZN = ZZ
      ZN1=ZN+ONE
      CALL PSI0(ZN1,PS)
      S1=PS+GE
      T=T-ZET3*S1
C
      DO 1 K=1,13
      L=K-1
      ZNK=ZN+DBLE(L)
      ZNK1=ZNK+ONE
      CALL PSI0(ZNK1,PS)
      S1=PS+GE
      R=ZN/ZNK
C
      T=T+CK2(K)*R*S1
1     CONTINUE
C
      DO 3 K=1,4
      L=K-1
      ZNK=ZN+DBLE(L)
      ZNK1=ZNK+ONE
      CALL PSI0(ZNK1,PS)
      CALL PSI1(ZNK1,PS1)
      S1=PS+GE
      S2=-PS1+ZET2
      R=ZN/ZNK
      T1=S1**2+S2
C
      T=T-P22(K)*R*T1
3     CONTINUE
C
      DO 4 K=1,13
      L=K-1
      ZNK=ZN+DBLE(L)
      ZNK1=ZNK+ONE
      CALL PSI0(ZNK1,PS)
      S1=PS+GE
C
      T=T-CK4(K)/ZNK*ZN/2.0D0
4     CONTINUE
C
      DO 5 K=1,5
      L=K-1
      ZNK=ZN+DBLE(L)
      ZNK1=ZNK+ONE
      CALL PSI0(ZNK1,PS)
      CALL PSI1(ZNK1,PS1)
      S1=PS+GE
      S2=-PS1+ZET2
      T1=-S1/ZNK
      T2=(S1**2+S2)/ZNK
C
      T=T-(P24(K)*T1+P34(K)*T2)*ZN/2.0D0
5     CONTINUE
C
C
      ACG20=T

C
      RETURN
      END


      COMPLEX*16 FUNCTION ACG22(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(X) LI2(X)/(X-1)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
!      IMPLICIT COMPLEX*16 (A-H,O-Z)
      IMPLICIT none
      COMPLEX*16 ZN,ZNK1,PS,PS1,PS2,T,TA
      REAL*8 CK1,P21,GE,ZERO,ONE,ZET2,PI
      INTEGER L
C
!      COMMON/ACLOG6/ CK1(12)
!      COMMON/POLY1 / P21(4)
!      COMMON/ACCON4/ GE
!      COMMON/VAL   / ZERO,ONE
C

      DIMENSION CK1(12),P21(4)
      PI = 3.141592653589793238462643D0
      ZET2 = PI**2/6.0D0
      DATA CK1/ 2.201218318731435D-8,!-0.283822933724932D+0,
     &          2.833327652357064D0,!0.999994319023731D+0,
     &         -1.8330909624101532D0,!-0.124975762907682D+1,
     &          0.7181879191200942,!0.607076808008983D+0,
     &         -0.280403220046588D-1,
     &         -0.181869786537805D+0,
     &          0.532318519269331D+0,
     &         -0.107281686995035D+1,
     &          0.138194913357518D+1,
     &         -0.111100841298484D+1,
     &          0.506649587198046D+0,
     &         -0.100672390783659D+0/
C----------------------------------------------------------------------
!      P11(1)=-49.0D0/36.0+ZET2
!      P11(2)=11.0D0/6.0D0
!      P11(3)=-7.0D0/12.0D0
!      P11(4)=1.0D0/9.0D0
C---------------------------------------------------------
      P21(1)=11.0D0/6.0
      P21(2)=-3.0D0
      P21(3)=3.0D0/2.0D0
      P21(4)=-1.0D0/3.0D0
C
!      CK1(1)=CK1(1)+P11(1)
!      CK1(2)=CK1(2)+P11(2)
!      CK1(3)=CK1(3)+P11(3)
!      CK1(4)=CK1(4)+P11(4)
!      write(*,*) CK1(1),CK1(2),CK1(3),CK1(4),CK1(5)

      ZERO=0.0D0
      ONE =1.0D0
      GE   = 0.57721566490153D+0

      T=DCMPLX(ZERO,ZERO)
C
      DO 1 L=1,12
      ZNK1=ZN+DBLE(L)
      CALL PSI1(ZNK1,PS1)
!      write(*,*)T,CK1(L)
      T=T+CK1(L)*PS1
1     CONTINUE
C
      DO 2 L=1,4
      ZNK1=ZN+DBLE(L)
      CALL PSI0(ZNK1,PS)
      CALL PSI1(ZNK1,PS1)
      CALL PSI2(ZNK1,PS2)
      TA=(PS+GE)*PS1-PS2/2.0D0
      T=T-P21(L)*TA
2     CONTINUE
C
      ACG22=T
C
      RETURN
      END

!------------------------------------------------------------------------


      SUBROUTINE GAMMAL(ZZ,RES)
C     -------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  LOG(GAMMA(Z)) FOR COMPLEX ARGUMENT
C
C************************************************************************
C
!      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT none
      COMPLEX*16 ZZ,Z,RES,T,T1,T2,X,X2,ONE
      REAL*8 PI,R
C
      Z=ZZ
      PI = 3.141592653589793238462643D0
C
      ONE=DCMPLX(1.0D0,0.0D0)
      T=DCMPLX(0.0D0,0.0D0)
2     R=SQRT(DREAL(Z)**2+DIMAG(Z)**2)
      IF(R.GT.10.0D0) GOTO 1
      T=T-LOG(Z)
      Z=Z+ONE
      GOTO 2
1     CONTINUE
C
      T1=Z*(LOG(Z)-1.0D0)+LOG(2.0D0*PI/Z)/2.0D0
C
      X=ONE/Z
      X2=X*X
      T2 = (1.D0/12.D0+(-1.D0/360.D0+(1.D0/1260.D0+(-1.D0/1680.D0+(1.D0/
     #1188.D0+(-691.D0/360360.D0+(1.D0/156.D0-3617.D0/122400.D0*X2)*X2
     #)*X2)*X2)*X2)*X2)*X2)*X
C
      RES=T1+T2+T
C
      RETURN
      END
      SUBROUTINE GAMMA(ZZ,RES)
C     ------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  GAMMA(Z) FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT none
      COMPLEX*16 ZZ,T,RES
C
      CALL GAMMAL(ZZ,T)
C
      RES=EXP(T)
C
      RETURN
      END
      SUBROUTINE BETA(AA,BB,RES)
C     --------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C---  BETA(A,B) FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT none
      COMPLEX*16 T1,T2,T3,T,RES,AA,BB
C
      CALL GAMMAL(AA,T1)
      CALL GAMMAL(BB,T2)
      CALL GAMMAL(AA+BB,T3)
      T=T1+T2-T3
C
      RES=EXP(T)
C
      RETURN
      END
      SUBROUTINE PSI0(ZZ,RES)
C     -----------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  PSI(Z) FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT none
!      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 Z,T,ONE,Y,T0,RES,ZZ,Y2
      REAL*8 R
C
      ONE=DCMPLX(1.0D0,0.0D0)
C
      Z=ZZ
      T=DCMPLX(0.0D0,0.0D0)
2     R=SQRT(DREAL(Z)**2+DIMAG(Z)**2)
      IF(R.GT.10.0D0) GOTO 1
      T=T-ONE/Z
      Z=Z+ONE
      GOTO 2
1     Y=ONE/Z
      Y2=Y*Y
!      write(*,*) Y2
      T0 = (-1.D0/2.D0+(-1.D0/12.D0+(1.D0/120.D0+(-1.D0/252.D0+(1.D0/240
     #.D0+(-1.D0/132.D0+(691.D0/32760.D0+(-ONE/12.0D0+ONE*3617.0D0
     #             /8160.D0*Y2
     #  )*Y2  )*Y2  )*Y2  )*Y2  )*Y2  )*Y2  )*Y)*Y-LOG(Y)
C
      RES=T+T0
C
      RETURN
      END
      SUBROUTINE PSI1(ZZ,RES)
C     -----------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  PSI'(Z) FOR COMPLEX ARGUMENT
C
C************************************************************************
C      IMPLICIT none
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 Z,T,ONE,Y,T0,RES,ZZ
C
      ONE=DCMPLX(1.0D0,0.0D0)
C
      Z=ZZ
      T=DCMPLX(0.0D0,0.0D0)
2     R=SQRT(DREAL(Z)**2+DIMAG(Z)**2)
      IF(R.GT.10.0D0) GOTO 1
      T=T+ONE/Z**2
      Z=Z+ONE
      GOTO 2
1     Y=ONE/Z
      Y2=Y*Y
      T0 = (1.D0+(1.D0/2.D0+(1.D0/6.D0+(-1.D0/30.D0+
     &(1.D0/42.D0+(-1.D0/30.D0+(5.D0/66.D0-691.D0/2730.D0*Y2)
     &*Y2)*Y2)*Y2)*Y2)*Y)*Y)*Y
C
      RES=T+T0
C
      RETURN
      END
      SUBROUTINE PSI2(ZZ,RES)
C     -----------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  PSI''(Z) FOR COMPLEX ARGUMENT
C
C************************************************************************
C      IMPLICIT none
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 Z,T,ONE,Y,T0,RES,ZZ
C
      ONE=DCMPLX(1.0D0,0.0D0)
      TWO=ONE*2.0D0
C
      Z=ZZ
      T=DCMPLX(0.0D0,0.0D0)
2     R=SQRT(DREAL(Z)**2+DIMAG(Z)**2)
      IF(R.GT.10.0D0) GOTO 1
      T=T-TWO/Z**3
      Z=Z+ONE
      GOTO 2
1     Y=ONE/Z
      Y2=Y*Y
      T0 =(-1.D0+(-1.D0+(-1.D0/2.D0+(1.D0/6.D0+(-1.D0/6.D0+(3.D0/
     &10.D0+(-5.D0/6.D0+691.D0/210.D0*Y2)*Y2)*Y2)*Y2)*Y2)*Y)*Y)*Y2
C
      RES=T+T0
C
      RETURN
      END
      SUBROUTINE PSI3(ZZ,RES)
C     -----------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  PSI'''(Z) FOR COMPLEX ARGUMENT
C
C************************************************************************
C      IMPLICIT none
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 Z,T,ONE,Y,T0,RES,ZZ
C
      ONE=DCMPLX(1.0D0,0.0D0)
      SIX=ONE*6.0D0
C
      Z=ZZ
      T=DCMPLX(0.0D0,0.0D0)
2     R=SQRT(DREAL(Z)**2+DIMAG(Z)**2)
      IF(R.GT.10.0D0) GOTO 1
      T=T+SIX/Z**4
      Z=Z+ONE
      GOTO 2
1     Y=ONE/Z
      Y2=Y*Y
C
      T0 = (2.D0+(3.D0+(2.D0+(-1.D0+(4.D0/3.D0+(-3.D0+(10.D0+(-691.D0/15
     #.D0+(280.D0-10851.D0/5.D0*Y2  )*Y2  )*Y2  )*Y2  )*Y2  )*Y2  )*Y2
     #)*Y)*Y)*Y2*Y
C
      RES=T+T0
C
      RETURN
      END
      SUBROUTINE BET(ZZ,RES)
C     ----------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  BET(Z) FOR COMPLEX ARGUMENT  $\beta(z)$
C
C************************************************************************
C      IMPLICIT none
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 Z,ONE,RES,V1,V2,Z1,Z2,ZZ
C
      ONE=DCMPLX(1.0D0,0.0D0)
C
      Z=ZZ
      Z1=(Z+ONE)/2.0D0
      Z2=Z/2.0D0
      CALL PSI0(Z1,V1)
      CALL PSI0(Z2,V2)
C
      RES=(V1-V2)/2.0D0
C
      RETURN
      END
      SUBROUTINE BET1(ZZ,RES)
C     -----------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  BET'(Z) FOR COMPLEX ARGUMENT  $\beta'(z)$
C
C************************************************************************
C      IMPLICIT none
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 Z,ONE,RES,V1,V2,Z1,Z2,ZZ
C
      ONE=DCMPLX(1.0D0,0.0D0)
C
      Z=ZZ
      Z1=(Z+ONE)/2.0D0
      Z2=Z/2.0D0
      CALL PSI1(Z1,V1)
      CALL PSI1(Z2,V2)
C
      RES=(V1-V2)/4.0D0
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FLI4(X)
C     -----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  LI4(X) FOR -1. LE . X . LE .+1
C
C************************************************************************
C      IMPLICIT none
      IMPLICIT REAL*8(A-H,O-Z)
C
      A=1D0
      F=0D0
      AN=0D0
      TCH=1D-16
1     AN=AN+1D0
      A=A*X
      B=A/AN**4
      F=F+B
      IF(ABS(B)-TCH)2,2,1
2     FLI4 = F
      END


