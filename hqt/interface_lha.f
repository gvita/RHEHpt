*****************
* LHAPDF version*
*****************

      subroutine lambda5qcd(lam,iset)
! Initializzation for native PDFs, dummy for LHAPFD
      implicit none
      integer iset
      real*8 lam
      end


      subroutine inizializzazPDF(pdfnamelen)
      integer flagalphas,mem,Nmem,pdfnamelen
      character*64 pdfname
      double precision alphasPDF,QMZ
      complex*16 alphas
      character(LEN = pdfnamelen) tryname
      common/flagalphas/flagalphas
      common/pdfname/pdfname
      common/mem/mem

      if(flagalphas.eq.0) then
		tryname = pdfname(1:pdfnamelen)
        call InitPDFsetByName(tryname)
        call numberPDF(Nmem)
        if(mem.gt.Nmem) then
          write(*,*)'input mem out of range!!!'
          stop
        endif

        call InitPDF(mem)
        flagalphas=1
      endif
      end


      function alphas(nq2)
c.....N.B.:Here you have in output alpha_s/(4*Pi)!!!!!
      implicit none
      double precision alphasPDF,xmu,pi
      complex*16 alphas,nq2
      real*8 f(-6:6)
      data pi/3.14159265358979323846264338d0/

      xmu=dreal(cdsqrt(nq2))
      alphas=DCMPLX(alphasPDF(xmu),0d0)
      alphas=alphas/4d0/Pi

      return
      end


      SUBROUTINE partons(xmu,X,FX,NF,ISET,IH)
! Unused variables NF,ISET
      implicit none
      integer Iprtn,ih,Irt,NF,ISET
      double precision fx(-NF:NF),x,xmu,sqrtxmu,fPDF(-6:6),temp
!      double precision UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU
c---  ih1=+1 proton
c---  ih1=-1 pbar

C---set to zero if x out of range
      if (x .ge. 1d0) then
          do Iprtn=-5,5
             fx(Iprtn)=0d0
          enddo
          return
      endif
      sqrtxmu=sqrt(xmu)
      call evolvePDF(x,sqrtxmu,fPDF)
      if (ih.eq.1) then
        do Iprtn=-5,5
          fx(Iprtn)=fPDF(Iprtn)/x
        enddo
      elseif(ih.eq.-1) then
        do Iprtn=-5,5
          fx(+Iprtn)=fPDF(-Iprtn)/x
        enddo
      endif

! With our convention the quarks u and d are exchanged
      temp=fx(1)
      fx(1)=fx(2)
      fx(2)=temp
      temp=fx(-1)
      fx(-1)=fx(-2)
      fx(-2)=temp

      return
      end
