      subroutine saveprofiles
       use constants , only :pi4, pme, pqe, c0, zero
       use spline
       use chebyshev
       use approximation       
       use plasma
       use rt_parameters
       use maxwell
      implicit none
      !real*8 fn
      !external fn
      integer i,k,j,klo,khi,ierr
      include 'for/parameter.inc'
      include 'for/const.inc'
      include 'for/status.inc'
      !integer ncheb
      integer im,ip
      real*8 p_in
  
      real*8 anz,apz,share

      real*8 pchm0
      real*8 fpol,dfmy

      !real*8 chebne,chebdne,chebddne

      real*8 xlogj

      !common /a0l3/ y2dn(501),y2tm(501),y2tmi(501)
      !common /a0l5/ y2zeff(501)

      real*8 znak_tor,znak_pol
      common/left/ znak_tor,znak_pol
      real*8 ynzmp(1001),pmp(1001),ynzmm(1001),pmm(1001)
      !common/ne_cheb/chebne(50),chebdne(50),chebddne(50),ncheb
      real*8 efld(100),r,vmax,fvt
      real*8 zff,zefff,fnr,fnrr
      real*8 pn,fn1,fn2,gst,dens,tmp,ft,vt,vclt
      real*8 znak
!      common/lh/vij(i0,100),fij0(i0,100,2),fij(i0,100,2),dfij(i0,100,2)
!     &,dij(i0,100,2),enorm(100),fst(100)
      real*8 calls
      common/firstcall/calls
      save share

      p_in=dble(QLH)    ! input LH power, MW

      if(calls.eq.zero) then
            call read_parameters('lhcd/ray_tracing.dat')
      end if

      call init_plasma(NA1,ABC,BTOR,RTOR,UPDWN,GP2,
     & AMETR,RHO,SHIF,ELON,TRIA,MU,NE,TE,TI,ZEF,UPL)

      do j=1,nr
       r=dble(j)/dble(nr+1)
       call lock(rh,nspl,r,klo,khi,ierr)
       if(ierr.eq.1) then
        write(*,*)'lock error in saveprofiles, Efield'
        write(*,*)'j=',j,' rh(j)=',rh(j),' r=',r
        pause
        stop
       end if
       call linf(rh,afld,r,efld(j),klo,khi)
       if(inew.eq.0) then !vardens
        pn=fn1(r,fnr)
       else
        pn=fn2(r,fnr,fnrr)
       end if
       vt=fvt(r)
       tmp=ft(r)/0.16d-8  !Te,  KeV
       dens=pn/1.d+13     !10^13 cm^-3
       xlogj=dlog(5.1527d7*tmp*16.d0*dsqrt(tmp)/dsqrt(dens))
       enorm(j)=(3.835d0/xlogj)*efld(j)*tmp/dens
       enorm(j)=enorm(j)*5.d0/(5.d0+zefff(r))
!!       fst(j)=pn*xlogj*c0**4/pi4/vt**3
       fst(j)=((5.d0+zefff(r))/5.d0)*pn*xlogj*c0**4/pi4/vt**3
       end do
 !     open(96, file='lhcd/out/difsave.dat', position='append')
      if(calls.eq.zero) then
       do j=1,nr
         r=dble(j)/dble(nr+1)
         vclt=3.d10/fvt(r)
         print *, vclt
         call init_vi(vclt, vij(:,j))
         call init_fmaxw_classic(vclt,enorm(j),fij(:,j,1),dfij(:,j,1))
         call init_fmaxw_ext(vclt,enorm(j),fij(:,j,2),dfij(:,j,2))     
       end do
       fij0(:,:,:)=fij(:,:,:)
       dij(:,:,:)=zero
      end if
      calls=1.d0
!      close(96)
      end

