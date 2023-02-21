      subroutine saveprofiles
       use constants , only :pi4, pme, pqe, c0, zero
       use spline
       use chebyshev
       use approximation       
       use plasma
       use rt_parameters
       use maxwell
      implicit none
      integer i,k,j,klo,khi,ierr
      include 'for/parameter.inc'
      include 'for/const.inc'
      include 'for/status.inc'
      integer im,ip
      real*8 p_in
       ! real*8 anz,apz,share
      !real*8 pchm0
      !real*8 fpol,dfmy
      !real*8 xlogj
      !real*8 znak_tor,znak_pol
      !common/left/ znak_tor,znak_pol
      !real*8 ynzmp(1001),pmp(1001),ynzmm(1001),pmm(1001)

      !real*8 efld(100),r,vmax
      !real*8 zff,fnr,fnrr
      !real*8 pn,gst,dens,tmp,vt,vclt
      !real*8 znak
      integer :: calls = 0
      !save share

      p_in=dble(QLH)    ! input LH power, MW

      if(calls.eq.0) then
            call read_parameters('lhcd/ray_tracing.dat')
      end if

      call init_plasma(NA1,ABC,BTOR,RTOR,UPDWN,GP2,
     & AMETR,RHO,SHIF,ELON,TRIA,MU,NE,TE,TI,ZEF,UPL)

      call calc_enorm
      
 !     open(96, file='lhcd/out/difsave.dat', position='append')
      if(calls.eq.0) then
            call init_maxwell
            calls=1
      end if
!      close(96)
      end

