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
      integer :: calls = 0

      if(calls.eq.0) then
            call read_parameters('lhcd/ray_tracing.dat')
      end if
      print *, 'saveprofiles'
      print *, AMETR(1), AMETR(2), AMETR(3)
      call init_plasma(NA1,ABC,BTOR,RTOR,UPDWN,GP2,
     & AMETR,RHO,SHIF,ELON,TRIA,MU,NE,TE,TI,ZEF,UPL)

      call calc_enorm
     
      if(calls.eq.0) then
            call init_maxwell
            calls=1
            pause 'after init maxwell'
      end if

      end

