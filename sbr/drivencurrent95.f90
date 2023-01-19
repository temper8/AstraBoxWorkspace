      subroutine drivencurrent95(outj, sigmaj, UPL, NRD, NA1, TIME, TAU, ROC, RTOR, GP2)
! ******************************************************************
!   outj(i)  = LH driven current density, MA/m^2
!   dndt(i)  = d^2Jr1/dt^2/E, MA/m^2/sec^2/(V/m), ~runaway d(el.density)/dt/E
!   djdt(i)  = dJr2/dt, time drivative of runaway current Jr2, MA/m^2/sec
!   outjrun(i)  = LH driven runaway current density, MA/m^2
!   outnerun(i) = runaway electron density/10^19 m^-3
! ******************************************************************
      implicit none
 ! variables imported from ASTRA
      integer NRD
      ! NRD     501 - Maximum size of the radial grid
      integer NA1
      ! NA1     Edge grid point number: ROC=RHO(NA1)      
      double precision TIME, TAU, RTOR, ROC, GP2
      double precision UPL(NRD)
!
      real*8 outj(NRD),sigmaj(NRD),afld(NRD),dtau
      integer i,inpt,ispectr
      real*8,dimension(:),allocatable:: outjp,outjm,ohjp,ohjm
      real*8 dt, cup,cup0,cum,cum0,cp,cm,cp0,cm0,aiint
      real*8, parameter :: zero=0.d0, eps=1.d-2 
!
      inpt=NA1
      allocate(outjp(inpt),outjm(inpt),ohjp(inpt),ohjm(inpt))
      do i=1,inpt
          afld(i)=UPL(i)/RTOR/GP2 !!variant
      end do
!
!!!!!!!!!!!!! starting LH current calculation !!!!!!!!!!!!!!!!!
      outj=zero
      outjp=zero
      outjm=zero
      ohjp=zero
      ohjm=zero
      cup=zero
      cum=zero
      cp=zero
      cm=zero
      cup0=zero
      cum0=zero
      cp0=zero
      cm0=zero
!
!!positive spectrum:
      ispectr=1
      call lhcurrent(outjp,ohjp,cup,cup0,inpt,ispectr)
      if(cup0.ne.zero) then
          cp0=aiint(ohjp,roc)
          if(cp0.ne.zero) then
          do i=1,inpt
              ohjp(i)=cup0*ohjp(i)/cp0
          end do
          end if
      end if
      if(cup.ne.zero) then
          cp=aiint(outjp,roc)
          if(cp.ne.zero) then
              do i=1,inpt
                  outjp(i)=cup*outjp(i)/cp
              end do
          end if
      end if

!!negative spectrum:
      ispectr=-1
      call lhcurrent(outjm,ohjm,cum,cum0,inpt,ispectr)
      if(cum0.ne.zero) then
          cm0=aiint(ohjm,roc)
          if(cm0.ne.zero) then
            do i=1,inpt
                ohjm(i)=cum0*ohjm(i)/cm0
          end do
          end if
      end if
      if(cum.ne.zero) then
          cm=aiint(outjm,roc)
          if(cm.ne.zero) then
              do i=1,inpt
                  outjm(i)=cum*outjm(i)/cm
              end do
          end if
      end if

      do i=1,inpt
          outj(i)=outjp(i)+outjm(i)
          sigmaj(i)=zero
          if(dabs(afld(i)).gt.eps) then
              sigmaj(i)=(ohjp(i)+ohjm(i))/afld(i)
          end if
!!!!       write(*,*) i,outj(i)
      end do
!

      write(*,*)'time=',time
      write(*,*)'cup=',cup,' cp=',cp
      write(*,*)'cum=',cum,' cm=',cm
      write(*,*)'cup0=',cup0,' cp0=',cp0
      write(*,*)'cum0=',cum0,' cm0=',cm0
      write(*,*)'sigma driven current, MA=',cp0+cm0
      write(*,*)'driven current, MA=',cup+cum
      write(*,*)
!
      call fokkerplanck_new(time, TAU)
!
      deallocate(outjp,outjm,ohjp,ohjm)
!
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine lhcurrent(outj,ohj,cuj,cujoh,inpt,ispectr)
!!      implicit real*8 (a-h,o-z)
      use plasma
      use maxwell
      use rt_parameters, only : nr
      implicit none
      real*8 outj(*),ohj(*),cuj,cujoh,curs,curs0,curdir
      real*8 currn,pqe,vt0,fvt,ccur,cfull,cfull0
      real*8 r,pn,fn1,fn2,fnr,fnrr,vt,vto!,rh1
      integer klo,khi,ierr,nrr,i,j,inpt,ispectr,ismthout
      !common /a0ab/ nr
      real*8 y2dn,y2tm,y2tmi
      common /a0l3/ y2dn(501),y2tm(501),y2tmi(501)
      integer inew
      common /cnew/ inew !est !sav2008
      real*8 zv1,zv2,sk,fout
      common/plosh/ zv1(100,2),zv2(100,2),sk(100)
      integer k
      !parameter(i0=1002)
      !real*8 vij,fij0,fij,dfij,dij,enorm,fst
      !common/lh/ vij(i0,100), fij0(i0,100,2), fij(i0,100,2), dfij(i0,100,2), dij(i0,100,2), enorm(100), fst(100)
      real*8,dimension(:),allocatable:: vj, fj, fj0, cur, cur0, currnt, rxx, wrk
      parameter(ismthout=1)
      interface 
      function currlhcd(v,f) result(curs)
            implicit none
            real*8 v(:),f(:)
            real*8 curs
      end function
      end interface

      allocate(vj(i0),fj(i0),fj0(i0),cur(nr),cur0(nr),currnt(nr+2),rxx(nr+2),wrk(nr+2))
!---------------------------------------------------
! initial constants
!---------------------------------------------------
      pqe=4.803e-10
      vt0=fvt(zero)
      ccur=pqe*vt0*0.333d-9
      curdir=-dble(ispectr)
!
      cfull=zero
      cfull0=zero
      k=(3-ispectr)/2
      do j=1,nr
          do i=1,i0
              vj(i)=vij(i,j) !Vpar/Vt
              fj0(i)=fij0(i,j,k)
              fj(i)=fij(i,j,k)-fij0(i,j,k)
          end do
          r=dble(j)/dble(nr+1)
          if(inew.eq.0) then !vardens
              pn=fn1(r,fnr)
          else
              pn=fn2(r,fnr,fnrr)
          end if
          vt=fvt(r)
          vto=vt/vt0
          curs  = currlhcd(vj,fj)
          cur(j)=curs*pn*ccur*curdir*vto  !Ampere/cm2
          cfull=cfull+cur(j)*sk(j)
          curs0 = currlhcd(vj,fj0)          
          cur0(j)=curs0*pn*ccur*curdir*vto  !Ampere/cm2
          cfull0=cfull0+cur0(j)*sk(j)
      end do
      cuj=cfull*1d-6   !driven current, MA
      cujoh=cfull0*1d-6   !driven current, MA
!!      write(*,*)
!!      write(*,*)'ccur',ccur,' curdir=',curdir,' nr=',nr
!!      write(*,*)'cu_out, MA=',cu_out,' cfull, A=',cfull
!!           close(111)
!      pause
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      currn=cur(1)                   ! Jstoped, A/cm^2
      currnt(1)=currn*1.d-2          ! Jstoped, MA/m^2
      rxx(1)=zero
      do j=1,nr
          rxx(j+1)=dble(j)/dble(nr+1)
          currn=cur(j)                   ! Jstopped, A/cm^2
          currnt(j+1)=currn*1.d-2        ! Jstoped, MA/m^2
      end do
      nrr=nr+2
      rxx(nrr)=1.d0
      currnt(nr+2)=zero
!
      if(ismthout.ne.0) then
          do i=1,nrr
              wrk(i)=currnt(i)
          end do
          call fsmoth4(rxx,wrk,nrr,currnt)
      end if
!
      rh(1)=rh1
      if(rh(inpt).gt.1d0) rh(inpt)=1.d0
      do j=1,inpt
          call lock2(rxx,nrr,rh(j),klo,khi,ierr)
          if(ierr.ne.0) then
               write(*,*)'lock2 error in current profile for ASTRA'
               write(*,*)'ierr=',ierr,' j=',j,' rh(j)=',rh(j)
               write(*,*)'rxx(1)=',rxx(1),' rxx(nrr)=',rxx(nrr)
               pause
          end if
          call linf(rxx,currnt,rh(j),fout,klo,khi)
          outj(j)=fout
      end do
      rh(1)=zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      currn=cur0(1)                   ! Jstoped, A/cm^2
      currnt(1)=currn*1.d-2          ! Jstoped, MA/m^2
      rxx(1)=zero
      do j=1,nr
          rxx(j+1)=dble(j)/dble(nr+1)
          currn=cur0(j)                   ! Jstopped, A/cm^2
          currnt(j+1)=currn*1.d-2        ! Jstoped, MA/m^2
      end do
      nrr=nr+2
      rxx(nrr)=1.d0
      currnt(nr+2)=zero
!
      if(ismthout.ne.0) then
          do i=1,nrr
              wrk(i)=currnt(i)
          end do
          call fsmoth4(rxx,wrk,nrr,currnt)
      end if
!
      rh(1)=rh1
      if(rh(inpt).gt.1d0) rh(inpt)=1.d0
      do j=1,inpt
          call lock2(rxx,nrr,rh(j),klo,khi,ierr)
          if(ierr.ne.0) then
              write(*,*)'#2 lock2 error in current profile for ASTRA'
              write(*,*)'ierr=',ierr,' j=',j,' rh(j)=',rh(j)
              write(*,*)'rxx(1)=',rxx(1),' rxx(nrr)=',rxx(nrr)
              pause
          end if
          call linf(rxx,currnt,rh(j),fout,klo,khi)
          ohj(j)=fout
      end do
      rh(1)=zero
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      deallocate(vj,fj,fj0,cur,cur0,currnt,rxx,wrk)
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function currlhcd(v,f) result(curs)
      implicit none
      real*8 v(:),f(:)
      real*8 curs      
      integer i0,k

      real*8 vl,vr,fl,fr
      curs=0.d0
      i0 = size(v)
      do k=1,i0-1
            vl=v(k)
            vr=v(k+1)
            fl=f(k)
            fr=f(k+1)
            curs=curs+(fl*vl+fr*vr)/2d0*(vr-vl)
      end do
end

subroutine currlhcd_old(i0,v,f,f0,curs,curs0)
      implicit none
      integer i0,k
      real*8 v(*),f(*),f0(*),curs,curs0
      real*8 vl,vr,fl,fr,zero
      parameter(zero=0.d0)
      curs=zero
      curs0=zero
      do k=1,i0-1
            vl=v(k)
            vr=v(k+1)
            fl=f(k)
            fr=f(k+1)
            curs=curs+(fl*vl+fr*vr)/2d0*(vr-vl)
            fl=f0(k)
            fr=f0(k+1)
            curs0=curs0+(fl*vl+fr*vr)/2d0*(vr-vl)
      end do
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ddc(diffusion)
      use maxwell
      implicit none
      common/testf/ tcur

      integer ntau,tc,koltoch,i,j,k,klo
      real*8 curtime,tcur
      real*16 tau0,spacing,curtime0
      parameter(tau0=3.000990745207882E-002)
      !common/lh/dij(1002,100,2)
      real*8 b,b1,b2,d,diffusion
!      real*8,dimension(:),allocatable:: diffusion
      integer i1,iunit6
      b1=0
      b2=60
      j=10
      k=1
!      allocate (diffusion(500))
! write(*,*)'time=',tcur
! do tc=1,koltoch
!       spacing=0.008/koltoch
!       curtime=tau0+spacing*tc!-0.0002
!      curtime0=curtime+0.000000001
!      if((tcur-curtime)*(tcur-curtime0).lt.zero) then
!       if((tcur-0.0301)*(tcur-0.0302).lt.zero) then
!       open(iunit6,file='lhcd/distribution/ddc.dat',position="append")
!      do i1=1,500
!      b=b1+(b2/500)*(i1-1)
!      diffusion(i1)=d(b)
!       do i=1,1001
!        write(iunit6,*) i, diffusion
!       end do
!      write(iunit6,*)
!       close(iunit6)
!       end if
! end do
! !     deallocate(diffusion)
      end














      real*8 function d(x)
      use maxwell
      implicit none
      !integer i0
      !parameter(i0=1002)
      !real*8 vij,fij0,fij,dfij,dij,enorm,fst
      !common/lh/vij(i0,100),fij0(i0,100,2),fij(i0,100,2),dfij(i0,100,2),dij(i0,100,2),enorm(100),fst(100)
      real*8,dimension(:),allocatable:: vvj,ddj
      integer klo,khi,ierr
      real*8 d0,x
      integer jindex,kindex,k,j,i
      common/dddql/ d0,jindex,kindex
      d=zero
      if(d0.eq.zero) return
      j=jindex
      if(x.ge.vij(i0,j)) return
      k=kindex
      allocate(vvj(i0),ddj(i0))
!      write(*,*)'function d(x): k=',k,' j=',j
      do i=1,i0
       vvj(i)=vij(i,j)
       ddj(i)=dij(i,j,k)
      end do
      call lock(vvj,i0,x,klo,khi,ierr)
      if(ierr.eq.1) then
       write(*,*)'lock error in finction d2(x)'
       write(*,*)'j=',j,' v=',x
       write(*,*)'vj(1)=',vvj(1),' vj(i0)=',vvj(i0)
       pause
       stop
      end if
      d=ddj(klo)
!!	call ddc(d)
!!!      call linf(vvj,ddj,x,d,klo,khi)
!
!      write(*,*)'klo=',klo,' khi=',khi
!      write(*,*)'vj(klo)=',vvj(klo),' vj(khi)=',vvj(khi)
!      write(*,*)'x=',x,' d=',d
!
      deallocate(vvj,ddj)
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      subroutine qromb(func,a,b,ss)
      implicit none
      integer jmax,jmaxp,k,km
      double precision a,b,func,ss,eps
      external func
      parameter (eps=1.d-6, jmax=200, jmaxp=jmax+1, k=5, km=k-1)
!u    uses polint,trapzd
      integer j
      double precision dss,h(jmaxp),s(jmaxp)
      h(1)=1.d0
      do 11 j=1,jmax
          call trapzd(func,a,b,s(j),j)
          if (j.ge.k) then
              call polint(h(j-km),s(j-km),k,0.d0,ss,dss)
              if (abs(dss).le.eps*abs(ss)) return
          endif
          s(j+1)=s(j)
          h(j+1)=0.25d0*h(j)
11    continue
      pause 'too many steps in qromb'
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine polint(xa,ya,n,x,y,dy)
      implicit none
      integer n,nmax
      double precision dy,x,y,xa(n),ya(n)
      parameter (nmax=10)
      integer i,m,ns
      double precision den,dif,dift,ho,hp,w,c(nmax),d(nmax)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
          dift=abs(x-xa(i))
          if (dift.lt.dif) then
              ns=i
              dif=dift
          endif
          c(i)=ya(i)
          d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
          do 12 i=1,n-m
              ho=xa(i)-x
              hp=xa(i+m)-x
              w=c(i+1)-d(i)
              den=ho-hp
              if(den.eq.0.d0)pause 'failure in polint'
              den=w/den
              d(i)=hp*den
              c(i)=ho*den
12        continue
          if (2*ns.lt.n-m)then
              dy=c(ns+1)
          else
              dy=d(ns)
              ns=ns-1
          endif
          y=y+dy
13    continue
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine trapzd(func,a,b,s,n)
      implicit none
      integer n
      double precision a,b,s,func
      external func
      integer it,j
      double precision del,sum,tnm,x
      if (n.eq.1) then
          s=0.5d0*(b-a)*(func(a)+func(b))
      else
          it=2**(n-2)
          tnm=it
          del=(b-a)/tnm
          x=a+0.5d0*del
          sum=0.d0
          do 11 j=1,it
             sum=sum+func(x)
              x=x+del
11        continue
          s=0.5d0*(s+(b-a)*sum/tnm)
      endif
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

