subroutine fokkerplanck(dtstep,time,nomer)
    implicit none
    real*8 t,dtstep,dtau,d5,check1,check2,check3
    integer nr,ispl,ip,im,ii,ibeg,nomer,vivod,iunit
    common /a0ab/ nr
    real*8 ynzm0,pm0,plaunp,plaunm,fmaxw,time,pachka
    common/grillspektr/ ynzm0(1001),pm0(1001),ispl, plaunp,plaunm,ip,im
    integer i0
    parameter(i0=1002)
    real*8 vij,fij0,fij,dfij,dij,enorm,fst
    common/lh/vij(i0,100),fij0(i0,100,2),fij(i0,100,2),dfij(i0,100,2), dij(i0,100,2),enorm(100),fst(100)
    integer n,i,j,it,nt,k
    real*8 xend,h,shift,ybeg,yend,tend,dt,dff
    real*8,dimension(:),allocatable:: y,x,xx,xxm,xxp,a,b,c,f
    real*8,dimension(:),allocatable:: vj,fj,dfj,d1,d2,d3,givi
    real*8 znak,alfa2,zero,dt0,h0,eps,r,fvt,fout1,fout2
    common/ef/ alfa2
    real*8 calls
    common/firstcall/calls
    real*8 d0
    integer jindex,kindex,klo,khi,ierr,klo1,khi1
    integer klo2,klo3,khi2,khi3,ierr1,ierr2,ierr3
    common/dddql/ d0,jindex,kindex
    parameter(zero=0.d0,dt0=0.1d0,h0=0.1d0,eps=1.d-7)
!
    do k=1,2
     kindex=k
     znak=2.d0*dble(k)-3.d0
!
!!       do j=1,nr
!!        fij0(i,j,k)=fmaxw(vij(i,j),znak*enorm(j),dff)
!!       end do
!!       if(1.gt.0) go to 1
!
     d0=zero
     do j=1,nr
      jindex=j
      dtau=dtstep*fst(j)
      nt=1
      if(dtau.gt.dt0) then
       nt=1+dtau/dt0
      end if
      dt=dtau/nt
      r=dble(j)/dble(nr+1)
      xend=3.d10/fvt(r)
!!        xend=vij(i0,j)
      n=xend/h0-1
      h=xend/dble(n+1)
      if(h.gt.h0) then
       n=n+1
       h=xend/dble(n+1)
      end if
      allocate(y(n),x(n+2),xx(n+1),a(n),b(n),c(n),f(n))
!!!!!! grid !!!!!!!!!
!!       shift=h*0.1d0 !0.01d0
      do i=1,n+2
       x(i)=h*dble(i-1) !+shift
      end do
      do i=1,n+1
       xx(i)=h/2.d0+h*dble(i-1) !+shift
      end do
      allocate(vj(i0),fj(i0))
      do i=1,i0
       vj(i)=vij(i,j)
       fj(i)=fij0(i,j,k)
      end do
      do i=1,n
       call lock(vj,i0,x(i+1),klo,khi,ierr)
       if(ierr.eq.1) then
        write(*,*)'lock error #1 in finction fokkerplanck'
        write(*,*)'j=',j,' v=',x(i+1)
        write(*,*)'vj(1)=',vj(1),' vj(i0)=',vj(i0)
        pause
        stop
       end if
       call linf(vj,fj,x(i+1),y(i),klo,khi)
      end do
      deallocate(fj)
      ybeg=fij0(1,j,k)  !boundary conditions
      yend=zero
      alfa2=znak*enorm(j)
!!!!!!!!!!!!!EVALUATING DIFFUSION!!!!!!!!!!!!!!!!!!

  allocate(d1(n+1),d2(n+1),d3(n+1))
  do i=1,n+1
  d1(i)=0d0
  d2(i)=0d0
  d3(i)=0d0
  end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do it=1,nt
       call abccoef(a,b,c,f,y,dt,n,ybeg,yend,xx,h,d1,d2,d3)
       call tridag(a,b,c,f,y,n)
!!         t=dt*dble(it)
      end do
      deallocate(d1,d2,d3)
  pachka=0.35
  if(time*nomer.gt.pachka)then
  open(iunit,file='lhcd/pressF.dat',position="append")
  do vivod=1,n
!!	write(iunit,*)time,vivod,f(vivod)
  end do
  end if
      allocate(fj(n+2))
      fj(1)=ybeg
      fj(n+2)=yend
      do i=1,n
       fj(i+1)=y(i)
      end do
      do i=2,i0-1
       if(vij(i,j).lt.xend) then
        call lock(x,n+2,vij(i,j),klo,khi,ierr)
        if(ierr.eq.1) then
         write(*,*)'lock error #2 in finction fokkerplanck'
         write(*,*)'j=',j,' vij=',vij(i,j)
         write(*,*)'x(1)=',x(1),' x(n+2)=',x(n+2)
         pause
         stop
         end if
        call linf(x,fj,vij(i,j),fij0(i,j,k),klo,khi)
       else
        fij0(i,j,k)=zero
       end if
      end do
      deallocate(fj)
      deallocate(y,x,xx,a,b,c,f)
      allocate(fj(i0),dfj(i0))
      do i=1,i0
       fj(i)=fij0(i,j,k)
       dfj(i)=zero
      end do
      do i=1,i0
       if(i.eq.1) then
        dfj(i)=zero
       else if(i.eq.i0) then
        dfj(i)=(fj(i)-fj(i-1))/vij(2,j)
       else
        dfj(i)=0.5d0*(fj(i+1)-fj(i-1))/vij(2,j)
       end if
      end do
      ii=0
      ibeg=0
      do i=i0-1,1,-1
       if(dfj(i).gt.zero) then
!          write(*,*) '#1 positive derivs'
!          write(*,*) '#1 df>0: i,j,k=',i,j,k
!          write(*,*) '#1 dfj(i),i,j,k=',dfj(i),i,j,k
!          write(*,*)
        fij0(i,j,k)=fij0(i+1,j,k)
        ii=i
       end if
       if(fij0(i,j,k).lt.fij0(i+1,j,k)) then 
        fij0(i,j,k)=fij0(i+1,j,k)
        ii=i
       end if
      end do
      ibeg=ii
!
      if(ibeg.gt.0) then
       call integral(ibeg,i0,vj,fj,fout1)
       do i=1,i0
        fj(i)=fij0(i,j,k)
       end do
       call integral(ibeg,i0,vj,fj,fout2)
       do i=ibeg,i0
        fij0(i,j,k)=fj(i)*fout1/fout2
       end do
!!         write(*,*)'#1 j,k,ibeg=',j,k,ibeg
!!         write(*,*)'#1 v(ibeg)=',vj(ibeg),' f1/f2=',fout1/fout2
      end if
!
      deallocate(vj,fj,dfj)

     end do
!!!!!!!!!!!!!!!!!
!
1      continue
!
     d0=1.d0
     do j=1,nr
      jindex=j
      dtau=dtstep*fst(j)
      nt=1
      if(dtau.gt.dt0) then
       nt=1+dtau/dt0
      end if
      dt=dtau/nt
      r=dble(j)/dble(nr+1)
      xend=3.d10/fvt(r)
!!        xend=vij(i0,j)
      n=xend/h0-1
      h=xend/dble(n+1)
      if(h.gt.h0) then
       n=n+1
       h=xend/dble(n+1)
      end if
      allocate(y(n),x(n+2),xx(n+1),a(n),b(n),c(n),f(n))
  allocate(xxm(n+1),xxp(n+1))
!!	write(*,*)'OSHIBKA2, j=',j
!!!!!! grid !!!!!!!!!
!!       shift=h*0.1d0 !0.01d0
      do i=1,n+2
       x(i)=h*dble(i-1) !+shift
      end do
      do i=1,n+1
       xx(i)=h/2.d0+h*dble(i-1) !+shift
!	 xxm(i)=xx(i)-h/2d0
!	 xxp(i)=xx(i)+h/2d0
      end do
  if (nomer.gt.9) then 
  open(iunit,file='lhcd/distribution/xxx.dat',position="append")
!	check1=xx(4)-h/2d0-xx(4)+h/2d0
!	check2=xx(5)+h/2d0-xxp(5)
!	check3=xx(5)-h/2d0-xxm(5)
!	do i=1,n
!	write(iunit,*)xx(1)-h/2d0-xxx(1),xx(1)+h/2d0-xxx(2)
!	write(iunit,*)xx(2)-h/2d0-xxx(2),xx(2)+h/2d0-xxx(4)
!	write(iunit,*)xx(2)-xxx(3),xx(1)-xxx(1)
!	write(iunit,*)check1,check2,check3
!	end do
  end if
  close(iunit)
      allocate(vj(i0),fj(i0))
      do i=1,i0
       vj(i)=vij(i,j)
       fj(i)=fij(i,j,k)
      end do
      do i=1,n
       call lock(vj,i0,x(i+1),klo,khi,ierr)
       if(ierr.eq.1) then
        write(*,*)'lock error #3 in finction fokkerplanck'
        write(*,*)'j=',j,' v=',x(i+1)
        write(*,*)'vj(1)=',vj(1),' vj(i0)=',vj(i0)
        pause
        stop
       end if
       call linf(vj,fj,x(i+1),y(i),klo,khi)
      end do
      deallocate(fj)
      ybeg=fij(1,j,k)  !boundary conditions
      yend=zero
      alfa2=znak*enorm(j)
!!!!!!!!!!!!!EVALUATING DIFFUSION!!!!!!!!!!!!!!!!!!
!!	write(*,*)'OSHIBKA, j=',j
  if (nomer.gt.9) then 
  open(iunit,file='lhcd/distribution/dddd.dat',position="append")
  end if
      allocate(d1(n+1),d2(n+1),d3(n+1))	
  do i=1,n+1
  call lock(vj,i0,xx(i),klo1,khi1,ierr1)
  call lock(vj,i0,xx(i)-h/2d0,klo2,khi2,ierr2)
  call lock(vj,i0,xx(i)+h/2d0,klo3,khi3,ierr3)
        if(ierr1.eq.1) then
  write(*,*)'lock error in finction d2(x)'
  write(*,*)'j=',j,' v=',xx(i)
  write(*,*)'klo1=',klo1,'khi1=',khi1,'i=',i
  write(*,*)'vj(1)=',vj(1),' vj(i0)=',vj(i0)
  pause
  stop
        end if
        if(ierr2.eq.1) then
  write(*,*)'lock error in finction d2(x)'
  write(*,*)'j=',j,' v=',xxm(i)
  write(*,*)'klo2=',klo2,'khi2=',khi2,'i=',i
  write(*,*)'vj(1)=',vj(1),' vj(i0)=',vj(i0)
  pause
  stop
        end if
        if(ierr3.eq.1) then
  write(*,*)'lock error in finction d2(x)'
  write(*,*)'j=',j,' v=',xxp(i)
  write(*,*)'klo3=',klo3,'khi3=',khi3,'i=',i
  write(*,*)'vj(1)=',vj(1),' vj(i0)=',vj(i0)
  pause
  stop
        end if
  d1(i)=dij(klo1,j,k)
  d2(i)=dij(klo2,j,k)
  d3(i)=dij(klo3,j,k)	
  if (nomer.gt.9) then 
!	write(iunit,*)d1(i),d2(i),d3(i)
  end if
  end do
  if (nomer.gt.9) then 
  close(iunit)
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!solve problem
      do it=1,nt
       call abccoef(a,b,c,f,y,dt,n,ybeg,yend,xx,h,d1,d2,d3)
       call tridag(a,b,c,f,y,n)
!!         t=dt*dble(it)
      end do
      deallocate(d1,d2,d3)
      allocate(fj(n+2))
      fj(1)=ybeg
      fj(n+2)=yend
      do i=1,n
       fj(i+1)=y(i)
      end do
  if (nomer.gt.9) then 
          !print *, "fk time=", time
          call write_distribution(fj, n, time)
!	open(iunit,file='lhcd/distribution/distribution.dat',position="append")
!	      do i=1,n
!      	      write(iunit,*)i,fj(i+1)
!	      end do
  end if
!	close(iunit)
      do i=2,i0-1
       if(vij(i,j).lt.xend) then
        call lock(x,n+2,vij(i,j),klo,khi,ierr)
        if(ierr.eq.1) then
         write(*,*)'lock error #4 in finction fokkerplanck'
         write(*,*)'j=',j,' vij=',vij(i,j)
         write(*,*)'x(1)=',x(1),' x(n+2)=',x(n+2)
         pause
         stop
        end if
        call linf(x,fj,vij(i,j),fij(i,j,k),klo,khi)
       else
        fij(i,j,k)=zero
       end if
      end do
      deallocate(fj)

      deallocate(y,x,xx,xxm,xxp,a,b,c,f)
      allocate(fj(i0),dfj(i0))
      do i=1,i0
       fj(i)=fij(i,j,k)
       dfj(i)=zero
      end do
      do i=1,i0
       if(i.eq.1) then
        dfj(i)=zero
       else if(i.eq.i0) then
        dfj(i)=(fj(i)-fj(i-1))/vij(2,j)
       else
        dfj(i)=0.5d0*(fj(i+1)-fj(i-1))/vij(2,j)
       end if
      end do
      ii=0
      ibeg=0
      do i=i0-1,1,-1
       if(dfj(i).gt.zero) then
!          write(*,*) '#2 positive derivs'
!c          write(*,*) '#2 df>0: i,j,k=',i,j,k
!c          write(*,*) '#2 dfj(i),i,j,k=',dfj(i),i,j,k
!c          write(*,*) '#2 fij=',fij(i,j,k)
!          write(*,*)
        fij(i,j,k)=fij(i+1,j,k)
        dfij(i,j,k)=dfij(i+1,j,k)
        ii=i
       end if
       if(fij(i,j,k).lt.fij(i+1,j,k)) then
        fij(i,j,k)=fij(i+1,j,k)
        dfij(i,j,k)=dfij(i+1,j,k)
        ii=i
       end if
      end do
      ibeg=ii
!
      if(ibeg.gt.0) then
       call integral(ibeg,i0,vj,fj,fout1)
       do i=1,i0
        fj(i)=fij(i,j,k)
        dfj(i)=dfij(i,j,k)
       end do
       call integral(ibeg,i0,vj,fj,fout2)
       do i=ibeg,i0
        fij(i,j,k)=fj(i)*fout1/fout2
        dfij(i,j,k)=dfj(i)*fout1/fout2
       end do
!!         write(*,*)'#2 j,k,ibeg=',j,k,ibeg
!!         write(*,*)'#2 v(ibeg)=',vj(ibeg),' f1/f2=',fout1/fout2
      end if
!
      deallocate(vj,fj,dfj)

     end do
    end do
!
    end