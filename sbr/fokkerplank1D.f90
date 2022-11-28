
subroutine fokkerplanck1D(vj, fj, dtstep)
    implicit none
    integer, parameter :: i0=1002
    real*8,dimension(:),allocatable:: y,x,xx,xxm,xxp,a,b,c,f
    real*8,dimension(:),allocatable:: vj,fj,dfj,d1,d2,d3,givi

    
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
!!      write(*,*)'#1 j,k,ibeg=',j,k,ibeg
!!      write(*,*)'#1 v(ibeg)=',vj(ibeg),' f1/f2=',fout1/fout2
    end if
!
    deallocate(vj,fj,dfj)

end subroutine fokkerplanck1D
