
subroutine fokkerplanck1D(h, n, nt, d1, d2, d3, vj, fj0)
    implicit none
    integer, intent(in) :: n, nt
    real*8, intent(in) :: h
    real*8, dimension(:), intent(in) :: vj, d1, d2, d3
    real*8, dimension(:), intent(inout) :: fj0
    integer, parameter :: i0=1002
    real*8,dimension(:),allocatable:: y,x,xx,xxm,xxp,a,b,c,f
    real*8,dimension(:),allocatable:: fj, dfj,  givi


    allocate(y(n),x(n+2),xx(n+1),a(n),b(n),c(n),f(n))
    !!!!!! grid !!!!!!!!!
    !!  shift=h*0.1d0 !0.01d0
    do i=1,n+2
        x(i)=h*dble(i-1) !+shift
    end do

    do i=1,n+1
        xx(i)=h/2.d0+h*dble(i-1) !+shift
    end do

    allocate(fj(i0))
    do i=1,i0
        fj(i)=fj0(i)
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
    ybeg=fj0(1)  !boundary conditions
    yend=zero
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do it=1,nt
        call abccoef(a,b,c,f,y,dt,n,ybeg,yend,xx,h,d1,d2,d3)
        call tridag(a,b,c,f,y,n)
    !!         t=dt*dble(it)
    end do
    
    allocate(fj(n+2))
    fj(1)=ybeg
    fj(n+2)=yend
    do i=1,n
        fj(i+1)=y(i)
    end do
    do i=2,i0-1
        if(vj(i).lt.xend) then
            call lock(x,n+2,vj(i),klo,khi,ierr)
            if(ierr.eq.1) then
                write(*,*)'lock error #2 in finction fokkerplanck'
                write(*,*)'j=', 123 ,' vij=',vj(i)
                write(*,*)'x(1)=',x(1),' x(n+2)=',x(n+2)
                pause
                stop   
            end if
            call linf(x,fj,vj(i),fj0(i),klo,khi)
        else
            fj0(i)=zero
        end if
    end do
    deallocate(fj)
    deallocate(y,x,xx,a,b,c,f)

    allocate(fj(i0),dfj(i0))
    do i=1,i0
        fj(i)=fj0(i)
        dfj(i)=zero
    end do
    do i=1,i0
        if(i.eq.1) then
            dfj(i)=zero
        else if(i.eq.i0) then
                dfj(i)=(fj(i)-fj(i-1))/vj(2)
            else
                dfj(i)=0.5d0*(fj(i+1)-fj(i-1))/vj(2)
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
            fj0(i)=fj0(i+1)
            ii=i
        end if
        if(fj0(i).lt.fj0(i+1)) then 
            fj0(i)=fj0(i+1)
            ii=i
        end if
    end do
    ibeg=ii
!
    if(ibeg.gt.0) then
        call integral(ibeg,i0,vj,fj,fout1)
        do i=1,i0
            fj(i)=fj0(i)
        end do
        call integral(ibeg,i0,vj,fj,fout2)
        do i=ibeg,i0
            fj0(i)=fj(i)*fout1/fout2
        end do
!!      write(*,*)'#1 j,k,ibeg=',j,k,ibeg
!!      write(*,*)'#1 v(ibeg)=',vj(ibeg),' f1/f2=',fout1/fout2
    end if
!
    deallocate(vj,fj,dfj)

end subroutine fokkerplanck1D
