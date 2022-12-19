
subroutine fokkerplanck1D(alfa2, h, n, dt, nt, xend, d1, d2, d3, vj, fj0, out_fj, dfj0)
    implicit none
    real*8, intent(in)  :: alfa2     
    integer, intent(in) :: n, nt
    real*8, intent(in) :: h, dt, xend
    real*8, intent(in) ::  d1(:), d2(:), d3(:), vj(:)
    real*8, intent(inout) :: fj0(:)
    real*8, intent(inout) :: out_fj(:)
    real*8, intent (inout), optional :: dfj0(:)

    integer :: i0
    real*8, parameter :: zero=0.d0
    real*8  y(n),x(n+2)
    real*8,dimension(:),allocatable:: fj, dfj,  givi
    integer i, ii, it, ibeg, klo, khi, ierr, klo1, khi1
    real*8 shift, ybeg, yend, tend, dff
    real*8 fout1,fout2
    interface 
    subroutine savelyev_solver(alfa2, nt, h, dt, n, ybeg, yend, d1,d2,d3, y)
        implicit none
        real*8, intent(in)  :: alfa2 
        integer, intent(in) :: nt, n
        real*8, intent(in)  :: h, dt
        real*8, intent(in)  :: ybeg, yend        
        real*8, intent(in)  :: d1(n+1),d2(n+1),d3(n+1)
        real*8, intent(inout) :: y(n)
    end subroutine
    subroutine burying_procedure(f0, dv, ibeg, df0)
        ! процедура закапывания
        implicit none
        real*8,  intent(inout)  :: f0(:)
        real*8,  intent(in)     :: dv
        real*8,  intent (inout), optional :: df0(:)
        integer, intent(out) :: ibeg
    end subroutine            
    end interface

    i0 = size(vj)

    !!!!!! grid !!!!!!!!!
    !!  shift=h*0.1d0 !0.01d0
    do i=1,n+2
        x(i)=h*dble(i-1) !+shift
    end do

    do i=1,n
        call lock(vj,i0,x(i+1),klo,khi,ierr)
        if(ierr.eq.1) then
            write(*,*)'lock error #1 in finction fokkerplanck'
            write(*,*)'j=', 123, ' v=', x(i+1)
            write(*,*)'vj(1)=',vj(1),' vj(i0)=',vj(i0)
            pause
            stop
        end if
        call linf(vj,fj0,x(i+1),y(i),klo,khi)
    end do

    ybeg=fj0(1)  !boundary conditions
    yend=zero
    !!!!!!!!!!!!   solve problem   !!!!!!!!!!!!!!!!!!!!!!!!!!
    call savelyev_solver(alfa2, nt, h, dt, n, ybeg, yend, d1,d2,d3, y)

    allocate(fj(n+2))
    fj(1)=ybeg
    fj(n+2)=yend
    do i=1,n
        fj(i+1)=y(i)
    end do
    out_fj(:) = fj(:)
    !call write_distribution(fj, n, 0.1999d0)
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



    if (present(dfj0)) then
        call burying_procedure(fj0, vj(2), ibeg, dfj0)
    else 
        call burying_procedure(fj0, vj(2), ibeg)
    end if

    allocate(fj(i0), dfj(i0))
    fj(:)=fj0(:)
!   нормировка распределения
    if(ibeg.gt.0) then
        call integral(ibeg,i0,vj,fj,fout1)
        do i=1,i0
            fj(i)=fj0(i)
            if (present(dfj0)) then
                dfj(i)=dfj0(i+1)
            end if
        end do
        call integral(ibeg,i0,vj,fj,fout2)
        do i=ibeg,i0
            fj0(i)=fj(i)*fout1/fout2
            if (present(dfj0)) then
                dfj0(i)=dfj(i)*fout1/fout2
            end if            
        end do
!!      write(*,*)'#1 j,k,ibeg=',j,k,ibeg
!!      write(*,*)'#1 v(ibeg)=',vj(ibeg),' f1/f2=',fout1/fout2
    end if
    deallocate(fj,dfj)

end subroutine fokkerplanck1D


subroutine burying_procedure(f0, dv, ibeg, df0)
    ! процедура закапывания
    implicit none
    real*8,  intent(inout)  :: f0(:)
    real*8,  intent(in)     :: dv
    real*8,  intent (inout), optional :: df0(:)
    integer, intent(out) :: ibeg

    integer i, ii,  i0
    real*8, allocatable  :: f(:), df(:)

    i0 = size(f0)
    allocate(f(i0), df(i0))
    
    f(:)=f0(:)
    df(:)=0d0

    do i=2, i0-1
        df(i)=0.5d0*(f(i+1)-f(i-1))/dv
    end do
    df(1)=0d0
    df(i0)=(f(i0)-f(i0-1))/dv

!   сдвиг расределения вправо. зачем-то ???
    ii=0
    ibeg=0
    do i=i0-1,1,-1
        if(df(i).gt.0d0) then
!          write(*,*) '#1 positive derivs'
!          write(*,*) '#1 df>0: i,j,k=',i,j,k
!          write(*,*) '#1 dfj(i),i,j,k=',dfj(i),i,j,k
!          write(*,*)
            f0(i)=f0(i+1)
            if (present(df0)) then
                df0(i)=df0(i+1)
            end if
            ii=i
        end if
        if(f0(i).lt.f0(i+1)) then 
            f0(i)=f0(i+1)
            if (present(df0)) then
                df0(i)=df0(i+1)
            end if            
            ii=i
        end if
    end do
    deallocate(f,df)
    ibeg=ii
end subroutine