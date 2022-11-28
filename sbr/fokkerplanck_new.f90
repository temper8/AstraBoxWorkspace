subroutine fokkerplanck_new(dtstep,time,nomer)
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
        do j=1,nr
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

            allocate(d1(n+1),d2(n+1),d3(n+1))

            d1(:)=0d0
            d2(:)=0d0
            d3(:)=0d0

            call fokkerplanck1D(h, n, nt, d1, d2, d3, vij(:,j), fij0(:,j,k))

            call init_diffusion(h, n, vij(:,j), dij(:,j,k), d1, d2, d3)

            call fokkerplanck1D(h, n, nt, d1, d2, d3, vij(:,j), fij0(:,j,k))

            deallocate(d1,d2,d3)

        end do
    end do

 end

 subroutine init_diffusion(h, n, vj, dj, d1, d2, d3)
    implicit none
    integer, intent(in) :: n
    real*8, intent(in) :: h
    real*8, dimension(:), intent(in) :: vj, dj
    real*8, dimension(:), intent(out) :: d1, d2, d3
    integer, parameter :: i0=1002

    do i=1,n+1
        xx(i)=h/2.d0+h*dble(i-1) !+shift
    end do

    do i=1,n+1
        call lock(vj, i0, xx(i), klo1, khi1, ierr1)
        call lock(vj, i0, xx(i)-h/2d0, klo2, khi2, ierr2)
        call lock(vj, i0, xx(i)+h/2d0, klo3, khi3, ierr3)
        if(ierr1.eq.1) then
            write(*,*)'lock error in finction d2(x)'
            write(*,*)'j=', 123,' v=', xx(i)
            write(*,*)'klo1=', klo1, 'khi1=', khi1, 'i=', i
            write(*,*)'vj(1)=', vj(1),' vj(i0)=', vj(i0)
            pause
            stop
        end if
        if(ierr2.eq.1) then
            write(*,*)'lock error in finction d2(x)'
            write(*,*)'j=', 123, ' v=', xxm(i)
            write(*,*)'klo2=', klo2, 'khi2=', khi2, 'i=',i
            write(*,*)'vj(1)=', vj(1), ' vj(i0)=', vj(i0)
            pause
            stop
        end if
        if(ierr3.eq.1) then
            write(*,*)'lock error in finction d2(x)'
            write(*,*)'j=', 123, ' v=', xxp(i)
            write(*,*)'klo3=', klo3, 'khi3=', khi3, 'i=',i
            write(*,*)'vj(1)=', vj(1), ' vj(i0)=', vj(i0)
            pause
            stop
        end if
        d1(i) = dj(klo1)
        d2(i) = dj(klo2)
        d3(i) = dj(klo3)	
    end do

end subroutine