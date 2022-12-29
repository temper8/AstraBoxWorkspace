!! calculation of distribution functions at time t1=t+dtau !!
subroutine fokkerplanck_new(time, TAU)
    use FokkerPlanck1D_mod
    use Utils
    implicit none

    type(FokkerPlanck1D) fp_test

    real*8, intent(in) :: time, TAU
    real*8 t, dtstep, dtau
    integer nr
    common /a0ab/ nr
    integer, parameter :: ntau = 10
    integer i0
    parameter(i0=1002)
    real*8 vij,fij0,fij,dfij,dij,enorm,fst
    common/lh/vij(i0,100),fij0(i0,100,2),fij(i0,100,2),dfij(i0,100,2), dij(i0,100,2),enorm(100),fst(100)
    integer n,i,j,it,nt,k
    real*8 xend,h,dt
    real*8 znak,alfa2,zero,dt0,h0,eps,r,fvt
    !common/ef/ alfa2
    
    real*8 d0
    integer jindex,kindex
    common/dddql/ d0,jindex,kindex
    parameter(zero=0.d0,dt0=0.1d0,h0=0.1d0,eps=1.d-7)
    real time1, time2

   
    interface 
    subroutine fokkerplanck1D_iter(alfa2, h, n, dt, nt, xend, d1, d2, d3, vj, fj0, out_fj, dfj0)
        implicit none
        real*8, intent(in)  :: alfa2           
        integer, intent(in) :: n, nt
        real*8, intent(in) :: h, dt, xend
        real*8, intent(in) ::  d1(:), d2(:), d3(:), vj(:)
        real*8, intent(inout) :: fj0(:)
        real*8, intent(inout) :: out_fj(:)        
        real*8, intent (inout), optional :: dfj0(:)
    end subroutine fokkerplanck1D_iter      
    
    subroutine init_diffusion(h, n, vj, dj, d1, d2, d3)
        implicit none
        integer, intent(in) :: n
        real*8, intent(in) :: h
        real*8, dimension(:), intent(in) :: vj, dj
        real*8, dimension(:), intent(out) :: d1, d2, d3
    end subroutine init_diffusion

    subroutine write_distribution(arr,N,time)
        implicit none
        real*8, intent(in) :: arr(*)
        integer, intent(in) :: N
        real*8, intent(in) :: time
    end subroutine write_distribution    

    subroutine write_matrix(arr,time, array_name)
        implicit none
        real*8, intent(in) :: arr(:,:)
        real*8, intent(in) :: time
        character(len=*), intent(in) :: array_name
    end subroutine write_matrix        

    subroutine write_v_array(v, a, time, array_name)
        implicit none
        real*8, intent(in) :: v(:,:)    
        real*8, intent(in) :: a(:,:,:)
        real*8, intent(in) :: time
        character(len=*), intent(in) :: array_name
    end subroutine    
    end interface 

    dtstep=TAU/dble(ntau) !seconds 

    print *, 'fokkerplanck_new'
    write(*,*)'time=',time,' dt=',dtstep

    time1 = sys_time()

    do j=1, nr
        jindex=j  ! common/dddql/ 
        dtau=dtstep*fst(j)
        nt=1
        if(dtau.gt.dt0) then
            nt=1+dtau/dt0
        end if
        dt=dtau/nt
        r=dble(j)/dble(nr+1)
        xend=3.d10/fvt(r)
        do k=1,2
            znak=2.d0*dble(k)-3.d0
            fp_test = FokkerPlanck1D(znak*enorm(j), xend, vij(:,j), fij0(:,j,k))
            call fp_test%init_zero_diffusion
            do i=1, ntau
                call fp_test%solve_time_step(dt, nt)
                !call fokkerplanck1D_iter(alfa2, h, n, dt, nt, xend, d1, d2, d3, vij(:,j), fij0(:,j,k), out_fj)
            end do
            fij0(:,j,k) = fp_test%f

            fp_test = FokkerPlanck1D(znak*enorm(j), xend, vij(:,j), fij(:,j,k))
            call fp_test%init_diffusion(dij(:,j,k))
            do i=1, ntau
                call fp_test%solve_time_step(dt, nt)
                !call fokkerplanck1D_iter(alfa2, h, n, dt, nt, xend, d1, d2, d3, vij(:,j), fij(:,j,k),out_fj, dfij(:,j,k))
            end do
            fij(:,j,k) = fp_test%f
        end do
   
        call write_distribution(fij0(:,j,k), i0, time)
        !call write_distribution(out_fj, n, time)
    end do

    write(*,*)'fokkerplanck nr= ',nr,' ntau =',ntau, 'nt =', nt


    call write_v_array(vij, fij0(:,1:nr,:), time, 'maxwell')
    call write_v_array(vij,  dij(:,1:nr,:), time, 'diffusion')
    !call write_matrix(dij(1:i0,1:nr,1), time, 'diffusion')
    time2 = sys_time() - time1
    print *, 'fokkerplanck_new eval time: ', time2    
 end


 subroutine init_diffusion(h, n, vj, dj, d1, d2, d3)
    ! инициализация диффузии для схемы савельева
    implicit none
    integer, intent(in) :: n
    real*8, intent(in) :: h
    real*8, dimension(:), intent(in) :: vj, dj
    real*8, dimension(:), intent(out) :: d1, d2, d3
    real*8, dimension(:), allocatable :: xx
    integer :: i0
    integer i, klo, khi, ierr, klo1, khi1
    integer klo2, klo3, khi2, khi3, ierr1, ierr2, ierr3

    i0 = size(vj)

    allocate(xx(n+1))
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
            write(*,*)'j=', 123, ' v=', xx(i)
            write(*,*)'klo2=', klo2, 'khi2=', khi2, 'i=',i
            write(*,*)'vj(1)=', vj(1), ' vj(i0)=', vj(i0)
            pause
            stop
        end if
        if(ierr3.eq.1) then
            write(*,*)'lock error in finction d2(x)'
            write(*,*)'j=', 123, ' v=', xx(i)
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