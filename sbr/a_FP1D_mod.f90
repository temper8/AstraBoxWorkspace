module FokkerPlanck1D_mod ! the module name defines the namespace
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    type FokkerPlanck1D 
        !- solver of FP eq
        integer          :: direction
        !- direction
        real(dp)         :: enorm
        !- электрическое поле
        real(dp)         :: v_lim
        !- верхняя граница скорости электронов
        real(dp), allocatable         :: v(:)
        !- сетка скоростей
        real(dp), allocatable         :: f(:)
        !- распределение
        integer         :: i0
        !- size of distribution grid
        real(dp)        :: alfa2  
        !- поле со знаком
        integer         :: n
        !- size of local grid        
        real(dp)        :: h  
        !- step of local grid 
        real(dp), allocatable ::  d1(:), d2(:), d3(:)
        !- диффузия
        contains
        procedure :: print => FokkerPlanck1D_print
        procedure :: solve_time_step => FokkerPlanck1D_solve_time_step
        procedure :: init_zero_diffusion => FokkerPlanck1D_init_zero_diffusion

    end type FokkerPlanck1D   
    interface FokkerPlanck1D
        procedure :: FokkerPlanck1D_constructor
    end interface FokkerPlanck1D

    contains
    function FokkerPlanck1D_constructor(dir, e, v_lim, v, f) result(this)
        !- конструктор для FokkerPlanck1D
        type(FokkerPlanck1D) :: this
        integer, value :: dir
        real(dp), value :: e, v_lim, v(:), f(:)
        integer  :: n
        real(dp) :: h
        real(dp), parameter :: h0 = 0.1d0
        !this%inst_field1 = cmplx(0.,0.) 
        this%direction = dir
        this%enorm     = e
        this%v_lim = v_lim
        this%v = v
        this%f = f
        this%i0 = size(v)
        this%alfa2=this%direction*this%enorm
        n = v_lim/h0-1
        h = v_lim/dble(n+1)
        if (h.gt.h0) then
            n = n+1
            h = v_lim/dble(n+1)
        end if
        this%n = n
        this%h = h
        print *, 'exe constructor FokkerPlanck1D'
    end function FokkerPlanck1D_constructor

    subroutine FokkerPlanck1D_print(this)
        class(FokkerPlanck1D), intent(in) :: this

        print *, 'direction = ', this%direction, 'e = ', this%enorm, 'i0 =', this%i0
    end subroutine FokkerPlanck1D_print

    subroutine FokkerPlanck1D_init_zero_diffusion(this)
      implicit none
      class(FokkerPlanck1D), intent(inout) :: this
      integer :: n
      n = this%n
      allocate(this%d1(n+1),this%d2(n+1),this%d3(n+1))
      this%d1(:)=0d0
      this%d2(:)=0d0
      this%d3(:)=0d0
    end subroutine FokkerPlanck1D_init_zero_diffusion

    subroutine FokkerPlanck1D_solve_time_step(this, dt, nt)
        implicit none
        class(FokkerPlanck1D), intent(inout) :: this
        
        integer, intent(in) :: nt
        real*8, intent(in)  :: dt


        !real*8, intent (inout), optional :: dfj0(:)

        real(dp), parameter :: zero=0.d0
        real(dp)  y(this%n+2),x(this%n+2)
        real(dp), dimension(:), allocatable:: fj, dfj,  givi
        integer i, ii, it, ibeg, klo, khi, ierr, klo1, khi1
        real(dp) shift, ybeg, yend, tend, dff

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
  
  
        subroutine burying_procedure(v, f0, df0)
            ! процедура закапывания
            implicit none
            real*8,  intent(in)     :: v(:)        
            real*8,  intent(inout)  :: f0(:)
            real*8,  intent (inout), optional :: df0(:)
        end subroutine            
        end interface
 
        !!!!!! grid !!!!!!!!!
        !!  shift=h*0.1d0 !0.01d0
        do i=1, this%n+2
            x(i)=this%h*dble(i-1) !+shift
        end do
  
        do i=1, this%n+1
            call lock(this%v,this%i0,x(i+1),klo,khi,ierr)
            if(ierr.eq.1) then
                write(*,*)'lock error #1 in finction fokkerplanck'
                write(*,*)'j=', 123, ' v=', x(i+1)
                write(*,*)'vj(1)=',this%v(1),' vj(i0)=',this%v(this%i0)
                pause
                stop
            end if
            call linf(this%v,this%f,x(i+1),y(i),klo,khi)
        end do
  
        ybeg = this%f(1)  !boundary conditions
        yend = this%f(this%i0) !zero
        !print *, ' yend =', yend 
        !!!!!!!!!!!!   solve problem   !!!!!!!!!!!!!!!!!!!!!!!!!!
        !call savelyev_solver(alfa2, nt, h, dt, n, ybeg, yend, d1,d2,d3, y)
        call teplova_khavin_solver(this%alfa2, nt, this%h, dt, this%n, ybeg, yend, this%d1,this%d2,this%d3, y)
        allocate(fj(this%n+2))
        fj(1)=ybeg
        fj(this%n+2)=yend
        do i=1,this%n
            fj(i+1)=y(i)
        end do

        do i=2,this%i0-1
          if(this%v(i).lt.this%v_lim) then
              call lock(x,this%n+2,this%v(i),klo,khi,ierr)
              if(ierr.eq.1) then
                  write(*,*)'lock error #2 in finction fokkerplanck'
                  write(*,*)'j=', 123 ,' vij=',this%v(i)
                  write(*,*)'x(1)=', x(1),' x(n+2)=',x(this%n+2)
                  pause
                  stop   
              end if
              call linf(x,fj,this%v(i),this%f(i),klo,khi)
          else
            this%f(i)=zero
          end if
      end do
      deallocate(fj)
  
      !if (present(dfj0)) then
         ! call burying_procedure(vj, fj0, dfj0)
      !else 
        !  call burying_procedure(vj, fj0)
      !end if
  
  
    end subroutine FokkerPlanck1D_solve_time_step


 end module FokkerPlanck1D_mod