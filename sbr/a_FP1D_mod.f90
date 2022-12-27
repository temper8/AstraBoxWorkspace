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
        !   complex         :: inst_field1
        contains
        procedure :: set   => set_e
        procedure :: print => e_print

    end type FokkerPlanck1D   
    interface FokkerPlanck1D
       procedure :: constructor
    end interface FokkerPlanck1D
 contains
    ! implement instance constructor with no arguments
    subroutine classname_ctor0(this)
      type(FokkerPlanck1D) :: this
      !this%inst_field1 = cmplx(0.,0.) 
    end subroutine classname_ctor0
    ! implement instance constructor with one argument
    function constructor( dir, e, v_lim, v, f) result(this)
      type(FokkerPlanck1D) :: this
      real, value :: dir, e, v_lim, v(:), f(:)
      !this%inst_field1 = cmplx(0.,0.) 
      this%direction = dir
      this%enorm     = e
      this%v_lim = v_lim
      this%v = v
      this%f = f
      print *, 'exe constructor FokkerPlanck1D'
    end function constructor
    subroutine set_e(this, value)
        class(FokkerPlanck1D), intent(in out) :: this
        integer, intent(in)          :: value
        this%enorm = value
    end subroutine set_e

    subroutine e_print(this)
        class(FokkerPlanck1D), intent(in) :: this

        print *, 'direction = ', this%direction, 'e = ', this%enorm
    end subroutine e_print


 end module FokkerPlanck1D_mod