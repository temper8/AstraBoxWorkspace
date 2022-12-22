module FokkerPlanck1D_mod ! the module name defines the namespace
    type FokkerPlanck1D ! classname is the class prototype name
        integer        :: direction
        real*8         :: enorm

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
    function constructor( dir, e) result(this)
      type(FokkerPlanck1D) :: this
      real,value :: dir, e
      !this%inst_field1 = cmplx(0.,0.) 
      this%direction = dir
      this%enorm     = e
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

    ! implement instance constructor with two arguments
    subroutine classname_ctor2(this,Value1,Value2)
      type(FokkerPlanck1D) :: this
      real,value :: Value1,Value2
      !this%inst_field1 = cmplx(Value1,Value2)
    end subroutine classname_ctor2
 end module FokkerPlanck1D_mod