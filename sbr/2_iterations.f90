module iterations
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64    
    implicit none
    type IterationResult

        integer :: number
        !! iteration number 'iteration=',iterat
        integer :: spectr_direction
        !! 'ispectr=',ispectr
        real(dp) :: P_launched
        !!P_launched, MW=',plaun
        real(dp) :: P_landau
        !!'P_landau, MW=',ol
        real(dp) :: P_coll
        !! 'P_coll, MW=',oc
        real(dp) :: P_alph
        !!'P_alph, MW=',oa 
        real(dp) :: alphas_power
        !!'Alphas power, MW=',fuspow
        real(dp) :: P_fast
        !write(*,*) 'P_fast (landau+coll), MW=',of
        real(dp) :: P_lost
        !write(*,*) 'P_lost, MW=',plost/xsgs
        real(dp) :: P_not_accounted
        !write(*,*) 'P_not accounted, MW=',pnab/xsgs
        real(dp) :: P_landau_strong_absorption
        !write(*,*) 'P_landau (strong absorption), MW=',ppv1/xsgs
        real(dp) :: P_landau_weak_absorption
        !write(*,*) 'P_landau (weak absorption), MW=',ppv2/xsgs
        real(dp) :: P_turns
        !write(*,*) 'P_turns, MW=', psum4/xsgs
        real(dp) :: efficiency
        !write(*,*) 'efficiency, I(MA)/P(MW)=',oi/plaun !sav2008
        !call integral(1,nspl,rh,con,avedens) !sav2010
        real(dp) :: avedens
        real(dp) :: r0
        !write (*,*) '<Ne>, m^-3=',avedens*1.d19,' R, m=',r0*1.d-2
        real(dp) :: eta_eff
        !eta_eff=1.d17*avedens*r0*oi/plaun
        !write (*,*) 'eta_eff=<Ne>*R*I/P, A/(W*m^2)=',eta_eff !sav2010
        real(dp) :: residual 
        !! невязка 'nevyazka=', pchg        
    
    contains
        procedure :: print => iteration_result_print
    
    end type IterationResult
contains
    subroutine iteration_result_print(this)
        class(IterationResult), intent(in) :: this
        print *, ' ---------'
        print *, 'ITERATION:'
        print *, 'iteration=', this%number
        print *, 'ispectr=', this%spectr_direction
        print *, 'P_launched, MW=', this%P_launched
        print *, 'P_landau, MW=', this%P_landau
        print *, 'P_coll, MW=', this%P_coll
        print *, 'P_alph, MW=', this%P_alph 
        print *, 'Alphas power, MW=', this%alphas_power
        print *, 'P_fast (landau+coll), MW=', this%P_fast
        print *, 'P_lost, MW=', this%P_lost
        print *, 'P_not accounted, MW=', this%P_not_accounted
        print *, 'P_landau (strong absorption), MW=', this%P_landau_strong_absorption
        print *, 'P_landau (weak absorption), MW=', this%P_landau_weak_absorption
        print *, 'P_turns, MW=', this%P_turns
        print *, 'efficiency, I(MA)/P(MW)=', this%efficiency 
        !call integral(1,nspl,rh,con,avedens) !sav2010
        print *, '<Ne>, m^-3=',this%avedens,' R, m=',this%r0
        !eta_eff=1.d17*avedens*r0*oi/plaun
        print *, 'eta_eff=<Ne>*R*I/P, A/(W*m^2)=',this%eta_eff
        print *, 'nevyazka=', this%residual
    end subroutine iteration_result_print
end module iterations