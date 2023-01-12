module spectrum1D
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none    
    integer :: ispl
    !+ size of spectrum
    real(dp) :: plaun
    !+ power of spectrum
    real(dp) :: ynzm0(1001)
    !+
    real(dp) :: pm0(1001)
    !+ 
    integer, parameter, private :: HEADER_LENGTH = 53

    contains
    subroutine read_positive_spectrum(file_name, p_in)
        implicit none
        character(*) file_name
        real(dp) :: p_in        
        integer, parameter :: iunit = 20        
        integer :: i, i1
        real(dp) :: anz,apz
        
        open(iunit, file= file_name)
        do i = 1, HEADER_LENGTH
            read(iunit,*)
        end do
        do i=1,10000
            read (iunit,*) anz,apz
            if(apz.eq.-88888.d0) then
                plaun=p_in*anz !input power in positive spectrum
                exit
            end if
            ynzm0(i)=anz
            pm0(i)=apz
            i1=i
        end do
        close(iunit)
        ispl=i1
    end subroutine read_positive_spectrum        

    subroutine read_negative_spectrum(file_name, p_in)
        implicit none
        character(*) file_name
        real(dp) :: p_in
        integer, parameter :: iunit = 20        
        integer :: i, i1
        real(dp) :: anz, apz

        open(iunit, file= file_name)        
        do i = 1, HEADER_LENGTH
            read(iunit,*)
        end do
        apz=0.d0
        do while(apz.ne.-88888.d0)
            read (iunit,*) anz,apz
        end do
        read(iunit,*)
        plaun=p_in*(1.d0-anz) !input power in negative spectrum
        if (plaun > 0.d0) then
            do i=1,10000
                read (iunit,*,end=10) ynzm0(i),pm0(i)
            i1=i
            end do        
        end if
10      close(iunit)
        ispl=i1        
    end subroutine read_negative_spectrum    

end module spectrum1D