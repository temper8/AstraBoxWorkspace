module spectrum_mod
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none    


    type spectrum_point
        real(dp) nz
        !! 
        real(dp) ny
        !!
        real(dp) power
        !! power

    contains
    end type spectrum_point

    type spectrum
        integer size
        !!
        real(dp) input_power
        !!
        real(dp) power_ratio
        !!
        real(dp) max_power
        !!
        real(dp) sum_power
        !! суммарная power        
        type(spectrum_point), allocatable ::  data(:)
        !! 
    contains

    end type spectrum

    interface spectrum
        module procedure :: spectrum_constructor
        !module procedure :: read_spectrum
    end interface spectrum    
contains
    function spectrum_constructor(size) result(this)
        !- конструктор для spectrum
        implicit none
        type(spectrum) :: this
        integer, value :: size
        this%size = size
        this%input_power = 0
        this%sum_power = 0
        allocate(this%data(size))
    end function spectrum_constructor     

    function read_spectrum(file_name) result(spectr)
        !- чтение spectrum из файла
        implicit none
        type(spectrum) :: spectr
        character (len = *), value :: file_name 
        logical                     :: res
        integer i,n,stat
        real(dp) sum_power
        !integer, value :: size
        print *, file_name      
        ! Check if the file exists
        inquire( file=trim(file_name), exist=res )        
        if (.not.res) then
            print *, 'spectrum file not exists'
            stop
        end if

        open(20,file=file_name)
        n=-1
        stat=0
        do while(stat == 0)
            n=n+1
            read (20,*,iostat=stat)
        enddo

        spectr%size = n
        spectr%input_power = 0
        spectr%sum_power = 0
        spectr%power_ratio = 1
        sum_power = 0
        allocate(spectr%data(n))        
        print *,'Spectrum size = ',  n
        rewind(20)
        do i=1,n
            read (20,*) spectr%data(i)%nz, spectr%data(i)%ny, spectr%data(i)%power
            sum_power = sum_power + spectr%data(i)%power
        enddo
        spectr%sum_power = sum_power
        close(20)

    end function read_spectrum         

    subroutine divide_spectrum(spectr, pos_spectr, neg_spectr)
        !! деление спектра на две части
        implicit none
        type(spectrum), intent(in)  :: spectr
        type(spectrum), intent(out) :: pos_spectr, neg_spectr
        type(spectrum):: tmp_spectr
        type(spectrum_point) :: p
        integer i, pos_n, neg_n

        pos_spectr = spectrum(spectr%size)
        tmp_spectr = spectrum(spectr%size)
        pos_n = 0
        neg_n = 0
        do i = 1, spectr%size
            p = spectr%data(i)

            if (p%nz>0) then
                pos_n = pos_n + 1                
                pos_spectr%data(pos_n) = p
                pos_spectr%sum_power = pos_spectr%sum_power + p%power
            end if
            if (p%nz<0) then
                neg_n = neg_n + 1                
                p%nz = -p%nz
                tmp_spectr%data(neg_n) = p
                tmp_spectr%sum_power = tmp_spectr%sum_power + p%power
            endif
        end do
        pos_spectr%size = pos_n

        neg_spectr = spectrum(neg_n)
        neg_spectr%sum_power = tmp_spectr%sum_power
        do i = 1, neg_n
            neg_spectr%data(i) = tmp_spectr%data(neg_n + 1 - i)
        end do
        neg_spectr%size = neg_n

        pos_spectr%power_ratio = pos_spectr%sum_power/spectr%sum_power
        neg_spectr%power_ratio = neg_spectr%sum_power/spectr%sum_power
        pos_spectr%input_power = pos_spectr%power_ratio * spectr%input_power
        neg_spectr%input_power = neg_spectr%power_ratio * spectr%input_power        
        print *, pos_n, neg_n        
        print *, 'sum_power ', spectr%sum_power, pos_spectr%sum_power, neg_spectr%sum_power
        print *, 'power_ratio ', pos_spectr%power_ratio, neg_spectr%power_ratio
        print *, 'input_power ', spectr%input_power, pos_spectr%input_power, neg_spectr%input_power
    end subroutine

end module spectrum_mod

module spectrum1D
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    use spectrum_mod
    implicit none    

    type(spectrum) :: full_spectrum
    type(spectrum) :: pos_spectr, neg_spectr
    integer :: ispl
    !! size of spectrum
    real(dp) :: plaun
    !! power of spectrum
    real(dp) :: ynzm0(1001)
    !+
    real(dp) :: pm0(1001)
    !+ 
    real(dp) :: ynzm(1001), pm(1001)
    !! бывший common /a0a1/ ynzm(1001),pm(1001)     
    real(dp) :: pabs
    !! бывший common /a0gh/ pabs

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

        if(ispl.gt.4001) stop 'too many points in spectrum'

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

        if(ispl.gt.4001) stop 'too many points in spectrum'

    end subroutine read_negative_spectrum    

    subroutine spectrum_approximation(ispectr)
    !! approximation of input LH spectrum
        use constants, only: zero, xsgs
        use spline
        use rt_parameters, only: nnz, ntet, pabs0
        implicit none
        integer, intent(in) :: ispectr
        real(dp) yn2z(1001),powinp(1001)
        integer innz, i
        real(dp) dxx, xx0, xx1, xx2, yy1, yy2, pinp
        real(dp) dpw, dpower, pwcurr, ptot, dynn
        real(dp) pmax, pnorm
        call splne(ynzm0,pm0,ispl,yn2z)
        innz=100*ispl
        dxx=(ynzm0(ispl)-ynzm0(1))/innz
        xx2=ynzm0(1)
        yy2=pm0(1)
        pinp=0d0
        do i=1,innz
              xx1=xx2
              yy1=yy2
              xx2=xx1+dxx
              call splnt(ynzm0,pm0,yn2z,ispl,xx2,yy2,dynn)
              dpw=.5d0*(yy2+yy1)*(xx2-xx1)
              pinp=pinp+dpw
        end do
  
        dpower=pinp/dble(nnz)
        xx2=ynzm0(1)
        yy2=pm0(1)
        pwcurr=zero
        ptot=zero
        do i=1,nnz-1
            xx0=xx2
  11        continue
            xx1=xx2
            yy1=yy2
            xx2=xx1+dxx
            call splnt(ynzm0,pm0,yn2z,ispl,xx2,yy2,dynn)
            dpw=.5d0*(yy2+yy1)*(xx2-xx1)
            if(pwcurr+dpw.gt.dpower) then
                xx2=xx1+dxx*(dpower-pwcurr)/dpw
                call splnt(ynzm0,pm0,yn2z,ispl,xx2,yy2,dynn)
                dpw=.5d0*(yy2+yy1)*(xx2-xx1)
                pwcurr=pwcurr+dpw
            else
                pwcurr=pwcurr+dpw
                go to 11
            end if
            ynzm(i)=.5d0*(xx2+xx0)
            pm(i)=pwcurr
            ptot=ptot+pwcurr
            pwcurr=zero
        end do
        ynzm(nnz)=.5d0*(ynzm0(ispl)+xx2)
        pm(nnz)=pinp-ptot
        pnorm=plaun*xsgs/(pinp*ntet)
        pmax=-1d+10
        do i=1,nnz
            call splnt(ynzm0,pm0,yn2z,ispl,ynzm(i),powinp(i),dynn)
            pm(i)=pm(i)*pnorm
            if (pm(i).gt.pmax) pmax=pm(i)
            ynzm(i)=dble(ispectr)*ynzm(i) !sav2009
        end do
        pabs=pabs0*pmax/1.d2
    end subroutine

    subroutine copy_to_spectrum_1D(spectr)
        use spectrum_mod
        implicit none
        type(spectrum) :: spectr
        type(spectrum_point) ::p
        integer i
        do i= 1, spectr%size
            p = spectr%data(i)
            ynzm0(i) = p%nz
            pm0(i) = p%power
        end do        
        plaun = spectr%input_power
        ispl = spectr%size
    end subroutine

    function create_spectrum() result(spectr)
        use spectrum_mod
        use rt_parameters, only: nnz
        implicit none
        type(spectrum) :: spectr
        type(spectrum_point) :: p
        integer i 
        real(dp) :: pmax
        pmax = 0
        spectr = spectrum(nnz)
        do i= 1, nnz
            p = spectrum_point(nz= ynzm(i), ny = 0, power = pm(i))
            if (pm(i) > pmax) pmax=pm(i)
            spectr%data(i) = p
        end do
        spectr%max_power = pmax
    end function

    subroutine write_spectrum(ispectr)
        implicit none
        integer, intent(in) :: ispectr        
        !       call get_unit(iunit)
        !       if(iunit.eq.0) then
        !        write(*,*)'no free units up to 299'
        !        pause
        !        stop
        !       end if
        !       if(ispectr.eq.1) then
        !        open(iunit,file='lhcd/out/used_spectrP.dat')
        !       else if(ispectr.eq.-1) then
        !        open(iunit,file='lhcd/out/used_spectrM.dat')
        !       end if
        !       do i=1,nnz
        !        write(iunit,1008) ynzm(i),powinp(i)
        !       end do
        !       write(iunit,*)
        !      close(iunit)
        !1008   format (1x,10(e14.7,3x))        
    end subroutine
end module spectrum1D