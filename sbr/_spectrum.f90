module spectrum1D
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none    
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
end module spectrum1D