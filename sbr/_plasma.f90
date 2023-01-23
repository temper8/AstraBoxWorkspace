module plasma
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64        
    implicit none
    integer ngrid, nspl
    !+ ASTRA radial grid number
    real(dp) tcur
    !+ время (придумать название для переменной получше)
    real(dp) rm
    !+ minor radius in mid-plane, cm
    real(dp) b_tor0, b_tor
    !+ временно нужно две переменных, тоже нужно исправить
    real(dp) r0
    real(dp) z0
    real(dp) rh1
    real(dp), dimension(:),allocatable:: con,tem,temi,zeff,afld
    real(dp), dimension(:),allocatable:: rh,rha,drhodr,delta,ell,gamm,amy

    integer, parameter :: ipsy = 5, ncoef = 5
    !+   ipsy = number of polinomial decomposition coefficients
    !+   used for interpolation of Zakharov's moments.
    real(dp), dimension(ipsy) :: cdl,cly,cgm,cmy,coeffs


    real(dp) y2dn(501),y2tm(501),y2tmi(501)
    !+ бывший common /a0l3/
    real(dp) y2zeff(501)
    !+ бывший common /a0l5/ 

    integer ncheb
    real(dp) chebne(50),chebdne(50),chebddne(50)    
    !+ бывший common/ne_cheb

contains
    subroutine init_plasma(NA1, ABC, BTOR, RTOR, UPDWN, GP2, AMETR, RHO, SHIF, ELON, TRIA,MU, NE, TE, TI, ZEF, UPL)
        use constants
        use approximation
        use rt_parameters
        use spline
        use chebyshev
        implicit none
        integer, intent(in)  :: NA1
        real(dp), intent(in) :: ABC, BTOR, RTOR, UPDWN, GP2
        real(dp), dimension(*) :: AMETR, RHO, SHIF, ELON, TRIA,MU,  NE, TE, TI, ZEF, UPL
        integer i, k
        integer, parameter :: N  = 501
        real(dp) :: znak_tor, znak_pol, fpol, dfmy
 
        ngrid = NA1
        nspl = ngrid
        if (.not. allocated(rh)) then
            allocate(rh(N),rha(N),drhodr(N),con(N),tem(N))
            allocate(temi(N),zeff(N), afld(N))
            allocate(delta(N),ell(N),gamm(N),amy(N))
        end if
        do i=1, ngrid
            rh(i)=AMETR(i)/ABC
            rha(i)=RHO(i)/ABC  !/ABC instead of /ROC is not a mistake!
            delta(i)=(SHIF(1)-SHIF(i))/ABC  !FRTC Shafr. shift. defin.
            ell(i)=ELON(i)
            gamm(i)=rh(i)*TRIA(i)
            con(i)=NE(i)
            tem(i)=TE(i)
            temi(i)=TI(i)
            zeff(i)=ZEF(i)
            afld(i)=UPL(i)/RTOR/GP2 !!variant
        end do
        !rh(ngrid)=1.d0
        rh1=rh(1)          !saving the first ASTRA radial grid element
        rh(1) = 0.0d0         !shifting the first element to zero
        rha(1) = 0.0d0        !shifting the first element to zero
        delta(1) = 0.0d0      !putting delta(rh=0.)=0.
        gamm(1) = 0.0d0       !putting gamm(rh=0.)=0.
  
        b_tor0=1.d4*BTOR*RTOR/(RTOR+SHIF(1)) !B_tor_(magnetic axis), Gauss
        rm=1.d2*ABC                       !minor radius in mid-plane, cm
        r0=1.d2*(RTOR+SHIF(1))     !x-coordinate of the magnetic axis, cm
        z0=1.d2*UPDWN              !z-coordinate of the magnetic axis, cm


    !   spline approximation of plasma profiles         
    !
    !   shift as a function of "minor radius":
        call approx(rh,delta,ngrid,polin1,ipsy-1,coeffs)
        cdl(1)=0.0d0
        do k=2,ipsy
            cdl(k)=coeffs(k-1)
        end do
 
    !   triangularity as a function of "minor radius":
        call approx(rh,gamm,ngrid,polin1,ipsy-1,coeffs)
        cgm(1)=0.0d0
        do k=2,ipsy
            cgm(k)=coeffs(k-1)
        end do
 
    !   ellipticity as a function of "minor radius":
        call approx(rh,ell,ngrid,polin,ipsy,cly)            

    !  "poloidal magnetic field":
        call diff(rh,rha,ngrid,drhodr)
 
        do i=2,ngrid
            amy(i)=1.d4*BTOR*MU(i)*rha(i)*drhodr(i)
            !print *, amy(i), BTOR, MU(i)
        end do
        !print *, '----------------'
        amy(1)=0.d0  
        
    !! amy=(btor/q)*rho*(drho/dr) is a function of "minor radius" r=rh(i).
    !! Poloidal magnetic field: B_pol=amy(r)*sqrt(g22/g), where g is
    !! determinant of 3D metric tensor and g22 is the (22) element of
    !! the tensor, normalized on ABC^4 and ABC^2, correspondingly.
    !!
    !!  Polinomial approximation of the amy(r):
    !    inpt2=ngrid-3
        call approx(rh,amy,ngrid-3,polin1,ipsy-1,coeffs)
        cmy(1)=0.d0
        do k=2,ipsy
         cmy(k)=coeffs(k-1)
        end do
  
        ! зачем-то меняет знак коэффициентов????
        znak_tor=dsign(1.d0,dble(itor))
        b_tor=znak_tor*dabs(b_tor0)
        fpol=fdf(1.d0,cmy,ncoef,dfmy)
        znak_pol=dsign(1.d0,dble(i_pol))*dsign(1.d0,fpol)
        do i=1,ncoef
         cmy(i)=znak_pol*cmy(i)
        end do

        
    !!!!!!!!!!!!!!! spline approximation of plasma profiles !!!!!!!!!!!!!!!!
        call splne(rh,con,nspl,y2dn)
        call splne(rh,tem,nspl,y2tm)
        call splne(rh,zeff,nspl,y2zeff)
        call splne(rh,temi,nspl,y2tmi)

        if(inew.ne.0) then
            ncheb=20
            call chebft1(zero,1.d0,chebne,ncheb,fn)
            call chder(zero,1.d0,chebne,chebdne,ncheb)
            call chder(zero,1.d0,chebdne,chebddne,ncheb)
        end if        
    end subroutine

    double precision  function fn(x)
    ! plasma  density,  cm^-3
        use spline      
    !      use plasma
        implicit real*8 (a-h,o-z)
        !common /a0l3/ y2dn(501),y2tm(501),y2tmi(501)
        !common /a0l4/ con(501),tem(501),temi(501),nspl
        parameter(zero=0.d0,alfa=4.d0,dr=.02d0)
        pa=dabs(x)
        if(pa.le.rh(nspl)) then
            call splnt(rh,con,y2dn,nspl,pa,y,dy)
        else
            r=pa-rh(nspl)
            y=con(nspl)*dexp(-alfa*(r/dr)**2)
        end if
        fn=y*1.d+13    !cm^-3
    end    
        
end module plasma
