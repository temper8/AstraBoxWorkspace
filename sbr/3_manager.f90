module manager_mod
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    real(dp) :: yn3
    !! common /abefo/ yn3
    integer :: iroot
    !!common /beo/ iroot
    integer :: ivar
    !!common /bdeo/ ivar    
    integer :: icall1,icall2
    !common /aef2/ icall1,icall2
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine manager(iterat,iw0, ntet, spectr)
        use constants            
        use plasma
        use rt_parameters, only : nr, ipri, iw, nmaxm, pabs0            
        !use spectrum1D, only: ynzm, pm
        use trajectory
        use spectrum_mod
        use iterator_mod,only: plost, pnab
        implicit none
        type (spectrum) spectr
        type (spectrum_point) point
        real*8 pabs
        integer iznzap(mpnt),iwzap(mpnt),irszap(mpnt)
        real*8 rzap(mpnt),tetzap(mpnt),xmzap(mpnt),yn3zap(mpnt)
        !common /a0a1/ ynzm(1001),pm(1001) 
        !common /a0a2/ tet1,tet2
        ! real*8 plost,pnab
        !common /a0a4/ plost,pnab
        real*8 rzz,tetzz,xmzz
        common /abc/ rzz,tetzz,xmzz,iznzz,iwzz,irszz
        common /abcd/ irs
        common /abcde/ izn
        common /abcdg/ iabsorp
        !common /abefo/ yn3
        real*8 pow
        common /acg/ pow
        !common /a0gh/ pabs
        !common /aef2/ icall1,icall2
        common /ag/ inak,lenstor,lfree
        !common/refl/nrefj(mpnt)
        integer lenstor, ntet, irs, iout, itr, inak, nnj,  n_it
        integer maxref, iterat, nmax0, ibad, itet, nref
        integer nbad1, nbad2, inz
        integer iw0, ifail, iabsirp, inak0,ib,ie,izn
        integer lfree, nmax, iabsorp, i, nb1,nb2
        integer iznzz, iwzz, irszz
        real(dp) htet, hr, yn, rin, xmin, rstart
        real(dp) xnr, powexit, dltpow,  pow1, pgamma, xm
        real(dp) tetin0, tetin, tet

        pabs = spectr%max_power*pabs0/1.d2
        ! print *, 'pabs =', pabs, spectr%max_power, pabs0
        lenstor = length
        htet = zero
        hr = 1.d0/dble(nr+1) !sav2008
        if (ntet.ne.1) htet = (tet2-tet1)/(ntet-1)
        irs = 1
        iout = 0
        mbeg(1) = 1
        itr = 0
        inak = 0
        nnj = 0
        do n_it = 0,3
            nnj = nnj+nmaxm(n_it+1)
        end do
        maxref = nnj
        if (iterat.lt.3) nmax0=nmaxm(iterat+1)
        if (iterat.ge.3) nmax0=nmaxm(4)
        if (ipri.gt.1) then
            write(*,1001) iterat+1
            write(*,1002)
        end if
        ibad = 0
        !--------------------------------------
        ! begin outer loop on teta
        !--------------------------------------
        do itet = 1,ntet
            nref = 0
            nbad1 = 0
            nbad2 = 0
            icall1 = 0
            icall2 = 0
            tetin = tet1+htet*(itet-1)
            !--------------------------------------
            ! begin inner loop on nz
            !--------------------------------------
            do inz = 1, spectr%size
                itr = itr+1
                !ipri          if(ipri.eq.4)  write(23,*)
                point = spectr%data(inz)
                if(iterat.eq.0) then
                    !-----------------------------------------
                    !    find initial radius for a trajectory
                    !    on the 1th iteration
                    !-----------------------------------------
                    yn = point%nz
                    pow = point%power
                    ! print *, 'init pow', pow
                    !yn=ynzm(inz) !sav2008, yn is introduced
                    !pow=pm(inz)
                    irs = 1
                    iw = iw0
                    rin = rini(xmin,tetin,xnr,point,hr,ifail)
                    if (ifail.eq.1) then
                        if (ipri.gt.1) write (*,*) 'error: no roots'
                        iabsorp = -1
                        inak0 = inak
                        go to 10
                    end if
                    rbeg(itr) = rin !sav2008
                    tetbeg(itr) = tetin !sav2008
                    xnrbeg(itr) = xnr !sav2008
                    xmbeg(itr) = xmin !sav2008
                    yn3beg(itr) = yn3 !sav2008
                else
                    if (mbad(itr).ne.0) then
                        plost = plost+point%power
                        go to 31
                    end if
                    ib = mbeg(itr)
                    ie = mend(itr)
                    powexit = point%power
                    dltpow = pabs
                    call dqliter(dltpow,ib,ie,hr,powexit,iout)
                    if (nmax0.eq.0) then
                        ib = mbeg(itr)
                        ie = mend(itr)
                        pow1 = powexit
                        pgamma = 1.d0-pow1/point%power
                        powexit = pow1/pgamma
                        dltpow = powexit-pow1+pabs
                        call dqliter(dltpow,ib,ie,hr,powexit,iout)
                        powexit = powexit-dltpow+pabs
                        if (powexit.lt.zero) powexit=zero
                        go to 30
                    end if
	                if (iout.eq.0) then
                        go to 30
                    else
                        tetin = tetzap(itr)
                        xmin = xmzap(itr)
                        rin = rzap(itr)
                        yn3 = yn3zap(itr)
                        pow = powexit
                        irs = irszap(itr)
                        iw = iwzap(itr)
                        izn = iznzap(itr)
                        jrad(ie+1) = 1
                        dland(ie+1) = lfree
                        inak = lfree-1
                    end if
                end if
                !---------------------------------------
                ! initial parameters for a trajectory
                !---------------------------------------
                xm = xmin
                rstart = rin !sav2008
                tet = tetin
                nmax = nmax0
                iabsorp = 0
                inak0 = inak
                ! print *, pow, nmax, icall1, iabsorp
                !-------------------------------------
                ! call ray tracing
                !-------------------------------------
                call traj(xm,tet,rstart,nmax,nb1,nb2,itet,inz, pabs) !sav2009
                nbad1 = nbad1+nb1
                nbad2 = nbad2+nb2
                nrefj(itr) = nrefj(itr)+nmax
                ! print *, pow, nmax, icall1, iabsorp
                powexit = pow
                nref = nref+nmax
10              if (iabsorp.lt.0) then
                    !-------------------------------------
                    !    encounted problems
                    !-------------------------------------
                    if (inak.eq.lenstor-1) then
                        write (*,*) 'fix maximal length'
                        nmax0 = 0
                        do i=1,4
                            nmaxm(i) = 0
                        end do
                        iout = 1
                        goto 20
                    end if
                    if (ipri.gt.1) then
                        tetin0=tet1+htet*(itet-1)
                        write (*,111) tetin0, point%nz

111                     format(1x,'traj. with tet0=',f10.5,1x,', Ninput=',f10.5,1x,'failed')
                    end if
                    mbad(itr) = 1
                    plost= plost+pow
                    inak = inak0
                    mend(itr) = inak-1
                    goto 30
                end if
                ! print *, 'plost', plost
                !---------------------------------------
                ! remember end point of trajectory
                !---------------------------------------
                rzap(itr) = rzz
                tetzap(itr) = tetzz
                xmzap(itr) = xmzz
                yn3zap(itr) = yn3
                iznzap(itr) = iznzz
                iwzap(itr) = iwzz
                irszap(itr) = irszz
                if (iterat.eq.0) then
                    if (itr.gt.1) mbeg(itr) = mend(itr-1)+2
                    mend(itr) = inak
                    jrad(mend(itr)+1) = 0
                    lfree = mend(itr)+2
                    inak = lfree-1
                end if
20              continue
                if(iout.ne.0) then
                    dcoll(ie+1) = inak
                    jrad(inak+1) = 0
                    lfree = inak+2
                end if
                if(nrefj(itr).gt.maxref.and.pow.gt.pabs) then !forced absorp
                    if(pow.ge.point%power) go to 30 !sav2008
                    ib = mbeg(itr)
                    ie = mend(itr)
                    pow1 = pow
                    pgamma = 1.d0-pow1/point%power
                    powexit = pow1/pgamma
                    dltpow = powexit-pow1+pabs
                    ! print *, '1 -- ', powexit, dltpow, pabs
                    call dqliter(dltpow,ib,ie,hr,powexit,iout)
                    powexit = powexit-dltpow+pabs
                    ! print *, '2 -- ', powexit, dltpow, pabs
                    if(powexit.lt.zero) powexit=zero
                end if
30              continue
                pnab = pnab+powexit
31              continue
            end do
            if(ipri.gt.1) write(*,1003)itet,icall1,icall2,nref,lfree-1,nbad1,nbad2
            ! print *,'itr =', itr
            ! print *, pnab,  dltpow
            !pause
        end do
1001    format (30x,i4,' iteration')
1002    format (6x,'n',5x,'call2',6x,'call4',6x,'nrefl',4x,'last',5x,'bad2',5x,'bad4')
1003    format (3x,i4,2(1x,i10),2x,i7,2x,i8,2(1x,i7),2(2x,i7))
1004    format(1x,i8)
1005    format(1x,i5)
1006    format (e14.7)
    end    


    real*8 function rini(xm,tet,xnr,point,hr,ifail) !sav2009
        use constants, only : zero
        use rt_parameters, only : inew
        use spectrum_mod
        implicit none
        type(spectrum_point) :: point
        real(dp) xm, tet,xnr,hr
        integer ifail, ntry
        real(dp) :: vgrp(3),vph(3)
        real(dp) :: ynz,ynpopq
        common /bcef/ ynz,ynpopq
        real(dp) :: g11,g12,g22,g33,gg,g,si,co
        common/metrika/g11,g12,g22,g33,gg,g,si,co !sav2009
        real(dp) :: pa, prt, prm
        real(dp) :: f1,f2
        real(dp),  parameter :: rhostart=1.d0
        integer, parameter :: ntry_max=5
        ifail = 1
        rini = zero
        ntry = 0
        pa = rhostart
        do while (ntry.lt.ntry_max.and.pa.ge.2d0*hr)
            pa = rhostart-hr*dble(ntry)-1.d-4
            ntry = ntry+1
            ivar = 1
            call disp2(pa,xm,tet,xnr,prt,prm)
            if (inew.gt.0) then !g' in ST and poloidal grill direction
                !yn3 = zero                 !Nfi=0
                !xm = yn*dsqrt(g22)/si      !given Npar at Nfi=0
                yn3 = point%ny**dsqrt(g33)/co     
                xm = point%nz*dsqrt(g22)/si
!!              xm=yn*dsqrt(g22)         !given yn=(N*jpol) at Nfi=0
            else !usual tokamak and toroidal grill direction
                !xm = zero               !N2=0
                !yn3 = yn*dsqrt(g33)/co  !if given Npar at Nteta=0
!!              yn3=yn*dsqrt(g33)       !if given Nfi at Nteta=0
                yn3 = point%nz*dsqrt(g33)/co    
                xm = point%ny*dsqrt(g22)/si
            end if
            ivar = 0
            iroot = 2
            call disp2(pa,xm,tet,xnr,f1,f2)
            if (f1.ge.zero.and.f2.ge.zero) then
                rini = pa
                ifail = 0
                return
            end if
        end do
      end      
end module manager_mod
