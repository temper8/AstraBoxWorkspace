module current
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    real(dp) :: dql(101,100)
    !!
    real(dp) :: pdl(100)
    !!
    real(dp) :: vzmin(100)
    !!
    real(dp) :: vzmax(100)
    !common /a0i3/ dql(101,100),pdl(100),vzmin(100),vzmax(100)
    real(dp) :: fcoll(100)
    real(dp) :: dens(100) 
    real(dp) :: eta(100)
    !common /a0i4/ fcoll(100),dens(100),eta(100)
    real(dp) :: dq1(101,100)
    real(dp) :: dq2(101,100)
    real(dp) :: pdc(100)
    real(dp) :: pda(100)
    real(dp) :: ppv1,ppv2
    !common/vvv1/dq1(101,100),dq2(101,100),pdc(100),pda(100),ppv1,ppv2
    real(dp) :: pdfast(100)
    !common /vvv3/ pdfast(100)
    real(dp) :: dqi0(50,100) 
    !common /alph/ dqi0(50,100)    
    real(dp) :: dncount(101,100)
    !common/findsigma/dncount(101,100)
contains

subroutine dfind(j,i,v,powpr,pil,pic,pia,df,decv,refr,vlf,vrt,ifast)
    use constants
    use plasma
    use rt_parameters
    implicit real*8 (a-h,o-z)
    integer i,j,ifast,k
     !common /a0i3/ dql(101,100),pdl(100),vzmin(100),vzmax(100)
     !common /a0i4/ fcoll(100),dens(100),eta(100)
     !common/vvv1/dq1(101,100),dq2(101,100),pdc(100),pda(100),ppv1,ppv2
     !common /vvv3/ pdfast(100)
     !common /alph/ dqi0(50,100)
     !common/findsigma/dncount(101,100)
     

    if(v.gt.cltn) return
    if(pil.gt.zero) then
        if(v.lt.vzmin(j)) vzmin(j)=v
        if(v.gt.vzmax(j)) vzmax(j)=v
    end if
    pchgl=zero
    pchgc=zero
    pchga=zero
    denom=pil+pic+pia
    powlandau=1.d0-dexp(-2.d0*pil)
    powdamped=1.d0-dexp(-2.d0*denom)
    domin=powpr*powdamped
    if(denom.ne.zero) then
!!       pchgl=powpr*(1.d0-dexp(-2d0*pil))
!!       pchgc=powpr*dexp(-2d0*pil)*dabs(-2d0*pic)
!!       pchga=powpr*dexp(-2d0*pil)*dabs(-2d0*pia)
        fff=domin/denom
        pchgl=dabs(pil*fff)
        pchgc=dabs(pic*fff)
        pchga=dabs(pia*fff)
    end if
    dd=zero
    if(pil.eq.zero) go to 1 !no Landau absorption
    if(powlandau.gt.pchm) then !strong absorption
        ppv1=ppv1+pchgl
        if(dabs(df).gt.tiny) then
            dd=dabs(-pchgl/vk(j)/(df*1.d10))
            dncount(i,j)=dncount(i,j)+1.d0
        else
            dd=zero
        end if
        dq1(i,j)=dq1(i,j)+dd
    else  ! weak absorption
        ppv2=ppv2+pchgl
        dd=dabs(2.d0*decv*powpr*1.d-10/vk(j))
        dncount(i,j)=dncount(i,j)+1.d0
        dq2(i,j)=dq2(i,j)+dd
    end if

1   continue
    dql(i,j)=dql(i,j)+dd
    pdl(j)=pdl(j)+pchgl
    pdc(j)=pdc(j)+pchgc
    pda(j)=pda(j)+pchga

    if(ifast.eq.-1) pdfast(j)=pdfast(j)+pchgl+pchgc+pchga
    if(itend0.gt.0) then
        parn=cltn/v
        dvz=vrt-vlf
        dnpar=cltn*dvz/v**2
        weight=(refr**2-eta(j))**2/(refr**2*parn**3)
!!!        adde=zze*(dd/dens(j))*weight
!!!        e2perp(i,j)=e2perp(i,j)+adde
        addd=zza*(dd/dens(j))*weight/fcoll(j)/refr**3
        arg=clt/(refr*valfa)
        do k=1,kv
            if(vperp(k,j).gt.arg) then
                hevis=dsqrt((vperp(k,j)-arg)*(vperp(k,j)+arg))
                adda=addd*hevis
                dqi0(k,j)=dqi0(k,j)+adda*dnpar
            end if
        end do
     end if
     return
     end    
end module current

module power
    implicit none
    
contains
    
end module power