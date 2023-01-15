module constants 
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64    
    implicit none
    real(dp), parameter :: pi=4.d0*datan(1.d0)
    real(dp), parameter :: pi2=2.d0*pi
    real(dp), parameter :: pi4=4.d0*pi
    real(dp), parameter :: piq=dsqrt(pi)

    real(dp), parameter :: talfa=3.5d0    
    !+ alpha particles' birth energy, MeV
    real(dp), parameter :: zalfa=2.d0     
    !+ alpha particles' electrical charge
    real(dp), parameter :: xmalfa=4.d0    
    !+ alpha particles' atomic mass
    real(dp), parameter :: tin=1d-7
    real(dp), parameter :: clt=3.0d+10
    real(dp), parameter :: pme=9.11e-28
    real(dp), parameter :: pqe=4.803e-10
    real(dp), parameter :: xlog=16.d0+dlog(16.d0)
    real(dp), parameter :: c0=dsqrt(pi4*pqe**2/pme)
    real(dp), parameter :: c1=pqe/pme/clt
    real(dp), parameter :: xsgs=1d+13
    real(dp), parameter :: xwtt=1d-7    
contains
    
end module constants 