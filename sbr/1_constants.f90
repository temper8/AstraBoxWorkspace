module constants 
    !! модуль с математическими и физическими константами
    use kind_module
    implicit none
    real(wp), parameter :: pi=4.d0*datan(1.d0)
    !! число Пи = 3.1415....
    real(wp), parameter :: pi2=2.d0*pi
    real(wp), parameter :: pi4=4.d0*pi
    real(wp), parameter :: piq=dsqrt(pi)

    real(wp), parameter :: zero=0.d0      
    real(wp), parameter :: one=1.d0
    real(wp), parameter :: two=2.d0
    real(wp), parameter :: tiny=1.d-100
    real(wp), parameter :: tin=1d-7

    real(wp), parameter :: talfa=3.5d0    
    !! alpha particles' birth energy, MeV
    real(wp), parameter :: zalfa=2.d0     
    !! alpha particles' electrical charge
    real(wp), parameter :: xmalfa=4.d0    
    !! alpha particles' atomic mass

    real(wp), parameter :: clt=3.0d+10
    !! скорость света
    real(wp), parameter :: pme=9.11e-28
    real(wp), parameter :: pqe=4.803e-10
    real(wp), parameter :: xlog=16.d0+dlog(16.d0)
    real(wp), parameter :: c0=dsqrt(pi4*pqe**2/pme)
    real(wp), parameter :: c1=pqe/pme/clt
    real(wp), parameter :: xsgs=1d+13
    real(wp), parameter :: xwtt=1d-7    

    real(wp), parameter :: cnst1=0.2965924106d-6
    !! cnst1=(m_e/m_p)**2, CGS
    real(wp), parameter :: cnst2=0.359680922d-35
    !! cnst2=(m_e/e)**2,  CGS
contains
    
end module constants 