module constants 
    !! модуль с математическими и физическими константами
    use kind_module
    implicit none

    real(wp), parameter :: zero = 0.0_wp 
    real(wp), parameter :: one  = 1.0_wp
    real(wp), parameter :: two  = 2.0_wp
    real(wp), parameter :: tiny = 1.e-100_wp
    real(wp), parameter :: tin  = 1e-7_wp

    real(wp), parameter :: pi = acos(-one) !! число Пи = 3.1415....
    real(wp), parameter :: pi2 = 2.d0*pi
    real(wp), parameter :: pi4 = 4.d0*pi
    real(wp), parameter :: piq = sqrt(pi)

    real(wp), parameter :: talfa  = 3.5_wp    !! alpha particles' birth energy, MeV
    real(wp), parameter :: zalfa  = 2.0_wp    !! alpha particles' electrical charge
    real(wp), parameter :: xmalfa = 4.0_wp  !! alpha particles' atomic mass

    real(wp), parameter :: clt = 3.0e+10_wp  !! скорость света

    real(wp), parameter :: pme  = 9.11e-28_wp
    real(wp), parameter :: pqe  = 4.803e-10_wp
    real(wp), parameter :: xlog = 16.0_wp + dlog(16.0_wp)
    real(wp), parameter :: c0 = sqrt(pi4*pqe**2/pme)
    real(wp), parameter :: c1 = pqe/pme/clt
    real(wp), parameter :: xsgs = 1e+13_wp
    real(wp), parameter :: xwtt = 1e-7_wp

    real(wp), parameter :: cnst1 = 0.2965924106e-6_wp
    !! cnst1=(m_e/m_p)**2, CGS
    real(wp), parameter :: cnst2 = 0.359680922e-35_wp
    !! cnst2=(m_e/e)**2,  CGS
contains
    
end module constants 