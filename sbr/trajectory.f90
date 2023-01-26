module trajectory
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64       
    implicit none
    integer, parameter :: length = 5000000
    integer, parameter :: mpnt = 10000

    real(dp) dland(length),dcoll(length),perpn(length),dalf(length)
    real(dp) vel(length),jrad(length),iww(length),tetai(length)
    real(dp) xnpar(length),izz(length)
    !! бывший common/agh/xnpar,vel,dland,dcoll,dalf,perpn,tetai,jrad,iww,izz
    real(dp) mbeg(mpnt),mend(mpnt),mbad(mpnt),rbeg(mpnt) !sav2008
    real(dp) tetbeg(mpnt),xnrbeg(mpnt),xmbeg(mpnt),yn3beg(mpnt)
    !! common/viewdat/mbeg,mend,mbad,rbeg,tetbeg,xnrbeg,xmbeg,yn3beg   
    data mbad /mpnt*0/
contains
    
end module trajectory