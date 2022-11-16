 
subroutine write_distribution(fj,N,time)
    implicit none

    integer i, N
    integer itime
    integer, parameter :: iunit = 20
    real*8 time
    real*8 fj(*)
    character(80) fname

    itime = INT(time*100000)
    if (MOD(itime, 10) == 0) then
        !print *, N, time, itime, MOD(itime, 10)
        write(fname,'("lhcd/distribution/dstr", f9.7,".dat")') time
        !print *, fname
        open(iunit,file=fname,position="append")
        do i=1,n
            write(iunit,*) i, fj(i)
        end do
        close(iunit)
    end if
end subroutine


