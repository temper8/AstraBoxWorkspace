 
subroutine write_array(arr, N, array_name)
    implicit none
    real*8, intent(in) :: arr(*)
    integer, intent(in) :: N
    character(len=*), intent(in) :: array_name
    
    integer i
    integer, parameter :: iunit = 21
    
    character(80) fname
    print *, 'write_array:', array_name, N
    write(fname,'("lhcd/distribution/", A,".dat")') array_name
    print *, fname
    open(iunit,file=fname,position="append")
        do i=1,n
            write(iunit,*) i, arr(i)
        end do
    close(iunit)
end subroutine

subroutine write_distribution(arr,N,time)
    implicit none
    real*8, intent(in) :: arr(*)
    integer, intent(in) :: N
    real*8, intent(in) :: time

    integer i
    integer itime
    integer, parameter :: iunit = 20
    character(80) fname

    itime = INT(time*100000)
    !print *, N, time, itime, MOD(itime, 10)
    if (MOD(itime, 10) == 0) then
        write(fname,'("lhcd/distribution/dstr", f9.7,".dat")') time
        !print *, fname
        open(iunit,file=fname,position="append")
        do i=1, n
            write(iunit,*) i, arr(i)
        end do
        close(iunit)
    end if
end subroutine


