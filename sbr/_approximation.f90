module approximation
    !+ polinomial approximation
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64      
    implicit none
    
contains
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function polin(k,x)
    implicit none
    integer k
    real(dp) x
    polin=1d0
    if(k.gt.1) polin=x**(k-1)
    return
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function polin1(k,x)
    implicit none
    integer k
    real(dp) x
    polin1=x**k
    return
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function polin2(k,x)
    implicit none
    integer k
    real(dp) x
    polin2=x**(k+1)
    return
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine approx(x,y,n,f,m,b)
!
!     y(i)=y(x(i))  the data to be approximated
!     n             number of points in the input data
!     m             number of coefficients of decomposition
!                   over base functions "f(k,x)" :
!                          y(x)=sum_1^m [b(k)*f(k,x)]
!     b(i)          found decomposition coefficients
!
    implicit real*8 (a-h,o-z)
    integer,  parameter :: np=20
    real(dp), parameter :: zero=0.d0
    real(dp) a(np,np),indx(np)
    real(dp) y(n),x(n),b(*)
    integer i,j,k,m,n
    if(m.gt.np) then
        write(*,*)'index error subroutine "approx"'
        return
    end if

    do j=1,m
        do k=1,j
            a(k,j)=zero
            do i=1,n
                a(k,j)=a(k,j)+f(j,x(i))*f(k,x(i))
            end do
        end do
    end do
    do k=2,m
        do j=1,k-1
            a(k,j)=a(j,k)
        end do
    end do

    do k=1,m
     b(k)=zero
      do i=1,n
       b(k)=b(k)+y(i)*f(k,x(i))
      end do
    end do

    call ludcmp(a,m,np,indx,d)
    call lubksb(a,m,np,indx,b)
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ludcmp(a,n,np,indx,d)
    implicit real*8 (a-h,o-z)
    integer,  parameter :: nmax=501
    real(dp), parameter :: tiny=1.d-20, zero=0.d0
    real(dp) a(np,np),indx(n),vv(nmax)
    integer i,j,k,m,n,np,imax
    d=1.d0
    do 12 i=1,n
    aamax=zero
    do 11 j=1,n
    if (dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
11      continue
    if (aamax.eq.zero) pause 'singular matrix.'
    vv(i)=1.d0/aamax
12    continue
    do 19 j=1,n
    if (j.gt.1) then
    do 14 i=1,j-1
      sum=a(i,j)
      if (i.gt.1)then
        do 13 k=1,i-1
          sum=sum-a(i,k)*a(k,j)
13            continue
        a(i,j)=sum
      endif
14        continue
    endif
    aamax=zero
    do 16 i=j,n
        sum=a(i,j)
        if (j.gt.1)then
            do 15 k=1,j-1
                sum=sum-a(i,k)*a(k,j)
15          continue
    a(i,j)=sum
    endif
    dum=vv(i)*dabs(sum)
    if (dum.ge.aamax) then
      imax=i
      aamax=dum
    endif
16      continue
    if (j.ne.imax)then
    do 17 k=1,n
      dum=a(imax,k)
      a(imax,k)=a(j,k)
      a(j,k)=dum
17        continue
    d=-d
    vv(imax)=vv(j)
    endif
    indx(j)=imax
    if(j.ne.n) then
    if(a(j,j).eq.zero) a(j,j)=tiny
    dum=1.d0/a(j,j)
    do 18 i=j+1,n
      a(i,j)=a(i,j)*dum
18        continue
    endif
19    continue
    if(a(n,n).eq.zero) a(n,n)=tiny
    return
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine lubksb(a,n,np,indx,b)
    implicit real*8 (a-h,o-z)
    real(dp), parameter :: zero=0.d0
    real(dp)  a(np,np),indx(n),b(n)
    integer i,j,ii,ll,n,np 
    ii=0
    do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
            do j=ii,i-1
                sum=sum-a(i,j)*b(j)
            end do
        else if (sum.ne.zero) then
            ii=i
        endif
        b(i)=sum
    end do
    do i=n,1,-1
        sum=b(i)
        if(i.lt.n)then
            do j=i+1,n
                sum=sum-a(i,j)*b(j)
            end do
        endif
        b(i)=sum/a(i,i)
    end do
    return
end    
end module approximation