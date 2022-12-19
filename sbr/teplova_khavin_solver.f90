subroutine teplova_khavin_solver(alfa2, nt, h, dt, n, ybeg, yend, d1,d2,d3, y)
    ! схема Ченга-Купера для уравнения Фоккера-Планка
    implicit none
    real*8, intent(in)  :: alfa2      
    integer, intent(in) :: nt, n
    real*8, intent(in)  :: h, dt
    real*8, intent(in)  :: ybeg, yend
    real*8, intent(in)  :: d1(n+1),d2(n+1),d3(n+1)
    real*8, intent(inout) :: y(n), y1(n)
    integer i, it
    real*8 xx(n+1), a(n),b(n),c(n),f(n)

    do i=1,n+1
          xx(i)=h/2.d0+h*dble(i-1) !+shift
    end do
    
    do i=1,n
        y1(i)=y(i+1)
    end do
    
    do it=1,nt
           !call ABCcoef(a,b,c,f,y1,dt,n,ybeg,yend,x,xx,h,D1)
        call TK_abcoef(alfa2, a,b,c,f,y1, dt, n, ybeg, yend, xx, h, d1)
        call tridag(a,b,c,f,y1,n)
                  
          !call tridag(a,b,c,f,y,n)   
        do i=1,n
            if (y1(i).lt.0.d0) then
                if (y1(i) > epsilon(y1(i))) then
                    y1(i)=0.d0
                else
                    write(*,*) 'y(i)=',y1(i),' lt negative epsilon=',epsilon(y1(i))
                    pause
                    stop
                endif
            endif
        enddo
    end do

    do i=1,n
        y(i+1)=y1(i)
    end do

end subroutine

! --
subroutine TK_abcoef(alfa2, A,B,C,f,Y,k,n,ybeg,yend,xx,h,df)
    implicit none
    real*8, intent(in)    :: alfa2
    real*8, intent(inout) :: a(n),b(n),c(n),f(n),y(n)
    real*8, intent(in)    :: dt
    integer, intent(in)   :: n
    real*8, intent(in)    :: ybeg, yend, h
    real*8, intent(in)    :: xx(n+1)
    real*8, intent(in)    :: df(n+1)

    integer i
    real*8 k,df(n)
    real*8 C1,B1,z,w,r,dlt,alfa2
    real*8 tmp1,tmp2,tmp3
    external C1,w,B1,dlt
    r=k/h
    do i=1,n
    
        tmp1=dlt(xx(i),h,df(i)) * B1(xx(i),alfa2)
        A(i)=-r*(C1(xx(i),df(i))/h-tmp1)

        tmp2=C1(xx(i+1),df(i+1))/h-dlt(xx(i+1),h,df(i+1), alfa2) * B1(xx(i+1),alfa2)
        tmp3=(1.d0-dlt(xx(i),h,df(i), alfa2)) * B1(xx(i),alfa2)
        B(i)=r*(tmp2+tmp3+C1(xx(i),df(i))/h)+1.d0
  
        tmp1=(1.d0-dlt(xx(i+1),h,df(i+1), alfa2)) * B1(xx(i+1),alfa2)
        C(i)=-r*(tmp1+C1(xx(i+1),df(i+1))/h)

        f(i)=Y(i)
    enddo
    f(1)=f(1)-A(1)*ybeg
    f(n)=f(n)-C(n)*yend !yend in either way=0 all the time
end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

real*8 function B1(xx, alfa2)
    implicit none
    real*8 xx,alfa2,beta
    B1=-alfa2+1.d0/(xx*xx)
end function

      
real*8 function C1(xx,dif)
    implicit none
    real*8 xx,dif
    C1=dif+1.d0/(xx*xx*xx)
end function
      
      
!real*8 function w(xx,h,dif)
!    implicit none
!    real*8 xx,h,B1,C1,dif
!    w=h*B1(xx)/C1(xx,dif)
!end function
      
      
function dlt(xx,h,dif, alfa2) result(res)
    implicit none
    real*8 res
    real*8 xx, h, dif, alfa2
    real*8 B1, C1, w
    w = h*B1(xx, alfa2)/C1(xx,dif)
    res = 1.d0/w-1.d0/(dexp(w)-1.d0)
end function
