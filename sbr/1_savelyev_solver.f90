module savelyev_solver_module
    use kind_module
    implicit none
    
    PRIVATE :: q, k
contains
    subroutine savelyev_solver(alfa2, nt, h, dt, n, ybeg, yend, d1,d2,d3, y)
        ! разностная схема Савельева для уравнения Фоккера-Планка
        implicit none
        real(wp), intent(in)  :: alfa2      
        integer, intent(in) :: nt, n
        real(wp), intent(in)  :: h, dt
        real(wp), intent(in)  :: ybeg, yend
        real(wp), intent(in)  :: d1(n+1),d2(n+1),d3(n+1)
        real(wp), intent(inout) :: y(n)
        integer i, it
        real(wp) xx(n+1), a(n),b(n),c(n),f(n)

        do i=1,n+1
            xx(i)=h/2.d0+h*dble(i-1) !+shift
        end do
    
        do it=1,nt
            call savelyev_abccoef(alfa2, a,b,c,f,y,dt,n,ybeg,yend,xx,h,d1,d2,d3)
            call tridag(a,b,c,f,y,n)      
        end do

    end subroutine

    !!!!!!! -- fill abc matrix
    subroutine savelyev_abccoef(alfa2, a,b,c, f, y, dt, n, ybeg, yend, xx, h, d1,d2,d3)
        implicit none
        real(wp), intent(in)    :: alfa2
        real(wp), intent(inout) :: a(n),b(n),c(n),f(n),y(n)
        real(wp), intent(in)    :: dt
        integer, intent(in)   :: n
        real(wp), intent(in)    :: ybeg, yend, h
        real(wp), intent(in)    :: xx(n+1)
        real(wp), intent(in)    :: d1(n+1),d2(n+1),d3(n+1)

        integer i,iunit,iunit2
        real(wp) a1(n),b1(n),c1(n),f1(n),a2(n),b2(n),c2(n),f2(n)
        real(wp) kinv,rs,rmink,rplusk,q,qf,r1,rmink2,rplusk2,kinv2
        real(wp) r,kappa,sum,bmin,bplus,sum2,sum3,sum4
        real(wp) dc,as(n+1),k,k2,d
        external kinv,rs,rmink,rplusk,q,kinv2,rmink2,rplusk2,d

        sum=(kinv(xx(1) - h/2d0, d2(1), alfa2) + kinv(xx(1) + h/2d0, d3(1), alfa2))*h/2d0
        as(1)=h/sum

        sum=(kinv(xx(2)-h/2d0, d2(2), alfa2)+kinv(xx(2)+h/2d0, d3(2), alfa2))*h/2d0
        as(2)=h/sum

        r=h/2d0*dabs(rs(xx(1)+h/2d0, d3(1), alfa2))/k(xx(1)+h/2d0, d3(1))
        kappa=1d0/(1d0+r)
        sum=(rmink(xx(1), d1(1), alfa2) + rmink(xx(2), d1(2), alfa2))*h/2d0
        bmin=sum/h

        sum = (rplusk(xx(1), d1(1), alfa2) + rplusk(xx(2), d1(2), alfa2))*h/2d0
        bplus = sum/h

        sum = qf(xx(2))-qf(xx(1))
        dc = sum/h

        a(1) = as(1)*(kappa/h**2 - bmin/h)
        c(1) = as(2)*(kappa/h**2 + bplus/h)
        b(1) = -(1d0/dt+a(1) + c(1) + dc)
        f(1) = -y(1)/dt-a(1)*ybeg
        do i=2,n
            sum = (kinv(xx(i+1) - h/2d0, d2(i+1)) + kinv(xx(i+1) + h/2d0, d3(i+1)))
            sum = sum*h/2d0
            as(i+1) = h/sum

            r = h/2d0*dabs(rs(xx(i) + h/2d0, alfa2))/k(xx(i) + h/2d0, d3(i))
            kappa=1d0/(1d0+r)
            sum = (rmink(xx(i), d1(i), alfa2) + rmink(xx(i+1), d1(i+1), alfa2))*h/2d0
            bmin = sum/h
            sum = (rplusk(xx(i), d1(i), alfa2) + rplusk(xx(i+1), d1(i+1), alfa2))*h/2d0
            bplus = sum/h
            sum = qf(xx(i+1)) - qf(xx(i))
            dc = sum/h
            
            a(i) = as(i)*(kappa/h**2 - bmin/h)
            c(i) = as(i+1)*(kappa/h**2 + bplus/h) 
            b(i) = -(1d0/dt + a(i) + c(i) + dc) 
            f(i) = -y(i)/dt
        end do
        f(n)=f(n)-c(n)*yend
        a(1)=0d0
        c(n)=0d0
    end

    real(wp) function rplusk(x,dif, alfa2)
        implicit none
        integer iunit
        real(wp) x,k,rs,dif,d,razn
        real(wp) alfa2      
        rplusk=0.5d0*(rs(x, alfa2)+dabs(rs(x, alfa2)))/k(x,dif)
    end

    real(wp) function rplusk2(x,dif, alfa2)
        implicit none
        integer iunit
        real(wp) x,k2,rs,dif,d,razn
        real(wp) alfa2      
        rplusk2=0.5d0*(rs(x, alfa2)+dabs(rs(x, alfa2)))/k2(x,dif)
    end

    real(wp) function rmink(x,dif, alfa2)
        implicit none
        integer iunit
        real(wp) x,k,rs,dif,d,razn
        real(wp) alfa2      
        rmink=0.5d0*(rs(x, alfa2)-dabs(rs(x, alfa2)))/k(x,dif)
    end

    real(wp) function rmink2(x,dif, alfa2)
        implicit none
        integer iunit
        real(wp) x,k2,rs,dif,d,razn
        real(wp) alfa2      
        rmink2=0.5d0*(rs(x, alfa2)-dabs(rs(x, alfa2)))/k2(x,dif)
    end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(wp) function rs(x, alfa2)
        implicit none
        real(wp) x
        real(wp) alfa2
        !common/ef/ alfa2
        rs=1d0/x**2-alfa2
    end

    real(wp) function q(x)
        implicit none
        real(wp) x
        q=2d0/x**3
    end

    real(wp) function qf(x)
        implicit none
        real(wp) x
        qf=-1d0/x**2
    end

    real(wp) function k(x,dif)
        implicit none
        integer iunit
        real(wp) x,dif,d,razn
        k=dif+1d0/x**3
    end

    real(wp) function k2(x,dif)
        implicit none
        integer iunit
        real(wp) x,dif,d,razn
        k2=d(x)+1d0/x**3
    end

    real(wp) function kinv(x,dif)
        implicit none
        integer iunit
        real(wp) x,dif,razn,d
        kinv=x**3/(dif*x**3+1d0)
    end

    real(wp) function kinv2(x,dif)
        implicit none
        integer iunit
        real(wp) x,dif,razn,d,kino
        kinv2=x**3/(d(x)*x**3+1d0)
    end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine tridag(a,b,c,r,u,n)
        implicit none
        integer, intent(in)    :: n
        real(wp),  intent(in)    :: a(n), b(n), c(n), r(n)
        real(wp),  intent(inout) :: u(n)
        integer, parameter :: nmax=1000000
        integer j
        real(wp) bet, gam(nmax)

        if(b(1).eq.0.d0) pause 'tridag: rewrite equations'
        bet=b(1)
        u(1)=r(1)/bet
        do j=2,n
            gam(j)=c(j-1)/bet
            bet=b(j)-a(j)*gam(j)
            if(bet.eq.0.d0) then
                    write(*,*)'b(j)=',b(j),'a(j)=',a(j),'gam(j)=',gam(j)
                    pause 'tridag failed'
            end if
            u(j)=(r(j)-a(j)*u(j-1))/bet
        end do
        do j=n-1,1,-1
            u(j)=u(j)-gam(j+1)*u(j+1)
        end do
    end subroutine
end module savelyev_solver_module