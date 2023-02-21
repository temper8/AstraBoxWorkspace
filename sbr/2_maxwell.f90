module maxwell
      use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64          
      use constants, only : zero
      implicit none
      integer, parameter :: i0 = 1002

      real(dp) vij(i0,100), fij0(i0,100,2), fij(i0,100,2)
      real(dp) dfij(i0,100,2), dij(i0,100,2)

contains
      subroutine init_vi(vclt, vi)
            real*8, intent(in) :: vclt
            real*8, intent(out) :: vi(i0)            
            real*8 vmax
            integer i
            vmax = 2.d0*vclt
            do i=1,i0
                  vi(i)=dble(i-1)*vmax/dble(i0-1)
            end do
      end subroutine      

      subroutine init_fmaxw_classic(vclt, enorm, fi, dfi)
            real*8, intent(in) :: vclt, enorm
            real*8, intent(out) :: fi(i0), dfi(i0)
            real*8 vi, vmax
            integer i
            vmax = 2.d0*vclt
            do i=1,i0
                  vi = dble(i-1)*vmax/dble(i0-1)
                  if(vi < vclt) then
                        fi(i) = fmaxw_classic(vi, enorm, dfi(i))
                  else
                        fi(i) = zero
                        dfi(i) = zero
                  end if
            end do
      end subroutine      

      subroutine init_fmaxw_ext(vclt, enorm, fi, dfi)
            real*8, intent(in) :: vclt, enorm
            real*8, intent(out) :: fi(i0), dfi(i0)
            real*8 vi, vmax
            integer i
            vmax=2.d0*vclt
            do i=1,i0
                  vi = dble(i-1)*vmax/dble(i0-1)
                  if(vi < vclt) then
                        fi(i) = fmaxw_ext(vi, enorm, dfi(i))
                  else
                        fi(i) = zero
                        dfi(i) = zero
                  end if
            end do
      end subroutine      

      double precision function funmaxwell(v,dfunmaxwell)
            implicit none
            real*8 v,dfunmaxwell,arg,pi2sqrt
            parameter(pi2sqrt=2.506628274631d0)

            arg=-0.5d0*v**2
            funmaxwell=dexp(arg)/pi2sqrt
            dfunmaxwell=-v*funmaxwell
      end      

      double precision function fmaxw_classic(v,alfa2,dfmaxw)
            implicit none
            real*8 v,alfa2,dfmaxw
            real*8 arg,alfa,api,b,psiq,f,df
            real*8 pi2sqrt,pisqrt,zero
            parameter(pi2sqrt=2.506628274631d0,pisqrt=1.77245385090552d0)
            parameter(zero=0.d0)
            
            arg=-0.5d0*v**2*(1.d0+0.5d0*alfa2*v**2)
            fmaxw_classic=dexp(arg)/pi2sqrt
            dfmaxw=-v*(1.d0+alfa2*v**2)*fmaxw_classic
      end

      double precision function fmaxw_ext(v,alfa2,dfmaxw)
            implicit none
            real*8 v,alfa2,dfmaxw
            real*8 arg,alfa,api,b,f,df
            real*8 pi2sqrt,pisqrt,zero
            parameter(pi2sqrt=2.506628274631d0,pisqrt=1.77245385090552d0)
            parameter(zero=0.d0)

            alfa=dsqrt(alfa2)
            api=2.d0*alfa*dexp(-0.25d0/alfa2)/pisqrt
            b=2.d0-erfcc(0.5d0/alfa)+api
            f=psiq(v,alfa2)
            fmaxw_ext=(f+api)/b/pi2sqrt
            df=-v*((1.d0-alfa2*v**2)*f+api)
            dfmaxw=df/b/pi2sqrt
      end      

      double precision function fmaxw(v,alfa2,dfmaxw)
            implicit none
            real*8 v,alfa2,dfmaxw
            real*8 arg,alfa,api,b,f,df
            real*8 pi2sqrt,pisqrt,zero
            parameter(pi2sqrt=2.506628274631d0,pisqrt=1.77245385090552d0)
            parameter(zero=0.d0)
            if(alfa2.le.zero) then
                  arg=-0.5d0*v**2*(1.d0-0.5d0*alfa2*v**2)
                  fmaxw=dexp(arg)/pi2sqrt
                  dfmaxw=-v*(1.d0-alfa2*v**2)*fmaxw
            else
                  alfa=dsqrt(alfa2)
                  api=2.d0*alfa*dexp(-0.25d0/alfa2)/pisqrt
                  b=2.d0-erfcc(0.5d0/alfa)+api
                  f=psiq(v,alfa2)
                  fmaxw=(f+api)/b/pi2sqrt
                  df=-v*((1.d0-alfa2*v**2)*f+api)
                  dfmaxw=df/b/pi2sqrt
            end if
      end      

      double precision function psiq(v,alfa2)
!!! psiq=exp(ksiV**2)*erfcc(ksiV)*exp(-0.25/alfa2)
            implicit none
            double precision v,alfa2,df
            double precision x,t,z,f,asymp,alfa,q,u
            double precision zero,zmax,pisqrt
            parameter(zero=0.d0,zmax=10.d0,pisqrt=1.77245385090552d0)
            
            alfa=dsqrt(alfa2)
            q=-0.25d0/alfa2
            x=0.5d0*(alfa*v**2-1.d0/alfa)
            z=abs(x)
            if(z.gt.zmax) then !asymptotics
                  f=dexp(q)*(1.d0-0.5d0/z**2+0.75d0/z**4-15.d0/8.d0/z**6)/z/pisqrt
            else
                  t=1.d0/(1.d0+0.5d0*z)
                  f=t*exp(q-1.26551223d0+t*(1.00002368d0+t*(.37409196d0+t*&
                  &(.09678418d0+t*(-.18628806d0+t*(.27886807d0+t*(-1.13520398d0+t*&
                  &(1.48851587d0+t*(-.82215223d0+t*.17087277d0)))))))))
            end if
            if(x.lt.zero) then
                  u=-0.5d0*v**2+0.25d0*alfa2*v**4 !u=x**2-0.25d0/alfa2
                  f=2.d0*dexp(u)-f
            end if
            psiq=f
            return
            end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function erfcc(x)
            implicit none
            double precision erfcc,x
            double precision t,z,zero,zmax,pisqrt
            parameter(zero=0.d0,zmax=10.d0,pisqrt=1.77245385090552d0)
            
            z=abs(x)
            if(z.gt.zmax) then !asymptotics
                  erfcc=(1.d0-0.5d0/z**2+0.75d0/z**4-15.d0/8.d0/z**6)/z/pisqrt
                  erfcc=exp(-z*z)*erfcc
            else
                  t=1.d0/(1.d0+0.5d0*z)
                  erfcc=t*exp(-z*z-1.26551223d0+t*(1.00002368d0+t*(.37409196d0+t*&
                  &(.09678418d0+t*(-.18628806d0+t*(.27886807d0+t*(-1.13520398d0+t*&
                  &(1.48851587d0+t*(-.82215223d0+t*.17087277d0)))))))))
            end if
            if(x.lt.zero) erfcc=2.d0-erfcc
            return
            end      
end module maxwell
