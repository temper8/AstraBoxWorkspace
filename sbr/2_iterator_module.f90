module iterator_mod
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64      
    implicit none
    
contains
    subroutine recalculate_f_for_a_new_mesh
        implicit none
        !-------------------------------------------
        !   recalculate f' for a new mesh
        !-------------------------------------------
        k=(3-ispectr)/2
        do j=1,nr
            r=hr*dble(j)
            vt=fvt(r)
            vto=vt/vt0
            if(iterat.gt.0) then
                v1=dmin1(vzmin(j),vz1(j))
                v2=dmax1(vzmax(j),vz2(j))
            else
                v1=vzmin(j)
                v2=vzmax(j)
            end if
            vmax=cltn/vto
            vp1=v1/vto
            vp2=v2/vto
            call gridvel(vp1,vp2,vmax,cdel,ni1,ni2,ipt1,kpt3,vrj)
            do i=1,i0
                vvj(i)=vij(i,j)
                vdfj(i)=dfij(i,j,k) !=dfundv(i,j)*vto**2
            end do
            do i=1,ipt
                call lock(vvj,i0,vrj(i),klo,khi,ierr)
                if(ierr.eq.1) then
            !!!         if(vrj(i).gt.vvj(i0)) exit
                    write(*,*)'lock error in new v-mesh'
                    write(*,*)'j=',j,' i0=',i0
                    write(*,*)'vvj(1)=',vvj(1),' vvj(i0)=',vvj(i0)
                    write(*,*)'i=',i,' vrj(i)=',vrj(i)
                    write(*,*)
                    pause'next key = stop'
                    stop
                end if
                call linf(vvj,vdfj,vrj(i),dfout,klo,khi)
                vgrid(i,j)=vrj(i)*vto
                dfundv(i,j)=dfout/vto**2
                if(dfundv(i,j).gt.zero) dfundv(i,j)=zero
            end do
            vz1(j)=v1
            vz2(j)=v2
        end do
    end subroutine
end module iterator_mod