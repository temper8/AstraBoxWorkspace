module FluxSurface_mod
    !+ все что связанно с магнитными поверхностями
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    type FluxSurface 
    !+ класс магнитной поверхности
        integer       :: index
        !+ номер магнитной поверхности
        real(dp) :: r
        !+ радиус
        real(dp) :: vmax
        !+ vmax=cltn/vto
        real(dp) :: vt
        !+ наверно тепловая скорость электронов????? vt=fvt(r)
        integer       :: ipt
        !+ размер vgrid
        real(dp), allocatable :: vgrid(:)
        !+ 
        real(dp), allocatable :: vr_grid(:)
        !+ бываший vrj
        real(dp), allocatable :: diffusion(:)
        !+ бывший dijk(i,j,k) или dj(i)
        !   complex         :: inst_field1
        contains
        !procedure :: set   => set_e
        !procedure :: print => e_print

    end type FluxSurface  

end module FluxSurface_mod