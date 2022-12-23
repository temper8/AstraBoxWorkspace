module FluxSurface_mod
    !+ все что связанно с магнитными поверхностями
    type FluxSurface 
    !+ класс магнитной поверхности
        integer       :: index
        !+ номер магнитной поверхности
        real (real64) :: r
        !+ радиус
        real (real64) :: vmax
        !+ vmax=cltn/vto
        real (real64) :: vt
        !+ наверно тепловая скорость электронов????? vt=fvt(r)
        integer       :: ipt
        !+ размер vgrid
        real (real64) :: vgrid(:)
        !+ 
        real (real64) :: vr_grid(:)
        !+ бываший vrj
        real (real64) :: diffusion(:)
        !+ бывший dijk(i,j,k) или dj(i)
        !   complex         :: inst_field1
        contains
        procedure :: set   => set_e
        procedure :: print => e_print

    end type FluxSurface  

end module FluxSurface_mod