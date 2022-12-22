module FluxSurface_mod ! the module name defines the namespace
    type FluxSurface ! classname is the class prototype name
        integer        :: direction
        real*8         :: enorm

        !   complex         :: inst_field1
        contains
        procedure :: set   => set_e
        procedure :: print => e_print

    end type FluxSurface  

end module FluxSurface_mod