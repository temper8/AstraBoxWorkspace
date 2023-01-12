module rt_parameters
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none
!   physical parameters 
    real(dp) :: freq
    !+ Freq,     RF frequency, GHz
    real(dp) :: xmi1
    !+  Mi1/Mp,  relative mass of ions 1
    real(dp) :: zi1
    !+ charge of ions 1
    real(dp) :: xmi2
    !+ Mi2/Mp,  relative mass of ions 2
    real(dp) :: zi2
    !+ charge of ions 2
    real(dp) :: dni2 
    !+  0.03   Ni2/Ni1, relative density of ions 2
    real(dp) :: xmi3
    !+  Mi3/Mp,  relative mass of ions 3
    real(dp) :: zi3
    !+  charge of ions 3
    real(dp) :: dni3
    !+  Ni3/Ni1, relative density of ions 3

!!!!!!!!!!!!!  parameters for alphas calculations !!!
    integer  ::  itend0
    !+ itend0,   if = 0, no alphas
    real(dp) ::  energy
    !+ energy,   max. perp. energy of alphas (MeV)
    real(dp) ::  factor   
    !+ factor,   factor in alpha source
    real(dp) ::  dra   
    !+ dra,      relative alpha source broadening (dr/a)
    integer  ::  kv    
    !+ kv,       V_perp  greed number    

!!!!!!!!!!!!! numerical parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer  ::  nr     
    !+ nr,  radial grid number  <= 505
    real(dp) ::  hmin1
    !+ hmin1, rel.(hr) min. step in the Fast comp. mode, <1.d0
    real(dp) ::  rrange
    !+ rrange,   rel.(hr) size of a 'turning' point region, <1.d0
    real(dp) ::  eps
    !+ eps,      accuracy
    real(dp) ::  hdrob
    !+ hdrob,    h4 correction,
    real(dp) ::  cleft
    !+ cleft,    left Vz plato border shift (<1)
    real(dp) ::  cright
    !+ cright,   right Vz plato border shift (>1)
    real(dp) ::  cdel
    !+ cdel,     (left part)/(Vz plato size)
    real(dp) ::  rbord 
    !+ rbord,    relative radius of reflection, <1.
    real(dp) ::  pchm
    !+ pchm,     threshold between 'strong' and weak' absorption, <1.
    real(dp) ::  pabs0
    !+ pabs,     part of remaining power interp. as absorption
    real(dp) ::  pgiter
    !+ pgiter,   relative accuracy to stop iterations
    integer ::   ni1     
    !+ ni1,      grid number in the left part of Vz plato
    integer ::   ni2     
    !+ ni2,      grid number in the right part of Vz plato
    integer ::   niterat     
    !+ niterat,  maximal number of iterations
    integer ::   nmaxm(4)
    !+ nmaxm(1), permitted reflections at 0 iteration
    !+ nmaxm(2), permitted reflections at 1 iteration
    !+ nmaxm(3), permitted reflections at 2 iteration
    !+ nmaxm(4), permitted reflections at 3 iteration
    integer ::   maxstep2  
    !+ maxstep2, maximal steps' number in Fast comp. mode
    integer ::   maxstep4    
    !+ maxstep4, maximal steps' number in Slow comp. mode    

!!!!!!!!!!!!!  options !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer ::  ipri
    !+ ipri, printing output monitoring: 0,1,2,3,4
    integer ::  iw
    !+ iw, initial mode (slow=1, fast=-1)
    integer ::  ismth
    !+ ismth, if=0, no smoothing in Ne(rho),Te(rho),Ti(rho)
    integer ::  ismthalf
    !+ ismthalf,  if=0, no smoothing in D_alpha(vperp)
    integer ::  ismthout
    !+ ismthout,  if=0, no smoothing in output profiles
    integer ::  inew
    !+ inew=0 for usual tokamak&Ntor_grill; 1 or 2 for g' in ST&Npol_grill
    integer ::  itor
    !+ itor,      +-1, Btor direction in right coord{drho,dteta,dfi}
    integer ::  ipol
    !+ ipol,      +-1, Bpol direction in right coord{drho,dteta,dfi}

  !!!!!!!!!!!!!  grill parameters and input LH spectrum !!!!!!!!!!!!
    real(dp) ::  zplus 
    !+ Zplus,    upper grill corner in centimeters
    real(dp) ::  zminus
    !+ Zminus,   lower grill corner in centimeters
    integer ::   ntet
    !+ ntet,     theta grid number
    integer ::   nnz
    !+ nnz,      N_phi grid number    

    contains      
    subroutine show_parameters()          
      print*, "Freq = ", freq          
      print*, "xmi1 = ", xmi1     
      print*, "zi1 = ",  zi1     
      print*, "xmi2 = ", xmi2     
      print*, "zi2 = ",  zi2     
      print*, "dni2 = ", dni2     

      print*, "---------- grill parameters and input LH spectrum "
      print*, "zplus = ",  zplus     
      print*, "zminus = ", zminus     
      print*, "ntet = ",  ntet     
      print*, "nnz = ", nnz           
    end subroutine show_parameters
    
    subroutine read_parameters(file_name)
        implicit none
        integer, parameter :: iunit = 20
        character(*) file_name
        print *, file_name

        open(iunit, file= file_name)
        !!!!!!!!!!!!!  read  physical parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               read(iunit,*)
               read(iunit,*) freq
               read(iunit,*) xmi1
               read(iunit,*) zi1
               read(iunit,*) xmi2
               read(iunit,*) zi2
               read(iunit,*) dni2
               read(iunit,*) xmi3
               read(iunit,*) zi3
               read(iunit,*) dni3
        !!!!!!!!!!!!!  read parameters for alphas calculation !!!!!!!!!!!!!!!!!!!
               read(iunit,*)
               read(iunit,*) itend0
               read(iunit,*) energy
               read(iunit,*) factor
               read(iunit,*) dra
               read(iunit,*) kv
        
        !!!!!!!!!!!!!  read  numerical parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               read(iunit,*)
               read(iunit,*) nr
               read(iunit,*) hmin1
               read(iunit,*) rrange
               read(iunit,*) eps
               read(iunit,*) hdrob
               read(iunit,*) cleft
               read(iunit,*) cright
               read(iunit,*) cdel
               read(iunit,*) rbord
               read(iunit,*) pchm
               read(iunit,*) pabs0
               read(iunit,*) pgiter
               read(iunit,*) ni1
               read(iunit,*) ni2
               read(iunit,*) niterat
               read(iunit,*) nmaxm(1)
               read(iunit,*) nmaxm(2)
               read(iunit,*) nmaxm(3)
               read(iunit,*) nmaxm(4)
               read(iunit,*) maxstep2
               read(iunit,*) maxstep4
        
        !!!!!!!!!!!!!  read  options !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               read(iunit,*)
               read(iunit,*) ipri
               read(iunit,*) iw
               read(iunit,*) ismth
               read(iunit,*) ismthalf
               read(iunit,*) ismthout
               read(iunit,*) inew
        
               read(iunit,*) itor     !Btor direction in right-hand {drho,dteta,dfi}
               read(iunit,*) ipol     !Bpol direction in right-hand {drho,dteta,dfi}
        
        !!!!!!!!!!!!!  read grill parameters and input LH spectrum !!!!!!!!!!!!
               read(iunit,*)
               read(iunit,*) zplus
               read(iunit,*) zminus
               read(iunit,*) ntet
               read(iunit,*) nnz
               read(iunit,*)
        close(iunit)        
        call show_parameters
    end subroutine read_parameters
end module rt_parameters