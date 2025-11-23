    ! module holding the global variables
    module global
    logical imigration              ! whether to consider the migration
    logical ibound                  ! whether to consider the boundary
    real gamma,c,g,rho,cd,ca,Pv     ! constants
    real miu, sigma                 ! viscosity and surface tension
    real pi                         
    real dt,t,tend,sdt              ! time increment, current time, end time
    real pamb,uamb(3)               ! ambient pressure and velocity
    real pressure_loc(3)            ! location where the pressure is to be output
    integer inc                     ! increment number
    character*1000 path             
    end module