    module global
    logical imigration
    logical ibound
    real gamma,c,g,rho,cd,ca,Pv
    real miu, sigma
    real pi
    real dt,t,tend,sdt
    real tdelay
    
    real pamb,uamb(3)
    
    real pressure_loc(3)
    !integer order
    integer inc
    character*1000 path
    
    real t_arrive   ! when the shock wave arrive at the pressure points
    real t_offset   ! time offset for bubble pressure
    real t_start    ! start time of bubble simulation
    real t_shock    ! end time of the shock wave
    logical iarrive ! whether shock wave arrive at the pressure points
    
    end module