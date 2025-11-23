    module model
        use fluid_eos
        implicit none
        real, parameter :: pi = 3.14159265358979323846
        real :: wcharge, hdepth
        real :: gaccel  ! gravitational acceleration
        real :: patm    ! atmospheric pressure
        real :: p0, rho0, c0, e0  ! initial pressure, density, sound speed, specific internal energy of water
        real :: alpha    ! reflection coefficient
        real :: rp, rz, dist  ! position of the pressure sensor from the initial bubble center
        type(eos) :: water, explosive


        ! parameters for shock wave and bubble
        real :: Rs, dRs ! shock wave radius
        real :: Rb0     ! initial bubble radius
        real :: Rb, dRb ! bubble radius
        real :: L, dL   ! computational domain length
        real :: pg      ! pressure at the left boundary
        real :: pg0     ! initial pressure inside the bubble
        real :: Ps      ! shock wave pressure
        real :: rhos     ! shock wave density
        real :: Us      ! shock wave material velocity
        real :: cs      ! shock wave sound speed


        ! global data
        real :: t, dt, tend, CFL
        real :: t_arrive  ! time when shock wave arrives at the sensor
        real :: t_far     ! time when the far-field simulation starts
        real :: t_bubble  ! time when bubble dynamics starts
        real :: t_shock_end ! time when shock wave passes the sensor
        real, parameter:: t_offset = 1e-5 ! time offset for sensor data


        ! used to store the influence of shock wave to bubble motion
        real, allocatable:: shock2bubble(:,:)
        integer :: n_s2b


        contains

        subroutine update_shock_front()
            ! calculate shock front velocity and material velocity
            ! based on Hugoniout conditions
            ! assuming Ps is known
            implicit none
            real :: M
            dRs = (Ps+water%Pw)*(water%gamma+1.0)/(P0+water%Pw) + &
                water%gamma-1.0
            M = sqrt(dRs/(2.0*water%gamma))
            dRs = M*c0
            us = 2*c0/(water%gamma+1.0)*(M**2-1.0)/M
            rhos = rho0*(water%gamma+1.0)*M**2/( &
                (water%gamma-1.0)*M**2+2.0)
            cs = water%sound_spd(rhos, Ps)
            return
        end subroutine update_shock_front
    end module