    subroutine initialize()
      use model
      use DG
      use near_field
      use bubble
      implicit none
      real :: mat_inp(10)

      ! print banner
      print*, '============================================'
      print*, '*       Undex-UBE Simulation Program       *'
      print*, '*         Developed by Yunlong Liu         *'
      print*, '*     in Harbin Engineering University     *'
      print*, '*       yunlong_liu@hrbeu.edu.cn           *'
      print*, '============================================'
      print*, ''
      print*, ''

      print*, '***        Initialization started        ***'



      open(unit=10, file='./input/case.inp',&
         status='old', action='read')
      read(10,*) wcharge
      read(10,*) hdepth
      read(10,*) tend
      read(10,*) rp
      read(10,*) rz
      read(10,*) gaccel
      read(10,*) patm
      read(10,*) alpha
      read(10,*) Ca
      read(10,*) Cd
      read(10,*) CFL

      close(10)
      
      dist = sqrt(rp**2 + rz**2)
      ! Initialize equations of state
      open(unit=10, file='./input/water.mat', &
        status='old', action='read')
      read(10,*) mat_inp(1)
      read(10,*) mat_inp(2)
      read(10,*) mat_inp(3)
      close(10)
      call water%init(1, mat_inp(1:3))
      open(unit=10, file='./input/jwl.mat', &
        status='old', action='read')
      read(10,*) mat_inp(1)
      read(10,*) mat_inp(2)
      read(10,*) mat_inp(3) 
      read(10,*) mat_inp(4)
      read(10,*) mat_inp(5)
      read(10,*) mat_inp(6)
      read(10,*) mat_inp(7)
      close(10)
      call explosive%init(2, mat_inp(1:8))

      ! calculate stastic pressure, sound speed
      rho0 = water%dens0
      p0 = patm + rho0 * gaccel * hdepth
      c0 = water%sound_spd(rho0, p0)
      e0 = water%internal_energy(p0, rho0)

      ! initialize time
      t = 0.0

      ! initialize DG basis
      call init_basis()

      ! initialize near-field mesh
      Rb0 = (3.0/(4.0*pi)*wcharge/explosive%dens0)**(1.0/3.0)
      Pg0 = explosive%press(explosive%dens0, explosive%e0)
      Pg = Pg0
      Ps = Pg
      call update_shock_front()
      ! initialize near-field condition
      call init_condi_near_field(50, 3)
      Rb = Rb0
      Rs = Rb*1.0001
      dRb = Us
      L = Rs-Rb
      dL = dRs - dRb


      ! open file for output
      open(200, file='output/sensor.dat', status='replace', action='write')
      write(200,"(2A15)") '# time', 'P(Pa)'

      open(201, file='output/bubble.dat', status='replace', action='write')
      write(201,"(4A15)") '# time', 'Rb(m)', 'r(m)', 'z(m)'

      !open(202, file='output/shockwave.dat', status='replace', action='write')
      !write(202,"(3A15)") '# time', 'Rs(m)', 'Pm(Pa)'


      allocate(shock2bubble(2,100000))
        shock2bubble = 0
        n_s2b = 0

      end subroutine initialize