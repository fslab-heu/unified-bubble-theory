    module shockwave
    
    ! variables definition
    
    ! global
    real pi
    real dt, dt0
    real time
    ! physical
    
    real Rb0, Rb, Rs, dRb, dRs, L, dL! gemoetry variables
    real Ma, us, ps, rhos, cs        ! shock front vars
    real p0,rho0,c0,pw,gamma         ! material of water
    real jwl_E,jwl_rho0,A,B,C,R1,R2,w      ! material of jwl
    
    real,parameter:: bias = 5.0
    integer,parameter:: ndg = 20
    real dRdg(ndg+1),Rdg(ndg+1),Ldg(ndg),dLdg(ndg)  ! 3 start points and lengths of the DG element
    real Rdg1(ndg+1)
    real flux(3,ndg+1)
    integer id_dg                       ! which element is currently 
    real,allocatable:: ans_dg(:,:,:),ans_dg1(:,:,:)
    ! backup variables for 2nd order marching
    real Rb1, Rs1, dRb1, dRs1, L1, dL1, Re1
    real,allocatable:: ans1(:,:), pans1(:)
    
    
    real pg                         ! internal gas pressure
    
    
    ! solution of the far-field wave equation
    real Re                         ! left end of the computational domain
    real,allocatable:: pans(:),prhs(:)
    real Lp,Kc,MKc
    
    !  numerical solution
    integer np, npmax
    real,allocatable:: mmatrix(:,:),imatrix(:,:), rhs(:,:),ans(:,:)
    ! dg solution
    !real,allocatable:: ans_dg(:,:,:)
    !real flux(2,3)
    
    ! saved varialbes
    integer s_nb,s_nb2,s_nb3
    
    
    contains
    include 'initialization.f90'
    !include 'shockwave_flux.f90'
    !include 'shockwave_buildmatrix.f90'
    include 'shockwave_march.f90'
    include 'basis.f90'
    include 'entropy.f90'
    include 'ini_pressure.f90'
    include 'far_field_march.f90'
    include 'universal_initial_condition.f90'
    include 'shockwave_DG.f90'
    
    
    function jwl(rho) result(p)
    implicit none
    real rho, p, rv
    
    rv = jwl_rho0/rho
    
    p = A*exp(-R1*rv) + B*exp(-R2*rv)+&
        C*rv**(-1-w)
    
    end function
    
    
    function pressure_jwl(E) result(p)
    implicit none
    real rho, p, rv, e, v
    
    rv = (Rb/Rb0)**3.0
    v = 4.0/3.0*pi*Rb**3
    
    p = A*(1-w/rv/R1)*exp(-R1*rv) &
        + B*(1-w/rv/R2)*exp(-R2*rv)+&
        w*E/v
    
    end function
    end module