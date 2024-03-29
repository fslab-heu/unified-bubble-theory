    subroutine initialization(Rp)
    use linear_equation
    use polynomial
    use fsintegration
    use math,only:length
    use global,only: gamma_b=>gamma,&
        rho_w=>rho,&
        csound=>c,&
        miu,sigma,tend,sdt,pressure_loc,pv,g,cd,ca,imigration,ibound,pamb,uamb
    use boundary
    use bubble,only:nbubble
    implicit none
    real wcharge,hdepth,Rp
    real entro,cs,alpha
    integer i,j
    
    print*,''
    print*,' ___________________________________________________'
    print*,'|_|_______________________________________________|_|'
    print*,'|_|                                               |_|'
    print*,'|_|    This code is designated to simulate the    |_|'
    print*,'|_|    underwater explosion pressure load using   |_|'
    print*,'|_|    the unified bubble theory.                 |_|'
    print*,'|_|    =======================================    |_|'
    print*,'|_|    The code is developed by Yun-Long Liu in   |_|'
    print*,'|_|    the  FSLAB  of  the  Harbin  Engineering   |_|'
    print*,'|_|    University.                                |_|'
    print*,'|_|    =======================================    |_|'
    print*,'|_|    Please contact the authors for more info:  |_|'
    print*,'|_|     yunlong_liu@hrbeu.edu.cn(Yun-Long Liu)    |_|'
    print*,'|_|      zhangaman@hrbeu.edu.cn (A-Man Zhang)     |_|'
    print*,'|_|_______________________________________________|_|'
    print*,'|_|_______________________________________________|_|'
    print*,''
    print*,''
    
    
    pi = asin(1.0)*2.0
    
    ! read jwl info
    open(100,file='jwl.in')
    read(100,*) A
    read(100,*) B
    read(100,*) C
    read(100,*) R1
    read(100,*) R2
    read(100,*) w
    read(100,*) jwl_rho0
    read(100,*) jwl_E
    close(100)
    ! read water info
    open(100,file='water.in')
    read(100,*) gamma
    read(100,*) pw
    read(100,*) rho0
    close(100)
    
    ! read bubble parameters
    open(100,file='bubble.in')
    read(100,*) miu
    read(100,*) sigma
    read(100,*) Pv
    read(100,*) g
    read(100,*) cd
    read(100,*) ca
    read(100,*) alpha
    close(100)
    
    
    
    
    ! calculate C if it is not input
    if(C==0)then
        C = -A*w/R1*exp(-R1)-B*w/R2*exp(-R2) + &
            w*jwl_E
    endif
    
    ! read case parameters
    open(100,file='case.in')
    read(100,*) wcharge
    read(100,*) hdepth
    read(100,*) pressure_loc
    read(100,*) tend
    read(100,*) sdt
    close(100)
    
    print*,''
    print*,''
    print*,'================================================='
    print*,'Case info:'
    print'(4A10)','charge','depth','tend','distance'
    print'(4F10.3)',wcharge,hdepth,tend,length(pressure_loc)
    print*,'================================================='
    print*,''
    sdt = sdt*0.1
    Rp = length(pressure_loc)    ! distance
    p0 = 1.03e5+Hdepth*g*rho0
    c0 = sqrt((p0+pw)/rho0*gamma)
    ! initialize bubble
    pamb = p0
    uamb = 0
    rho_w = rho0
    gamma_b = w+1
    csound = c0
    imigration = .true.
    ibound = .false.
    nbubble = 1
    call bound1%ini_bound([0.0,0.0,1.0],[0.0,0.0,hdepth],alpha)
    
    ! number of orders
    
    npmax = 20
    
    ! initiate the legendre basis
    call initiate_legendre(npmax,0.0,1.0)
    
    np = 1
    
    ! set geometry
    Rb = (wcharge/jwl_rho0*0.75/pi)**(1/3.0)
    Rb0 = Rb
    Rs = Rb*1.01
    L = Rs - Rb
    
    ! solve the RH condition to obtain the initial velocity
    pg = jwl(jwl_rho0)
    Ma = ((pg+pw)/(p0+pw)*(gamma+1)+gamma-1)/2/gamma
    Ma = sqrt(Ma)
    cs = c0/Ma/(gamma+1)*&                           ! sound speed
        sqrt((2*gamma*Ma**2-gamma+1)*&
        (gamma*Ma**2-Ma**2+2.0))
    rhos = rho0*(gamma+1)*Ma**2/((gamma-1)*Ma**2+2.0)
    dRb = 2*c0/(gamma+1)*(Ma**2-1)/Ma                ! velocity
    dRs = Ma*c0
    dL = dRs - dRb
    !(pg+pw)/(rhos/1000)**gamma              ! divide rho by 1000 to normalize entropy
    ! set the initial condition of computational domain
    allocate(mmatrix(npmax,npmax),imatrix(npmax,npmax),ans(npmax,3),rhs(npmax,3))
    allocate(ans1(npmax,3))
    mmatrix = 0
    ans = 0
    ans(1,1) = rhos
    ans(1,2) = rhos*dRb
    ans(1,3) = p2e(pg) + 0.5*rhos*dRb**2     ! total energy

    
    ! calculate the mass matrix and its invert
    do i = 1,npmax
        do j = 1,i
            s_nb = i
            s_nb2 = j
            mmatrix(i,j) = adaptive_integrate(&
                0.0,1.0,basisn)
            if(i.ne.j)then
                mmatrix(j,i)=mmatrix(i,j)
            endif
        enddo
    enddo
    
    call inv_lie(mmatrix(1:np,1:np),imatrix(1:np,1:np),0)
    
    ! update dt
    dt = L/np/5.0/(dRb+cs)
    end subroutine