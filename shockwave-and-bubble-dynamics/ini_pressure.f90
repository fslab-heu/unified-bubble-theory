    subroutine ini_pressure()
    use fsintegration
    use linear_equation
    implicit none
    integer i,j,k,npwave
    
    npwave = 16
    
    Lp = 0.5*Rs
    Re = Rs - Lp
    Kc = 0.5*(gamma-1)/sqrt(gamma*rho0*(p0+pw))
    MKc = Kc! + 1.0/(rho0*c0)
    allocate(pans(npmax),prhs(npmax),pans1(npmax))
    pans = 0
    prhs = 0
    do i = 1,npwave
        s_nb = i
        prhs(i) = adaptive_integrate_gl(0.0,1.0,pressure_L2,&
            debug = .false.,np = 5)
    enddo
    np = npwave
    call inv_lie(mmatrix(1:np,1:np),imatrix(1:np,1:np),0)
    pans(1:np) = matmul(imatrix(1:np,1:np),prhs(1:np))
    end subroutine
    
    function pressure_L2(xp) result(res)
    implicit none
    real xp
    real res
    ! local
    real xe,r,cv(3),press
    integer i
    r = Re + xp*Lp
    
    do i=1,ndg
        if(r>Rdg(ndg+1-i))then
            exit
        endif
    enddo
    id_dg = ndg+1-i
    xe = (r-Rdg(id_dg))/Ldg(id_dg)
    cv = 0
    do i=1,np
        cv = cv + ans_dg(i,:,id_dg)*basis(xe,i)
    enddo
    press = e2p(cv(3) - 0.5*cv(2)**2/cv(1))
    press = press - p0 ! use dynamic pressure in wave equation
    res = press * basis(xp,s_nb)
    return
    end function