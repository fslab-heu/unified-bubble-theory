    function int_rhs_far(x) result(res)
    use fsintegration
    implicit none
    ! input
    real x
    ! output
    real res
    ! local
    real press,ch,csound,px,r
    integer i
    ! calculate conservative vectors
    r = Re+x*Lp
    press = 0
    do i=1,np
        press = press + basis(x,i)*pans(i)
    enddo
    csound = c0 + Kc*press
    ch = csound - dRs
    res = dbasis(x,s_nb)*far_field_flux(press)/Lp    ! flux
    res = res + basis(x,s_nb)*far_field_source(press,r)          ! source

    return
    end function
    
    function far_field_flux(press) result(res)
    implicit none
    real press, res
    res = (c0-dRs)*press + 0.75*MKc*press**2
    return
    end function
    
    function far_field_source(press,r) result(res)
    implicit none
    real press, res,r
    res = -(c0+MKc*press)*press/r
    return
    end function
    
    function far_field_march() result(dans)
    use fsintegration
    implicit none
    real dans(npmax)
    integer i
    real press,csound,ch
    real eps,rho1                ! energy loss
    
    dans = 0
    press = 0
    do i=1,np
        press = press + basis(1.0,i)*pans(i)
    enddo
    !csound = c0 + Kc*press
    
    
    ! solve the shock front with the Hugoniot condition
    ! to find the energy loss
    ! Note: press  =  ps - p0
    
    Ma = ((press+pw+p0)/(p0+pw)*(gamma+1)+gamma-1)/2/gamma
    Ma = sqrt(Ma)
    dRs = Ma*c0
    rhos = rho0*(gamma+1)*Ma**2/((gamma-1)*Ma**2+2.0)
    rho1 = rhos*((p0+pw)/(press+pw+p0))**(1.0/gamma)
    eps = 0.5*(2*p0+press)*(1/rho0-1/rhos)+&
        pw*(1/rho1-1/rhos)+&
        (press+p0+pw)/(gamma-1)*(rho1**(gamma-1)-rhos**(gamma-1))/rhos**gamma
    
    eps = eps*dRs/2.0/(press+p0)*c0**2*rho0**2
    !eps = 0
    
    
    ! integrating the flux and source term
    prhs = 0
    do i=1,np
        s_nb = i
        prhs(i) = adaptive_integrate_gl(0.0,1.0,int_rhs_far,np=5)
    enddo

    ! boundary conditions
    
    ! right end
    
    do i=1,np
        prhs(i) = prhs(i) - basis(1.0,i)/Lp*far_field_flux(press)   ! flux
        !prhs(i) = prhs(i) - basis(1.0,i)*eps/Lp
    enddo
    
    ! left end
    press = 0
    do i=1,np
        press = press + basis(0.0,i)*pans(i)
    enddo
    do i=1,np
        prhs(i) = prhs(i) + basis(0.0,i)/Lp*far_field_flux(press)
    enddo
    
    ! solve the equation
    
    dans = matmul(imatrix(1:np,1:np),prhs(1:np))
    end function