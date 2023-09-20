    function energy_k(x) result(res)
    implicit none
    real x ,res
    integer i
    real cv(3),r
    r = x*L+Rb
    cv = 0
    do i=1,np
        cv = cv + ans(i,:)*basis(x,i)
    enddo
    
    res = 0.5*cv(2)**2/cv(1)*4*pi*r**2
    return
    end function
    
    function equivalent_vel() result(res)
    use FSintegration
    implicit none
    real res
    real ek
    !real, external:: energy_k
    ek = adaptive_integrate_gl(0.0,0.65,energy_k,np=5)*L
    res = ek/(2.0*pi*Rb**3*rho0)
    res = sqrt(res)
    return
    end function
    
    function energy_k_dg(x) result(res)
    implicit none
    real x ,res
    integer i
    real cv(3),r
    r = x*Ldg(id_dg)+Rdg(id_dg)
    cv = 0
    do i=1,np
        cv = cv + ans_dg(i,:,id_dg)*basis(x,i)
    enddo
    
    res = 0.5*cv(2)**2/cv(1)*4*pi*r**2
    return
    end function
    
    function equivalent_vel_dg() result(res)
    use FSintegration
    implicit none
    real res
    real ek,xend,x1,x2
    integer i!,id_dg
    logical lfin
    !real, external:: energy_k
    xend = 0.65
    xend = xend*L
    ek = 0
    
    do id_dg=1,ndg
        x1 = 0
        lfin = .false.
        if(xend>Ldg(id_dg))then
            x2 = 1.0
        else
            x2 = xend/Ldg(id_dg)
            lfin = .true.
        endif
        xend = xend - Ldg(id_dg)
        ek = ek+&
            adaptive_integrate_gl(x1,x2,energy_k_dg,np=5)*Ldg(id_dg)
        if(lfin)exit
    enddo
    res = ek/(2.0*pi*Rb**3*rho0)
    res = sqrt(res)
    return
    end function