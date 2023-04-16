
    
    function int_rhs(x) result(res)
    use fsintegration
    implicit none
    real x
    real res,cv(3),cf,pp,uhx,ux,rhoux,rhox
    integer i
    real uh,u,xx,ein
    ! calculate conservative vectors
    cv = 0
    do i=1,np
        cv = cv + basis(x,i)*ans(i,:)
    enddo
    u = cv(2)/cv(1)
    uh = (u-x*dL-dRb)/L
    xx = x*L+Rb
    res = 0
    if(s_nb2==1)then            ! rho
        res = dbasis(x,s_nb)*uh*cv(1)               ! f
        res = res - basis(x,s_nb)*(2*cv(2)/xx + dL/L*cv(1))          ! s
    elseif(s_nb2==2)then        ! rho*u
        ! calculate pressure
        !pp = (cv(1)/1000)**gamma*cv(3) - pw
        ein = cv(3) - 0.5*cv(2)**2/cv(1)
        pp = e2p(ein)
        !pp = entro2pressure(cv(3),cv(1))
        
        res = dbasis(x,s_nb)*(uh*cv(2)+pp/L)                        ! f
        res = res - basis(x,s_nb)*(2*cv(1)*u**2/xx + dL/L*cv(2))    ! s
    elseif(s_nb2==3)then        ! s
        ein = cv(3) - 0.5*cv(2)**2/cv(1)
        pp = e2p(ein)
        res = dbasis(x,s_nb)*(uh*cv(3) + u*pp/L)                    ! f
        res = res - basis(x,s_nb)*(2*u*(cv(3)+pp)/xx + dL/L*cv(3))  ! s
    endif

    return
    end function
    
    function march() result(dans)
    use fsintegration
    !use RiemannSolver
    implicit none
    real a1,a2
    !real, external:: basis
    real  rhog
    integer i,j,k
    real dans(npmax,3),px,ux,eta,dma,cvb(3),dcvb(3),cvs(3),dcvs(3)
    real uh,e0
    
    rhog = jwl_rho0*(Rb0/Rb)**3.0
    pg = jwl(rhog)
    
    ! bubble surface
    cvb = 0
    do i=1,np
        cvb = cvb + ans(i,:)*basis(0.0,i)
    enddo
    dRb = cvb(2)/cvb(1)
    
    ! shock front
    cvs = 0
    do i=1,np
        cvs = cvs + ans(i,:)*basis(1.0,i)
    enddo
    rhos = cvs(1)
    us = cvs(2)/rhos
    ps = e2p(cvs(3)-0.5*rhos*us**2)
    Ma = ((ps+pw)/(p0+pw)*(gamma+1)+gamma-1)/2/gamma
    Ma = sqrt(Ma)
    dRs = Ma*c0
    
    dL = dRs - dRb
    L = Rs - Rb
    
    ! calculate right-hand side integration
    rhs = 0
    do i=1,np
        s_nb = i
        do j=1,3
            s_nb2 = j
            rhs(i,j) = adaptive_integrate_gl(0.0,1.0,int_rhs,np=5)
        enddo
    enddo
    ! calculate the flux
    ! left end ( bubble surface)
    do i = 1,np
        rhs(i,2) = rhs(i,2) + basis(0.0,i)*pg/L
        rhs(i,3) = rhs(i,3) + basis(0.0,i)*pg*dRb/L
    enddo
    
    
    ! right end ( shock front)
    
    e0 = p2e(p0)

    do i = 1,np
        rhs(i,1) = rhs(i,1) + basis(1.0,i)*dRs*rho0/L
        rhs(i,2) = rhs(i,2) - basis(1.0,i)*(p0/L)
        rhs(i,3) = rhs(i,3) + basis(1.0,i)*e0*dRs/L
    enddo
    
    ! solve system
    
    dans = 0
    dans(1:np,:) = matmul(imatrix(1:np,1:np),rhs(1:np,:))

    end function