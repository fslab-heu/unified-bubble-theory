    subroutine update_dg()
    implicit none
    real cv(3)
    integer i
    real rhog,weight,ctemp
    
    L = Rs-Rb
    Ldg = Rdg(2:ndg+1)-Rdg(1:ndg)
    ! bubble surface
    
    rhog = jwl_rho0*(Rb0/Rb)**3.0
    pg = jwl(rhog)
    cv = 0
    do i=1,np
        cv = cv + ans_dg(i,:,1)*basis(0.0,i)
    enddo
    dRb = cv(2)/cv(1)
    dRdg(1)=dRb
    ctemp = soundspeed(pg,cv(1))
    dt0 = min(dt0,Ldg(1)/(dRb+ctemp))
    ! shock front
    cv = 0
    do i=1,np
        cv = cv + ans_dg(i,:,ndg)*basis(1.0,i)
    enddo
    rhos = cv(1)
    us = cv(2)/rhos
    ps = e2p(cv(3)-0.5*rhos*us**2)
    ctemp = soundspeed(ps,rhos)
    Ma = ((ps+pw)/(p0+pw)*(gamma+1)+gamma-1)/2/gamma
    Ma = sqrt(Ma)
    Ma = min(Ma,&
        0.25/c0*(us*(gamma+1)+sqrt(gamma**2*us**2+2*gamma*us**2+16*c0**2+us**2)))
    Ma = min(Ma,&
        -sqrt(-2*rhos*(gamma*(rhos-rho0)-rhos-rho0))/(gamma*(rhos-rho0)-rhos-rho0))

    dRs = Ma*c0
    !dRs = HLLC_SR(cv,[rho0,0.0,p2e(p0)])
    
    dt0 = min(dt0,Ldg(ndg)/(dRdg(ndg+1)+us+ctemp))
    
    do i=1,ndg+1
        weight = (Rdg(i)-Rb)/L
        dRdg(i) = dRb*(1-weight) + dRs*weight
    enddo
    
    dL = dRs-dRb
    
    dLdg = dRdg(2:ndg+1)-dRdg(1:ndg)
    
    end subroutine
    
    
    subroutine DG_init()
    use fsintegration
    implicit none
    integer i,npdg
    
    npdg=4
    
    allocate(ans_dg(npmax,3,ndg),ans_dg1(npmax,3,ndg))
    ans_dg=0
    
    ! initialize geometry variables
    ! we use 3 elements
    L = Rs - Rb
    dL = dRs-dRb
    
    
    Rdg(1) = 0
    do i=1,ndg
        Ldg(i) = (bias-1)*(1-abs(i-1-(ndg-1)/2.0)*2/(ndg-1))+1.0
        Rdg(i+1) = Rdg(i) + Ldg(i)
    enddo
    Rdg = Rdg / sum(Ldg)*L+Rb
    Ldg = Ldg / sum(Ldg)*L
    

    
    ! L2 projection for the DG elements

    do s_nb=1,ndg     ! loop for element
        rhs = 0
        do s_nb2=1,3    ! loop for var
            do s_nb3=1,npdg   ! loop for polynomial
                rhs(s_nb3,s_nb2) = adaptive_integrate_gl(&
                    0.0,1.0,shock_dg_L2,debug = .false.,np = 5)
            enddo
        enddo
        ans_dg(1:npdg,1:3,s_nb) = matmul(imatrix(1:npdg,1:npdg),rhs(1:npdg,:))
    enddo
    np=npdg
    call update_dg
    ! evaluate dt
    end subroutine
    
    
    
    function shock_dg_L2(xp) result(res)
    implicit none
    real xp
    real res
    ! local
    real xe,r
    integer i
    
    ! s_nb element
    ! s_nb2 var
    ! s_nb3 polynomial
    
    r = Rdg(s_nb) + xp*Ldg(s_nb)
    xe = (r-Rb)/L

    res = 0
    do i=1,np
        res = res + ans(i,s_nb2)*basis(xe,i)
    enddo
    res = res*basis(xp,s_nb3)
    return
    end function
    
    
    
    
    
    
    
    
    
    
    
    
    
    function int_rhs_dg(x) result(res)
    use fsintegration
    implicit none
    real x
    real res,cv(3),cf,pp,uhx,ux,rhoux,rhox
    integer i
    real uh,u,xx,ein,sc
    ! calculate conservative vectors
    cv = 0
    do i=1,np
        cv = cv + basis(x,i)*ans_dg(i,:,id_dg)
    enddo
    u = cv(2)/cv(1)
    xx = x*Ldg(id_dg)+Rdg(id_dg)
    uh = (u-x*dLdg(id_dg)-dRdg(id_dg))/Ldg(id_dg)
    sc = dLdg(id_dg)/Ldg(id_dg)
    res = 0
    if(s_nb2==1)then            ! rho
        res = dbasis(x,s_nb)*uh*cv(1)               ! f
        res = res - basis(x,s_nb)*(2*cv(2)/xx + sc*cv(1))          ! s
    elseif(s_nb2==2)then        ! rho*u
        ein = cv(3) - 0.5*cv(2)**2/cv(1)
        pp = e2p(ein)
        res = dbasis(x,s_nb)*(uh*cv(2)+pp/Ldg(id_dg))                      ! f
        res = res - basis(x,s_nb)*(2*cv(1)*u**2/xx + sc*cv(2))    ! s
    elseif(s_nb2==3)then        ! s
        ein = cv(3) - 0.5*cv(2)**2/cv(1)
        pp = e2p(ein)
        res = dbasis(x,s_nb)*(uh*cv(3) + u*pp/Ldg(id_dg))                    ! f
        res = res - basis(x,s_nb)*(2*u*(cv(3)+pp)/xx + sc*cv(3))  ! s
    endif

    return
    end function
    
    function march_dg() result(dans)
    use fsintegration
    !use RiemannSolver
    implicit none
    real a1,a2
    !real, external:: basis
    real  rhog
    integer i,j,k
    real dans(npmax,3),px,ux,eta,dma,cvb(3),dcvb(3),cvs(3),dcvs(3)
    real uh,e0
    
    ! calculate right-hand side integration
    rhs = 0
    do i=1,np
        s_nb = i
        do j=1,3
            s_nb2 = j
            !rhs(i,j) = adaptive_integrate_gl(0.0,1.0,int_rhs_dg,np=5)
            rhs(i,j) = quadrature_legendre(0.0,1.0,int_rhs_dg,np=5)
        enddo
    enddo
    ! flux
    do i = 1,np
        rhs(i,:) = rhs(i,:) + basis(0.0,i)*flux(:,id_dg)/Ldg(id_dg)   ! left end
        rhs(i,:) = rhs(i,:) - basis(1.0,i)*flux(:,id_dg+1)/Ldg(id_dg) ! right end
    enddo
    
    ! solve system
    dans = 0
    dans(1:np,:) = matmul(imatrix(1:np,1:np),rhs(1:np,:))

    end function
    
    subroutine calculate_flux()
    
    implicit none
    real cv(3),cvs(3),UL(3),UR(3)
    real rhog,e0,ctemp
    integer i,j
    dt0 = 1e10
    
    call update_dg()
    
    ! left end
    
    
    flux(:,1) = [0.0,pg,dRb*pg]

    ! right end
    
    e0 = p2e(p0)
    flux(:,ndg+1) = [-rho0*dRs,p0,-e0*dRs]
    
    
    ! middle
    do j=2,ndg
        UL = 0
        UR = 0
        do i=1,np
            UL = UL+ ans_dg(i,:,j-1)*basis(1.0,i)
            UR = UR+ ans_dg(i,:,j)*basis(0.0,i)
        enddo
        ctemp = soundspeed(e2p(UL(3)-0.5*UL(2)**2/UL(1)),UL(1))
        dt0 = min(dt0,Ldg(j-1)/(dRdg(j)+UL(2)/UL(1)+ctemp))
        ctemp = soundspeed(e2p(UR(3)-0.5*UR(2)**2/UR(1)),UR(1))
        dt0 = min(dt0,Ldg(j)/(dRdg(j)+UR(2)/UR(1)+ctemp))
        flux(:,j) = HLLC_flux(UL,UR,dRdg(j))
    enddo
    end subroutine
    
    
    function HLLC_flux(El,Er,ui) result(res)
    implicit none
    real el(3),er(3)
    real res(3)
    real ul,ur,einl,einr,pl,pr,al,ar
    real hl,hr,u0,h0,a0,sl,sr
    real fl(3),fr(3),eint(3)
    real pvrs,ui
    real ustar,pstar

    Ul=El(2)/El(1);Ur=Er(2)/Er(1);
    Einl=El(3)/El(1)-0.5*(Ul**2);
    Einr=Er(3)/Er(1)-0.5*(Ur**2);
    pl = e2p(Einl*El(1))
    pr = e2p(Einr*Er(1))

    al = soundspeed(Pl,El(1))
    ar = soundspeed(Pr,Er(1))

    Pvrs = Pl+Pr - (Ur-Ul)*(al+ar)*(El(1)+Er(1))/4.0
    Pvrs = Pvrs/2.0
    if(Pvrs>Pl)then
        sl = Ul - al*sqrt(1+(1+gamma)/2/gamma*((Pvrs+pw)/(Pl+pw)-1))
    else
        sl = Ul - al
    endif
    if(Pvrs>Pr)then
        sr = Ur + ar*sqrt(1+(1+gamma)/2/gamma*((Pvrs+pw)/(Pr+pw)-1))
    else
        sr = Ur + ar
    endif

    Ustar = (Pr-Pl+El(1)*Ul*(Sl-Ul)-Er(1)*Ur*(Sr-Ur))/(El(1)*(Sl-Ul)-Er(1)*(Sr-Ur));
    Pstar = Pl+El(1)*(Sl-Ul)*(Ustar-Ul);
    if(Ustar-ui>0)then
        if(Sl-ui>0)then
            Eint = El
        else
            Eint(:) = El(1)*(Sl-Ul)/(Sl-Ustar);
            Eint(2)=Eint(2)*Ustar;
        endif
    else
        if(Sr-ui<0)then
            Eint = Er
        else
            Eint(:) = Er(1)*(Sr-Ur)/(Sr-Ustar);
            Eint(2)=Eint(2)*Ustar;
            Eint(3)=Eint(3)*(Er(3)/Er(1)+(Ustar-Ur)*(Ustar+Pr/Er(1)/(Sr-Ur)));
        endif
    endif
    
    res = U2F(Eint,ui)
    
    
    end function
    
    function U2F(U,ui) result(res)
    implicit none
    real U(3),ui,res(3)
    real ein,vel,pre,uh
    vel = U(2)/U(1)
    ein = U(3) - 0.5*U(1)*vel**2
    pre = e2p(ein)
    uh = vel-ui
    res(1) = uh*U(1)
    res(2) = U(2)*uh + pre
    res(3) = U(3)*uh + vel*pre
    return
    end function
    
    
    function HLLC_sr(El,Er) result(sr)
    implicit none
    real el(3),er(3)
    real res(3)
    real ul,ur,einl,einr,pl,pr,al,ar
    real hl,hr,u0,h0,a0,sl,sr
    real fl(3),fr(3),eint(3)
    real pvrs,ui

    Ul=El(2)/El(1);Ur=Er(2)/Er(1);
    Einl=El(3)/El(1)-0.5*(Ul**2);
    Einr=Er(3)/Er(1)-0.5*(Ur**2);
    pl = e2p(Einl*El(1))
    pr = e2p(Einr*Er(1))

    al = soundspeed(Pl,El(1))
    ar = soundspeed(Pr,Er(1))

    Pvrs = Pl+Pr - (Ur-Ul)*(al+ar)*(El(1)+Er(1))/4.0
    Pvrs = Pvrs/2.0
    if(Pvrs>Pl)then
        sl = Ul - al*sqrt(1+(1+gamma)/2/gamma*((Pvrs+pw)/(Pl+pw)-1))
    else
        sl = Ul - al
    endif
    if(Pvrs>Pr)then
        sr = Ur + ar*sqrt(1+(1+gamma)/2/gamma*((Pvrs+pw)/(Pr+pw)-1))
    else
        sr = Ur + ar
    endif
    end function