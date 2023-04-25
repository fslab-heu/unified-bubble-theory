    ! module to hold the variables and functions of bubble
    module bubble
        private:: length
        integer nbubble
        type sbubble
            ! static variables
            integer id,fid
            ! current variables
            real r,dr,ddr,dddr,p,vel(3),center(3),h,dh
            ! ambient fluid conditions
            real pamb,dpamb,grad_pamb(3),uamb(3),duamb(3)
            ! initial conditions
            real r0,p0,rm,tref,center0(3),ts
            logical,allocatable:: impacted(:,:)
            logical started
            
            ! history for interpolation
            ! 1-time 2-R 3-dR 4-ddR 5-P 6-H 7-dH
            real,allocatable:: bhis(:,:)
        contains
        procedure ini_bubble
        procedure collect_ambient
        procedure advance
        procedure print_bubble  
        procedure update
        end type
        type(sbubble):: bubbles(100)
    contains
    
    subroutine update(this)
    ! update variables at a new time increment
    use global
    implicit none
    class(sbubble):: this
    
    this%p = pressure(this,this%r)
    this%h = enthalpy(this,this%r)
    this%dh = denthalpy(this,this%r,this%dr)
    
    this%bhis(inc+1,1:7) = [t+dt,&
    this%r,this%dr,this%ddr,&
    this%p,this%h,this%dh]
    
    end subroutine
    
    subroutine ini_bubble(this,ts,pos,r,dr,p0,id)
    ! initialize a bubble
    use global
    class(sbubble):: this
    real pos(3),r,dr,p0
    integer id
    real sol1(4),dat1(7)
    character*100 name
    this%ts=ts
    this%r0=r
    this%p=p0
    this%center=pos
    this%center0=pos
    this%r=r
    this%dr=dr
    this%p0=p0
    this%vel=0
    this%p=p0
    this%pamb=pamb - rho*g*this%center(3)
    this%dpamb=0
    this%uamb=uamb
    this%id=id
    this%fid=id+100
    this%started=.false.
    allocate(this%bhis(1000000,8))
    this%bhis=0
    
    allocate(this%impacted(nbubble,2))
    this%impacted = .false.
    
    
    
    ! calculate rm with the initial conditions
    this%rm = 2*pi*this%r0**3*rho*this%dr**2 &   ! kinetic energy
        + 4.0/3.0*this%p0/(gamma-1)*this%r0**3 & ! internal energy
        + 4.0/3.0*pi*this%r0**3.0*this%pamb
    this%rm = (0.75*this%rm / this%pamb /pi)**0.33333
    
    this%tref = this%rm*sqrt(rho/this%pamb)

    
    call this%update()
    
    call this%collect_ambient(1)
    dat1 = [this%r,this%dr,&
        enthalpy(this,this%r),denthalpy(this,this%r,this%dr),&
        this%vel(1),this%vel(2),this%vel(3)]
    
    sol1 = solve(this,dat1)
    this%ddr = sol1(1)
    
    this%bhis(1,4) = this%ddr
    
    ! open file for output
    write(name,*) id
    open(this%fid,file=trim(path)//'/bubble_'//trim(adjustl(name))//'.dat')
    write(this%fid,"(10A15)")'#time','R','pressure','-rho*dphi*R','ux','uy','uz','vx','vy','vz'
    return
    end subroutine
    
    subroutine collect_ambient(this,stg)
    ! calculate the aimbient variables of the bubble center
    use global
    use boundary
    implicit none
    class(sbubble):: this
    
    integer i,stg
    real pos(3)
    real vel(3),dphi,ddphi
    real time
    real grad_dphi(3), duamb(3)
    pos = this%center
    this%pamb = pamb - rho*g*this%center(3)
    this%uamb = uamb
    this%dpamb = 0
    this%grad_pamb = [0.0,0.0,-rho*g]
    this%duamb = 0

    
    do i=1,nbubble
        if(i==this%id)cycle
        if(bubbles(i)%started)then
        call induced_flow(bubbles(i),pos,&
            vel,duamb,dphi,ddphi,grad_dphi,stg,this%impacted(i,:))
        this%uamb = this%uamb+vel
        this%pamb = this%pamb - rho*dphi
        this%dpamb = this%dpamb - rho*ddphi
        this%grad_pamb = this%grad_pamb - rho *grad_dphi
        this%duamb = this%duamb + duamb
        endif
    enddo

    if(ibound)then
        do i=1,nbubble
            if(bubbles(i)%started)then
            call induced_flow_bound(bubbles(i),pos,&
                vel,duamb,dphi,ddphi,grad_dphi,stg)
            this%uamb = this%uamb+vel
            this%pamb = this%pamb - rho*dphi
            this%dpamb = this%dpamb - rho*ddphi
            this%grad_pamb = this%grad_pamb - rho *grad_dphi
            this%duamb = this%duamb + duamb
            endif
        enddo
    endif

    this%pamb = this%pamb - 0.5*rho*(length(this%uamb))**2.0
    this%dpamb = this%dpamb - rho*dot_product(this%uamb,this%duamb)
    
    end subroutine
    
    
    function collect_induced_field(pos,time) result(output)
    ! calculate the pressure at a given position
    use global
    use boundary
    implicit none
    integer i
    real pos(3)
    real vel(3),dphi,ddphi
    real time
    real grad_dphi(3), duamb(3)
    real output(11)
    
    output(1) = pamb - rho*g*pos(3) ! pressure
    output(2) = 0                   ! dp
    output(3:5) = uamb              ! uamb
    output(6:8) = [0.0,0.0,-rho*g]  ! nabla p
    output(9:11) = 0                ! duamb
    
    do i=1,nbubble
        if(bubbles(i)%started)then
        call induced_flow(bubbles(i),pos,vel,duamb,dphi,ddphi,grad_dphi,1)
        output(3:5) = output(3:5)+vel
        output(1) = output(1) - rho*dphi
        output(2) = output(2) - rho*ddphi
        output(6:8) = output(6:8) - rho *grad_dphi
        output(9:11) = output(9:11) + duamb
        endif
    enddo

    if(ibound)then
        do i=1,nbubble
            if(bubbles(i)%started)then
            call induced_flow_bound(bubbles(i),pos,vel,duamb,dphi,ddphi,grad_dphi,1)
            output(3:5) = output(3:5)+vel
            output(1) = output(1) - rho*dphi
            output(2) = output(2) - rho*ddphi
            output(6:8) = output(6:8) - rho *grad_dphi
            output(9:11) = output(9:11) + duamb
            endif
        enddo
    endif

    output(1) = output(1) - 0.5*rho*(length(output(3:5)))**2.0
    output(2) = output(2) - rho*dot_product(output(3:5),output(9:11))
    
    end function
    
    
    
    
    subroutine induced_flow(this,pos,vel,dvel,&
        dphi,ddphi,grad_dphi,stg,impacted)
    use global
    implicit none
    class(sbubble):: this
    real vel(3),dphi,pos(3),ddphi
    real center(3),sub(3),dis
    real tint,time
    integer i,m,n,stg
    real wm,wn,tinc
    real r, dr, h, ddr,dh
    real dvel(3),grad_dphi(3),dir(3)
    logical,optional:: impacted(:)
    
    if(stg==1)then
        time=t
    elseif(stg==2)then
        time=t+0.5*dt
    elseif(stg==3)then
        time=t+dt
    endif
    
    vel = 0
    dphi = 0
    ddphi = 0
    dvel = 0
    grad_dphi =0
    
    center=this%center
    sub = center-pos
    dis = length(sub)
    dir = -sub/dis
    tint = time - (dis-this%r)/c
    if(tint<this%ts)then
        m=1;n=2;wm=0;wn=0;
        return
    else
        do i=1,inc
            if(tint<=this%bhis(i,1))exit
        enddo
        if(i>inc)then
            print*,'Error in interpolation!'
            print*,'Check for the distance between the bubbles!'
            stop
        endif
        m=i-1;n=i;
        tinc=this%bhis(n,1)-this%bhis(m,1)
        wn = (tint-this%bhis(m,1))/tinc
        wm = 1-wn

    endif
    
    r = wm*this%bhis(m,2)+wn*this%bhis(n,2)
    dr = wm*this%bhis(m,3)+wn*this%bhis(n,3)
    ddr = wm*this%bhis(m,4)+wn*this%bhis(n,4)
    
    ! h and dh must be intepolated from existing data
    ! since the ambient pressure are not stored for history
    h = wm*this%bhis(m,6)+wn*this%bhis(n,6)
    dh = wm*this%bhis(m,7)+wn*this%bhis(n,7)

    vel = dir*r/dis**2*(r*dr+(dis-r)/c*(h+0.5*dr**2))
    dphi = -r/dis*(h+0.5*dr**2.0)
    ddphi = -(dr*(h+0.5*dr**2.0) + r*(dh+dr*ddr))/dis
    ! calculate grad_dphi
    grad_dphi = -dphi*dir/dis - &
        dir*ddphi/c
    ! calculate dvel
    dvel = - sub / dis**3.0*(&
        this%bhis(n,2)/c*(dis-this%bhis(n,2))*(         &
        this%bhis(n,6)+0.5*this%bhis(n,3)**2.0)+        &
        this%bhis(n,2)**2.0*this%bhis(n,3)) +           &
        sub / dis**3.0*(this%bhis(m,2)/c*(dis-          &
        this%bhis(m,2))*(this%bhis(m,6)+0.5*            &
        this%bhis(m,3)**2.0)+this%bhis(m,2)**2.0*       &
        this%bhis(m,3))
    dvel =dvel/tinc
    !
    ! correction for the initial shock
    
    if(present(impacted))then
        if(.not.all(impacted).and.(wm+wn)>0.5)then
            if((stg==1.and.(.not.impacted(1))).or.      &
                (stg==2.and.impacted(1)))then
                dh = dh + h/dt
                ddphi = ddphi + dphi/dt
                grad_dphi = grad_dphi - dir*dphi/(dt*c)
                dvel = dvel + vel/dt
                if(stg==1)impacted(1)=.true.
                if(stg==2)impacted(2)=.true.
            endif
        endif
    endif
    end subroutine
    
    subroutine induced_flow_bound(this,pos,vel,dvel,    &
        dphi,ddphi,grad_dphi,stg,impacted)
    use global
    use boundary
    implicit none
    class(sbubble):: this
    real vel(3),dphi,pos(3),ddphi
    real center(3),sub(3),dis
    real tint,time
    integer i,m,n,stg
    real wm,wn,tinc
    real r, dr, h, ddr,dh
    real dvel(3),grad_dphi(3),alpha,dir(3)
    logical,optional:: impacted(2)
    if(stg==1)then
        time=t
    elseif(stg==2)then
        time=t+0.5*dt
    elseif(stg==3)then
        time=t+dt
    endif
    
    vel = 0
    dphi = 0
    ddphi = 0
    dvel = 0
    grad_dphi =0
    
    
    alpha = bound1%alpha
    center=this%center
    center = center-(2*dot_product(bound1%norm,         &
        center+bound1%b))*bound1%norm
    sub = center-pos
    dis = length(sub)
    dir = - sub / dis
    tint = time - (dis-this%r)/c
    if(tint<this%ts)then
        m=1;n=2;wm=0;wn=0;
        return
    else
        do i=1,inc
            if(tint<this%bhis(i,1))exit
        enddo
        if(i>inc)then
            print*,'Error in interpolation for image bubble!'
            stop
        endif
        
        m=i-1;n=i;
        tinc=this%bhis(n,1)-this%bhis(m,1)
        wn = (tint-this%bhis(m,1))/tinc
        wm = 1-wn
    endif
    
    r = wm*this%bhis(m,2)+wn*this%bhis(n,2)
    dr = wm*this%bhis(m,3)+wn*this%bhis(n,3)
    ddr = wm*this%bhis(m,4)+wn*this%bhis(n,4)
    
    ! h = enthalpy(this,r,3)
    ! dh = denthalpy(this,r,dr,3)
    ! h and dh must be intepolated from existing data
    ! since the ambient pressure are not stored for history
    h = wm*this%bhis(m,6)+wn*this%bhis(n,6)
    dh = wm*this%bhis(m,7)+wn*this%bhis(n,7)
    

    vel = dir*r/dis**2*(r*dr+(dis-r)/c*(h+0.5*dr**2))
    dphi = -r/dis*(h+0.5*dr**2.0)
    ddphi = -(dr*(h+0.5*dr**2.0) + r*(dh+dr*ddr))/dis
    ! calculate grad_dphi
    grad_dphi = -dphi*dir/dis - &
        dir*ddphi/c
    ! calculate dvel
    dvel = - sub / dis**3.0*(&
        this%bhis(n,2)/c*(dis-this%bhis(n,2))*(         &
        this%bhis(n,6)+0.5*this%bhis(n,3)**2.0)+        &
        this%bhis(n,2)**2.0*this%bhis(n,3)) +           &
        sub / dis**3.0*(this%bhis(m,2)/c*(dis-          &
        this%bhis(m,2))*(this%bhis(m,6)+0.5*            &
        this%bhis(m,3)**2.0)+this%bhis(m,2)**2.0        &
        *this%bhis(m,3))
    dvel =dvel/tinc
    ! correction for the initial shock
    
    if(present(impacted))then
        if(.not.all(impacted).and.(wm+wn)>0.5)then
            if((stg==1.and.(.not.impacted(1))).or.&
                (stg==2.and.impacted(1)))then
                dh = dh + h/dt
                ddphi = ddphi + dphi/dt
                grad_dphi = grad_dphi - dir*dphi/(dt*c)
                dvel = dvel + vel/dt
                if(stg==1)impacted(1)=.true.
                if(stg==2)impacted(2)=.true.
            endif
        endif
    endif
    vel=vel*alpha
    dvel=dvel*alpha
    dphi=dphi*alpha
    ddphi=ddphi*alpha
    grad_dphi=grad_dphi*alpha
    end subroutine
    
    real function length(vect)
        real vect(:)
        integer N
        N=size(vect)
        length=sqrt(sum(vect**2))
        return
    end function
    
    
    
    subroutine advance(this)
    use global
    implicit none
    class(sbubble):: this
    real dat1(7),dat2(7),dat3(7),dat4(7)
    real sol1(4),sol2(4),sol3(4),sol4(4)
    
    if(t+dt>this%ts)then
        this%started=.true.
        call this%collect_ambient(1)
        dat1 = [this%r,this%dr,enthalpy(this,this%r),       &
            denthalpy(this,this%r,this%dr),&
            this%vel(1),this%vel(2),this%vel(3)]
    
        sol1 = solve(this,dat1)
    
        call this%collect_ambient(2)
        dat2(1:2) = dat1(1:2)+0.5*dt*[dat1(2),sol1(1)]
        dat2(3:4) = [enthalpy(this,dat2(1)),                &
            denthalpy(this,dat2(1),dat2(2))]
    
        dat2(5:7) = this%vel+0.5*dt*sol1(2:4)

        sol2 = solve(this,dat2)
        this%r = this%r+dt*dat2(2)
        this%dr = this%dr+dt*sol2(1)
        this%ddr = sol1(1)+2*(sol2(1)-sol1(1))
    
    
        this%center = this%center+dt*dat2(5:7)
        this%vel = this%vel+ dt*sol2(2:4)
    endif
    call this%collect_ambient(3)
    
    call this%update()
    end subroutine
    
    function rk(this,input,slop,step)
    class(sbubble):: this
    real rk(4)
    real input(4),slop(2)
    integer step
    rk = 0
    end function
    
    
    subroutine print_bubble(this)
    use global
    class(sbubble):: this
    write(this%fid,"(10E15.7)") t,this%r,&
        this%p,this%R*rho*(this%h+this%dr**2),&
        this%center-this%center0, this%vel
    end subroutine
    
    real function enthalpy(this,r)
    use global
    class(sbubble):: this
    real r
    real pb,pa
    real e
    pb = pressure(this,r) 
    pa = this%pamb
    e = (pb-pa)/rho/c**2
    enthalpy = e-0.5*e**2
    enthalpy = enthalpy*c**2
    
    end function
    
    real function pressure(this,r)
    use global
    implicit none
    class(sbubble):: this
    real r
    pressure = this%p0 *(this%r0/r)**(3*gamma)-     &
        2*sigma/this%R -                            &    ! surface tension
        4*miu*this%dR/this%R                        &    ! viscosity
        + Pv                                             ! Pv
    end function
    
    real function dpressure(this,r,dr)
    use global
    implicit none
    class(sbubble):: this
    real r,dr
    dpressure = -3*gamma*pressure(this,r)/r*dr+     &
        2*sigma/this%R**2*this%dR +                 &   ! surface tension
        4*miu*this%dR**2/this%R**2                      ! viscosity
    end function
    
    
    real function denthalpy(this,r,dr)
    use global
    implicit none
    class(sbubble):: this
    real r,dr
    real dpa,dpb,e,de,pa,pb
    
    pb = pressure(this,r)
    pa = this%pamb
    e= (pb-pa)/rho/c**2
    dpa = this%dpamb
    dpb = dpressure(this,r,dr)
    de = (dpb-dpa)/rho/c**2
    denthalpy = de- e*de+ e**2*de;
    denthalpy = denthalpy*c**2
    
    end function
    
    
    function solve(this,input)
    use global
    implicit none
    class(sbubble):: this
    real,dimension(4):: solve
    real input(7)
    real r,dr,p,h
    real A,u2,uref(3)
    real du2, dh,du(3),u(3)
    real ss(3)
    real Pa,dPa, NPa(3)
    r=input(1)
    dr=input(2)
    h=input(3)
    dh = input(4)
    u = input(5:7)
    uref = u-this%uamb
    u2 = dot_product(uref,uref)
    
    Pa = this%pamb
    dPa = this%dpamb
    NPa = this%grad_pamb
    
    ! solution for vel
    ss = u-this%uamb
    ss = ss*length(ss)
    solve = 0
    if(imigration)then
        solve(2:4) = - ( 24.0*1e9*Ca*dR*rho*uref +      &
            8.0*1e9*R*0*rho*uref                        &
                + 3.0*1e9*Cd*rho*ss                     &
               -12.0*dR**2.0*rho*uref                   &
                + 8.0*1e9*NPA*R                         &
                 )/4.0/R/rho/(2*Ca*1e9-dR)
        
    endif
    
    solve(1) = ((-6.0 * dR ** 2 * C + C * u2 + 2.0    &
     * dot_product(uref,solve(2:4)) * R + 2.0 * dR ** 3   &
      + dR * u2 + 4.0 * C * H + 4.0 * dH * R +        &
      4.0 * H * dR) / R / (C - dR)) / 4.0
    solve(2:4) = solve(2:4) + this%duamb 
    
    
    return
    
    end function
    
    subroutine collect_dt()
    use global
    implicit none
    integer i,j
    real dis
    dt = 1e10
    do i=1,nbubble
        dt = min(dt,&
            0.01*bubbles(i)%tref*(bubbles(i)%r/bubbles(i)%rm)**2.0)
        do j=1,i-1
            dis=length(bubbles(i)%center-bubbles(j)%center)
            dt = min(dt,dis/c*0.3)
        enddo
        
    enddo
    dt = sdt*dt
    end subroutine
    
    function pressure_at_position(pos)
    use global
    use boundary
    implicit none
    real pos(3)
    real pressure_at_position
    integer i
    real dphi,vel(3),ddphi,dphi1,vel1(3)
    real grad_phi(3),dvel(3)
    dphi=0;
    vel=uamb
    
    do i=1,nbubble
        call induced_flow(bubbles(i),pos,vel1,dvel,     &
            dphi1,ddphi,grad_phi,1)
        vel=vel+vel1
        dphi=dphi+dphi1
    enddo
    if(ibound)then
        do i=1,nbubble
            call induced_flow_bound(bubbles(i),pos,     &
                vel1,dvel,dphi1,ddphi,grad_phi,1)
            vel=vel+vel1
            dphi=dphi+dphi1
        enddo
    endif
    pressure_at_position = -rho*(dphi +                 &
        0.5*dot_product(vel,vel)) &
        + rho * 0.5*dot_product(uamb,uamb)
    end function
    
    
    end module