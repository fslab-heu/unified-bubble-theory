!  main.f90
!
!

!****************************************************************************
!
!  PROGRAM: undex_ubt
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program undex_ubt
    use math
    !use char_op
    use linear_equation
    use shockwave
    ! modules for bubble
    use global,only: t,csound=>c,tend,dtb=>dt,&
        inc,pressure_loc,iarrive,t_arrive,t_start,&
        t_shock,pamb,tdelay,fpath
    use bubble
    implicit none
    real rp
    integer n,n_arrive
    logical iexit
    real pdans(100)
    real dans(100,3),t_offset,pout,drb0,R0,Pg0
    real state(11)
    character*(100) line
    
    call initialization(rp)
    open(104,file=trim(fpath)//'output/pressure.dat')
    open(105,file=trim(fpath)//'output/bubble.dat')
    write(105,'(3A15)') '#time','radius','migration'
    time = 0
    write(104,'(2E15.6)') 0,0
    call output(Rp)
    n = 0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !             first step, 1 DG cell, high order
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    print*,'-------------------------------------------------'
    print*,'Start near-field solver!'
    iarrive = .false.
    print*,'-------------------------------------------------'
    print'(A10,3A13)','increment','time','Rs','Rb'
    print*,'-------------------------------------------------'
    do while(.true.)
        ans1 = ans
        Rb1 = Rb
        Rs1 = Rs
        dans(1:npmax,:) = march()
        ans = ans + 0.5*dt*dans(1:npmax,:)
        Rb = Rb + 0.5*dt*dRb
        Rs = Rs + 0.5*dt*dRs
        dans(1:npmax,:) = march()
        ans = ans1 + dt*dans(1:npmax,:)
        Rb = Rb1 + dt*dRb
        Rs = Rs1 + dt*dRs
        time = time+dt
        n = n+1
        call output(Rp)
        
        if (L>2*Rb)then
            exit
        endif

        if(L>1.2*Rb.and.np==8)then
            np = 9
            call inv_lie(mmatrix(1:np,1:np),imatrix(1:np,1:np),0)
            call inv_lie(mmatrix(1:np,1:np),imatrix(1:np,1:np),0)
        elseif(L>0.8*Rb.and.np==7)then
            np = 8
            call inv_lie(mmatrix(1:np,1:np),imatrix(1:np,1:np),0)
        elseif(L>0.5*Rb.and.np==6)then
            np = 7
            call inv_lie(mmatrix(1:np,1:np),imatrix(1:np,1:np),0)
        elseif(L>0.3*Rb.and.np==5)then
            np = 6
            call inv_lie(mmatrix(1:np,1:np),imatrix(1:np,1:np),0)
        elseif(L>0.2*Rb.and.np==4)then
            np = 5
            call inv_lie(mmatrix(1:np,1:np),imatrix(1:np,1:np),0)
        elseif(L>0.1*Rb.and.np==3)then
            np = 4
            call inv_lie(mmatrix(1:np,1:np),imatrix(1:np,1:np),0)
        elseif(L>0.05*Rb.and.np==2)then
            np = 3
            call inv_lie(mmatrix(1:np,1:np),imatrix(1:np,1:np),0)
        elseif(L>0.02*Rb.and.np==1)then
            np = 2
            call inv_lie(mmatrix(1:np,1:np),imatrix(1:np,1:np),0)
        endif
        ! update dt
        dt = L/(np+1)**2/(dRs+cs)
        if(mod(n,10)==0)   write(105,'(3E15.5)') time,Rb,0.0
        if(mod(n,1000)==0)then
                print*,'-------------------------------------------------'
                print'(A10,3A13)','increment','time','Rs','Rb'
                print*,'-------------------------------------------------'
        endif
        if(dt<0)stop
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !             second step, 20 DG cell, low order
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call DG_init()
    dt= 0.5*dt
    call output_dg(Rp)

    iexit = .false.
    n_arrive = 0
    do while(.true.)
        time=time+dt
        
        ! backup data
        ans_dg1=ans_dg
        Rb1=Rb
        Rs1=Rs
        Rdg1=Rdg
        ! flux
        call pp_limiter()
        call calculate_flux()
        
        
        do id_dg=1,ndg
            dans = 0
            dans(1:npmax,:) = march_dg()
            ans_dg(:,:,id_dg) = ans_dg(:,:,id_dg) + &
                0.5*dt*dans(1:npmax,:)
        enddo
        ans_dg(3:np,:,ndg)=0
        Rb = Rb+  0.5*dt*dRb
        Rs = Rs + 0.5*dt*dRs
        Rdg = Rdg + 0.5*dt*dRdg
        ! update geometry variables
        call calculate_flux()
        do id_dg=1,ndg
            dans = 0
            dans(1:npmax,:) = march_dg()
            ans_dg(:,:,id_dg) = ans_dg1(:,:,id_dg) + &
                dt*dans(1:npmax,:)
        enddo
        ans_dg(3:np,:,ndg)=0
        Rb = Rb1+  dt*dRb
        Rs = Rs1 + dt*dRs
        Rdg = Rdg1 + dt*dRdg
        dt = dt0/(np+1)*0.2
        n=n+1
        
        if(Rs>Rp.and.(.not.iarrive))then
            iarrive  = .true.
            t_arrive = time
            n_arrive = n
        endif
        
        if(mod(n - n_arrive,10)==0)then
            call output_dg(Rp)
        endif
        if(Rs>6*Rb.and..not.iexit)then
            dRb0 = equivalent_vel_dg()
            R0 = Rb
            pg0=pg
            t_start = time
            iexit = .true.
            !exit
        endif
        if(time>2*t_start.and.iexit)then
            exit
        endif
        if(mod(n,100)==0)then
            print'(I10,3E13.3)',n,time,Rs,Rb
        endif
        if(mod(n,10)==0.and..not.iexit)&
            write(105,'(3E15.5)') time,Rb,0.0
        if(mod(n,2000)==0)then
            print*,'-------------------------------------------------'
            print'(A10,3A13)','increment','time','Rs','Rb'
            print*,'-------------------------------------------------'
        endif
    enddo
    

    print*,'-------------------------------------------------'
    print*,'-------------------------------------------------'
    print*,'Finished near-field solver!'
    print*,'-------------------------------------------------'
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !             third step, farfield shock wave
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    print*,'-------------------------------------------------'
    print*,'Start the solution of far-field solver!'
    print*,'-------------------------------------------------'
    print*,'-------------------------------------------------'
    ! L2 project to find the initial solution of p
    call ini_pressure()
    ! solve the transportation equation of shock wave energy
    iexit = .false.
    do while(Re<Rp)
        if(Rs<Rp)then
            dt = min(Rs / 100.0 / c0,&
                (Rp+0.01*Lp-Rs)/c0)
        else
            dt = Lp / 100.0 / c0
        endif
        
        if(Re+dt*c0>Rp) then
            dt = (Rp-Re)/c0
            iexit = .true.
        endif
        pans1 = pans
        Re1= Re
        Rs1 = Rs
    
        pdans(1:npmax) = far_field_march()
        pans = pans + 0.5*dt*pdans(1:npmax)
        Re = Re + 0.5*dt*dRs
        Rs = Rs + 0.5*dt*dRs
        pdans(1:npmax) = far_field_march()
        Re = Re1 + dt*dRs
        Rs = Rs1 + dt*dRs
        pans = pans1 + dt*pdans(1:npmax)
        time = time+dt
        
        if(Rs>Rp.and..not.iarrive)then
            iarrive = .true.
            t_arrive = time
        endif
        
        if(iarrive)then
            call far_field_output(Rp)
        endif
        if(iexit)exit
        n = n+1
    enddo
    t_offset = time - t_arrive
    
    print*,'-------------------------------------------------'
    print*,'-------------------------------------------------'
    print*,'Finished far-field solver!'
    print*,'-------------------------------------------------'

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !             4th step, bubble solving
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    print*,'-------------------------------------------------'
    print*,'Start unified bubble dynamics solver!'
    print*,'-------------------------------------------------'
    print*,'-------------------------------------------------'
    call initiation(R0,dRb0,pg0)
    inc = 0
    t = 0
    t_shock = time-t_arrive+tdelay  ! end time of shock wave
    print*,'-------------------------------------------------'
    print'(A10,3A13)','increment','time','Rb','Zb'
    print*,'-------------------------------------------------'
    do while(t + t_start <tend)
        call collect_dt()
        call bubbles(1)%advance
        time = t + (Rp-bubbles(1)%R)/c0
        if(mod(inc,20)==0 .and. time+t_start - t_arrive > t_shock)then
            !state = collect_induced_field(pressure_loc,time)
            pout = rho0*(2*bubbles(1)%R*bubbles(1)%dR**2+&
                bubbles(1)%R**2*bubbles(1)%ddR)/(rp-bubbles(1)%R)
            write(104,'(2E15.6)') time + t_start - t_arrive + tdelay, &
                pout
        endif
        if(mod(inc,1000)==0)then
            print'(I10,3E13.3)',inc,t,bubbles(1)%R,bubbles(1)%center(3)
        endif
        if(mod(inc,20000)==0)then
            print*,'-------------------------------------------------'
            print'(A10,3A13)','increment','time','Rb','Zb'
            print*,'-------------------------------------------------'
        endif
        if(mod(inc,20)==0)   write(105,'(3E15.5)') t+t_start,bubbles(1)%R,bubbles(1)%center(3)
        inc=inc+1
        t=t+dtb
    enddo
    close(104)
    close(105)
    
    
    print*,''
    print*,''
    print*,''
    print*,'-------------------------------------------------'
    print*,'Simulation Completed!'
    print*,'Results can be found in the folder of output.'
    end program undex_ubt

