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
        inc,pressure_loc,iarrive,t_arrive,t_start,t_shock
    use bubble
    implicit none
    real rp
    integer n
    logical iexit
    real pdans(100)
    real dans(100,3),t_offset,pout,drb0,R0,Pg0
    character*(100) line
    
    call initialization(rp)
    open(100,file='./output/pressure.dat')
    open(101,file='./output/xloc.dat')
    open(102,file='./output/history.dat')
    open(103,file='./output/energy.dat')
    open(104,file='./output/p.dat')
    open(105,file='./output/bubble.dat')
    write(105,'(3A15)') '#time','radius','migration'
    write(102,'(8A15)') '#time', 'Rb', 'Rs', 'dRb', 'Modified dRb', 'Pm','Pb','Pr'
    time = 0
    write(104,'(2E15.6)') 0,0
    call output(Rp)
    !dt = 1e-5
    n = 0
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
            !print*,'Rb',np
            call inv_lie(mmatrix(1:np,1:np),imatrix(1:np,1:np),0)
        elseif(L>0.8*Rb.and.np==7)then
            np = 8
            !print*,'Rb',np
            call inv_lie(mmatrix(1:np,1:np),imatrix(1:np,1:np),0)
        elseif(L>0.5*Rb.and.np==6)then
            np = 7
            !print*,'Rb',np
            call inv_lie(mmatrix(1:np,1:np),imatrix(1:np,1:np),0)
        elseif(L>0.3*Rb.and.np==5)then
            np = 6
            !print*,'Rb',np
            call inv_lie(mmatrix(1:np,1:np),imatrix(1:np,1:np),0)
        elseif(L>0.2*Rb.and.np==4)then
            np = 5
            !print*,'Rb',np
            call inv_lie(mmatrix(1:np,1:np),imatrix(1:np,1:np),0)
        elseif(L>0.1*Rb.and.np==3)then
            np = 4
            !print*,'Rb',np
            call inv_lie(mmatrix(1:np,1:np),imatrix(1:np,1:np),0)
        elseif(L>0.05*Rb.and.np==2)then
            np = 3
            !print*,'Rb',np
            call inv_lie(mmatrix(1:np,1:np),imatrix(1:np,1:np),0)
        elseif(L>0.02*Rb.and.np==1)then
            np = 2
            !print*,'Rb',np
            call inv_lie(mmatrix(1:np,1:np),imatrix(1:np,1:np),0)
        endif
        ! update dt
        dt = L/(np+1)**2/(dRs+cs)
        if(mod(n,50)==0)   write(105,'(3E15.5)') time,Rb,0.0
        if(mod(n,1000)==0)then
                print*,'-------------------------------------------------'
                print'(A10,3A13)','increment','time','Rs','Rb'
                print*,'-------------------------------------------------'
        endif
        if(dt<0)stop
    enddo
    ! start DG simulation
    call DG_init()
    dt= 0.5*dt
    call output_dg(Rp)

    open(108,file='./output/dg.dat')
    write(108,'(5A15)') '#time', 'Rb', 'R1', 'R2', 'Rs'
    iexit = .false.
    do while(.true.)
        time=time+dt
        
        ! backup data
        ans_dg1=ans_dg
        Rb1=Rb
        Rs1=Rs
        Rdg1=Rdg
        ! flux
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
        endif
        
        if(mod(n,10)==0)then
            !print'(4E10.3,I3)',time,dt,Rb,Rs,np
            call output_dg(Rp)
        endif
        if(Rs>6*Rb.and..not.iexit)then
            dRb0 = equivalent_vel_dg()
            R0 = Rb
            pg0=pg
            t_start = time
            iexit = .true.
        endif
        if(time>2*t_start.and.iexit)then
            exit
        endif
        if(mod(n,100)==0)then
            print'(I10,3E13.3)',n,time,Rs,Rb
        endif
        if(mod(n,2000)==0)then
            print*,'-------------------------------------------------'
            print'(A10,3A13)','increment','time','Rs','Rb'
            print*,'-------------------------------------------------'
        endif
        if(mod(n,200)==0.and..not.iexit)   write(105,'(3E15.5)') time,Rb,0.0
    enddo
    
    close(100)
    close(101)
    close(102)
    close(103)
    !stop
    ! start far-field pressure prediction
    print*,'-------------------------------------------------'
    print*,'-------------------------------------------------'
    print*,'Finished near-field solver!'
    print*,'-------------------------------------------------'
    ! record the time of the end of the nearfield solver
    
    

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

    print*,'-------------------------------------------------'
    print*,'Start unified bubble dynamics solver!'
    print*,'-------------------------------------------------'
    print*,'-------------------------------------------------'
    call initiation(R0,dRb0,pg0)
    inc = 0
    t = 0
    t_shock = time - t_arrive
    print*,'-------------------------------------------------'
    print'(A10,3A13)','increment','time','Rb','Zb'
    print*,'-------------------------------------------------'
    do while(t + t_start - t_arrive<tend)
        call collect_dt()
        call bubbles(1)%advance
        !if(mod(inc,1)==0)&
        !    call bubbles(1)%print_bubble
        if(mod(inc,20)==0 .and. t + t_start - t_arrive > t_offset .and.&
            t > Rp/csound)then
            write(104,'(2E15.6)') t + t_start - t_arrive, &
                collect_induced_field(pressure_loc,t)
        endif
        if(mod(inc,500)==0)then
            print'(I10,3E13.3)',inc,t,bubbles(1)%R,bubbles(1)%center(3)
        endif
        if(mod(inc,20)==0)   write(105,'(3E15.5)') t+t_start,bubbles(1)%R,bubbles(1)%center(3)
        inc=inc+1
        t=t+dtb
    enddo
    close(104)
    close(105)
    print*,'Simulation Completed!'
    pause
    
    end program undex_ubt

