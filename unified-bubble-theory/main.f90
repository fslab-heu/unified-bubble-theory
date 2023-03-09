    program main
    use bubble
    use global
    implicit none
    integer i
    
    path='output'
    call initiation()
    inc=1
    
    open(100,file=trim(path)//'/pressure.dat')
    write(100,'(12A15)') '#time' , 'pressure',      &
        'dp','vx','vy','vz','px','py','pz','dvx',   &
        'dvy','dvz'
    
    do while(t<tend)
        call collect_dt()

        do i=1,nbubble
            call bubbles(i)%advance
            if(mod(inc,10)==0)&
                call bubbles(i)%print_bubble
        enddo
        if(mod(inc,100)==0)then
            
            write(100,'(12E15.6)') t,               &
                collect_induced_field(pressure_loc,t)
        endif
        if(mod(inc,500)==0)then
            print'(I8,E12.4)',inc,t
        endif
        
        if(mod(inc,10000)==1)then
            print*,'-------------------------'
            print'(A8,A12)','inc','time'
            print*,'-------------------------'
        endif
        
        inc=inc+1
        t=t+dt
    enddo
    print*,'Simulation Completed!'
    end program