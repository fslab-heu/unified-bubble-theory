    subroutine output_DG(Rp)
    use shockwave
    !use math
    use global,only: iarrive,t_arrive,tdelay
    implicit none
    integer n,i,j
    real dx,s,rho,ein,mom
    real,allocatable:: x(:),pressure(:),entropy(:)
    real,external:: energy
    real pr,Rp,cv(3),zetar,ur,zeta
    character*100 fmt
    logical, save:: touched = .false.
    n = 100
    dx = L/n
    write(fmt, "('(E13.5,',I3,'E13.5)')") n+1
    allocate(x(n+1),pressure(n+1),entropy(n+1))
    do i=1,n+1
        x(i)=(i-1)*dx+Rb
        do j=1,ndg
            if(x(i)>=Rdg(ndg+1-j))then
                exit
            endif
        enddo
        s_nb = ndg-j+1    ! element id
        zeta = (x(i)-Rdg(s_nb))/Ldg(s_nb)
        
        
        rho = 0
        s = 0
        mom =0
        do j=1,np
            s =   s +   ans_dg(j,3,s_nb)*basis(zeta,j)
            rho = rho + ans_dg(j,1,s_nb)*basis(zeta,j)
            mom = mom + ans_dg(j,2,s_nb)*basis(zeta,j)
        enddo
        pressure(i) = e2p(s-0.5*mom**2/rho)
    enddo
    
    if(Rp<Rb)then
        pr = pg
        ur = 0
    elseif(Rp>Rs)then
        pr = p0
        ur = 0
    else
        !if(.not.touched)then
        !    write(102,'(8E15.7)') time, Rb, Rs, dRb, equivalent_vel_dg(), pressure(n+1),pg,p0
        !    touched = .true.
        !endif
        
        cv = 0
        do j=1,ndg
            if(Rp>=Rdg(ndg+1-j))then
                exit
            endif
        enddo
        s_nb = ndg-j+1    ! element id
        zetar = (Rp-Rdg(s_nb))/Ldg(s_nb)
        
        do j=1,np
            cv = cv + ans_dg(j,:,s_nb)*basis(zetar,j)
        enddo
        ur = cv(2)/cv(1)
        pr = e2p(cv(3)-0.5*cv(2)**2/cv(1))
    endif
    if(iarrive)then
        write(104,'(2E15.6)') time-t_arrive+tdelay,pr - p0
    endif
    end subroutine