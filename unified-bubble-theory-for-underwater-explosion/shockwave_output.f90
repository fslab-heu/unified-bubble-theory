    subroutine output(Rp)
    use shockwave
    use math
    implicit none
    integer n,i,j
    real dx,s,rho,ein,mom
    real,allocatable:: x(:),zeta(:),pressure(:),entropy(:)
    real,external:: energy
    real pr,Rp,cv(3),zetar,ur
    character*100 fmt
    logical, save:: touched = .false.
    n = 100
    dx = 1.0/n
    write(fmt, "('(E13.5,',I3,'E13.5)')") n+1
    allocate(x(n+1),zeta(n+1),pressure(n+1),entropy(n+1))
    do i=1,n+1
        zeta(i) = (i-1)*dx
        rho = 0
        s = 0
        mom =0
        do j=1,np
            s = s + ans(j,3)*basis(zeta(i),j)
            rho = rho + ans(j,1)*basis(zeta(i),j)
            mom = mom + ans(j,2)*basis(zeta(i),j)
        enddo
        pressure(i) = e2p(s-0.5*mom**2/rho)
        x(i) = zeta(i)*L+Rb
        entropy(i)=s
    enddo
    continue
    
    if(Rp<Rb)then
        pr = pg
        ur = 0
    elseif(Rp>Rs)then
        pr = p0
        ur = 0
    else
        if(.not.touched)then
            !write(102,'(8E15.7)') time, Rb, Rs, dRb, equivalent_vel_dg(), pressure(n+1),pg,p0
            touched = .true.
        endif
        cv = 0
        zetar = (Rp - Rb)/L
        do j=1,np
            cv = cv + ans(j,:)*basis(zetar,j)
        enddo
        ur = cv(2)/cv(1)
        pr = e2p(cv(3)-0.5*cv(2)**2/cv(1))
    endif
    
    !write(100,trim(fmt)) time,pressure
    !write(101,trim(fmt)) time,x
    !write(102,'(8E15.7)') time, Rb, Rs, dRb, equivalent_vel(), pressure(n+1),pg,Pr
    end subroutine
    
    
    function energy(x) result(res)
    use shockwave
    implicit none    
    real r,res,x
    integer i
    real cv
    r = x*L+Rb
    cv = 0
    do i=1,np
        cv = cv + ans(i,1)*basis(x,i)
    enddo
    res = cv*4*pi*r**2*L
    
    end function