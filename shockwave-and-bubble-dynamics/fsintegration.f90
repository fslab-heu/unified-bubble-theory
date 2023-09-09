    module fsintegration
    
    type gauss_lobatto
        integer npoints
        real,allocatable:: pos(:)
        real,allocatable:: weight(:)
    contains
    procedure gl_init 
    end type
    
    type gauss_legendre
        integer npoints
        real,allocatable:: pos(:)
        real,allocatable:: weight(:)
    contains
    procedure gl_init_legendre
    end type
    
    type(gauss_legendre) gl0
    
    contains
    
    subroutine gl_init(this,npoints)
    implicit none
    class(gauss_lobatto)::this
    ! input
    integer, intent(in):: npoints
    this%npoints=npoints
    allocate(this%pos(npoints),this%weight(npoints))
    if(npoints==3)then
        this%pos = [-1,0,1]
        this%weight = [1,4,1]/3.0
    elseif(npoints==4)then
        this%pos = [-1.0,-0.2*sqrt(5.0),0.2*sqrt(5.0),1.0]
        this%weight = [1,5,5,1]/6.0
    elseif(npoints==5)then
        this%pos = [-1.0,-sqrt(21.0)/7.0,0.0,sqrt(21.0)/7.0,1.0]
        this%weight = [0.1,49.0/90.0,32.0/45.0,49.0/90.0,0.1]
    else
        print*,'Error: Invalid npoints for Gauss Lobatto quadrature!'
        stop
    endif

    return
    end subroutine
    
    function adaptive_integrate(xstart,xend,func,ecr,debug) result(res)
    implicit none
    ! input
    real xstart,xend
    real,external:: func
    real,optional:: ecr
    logical,optional:: debug
    ! output
    real res
    ! local
    real ecr0,dx,x(5),y(5),yend,err
    real area,area_sub,area_sup,area_bak
    integer n,m,i,k
    
    ! initialize parameters
    if(present(ecr))then
        ecr0 = ecr
    else
        ecr0 = 1e-9
    endif
    
    
    dx =(xend - xstart)/5;
    
    x = xstart + dx*&
        [0.0,0.25,0.5,0.75,1.0]
    do i=1,5
        y(i)=func(x(i))
    enddo
    yend = func(xend)
    res=0
    n = 0
    m = 0
    area = 0
    if(present(debug))then
        if(debug) print'(A6,A6,A10,A10)','inc','iter','dx','err'
    endif
    
    do while(.true.)
        area_bak = area
        k=1
        do while(.true.)
            area = dx*(y(1) + 4*y(3) + y(5))/6.0
            area_sub = dx*(y(1) + 4*y(2) + y(3))/12.0 + &
                dx*(y(3) + 4*y(4) + y(5))/12.0
            err = abs(area - area_sub)
            if(err>max(ecr0,ecr0*abs(area)))then
                dx = dx/2.0
                x = x(1) + dx*&
                    [0.0,0.25,0.5,0.75,1.0]
                do i=1,5
                    y(i)=func(x(i))
                enddo
                m = 0
                k = k + 1
            else
                m = m + 1
                exit
            endif
        enddo
        n = n + 1
        res = res + area
        if(present(debug))then
            if(debug)print'(I6,I6,E10.3,E10.3)',n,k,dx,err
        endif
        
        x(1) = x(5)
        ! determine whether to enlarge the the next step
        if(m >= 2)then
            area_sup = 2*dx*(func(x(1)-2*dx)+4.0*y(1)+y(5))/6.0
            if(abs(area_sup - area - area_bak)<=max(ecr0,ecr0*abs(area_sup)))then
                dx = dx*2.0
            endif
        endif
        
        ! if the next step exceeds xend, limit the dx
        if(abs((xend-x(1))*yend)<max(ecr0,ecr0*abs(res)))then
            exit
        elseif(x(1)+dx>xend)then
            dx = xend - x(1)
        endif
        x = x(1) + dx*&
            [0.0,0.25,0.5,0.75,1.0]
        do i=1,5
            y(i)=func(x(i))
        enddo
    enddo
    end function
    
    
    function adaptive_integrate_gl(xstart,xend,func,ecr,debug,np) result(res)
    ! use 5 points Gauss Labatoo integration for each segment
    implicit none
    ! input
    real xstart,xend
    real,external:: func
    real,optional:: ecr
    logical,optional:: debug
    integer,optional:: np
    ! output
    real res
    ! local
    real ecr0,dx,x,y,yend,err
    real area,area_sub,area_sup,area_bak
    integer n,m,i,k,np0

    ! initialize parameters
    if(present(ecr))then
        ecr0 = ecr
    else
        ecr0 = 1e-9
    endif
    
    if(present(np))then
        np0 = np
    else
        np0 = 5
    endif
    
       
    dx =(xend - xstart)/5.0;
    
    x= xstart 
    
    yend = func(xend)
    res=0
    n = 0
    m = 0
    area = 0
    if(present(debug))then
        if(debug) print'(A6,A6,A10,A10)','inc','iter','dx','err'
    endif
    
    do while(.true.)
        area_bak = area
        k=1
        do while(.true.)
            area = quadrature_gl(x,x+dx,func,np0)
            area_sub = quadrature_gl(x,x+0.5*dx,func,np0) + quadrature_gl(x+0.5*dx,x+dx,func,np0)
            err = abs(area - area_sub)
            if(err>max(ecr0,ecr0*abs(area)))then
                dx = dx/2.0
                m = 0
                k = k + 1
            else
                m = m + 1
                exit
            endif
        enddo
        n = n + 1
        res = res + area
        if(present(debug))then
            if(debug)print'(I6,I6,E10.3,E10.3)',n,k,dx,err
        endif
        
        x=x+dx
        ! determine whether to enlarge the the next step
        if(m >= 2)then
            area_sup = quadrature_gl(x-2*dx,x,func,np0)
            if(abs(area_sup - area - area_bak)<=max(ecr0,ecr0*abs(area_sup)))then
                dx = dx*2.0
            endif
        endif
        
        ! if the next step exceeds xend, limit the dx
        if(abs((xend-x)*yend)<max(ecr0,ecr0*abs(res)))then
            exit
        elseif(x+dx>xend)then
            dx = xend - x
        endif
    enddo
    end function
    
    function quadrature_gl(x1,x2,func,np) result(res)
    implicit none
    real x1,x2
    real,external:: func
    integer np
    real res
    ! local
    type(gauss_lobatto)::this
    real L,x,xc
    integer i
    real,allocatable,dimension(:):: xg,wg
    
    allocate (xg(np))
    allocate (wg(np))
    call this%gl_init(np)
    xg=this%pos
    wg=this%weight
       
    L = x2 - x1                
    xc = (x1+x2)/2.0
    
    res = 0
    do i=1,np
        x = xc+ xg(i)*L/2.0
        res = res + wg(i)*func(x)
    enddo
    res = res*L/2.0
    deallocate(xg,wg)
    return
    end function
    
    
    subroutine gl_init_legendre(this,npoints)
    implicit none
    class(gauss_legendre)::this
    ! input
    integer, intent(in):: npoints
    this%npoints=npoints
    allocate(this%pos(npoints),this%weight(npoints))
    if(npoints==3)then
        this%pos = [-sqrt(3.0/5.0),sqrt(3.0/5.0),0.0]
        this%weight = [5.0/9.0,5.0/9.0,8.0/9.0]
    elseif(npoints==4)then
        this%pos = [sqrt(3.0/7.0-2.0/7.0*sqrt(6.0/5.0)),sqrt(3.0/7.0+2.0/7.0*sqrt(6.0/5.0)),&
                   -sqrt(3.0/7.0-2.0/7.0*sqrt(6.0/5.0)),-sqrt(3.0/7.0+2.0/7.0*sqrt(6.0/5.0))]
        this%weight = [0.5+sqrt(30.0)/36.0,0.5-sqrt(30.0)/36.0,&
                       0.5+sqrt(30.0)/36.0,0.5-sqrt(30.0)/36.0]
    elseif(npoints==5)then
        this%pos = [0.0,sqrt(245.0-14.0*sqrt(70.0))/21.0,-sqrt(245.0-14.0*sqrt(70.0))/21.0,&
                    sqrt(245.0+14.0*sqrt(70.0))/21.0,-sqrt(245.0+14.0*sqrt(70.0))/21.0]
        this%weight = [128.0/225.0,(322.0+13.0*sqrt(70.0))/900.0,(322.0+13.0*sqrt(70.0))/900.0,&
                      (322.0-13.0*sqrt(70.0))/900.0,(322.0-13.0*sqrt(70.0))/900.0]
    elseif(npoints==2)then
        this%pos = [1.0/sqrt(3.0),-1.0/sqrt(3.0)]
        this%weight = [1.0,1.0]
    else
        print*,'Error: Invalid npoints for Gauss legendre quadrature!'
        stop
    endif
    
    end  subroutine
    
    
    function adaptive_integrate_legendre(xstart,xend,func,ecr,debug,np) result(res)
    ! use 5 points Gauss Labatoo integration for each segment
    implicit none
    ! input
    real xstart,xend
    real,external:: func
    real,optional:: ecr
    logical,optional:: debug
    integer,optional:: np
    ! output
    real res
    ! local
    real ecr0,dx,x,y,yend,err
    real area,area_sub,area_sup,area_bak
    integer n,m,i,k,np0

    ! initialize parameters
    if(present(ecr))then
        ecr0 = ecr
    else
        ecr0 = 1e-9
    endif

    if(present(np))then
        np0 = np
    else
        np0 = 5
    endif
    
    
    dx =(xend - xstart)/5.0;
    
    x= xstart 
    
    yend = func(xend)
    res=0
    n = 0
    m = 0
    area = 0
    if(present(debug))then
        if(debug) print'(A6,A6,A10,A10)','inc','iter','dx','err'
    endif
    
    do while(.true.)
        area_bak = area
        k=1
        do while(.true.)
            area = quadrature_legendre(x,x+dx,func,np0)
            area_sub = quadrature_legendre(x,x+0.5*dx,func,np0) +quadrature_legendre(x+0.5*dx,x+dx,func,np0)
            err = abs(area - area_sub)
            if(err>max(ecr0,ecr0*abs(area)))then
                dx = dx/2.0
                m = 0
                k = k + 1
            else
                m = m + 1
                exit
            endif
        enddo
        n = n + 1
        res = res + area
        if(present(debug))then
            if(debug)print'(I6,I6,E10.3,E10.3)',n,k,dx,err
        endif
        
        x=x+dx
        ! determine whether to enlarge the the next step
        if(m >= 2)then
            area_sup = quadrature_legendre(x-2*dx,x,func,np0)
            if(abs(area_sup - area - area_bak)<=max(ecr0,ecr0*abs(area_sup)))then
                dx = dx*2.0
            endif
        endif
        
        ! if the next step exceeds xend, limit the dx
        if(abs((xend-x)*yend)<max(ecr0,ecr0*abs(res)))then
            exit
        elseif(x+dx>xend)then
            dx = xend - x
        endif
    enddo
    end function
    
    function quadrature_legendre(x1,x2,func,np) result(res)
    implicit none
    real x1,x2
    real,external:: func
    integer np
    real res
    ! local
    type(gauss_legendre)::this
    real L,x,xc
    integer i
    real,allocatable,dimension(:):: xg,wg
    
    allocate (xg(np))
    allocate (wg(np))
    call this%gl_init_legendre(np)                                                                                         
    xg=this%pos
    wg=this%weight
       
    L = x2 - x1
    xc = (x1+x2)/2.0
    
    res = 0
    do i=1,np
        x = xc+ xg(i)*L/2.0
        res = res + wg(i)*func(x)
    enddo
    res = res*L/2.0
    deallocate(xg,wg)
    return
    end function
    
    end module