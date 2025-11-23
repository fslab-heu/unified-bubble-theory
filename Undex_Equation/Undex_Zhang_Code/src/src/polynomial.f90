    module polynomial
    ! in this module, the legendre polynomial with arbitrary 
    ! order is implemented.
    ! The initiation subroutine should be firstly called
    ! then, the instance leg_ply can be used to evaluate the
    ! polynomial and its derivative
    
    ! call initiate_legendre(m)
    ! res = leg_ply(m)%evaluate(x)
    ! res = leg_ply(m)%derivative(x)
    integer nmax
    real,private:: xs, xe, length
    real,private, allocatable:: coeffs(:,:)
    
    logical,private:: initiated
    
    type legendre
        integer order
        real,allocatable:: ia(:), &     ! integration coefficients, from -1 to zeta
            a(:),&                      ! coefficients
            da(:),dda(:),ddda(:)        ! 1st, 2nd and 3rd derivative
    contains
    procedure :: initiate => basis_initiate
    procedure :: integrate
    procedure :: evaluate
    procedure :: derivative
    procedure :: second_derivative
    procedure :: third_derivative
    end type
    
    class(legendre),pointer:: leg_ply(:)
    contains
    
    subroutine initiate_legendre(m,x1,x2)
    implicit none
    integer,intent(in):: m
    real,optional:: x1,x2   ! these are the range of the ply
                            ! default value is -1:1
    ! local
    integer i,n
    nmax=m
    initiated = .true.
    allocate(coeffs(0:max(m,2),0:max(m,2)))
    coeffs = 0
    coeffs(0,0) = 1
    coeffs(0:1,1) = [0.0,1.0]
    coeffs(0:2,2) = [-0.5,0.0,1.5]
    
    ! define the mapping function
    if(present(x1) .neqv. present(x2))then
        print*,'Error: x1 and x2 must be given at the same time!'
        stop
    endif
    
    if(present(x1).and.present(x2))then
        xs = x1
        xe = x2
    else
        xs = -1
        xe = 1
    endif
    
    length = xe - xs
    
    do n=2,m-1
        coeffs(1:n+1,n+1) = (2*n+1)*coeffs(0:n,n)
        coeffs(0:n-1,n+1)=coeffs(0:n-1,n+1)-&
            n*coeffs(0:n-1,n-1)
        coeffs(:,n+1) = coeffs(:,n+1)/(n+1)
    enddo
    
    allocate(leg_ply(0:m))
    do i=0,m
        call leg_ply(i)%initiate(i)    
    enddo
    
    
    end subroutine
    subroutine basis_initiate(this,n)
    implicit none
    class(legendre):: this
    integer n
    ! local
    integer i
    ! code
    if(n>nmax)then
        print*,'Error: n must be smaller than nmax'
        stop
    endif
    this%order = n
    allocate(this%ia(0:n+1),this%a(0:n),this%da(0:n-1),this%dda(0:n-2),this%ddda(0:n-3))
    this%a(0:n) = coeffs(0:n,n)
    do i=1,n+1
        this%ia(i) = this%a(i-1)/(i+2)
    enddo
    do i=0,n-1
        this%da(i) = this%a(i+1)*(i+1)
    enddo
    do i=0,n-2
        this%dda(i) = this%da(i+1)*(i+1)
    enddo
    do i=0,n-3
        this%ddda(i) = this%dda(i+1)*(i+1)
    enddo
    end subroutine
    
    function evaluate(this,x0) result (y)
    implicit none
    class(legendre):: this
    real, intent(in):: x0
    real y
    ! local
    integer i
    real x
    if(.not.initiated)then
        print*,'Error: the polynomial module used without initiation!'
        print*,'Please call basis_initiate fistly!'
        stop
    endif
    
    x = (x0 - xs)/length*2.0 - 1.0
    y = this%a(0)
    do i=1,this%order
        y = y + this%a(i)*x**i
    enddo
    return
    end function

    function integrate(this,x0) result (y)
    implicit none
    class(legendre):: this
    real, intent(in):: x0
    real y
    ! local
    integer i
    real x
    if(.not.initiated)then
        print*,'Error: the polynomial module used without initiation!'
        print*,'Please call basis_initiate fistly!'
        stop
    endif
    
    x = (x0 - xs)/length*2.0 - 1.0
    y = this%ia(0)
    do i=1,this%order
        y = y + this%ia(i)*x**i
    enddo
    return
    end function

    function derivative(this,x0) result(dy)
    implicit none
    class(legendre):: this
    real, intent(in):: x0
    real dy
    ! local
    integer i
    real x
    if(.not.initiated)then
        print*,'Error: the polynomial module used without initiation!'
        print*,'Please call basis_initiate fistly!'
        stop
    endif
    x = (x0 - xs)/length*2.0 - 1.0
    
    if(this%order==0)then
        dy=0
        return
    endif
    
    dy = this%da(0)
    do i=1,this%order-1
        dy = dy+ this%da(i)*x**i
    enddo
    dy = dy/length*2.0
    return
    end function
    
    function second_derivative(this,x0) result(ddy)
    implicit none
    class(legendre):: this
    real, intent(in):: x0
    real ddy
    ! local
    integer i
    real x
    if(.not.initiated)then
        print*,'Error: the polynomial module used without initiation!'
        print*,'Please call basis_initiate fistly!'
        stop
    endif
    x = (x0 - xs)/length*2.0 - 1.0
    
    if(this%order<=1)then
        ddy=0
        return
    endif
    
    ddy = this%dda(0)
    do i=1,this%order-2
        ddy = ddy+ this%dda(i)*x**i
    enddo
    ddy = ddy/(length**2)*(2.0**2)
    return
    end function
    
    function third_derivative(this,x0) result(dddy)
    implicit none
    class(legendre):: this
    real, intent(in):: x0
    real dddy
    ! local
    integer i
    real x
    if(.not.initiated)then
        print*,'Error: the polynomial module used without initiation!'
        print*,'Please call basis_initiate fistly!'
        stop
    endif
    x = (x0 - xs)/length*2.0 - 1.0
    
    if(this%order<=2)then
        dddy=0
        return
    endif
    
    dddy = this%ddda(0)
    do i=1,this%order-3
        dddy = dddy+ this%ddda(i)*x**i
    enddo
    dddy = dddy/(length**3)*(2.0**3)
    return
    end function
    
    end module