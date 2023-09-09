    ! in this module, the legendre polynomial with arbitrary 
    ! order is implemented.
    
    module legendre
    integer nmax
    integer,allocatable:: coeffs(:,:)
    
    type basis
        integer order
        real,allocatable:: a(:),da(:)   ! coefficients for the 
                                        ! polynomial and its derivative
    contains
    procedure initiate => basis_initiate
    procedure evaluate
    procedure derivative
    end type

    contains
    
    subroutine initiate(n)
    implicit none
    ! local
    integer i,n
    nmax=n
    allocate(coeffs(n,n))
    coeffs = 0
    coeffs(1,1) = 1
    coeffs(1:2,2) = [0.0,1.0]
    coeffs(1:3,3) = [0.5,0.0,1.5]
    
    do i=4,nmax
        n = i -1
        coeffs(1:n-1,i) = - n*coeffs(1:n-1,n-1)
        coeffs(2:i,i) = coeffs(2:i,i) + &
            (2*n+1)*coeffs(1:n,n)
    enddo
    end subroutine
    subroutine basis_initiate(this,n)
    implicit none
    class(basis):: this
    integer n
    ! local
    integer i
    ! code
    if(n>nmax)then
        print*,'Error: n must be smaller than nmax'
        stop
    endif
    
    allocate(this%a(n))
    this%a(1:n) = coeffs(1:n,n)
    do i=1,n-1
        this%da(i) = this%a(i+1)*i
    enddo
    end subroutine
    
    function evaluate(this,x) result (y)
    implicit none
    class(basis):: this
    real, intent(in):: x
    real y
    ! local
    integer i
    y = this%a(1)
    do i=2,this%order
        y = y + this%a(i)*x**i
    enddo
    return
    end function
    function derivative(this,x) result(dy)
    implicit none
    class(basis):: this
    real, intent(in):: x
    real dy
    ! local
    integer i
    dy = this%da(1)
    do i=2,this%order-1
        dy = dy+ this%da(i)*x**i
    enddo
    return
    end function
    end module