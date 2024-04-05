    ! module to hold the boundary information
    ! in this code, it refers to the freesurface
    module boundary
    private :: length
    type bound
        real norm(3)
        real b
        real alpha
    contains
    procedure:: ini_bound
    end type
    
    type(bound):: bound1
    contains
    subroutine ini_bound(this,norm,pos,alpha)
    implicit none
    class(bound):: this
    real norm(3),pos(3),alpha
    this%alpha = alpha
    this%norm = norm/length(norm)
    this%b = -dot_product(this%norm,pos)
    return
    end subroutine
    real function length(vect)
        real vect(:)
        integer N
        N=size(vect)
        length=sqrt(sum(vect**2))
        return
    end function
    
    end module