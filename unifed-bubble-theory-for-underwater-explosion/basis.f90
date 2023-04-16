!   wrapper functions for basis of RKDGM
    
    function basis(x,n) result(res)
    use polynomial
    implicit none
    real res,x
    integer n
    res = leg_ply(n-1)%evaluate(x)
    return
    end function
    
    function dbasis(x,n) result(res)
    use polynomial
    implicit none
    real res,x
    integer n
    res = leg_ply(n-1)%derivative(x)
    return
    end function
    
    function basisn(x) result(res)
    use polynomial
    implicit none
    real x
    real res
    res = basis(x,s_nb)*basis(x,s_nb2)
    return
    end function