    subroutine far_field_output(Rp)
    use polynomial
    use shockwave
    use global,only:t_arrive,tdelay
    implicit none
    integer i,n,j
    real press,dx,x,t,Rp
    
    x = (Rp-Re)/Lp
    press = 0
    do j=1,np
        press = press + basis(x,j)*pans(j)
    enddo
    write(104,'(2E15.6)'),time-t_arrive+tdelay,press
    end subroutine