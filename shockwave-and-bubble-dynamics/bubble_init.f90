    subroutine initiation(r,dr,p)
    ! initiate the bubble simulation
    use bubble
    use global
    use boundary
    implicit none
    integer i
    ! input
    real r,dr,p
    ! code
    
    pi = 2.0*acos(0.0)
    
    ! initiate the bubble
    call bubbles(1)%ini_bubble(0.0,[0.0,0.0,0.0],r,dr,p,1)
    ! calculate dt
    call collect_dt()
    
    end subroutine