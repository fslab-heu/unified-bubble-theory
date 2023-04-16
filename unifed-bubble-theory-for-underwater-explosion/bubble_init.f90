    subroutine initiation(r,dr,p)
    use bubble
    use global
    use boundary
    implicit none
    integer i
    real pos(3),r,dr,p,ts
    
    pi = 2.0*acos(0.0)
    ts = 0
    pos = 0
    call bubbles(1)%ini_bubble(0.0,pos,r,dr,p,1)
    call collect_dt()
    
    end subroutine