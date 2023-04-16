    ! read case parameters and initialize the programe
    subroutine initiation()
    use bubble
    use global
    use boundary
    implicit none
    integer i
    character*100 fname
    real pos(3),r,dr,p,norm(3),alpha,ts
    
    pi = 2.0*acos(0.0)
    
    ! constants
    fname = './constants.in'
    
    open(100,file=fname,status = 'old')
    read(100,*) gamma
    read(100,*) c
    read(100,*) miu
    read(100,*) sigma
    read(100,*) Pv
    read(100,*) g
    read(100,*) rho
    read(100,*) cd
    read(100,*) ca
    read(100,*) pamb
    read(100,*) uamb
    
    close(100)
    
    
    ! case parameters
    fname ='./case.in'
    open(100,file=fname,status = 'old')
    read(100,*) tend
    read(100,*) sdt
    read(100,*) pressure_loc
    read(100,*) ibound
    read(100,*) pos
    read(100,*) norm
    read(100,*) alpha
    
    call bound1%ini_bound(norm,pos,alpha)
    read(100,*) imigration
    read(100,*) nbubble
    do i=1,nbubble
        read(100,*) ts,pos,r,dr,p
        call bubbles(i)%ini_bubble(ts,pos,r,dr,p,i)
    enddo
    close(100)
    
    ! set the initial time increment
    call collect_dt()
    
    end subroutine