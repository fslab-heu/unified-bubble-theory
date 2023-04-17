    subroutine write_restart()
    use shockwave
    implicit none
    open(1000,file='nearfield.res',status='new',form='unformatted')
    write(1000) time,dt,np,Rb,Rs,dRb,dRs,pg
    write(1000) ans,imatrix
    close(1000)
    
    
    end subroutine
    
    subroutine read_restart()
    use shockwave
    implicit none
    open(1000,file='nearfield.res',status='old',form='unformatted')
    read(1000) time,dt,np,Rb,Rs,dRb,dRs,pg
    read(1000) ans,imatrix
    close(1000)
    
    end subroutine