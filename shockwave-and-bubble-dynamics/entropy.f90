    function e2p(e) result (pressure)
    implicit none
    real e,  pressure
    
    pressure = e*(gamma-1) - gamma*pw
    
    return
    end function
    
    function p2e(pressure) result (e)
    implicit none
    real e,  pressure
    
    e = (pressure+ gamma*pw)/(gamma-1)
    
    return
    end function
    
    function soundspeed(pressure,rho) result (c)
    implicit none
    real rho,  pressure,c
    real root
    root = (pressure+pw)/rho*gamma
    if(root<=0)then
        c=-1
    else
        c = sqrt(root)
    endif
    

    return
    end function
    
    function entro2pressure(entro, rho) result(pressure)
    implicit none
    real entro, rho
    real pressure
    
    pressure = entro*(rho/1000)**gamma -pw/1e5
    pressure = pressure*1e5
    end function
    
    function entro2rho(entro, pressure) result(rho)
    implicit none
    real entro, rho
    real pressure
    
    rho = ((pressure+pw)/1e5/entro)**(1/gamma)
    rho = rho * 1000
    end function
    
    function entroCal(rho, pressure) result(entro)
    implicit none
    real rho, pressure, entro
    entro = (pressure+pw)/1e5/(rho/1000)**gamma
    
    end function
    
    function entroDp(entro, entrox, rho, rhox) result(dp)
    implicit none
    real entro, entrox, rho, rhox
    real dp
    dp = entrox*(rho/1000)**gamma + &
        entro*(rho/1000)**(gamma-1)*rhox/1000
    dp = dp*1e5
    
    end function
    