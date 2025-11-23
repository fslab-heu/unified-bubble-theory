    submodule (fluid_eos) func_entropy_sub
    implicit none
    contains
    module function entropy(this,p,dens)
    implicit none
    class(eos):: this
    ! input
    real p,dens
    ! output
    real entro2dens
    real temp,rho0
    real entropy

    rho0 = this%dens0/dens
    if(this%eos_type==0)then
        entropy=0
    elseif(this%eos_type==1)then
        entropy = (p + this%pw)/dens**this%gamma
    else        !---not verified JWL entropy
        entropy = this%e0 * this%w * this%dens0**(-this%w)
        write(*,*) 'Unknown material'
    endif
    return
    end function

    module function dens_from_entropy(this,entro, p)
    implicit none
    class(eos):: this
    ! input
    real entro, p
    ! output
    real dens_from_entropy

    if(this%eos_type==0)then
        dens_from_entropy=0
    elseif(this%eos_type==1)then
        dens_from_entropy = ((p+this%pw)/(entro))**(1.0/this%gamma)
    else        !---not verified JWL entropy
        ! entropy = this%e0 * this%w * this%dens0**(-this%w)
        ! write(*,*) 'Unknown material'
    endif
    return
    end function
    end submodule