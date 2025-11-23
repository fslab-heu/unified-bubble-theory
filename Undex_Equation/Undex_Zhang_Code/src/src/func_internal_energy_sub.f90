submodule(fluid_eos) fluid_eos_interal_energy
contains
    module function eos_internal_energy(this,pressure,density)! result(ein)
    implicit none
    class(eos):: this
    ! input
    real pressure, density
    ! output
    real eos_internal_energy
    real rho0,e1,e2
    
    if(this%eos_type==0)then
        eos_internal_energy=0
    elseif(this%eos_type==1)then
        eos_internal_energy=(pressure+this%gamma*this%Pw)/(this%gamma-1)
    elseif(this%eos_type==2)then
        rho0=this%dens0/density
        e1=exp(-rho0*this%R1);e2=exp(-rho0*this%R2)
        eos_internal_energy=(pressure-this%A*(1-this%w/rho0/this%R1)*e1-     &
            this%B*(1-this%w/rho0/this%R2)*e2)/this%w
    else
        write(*,*) 'Unknown material Report from energe'
        stop
    endif    
    
    return
    end function
end submodule
