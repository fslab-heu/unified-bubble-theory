submodule(fluid_eos) fluid_eos_cspeed
contains
   module function eos_sound_speed(this,dens,p) result(c)
      implicit none
      class(eos):: this
      ! input
      real dens,p
      ! output
      real c
      ! local
      real Ce,rho0,e1,e2
      real sc
      sc = 1e-12
      c=-1
      if(this%eos_type==0)then
         c=0
      elseif(this%eos_type==1)then
         if(dens>sc)then
            c=(P+this%Pw)*this%gamma/dens
         endif
      elseif(this%eos_type==2)then
         rho0=this%dens0/dens
         e1=exp(-rho0*this%R1)
         e2=exp(-rho0*this%R2)
         Ce=(P-this%A*e1-this%B*e2)*rho0**(1+this%w)
         c = Ce/rho0**this%w*(1+this%w)/this%dens0+&
            (this%A*this%R1*e1+this%B*this%R2*e2)*rho0**2/this%dens0
      else
         write(*,*) 'Unknown material Report from csound'
      endif
      if(c<0)then
         c=-1
      else
         c=sqrt(c)
      endif
   end function


   module function eos_sound_speed_2(this,dens,p) result(c2)
      ! returns the square of the sound speed
      ! used to check for parabolicity, when c^2 < 0
      implicit none
      class(eos):: this
      ! input
      real dens,p
      ! output
      real c2
      ! local
      real Ce,rho0,e1,e2
      real sc
      sc = 1e-12
      if(this%eos_type==0)then   ! void
         c2=0
      elseif(this%eos_type==1)then  
         if(dens>sc)then
            c2=(P+this%Pw)*this%gamma/dens
         endif
      elseif(this%eos_type==2)then
         rho0=this%dens0/dens
         e1=exp(-rho0*this%R1)
         e2=exp(-rho0*this%R2)
         Ce=(P-this%A*e1-this%B*e2)*rho0**(1+this%w)
         c2 = Ce/rho0**this%w*(1+this%w)/this%dens0+&
            (this%A*this%R1*e1+this%B*this%R2*e2)*rho0**2/this%dens0
      else
         write(*,*) 'Unknown material Report from csound'
      endif
   end function

end submodule

