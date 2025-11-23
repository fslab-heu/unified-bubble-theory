submodule(fluid_eos)  fluid_eos_pressure
contains
   module function eos_pressure(this,dens,En,Fb) result(p)
      implicit none
      class(eos):: this
      ! input
      real dens, en
      real,optional:: Fb
      ! output
      real p
      real rho0
      !global sc
      real sc!,pstatic
      sc = 1e-12

      p = 0
      if(this%eos_type==0)then
         p=0
      elseif(this%eos_type==1)then
         p=En*(this%gamma-1)            &
            - this%gamma*this%Pw
      elseif(this%eos_type==2)then
         rho0=this%dens0/dens          !!!!rho0=1/��
         p=this%A*(1-this%W/this%R1/rho0)*exp(-rho0*this%R1)+&
            this%B*(1-this%W/this%R2/rho0)*exp(-rho0*this%R2)+&
            this%W*En
         ! detonation p
         if(this%iexpo.and.present(Fb))then
            if(Fb.le.sc)then
               ! explosive p equal to enviorment p
               p = this%pstatic
            else
               ! explosive p * burn fraction
               p = p * Fb
            endif
         endif
      else
         write(*,*) 'Unknown material Report from p'
      endif
      !if(p<sc) p=sc

      return
   end function

   module function eos_pressure_adb(this,P,dens_temp,dens_ifluid)
      implicit none
      class(eos):: this
      ! input
      real p,dens_temp,dens_ifluid
      ! output
      real eos_pressure_adb
      ! local
      real r_rho,e1,e2,c
      real sc
      sc = 1e-12
      
      if(this%eos_type==0)then
          eos_pressure_adb=0
      else
          r_rho=dens_ifluid/this%dens0
          if(this%eos_type==1)then
              eos_pressure_adb=(P+this%Pw)*r_rho**this%gamma-this%Pw
          elseif(this%eos_type==2)then
              e1=exp(-this%dens0/dens_temp*this%R1)
              e2=exp(-this%dens0/dens_temp*this%R2)
              C=(p-this%A*e1-this%B*e2)*(this%dens0/dens_temp)**(1+this%w)
              e1=exp(-this%dens0/dens_ifluid*this%R1)
              e2=exp(-this%dens0/dens_ifluid*this%R2)
              eos_pressure_adb=this%A*e1+this%B*e2           &
                 + C/(this%dens0/dens_ifluid)**(1+this%w)
          else
              write(*,*) 'Unknown material-adbp'
              stop
          endif
      endif
      return
      end function
end submodule
