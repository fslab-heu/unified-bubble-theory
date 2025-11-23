submodule(fluid_eos)  fluid_eos_init
contains
   module subroutine eos_init(this,eos_type,data_in,iexpo0, thermal_in)
    ! eos_type: 1 for tammann , 2 for JWL
    ! data_in : tammann: gamma, pw, dens0
    !           JWL: A, B, R1, R2, w, dens0, e0
      implicit none
      class(eos):: this
      ! input
      integer eos_type
      real data_in(:)
      ! local
      real e1, e2
      logical,optional::iexpo0
      real, optional:: thermal_in(:)


      !------- code
      this%eos_type = eos_type

      if(eos_type==1)then             ! tamman eos
         this%gamma = data_in(1)
         this%pw = data_in(2)
         this%dens0 = data_in(3)
         if(present(thermal_in)) then
            this%Cv = thermal_in(1)
            this%T0 = thermal_in(2)
         else
            this%Cv = 0.0
            this%T0 = 0.0
         end if
         !this%entropy = (1.0+this%B)/this%R1**this%A
      elseif(eos_type==2)then         ! jwl eos
         this%A =    data_in(1) !* k_jwl
         this%B =    data_in(2) !* k_jwl
         this%R1 =   data_in(3)
         this%R2 =   data_in(4)
         this%w =    data_in(5)
         this%dens0 =data_in(6)
         this%e0 =   data_in(7)

         e1 = exp(-this%R1)
         e2 = exp(-this%R2)

         this%C = -this%w/this%R1*this%A*e1 - &
            this%w/this%R2*this%B*e2+this%w*this%e0
         !this%entropy = this%e0 * this%w * this%dens0**(-this%w)

         ! for jwl, add detonation switch: iexpo (default set to false)
         this%iexpo = .false.!
         if(present(iexpo0)) this%iexpo = iexpo0

      else
         print*,'Error: unknown fluid type:',eos_type
         stop
      endif
   end subroutine
end submodule
