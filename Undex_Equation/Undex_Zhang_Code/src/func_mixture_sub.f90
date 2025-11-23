submodule(fluid_eos) fluid_eos_mixture
contains
    module subroutine eos_mix_init(this,eos1,eos2)
        implicit none
        class(eos_mix):: this
        ! input
        class(eos), target:: eos1, eos2

        this%eos1 => eos1
        this%eos2 => eos2
        ! check for thermal parameters initiation
        if(this%eos1%Cv==0.0 .or. this%eos2%Cv==0.0) then
            write(*,*) 'Error: One of the mixture components has Cv=0.0'
            error stop
        end if

    end subroutine

    module function eos_mix_pressure(this,g1, g2, en) result(p)
        implicit none
        class(eos_mix):: this
        ! input
        real g1, g2, en
        ! output
        real p
        ! local
        real R1, R2
        real e1, e2
        real D1, DT, D2
        real P1, P2
        ! code
        P1 = this%eos1%pw*this%eos1%gamma
        P2 = this%eos2%pw*this%eos2%gamma
        D1 = g1*this%eos1%Cv+g2*this%eos2%Cv
        DT = this%eos1%Cv*this%eos1%T0 - this%eos2%Cv*this%eos2%T0
        e1 = (g1+g2)*en*this%eos1%Cv/D1 - g2*DT*this%eos2%Cv/D1
        e2 = (g1+g2)*en*this%eos2%Cv/D1 + g1*DT*this%eos1%Cv/D1
        R1 = g1*e1*(this%eos1%gamma - 1.0)
        R2 = g2*e2*(this%eos2%gamma - 1.0)
        D2 = sqrt((R1-R2 - P1+P2)**2 + 4.0*R1*R2)
        p = 0.5*(R1 + R2 - P1 - P2 + D2)
    end function

    module function eos_mix_sound_speed(this,g1, g2, en) result(c)
        implicit none
        class(eos_mix):: this
        ! input
        real g1, g2, en
        ! output
        real c
        ! local
        real R1, R2
        real e1, e2
        real D1, DT, D2
        real P1, P2, p
        ! code
        P1 = this%eos1%pw*this%eos1%gamma
        P2 = this%eos2%pw*this%eos2%gamma
        D1 = g1*this%eos1%Cv+g2*this%eos2%Cv
        DT = this%eos1%Cv*this%eos1%T0 - this%eos2%Cv*this%eos2%T0
        e1 = (g1+g2)*en*this%eos1%Cv/D1 - g2*DT*this%eos2%Cv/D1
        e2 = (g1+g2)*en*this%eos2%Cv/D1 + g1*DT*this%eos1%Cv/D1
        R1 = g1*e1*(this%eos1%gamma - 1.0)
        R2 = g2*e2*(this%eos2%gamma - 1.0)
        D2 = sqrt((R1-R2 - P1+P2)**2 + 4.0*R1*R2)
        p = 0.5*(R1 + R2 - P1 - P2 + D2)
        
        c = g1*(this%eos1%gamma-1)*(p+P2)*(e1+ this%eos1%Cv/D1*p) + &
            g2*(this%eos2%gamma-1)*(p+P1)*(e2+ this%eos2%Cv/D1*p)
        if(c.lt.0.0) then
            write(*,*) 'Error: Negative sound speed squared in eos_mix_sound_speed'
            error stop
        end if
        c = sqrt(c/(g1+g2)/D1)
        
    end function
end submodule
