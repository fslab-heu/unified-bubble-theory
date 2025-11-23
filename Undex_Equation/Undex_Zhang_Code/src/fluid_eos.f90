    module fluid_eos
    
    type:: eos
        integer eos_type        ! 1 for tamman eos, 2 for JWL
        ! parameters for ideal gas eos
        real gamma, pw          ! parameters for tamman eos
        ! parameters for explosive
        real A, B, C, R1, R2, w ! parameters for JWL
        real e0,dens0          ! initial data for jwl
        logical iexpo       !
        integer deton_type      !   1 for point-detonation, 2 for plane_detonation
        real deton_pos(3), Vexp(3), Pstatic     !      detonation position, wave velocity, pressure for unexploded part
        real deton_dir(3)           !   direction for plane_detonation
        real ini_time                   !   define intial stage time to control minor dt

        ! thermal parameters
        real Cv, T0, p0         ! specific heat, ref temperature, ref pressure
    contains
    procedure:: init=>eos_init
    procedure:: press=>eos_pressure
    procedure:: sound_spd=>eos_sound_speed
    procedure:: sound_spd2=> eos_sound_speed_2  ! square of sound speed, used to check
                                                ! parabolicity
    procedure:: internal_energy => eos_internal_energy
    ! procedure:: entro2dens
    procedure:: press_adb => eos_pressure_adb       ! calculate adiabatic pressure
    procedure:: entropy
    procedure:: dens_from_entropy
    end type

    type:: eos_mix
        type(eos), pointer:: eos1 => null(), eos2 => null()
        real R1, R2
        contains
        procedure:: init=>eos_mix_init
        procedure:: press=>eos_mix_pressure
        procedure:: sound_spd=>eos_mix_sound_speed
    end type

    
    type(eos)::eos0(10)
    type(eos_mix):: eosm0(10)


    ! interface definition
    interface
        module subroutine eos_init(this,eos_type,data_in,iexpo0, thermal_in)
            implicit none
            class(eos):: this
            ! input
            integer eos_type
            real data_in(:)
            logical,optional::iexpo0
            real, optional:: thermal_in(:)
        end subroutine
    end interface

    interface
        module function eos_pressure(this,dens,En,Fb) result(p)
            implicit none
            class(eos):: this
            ! input
            real dens, en
            real,optional:: Fb
            ! output
            real p
        end function
    end interface

    interface
        module function entropy(this,p,dens)
            implicit none
            class(eos):: this
            ! input
            real p,dens
            ! output
            real entropy
        end function
        module function dens_from_entropy(this,entro, p)
            implicit none
            class(eos):: this
            ! input
            real entro, p
            ! output
            real dens_from_entropy
        end function
    end interface

    interface
        module function eos_sound_speed(this,dens,p) result(c)
            implicit none
            class(eos):: this
            ! input
            real dens,p
            ! output
            real c
        end function
    end interface

    interface
        module function eos_sound_speed_2(this,dens,p) result(c2)
            implicit none
            class(eos):: this
            ! input
            real dens,p
            ! output
            real c2
        end function
    end interface

    interface
    module function eos_internal_energy(this,pressure,density)
        implicit none
        class(eos):: this
        ! input
        real pressure, density
        ! output
        real eos_internal_energy
    end function
    end interface
    interface
    module function eos_pressure_adb(this,P,dens_temp,dens_ifluid)! dens_temp,dens_fluid
        implicit none
        class(eos):: this
        ! input
        real p,dens_temp,dens_ifluid
        ! output
        real eos_pressure_adb
    end function
    end interface

    interface
        module subroutine eos_mix_init(this,eos1,eos2)
            implicit none
            class(eos_mix):: this
            ! input
            class(eos), target:: eos1, eos2
        end subroutine
    end interface
    interface
        module function eos_mix_pressure(this,g1, g2, en) result(p)
            implicit none
            class(eos_mix):: this
            ! input
            real g1, g2, en
            ! output
            real p
        end function
    end interface
    interface
        module function eos_mix_sound_speed(this,g1, g2, en) result(c)
            implicit none
            class(eos_mix):: this
            ! input
            real g1, g2, en
            ! output
            real c
        end function
    end interface

    contains
    
    real function clamp(L,R,x,increase)
        real L,R,X
        integer increase
        if(x<L)then
            clamp=0
        elseif(x>R)then
            clamp=1
        else
            clamp=(x-L)/(R-L)
        endif
        if (increase==0)then
            clamp=1-clamp
        endif
        return
    end function
    end module