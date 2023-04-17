platform := win
platform := linux

fcc := ifort
FLAGS := -O3  -i4 -fpp  

# basic libs
sources := linear_equation.f90 \
	math.f90 \
	polynomial.f90 \
	fsintegration.f90

# bubble solver
sources += global.f90 \
    boundary.f90 \
	bubble_module.f90 \
	bubble_init.f90

sources += shock_wave_module.f90 \
	far_field_output.f90 \
    restart.f90 \
    shockwave_DG_output.f90 \
    shockwave_output.f90 \
    update_dt.f90 \
	main.f90
    
    
	
			

ifeq ("$(platform)","win")
	suffix := obj
	rm := del
    FLAGS += /real-size:64  /NODEFAULTLIB:MSVCRT.lib
else
	suffix := o
	rm := rm
	FLAGS += -r8
endif



%.$(suffix):%.f90
	$(fcc) $(FLAGS) -c $< 


objects := $(sources:.f90=.$(suffix))


all:$(objects)
	$(fcc) -o main.exe $(FLAGS)  $(objects) $(LIBS)

.PHONY: clean
clean:
	$(rm) *.exe
	$(rm) *.mod
	$(rm) *.$(suffix)
	$(rm) *.exp
	$(rm) *.a