#platform := win
platform := linux

fcc := ifort
FLAGS := -O3  -r8 -i4 -fpp    
sources := global.f90 \
    boundary.f90 \
	bubble_module.f90 \
	initiation.f90 \
	main.f90 
	
			

ifeq ("$(platform)","win")
	suffix := obj
	rm := del
else
	suffix := o
	rm := rm
endif



%.$(suffix):%.f90
	$(fcc) -c $<


objects := $(sources:.f90=.$(suffix))

all:$(objects)
	$(fcc) -o main.exe $(objects)

.PHONY: clean
clean:
	$(rm) *.exe
	$(rm) *.mod
	$(rm) *.$(suffix)