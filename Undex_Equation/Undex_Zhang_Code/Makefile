# Fortran compiler
FC = ifort

EXE = main
FLAGS = -g
FLAGS += -r8
FLAGS += -fpe0 -traceback

# basic fsmath source files
SRC = linear_equations.f90 \
	math.f90 \
	polynomial.f90 \
	integration.f90

# fluid eos source files
SRC += fluid_eos.f90 \
	func_pressure_sub.f90 \
	initiation_sub.f90 \
	func_sound_speed_sub.f90 \
	func_internal_energy_sub.f90 \
	func_mixture_sub.f90 \
	func_entropy.f90

# undex_ube source files
SRC += model.f90 \
	DG.f90 \
	near_field.f90 \
	far_field.f90 \
	bubble.f90 \
	advance.f90 \
	initialize.f90 \
	DG_sub.f90 \
	main.f90

OBJ = $(SRC:.f90=.o)
all: $(EXE)

$(EXE): $(OBJ)
	$(FC) -o $@ $^ $(FLAGS)

%.o: src/%.f90
	$(FC) -c $< $(FLAGS)

clean:
	rm -f $(OBJ) $(EXE) *.mod *.smod src/*.mod src/*.smod src/*.o
