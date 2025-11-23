# Fortran compiler
FC = ifort

EXE = main
FLAGS = -g -I./lib -L./lib -lfsmath -lfseos
FLAGS += -r8
FLAGS += -fpe0 -traceback
SRC = model.f90 \
	DG.f90 \
	far_field.f90 \
	near_field.f90 \
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
	rm -f $(OBJ) $(EXE)
