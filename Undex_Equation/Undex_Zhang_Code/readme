### 1. Purpose
This is the program to simulate the underwater explosion bubble based on the unified bubble theory. There are 3 solvers used in this code:



```
# Underwater Explosion Bubble Simulator

## Purpose

This program simulates underwater explosion bubbles using a unified-bubble approach. It contains three coupled solvers:

- Near-field solver: discontinuous Galerkin solution of the nonlinear Euler equations in a moving coordinate system to capture shock formation and early high-Mach bubble expansion.
- Far-field solver: nonlinear transport-type solver that propagates the shock to larger distances using initial data from the near-field solver.
- Bubble dynamics solver: unified bubble model that computes the later-stage bubble motion and pressure history.

## Build

The project provides a `Makefile` in the project root. By default the `Makefile` sets the Fortran compiler in the `FC` variable. To build:

```sh
make        # uses FC from the Makefile (edit FC to ifort/gfortran/ifx as needed)
```

If you wish to use a different compiler, change `FC` and the correspoinding flags as you like.

On windows, you can use Visual Studio with intel oneapi to compile and run the code. 

## Run

Place the required input files (see below) in the project root and run the produced executable (default `main`):

```sh
./main
```

The program prints stage markers and periodic progress lines to standard output. Progress messages and diagnostic frequency are controlled in `src/advance.f90` and the top-level messages in `src/main.f90`.

## Required input files

Put these files in the project root (same folder as the executable):

- `case.inp` — primary run configuration. Typical items included here are simulation time settings (`tend`, `t_arrive`), time step controls, mesh/configuration choices, bubble initial parameters (initial radius `Rb`, initial pressure, etc.) and output frequency. The exact keyword names follow the parsing code in `initialize.f90` (look for reads of `case.inp`).

- `water.mat` — material/EOS data for the surrounding fluid (water). This file contains material parameters used by the equation-of-state routines. Keep it in the project root so the binary can read material properties at startup.

- `jwl.mat` — JWL (Jones–Wilkins–Lee) explosive parameters used when modeling explosive loading. If your case involves an explosive charge, place the corresponding JWL parameter file here.

Notes:
- The code expects these names exactly unless you change the I/O code in `initialize.f90` or other initialization routines.

## Output / result files

The code writes several result files during a run. Typical files and their meanings:

- Unit 101..103 (near-field outputs): time series or snapshot files produced by the near-field solver. Check `src/advance.f90` and `src/near_field.f90` for exact formats and file unit numbers used by `output_near_field`.
- Unit 200, 201 (bubble outputs): bubble-derived time histories (pressure pulses, R(t), center position, etc.). The code uses formatted `write` calls to these units in the bubble dynamics routine.
- Far-field output files: produced by the far-field solver via `output_far_field` (see `src/far_field.f90` for filenames and formats).
- `main` (or configured executable): the compiled program. Standard output contains progress logs printed during the simulation.

If you need specific column formats or a description of each output file, I can extract the exact `write`/`open` statements from the source and document them here.

## References
```
1. Zhang, A.M., et al., A unified theory for bubble dynamics. Physics of Fluids, 2023. 35(3): p. 033323.
2. Zhang, A.M., et al., Theoretical study on bubble dynamics under hybrid-boundary and multi-bubble conditions using the unified equation. Science China Physics, Mechanics & Astronomy, 2023. 66(12).
3. Liu, Y.L., Zhang, A.M., et al., A novel model for full-range underwater explosion shock wave formation and propagation. Under review.
```

## Contact / further help

Contace the author for further asistance:
yunlong_liu@hrbeu.edu.cn


