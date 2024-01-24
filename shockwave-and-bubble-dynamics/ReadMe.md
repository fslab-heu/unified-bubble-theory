### 1. Purpose
This is the program to simulate the underwater explosion bubble based on the unified bubble theory. 
The early stage of the bubble expansion with high Mach number is solved with the discontinuous galerkin method for the nonlinear Euler equation.
The later stage of the bubble motion is solved with the unified bubble theory. The progame can give the motion of the bubble, and the pressure history in the fluid field.

### 2. Usage
The input parameters are read from 3 input files, i.e. 'bubble.in', 'jwl.in' and 'water.in'.
Please see the example input files for instruction.

### 3. Output
You can find simulation results in the `output` folder. The `bubble.dat` holds the radius and migration histories of the bubble.
The `pressure.dat` holds the pressure history at the measuring point.