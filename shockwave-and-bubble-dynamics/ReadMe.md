### 1. Purpose
This is the program to simulate the underwater explosion bubble based on the unified bubble theory. There are 3 solvers used in this code:
a. The near-field shock wave formation and the early stage of the bubble expansion with high Mach number is solved with the discontinuous galerkin method for the nonlinear Euler equation in the moving coordinate system.
b. The far-field shock wave is calculated with the nonlinear transportation equation with the initial condition given by the near-field shock wave solver.
c. The later stage of the bubble motion is solved with the unified bubble theory. The progame can give the motion of the bubble, and the pressure history in the fluid field.

### 2. Usage
The excutable file locates in the folder named `binary`. Besides, some dependencies are also included.You may run `undex_zhang_et_al.exe` to start the simulation.
The input parameters are read from 3 input files, i.e. 'bubble.in', 'jwl.in' and 'water.in'.
Please see the example input files for instruction.

Examples:

`case.in` defines the basic case parameters which should be altered almost in each new case:
```
1			! charge weight in kg
5			! charge depth in meters
0,20,0		! coordinate (x,y,z) of pressure point with respect to charge center
0.3			! total solution time
1		! dt scale, default value: 1.0, which is used to scale the time step of the bubble solver.
```
`jwl.in` defines the material constants for the JWL EoS of explosion products. These parameters are used in the near-field shock wave solver.
Here in an example for the TNT explosive:
```
371.2e9			! A
3.231e9			! B
0  1.28e9		! C, if input 0, then calculate with E0
4.15			! R1
0.95			! R2
0.3				! w
1630			! rho0
7e9				! E0*rho0, specific energy per unit volume
```
`water.in` defines the material constants of the Tammann EoS for the surrounding water.
The sound speed calculated from these parameters are also used in the bubble simulation.
Here is an example for the water often used in papers, and the corresponding sound speed is 1427m/s at the atmosphere:
```
7.15		! gamma
3.309e8		! pw
1000		! rho0
```
`bubble.in` defines the constants for the bubble simulation.
```
1e-3	! dynamic viscosity of water
0.075	! surface tension of water
2338	! cavitation limit of water
9.8		! gravity
1.0  	! drag coefficient of the bubble
0.5		! added-mass coefficient of the bubble
-1.0	! reflection coefficient of boudnary, -1 for ideal free surface and 1 for ideal solid wall
```
Note that the boundary included in this code is only the free surface which is often used in UNDEX simulation.
The reflection coefficient should be ranging from -1 to 0. Typically, -1 is a good choice.

An optional argument can be passed in when runing the excutable in the command line terminal. 
The argument should be the working directory enclosed by double quote and ended by `/`. For example:
```
undex_ubt.exe "C:/user/Administrator/fsundex/"
```
will lead the code to read ini files from the given directory and put the results there too.
### 3. Output
You can find simulation results in the `output` folder. The `bubble.dat` holds the radius and migration histories of the bubble.
The `pressure.dat` holds the pressure history at the measuring point.
Note that even the free surface is included, the pressure at the measuring point doesnt consider the reflected wave from the free surface.
The reason is that the reflected wave is usually considered in the far-field shock analysis by adding another source.
The mixed pressure will interfere with the solution procedures.


