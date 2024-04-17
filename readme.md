### Unified bubble theory
#### 1. Purpose
This is the fortran code for the unified bubble theory which is used to simulate bubble dynamics in different circumstances. 

#### 2. Usage
The code is written with modern Fortran. Thus, a Fortran compiler is required. We recommend the intel oneapi HPC toolkit.
On windows, the visual studio is used as the IDE. 
On linux, it can be compiled with make utility.

#### 3. Contents
There are two folders here:
1. unified-bubble-theory: the bubble dynamics model to simulate the bubble motion in different circumstances.
2. unified-bubble-theory-for-underwater-explosion: the complete code for the underwater explosion simulation, including the nearfield solver, far-field solver and bubble solver.

The instructions for both can be found in the `readme.md` file in each folder, respectively.

#### 3. References
Users can refer to the following paper for detailed theory and equations:

[A-Man Zhang,  Shi-Min Li,  Pu Cui,  Shuai Li, and  Yun-Long Liu, A unified theory for bubble dynamics. Physics of Fluids, 2023. 35(3): 033323](https://aip.scitation.org/doi/10.1063/5.0145415)
