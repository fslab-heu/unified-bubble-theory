### Unified bubble theory
#### 1. Purpose
This is the fortran code for the unified bubble theory which is used to simulate
bubble dynamics in different circumstances. 

#### 2. Usage
The code is written with modern fortran. Thus, a fortran compiler is required.
On windows, the visual studio is used as the IDE. 
On linux, it can be compiled with make ultility.

#### 3. Examples
Here we present an example for the 'constants.dat':
```
1.4							! gamma
1500						! sound speed
1e-3						! dynamic viscosity
0.075						! surface tension
2338						! cavitition limit
9.8							! gravity 
1000						! density of surrounding fluid
0.75						! drag coefficient
1.0							! added-mass coefficient
1.0e5						! ambient pressure at z = 0
0,0,0						! current velocity
```


Here we present an example for the 'case.dat':

```
0.0055									! solution time
0.010000								! dt scale
0, 0.7, -1.5							! pressure investigation
.false.									! whether to calculate boundary
0.000000, 0.000000, 0.000000			! boundary center location
0.000000, 0.000000, -1.000000			! boundary normal direction
1.000000								! alpha for boundary
.true.									! migration
2										! number of bubbles
0.0, -93e-3,	0,	-0.1,	0.00263,	120,	1.2e6	! info for bubble 1: delaytime,cx,cy,cz,r0,dr0,p0
0.0,	  0,    0,  -0.1,   0.00242,    120,    1.2e6	! info for bubble 2: delaytime,cx,cy,cz,r0,dr0,p0
```

#### 3. References
Users can refer to the following paper for detailed theory and equations:

```
A.M. Zhang, S.M. Li, S. Li, P. Cui, Y.L. Liu.A unified theory for bubble dynamics. arXiv preprint arXiv:2301.13698, 2023
```