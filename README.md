# Approx-ViscoElastic-CoR
This MATLAB toolbox contains implementations of several methods for computing the Coefficient of Restitution for a spherical bead impact based on a general nonlinear viscoelastic contact model, while including the effects of any external load on the bead. The toolbox provides implementations of the methods discussed in [1,2] for both numerical and approximate analytical Coefficient of Restitution, according to the specified viscoelastic model. The set of parameters that define the dynamics of the bead during impact, according to the viscoelastic model, include:
 1) The mass of the bead (must be > 0)
 2) Stiffness constant for the viscoelastic model (must be >=0)
 3) Damping Ratio for the viscoelastic model (must be >= 0)
 4) Exponent on the nonlinear stiffness term (must be >=1)
 5) Exponent on the nonlinear damping term  (must be >=1)
 6) Pre-impact velocity (must be > 0, the sign convention used in the method considers the velocity positive towards the direction of the impact)
 7) Gravitational acceleration constant (must be >=0, nominal = 9.8)
 8) Magnitude of the compressive external force (must be >= 0)
See the references [1,2] for more details on the method of computation.

## Installation
Simply download or cloan this repository and set the path to the folder on MATLAB.

## Instructions on usage
This toolbox allows user to either directly work with the physical (dimensional) parameters or scaled (non-dimensional) parameters. To work dirctly with the physical (dimensional) quantities, one simply needs to call the function `getCOR()`. The `getCOR()` function can be used to compute the Coefficient of Restitution using any of the methods described in [1,2]. The general synatax for this function is as follows:
```
getCOR(m,k,gam0,v0,alp,bet,g,F,method)
```
Here the arguments 
 - `m` accepts the value of the mass (must be > 0)
 - `k` accepts the value of the stiffness constant (must be >=0)
 - `gam0` accepts the value of the damping ratio (must be >=0)
 - `v0` accepts the value of pre-impact velocity (must be >0)
 - `alp` accepts the value of the exponent on the nonlinear stiffness term (must be >=1)
 - `bet` accepts the value of the exponent on the nonlinear damping term (must be >=1)
 - `g` accepts the value for gravitational acceleration constant (must be >=0)
 - `F` accepts the value for the magnitude of the compressive force (must be >=0)
 - `method` accepts a character array to set the method used for the computation. A list of valid inputs for the `method` arguement and their respective descriptions are given below:
    - `method = 'ord2approx'` computes the Coefficient of Restitution using a second order approximation, while approximating

