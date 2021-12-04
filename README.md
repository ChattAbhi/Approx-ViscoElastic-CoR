# Approx-ViscoElastic-CoR
This MATLAB toolbox contains implementations of several methods for computing the Coefficient of Restitution for a spherical bead impact based on a general nonlinear viscoelastic contact model, while including the effects of any external load on the bead. The toolbox provides implementations of the methods discussed in [1,2] for both numerical and approximate analytical Coefficient of Restitution, according to the specified viscoelastic model. The set of parameters that define the dynamics of the bead during impact, according to the viscoelastic model, include:
 1) The mass of the bead (must be > 0)
 2) Stiffness constant for the viscoelastic model (must be >=0)
 3) Damping Ratio for the viscoelastic model (must be >= 0)
 4) Exponent on the nonlinear stiffness term (must be >=1)
 5) Exponent on the nonlinear damping term  (must be >=1)
 6) Pre-impact velocity (must be > 0, the sign convention used in the method considers the velocity positive towards the direction of the impact)
 7) Gravitational acceleration constant (must be >=0, nominal = 9.8 m/s (if using SI units))
 8) Magnitude of the compressive external force (must be >= 0)
See the references [1,2] for more details on the method of computation.

## Installation
Simply download or cloan this repository and set path to the folder on MATLAB, using the `addpath` command.

## Instructions on usage
To see the description of any MATLAB function file provieded as part of this toolbox, type `>> help function_name` in your MATLAB command window. `function_name` being the name of any file provided in this toolbox. User can either directly work with the unscaled physical (dimensional) parameters or scaled (dimentionless) parameters. 

### Working with unscaled physical quantitues 
To work dirctly with the physical (dimensional) quantities, one simply needs to call the function `getCOR()`. The `getCOR()` function can be used to compute the Coefficient of Restitution using any of the methods described in [1,2]. The general syntax for this function is as follows:
```
getCOR(m,k,gam0,v0,alp,bet,g,F,method)
```
The definitions of the arguments are listed below:
 - `m` is the value of the mass in appropriate units (must be > 0)
 - `k` is the value of the stiffness constant in appropriate units (must be >=0)
 - `gam0` is the value of the damping ratio in appropriate units (must be >=0)
 - `v0` is the value of pre-impact velocity in appropriate units (must be >0)
 - `alp` is the value of the exponent on the nonlinear stiffness term (must be >=1)
 - `bet` is the value of the exponent on the nonlinear damping term (must be >=1)
 - `g` is the value for gravitational acceleration constant in appropriate units (must be >=0)
 - `F` is the value for the magnitude of the compressive force in appropriate units (must be >=0)
 - `method` is a character array or a string input to set the method used for the computation. A list of valid inputs for the `method` arguement and their respective descriptions are given below:
    - `method = 'ord2approx'` (Default) computes the approximate value for the Coefficient of Restitution based on a second-order expansion, while approximating the coefficient integral functions. The accuracy of this approximation method is higher when the external force on the system is either very small or very large. 
    - `method = 'ord1approx'` computes the approximate value for the Coefficient of Restitution using the first-order expansion. This approximation method is very accurate when the external load on the system is small, and becomes less accurate as the magnitude of the external load increases.  
    - `method = 'ord2schw-pos'` computes the approximate value for the Coefficient of Restitution when `bet=3/2`, using a second-order expansion that is based on the coefficient values reported by Schwager and Poschel in [3]. This approximation method is very accurate when the external load on the system is small, and becomes less accurate as the magnitude of the external load increases.  
    - `method = 'ord2numint'` computes the approximate value for the Coefficient of Restitution using the second-order exapansion, while numerically computing the integral coefficient functions. The accuracy of this approximation method is higher when the external force on the system is either very small or very large. 
    - `method = 'dirnumint'` computes the value for the Coefficient of Restitution by numerically integrating the viscoelastic contact law-based equation of motion. This method is accurate up to the specified absolute and relative error tolerences for the numerical integration.
Additional settings for the `getCOR()` function can be specified if a numerical integration based method is selected (i.e. `method='ord2numint'` or `method=dirnumint`). These settings are passed in as Name-Value pairs to the function call,
```
getCOR(m,k,gam0,v0,alp,bet,g,F,method,Name,Value,_)
```
where `Name` is a string or character array input representing the setting, and `Value` is the numerical value corresponding to the setting. The list of acceptable `Name` keywords are listed below:
- `Name = 'AbsTol'` used to set the absolute error tolerence for the numerical integration. The default value is `1E-8` for `method = ord2numint` and `1E-14` for `'method = dirnumint'`.
- `Name = 'RelTol'` used to set the relative error tolerence for the numerical integration. The default value is `1E-7` for `method = 'ord2numint'` and `1E-10` for `method = dirnumint`. 
- `Name = 'MaxIter'` used to set the maximum number of allowed iterations for the numerical integrator. The default value is `100`.
- `Name = 'EquTol'` used to set a tolerence for detecting the equilibrium point during the numerical integration. This setting is only valid for `method = 'dirnumint'`. The default value is `1E-8`.

An example of a bouncing steel ball simulation is provided in the file `ExampleBouncingBall\ExampleBouncingBall.m`. Running this file generates plots for the height, velocity, and the Coefficient of Restitution values of the ball, as shown below,
<p align="center">
  <img src="https://github.com/ChattAbhi/Approx-ViscoElastic-CoR/blob/main/ExampleBouncingBall/bounce_height.png" width="450">
  <img src="https://github.com/ChattAbhi/Approx-ViscoElastic-CoR/blob/main/ExampleBouncingBall/bounce_velocity.png" width="450">
  <img src="https://github.com/ChattAbhi/Approx-ViscoElastic-CoR/blob/main/ExampleBouncingBall/bounce_CoR.png" width="450">
</p>

For additional help with the `getCOR()` function, type `>> help getCOR` in the MATLAB command window.

### Working with scaled dimensionless parameters 
To work with the scaled (dimensionless) parameters, one first needs to convert the physical quantities to scaled parameter values. This can be done using the command:
```
[gam,gtil]=ScaledParams(m,k,gam0,v0,alp,bet,g,F)
```
This will compute the values of the two scaled parameters required for the computations in the variables `gam` and `gtil`. The definitions of the inputs `m`, `k`, `gam0`, `v0`, `alp`, `bet`, `g`, and `F` are same as described above. For more details on the scaling, see [1,2].

Once the scaled paramters `gam` and `gtil` is known, the function `analyticCOR()` function can be used to calculate the Coefficient of Restitution using the second-order expansion, as described in [1,2]. This function can be called as,
```
analyticCOR(alp,bet,gam,gtil,choice)
```
The inputs `alp`, `bet`, `gam`, and `gtil` are known, and the input `choice` is a string or character array used to set the computation method. The two valid inputs for `choice` are `'approx'` (Default) and `'numint'`. When `choice = 'approx'` the function uses approximate coefficient integral values to evaluate the second-order expansion, and when `choice = 'numint'` the coefficient integrals are numerically evaluated. In case of `choice = 'numint'`, additional arguements related to numerical integration settings can be passed in as Name-Value pairs. The valid Name keywords are `AbsTol`, `RelTol`, and `MaxIter`, definitions are identical to what is mentioned in the previous section.

To compute the Coefficient of Restitution using the first-order exapansion, as described in [1,2], using the scaled parameters `gam` and `gtil`, one needs to first compute the constant coefficients for the first order expansion. This can be done using the command:
```
[C0,C1] = ConstantCORCoeffs(alp,bet) 
```
the outputs `C0` and `C1` are the coefficients needed for the first-order approximation. So the value of the Coefficient of Restitution, stored is some variable `e` can be computed as,
```
e = 1 - gam*C0 - gam*gtil*C1
```
Note, that with this computation the Coefficient of Restitution will not be bounded in the set [0,1], to keep this value bounded one can alternatively compute: 
```
e = max([1 - gam*C0 - gam*gtil*C1,0]);
```

When `bet=3/2`, the function `ConstantCORCoeffs` can also be used to compute the second-order constant coefficient using the value reported by Schwager and Poschel in [3]. This can be done by calling the function as,
```
[C0,C1,C2] = ConstantCORCoeffs(alp,3/2) 
```
the output `C2` is the second-order coefficient. This can be used to compute the Coefficient of Restitution value `e` as,
```
e = 1 - gam*C0 - gam*gtil*C1 - C2*gam^2
```
Again to keep the value of `e` bounded within the set [0,1], one can alternatively compute:
```
e = max([1 - gam*C0 - gam*gtil*C1 - C2*gam^2,0]);
```

Lastly, to compute the value of the Coefficient of Restitution via a direct numerical integration approach, one can use the function `numericCOR()`. This function can be called as,
```
numericCOR(alp,bet,gam,gtil)
```
to output the value of the numerically integrated value of the Coefficient of Restitution. Additional settings can be passed in to the function call as Name-Value pair arguments. The absolute error, relative error, and equilibrium tolerences can be set using the Name keywords `AbsTol`, `RelTol`, and `EquTol`, respectively. To view plots the numerically integrated state-trajectories (scaled displacement and velocities), one can use the Name keyword `'ShowPlots'` and set the corresponding boolean value to `true`. These plots can be automatically saved as .fig or .pdf files using the Name keywords `'figsave'` and `'pdfsave'` respectively. One can also obtain the state-trajectory data using the name keyword `'OutDat'`, using the following syntax,
```
[e,T,U,dU]=numericCOR(alp,bet,gam,gtil,'OutDat',true)
```
which yields the outputs `e` as the value of the Coefficient of Restitution, `T` as the time date, `U` as the scaled-displacement data, and `dU` as the scaled-velocity data. See `>> help numericCOR` for more details. 

While working with the scaled parmeters `gam` and `gtil`, the references [1,2] also discuss about a critical value `gam_c`, such that for any value of `gam >= gam_c`, the Coefficient of Restitution value must be 0. This toolbox also provides a function `gamma_c` to calculate this value of `gam_c`, given the value of `gtil`. This function can be used as follows:
```
gam_c = gamma_c(alp,bet,gtil)
```
Hence, `gam_c` may be used to determine whether or not the spherical object will rebound after the impact. See references [1,2] for details.

The simulation results in references [1,2] were generated using the files in the folder `PaperResults`. The file `PaperResults/GenFigs.m` can be used to generate plots of all simulation results in [1,2], which are then stored in various sub-folder within `PaperResults/Figures`. If the user finds this toolbox useful, the authors would appreciate a citation to [1,2].

## References 
* [1] Chatterjee, A., James, G., Brogliato, B., "Approximate coefficient of restitution for nonlinear viscoelastic contact with external load" (Under Review). Granular Matter. Preprint: https://hal.inria.fr/hal-03463883
* [2] Chatterjee, A., James, G., Brogliato, B., "Approximate analytical coefficient of restitution formulation for single bead impact with external load, using nonlin- ear visco-elastic models". Research Report hal-03462750, INRIA (December 2021). https://hal.inria.fr/hal-03462750
* [3] Schwager, T., Po ̈schel, T., "Coefficient of resti- tution for viscoelastic spheres: The effect of delayed recovery". Physical Review E 78(5), 051304 (2008)
