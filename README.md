# Toro

This package implements the exact Riemann solver for the Euler equations as described in

Toro, E. F. "Riemann Solvers and Numerical Methods for Fluid Dynamics" 2nd ed. Ch 4, pp 115-162. Springer. 1999 ISBN: 3-540-65966-8.

# Usage

There are currently three main functions to use in this package:

The first is:
```julia
CenterPressure( gamma_L, P_L, rho_L, u_L, gamma_R, P_R, rho_R, u_R )
```
This returns an estimate of the pressure in the "center region" (also called the "star region") that satisfies Eq. 4.5. The iterative procedure required to find this value is currently performed with a bisection method.
Here, 
* gamma is the ratio of specific heats
* P is the pressure
* rho is the density
* u is the velocity
* The subscripts \_L or \_R indicate whether the variable refers to properties to the (L)eft or (R)ight of the interface.

The velocity in the "center region"/"star region" is found using
```julia
CenterVelocity( P_0, gamma_L, P_L, rho_L, u_L, gamma_R, P_R, rho_R, u_R )
```
The arguments here have the same definitions as `CenterPressure`.
The new argument `P_0` refers to the pressure in the center/star region, as found using `CenterPressure`.

The density on either edge of the volume after the interaction of any shock/rarefaction from the interface is found using
```julia
EdgeDensity( P_0, gamma, P, rho, u )
```
where these variables have the same definitions as above, but the subscripts are omitted as they refer to the initial state (before the interaction of any shock/rarefaction wave) of the gas the wave is moving in to.
In other words, the initial state of the left or right gas.

As always, information on the functions in this package, including equations and implementation details, can be found using
```julia
?> Toro
?> CenterPressure
?> CenterVelocity
?> EdgeDensity
```

# Notes

While care has been taken to ensure the results are accurate, bugs can and do creep in and so you should be sure the results you obtain make sense.
If you find anything not behaving as expected, please open an issue to report it or a pull request with a fix.
Tests on the validity of the solutions found from this package are performed against the tabulated values in Table 4.3 of the book.
If there are additional tests that would be good to have, feel free to open a pull request.

# Roadmap

There are a few more features I would like to add to this package:

* Implement a Gas (or similar) type that contains information on density, pressure, velocity, and ratio of specific heats. This will allow overloaded function definitions that take in far fewer arguments than the current implementation requires, thus making the package easier to use.
* Improve the iteration method in CenterPressure to use Newton-Raphson or similar for faster convergence.
