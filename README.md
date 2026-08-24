# Quasi-1D Nozzle Flow Solver

A numerical solver for transient compressible flow through
converging-diverging nozzles using the MacCormack predictor-corrector
method.

The solver integrates the quasi-one-dimensional Euler equations in
conservative form and marches the solution in time until a steady-state
flow field is obtained. Geometric source terms account for variations in
nozzle area, while artificial viscosity enhances numerical stability in
transonic and shock-containing cases.

------------------------------------------------------------------------

## Features

-   Quasi-one-dimensional compressible flow formulation
-   Conservative Euler equations
-   MacCormack predictor-corrector time integration
-   CFL-based adaptive timestep
-   Pressure-based artificial viscosity
-   Subsonic and supersonic outlet boundary conditions
-   Custom nozzle area distributions
-   Computation of:
    -   Pressure
    -   Temperature
    -   Density
    -   Velocity
    -   Mach number
    -   Mass flow rate

------------------------------------------------------------------------

## Governing Equations

The solver advances the following conservative variables:

$$U_1 = \rho A$$
$$U_2 = \rho u A$$
$$U_3 = \rho A \left(\frac{T}{\gamma - 1}+\frac{\gamma}{2}u^2\right)$$

The quasi-one-dimensional Euler equations are solved with geometric
source terms arising from nozzle area variation.

------------------------------------------------------------------------

## Numerical Method

The solution is obtained using MacCormack's two-step predictor-corrector
scheme:

1.  Forward-difference predictor step
2.  Backward-difference corrector step
3.  Averaging of the predicted and corrected derivatives

The timestep is computed dynamically from a CFL stability condition.

Artificial viscosity is applied to suppress non-physical oscillations
and improve convergence in shock-capturing cases.

------------------------------------------------------------------------

## Default Geometry

The default converging-diverging nozzle is defined by

$$\[ A(x)=1+2.2(x-1.5)\^2 \]$$

with the area normalized by the throat area.

------------------------------------------------------------------------

## Repository Structure

``` text
main.py                    # Solver driver

mesh.py                    # Mesh generation and nozzle geometry
initialization.py          # Initial flow field
physics.py                 # Variable conversions and flux evaluation
solver.py                  # MacCormack solver implementation
boundary_conditions.py     # Boundary treatments
postprocess.py             # Visualization
```

------------------------------------------------------------------------

## Convergence History

The solver was marched for approximately 2500 iterations. The residuals
for continuity, momentum, and energy decrease by several orders of
magnitude, reaching the prescribed convergence tolerance of $$(10\^{-8})$$.

The oscillatory behavior during convergence is characteristic of the
transient approach to the steady solution, while the overall downward
trend indicates stable convergence.

<img width="1669" height="831" alt="Image" src="https://github.com/user-attachments/assets/70707cff-bf36-41e3-8d30-511eadb6da7f" />

------------------------------------------------------------------------

## Numerical Validation

The computed Mach-number distribution was compared against the
analytical isentropic solution for the same nozzle geometry.

The numerical and analytical solutions show close agreement over the
nozzle domain, demonstrating that the solver correctly captures the
expected subsonic-to-supersonic acceleration through the
converging-diverging nozzle.

<img width="1669" height="831" alt="Image" src="https://github.com/user-attachments/assets/63db4fda-c52d-45aa-acd1-a3fe9af40fd9" />

### Absolute Error

The absolute difference between the numerical and analytical Mach
numbers remains small across the domain. The largest deviations occur
near the throat and in regions where the numerical solution transitions
between flow regimes.

### Relative Error

The relative error decreases substantially through the nozzle, reaching
values close to zero over much of the domain. The larger relative error
near the inlet is expected because the Mach number itself is small
there, making the relative error more sensitive to small absolute
differences.

------------------------------------------------------------------------

## Example Output

The solver produces steady-state distributions of:

-   Pressure
-   Temperature
-   Mach number
-   Mass flow rate

along the nozzle centerline.

The convergence history and numerical-validation plots above provide
additional evidence of numerical stability and agreement with the
analytical solution.

------------------------------------------------------------------------

## Current Limitations

-   Designed primarily as an educational and research-oriented CFD
    implementation
-   Sensitivity to initialization strategy
-   Sensitivity to nozzle area distributions
-   Numerical accuracy and convergence can depend on the treatment of
    boundary conditions and artificial viscosity

------------------------------------------------------------------------

## Future Work

-   Improve robustness for arbitrary nozzle geometries and area
    distributions.
-   Develop a generalized initialization strategy that converges
    reliably across a wider range of flow regimes.
-   Investigate the sensitivity of convergence and stability to boundary
    condition treatment.
-   Improve code modularity and support for user-defined nozzle
    geometries.
-   Extend validation to cases involving shocks and stronger transonic
    effects.

------------------------------------------------------------------------
