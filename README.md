# Quasi-1D Nozzle Flow Solver

A numerical solver for transient compressible flow through converging-diverging nozzles using MacCormack's predictor-corrector method.

The solver integrates the quasi-one-dimensional Euler equations in conservative form and marches the solution in time until a steady-state flow field is obtained. Geometric source terms account for variations in nozzle area, while artificial viscosity enhances numerical stability in transonic regions.

---

## Features

* Quasi-one-dimensional compressible flow formulation
* Conservative Euler equations
* MacCormack predictor-corrector time integration
* CFL-based adaptive timestep
* Pressure-based artificial viscosity
* Subsonic and supersonic outlet boundary conditions
* Custom nozzle area distributions
* Computation of:

  * Pressure
  * Temperature
  * Density
  * Velocity
  * Mach number
  * Mass flow rate

---

## Governing Equations

The solver advances the conservative variables

The conservative variables are:

U1 = ρA

U2 = ρuA

U3 = ρA[T/(γ−1) + (γ/2)u²]

using the quasi-one-dimensional Euler equations with geometric source terms arising from nozzle area variation.

---

## Numerical Method

The solution is obtained using MacCormack's two-step predictor-corrector scheme:

1. Forward-difference predictor step
2. Backward-difference corrector step
3. Averaging of predicted and corrected derivatives

The timestep is computed dynamically from a CFL stability condition.

Artificial viscosity is applied to suppress non-physical oscillations and improve convergence in shock-capturing cases.

---

## Default Geometry

The default converging-diverging nozzle is

[
A(x)=1+2.2(x-1.5)^2
]

normalized by the throat area.

---

## Repository Structure

```text
main.py                    # Solver driver

mesh.py                    # Mesh generation and nozzle geometry
initialization.py          # Initial flow field
physics.py                 # Variable conversions and flux evaluation
solver.py                  # MacCormack solver implementation
boundary_conditions.py     # Boundary treatments
postprocess.py             # Visualization
```

---

## Example Output

The solver produces steady-state distributions of:

* Pressure
* Temperature
* Mach number
* Mass flow rate

along the nozzle centerline.

---

## Current Limitations

* Designed primarily as an educational and research-oriented CFD implementation
* Sensitivity to initialization strategy
* Sensitivity to nozzle area distributions

---

## Future Work

- Improve robustness for arbitrary nozzle geometries and area distributions.
- Develop a generalized initialization strategy that converges reliably across a wider range of flow regimes.
- Investigate the sensitivity of convergence and stability to boundary condition treatment.
- Improve code modularity and support for user-defined nozzle geometries.

  ---
  ---
