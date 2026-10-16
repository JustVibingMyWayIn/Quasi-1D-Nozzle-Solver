# Quasi-1D Nozzle Flow Solver

**A finite-difference analysis of compressible flow through converging–diverging nozzles, using the quasi-one-dimensional Euler equations and the MacCormack predictor–corrector scheme.**

---

## Objective

The objective of this project is to perform a **finite-difference method (FDM) analysis of the quasi-one-dimensional (quasi-1D) nozzle flow model** for compressible, inviscid, adiabatic flow through a converging–diverging duct of arbitrary cross-sectional area distribution `A(x)`.

Concretely, the solver:

- Integrates the time-dependent, conservative-form quasi-1D Euler equations to a steady state using MacCormack's explicit predictor–corrector scheme.
- Supports arbitrary smooth nozzle area distributions, not just a single fixed shape.
- Initializes and validates the numerical steady state against a closed-form analytical solution (isentropic, and isentropic-with-a-normal-shock) of the same quasi-1D equations.
- Reports convergence behavior and discretization error so that the fidelity of the FDM solution can be judged quantitatively.

This is a from-scratch numerical methods project, not a wrapper around an existing CFD package — the flux evaluation, time-stepping, artificial viscosity, boundary conditions, and the analytical reference solution are all implemented directly.

---

## Theory and Governing Equations

The flow is modeled as inviscid, adiabatic, and calorically perfect, with properties uniform across each cross-section (the quasi-1D assumption). The solver advances the conservative variable vector

$$U_1 = \rho A, \qquad U_2 = \rho u A, \qquad U_3 = \rho A\left(\frac{T}{\gamma - 1} + \frac{\gamma}{2}u^2\right)$$

where all quantities are non-dimensionalized by reservoir stagnation conditions ($p_0 = \rho_0 = T_0 = 1$) and the local speed of sound at the reservoir. The quasi-one-dimensional Euler equations are solved in conservative form, with a geometric source term arising from the area gradient `dA/dx` (area change acts like a source/sink of momentum flux, rather than being folded into the convective flux itself).

The local Mach number, pressure, density, and temperature are recovered from the conservative variables through the isentropic/perfect-gas relations implemented in `physics.py`.

---

## Numerical Method

### MacCormack Predictor–Corrector Scheme

The governing equations are integrated in pseudo-time using MacCormack's explicit two-step scheme:

1. **Predictor step** — forward-difference spatial derivatives advance the solution by one time step.
2. **Corrector step** — backward-difference spatial derivatives are applied to the predicted state.
3. The final update is the average of the predicted and corrected time derivatives, giving the scheme second-order accuracy in space while remaining explicit.

### Time Stepping

The time step is recomputed every iteration from a CFL stability condition based on the local flow speed and speed of sound (`physics.compute_time_step`), so the march automatically adapts to the local flow state rather than using a fixed `dt`.

### Artificial Viscosity

A pressure-gradient-switched artificial viscosity term (`solver.artificial_viscosity`) is added to each conservative variable at every step. It is proportional to the local second difference of the pressure field, scaled by a pressure-based sensor so that it activates near sharp gradients (e.g. a captured shock) and stays small in smooth regions. This suppresses non-physical oscillations without requiring a different scheme in shock-containing cases.

### Convergence Monitoring

At every iteration, the RMS change in each conservative variable between successive steps is recorded (`convergence.py`). The run is considered converged once the residual for continuity, momentum, *and* energy all drop below a prescribed tolerance (10⁻⁸ in the default configuration), and the residual history can be plotted on a semi-log axis to inspect the convergence trend.

### Physical-Validity Guard

`physics.check_physical` checks every field after each step for non-positive density/temperature or non-finite values and raises an immediate, descriptive error (rather than silently propagating NaNs) if the solution becomes unphysical — pointing to the likely remedy (lower CFL, higher artificial viscosity coefficient, or a finer grid).

---

## Nozzle Geometry Module

Nozzle shapes are not hard-coded. `geometry.py` defines a `Nozzle` class that wraps an arbitrary, user-suppliable area function `A(x)` so the rest of the solver is geometry-agnostic. A nozzle can be constructed:

- **Analytically**, from a closed-form area function or its exact derivative (`Nozzle.from_function`)
- **From a parabolic area law**, `A(x) = 1 + k(x - L/2)^2` (`Nozzle.parabolic`, the default case)
- **From a smooth cosine-blended converging–diverging contour** with independently specified inlet and exit area ratios (`Nozzle.cosine`)
- **From a radius distribution**, `A = \pi r(x)^2` (`Nozzle.from_radius`)
- **From discrete (x, A) data points**, fitted with a cubic spline (`Nozzle.from_points`)
- **From a CSV file** of measured or design coordinates (`Nozzle.from_csv`)

The class also locates the throat (global area minimum) precisely via bounded scalar minimization rather than just the nearest grid point, and exposes a `validate()` method that:

- Confirms the area distribution is positive and finite everywhere.
- Warns if more than one local minimum exists (multiple throats, which the steady analytical comparison does not support).
- Warns if the minimum lies at a domain boundary (no true converging–diverging throat).
- Estimates the equivalent wall half-angle of the duct and warns if it exceeds a configurable threshold (default 45°), flagging where the quasi-1D assumption itself becomes questionable.
- Checks whether the chosen grid resolution adequately resolves the area variation.

---

## Initial Conditions and Analytical Reference Solution

`Analytical.py` implements the closed-form quasi-1D isentropic flow solution (and the isentropic-flow-with-a-normal-shock solution) for a given area distribution and exit pressure, using the area–Mach relation solved numerically with Brent's method. It determines the flow regime for a given back pressure — subsonic throughout, shock-free supersonic (design condition), a normal shock standing inside the duct, or over-expanded with the shock pushed outside the duct — and, where a shock stands inside the duct, locates it by root-finding on the downstream exit-pressure match.

This analytical solution serves two roles in the project:

1. **Initialization** (`initialization.py`): the FDM march starts from the analytical solution for the given geometry rather than an arbitrary guess, improving robustness and reducing the number of iterations needed to reach steady state.
2. **Validation**: once the MacCormack march converges, the numerical Mach-number distribution is compared directly against this analytical solution to quantify discretization error.

---

## Boundary Conditions

`boundary_conditions.py` implements:

- **Subsonic reservoir inlet**: stagnation pressure and temperature are held fixed at the reservoir values, while the inlet velocity is linearly extrapolated from the interior; static density and temperature then follow isentropically. This remains valid even when the inlet Mach number is not small.
- **Subsonic outlet**: a specified static back pressure is imposed, with density and velocity extrapolated from the interior — producing either a fully subsonic profile or a normal shock inside the duct, depending on the pressure ratio.
- **Supersonic outlet**: all variables are linearly extrapolated from the interior, appropriate when the exit flow is supersonic and no physical boundary condition should be imposed there.

---

## Repository Structure

```
main.py                    # Driver: sets up geometry, mesh, runs the MacCormack march, plots results
geometry.py                 # Nozzle class: arbitrary area distributions, throat detection, validation
mesh.py                     # Grid generation, throat-normalized area array
initialization.py           # Analytical initial flow field
Analytical.py                # Closed-form isentropic / normal-shock quasi-1D solution
physics.py                  # Conservative <-> primitive conversions, flux evaluation, CFL time step
solver.py                    # MacCormack predictor-corrector step, artificial viscosity
boundary_conditions.py      # Inlet/outlet boundary treatments
convergence.py               # Residual tracking, convergence check, convergence plot
postprocess.py               # Result and numerical-vs-analytical comparison plots
requirements.txt             # numpy, matplotlib
```

---

## Getting Started

```bash
pip install -r requirements.txt
python main.py
```

`main.py` is the single entry point. The nozzle geometry, grid size `N`, CFL number, artificial viscosity coefficient `CX`, exit pressure (or `None` for a fully supersonic/shock-free exit), and convergence tolerance are all set as constants at the top of the file. Several alternative nozzle definitions (cosine contour, custom function, radius-based, point-based, CSV-based) are included as commented-out examples, so swapping geometries does not require touching any solver code.

Running `main.py` prints the detected throat location and equivalent wall half-angle, reports whether the run converged (and after how many iterations), and then produces three plots: the converged flow-field distributions, the numerical-vs-analytical Mach comparison with absolute and relative error, and the residual convergence history — along with the max/mean/RMS Mach-number error against the analytical solution for that run.

---

## Results

### Convergence History

The solver was marched until the RMS residual for continuity, momentum, and energy all fell below the prescribed tolerance of 10⁻⁸. The residuals decrease by several orders of magnitude; the oscillatory behavior visible along the way is characteristic of the transient approach to the steady solution, and the overall downward trend indicates stable convergence.

<img width="1669" height="831" alt="Convergence history for continuity, momentum, and energy residuals" src="https://github.com/user-attachments/assets/70707cff-bf36-41e3-8d30-511eadb6da7f" />

### Numerical vs. Analytical Validation

The converged Mach-number distribution is compared against the closed-form isentropic solution for the same area distribution. The two agree closely over the nozzle domain, confirming that the solver captures the expected subsonic-to-supersonic acceleration through the converging–diverging duct.

<img width="1669" height="831" alt="Numerical vs analytical Mach number, with absolute and relative error" src="https://github.com/user-attachments/assets/63db4fda-c52d-45aa-acd1-a3fe9af40fd9" />

The absolute Mach-number error remains small across the domain, with the largest deviations near the throat and in regions of rapid flow transition. The relative error is largest near the inlet, which is expected since the Mach number itself is small there, making the relative metric more sensitive to small absolute differences — it decreases substantially through the rest of the duct.

---

## Applications

A validated quasi-1D FDM nozzle model like this one is primarily a tool for fast, physically grounded **engineering judgment** — the kind of first-pass analysis used to size a design or sanity-check a more expensive one before committing resources to higher-fidelity work. Representative uses include:

- **Preliminary nozzle sizing and trade studies.** Rapidly sweeping throat area, expansion (exit-to-throat area) ratio, and overall contour for rocket, turbojet/turbofan, or ramjet exhaust nozzles, well before a 3D CAD model or mesh exists.
- **Choked mass flow rate estimation.** Computing the critical pressure ratio and resulting mass flow through a throat, relevant to flow metering devices (venturi/flow nozzles) and to sizing propulsion feed systems.
- **Design-condition and off-design screening.** Determining, for a given area ratio, the back pressure at which a nozzle runs shock-free (subsonic or supersonic design) versus the pressure range over which a normal shock sits inside the duct — informing where a design will under- or over-expand.
- **Operating-envelope / altitude-sensitivity studies.** For rocket nozzles, approximating how exit flow behavior (and hence a first-order estimate of thrust coefficient) shifts as ambient back pressure changes with altitude.
- **Supersonic wind-tunnel nozzle contour checks.** Verifying that a candidate area distribution reaches the intended test-section Mach number shock-free before refining the contour with a full method-of-characteristics design.
- **A fast baseline against 3D CFD.** Serving as a cheap, independent check on mass flow, throat Mach number, or exit conditions predicted by a full 2D/3D Navier–Stokes simulation of the same nozzle, to catch gross modeling or setup errors early.
- **Teaching and verification.** A compact, fully transparent implementation of a classic compressible-flow CFD problem (this exact setup follows the quasi-1D nozzle case study in Anderson's *Computational Fluid Dynamics: The Basics with Applications*), useful for learning finite-difference methods, artificial viscosity, and convergence behavior.

---

## Limitations of the Quasi-1D Approach

These are limitations of the **quasi-1D physical model itself**, not just this implementation — no amount of additional coding effort on this codebase removes them, because they follow from the governing assumptions.

- **No cross-sectional variation.** Flow properties are assumed uniform across each cross-section. The model cannot capture boundary-layer growth along the wall, flow separation, or non-uniform velocity/pressure profiles across the duct — all of which require at least a 2D (axisymmetric) treatment.
- **Inviscid and adiabatic.** There is no viscosity or heat transfer in the governing equations, so wall friction losses, viscous boundary-layer blockage of the effective flow area, and real discharge coefficients are not predicted. These must come from empirical correlations, boundary-layer corrections, or higher-fidelity viscous CFD.
- **Breaks down for strongly divergent geometries.** The quasi-1D assumption is only accurate for small equivalent wall half-angles; the geometry module itself flags nozzles whose wall angle grows large as geometries where "the quasi-1D assumption is questionable." Short, abruptly diverging, or non-smooth nozzle shapes fall outside the model's valid range.
- **Only a single, normal, duct-axis shock is resolved.** The analytical and numerical machinery can capture one normal shock standing perpendicular to the axis. It cannot predict oblique shocks, shock trains, Mach-disk structures in under-expanded jets, or shock/boundary-layer interaction — the analytical model explicitly labels over-expanded exhaust as "not resolved by 1D theory" once the shock structure moves outside the duct.
- **No true 3D nozzle-shape effects.** Asymmetric, scarfed, plug, or rectangular nozzle geometries, and any secondary or swirling flow structures, are outside what a 1D area distribution can represent.
- **Ideal, calorically perfect gas.** A constant ratio of specific heats is assumed throughout. Real-gas effects, combustion chemistry, variable specific heats with temperature, and two-phase (e.g. particulate-laden rocket exhaust) flows are not modeled.
- **Steady-state design tool, not a transient simulator.** The time-marching here is a pseudo-time method used to reach a steady solution; it does not represent genuine unsteady phenomena such as engine start-up/shutdown transients, flow oscillations, or combustion instability.
- **No coupling to the external flow field.** Once flow leaves the nozzle exit plane, plume expansion, base flow, and interaction with an external freestream are outside the scope of an internal-duct quasi-1D model.
- **Purely aerodynamic/thermodynamic.** The model says nothing about structural loads, thermal stresses, wall cooling requirements, or erosion — all of which are necessary for an actual engineering design and must be assessed separately.

In short: this tool is well suited to fast, first-pass sizing and engineering judgment, and to verifying that higher-fidelity tools are behaving sensibly — but it is not a substitute for viscous 2D/3D CFD, structural analysis, or experimental testing in a final design.

---

## Current Code Limitations and Future Work

Separate from the physical limitations above, there are implementation-specific items that could be extended:

- Validation is currently limited to shock-free and single-normal-shock cases; extending systematic validation to stronger transonic and multi-shock scenarios would be useful.
- Solution robustness still depends on the choice of CFL number, artificial viscosity coefficient, and grid resolution near the throat; a more automated robustness/sensitivity study across these parameters would strengthen the tool.
- The solver is single-threaded NumPy; for very fine grids or parameter sweeps, a vectorized batch-run mode would speed up trade studies.
- No automated test suite currently exercises the geometry, analytical, and solver modules against known closed-form cases.

---

## Reference

The quasi-1D nozzle formulation, MacCormack scheme, and artificial viscosity treatment used here follow the classic compressible-flow CFD case study in John D. Anderson, *Computational Fluid Dynamics: The Basics with Applications*.
