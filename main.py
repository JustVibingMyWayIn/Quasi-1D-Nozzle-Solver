import numpy as np

from geometry import Nozzle
from mesh import generate_mesh
from initialization import initialize_flow
from solver import maccormack_step
from postprocess import plot_results, plot_mach_comparison
from physics import conservative_to_primitive, check_physical
from Analytical import Analytical
from convergence import new_residual_history, update_residual_history, has_converged, plot_convergence

# constants
GAMMA = 1.4
N = 61
CFL = 0.5
CX = 0.2
ITERATIONS = 5000
EXIT_PRESSURE = None
CONVERGENCE_TOL = 1e-8  # stop early once the RMS residual drops below this

# geometry: any positive, smooth A(x) can be supplied here
NOZZLE = Nozzle.parabolic(length=3.0)
# NOZZLE = Nozzle.cosine(length=3.0, inlet_ratio=4.0, exit_ratio=3.0, throat_fraction=0.4)
# NOZZLE = Nozzle.from_function(lambda x: 1 + 3.0*np.sin(np.pi*(x - 1.0)/2.0)**2 * (x < 1.0) + 2.0*np.sin(np.pi*(x - 1.0)/4.0)**2 * (x >= 1.0), 0.0, 3.0)
# NOZZLE = Nozzle.from_radius(lambda x: 0.6 + 0.25*(x - 1.2)**2, 0.0, 3.0)
# NOZZLE = Nozzle.from_points([0.0, 0.75, 1.5, 2.25, 3.0], [5.0, 2.0, 1.0, 2.0, 3.0])
# NOZZLE = Nozzle.from_csv("my_nozzle.csv")  # columns: x, A

geometry_info = NOZZLE.validate(n=N)
print(f"Nozzle: {NOZZLE.name}, throat at x = {geometry_info['throat'][0]:.4f}, "
      f"max equivalent wall half-angle = {geometry_info['max_half_angle_deg']:.1f} deg")

# mesh
x, dx, area = generate_mesh(NOZZLE, N)

# initialize
U1, U2, U3 = initialize_flow(x, area, GAMMA, EXIT_PRESSURE, nozzle=NOZZLE)

# solve
residual_history = new_residual_history()

for i in range(ITERATIONS):
    U1_old, U2_old, U3_old = U1.copy(), U2.copy(), U3.copy()

    U1, U2, U3 = maccormack_step(U1,U2,U3,area,dx,GAMMA,CFL,CX,EXIT_PRESSURE)

    check_physical(U1, U2, U3, area, GAMMA, x, i + 1)

    update_residual_history(residual_history, U1_old, U1, U2_old, U2, U3_old, U3)

    if has_converged(residual_history, CONVERGENCE_TOL):
        print(f"Converged after {i + 1} iterations (residual < {CONVERGENCE_TOL:.0e})")
        break
else:
    print(f"Reached max iterations ({ITERATIONS}) without meeting tolerance {CONVERGENCE_TOL:.0e}")

# extract final numerical Mach profile
rho, velocity, temperature, pressure, mach_numerical = conservative_to_primitive(U1,U2,U3,area,GAMMA)

# analytical (isentropic, shock-free) Mach profile for comparison
analytical = Analytical(EXIT_PRESSURE, area, GAMMA, x=x, nozzle=NOZZLE)
mach_analytical = analytical.mach_number(EXIT_PRESSURE, area)

# plots
plot_results(x,area,U1,U2,U3,GAMMA)
plot_mach_comparison(x, mach_numerical, mach_analytical)
plot_convergence(residual_history, tol=CONVERGENCE_TOL)
