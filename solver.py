import numpy as np

from physics import (conservative_to_primitive,compute_fluxes,compute_time_step)
from boundary_conditions import apply_boundary_conditions

def artificial_viscosity(U, pressure, Cx):
    pressure_diff = (pressure[2:]-2*pressure[1:-1]+pressure[:-2])
    pressure_ref = pressure_diff + 4 * pressure[1:-1]

    coeff = (Cx*np.abs(pressure_diff)/pressure_ref)

    viscosity = coeff * (U[2:]-2*U[1:-1]+U[:-2])

    return viscosity

def maccormack_step(U1,U2,U3,area,dx,gamma,CFL,Cx,exit_pressure):

    rho, velocity, temperature, pressure, mach = (conservative_to_primitive(U1,U2,U3,area,gamma))

    dt = compute_time_step(dx,velocity,temperature,CFL)
    log_area = np.log(area)
    F1, F2, F3 = compute_fluxes(U1,U2,U3,gamma)

    J2 = F2 - (U2**2 / U1)

    # predictor arrays
    U1_p = U1.copy()
    U2_p = U2.copy()
    U3_p = U3.copy()

    # forward difference
    dF1dx = (F1[2:] - F1[1:-1]) / dx
    dF2dx = (F2[2:] - F2[1:-1]) / dx
    dF3dx = (F3[2:] - F3[1:-1]) / dx

    dlogAdx = ((log_area[2:]-log_area[1:-1])/dx)

    dU1dt = -dF1dx
    dU2dt = -(dF2dx-J2[1:-1]*dlogAdx)
    dU3dt = -dF3dx

    S1 = artificial_viscosity(U1, pressure, Cx)
    S2 = artificial_viscosity(U2, pressure, Cx)
    S3 = artificial_viscosity(U3, pressure, Cx)

    U1_p[1:-1] = U1[1:-1] + dU1dt * dt + S1
    U2_p[1:-1] = U2[1:-1] + dU2dt * dt + S2
    U3_p[1:-1] = U3[1:-1] + dU3dt * dt + S3

    U1_p, U2_p, U3_p = apply_boundary_conditions(U1_p,U2_p,U3_p,area,gamma,exit_pressure)

    # predictor primitive variables
    rho_p, velocity_p, temperature_p, pressure_p, mach_p = (conservative_to_primitive(U1_p,U2_p,U3_p,area,gamma))

    F1_p, F2_p, F3_p = compute_fluxes(U1_p,U2_p,U3_p,gamma)
    J2_p = F2_p - (U2_p**2 / U1_p)

    # backward difference
    dF1dx_p = (F1_p[1:-1] - F1_p[:-2]) / dx
    dF2dx_p = (F2_p[1:-1] - F2_p[:-2]) / dx
    dF3dx_p = (F3_p[1:-1] - F3_p[:-2]) / dx

    dlogAdx_p = ((log_area[1:-1]-log_area[:-2])/dx)

    dU1dt_p = -dF1dx_p
    dU2dt_p = -(dF2dx_p-J2_p[1:-1]*dlogAdx_p)
    dU3dt_p = -dF3dx_p

    S1_p = artificial_viscosity(U1_p, pressure_p, Cx)
    S2_p = artificial_viscosity(U2_p, pressure_p, Cx)
    S3_p = artificial_viscosity(U3_p, pressure_p, Cx)

    U1[1:-1] += (0.5*(dU1dt+dU1dt_p)*dt+S1_p)
    U2[1:-1] += (0.5*(dU2dt+dU2dt_p)*dt+S2_p)
    U3[1:-1] += (0.5*(dU3dt+dU3dt_p)*dt+S3_p)

    U1, U2, U3 = apply_boundary_conditions(U1,U2,U3,area,gamma,exit_pressure)

    return U1, U2, U3

'''print(f"Exit pressure BC = {EXIT_PRESSURE}")

if EXIT_PRESSURE is None:
    print("Running supersonic outlet case")
else:
    print(
        f"Running pressure outlet case: {EXIT_PRESSURE}"
    )'''