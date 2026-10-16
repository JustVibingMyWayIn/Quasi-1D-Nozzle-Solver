import numpy as np


def conservative_to_primitive(U1, U2, U3, area, gamma):
    rho = U1 / area
    velocity = U2 / U1
    temperature = ((gamma - 1)*((U3/U1)-(gamma/2)*velocity**2))
    pressure = rho * temperature
    mach = velocity / np.sqrt(temperature)
    return rho, velocity, temperature, pressure, mach

def primitive_to_conservative(rho, velocity, temperature, area, gamma):
    U1 = rho * area
    U2 = rho * area * velocity
    U3 = rho*area*(temperature/(gamma-1)+(gamma/2)*velocity**2)
    return U1, U2, U3

def compute_fluxes(U1, U2, U3, gamma):
    F1 = U2
    F2 = ((U2**2)/U1+((gamma-1)/gamma)*(U3-(gamma/2)*(U2**2/U1)))
    F3 = (gamma*U2*U3/U1-(gamma*(gamma-1)/2)*U2*(U2/U1)**2)
    return F1, F2, F3

def compute_time_step(dx, velocity, temperature, CFL):
    dt = np.min(CFL*dx/(np.sqrt(temperature)+np.abs(velocity)))
    return dt

def check_physical(U1, U2, U3, area, gamma, x, iteration):
    #Stop with a clear message instead of propagating NaNs when the solution becomes unphysical
    rho, velocity, temperature, pressure, mach = conservative_to_primitive(U1, U2, U3, area, gamma)
    bad = (rho <= 0) | (temperature <= 0) | ~np.isfinite(U1) | ~np.isfinite(U2) | ~np.isfinite(U3)
    if np.any(bad):
        i = int(np.argmax(bad))
        raise FloatingPointError(
            f"Unphysical state at iteration {iteration}, x = {x[i]:.4f}: "
            f"rho = {rho[i]:.4g}, T = {temperature[i]:.4g}, u = {velocity[i]:.4g}. "
            f"Try a lower CFL, a higher CX, or a finer grid.")