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
    dt = np.min(CFL*dx/(np.sqrt(temperature)+velocity))
    return dt