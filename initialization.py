import numpy as np
from physics import primitive_to_conservative

def initialize_flow(x, area, gamma):
    rho = np.zeros_like(x)
    temperature = np.zeros_like(x)

    mask1 = x <= 0.5
    mask2 = (x > 0.5) & (x <= 1.5)
    mask3 = (x > 1.5) & (x <= 2.1)
    mask4 = x > 2.1

    rho[mask1] = 1.0
    temperature[mask1] = 1.0
    rho[mask2] = 1 - 0.366 * (x[mask2] - 0.5)
    temperature[mask2] = 1 - 0.167 * (x[mask2] - 0.5)

    rho[mask3] = 0.634 - 0.702 * (x[mask3] - 1.5)
    temperature[mask3] = 0.833 - 0.490 * (x[mask3] - 1.5)

    rho[mask4] = 0.589 + 0.102 * (x[mask4] - 2.1)
    temperature[mask4] = 0.94 + 0.062 * (x[mask4] - 2.1)

    velocity = 0.58 / (rho * area)

    U1, U2, U3 = primitive_to_conservative(rho, velocity, temperature, area, gamma)

    return U1, U2, U3