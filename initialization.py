import numpy as np
import matplotlib.pyplot as plt
from physics import primitive_to_conservative
from Analytical import Analytical

def initialize_flow(x, area, gamma, exit_pressure=None, nozzle=None, show_plot=False):
    #Geometry independent initial condition: the quasi-1D isentropic (or normal shock) solution
    #for the supplied area distribution, with reservoir conditions p0 = rho0 = T0 = 1
    analytical = Analytical(exit_pressure, area, gamma, x=x, nozzle=nozzle)
    solution = analytical.solve()

    mach = solution["mach"]
    temperature = solution["temperature"]
    rho = solution["density"]
    velocity = mach * np.sqrt(temperature)

    if show_plot:
        plt.plot(x, mach, label='INITIAL MACH NUMBER')
        plt.xlabel('Position')
        plt.ylabel('Mach Number')
        plt.title('Initial Flow Field')
        plt.legend()
        plt.grid(True)
        plt.show()

    U1, U2, U3 = primitive_to_conservative(rho, velocity, temperature, area, gamma)

    return U1, U2, U3
