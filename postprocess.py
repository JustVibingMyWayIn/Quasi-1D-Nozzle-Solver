import matplotlib.pyplot as plt
from physics import conservative_to_primitive

def plot_results(x,area,U1,U2,U3,gamma):

    rho, velocity, temperature, pressure, mach = (conservative_to_primitive(U1,U2,U3,area,gamma))

    plt.figure(figsize=(10, 6))

    plt.plot(x, pressure, label="Pressure")
    plt.plot(x, temperature, label="Temperature")
    plt.plot(x, mach, label="Mach")
    plt.plot(x, U2, label="Mass Flow")
    plt.plot(x, area, label="Area")

    plt.xlabel("x")
    plt.legend()
    plt.grid()

    plt.show()