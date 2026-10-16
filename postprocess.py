import numpy as np
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

def plot_mach_comparison(x, mach_numerical, mach_analytical):
    error = np.abs(mach_numerical - mach_analytical)
    relerror = (error / (mach_analytical))*100

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

    ax1.plot(x, mach_numerical, 'o-', label="Numerical (MacCormack)", markersize=4)
    ax1.plot(x, mach_analytical, '--', label="Analytical (Isentropic)", linewidth=2)
    ax1.set_ylabel("Mach Number")
    ax1.set_title("Numerical vs Analytical Mach Number")
    ax1.legend()
    ax1.grid()

    ax2.plot(x, error, color='red', label="|Numerical - Analytical|")
    ax2.set_xlabel("x")
    ax2.set_ylabel("Absolute Error")
    ax2.legend()
    ax2.grid()

    ax3.plot(x, relerror, 'o-', label = "relative error in calculation", markersize=4, linewidth=2)
    ax3.set_xlabel("x")
    ax3.set_ylabel("Relative error")
    ax3.legend()
    ax3.grid()

    plt.tight_layout()
    plt.show()

    print(f"Max Mach error:  {np.max(error):.6f}")
    print(f"Mean Mach error: {np.mean(error):.6f}")
    print(f"RMS Mach error:  {np.sqrt(np.mean(error**2)):.6f}")