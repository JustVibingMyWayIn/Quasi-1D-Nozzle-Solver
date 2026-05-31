def apply_boundary_conditions(U1, U2, U3, area, gamma, exit_pressure):
    if exit_pressure is None:
        return apply_supersonic_boundary_conditions(U1, U2, U3, area, gamma)
    return apply_subsonic_boundary_conditions(U1, U2, U3, area, gamma, exit_pressure)

def apply_subsonic_boundary_conditions(U1, U2, U3, area, gamma, exit_pressure):
    # inlet
    U1[0] = area[0]
    U2[0] = 2 * U2[1] - U2[2]
    U3[0] = area[0] * (1/(gamma-1)+(gamma/2)*(U2[0]/U1[0])**2)

    # outlet
    U1[-1] = 2 * U1[-2] - U1[-3]
    U2[-1] = 2 * U2[-2] - U2[-3]
    U3[-1] = (exit_pressure*area[-1]/(gamma-1)+(gamma/2)*area[-1]*(U2[-1]/U1[-1])**2)

    return U1, U2, U3

def apply_supersonic_boundary_conditions(U1, U2, U3, area, gamma):
    U1[0] = area[0]
    U2[0] = 2* U2[1] - U2[2]
    U3[0] = area[0] * (1/(gamma-1)+(gamma/2)*(U2[0]/U1[0])**2)

    U1[-1] = 2 * U1[-2] - U1[-3]
    U2[-1] = 2 * U2[-2] - U2[-3]
    U3[-1] = 2 * U3[-2] - U3[-3]

    return U1, U2, U3