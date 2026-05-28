def apply_boundary_conditions(U1, U2, U3, area, gamma, exit_pressure):
    # inlet
    U1[0] = area[0]
    U2[0] = 2 * U2[1] - U2[2]
    U3[0] = area[0] * (1/(gamma-1)+(gamma/2)*(U2[0]/U1[0])**2)

    # outlet
    U1[-1] = 2 * U1[-2] - U1[-3]
    U2[-1] = 2 * U2[-2] - U2[-3]
    U3[-1] = (exit_pressure*area[-1]/(gamma-1)+(gamma/2)*area[-1]*(U2[-1]/U1[-1])**2)

    return U1, U2, U3