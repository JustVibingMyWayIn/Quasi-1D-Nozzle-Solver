def apply_boundary_conditions(U1, U2, U3, area, gamma, exit_pressure):
    if exit_pressure is None:
        return apply_supersonic_boundary_conditions(U1, U2, U3, area, gamma)
    return apply_subsonic_boundary_conditions(U1, U2, U3, area, gamma, exit_pressure)

def apply_inlet_boundary_condition(U1, U2, U3, area, gamma):
    #Subsonic inlet fed from a reservoir: p0 = T0 = 1 are held fixed and velocity is extrapolated.
    #Static T and rho then follow isentropically, which stays valid when the inlet Mach number is not small.
    u1 = U2[1] / U1[1]
    u2 = U2[2] / U1[2]
    u0 = 2 * u1 - u2

    T0 = 1 - 0.5 * (gamma - 1) * u0**2
    rho0 = T0 ** (1 / (gamma - 1))

    U1[0] = rho0 * area[0]
    U2[0] = U1[0] * u0
    U3[0] = U1[0] * (T0/(gamma-1)+(gamma/2)*u0**2)

def apply_subsonic_boundary_conditions(U1, U2, U3, area, gamma, exit_pressure):
    # inlet
    apply_inlet_boundary_condition(U1, U2, U3, area, gamma)

    # outlet
    U1[-1] = 2 * U1[-2] - U1[-3]
    U2[-1] = 2 * U2[-2] - U2[-3]
    U3[-1] = (exit_pressure*area[-1]/(gamma-1)+(gamma/2)*U2[-1]*(U2[-1]/U1[-1]))

    return U1, U2, U3

def apply_supersonic_boundary_conditions(U1, U2, U3, area, gamma):
    apply_inlet_boundary_condition(U1, U2, U3, area, gamma)

    U1[-1] = 2 * U1[-2] - U1[-3]
    U2[-1] = 2 * U2[-2] - U2[-3]
    U3[-1] = 2 * U3[-2] - U3[-3]

    return U1, U2, U3
