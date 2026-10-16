#analytical solution
#Takes the mesh generated for this particular area ratio and provides theoretical Mach Number, Pressure, and Temperature Profiles

import numpy as np
from scipy.optimize import brentq

_UNSET = object()


class Analytical:
    def __init__(self, P_exit, area, gamma, x=None, nozzle=None):
        self.gamma = gamma
        self.P_critical = (2 / (self.gamma + 1)) ** (self.gamma / (self.gamma - 1))
        self.x = None if x is None else np.asarray(x, dtype=float) #Grid, needed to use the continuous geometry
        self.nozzle = nozzle #Optional geometry.Nozzle, gives the exact throat instead of the nearest grid point
        self._set_area(area)
        self.P_exit = P_exit

    def _uses_geometry(self):
        return self.nozzle is not None and self.x is not None

    def _set_area(self, area): #Takes the area distribution, and finds all important parameters
        self.area = np.asarray(area, dtype=float) #Area distribution along Nozzle
        self.n = len(self.area) #Number of points taken for computation
        self.throat_idx = int(np.argmin(self.area)) #Index of the grid point nearest the throat
        self.throat_pos = float(self.throat_idx) #Fractional index of the throat
        self.A_star = self.area[self.throat_idx] #Area of Throat

        if self._uses_geometry():
            #Exact throat from the continuous geometry, scaled into the same units as the supplied area array
            x_t, A_t = self.nozzle.throat()
            self._area_scale = self.area[0] / float(self.nozzle.area(self.x[0]))
            self.A_star = A_t * self._area_scale
            self.throat_pos = float(np.interp(x_t, self.x, np.arange(self.n)))

        self.P_sub_design, self.P_shock_exit, self.P_super_design = self._design_pressures() #Exit pressures for each flow regime

    def _area_at(self, s):
        if self._uses_geometry():
            x_s = np.interp(s, np.arange(self.n), self.x)
            return float(self.nozzle.area(x_s)) * self._area_scale #Exact area at the shock location
        return np.interp(s, np.arange(self.n), self.area) #Perform linear interpolation to find shock location, if flow regime is appropriate

    def area_ratio_from_mach(self, mach):
        term = (2 / (self.gamma + 1)) * (1 + (self.gamma - 1) / 2 * mach ** 2)
        return (1 / mach) * term ** ((self.gamma + 1) / (2 * (self.gamma - 1))) #Pretty self explanatory

    def mach_from_area(self, area_ratio, supersonic):
        if area_ratio <= 1.0 + 1e-12:
            return 1.0

        f = lambda M: self.area_ratio_from_mach(M) - area_ratio #Solving for Mach number known area ratio by Brent's method, which requires equailty with 0

        if supersonic:
            lo, hi = 1.0 + 1e-9, 50.0
            while f(hi) < 0:
                hi *= 2
            return brentq(f, lo, hi, xtol=1e-10)
        else:
            return brentq(f, 1e-9, 1.0 - 1e-9, xtol=1e-10)

    #Known gamma and mach, all properties can be ascertained from isentropic relations
    def pressure_ratio(self, mach):
        return (1 + (self.gamma - 1) / 2 * mach ** 2) ** (-self.gamma / (self.gamma - 1))

    def temperature_ratio(self, mach):
        return 1 / (1 + (self.gamma - 1) / 2 * mach ** 2)

    def density_ratio(self, mach):
        return (1 + (self.gamma - 1) / 2 * mach ** 2) ** (-1 / (self.gamma - 1))

    def normal_shock(self, M1):
        gamma = self.gamma
        M2 = np.sqrt((1 + (gamma - 1) / 2 * M1 ** 2) / (gamma * M1 ** 2 - (gamma - 1) / 2))
        p0_ratio = (
            ((gamma + 1) * M1 ** 2 / ((gamma - 1) * M1 ** 2 + 2)) ** (gamma / (gamma - 1))
            * ((gamma + 1) / (2 * gamma * M1 ** 2 - (gamma - 1))) ** (1 / (gamma - 1))
        )
        return M2, p0_ratio

    def _design_pressures(self):
        AR_exit = self.area[-1] / self.A_star

        M_e1 = self.mach_from_area(AR_exit, supersonic=False)
        P_sub_design = self.pressure_ratio(M_e1)

        M_e3 = self.mach_from_area(AR_exit, supersonic=True)
        P_super_design = self.pressure_ratio(M_e3)

        M2, p0_ratio = self.normal_shock(M_e3)
        P_shock_exit = self.pressure_ratio(M2) * p0_ratio

        return P_sub_design, P_shock_exit, P_super_design

    def flow_regime(self, P_exit_ratio):
        if P_exit_ratio is None:
            return "supersonic design (shock-free)"
        eps = 1e-9
        if P_exit_ratio >= self.P_sub_design - eps:
            return "subsonic (no shock)"
        if P_exit_ratio > self.P_shock_exit + eps:
            return "normal shock in duct"
        if P_exit_ratio > self.P_super_design + eps:
            return "overexpanded (oblique shocks outside duct, not resolved by 1D theory)"
        return "supersonic (shock-free; under-expanded if below design)"

    def _exit_pressure_for_shock_at(self, s):
        A_shock = self._area_at(s)
        M1 = self.mach_from_area(A_shock / self.A_star, supersonic=True)
        M2, p0_ratio = self.normal_shock(M1)
        A_star_downstream = A_shock / self.area_ratio_from_mach(M2)

        AR_exit_downstream = self.area[-1] / A_star_downstream
        M_exit = self.mach_from_area(AR_exit_downstream, supersonic=False)
        return self.pressure_ratio(M_exit) * p0_ratio

    def _find_shock_location(self, P_exit_ratio):
        f = lambda s: self._exit_pressure_for_shock_at(s) - P_exit_ratio
        lo = self.throat_pos + 1e-3
        hi = (self.n - 1) - 1e-3
        return brentq(f, lo, hi, xtol=1e-8)

    def solve(self, P_exit=_UNSET, area=None):
        if area is not None: #If area is provided, use the provided area
            self._set_area(area)
        if P_exit is not _UNSET:
            self.P_exit = P_exit
        '''
        Area is not provided by default, neither is exit pressure. 
        '''
        n = self.n 
        mach = np.zeros(n)
        P0_local = np.ones(n)
        shock_index = None
        regime = self.flow_regime(self.P_exit)

        if self.P_exit is None or regime != "normal shock in duct":
            if regime == "subsonic (no shock)": 
                f = lambda M: self.pressure_ratio(M) - self.P_exit
                M_exit = brentq(f, 1e-9, 1.0 - 1e-9, xtol=1e-10)
                A_star_eff = self.area[-1] / self.area_ratio_from_mach(M_exit) #A/(A/A*) = A*

                for i in range(n):
                    AR = self.area[i] / A_star_eff
                    mach[i] = self.mach_from_area(AR, supersonic=False)
            else:
                for i in range(n):
                    AR = self.area[i] / self.A_star
                    mach[i] = self.mach_from_area(AR, supersonic=(i > self.throat_pos))
        else: #Normal shock finder
            s_shock = self._find_shock_location(self.P_exit)
            shock_index = s_shock
            shock_cell = int(np.ceil(s_shock))

            A_shock = self._area_at(s_shock)
            M1 = self.mach_from_area(A_shock / self.A_star, supersonic=True)
            M2, p0_ratio = self.normal_shock(M1)
            A_star_downstream = A_shock / self.area_ratio_from_mach(M2)

            for i in range(n):
                if i <= self.throat_pos:
                    AR = self.area[i] / self.A_star
                    mach[i] = self.mach_from_area(AR, supersonic=False)
                elif i < shock_cell:
                    AR = self.area[i] / self.A_star
                    mach[i] = self.mach_from_area(AR, supersonic=True)
                else:
                    AR = self.area[i] / A_star_downstream
                    mach[i] = self.mach_from_area(AR, supersonic=False)
                    P0_local[i] = p0_ratio

        pressure = P0_local * self.pressure_ratio(mach)
        temperature = self.temperature_ratio(mach)
        density = P0_local * self.density_ratio(mach)

        return {
            "mach": mach,
            "pressure": pressure,
            "temperature": temperature,
            "density": density,
            "regime": regime,
            "shock_index": shock_index,
        }

    def mach_number(self, P_exit, area):
        return self.solve(P_exit=P_exit, area=area)["mach"]

    def pressure_distribution(self, mach):
        return self.pressure_ratio(mach)

    def temperature_distribution(self, mach):
        return self.temperature_ratio(mach)

    def density_distribution(self, mach, pressure=None, temperature=None, gamma=None):
        return self.density_ratio(mach)