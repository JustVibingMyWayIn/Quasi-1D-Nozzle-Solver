#Nozzle geometry definition
#Wraps an arbitrary area distribution A(x) so the rest of the solver does not depend on a particular shape

import warnings

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar


class Nozzle:
    def __init__(self, area_func, x_start, x_end, dadx_func=None, name="custom"):
        if x_end <= x_start:
            raise ValueError("x_end must be greater than x_start")
        self._area_func = area_func
        self._dadx_func = dadx_func #Optional exact derivative, otherwise central differences are used
        self.x_start = float(x_start)
        self.x_end = float(x_end)
        self.length = self.x_end - self.x_start
        self.name = name
        self._throats = None

    # ---------- constructors ----------

    @classmethod
    def from_function(cls, area_func, x_start, x_end, dadx_func=None, name="custom"):
        return cls(area_func, x_start, x_end, dadx_func, name)

    @classmethod
    def from_radius(cls, radius_func, x_start, x_end, name="custom (radius)"):
        return cls(lambda x: np.pi * radius_func(x) ** 2, x_start, x_end, name=name)

    @classmethod
    def from_points(cls, x_points, area_points, name="custom (points)"):
        x_points = np.asarray(x_points, dtype=float)
        area_points = np.asarray(area_points, dtype=float)
        spline = CubicSpline(x_points, area_points) #C2 continuous fit through the supplied points
        return cls(spline, x_points[0], x_points[-1], spline.derivative(), name)

    @classmethod
    def from_csv(cls, path, delimiter=",", skiprows=0, name=None):
        data = np.loadtxt(path, delimiter=delimiter, skiprows=skiprows)
        return cls.from_points(data[:, 0], data[:, 1], name or f"csv: {path}")

    @classmethod
    def parabolic(cls, length=3.0, k=2.2):
        #Anderson's nozzle: A = 1 + k(x - L/2)^2
        xt = 0.5 * length
        return cls(lambda x: 1 + k * (x - xt) ** 2, 0.0, length,
                   lambda x: 2 * k * (x - xt), name="parabolic")

    @classmethod
    def cosine(cls, length=3.0, inlet_ratio=4.0, exit_ratio=3.0, throat_fraction=0.4):
        #Smooth converging-diverging nozzle built from two half cosine waves meeting at the throat
        xt = throat_fraction * length

        def area(x):
            x = np.asarray(x, dtype=float)
            conv = 1 + 0.5 * (inlet_ratio - 1) * (1 + np.cos(np.pi * x / xt))
            div = 1 + 0.5 * (exit_ratio - 1) * (1 - np.cos(np.pi * (x - xt) / (length - xt)))
            return np.where(x <= xt, conv, div)

        return cls(area, 0.0, length, name="cosine")

    # ---------- evaluation ----------

    def area(self, x):
        return np.asarray(self._area_func(np.asarray(x, dtype=float)), dtype=float)

    def dadx(self, x):
        x = np.asarray(x, dtype=float)
        if self._dadx_func is not None:
            return np.asarray(self._dadx_func(x), dtype=float)
        h = 1e-6 * self.length
        return (self.area(x + h) - self.area(x - h)) / (2 * h)

    # ---------- throat detection ----------

    def throats(self, samples=4001):
        #Returns every local minimum of A(x) as (x, A), refined with a bounded scalar minimisation
        if self._throats is not None:
            return self._throats

        xs = np.linspace(self.x_start, self.x_end, samples)
        a = self.area(xs)
        found = []

        for i in range(1, samples - 1):
            if a[i] <= a[i - 1] and a[i] < a[i + 1]:
                res = minimize_scalar(lambda s: float(self.area(s)),
                                      bounds=(xs[i - 1], xs[i + 1]), method="bounded",
                                      options={"xatol": 1e-12})
                found.append((float(res.x), float(res.fun)))

        #Minimum located at an end of the domain (purely converging or diverging nozzle)
        if a[0] < a[1]:
            found.insert(0, (self.x_start, float(a[0])))
        if a[-1] < a[-2]:
            found.append((self.x_end, float(a[-1])))

        if not found:
            i = int(np.argmin(a))
            found.append((float(xs[i]), float(a[i])))

        self._throats = found
        return found

    def throat(self):
        #Global minimum, i.e. the sonic throat for choked flow
        return min(self.throats(), key=lambda t: t[1])

    # ---------- checks ----------

    def validate(self, n=None, max_half_angle_deg=45.0, interp_tol=0.01, samples=4001):
        xs = np.linspace(self.x_start, self.x_end, samples)
        a = self.area(xs)

        if not np.all(np.isfinite(a)):
            raise ValueError(f"[{self.name}] area function returned non-finite values")
        if np.any(a <= 0):
            bad = xs[np.argmax(a <= 0)]
            raise ValueError(f"[{self.name}] area must be positive everywhere (A <= 0 at x = {bad:.4f})")

        throats = self.throats(samples)
        x_t, _ = self.throat()

        interior = [t for t in throats if self.x_start < t[0] < self.x_end]
        if len(interior) > 1:
            locs = ", ".join(f"x = {t[0]:.4f}" for t in interior)
            warnings.warn(f"[{self.name}] multiple throats found ({locs}); the steady solution may not be "
                          f"unique and the analytical comparison assumes a single throat")
        if x_t in (self.x_start, self.x_end):
            warnings.warn(f"[{self.name}] minimum area lies at the domain boundary (x = {x_t:.4f}); "
                          f"no converging-diverging throat, choked/supersonic boundary conditions do not apply")

        #Wall angle of the equivalent circular duct: dr/dx = A' / (2 sqrt(pi A))
        slope = np.abs(self.dadx(xs)) / (2 * np.sqrt(np.pi * a))
        half_angle = np.degrees(np.arctan(np.max(slope)))
        if half_angle > max_half_angle_deg:
            warnings.warn(f"[{self.name}] equivalent wall half-angle reaches {half_angle:.1f} deg; "
                          f"the quasi-1D assumption is questionable")

        if n is not None:
            #Check the grid resolves the geometry: linear interpolation error at cell midpoints
            x = np.linspace(self.x_start, self.x_end, n)
            ax = self.area(x)
            xm = 0.5 * (x[1:] + x[:-1])
            err = np.max(np.abs(self.area(xm) - 0.5 * (ax[1:] + ax[:-1])) / self.area(xm))
            if err > interp_tol:
                warnings.warn(f"[{self.name}] N = {n} under-resolves the geometry (max midpoint "
                              f"interpolation error {100 * err:.2f}%); consider increasing N")

        return {"throats": throats, "throat": (x_t, self.throat()[1]), "max_half_angle_deg": half_angle}
