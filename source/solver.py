import numpy as np
from scipy.optimize import bisect
from scipy.integrate import simpson

# --- GRID ---
class RadialGrid:
    def __init__(self, r_min=1e-6, r_max=10.0, h=1e-3):
        self.r_min = r_min
        self.r_max = r_max
        self.h = h
        self.r = np.arange(r_min, r_max, h)
        self.N = len(self.r)
        self.r_rev = np.flip(self.r)


# --- NUMERICS ---
def verlet_integrate_1D(fun, y_0, y_1, line_grid):
    h = np.abs(line_grid[1] - line_grid[0])
    y_integrated = np.zeros(len(line_grid))
    y_integrated[0] = y_0
    y_integrated[1] = y_1

    h2 = h**2
    for i in range(2, len(line_grid)):
        y_prev = y_integrated[i - 1]
        y_pprev = y_integrated[i - 2]
        f_prev = fun(line_grid[i - 1], y_prev, i - 1)
        y_next = 2 * y_prev - y_pprev + f_prev * h2
        y_integrated[i] = y_next

    return y_integrated


def get_norm_factor(u_r, line_grid):
    norm = simpson(u_r**2, line_grid)
    return np.sqrt(norm)


# --- SHRÖDINGER SOLVER ---
def F(E, Z, V_eff_r, r, u_r):
    return 2 * (- Z / r - E + V_eff_r) * u_r


def solve_shrodinger(grid, Z, V_eff, E_bounds, E_rough_step):
    h = grid.h
    u_rmax = grid.r_max * np.exp(-Z * grid.r_max)
    u_rmax_h = (grid.r_max - grid.h) * np.exp(-Z * grid.r_max + Z * h)

    # Flip V_eff to match the reversed grid for backward integration
    V_eff_rev = np.flip(V_eff)

    # Look for two points where u(0) chenges sign for the bisect method
    u_0 = verlet_integrate_1D(
        lambda r, ur, idx: F(E_bounds[0], Z, V_eff_rev[idx], r, ur),
        u_rmax,
        u_rmax_h,
        grid.r_rev,
    )
    u_prev = u_0

    E_bisect_range = None
    for E_test in np.arange(E_bounds[0] + E_rough_step, E_bounds[1], E_rough_step):
        u_i = verlet_integrate_1D(
            lambda r, ur, idx: F(E_test, Z, V_eff_rev[idx], r, ur),
            u_rmax,
            u_rmax_h,
            grid.r_rev,
        )
        if np.sign(u_prev[-1]) != np.sign(u_i[-1]):
            E_bisect_range = np.array([E_test, E_test - E_rough_step])
            break
        else:
            u_prev = u_i

    if E_bisect_range is None:
        raise ValueError(
            f"No sign change found in energy range {E_bounds}. Try expanding the search range or reducing E_rough_step."
        )

    # Proper bisect method
    E_root = bisect(
        lambda e: verlet_integrate_1D(
            lambda r, ur, idx: F(e, Z, V_eff_rev[idx], r, ur),
            u_rmax,
            u_rmax_h,
            grid.r_rev,
        )[-1],
        E_bisect_range[0],
        E_bisect_range[1],
        xtol=grid.h,
    )

    # NOTE: np.flip is necessary since we are integrating backwords r_max -> 0
    u_int = np.flip(
        verlet_integrate_1D(
            lambda r, ur, idx: F(E_root, Z, V_eff_rev[idx], r, ur),
            u_rmax,
            u_rmax_h,
            grid.r_rev,
        )
    )

    norm = get_norm_factor(u_int, grid.r)

    return u_int / norm, E_root