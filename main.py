#   ____  _____ _____   _   _      _ _                      _   _                   #
#  |  _ \|  ___|_   _| | | | | ___| (_)_   _ _ __ ___      / \ | |_ ___  _ __ ___   #
#  | | | | |_    | |   | |_| |/ _ \ | | | | | '_ ` _ \    / _ \| __/ _ \| '_ ` _ \  #
#  | |_| |  _|   | |   |  _  |  __/ | | |_| | | | | | |  / ___ \ || (_) | | | | | | #
#  |____/|_|     |_|   |_| |_|\___|_|_|\__,_|_| |_| |_| /_/   \_\__\___/|_| |_| |_| #
#                                                                                   #


import numpy as np
from tqdm import tqdm
from scipy.optimize import bisect
from scipy.integrate import simpson
import matplotlib.pyplot as plt


def F(E, r, u_r, V_eff=0):
    return (-4 / r - 2 * E + 2 * V_eff) * u_r


def verlet_integrate_1D(fun, y_0, y_1, line_grid):
    h = np.abs(line_grid[1] - line_grid[0])
    y_integrated = np.zeros(len(line_grid))
    y_integrated[0] = y_0
    y_integrated[1] = y_1

    h2 = h**2
    for i in range(2, len(line_grid)):
        y_prev = y_integrated[i - 1]
        y_pprev = y_integrated[i - 2]
        y_next = 2 * y_prev - y_pprev + fun(line_grid[i - 1], y_prev) * h2
        y_integrated[i] = y_next

    return y_integrated


def get_norm_factor(u_r, line_grid):
    norm = simpson(u_r**2, line_grid) * 4 * np.pi
    return np.sqrt(norm)


def poisson_integrate(u_r, line_grid):
    h = np.abs(line_grid[1] - line_grid[0])

    def U_2_poisson(r, ur):
        return -(ur**2) / r

    U_part = verlet_integrate_1D(U_2_poisson, 0, h, line_grid)
    alpha = (2 - U_part[-1]) / line_grid[-1]

    V_hartree = U_part / line_grid + alpha
    return V_hartree


class RadialGrid:
    def __init__(self, r_min=1e-6, r_max=10.0, h=1e-3):
        self.r_min = r_min
        self.r_max = r_max
        self.h = h
        self.r = np.arange(r_min, r_max, h)
        self.N = len(self.r)
        self.r_rev = np.flip(self.r)


# --- PARAMETERS ---
E_thr = 1e-3

# --- GRID ---
grid = RadialGrid(r_min=1e-6, r_max=10.0, h=1e-3)

# --- PARAMETERS FOR BISECT-PRELIMINARY TO BISECT
E_range_search = (-3, 0)
E_step = 0.05

u_rmax = grid.r_max * np.exp(-2 * grid.r_max)
u_rmax_h = (grid.r_max - grid.h) * np.exp(-2 * grid.r_max + grid.h)

# Look for two points where u_int chenges sign since we want the value of E
# for which u_int=0

u_0 = verlet_integrate_1D(
    lambda r, ur: F(E_range_search[0], r, ur), u_rmax, u_rmax_h, grid.r_rev
)
u_prev = u_0

for E_test in np.arange(E_range_search[0] + E_step, E_range_search[1], E_step):
    u_i = verlet_integrate_1D(
        lambda r, ur: F(E_test, r, ur), u_rmax, u_rmax_h, grid.r_rev
    )
    if np.sign(u_prev[-1]) != np.sign(u_i[-1]):
        E_bisect_range = np.array([E_test, E_test - E_step])
        break
    else:
        u_prev = u_i

E_root = bisect(
    lambda e: verlet_integrate_1D(
        lambda r, ur: F(e, r, ur), u_rmax, u_rmax_h, grid.r_rev
    )[-1],
    E_bisect_range[0],
    E_bisect_range[1],
    xtol=grid.h,
)

# NOTE: np.flip is necessary since we are integrating backwords r_max -> 0
u_int = np.flip(
    verlet_integrate_1D(lambda r, ur: F(E_root, r, ur), u_rmax, u_rmax_h, grid.r_rev)
)

print(f"E_root: {E_root}")

plt.figure(figsize=(10, 6))
plt.plot(grid.r, u_int)
plt.xlim(-0.2, 6)
plt.xlabel("r")
plt.ylabel("u(r)")
plt.title("Integrated u(r)")
plt.grid(True)
plt.show()
