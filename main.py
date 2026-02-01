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


class RadialGrid:
    def __init__(self, r_min=1e-6, r_max=10.0, h=1e-3):
        self.r_min = r_min
        self.r_max = r_max
        self.h = h
        self.r = np.arange(r_min, r_max, h)
        self.N = len(self.r)
        self.r_rev = np.flip(self.r)


def F(E, V_eff_r, r, u_r):
    return (-4 / r - 2 * E + 2 * V_eff_r) * u_r


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


def poisson_integrate(u_r, line_grid):
    h = np.abs(line_grid[1] - line_grid[0])

    def U_2_poisson(r, _, idx):
        return -(u_r[idx] ** 2) / r

    U_part = verlet_integrate_1D(U_2_poisson, 0, h, line_grid)
    alpha = (1 - U_part[-1]) / line_grid[-1]

    V_hartree = (U_part / line_grid) + alpha
    V_hartree[0] = V_hartree[1]
    return V_hartree


def solve_shrodinger(grid, Z, V_eff, E_bounds, E_rough_step):
    h = grid.h
    u_rmax = grid.r_max * np.exp(-Z * grid.r_max)
    u_rmax_h = (grid.r_max - grid.h) * np.exp(-Z * grid.r_max + Z * h)

    # Flip V_eff to match the reversed grid for backward integration
    V_eff_rev = np.flip(V_eff)

    # Look for two points where u(0) chenges sign for the bisect method
    u_0 = verlet_integrate_1D(
        lambda r, ur, idx: F(E_bounds[0], V_eff_rev[idx], r, ur),
        u_rmax,
        u_rmax_h,
        grid.r_rev,
    )
    u_prev = u_0

    E_bisect_range = None
    for E_test in np.arange(E_bounds[0] + E_rough_step, E_bounds[1], E_rough_step):
        u_i = verlet_integrate_1D(
            lambda r, ur, idx: F(E_test, V_eff_rev[idx], r, ur),
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
            lambda r, ur, idx: F(e, V_eff_rev[idx], r, ur),
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
            lambda r, ur, idx: F(E_root, V_eff_rev[idx], r, ur),
            u_rmax,
            u_rmax_h,
            grid.r_rev,
        )
    )

    norm = get_norm_factor(u_int, grid.r)

    return u_int / norm, E_root


# --- GENERAL PARAMETERS --- #
Z = 2
convergence_threshold = 1e-4

# --- GRID --- #
grid = RadialGrid(r_min=1e-6, r_max=10.0, h=1e-3)

# --- PARAMETERS FOR BISECT-PRELIMINARY TO BISECT --- #
E_search_range = (-3, 0)
E_rough_step = 0.1

# --- FIRST CALCULATION WITH V_eff=0 --- #
u_ind, E_root = solve_shrodinger(
    grid, Z, np.zeros(len(grid.r)), E_search_range, E_rough_step
)

print(f"Initial single electron eigenvalue: {E_root:.4f}")

# --- SELF-CONSISTENT PART ---#
u_old = u_ind
E_old = E_root
V_eff = np.zeros_like(grid.r)

print("\n=== ENTERING SELF CONSISTENT LOOP ===")
iteration = 1
while True:
    V_hartree_old = poisson_integrate(u_old, grid.r)
    TOTEN_old = 2 * E_old - simpson(V_hartree_old * u_old**2, grid.r)
    V_eff = V_hartree_old

    u_new, E_new = solve_shrodinger(grid, Z, V_eff, E_search_range, E_rough_step)
    V_hartree_new = poisson_integrate(u_new, grid.r)
    TOTEN_new = 2 * E_new - simpson(V_hartree_new * u_new**2, grid.r)

    E_diff = np.abs(TOTEN_new - TOTEN_old)
    print(f"Iteration {iteration}: ΔE = {E_diff}")

    if E_diff < convergence_threshold:
        break

    u_old = u_new
    E_old = E_new

    iteration += 1

u_sol = u_new
E_1_sol = E_new
TOTEN_sol = TOTEN_new

print("\n=== RESULTS WITH POTENTIALS ===")
print(f"Single electron eigenvalue E_1: {E_1_sol:.4f}")
print(f"Total energy TOTEN: {TOTEN_sol:.4f}")

plt.figure(figsize=(10, 6))
plt.plot(grid.r, u_sol, label="u_sol")
plt.plot(grid.r, u_ind, label="u_ind")
plt.legend()
plt.xlim(-0.2, 6)
plt.xlabel("r")
plt.ylabel("u(r)")
plt.title("Integrated u(r)")
plt.grid(True)
plt.show()
