import numpy as np
import matplotlib.pyplot as plt

from source.solver import RadialGrid, solve_shrodinger
from source.dft_potentials import get_V_h, get_V_x, get_V_c, get_TOTEN

print(
    r"""
╔════════════════════════════════════════════════════════════════════════════════════╗
║   ____  _____ _____   _   _      _ _                      _   _                    ║
║  |  _ \|  ___|_   _| | | | | ___| (_)_   _ _ __ ___      / \ | |_ ___  _ __ ___    ║
║  | | | | |_    | |   | |_| |/ _ \ | | | | | '_ ` _ \    / _ \| __/ _ \| '_ ` _ \   ║
║  | |_| |  _|   | |   |  _  |  __/ | | |_| | | | | | |  / ___ \ || (_) | | | | | |  ║
║  |____/|_|     |_|   |_| |_|\___|_|_|\__,_|_| |_| |_| /_/   \_\__\___/|_| |_| |_|  ║
║                                                                                    ║
╚════════════════════════════════════════════════════════════════════════════════════╝
"""
)

# --- GENERAL PARAMETERS --- #
Z = 2
convergence_threshold = 1e-4
mixing_alpha = 0.1

# --- GRID --- #
grid = RadialGrid(r_min=1e-6, r_max=10.0, h=1e-3)

# --- PARAMETERS FOR BISECT-PRELIMINARY TO BISECT --- #
E_search_range = (-3, -0.01)
E_rough_step = 0.1

# --- FIRST CALCULATION WITH V_eff=0 --- #
u_ind, E_root = solve_shrodinger(
    grid, Z, np.zeros(len(grid.r)), E_search_range, E_rough_step
)

print(f"Initial single electron eigenvalue: {E_root:.4f}")

# -------------------------------
#  --- SELF-CONSISTENT PART ---
# -------------------------------

V_H_calc = 2 * get_V_h(u_ind, grid.r)
V_X_calc = get_V_x(u_ind, grid.r)
V_C_calc, ec_calc = get_V_c(u_ind, grid.r)
V_eff_input = V_H_calc + V_X_calc + V_C_calc
TOTEN_old = get_TOTEN(E_root, u_ind, grid.r, V_H_calc, V_X_calc, V_C_calc, ec_calc)

print("\n═══ ENTERING SELF CONSISTENT LOOP ═══")
iteration = 1
while True:

    u_new, E_new = solve_shrodinger(grid, Z, V_eff_input, E_search_range, E_rough_step)
    V_H_new = 2 * get_V_h(u_new, grid.r)
    V_X_new = get_V_x(u_new, grid.r)
    V_C_new, ec_new = get_V_c(u_new, grid.r)
    V_eff_output = V_H_new + V_X_new + V_C_new

    TOTEN_new = get_TOTEN(E_new, u_new, grid.r, V_H_new, V_X_new, V_C_new, ec_new)

    E_diff = np.abs(TOTEN_new - TOTEN_old)
    print(f"Iteration {iteration:3d}: ΔE = {E_diff:.6e} | E_tot = {TOTEN_new:.6f}")

    if E_diff < convergence_threshold:
        break

    TOTEN_old = TOTEN_new
    V_eff_input = mixing_alpha * V_eff_output + (1 - mixing_alpha) * V_eff_input
    iteration += 1

u_sol = u_new
eigen_sol = E_new
TOTEN_sol = TOTEN_new

print("\n═══ RESULTS ═══")
print(f"Single electron eigenvalue E_1: {eigen_sol:.4f}")
print(f"Total energy TOTEN: {TOTEN_sol:.4f}")

print("\n" + "═" * 87)

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

# TODO: Make F depend on Z, so it is more elastic
# TODO: Make an yaml input file
# TODO: Make possible choso which potentials include (Note the facor 2 in front of hartree potential change)