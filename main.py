import numpy as np
import matplotlib.pyplot as plt

from source.solver import RadialGrid, solve_shrodinger
from source.dft_potentials import get_V_h, get_V_x, get_V_c, get_V_eff, get_TOTEN

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
max_iterations = 70
use_exchange = False
use_correlation = False

# --- GRID --- #
grid = RadialGrid(r_min=1e-12, r_max=10.0, h=1e-3)

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

V_eff_input, V_H_0, V_X_0, V_C_0, ec_0 = get_V_eff(
    u_ind, grid.r, use_exchange, use_correlation
)
TOTEN_old = get_TOTEN(E_root, u_ind, grid.r, V_H_0, V_X_0, V_C_0, ec_0)

print("\n═══ ENTERING SELF CONSISTENT LOOP ═══")
iteration = 1
while iteration < max_iterations:

    u_new, E_new = solve_shrodinger(grid, Z, V_eff_input, E_search_range, E_rough_step)

    V_eff_output, V_H_new, V_X_new, V_C_new, ec_new = get_V_eff(
        u_new, grid.r, use_exchange, use_correlation
    )

    TOTEN_new = get_TOTEN(E_new, u_new, grid.r, V_H_new, V_X_new, V_C_new, ec_new)

    E_diff = np.abs(TOTEN_new - TOTEN_old)
    print(f"Iteration {iteration:3d}: ΔE= {E_diff:.6e} | E_1= {E_new:.6f} | E_tot= {TOTEN_new:.6f}")

    if E_diff < convergence_threshold:
        break

    TOTEN_old = TOTEN_new
    V_eff_input = mixing_alpha * V_eff_output + (1 - mixing_alpha) * V_eff_input
    iteration += 1


print("\n═══ RESULTS ═══")
print(f"Single electron eigenvalue E_1: {E_new:.6f}")
print(f"Total energy TOTEN: {TOTEN_new:.6f}")

print("\n" + "═" * 87)

plt.figure(figsize=(10, 6))
plt.plot(grid.r, u_new, label="u_sol")
plt.plot(grid.r, u_ind, label="u_ind")
plt.legend()
plt.xlim(-0.2, 6)
plt.xlabel("r")
plt.ylabel("u(r)")
plt.title("Integrated u(r)")
plt.grid(True)
plt.show()

# TODO: Make an yaml input file