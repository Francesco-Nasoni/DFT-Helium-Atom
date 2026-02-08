import numpy as np
import matplotlib.pyplot as plt

from source.solver import RadialGrid, solve_shrodinger
from source.dft_potentials import get_V_eff, get_TOTEN
from source.io_utils import load_config_yaml

def main():
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

    cfg = load_config_yaml("config.yaml")
    Z = 2

    # --- GRID --- #
    r_min = cfg["grid"]["r_min"]
    r_max = cfg["grid"]["r_max"]
    h = cfg["grid"]["h"]

    # --- BISECT --- #
    E_search_range = (cfg["bisect"]["E_min"], cfg["bisect"]["E_max"])
    E_rough_step = cfg["bisect"]["rough_step"]

    # --- SCF PARAMETERS --- #
    convergence_threshold = cfg["scf"]["TOTEN_threshold"]
    mixing_alpha = cfg["scf"]["mix_alpha"]
    max_iterations = cfg["scf"]["max_iter"]
    use_exchange = cfg["model"]["use_exchange"]
    use_correlation = cfg["model"]["use_correlation"]


    grid = RadialGrid(r_min, r_max, h)

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

if __name__ == "__main__":
    main()