import numpy as np
from pathlib import Path

from source.solver import RadialGrid, solve_shrodinger
from source.dft_potentials import get_V_eff, get_TOTEN
from source.io_utils import load_config_yaml, append_csv, save_profiles

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

    # --- GRID PARAMETERS--- #
    r_min = cfg["grid"]["r_min"]
    r_max = cfg["grid"]["r_max"]
    h = cfg["grid"]["h"]

    # --- BISECT PARAMETERS--- #
    E_search_range = (cfg["bisect"]["E_min"], cfg["bisect"]["E_max"])
    E_rough_step = cfg["bisect"]["rough_step"]

    # --- SCF PARAMETERS --- #
    convergence_threshold = cfg["scf"]["TOTEN_threshold"]
    mixing_alpha = cfg["scf"]["mix_alpha"]
    max_iterations = cfg["scf"]["max_iter"]
    use_exchange = cfg["model"]["use_exchange"]
    use_correlation = cfg["model"]["use_correlation"]

    # --- OUTPUT --- #
    out_dir = Path(cfg["output"]["out_dir"])
    scf_log_file_name = cfg["output"]["scf_log_csv"]
    profiles_save_file_name = cfg["output"]["profiles_dat"]

    # Delete existing log and profiles files if they exist
    (out_dir / scf_log_file_name).unlink(missing_ok=True)
    (out_dir / profiles_save_file_name).unlink(missing_ok=True)
    
    # Create output_dir if it does't exist
    out_dir.mkdir(exist_ok=True)

    # Define the grid
    grid = RadialGrid(r_min, r_max, h)

    # --- FIRST CALCULATION (V_eff=0) --- #
    V_eff_0 = np.zeros(len(grid.r))
    u_ind, E_root = solve_shrodinger(
        grid, Z, V_eff_0, E_search_range, E_rough_step
    )

    print("\n═══ INDEPENDENT ELECTRON RESULTS ═══")
    print(f"Single electron eigenvalue E_1: {E_root:.3f}")
    print(f"Total energy TOTEN: {2*E_root:.3f}")

    # ------------------------------- #
    #  --- SELF-CONSISTENT PART ---
    # ------------------------------- #

    V_eff_input, V_H_0, V_X_0, V_C_0, ec_0 = get_V_eff(
        u_ind, grid.r, use_exchange, use_correlation
    )
    TOTEN_old = get_TOTEN(E_root, u_ind, grid.r, V_H_0, V_X_0, V_C_0, ec_0)

    append_csv(out_dir/scf_log_file_name, 0, E_root, 2*E_root, 0)

    print("\n═══ ENTERING SELF CONSISTENT LOOP ═══")
    it = 1
    while it < max_iterations:

        u_new, E_new = solve_shrodinger(grid, Z, V_eff_input, E_search_range, E_rough_step)
        V_eff_output, V_H_new, V_X_new, V_C_new, ec_new = get_V_eff(
            u_new, grid.r, use_exchange, use_correlation
        )
        TOTEN_new = get_TOTEN(E_new, u_new, grid.r, V_H_new, V_X_new, V_C_new, ec_new)

        E_diff = np.abs(TOTEN_new - TOTEN_old)
        print(f"Iteration {it:3d}: ΔE= {E_diff:.6e} | E_1= {E_new:.6f} | E_tot= {TOTEN_new:.6f}")
        append_csv(out_dir/scf_log_file_name, it, E_new, TOTEN_new, E_diff)

        if E_diff < convergence_threshold:
            break

        TOTEN_old = TOTEN_new
        V_eff_input = mixing_alpha * V_eff_output + (1 - mixing_alpha) * V_eff_input
        it += 1

    save_profiles(out_dir/profiles_save_file_name, grid.r, u=u_new, V_H=V_H_new, V_x=V_X_new, V_c=V_C_new, V_eff=V_eff_output)
    print("\n═══ FINAL RESULTS ═══")
    print(f"Single electron eigenvalue E_1: {E_new:.6f}")
    print(f"Total energy TOTEN: {TOTEN_new:.6f}")
    print(f"\nPotential and radial probability amplitude were saved in '{out_dir}' directory")
    print("═" * 86)


if __name__ == "__main__":
    main()