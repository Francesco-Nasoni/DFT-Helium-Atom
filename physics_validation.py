import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.integrate import simpson

from source.io_utils import load_config_yaml
from source.solver import RadialGrid, solve_schrodinger
from source.dft_potentials import get_V_h


if __name__ == "__main__":
    
    # ================= PLOT SETTINGS ==================
    LABEL_SIZE = 14
    LINE_WIDTH = 2
    # ==================================================

    SCRIPT_DIR = Path(__file__).parent

    cfg = load_config_yaml("config.yaml")
    
    # --- GRID PARAMETERS--- #
    r_min = cfg["grid"]["r_min"]
    r_max = cfg["grid"]["r_max"]
    h = cfg["grid"]["h"]

    # --- BISECT PARAMETERS--- #
    E_search_range = (cfg["bisect"]["E_min"], cfg["bisect"]["E_max"])
    E_rough_step = cfg["bisect"]["rough_step"]
    
    grid = RadialGrid(r_min, r_max, h)
    V_eff_0 = np.zeros(len(grid.r))
    
    # Solve for hydrogenic helium (Z=2, independent electrons)
    u_ind, _ = solve_schrodinger(grid, 2, V_eff_0, E_search_range, E_rough_step)
    
    # Solve for hydrogen (Z=1)
    u_hyd, _ = solve_schrodinger(grid, 1, V_eff_0, [-1, 0], 0.1)
    
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # ---- LEFT PANEL: HARTREE POTENTIAL COMPARISON ----
    V_H_hyd = get_V_h(u_hyd, grid.r)
    
    V_H_analytical = (-(grid.r + 1) * np.exp(-2 * grid.r) + 1)/grid.r
    diff_hartree = np.max(np.abs(V_H_hyd - V_H_analytical))
    print(f"Hartree potential - Maximum absolute difference: {diff_hartree:.3e}")
    
    ax1.plot(
        grid.r, 
        V_H_hyd, 
        label="Numerical Hartree Potential",
        color="tab:blue",
        lw=LINE_WIDTH
    )
    ax1.plot(
        grid.r, 
        V_H_analytical, 
        label=r"Analytical: $(-(r+1)e^{-2r} + 1)/r$",
        color="tab:orange",
        linestyle="--",
        lw=LINE_WIDTH
    )
    ax1.set_xlabel("r (Bohr)", fontsize=LABEL_SIZE)
    ax1.set_ylabel(r"$V_H(r)$", fontsize=LABEL_SIZE)
    ax1.set_title(f"Hartree Potential Validation for Hydrogen\n(Max Abs. Diff: {diff_hartree:.3e})", 
                  fontsize=LABEL_SIZE * 1.1, pad=15)
    ax1.legend(fontsize=LABEL_SIZE * 0.8, loc="best")
    ax1.grid(True, alpha=0.3)
    ax1.tick_params(labelsize=LABEL_SIZE * 0.8)
    ax1.set_xlim([-0.1, 10])

    # ---- RIGHT PANEL: WAVEFUNCTION COMPARISON ----
    u_analytical = grid.r * np.exp(-2 * grid.r)
    Norm = simpson(u_analytical**2, grid.r)
    u_analytical = u_analytical / np.sqrt(Norm)
    
    diff_wavefunction = np.max(np.abs(u_ind - u_analytical))
    print(f"Wavefunction - Maximum absolute difference: {diff_wavefunction:.3e}")
    
    ax2.plot(
        grid.r, 
        u_ind, 
        label="Numerical Solution (Z=2)",
        color="tab:blue",
        lw=LINE_WIDTH
    )
    ax2.plot(
        grid.r, 
        u_analytical, 
        label=r"Analytical: $\frac{1}{\sqrt{N}}r e^{-2r}$",
        color="tab:orange",
        linestyle="--",
        lw=LINE_WIDTH
    )
    ax2.set_xlabel("r (Bohr)", fontsize=LABEL_SIZE)
    ax2.set_ylabel(r"$u(r)$", fontsize=LABEL_SIZE)
    ax2.set_title(f"Radial Wavefunction Comparison\n(Max Abs. Diff: {diff_wavefunction:.3e})", 
                  fontsize=LABEL_SIZE * 1.1, pad=15)
    ax2.legend(fontsize=LABEL_SIZE * 0.8, loc="best")
    ax2.grid(True, alpha=0.3)
    ax2.tick_params(labelsize=LABEL_SIZE * 0.8)
    ax2.set_xlim([-0.1, 5])
    
    # Save figure
    fig.tight_layout()
    save_path = SCRIPT_DIR / "physics_validation.png"
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"\nSaved combined validation plot to: {save_path}")
    
    plt.show()