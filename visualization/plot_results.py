import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def plot_scf_convergence(
    csv_path, label_size, line_width, marker_size, figsize, save_path=None
):
    df = pd.read_csv(csv_path)

    # Filter out iteration 0 for nicer plot
    df = df[df["dE"] > 0]

    fig, ax1 = plt.subplots(figsize=figsize)
    
    color = "tab:blue"
    ax1.set_xlabel("Iteration", fontsize=label_size)
    ax1.set_ylabel("Total Energy (Ha)", color=color, fontsize=label_size)
    ax1.plot(
        df["iter"],
        df["E_tot"],
        marker="o",
        color=color,
        label="$E_{tot}$",
        linewidth=line_width,
        markersize=marker_size,
    )
    ax1.tick_params(axis="y", labelcolor=color, labelsize=label_size * 0.8)
    ax1.tick_params(axis="x", labelsize=label_size * 0.8)
    ax1.grid(True, linestyle="--", alpha=0.6)


    ax2 = ax1.twinx()
    color = "tab:red"
    ax2.set_ylabel(r"$\Delta E$ (Ha)", color=color, fontsize=label_size)
    ax2.semilogy(
        df["iter"],
        df["dE"],
        marker="x",
        color=color,
        label=r"$\Delta E$",
        linewidth=line_width,
        markersize=marker_size,
    )
    ax2.tick_params(axis="y", labelcolor=color, labelsize=label_size * 0.8)

    plt.title("Self-Consistent Energy Convergence", fontsize=label_size * 1.2)
    fig.tight_layout()
    if save_path:
        plt.savefig(save_path, format="pdf", bbox_inches="tight")


def plot_radial_density(df, label_size, line_width, figsize, save_path=None):
    fig, ax = plt.subplots(figsize=figsize)

    r_plot = df["r"][1:]
    u_plot = df["u"][1:]

    ax.plot(
        r_plot,
        u_plot**2,
        label=r"$|u(r)|^2$ Probability Density",
        color="black",
        lw=line_width,
    )
    ax.plot(
        r_plot,
        u_plot,
        label=r"$u(r)$ Probability Amplitude",
        color="tab:gray",
        linestyle="--",
        alpha=0.8,
        lw=line_width,
    )

    ax.set_xlabel("r (Bohr)", fontsize=label_size)
    ax.set_ylabel(r"Radial Profile", fontsize=label_size)
    ax.set_title("Radial Probability Density and Amplitude", fontsize=label_size * 1.2)
    ax.grid(True, alpha=0.3)
    ax.set_xlim([-0.1, 5])
    ax.legend(fontsize=label_size * 0.8)
    ax.tick_params(labelsize=label_size * 0.8)

    fig.tight_layout()
    if save_path:
        plt.savefig(save_path, format="pdf", bbox_inches="tight")


def plot_potentials(df, label_size, line_width, figsize, save_path=None):
    fig, ax = plt.subplots(figsize=figsize)

    r_plot = df["r"][1:]

    ax.plot(r_plot, df["V_eff"][1:], label=r"$V_{eff}$", lw=line_width)
    ax.plot(r_plot, df["V_H"][1:], label=r"$V_{H}$", linestyle="--", lw=line_width)
    ax.plot(r_plot, df["V_x"][1:], label=r"$V_{x}$", linestyle=":", lw=line_width)
    ax.plot(r_plot, df["V_c"][1:], label=r"$V_{c}$", linestyle="-.", lw=line_width)

    # Also plot the bare nuclear potential
    Z = 2 
    nuc_pot = -Z / r_plot
    ax.plot(r_plot, nuc_pot, label="$-Z/r$", color="gray", alpha=0.5, lw=line_width)

    ax.set_xlabel("r (Bohr)", fontsize=label_size)
    ax.set_ylabel("Potential (Ha)", fontsize=label_size)
    ax.set_title("Effective Potential Components", fontsize=label_size * 1.2)
    ax.set_ylim([-4, 3.5])
    ax.set_xlim([-0.1, 5])
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=label_size * 0.8)
    ax.tick_params(labelsize=label_size * 0.8)

    fig.tight_layout()
    if save_path:
        plt.savefig(save_path, format="pdf", bbox_inches="tight")


if __name__ == "__main__":

    # ---------- SETTINGS ----------
    FIGSIZE_1 = (9, 6)
    FIGSIZE_2 = (8.2, 6)
    LABEL_SIZE = 14
    LINE_WIDTH = 2
    MARKER_SIZE = 6.0
    SAVE_PDF = True
    # ------------------------------

    script_dir = Path(__file__).parent
    base_dir = script_dir.parent
    outputs_dir = base_dir / "outputs"

    scf_log = outputs_dir / "scf_log.csv"
    profiles_dat = outputs_dir / "profiles_final.dat"

    if scf_log.exists():
        save_path = script_dir / "scf_convergence.pdf" if SAVE_PDF else None
        plot_scf_convergence(
            scf_log, LABEL_SIZE, LINE_WIDTH, MARKER_SIZE, FIGSIZE_1, save_path
        )
    else:
        print(f"File not found: {scf_log}")

    if profiles_dat.exists():
        with open(profiles_dat, "r") as f:
            header = f.readline().strip().replace("#", "").split()
        data = np.loadtxt(profiles_dat)
        df = pd.DataFrame(data, columns=header)

        save_path_density = script_dir / "radial_density.pdf" if SAVE_PDF else None
        save_path_potentials = script_dir / "potentials.pdf" if SAVE_PDF else None

        plot_radial_density(df, LABEL_SIZE, LINE_WIDTH, FIGSIZE_2, save_path_density)
        plot_potentials(df, LABEL_SIZE, LINE_WIDTH, FIGSIZE_2, save_path_potentials)
    else:
        print(f"File not found: {profiles_dat}")

    plt.show()
