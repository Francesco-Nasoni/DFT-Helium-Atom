# Density Functional Theory for Helium Atom

The goal of this project is to implement a **Self-Consistent Field (SCF)** solver to determine the ground state properties of the Helium atom using **Density Functional Theory (DFT)** by employing different approximation levels for the effective potential (Hartree only, Hartree + exchange, Hartree + exchange + correlation)

## Table Of Contents

1. [Dependencies](#dependencies)  
2. [Usage](#usage)  
3. [Configuration](#configuration)
4. [Outputs](#outputs)
5. [Background theory](#background-theory)

## Dependencies

### Core (required)

- Python 3.x
- numpy
- scipy
- pyyaml

### Optional (only for the provided plotting script)

- pandas
- matplotlib

## Usage

### 1) Configure parameters
Edit `config.yaml` to choose grid parameters, SCF thresholds/mixing, and whether to include exchange/correlation.

### 2) Run the solver
From the repository root:

```bash
python main.py
```

The program performs an initial “independent-electron” solve and then enters an SCF loop, writing outputs into the configured output directory.

### 3) (Optional) Quick plots
After a successful run, generate plots with the convenience script:

```bash
python visualization/plot_results.py
```

This reads the output files and produces figures (SCF convergence, radial density/amplitude, potential profiles). You can also analyze the generated files with any other tool (MATLAB, Julia, gnuplot, your own scripts, etc.).

## Configuration

All runtime parameters are controlled via `config.yaml`.

The parameters are:

### Grid
- `r_min`, `r_max`: radial domain
- `h`: grid step

### Eigenvalue search
- `E_min`, `E_max`: energy search interval
- `rough_step`: rough scan step to locate a sign change before bisection

### SCF
- `max_iter`: maximum iterations for the self consistent cycle
- `TOTEN_threshold`: convergence threshold based on `|ΔE_total|`
- `mix_alpha`: linear mixing parameter for the effective potential

### Model toggles
- `use_exchange`: enable LDA exchange potential
- `use_correlation`: enable LDA correlation potential (using Ceperley-Alder parameterization)

### Outputs
- `out_dir`: output directory
- `scf_log_csv`: SCF log filename
- `profiles_dat`: radial profiles filename

Implementation note: if `r_min` is set to `0`, it is internally replaced by a small positive value to avoid divisions by zero at the origin.

## Outputs

Two main files are produced under `out_dir` (default `outputs/`):

### 1) `scf_log.csv`
A CSV log of the SCF iterations, including quantities such as:
- `iter`
- `eps_1s` (single-particle eigenvalue)
- `E_tot` (total energy)
- `dE` (difference from previous iteration)

### 2) `profiles_final.dat`
A whitespace-separated table of final radial profiles, typically including columns such as:
- `r` (grid)
- `u` (radial amplitude)
- `V_H` (Hartree potential)
- `V_x` (exchange, if enabled)
- `V_c` (correlation, if enabled)
- `V_eff` (final effective potential)

Note: the names are set in the confing.yaml file. By default, they are `scf_log.csv` `profiles_final.dat`.

---

## Background Theory  
