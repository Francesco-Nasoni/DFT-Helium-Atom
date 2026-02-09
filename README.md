# Density Functional Theory for Helium Atom

The goal of this project is to implement a **Self-Consistent Field (SCF)** solver to determine the ground state properties of the Helium atom using **Density Functional Theory (DFT)** by employing different approximation levels for the effective potential (Hartree only, Hartree + exchange, Hartree + exchange + correlation)

## Table Of Contents

1. [Dependencies](#dependencies)  
2. [Usage](#usage)  
3. [Configuration](#configuration)
4. [Outputs](#outputs)
5. [Background Theory](#background-theory)

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

## Output

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

When a system is studied using DFT, one introduces a fictitious system of non-interacting electrons that reproduces the exact interacting density. This system is governed by the Kohn–Sham (KS) equation which reads:

$$
\left[-\frac{1}{2}\nabla^2 + V_\mathrm{eff}(\mathbf{r})\right]\phi_i(\mathbf{r}) = \varepsilon_i \phi_i(\mathbf{r}),
$$

with the density:

$$
n(\mathbf{r}) = \sum_i f_i|\phi_i(\mathbf{r})|^2,
$$

where $f_i$ are occupations.

The effective potential is decomposed as:

$$
V_{\mathrm{eff}}(\mathbf{r}) = V_{\mathrm{ext}}(\mathbf{r}) + V_{H}[n] (\mathbf{r}) + V_{xc}[n] (\mathbf{r}),
$$

where:

- $V_\mathrm{ext}$ is the nuclear attraction (for helium: $V_\mathrm{ext}(r) = -Z/r$),
- $V_H$ is the Hartree (classical Coulomb) potential,
- $V_{xc}$ is the exchange-correlation potential (all many-body effects beyond Hartree).

### Hartree potential

The Hartree potential is:

$$
V_H(\mathbf{r}) = \int \frac{n(\mathbf{r}')}{|\mathbf{r} - \mathbf{r}'|}\,d^3r'.
$$

It can also be obtained by solving Poisson’s equation:

$$
\nabla^2 V_H(\mathbf{r}) = -4\pi n(\mathbf{r})
$$

<!-- with boundary condition $V_H(r\to\infty)\to 0$ (for finite systems). -->

### Exchange-correlation

The exchange-correlation energy is defined by:

$$
E_{xc}[n] = E[n] - T_s[n] - \int V_\mathrm{ext}(\mathbf{r})n(\mathbf{r})\,d^3r - E_H[n],
$$

where $T_s[n]$ is the kinetic energy of the non-interacting KS system and the Hartree energy is:

$$
E_H[n] = \frac{1}{2}\iint \frac{n(\mathbf{r})n(\mathbf{r}')}{|\mathbf{r}-\mathbf{r}'|}\,d^3r\,d^3r'.
$$

Given an approximation for $E_{xc}[n]$, the KS exchange-correlation potential is obtained by functional derivative:

$$
V_{xc}(\mathbf{r}) = \frac{\delta E_{xc}[n]}{\delta n(\mathbf{r})}.
$$

### Local Density Approximation (LDA) for exchange and correlation

The exact form of the exchange-correlation potential $V_{xc}(\mathbf{r})$ is unknown, it must be expressed using models and approximations. In general, for a non-homogeneous system, it is a function of the density and its gradients

$$V_{xc}(\mathbf{r}) = V_{xc}(n, \nabla n, \nabla^2 n, \ldots)$$

The Local Density Approximation (LDA) assumes that, at each point in space $\mathbf{r}$, the XC potential is well approximated by that of a uniform electron gas (UEG) evaluated at the local density $n(\mathbf{r})$.

The XC potential can be separated into exchange and correlation contributions:

$$
V_{xc}(r) = V_x(r) + V_c(r).
$$

### LDA exchange

Within LDA, the exchange potential is taken from the UEG result:

$$
V_x(r) = -\left(\frac{3}{\pi}\right)^{1/3} n^{1/3}(r).$$

For the helium atom, the density entering the Hartree and XC terms is the **full density** (two electrons in the 1s orbital). Using the radial representation in terms of the reduced radial function $u(r)$, this can be rewritten explicitly as a function of $u(r)$ and $r$:

$$
V_x(r) = -\left[\frac{3u^2(r)}{2\pi^2 r^2}\right]^{1/3}.
$$

### LDA correlation (Ceperley–Alder / Perdew–Zunger parameterization)

Correlation is included through a standard UEG-based parameterization derived from Quantum Monte Carlo data (Ceperley–Alder) and fit by Perdew–Zunger. Define the Wigner–Seitz radius $r_s$ by:

$$
n = \frac{3}{4\pi r_s^3}.
$$

The correlation potential is expressed in terms of the correlation energy per particle $\varepsilon_c(r_s)$ as:

$$
V_c(r_s) = \left(1-\frac{r_s}{3}\frac{d}{dr_s}\right)\varepsilon_c(r_s).
$$

where $\varepsilon_c$ is the correlation energy parameter defined by

$$
E_c = \int d^3r\, \varepsilon_c [n(\mathbf{r})] (\mathbf{r}),
$$

a commonly used parameterization for $\varepsilon_c$ is:

- For $r_s \ge 1$:

$$
\varepsilon_c = \frac{\gamma}{1+\beta_1\sqrt{r_s}+\beta_2 r_s},
$$

which leads to:

$$
V_c(r_s)= \varepsilon_c
\frac{1+\frac{7}{6}\beta_1\sqrt{r_s}+\frac{4}{3}\beta_2 r_s}{1+\beta_1\sqrt{r_s}+\beta_2 r_s}.
$$

- For $r_s < 1$:
  
$$
\varepsilon_c = A\ln r_s + B + C r_s \ln r_s + D r_s,
$$

which leads to:

$$
V_c(r_s)=A\ln r_s + B - \frac{A}{3} + \frac{2}{3}C r_s \ln r_s + \frac{(2D-C)}{3}r_s.
$$

where $\beta_1$, $\beta_2$, $\gamma$, A, B, C, D are tabulated parameters.

### Total Energy calculation

The LDA $V_x(r)$ and $V_c(r)$ enter the self-consistent Kohn–Sham effective potential

$$
V_\mathrm{eff}(r)=V_\mathrm{ext}(r)+V_H(r)+V_x(r)+V_c(r).
$$

Once the KS equations are solved for each electron, the total energy of the system is obtained by summing up the single electron eigenvalues $\varepsilon_i$ and removing double counted interactions.

For the Helium atom, for which $\sum_i\varepsilon_i=2\varepsilon$, the energy of the system is given by

$$
E = 2\varepsilon - \int dr\, V_H(r)\,u^2(r) - \frac{1}{2}\int dr\, u^2(r)\,V_x(r) +\int dr\, 2u^2(r)(\varepsilon_c(r)-V_c(r))
$$

