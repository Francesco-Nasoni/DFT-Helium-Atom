import numpy as np
from scipy.integrate import simpson
from source.solver import verlet_integrate_1D


def get_V_h(u_r, line_grid):
    h = np.abs(line_grid[1] - line_grid[0])

    def U_2_poisson(r, _, idx):
        r_eff = max(r, 1e-12)
        return -(u_r[idx] ** 2) / r_eff

    # We solve the particular solution by integrating using Verlet and boundary
    # conditions U(0)=0 and U(h)=0 so that both U and U' are zero in the origin
    U_part = verlet_integrate_1D(U_2_poisson, 0, 0, line_grid)
    qmax = simpson(u_r**2, line_grid)
    alpha = (qmax - U_part[-1]) / line_grid[-1]

    V_hartree = (U_part / line_grid) + alpha
    V_hartree[0] = V_hartree[1]
    return V_hartree


def get_V_x(u_r, line_grid):
    vx = -((3 * (u_r**2) / 2 / np.pi / np.pi / line_grid / line_grid) ** (1 / 3))
    return vx


def get_V_c(u_r, line_grid):

    # PZ81 coefficients
    gamma = -0.1423
    beta1 = 1.0529
    beta2 = 0.3334
    A = 0.0311
    B = -0.048
    C = 0.0020
    D = -0.0116

    r_s = (3 * line_grid**2 / (2 * u_r**2)) ** (1 / 3)
    mask_rs_low = r_s > 1

    sqrt_rs = np.sqrt(r_s)
    log_rs = np.log(r_s)

    # High density
    ec_low = gamma / (1 + beta1 * sqrt_rs + beta2 * r_s)
    numerator = 1 + (7 / 6) * beta1 * sqrt_rs + (4 / 3) * beta2 * r_s
    denominator = 1 + beta1 * sqrt_rs + beta2 * r_s
    vc_low = ec_low * (numerator / denominator)

    # Low density
    ec_high = A * log_rs + B + C * r_s * log_rs + D * r_s
    vc_high = A * log_rs + B - A / 3 + (2 / 3) * C * r_s * log_rs + (2 * D - C) * r_s / 3

    v_c = np.where(mask_rs_low, vc_low, vc_high)
    ec_r = np.where(mask_rs_low, ec_low, ec_high)

    return v_c, ec_r


def get_V_eff(u_r, r, use_exchange, use_correlation):

    V_eff = np.zeros(r.shape)
    V_X = None
    V_C = None
    ec = None
    if use_exchange or use_correlation:
        V_H = 2 * get_V_h(u_r, r)
        if use_exchange:
            V_X = get_V_x(u_r, r)
            V_eff += V_X

        if use_correlation:
            V_C, ec = get_V_c(u_r, r)
            V_eff += V_C
    else:
        V_H = get_V_h(u_r, r)

    V_eff += V_H

    return V_eff, V_H, V_X, V_C, ec


def get_TOTEN(E, u_r, line_grid, v_h, v_x=None, v_c=None, ec_c=None):
    toten = 2 * E - simpson(v_h * u_r**2, line_grid)
    if v_x is not None:
        toten += -0.5 * simpson(v_x * u_r**2, line_grid)
    if v_c is not None:
        toten += simpson((ec_c - v_c) * 2 * u_r**2, line_grid)

    return toten
