import numpy as np
from scipy.optimize import root_scalar

def generate_catenary_nodes(P1, P2, g_vec, length_factor=1.5, n=10):
    P1 = np.array(P1, dtype=float)
    P2 = np.array(P2, dtype=float)
    g_vec = np.array(g_vec, dtype=float)
    
    # Normalize gravity and compute local frame
    L = P2 - P1
    d = np.linalg.norm(L)
    s = length_factor * d
    
    e_x = L / d
    g_proj = g_vec - np.dot(g_vec, e_x) * e_x
    g_proj_norm = np.linalg.norm(g_proj)
    if g_proj_norm < 1e-8:
        raise ValueError("Gravity vector is nearly parallel to the rope.")
    e_y = g_proj / g_proj_norm
    e_z = np.cross(e_x, e_y)
    R = np.vstack((e_x, e_y, e_z)).T

    # Coordinates in local frame
    dx = d
    dy = np.dot(P2 - P1, e_y)
    dz = np.dot(P2 - P1, e_z)

    # Horizontal 2D catenary in x-y plane
    def arc_length_eq(a):
        return 2 * a * np.sinh(dx / (2 * a)) - s

    sol = root_scalar(arc_length_eq, bracket=[1e-6, 1e6], method='brentq')
    if not sol.converged:
        raise RuntimeError("Catenary solver did not converge")
    a = sol.root

    x_vals = np.linspace(0, dx, n)
    x_mid = dx / 2
    y0 = dy / 2 - a * np.cosh((x_mid - x_mid) / a) + a * np.cosh((x_vals - x_mid) / a)
    y_vals = y0

    # z component interpolated linearly
    z_vals = np.linspace(0, dz, n)

    # Compose local positions
    local_pts = np.vstack((x_vals, y_vals, z_vals)).T
    world_pts = (R @ local_pts.T).T + P1

    return world_pts
