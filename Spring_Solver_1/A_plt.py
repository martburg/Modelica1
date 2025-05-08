import numpy as np
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt

def generate_catenary_nodes(P1, P2, g_vec, length_factor=1.5, n=10):
    P1 = np.array(P1, dtype=float)
    P2 = np.array(P2, dtype=float)
    g_vec = np.array(g_vec, dtype=float)
    
    # Direction vector from P1 to P2
    L = P2 - P1
    dx = np.linalg.norm(L)
    s = length_factor * dx

    # Define local frame
    e_x = L / dx
    g_proj = g_vec - np.dot(g_vec, e_x) * e_x
    g_proj_norm = np.linalg.norm(g_proj)
    if g_proj_norm < 1e-8:
        raise ValueError("Gravity vector is nearly parallel to the rope.")
    e_y = g_proj / g_proj_norm
    e_z = np.cross(e_x, e_y)
    R = np.vstack((e_x, e_y, e_z)).T

    # Relative vector in local frame
    rel = P2 - P1
    dx_local = np.dot(rel, e_x)
    dy_local = np.dot(rel, e_y)
    dz_local = np.dot(rel, e_z)

    # === Dimensionless catenary solve ===
    s_scaled = s / dx_local

    def arc_length_eq_scaled(a_scaled):
        return 2 * a_scaled * np.sinh(1 / (2 * a_scaled)) - s_scaled

    sol = root_scalar(arc_length_eq_scaled, bracket=[1e-3, 1e3], method='brentq')
    if not sol.converged:
        raise RuntimeError("Catenary solver did not converge")

    a = sol.root * dx_local

    # Generate local coordinates
    x_vals = np.linspace(0, dx_local, n)
    x_mid = dx_local / 2
    y_vals = -a * (np.cosh((x_vals - x_mid) / a) - np.cosh(0 / a))
    y_vals += dy_local / dx_local * x_vals  # linear offset to match vertical drop
    z_vals = np.linspace(0, dz_local, n)

    local_pts = np.vstack((x_vals, y_vals, z_vals)).T
    world_pts = (R @ local_pts.T).T + P1
    return world_pts

# === Example usage and plot
if __name__ == "__main__":
    P1 = np.array([0.0, 0.0,10.0])
    P2 = np.array([10.0, 0.0, 0.0])
    g_vec = np.array([0.0, 0.0, -9.81])
    length_factor = 1.5
    n = 10

    rope_xyz = generate_catenary_nodes(P1, P2, g_vec, length_factor, n)

    print("Rope node coordinates:")
    print(rope_xyz)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(rope_xyz[:, 0], rope_xyz[:, 1], rope_xyz[:, 2], '-o', label='Catenary Rope')
    ax.quiver(P1[0], P1[1], P1[2], g_vec[0], g_vec[1], g_vec[2], length=10, color='r', label='Gravity')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Stable 3D Catenary Curve')
    ax.legend()
    plt.show()
