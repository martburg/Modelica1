import numpy as np
from scipy.optimize import root
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def compute_spring_force(p1, p2, s0, c):
    delta = p2 - p1
    length = np.linalg.norm(delta)
    force = - (length - s0) * c * delta / length  # Spring force
    return force, length  # Return the spring force and length for debugging

def equilibrium_residuals(positions_flat, P1, P2, n, s0_list, total_mass, c_real, g_scaled):
    positions = positions_flat.reshape((n - 1, 3))
    all_points = [P1] + list(positions) + [P2]

    residuals = []
    m = total_mass / (n - 1)

    for i in range(1, n):
        left = all_points[i - 1]
        mid = all_points[i]
        right = all_points[i + 1]

        F_left, length_left = compute_spring_force(left, mid, s0_list[i - 1], c_real)
        F_right, length_right = compute_spring_force(right, mid, s0_list[i], c_real)
        F_gravity = np.array([0, 0, -m * g_scaled])
        net_force = F_left + F_right + F_gravity

        residuals.append(net_force)

    return np.concatenate(residuals)

def refined_initial_guess(P1, P2, n, max_sag=0.5):
    # More gradual initial guess: cubic interpolation between P1 and P2
    P1 = np.array(P1)
    P2 = np.array(P2)
    guesses = []
    for i in range(1, n):
        t = i / (n - 1)
        point = (1 - t) * P1 + t * P2
        # Gradually increase the sag with cubic function, more natural curve
        sag = -max_sag * (1 - 4 * t * (1 - t))  # Cubic curve sag
        point[2] += sag
        guesses.append(point)
    return np.array(guesses)

def solve_mass_positions(P1, P2, n, s0_list, c_real, total_mass_real, g=9.81, max_sag=0.5):
    P1_norm = np.array(P1)
    P2_norm = np.array(P2)
    s0_norm = s0_list
    g_scaled = g  # Keep real gravity, no scaling


    # Better initial guess with gradual sag
    guesses = refined_initial_guess(P1_norm, P2_norm, n, max_sag)
    positions_flat = guesses.flatten()

    result = root(
        equilibrium_residuals,
        positions_flat,
        args=(P1_norm, P2_norm, n, s0_norm, total_mass_real, c_real, g_scaled),
        method='hybr',
        options={'maxfev': 100000}
    )

    if not result.success:
        raise RuntimeError("Solver failed: " + result.message)

    positions_norm = result.x.reshape((n - 1, 3))
    return positions_norm #* L   for i in range(1, n):
 


if __name__ == "__main__":
    P1 = np.array([15, -7, 0])
    P2 = np.array([-15, 70, 0])
    n = 50
    total_length = np.linalg.norm(P2 - P1)
    s0_list = [total_length / n] * n
    c = 100
    total_mass = 20

    positions = solve_mass_positions(P1, P2, n, s0_list, c, total_mass, max_sag=4.5)
