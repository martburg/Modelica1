import numpy as np
from scipy.optimize import root

def compute_spring_force(p1, p2, s0, c):
    delta = p2 - p1
    length = np.linalg.norm(delta)
    force = - (length - s0) * c * delta / length
    return force, length

def compute_spring_lengths(positions, P1, P2):
    all_points = [P1] + list(positions) + [P2]
    return [np.linalg.norm(all_points[i + 1] - all_points[i]) for i in range(len(all_points) - 1)]

def equilibrium_residuals(positions_flat, P1, P2, n, s0_list, total_mass, c_real, g_scaled):
    positions = positions_flat.reshape((n - 1, 3))
    all_points = [P1] + list(positions) + [P2]
    residuals = []
    m = total_mass / (n - 1)

    for i in range(1, n):
        left, mid, right = all_points[i - 1], all_points[i], all_points[i + 1]
        F_left, _ = compute_spring_force(left, mid, s0_list[i - 1], c_real)
        F_right, _ = compute_spring_force(right, mid, s0_list[i], c_real)
        F_gravity = np.array([0, 0, -m * g_scaled])
        residuals.append(F_left + F_right + F_gravity)

    return np.concatenate(residuals)

def refined_initial_guess(P1, P2, n, max_sag=0.5):
    P1 = np.array(P1)
    P2 = np.array(P2)
    guesses = []
    for i in range(1, n):
        t = i / (n - 1)
        point = (1 - t) * P1 + t * P2
        sag = -max_sag * (1 - 4 * t * (1 - t))  # Cubic curve for sag
        point[2] += sag
        guesses.append(point)
    return np.array(guesses)

def solve_mass_positions(P1, P2, n, s0_list, c_real, total_mass_real, g=9.81, max_sag=0.5):
    P1 = np.array(P1, dtype=np.float64)
    P2 = np.array(P2, dtype=np.float64)
    s0_list = np.array(s0_list, dtype=np.float64)
    guesses = refined_initial_guess(P1, P2, n, max_sag)
    positions_flat = guesses.flatten()

    result = root(
        equilibrium_residuals,
        positions_flat,
        args=(P1, P2, n, s0_list, total_mass_real, c_real, g),
        method='hybr',
        options={'maxfev': 100000}
    )

    if not result.success:
        raise RuntimeError("Solver failed: " + result.message)

    return result.x.reshape((n - 1, 3))
