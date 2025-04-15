import numpy as np
from scipy.optimize import root
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def compute_spring_force(p1, p2, s0, c=1.0):
    delta = p2 - p1
    length = np.linalg.norm(delta)
    force = - (length - s0) * c * delta / length
    return force, length

def compute_spring_lengths(positions, P1, P2):
    all_points = [P1] + list(positions) + [P2]
    spring_lengths = []
    for i in range(len(all_points) - 1):
        length = np.linalg.norm(all_points[i + 1] - all_points[i])
        spring_lengths.append(length)
    return spring_lengths

def equilibrium_residuals(positions_flat, P1, P2, n, s0_list, total_mass, g_scaled):
    positions = positions_flat.reshape((n - 1, 3))
    all_points = [P1] + list(positions) + [P2]

    residuals = []
    m = total_mass / (n - 1)

    for i in range(1, n):
        left = all_points[i - 1]
        mid = all_points[i]
        right = all_points[i + 1]

        F_left, _ = compute_spring_force(left, mid, s0_list[i - 1])
        F_right, _ = compute_spring_force(right, mid, s0_list[i])
        F_gravity = np.array([0, 0, -m * g_scaled])
        net_force = F_left + F_right + F_gravity
        residuals.append(net_force)

    return np.concatenate(residuals)

def refined_initial_guess(P1, P2, n, max_sag=0.5):
    P1 = np.array(P1)
    P2 = np.array(P2)
    guesses = []
    for i in range(1, n):
        t = i / (n - 1)
        point = (1 - t) * P1 + t * P2
        sag = -max_sag * (1 - 4 * t * (1 - t))
        point[2] += sag
        guesses.append(point)
    return np.array(guesses)

def solve_mass_positions(P1, P2, n, s0_list, c_real, total_mass_real, g=9.81, max_sag=0.5):
    L = np.linalg.norm(P2 - P1)
    P1_norm = np.array(P1) / L
    P2_norm = np.array(P2) / L
    s0_norm = [s / L for s in s0_list]
    g_scaled = g / (c_real * L)

    guesses = refined_initial_guess(P1_norm, P2_norm, n, max_sag)
    positions_flat = guesses.flatten()

    result = root(
        equilibrium_residuals,
        positions_flat,
        args=(P1_norm, P2_norm, n, s0_norm, total_mass_real, g_scaled),
        method='hybr',
        options={'maxfev': 10000}
    )

    if not result.success:
        raise RuntimeError("Solver failed: " + result.message)

    positions_norm = result.x.reshape((n - 1, 3))
    return positions_norm * L

def verify_equilibrium_real_units(P1, P2, positions, s0_list, c, total_mass, g=9.81, tolerance=1e-6):
    n = len(s0_list)
    m = total_mass / (n - 1)
    all_points = [np.array(P1)] + list(positions) + [np.array(P2)]
    passed = True

    print("\n=== Equilibrium Check Per Mass (Real Units) ===")
    for i in range(1, n):
        left = all_points[i - 1]
        mid = all_points[i]
        right = all_points[i + 1]

        F_left, _ = compute_spring_force(left, mid, s0_list[i - 1], c)
        F_right, _ = compute_spring_force(right, mid, s0_list[i], c)
        F_gravity = np.array([0, 0, -m * g])

        net_force = F_left + F_right + F_gravity
        norm = np.linalg.norm(net_force)

        print(f"Mass {i}: Net Force = {net_force} N, |F| = {norm:.3e} N")
        if norm > tolerance:
            passed = False

    if passed:
        print("\u2705 All masses are in force equilibrium (within tolerance).")
    else:
        print("\u274C Some masses are NOT in force equilibrium!")

def plot_spring_mass_system(P1, P2, mass_positions, spring_lengths):
    fig = plt.figure(figsize=(12, 6))
    ax1 = fig.add_subplot(121, projection='3d')
    all_points = [P1] + list(mass_positions) + [P2]
    all_points = np.array(all_points)

    for i in range(len(all_points) - 1):
        xs = [all_points[i][0], all_points[i+1][0]]
        ys = [all_points[i][1], all_points[i+1][1]]
        zs = [all_points[i][2], all_points[i+1][2]]
        ax1.plot(xs, ys, zs, 'b-', linewidth=2)

    ax1.scatter(mass_positions[:, 0], mass_positions[:, 1], mass_positions[:, 2], 
                color='red', s=50, label='Masses')
    ax1.scatter([P1[0], P2[0]], [P1[1], P2[1]], [P1[2], P2[2]], 
                color='black', s=50, label='Fixed Points')

    y_min, y_max = np.min(all_points[:,1]), np.max(all_points[:,1])
    if abs(y_max - y_min) < 1e-6:
        ax1.set_ylim(-0.5, 0.5)

    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Z (sag downward)')
    ax1.set_title('Mass-Spring Chain in 3D')
    ax1.legend()
    ax1.view_init(elev=20, azim=-60)

    ax2 = fig.add_subplot(122)
    ax2.plot(spring_lengths, marker='o', linestyle='-', color='b', label='Spring Lengths')
    ax2.set_xlabel('Spring Index')
    ax2.set_ylabel('Spring Length (m)')
    ax2.set_title('Spring Lengths in the System')
    ax2.legend()

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    P1 = np.array([-2, 0, 0])
    P2 = np.array([2, 0, 0])
    n = 4
    total_length = np.linalg.norm(P2 - P1)
    s0_list = [total_length / n] * n
    c = 100  # real spring constant in N/m
    total_mass = 20  # real total mass in kg

    positions = solve_mass_positions(P1, P2, n, s0_list, c, total_mass, max_sag=4.5)
    verify_equilibrium_real_units(P1, P2, positions, s0_list, c, total_mass)

    for i, pos in enumerate(positions, start=1):
        print(f"Mass {i}: x = {pos[0]:.3f}, y = {pos[1]:.3f}, z = {pos[2]:.3f}")

    spring_lengths = compute_spring_lengths(positions, P1, P2)
    print("Spring lengths:", spring_lengths)
    plot_spring_mass_system(P1, P2, positions, spring_lengths)