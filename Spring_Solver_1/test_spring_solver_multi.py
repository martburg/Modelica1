import ctypes
import numpy as np
import matplotlib.pyplot as plt

from numpy.ctypeslib import ndpointer

# Load the DLL
dll_path = '.\\Resources\\Library\\spring_solver.dll'  # Change path as needed
#dll_path = '.\\Resources\\Library\\spring_solver_ref_ref.dll'  # Change path as needed
#dll_path = './Resources/Library/solve_spring_mass_v22.so'

lib = ctypes.CDLL(dll_path)


# Setup the prototype
lib.solve_spring_mass_c.argtypes = [
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # P1
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # P2
    ctypes.c_int,                                       # n
    ctypes.c_double,                                    # c
    ctypes.c_double,                                    # total_mass
    ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),   # g
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")    # out_positions
]
lib.solve_spring_mass_c.restype = ctypes.c_int

def set_axes_equal(ax):
    """Set 3D plot axes to equal scale."""
    limits = np.array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        ax.get_zlim3d(),
    ])
    spans = limits[:, 1] - limits[:, 0]
    centers = np.mean(limits, axis=1)
    max_span = max(spans)

    new_limits = np.array([
        centers - max_span / 2,
        centers + max_span / 2
    ]).T

    ax.set_xlim3d(new_limits[0])
    ax.set_ylim3d(new_limits[1])
    ax.set_zlim3d(new_limits[2])

def run_and_plot(P1, P2, n, c, total_mass, g):
    P1 = np.array(P1, dtype=np.float64)
    P2 = np.array(P2, dtype=np.float64)
    out_positions = np.zeros((n - 1) * 3, dtype=np.float64)

    g_vec = np.array([0.0, 0., -g], dtype=np.float64)

    print(f"\nRunning solver with n={n}, c={c}, mass={total_mass}")
    res = lib.solve_spring_mass_c(P1, P2, n, c, total_mass, g_vec, out_positions)
    if res != 0:
        print("❌ Solver failed")
        return

    positions = out_positions.reshape((n - 1, 3))
    all_points = np.vstack([P1, positions, P2])
    spring_lengths = np.linalg.norm(np.diff(all_points, axis=0), axis=1)

    fig = plt.figure(figsize=(12, 6))
    ax1 = fig.add_subplot(121, projection='3d')
    ax1.plot(all_points[:, 0], all_points[:, 1], all_points[:, 2], '-o', label='Springs')
    ax1.scatter(P1[0], P1[1], P1[2], color='black', label='Fixed P1')
    ax1.scatter(P2[0], P2[1], P2[2], color='black', label='Fixed P2')
    ax1.set_title("3D View of Spring-Mass System")
    ax1.set_xlabel("X")
    ax1.set_ylabel("Y")
    ax1.set_zlabel("Z")
    ax1.legend()
    ax1.view_init(elev=20, azim=-60)
    set_axes_equal(ax1)  # <-- Enforce equal scaling for all axes

    ax2 = fig.add_subplot(122)
    ax2.plot(spring_lengths, marker='o')
    ax2.set_title("Spring Lengths")
    ax2.set_xlabel("Spring Index")
    ax2.set_ylabel("Length [m]")
    plt.tight_layout()
    plt.show()

# Example test cases
#run_and_plot(P1=[0, 0, 0], P2=[0, 100, 10], n=20, c=39250, total_mass=0.9, g=9.81)
#run_and_plot(P1=[0, 0, 0], P2=[0, 100, 10], n=20, c=15700, total_mass=1.1, g=9.81)
#run_and_plot(P1=[0, 0, 0], P2=[0, 100, 10], n=20, c=188400, total_mass=0.7, g=9.81)
#run_and_plot(P1=[0, 0, 0], P2=[0, 100, 10], n=20, c=2355000, total_mass=0.76, g=9.81)
#run_and_plot(P1=[0, 0, 0], P2=[0, 100, 10], n=20, c=12560000, total_mass=1.0, g=9.81)
#run_and_plot(P1=[0, 0, 0], P2=[0, 100, 10], n=20, c=31400000, total_mass=6.12, g=9.81)
run_and_plot(P1=[0, 0, 0], P2=[0, 1,10], n=20, c=39250, total_mass=0.09, g=9.81)
run_and_plot(P1=[0, 0, 0], P2=[0, 1, 10], n=20, c=15700, total_mass=0.11, g=9.81)
run_and_plot(P1=[0, 0, 0], P2=[0, 1, 10], n=20, c=188400, total_mass=0.07, g=9.81)
run_and_plot(P1=[0, 0, 0], P2=[0, 1, 10], n=20, c=2355000, total_mass=0.076, g=9.81)
run_and_plot(P1=[0, 0, 0], P2=[0, 1, 10], n=20, c=12560000, total_mass=0.10, g=9.81)
run_and_plot(P1=[0, 0, 0], P2=[0, 1, 10], n=20, c=31400000, total_mass=0.612, g=9.81)
#run_and_plot(P1=[0, 0, 0], P2=[0, 7, 0], n=10, c=100, total_mass=20, g=9.81)
#run_and_plot(P1=[0, 0, 0], P2=[0, 6, 0], n=10, c=100, total_mass=20, g=9.81)
#run_and_plot(P1=[0, 0, 0], P2=[0, 5, 0], n=10, c=100, total_mass=20, g=9.81)
#run_and_plot(P1=[0, 0, 0], P2=[0, 4, 0], n=10, c=100, total_mass=20, g=9.81)
#run_and_plot(P1=[0, 0, 0], P2=[0, 3, 0], n=10, c=100, total_mass=20, g=9.81)
#run_and_plot(P1=[0, 0, 0], P2=[0, 2, 0], n=10, c=100, total_mass=20, g=9.81)
#run_and_plot(P1=[0, 0, 0], P2=[0, 1.7, 0], n=10, c=100, total_mass=20, g=9.81)
#run_and_plot(P1=[0, 0, 0], P2=[0, 1.6, 0], n=10, c=100, total_mass=20, g=9.81)
#run_and_plot(P1=[0, 0, 0], P2=[0, 1.5, 0], n=10, c=100, total_mass=20, g=9.81)
#run_and_plot(P1=[0, 0, 0], P2=[0, 1.4, 0], n=10, c=100, total_mass=20, g=9.81)
#run_and_plot(P1=[0, 0, 0], P2=[0, 1.3, 10], n=10, c=100, total_mass=20, g=9.81)
#run_and_plot(P1=[0, 0, 0], P2=[0, 1.2, 0], n=10, c=100, total_mass=20, g=9.81)
#run_and_plot(P1=[0, 0, 0], P2=[0, 1.1, 0], n=10, c=100, total_mass=20, g=9.81)
#run_and_plot(P1=[0, 0, 0], P2=[0, 0.9, 0], n=10, c=100, total_mass=20, g=9.81)
#run_and_plot(P1=[0, 0, 0], P2=[0, 0.8, 10], n=10, c=100, total_mass=20, g=9.81)
#run_and_plot(P1=[0, 0, 0], P2=[0, 0.7, 0], n=10, c=100, total_mass=20, g=9.81)
#run_and_plot(P1=[0, 0, 0], P2=[0, 0.6, 10], n=10, c=100, total_mass=20, g=9.81)
#run_and_plot(P1=[0, 0, 0], P2=[0, 0.5, 0], n=10, c=100, total_mass=20, g=9.81)
#run_and_plot(P1=[0, 0, 0], P2=[0, 0.4, 0], n=10, c=1000, total_mass=20, g=9.81)
#run_and_plot(P1=[0, 0, 0], P2=[0, 0.3, 0], n=10, c=10000, total_mass=20, g=9.81)
#run_and_plot(P1=[0, 0, 0], P2=[0, 0.2, 0], n=10, c=100000, total_mass=20, g=9.81)
run_and_plot(P1=[0, 0, 0], P2=[0, 0.1, 10], n=20, c=39250, total_mass=0.9, g=9.81)
