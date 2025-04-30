import ctypes
import numpy as np
import matplotlib.pyplot as plt

from numpy.ctypeslib import ndpointer

# Load the DLL
#dll_path = '.\\Resources\\Library\\spring_solver.dll'  # Change path as needed
dll_path = './Resources/Library/solve_spring_mass_v22.so'

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

def run_and_plot(P1, P2, n, c, total_mass, g):
    P1 = np.array(P1, dtype=np.float64)
    P2 = np.array(P2, dtype=np.float64)
    out_positions = np.zeros((n - 1) * 3, dtype=np.float64)

    g_vec = np.array([0.0, 0,-g], dtype=np.float64)

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

    ax2 = fig.add_subplot(122)
    ax2.plot(spring_lengths, marker='o')
    ax2.set_title("Spring Lengths")
    ax2.set_xlabel("Spring Index")
    ax2.set_ylabel("Length [m]")
    plt.tight_layout()
    plt.show()

# Example test cases
run_and_plot(P1=[0, 0, 0], P2=[0, 1.9, 10], n=10, c=100, total_mass=20, g=9.81)
run_and_plot(P1=[0, 0, 0], P2=[0, 1.8, 10], n=10, c=100, total_mass=20, g=9.81)
run_and_plot(P1=[0, 0, 0], P2=[0, 1.7, 10], n=10, c=10, total_mass=20, g=9.81)
run_and_plot(P1=[0, 0, 0], P2=[0, 1.6, 10], n=10, c=100, total_mass=20, g=9.81)
run_and_plot(P1=[0, 0, 0], P2=[0, 1.5, 10], n=11, c=100, total_mass=20, g=9.81)
run_and_plot(P1=[0, 0, 0], P2=[0, 1.4, 10], n=10, c=100, total_mass=20, g=9.81)
run_and_plot(P1=[0, 0, 0], P2=[0, 1.3, 10], n=10, c=100, total_mass=20, g=9.81)
run_and_plot(P1=[0, 0, 0], P2=[0, 1.2, 10], n=10, c=100, total_mass=20, g=9.81)
run_and_plot(P1=[0, 0, 0], P2=[0, 1.1, 10], n=10, c=100, total_mass=20, g=9.81)
run_and_plot(P1=[0, 0, 0], P2=[0, 1.0, 10], n=10, c=100, total_mass=20, g=9.81)
run_and_plot(P1=[0, 0, 0], P2=[0, 0.9, 10], n=10, c=100, total_mass=20, g=9.81)
run_and_plot(P1=[0, 0, 0], P2=[0, 0.8, 10], n=10, c=100, total_mass=20, g=9.81)
run_and_plot(P1=[0, 0, 0], P2=[0, 0.7, 10], n=10, c=100, total_mass=20, g=9.81)
run_and_plot(P1=[0, 0, 0], P2=[0, 0.6, 10], n=10, c=100, total_mass=20, g=9.81)
run_and_plot(P1=[0, 0, 0], P2=[0, 0.5, 10], n=11, c=100, total_mass=20, g=9.81)
run_and_plot(P1=[0, 0, 0], P2=[0, 0.4, 10], n=10, c=100, total_mass=20, g=9.81)
run_and_plot(P1=[0, 0, 0], P2=[0, 0.3, 10], n=10, c=100, total_mass=20, g=9.81)
run_and_plot(P1=[0, 0, 0], P2=[0, 0.2, 10], n=10, c=100, total_mass=20, g=9.81)
run_and_plot(P1=[0, 0, 0], P2=[0, 0.1, 10], n=10, c=100, total_mass=20, g=9.81)