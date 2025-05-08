import ctypes
import numpy as np
from numpy.linalg import norm
from numpy.ctypeslib import ndpointer

# === Load shared library ===
dll_path = '.\\Resources\\Library\\spring_solver.dll'  # Update path if needed
lib = ctypes.CDLL(dll_path)

# === Define C function signature ===
lib.solve_spring_mass_c.argtypes = [
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # P1
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # P2
    ctypes.c_int,                                       # n
    ctypes.c_double,                                    # total_mass
    ctypes.c_double,                                    # length_factor
    ctypes.c_double,                                    # rope_diameter
    ctypes.c_double,                                    # youngs_modulus
    ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),   # g
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")    # out_positions
]
lib.solve_spring_mass_c.restype = ctypes.c_int

# === Solver function to compute required length factor ===
def solve_required_length_factor(P1, P2, total_mass, rope_diameter, youngs_modulus,
                                  g_vec, desired_sag, n=20, tol=1e-4, max_iter=50):
    P1 = np.ascontiguousarray(P1, dtype=np.float64)
    P2 = np.ascontiguousarray(P2, dtype=np.float64)
    g_vec = np.ascontiguousarray(g_vec, dtype=np.float64)
    assert P1.shape == (3,) and P2.shape == (3,) and g_vec.shape == (3,)

    dof = (n - 1) * 3
    out_positions = np.zeros(dof, dtype=np.float64)

    # Bisection search
    lf_low = 1.0
    lf_high = 1.05

    for iteration in range(max_iter):
        lf_mid = 0.5 * (lf_low + lf_high)

        ret = lib.solve_spring_mass_c(
            P1, P2,
            ctypes.c_int(n),
            ctypes.c_double(total_mass),
            ctypes.c_double(lf_mid),
            ctypes.c_double(rope_diameter),
            ctypes.c_double(youngs_modulus),
            g_vec,
            out_positions
        )

        if ret != 0:
            raise RuntimeError(f"Solver failed at iteration {iteration}, length_factor={lf_mid}")

        midpoint = out_positions[(n//2 - 1)*3:(n//2)*3]
        line_mid = 0.5 * (P1 + P2)
        sag = line_mid[2] - midpoint[2]

        if abs(sag - desired_sag) < tol:
            return lf_mid, out_positions.copy()

        if sag < desired_sag:
            lf_low = lf_mid
        else:
            lf_high = lf_mid

    raise RuntimeError(f"Failed to converge to desired sag {desired_sag} after {max_iter} iterations")

# === Main usage example ===
if __name__ == "__main__":
    P1 = np.array([0, 0, 0], dtype=np.float64)
    P2 = np.array([0, 100, 0], dtype=np.float64)
    total_mass = 1.0
    rope_diameter = 0.01
    youngs_modulus = 1e10
    g_vec = np.array([0, 0, -9.81], dtype=np.float64)
    desired_sag = 4.0
    n = 20

    print(f"Finding required rope length factor for sag = {desired_sag} meters...")

    length_factor, positions = solve_required_length_factor(
        P1, P2, total_mass, rope_diameter, youngs_modulus, g_vec, desired_sag, n
    )

    print(f"\n✅ Required length factor: {length_factor:.6f}\n")
    print("=== Rope Node Positions ===")
    for i in range(n - 1):
        x, y, z = positions[i*3:i*3+3]
        print(f"Node {i + 1:2d}: x = {x: .3f}, y = {y: .3f}, z = {z: .3f}")
