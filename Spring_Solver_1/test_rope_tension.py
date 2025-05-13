import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer

# Load the DLL
lib = ctypes.CDLL('./Resources/Library/solve_rope_tension_lapak.dll')

# Define function signature (must match the .h/.c exactly!)
lib.solve_rope_tension.argtypes = [
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # P1
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # P2
    ctypes.c_int,                                      # n
    ctypes.c_double,                                   # total_mass
    ctypes.c_double,                                   # rope_diameter
    ctypes.c_double,                                   # youngs_modulus
    ctypes.c_double,                                   # F_target
    ctypes.POINTER(ctypes.c_double),                   # L_out    
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # g_vec
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # out_positions
    ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),  # F_P1_n
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # F_P2_n
    ctypes.POINTER(ctypes.c_int),                      # Status_dynamic
    ctypes.POINTER(ctypes.c_int)                       # Status_newton    
]
lib.solve_rope_tension.restype = ctypes.c_int  # Returns SolveStatus

# -------------------------------
# Example Input
# -------------------------------

P1 = np.array([0.0, 0.0, 0.0], dtype=np.float64)
P2 = np.array([10.0, 0.0, 0.0], dtype=np.float64)
tension_mag = 30.0  # Newtons

n = 10
total_mass = 1.0
rope_diameter = 0.02
youngs_modulus = 1e9
g_vec = np.array([0.0, 0.0, -9.81], dtype=np.float64)

out_positions = np.zeros((n - 1) * 3, dtype=np.float64)
F_P1_n = np.zeros(3, dtype=np.float64)
F_P2_n = np.zeros(3, dtype=np.float64)
L_out = ctypes.c_double()
status_dynamic = ctypes.c_int()
status_newton = ctypes.c_int()

# -------------------------------
# Call the Solver
# -------------------------------

result = lib.solve_rope_tension(
    P1, P2,
    n,
    total_mass,
    rope_diameter,
    youngs_modulus,
    tension_mag,
    ctypes.byref(L_out),
    g_vec,
    out_positions,
    F_P1_n,
    F_P2_n,
    ctypes.byref(status_dynamic),
    ctypes.byref(status_newton),
)

# -------------------------------
# Output
# -------------------------------

print("Solver return code:", result)
print("Status (dynamic):", status_dynamic.value)
print("Status (newton):", status_newton.value)
print("Rope length:", L_out.value)
print("Force on P1 (Newton-based):", F_P1_n)
print("Force on P2 (Newton-based):", F_P2_n)
print("Node positions:")
print(out_positions.reshape(-1, 3))
