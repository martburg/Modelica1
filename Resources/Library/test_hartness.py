import ctypes
import numpy as np
import os

print("Start")

# Load your DLL
dll_path = os.path.abspath("spring_solver.dll")  # or .pyd if that's what you've got
print(dll_path)
solver = ctypes.CDLL(dll_path)

# Define function prototype
solver.solve_spring_mass_c.argtypes = [
    ctypes.POINTER(ctypes.c_double),  # P1
    ctypes.POINTER(ctypes.c_double),  # P2
    ctypes.c_int,                     # n
    ctypes.POINTER(ctypes.c_double),  # s0_list
    ctypes.c_double,                 # c
    ctypes.c_double,                 # total_mass
    ctypes.c_double,                 # g
    ctypes.c_double,                 # max_sag
    ctypes.POINTER(ctypes.c_double)  # out_positions
]
solver.solve_spring_mass_c.restype = ctypes.c_int

# Prepare data
n = 10
P1 = np.array([0.0, 0.0, 0.0], dtype=np.float64)
P2 = np.array([9.0, 0.0, 0.0], dtype=np.float64)
s0_list = np.ones(n, dtype=np.float64)  # Just dummy data
out_positions = np.zeros((n - 1, 3), dtype=np.float64)

res = solver.solve_spring_mass_c(
    P1.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    P2.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    n,
    s0_list.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    ctypes.c_double(100.0),
    ctypes.c_double(20.0),
    ctypes.c_double(9.81),
    ctypes.c_double(10.0),
    out_positions.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
)

print(f"Returned code: {res}")
print("Output positions:")
print(out_positions)
