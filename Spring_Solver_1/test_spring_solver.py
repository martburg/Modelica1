import ctypes
import numpy as np

# Load the compiled DLL
dll_path = ".\\Spring_Solver_1\\spring_solver.dll"
lib = ctypes.CDLL(dll_path)

# Define the argument and return types for the C function
lib.solve_spring_mass_c.argtypes = [
    ctypes.POINTER(ctypes.c_double),  # P1
    ctypes.POINTER(ctypes.c_double),  # P2
    ctypes.c_int,                     # n
    ctypes.c_double,                  # c_real
    ctypes.c_double,                  # total_mass_real
    ctypes.c_double,                  # g
    ctypes.c_double,                  # max_sag
    ctypes.POINTER(ctypes.c_double)   # out_positions
]
lib.solve_spring_mass_c.restype = ctypes.c_int

def test_spring_solver():
    n = 5  # number of masses
    P1 = np.array([-10.0, 0.0, 0.0], dtype=np.float64)
    P2 = np.array([10.0, 0.0, 0.0], dtype=np.float64)
    c = 100.0
    total_mass = 50.0
    g = 9.81
    max_sag = 1.5

    # Allocate output buffer (n-1 masses × 3D)
    out_positions = np.zeros((n - 1) * 3, dtype=np.float64)

    result = lib.solve_spring_mass_c(
        P1.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        P2.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_int(n),
        ctypes.c_double(c),
        ctypes.c_double(total_mass),
        ctypes.c_double(g),
        ctypes.c_double(max_sag),
        out_positions.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    )

    if result != 0:
        print(f"Solver failed with code: {result}")
    else:
        print("Solver succeeded.")
        positions = out_positions.reshape((n - 1, 3))
        for i, pos in enumerate(positions):
            print(f"Mass {i+1}: x={pos[0]:.3f}, y={pos[1]:.3f}, z={pos[2]:.3f}")

if __name__ == "__main__":
    test_spring_solver()
