import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer

np.random.seed(42)  # Fixed seed for reproducibility

# Load DLL
dll_path = './Resources/Library/spring_solver_lapak.dll'
lib = ctypes.CDLL(dll_path)

lib.solve_spring_mass_c.argtypes = [
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # P1
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # P2
    ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),  # g_vec
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # out_positions
    ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),  # F_P1_n
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # F_P2_n
    ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),  # F_P1_w
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # F_P2_w
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
]
lib.solve_spring_mass_c.restype = ctypes.c_int

def force_test_case():
    P1 = [0,0,0]
    P2 = [0,100,0]
    P1 = np.array(P1)
    P2 = np.array(P2)
    g_vec = [-9.81,0,0]


    return {
        "P1": P1, "P2": P2,
        "n": 20, "total_mass": 1,
        "length_factor": 1.5,
        "rope_diameter": 0.01,
        "youngs_modulus": 1e9,
        "g_vec": g_vec
    }

def generate_test_case():
    # Random 3D endpoints within ±100 m cube
    P1 = np.random.uniform(-100, 100, 3).tolist()
    P2 = np.random.uniform(-100, 100, 3).tolist()

    # Number of mass points
    n = np.random.randint(10, 100)

    # Length factor: rope is longer than straight-line distance
    length_factor = np.random.uniform(1.1, 2.5)

    # Rope diameter: 1–8 cm
    rope_diameter = np.random.uniform(0.01, 0.08)

    # Density: 500–8000 kg/m³ (nylon to steel range)
    density = np.random.uniform(500, 8000)

    # Young’s modulus: 1e9–2e11 Pa
    youngs_modulus = np.random.uniform(1e9, 2e11)

    # Gravity vector
    g_vec = np.random.uniform(-15, 15, 3).tolist()

    # Compute length and mass from geometry
    P1_np = np.array(P1)
    P2_np = np.array(P2)
    L_straight = np.linalg.norm(P2_np - P1_np)
    L0 = L_straight * length_factor
    area = np.pi * (rope_diameter ** 2) / 4.0
    volume = area * L0
    total_mass = density * volume

    return {
        "P1": P1, "P2": P2,
        "n": n, "total_mass": total_mass,
        "length_factor": length_factor,
        "rope_diameter": rope_diameter,
        "youngs_modulus": youngs_modulus,
        "g_vec": g_vec
    }

def run_and_test(P1, P2, n, total_mass, length_factor, rope_diameter,
                 youngs_modulus, g_vec, tol_percent=2.5, force_tol=0.02, verbose=True):
    P1 = np.array(P1, dtype=np.float64)
    P2 = np.array(P2, dtype=np.float64)
    g_vec = np.array(g_vec, dtype=np.float64)

    out_positions = np.zeros((n - 1) * 3, dtype=np.float64)
    F_P1_n = np.zeros(3, dtype=np.float64)
    F_P2_n = np.zeros(3, dtype=np.float64)
    F_P1_w = np.zeros(3, dtype=np.float64)
    F_P2_w = np.zeros(3, dtype=np.float64)

    Length_init    = ctypes.c_double()
    Length_dyn     = ctypes.c_double()
    Length_newton  = ctypes.c_double()
    Status_dyn     = ctypes.c_int()
    Status_newton  = ctypes.c_int()

    res = lib.solve_spring_mass_c(
        P1, P2, n, total_mass, length_factor, rope_diameter, youngs_modulus,
        g_vec, out_positions, F_P1_n, F_P2_n, F_P1_w, F_P2_w,
        ctypes.byref(Length_init), ctypes.byref(Length_dyn), ctypes.byref(Length_newton),
        ctypes.byref(Status_dyn), ctypes.byref(Status_newton)
    )

    if Status_dyn.value != 0 or Status_newton.value != 0:
        if verbose:
            print("❌ Solver failed: Status_Dynamic =", Status_dyn.value,
                  "Status_Newton =", Status_newton.value)
        return False

    # Segment length deviation check
    points = np.vstack([P1, out_positions.reshape(-1, 3), P2])
    lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
    deviation = np.max(np.abs(lengths - np.mean(lengths)) / np.mean(lengths)) * 100

    if deviation > tol_percent:
        if verbose:
            print(f"❌ Length deviation too high: {deviation:.2f}%")
        return False

    # Gravity-aligned force comparison
    g_norm = np.linalg.norm(g_vec)
    if g_norm < 1e-6:
        if verbose:
            print("⚠️ Skipping force check: gravity vector too small.")
        return True  # Can't compare directionally aligned force
    g_hat = g_vec / g_norm
    expected_force = total_mass * g_norm
    allowed_error = force_tol * expected_force

    # Projected Newton-based force
    F_combined_n = F_P1_n + F_P2_n
    F_proj_n = np.dot(F_combined_n, g_hat)
    F_error_n = abs(F_proj_n - expected_force)

    # Projected Weight-based force
    F_combined_w = F_P1_w + F_P2_w
    F_proj_w = np.dot(F_combined_w, g_hat)
    F_error_w = abs(F_proj_w - expected_force)

    if F_error_n > allowed_error or F_error_w > allowed_error:
        if verbose:
            print("❌ Gravity-aligned force imbalance too high:")
            print(f"   Newton-based: {F_proj_n:.2f} N vs {expected_force:.2f} N "
                  f"(error: {F_error_n:.2e} N)")
            print(f"   Weight-based: {F_proj_w:.2f} N vs {expected_force:.2f} N "
                  f"(error: {F_error_w:.2e} N)")
        return False

    if verbose:
        print(f"✅ Success. Max segment deviation = {deviation:.2f}%, "
              f"Force errors: Newton = {F_error_n:.2e}, Weight = {F_error_w:.2e}")

    return True


# Run semi-random test batch
num_tests = 100
failures = []

for i in range(num_tests):
    print(f"\n--- Running Test Case {i + 1} ---")
    test = generate_test_case()
    success = run_and_test(**test, verbose=True)
    if not success:
        failures.append(test)

print(f"\nSummary: {num_tests - len(failures)} passed / {num_tests} total.")
if failures:
    print("\nFailed test inputs:")
    for f in failures:
        for k, v in f.items():
            print(f"{k} = {v}")
        print("---")
