import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer

np.random.seed(42)  # Fixed seed for reproducibility

# Load DLL
dll_path = './Resources/Library/solve_rope_length_lapak.dll'
lib = ctypes.CDLL(dll_path)

lib.solve_rope_length.argtypes = [
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # P1
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # P2
    ctypes.c_int,                                      # n
    ctypes.c_double,                                   # total_mass
    ctypes.c_double,                                   # length_factor
    ctypes.c_double,                                   # rope_diameter
    ctypes.c_double,                                   # youngs_modulus
    ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),  # g_vec
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # out_positions
    ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),  # F_P1_n
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # F_P2_n
    ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),  # F_P1_w
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # F_P2_w
    ctypes.POINTER(ctypes.c_double),                   # Length_initial
    ctypes.POINTER(ctypes.c_double),                   # Length_cat
    ctypes.POINTER(ctypes.c_double),                   # Length_dynamic
    ctypes.POINTER(ctypes.c_double),                   # Length_newton
    ctypes.POINTER(ctypes.c_int),                      # Status_dynamic
    ctypes.POINTER(ctypes.c_int),                      # Status_newton
    ctypes.c_int,                                      # debug_level
]
lib.solve_rope_length.restype = ctypes.c_int

# Define ctypes prototypes
lib.solve_rope_tension.argtypes = [
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # P1
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # P2
    ctypes.c_int,                                      # n
    ctypes.c_double,                                   # total_mass
    ctypes.c_double,                                   # rope_diameter
    ctypes.c_double,                                   # youngs_modulus
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # g_vec
    ctypes.c_double,                                   # F_target
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # out_positions
    ctypes.POINTER(ctypes.c_double),                   # out_length_factor
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # F_P1_out
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # F_P2_out
    ctypes.POINTER(ctypes.c_int),                      # Status_dynamic
    ctypes.POINTER(ctypes.c_int),                      # Status_newton
    ctypes.c_int                                        # debug_level
]
lib.solve_rope_tension.restype = ctypes.c_int

def test_case_tension():
    P1 = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    P2 = np.array([0.0, 10.0, 0.0], dtype=np.float64)
    n = 20
    total_mass = 1.0
    rope_diameter = 0.01
    youngs_modulus = 1e9
    g_vec = np.array([0.0, 0.0, -9.81], dtype=np.float64)
    F_target = 5.0  # Newtons

    return {
        "P1": P1, "P2": P2,
        "n": n, "total_mass": total_mass,
        "rope_diameter": rope_diameter,
        "youngs_modulus": youngs_modulus,
        "g_vec": g_vec,
        "F_target":F_target
    }

def run_and_test_tension(P1, P2, n, total_mass, rope_diameter,
                 youngs_modulus, g_vec, F_target):

    # === Outputs ===
    out_positions = np.zeros((n - 1) * 3, dtype=np.float64)
    out_length_factor = ctypes.c_double()
    F_P1_out = np.zeros(3, dtype=np.float64)
    F_P2_out = np.zeros(3, dtype=np.float64)
    Status_dynamic = ctypes.c_int()
    Status_newton = ctypes.c_int()

    # === Call DLL Function ===
    result = lib.solve_rope_tension(
        P1, P2,
        n, total_mass,
        rope_diameter, youngs_modulus,
        g_vec, F_target,
        out_positions,
        ctypes.byref(out_length_factor),
        F_P1_out, F_P2_out,
        ctypes.byref(Status_dynamic),
        ctypes.byref(Status_newton),
        3  # debug_level: 0–5
    )

    # === Output ===
    print("Return code:", result)
    print("Computed Length Factor:", out_length_factor.value)
    print("Force at P1:", F_P1_out)
    print("Force at P2:", F_P2_out)
    print("Dynamic relaxation status:", Status_dynamic.value)
    print("Newton solver status:", Status_newton.value)
    print("Node Positions:\n", out_positions.reshape(-1, 3))

def generate_weight_test():

    P1 =[0,0,0]
    P2 =[0,70,0]
    n = 20
    total_mass = 1
    length_factor = 1.3
    rope_diameter = 0.01
    youngs_modulus = 1e9
    g_vec = [0,0,-9.81]

    return {
        "P1": P1, "P2": P2,
        "n": n, "total_mass": total_mass,
        "length_factor": length_factor,
        "rope_diameter": rope_diameter,
        "youngs_modulus": youngs_modulus,
        "g_vec": g_vec
    }

def generate_test_case():
    P1 = np.random.uniform(-100, 100, 3).tolist()
    P2 = np.random.uniform(-100, 100, 3).tolist()
    n = np.random.randint(10, 100)
    length_factor = np.random.uniform(1.1, 2.5)
    rope_diameter = np.random.uniform(0.01, 0.08)
    density = np.random.uniform(500, 8000)
    youngs_modulus = np.random.uniform(1e9, 2e11)
    g_vec = np.random.uniform(-15, 15, 3).tolist()

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

def run_and_test_length(P1, P2, n, total_mass, length_factor, rope_diameter,
                 youngs_modulus, g_vec, tol_percent=2.5, force_tol=0.02, verbose=True):
    length_tol_percent = 2.5  # max allowed % difference between arc length reports

    P1 = np.array(P1, dtype=np.float64)
    P2 = np.array(P2, dtype=np.float64)
    g_vec = np.array(g_vec, dtype=np.float64)

    out_positions = np.zeros((n - 1) * 3, dtype=np.float64)
    F_P1_n = np.zeros(3, dtype=np.float64)
    F_P2_n = np.zeros(3, dtype=np.float64)
    F_P1_w = np.zeros(3, dtype=np.float64)
    F_P2_w = np.zeros(3, dtype=np.float64)

    Length_init = ctypes.c_double()
    Length_cat = ctypes.c_double()
    Length_dyn = ctypes.c_double()
    Length_newton = ctypes.c_double()
    Status_dyn = ctypes.c_int()
    Status_newton = ctypes.c_int()

    res = lib.solve_rope_length(
        P1, P2, n, total_mass, length_factor, rope_diameter, youngs_modulus,
        g_vec, out_positions, F_P1_n, F_P2_n, F_P1_w, F_P2_w,
        ctypes.byref(Length_init), ctypes.byref(Length_cat), ctypes.byref(Length_dyn), ctypes.byref(Length_newton),
        ctypes.byref(Status_dyn), ctypes.byref(Status_newton), debug_level
    )

    if Status_dyn.value != 0 or Status_newton.value != 0:
        reason = f"Solver failed: Status_Dynamic={Status_dyn.value}, Status_Newton={Status_newton.value}"
        if verbose:
            print(f"❌ {reason}")
        return False, reason

    # Check arc length consistency
    lengths = [Length_init.value, Length_dyn.value, Length_newton.value]
    mean_len = np.mean(lengths)
    max_dev = max(abs(l - mean_len) / mean_len * 100 for l in lengths)
    if max_dev > length_tol_percent:
        reason = (f"Arc length mismatch: Init = {Length_init.value:.3f}, "
                  f"Dynamic = {Length_dyn.value:.3f}, Newton = {Length_newton.value:.3f} "
                  f"(max deviation: {max_dev:.2f}%)")
        if verbose:
            print(f"❌ {reason}")
        return False, reason

    # Segment length deviation check
    points = np.vstack([P1, out_positions.reshape(-1, 3), P2])
    segment_lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
    deviation = np.max(np.abs(segment_lengths - np.mean(segment_lengths)) / np.mean(segment_lengths)) * 100
    if deviation > tol_percent:
        reason = f"Length deviation too high: {deviation:.2f}%"
        if verbose:
            print(f"❌ {reason}")
        return False, reason

    # Gravity-aligned force comparison
    g_norm = np.linalg.norm(g_vec)
    if g_norm < 1e-6:
        if verbose:
            print("⚠️ Skipping force check: gravity vector too small.")
        return True, ""

    g_hat = g_vec / g_norm
    expected_force = total_mass * g_norm
    allowed_error = force_tol * expected_force

    F_combined_n = F_P1_n + F_P2_n
    F_proj_n = np.dot(F_combined_n, g_hat)
    F_error_n = abs(F_proj_n - expected_force)

    F_combined_w = F_P1_w + F_P2_w
    F_proj_w = np.dot(F_combined_w, g_hat)
    F_error_w = abs(F_proj_w - expected_force)

    if F_error_n > allowed_error or F_error_w > allowed_error:
        reason = (f"Gravity-aligned force imbalance too high:\n"
                  f"   Newton-based: {F_proj_n:.2f} N vs {expected_force:.2f} N "
                  f"(error: {F_error_n:.2e} N)\n"
                  f"   Weight-based: {F_proj_w:.2f} N vs {expected_force:.2f} N "
                  f"(error: {F_error_w:.2e} N)")
        if verbose:
            print(f"❌ {reason}")
        return False, reason

    if verbose:
        print(f"✅ Success. Max segment deviation = {deviation:.2f}%, "
              f"Force errors: Weight = {F_error_w:.2e}, Springs = {F_error_n:.2e}")
    return True, ""

# Run semi-random test batch
num_tests = 1
failures = []
debug_level = 5

for i in range(num_tests):
    print(f"\n--- Running Test Case {i + 1} ---")
    test = test_case_tension()
    #test = generate_weight_test()
    #test = generate_test_case()
    #run_and_test_length(**test)
    run_and_test_tension(**test)
    #success, reason = run_and_test_length(**test, verbose=True)
    #if not success:
    #    test["reason"] = reason
    #    failures.append(test)

# Summary
print(f"\nSummary: {num_tests - len(failures)} passed / {num_tests} total.")
if failures:
    print("\nFailed test inputs:")
    for f in failures:
        for k, v in f.items():
            if k != "reason":
                print(f"{k} = {v}")
        print("Reason:", f["reason"])
        print("---")
