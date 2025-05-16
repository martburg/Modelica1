import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D



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
    P2 = np.array([0.0, 70.0, 0.0], dtype=np.float64)
    n = 20
    total_mass = 1.0
    rope_diameter = 0.01
    youngs_modulus = 1e9
    g_vec = np.array([0.0, 0.0, -9.81], dtype=np.float64)
    F_target = 50.0  # Newtons

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
    P2 =[70,0,0]
    n = 20
    total_mass = 5
    length_factor = 1.6
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
        return False, reason, out_positions

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
        return False, reason, out_positions

    # Segment length deviation check
    points = np.vstack([P1, out_positions.reshape(-1, 3), P2])
    segment_lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
    deviation = np.max(np.abs(segment_lengths - np.mean(segment_lengths)) / np.mean(segment_lengths)) * 100
    if deviation > tol_percent:
        reason = f"Length deviation too high: {deviation:.2f}%"
        if verbose:
            print(f"❌ {reason}")
        return False, reason, out_positions

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
        return False, reason, out_positions

    if verbose:
        print(f"✅ Success. Max segment deviation = {deviation:.2f}%, "
              f"Force errors: Weight = {F_error_w:.2e}, Springs = {F_error_n:.2e}")
    return True, "", out_positions

def plot_rope_colored(P1, P2, positions, g_vec=np.array([0, 0, -1])):
    P1 = np.asarray(P1)
    P2 = np.asarray(P2)
    g_vec = np.asarray(g_vec)
    nodes = np.vstack([P1, positions.reshape(-1, 3), P2])

    # --- Define gravity-aligned rotation matrix ---
    def get_gravity_aligned_rotation(P1, P2, g_vec):
        z_unit = -g_vec / np.linalg.norm(g_vec)
        x_raw = P2 - P1
        x_unit = x_raw / np.linalg.norm(x_raw)

        # Make x_unit orthogonal to z_unit
        x_unit = x_unit - np.dot(x_unit, z_unit) * z_unit
        x_unit /= np.linalg.norm(x_unit)

        y_unit = np.cross(z_unit, x_unit)
        y_unit /= np.linalg.norm(y_unit)

        R = np.column_stack([x_unit, y_unit, z_unit])
        return R

    # --- Apply rotation to align gravity with Z′ ---
    R = get_gravity_aligned_rotation(P1, P2, g_vec)
    nodes_rot = (R.T @ nodes.T).T  # Apply active rotation

    # === Print transformed positions ===
    print("\n[INFO] Transformed (gravity-aligned) node positions:")
    for i, pt in enumerate(nodes_rot):
        print(f"Node {i:2d}: [{pt[0]: .3f}, {pt[1]: .3f}, {pt[2]: .3f}]")

    # === Check planarity ===
    y_dev = np.std(nodes_rot[:, 1])
    print(f"[DEBUG] Std deviation in rotated Y (should be near 0): {y_dev:.6f}")

    # === Stretch coloring ===
    lengths = np.linalg.norm(np.diff(nodes_rot, axis=0), axis=1)
    mean_len = np.mean(lengths)
    stretch = (lengths - mean_len) / mean_len
    norm = plt.Normalize(vmin=np.min(stretch), vmax=np.max(stretch))
    colors = cm.viridis(norm(stretch))

    # === Plot ===
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection='3d')
    for i in range(len(nodes_rot) - 1):
        seg = np.vstack([nodes_rot[i], nodes_rot[i + 1]])
        ax.plot(seg[:, 0], seg[:, 1], seg[:, 2], color=colors[i])

    # Show original straight rope (rotated)
    ref_line = np.vstack([P1, P2]) @ R
    ax.plot(ref_line[:, 0], ref_line[:, 1], ref_line[:, 2], '--', color='gray', label='Straight Line')

    ax.set_xlabel("X′ (rope span)")
    ax.set_ylabel("Y′ (planarity deviation)")
    ax.set_zlabel("Z′ (aligned with gravity ↓)")
    ax.set_title('Rope with Gravity-Aligned Z′ Axis and Stretch Coloring')

    mappable = cm.ScalarMappable(cmap='viridis', norm=norm)
    mappable.set_array([])
    cbar = fig.colorbar(mappable, ax=ax, shrink=0.6)
    cbar.set_label('Relative Stretch')
    ax.legend()
    plt.tight_layout()
    plt.pause(0.1)

    plt.show()



# Run semi-random test batch
np.random.seed(42)  # Fixed seed for reproducibility
num_tests = 10
failures = []
debug_level = 3
#debug_mode = 0  # no debug
#debug_mode = 1  #weight case
#debug_mode = 2  #length test silent
debug_mode = 3  #length test with plot
#debug_mode = 4  #tension test 
#debug_mode = 3  #length test with plot

for i in range(num_tests):
    print(f"\n--- Running Test Case {i + 1} ---")
    if debug_mode == 0:
        debug_level = 0
    elif debug_mode == 1:
        test = generate_weight_test()
        success, reason, positions = run_and_test_length(**test, verbose=True)
        if not success:
            test["reason"] = reason
            failures.append(test)
    elif debug_mode == 2:
        test = generate_test_case()
        success, reason, positions = run_and_test_length(**test, verbose=True)
        if not success:
            test["reason"] = reason
            failures.append(test)
    elif debug_mode == 3:
        test = generate_test_case()
        success, reason, positions = run_and_test_length(**test, verbose=True)
        #### inset plot code 
        if not success:
            test["reason"] = reason
            test["positions"] = positions.copy()
            failures.append(test)        #        
    elif debug_mode == 4:
        test = test_case_tension()
        run_and_test_tension()

# Summary
print(f"\nSummary: {num_tests - len(failures)} passed / {num_tests} total.")
if failures:
    print("\nFailed test inputs:")
    for f in failures:
        for k, v in f.items():
            if k not in ("reason", "positions"):
                print(f"{k} = {v}")
        print("Reason:", f["reason"])
        print("---")
    for f in failures:
        if "positions" in f:
            plot_rope_colored(np.array(f["P1"]), np.array(f["P2"]), f["positions"])
