import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer
import matplotlib.pyplot as plt

from matplotlib.patches import Patch  # for proxy legend entry




from matplotlib.collections import LineCollection
from matplotlib.patches import Patch


# Load DLL
dll_path = './Resources/Library/solve_rope_length_lapak.so'
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

lib.solve_rope_tension_arc.argtypes = [
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
    ctypes.POINTER(ctypes.c_int),                      # Status_newton
    ctypes.c_int                                        # debug_level
]
lib.solve_rope_tension_arc.restype = ctypes.c_int

lib.generate_lambda_force_table.argtypes = [
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # P1
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # P2
    ctypes.c_int,                                      # n
    ctypes.c_double,                                   # total_mass
    ctypes.c_double,                                   # rope_diameter
    ctypes.c_double,                                   # youngs_modulus
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # g_vec
    ctypes.c_double,                                   # lambda_min
    ctypes.c_double,                                   # lambda_max
    ctypes.c_int,                                      # num_samples
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # lambda_out
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # F_P1_magnitude_out
    ctypes.c_int                                       # debug_level
]
lib.generate_lambda_force_table.restype = ctypes.c_int

def test_case_tension():
    P1 = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    P2 = np.array([0.0, 70.0, 0.0], dtype=np.float64)
    n = 10
    total_mass = 7.0
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
                         youngs_modulus, g_vec, F_target,
                         force_tol=0.02, tol_percent=2.5, verbose=True, debug_level=3):
    """
    Executes solve_rope_tension and evaluates its correctness.
    Checks: return code, force error, planarity.
    """
    P1 = np.asarray(P1, dtype=np.float64)
    P2 = np.asarray(P2, dtype=np.float64)
    g_vec = np.asarray(g_vec, dtype=np.float64)

    # === Outputs ===
    out_positions = np.zeros((n - 1) * 3, dtype=np.float64)
    out_length_factor = ctypes.c_double()
    F_P1_out = np.zeros(3, dtype=np.float64)
    F_P2_out = np.zeros(3, dtype=np.float64)
    Status_dynamic = ctypes.c_int()
    Status_newton = ctypes.c_int()

    # === DLL Call ===
    result = lib.solve_rope_tension(
        P1, P2, n, total_mass,
        rope_diameter, youngs_modulus,
        g_vec, F_target,
        out_positions,
        ctypes.byref(out_length_factor),
        F_P1_out, F_P2_out,
        ctypes.byref(Status_dynamic),
        ctypes.byref(Status_newton),
        debug_level
    )

    # === Early failure check ===
    if result != 0 or Status_dynamic.value != 0 or Status_newton.value != 0:
        reason = f"Solver failed: ret={result}, dyn={Status_dynamic.value}, newton={Status_newton.value}"
        if verbose:
            print(f"❌ {reason}")
        return False, reason, out_positions, Status_dynamic.value, Status_newton.value

    # === Force check ===
    F_proj = np.dot(F_P1_out, g_vec) / np.linalg.norm(g_vec)
    error_force = abs(F_proj - F_target)
    allowed_error = force_tol * F_target

    if error_force > allowed_error:
        reason = (f"Force mismatch: projected P1 = {F_proj:.2f} N, "
                  f"expected = {F_target:.2f} N (error = {error_force:.2e} N)")
        if verbose:
            print(f"❌ {reason}")
        return False, reason, out_positions, Status_dynamic.value, Status_newton.value

    # === Planarity check ===
    cord_vec = P2 - P1
    cord_vec /= np.linalg.norm(cord_vec)
    g_proj = g_vec - np.dot(g_vec, cord_vec) * cord_vec

    if np.linalg.norm(g_proj) < 1e-8:
        planarity_error = None
        if verbose:
            print("⚠️ Skipping planarity check (g ‖ cord).")
    else:
        plane_normal = np.cross(cord_vec, g_proj)
        plane_normal /= np.linalg.norm(plane_normal)

        all_nodes = np.vstack([P1, out_positions.reshape(-1, 3), P2])
        distances = np.dot(all_nodes - P1, plane_normal)
        planarity_error = np.max(np.abs(distances))

        allowed_deviation = tol_percent / 100 * np.linalg.norm(P2 - P1)
        if planarity_error > allowed_deviation:
            reason = f"Planarity deviation too high: {planarity_error:.2e} m"
            if verbose:
                print(f"❌ {reason}")
            return False, reason, out_positions, Status_dynamic.value, Status_newton.value

    # === Success ===
    if verbose:
        print(f"✅ Success. Force error = {error_force:.2e}, Length Factor = {out_length_factor.value:.5f}")
        if planarity_error is not None:
            print(f"   Max Planarity Deviation = {planarity_error:.2e} m")

    return True, "", out_positions, Status_dynamic.value, Status_newton.value

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
    length_tol_percent = 2.5

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
        return False, reason, out_positions, Status_dyn.value, Status_newton.value

    lengths = [Length_init.value, Length_dyn.value, Length_newton.value]
    mean_len = np.mean(lengths)
    max_dev = max(abs(l - mean_len) / mean_len * 100 for l in lengths)
    if max_dev > length_tol_percent:
        reason = (f"Arc length mismatch: Init = {Length_init.value:.3f}, "
                  f"Dynamic = {Length_dyn.value:.3f}, Newton = {Length_newton.value:.3f} "
                  f"(max deviation: {max_dev:.2f}%)")
        if verbose:
            print(f"❌ {reason}")
        return False, reason, out_positions, Status_dyn.value, Status_newton.value

    points = np.vstack([P1, out_positions.reshape(-1, 3), P2])
    segment_lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
    deviation = np.max(np.abs(segment_lengths - np.mean(segment_lengths)) / np.mean(segment_lengths)) * 100
    if deviation > tol_percent:
        reason = f"Length deviation too high: {deviation:.2f}%"
        if verbose:
            print(f"❌ {reason}")
        return False, reason, out_positions, Status_dyn.value, Status_newton.value

    g_norm = np.linalg.norm(g_vec)
    if g_norm < 1e-6:
        if verbose:
            print("⚠️ Skipping force check: gravity vector too small.")
        return True, "", out_positions, Status_dyn.value, Status_newton.value

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
        return False, reason, out_positions, Status_dyn.value, Status_newton.value

    cord_vec = P2 - P1
    cord_vec /= np.linalg.norm(cord_vec)
    g_proj = g_vec - np.dot(g_vec, cord_vec) * cord_vec

    if np.linalg.norm(g_proj) < 1e-8:
        if verbose:
            print("⚠️ Skipping plane deviation check: g_vec parallel to cord.")
        max_plane_deviation = None
    else:
        plane_normal = np.cross(cord_vec, g_proj)
        plane_normal /= np.linalg.norm(plane_normal)

        all_nodes = np.vstack([P1, out_positions.reshape(-1, 3), P2])
        distances = np.dot(all_nodes - P1, plane_normal)
        max_plane_deviation = np.max(np.abs(distances))

        rope_span = np.linalg.norm(P2 - P1)
        allowed_plane_deviation = 1e-4 * rope_span

        if max_plane_deviation > allowed_plane_deviation:
            reason = f"Node deviation from rope–gravity plane too high: {max_plane_deviation:.2e} m"
            if verbose:
                print(f"❌ {reason}")
            return False, reason, out_positions, Status_dyn.value, Status_newton.value

    if verbose:
        extra = ""
        if max_plane_deviation is not None:
            extra = f", Max plane deviation = {max_plane_deviation:.2e} m"
        print(f"✅ Success. Max segment deviation = {deviation:.2f}%, "
              f"Force errors: Weight = {F_error_w:.2e}, Springs = {F_error_n:.2e}{extra}")
    return True, "", out_positions, Status_dyn.value, Status_newton.value

def plot_rope_2d(P1, P2, x, g_vec,
                                                     total_mass=None, length_factor=None,
                                                     rope_diameter=None, youngs_modulus=None,
                                                     status_tuple=None, fail_reason=None):
    """
    2D rope plot with diagnostics and failure reason shown in the subtitle (title second line).
    Includes sag, span-to-sag ratio, planarity error, and segment coloring.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.patches import Patch

    x = np.asarray(x).reshape(-1, 3)
    n = len(x)
    P1 = np.asarray(P1)
    P2 = np.asarray(P2)
    g_vec = np.asarray(g_vec)
    g_mag = np.linalg.norm(g_vec)

    if status_tuple:
        print(f"Solver status: dyn={status_tuple[0]}, newton={status_tuple[1]}")
        success = all(s == 0 for s in status_tuple)
        print("Success" if success else "Solver failed")

    # Frame setup
    ez = g_vec / g_mag
    ex = P2 - P1
    ex -= np.dot(ex, ez) * ez
    ex /= np.linalg.norm(ex)
    ey = np.cross(ez, ex)
    R = np.vstack([ex, ey, ez])

    x_rel = x - x[0]
    x_proj = x_rel @ R.T
    x2d = x_proj[:, [0, 2]]
    x2d[:, 1] *= -1

    seg_lengths = np.linalg.norm(np.diff(x, axis=0), axis=1)
    arc_length = np.sum(seg_lengths)
    cord_length = np.linalg.norm(P2 - P1)
    area = np.pi * (rope_diameter ** 2) / 4.0 if rope_diameter else None
    volume = area * (cord_length * length_factor) if area and length_factor else None
    density = total_mass / volume if volume else None

    all_nodes = np.vstack([P1, x, P2])
    sag = np.max(all_nodes[:, 2]) - np.min(all_nodes[:, 2])
    span_to_sag = cord_length / sag if sag > 1e-6 else float('inf')

    u = (P2 - P1) / np.linalg.norm(P2 - P1)
    g_proj = g_vec - np.dot(g_vec, u) * u
    if np.linalg.norm(g_proj) < 1e-8:
        planarity_error = None
    else:
        normal = np.cross(u, g_proj)
        normal /= np.linalg.norm(normal)
        distances = np.dot(all_nodes - P1, normal)
        planarity_error = np.max(np.abs(distances))

    # Plot
    fig, ax2d = plt.subplots(figsize=(10, 5))
    norm = plt.Normalize(vmin=np.min(seg_lengths), vmax=np.max(seg_lengths))
    cmap = plt.get_cmap("viridis")
    colors = cmap(norm(seg_lengths))
    segments = [[x2d[i], x2d[i + 1]] for i in range(n - 1)]
    lc = LineCollection(segments, colors=colors, linewidths=2)
    ax2d.add_collection(lc)
    ax2d.scatter(*x2d.T, c="k", s=2, zorder=3)

    title = "Rope Projection into Gravity-Aligned Plane"
    subtitle = f"Reason: {fail_reason}" if fail_reason else None
    ax2d.set_title(f"{title}\n{subtitle}" if subtitle else title)

    ax2d.set_xlabel("Horizontal along rope")
    ax2d.set_ylabel("Vertical (gravity down)")
    ax2d.axis('equal')

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax2d, orientation='vertical', pad=0.02)
    cbar.set_label("Segment Length [m]")
    def small_patch(label):
            return Patch(facecolor='none', edgecolor='k', label=label, linewidth=0.5)

    # Legend
    legend_elements = []
    if total_mass is not None:
        legend_elements.append(small_patch(f"Mass: {total_mass:.2f} kg"))
    if density:
        legend_elements.append(small_patch(f"Density: {density:.1f} kg/m³"))
    if length_factor:
        legend_elements.append(small_patch(f"Length Factor: {length_factor:.4f}"))
    legend_elements.append(small_patch(f"Arc Length: {arc_length:.3f} m"))
    legend_elements.append(small_patch(f"Cord Length: {cord_length:.3f} m"))
    legend_elements.append(small_patch(f"Span-to-Sag Ratio: {span_to_sag:.2f}"))
    legend_elements.append(small_patch(f"Sag Depth: {sag:.3f} m"))
    if rope_diameter:
        legend_elements.append(small_patch(f"Diameter: {rope_diameter*1000:.3f} mm"))
    if youngs_modulus:
        legend_elements.append(small_patch(f"Young's Modulus: {youngs_modulus/1e9:.3f} GPa"))
    legend_elements.append(small_patch(f"Gravity Magnitude: {g_mag:.2f} m/s²"))
    if planarity_error is not None:
        legend_elements.append(small_patch(f"Planarity Error: {planarity_error:.2e} m"))
    legend_elements.append(small_patch(f"Nodes: {n}"))
    if status_tuple:
        legend_elements.append(small_patch(f"Status: dyn={status_tuple[0]}, newton={status_tuple[1]}"))

    ax2d.legend(handles=legend_elements, loc='lower left', frameon=True, handlelength=1.0, borderpad=0.5, labelspacing=0.3)
    plt.tight_layout()
    plt.show()

def run_and_test_tension_arc(P1, P2, n, total_mass, rope_diameter,
                              youngs_modulus, g_vec, F_target,
                              force_tol=0.02, tol_percent=2.5,
                              verbose=True, debug_level=3):
    """
    Executes solve_rope_tension_arc and evaluates correctness.
    """
    P1 = np.asarray(P1, dtype=np.float64)
    P2 = np.asarray(P2, dtype=np.float64)
    g_vec = np.asarray(g_vec, dtype=np.float64)

    # --- Compute initial guess for s0 ---
    L_straight = np.linalg.norm(P2 - P1)*1.5
    s0 = np.full(n, L_straight / n, dtype=np.float64)

    # --- Outputs ---
    out_positions = np.zeros(n * 3, dtype=np.float64)
    out_length_factor = ctypes.c_double()
    F_P1_out = np.zeros(3, dtype=np.float64)
    F_P2_out = np.zeros(3, dtype=np.float64)
    Status_newton = ctypes.c_int()

    # --- DLL call ---
    result = lib.solve_rope_tension_arc(
                                        P1, P2, n, total_mass,
                                        rope_diameter, youngs_modulus,
                                        g_vec, F_target,
                                        out_positions,
                                        ctypes.byref(out_length_factor),
                                        F_P1_out, F_P2_out,
                                        ctypes.byref(Status_newton),
                                        debug_level )

    if result != 0 or Status_newton.value != 0:
        reason = f"Solver failed: ret={result}, newton={Status_newton.value}"
        if verbose:
            print(f"❌ {reason}")
        return False, reason, out_positions, 0, Status_newton.value

    # --- Force check ---
    g_norm = np.linalg.norm(g_vec)
    g_unit = g_vec / g_norm
    force_proj = np.dot(F_P1_out, g_unit)
    error_force = abs(force_proj - F_target)
    allowed_error = force_tol * F_target

    if error_force > allowed_error:
        reason = (f"Force mismatch: projected P1 = {force_proj:.2f} N, "
                  f"expected = {F_target:.2f} N (error = {error_force:.2e} N)")
        if verbose:
            print(f"❌ {reason}")
        return False, reason, out_positions, 0, Status_newton.value

    # --- Planarity check ---
    cord_vec = P2 - P1
    cord_vec /= np.linalg.norm(cord_vec)
    g_proj = g_vec - np.dot(g_vec, cord_vec) * cord_vec
    if np.linalg.norm(g_proj) < 1e-8:
        planarity_error = None
        if verbose:
            print("⚠️ Skipping planarity check (g ‖ cord).")
    else:
        plane_normal = np.cross(cord_vec, g_proj)
        plane_normal /= np.linalg.norm(plane_normal)
        all_nodes = out_positions.reshape(-1, 3)
        distances = np.dot(all_nodes - P1, plane_normal)
        planarity_error = np.max(np.abs(distances))
        allowed_deviation = tol_percent / 100 * np.linalg.norm(P2 - P1)
        if planarity_error > allowed_deviation:
            reason = f"Planarity deviation too high: {planarity_error:.2e} m"
            if verbose:
                print(f"❌ {reason}")
            return False, reason, out_positions, 0, Status_newton.value

    # --- Success ---
    if verbose:
        print(f"✅ Success. Force error = {error_force:.2e}, Length Factor = {out_length_factor.value:.5f}")
        if planarity_error is not None:
            print(f"   Max Planarity Deviation = {planarity_error:.2e} m")

    return True, "", out_positions, 0, Status_newton.value

def run_lambda_force_table(P1, P2, n, total_mass, rope_diameter, youngs_modulus, g_vec,
                           lambda_min=1.0, lambda_max=3.0, num_samples=100):
    P1 = np.array(P1, dtype=np.float64)
    P2 = np.array(P2, dtype=np.float64)
    g_vec = np.array(g_vec, dtype=np.float64)

    lambda_out = np.zeros(num_samples, dtype=np.float64)
    F_P1_mags = np.zeros(num_samples, dtype=np.float64)

    ret = lib.generate_lambda_force_table(
        P1, P2, n, total_mass, rope_diameter, youngs_modulus, g_vec,
        lambda_min, lambda_max, num_samples, lambda_out, F_P1_mags, debug_level
    )

    if ret != 0:
        return None, f"generate_lambda_force_table failed with code {ret}"

    return lambda_out, F_P1_mags


# Run semi-random test batch
np.random.seed(42)  # Fixed seed for reproducibility
num_tests = 1
failures = []
debug_level = 3
#debug_mode = 0  # no debug
#debug_mode = 1  #weight case
#debug_mode = 2  #length test silent
#debug_mode = 3  #length test with plot
#debug_mode = 4  #tension test 
#debug_mode = 5   # arc test
debug_mode = 6   # arc test


for i in range(num_tests):
    print(f"\n--- Running Test Case {i + 1} ---")
    if debug_mode == 0:
        debug_level = 0
    elif debug_mode == 1:
        test = generate_weight_test()
        success, reason, positions, status_dyn, status_newton = run_and_test_length(**test, verbose=True)
        if not success:
            test["reason"] = reason
            failures.append(test)
    elif debug_mode == 2:
        test = generate_test_case()
        success, reason, positions, status_dyn, status_newton = run_and_test_length(**test, verbose=True)
        if not success:
            test["reason"] = reason
            test["positions"] = positions.copy()
            test["status_dyn"] = status_dyn
            test["status_newton"] = status_newton
            failures.append(test)
    elif debug_mode == 3:
        test = generate_test_case()
        success, reason, positions, status_dyn, status_newton = run_and_test_length(**test, verbose=True)
        if not success:
            test["reason"] = reason
            test["positions"] = positions.copy()
            failures.append(test)        #        
    elif debug_mode == 4:
        test = test_case_tension()
        success, reason, positions, status_dyn, status_newton = run_and_test_tension(**test, verbose=True, debug_level=debug_level)
        if not success:
            test["reason"] = reason
            test["positions"] = positions.copy()
            test["status_dyn"] = status_dyn
            test["status_newton"] = status_newton
            failures.append(test)
    elif debug_mode == 5:
        test = test_case_tension()
        success, reason, positions, status_dyn, status_newton = run_and_test_tension_arc(**test, verbose=True, debug_level=debug_level)
        if not success:
            test["reason"] = reason
            test["positions"] = positions.copy()
            test["status_dyn"] = status_dyn
            test["status_newton"] = status_newton
            failures.append(test)
    elif debug_mode == 6:
        test = test_case_tension()
        lambda_vals, F_vals = run_lambda_force_table(
            test["P1"], test["P2"], test["n"], test["total_mass"],
            test["rope_diameter"], test["youngs_modulus"], test["g_vec"],
            lambda_min=1.003, lambda_max=4.5, num_samples=200 )

        if lambda_vals is not None:
            plt.plot(lambda_vals, F_vals, label="||F_P1|| vs λ")
            plt.xlabel("Length Factor λ")
            plt.ylabel("Force Magnitude at P1 [N]")
            plt.title("Force vs Rope Length Factor")
            plt.grid(True)
            plt.legend()
            plt.tight_layout()
            plt.show()
        else:
            print(F_vals)




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
        plot_rope_2d(
            f["P1"], f["P2"], f["positions"],
            g_vec=np.array(f["g_vec"]),
            total_mass=f.get("total_mass"),
            length_factor=f.get("length_factor"),
            rope_diameter=f.get("rope_diameter"),
            youngs_modulus=f.get("youngs_modulus"),
            status_tuple=(f.get("status_dyn", -99), f.get("status_newton", -99)),
            fail_reason=f.get("reason", "")
        )