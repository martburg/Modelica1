# -*- coding: utf-8 -*-
"""
@file    inverse_force_model.py
@brief   Wrapper and inverse solver for rope length vs. tension using solve_rope_length_lapak.dll

This script loads the shared rope simulation DLL, samples tension values over stretch factors (λ),
builds a composite force model, and inverts it to find the rope length(s) needed to achieve a
desired tension.

@details
The workflow includes:
  - Running the DLL to generate (λ, force) tables for a given physical rope configuration.
  - Filtering and blending Newton- and weight-based endpoint forces.
  - Constructing a smooth composite force model using splines, Hermite patches, and physical tails.
  - Inverting the force model using bisection to determine stretch factor(s) that match a target force.
  - Plotting both the forward and inverse mappings.
  - Performing batch testing using random rope parameters.

@capabilities
  - Automatic sampling and diagnostics via multiprocessing.
  - Support for failure detection and timeouts.
  - Filtering of poor quality samples via agreement and monotonicity.
  - Generates diagnostic plots of force vs λ and inverted force vs length.

@requires
  - Numpy
  - Matplotlib
  - Scipy
  - solve_rope_length_lapak.dll

@note
  This script is suitable for validation of force models, inversion testing, and fitting
  approximate rope behavior in parameter studies.

@author
  Martin Burger

@date
  2025-05-20
"""

import ctypes
import numpy as np
import matplotlib.pyplot as plt
from numpy.ctypeslib import ndpointer
from scipy.interpolate import UnivariateSpline
import multiprocessing as mp
import traceback
import time

# === Load DLL ===
dll_path = "./Resources/Library/solve_rope_length_lapak.dll"
lib = ctypes.CDLL(dll_path)

# === Define DLL function signature ===
lib.generate_lambda_force_table.argtypes = [
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),                      # P1
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),                      # P2
    ctypes.c_int,                                                          # n
    ctypes.c_double,                                                       # total_mass  
    ctypes.c_double,                                                       # rope_diameter
    ctypes.c_double,                                                       # youngs_modulus
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),                      # g_vec
    ctypes.c_double,                                                       # lambda_start
    ctypes.c_double,                                                       # lambda_end
    ctypes.c_int,                                                          # num_samples
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),                      # lambda_out  
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),                      # F_P1_n_out
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),                      # F_P2_n_out
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),                      # F_P1_w_out
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),                      # F_P2_w_out
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),                      # L_init_out
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),                      # L_cat_out
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),                      # L_dyn_out
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),                      # L_newton_out
    ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),                         # status_dynamic_out
    ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),                         # status_newton_out
    ctypes.c_int                                                           # debug_level    
]

lib.generate_lambda_force_table.restype = ctypes.c_int

def invert_force_model(F_model, lambda_min=1.0, lambda_max=3.5, F_target=100.0, tol=1e-3, num_points=1000):
   """
    @brief Invert a force model to find λ such that F_model(λ) ≈ F_target.

    @param[in] F_model     Function mapping λ to force
    @param[in] lambda_min  Minimum λ to search
    @param[in] lambda_max  Maximum λ to search
    @param[in] F_target    Target force magnitude
    @param[in] tol         Tolerance for force matching
    @param[in] num_points  Resolution of λ sampling grid

    @returns List of λ values where F_model(λ) ≈ F_target (zero-crossings of F - F_target)
    """
    λ_space = np.linspace(lambda_min, lambda_max, num_points)
    F_space = F_model(λ_space)
    delta = F_space - F_target

    roots = []
    for i in range(len(delta) - 1):
        if delta[i] * delta[i + 1] < 0:
            λ0, λ1 = λ_space[i], λ_space[i + 1]
            for _ in range(20):
                λ_mid = 0.5 * (λ0 + λ1)
                f_mid = F_model(λ_mid) - F_target
                if abs(f_mid) < tol:
                    break
                if (F_model(λ0) - F_target) * f_mid < 0:
                    λ1 = λ_mid
                else:
                    λ0 = λ_mid
            roots.append(λ_mid)
    return sorted(roots)

def filter_monotonic_decreasing(lambda_out, force_array, base_mask, lower=1.0, upper=1.5):
    """
    @brief Filters a force array to retain only monotonic decreasing behavior in [lower, upper].

    @param[in] lambda_out   Array of λ values
    @param[in] force_array  Corresponding force magnitudes
    @param[in] base_mask    Initial boolean mask of valid points
    @param[in] lower        Lower λ bound for enforcing monotonicity
    @param[in] upper        Upper λ bound for enforcing monotonicity

    @returns New mask with non-monotonic entries in [lower, upper] removed
    """
    new_mask = base_mask.copy()
    in_range = (lambda_out >= lower) & (lambda_out <= upper) & base_mask
    indices = np.where(in_range)[0]

    for i in range(1, len(indices)):
        idx_prev = indices[i - 1]
        idx_curr = indices[i]
        if force_array[idx_curr] > force_array[idx_prev]:
            new_mask[idx_curr] = False  # violates monotonicity

    return new_mask

def build_composite_force_model(
    lambda_out, F_P1_n_out, F_P1_w_out,
    status_dynamic_out, status_newton_out,
    P1, P2, total_mass, rope_diameter,
    lambda_switch=1.6, lambda_tail_start=3, lambda_max=3.5,
    tolerance=0.05, clip_above_start=True, g_mag=9.81):

    """
    @brief Construct a composite force model F(λ) from simulation samples.

    @details
      Combines spline fitting in [1, λ_switch], Hermite interpolation in [λ_switch, λ_tail_start],
      and a physically derived linear tail beyond λ_tail_start.

    @param[in] lambda_out         Sampled λ values
    @param[in] F_P1_n_out         Newton-based force vectors at P1
    @param[in] F_P1_w_out         Weight-based force vectors at P1
    @param[in] status_dynamic_out Status codes from dynamic solver
    @param[in] status_newton_out  Status codes from Newton solver
    @param[in] P1, P2             Endpoint coordinates
    @param[in] total_mass         Total rope mass
    @param[in] rope_diameter      Rope diameter [m]
    @param[in] lambda_switch      Switch point between spline and Hermite patch
    @param[in] lambda_tail_start  End of Hermite region, start of physical tail
    @param[in] lambda_max         Maximum λ considered
    @param[in] tolerance          Max relative difference allowed between Newton/Weight forces
    @param[in] clip_above_start   Filter points with force > 1.05 × initial force
    @param[in] g_mag              Gravitational constant (default 9.81)

    @returns
      - F_model(λ): callable composite force model
      - span: endpoint distance
      - f_min, f_max: min and max valid force values
    """


    F_P1_n_mag = np.linalg.norm(F_P1_n_out.reshape(-1, 3), axis=1)
    F_P1_w_mag = np.linalg.norm(F_P1_w_out.reshape(-1, 3), axis=1)
    avg_force = 0.5 * (F_P1_n_mag + F_P1_w_mag)
    rel_diff = np.abs(F_P1_n_mag - F_P1_w_mag) / np.maximum(F_P1_n_mag, 1e-6)
    valid = (np.array(status_dynamic_out) == 0) & (np.array(status_newton_out) == 0)
    agree = rel_diff < tolerance
    below_start = avg_force < avg_force[0] * 1.05 if clip_above_start else np.ones_like(avg_force, dtype=bool)
    good = valid & agree & below_start
    # Apply monotonic filter in [1.0, 1.5]
    good = filter_monotonic_decreasing(lambda_out, avg_force, good, 1.0, 1.5)

    # Final data for spline
    λ_good = lambda_out[good]
    F_good = avg_force[good]

        # Diagnostic summary
    total_points = len(lambda_out)
    num_valid = np.sum(valid)
    num_agree = np.sum(valid & agree)
    num_below = np.sum(valid & agree & below_start)
    num_final = np.sum(good)

    print("=== Force Model Sample Filtering Summary ===")
    print(f"Total samples       : {total_points}")
    print(f"Valid solver output : {num_valid}")
    print(f"Newton ≈ Weight     : {num_agree}")
    print(f"Below rise threshold: {num_below}")
    print(f"Final used (incl. monotonic): {num_final}")
    print("============================================")

    spline_mask = λ_good <= lambda_switch
    spline = UnivariateSpline(λ_good[spline_mask], F_good[spline_mask], s=0)
    λ0 = lambda_switch
    λ1 = lambda_tail_start
    f0 = spline(λ0)
    df0 = spline.derivative()(λ0)
    in_tail_range = (lambda_out >= λ0) & (lambda_out <= λ1) & good
    λ_tail_fit = lambda_out[in_tail_range]
    F_tail_fit = avg_force[in_tail_range]
    slope, intercept = np.polyfit(λ_tail_fit, F_tail_fit, deg=1)
    f1 = slope * λ1 + intercept
    df1 = slope

    span = np.linalg.norm(np.array(P2) - np.array(P1))
    A = np.pi * rope_diameter**2 / 4
    total_length_max = lambda_max * span
    density = total_mass / (total_length_max * A)
    tail_slope_physical = g_mag * density * A * span * 0.5

    def hermite_patch(λ):
        t = (λ - λ0) / (λ1 - λ0)
        h00 = 2 * t**3 - 3 * t**2 + 1
        h10 = t**3 - 2 * t**2 + t
        h01 = -2 * t**3 + 3 * t**2
        h11 = t**3 - t**2
        return f0 * h00 + df0 * (λ1 - λ0) * h10 + f1 * h01 + df1 * (λ1 - λ0) * h11

    def tail_force(λ):
        return f1 + tail_slope_physical * (λ - λ1)

    def F_model(λ):
        λ = np.asarray(λ)
        return np.piecewise(
            λ,
            [λ <= λ0, (λ > λ0) & (λ <= λ1), λ > λ1],
            [lambda x: spline(x), lambda x: hermite_patch(x), lambda x: tail_force(x)]
        )

    λ_dense = np.linspace(min(lambda_out), lambda_max, 600)
    F_dense = F_model(λ_dense)

    plt.figure(figsize=(10, 5))
    plt.plot(lambda_out[good], avg_force[good], 'ko', label="Valid Data (mean of Newton & Weight)")
    plt.plot(λ_dense, F_dense, 'r-', label="Composite Force Model", lw=2)
    plt.axvline(λ0, color='gray', linestyle='--', label="λ_switch")
    plt.axvline(λ1, color='lightgray', linestyle='--', label="λ_tail_start")
    plt.xlabel("Length Factor λ")
    plt.ylabel("‖F_P1‖ [N]")
    plt.title("Force Model | span=%.1f m, mass=%.1f kg, dia=%.3f m" % (
        np.linalg.norm(np.array(P2) - np.array(P1)), total_mass, rope_diameter))
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    return F_model, span, min(avg_force[good]), max(avg_force[good])

def _dll_worker(args, return_dict, debug_level, log_queue=None):
    """
    @brief Internal multiprocessing worker that wraps a call to the DLL function.

    @param[in] args         Tuple of all solver inputs
    @param[out] return_dict Multiprocessing dictionary with results and status
    @param[in] debug_level  Logging verbosity level (0 = silent)
    @param[in,out] log_queue Optional queue to stream logs back to parent process
    """
    try:
        if log_queue:
            log_queue.put("[Worker] ✅ Reached beginning of worker function")

        (P1, P2, n, total_mass, rope_diameter, youngs_modulus, g_vec,
         lambda_start, lambda_end, num_samples) = args

        lambda_out = np.zeros(num_samples)
        F_P1_n_out = np.zeros(num_samples * 3)
        F_P1_w_out = np.zeros(num_samples * 3)
        F_P2_n_out = np.zeros(num_samples * 3)
        F_P2_w_out = np.zeros(num_samples * 3)
        L_init_out = np.zeros(num_samples)
        L_cat_out = np.zeros(num_samples)
        L_dyn_out = np.zeros(num_samples)
        L_newton_out = np.zeros(num_samples)
        status_dyn_out = np.zeros(num_samples, dtype=np.int32)
        status_newton_out = np.zeros(num_samples, dtype=np.int32)

        log_queue.put("[Worker] ✅ Finished allocating NumPy arrays")

        ret = lib.generate_lambda_force_table(
            P1, P2, n, total_mass, rope_diameter, youngs_modulus, g_vec,
            lambda_start, lambda_end, num_samples,
            lambda_out,
            F_P1_n_out, F_P2_n_out, F_P1_w_out, F_P2_w_out,
            L_init_out, L_cat_out, L_dyn_out, L_newton_out,
            status_dyn_out, status_newton_out, debug_level
        )

        log_queue.put(f"[Worker] ✅ DLL call returned: {ret}")

        return_dict['ret'] = ret
        return_dict['lambda_out'] = lambda_out
        return_dict['F_P1_n_out'] = F_P1_n_out
        return_dict['F_P1_w_out'] = F_P1_w_out
        return_dict['status_dyn_out'] = status_dyn_out
        return_dict['status_newton_out'] = status_newton_out

        log_queue.put("[Worker] ✅ Finished processing.")

    except Exception as e:
        if log_queue:
            log_queue.put(f"[Worker] ❌ Exception occurred: {e}")
            log_queue.put(traceback.format_exc())

def run_force_model_and_plot_inverse(P1, P2,
                                     n, total_mass, rope_diameter, youngs_modulus,
                                    g_vec, timeout_sec=20, debug_level=0):
    """
    @brief Runs a full forward solve and builds an inverse model of rope length vs. force.

    @details
      - Executes a forward rope simulation via DLL
      - Builds a composite force model F(λ)
      - Inverts F(λ) to find λ values for a sweep of target forces
      - Plots rope length vs. force using the span-scaled λ

    @param[in] P1, P2         Endpoint coordinates
    @param[in] n              Number of segments
    @param[in] total_mass     Rope mass
    @param[in] rope_diameter  Rope diameter
    @param[in] youngs_modulus Elastic modulus
    @param[in] g_vec          Gravity vector
    @param[in] timeout_sec    Timeout in seconds for each solve
    @param[in] debug_level    Logging verbosity

    @returns "success", "timeout", or "error"
    """
    
    args = (P1, P2, n, total_mass, rope_diameter, youngs_modulus, g_vec, 1.00001, 3.5, 100)
    mgr = mp.Manager()
    return_dict = mgr.dict()
    log_queue = mp.Manager().Queue()
    proc = mp.Process(target=_dll_worker, args=(args, return_dict, debug_level, log_queue))
    proc.start()
    # Optional: flush any early logs while still running

    time.sleep(0.2)  # allow buffer to fill
    while not log_queue.empty():
        try:
            msg = log_queue.get_nowait()
            print(msg)
        except:
            break

    proc.join(timeout=timeout_sec)
    # Drain log queue after join
    while not log_queue.empty():
        try:
            msg = log_queue.get_nowait()
            print(msg)
        except:
            break

    if proc.is_alive():
        proc.terminate()
        proc.join()
        print("⏰ Test timed out after", timeout_sec, "seconds.")
        return "timeout"

    if return_dict.get("ret", -999) != 0:
        print("❌ Solver returned error code:", return_dict["ret"])
        return "error"

    # Build force model and get span, f_min, f_max
    F_model, span, f_min, f_max = build_composite_force_model(
        return_dict["lambda_out"],
        return_dict["F_P1_n_out"],
        return_dict["F_P1_w_out"],
        return_dict["status_dyn_out"],
        return_dict["status_newton_out"],
        P1, P2, total_mass, rope_diameter
    )

    # Sweep force from f_max to f_min
    F_sweep = np.linspace(f_max, f_min, 200)
    L_all = []
    for F in F_sweep:
        roots = invert_force_model(F_model, 1.0, 3.5, F_target=F)
        # L_all.append([r * span for r in roots]) resale to span
        L_all.append([r  for r in roots])

    # Plot: force vs rope length (may have multiple solutions per force)
    plt.figure(figsize=(10, 5))
    for i, L_vals in enumerate(L_all):
        for L in L_vals:
            plt.plot(L, F_sweep[i], 'bo', markersize=3)
    plt.xlabel("Rope Length [m]")
    plt.ylabel("Force at P1 [N]")
    plt.title(f"Inverse Force Model | span={span:.1f} m, mass={total_mass:.1f} kg, dia={rope_diameter:.3f} m")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    return "success"

def random_unit_vector():
    v = np.random.normal(size=3)
    return v / np.linalg.norm(v)

if __name__ == '__main__':

    total_tests = 0
    num_success = 0
    num_timeout = 0
    num_error = 0
    debug_level = 0

    # Run multiple tests with random parameters

    for i in range(20):
        P1 = np.array([0.0, 0.0, 0.0], dtype=np.float64)
        direction = random_unit_vector()
        span = np.random.uniform(30.0, 100.0)
        P2 = P1 + direction * span
        total_mass = np.random.uniform(5.0, 150.0)
        rope_diameter = np.random.uniform(0.01, 0.05)
        g_vec = np.array([0.0, 0.0, -9.81], dtype=np.float64)

        print(f"\n=== Running test {i+1} ===")
        result = run_force_model_and_plot_inverse(
            P1, P2, 20, total_mass, rope_diameter, 2e9, g_vec,
            timeout_sec=20, debug_level=debug_level)

        total_tests += 1
        if result == "success":
            num_success += 1
        elif result == "timeout":
            num_timeout += 1
        else:
            num_error += 1

        print(f"Test {i + 1} result: {result}")

    # Final summary
    print("\n=== Overall Test Summary ===")
    print(f"Total tests run     : {total_tests}")
    print(f"Successful          : {num_success}")
    print(f"Timed out (>20 sec) : {num_timeout}")
    print(f"Solver errors       : {num_error}")
    print("============================")