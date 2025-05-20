# -*- coding: utf-8 -*-
"""
@file    run_multiple_solve_rope_lenght_tests.py
@brief   Python ctypes wrapper, test harness, and visualizer for solve_rope_length_lapak.dll

This script wraps the `solve_rope_length_lapak.dll` shared library using Python's `ctypes` 
interface and provides tools to:

- Generate rope force vs. length factor tables via `generate_lambda_force_table`.
- Plot Newton- and weight-based endpoint forces as a function of length factor λ.
- Visualize rope length outcomes from different solver stages.
- Run randomized test cases in parallel using `multiprocessing`.
- Automatically detect and report failed or timeout cases.
- Optionally visualize results for each test case.
- Build inverse models from (λ, force) samples.

@details
The wrapped DLL implements a 3D spring-mass rope simulation using:
  - Initial position estimation via catenary or cubic parabola
  - Dynamic relaxation for coarse equilibrium
  - Newton refinement for final convergence
  - Endpoint force reporting (Newton-based and weight-based)
  - Full diagnostic output and solver status flags

@requires
  - Python 3.6+
  - Numpy
  - Matplotlib
  - Scipy (optional, for spline fitting or inversion)
  - solve_rope_length_lapak.dll (compiled from C source)

@usage
  Run this script directly with command-line options:
    > python solve_rope_length-py.py -n 10 --segments 30 --debug --plot

@command_line_args
  -n / --num_tests      : Number of randomized test cases to run
  --segments            : Number of discrete rope segments (default = 20)
  --debug               : Enable verbose logging during DLL call
  --plot                : Show plots for all successful test cases

@functions
  - run_single_test         : Run one rope simulation with given input and capture results
  - run_parallel_tests      : Run N randomized simulations and collect statistics
  - plot_force_model        : Plot ‖F_P1‖ vs λ using Newton and weight-based decompositions
  - plot_inverse_mapping    : Visualize inverted (force → rope length) model
  - generate_random_case    : Generate a randomized test scenario (endpoints, mass, etc.)

@dll_interface
  Function: generate_lambda_force_table
    Inputs:
      - Endpoints: P1, P2 (3D)
      - Discretization count: n
      - Physical parameters: total_mass, rope_diameter, youngs_modulus
      - Gravity vector: g_vec
      - Sampling range: lambda_start, lambda_end, num_samples
      - Debug level
    Outputs (via pointers):
      - λ samples
      - Forces: Newton and weight-based (P1 and P2)
      - Lengths: Init, Cat, Dynamic, Newton
      - Status codes: dynamic and Newton phases

@outputs
  For each test, the script returns or plots:
    - λ values
    - Force magnitudes at P1/P2 via two decomposition methods
    - Solver success/failure status
    - Dynamic and Newton rope lengths
    - Optional plots for diagnostics

@limitations
  - Large n (>200) may increase runtime and risk instability.
  - Requires compiled DLL compatible with Python calling conventions.

@authors
  - Martin Burger
  - Generated with support from ChatGPT

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
dll_path = "../Resources/Library/solve_rope_length_lapak.dll"
lib = ctypes.CDLL(dll_path)

# === Define DLL function signature ===
lib.generate_lambda_force_table.argtypes = [
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # P1
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # P2
    ctypes.c_int,                                      # n
    ctypes.c_double,                                   # total_mass
    ctypes.c_double,                                   # rope_diameter
    ctypes.c_double,                                   # youngs_modulus
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # g_vec
    ctypes.c_double,                                   # lambda_start
    ctypes.c_double,                                   # lambda_end
    ctypes.c_int,                                      # num_samples
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # lambda_out
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # F_P1_n_out
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # F_P2_n_out
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # F_P1_w_out
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # F_P2_w_out
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # L_init_out
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # L_cat_out
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # L_dyn_out
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # L_newton_out
    ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),     # status_dynamic_out
    ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),     # status_newton_out
    ctypes.c_int                                       # debug_level
]

lib.generate_lambda_force_table.restype = ctypes.c_int

#############################################
#           Worker and Solver Logic         #
#############################################

def _dll_worker(args, return_dict, debug_level, log_queue=None):
    try:
        if log_queue:
            log_queue.put("[Worker] Started")

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

        if log_queue:
            log_queue.put(f"[Worker] Inputs:")
            log_queue.put(f"  span = {np.linalg.norm(P2 - P1):.2f}")
            log_queue.put(f"  total_mass = {total_mass:.2f}")
            log_queue.put(f"  rope_diameter = {rope_diameter:.4f}, n = {n}")
            log_queue.put(f"  P1 = {P1.tolist()}")
            log_queue.put(f"  P2 = {P2.tolist()}")
            log_queue.put(f"  g_vec = {g_vec.tolist()}")

        ret = lib.generate_lambda_force_table(
            P1, P2, n, total_mass, rope_diameter, youngs_modulus, g_vec,
            lambda_start, lambda_end, num_samples,
            lambda_out,
            F_P1_n_out, F_P2_n_out, F_P1_w_out, F_P2_w_out,
            L_init_out, L_cat_out, L_dyn_out, L_newton_out,
            status_dyn_out, status_newton_out, debug_level
        )

        return_dict['ret'] = ret
        return_dict['lambda_out'] = lambda_out
        return_dict['F_P1_n_out'] = F_P1_n_out
        return_dict['F_P1_w_out'] = F_P1_w_out
        return_dict['status_dyn_out'] = status_dyn_out
        return_dict['status_newton_out'] = status_newton_out
        return_dict['P1'] = P1
        return_dict['P2'] = P2

        if log_queue:
            log_queue.put(f"[Worker] DLL call returned: {ret}")
            log_queue.put("[Worker] Done")

    except Exception as e:
        if log_queue:
            log_queue.put(f"[Worker Error] {e}")
            log_queue.put(traceback.format_exc())

def run_single_test(P1, P2, n, total_mass, rope_diameter, youngs_modulus, g_vec,
                    lambda_start=1.00001, lambda_end=3.5, num_samples=100, timeout_sec=20, debug_level=0):
    args = (P1, P2, n, total_mass, rope_diameter, youngs_modulus, g_vec, lambda_start, lambda_end, num_samples)
    mgr = mp.Manager()
    return_dict = mgr.dict()
    log_queue = mgr.Queue()

    proc = mp.Process(target=_dll_worker, args=(args, return_dict, debug_level, log_queue))
    proc.start()
    proc.join(timeout=timeout_sec)

    logs = []
    while not log_queue.empty():
        try:
            logs.append(log_queue.get_nowait())
        except:
            break

    if proc.is_alive():
        proc.terminate()
        proc.join()
        logs.append("⏰ Timeout occurred.")
        return "timeout", logs

    if return_dict.get("ret", -999) != 0:
        logs.append(f"❌ DLL returned non-zero: {return_dict.get('ret')}")
        return "error", logs

    logs.append("✅ DLL call succeeded.")
    return "success", return_dict

##############################################
#           Plotting                         #
##############################################

def plot_force_model(lambda_out, F_P1_n_out, F_P1_w_out, status_dyn_out, status_newton_out, P1, P2, total_mass, rope_diameter):
    """
    Plots the average force vs. lambda using the Newton and weight-based force outputs.
    Only valid and agreeing results are shown.
    """
    F_P1_n_mag = np.linalg.norm(F_P1_n_out.reshape(-1, 3), axis=1)
    F_P1_w_mag = np.linalg.norm(F_P1_w_out.reshape(-1, 3), axis=1)
    avg_force = 0.5 * (F_P1_n_mag + F_P1_w_mag)
    rel_diff = np.abs(F_P1_n_mag - F_P1_w_mag) / np.maximum(F_P1_n_mag, 1e-6)

    valid = (np.array(status_dyn_out) == 0) & (np.array(status_newton_out) == 0)
    agree = rel_diff < 0.05
    mask = valid & agree

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(lambda_out[mask], avg_force[mask], 'bo', label="Valid & agreeing (mean of Newton & weight)")
    ax.set_xlabel("Length Factor λ")
    ax.set_ylabel("‖F_P1‖ [N]")
    ax.set_title(f"Force vs. λ | span={np.linalg.norm(P2 - P1):.1f} m, mass={total_mass:.1f} kg, dia={rope_diameter:.3f} m")
    ax.grid(True)
    ax.legend()
    plt.tight_layout()
    return fig

def plot_inverse_mapping(L_all, F_sweep):
    """
    Plots the inverted force-to-length mapping results.
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    for i, L_vals in enumerate(L_all):
        for L in L_vals:
            ax.plot(L, F_sweep[i], 'bo', markersize=3)
    ax.set_xlabel("Rope Length [m]")
    ax.set_ylabel("Force at P1 [N]")
    ax.set_title("Inverse Force Model: Rope Length vs. Force")
    ax.grid(True)
    plt.tight_layout()
    return fig

#############################################
#           Parallel Test Harness           #
#############################################

def random_unit_vector():
    v = np.random.normal(size=3)
    return v / np.linalg.norm(v)

def generate_random_case():
    P1 = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    direction = random_unit_vector()
    span = np.random.uniform(30.0, 100.0)
    P2 = P1 + direction * span
    total_mass = np.random.uniform(5.0, 150.0)
    rope_diameter = np.random.uniform(0.01, 0.05)
    g_vec = np.array([0.0, 0.0, -9.81], dtype=np.float64)
    return P1, P2, span, total_mass, rope_diameter, g_vec

def run_parallel_tests(n_cases=20, n=20, youngs_modulus=2e9, debug_level=0):
    results = []
    for i in range(n_cases):
        P1, P2, span, total_mass, rope_diameter, g_vec = generate_random_case()
        status, payload = run_single_test(
            P1, P2, n, total_mass, rope_diameter, youngs_modulus, g_vec, debug_level=debug_level)
        results.append((status, span, total_mass, rope_diameter, payload))
        print(f"[{i+1}/{n_cases}] Result: {status}")
        if status != "success":
             for log in payload:
                 print("  " + log)
    return results

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Run random rope solver tests and visualize results.")
    parser.add_argument('-n', '--num_tests', type=int, default=5, help='Number of random test cases to run')
    parser.add_argument('--segments', type=int, default=20, help='Number of rope segments')
    parser.add_argument('--debug', action='store_true', help='Enable debug logging')
    parser.add_argument('--plot', action='store_true', help='Show plots for each test')
    args = parser.parse_args()

    debug_level = 3 if args.debug else 0

    results = run_parallel_tests(
        n_cases=args.num_tests,
        n=args.segments,
        youngs_modulus=2e9,
        debug_level=debug_level
    )

    num_success = sum(1 for r in results if r[0] == "success")
    num_error = sum(1 for r in results if r[0] == "error")
    num_timeout = sum(1 for r in results if r[0] == "timeout")

    print("\n=== Summary ===")
    print(f"Success  : {num_success}")
    print(f"Error    : {num_error}")
    print(f"Timeouts : {num_timeout}")

    if args.plot:
        for i, (status, span, total_mass, rope_diameter, data) in enumerate(results):
            if status != "success":
                continue

            fig = plot_force_model(
                data['lambda_out'],
                data['F_P1_n_out'],
                data['F_P1_w_out'],
                data['status_dyn_out'],
                data['status_newton_out'],
                data.get('P1', np.array([0.0, 0.0, 0.0])),  # Fallback if missing
                data.get('P2', np.array([0.0, 0.0, -1.0])),
                total_mass,
                rope_diameter
            )
            fig.suptitle(f"Test #{i+1} | span={span:.1f} m")
            plt.show()
