import ctypes
import numpy as np
import matplotlib.pyplot as plt

# === Load DLL ===
dll_path = "./Resources/Library/solve_rope_length_lapak.dll"
lib = ctypes.CDLL(dll_path)

# === Define DLL function signature ===
lib.generate_lambda_force_table.argtypes = [
    np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # P1
    np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # P2
    ctypes.c_int,            # n
    ctypes.c_double,         # total_mass
    ctypes.c_double,         # rope_diameter
    ctypes.c_double,         # youngs_modulus
    np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # g_vec
    ctypes.c_double,         # lambda_start
    ctypes.c_double,         # lambda_end
    ctypes.c_int,            # num_samples
    np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # lambda_out
    np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # F_P1_n_out
    np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # F_P2_n_out
    np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # F_P1_w_out
    np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # F_P2_w_out
    np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # L_init_out
    np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # L_cat_out
    np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # L_dyn_out
    np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # L_newton_out
    np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),     # status_dynamic_out
    np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),     # status_newton_out
]
lib.generate_lambda_force_table.restype = ctypes.c_int

# === Input parameters ===
P1 = np.array([0.0, 0.0, 0.0], dtype=np.float64)
P2 = np.array([70.0, 0.0, 70.0], dtype=np.float64)
n = 20
total_mass = 500.0
rope_diameter = 0.02
youngs_modulus = 2e9
g_vec = np.array([0.0, 0.0, -9.81], dtype=np.float64)

lambda_start = 1.00001
lambda_end = 3
num_samples = 100

# === Allocate output arrays ===
lambda_out = np.zeros(num_samples, dtype=np.float64)

F_P1_n_out = np.zeros(num_samples * 3, dtype=np.float64)
F_P2_n_out = np.zeros(num_samples * 3, dtype=np.float64)
F_P1_w_out = np.zeros(num_samples * 3, dtype=np.float64)
F_P2_w_out = np.zeros(num_samples * 3, dtype=np.float64)

L_init_out   = np.zeros(num_samples, dtype=np.float64)
L_cat_out    = np.zeros(num_samples, dtype=np.float64)
L_dyn_out    = np.zeros(num_samples, dtype=np.float64)
L_newton_out = np.zeros(num_samples, dtype=np.float64)

status_dynamic_out = np.zeros(num_samples, dtype=np.int32)
status_newton_out  = np.zeros(num_samples, dtype=np.int32)

# === Call the DLL ===
ret_code = lib.generate_lambda_force_table(
    P1, P2, n, total_mass, rope_diameter, youngs_modulus, g_vec,
    lambda_start, lambda_end, num_samples,
    lambda_out,
    F_P1_n_out, F_P2_n_out,
    F_P1_w_out, F_P2_w_out,
    L_init_out, L_cat_out, L_dyn_out, L_newton_out,
    status_dynamic_out, status_newton_out
)

if ret_code != 0:
    raise RuntimeError(f"generate_lambda_force_table failed with code {ret_code}")

# === Plotting ===
def plot_lambda_force_table(
    lambda_out,
    F_P1_n_out, F_P2_n_out,
    F_P1_w_out, F_P2_w_out,
    L_dyn_out, L_newton_out,
    status_dynamic_out, status_newton_out,
    force_clip_threshold=None  # New: optional clipping threshold
):
    num_samples = len(lambda_out)

    # === Compute magnitudes ===
    F_P1_n_mag = np.linalg.norm(F_P1_n_out.reshape(-1, 3), axis=1)
    F_P2_n_mag = np.linalg.norm(F_P2_n_out.reshape(-1, 3), axis=1)
    F_P1_w_mag = np.linalg.norm(F_P1_w_out.reshape(-1, 3), axis=1)
    F_P2_w_mag = np.linalg.norm(F_P2_w_out.reshape(-1, 3), axis=1)

    # === Clip forces if threshold is provided ===
    if force_clip_threshold is not None:
        F_P1_n_mag = np.clip(F_P1_n_mag, 0, force_clip_threshold)
        F_P2_n_mag = np.clip(F_P2_n_mag, 0, force_clip_threshold)
        F_P1_w_mag = np.clip(F_P1_w_mag, 0, force_clip_threshold)
        F_P2_w_mag = np.clip(F_P2_w_mag, 0, force_clip_threshold)

    failed = (np.array(status_dynamic_out) != 0) | (np.array(status_newton_out) != 0)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True, gridspec_kw={'hspace': 0.3})

    # === Force plot ===
    ax1.plot(lambda_out, F_P1_n_mag, label="‖F_P1‖ (Newton)", lw=2)
    ax1.plot(lambda_out, F_P2_n_mag, label="‖F_P2‖ (Newton)", lw=2, ls='--')
    ax1.plot(lambda_out, F_P1_w_mag, label="‖F_P1‖ (Weight)", lw=2)
    ax1.plot(lambda_out, F_P2_w_mag, label="‖F_P2‖ (Weight)", lw=2, ls='--')

    ax1.scatter(lambda_out[failed], F_P1_n_mag[failed], color='red', label='Fail (F_P1)', marker='x')
    ax1.scatter(lambda_out[failed], F_P2_n_mag[failed], color='black', label='Fail (F_P2)', marker='x')

    ax1.set_ylabel("Force Magnitude [N]")
    ax1.set_title("Force Magnitudes vs Length Factor λ")
    ax1.grid(True)
    ax1.legend()

    # === Length plot ===
    ax2.plot(lambda_out, L_dyn_out, label="L_dynamic", lw=2)
    ax2.plot(lambda_out, L_newton_out, label="L_newton", lw=2, ls='--')
    ax2.scatter(lambda_out[failed], L_dyn_out[failed], color='red', marker='x', label='Fail (L_dyn)')
    ax2.scatter(lambda_out[failed], L_newton_out[failed], color='black', marker='x', label='Fail (L_newton)')

    ax2.set_xlabel("Length Factor λ")
    ax2.set_ylabel("Rope Length [m]")
    ax2.set_title("Rope Lengths vs λ")
    ax2.grid(True)
    ax2.legend()

    plt.tight_layout()
    plt.show()


# === Execute plot ===
plot_lambda_force_table(
    lambda_out,
    F_P1_n_out, F_P2_n_out,
    F_P1_w_out, F_P2_w_out,
    L_dyn_out, L_newton_out,
    status_dynamic_out, status_newton_out,
    force_clip_threshold=10000# clip any force > 100 N
)