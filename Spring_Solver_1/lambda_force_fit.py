import ctypes
from scipy.interpolate import UnivariateSpline
from scipy.stats import linregress
import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import UnivariateSpline
import numpy as np
import matplotlib.pyplot as plt

def build_composite_force_model(
    lambda_out,
    F_P1_n_out, F_P1_w_out,
    status_dynamic_out, status_newton_out,
    P1, P2,
    total_mass, rope_diameter,
    lambda_switch=1.8,
    lambda_tail_start=2.6,
    lambda_max=3.0,
    tolerance=0.05,
    clip_above_start=True,
    g_mag=9.81
):
    """
    Builds a C¹-smooth force model:
    - spline up to lambda_switch
    - Hermite patch to lambda_tail_start
    - analytical linear tail from rope weight
    """
    F_P1_n_mag = np.linalg.norm(F_P1_n_out.reshape(-1, 3), axis=1)
    F_P1_w_mag = np.linalg.norm(F_P1_w_out.reshape(-1, 3), axis=1)

    avg_force = 0.5 * (F_P1_n_mag + F_P1_w_mag)
    rel_diff = np.abs(F_P1_n_mag - F_P1_w_mag) / np.maximum(F_P1_n_mag, 1e-6)

    valid = (np.array(status_dynamic_out) == 0) & (np.array(status_newton_out) == 0)
    agree = rel_diff < tolerance
    below_start = avg_force < avg_force[0] * 1.05 if clip_above_start else np.ones_like(avg_force, dtype=bool)

    good = valid & agree & below_start
    λ_good = lambda_out[good]
    F_good = avg_force[good]

    # Spline
    spline_mask = λ_good <= lambda_switch
    spline = UnivariateSpline(λ_good[spline_mask], F_good[spline_mask], s=0)
    λ0 = lambda_switch
    λ1 = lambda_tail_start
    f0 = spline(λ0)
    df0 = spline.derivative()(λ0)

    # Linear extrapolation from valid points in [λ0, λ1] to estimate f1
    in_tail_range = (lambda_out >= λ0) & (lambda_out <= λ1) & good
    λ_tail_fit = lambda_out[in_tail_range]
    F_tail_fit = avg_force[in_tail_range]
    slope, intercept = np.polyfit(λ_tail_fit, F_tail_fit, deg=1)
    f1 = slope * λ1 + intercept
    df1 = slope

    # Geometry and material
    span = np.linalg.norm(np.array(P2) - np.array(P1))
    A = np.pi * rope_diameter**2 / 4
    total_length_max = lambda_max * span
    density = total_mass / (total_length_max * A)

    # Tail slope from additional mass
    tail_slope_physical = g_mag * density * A * span * 0.5

    # Hermite patch
    def hermite_patch(λ):
        t = (λ - λ0) / (λ1 - λ0)
        h00 = 2 * t**3 - 3 * t**2 + 1
        h10 = t**3 - 2 * t**2 + t
        h01 = -2 * t**3 + 3 * t**2
        h11 = t**3 - t**2
        return (
            f0 * h00 + df0 * (λ1 - λ0) * h10 +
            f1 * h01 + df1 * (λ1 - λ0) * h11
        )

    # Linear tail
    def tail_force(λ):
        return f1 + tail_slope_physical * (λ - λ1)

    # Piecewise force model
    def F_model(λ):
        λ = np.asarray(λ)
        return np.piecewise(
            λ,
            [λ <= λ0, (λ > λ0) & (λ <= λ1), λ > λ1],
            [lambda x: spline(x),
             lambda x: hermite_patch(x),
             lambda x: tail_force(x)]
        )

    # Plot
    λ_dense = np.linspace(min(lambda_out), lambda_max, 600)
    F_dense = F_model(λ_dense)

    plt.figure(figsize=(10, 5))
    plt.plot(lambda_out[good], avg_force[good], 'ko', label="Valid Data (mean of Newton & Weight)")
    plt.plot(λ_dense, F_dense, 'r-', label="Composite Force Model", lw=2)
    plt.axvline(λ0, color='gray', linestyle='--', label="λ_switch")
    plt.axvline(λ1, color='lightgray', linestyle='--', label="λ_tail_start")
    plt.xlabel("Length Factor λ")
    plt.ylabel("‖F_P1‖ [N]")
    plt.title("Spline + Hermite + Linear Force Model")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    return F_model

# === Load DLL ===
dll_path = "./Resources/Library/solve_rope_length_lapak.so"
lib = ctypes.CDLL(dll_path)

# === Define DLL function signature ===
lib.generate_lambda_force_table.argtypes = [
    np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # P1
    np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # P2
    ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double,
    np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # g_vec
    ctypes.c_double, ctypes.c_double, ctypes.c_int,
    np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
    *[np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS") for _ in range(8)],
    *[np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS") for _ in range(2)],
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
lambda_end = 3.0
num_samples = 100

# === Output arrays ===
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

# === Call DLL ===
ret = lib.generate_lambda_force_table(
    P1, P2, n, total_mass, rope_diameter, youngs_modulus, g_vec,
    lambda_start, lambda_end, num_samples,
    lambda_out,
    F_P1_n_out, F_P2_n_out, F_P1_w_out, F_P2_w_out,
    L_init_out, L_cat_out, L_dyn_out, L_newton_out,
    status_dyn_out, status_newton_out
)
if ret != 0:
    raise RuntimeError(f"DLL call failed: {ret}")

F_model = build_composite_force_model(
    lambda_out,
    F_P1_n_out, F_P1_w_out,
    status_dyn_out, status_newton_out,
    P1, P2,
    total_mass, rope_diameter,
    lambda_switch=1.8,
    lambda_tail_start=3,
    lambda_max=3.0
)

# === Postprocess ===
F_P1_n = F_P1_n_out.reshape(-1, 3)
F_P1_w = F_P1_w_out.reshape(-1, 3)
F_P1_n_mag = np.linalg.norm(F_P1_n, axis=1)
F_P1_w_mag = np.linalg.norm(F_P1_w, axis=1)

# Mean force, filtering extreme differences
valid = (np.abs(F_P1_n_mag - F_P1_w_mag) / np.maximum(F_P1_w_mag, 1e-6)) < 0.05
F_mean = 0.5 * (F_P1_n_mag + F_P1_w_mag)

# Clip out large forces > F_start
F_start = F_mean[0]
valid &= F_mean < F_start * 1.01  # allow 1% margin

# Region selection
region_I = lambda_out < 1.25
region_III = lambda_out > 1.75
valid_I = valid & region_I
valid_III = valid & region_III

# === Fit functions ===
def catenary_like(lmbd, a, b, c):
    return a * np.cosh(b * (lmbd - 1)) + c

def linear_rise(lmbd, a, b):
    return a + b * (lmbd - 1.75)

# === Fit and plot ===
popt_I, _ = curve_fit(catenary_like, lambda_out[valid_I], F_mean[valid_I], maxfev=10000)
popt_III, _ = curve_fit(linear_rise, lambda_out[valid_III], F_mean[valid_III])

lambda_fit = np.linspace(lambda_start, lambda_end, 300)
F_fit_I = catenary_like(lambda_fit, *popt_I)
F_fit_III = linear_rise(lambda_fit, *popt_III)

# === Plot ===
plt.figure(figsize=(10, 6))
plt.plot(lambda_out, F_P1_n_mag, label="‖F_P1‖ Newton", lw=1)
plt.plot(lambda_out, F_P1_w_mag, label="‖F_P1‖ Weight", lw=1)
plt.plot(lambda_out[valid], F_mean[valid], 'k.', label="Mean Force (valid)")

plt.plot(lambda_fit, F_fit_I, 'r--', lw=2, label="Fit Region I (λ < 1.25)")
plt.plot(lambda_fit, F_fit_III, 'b--', lw=2, label="Fit Region III (λ > 1.75)")

plt.axvline(1.25, color='gray', linestyle=':', lw=1)
plt.axvline(1.75, color='gray', linestyle=':', lw=1)

plt.xlabel("Length Factor λ")
plt.ylabel("Force Magnitude [N]")
plt.title("Piecewise Fit of Force vs Length Factor")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()


