import numpy as np
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def solve_catenary_parameter(D, sag):
    def equation(a):
        return a * np.cosh(D / (2 * a)) - a - sag
    result = root_scalar(equation, bracket=[1e-3, 1e6], method='brentq')
    if not result.converged:
        raise RuntimeError("Failed to solve for catenary parameter a")
    return result.root

def compute_catenary_arc_length(D, a):
    return 2 * a * np.sinh(D / (2 * a))

def compute_vertical_force(total_mass, g, arc_length):
    return total_mass * g / arc_length  # weight per unit length

def compute_tension_components_corrected(a, w, arc_length):
    T_horiz = a * w
    T_vert = w * arc_length / 2
    T_total = np.hypot(T_horiz, T_vert)
    return T_total, T_vert, T_horiz

def generate_catenary_points(D, a, num_points=100):
    y = np.linspace(0, D, num_points)
    z = a * np.cosh((y - D/2) / a) - a
    return y, z

# Optional: Your rope node positions here
# For example (y_i, z_i) from your rope model (insert real values):
# Rope nodes from final solver output
rope_y = np.array([
    2.976355, 6.239358, 9.841626, 13.847721, 18.335048, 23.390943,
    29.099551, 35.506746, 42.553013, 50.000000, 57.446987, 64.493254,
    70.900449, 76.609057, 81.664952, 86.152279, 90.158374, 93.760642,
    97.023645
])
rope_z = np.array([
   -6.889731, -13.647919, -20.231030, -26.575985, -32.589720, -38.133480,
   -43.001946, -46.904976, -49.480373, -50.387660, -49.480373, -46.904976,
   -43.001946, -38.133480, -32.589720, -26.575985, -20.231030, -13.647919,
   -6.889731
])

# Parameters
D = 100.0            # Horizontal span
sag = 50.39          # Depth of sag
total_mass = 1.0     # kg
g = 9.81             # m/s²

# Solve catenary
a = solve_catenary_parameter(D, sag)
arc_len = compute_catenary_arc_length(D, a)
w = compute_vertical_force(total_mass, g, arc_len)
T, T_vert, T_horiz = compute_tension_components_corrected(a, w, arc_len)

# Print summary
print(f"Catenary parameter a:        {a:.6f} m")
print(f"Total arc length:            {arc_len:.6f} m")
print(f"Weight per unit length:      {w:.6f} N/m")
print(f"Total weight:                {total_mass * g:.6f} N")
print(f"Endpoint tension:            {T:.6f} N")
print(f"  Vertical component:        {T_vert:.6f} N")
print(f"  Horizontal component:      {T_horiz:.6f} N")

# Generate catenary curve
cat_y, cat_z = generate_catenary_points(D, a)

cat_z = cat_z

# Shift catenary to align with rope sag
rope_sag_min = np.min(rope_z)
cat_sag_min = np.min(cat_z)
z_shift = rope_sag_min - cat_sag_min
cat_z += z_shift


# Plot
plt.figure(figsize=(10, 6))
plt.plot(cat_y, cat_z, label="Catenary (theoretical)", linewidth=2)
plt.plot(rope_y, rope_z, 'o-', label="Spring-mass rope (simulated)", color='orange')
plt.xlabel("Horizontal position (m)")
plt.ylabel("Vertical position (m)")
plt.title("Rope Shape Comparison: Catenary vs Spring-Mass Model")
plt.grid(True)
plt.legend()
#plt.gca().invert_yaxis()
plt.axis("equal")
plt.tight_layout()
plt.show()
