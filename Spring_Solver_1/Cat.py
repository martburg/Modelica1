import numpy as np
import matplotlib.pyplot as plt

# Inputs
P1 = np.array([0.0, 0.0, 0.0])
P2 = np.array([0.0, 0.1, 10.0])
n = 10
c = 100.0
total_mass = 20.0
m = total_mass / n
g_vec = np.array([0.0, 0.0, -9.81])
g_norm = np.linalg.norm(g_vec)
e_g = g_vec / g_norm

# Compute rope direction
dx_line = P2 - P1
L = np.linalg.norm(dx_line)
e_rope = dx_line / L

# Project gravity onto rope direction
dot_g_rope = np.dot(e_g, e_rope)
g_perp = e_g - dot_g_rope * e_rope
g_perp_norm = np.linalg.norm(g_perp)
if g_perp_norm > 1e-12:
    g_perp /= g_perp_norm
else:
    g_perp[:] = 0.0

# Adapted rest lengths
s0_adapted = np.array([(L / n) + (m * g_norm * (n - i)) / c for i in range(n)])
s0 = L / n  # uniform reference length
T0 = c * s0
a = (T0 / (m * g_norm)) if (T0 > 0 and g_norm > 0) else 1.0
print(a)

# Compute initial relaxed positions
positions = []
for i in range(n):
    t = (i + 1) / n
    base = (1 - t) * P1 + t * P2
    print(base)
    s = t * L
    sag_amount = a * (np.cosh(s / a) - 1.0)
    print(sag_amount)
    pos = base + sag_amount * e_g  # Apply sag in gravity direction

    positions.append(pos)

positions = np.array(positions)

# Plot
y = positions[:, 1]
z = positions[:, 2]

plt.figure(figsize=(6, 8))
plt.plot(y, z, marker='o', linestyle='-', label="Catenary-initialized Positions (Gravity-Aligned)")
plt.xlabel("Y")
plt.ylabel("Z")
plt.title("Initial Relaxed Positions (Aligned with Gravity)")
plt.grid(True)
plt.axis('equal')
plt.legend()
plt.tight_layout()
plt.show()
