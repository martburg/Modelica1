import numpy as np
import matplotlib.pyplot as plt

# Given endpoints
P1 = np.array([0.0, 0.0, 0.0])
P2 = np.array([0.0, 0.167595, 0.985856])

# Gravity vector
g_vec = np.array([0.0, 0.0, -9.81])
e_g = g_vec / np.linalg.norm(g_vec)

# Number of internal nodes
n = 10
L = np.linalg.norm(P2 - P1)
sag_depth = 0.05 * L  # 5% of rope length

# Generate cubic parabola along gravity
points = []
for i in range(1, n):
    t = i / n
    base = (1 - t) * P1 + t * P2
    sag = -140.0 * sag_depth * t * (1.0 - t)
    pos = base - sag * e_g
    points.append(pos)

points = np.array(points)

# Add endpoints for full curve
curve = np.vstack([P1, points, P2])

# Plot the curve
fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111, projection='3d')
ax.plot(curve[:, 0], curve[:, 1], curve[:, 2], marker='o', label='Cubic parabola (gravity sag)')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Cubic Parabola from P1 to P2 with Sag along Gravity')
ax.legend()
ax.view_init(elev=20, azim=120)
plt.tight_layout()
plt.show()
