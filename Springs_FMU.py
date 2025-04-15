import numpy as np
from scipy.optimize import root
from fmpy import *
from fmpy.model_description import read_model_description


def compute_spring_force(p1, p2, s0, c):
    delta = p2 - p1
    length = np.linalg.norm(delta)
    force = - (length - s0) * c * delta / length  # Spring force
    return force, length  # Return the spring force and length for debugging

def compute_spring_lengths(positions, P1, P2):
    # Compute the spring lengths for all springs in the system
    all_points = [P1] + list(positions) + [P2]
    spring_lengths = []
    
    for i in range(len(all_points) - 1):
        length = np.linalg.norm(all_points[i + 1] - all_points[i])
        spring_lengths.append(length)
    
    return spring_lengths

def equilibrium_residuals(positions_flat, P1, P2, n, s0_list, total_mass, c_real, g_scaled):
    positions = positions_flat.reshape((n - 1, 3))
    all_points = [P1] + list(positions) + [P2]

    residuals = []
    m = total_mass / (n - 1)

    for i in range(1, n):
        left = all_points[i - 1]
        mid = all_points[i]
        right = all_points[i + 1]

        F_left, length_left = compute_spring_force(left, mid, s0_list[i - 1], c_real)
        F_right, length_right = compute_spring_force(right, mid, s0_list[i], c_real)
        F_gravity = np.array([0, 0, -m * g_scaled])
        net_force = F_left + F_right + F_gravity
        residuals.append(net_force)

    return np.concatenate(residuals)

def refined_initial_guess(P1, P2, n, max_sag=0.5):
    # More gradual initial guess: cubic interpolation between P1 and P2
    P1 = np.array(P1)
    P2 = np.array(P2)
    guesses = []
    for i in range(1, n):
        t = i / (n - 1)
        point = (1 - t) * P1 + t * P2
        # Gradually increase the sag with cubic function, more natural curve
        sag = -max_sag * (1 - 4 * t * (1 - t))  # Cubic curve sag
        point[2] += sag
        guesses.append(point)
    return np.array(guesses)

def solve_mass_positions(P1, P2, n, s0_list, c_real, total_mass_real, g=9.81, max_sag=0.5):
    L = np.linalg.norm(P2 - P1)
    P1_norm = np.array(P1) / L
    P2_norm = np.array(P2) / L
    s0_norm = [s / L for s in s0_list]
    g_scaled = g / (c_real * L)

    # Better initial guess with gradual sag
    guesses = refined_initial_guess(P1_norm, P2_norm, n, max_sag)
    positions_flat = guesses.flatten()

    result = root(
        equilibrium_residuals,
        positions_flat,
        args=(P1_norm, P2_norm, n, s0_norm, total_mass_real, c_real, g_scaled),
        method='hybr',
        options={'maxfev': 10000}
    )

    if not result.success:
        raise RuntimeError("Solver failed: " + result.message)

    positions_norm = result.x.reshape((n - 1, 3))
    return positions_norm * L

fmu_filename = 'Sag.fmu'
unzipdir = extract(fmu_filename)

model_description = read_model_description(fmu_filename)


# Filter parameters
#parameters = [v for v in model_description.modelVariables if v.causality== 'parameter']
#print(parameters)
# Print all parameter names and their types

target_vars = ['f1_r[1]','f1_r[2]', 'f1_r[3]',"f2_r[1]","f2_r[2]","f2_r[3]"]

vars = {}
for variable in model_description.modelVariables:
    if variable.name in target_vars and variable.start is not None:
        vars.
        print(f"{variable.name} = {variable.start}")

print(model_description.modelVariables.f1_r[1].start)       
#P1 = [v for v in parameters if v.name == "f1_r[1]"]
#P2 = [v for v in parameters if v.name == "f1_r[2]"]     
#print(P1,P2)


body_r_0 = positions = solve_mass_positions(P1, P2, n, s0_list, c, total_mass, max_sag=0.5)

body_r_0_flat = body_r_0.flatten()



model_description = read_model_description(fmu_filename)
model = instantiate_fmu(model_description, unzipdir, fmi_type='CoSimulation')

# Map parameter names
parameter_names = [f'body_r_0[{i+1},{j+1}]' for i in range(body_r_0.shape[0]) for j in range(3)]

# Set initial parameters
model.setupExperiment(startTime=0.0)
model.enterInitializationMode()
for name, value in zip(parameter_names, body_r_0_flat):
    model.setReal([model.getVariableByName(name).valueReference], [value])
model.exitInitializationMode()

# Simulate
results = simulate_fmu(fmu_filename, start_time=0.0, stop_time=10.0,
                       input=None, fmi_type='CoSimulation',
                       output=['bodies[1].r_0[1]', 'bodies[1].r_0[2]'],  # just an example
                       start_values={name: val for name, val in zip(parameter_names, body_r_0_flat)})
