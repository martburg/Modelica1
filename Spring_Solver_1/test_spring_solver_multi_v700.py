import ctypes
import os
import numpy as np
import matplotlib.pyplot as plt
from numpy.ctypeslib import ndpointer

dll_path = '.\\Resources\\Library\\spring_solver_lapak.dll'  # Change path as needed
lib = ctypes.CDLL(dll_path)

# Setup the prototype
lib.solve_spring_mass_c.argtypes = [
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # P1
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # P2
    ctypes.c_int,                                       # n
    ctypes.c_double,                                    # total_mass
    ctypes.c_double,                                    # length_factor
    ctypes.c_double,                                    # rope_diameter
    ctypes.c_double,                                    # youngs_modulus
    ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),   # g
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),   # out_positions
    ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),   # F_P1_n
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),   # F_P2_n
    ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),   # F_p1_w
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),   # F_P2_w
    ctypes.POINTER(ctypes.c_double),                    # Length Init
    ctypes.POINTER(ctypes.c_double),                    # Length Dynamic
    ctypes.POINTER(ctypes.c_double),                    # Length Newton
    ctypes.POINTER(ctypes.c_double),                    # Length_Energy
    ctypes.POINTER(ctypes.c_int),                       # Status Dynamic
    ctypes.POINTER(ctypes.c_int),                       # Status Newton
    ctypes.POINTER(ctypes.c_int),                       # Status Energy
]
lib.solve_spring_mass_c.restype = ctypes.c_int

def run_and_test(P1, P2, n, total_mass, length_factor, rope_diameter,
                 youngs_modulus, g_vec, tol_percent=2.5, verbose=False):
    P1 = np.array(P1, dtype=np.float64)
    P2 = np.array(P2, dtype=np.float64)
    out_positions = np.zeros((n - 1) * 3, dtype=np.float64)

    F_P1_n = np.zeros(3, dtype=np.float64)
    F_P2_n = np.zeros(3, dtype=np.float64)
    F_P1_w = np.zeros(3, dtype=np.float64)
    F_P2_w = np.zeros(3, dtype=np.float64)

    Length_init     = ctypes.c_double(0.0)
    Length_dynamic  = ctypes.c_double(0.0)
    Length_Newton   = ctypes.c_double(0.0)
    Length_energy   = ctypes.c_double(0.0)

    Status_Dynamic  = ctypes.c_int(0)
    Status_Newton   = ctypes.c_int(0)
    Status_Energy   = ctypes.c_int(0)

    g_vec = np.array(g_vec, dtype=np.float64)

    res = lib.solve_spring_mass_c(
        P1, P2, n, total_mass, length_factor, rope_diameter, youngs_modulus, g_vec, out_positions,
        F_P1_n, F_P2_n, F_P1_w, F_P2_w,
        ctypes.byref(Length_init), ctypes.byref(Length_dynamic),
        ctypes.byref(Length_Newton), ctypes.byref(Length_energy),
        ctypes.byref(Status_Dynamic), ctypes.byref(Status_Newton), ctypes.byref(Status_Energy)
    )
    print(Status_Energy.value)
    # Check solver statuses
    if Status_Dynamic.value != 0 or Status_Newton.value != 0 or Status_Energy.value not in [0, -6]:
        if verbose:
            print("❌ Solver status failed.")
            print(f"Status_Dynamic = {Status_Dynamic.value}, Status_Newton = {Status_Newton.value}, Status_Energy = {Status_Energy.value}")
        return False

    positions = out_positions.reshape((n - 1, 3))
    all_points = np.vstack([P1, positions, P2])
    spring_lengths = np.linalg.norm(np.diff(all_points, axis=0), axis=1)

    mean_length = np.mean(spring_lengths)
    deviations = np.abs(spring_lengths - mean_length) / mean_length * 100.0
    max_deviation = np.max(deviations)


    if res == -3:
        print("Implausible Density")
        return False
    if res == -4:
        print("Gravity to low")
        return False
    if res == -5:
        print("Gravity to high")
        return False

    if verbose:
        print(f"Max deviation = {max_deviation:.2f}%")

    if max_deviation > tol_percent:
        if verbose:
            print(f"❌ Spring lengths differ by more than {tol_percent}%")
        return False

    if verbose:
        print("✅ Test passed.")
    return True

test_cases = [

{
    "P1": [0.48, 94.23, -48.23],
    "P2": [0.91, 71.91, -10.35],
    "n": 93,
    "total_mass": 189.64,
    "length_factor": 2.47,
    "rope_diameter": 0.04715,
    "youngs_modulus": 37304794118.23,
    "g_vec": [-8.82, 14.22, -4.4]
},

{
    "P1": [-39.72, -93.79, -21.63],
    "P2": [-59.96, 36.01, -63.74],
    "n": 51,
    "total_mass": 269.55,
    "length_factor": 2.1,
    "rope_diameter": 0.03442,
    "youngs_modulus": 72266166781.55,
    "g_vec": [6.99, -9.46, -5.79]
},

{
    "P1": [70.52, 14.39, -57.23],
    "P2": [80.02, -31.47, -92.2],
    "n": 51,
    "total_mass": 228.06,
    "length_factor": 1.28,
    "rope_diameter": 0.0623,
    "youngs_modulus": 161390775428.56,
    "g_vec": [-0.22, -15.04, -10.46]
},

{
    "P1": [87.38, 21.31, 42.07],
    "P2": [42.72, 98.56, -31.69],
    "n": 41,
    "total_mass": 89.83,
    "length_factor": 1.96,
    "rope_diameter": 0.02245,
    "youngs_modulus": 185534403806.82,
    "g_vec": [-13.49, 0.77, 1.33]
},

{
    "P1": [-47.09, -92.09, -78.35],
    "P2": [-51.95, -45.22, -93.44],
    "n": 11,
    "total_mass": 83.34,
    "length_factor": 1.57,
    "rope_diameter": 0.03696,
    "youngs_modulus": 65591497467.31,
    "g_vec": [-8.83, -12.28, -4.86]
},

{
    "P1": [-82.05, -94.96, -59.58],
    "P2": [51.33, 98.16, 85.71],
    "n": 76,
    "total_mass": 1937.8,
    "length_factor": 2.17,
    "rope_diameter": 0.06418,
    "youngs_modulus": 179101960393.6,
    "g_vec": [0.17, 1.27, -7.42]
},

{
    "P1": [-13.13, 53.38, 61.28],
    "P2": [-38.31, 85.85, -81.13],
    "n": 85,
    "total_mass": 418.93,
    "length_factor": 1.15,
    "rope_diameter": 0.05594,
    "youngs_modulus": 154039688019.34,
    "g_vec": [19.77, 10.8, -16.27]
},

{
    "P1": [59.84, 20.48, 71.37],
    "P2": [-80.42, -44.82, -17.16],
    "n": 14,
    "total_mass": 334.53,
    "length_factor": 1.7,
    "rope_diameter": 0.03749,
    "youngs_modulus": 54075490695.82,
    "g_vec": [-18.3, -2.07, -9.39]
},

{
    "P1": [78.75, -14.81, -34.44],
    "P2": [0.5, -16.45, -48.53],
    "n": 89,
    "total_mass": 258.24,
    "length_factor": 1.68,
    "rope_diameter": 0.04961,
    "youngs_modulus": 164784233086.55,
    "g_vec": [0.11, -13.93, -5.92]
},

{
    "P1": [-48.66, 93.95, 56.94],
    "P2": [93.52, 67.63, -71.58],
    "n": 24,
    "total_mass": 578.43,
    "length_factor": 1.47,
    "rope_diameter": 0.05089,
    "youngs_modulus": 23046058802.32,
    "g_vec": [1.06, -13.2, 12.73]
},

{
    "P1": [12.77, 44.21, -88.23],
    "P2": [19.57, 63.18, -86.49],
    "n": 93,
    "total_mass": 12.84,
    "length_factor": 2.41,
    "rope_diameter": 0.01831,
    "youngs_modulus": 178946769057.13,
    "g_vec": [0.49, -9.66, 6.52]
},

{
    "P1": [-62.81, -35.35, 17.53],
    "P2": [43.61, -11.47, -68.74],
    "n": 79,
    "total_mass": 982.2,
    "length_factor": 1.78,
    "rope_diameter": 0.07108,
    "youngs_modulus": 17560998960.48,
    "g_vec": [-15.83, 19.99, -1.79]
},

{
    "P1": [-77.2, 85.93, 89.5],
    "P2": [-37.06, -28.97, -9.51],
    "n": 93,
    "total_mass": 198.03,
    "length_factor": 2.09,
    "rope_diameter": 0.02773,
    "youngs_modulus": 13440463990.22,
    "g_vec": [-3.35, -12.4, 15.94]
},

{
    "P1": [-92.25, 30.75, -32.2],
    "P2": [31.09, -81.2, 89.97],
    "n": 88,
    "total_mass": 276.07,
    "length_factor": 1.62,
    "rope_diameter": 0.03241,
    "youngs_modulus": 59380728416.63,
    "g_vec": [9.94, 15.05, -14.54]
},

{
    "P1": [-87.87, -26.22, -49.09],
    "P2": [-22.82, -60.61, -47.62],
    "n": 58,
    "total_mass": 282.77,
    "length_factor": 2.18,
    "rope_diameter": 0.04737,
    "youngs_modulus": 75779802590.03,
    "g_vec": [5.08, -4.47, -6.18]
},

{
    "P1": [17.5, 93.24, 4.4],
    "P2": [81.96, 58.73, -88.26],
    "n": 49,
    "total_mass": 201.99,
    "length_factor": 1.23,
    "rope_diameter": 0.04209,
    "youngs_modulus": 25321492786.87,
    "g_vec": [-2.84, 7.93, 19.47]
},

{
    "P1": [50.3, 1.54, 44.69],
    "P2": [73.87, 25.79, 22.49],
    "n": 60,
    "total_mass": 9.41,
    "length_factor": 1.88,
    "rope_diameter": 0.01255,
    "youngs_modulus": 43261699266.41,
    "g_vec": [10.18, 13.7, 14.76]
},

{
    "P1": [72.07, 29.68, -15.96],
    "P2": [-62.38, -32.2, 40.76],
    "n": 34,
    "total_mass": 692.47,
    "length_factor": 1.38,
    "rope_diameter": 0.06349,
    "youngs_modulus": 194414314724.41,
    "g_vec": [-11.83, -16.26, 15.52]
},

{
    "P1": [41.64, 22.1, -38.47],
    "P2": [-40.17, -92.1, 45.7],
    "n": 36,
    "total_mass": 1622.49,
    "length_factor": 2.16,
    "rope_diameter": 0.07642,
    "youngs_modulus": 94333321722.97,
    "g_vec": [-18.27, -19.81, 9.18]
},

{
    "P1": [-2.79, -12.77, -35.07],
    "P2": [4.34, -13.77, 19.63],
    "n": 78,
    "total_mass": 79.79,
    "length_factor": 1.35,
    "rope_diameter": 0.03693,
    "youngs_modulus": 57501728979.5,
    "g_vec": [15.92, 19.8, -16.85]
},

{
    "P1": [-65.09, 52.82, -0.74],
    "P2": [83.85, -93.48, -86.6],
    "n": 41,
    "total_mass": 1538.82,
    "length_factor": 1.39,
    "rope_diameter": 0.07902,
    "youngs_modulus": 39631567105.0,
    "g_vec": [-3.87, -16.79, -9.49]
},

{
    "P1": [14.73, -0.4, 16.09],
    "P2": [-24.65, 65.35, -72.28],
    "n": 68,
    "total_mass": 223.42,
    "length_factor": 1.89,
    "rope_diameter": 0.03587,
    "youngs_modulus": 34072114431.48,
    "g_vec": [-2.86, 7.59, 14.73]
},

{
    "P1": [-86.17, 45.58, -42.39],
    "P2": [44.75, -86.72, 72.51],
    "n": 16,
    "total_mass": 76.78,
    "length_factor": 2.3,
    "rope_diameter": 0.01394,
    "youngs_modulus": 46948350775.49,
    "g_vec": [13.42, -16.44, -19.05]
},

{
    "P1": [9.65, -14.75, 59.11],
    "P2": [-2.22, 0.13, 97.82],
    "n": 20,
    "total_mass": 149.62,
    "length_factor": 2.34,
    "rope_diameter": 0.04344,
    "youngs_modulus": 137046071727.03,
    "g_vec": [-13.99, 13.23, -17.93]
},

{
    "P1": [-47.34, -16.29, -26.68],
    "P2": [4.29, -7.55, -34.62],
    "n": 69,
    "total_mass": 40.55,
    "length_factor": 2.32,
    "rope_diameter": 0.0205,
    "youngs_modulus": 120427857172.77,
    "g_vec": [-7.78, -9.97, -5.4]
},

{
    "P1": [-73.24, 38.77, 69.53],
    "P2": [-52.49, -69.88, -98.28],
    "n": 19,
    "total_mass": 58.6,
    "length_factor": 2.15,
    "rope_diameter": 0.01314,
    "youngs_modulus": 118430643652.48,
    "g_vec": [-18.11, 10.49, 19.49]
},

{
    "P1": [51.22, -28.16, -97.81],
    "P2": [19.04, 98.25, -48.48],
    "n": 16,
    "total_mass": 648.25,
    "length_factor": 2.21,
    "rope_diameter": 0.05175,
    "youngs_modulus": 150617123991.37,
    "g_vec": [-12.42, -15.04, 2.21]
},

{
    "P1": [78.51, 61.15, -50.55],
    "P2": [11.18, -40.52, -34.04],
    "n": 92,
    "total_mass": 274.47,
    "length_factor": 1.85,
    "rope_diameter": 0.03918,
    "youngs_modulus": 21466338634.17,
    "g_vec": [17.28, 7.28, 18.98]
},

{
    "P1": [-52.66, 5.53, -86.53],
    "P2": [86.29, -75.23, -60.4],
    "n": 53,
    "total_mass": 567.28,
    "length_factor": 1.2,
    "rope_diameter": 0.0608,
    "youngs_modulus": 138372493786.52,
    "g_vec": [10.61, -16.12, 16.01]
},

{
    "P1": [26.33, 51.47, 65.64],
    "P2": [43.67, 17.24, -35.51],
    "n": 70,
    "total_mass": 463.81,
    "length_factor": 1.71,
    "rope_diameter": 0.0565,
    "youngs_modulus": 199533501793.8,
    "g_vec": [7.89, 18.28, 3.5]
},

{
    "P1": [-67.92, -59.52, 7.89],
    "P2": [-3.72, 14.31, -90.19],
    "n": 78,
    "total_mass": 720.2,
    "length_factor": 2.11,
    "rope_diameter": 0.05601,
    "youngs_modulus": 128313025196.14,
    "g_vec": [-3.83, 3.32, 3.41]
},

{
    "P1": [58.05, 24.95, 7.5],
    "P2": [21.24, -79.76, -45.37],
    "n": 26,
    "total_mass": 49.7,
    "length_factor": 1.93,
    "rope_diameter": 0.01633,
    "youngs_modulus": 146805971792.15,
    "g_vec": [1.1, 14.92, 7.16]
},

{
    "P1": [65.49, 18.69, 18.27],
    "P2": [-32.28, -10.75, 72.82],
    "n": 29,
    "total_mass": 21.64,
    "length_factor": 1.3,
    "rope_diameter": 0.01353,
    "youngs_modulus": 38613107655.98,
    "g_vec": [8.93, 0.79, -5.27]
},

{
    "P1": [-25.32, 88.83, 1.03],
    "P2": [-53.03, -78.85, -32.34],
    "n": 93,
    "total_mass": 750.91,
    "length_factor": 1.26,
    "rope_diameter": 0.06619,
    "youngs_modulus": 191819621713.59,
    "g_vec": [12.63, -4.67, 10.9]
},

{
    "P1": [-19.57, -89.81, 84.66],
    "P2": [76.73, 4.43, 10.45],
    "n": 64,
    "total_mass": 38.08,
    "length_factor": 2.06,
    "rope_diameter": 0.01237,
    "youngs_modulus": 190099748997.84,
    "g_vec": [-11.53, 12.78, -0.06]
},

{
    "P1": [65.83, -44.56, -89.33],
    "P2": [88.07, 87.41, 6.48],
    "n": 73,
    "total_mass": 54.04,
    "length_factor": 1.24,
    "rope_diameter": 0.01836,
    "youngs_modulus": 167329782873.84,
    "g_vec": [5.69, -1.5, -3.21]
},

{
    "P1": [-70.83, -54.66, 9.72],
    "P2": [-56.97, -67.36, 12.66],
    "n": 30,
    "total_mass": 15.05,
    "length_factor": 2.38,
    "rope_diameter": 0.02057,
    "youngs_modulus": 180498968584.91,
    "g_vec": [-5.28, 5.05, 16.49]
},

{
    "P1": [65.82, 26.26, 24.03],
    "P2": [27.84, 38.1, -96.7],
    "n": 82,
    "total_mass": 148.24,
    "length_factor": 1.12,
    "rope_diameter": 0.03641,
    "youngs_modulus": 101442951961.99,
    "g_vec": [18.32, 3.42, 11.26]
},

{
    "P1": [-12.6, -32.06, -50.03],
    "P2": [5.51, -24.94, 93.99],
    "n": 55,
    "total_mass": 473.96,
    "length_factor": 1.8,
    "rope_diameter": 0.04803,
    "youngs_modulus": 128751655680.2,
    "g_vec": [-13.08, -7.21, -8.24]
},

{
    "P1": [90.0, -42.96, 38.55],
    "P2": [-94.99, -47.52, -7.26],
    "n": 70,
    "total_mass": 1056.63,
    "length_factor": 1.82,
    "rope_diameter": 0.06227,
    "youngs_modulus": 79219089603.18,
    "g_vec": [-18.76, 4.71, 7.13]
},

{
    "P1": [-91.75, -80.01, 31.0],
    "P2": [99.11, -96.84, 11.52],
    "n": 64,
    "total_mass": 896.82,
    "length_factor": 2.08,
    "rope_diameter": 0.05339,
    "youngs_modulus": 171794660276.19,
    "g_vec": [10.04, 13.07, 1.81]
},

{
    "P1": [-75.54, 85.24, -94.8],
    "P2": [67.24, -84.65, -61.07],
    "n": 17,
    "total_mass": 142.56,
    "length_factor": 1.73,
    "rope_diameter": 0.02162,
    "youngs_modulus": 153534529929.42,
    "g_vec": [-19.84, -2.19, -6.63]
},

{
    "P1": [78.63, 30.51, -75.12],
    "P2": [-40.92, 82.66, -82.25],
    "n": 46,
    "total_mass": 313.2,
    "length_factor": 2.03,
    "rope_diameter": 0.03878,
    "youngs_modulus": 91354552214.82,
    "g_vec": [-3.73, 17.7, 18.03]
},

{
    "P1": [50.69, -93.68, -71.17],
    "P2": [22.84, -53.49, 33.21],
    "n": 17,
    "total_mass": 390.81,
    "length_factor": 2.18,
    "rope_diameter": 0.0445,
    "youngs_modulus": 98857012764.6,
    "g_vec": [-5.95, 15.47, 12.77]
},

{
    "P1": [-60.5, 72.18, 64.61],
    "P2": [80.43, 75.5, 91.92],
    "n": 34,
    "total_mass": 331.98,
    "length_factor": 1.88,
    "rope_diameter": 0.03957,
    "youngs_modulus": 68527045306.49,
    "g_vec": [-10.87, 17.27, 2.85]
},

{
    "P1": [-42.29, -51.74, -83.51],
    "P2": [-0.69, 77.55, 26.03],
    "n": 43,
    "total_mass": 187.36,
    "length_factor": 2.25,
    "rope_diameter": 0.02465,
    "youngs_modulus": 143327385117.9,
    "g_vec": [0.84, -2.86, -1.42]
},

{
    "P1": [52.0, 14.97, 72.59],
    "P2": [-41.79, 76.32, 45.88],
    "n": 27,
    "total_mass": 966.97,
    "length_factor": 1.69,
    "rope_diameter": 0.07952,
    "youngs_modulus": 162784997585.08,
    "g_vec": [11.95, 6.58, 13.85]
},

{
    "P1": [48.92, -9.43, 85.39],
    "P2": [-21.11, 13.09, -78.81],
    "n": 26,
    "total_mass": 287.54,
    "length_factor": 1.44,
    "rope_diameter": 0.03759,
    "youngs_modulus": 161165334100.51,
    "g_vec": [-10.97, 13.63, -4.14]
},

{
    "P1": [79.03, -85.98, 14.96],
    "P2": [-6.55, 22.86, -10.54],
    "n": 54,
    "total_mass": 395.29,
    "length_factor": 1.76,
    "rope_diameter": 0.04507,
    "youngs_modulus": 163568883638.52,
    "g_vec": [-4.05, -1.12, 3.28]
},

{
    "P1": [49.13, -31.72, 60.12],
    "P2": [70.73, -69.3, 16.54],
    "n": 53,
    "total_mass": 228.99,
    "length_factor": 1.46,
    "rope_diameter": 0.057,
    "youngs_modulus": 189729807118.56,
    "g_vec": [-8.59, 9.57, -10.45]
},

{
    "P1": [-23.74, 75.15, 97.39],
    "P2": [74.97, -66.83, 63.93],
    "n": 56,
    "total_mass": 867.11,
    "length_factor": 1.85,
    "rope_diameter": 0.05821,
    "youngs_modulus": 77578926376.27,
    "g_vec": [-18.71, -4.27, 14.83]
},

{
    "P1": [62.46, 41.91, -6.32],
    "P2": [78.5, -54.01, -65.92],
    "n": 12,
    "total_mass": 87.95,
    "length_factor": 1.3,
    "rope_diameter": 0.02748,
    "youngs_modulus": 45287244106.86,
    "g_vec": [5.71, 18.88, -10.11]
},

{
    "P1": [0.64, 98.17, -97.57],
    "P2": [11.3, 2.46, 61.49],
    "n": 87,
    "total_mass": 1229.47,
    "length_factor": 2.2,
    "rope_diameter": 0.06186,
    "youngs_modulus": 27341297522.56,
    "g_vec": [-5.79, -3.59, -13.28]
},

{
    "P1": [-21.31, 78.58, -99.46],
    "P2": [-40.25, 11.01, -33.18],
    "n": 18,
    "total_mass": 351.05,
    "length_factor": 2.09,
    "rope_diameter": 0.04707,
    "youngs_modulus": 18630935078.4,
    "g_vec": [-9.12, -1.6, -12.74]
},

{
    "P1": [68.88, -65.0, 52.96],
    "P2": [-27.89, 56.55, -58.75],
    "n": 61,
    "total_mass": 210.6,
    "length_factor": 2.11,
    "rope_diameter": 0.02577,
    "youngs_modulus": 167478321826.3,
    "g_vec": [-12.53, -0.12, -18.68]
},

{
    "P1": [-83.84, 76.42, 59.5],
    "P2": [73.94, 39.69, 34.83],
    "n": 12,
    "total_mass": 57.55,
    "length_factor": 1.54,
    "rope_diameter": 0.01704,
    "youngs_modulus": 188961586337.09,
    "g_vec": [1.0, 6.11, 16.76]
},

{
    "P1": [-10.95, 35.07, 68.22],
    "P2": [98.2, -67.64, 59.9],
    "n": 16,
    "total_mass": 287.47,
    "length_factor": 2.01,
    "rope_diameter": 0.03483,
    "youngs_modulus": 102814498725.81,
    "g_vec": [15.5, -18.04, 10.16]
},

{
    "P1": [-48.05, -72.96, -76.61],
    "P2": [-79.34, 27.56, -29.72],
    "n": 50,
    "total_mass": 497.14,
    "length_factor": 2.41,
    "rope_diameter": 0.04774,
    "youngs_modulus": 167174235624.15,
    "g_vec": [-11.33, 7.68, 7.47]
},

{
    "P1": [-84.97, 47.97, 97.69],
    "P2": [69.32, 15.93, 25.35],
    "n": 44,
    "total_mass": 1555.42,
    "length_factor": 1.79,
    "rope_diameter": 0.07988,
    "youngs_modulus": 12270200114.03,
    "g_vec": [-5.87, 7.18, 10.82]
},

{
    "P1": [19.17, -76.2, -51.05],
    "P2": [-34.24, 46.21, -63.24],
    "n": 49,
    "total_mass": 418.84,
    "length_factor": 1.16,
    "rope_diameter": 0.05855,
    "youngs_modulus": 111739272552.05,
    "g_vec": [-13.37, -18.2, 14.93]
},

{
    "P1": [-43.49, 19.35, -13.08],
    "P2": [84.66, -42.44, 35.83],
    "n": 12,
    "total_mass": 788.19,
    "length_factor": 1.75,
    "rope_diameter": 0.06174,
    "youngs_modulus": 44006487931.06,
    "g_vec": [-17.84, 18.23, 19.23]
},

{
    "P1": [62.5, 47.27, -55.61],
    "P2": [-9.76, -7.44, 65.28],
    "n": 65,
    "total_mass": 160.9,
    "length_factor": 1.62,
    "rope_diameter": 0.02893,
    "youngs_modulus": 157390357751.07,
    "g_vec": [7.35, 5.73, 5.99]
},

{
    "P1": [-13.4, -12.35, -50.93],
    "P2": [-58.11, 17.56, -53.49],
    "n": 14,
    "total_mass": 96.44,
    "length_factor": 1.7,
    "rope_diameter": 0.03662,
    "youngs_modulus": 163010814233.38,
    "g_vec": [-17.43, 10.7, -4.94]
},

{
    "P1": [-35.83, 51.61, 87.91],
    "P2": [32.17, -13.66, 49.56],
    "n": 47,
    "total_mass": 103.84,
    "length_factor": 1.25,
    "rope_diameter": 0.03224,
    "youngs_modulus": 24780137657.05,
    "g_vec": [-1.68, -3.42, 1.01]
},

{
    "P1": [-17.91, 11.39, 86.14],
    "P2": [-92.76, 24.29, -8.45],
    "n": 47,
    "total_mass": 89.41,
    "length_factor": 1.48,
    "rope_diameter": 0.02518,
    "youngs_modulus": 80267254013.09,
    "g_vec": [-19.19, -10.02, -15.26]
},

{
    "P1": [5.81, -99.8, -67.63],
    "P2": [-81.29, 9.93, -34.19],
    "n": 73,
    "total_mass": 814.61,
    "length_factor": 1.25,
    "rope_diameter": 0.0759,
    "youngs_modulus": 29317899833.76,
    "g_vec": [1.94, -11.13, -19.03]
},

{
    "P1": [96.26, 23.53, 22.06],
    "P2": [97.34, -82.31, 94.9],
    "n": 21,
    "total_mass": 329.43,
    "length_factor": 2.4,
    "rope_diameter": 0.03688,
    "youngs_modulus": 29854053410.77,
    "g_vec": [-11.84, 4.92, -3.2]
},

{
    "P1": [99.41, -24.17, 29.75],
    "P2": [93.58, -88.48, 84.57],
    "n": 20,
    "total_mass": 304.13,
    "length_factor": 2.04,
    "rope_diameter": 0.04734,
    "youngs_modulus": 31005485769.43,
    "g_vec": [2.71, 8.18, -14.03]
},

{
    "P1": [-51.18, 70.88, -93.11],
    "P2": [41.77, 99.88, -98.5],
    "n": 51,
    "total_mass": 415.05,
    "length_factor": 2.06,
    "rope_diameter": 0.05129,
    "youngs_modulus": 117633388159.33,
    "g_vec": [-13.56, 6.97, -17.6]
},

{
    "P1": [-26.96, 39.46, 18.86],
    "P2": [-58.47, 12.63, 76.38],
    "n": 69,
    "total_mass": 20.54,
    "length_factor": 1.44,
    "rope_diameter": 0.01601,
    "youngs_modulus": 197325317245.82,
    "g_vec": [6.67, -1.13, -7.14]
},

{
    "P1": [-55.5, 62.98, -66.27],
    "P2": [-43.28, -16.34, -35.12],
    "n": 11,
    "total_mass": 63.76,
    "length_factor": 2.19,
    "rope_diameter": 0.02075,
    "youngs_modulus": 29597903191.81,
    "g_vec": [11.93, 18.31, 15.01]
},

{
    "P1": [82.68, 3.58, -32.72],
    "P2": [-68.59, -97.64, -50.41],
    "n": 20,
    "total_mass": 80.58,
    "length_factor": 1.93,
    "rope_diameter": 0.01705,
    "youngs_modulus": 108592911737.42,
    "g_vec": [16.56, 11.89, 18.79]
},

{
    "P1": [-35.91, -43.31, 14.77],
    "P2": [-63.34, -6.17, 30.16],
    "n": 33,
    "total_mass": 350.41,
    "length_factor": 2.02,
    "rope_diameter": 0.06737,
    "youngs_modulus": 142003088235.1,
    "g_vec": [-14.87, 17.22, 0.23]
},

{
    "P1": [43.23, 75.02, 22.72],
    "P2": [-44.36, -80.1, 21.25],
    "n": 28,
    "total_mass": 432.99,
    "length_factor": 1.32,
    "rope_diameter": 0.04842,
    "youngs_modulus": 158175117521.34,
    "g_vec": [-17.35, -10.35, -17.57]
},

{
    "P1": [64.38, -38.09, 9.79],
    "P2": [12.6, -5.42, 31.02],
    "n": 64,
    "total_mass": 454.84,
    "length_factor": 1.93,
    "rope_diameter": 0.06805,
    "youngs_modulus": 48978500891.27,
    "g_vec": [-11.6, 19.17, -16.08]
},

{
    "P1": [49.14, -83.1, 54.88],
    "P2": [18.0, -77.85, -84.1],
    "n": 57,
    "total_mass": 934.41,
    "length_factor": 2.08,
    "rope_diameter": 0.06335,
    "youngs_modulus": 141791814359.13,
    "g_vec": [10.03, -8.28, 15.57]
},

{
    "P1": [89.63, -81.69, -34.77],
    "P2": [-53.07, -21.35, 48.47],
    "n": 35,
    "total_mass": 80.33,
    "length_factor": 1.47,
    "rope_diameter": 0.01989,
    "youngs_modulus": 171937463445.38,
    "g_vec": [-1.03, 2.83, -7.43]
},

{
    "P1": [20.48, -24.4, 57.98],
    "P2": [-87.11, -49.22, -53.98],
    "n": 61,
    "total_mass": 17.99,
    "length_factor": 1.36,
    "rope_diameter": 0.01035,
    "youngs_modulus": 75950068103.07,
    "g_vec": [19.41, 8.24, 4.71]
},

{
    "P1": [84.86, 52.94, -10.52],
    "P2": [49.9, 95.08, -7.8],
    "n": 67,
    "total_mass": 277.2,
    "length_factor": 1.24,
    "rope_diameter": 0.07206,
    "youngs_modulus": 129693149635.82,
    "g_vec": [5.7, 16.75, -10.05]
},

{
    "P1": [56.57, 22.92, -0.89],
    "P2": [23.47, -69.67, -56.89],
    "n": 38,
    "total_mass": 152.84,
    "length_factor": 2.34,
    "rope_diameter": 0.02711,
    "youngs_modulus": 89211883528.79,
    "g_vec": [-9.2, 3.91, -2.32]
},

{
    "P1": [-18.06, -63.61, 46.33],
    "P2": [-10.31, -96.75, 92.8],
    "n": 43,
    "total_mass": 383.52,
    "length_factor": 1.57,
    "rope_diameter": 0.07348,
    "youngs_modulus": 122195626096.56,
    "g_vec": [12.63, 8.94, 14.76]
},

{
    "P1": [21.76, -46.54, -55.61],
    "P2": [-77.22, 31.86, 56.18],
    "n": 81,
    "total_mass": 431.74,
    "length_factor": 1.87,
    "rope_diameter": 0.04175,
    "youngs_modulus": 174965882623.27,
    "g_vec": [10.44, -7.35, -11.79]
},

{
    "P1": [83.21, 70.05, 92.92],
    "P2": [-6.9, -90.51, -83.29],
    "n": 92,
    "total_mass": 1614.01,
    "length_factor": 1.37,
    "rope_diameter": 0.07672,
    "youngs_modulus": 191876260670.32,
    "g_vec": [-4.39, -9.69, -13.2]
},

{
    "P1": [52.17, 28.19, 35.2],
    "P2": [-33.71, 29.29, -62.66],
    "n": 10,
    "total_mass": 738.2,
    "length_factor": 2.02,
    "rope_diameter": 0.05978,
    "youngs_modulus": 126492922475.81,
    "g_vec": [17.51, -7.93, -9.79]
},

{
    "P1": [-68.1, 34.67, -45.88],
    "P2": [-45.43, -67.56, -33.52],
    "n": 50,
    "total_mass": 85.14,
    "length_factor": 1.86,
    "rope_diameter": 0.02351,
    "youngs_modulus": 25784192160.05,
    "g_vec": [-6.75, -1.83, -12.38]
},

{
    "P1": [37.78, 88.79, -14.64],
    "P2": [80.22, -80.34, -91.41],
    "n": 13,
    "total_mass": 64.05,
    "length_factor": 2.19,
    "rope_diameter": 0.01398,
    "youngs_modulus": 60674493160.11,
    "g_vec": [8.85, 10.44, 13.72]
},

{
    "P1": [-65.4, 99.2, -87.27],
    "P2": [65.4, 27.9, 47.02],
    "n": 43,
    "total_mass": 30.93,
    "length_factor": 1.94,
    "rope_diameter": 0.01006,
    "youngs_modulus": 70778189278.27,
    "g_vec": [8.47, 14.34, -12.65]
},

{
    "P1": [-75.24, -62.77, -33.65],
    "P2": [-85.95, -88.62, 83.21],
    "n": 16,
    "total_mass": 79.11,
    "length_factor": 2.1,
    "rope_diameter": 0.01998,
    "youngs_modulus": 138132508150.3,
    "g_vec": [-7.52, 19.29, 18.57]
},

{
    "P1": [-80.76, -62.29, -84.67],
    "P2": [20.83, 97.32, 34.4],
    "n": 64,
    "total_mass": 66.55,
    "length_factor": 2.31,
    "rope_diameter": 0.01281,
    "youngs_modulus": 72942011292.86,
    "g_vec": [8.43, 19.38, 7.12]
},

{
    "P1": [-99.13, 1.45, 87.4],
    "P2": [-48.53, -1.7, -65.24],
    "n": 61,
    "total_mass": 816.01,
    "length_factor": 1.3,
    "rope_diameter": 0.07049,
    "youngs_modulus": 39013943234.29,
    "g_vec": [-17.69, -16.78, 9.08]
},

{
    "P1": [-94.63, 10.69, 39.98],
    "P2": [29.58, -22.39, 11.0],
    "n": 35,
    "total_mass": 991.32,
    "length_factor": 2.39,
    "rope_diameter": 0.06331,
    "youngs_modulus": 183336890097.03,
    "g_vec": [-16.24, -13.83, 17.27]
},

{
    "P1": [29.88, 85.65, -16.66],
    "P2": [-99.28, -87.49, -78.1],
    "n": 70,
    "total_mass": 1453.75,
    "length_factor": 1.93,
    "rope_diameter": 0.06535,
    "youngs_modulus": 149582841413.61,
    "g_vec": [-16.19, -5.93, 8.05]
},

{
    "P1": [-29.04, -97.4, -43.88],
    "P2": [-33.88, -45.15, 62.37],
    "n": 23,
    "total_mass": 845.05,
    "length_factor": 2.05,
    "rope_diameter": 0.06655,
    "youngs_modulus": 38975588919.41,
    "g_vec": [10.4, 7.82, -14.44]
},

{
    "P1": [29.62, -3.3, -12.13],
    "P2": [-77.58, 3.66, -80.01],
    "n": 74,
    "total_mass": 143.14,
    "length_factor": 1.6,
    "rope_diameter": 0.02994,
    "youngs_modulus": 179651106531.18,
    "g_vec": [-13.43, -14.7, 18.16]
},

{
    "P1": [12.59, 83.83, -40.02],
    "P2": [-77.03, 11.24, -30.89],
    "n": 38,
    "total_mass": 1033.79,
    "length_factor": 1.87,
    "rope_diameter": 0.078,
    "youngs_modulus": 54191369473.33,
    "g_vec": [12.08, 15.18, -1.1]
},

{
    "P1": [-2.2, -92.54, 60.8],
    "P2": [-79.71, -38.25, 13.06],
    "n": 85,
    "total_mass": 33.24,
    "length_factor": 1.77,
    "rope_diameter": 0.01502,
    "youngs_modulus": 45623600797.2,
    "g_vec": [-18.28, -16.83, 7.84]
},

{
    "P1": [-62.24, 52.68, 90.31],
    "P2": [-71.96, 91.94, 18.95],
    "n": 40,
    "total_mass": 520.73,
    "length_factor": 1.78,
    "rope_diameter": 0.06739,
    "youngs_modulus": 106676383917.07,
    "g_vec": [5.67, 8.4, 8.16]
},

{
    "P1": [-35.23, -78.95, 32.93],
    "P2": [-93.28, -44.87, -81.15],
    "n": 21,
    "total_mass": 199.34,
    "length_factor": 1.23,
    "rope_diameter": 0.03947,
    "youngs_modulus": 70555446922.94,
    "g_vec": [-12.37, 5.55, 5.14]
},

{
    "P1": [44.09, 15.81, -1.69],
    "P2": [-2.6, 93.26, 51.26],
    "n": 49,
    "total_mass": 736.51,
    "length_factor": 1.55,
    "rope_diameter": 0.07598,
    "youngs_modulus": 51195256985.33,
    "g_vec": [-0.12, -0.83, -18.94]
},

{
    "P1": [99.88, 43.74, -90.53],
    "P2": [-87.62, -94.9, 19.92],
    "n": 49,
    "total_mass": 923.56,
    "length_factor": 1.96,
    "rope_diameter": 0.04822,
    "youngs_modulus": 12803868756.23,
    "g_vec": [-17.14, 17.34, 17.48]
}
]

failed_test_logs = []
for i, case in enumerate(test_cases):
    print(f"\n--- Running Test Case {i + 1} ---")
    result = run_and_test(**case)
    #print(f"Result: {'PASS' if result else 'FAIL'}")
    if not result:
        print("❌ %d Test failed with input:"%i)
        for key, value in case.items():
            print(f"    {key} = {value}")
