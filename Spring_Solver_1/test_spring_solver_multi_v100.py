import ctypes
import os
import numpy as np
import matplotlib.pyplot as plt
from numpy.ctypeslib import ndpointer

# --- Step 1: Ensure mingw64/bin is in PATH ---
#mingw_bin = r"C:\msys64\mingw64\bin"
#if mingw_bin not in os.environ["PATH"]:
#    print(f"🛠 Adding {mingw_bin} to PATH")
##    os.environ["PATH"] = mingw_bin + ";" + os.environ["PATH"]
#else:
#    print(f"✅ {mingw_bin} is already in PATH")

# --- Step 2: Add your DLL folder to PATH ---
#dll_dir = r"C:\Users\Martin\Documents\MapleSim\Modelica\Modelica1\Resources\Library"
#if dll_dir not in os.environ["PATH"]:
#    print(f"🛠 Adding {dll_dir} to PATH")
#    os.environ["PATH"] = dll_dir + ";" + os.environ["PATH"]

# --- Step 3: Define and check DLL path ---
#dll_path = os.path.join(dll_dir, "spring_solver_lapak.dll")
#print(f"🔍 Loading DLL from: {dll_path}")
#assert os.path.exists(dll_path), "❌ DLL file does not exist!"



# Load the DLL
#dll_path = '.\\Resources\\Library\\spring_solver.dll'  # Change path as needed
dll_path = '.\\Resources\\Library\\spring_solver_lapak.dll'  # Change path as needed
#dll_path = '.\\Resources\\Library\\spring_solver_ref_ref.dll'  # Change path as needed
#dll_path = './Resources/Library/solve_spring_mass_v700.so'

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

def set_axes_equal(ax):
    """Set 3D plot axes to equal scale."""
    limits = np.array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        ax.get_zlim3d(),
    ])
    spans = limits[:, 1] - limits[:, 0]
    centers = np.mean(limits, axis=1)
    max_span = max(spans)

    new_limits = np.array([
        centers - max_span / 2,
        centers + max_span / 2
    ]).T

    ax.set_xlim3d(new_limits[0])
    ax.set_ylim3d(new_limits[1])
    ax.set_zlim3d(new_limits[2])

def run_and_plot(P1, P2, n, total_mass, length_factor, rope_diameter,
                  youngs_modulus,g_vec):
    P1 = np.array(P1, dtype=np.float64)
    P2 = np.array(P2, dtype=np.float64)
    out_positions = np.zeros((n - 1) * 3, dtype=np.float64)

    F_P1_n= np.zeros( 3, dtype=np.float64)
    F_P2_n= np.zeros( 3, dtype=np.float64)
    F_P1_w= np.zeros( 3, dtype=np.float64)
    F_P2_w= np.zeros( 3, dtype=np.float64)
    Length_init     = ctypes.c_double(0.0)
    Length_dynamic  = ctypes.c_double(0.0)
    Length_Newton   = ctypes.c_double(0.0)
    Length_energy   = ctypes.c_double(0.0)

    Status_Dynamic  = ctypes.c_int(0)
    Status_Newton   = ctypes.c_int(0)
    Status_Energy   = ctypes.c_int(0)

    g = 9.81
    g_vec = np.array([0.0, 0., -9.81], dtype=np.float64)

    #print(f"\nRunning solver with n={n}, length_factor={length_factor}, mass={total_mass}")
    res = lib.solve_spring_mass_c(
        P1, P2, n, total_mass, length_factor, rope_diameter, youngs_modulus, g_vec, out_positions,
        F_P1_n, F_P2_n, F_P1_w, F_P2_w,
        ctypes.byref(Length_init), ctypes.byref(Length_dynamic),
        ctypes.byref(Length_Newton), ctypes.byref(Length_energy),
        ctypes.byref(Status_Dynamic), ctypes.byref(Status_Newton), ctypes.byref(Status_Energy)
)
    F_P1_n= np.array(F_P1_n, dtype=np.float64)
    F_P2_n= np.array(F_P2_n, dtype=np.float64)
    F_P1_w= np.array(F_P1_w, dtype=np.float64)
    F_P2_w= np.array(F_P2_w, dtype=np.float64)

    print("\n")

    print("F P1 Newton = ",F_P1_n)   
    print("F P2 Newton = ",F_P2_n) 
    print("F P1 Weight = ",F_P1_w)  
    print("F P2 Weight = ",F_P2_w) 
    print(f"Length Init    = {Length_init.value:.6f}")
    print(f"Length Dynamic = {Length_dynamic.value:.6f}")
    print(f"Length Newton  = {Length_Newton.value:.6f}")
    print(f"Length Energy  = {Length_energy.value:.6f}")
    print(f"Status Dynamic = {Status_Dynamic.value}")
    print(f"Status Newton  = {Status_Newton.value}")
    print(f"Status Energy  = {Status_Energy.value}")

    if res != 0:
        print("❌ Solver failed")
        return

    positions = out_positions.reshape((n - 1, 3))
    all_points = np.vstack([P1, positions, P2])
    spring_lengths = np.linalg.norm(np.diff(all_points, axis=0), axis=1)

    fig = plt.figure(figsize=(12, 6))
    ax1 = fig.add_subplot(121, projection='3d')
    ax1.plot(all_points[:, 0], all_points[:, 1], all_points[:, 2], '-o', label='Springs')
    ax1.scatter(P1[0], P1[1], P1[2], color='black', label='Fixed P1')
    ax1.scatter(P2[0], P2[1], P2[2], color='black', label='Fixed P2')
    ax1.set_title("3D View of Spring-Mass System")
    ax1.set_xlabel("X")
    ax1.set_ylabel("Y")
    ax1.set_zlabel("Z")
    ax1.legend()
    ax1.view_init(elev=20, azim=-60)
    set_axes_equal(ax1)  # <-- Enforce equal scaling for all axes

    ax2 = fig.add_subplot(122)
    ax2.plot(spring_lengths, marker='o')
    ax2.set_title("Spring Lengths")
    ax2.set_xlabel("Spring Index")
    ax2.set_ylabel("Length [m]")
    plt.tight_layout()
    plt.show()

# Example test cases
run_and_plot(P1=[20, 0, 0], P2=[50, 70, 10], n=20,  total_mass=20, rope_diameter= 0.01, length_factor= 1.5,
             youngs_modulus= 1e10, g_vec=[0,0,-9.81])
