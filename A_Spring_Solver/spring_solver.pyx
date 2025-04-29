#spring_solver.pyx

from libc.stdlib cimport malloc, free
import numpy as np
cimport numpy as np
from scipy.optimize import root
cimport cython
from libc.stdio cimport printf
from numpy cimport import_array
from numpy cimport npy_intp
from numpy cimport PyArray_DATA, PyArrayObject

cdef extern from "Python.h":
    int Py_IsInitialized()
    void Py_Initialize()

cdef int import_array_wrapper():
    printf(">>> [DEBUG] import_array_wrapper called\n") 
    if import_array() < 0:
        printf(">>> [ERROR] import_array failed\n")
        return -1
    printf(">>> [DEBUG] import array done\n") 
    return 0

# Expose C entry point for Modelica
cdef public int solve_spring_mass_c(double* P1, double* P2, int n,
                                    double* s0_list, double c,
                                    double total_mass, double g,
                                    double max_sag,
                                    double* out_positions):
    printf(">>> [DEBUG] Starting solve_spring_mass_c\n")
    printf("    n = %d, c = %f, total_mass = %f, g = %f, max_sag = %f\n", n, c, total_mass, g, max_sag)

    if n < 2:
        printf(">>> [ERROR] Invalid value of n: must be >= 2\n")
        return 2

    if not Py_IsInitialized():
        printf(">>> [DEBUG] Python not initialized. Initializing...\n")
        Py_Initialize()
    else:
        printf(">>> [DEBUG] Python already initialized.\n")

    if import_array_wrapper() < 0:
        printf(">>> [ERROR] NumPy import_array() failed\n")
        return 1

    printf(">>> [DEBUG] NumPy array API initialized\n")

    # Proper dimension declarations

    #cdef npy_intp* dim3 = <npy_intp*> malloc(sizeof(npy_intp))
    #dim3[0] = 3
    #cdef np.ndarray[np.float64_t, ndim=1] np_P1
    #np_P1 = np.PyArray_SimpleNewFromData(1, dim3, np.NPY_FLOAT64, <void*>P1)
    #free(dim3)
    #cdef npy_intp* dim_3 = <npy_intp*> malloc(sizeof(npy_intp))
    #dim_3[0] = 3
    #cdef np.ndarray[np.float64_t, ndim=1] np_P2
    #np_P2 = np.PyArray_SimpleNewFromData(1, dim_3, np.NPY_FLOAT64, <void*>P2)
    #free(dim_3)
    #cdef npy_intp* dim__3 = <npy_intp*> malloc(sizeof(npy_intp))
    #dim__3[0] = n
    #cdef np.ndarray[np.float64_t, ndim=1] np_s0
    #np_s0 = np.PyArray_SimpleNewFromData(1, dim__3, np.NPY_FLOAT64, <void*>s0_list)
    #free(dim__3)
    #cdef np.ndarray[np.float64_t, ndim=2] np_out = np.empty((n - 1, 3), dtype=np.float64)

    printf(">>> [DEBUG] Input arrays created\n")
    #cdef double* pP1 = <double*> PyArray_DATA(np_P1)
    #printf("    P1 = [%f, %f, %f]\n", pP1[0], pP1[1], pP1[2])

    #cdef double* pP2 = <double*> PyArray_DATA(np_P2)
    #cdef double* ps0 = <double*> PyArray_DATA(np_s0)

    #printf("    P2 = [%f, %f, %f]\n", pP2[0], pP2[1], pP2[2])
    #printf("    s0[0] = %f, s0[%d] = %f\n", ps0[0], n - 1, ps0[n - 1])

    printf(">>> [DEBUG] Ready to call solve_spring_mass\n")
    #solve_spring_mass(np_P1, np_P2, n, np_s0, c, total_mass, g, max_sag, np_out)

    printf(">>> [DEBUG] Copying result back to out_positions\n")
    #cdef int i, j
    #for i in range(n - 1):
    #    for j in range(3):
    #        out_positions[i * 3 + j] = np_out[i, j]

    printf(">>> [DEBUG] Finished successfully\n")
    return 0


@cython.boundscheck(False)
@cython.wraparound(False)
cdef void compute_spring_force(np.ndarray[np.float64_t, ndim=1] p1,
                                np.ndarray[np.float64_t, ndim=1] p2,
                                double s0, double c,
                                np.ndarray[np.float64_t, ndim=1] force_out,
                                double* length_out):
    cdef np.ndarray[np.float64_t, ndim=1] delta = p2 - p1
    cdef double length = np.linalg.norm(delta)
    force_out[:] = - (length - s0) * c * delta / length
    length_out[0] = length

@cython.boundscheck(False)
@cython.wraparound(False)
def equilibrium_residuals(np.ndarray[np.float64_t, ndim=1] positions_flat,
                          np.ndarray[np.float64_t, ndim=1] P1,
                          np.ndarray[np.float64_t, ndim=1] P2,
                          int n,
                          np.ndarray[np.float64_t, ndim=1] s0_list,
                          double total_mass,
                          double c_real,
                          double g_scaled):
    cdef np.ndarray[np.float64_t, ndim=2] positions = positions_flat.reshape((n - 1, 3))
    cdef list all_points = [P1] + [positions[i] for i in range(n - 1)] + [P2]
    cdef list residuals = []
    cdef int i
    cdef double m = total_mass / (n - 1)
    cdef np.ndarray[np.float64_t, ndim=1] F_left = np.zeros(3)
    cdef np.ndarray[np.float64_t, ndim=1] F_right = np.zeros(3)
    cdef double dummy_len

    for i in range(1, n):
        compute_spring_force(all_points[i - 1], all_points[i], s0_list[i - 1], c_real, F_left, &dummy_len)
        compute_spring_force(all_points[i + 1], all_points[i], s0_list[i], c_real, F_right, &dummy_len)
        F_gravity = np.array([0, 0, -m * g_scaled], dtype=np.float64)
        net_force = F_left + F_right + F_gravity
        residuals.append(net_force)

    return np.concatenate(residuals)

@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.ndarray[np.float64_t, ndim=2] refined_initial_guess(np.ndarray[np.float64_t, ndim=1] P1,
                                                             np.ndarray[np.float64_t, ndim=1] P2,
                                                             int n,
                                                             double max_sag):
    cdef list guesses = []
    cdef int i
    cdef double t, sag
    cdef np.ndarray[np.float64_t, ndim=1] point
    for i in range(1, n):
        t = i / (n - 1)
        point = (1 - t) * P1 + t * P2
        sag = -max_sag * (1 - 4 * t * (1 - t))
        point = point.copy()
        point[2] += sag
        guesses.append(point)
    return np.array(guesses, dtype=np.float64)

@cython.boundscheck(False)
@cython.wraparound(False)
def solve_spring_mass(np.ndarray[np.float64_t, ndim=1] P1,
                      np.ndarray[np.float64_t, ndim=1] P2,
                      int n,
                      np.ndarray[np.float64_t, ndim=1] s0_list,
                      double c,
                      double total_mass,
                      double g,
                      double max_sag,
                      np.ndarray[np.float64_t, ndim=2] out_positions):
    cdef np.ndarray[np.float64_t, ndim=2] guesses = refined_initial_guess(P1, P2, n, max_sag)
    cdef np.ndarray[np.float64_t, ndim=1] flat_guess = guesses.flatten()

    result = root(equilibrium_residuals, flat_guess,
                  args=(P1, P2, n, s0_list, total_mass, c, g),
                  method='hybr', options={'maxfev': 100000})

    if not result.success:
        raise RuntimeError("Solver failed with or : " + result.message)

    cdef np.ndarray[np.float64_t, ndim=2] positions = result.x.reshape((n - 1, 3))
    cdef int i, j
    for i in range(n - 1):
        for j in range(3):
            out_positions[i, j] = positions[i, j]

#