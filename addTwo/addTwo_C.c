# spring_solver.pyx
#ifdef _WIN32
__declspec(dllexport)
#endif


from libc.stdio cimport fopen, fputs, fclose, FILE

cdef extern from "windows.h":
    void OutputDebugStringA(const char* lpOutputString)

cdef void dbgprint(const char* msg):
    OutputDebugStringA(msg)

from libc.stdlib cimport malloc, free
import numpy as np
cimport numpy as np
from scipy.optimize import root
cimport cython
from libc.stdio cimport printf
from numpy cimport import_array

# Tell Cython about numpy array initialization
np.import_array()

cdef void dbgfile(const char* msg):
    cdef FILE* f = fopen("debug_log.txt", "a")
    if f != NULL:
        fputs(msg, f)
        fclose(f)

# Expose C entry point for Modelica
cdef public int solve_spring_mass_c(double* P1, double* P2, int n,
                                    double* s0_list, double c,
                                    double total_mass, double g,
                                    double max_sag,
                                    double* out_positions) with gil:

    dbgprint(b"Entered solve_spring_mass_c\n")
    dbgfile(b"Starting solver\n")

    printf(">>> [DEBUG] Starting solve_spring_mass_c\n")
    printf("    n = %d, c = %f, total_mass = %f, g = %f, max_sag = %f\n", n, c, total_mass, g, max_sag)
#
#    np.import_array()
#
#    cdef int i, j
#    cdef np.ndarray[np.float64_t, ndim=1] np_P1 = np.PyArray_SimpleNewFromData(1, [3], np.NPY_FLOAT64, <void*>P1)
#    cdef np.ndarray[np.float64_t, ndim=1] np_P2 = np.PyArray_SimpleNewFromData(1, [3], np.NPY_FLOAT64, <void*>P2)
#    cdef np.ndarray[np.float64_t, ndim=1] np_s0 = np.PyArray_SimpleNewFromData(1, [n], np.NPY_FLOAT64, <void*>s0_list)
#    cdef np.ndarray[np.float64_t, ndim=2] np_out = np.empty((n - 1, 3), dtype=np.float64)
#
#    solve_spring_mass(np_P1, np_P2, n, np_s0, c, total_mass, g, max_sag, np_out)
#
#    for i in range(n - 1):
#        for j in range(3):
#            out_positions[i * 3 + j] = np_out[i, j]

    return 0
double addTwo_c(double a, double b) {
    return a + b;
}
