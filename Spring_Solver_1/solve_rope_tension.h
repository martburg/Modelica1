
typedef enum {
    SOLVE_SUCCESS = 0,
    SOLVE_ERROR_NAN = -1,
    SOLVE_ERROR_LINE_SEARCH_FAILED = -2,
    SOLVE_ERROR_JACOBIAN_FAILED = -3,
    SOLVE_ERROR_MAX_ITER = -4,
    SOLVE_ERROR_ALLOC = -5,
    SOLVE_ERROR = -6,
    INVALD_INPUT = -7
} SolveStatus;

#ifdef _WIN32
#define DLL_EXPORT __declspec(dllexport)
#else
#define DLL_EXPORT
#endif


#ifndef ROPE_SOLVER_TENSION_H
#define ROPE_SOLVER_TENSION_H

DLL_EXPORT int solve_rope_tension(
    const double *P1, const double *P2,
    int n, double total_mass, double rope_diameter,
    double youngs_modulus,
    const double F_target,
    double *L_out,
    const double *g_vec,
    double *x_out,
    double *F_P1_reported,
    double *F_P2_reported,
    int* Status_dynamic,
    int* Status_newton
) ;

#endif