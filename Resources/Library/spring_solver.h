// spring_solver.h
#ifndef SPRING_SOLVER_H
#define SPRING_SOLVER_H

#ifdef _WIN32
#define DLL_EXPORT __declspec(dllexport)
#else
#define DLL_EXPORT
#endif

#ifdef __cplusplus
extern "C" {
#endif

DLL_EXPORT int solve_spring_mass_c(
    double* P1, double* P2, int n,
    double c,
    double total_mass, double* g_vec,
    double* out_positions
);

#ifdef __cplusplus
}
#endif

#endif // SPRING_SOLVER_H
