/*
 * spring_mass_solver_final.h
 * Header for spring_mass_solver_final.c
 */

 #ifndef SPRING_MASS_SOLVER_FINAL_H
 #define SPRING_MASS_SOLVER_FINAL_H
 
 #ifdef __cplusplus
 extern "C" {
 #endif
 
 #ifdef _WIN32
 #define DLL_EXPORT __declspec(dllexport)
 #else
 #define DLL_EXPORT
 #endif
 
 /* Error codes */
 #define ERR_ALLOC_FAILED -1
 #define ERR_SOLVE_FAILED -2
 #define ERR_NEWTON_DID_NOT_CONVERGE -3
 
 /**
  * @brief Solve the static equilibrium of a spring-mass system under gravity.
  *
  * @param P1 3D coordinates of the first fixed point.
  * @param P2 3D coordinates of the second fixed point.
  * @param n Number of segments (n-1 intermediate masses).
  * @param c Spring constant for each segment.
  * @param total_mass Total mass distributed evenly across all segments.
  * @param g_vec Gravity vector (array of 3 doubles).
  * @param out_positions Output array of size (n-1)*3, filled with positions of internal nodes.
  * @return 0 if successful, otherwise a negative error code.
  */
 DLL_EXPORT int solve_spring_mass_c(
     double* P1,
     double* P2,
     int n,
     double c,
     double total_mass,
     double* g_vec,
     double* out_positions
 );
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* SPRING_MASS_SOLVER_FINAL_H */
 