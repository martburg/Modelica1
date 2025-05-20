/**
 * @file solve_rope_length.h
 * @brief Header for the 3D spring-mass rope solver with Newton and dynamic relaxation.
 *
 * This solver simulates a rope suspended between two endpoints in 3D space, modeled as
 * a chain of masses and linear springs. The implementation supports:
 * - Arc-length parameterized initialization using a cubic parabola.
 * - Dynamic relaxation for rough equilibrium.
 * - Newton-based refinement using analytic Jacobian.
 * - Force decomposition (spring-based and weight-based).
 * - Tension scanning over a range of rope lengths (lambda sweep).
 *
 * Compatible with C/C++, Python (via ctypes), and Modelica.
 */

 #ifndef SOLVE_ROPE_LENGTH_H
 #define SOLVE_ROPE_LENGTH_H
 
 #ifdef _WIN32
     #define DLL_EXPORT __declspec(dllexport)
 #else
     #define DLL_EXPORT
 #endif
 
 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /**
  * @brief Status and error codes returned by solver functions.
  */
 enum SolveStatus {
     SOLVE_SUCCESS                = 0,  ///< Success
     SOLVE_ERROR_ALLOC           = -1, ///< Memory allocation failed
     SOLVE_ERROR_JACOBIAN_FAILED = -2, ///< Newton Jacobian solve failed
     SOLVE_ERROR_MAX_ITER        = -3, ///< Solver did not converge
     INVALD_INPUT                = -4, ///< Input validation failed
     SOLVE_ERROR_TOO_LOW_GRAVITY = -5, ///< Gravity magnitude too low
     SOLVE_ERROR_TOO_HIGH_GRAVITY= -6  ///< Gravity magnitude too high
 };
 
 /**
  * @brief Solve the static equilibrium of a 3D rope under gravity.
  *
  * @param[in]  P1                3D coordinates of the first endpoint.
  * @param[in]  P2                3D coordinates of the second endpoint.
  * @param[in]  n                 Number of masses (n-1 segments).
  * @param[in]  total_mass        Total mass of the rope.
  * @param[in]  length_factor     Unstretched rope length = factor * straight-line distance.
  * @param[in]  rope_diameter     Rope diameter (used for stiffness).
  * @param[in]  youngs_modulus    Young’s modulus of the rope material.
  * @param[in]  g_vec             Gravity vector [3].
  * @param[out] out_positions     Output node positions, flattened [(n-1) * 3].
  * @param[out] F_P1_out_w        Endpoint force at P1 (weight-based method).
  * @param[out] F_P2_out_w        Endpoint force at P2 (weight-based method).
  * @param[out] F_P1_out_n        Endpoint force at P1 (spring-based Newton).
  * @param[out] F_P2_out_n        Endpoint force at P2 (spring-based Newton).
  * @param[out] Length_initial    Arc length after cubic parabola initialization.
  * @param[out] Length_cat        Reserved (currently set to 0.0).
  * @param[out] Length_dynamic    Arc length after dynamic relaxation.
  * @param[out] Length_newton     Arc length after Newton refinement.
  * @param[out] Status_dynamic    Status of dynamic relaxation (0 = success).
  * @param[out] Status_newton     Status of Newton solver (0 = success).
  * @param[in]  debug_level       Verbosity level (0 = silent, 5 = debug).
  *
  * @return SolveStatus code: 0 = success, <0 = error.
  */
 DLL_EXPORT int solve_rope_length(
     double* P1, double* P2,
     int n, double total_mass, double length_factor,
     double rope_diameter, double youngs_modulus,
     double* g_vec, double* out_positions,
     double* F_P1_out_w,
     double* F_P2_out_w,
     double* F_P1_out_n,
     double* F_P2_out_n,
     double* Length_initial,
     double* Length_cat,
     double* Length_dynamic,
     double* Length_newton,
     int* Status_dynamic,
     int* Status_newton,
     int debug_level
 );
 
 /**
  * @brief Generate a force vs. rope length curve by sweeping the length factor λ.
  *
  * This function runs `solve_rope_length` multiple times with different λ values
  * between `lambda_start` and `lambda_end`, storing results for each run.
  *
  * @param[in]  P1, P2            Endpoints (3D coordinates).
  * @param[in]  n                 Number of masses.
  * @param[in]  total_mass        Total rope mass.
  * @param[in]  rope_diameter     Rope diameter.
  * @param[in]  youngs_modulus    Young’s modulus.
  * @param[in]  g_vec             Gravity vector [3].
  * @param[in]  lambda_start      Minimum length factor.
  * @param[in]  lambda_end        Maximum length factor.
  * @param[in]  num_samples       Number of length factors to sample.
  * @param[out] lambda_out        λ values used [num_samples].
  * @param[out] F_P1_n_out        Newton-based F_P1 forces [3 * num_samples].
  * @param[out] F_P2_n_out        Newton-based F_P2 forces [3 * num_samples].
  * @param[out] F_P1_w_out        Weight-based F_P1 forces [3 * num_samples].
  * @param[out] F_P2_w_out        Weight-based F_P2 forces [3 * num_samples].
  * @param[out] L_init_out        Initial lengths [num_samples].
  * @param[out] L_cat_out         Reserved (set to 0) [num_samples].
  * @param[out] L_dyn_out         Dynamic lengths [num_samples].
  * @param[out] L_newton_out      Newton lengths [num_samples].
  * @param[out] status_dynamic_out Dynamic relaxation status codes [num_samples].
  * @param[out] status_newton_out  Newton solver status codes [num_samples].
  * @param[in]  debug_level       Verbosity level (0 = silent).
  *
  * @return SolveStatus code: 0 = success, <0 = error.
  */
 DLL_EXPORT int generate_lambda_force_table(
     const double* P1, const double* P2,
     int n, double total_mass,
     double rope_diameter, double youngs_modulus,
     const double* g_vec,
     double lambda_start, double lambda_end, int num_samples,
     double* lambda_out,
     double* F_P1_n_out,
     double* F_P2_n_out,
     double* F_P1_w_out,
     double* F_P2_w_out,
     double* L_init_out,
     double* L_cat_out,
     double* L_dyn_out,
     double* L_newton_out,
     int* status_dynamic_out,
     int* status_newton_out,
     int debug_level
 );
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif // SOLVE_ROPE_LENGTH_H
// End of solve_rope_length.h
// This file is part of the 3D spring-mass rope solver. 