/**
 * @file solve_rope_length_v09.c
 * @brief Implementation of the 3D spring-mass rope solver.
 *
 * This source file defines the full implementation of the `solve_rope_length` interface
 * declared in `solve_rope_length.h`, along with all internal numerical routines used
 * for initialization, dynamic relaxation, Newton refinement, and force decomposition.
 *
 * The module supports dynamic linking as a shared library (DLL or .so) and is compatible
 * with external environments such as C/C++, Python (via ctypes), or Modelica.
 */

#include "solve_rope_length_v09.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

enum LogLevel {
    LOG_NONE = 0,
    LOG_ERROR = 1,
    LOG_WARN  = 2,
    LOG_INFO  = 3,
    LOG_MOREINFO = 4,
    LOG_DEBUG = 5
};

static int CURRENT_LOG_LEVEL = 5;  // Default: silent

#define log_debug(fmt, ...) if (CURRENT_LOG_LEVEL >= LOG_DEBUG) fprintf(stderr, "[DEBUG] " fmt, ##__VA_ARGS__)
#define log_moreinfo(fmt, ...) if (CURRENT_LOG_LEVEL >= LOG_MOREINFO) fprintf(stderr, "[MOREINFO] " fmt, ##__VA_ARGS__)
#define log_info(fmt, ...)  if (CURRENT_LOG_LEVEL >= LOG_INFO)  fprintf(stderr, "[INFO] "  fmt, ##__VA_ARGS__)
#define log_warn(fmt, ...)  if (CURRENT_LOG_LEVEL >= LOG_WARN)  fprintf(stderr, "[WARN] "  fmt, ##__VA_ARGS__)
#define log_error(fmt, ...) if (CURRENT_LOG_LEVEL >= LOG_ERROR) fprintf(stderr, "[ERROR] " fmt, ##__VA_ARGS__)


#define IS_INVALID(x) (!(isfinite(x)))
#define ARC_STEPS_FINE 1000
#define NEWTON_TOL     1e-10
#define NEWTON_MAX     20

#define MAX_ITER_CAT 100
#define TOL_CAT 1e-8
#define SIGN(x) ((x) >= 0 ? 1 : -1)

#define MAX_ITER_ARC 100
#define TOL_ARC 1e-8

extern void dgesv_(int* N, int* NRHS, double* A, int* LDA,
    int* IPIV, double* B, int* LDB, int* INFO);

// Recast of spring-mass rope solver as energy minimization (C code)
// Computes total energy and its gradient (i.e., force residual) using BFGS optimizer

static inline double norm3(const double v[3]) {
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

double vec3_dot(const double v1[3], const double v2[3]) {
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

// Helper to compute projection
void project_parallel_and_perpendicular(const double* F, const double* g_unit, double* F_parallel, double* F_perp) {
    double dot = vec3_dot(F, g_unit);
    for (int j = 0; j < 3; ++j) {
        F_parallel[j] = dot * g_unit[j];
        F_perp[j] = F[j] - F_parallel[j];
    }
}

static void compute_spring_force(const double *p1, const double *p2, double s0, double c, double *F_out) {
    double delta[3];
    double length_sq = 0.0;
    for (int i = 0; i < 3; ++i) {
        delta[i] = p2[i] - p1[i];
        length_sq += delta[i] * delta[i];
    }
    double length = sqrt(length_sq);

    if (length < 1e-8) length = 1e-8; // prevent division by zero
    double coeff = -(length - s0) * c / length;
    for (int i = 0; i < 3; ++i)
        F_out[i] = coeff * delta[i];
}

static void compute_residuals(
    const double *positions_flat, const double *P1, const double *P2, int n,
    const double* s0, double c, double m, const double* g_vec,  double *residuals_out)
{
    const double *prev = P1;

    for (int i = 0; i < n-1; ++i) {
        const double *mid = &positions_flat[i*3];
        const double *next = (i == n-2) ? P2 : &positions_flat[(i+1)*3];
        double F_l[3], F_r[3];

        compute_spring_force(prev, mid, s0[i], c, F_l);
        compute_spring_force(next, mid, s0[i+1], c, F_r);

        for (int j = 0; j < 3; ++j)
            residuals_out[i*3+j] = F_l[j] + F_r[j] + m * g_vec[j];

        prev = mid;
    }
}

static void cubic_parabola_point(double* out, const double* P1, const double* P2, const double* dir, double a, double t) {
    double sag = -4.0 * a * t * (1.0 - t);
    for (int j = 0; j < 3; ++j)
        out[j] = (1 - t) * P1[j] + t * P2[j] + sag * dir[j];
}

static void add_spring_jacobian(double *J, int i, int j, int dof, const double *xi, const double *xj, double rest_len, double k, int sign) {
    double d[3];
    for (int k = 0; k < 3; ++k)
        d[k] = xj[k] - xi[k];

    double L = norm3(d);
    if (L < 1e-8) return;  // prevent division by zero

    // Compute coefficients based on the spring force derivative
    double coeff1 = -k * (1.0 - rest_len / L);                // transverse
    double coeff2 = -k * rest_len / (L * L * L);              // longitudinal

    // Compute Jacobian block
    double block[9];
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            double I_rc = (r == c) ? 1.0 : 0.0;
            block[r * 3 + c] = coeff1 * I_rc + coeff2 * d[r] * d[c];
        }
    }

    // Insert block into global Jacobian
    int row = i * 3;
    int col = j * 3;
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            J[(row + r) * dof + (col + c)] += sign * block[r * 3 + c];
        }
    }
}

static void report_endpoint_forces_springs(
    const double* P1, const double* P2, const double* x, int n,
    const double* s0_rest, double c, const double* g_vec,
    double m, double scale_force, double total_mass, double* F_P1_out,double* F_P2_out)
{

    double F_P1[3] = {0}, F_P2[3] = {0};
    compute_spring_force(&x[0], P1, s0_rest[0], c, F_P1); // Force ON P1 from x[0]
    compute_spring_force(&x[(n - 2) * 3], P2, s0_rest[n - 1], c, F_P2); // Force ON P2 from x[n-2]

    // Gravity unit vector
    double g_norm = norm3(g_vec);
    double g_unit[3] = {
        g_vec[0] / g_norm,
        g_vec[1] / g_norm,
        g_vec[2] / g_norm
    };

    // Rope direction unit vector
    double d[3] = {
        P2[0] - P1[0],
        P2[1] - P1[1],
        P2[2] - P1[2]
    };
    double d_norm = norm3(d);
    double d_unit[3] = {
        d[0] / d_norm,
        d[1] / d_norm,
        d[2] / d_norm
    };

    // Plane vector perpendicular to gravity
    double d_dot_g = vec3_dot(d_unit, g_unit);
    double e_plane[3] = {
        d_unit[0] - d_dot_g * g_unit[0],
        d_unit[1] - d_dot_g * g_unit[1],
        d_unit[2] - d_dot_g * g_unit[2]
    };
    double e_plane_norm = norm3(e_plane);
    if (e_plane_norm > 1e-12) {
        for (int i = 0; i < 3; ++i)
            e_plane[i] /= e_plane_norm;
    } else {
        e_plane[0] = e_plane[1] = e_plane[2] = 0.0;
    }
     // Decompose
    double f1_parallel = vec3_dot(F_P1, g_unit);
    double f1_plane    = vec3_dot(F_P1, e_plane);
    double f2_parallel = vec3_dot(F_P2, g_unit);
    double f2_plane    = vec3_dot(F_P2, e_plane);

    log_info("=== Endpoint Force Decomposition (Spring elongation Newton) ===\n");
    log_info("F_P1: [%f, %f, %f]\n", F_P1[0], F_P1[1], F_P1[2]);
    log_info("||F_P1_perp|| = %.6f, ||F_P1_para|| = %.6f\n",fabs(f1_plane),fabs(f1_parallel));
    //log_info("F_P1: [%f, %f, %f]\n", F_P1[0], F_P1[1], F_P1[2]);
    log_info("F_P2: [%f, %f, %f]\n", F_P2[0], F_P2[1], F_P2[2]);
    //log_info("||P2 F|| = %.6f\n",norm3(F_P2));
    log_info("||F_P2_perp|| = %.6f, ||F_P2_para|| = %.6f\n\n",fabs(f2_plane),fabs(f2_parallel));
    
   
/*    double F_net[3] = {F_P1[0] + F_P2[0], F_P1[1] + F_P2[1], F_P1[2] + F_P2[2]};
    double G[3] = {
        total_mass * g_vec[0],
        total_mass * g_vec[1],
        total_mass * g_vec[2],
    };
    log_info("Delta F (F_net - m*g) = %e\n", norm(G) -norm3(F_net));*/
    
    F_P1_out[0] = F_P1[0];F_P1_out[1] = F_P1[1];F_P1_out[2] = F_P1[2];
    F_P2_out[0] = F_P2[0];F_P2_out[1] = F_P2[1];F_P2_out[2] = F_P2[2];  
}

void report_endpoint_forces_weight(
    const double* P1, const double* P2, const double* x, int n,
    const double* g_vec, double total_mass,
    double* F_P1_out, double* F_P2_out)
{
    log_debug("P1:[%f,%f,%f], P2:[%f,%f,%f]\n",P1[0],P1[1],P1[2],P2[0],P2[1],P2[2]);
    log_debug("n:%d, g:[%f,%f.%f], m_%f\n",n,g_vec[0],g_vec[1],g_vec[2],total_mass);
    log_debug("Positions input to report weight\n");
    for (int i = 0; i < n - 1; ++i) {
        log_debug("Node %-3d: [%12.6f, %12.6f, %12.6f]\n", i + 1,
            x[i*3+0], x[i*3+1],x[i*3+2] );
    }

    // Normalize gravity vector
    double g_unit[3];
    double g_mag = norm3(g_vec);
    for (int j = 0; j < 3; ++j) g_unit[j] = g_vec[j] / g_mag;

    double G[3] = {
         total_mass * g_vec[0],
         total_mass * g_vec[1],
         total_mass * g_vec[2],
    };
    
    // Tangent vectors at endpoints
    double t1[3] = {
        x[0] - P1[0],
        x[1] - P1[1],
        x[2] - P1[2]
    };
    double t2[3] = {
        P2[0] - x[(n - 2) * 3 + 0],
        P2[1] - x[(n - 2) * 3 + 1],
        P2[2] - x[(n - 2) * 3 + 2]
    };

    double t1_norm = norm3(t1);
    double t2_norm = norm3(t2);
    for (int j = 0; j < 3; ++j) {
        t1[j] /= t1_norm;
        t2[j] /= t2_norm;
    }

    // Construct 3x2 matrix T = [t1 t2]

    // Compute T^T T (2x2)
    double A[4] = {
        vec3_dot(t1, t1), vec3_dot(t1, t2),
        vec3_dot(t2, t1), vec3_dot(t2, t2)
    };

    // Compute T^T G (2x1)
    double b[2] = {
        vec3_dot(t1, G),
        vec3_dot(t2, G)
    };

    // Solve A * lambda = b using 2x2 linear solver
    double det = A[0]*A[3] - A[1]*A[2];
    double lambda1 = ( A[3]*b[0] - A[1]*b[1]) / det;
    double lambda2 = (-A[2]*b[0] + A[0]*b[1]) / det;

    // Compute F1 = lambda1 * t1, F2 = lambda2 * t2
    double F_P1[3];
    double F_P2[3];
    for (int i = 0; i < 3; ++i) {
        F_P1[i] = lambda1 * t1[i];
        F_P2[i] = lambda2 * t2[i];
    }

    // Split forces into parallel/perpendicular components
    double F_P1_parallel[3], F_P1_perp[3];
    double F_P2_parallel[3], F_P2_perp[3];
    project_parallel_and_perpendicular(F_P1, g_unit, F_P1_parallel, F_P1_perp);
    project_parallel_and_perpendicular(F_P2, g_unit, F_P2_parallel, F_P2_perp);

    log_info("=== Endpoint Force Decomposition (Weight) ===\n");
    log_info("F_P1: [%f, %f, %f]\n", F_P1[0], F_P1[1], F_P1[2]);
    log_info("||F_P1_perp|| = %.6f, ||F_P1_para|| = %.6f\n", fabs(norm3(F_P1_perp)),fabs( norm3(F_P1_parallel)));
    log_info("F_P2: [%f, %f, %f]\n", F_P2[0], F_P2[1], F_P2[2]);
    log_info("||F_P2_perp|| = %.6f, ||F_P2_para|| = %.6f\n", fabs(norm3(F_P2_perp)),fabs( norm3(F_P2_parallel)));

    F_P1_out[0] = F_P1[0]; F_P1_out[1] = F_P1[1]; F_P1_out[2] = F_P1[2];
    F_P2_out[0] = F_P2[0]; F_P2_out[1] = F_P2[1]; F_P2_out[2] = F_P2[2];
}

// Arc length from t=0 to t=t_target
static double arc_length_to_t(const double* P1, const double* P2, const double* dir, double a, double t_target, int steps) {
    double sum = 0.0;
    double prev[3], curr[3];
    cubic_parabola_point(prev, P1, P2, dir, a, 0.0);
    for (int i = 1; i <= steps; ++i) {
        double t = t_target * ((double)i / steps);
        cubic_parabola_point(curr, P1, P2, dir, a, t);
        double dx = curr[0] - prev[0];
        double dy = curr[1] - prev[1];
        double dz = curr[2] - prev[2];
        sum += sqrt(dx*dx + dy*dy + dz*dz);
        for (int j = 0; j < 3; ++j) prev[j] = curr[j];
    }
    return sum;
}

double init_dynamic_relaxation(
    double* x, const double* P1, const double* P2,
    int n, const double* g_vec,
    double* s0_post_init, double L0_scaled, double scale_pos)
{
    //double L0_scaled = L0 / scale_pos;

    // Normalize gravity vector
    double g_norm = sqrt(g_vec[0]*g_vec[0] + g_vec[1]*g_vec[1] + g_vec[2]*g_vec[2]);
    double e_g[3] = {g_vec[0]/g_norm, g_vec[1]/g_norm, g_vec[2]/g_norm};
    double dir[3] = {-e_g[0], -e_g[1], -e_g[2]};  // sag direction

    // Fit parabola to match L0
    double a_low = 0.0, a_high = 1.0, a_mid = 0.0;
    for (int iter = 0; iter < 100; ++iter) {
        a_mid = 0.5 * (a_low + a_high);
        double arc_len = arc_length_to_t(P1, P2, dir, a_mid, 1.0, ARC_STEPS_FINE);
        if (fabs(arc_len - L0_scaled) < 1e-6) break;
        if (arc_len < L0_scaled) a_low = a_mid;
        else a_high = a_mid;
    }

    double total_arc = arc_length_to_t(P1, P2, dir, a_mid, 1.0, ARC_STEPS_FINE);

    // Internal points via Newton-Raphson arc-length solve
    for (int i = 0; i < n - 1; ++i) {
        double s_target = total_arc * (double)(i + 1) / n;

        // Initial guess: linear t
        double t = (double)(i + 1) / n;

        for (int it = 0; it < NEWTON_MAX; ++it) {
            double s_val = arc_length_to_t(P1, P2, dir, a_mid, t, ARC_STEPS_FINE);
            double s_val_eps = arc_length_to_t(P1, P2, dir, a_mid, t + 1e-6, ARC_STEPS_FINE);
            double dsdT = (s_val_eps - s_val) / 1e-6;

            double delta = (s_val - s_target) / dsdT;
            t -= delta;
            if (fabs(delta) < NEWTON_TOL) break;
            if (t < 0.0) t = 0.0;
            if (t > 1.0) t = 1.0;
        }

        // Compute point at converged t
        cubic_parabola_point(&x[i * 3], P1, P2, dir, a_mid, t);
    }

    // Compute s0_post_init
    double s0_post_init_sum = 0.0;
    for (int i = 0; i < n; ++i) {
        double xi[3], xi1[3];
        if (i == 0) {
            for (int j = 0; j < 3; ++j) {
                xi[j] = P1[j];
                xi1[j] = x[j];
            }
        } else if (i == n - 1) {
            for (int j = 0; j < 3; ++j) {
                xi[j] = x[(n - 2) * 3 + j];
                xi1[j] = P2[j];
            }
        } else {
            for (int j = 0; j < 3; ++j) {
                xi[j] = x[(i - 1) * 3 + j];
                xi1[j] = x[i * 3 + j];
            }
        }
        double dx = xi1[0] - xi[0];
        double dy = xi1[1] - xi[1];
        double dz = xi1[2] - xi[2];
        s0_post_init[i] = sqrt(dx*dx + dy*dy + dz*dz);
        s0_post_init_sum += s0_post_init[i];
    }

    return s0_post_init_sum;
}

int dynamic_relaxation(
    double* x, double *x_out, const double* P1, const double* P2, int n,
    double* s0_post_init, double c, double m, const double* g_vec,
    double dt_given, int max_steps, double scale_pos)
{
    int dof = (n - 1) * 3;

    memcpy(x_out, x, sizeof(double) * dof);

    double* v = calloc(dof, sizeof(double));
    double* F = calloc(dof, sizeof(double));
    int converged = 0;

    double omega = sqrt(c / m);
    double dt = dt_given > 0 ? fmin(dt_given, 0.05 / omega) : 0.05 / omega;
    double damping = sqrt(c * m);  // critical damping
    double damping_target = 1.2 * damping;
    double prev_v_norm = 1e10;

    for (int step = 0; step < max_steps; ++step) {
        memset(F, 0, dof * sizeof(double));

        for (int i = 0; i < n - 1; ++i) {
            const double* prev = (i == 0) ? P1 : &x_out[(i - 1) * 3];
            const double* mid  = &x_out[i * 3];
            const double* next = (i == n - 2) ? P2 : &x_out[(i + 1) * 3];

            double dx1[3] = {mid[0] - prev[0], mid[1] - prev[1], mid[2] - prev[2]};
            double dx2[3] = {mid[0] - next[0], mid[1] - next[1], mid[2] - next[2]};

            double L1 = sqrt(dx1[0]*dx1[0] + dx1[1]*dx1[1] + dx1[2]*dx1[2]);
            double L2 = sqrt(dx2[0]*dx2[0] + dx2[1]*dx2[1] + dx2[2]*dx2[2]);

            for (int j = 0; j < 3; ++j) {
                double F1 = -c * (L1 - s0_post_init[i]) / (L1 + 1e-12) * dx1[j];
                double F2 = -c * (L2 - s0_post_init[i + 1]) / (L2 + 1e-12) * dx2[j];
                F[i * 3 + j] += F1 + F2 + m * g_vec[j] - damping * v[i * 3 + j];
            }
        }

        double v_norm_sq = 0.0;
        for (int i = 0; i < dof; ++i) {
            double a = F[i] / m;
            v[i] += a * dt;
            x_out[i] += v[i] * dt;
            v_norm_sq += v[i] * v[i];
        }

        if (step % 100 == 0) {
            if (damping > damping_target)
                damping *= 0.95;
            if (v_norm_sq < prev_v_norm && dt < 0.1 / omega)
                dt *= 1.02;
            prev_v_norm = v_norm_sq;
        }

        if (v_norm_sq < 1e-7 && step > 1000) {
            damping *= 1.1;
            dt *= 0.9;
        }

        if (v_norm_sq < 1e-8 && step > 5000) {
            memset(v, 0, dof * sizeof(double));
        }
        if (step % 10000 == 0){
            log_debug("v_norm_sq %.6e\n",v_norm_sq);
        } 

        if (v_norm_sq < 1e-8) {
            log_moreinfo("Dynamic relaxation CONVERGED at step %d with v_norm_sq %.3e\n", step, v_norm_sq);
            converged = 1;
            break;
        }
    }
    free(v);
    free(F);
    if (converged == 0){
        log_error("\nDynamic relaxation didn NOT CONVERGE after %d iterations\n", max_steps);
        return SOLVE_ERROR_MAX_ITER;
    }; 
    return 0;
}

int analytic_newton_solver_3d(
    double *x, double *x_out, const double *P1, const double *P2, int n,
    const double *s0, double k, double m, const double *g_vec, double step_size, int max_iter, double tol)
{
    int dof = 3 * (n - 1);
    memcpy(x_out, x, sizeof(double) * dof);

    double *res = calloc(dof, sizeof(double));
    double *J = calloc(dof * dof, sizeof(double));
    int *ipiv = calloc(dof, sizeof(int));
    double *dx = calloc(dof, sizeof(double));
    if (!res || !J || !ipiv || !dx){
        free(res); free(J); free(ipiv); free(dx);
        return -1;}

    const double c1 = 1e-4;
    const double alpha_min = 1e-10;

    for (int iter = 0; iter < max_iter; ++iter) {
        memset(res, 0, sizeof(double) * dof);
        memset(J, 0, sizeof(double) * dof * dof);

        // Compute residuals and Jacobian
        for (int i = 0; i < n - 1; ++i) {
            const double *xi = &x_out[i * 3];
            const double *x_prev = (i == 0) ? P1 : &x_out[(i - 1) * 3];
            const double *x_next = (i == n - 2) ? P2 : &x_out[(i + 1) * 3];

            double F_prev[3], F_next[3];
            compute_spring_force(x_prev, xi, s0[i], k, F_prev);
            compute_spring_force(x_next, xi, s0[i + 1], k, F_next);

            for (int j = 0; j < 3; ++j)
                res[i * 3 + j] = F_prev[j] + F_next[j] + m * g_vec[j];

            add_spring_jacobian(J, i, i, dof, xi, x_prev, s0[i], k, +1);
            add_spring_jacobian(J, i, i, dof, xi, x_next, s0[i + 1], k, +1);
            if (i > 0)
                add_spring_jacobian(J, i, i - 1, dof, xi, x_prev, s0[i], k, -1);
            if (i < n - 2)
                add_spring_jacobian(J, i, i + 1, dof, xi, x_next, s0[i + 1], k, -1);
        }

        // Regularize diagonal
        for (int i = 0; i < dof; ++i)
            J[i * dof + i] += 1e-8;

        // Solve J dx = -res
        for (int i = 0; i < dof; ++i) dx[i] = -res[i];
        int info, nrhs = 1, lda = dof, ldb = dof;
        dgesv_(&dof, &nrhs, J, &lda, ipiv, dx, &ldb, &info);
        if (info < 0) {
            log_error("Argument %d had an illegal value\n", -info);
            free(res); free(J); free(ipiv); free(dx);
            return SOLVE_ERROR_JACOBIAN_FAILED;
        } else if (info > 0) {
            log_error("Matrix is singular: U(%d,%d) is exactly zero\n", info, info);
            free(res); free(J); free(ipiv); free(dx);
            return SOLVE_ERROR_JACOBIAN_FAILED;
        } else {
            log_info("Solution computed successfully\n");
        }
        
        // --- Armijo line search ---
        double orig_res_norm2 = 0.0;
        for (int i = 0; i < dof; ++i)
            orig_res_norm2 += res[i] * res[i];

        double alpha = 1.0;
        int wolfe_ok = 0;

        double *x_trial = malloc(dof * sizeof(double));
        double *res_trial = malloc(dof * sizeof(double));
        for (int trial = 0; trial < 12; ++trial) {
            for (int i = 0; i < dof; ++i)
                x_trial[i] = x_out[i] + alpha * dx[i];

            compute_residuals(x_trial, P1, P2, n, s0, k, m, g_vec, res_trial);

            double trial_res_norm2 = 0.0;
            for (int i = 0; i < dof; ++i)
                trial_res_norm2 += res_trial[i] * res_trial[i];

            if (trial_res_norm2 <= orig_res_norm2 + c1 * alpha * (-orig_res_norm2)) {
                memcpy(x_out, x_trial, sizeof(double) * dof);
                wolfe_ok = 1;
                break;
            }
            alpha *= 0.5;
            if (alpha < alpha_min) break;
        }
        free(x_trial);
        free(res_trial);

        if (!wolfe_ok) {
            log_debug("Wolfe line search failed; using fallback damping\n");
            for (int i = 0; i < dof; ++i)
            x_out[i] += 0.25 * dx[i];  // damped step
        }

        // Compute step size
        double norm_dx = 0.0;
        for (int i = 0; i < dof; ++i)
            norm_dx += dx[i] * dx[i];
        norm_dx = sqrt(norm_dx);

        if (norm_dx < tol) {
            log_moreinfo("3D Newton CONVERGED in %d steps with step norm %.3e\n", iter, norm_dx);
            free(res); free(J); free(ipiv); free(dx);
            return 0;
        }
    }

    log_error("3D Newton did NOT CONVERGE after %d iterations\n", max_iter);
    free(res); free(J); free(ipiv); free(dx);
    return SOLVE_ERROR_MAX_ITER;
}

double Delta_F(double total_mass, double *g_vec, double F_P1_out_n[3], double F_P2_out_n[3])
{
    double TM[3];
    TM[0] = total_mass * g_vec[0];
    TM[1] = total_mass * g_vec[1];
    TM[2] = total_mass * g_vec[2];
    double TMG = norm3(TM);

    double F_net[3] = {
        F_P1_out_n[0] + F_P2_out_n[0],
        F_P1_out_n[1] + F_P2_out_n[1],
        F_P1_out_n[2] + F_P2_out_n[2]};

    return TMG - norm3(F_net);
}

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
    int debug_level)
{
    CURRENT_LOG_LEVEL = debug_level;
    
    log_info("P1 = [%f, %f, %f], P2 = [%f, %f, %f]\n", P1[0], P1[1], P1[2], P2[0], P2[1], P2[2]);
    log_info("n = %i , Total Mass = %f , Length Factor = %f, Rope Diameter = %f, Young's MOdulus = %f \n", n, total_mass, length_factor, rope_diameter, youngs_modulus);
    log_info("g_vec = [%f,%f,%f]\n\n", g_vec[0], g_vec[1], g_vec[2]);

    int dof = (n - 1) * 3;
    double *x = malloc(dof * sizeof(double));
    double *x_init = malloc(dof * sizeof(double));
    double *x_cat = malloc((n - 2) * 3 * sizeof(double));
    double *x_relaxed = malloc(dof * sizeof(double));
    double *x_newton_3d = malloc(dof * sizeof(double));
    double *s0_init = malloc(n * sizeof(double));
    double *s0_cat = malloc(n * sizeof(double));
    double *s0_init_scaled = malloc(n * sizeof(double));
    double *s0_post_init = malloc(n * sizeof(double));
    double *s0_post_relax = malloc(n * sizeof(double));
    double *s0_post_newton_3d = malloc(n * sizeof(double));

    if (!x || !x_init || !x_cat || !x_relaxed || !x_newton_3d || !s0_init || !s0_cat || !s0_init_scaled ||
        !s0_post_init || !s0_post_relax || !s0_post_newton_3d ){
        free(x);free(x_relaxed);free(x_newton_3d);free(s0_init);free(s0_init_scaled);
        free(s0_post_init);free(s0_post_relax);free(s0_post_newton_3d);
        return -1;}

    if ((n < 2) || rope_diameter <= 0.001 || youngs_modulus <1000 || length_factor < 1 ) return INVALD_INPUT;    

    double dx_line[3] = {P2[0] - P1[0], P2[1] - P1[1], P2[2] - P1[2]};
    double L_straight = sqrt(dx_line[0]*dx_line[0] + dx_line[1]*dx_line[1] + dx_line[2]*dx_line[2]);
    double L0 = L_straight * length_factor;

/*    
    double Ar = 3.1415 * rope_diameter * rope_diameter / 4.0;
    double rope_volume = Ar * L0;
    double density = total_mass / rope_volume;
*/
    // Check if rope is nearly aligned with gravity
    double rope_dir[3] = {dx_line[0], dx_line[1], dx_line[2]};
    double rope_len = L_straight;
    for (int i = 0; i < 3; ++i) rope_dir[i] /= rope_len;

    double gravity_unit[3] = {g_vec[0], g_vec[1], g_vec[2]};
    double g_len = norm3(gravity_unit);
    for (int i = 0; i < 3; ++i) gravity_unit[i] /= g_len;

    if (g_len < 1.62) {
        log_error("Gravity to low (%.2f m/s²)\n", g_len);
        return -4;
    } else if (g_len > 49.05) {
        log_error("Gravity to high (%.2f m/s²)\n", g_len);
        return -5;
    } 

    double cos_theta = fabs(vec3_dot(rope_dir, gravity_unit));
    if (cos_theta > 0.999) {
        log_warn("Rope is nearly aligned with gravity (cosθ = %.5f). Numerical issues may occur.\n", cos_theta);
    }

    for (int i = 0; i < n; ++i) {
        s0_init[i] = L0 / n;
        s0_init_scaled[i] = length_factor / n;
    }

    // Cross-sectional area A = pi * d^2 / 4
    double A = 3.1415926 * rope_diameter * rope_diameter / 4.0;
    double L_seg = L0 / n;
    double c = (youngs_modulus * A) / L_seg;
    double m = total_mass / (n - 1);
    double g_norm = sqrt(g_vec[0]*g_vec[0] + g_vec[1]*g_vec[1] + g_vec[2]*g_vec[2]);

    // Scaling factors
    double scale_pos = L_straight;
    double scale_mass = total_mass;
    double scale_force = total_mass * g_norm;

    // Scaled problem inputs
    double P1_scaled[3] = {P1[0] / scale_pos, P1[1] / scale_pos, P1[2] / scale_pos};
    double P2_scaled[3] = {P2[0] / scale_pos, P2[1] / scale_pos, P2[2] / scale_pos};
    double m_scaled = m / scale_mass;
    double g_vec_scaled[3] = {g_vec[0] / scale_force, g_vec[1] / scale_force, g_vec[2] / scale_force};
    double c_scaled = c / scale_force;
    double L0_scaled = L0 / scale_pos;

    // --- Initialize rope with sag ---
    double s0_post_init_sum = init_dynamic_relaxation(
        x, P1_scaled, P2_scaled, n, g_vec_scaled, s0_post_init, L0_scaled, scale_pos);

    for (int i = 0; i < dof; ++i) {
        x_init[i] = x[i] * scale_pos;
        }
    log_moreinfo("Positions initiated from cubic parabola:\n");
    for (int i = 0; i < n - 1; ++i) {
        log_moreinfo("s0_post_init %-2d : %-10.6f    ", i, s0_post_init[i] * scale_pos);
        log_moreinfo("Node %-3d: [%-12.6f, %12.6f, %12.6f]\n", i + 1,
            x[i*3+0] * scale_pos, x[i*3+1] * scale_pos, x[i*3+2] * scale_pos);
    }
    log_moreinfo("s0_post_init %-3d : %f    \n",n-1,s0_post_init[n-1] * scale_pos);    
    log_info("Initial Length   = %f, Delta Length = %f]\n", s0_post_init_sum * scale_pos, s0_post_init_sum * scale_pos - L0);
    *Length_initial = s0_post_init_sum * scale_pos;

    // --- Relax dynamically ---

    *Status_dynamic = dynamic_relaxation(
        x, x_relaxed, P1_scaled, P2_scaled, n, s0_init_scaled, c_scaled, m_scaled, g_vec_scaled,0.001, 10000000, scale_pos);
    if (*Status_dynamic != 0) {
        log_error("Dynamic relaxation failed\n");
    }
    for (int i = 0; i < dof; ++i) {
        x_relaxed[i] *= scale_pos;
    }
    // Compute Dynamic Relaxation segment lengths
    double s0_post_relax_sum = 0;
    for (int i = 0; i < n; ++i) {
        double xi[3], xi1[3];
        if (i == 0) {
            for (int j = 0; j < 3; ++j) {
                xi[j]  = P1[j];
                xi1[j] = x_relaxed[j];
            }
        } else if (i == n - 1) {
            for (int j = 0; j < 3; ++j) {
                xi[j]  = x_relaxed[(n - 2) * 3 + j];
                xi1[j] = P2[j];
            }
        } else {
            for (int j = 0; j < 3; ++j) {
                xi[j]  = x_relaxed[(i - 1) * 3 + j];
                xi1[j] = x_relaxed[i * 3 + j];
            }
        }
        double dx[3] = {
            xi1[0] - xi[0], xi1[1] - xi[1], xi1[2] - xi[2]
        };
        s0_post_relax[i] = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
        s0_post_relax_sum += s0_post_relax[i];
    } 
    log_moreinfo("Positions after dynamic relaxation:\n");
    for (int i = 0; i < n - 1; ++i) {
        log_moreinfo("s0_post_relax %-3d : %-10.6f    ", i, s0_post_relax[i]);
        log_moreinfo("Node %-3d: [%12.6f, %12.6f, %12.6f]\n", i + 1,
            x_relaxed[i*3+0], x_relaxed[i*3+1], x_relaxed[i*3+2]);
    }
    log_moreinfo("s0_post_relax %d : %f    \n",n-1,s0_post_relax[n-1]);   
    log_info("Relaxed Length   = %f, Delta Length = %f]\n", s0_post_relax_sum, s0_post_relax_sum - L0);
    *Length_dynamic = s0_post_relax_sum;

    // --- Newton 3D solver ---

    *Status_newton = analytic_newton_solver_3d(x_relaxed, x_newton_3d, P1, P2, n, s0_post_relax, c, m, g_vec, 1, 20000, 1e-6);
    if (*Status_newton != 0) {
        log_error("3D Newton solver failed\n");
    }
    // Compute 3D Newton segment lengths
    double s0_newton_3d_sum = 0;
    for (int i = 0; i < n; ++i) {
        double xi[3], xi1[3];
        if (i == 0) {
            for (int j = 0; j < 3; ++j) {
                xi[j]  = P1[j];
                xi1[j] = x_newton_3d[j];
            }
        } else if (i == n - 1) {
            for (int j = 0; j < 3; ++j) {
                xi[j]  = x_newton_3d[(n - 2) * 3 + j];
                xi1[j] = P2[j];
            }
        } else {
            for (int j = 0; j < 3; ++j) {
                xi[j]  = x_newton_3d[(i - 1) * 3 + j];
                xi1[j] = x_newton_3d[i * 3 + j];
            }
        }
        double dx[3] = {
            xi1[0] - xi[0], xi1[1] - xi[1], xi1[2] - xi[2]
        };
        s0_post_newton_3d[i] = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
        s0_newton_3d_sum += s0_post_newton_3d[i];
    }    
    log_moreinfo("Positions after 3D Newton :\n");
    for (int i = 0; i < n - 1; ++i) {
        log_moreinfo("s0_post_3D_Newton %-3d : %-10.6f    ", i, s0_post_newton_3d[i] );
        log_moreinfo("Node %-3d: [%12.6f, %12.6f, %12.6f]\n", i + 1,
            x_newton_3d[i*3+0], x_newton_3d[i*3+1], x_newton_3d[i*3+2] );
    }
    log_moreinfo("s0_post_3D_Newton %-3d : %-10.6f    \n", n - 1, s0_post_newton_3d[n - 1]);
    log_info("3D Newton Length = %f, Delta Length = %f]\n\n", s0_newton_3d_sum, s0_newton_3d_sum - L0); 
    *Length_newton = s0_newton_3d_sum;

    // --- chose output ---

    memcpy(out_positions, x_newton_3d, dof * sizeof(double));

    report_endpoint_forces_springs(P1, P2, x_newton_3d, n, s0_post_relax, c, g_vec, m, scale_force,total_mass, F_P1_out_n, F_P2_out_n);
    log_info("||P1 F|| = %.6f\n",norm3(F_P1_out_n));
    log_info("||P2 F|| = %.6f\n",norm3(F_P2_out_n));

    log_warn( "Delta F - m*g (Springs) = %.6e\n\n",Delta_F(total_mass, g_vec, F_P1_out_n, F_P2_out_n));

    report_endpoint_forces_weight(P1,P2,x_newton_3d,n,g_vec,total_mass, F_P1_out_w, F_P2_out_w);
    log_info("||P1 F|| = %.6f\n",norm3(F_P1_out_w));
    log_info("||P2 F|| = %.6f\n",norm3(F_P2_out_w));

    log_warn( "Delta F - m*g (Weights) = %.6e\n\n",Delta_F(total_mass, g_vec, F_P1_out_w, F_P2_out_w)); 
    
    // --- Max sag log ---
    double max_sag = -1.0;
    int max_sag_idx = -1;
    double sag_dir[3] = {P2[0] - P1[0], P2[1] - P1[1], P2[2] - P1[2]};
    double sag_dir_len = norm3(sag_dir);
    for (int i = 0; i < 3; ++i) sag_dir[i] /= sag_dir_len;

    for (int i = 0; i < n - 1; ++i) {
        double* xi = &x_newton_3d[i * 3];
        double t;
        double proj[3];
        double dx[3] = {xi[0] - P1[0], xi[1] - P1[1], xi[2] - P1[2]};
        t = vec3_dot(dx, sag_dir);
        for (int j = 0; j < 3; ++j) proj[j] = P1[j] + t * sag_dir[j];
        double d[3] = {xi[0] - proj[0], xi[1] - proj[1], xi[2] - proj[2]};
        double dist = norm3(d);
        if (dist > max_sag) {
            max_sag = dist;
            max_sag_idx = i;
        }
    }
    log_info("Max sag at node %d: sag = %.6f, position = [%.6f, %.6f, %.6f]\n",
        max_sag_idx + 1, max_sag,
        x_newton_3d[max_sag_idx * 3 + 0],
        x_newton_3d[max_sag_idx * 3 + 1],
        x_newton_3d[max_sag_idx * 3 + 2]);

    log_moreinfo("Status Dynamc relaxation %d Status Newton %d\n",*Status_dynamic, *Status_newton);

    free(x); free(x_init); free(x_relaxed);
    free(x_newton_3d); free(s0_init);free(x_cat); free(s0_cat); free(s0_init_scaled); free(s0_post_init); free(s0_post_relax);
    free(s0_post_newton_3d);

    if (*Status_dynamic != 0 || *Status_newton != 0) {
        return SOLVE_ERROR;  // or a better code
    }
    return 0;
}


DLL_EXPORT int generate_lambda_force_table(
    const double* P1, const double* P2,
    int n, double total_mass,
    double rope_diameter, double youngs_modulus,
    const double* g_vec,
    double lambda_start, double lambda_end, int num_samples,
    double* lambda_out,           // Output: array of lambda values [num_samples]
    double* F_P1_n_out,          // Output: Newton-based F_P1 (3 * num_samples)
    double* F_P2_n_out,          // Output: Newton-based F_P2 (3 * num_samples)
    double* F_P1_w_out,          // Output: Weight-based F_P1 (3 * num_samples)
    double* F_P2_w_out,          // Output: Weight-based F_P2 (3 * num_samples)
    double* L_init_out,          // Output: Initial length [num_samples]
    double* L_cat_out,           // Output: Cat length [num_samples]
    double* L_dyn_out,           // Output: Dynamic length [num_samples]
    double* L_newton_out,        // Output: Newton length [num_samples]
    int* status_dynamic_out,     // Output: dynamic status [num_samples]
    int* status_newton_out       // Output: newton status [num_samples]
)
{
    if (!P1 || !P2 || !g_vec || !lambda_out || !F_P1_n_out || !F_P2_n_out || !F_P1_w_out || !F_P2_w_out ||
        !L_init_out || !L_cat_out || !L_dyn_out || !L_newton_out || !status_dynamic_out || !status_newton_out)
        return SOLVE_ERROR_ALLOC;

    if (num_samples < 2 || lambda_end <= lambda_start)
        return INVALD_INPUT;

    double* positions = malloc((n - 1) * 3 * sizeof(double));
    if (!positions)
        return SOLVE_ERROR_ALLOC;

    const double bias_strength = 4.0;       // steepness of transition near t = 0.5
    const double stretch_exponent = 3.0;    // how much to lift tail

    for (int i = 0; i < num_samples; ++i){ 
        double t = (double)i / (num_samples - 1);
        double lambda;

        // Parameters
        const double switch_point = 0.65;      // [0–1], where to switch to tail boost
        const double bias_strength = 5.0;
        const double front_exponent = 2.0;
        const double tail_boost = 1.4;         // >1 → denser upper tail

        if (t < switch_point) {
            double s = (tanh(bias_strength * (t - 0.5)) + 1.0) / 2.0;
            double stretched = pow(s, front_exponent);
            lambda = lambda_start + (lambda_end - lambda_start) * (switch_point * stretched);
        } else {
            // Linear taper in [switch_point, 1]
            double u = (t - switch_point) / (1.0 - switch_point);  // ∈ [0, 1]
            double linear_tail = switch_point + (1.0 - switch_point) * pow(u, tail_boost);
            lambda = lambda_start + (lambda_end - lambda_start) * linear_tail;
        }
        lambda_out[i] = lambda;
        double F_P1_n[3], F_P2_n[3], F_P1_w[3], F_P2_w[3];
        double Length_initial = NAN, Length_cat = NAN, Length_dynamic = NAN, Length_newton = NAN;
        int status_dynamic = -999, status_newton = -999;

        int status = solve_rope_length(
            P1, P2, n, total_mass, lambda,
            rope_diameter, youngs_modulus,
            g_vec, positions,
            F_P1_n, F_P2_n,
            F_P1_w, F_P2_w,
            &Length_initial, &Length_cat, &Length_dynamic, &Length_newton,
            &status_dynamic, &status_newton,
            0  // debug_level
        );

        for (int j = 0; j < 3; ++j) {
            F_P1_n_out[i * 3 + j] = F_P1_n[j];
            F_P2_n_out[i * 3 + j] = F_P2_n[j];
            F_P1_w_out[i * 3 + j] = F_P1_w[j];
            F_P2_w_out[i * 3 + j] = F_P2_w[j];
        }

        L_init_out[i]   = Length_initial;
        L_cat_out[i]    = Length_cat;
        L_dyn_out[i]    = Length_dynamic;
        L_newton_out[i] = Length_newton;

        status_dynamic_out[i] = status_dynamic;
        status_newton_out[i]  = status_newton;
}

    free(positions);
    return 0;
}