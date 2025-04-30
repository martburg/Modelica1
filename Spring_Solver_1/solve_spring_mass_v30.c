// spring_solver_cleaned.c

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifdef _WIN32
#define DLL_EXPORT __declspec(dllexport)
#else
#define DLL_EXPORT
#endif

#define MAX_ITER 100
#define TOL 1e-6
#define LS_REDUCTION 0.5
#define LS_MAX_TRIALS 8

// ==================== Utility Functions ====================

static double norm3(const double* v) {
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

static void normalize3(double* v) {
    double n = norm3(v);
    if (n > 1e-12) {
        v[0] /= n;
        v[1] /= n;
        v[2] /= n;
    }
}

static double dot3(const double* a, const double* b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// ==================== Spring Force ====================

static void compute_spring_force(const double* p_other, const double* p_self,
                                 double s0_local, double c, double* F) {
    double dx[3] = {
        p_other[0] - p_self[0],
        p_other[1] - p_self[1],
        p_other[2] - p_self[2]
    };
    double L = norm3(dx);
    if (L > 1e-12) {
        double force_magnitude = c * (L - s0_local) / L;
        F[0] = force_magnitude * dx[0];
        F[1] = force_magnitude * dx[1];
        F[2] = force_magnitude * dx[2];
    } else {
        F[0] = F[1] = F[2] = 0.0;
    }
}

// ==================== Dynamic Relaxation ====================

static void dynamic_relaxation(
    double* x, const double* P1, const double* P2, int n,
    double s0, double c, double m, const double* g_vec,
    double dt_given, int max_steps)
{
    int dof = (n-1) * 3;
    double* v = calloc(dof, sizeof(double));
    double* F = calloc(dof, sizeof(double));
    double* s0_adapted = calloc(n, sizeof(double)); // NEW: adapted spring lengths

    // --- Compute rope and gravity properties ---
    double dx_line[3] = {P2[0] - P1[0], P2[1] - P1[1], P2[2] - P1[2]};
    double L = sqrt(dx_line[0]*dx_line[0] + dx_line[1]*dx_line[1] + dx_line[2]*dx_line[2]);
    double e_rope[3] = {dx_line[0]/L, dx_line[1]/L, dx_line[2]/L};

    double g_norm = sqrt(g_vec[0]*g_vec[0] + g_vec[1]*g_vec[1] + g_vec[2]*g_vec[2]);
    double e_g[3] = {g_vec[0]/g_norm, g_vec[1]/g_norm, g_vec[2]/g_norm};

    // --- Compute plane normal (spanned by rope and gravity) ---
    double plane_normal[3] = {
        e_rope[1]*e_g[2] - e_rope[2]*e_g[1],
        e_rope[2]*e_g[0] - e_rope[0]*e_g[2],
        e_rope[0]*e_g[1] - e_rope[1]*e_g[0]
    };
    double plane_normal_norm = sqrt(
        plane_normal[0]*plane_normal[0] +
        plane_normal[1]*plane_normal[1] +
        plane_normal[2]*plane_normal[2]
    );
    if (plane_normal_norm < 1e-12) {
        if (fabs(e_rope[0]) < fabs(e_rope[1]) && fabs(e_rope[0]) < fabs(e_rope[2])) {
            plane_normal[0] = 0.0; plane_normal[1] = -e_rope[2]; plane_normal[2] = e_rope[1];
        } else if (fabs(e_rope[1]) < fabs(e_rope[2])) {
            plane_normal[0] = -e_rope[2]; plane_normal[1] = 0.0; plane_normal[2] = e_rope[0];
        } else {
            plane_normal[0] = -e_rope[1]; plane_normal[1] = e_rope[0]; plane_normal[2] = 0.0;
        }
        plane_normal_norm = sqrt(
            plane_normal[0]*plane_normal[0] +
            plane_normal[1]*plane_normal[1] +
            plane_normal[2]*plane_normal[2]
        );
    }
    plane_normal[0] /= plane_normal_norm;
    plane_normal[1] /= plane_normal_norm;
    plane_normal[2] /= plane_normal_norm;

    // --- Precompute adapted s0 ---
    for (int i = 0; i < n; ++i) {
        int masses_below = n - i;
        s0_adapted[i] = (L / n) + (m * g_norm * masses_below) / c;
    }

    // --- Initialize straight line shape ---
    for (int i = 0; i < n-1; ++i) {
        double t = (double)(i+1) / (double)n;
        for (int j = 0; j < 3; ++j)
            x[i*3+j] = (1.0 - t) * P1[j] + t * P2[j];

    } 

    // --- Relaxation loop ---
    double omega = sqrt(c / m);
    double dt = fmin(dt_given, 0.2 / omega);
    double damping = 10.5;
    double stop_threshold = 1e-4;

    for (int step = 0; step < max_steps; ++step) {
        memset(F, 0, dof * sizeof(double));

        for (int i = 0; i < n-1; ++i) {
            const double* mid = &x[i*3];
            const double* prev = (i == 0) ? P1 : &x[(i-1)*3];
            const double* next = (i == n-2) ? P2 : &x[(i+1)*3];
            double F_l[3], F_r[3];

            compute_spring_force(prev, mid, s0_adapted[i], c, F_l);
            compute_spring_force(next, mid, s0_adapted[i+1], c, F_r);

            for (int j = 0; j < 3; ++j)
                F[i*3+j] = F_l[j] + F_r[j] + m * g_vec[j] - damping * v[i*3+j];

            // --- Plane correction ---
            double p_minus_p1[3] = {
                mid[0] - P1[0],
                mid[1] - P1[1],
                mid[2] - P1[2]
            };
            double distance_to_plane =
                p_minus_p1[0]*plane_normal[0] +
                p_minus_p1[1]*plane_normal[1] +
                p_minus_p1[2]*plane_normal[2];

            double k_plane_correction = 2.0 * c;
            for (int j = 0; j < 3; ++j)
                F[i*3+j] += -k_plane_correction * distance_to_plane * plane_normal[j];

            // --- Clamp at ends ---
            double t_proj = (p_minus_p1[0]*e_rope[0] + p_minus_p1[1]*e_rope[1] + p_minus_p1[2]*e_rope[2]) / L;
            if (t_proj < 0.0) {
                for (int j = 0; j < 3; ++j)
                    F[i*3+j] += -2.0*c * (mid[j] - P1[j]);
            } else if (t_proj > 1.0) {
                for (int j = 0; j < 3; ++j)
                    F[i*3+j] += -2.0*c * (mid[j] - P2[j]);
            }
        }

        // Update positions
        double v_norm_sq = 0.0;
        for (int i = 0; i < dof; ++i) {
            double a = F[i] / m;
            v[i] += a * dt;
            x[i] += v[i] * dt;
            v_norm_sq += v[i] * v[i];
        }

        if (sqrt(v_norm_sq) < stop_threshold) {
            printf("Dynamic relaxation converged at step %d\n", step);
            break;
        }
    }

    free(v);
    free(F);
    free(s0_adapted);
}

// ==================== Residual Computation ====================

static void compute_residuals(const double* x, const double* P1, const double* P2, int n,
                              const double* s0_adapted, double c, double m, const double* g_vec,
                              double* R_out) {
    int dof = (n - 1) * 3;
    memset(R_out, 0, dof * sizeof(double));

    for (int i = 0; i < n - 1; ++i) {
        const double* mid = &x[i * 3];
        const double* prev = (i == 0) ? P1 : &x[(i - 1) * 3];
        const double* next = (i == n - 2) ? P2 : &x[(i + 1) * 3];
        double F_l[3], F_r[3];

        compute_spring_force(prev, mid, s0_adapted[i], c, F_l);
        compute_spring_force(next, mid, s0_adapted[i + 1], c, F_r);

        for (int j = 0; j < 3; ++j)
            R_out[i * 3 + j] = F_l[j] + F_r[j] + m * g_vec[j];
    }
}

// ==================== Linear Solver ====================

static int solve_linear_system(double* A, double* b, int n) {
    for (int k = 0; k < n; ++k) {
        int max_row = k;
        for (int i = k + 1; i < n; ++i)
            if (fabs(A[i * n + k]) > fabs(A[max_row * n + k]))
                max_row = i;

        if (fabs(A[max_row * n + k]) < 1e-12) return -1;

        if (max_row != k) {
            for (int j = 0; j < n; ++j) {
                double tmp = A[k * n + j];
                A[k * n + j] = A[max_row * n + j];
                A[max_row * n + j] = tmp;
            }
            double tmp = b[k];
            b[k] = b[max_row];
            b[max_row] = tmp;
        }

        for (int i = k + 1; i < n; ++i) {
            double f = A[i * n + k] / A[k * n + k];
            for (int j = k; j < n; ++j)
                A[i * n + j] -= f * A[k * n + j];
            b[i] -= f * b[k];
        }
    }

    for (int i = n - 1; i >= 0; --i) {
        for (int j = i + 1; j < n; ++j)
            b[i] -= A[i * n + j] * b[j];
        b[i] /= A[i * n + i];
    }
    return 0;
}

// ==================== Main Entry Point ====================

DLL_EXPORT int solve_spring_mass_c(
    double* P1, double* P2, int n,
    double c, double total_mass, double* g_vec,
    double* out_positions)
{
    int dof = (n - 1) * 3;
    double *x = malloc(dof * sizeof(double));
    double *res = malloc(dof * sizeof(double));
    double *J = malloc(dof * dof * sizeof(double));
    double *dx = malloc(dof * sizeof(double));
    double *s0_adapted = calloc(n, sizeof(double));
    if (!x || !res || !J || !dx || !s0_adapted) return -1;

    double dx_line[3] = {P2[0] - P1[0], P2[1] - P1[1], P2[2] - P1[2]};
    double L = norm3(dx_line);
    double s0 = L / n;

    double scale_pos = L;
    double g_norm = norm3(g_vec);
    double scale_force = total_mass * g_norm;

    double P1_scaled[3] = {P1[0]/scale_pos, P1[1]/scale_pos, P1[2]/scale_pos};
    double P2_scaled[3] = {P2[0]/scale_pos, P2[1]/scale_pos, P2[2]/scale_pos};
    double g_vec_scaled[3] = {g_vec[0]/scale_force, g_vec[1]/scale_force, g_vec[2]/scale_force};
    double c_scaled = c / scale_force;
    double m_scaled = total_mass / ((n-1) * scale_force);

    for (int i = 0; i < n - 1; ++i) {
        double t = (double)(i + 1) / (double)n;
        for (int j = 0; j < 3; ++j)
            x[i * 3 + j] = (1.0 - t) * P1_scaled[j] + t * P2_scaled[j];
    }

    dynamic_relaxation(x, P1_scaled, P2_scaled, n, s0 / scale_pos, c_scaled, m_scaled, g_vec_scaled, 0.001, 100000);

    for (int i = 0; i < n; ++i)
        s0_adapted[i] = s0 / scale_pos;

    for (int iter = 0; iter < MAX_ITER; ++iter) {
        compute_residuals(x, P1_scaled, P2_scaled, n, s0_adapted, c_scaled, m_scaled, g_vec_scaled, res);
        double res_norm = 0.0;
        for (int i = 0; i < dof; ++i)
            res_norm += res[i] * res[i];
        res_norm = sqrt(res_norm);
        if (res_norm < TOL) break;

        for (int i = 0; i < dof; ++i) {
            double eps = fmax(1e-6, 1e-4 * fabs(x[i]));
            double orig = x[i];
            x[i] = orig + eps;
            double* res_plus = malloc(dof * sizeof(double));
            compute_residuals(x, P1_scaled, P2_scaled, n, s0_adapted, c_scaled, m_scaled, g_vec_scaled, res_plus);
            x[i] = orig - eps;
            double* res_minus = malloc(dof * sizeof(double));
            compute_residuals(x, P1_scaled, P2_scaled, n, s0_adapted, c_scaled, m_scaled, g_vec_scaled, res_minus);
            x[i] = orig;
            for (int j = 0; j < dof; ++j)
                J[j * dof + i] = (res_plus[j] - res_minus[j]) / (2 * eps);
            free(res_plus);
            free(res_minus);
        }

        for (int i = 0; i < dof; ++i)
            dx[i] = -res[i];
        if (solve_linear_system(J, dx, dof) != 0) return -2;

        double alpha = 1.0;
        for (int trial = 0; trial < LS_MAX_TRIALS; ++trial) {
            double* x_trial = malloc(dof * sizeof(double));
            for (int i = 0; i < dof; ++i)
                x_trial[i] = x[i] + alpha * dx[i];
            double* res_trial = malloc(dof * sizeof(double));
            compute_residuals(x_trial, P1_scaled, P2_scaled, n, s0_adapted, c_scaled, m_scaled, g_vec_scaled, res_trial);
            double res_trial_norm = 0.0;
            for (int i = 0; i < dof; ++i)
                res_trial_norm += res_trial[i] * res_trial[i];
            res_trial_norm = sqrt(res_trial_norm);
            free(res_trial);

            if (res_trial_norm < res_norm) {
                memcpy(x, x_trial, dof * sizeof(double));
                free(x_trial);
                break;
            } else {
                alpha *= LS_REDUCTION;
                free(x_trial);
            }
        }
    }

    for (int i = 0; i < dof; ++i)
        out_positions[i] = x[i] * scale_pos;

    free(x);
    free(res);
    free(J);
    free(dx);
    free(s0_adapted);
    return 0;
}
