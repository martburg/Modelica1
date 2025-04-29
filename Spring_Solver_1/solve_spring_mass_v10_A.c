#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define MAX_ITER 100
#define TOL 1e-6
#define LS_REDUCTION 0.5
#define LS_MAX_TRIALS 8

static void compute_spring_force(const double *p1, const double *p2, double s0, double c, double *F_out) {
    double delta[3];
    double length_sq = 0.0;
    for (int i = 0; i < 3; ++i) {
        delta[i] = p2[i] - p1[i];
        length_sq += delta[i] * delta[i];
    }
    double length = sqrt(length_sq);
    if (length < 1e-12) length = 1e-12; // prevent division by zero
    double coeff = -(length - s0) * c / length;
    for (int i = 0; i < 3; ++i)
        F_out[i] = coeff * delta[i];
}

static void compute_residuals(
    const double *positions_flat, const double *P1, const double *P2, int n,
    double s0, double c, double m, const double* g_vec, double scale_pos, double *residuals_out)
{
    const double *prev = P1;

    // Normalize g_vec
    double g_vec_unit[3];
    double g_norm = sqrt(g_vec[0]*g_vec[0] + g_vec[1]*g_vec[1] + g_vec[2]*g_vec[2]);
    g_vec_unit[0] = g_vec[0] / g_norm;
    g_vec_unit[1] = g_vec[1] / g_norm;
    g_vec_unit[2] = g_vec[2] / g_norm;

    // Vector P1 -> P2
    double v[3] = {P2[0]-P1[0], P2[1]-P1[1], P2[2]-P1[2]};
    double v_norm_sq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    double v_norm = sqrt(v_norm_sq);

    // Alignment
    double dot_vg = v[0]*g_vec_unit[0] + v[1]*g_vec_unit[1] + v[2]*g_vec_unit[2];
    double cos_theta = dot_vg / v_norm;
    double align = fabs(cos_theta);

    // Dynamic penalty scaling
    double penalty_scale = 1.0 - align;
    double k_soft_base = 0.1 * c;
    double k_soft = k_soft_base * (0.01 + 0.99 * penalty_scale);

    // Lateral correction stiffness
    double k_lat_base = 0.1 * c; // or 0.2 * c if you want stronger pullback
    double k_lat = k_lat_base * (0.01 + 0.99 * penalty_scale); // optionally also scaled

    for (int i = 0; i < n-1; ++i) {
        const double *mid = &positions_flat[i*3];
        const double *next = (i == n-2) ? P2 : &positions_flat[(i+1)*3];
        double F_l[3], F_r[3];

        compute_spring_force(prev, mid, s0, c, F_l);
        compute_spring_force(next, mid, s0, c, F_r);

        for (int j = 0; j < 3; ++j)
            residuals_out[i*3+j] = F_l[j] + F_r[j] + m * g_vec[j];

        // === Penalty spring to rope line ===

        double p_minus_p0[3] = {mid[0] - P1[0], mid[1] - P1[1], mid[2] - P1[2]};
        double t_proj = (p_minus_p0[0]*v[0] + p_minus_p0[1]*v[1] + p_minus_p0[2]*v[2]) / v_norm_sq;
        double projection[3] = {P1[0] + t_proj*v[0], P1[1] + t_proj*v[1], P1[2] + t_proj*v[2]};

        for (int j = 0; j < 3; ++j)
            residuals_out[i*3+j] += -k_soft * (mid[j] - projection[j]) / scale_pos;

        // === New: Lateral deviation correction ===

        double deviation[3] = {mid[0] - projection[0], mid[1] - projection[1], mid[2] - projection[2]};
        double t_lateral = (deviation[0]*v[0] + deviation[1]*v[1] + deviation[2]*v[2]) / v_norm_sq;
        double deviation_along[3] = { t_lateral*v[0], t_lateral*v[1], t_lateral*v[2] };
        double deviation_lateral[3] = {
            deviation[0] - deviation_along[0],
            deviation[1] - deviation_along[1],
            deviation[2] - deviation_along[2]
        };

        for (int j = 0; j < 3; ++j)
            residuals_out[i*3+j] += -k_lat * deviation_lateral[j] / scale_pos;

        prev = mid;
    }
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Assume you already have:
// void compute_spring_force(const double* p1, const double* p2, double s0, double c, double* F_out);

static void dynamic_relaxation(
    double* x, const double* P1, const double* P2, int n,
    double s0, double c_initial, double m, const double* g_vec,
    double dt_given, int max_steps)
{
    int dof = (n-1) * 3;
    double* v = calloc(dof, sizeof(double));
    double* F = calloc(dof, sizeof(double));
    double* s0_adapted = calloc(n, sizeof(double));

    // Compute rope and gravity properties
    double dx_line[3] = {P2[0] - P1[0], P2[1] - P1[1], P2[2] - P1[2]};
    double L = sqrt(dx_line[0]*dx_line[0] + dx_line[1]*dx_line[1] + dx_line[2]*dx_line[2]);
    double e_rope[3] = {dx_line[0]/L, dx_line[1]/L, dx_line[2]/L};

    double g_norm = sqrt(g_vec[0]*g_vec[0] + g_vec[1]*g_vec[1] + g_vec[2]*g_vec[2]);
    double e_g[3] = {g_vec[0]/g_norm, g_vec[1]/g_norm, g_vec[2]/g_norm};

    // --- Gravity component orthogonal to rope ---
    double dot_g_rope = e_g[0]*e_rope[0] + e_g[1]*e_rope[1] + e_g[2]*e_rope[2];
    double g_perp[3] = {
        e_g[0] - dot_g_rope * e_rope[0],
        e_g[1] - dot_g_rope * e_rope[1],
        e_g[2] - dot_g_rope * e_rope[2]
    };
    double g_perp_norm = sqrt(g_perp[0]*g_perp[0] + g_perp[1]*g_perp[1] + g_perp[2]*g_perp[2]);
    if (g_perp_norm > 1e-12) {
        g_perp[0] /= g_perp_norm;
        g_perp[1] /= g_perp_norm;
        g_perp[2] /= g_perp_norm;
    } else {
        g_perp[0] = g_perp[1] = g_perp[2] = 0.0;
    }

    // --- Precompute adapted s0 ---
    for (int i = 0; i < n; ++i) {
        int masses_below = n - i;
        s0_adapted[i] = (L / n) + (m * g_norm * masses_below) / c_initial;
    }

    // --- Initialize shape + sag ---
    double T0 = c_initial * s0;
    double a = (T0 > 0 && g_norm > 0) ? (T0 / (m * g_norm)) : 1.0;

    for (int i = 0; i < n-1; ++i) {
        double t = (double)(i+1) / (double)n;
        for (int j = 0; j < 3; ++j)
            x[i*3+j] = (1.0 - t) * P1[j] + t * P2[j];

        double s = t * L;
        double sag_amount = a * (cosh(s/a) - 1.0);

        for (int j = 0; j < 3; ++j)
            x[i*3+j] += sag_amount * g_perp[j];
    }

    printf("\nInitial relaxed positions before relaxation:\n");
    for (int i = 0; i < n-1; ++i)
        printf("Node %d: [%.6f, %.6f, %.6f]\n", i+1, x[i*3+0], x[i*3+1], x[i*3+2]);

    // --- Relaxation loop ---
    double c = c_initial;
    double damping_factor = 2.5; // fixed
    double stop_threshold = 1e-4;
    double k_end_correction = 2.0; // fixed
    double v_norm_window = 0.0;
    int window_size = 5000; // check every 5000 steps
    double relax_factor = 0.95;
    double stiffen_factor = 100.0;

    // For monitoring
    double last_vnorm = 1e10;
    double moving_average_vnorm = 10;
    const double alpha = 0.5;  // moving average update factor

    double dt, damping;
    dt = fmin(dt_given, 0.05 / sqrt(c / m));
    damping = damping_factor * sqrt(c / m);

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

            double p_minus_p1[3] = {
                mid[0] - P1[0],
                mid[1] - P1[1],
                mid[2] - P1[2]
            };
            double t_proj = (p_minus_p1[0]*e_rope[0] + p_minus_p1[1]*e_rope[1] + p_minus_p1[2]*e_rope[2]) / L;

            if (t_proj < 0.0) {
                for (int j = 0; j < 3; ++j)
                    F[i*3+j] += -k_end_correction * c * (mid[j] - P1[j]);
            } else if (t_proj > 1.0) {
                for (int j = 0; j < 3; ++j)
                    F[i*3+j] += -k_end_correction * c * (mid[j] - P2[j]);
            }
        }

        // Update velocities and positions
        double v_norm_sq = 0.0;
        for (int i = 0; i < dof; ++i) {
            double a = F[i] / m;
            v[i] += a * dt;
            x[i] += v[i] * dt;
            v_norm_sq += v[i] * v[i];
        }

        
        double v_norm = sqrt(v_norm_sq);
        printf("v_norm %e \n", v_norm);

        // Moving average
        moving_average_vnorm = (1.0 - alpha) * moving_average_vnorm + alpha * v_norm;

        if (step % 1 == 0)
        printf("Moving Average %e \n", moving_average_vnorm);


        // Monitoring every window_size steps
        if (step % window_size == 0 && step > 0) {
            printf("Relaxation step %d: v_norm = %.6e (c = %.3e)\n", step, v_norm, c);
        }
        if (moving_average_vnorm > stop_threshold) {
            // Oscillation growing: stiffen                
            c *= stiffen_factor;
            printf("Oscillation detected, stiffening c = %e\n",c);
            dt = fmin(dt_given, 0.05 / sqrt(c / m));
            damping = damping_factor * sqrt(c / m);
            // Recompute s0_adapted
            for (int i = 0; i < n; ++i) {
                int masses_below = n - i;
                s0_adapted[i] = (L / n) + (m * g_norm * masses_below) / c;
            }
        } else {
            // Converging: relax a bit
            c *= relax_factor;
            c= 1000;
            printf("convergence detected, relaxing c = %e\n",c);
            dt = fmin(dt_given, 0.05 / sqrt(c / m));
            damping = damping_factor * sqrt(c / m);
            for (int i = 0; i < n; ++i) {
                int masses_below = n - i;
                s0_adapted[i] = (L / n) + (m * g_norm * masses_below) / c;
            }
        }
        v_norm_window = v_norm;
        

        if (v_norm < stop_threshold) {
            printf("Dynamic relaxation converged at step %d\n", step);
            break;
        }
    }

    free(v);
    free(F);
    free(s0_adapted);
}



static int solve_linear_system(double *A, double *b, int n) {
    for (int k = 0; k < n; ++k) {
        int max_row = k;
        for (int i = k + 1; i < n; ++i)
            if (fabs(A[i * n + k]) > fabs(A[max_row * n + k])) max_row = i;

        if (fabs(A[max_row * n + k]) < 1e-12) return -1;

        if (max_row != k) {
            for (int j = 0; j < n; ++j) {
                double tmp = A[k * n + j];
                A[k * n + j] = A[max_row * n + j];
                A[max_row * n + j] = tmp;
            }
            double tmp = b[k]; b[k] = b[max_row]; b[max_row] = tmp;
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

#ifdef _WIN32
#define DLL_EXPORT __declspec(dllexport)
#else
#define DLL_EXPORT
#endif

DLL_EXPORT int solve_spring_mass_c(
    double* P1, double* P2, int n,
    double c, double total_mass, double* g_vec,
    double* out_positions)
{
    printf("P1 = [%f, %f, %f], P2 = [%f, %f, %f]\n", P1[0], P1[1], P1[2], P2[0], P2[1], P2[2]);
    printf("n = %i , c = %f , TM = %f ",n,c,total_mass);
    printf("g_vec = [%f,%f,%f]\n",g_vec[0],g_vec[1],g_vec[2]);

    int dof = (n - 1) * 3;
    double *x = malloc(dof * sizeof(double));
    double *res = malloc(dof * sizeof(double));
    double *J = malloc(dof * dof * sizeof(double));
    double *dx = malloc(dof * sizeof(double));
    double m = total_mass / (n - 1);
    if (!x || !res || !J || !dx) return -1;

    double dx_line[3] = {P2[0] - P1[0], P2[1] - P1[1], P2[2] - P1[2]};
    double L = sqrt(dx_line[0]*dx_line[0] + dx_line[1]*dx_line[1] + dx_line[2]*dx_line[2]);
    double s0 = L / n;
    double s0_scaled = 1.0 / n;

    double scale_pos = L;
    double g_norm = sqrt(g_vec[0]*g_vec[0] + g_vec[1]*g_vec[1] + g_vec[2]*g_vec[2]);
    double scale_force = m * g_norm;

    double P1_scaled[3] = {P1[0]/scale_pos, P1[1]/scale_pos, P1[2]/scale_pos};
    double P2_scaled[3] = {P2[0]/scale_pos, P2[1]/scale_pos, P2[2]/scale_pos};
    double g_vec_scaled[3] = {g_vec[0]/scale_force, g_vec[1]/scale_force, g_vec[2]/scale_force};
    double c_scaled = c / scale_force;
    double m_scaled = m / scale_force;

    // Simple straight line init
    for (int i = 0; i < n-1; ++i) {
        double t = (double)(i+1) / (double)n;
        for (int j = 0; j < 3; ++j)
            x[i*3+j] = (1.0 - t) * P1_scaled[j] + t * P2_scaled[j];
    }

    dynamic_relaxation(x, P1_scaled, P2_scaled, n, s0_scaled, c_scaled, m_scaled, g_vec_scaled, 0.001, 1000000);

    printf("\nInitial relaxed positions after dynamic relaxation:\n");
    for (int i = 0; i < n - 1; ++i) {
        printf("Node %d: [%f, %f, %f]\n", i + 1, x[i*3+0] * scale_pos, x[i*3+1] * scale_pos, x[i*3+2] * scale_pos);
    }

    for (int iter = 0; iter < MAX_ITER; ++iter) {
        compute_residuals(x, P1_scaled, P2_scaled, n, s0_scaled, c_scaled, m_scaled, g_vec_scaled, scale_pos, res);
        double res_norm = 0.0;
        for (int i = 0; i < dof; ++i) res_norm += res[i] * res[i];
        res_norm = sqrt(res_norm);

        if (res_norm < TOL) break;

        printf("Iteration %d: residual norm = %.6e, max residual component = %.6e\n", iter, res_norm, TOL);

        if (res_norm < TOL) {
            printf("Converged at iteration %d!\n", iter);
            break;
        }

        // FD Jacobian
        for (int i = 0; i < dof; ++i) {
            double eps = fmax(1e-6, 1e-4 * fabs(x[i]));
            double orig = x[i];
            x[i] = orig + eps;
            double *res_plus = malloc(dof * sizeof(double));
            compute_residuals(x, P1_scaled, P2_scaled, n, s0_scaled, c_scaled, m_scaled, g_vec_scaled, scale_pos, res_plus);
            x[i] = orig - eps;
            double *res_minus = malloc(dof * sizeof(double));
            compute_residuals(x, P1_scaled, P2_scaled, n, s0_scaled, c_scaled, m_scaled, g_vec_scaled, scale_pos, res_minus);
            x[i] = orig;
            for (int j = 0; j < dof; ++j)
                J[j*dof+i] = (res_plus[j] - res_minus[j]) / (2*eps);
            free(res_plus);
            free(res_minus);
        }

        for (int i = 0; i < dof; ++i) dx[i] = -res[i];
        if (solve_linear_system(J, dx, dof) != 0) return -2;

        // Line search
        double alpha = 1.0;
        for (int trial = 0; trial < LS_MAX_TRIALS; ++trial) {
            double *x_trial = malloc(dof * sizeof(double));
            for (int i = 0; i < dof; ++i)
                x_trial[i] = x[i] + alpha * dx[i];
            double *res_trial = malloc(dof * sizeof(double));
            compute_residuals(x_trial, P1_scaled, P2_scaled, n, s0_scaled, c_scaled, m_scaled, g_vec_scaled, scale_pos, res_trial);
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
        x[i] *= scale_pos;

    memcpy(out_positions, x, dof * sizeof(double));

    printf("Final output positions:\n");

    for (int i = 0; i < n - 1; ++i) {
        printf("Node %d: [%f, %f, %f]\n", i + 1, out_positions[i*3+0], out_positions[i*3+1], out_positions[i*3+2]);
    }


    free(x); free(res); free(J); free(dx);
    return 0;
}
