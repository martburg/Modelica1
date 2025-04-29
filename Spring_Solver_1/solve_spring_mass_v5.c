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
    double k_soft = 0.1 * c;  // 10% penalty spring

    double g_vec_unit[3];
    double g_norm = sqrt(g_vec[0]*g_vec[0] + g_vec[1]*g_vec[1] + g_vec[2]*g_vec[2]);
    g_vec_unit[0] = g_vec[0] / g_norm;
    g_vec_unit[1] = g_vec[1] / g_norm;
    g_vec_unit[2] = g_vec[2] / g_norm;

    // Vector P1->P2
    double v[3] = {P2[0]-P1[0], P2[1]-P1[1], P2[2]-P1[2]};
    double v_norm_sq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    double v_norm = sqrt(v_norm_sq);

    // Alignment between rope and gravity
    double dot_vg = v[0]*g_vec_unit[0] + v[1]*g_vec_unit[1] + v[2]*g_vec_unit[2];
    double cos_theta = dot_vg / sqrt(v_norm_sq);
    double align = fabs(cos_theta); // 1 = fully aligned

    // Dynamic penalty scaling
    double penalty_scale = 1.0 - align; // 0 when parallel, 1 when orthogonal
    double k_soft_base = 0.1 * c; // your original penalty base
    double k_soft_1 = k_soft_base * (0.01 + 0.99 * penalty_scale);

    int nearly_parallel = fabs(fabs(cos_theta) - 1.0) < 1e-3;

    for (int i = 0; i < n-1; ++i) {
        const double *mid = &positions_flat[i*3];
        const double *next = (i == n-2) ? P2 : &positions_flat[(i+1)*3];
        double F_l[3], F_r[3];

        compute_spring_force(prev, mid, s0, c, F_l);
        compute_spring_force(next, mid, s0, c, F_r);

        for (int j = 0; j < 3; ++j)
            residuals_out[i*3+j] = F_l[j] + F_r[j] + m * g_vec[j];

        // Add penalty force if not nearly parallel
        double p_minus_p0[3] = {mid[0] - P1[0], mid[1] - P1[1], mid[2] - P1[2]};
        double t_proj = (p_minus_p0[0]*v[0] + p_minus_p0[1]*v[1] + p_minus_p0[2]*v[2]) / v_norm_sq;
        double projection[3] = {P1[0] + t_proj*v[0], P1[1] + t_proj*v[1], P1[2] + t_proj*v[2]};

        for (int j = 0; j < 3; ++j)
            residuals_out[i*3+j] += -k_soft_1 * (mid[j] - projection[j]) / scale_pos;
        

        prev = mid;
    }
}


static void dynamic_relaxation(
    double* x, const double* P1, const double* P2, int n,
    double s0, double c, double m, const double* g_vec,
    double dt_given, int max_steps)
{
    int dof = (n-1) * 3;
    double* v = calloc(dof, sizeof(double));
    double* F = calloc(dof, sizeof(double));
    double omega = sqrt(c / m);
    double dt = fmin(dt_given, 0.2 / omega);
    double damping = 1.5; // stronger damping
    double stop_threshold = 1e-4;

    for (int step = 0; step < max_steps; ++step) {
        for (int i = 0; i < dof; ++i) F[i] = 0.0;

        for (int i = 0; i < n-1; ++i) {
            const double *mid = &x[i*3];
            const double *prev = (i == 0) ? P1 : &x[(i-1)*3];
            const double *next = (i == n-2) ? P2 : &x[(i+1)*3];
            double F_l[3], F_r[3];
            compute_spring_force(prev, mid, s0, c, F_l);
            compute_spring_force(next, mid, s0, c, F_r);
            for (int j = 0; j < 3; ++j)
                F[i*3+j] += F_l[j] + F_r[j] + m * g_vec[j] - damping * v[i*3+j];
        }

        double v_norm_sq = 0.0;
        for (int i = 0; i < dof; ++i) {
            double a = F[i] / m;
            v[i] += a * dt;
            x[i] += v[i] * dt;
            v_norm_sq += v[i] * v[i];
        }
        if (sqrt(v_norm_sq) < stop_threshold) {
            break;
        }
    }

    free(v);
    free(F);
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

    dynamic_relaxation(x, P1_scaled, P2_scaled, n, s0_scaled, c_scaled, m_scaled, g_vec_scaled, 0.01, 1000);

    for (int iter = 0; iter < MAX_ITER; ++iter) {
        compute_residuals(x, P1_scaled, P2_scaled, n, s0_scaled, c_scaled, m_scaled, g_vec_scaled, scale_pos, res);
        double res_norm = 0.0;
        for (int i = 0; i < dof; ++i) res_norm += res[i] * res[i];
        res_norm = sqrt(res_norm);

        if (res_norm < TOL) break;

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
