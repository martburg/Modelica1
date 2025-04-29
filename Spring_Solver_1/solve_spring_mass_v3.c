// solve_spring_mass_newton.c - Newton-based solver with finite-difference Jacobian and rescaling
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define MAX_ITER 100
#define TOL 1e-6
#define FD_EPS 1e-6

static void compute_spring_force(const double *p1, const double *p2, double s0, double c, double *F_out) {
    double delta[3], length = 0.0;
    for (int i = 0; i < 3; ++i) {
        delta[i] = p2[i] - p1[i];
        length += delta[i] * delta[i];
    }
    length = sqrt(length);
    double coeff = - (length - s0) * c / length;
    for (int i = 0; i < 3; ++i)
        F_out[i] = coeff * delta[i];
}

static void dynamic_relaxation(
    double* x, const double* P1, const double* P2, int n,
    double s0, double c, double m, const double* g_vec,
    double dt_given, int max_steps)
{
    int dof = (n-1) * 3;
    double* v = calloc(dof, sizeof(double));
    double* F = calloc(dof, sizeof(double));

    const double damping = 0.5;
    const double stop_threshold = 1e-4;

    // Adaptive dt based on stiffness and mass
    double omega = sqrt(c / m);
    double dt = fmin(dt_given, 0.2 / omega);

    //printf("Dynamic relaxation using dt = %.6e\n", dt);

    for (int step = 0; step < max_steps; ++step) {
        for (int i = 0; i < dof; ++i) F[i] = 0.0;

        for (int i = 0; i < n-1; ++i) {
            const double *mid = &x[i*3];
            const double *prev = (i == 0) ? P1 : &x[(i-1)*3];
            const double *next = (i == n-2) ? P2 : &x[(i+1)*3];
            double F_l[3], F_r[3];

            compute_spring_force(prev, mid, s0, c, F_l);
            compute_spring_force(next, mid, s0, c, F_r);

            for (int j = 0; j < 3; ++j) {
                F[i*3+j] += F_l[j] + F_r[j] + m * g_vec[j] - damping * v[i*3+j];
            }
        }

        double v_norm_sq = 0.0;

        for (int i = 0; i < dof; ++i) {
            double a = F[i] / m;
            v[i] += a * dt;
            x[i] += v[i] * dt;
            v_norm_sq += v[i] * v[i];
        }

        //if (sqrt(v_norm_sq) < stop_threshold) {
        //    printf("Relaxation converged after %d steps (v_norm = %.6e)\n", step + 1, sqrt(v_norm_sq));
        //    break;
        //}
    }

    free(v);
    free(F);
}


static void compute_residuals(
    const double *positions_flat, const double *P1, const double *P2, int n,
    double s0, double c, double m, double* g_vec, double *residuals_out)
{
    const double *prev = P1, *next;
    double F_l[3], F_r[3];
    for (int i = 0; i < n - 1; ++i) {
        const double *mid = &positions_flat[i * 3];
        next = (i == n - 2) ? P2 : &positions_flat[(i + 1) * 3];
        compute_spring_force(prev, mid, s0, c, F_l);
        compute_spring_force(next, mid, s0, c, F_r);
        for (int j = 0; j < 3; ++j)
            residuals_out[i * 3 + j] = F_l[j] + F_r[j] + m * g_vec[j];
        prev = mid;
    }
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

    double preload_force = 0.5 * m * g_norm;
    double l0 = L / n;
    double l = l0 + preload_force / c;
    double d = l0;
    double sag = sqrt((3.0 / 8.0) * d * d * (l / d - 1.0));

    double P1_scaled[3] = {P1[0]/scale_pos, P1[1]/scale_pos, P1[2]/scale_pos};
    double P2_scaled[3] = {P2[0]/scale_pos, P2[1]/scale_pos, P2[2]/scale_pos};
    double g_vec_scaled[3] = {g_vec[0]/scale_force, g_vec[1]/scale_force, g_vec[2]/scale_force};
    double c_scaled = c / scale_force;
    double m_scaled = m / scale_force;

    // Create local frame: e_x along P1->P2, e_y perpendicular
    double e_x[3] = {dx_line[0]/L, dx_line[1]/L, dx_line[2]/L};
    double tmp[3] = {0, 0, 1};
    if (fabs(e_x[2]) > 0.9) { tmp[0] = 1; tmp[1] = 0; tmp[2] = 0; }
    double e_y[3] = {
        e_x[1]*tmp[2] - e_x[2]*tmp[1],
        e_x[2]*tmp[0] - e_x[0]*tmp[2],
        e_x[0]*tmp[1] - e_x[1]*tmp[0]
    };
    double norm_e_y = sqrt(e_y[0]*e_y[0] + e_y[1]*e_y[1] + e_y[2]*e_y[2]);
    for (int j = 0; j < 3; ++j)
        e_y[j] /= norm_e_y;

    // Calculate angle between cable and gravity
    double cos_theta = fabs(
        (dx_line[0]*g_vec[0] + dx_line[1]*g_vec[1] + dx_line[2]*g_vec[2]) / 
        (L * g_norm)
    );

    if (cos_theta > 0.9) {
        printf("[Init] Cable is nearly vertical (cos(theta) = %.6f). Nonuniform gravity stretch initialization.\n", cos_theta);

        double unstretched_length = L / n;
        double z_curr = 0.0;
        double total_stretched_length = 0.0;

        // First pass: calculate total stretched length
        for (int i = 0; i < n; ++i) {
            double weight_below = (n - i) * m * g_norm;
            double weight_next = (n - i - 1) * m * g_norm;
            double avg_weight = (weight_below + weight_next) * 0.5;
            double stretch = avg_weight / c;
            total_stretched_length += unstretched_length + stretch;
        }

        // Second pass: build node positions
        z_curr = 0.0;
        for (int i = 0; i < n-1; ++i) {
            double weight_below = (n - i) * m * g_norm;
            double weight_next = (n - i - 1) * m * g_norm;
            double avg_weight = (weight_below + weight_next) * 0.5;
            double stretch = avg_weight / c;
            double segment_length = unstretched_length + stretch;
            z_curr += segment_length;

            double t = z_curr / total_stretched_length; // normalized position
            for (int j = 0; j < 3; ++j)
                x[i*3 + j] = (1.0 - t) * P1_scaled[j] + t * P2_scaled[j];
        }
    } else {
        printf("[Init] Cable is tilted/horizontal (cos(theta) = %.6f). Sag applied.\n", cos_theta);

        for (int i = 0; i < n-1; ++i) {
            double t = (double)(i+1) / (double)n;
            double r_linear[3];
            for (int j = 0; j < 3; ++j)
                r_linear[j] = (1.0 - t) * P1_scaled[j] + t * P2_scaled[j];
            double sag_profile = -4.0 * (sag/scale_pos) * t * (1.0 - t);
            for (int j = 0; j < 3; ++j)
                x[i*3 + j] = r_linear[j] + sag_profile * e_y[j];
        }
    }

    dynamic_relaxation(x, P1_scaled, P2_scaled, n, s0_scaled, c_scaled, m_scaled, g_vec_scaled, 0.01, 1000);

    for (int iter = 0; iter < MAX_ITER; ++iter) {
        compute_residuals(x, P1_scaled, P2_scaled, n, s0_scaled, c_scaled, m_scaled, g_vec_scaled, res);

        double res_norm = 0.0;
        for (int i = 0; i < dof; ++i)
            res_norm += res[i] * res[i];
        res_norm = sqrt(res_norm);

        for (int i = 0; i < dof; ++i) {
            double orig = x[i];
            x[i] = orig + FD_EPS;
            double *res_plus = malloc(dof * sizeof(double));
            compute_residuals(x, P1_scaled, P2_scaled, n, s0_scaled, c_scaled, m_scaled, g_vec_scaled, res_plus);
            x[i] = orig - FD_EPS;
            double *res_minus = malloc(dof * sizeof(double));
            compute_residuals(x, P1_scaled, P2_scaled, n, s0_scaled, c_scaled, m_scaled, g_vec_scaled, res_minus);
            x[i] = orig;
            for (int j = 0; j < dof; ++j)
                J[j * dof + i] = (res_plus[j] - res_minus[j]) / (2 * FD_EPS);
            free(res_plus);
            free(res_minus);
        }

        for (int i = 0; i < dof; ++i) dx[i] = -res[i];
        if (solve_linear_system(J, dx, dof) != 0) return -2;
        double dx_norm = 0.0;
        for (int i = 0; i < dof; ++i) {
            x[i] += dx[i];
            dx_norm += dx[i] * dx[i];
        }
        dx_norm = sqrt(dx_norm);
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