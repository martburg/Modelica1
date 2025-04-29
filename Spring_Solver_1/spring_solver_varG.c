// solve_spring_mass_newton_debug.c - Newton-based solver with debug output
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define MAX_ITER 100
#define TOL 1e-6
#define FD_EPS 1e-6

// Compute spring force between two 3D points
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

// Compute residuals
static void compute_residuals(
    const double *positions_flat, const double *P1, const double *P2, int n,
    double s0, double c, double m, double* g_vec, double *residuals_out)
{
    const double *prev = P1, *next;
    double F_l[3], F_r[3], F_g[3];
    for (int i = 0; i < n - 1; ++i) {
        const double *mid = &positions_flat[i * 3];
        next = (i == n - 2) ? P2 : &positions_flat[(i + 1) * 3];
        compute_spring_force(prev, mid, s0, c, F_l);
        compute_spring_force(next, mid, s0, c, F_r);
        F_g[0] = -m * g_vec[0];
        F_g[1] = -m * g_vec[1];
        F_g[2] = -m * g_vec[2];
        for (int j = 0; j < 3; ++j)
            residuals_out[i * 3 + j] = F_l[j] + F_r[j] + F_g[j];
        prev = mid;
    }
}

// Solve linear system
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
    double c, double total_mass, double* g_vec, double max_sag,
    double* out_positions)
{
    printf("Solve_spring_mass_c called\n");
    int dof = (n - 1) * 3;
    double *x = malloc(dof * sizeof(double));
    double *res = malloc(dof * sizeof(double));
    double *J = malloc(dof * dof * sizeof(double));
    double *dx = malloc(dof * sizeof(double));
    double m = total_mass / (n - 1);
    if (!x || !res || !J || !dx) return -1;

    double dx_line[3] = {P2[0] - P1[0], P2[1] - P1[1], P2[2] - P1[2]};
    double L = sqrt(dx_line[0]*dx_line[0] + dx_line[1]*dx_line[1] + dx_line[2]*dx_line[2]);
    double s0 = (L / n) * 0.9;
    double g_dir[3];

    double g_norm = sqrt(g_vec[0]*g_vec[0] + g_vec[1]*g_vec[1] + g_vec[2]*g_vec[2]);
    for (int j = 0; j < 3; ++j)
        g_dir[j] = g_vec[j] / g_norm;

    printf("=== Initial Guess ===\n");
    for (int i = 1; i < n; ++i) {
        double t = (double)i / n;
        double base[3];
        for (int j = 0; j < 3; ++j)
            base[j] = (1 - t) * P1[j] + t * P2[j];

        double sag = -max_sag * (1 - 4 * t * (1 - t));
        for (int j = 0; j < 3; ++j) {
            x[(i - 1) * 3 + j] = base[j] + sag * g_dir[j];
            printf("Init Mass %d Axis %d: %f\n", i, j, x[(i - 1) * 3 + j]);
        }
    }

    for (int iter = 0; iter < MAX_ITER; ++iter) {
        compute_residuals(x, P1, P2, n, s0, c, m, g_vec, res);

        for (int i = 0; i < dof; ++i) {
            double orig = x[i];
            x[i] = orig + FD_EPS;
            double *res_plus = malloc(dof * sizeof(double));
            compute_residuals(x, P1, P2, n, s0, c, m, g_vec, res_plus);

            x[i] = orig - FD_EPS;
            double *res_minus = malloc(dof * sizeof(double));
            compute_residuals(x, P1, P2, n, s0, c, m, g_vec, res_minus);

            x[i] = orig;
            for (int j = 0; j < dof; ++j)
                J[j * dof + i] = (res_plus[j] - res_minus[j]) / (2 * FD_EPS);

            free(res_plus);
            free(res_minus);
        }

        for (int i = 0; i < dof; ++i)
            dx[i] = -res[i];

        if (solve_linear_system(J, dx, dof) != 0) {
            printf("❌ Linear solver failed!\n");
            return -2;
        }

        double norm = 0.0;
        for (int i = 0; i < dof; ++i) {
            x[i] += dx[i];
            norm += dx[i] * dx[i];
        }

        printf("Iter %3d | Residual Norm = %.8f\n", iter, sqrt(norm));
        if (sqrt(norm) < TOL) {
            printf("✅ Converged after %d iterations.\n", iter);
            break;
        }
    }

    memcpy(out_positions, x, dof * sizeof(double));
    free(x); free(res); free(J); free(dx);
    return 0;
}
