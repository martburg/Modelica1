// solve_spring_mass_newton.c - Newton-based solver with finite-difference Jacobian
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

// Residual calculation for all internal nodes
static void compute_residuals(
    const double *positions_flat, const double *P1, const double *P2, int n,
    double s0, double c, double m, double g, double *residuals_out)
{
    const double *prev = P1, *next;
    double F_l[3], F_r[3], F_g[3];
    for (int i = 0; i < n - 1; ++i) {
        const double *mid = &positions_flat[i * 3];
        next = (i == n - 2) ? P2 : &positions_flat[(i + 1) * 3];
        compute_spring_force(prev, mid, s0, c, F_l);
        compute_spring_force(next, mid, s0, c, F_r);
        F_g[0] = 0.0; F_g[1] = 0.0; F_g[2] = -m * g;
        for (int j = 0; j < 3; ++j)
            residuals_out[i * 3 + j] = F_l[j] + F_r[j] + F_g[j];
        prev = mid;
    }
}

// Basic dense linear solver using Gaussian elimination with partial pivoting
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
    double c, double total_mass, double g, double max_sag,
    double* out_positions)
{
    int dof = (n - 1) * 3;
    double *x = malloc(dof * sizeof(double));
    double *res = malloc(dof * sizeof(double));
    double *J = malloc(dof * dof * sizeof(double));
    double *dx = malloc(dof * sizeof(double));
    double m = total_mass / (n - 1);
    if (!x || !res || !J || !dx) return -1;

    double dx_line[3] = {P2[0] - P1[0], P2[1] - P1[1], P2[2] - P1[2]};
    double L = sqrt(dx_line[0]*dx_line[0] + dx_line[1]*dx_line[1] + dx_line[2]*dx_line[2]);
    double s0 = (L / n) *0.9;

// Initial guess using cubic sag in local coordinate system
for (int i = 1; i < n; ++i) {
    double t = (double)i / n;

    // Base position along the line
    double base[3];
    for (int j = 0; j < 3; ++j)
        base[j] = (1 - t) * P1[j] + t * P2[j];

    // Build local frame: e_z = -gravity, e_x = direction, e_y = e_z × e_x
    double ez[3] = {0.0, 0.0, -1.0};
    double ex[3] = {dx_line[0] / L, dx_line[1] / L, dx_line[2] / L};

    // e_y = e_z × e_x
    double ey[3] = {
        ez[1] * ex[2] - ez[2] * ex[1],
        ez[2] * ex[0] - ez[0] * ex[2],
        ez[0] * ex[1] - ez[1] * ex[0]
    };

    // Normalize e_y
    double ey_norm = sqrt(ey[0]*ey[0] + ey[1]*ey[1] + ey[2]*ey[2]);
    for (int j = 0; j < 3; ++j) ey[j] /= ey_norm;

    // Recompute e_z = e_x × e_y to ensure orthogonality
    ez[0] = ex[1] * ey[2] - ex[2] * ey[1];
    ez[1] = ex[2] * ey[0] - ex[0] * ey[2];
    ez[2] = ex[0] * ey[1] - ex[1] * ey[0];

    // Apply sag in local -Z (global direction of gravity)
    double sag = -max_sag * (1 - 4 * t * (1 - t));
    for (int j = 0; j < 3; ++j)
        x[(i - 1) * 3 + j] = base[j] + sag * ez[j];
}

    for (int iter = 0; iter < MAX_ITER; ++iter) {
        compute_residuals(x, P1, P2, n, s0, c, m, g, res);

        // Jacobian via central difference
        for (int i = 0; i < dof; ++i) {
            double orig = x[i];
            x[i] = orig + FD_EPS;
            double *res_plus = malloc(dof * sizeof(double));
            compute_residuals(x, P1, P2, n, s0, c, m, g, res_plus);
            x[i] = orig - FD_EPS;
            double *res_minus = malloc(dof * sizeof(double));
            compute_residuals(x, P1, P2, n, s0, c, m, g, res_minus);
            x[i] = orig;
            for (int j = 0; j < dof; ++j)
                J[j * dof + i] = (res_plus[j] - res_minus[j]) / (2 * FD_EPS);
            free(res_plus); free(res_minus);
        }

        for (int i = 0; i < dof; ++i) dx[i] = -res[i];
        if (solve_linear_system(J, dx, dof) != 0) return -2;
        double norm = 0.0;
        for (int i = 0; i < dof; ++i) {
            x[i] += dx[i];
            norm += dx[i] * dx[i];
        }
        if (sqrt(norm) < TOL) break;
    }

    memcpy(out_positions, x, dof * sizeof(double));
    free(x); free(res); free(J); free(dx);
    return 0;
}