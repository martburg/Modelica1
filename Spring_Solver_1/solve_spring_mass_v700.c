#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define LEGACY_MAX_ITER 100
#define TOL 1e-6
#define LS_REDUCTION 0.5
#define LS_MAX_TRIALS 8
#ifdef _WIN32
#define DLL_EXPORT __declspec(dllexport)
#else
#define DLL_EXPORT
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define IS_INVALID(x) (!(isfinite(x)))
#define ARC_STEPS_FINE 1000
#define NEWTON_TOL     1e-10
#define NEWTON_MAX     20

typedef enum {
    SOLVE_SUCCESS = 0,
    SOLVE_ERROR_NAN = -1,
    SOLVE_ERROR_LINE_SEARCH_FAILED = -2,
    SOLVE_ERROR_JACOBIAN_FAILED = -3,
    SOLVE_ERROR_MAX_ITER = -4,
    SOLVE_ERROR_ALLOC = -5
} SolveStatus;

extern void dgesv_(int* N, int* NRHS, double* A, int* LDA,
    int* IPIV, double* B, int* LDB, int* INFO);

// Recast of spring-mass rope solver as energy minimization (C code)
// Computes total energy and its gradient (i.e., force residual) using BFGS optimizer

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

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

// Compute total potential energy: springs + gravity
double compute_total_energy(
    const double *x, const double *P1, const double *P2, int n,
    const double *s0, double k, double m, const double *g_vec)
{
    double E = 0.0;
    const double *prev = P1;
    for (int i = 0; i < n; ++i) {
        const double *next = (i == n - 1) ? P2 : &x[i * 3];

        double d[3] = { next[0] - prev[0], next[1] - prev[1], next[2] - prev[2] };
        double L = norm3(d);
        double delta = L - s0[i];
        E += 0.5 * k * delta * delta;

        prev = next;
    }
    for (int i = 0; i < n - 1; ++i) {
        const double *xi = &x[i * 3];
        E += -m * (g_vec[0]*xi[0] + g_vec[1]*xi[1] + g_vec[2]*xi[2]);
    }
    return E;
}

// Compute gradient (i.e., force residual) of the energy
void compute_energy_gradient(
    const double *x, const double *P1, const double *P2, int n,
    const double *s0, double k, double m, const double *g_vec,
    double *grad)
{
    memset(grad, 0, sizeof(double) * 3 * (n - 1));

    // First segment: P1 to x[0]
    {
        const double *x0 = &x[0];
        double d[3] = { x0[0] - P1[0], x0[1] - P1[1], x0[2] - P1[2] };
        double L = norm3(d);
        if (L < 1e-8) L = 1e-8;
        double coeff = k * (L - s0[0]) / L;
        for (int j = 0; j < 3; ++j)
            grad[0 * 3 + j] += coeff * d[j];
    }

    // Internal segments
    for (int i = 1; i < n - 1; ++i) {
        const double *xi_prev = &x[(i - 1) * 3];
        const double *xi      = &x[i * 3];

        double d[3] = { xi[0] - xi_prev[0], xi[1] - xi_prev[1], xi[2] - xi_prev[2] };
        double L = norm3(d);
        if (L < 1e-8) L = 1e-8;
        double coeff = k * (L - s0[i]) / L;

        for (int j = 0; j < 3; ++j) {
            grad[(i - 1) * 3 + j] -= coeff * d[j];
            grad[i * 3 + j]       += coeff * d[j];
        }
    }

    // Last segment: x[n-2] to P2
    {
        const double *xn_1 = &x[(n - 2) * 3];
        double d[3] = { P2[0] - xn_1[0], P2[1] - xn_1[1], P2[2] - xn_1[2] };
        double L = norm3(d);
        if (L < 1e-8) L = 1e-8;
        double coeff = k * (L - s0[n - 1]) / L;
        for (int j = 0; j < 3; ++j)
            grad[(n - 2) * 3 + j] -= coeff * d[j];
    }

    // Gravity
    for (int i = 0; i < n - 1; ++i) {
        grad[i * 3 + 0] += m * g_vec[0];
        grad[i * 3 + 1] += m * g_vec[1];
        grad[i * 3 + 2] += m * g_vec[2];
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

    // Normalize g_vec
    double g_vec_unit[3];
    double g_norm = sqrt(g_vec[0]*g_vec[0] + g_vec[1]*g_vec[1] + g_vec[2]*g_vec[2]);
    g_vec_unit[0] = g_vec[0] / g_norm;
    g_vec_unit[1] = g_vec[1] / g_norm;
    g_vec_unit[2] = g_vec[2] / g_norm;

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

// Full BFGS optimizer for energy minimization
int energy_minimize(
    double *x, double *x_out, const double *P1, const double *P2, int n,
    const double *s0, double k, double m, const double *g_vec,
    double step_size, int max_iter, double tol)
{
    int dof = 3 * (n - 1);
    memcpy(x_out, x, sizeof(double) * dof);

    double *grad      = malloc(dof * sizeof(double));
    double *grad_prev = malloc(dof * sizeof(double));
    double *s         = malloc(dof * sizeof(double));
    double *y         = malloc(dof * sizeof(double));
    double *H         = calloc(dof * dof, sizeof(double));
    double *p         = malloc(dof * sizeof(double));
    double *x_prev    = malloc(dof * sizeof(double));
    double *Hy        = malloc(dof * sizeof(double));
    double *grad_trial = malloc(sizeof(double) * dof);
    if (!grad || !grad_prev || !s || !y || !H || !p || !x_prev || !Hy)
        return -1;

    for (int i = 0; i < dof; ++i)
        H[i * dof + i] = 1.0;

    compute_energy_gradient(x_out, P1, P2, n, s0, k, m, g_vec, grad);
    double prev_energy = compute_total_energy(x_out, P1, P2, n, s0, k, m, g_vec);

    for (int iter = 0; iter < max_iter; ++iter) {
        double grad_norm = 0.0;
        for (int i = 0; i < dof; ++i) {
            if (!isfinite(grad[i])) {
                printf("NaN in gradient at iter %d, index %d\n", iter, i);
                free(grad); free(grad_prev); free(s); free(y); free(H);
                free(p); free(x_prev); free(Hy);
                return SOLVE_ERROR_NAN;
            }
            grad_norm += grad[i] * grad[i];
        }
        grad_norm = sqrt(grad_norm);
    
        if (grad_norm < tol) {
            printf("\nBFGS CONVERGED 3D in %d steps with grad norm %.3e\n", iter, grad_norm);
            free(grad); free(grad_prev); free(s); free(y); free(H);
            free(p); free(x_prev); free(Hy);
            return 0;
        }
    
        // Compute p = -H * grad
        for (int i = 0; i < dof; ++i) {
            p[i] = 0.0;
            for (int j = 0; j < dof; ++j)
                p[i] -= H[i * dof + j] * grad[j];
        }
    
        // Save previous state
        memcpy(x_prev, x_out, sizeof(double) * dof);
        memcpy(grad_prev, grad, sizeof(double) * dof);
    
        // Compute x_trial = x_out + step_size * p
        double *x_trial = malloc(sizeof(double) * dof);
        for (int i = 0; i < dof; ++i) {
            x_trial[i] = x_out[i] + step_size * p[i];
            if (!isfinite(x_trial[i])) {
                //printf("NaN in x_trial at iter %d, index %d\n", iter, i);
                free(grad); free(grad_prev); free(s); free(y); free(H);
                free(p); free(x_prev); free(Hy);
                return SOLVE_ERROR_NAN;
            }
        }
    
        // Compute trial gradient and energy
        compute_energy_gradient(x_trial, P1, P2, n, s0, k, m, g_vec, grad_trial);
        double energy = compute_total_energy(x_trial, P1, P2, n, s0, k, m, g_vec);
        
        if (!isfinite(energy)) {
            //printf("NaN in energy at iter %d\n", iter);
            free(x_trial);
            free(grad); free(grad_prev); free(s); free(y); free(H);
            free(p); free(x_prev); free(Hy);
            return SOLVE_ERROR_NAN;
        }
        for (int i = 0; i < dof; ++i) {
            if (!isfinite(grad_trial[i])) {
                //printf("NaN in trial gradient at iter %d, index %d\n", iter, i);
                free(x_trial);
                free(grad); free(grad_prev); free(s); free(y); free(H);
                free(p); free(x_prev); free(Hy);
                return SOLVE_ERROR_NAN;
            }
        }
    
        if (fabs(energy - prev_energy) < tol * 1e-2 || energy < prev_energy) {
            prev_energy = energy;
            memcpy(x_out, x_trial, sizeof(double) * dof);
            memcpy(grad, grad_trial, sizeof(double) * dof);
            free(x_trial);
        } else {
            //printf("Rejected step at iter %d: energy increased\n", iter);
            free(grad); free(grad_prev); free(s); free(y); free(H);
            free(p); free(x_prev); free(Hy);
            return 0;
        }
    
        prev_energy = energy;
        memcpy(x_out, x_trial, sizeof(double) * dof);
        free(x_trial);
    
        // Update s, y
        for (int i = 0; i < dof; ++i) {
            s[i] = x_out[i] - x_prev[i];
            y[i] = grad[i] - grad_prev[i];
        }
    
        double sTy = 0.0;
        for (int i = 0; i < dof; ++i) {
            if (!isfinite(s[i]) || !isfinite(y[i])) {
                printf("NaN in s or y at iter %d, index %d\n", iter, i);
                free(grad); free(grad_prev); free(s); free(y); free(H);
                free(p); free(x_prev); free(Hy);
                return SOLVE_ERROR_NAN;
            }
            sTy += s[i] * y[i];
        }
    
        if (fabs(sTy) < 1e-12) {
            for (int i = 0; i < dof; ++i)
                x_out[i] -= step_size * grad[i];  // fallback
            compute_energy_gradient(x_out, P1, P2, n, s0, k, m, g_vec, grad);
            continue;
        }
    
        double rho = 1.0 / sTy;
    
        // Hy = H * y
        for (int i = 0; i < dof; ++i) {
            Hy[i] = 0.0;
            for (int j = 0; j < dof; ++j)
                Hy[i] += H[i * dof + j] * y[j];
        }
    
        for (int i = 0; i < dof; ++i)
            for (int j = 0; j < dof; ++j) {
                double ssT = s[i] * s[j];
                double Hy_sT = Hy[i] * s[j];
                double s_HyT = s[i] * Hy[j];
                H[i * dof + j] += rho * (ssT - Hy_sT - s_HyT);
            }
    }

    printf("\nBFGS did NOT CONVERGE after %d iterations\n", max_iter);

cleanup:
    free(grad); free(grad_prev); free(s); free(y); free(H);
    free(p); free(x_prev); free(Hy);
    return SOLVE_ERROR_MAX_ITER;
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

static void report_endpoint_forces(
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

    printf("\n=== Endpoint Force Decomposition (Spring elongation Newton) ===\n");
    printf("|P1 F| = %.6f, |F_hor| = %.6f, |F_ver| = %.6f, F_P1: [%f, %f, %f]\n",norm3(F_P1),f1_parallel,f1_plane, F_P1[0], F_P1[1], F_P1[2]);
    printf("|P2 F| = %.6f, |F_hor| = %.6f, |F_ver| = %.6f, F_P2: [%f, %f, %f]\n",norm3(F_P2),f2_parallel,f2_plane, F_P2[0], F_P2[1], F_P2[2]);

    // Net external force (in Newtons now)
    double F_net[3] = {F_P1[0] + F_P2[0], F_P1[1] + F_P2[1], F_P1[2] + F_P2[2]};

    double TM[3];
    TM[0] = total_mass * g_vec[0];
    TM[1] = total_mass * g_vec[1];
    TM[2] = total_mass * g_vec[2];
    double TMG = norm3(TM);

    printf("Delta F (F_net - m*g = %e\n", TMG -norm3(F_net)); 
    
    F_P1_out = F_P1;
    F_P2_out = F_P2;   
}

void report_endpoint_forces_weight_partitioned(
    const double* P1, const double* P2, const double* x, int n,
    const double* g_vec, double total_mass,
    double* F_P1_out, double* F_P2_out)
{
    int lowest_index = -1;
    double min_z = 1e20;

    // Find the node with the lowest height (assume gravity is along z)
    for (int i = 0; i < n - 1; ++i) {
        if (x[i * 3 + 2] < min_z) {
            min_z = x[i * 3 + 2];
            lowest_index = i;
        }
    }

    double m = total_mass / (n - 1);
    double g_mag = norm3(g_vec);
    double weight_P1 = m * (lowest_index + 1);        // Includes up to and including lowest node
    double weight_P2 = m * ((n - 1) - (lowest_index)); // Remaining nodes

    double total = weight_P1 + weight_P2;
    double f1_frac = weight_P1 / total;
    double f2_frac = weight_P2 / total;

    // Compute force vectors
    double P1_g_vec[3], P2_g_vec[3];
    for (int j = 0; j < 3; ++j) {
        P1_g_vec[j] = f1_frac * total_mass * g_vec[j];
        P2_g_vec[j] = f2_frac * total_mass * g_vec[j];
    }

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
    // Flip tangent direction if it's not pointing along g_vec
    double g_unit[3] = {
    g_vec[0] / g_mag,
    g_vec[1] / g_mag,
    g_vec[2] / g_mag
    };
    if (vec3_dot(t1, g_unit) < 0) {
        for (int j = 0; j < 3; ++j) t1[j] *= -1;
    }
    if (vec3_dot(t2, g_unit) < 0) {
        for (int j = 0; j < 3; ++j) t2[j] *= -1;
    }

    // Compute cosine of angle between g_vec and tangent
    double cos1 = vec3_dot(t1, g_vec) / (norm3(g_vec));
    double cos2 = vec3_dot(t2, g_vec) / (norm3(g_vec));

    // Scale to get full force vector along tangent
    double F1_mag = norm3(P1_g_vec) / fabs(cos1);
    double F2_mag = norm3(P2_g_vec) / fabs(cos2);

    double F_P1[3], F_P2[3];
    for (int j = 0; j < 3; ++j) {
        F_P1[j] = F1_mag * t1[j];
        F_P2[j] = F2_mag * t2[j];
    }
    // Normalize gravity vector
    double g_unit_vec[3] = {
        g_vec[0] / g_mag,
        g_vec[1] / g_mag,
        g_vec[2] / g_mag
    };

    // Split F_P1
    double F_P1_parallel[3], F_P1_perp[3];
    project_parallel_and_perpendicular(F_P1, g_unit_vec, F_P1_parallel, F_P1_perp);

    // Split F_P2
    double F_P2_parallel[3], F_P2_perp[3];
    project_parallel_and_perpendicular(F_P2, g_unit_vec, F_P2_parallel, F_P2_perp);

    printf("\n=== Endpoint Force Decomposition (Weight) ===\n");
    printf("|P1 F| = %.6f, |F_hor| = %.6f, |F_ver| = %.6f, F_P1: [%f, %f, %f]\n",norm3(F_P1),norm3(F_P1_parallel),norm3(F_P1_perp), F_P1[0], F_P1[1], F_P1[2]);
    printf("|P2 F| = %.6f, |F_hor| = %.6f, |F_ver| = %.6f, F_P2: [%f, %f, %f]\n",norm3(F_P2),norm3(F_P2_parallel),norm3(F_P2_perp), F_P2[0], F_P2[1], F_P2[2]);

    double F_net[3] = {
        F_P1[0] + F_P2[0],
        F_P1[1] + F_P2[1],
        F_P1[2] + F_P2[2]
    };
    double TM[3];
        TM[0] = total_mass * g_vec[0];
        TM[1] = total_mass * g_vec[1];
        TM[2] = total_mass * g_vec[2];
    double TMG = norm3(TM);

    printf("Delta F (F_net - m*g = %e\n", TMG -norm3(F_net));

    F_P1_out = F_P1;
    F_P2_out = F_P2;
}

void lin_comb3(double out[3], const double a, const double v1[3], const double b, const double v2[3]) {
    for (int i = 0; i < 3; ++i)
        out[i] = a * v1[i] + b * v2[i];
}

double arc_length_sagged(
    const double P1[3], const double P2[3], const double e_g[3],
    double a, int n)
{
    double length = 0.0;
    for (int i = 0; i < n; ++i) {
        double t0 = (double)i / n;
        double t1 = (double)(i + 1) / n;

        double base0[3], base1[3], sagged0[3], sagged1[3];

        // base points
        lin_comb3(base0, 1.0 - t0, P1, t0, P2);
        lin_comb3(base1, 1.0 - t1, P1, t1, P2);

        // sag values
        double sag0 = -4.0 * a * t0 * (1.0 - t0);
        double sag1 = -4.0 * a * t1 * (1.0 - t1);

        // sagged positions
        for (int j = 0; j < 3; ++j) {
            sagged0[j] = base0[j] + sag0 * e_g[j];
            sagged1[j] = base1[j] + sag1 * e_g[j];
        }

        // segment length
        double seg[3] = {
            sagged1[0] - sagged0[0],
            sagged1[1] - sagged0[1],
            sagged1[2] - sagged0[2]
        };
        length += norm3(seg);
    }
    return length;
}

void vec_add3(double out[3], const double v1[3], const double v2[3]) {
    for (int i = 0; i < 3; ++i)
        out[i] = v1[i] + v2[i];
}

double find_sag_depth(
    const double P1[3], const double P2[3], const double g_vec[3],
    double s0_stretched, int n, double tol)
{
    // Normalize gravity
    double g_norm = norm3(g_vec);
    double e_g[3] = { g_vec[0]/g_norm, g_vec[1]/g_norm, g_vec[2]/g_norm };

    // Straight-line length
    double dx[3] = { P2[0] - P1[0], P2[1] - P1[1], P2[2] - P1[2] };
    double L = norm3(dx);

    // Bisection bounds
    double a_low = 0.0;
    double a_high = 0.5 * L;
    double a_mid;

    for (int iter = 0; iter < 100; ++iter) {
        a_mid = 0.5 * (a_low + a_high);
        double len = arc_length_sagged(P1, P2, e_g, a_mid, n);

        if (fabs(len - s0_stretched) < tol)
            break;
        if (len < s0_stretched)
            a_low = a_mid;
        else
            a_high = a_mid;
    }

    return a_mid;
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
    double* s0_post_init, double L0, double scale_pos)
{
    double L0_scaled = L0 / scale_pos;

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

double dynamic_relaxation(
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
        //if (step % 10000 == 0){
        //    printf("v_norm_sq %.6e\n",v_norm_sq);
        //} 

        if (v_norm_sq < 1e-8) {
            printf("\nDynamic relaxation CONVERGED at step %d with v_norm_sq %.3e\n", step, v_norm_sq);
            converged = 1;
            break;
        }
    }
    free(v);
    free(F);
    if (converged == 0){
        printf("\nDynamic relaxation didn NOT CONVERGE after %d iterations\n", max_steps);
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
    if (!res || !J || !ipiv || !dx) return -1;

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
        if (info != 0) {
            //printf("dgesv failed (info = %d)\n", info);
            free(res); free(J); free(ipiv); free(dx);
            return SOLVE_ERROR_JACOBIAN_FAILED;
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
            //printf("Wolfe line search failed; using fallback damping\n");
            for (int i = 0; i < dof; ++i)
            x_out[i] += 0.25 * dx[i];  // damped step
        }

        // Compute step size
        double norm_dx = 0.0;
        for (int i = 0; i < dof; ++i)
            norm_dx += dx[i] * dx[i];
        norm_dx = sqrt(norm_dx);

        if (norm_dx < tol) {
            printf("3D Newton CONVERGED in %d steps with step norm %.3e\n", iter, norm_dx);
            free(res); free(J); free(ipiv); free(dx);
            return 0;
        }
    }

    //printf("3D Newton did NOT CONVERGE after %d iterations\n", max_iter);
    free(res); free(J); free(ipiv); free(dx);
    return SOLVE_ERROR_MAX_ITER;
}

DLL_EXPORT int solve_spring_mass_c(
    double* P1, double* P2,
    int n, double total_mass, double length_factor,
    double rope_diameter, double youngs_modulus,
    double* g_vec, double* out_positions)
{
    printf("P1 = [%f, %f, %f], P2 = [%f, %f, %f]\n", P1[0], P1[1], P1[2], P2[0], P2[1], P2[2]);
    printf("n = %i , Total Mass = %f , Length Factor = %f, Rope Diameter = %f, Young's MOdulus = %f ", n, total_mass, length_factor, rope_diameter, youngs_modulus);
    printf("g_vec = [%f,%f,%f]\n\n", g_vec[0], g_vec[1], g_vec[2]);

    int dof = (n - 1) * 3;
    int status;
    double *x = malloc(dof * sizeof(double));
    double *x_init = malloc(dof * sizeof(double));
    double *x_relaxed = malloc(dof * sizeof(double));
    double *x_energy = malloc(dof * sizeof(double));
    double *x_newton_3d = malloc(dof * sizeof(double));
    double *s0_init = malloc(n * sizeof(double));
    double *s0_init_scaled = malloc(n * sizeof(double));
    double *s0_post_init = malloc(n * sizeof(double));
    double *s0_post_relax = malloc(n * sizeof(double));
    double *s0_post_energy = malloc(n * sizeof(double));
    double *s0_post_newton_3d = malloc(n * sizeof(double));
    double *dx = malloc(dof * sizeof(double));
    double F_P1_out_w[3];
    double F_P2_out_w[3];
    double F_P1_out_n[3];
    double F_P2_out_n[3];


    if (!x || !x_init || !x_relaxed || !x_energy || !x_newton_3d ||!s0_init || !s0_init_scaled ||
         !s0_post_init || !s0_post_relax || !s0_post_newton_3d || !s0_post_energy || !dx )
        return -1;

    double dx_line[3] = {P2[0] - P1[0], P2[1] - P1[1], P2[2] - P1[2]};
    double L_straight = sqrt(dx_line[0]*dx_line[0] + dx_line[1]*dx_line[1] + dx_line[2]*dx_line[2]);
    double L0 = L_straight * length_factor;

    // Check if rope is nearly aligned with gravity
    double rope_dir[3] = {dx_line[0], dx_line[1], dx_line[2]};
    double rope_len = L_straight;
    for (int i = 0; i < 3; ++i) rope_dir[i] /= rope_len;

    double gravity_unit[3] = {g_vec[0], g_vec[1], g_vec[2]};
    double g_len = norm3(gravity_unit);
    for (int i = 0; i < 3; ++i) gravity_unit[i] /= g_len;

    double cos_theta = fabs(vec3_dot(rope_dir, gravity_unit));
    if (cos_theta > 0.999) {
        printf("Rope is nearly aligned with gravity (cosθ = %.5f). Numerical issues may occur.\n", cos_theta);
}

    for (int i = 0; i < n; ++i) {
        s0_init[i] = L0 / n;
        s0_init_scaled[i] = length_factor / n;
    }

    // Cross-sectional area A = pi * d^2 / 4
    double A = 3.1515926 * rope_diameter * rope_diameter / 4.0;
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

    // --- Initialize rope with sag ---
    double s0_post_init_sum = init_dynamic_relaxation(
        x, P1_scaled, P2_scaled, n, g_vec_scaled, s0_post_init, L0, scale_pos);

    for (int i = 0; i < dof; ++i) {
        x_init[i] = x[i] * scale_pos;
        }
    printf("\nPositions initiated from parabola:\n");
    for (int i = 0; i < n - 1; ++i) {
        printf("s0_post_init %-2d : %-10.6f    ", i, s0_post_init[i] * scale_pos);
        printf("Node %-3d: [%-12.6f, %12.6f, %12.6f]\n", i + 1,
            x[i*3+0] * scale_pos, x[i*3+1] * scale_pos, x[i*3+2] * scale_pos);
    }
    printf("s0_post_init %-3d : %f    \n",n-1,s0_post_init[n-1] * scale_pos);    
    printf("Initial Length = %f, Delta Length = %f]\n", s0_post_init_sum * scale_pos, s0_post_init_sum * scale_pos - L0);

    // --- Relax dynamically ---

    status = dynamic_relaxation(
        x, x_relaxed, P1_scaled, P2_scaled, n, s0_init_scaled, c_scaled, m_scaled, g_vec_scaled,0.001, 10000000, scale_pos);
    if (status != 0) {
        printf("\n Dynamic relaxation failed\n");
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
    printf("\nPositions after dynamic relaxation:\n");
    for (int i = 0; i < n - 1; ++i) {
        printf("s0_post_relax %-3d : %-10.6f    ", i, s0_post_relax[i]);
        printf("Node %-3d: [%12.6f, %12.6f, %12.6f]\n", i + 1,
            x_relaxed[i*3+0], x_relaxed[i*3+1], x_relaxed[i*3+2]);
    }
    printf("s0_post_relax %d : %f    \n",n-1,s0_post_relax[n-1]);   
    printf("Relaxed Length = %f, Delta Length = %f]\n\n", s0_post_relax_sum, s0_post_relax_sum - L0);

    // --- Newton 3D solver ---

    status = analytic_newton_solver_3d(x_relaxed, x_newton_3d, P1, P2, n, s0_post_relax, c, m, g_vec, 1, 20000, 1e-6);
    if (status != 0) {
        printf("\n 3D Newton solver failed\n");
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
    printf("\nPositions after 3D Newton :\n");
    for (int i = 0; i < n - 1; ++i) {
        printf("s0_post_3D_Newton %-3d : %-10.6f    ", i, s0_post_newton_3d[i] );
        printf("Node %-3d: [%12.6f, %12.6f, %12.6f]\n", i + 1,
            x_newton_3d[i*3+0], x_newton_3d[i*3+1], x_newton_3d[i*3+2] );
    }
    printf("s0_post_3D_Newton %-3d : %-10.6f    \n", n - 1, s0_post_newton_3d[n - 1]);

    printf("3D Newton Length = %f, Delta Length = %f]\n", s0_newton_3d_sum, s0_newton_3d_sum - L0); 
    
    // --- ENERGY solver ---    

    status = energy_minimize(x_relaxed, x_energy, P1, P2, n,s0_post_relax, c,m, g_vec,1, 20000, 1e-4);
    if (status != 0) {
        printf("\n Energy solver failed\n");
    }
    // Compute Energy segment lengths
    double s0_energy_sum = 0;
    for (int i = 0; i < n; ++i) {
        double xi[3], xi1[3];
        if (i == 0) {
            for (int j = 0; j < 3; ++j) {
                xi[j]  = P1[j];
                xi1[j] = x_energy[j];
            }
        } else if (i == n - 1) {
            for (int j = 0; j < 3; ++j) {
                xi[j]  = x_energy[(n - 2) * 3 + j];
                xi1[j] = P2[j];
            }
        } else {
            for (int j = 0; j < 3; ++j) {
                xi[j]  = x_energy[(i - 1) * 3 + j];
                xi1[j] = x_energy[i * 3 + j];
            }
        }
        double dx[3] = {
            xi1[0] - xi[0], xi1[1] - xi[1], xi1[2] - xi[2]
        };
        s0_post_energy[i] = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
        s0_energy_sum += s0_post_energy[i];
    }    
    printf("\nPositions after Energy):\n");
    for (int i = 0; i < n - 1; ++i) {
        printf("s0_post_Energy %-3d : %-10.6f    ", i, s0_post_energy[i]);
        printf("Node %-3d: [%12.6f, %12.6f, %12.6f]\n", i + 1,
            x_energy[i*3+0], x_energy[i*3+1], x_energy[i*3+2]);
    }
    printf("s0_post_Energy %-3d : %-10.6f\n", n - 1, s0_post_energy[n - 1]);
   
    printf("Energy Length = %f, Delta Length = %f\n\n", s0_energy_sum, s0_energy_sum - L0);

    // --- chose output ---

    memcpy(out_positions, x_energy, dof * sizeof(double));

    report_endpoint_forces(P1, P2, x_newton_3d, n, s0_post_relax, c, g_vec, m, scale_force,total_mass, F_P1_out_n, F_P2_out_n);

    report_endpoint_forces_weight_partitioned(P1,P2,x_energy,n,g_vec,total_mass, F_P1_out_w, F_P2_out_w);

    free(x); free(x_init); free(x_relaxed); free(x_energy);
    free(x_newton_3d); free(s0_init); free(s0_init_scaled); free(s0_post_init); free(s0_post_relax);
    free(s0_post_newton_3d); free(s0_post_energy);free(dx);
    return 0;
}
