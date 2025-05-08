#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define MAX_ITER 100
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

double norm3(const double v[3]) {
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

// Cubic parabola arc length calculator
static double arc_length_cubic_parabola(const double P1[3], const double P2[3], const double dir[3], double a, int n) {
    double length = 0.0;
    for (int i = 0; i < n; ++i) {
        double t0 = (double)i / n;
        double t1 = (double)(i + 1) / n;

        double base0[3], base1[3], sagged0[3], sagged1[3];
        for (int j = 0; j < 3; ++j) {
            base0[j] = (1 - t0) * P1[j] + t0 * P2[j];
            base1[j] = (1 - t1) * P1[j] + t1 * P2[j];
        }
        double s0 = -4 * a * t0 * (1 - t0);
        double s1 = -4 * a * t1 * (1 - t1);
        for (int j = 0; j < 3; ++j) {
            sagged0[j] = base0[j] + s0 * dir[j];
            sagged1[j] = base1[j] + s1 * dir[j];
        }
        double d[3] = {
            sagged1[0] - sagged0[0],
            sagged1[1] - sagged0[1],
            sagged1[2] - sagged0[2]
        };
        length += sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
    }
    return length;
}

double vec3_dot(const double v1[3], const double v2[3]) {
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

double estimate_sag_depth(double total_mass, const double g[3], double L, double EA) {
    double weight = total_mass * norm3(g);
    double T0 = EA / L;  // Rough tension scale
    return weight * L / (8 * T0);  // Parabolic sag estimate
}

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

static void report_endpoint_forces(
    const double* P1, const double* P2, const double* x, int n,
    const double* s0_rest, double c, const double* g_vec,
    double m, double scale_force)
{
    double F_P1[3] = {0}, F_P2[3] = {0};

    compute_spring_force(&x[0], P1, s0_rest[0], c, F_P1); // Force ON P1 from x[0]
    compute_spring_force(&x[(n - 2) * 3], P2, s0_rest[n - 1], c, F_P2); // Force ON P2 from x[n-2]

    // 🔁 RESCALE FORCES TO PHYSICAL UNITS
    for (int j = 0; j < 3; ++j) {
        F_P1[j] *= scale_force;
        F_P2[j] *= scale_force;
    }

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

    printf("\n=== Endpoint Force Decomposition ===\n");
    printf("Force on P1: [%f, %f, %f]\n", F_P1[0], F_P1[1], F_P1[2]);
    printf("  | Parallel to gravity:   %+f\n", f1_parallel);
    printf("  | In-plane perpendicular: %+f\n", f1_plane);

    printf("Force on P2: [%f, %f, %f]\n", F_P2[0], F_P2[1], F_P2[2]);
    printf("  | Parallel to gravity:   %+f\n", f2_parallel);
    printf("  | In-plane perpendicular: %+f\n", f2_plane);

    // Net external force (in Newtons now)
    double F_net[3] = {F_P1[0] + F_P2[0], F_P1[1] + F_P2[1], F_P1[2] + F_P2[2]};

    printf("\nNet external force on system (should be near zero): [%e, %e, %e]\n", F_net[0], F_net[1], F_net[2]);
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

void init_dynamic_relaxation(
    double* x, const double* P1, const double* P2,
    int n, double c_scaled,  double m_scaled, const double* g_vec,
    double* s0_relaxed, double L0, double scale_pos)
{
    double L0_scaled = L0 / scale_pos;
    double d[3] = {P2[0] - P1[0], P2[1] - P1[1], P2[2] - P1[2]};
    double L_straight = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);

    double g_norm = sqrt(g_vec[0]*g_vec[0] + g_vec[1]*g_vec[1] + g_vec[2]*g_vec[2]);
    double e_g[3] = {g_vec[0]/g_norm, g_vec[1]/g_norm, g_vec[2]/g_norm};

    // Sag direction is always opposite gravity
    double dir[3] = {-e_g[0], -e_g[1], -e_g[2]};

    // Fit parabola to match L0 using arc length
    double a_low = 0.0, a_high = 1.0, a_mid = 0.0;
    for (int iter = 0; iter < 100; ++iter) {
        a_mid = 0.5 * (a_low + a_high);
        double arc_len = arc_length_cubic_parabola(P1, P2, dir, a_mid, n);
        if (fabs(arc_len - L0_scaled) < 1e-6) break;
        if (arc_len < L0_scaled) a_low = a_mid;
        else a_high = a_mid;
    }

    // Apply parabola sag in gravity direction
    for (int i = 0; i < n - 1; ++i) {
        double t = (double)(i + 1) / n;
        double sag = -4.0 * a_mid * t * (1.0 - t);
        for (int j = 0; j < 3; ++j)
            x[i * 3 + j] = (1 - t) * P1[j] + t * P2[j] + sag * dir[j];
    }
/*
    // Compute s0_relaxed from actual distances in x[] and endpoints
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
        double dx[3] = {
            xi1[0] - xi[0],
            xi1[1] - xi[1],
            xi1[2] - xi[2]
        };
        s0_relaxed[i] = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]) ;
    }*/
    double s0_uniform = L0 / n;
    for (int i = 0; i < n; ++i)
        s0_relaxed[i] = s0_uniform;   
}

static void compute_residuals(
    const double *positions_flat, const double *P1, const double *P2, int n,
    double* s0, double c, double m, const double* g_vec, const double scale_pos, double *residuals_out)
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

        compute_spring_force(prev, mid, s0[i], c, F_l);
        compute_spring_force(next, mid, s0[i+1], c, F_r);

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

void dynamic_relaxation(
    double* x, const double* P1, const double* P2, int n,
    double* s0_relaxed, double c, double m, const double* g_vec,
    double dt_given, int max_steps, const double scale_pos)
{
    int dof = (n - 1) * 3;
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
            const double* prev = (i == 0) ? P1 : &x[(i - 1) * 3];
            const double* mid  = &x[i * 3];
            const double* next = (i == n - 2) ? P2 : &x[(i + 1) * 3];

            double dx1[3] = {mid[0] - prev[0], mid[1] - prev[1], mid[2] - prev[2]};
            double dx2[3] = {mid[0] - next[0], mid[1] - next[1], mid[2] - next[2]};

            double L1 = sqrt(dx1[0]*dx1[0] + dx1[1]*dx1[1] + dx1[2]*dx1[2]);
            double L2 = sqrt(dx2[0]*dx2[0] + dx2[1]*dx2[1] + dx2[2]*dx2[2]);

            for (int j = 0; j < 3; ++j) {
                double F1 = -c * (L1 - s0_relaxed[i]) / (L1 + 1e-12) * dx1[j];
                double F2 = -c * (L2 - s0_relaxed[i + 1]) / (L2 + 1e-12) * dx2[j];
                F[i * 3 + j] += F1 + F2 + m * g_vec[j] - damping * v[i * 3 + j];
            }
        }

        double v_norm_sq = 0.0;
        for (int i = 0; i < dof; ++i) {
            double a = F[i] / m;
            v[i] += a * dt;
            x[i] += v[i] * dt;
            v_norm_sq += v[i] * v[i];
        }

        if (step % 100 == 0) {
            if (damping > damping_target)
                damping *= 0.95;
            if (v_norm_sq < prev_v_norm && dt < 0.1 / omega)
                dt *= 1.02;
            prev_v_norm = v_norm_sq;
        }

        if (v_norm_sq < 1e-10 && step > 1000) {
            damping *= 1.1;
            dt *= 0.9;
        }

        if (v_norm_sq < 1e-8 && step > 5000) {
            memset(v, 0, dof * sizeof(double));
        }

        if (v_norm_sq < 1e-12) {
            printf("\nDynamic relaxation converged at step %d\n\n", step);
            converged = 1;
            break;
        }
    }

    if (!converged)
        printf("\nDynamic relaxation did NOT converge\n\n");
/*
    // Recompute s0_relaxed from final positions
    for (int i = 0; i < n; ++i) {
        double xi[3], xi1[3];
        if (i == 0) {
            memcpy(xi, P1, sizeof(double) * 3);
            memcpy(xi1, &x[0], sizeof(double) * 3);
        } else if (i == n - 1) {
            memcpy(xi, &x[(n - 2) * 3], sizeof(double) * 3);
            memcpy(xi1, P2, sizeof(double) * 3);
        } else {
            memcpy(xi, &x[(i - 1) * 3], sizeof(double) * 3);
            memcpy(xi1, &x[i * 3], sizeof(double) * 3);
        }
        double dx[3] = {
            xi1[0] - xi[0],
            xi1[1] - xi[1],
            xi1[2] - xi[2]
        };
        s0_relaxed[i] = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
    }
*/
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

int newton_refine(
    double* x, const double* P1, const double* P2, int n,
    double* s0, double c, double m, const double* g_vec, double scale_pos)
{
    int dof = (n - 1) * 3;
    double* res = malloc(dof * sizeof(double));
    double* dx = malloc(dof * sizeof(double));
    double* J = malloc(dof * dof * sizeof(double));
    if (!res || !dx || !J) return -1;

    for (int iter = 0; iter < MAX_ITER; ++iter) {
        compute_residuals(x, P1, P2, n, s0, c, m, g_vec, scale_pos, res);

        double res_norm = 0.0;
        for (int i = 0; i < dof; ++i) res_norm += res[i] * res[i];
        res_norm = sqrt(res_norm);

        if (res_norm < TOL) {
            printf("\nNewton converged at iteration %d!\n", iter);
            break;
        }

        // Finite difference Jacobian
        for (int i = 0; i < dof; ++i) {
            double eps = fmax(1e-6, 1e-4 * fabs(x[i]));
            double orig = x[i];
            x[i] = orig + eps;
            double* res_plus = malloc(dof * sizeof(double));
            compute_residuals(x, P1, P2, n, s0, c, m, g_vec, scale_pos,res_plus);
            x[i] = orig - eps;
            double* res_minus = malloc(dof * sizeof(double));
            compute_residuals(x, P1, P2, n, s0, c, m, g_vec,scale_pos, res_minus);
            x[i] = orig;

            for (int j = 0; j < dof; ++j)
                J[j * dof + i] = (res_plus[j] - res_minus[j]) / (2 * eps);

            free(res_plus);
            free(res_minus);
        }

        for (int i = 0; i < dof; ++i) dx[i] = -res[i];

        // Solve J dx = -res using simple Gauss elimination
        // Here you would use solve_linear_system(J, dx, dof)
        // Assume already implemented elsewhere

        // Basic damping
        for (int i = 0; i < dof; ++i)
            dx[i] *= 0.5;

        // Line search
        double alpha = 1.0;
        for (int trial = 0; trial < LS_MAX_TRIALS; ++trial) {
            double* x_trial = malloc(dof * sizeof(double));
            for (int i = 0; i < dof; ++i)
                x_trial[i] = x[i] + alpha * dx[i];

            double* res_trial = malloc(dof * sizeof(double));
            compute_residuals(x_trial, P1, P2, n, s0, c, m, g_vec,scale_pos, res_trial);

            double res_trial_norm = 0.0;
            for (int i = 0; i < dof; ++i)
                res_trial_norm += res_trial[i] * res_trial[i];
            res_trial_norm = sqrt(res_trial_norm);

            if (res_trial_norm < res_norm) {
                memcpy(x, x_trial, dof * sizeof(double));
                free(x_trial);
                free(res_trial);
                break;
            } else {
                alpha *= LS_REDUCTION;
                free(x_trial);
                free(res_trial);
            }
        }
    }

    free(res);
    free(dx);
    free(J);
    return 0;
}

DLL_EXPORT int solve_spring_mass_c(
    double* P1, double* P2,
    int n, double total_mass, double length_factor,
    double rope_diameter, double youngs_modulus,
    double* g_vec, double* out_positions)
{
    printf("P1 = [%f, %f, %f], P2 = [%f, %f, %f]\n", P1[0], P1[1], P1[2], P2[0], P2[1], P2[2]);
    printf("n = %i , Total Mass = %f , Length Factor = %f, Rope Diameter = %f, Young's MOdulus = %f ",n,total_mass, length_factor,rope_diameter,youngs_modulus);
    printf("g_vec = [%f,%f,%f]\n",g_vec[0],g_vec[1],g_vec[2]);

    int dof = (n - 1) * 3;
    double *x = malloc(dof * sizeof(double));
    double *res = malloc(dof * sizeof(double));
    double *J = malloc(dof * dof * sizeof(double));
    double *s0_relaxed = malloc(n * sizeof(double));
    double *s1 = malloc(n * sizeof(double));
    double *dx = malloc(dof * sizeof(double));

    if (!x || !res || !J || !s0_relaxed || !s1 || !dx)
        return -1;

    double dx_line[3] = {P2[0] - P1[0], P2[1] - P1[1], P2[2] - P1[2]};
    double L_straight = sqrt(dx_line[0]*dx_line[0] + dx_line[1]*dx_line[1] + dx_line[2]*dx_line[2]);
    double L0 = L_straight * length_factor;

    // Cross-sectional area A = pi * d^2 / 4
    double A = 3.1515926 * rope_diameter * rope_diameter / 4.0;
    // Stiffness per segment: c = EA / L_seg
    double L_seg = L0 / n;
    double c = (youngs_modulus * A) / L_seg;

    double m = total_mass / (n - 1);

    double g_norm = sqrt(g_vec[0]*g_vec[0] + g_vec[1]*g_vec[1] + g_vec[2]*g_vec[2]);


        // --- Scaling factors ---
    double scale_pos = L_straight;
    double scale_mass = total_mass;
    double scale_force = scale_mass * g_norm;  // = total_mass * |g|

    // --- Apply scaling ---
    double P1_scaled[3] = {P1[0]/scale_pos, P1[1]/scale_pos, P1[2]/scale_pos};
    double P2_scaled[3] = {P2[0]/scale_pos, P2[1]/scale_pos, P2[2]/scale_pos};

    double m_scaled = m / scale_mass;
    double g_vec_scaled[3] = {
        g_vec[0] / scale_force,
        g_vec[1] / scale_force,
        g_vec[2] / scale_force
    };

double c_scaled = c / scale_force;

    init_dynamic_relaxation( x, P1_scaled, P2_scaled, n, c_scaled, m_scaled, g_vec_scaled, s0_relaxed, L0, scale_pos); 

    /*
    printf("\ns0_relaxed BEFORE dynamic relaxation:\n");
    for (int i = 0; i < n ; ++i) {
        printf("Spring %d: %f]\n", i + 1,s0_relaxed[i]);
    }

    printf("\nInitial relaxed positions BEFORE dynamic relaxation:\n");
    for (int i = 0; i < n - 1; ++i) {
        printf("Node %d: [%f, %f, %f]\n", i + 1, x[i*3+0] * scale_pos, x[i*3+1] * scale_pos, x[i*3+2] * scale_pos);
    }
*/

    //dynamic_relaxation(x, P1_scaled, P2_scaled, n, s0_relaxed, c_scaled, m_scaled, g_vec_scaled, 0.001, 1000000,scale_pos);

    /*
    printf("\ns0_relaxed AFTER dynamic relaxation:\n");
    for (int i = 0; i < n ; ++i) {
        printf("Spring %d: %f]\n", i + 1,s0_relaxed[i]);
    }

    printf("\nInitial relaxed positions AFTER dynamic relaxation:\n");

    for (int i = 0; i < n - 1; ++i) {
        printf("Node %d: [%f, %f, %f]\n", i + 1, x[i*3+0] * scale_pos, x[i*3+1] * scale_pos, x[i*3+2] * scale_pos);
    }
/*
    for (int i = 0; i < n ; ++i) {
        s0_relaxed[i] *= scale_pos;
    }
*/

    // rescale to original size
    for (int i = 0; i < dof; ++i)
        x[i] *= scale_pos;
/*
    for (int i = 0; i < n; ++i)
        s0_relaxed[i] *= scale_pos;
*/
    newton_refine(x,P1,P2,n,s0_relaxed,c,m,g_vec,scale_pos );

    // Compute segment lengths and store in s1
    double s1_sum = 0;
    for (int i = 0; i < n; ++i) {
        double xi[3], xi1[3];

        if (i == 0) {
            // First segment: P1 to x[0]
            for (int j = 0; j < 3; ++j) {
                xi[j]  = P1[j];
                xi1[j] = x[j];
            }
        } else if (i == n - 1) {
            // Last segment: x[n-2] to P2
            for (int j = 0; j < 3; ++j) {
                xi[j]  = x[(n - 2) * 3 + j];
                xi1[j] = P2[j];
            }
        } else {
            // Middle segments: x[i-1] to x[i]
            for (int j = 0; j < 3; ++j) {
                xi[j]  = x[(i - 1) * 3 + j];
                xi1[j] = x[i * 3 + j];
            }
        }

        double dx[3] = {
            xi1[0] - xi[0],
            xi1[1] - xi[1],
            xi1[2] - xi[2]
        };
        s1[i] = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
        s1_sum += s1[i];
    }

    memcpy(out_positions, x, dof * sizeof(double));

    printf("Final Positions after init_dynamic_relaxation an newton solver:\n");

    for (int i = 0; i < n - 1; ++i) {
        printf("s0_relaxed %d : %f    ",i,s0_relaxed[i]);
        printf("Node %d: [%f, %f, %f]\n", i + 1, out_positions[i*3+0], out_positions[i*3+1], out_positions[i*3+2]);
    }
    printf("s0_relaxed %d : %f    ",n-1,s0_relaxed[n-1]);

    printf("\nLength = %f, Delta Length = %f]\n", s1_sum, s1_sum - L0);

    report_endpoint_forces(P1, P2, x, n, s0_relaxed, c, g_vec, m, scale_force);

    free(x); free(res); free(J);free(s0_relaxed);free(dx),free(s1);
    return 0;
}
