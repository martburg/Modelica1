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

double init_dynamic_relaxation( double* x,
     const double* P1, const double* P2,
     int n, double c, double m, const double* g_vec, double scale_pos, double* s0_relaxed)
{
    int dof = (n-1) * 3;
    double* v = calloc(dof, sizeof(double));
    double* F = calloc(dof, sizeof(double));
    
    double dx_line[3] = {P2[0] - P1[0], P2[1] - P1[1], P2[2] - P1[2]};
    double L = sqrt(dx_line[0]*dx_line[0] + dx_line[1]*dx_line[1] + dx_line[2]*dx_line[2]);
    double e_rope[3] = {dx_line[0]/L, dx_line[1]/L, dx_line[2]/L};

    double g_norm = sqrt(g_vec[0]*g_vec[0] + g_vec[1]*g_vec[1] + g_vec[2]*g_vec[2]);
    double e_g[3] = {g_vec[0]/g_norm, g_vec[1]/g_norm, g_vec[2]/g_norm};

    double dot_g_rope = e_g[0]*e_rope[0] + e_g[1]*e_rope[1] + e_g[2]*e_rope[2];
    double g_perp[3] = {
        e_g[0] - dot_g_rope * e_rope[0],
        e_g[1] - dot_g_rope * e_rope[1],
        e_g[2] - dot_g_rope * e_rope[2]
    };
    // Project gravity onto plane perpendicular to rope
    double g_perp_norm = sqrt(g_perp[0]*g_perp[0] + g_perp[1]*g_perp[1] + g_perp[2]*g_perp[2]);
    if (g_perp_norm > 1e-12) {
        g_perp[0] /= g_perp_norm;
        g_perp[1] /= g_perp_norm;
        g_perp[2] /= g_perp_norm;
    } else {
        g_perp[0] = g_perp[1] = g_perp[2] = 0.0;
    }

    double s0_stretched = 0;
    for (int i = 0; i < n; ++i) {
        int masses_below = n - i;
        s0_relaxed[i] = (L / n) + (m * g_norm * masses_below) / c;
        s0_stretched += s0_relaxed[i];
    }

    // Compute cosine of angle between e_rope and e_g
    double dot = e_rope[0]*e_g[0] + e_rope[1]*e_g[1] + e_rope[2]*e_g[2];
    dot = fmax(-1.0, fmin(1.0, dot));  // Clamp to valid acos domain

    double angle_rad = acos(dot);
    double angle_deg = angle_rad * 180.0 / 3.1415;

    double a_fit = find_sag_depth(P1, P2, g_vec, s0_stretched, n, 1e-6);
    //printf("a_fit = %.2f : \n", a_fit);
    double sag_factor = 10;

    if (angle_deg < 10.0 || angle_deg > 170.0) {
        //printf("Near-vertical case (%.2f): initializing cubic parabola along gravity\n", angle_deg);
    
        for (int i = 0; i < n - 1; ++i) {
            double t = (double)(i + 1) / (double)n;
    
            // Straight-line base interpolation
            double base[3] = {
                (1 - t) * P1[0] + t * P2[0],
                (1 - t) * P1[1] + t * P2[1],
                (1 - t) * P1[2] + t * P2[2]
            };
    
            // Parabolic sag in gravity direction (s(t) = -4a t(1 - t))
            double sag = a_fit * sag_factor * t * (1.0 - t);
    
            for (int j = 0; j < 3; ++j)
                x[i*3 + j] = base[j] + sag * e_g[j];
        }

        // Compute segment lengths and store in s0_adapted
        double s0_sum = 0;
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
            s0_relaxed[i] = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
            s0_sum += s0_relaxed[i];
        }
/*        double s0_sum1 = 0;
        for (int i = 0; i < n; ++i) {
            s0_relaxed[i] /= s0_sum;
            s0_sum1 += s0_relaxed[i];
        }*/
    }
 else {
        for (int i = 0; i < n; ++i) {
            int masses_below = n - i;
            s0_relaxed[i] = (L / n) ;
        }
        double a = 0.2;

        for (int i = 0; i < n-1; ++i) {
            double t = (double)(i+1) / (double)n;
            for (int j = 0; j < 3; ++j)
                x[i*3+j] = (1.0 - t) * P1[j] + t * P2[j];

            double s = t * L;
            double sag_amount = a * (cosh(s/a) - 1.0);
            //printf("sag_amount (%.6f)\n", sag_amount);
            for (int j = 0; j < 3; ++j)
                x[i*3+j] += sag_amount /100 * g_perp[j];
        }
    }
    free(v);
    free(F);
    return 1;
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

static double compute_energy(
    const double* x, const double* P1, const double* P2, int n,
    double s0, double c, double m, const double* g_vec)
{
    double energy = 0.0;
    const double *prev = P1;

    for (int i = 0; i < n-1; ++i) {
        const double *mid = &x[i*3];
        const double *next = (i == n-2) ? P2 : &x[(i+1)*3];

        // prev -> mid spring
        double delta_l[3] = {mid[0] - prev[0], mid[1] - prev[1], mid[2] - prev[2]};
        double l_l = sqrt(delta_l[0]*delta_l[0] + delta_l[1]*delta_l[1] + delta_l[2]*delta_l[2]);
        if (l_l > s0)
            energy += 0.5 * c * (l_l - s0) * (l_l - s0);

        // Gravity potential energy (mass m at each node)
        energy += m * (-g_vec[0] * mid[0] - g_vec[1] * mid[1] - g_vec[2] * mid[2]);

        prev = mid;
    }

    return energy;
}

static void compute_energy_gradient(
    const double* x, const double* P1, const double* P2, int n,
    double s0, double c, double m, const double* g_vec,
    double* grad_out)
{
    memset(grad_out, 0, (n-1)*3*sizeof(double));

    const double *prev = P1;

    for (int i = 0; i < n-1; ++i) {
        const double *mid = &x[i*3];
        const double *next = (i == n-2) ? P2 : &x[(i+1)*3];

        // prev -> mid spring
        double delta_l[3] = {mid[0] - prev[0], mid[1] - prev[1], mid[2] - prev[2]};
        double l_l = sqrt(delta_l[0]*delta_l[0] + delta_l[1]*delta_l[1] + delta_l[2]*delta_l[2]);
        if (l_l > s0) {
            double coeff_l = c * (1.0 - s0 / l_l);
            for (int j = 0; j < 3; ++j)
                grad_out[i*3+j] += coeff_l * delta_l[j];
        }

        // mid -> next spring
        double delta_r[3] = {next[0] - mid[0], next[1] - mid[1], next[2] - mid[2]};
        double l_r = sqrt(delta_r[0]*delta_r[0] + delta_r[1]*delta_r[1] + delta_r[2]*delta_r[2]);
        if (l_r > s0) {
            double coeff_r = c * (1.0 - s0 / l_r);
            for (int j = 0; j < 3; ++j)
                grad_out[i*3+j] -= coeff_r * delta_r[j];
        }

        // Gravity contribution
        for (int j = 0; j < 3; ++j)
            grad_out[i*3+j] += m * g_vec[j];

        prev = mid;
    }
}

static void dynamic_relaxation(
    double* x, const double* P1, const double* P2, int n,
    double* s0_relaxed, double c, double m, const double* g_vec,
    double dt_given, int max_steps, double scale_pos)
{
    int dof = (n-1) * 3;
    double* v = calloc(dof, sizeof(double));
    double* F = calloc(dof, sizeof(double));
    
    double omega = sqrt(c / m);
    double dt = fmin(dt_given, 0.05 / omega);
    double stop_threshold = 1e-4;

    double damping = 3* sqrt(c / m);

    for (int step = 0; step < max_steps; ++step) {
        memset(F, 0, dof * sizeof(double));

        for (int i = 0; i < n-1; ++i) {
            const double* prev = (i == 0) ? P1 : &x[(i-1)*3];
            const double* mid = &x[i*3];
            const double* next = (i == n-2) ? P2 : &x[(i+1)*3];

            double F_l[3];
            compute_spring_force(prev, mid, s0_relaxed[i], c, F_l);
            for (int j = 0; j < 3; ++j)
                F[i * 3 + j] += F_l[j];

            if (i < n-1){
                double F_r[3]; 
                compute_spring_force(next, mid, s0_relaxed[i+1], c, F_r);
                for (int j = 0; j < 3; ++j)
                    F[i * 3 + j] += F_r[j];
            }            

            for (int j = 0; j < 3; ++j)
                F[i*3+j] +=  m * g_vec[j] - damping * v[i*3+j];
        }
        double v_norm_sq = 0.0;
        for (int i = 0; i < dof; ++i) {
            double a = F[i] / m;
            v[i] += a * dt;
            x[i] += v[i] * dt;
            v_norm_sq += v[i] * v[i];
        }
        if (step % 5000 == 0)
            printf("Relaxation step %d: v_norm = %.6e\n", step, sqrt(v_norm_sq));
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
    double* s0_relaxed = malloc((n) * sizeof(double));
    //double* s0_relaxed = malloc((n) * sizeof(double));
    double m = total_mass / (n - 1);
    if (!x || !res || !J ) return -1;

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

    init_dynamic_relaxation( x, P1_scaled, P2_scaled, n, c_scaled, m_scaled, g_vec_scaled, scale_pos, s0_relaxed); 

    printf("\ns0_relaxed BEFORE dynamic relaxation:\n");
    for (int i = 0; i < n ; ++i) {
        printf("Spring %d: %f]\n", i + 1,s0_relaxed[i]);
    }

    printf("\nInitial relaxed positions BEFORE dynamic relaxation:\n");
    for (int i = 0; i < n - 1; ++i) {
        printf("Node %d: [%f, %f, %f]\n", i + 1, x[i*3+0] * scale_pos, x[i*3+1] * scale_pos, x[i*3+2] * scale_pos);
    }

    dynamic_relaxation(x, P1_scaled, P2_scaled, n, s0_relaxed, c_scaled, m_scaled, g_vec_scaled, 0.001, 100000, scale_pos);

    printf("\ns0_relaxed AFTER dynamic relaxation:\n");
    for (int i = 0; i < n ; ++i) {
        printf("Spring %d: %f]\n", i + 1,s0_relaxed[i]);
    }

    printf("\nInitial relaxed positions AFTER dynamic relaxation:\n");
    for (int i = 0; i < n - 1; ++i) {
        printf("Node %d: [%f, %f, %f]\n", i + 1, x[i*3+0] * scale_pos, x[i*3+1] * scale_pos, x[i*3+2] * scale_pos);
    }
/*
    for (int i = 0; i < n - 1; ++i) {
        double* p1 = (i == 0)     ? P1_scaled           : &x[(i - 1) * 3];
        double* p2 = (i == n - 1) ? (double*)P2_scaled  : &x[i * 3];
    
        double dx[3] = {
            p2[0] - p1[0],
            p2[1] - p1[1],
            p2[2] - p1[2]
        };
        s0_adapted_scaled[i] = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
    }    


    for (int iter = 0; iter < MAX_ITER; ++iter) {
        compute_residuals(x, P1_scaled, P2_scaled, n, s0_scaled, c_scaled, m_scaled, g_vec_scaled, scale_pos, res);
        double res_norm = 0.0;
        for (int i = 0; i < dof; ++i) res_norm += res[i] * res[i];
        res_norm = sqrt(res_norm);

        if (res_norm < TOL) break;

        //printf("Iteration %d: residual norm = %.6e, max residual component = %.6e\n", iter, res_norm, TOL);

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
*/
    for (int i = 0; i < dof; ++i)
        x[i] *= scale_pos;

    memcpy(out_positions, x, dof * sizeof(double));

    printf("Final output positions:\n");

    for (int i = 0; i < n - 1; ++i) {
        printf("Node %d: [%f, %f, %f]\n", i + 1, out_positions[i*3+0], out_positions[i*3+1], out_positions[i*3+2]);
    }


    free(x); free(res); free(J);free(s0_relaxed);
    return 0;
}
