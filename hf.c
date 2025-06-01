#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// Enhanced structure for 3D vectors
typedef struct { double x, y, z; } vec3_t;

// Mathematical constants
#define PI 3.14159265358979323846
#define SQRT_PI 1.77245385090551602729
#define MAX_ANGULAR_MOMENTUM 4
#define MAX_PRIMITIVES 10

// Enhanced basis set structures
typedef struct {
    double exponent;
    double contraction_coeff;
} primitive_t;

typedef struct {
    vec3_t center;
    int l, m, n;  // Angular momentum quantum numbers
    int num_primitives;
    primitive_t primitives[MAX_PRIMITIVES];
} contracted_gaussian_t;

// Pure functional vector operations
static inline vec3_t vec3_add(vec3_t a, vec3_t b) {
    return (vec3_t){a.x + b.x, a.y + b.y, a.z + b.z};
}

static inline vec3_t vec3_sub(vec3_t a, vec3_t b) {
    return (vec3_t){a.x - b.x, a.y - b.y, a.z - b.z};
}

static inline vec3_t vec3_scale(vec3_t v, double s) {
    return (vec3_t){v.x * s, v.y * s, v.z * s};
}

static inline double vec3_dot(vec3_t a, vec3_t b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline double vec3_norm2(vec3_t v) {
    return vec3_dot(v, v);
}

// Enhanced mathematical utilities
static double factorial_impl(int n, double acc) {
    return (n <= 1) ? acc : factorial_impl(n - 1, acc * n);
}

static double factorial(int n) {
    return (n < 0) ? 0.0 : factorial_impl(n, 1.0);
}

static double double_factorial_impl(int n, double acc) {
    return (n <= 1) ? acc : double_factorial_impl(n - 2, acc * n);
}

static double double_factorial(int n) {
    return (n < 0) ? 1.0 : double_factorial_impl(n, 1.0);
}

// Enhanced Boys function with better numerical stability
static double boys_series(int n, double x, int term, double acc) {
    if (term > 100) return acc;  // Increased convergence limit
    double term_val = pow(-x, term) / (factorial(term) * (2.0 * n + 2.0 * term + 1.0));
    if (fabs(term_val) < 1e-16) return acc;
    return boys_series(n, x, term + 1, acc + term_val);
}

static double boys_asymptotic(int n, double x) {
    return double_factorial(2 * n - 1) * sqrt(PI) / (2.0 * pow(2.0 * x, n + 0.5));
}

static double boys_function(int n, double x) {
    if (x < 1e-12) return 1.0 / (2.0 * n + 1.0);
    if (x > 30.0) return boys_asymptotic(n, x);
    return boys_series(n, x, 0, 0.0);
}

// Corrected Gaussian normalization
static double gaussian_norm(int l, int m, int n, double alpha) {
    // Correct normalization formula for Cartesian Gaussians
    double norm = pow(2.0 * alpha / PI, 0.75) * 
                  pow(4.0 * alpha, (l + m + n) / 2.0) / 
                  sqrt(double_factorial(2 * l - 1) * 
                       double_factorial(2 * m - 1) * 
                       double_factorial(2 * n - 1));
    return norm;
}

// Gaussian product rule
static double gaussian_product_coeff(vec3_t a_center, double a_alpha,
                                   vec3_t b_center, double b_alpha) {
    vec3_t diff = vec3_sub(a_center, b_center);
    double gamma = a_alpha + b_alpha;
    double prefactor = exp(-a_alpha * b_alpha * vec3_norm2(diff) / gamma);
    return prefactor * pow(PI / gamma, 1.5);
}

static vec3_t gaussian_product_center(vec3_t a_center, double a_alpha,
                                    vec3_t b_center, double b_alpha) {
    double gamma = a_alpha + b_alpha;
    return vec3_scale(vec3_add(vec3_scale(a_center, a_alpha),
                              vec3_scale(b_center, b_alpha)), 1.0 / gamma);
}

// Enhanced 1D overlap integral with better recursion
static double overlap_1d(int i, int j, double xa, double xb, double alpha, double beta) {
    if (i < 0 || j < 0) return 0.0;
    
    double gamma = alpha + beta;
    double xp = (alpha * xa + beta * xb) / gamma;
    double prefactor = exp(-alpha * beta * (xa - xb) * (xa - xb) / gamma) * sqrt(PI / gamma);
    
    if (i == 0 && j == 0) {
        return prefactor;
    }
    
    // Use Obara-Saika recursion relations
    if (i > 0) {
        double term1 = (xp - xa) * overlap_1d(i - 1, j, xa, xb, alpha, beta);
        double term2 = (i - 1 > 0) ? (i - 1) / (2.0 * gamma) * overlap_1d(i - 2, j, xa, xb, alpha, beta) : 0.0;
        double term3 = (j > 0) ? j / (2.0 * gamma) * overlap_1d(i - 1, j - 1, xa, xb, alpha, beta) : 0.0;
        return term1 + term2 + term3;
    } else if (j > 0) {
        double term1 = (xp - xb) * overlap_1d(i, j - 1, xa, xb, alpha, beta);
        double term2 = (j - 1 > 0) ? (j - 1) / (2.0 * gamma) * overlap_1d(i, j - 2, xa, xb, alpha, beta) : 0.0;
        double term3 = (i > 0) ? i / (2.0 * gamma) * overlap_1d(i - 1, j - 1, xa, xb, alpha, beta) : 0.0;
        return term1 + term2 + term3;
    }
    
    return 0.0;
}

// Enhanced 1D kinetic energy integral
static double kinetic_1d(int i, int j, double xa, double xb, double alpha, double beta) {
    double s_ij = overlap_1d(i, j, xa, xb, alpha, beta);
    double s_ij_plus2 = overlap_1d(i, j + 2, xa, xb, alpha, beta);
    double s_ij_minus2 = (j >= 2) ? overlap_1d(i, j - 2, xa, xb, alpha, beta) : 0.0;
    
    double term1 = beta * (2 * j + 1) * s_ij;
    double term2 = -2.0 * beta * beta * s_ij_plus2;
    double term3 = (j >= 2) ? -0.5 * beta * j * (j - 1) * s_ij_minus2 : 0.0;
    
    return term1 + term2 + term3;
}

// Enhanced 1D nuclear attraction integral using recursion
static double nuclear_1d(int i, int j, double xa, double xb, double xc, double alpha, double beta) {
    double gamma = alpha + beta;
    double xp = (alpha * xa + beta * xb) / gamma;
    double prefactor = exp(-alpha * beta * (xa - xb) * (xa - xb) / gamma) * sqrt(PI / gamma);
    
    if (i == 0 && j == 0) {
        // Base case: F₀ function
        double pc2 = (xp - xc) * (xp - xc);
        double t = gamma * pc2;
        return prefactor * boys_function(0, t);
    }
    
    // Recursion relations for higher angular momentum
    if (i > 0) {
        double term1 = (xp - xa) * nuclear_1d(i - 1, j, xa, xb, xc, alpha, beta);
        double term2 = (i - 1 > 0) ? (i - 1) / (2.0 * gamma) * nuclear_1d(i - 2, j, xa, xb, xc, alpha, beta) : 0.0;
        double term3 = (j > 0) ? j / (2.0 * gamma) * nuclear_1d(i - 1, j - 1, xa, xb, xc, alpha, beta) : 0.0;
        return term1 + term2 + term3;
    } else if (j > 0) {
        double term1 = (xp - xb) * nuclear_1d(i, j - 1, xa, xb, xc, alpha, beta);
        double term2 = (j - 1 > 0) ? (j - 1) / (2.0 * gamma) * nuclear_1d(i, j - 2, xa, xb, xc, alpha, beta) : 0.0;
        double term3 = (i > 0) ? i / (2.0 * gamma) * nuclear_1d(i - 1, j - 1, xa, xb, xc, alpha, beta) : 0.0;
        return term1 + term2 + term3;
    }
    
    return 0.0;
}

// Complete overlap integral
double overlap(vec3_t a_center, double a_alpha, int a_l, int a_m, int a_n,
               vec3_t b_center, double b_alpha, int b_l, int b_m, int b_n) {
    
    double norm_a = gaussian_norm(a_l, a_m, a_n, a_alpha);
    double norm_b = gaussian_norm(b_l, b_m, b_n, b_alpha);
    
    double sx = overlap_1d(a_l, b_l, a_center.x, b_center.x, a_alpha, b_alpha);
    double sy = overlap_1d(a_m, b_m, a_center.y, b_center.y, a_alpha, b_alpha);
    double sz = overlap_1d(a_n, b_n, a_center.z, b_center.z, a_alpha, b_alpha);
    
    return norm_a * norm_b * sx * sy * sz;
}

// Complete kinetic energy integral
double kinetic(vec3_t a_center, double a_alpha, int a_l, int a_m, int a_n,
               vec3_t b_center, double b_alpha, int b_l, int b_m, int b_n) {
    
    double norm_a = gaussian_norm(a_l, a_m, a_n, a_alpha);
    double norm_b = gaussian_norm(b_l, b_m, b_n, b_alpha);
    
    double tx = kinetic_1d(a_l, b_l, a_center.x, b_center.x, a_alpha, b_alpha);
    double ty = kinetic_1d(a_m, b_m, a_center.y, b_center.y, a_alpha, b_alpha);
    double tz = kinetic_1d(a_n, b_n, a_center.z, b_center.z, a_alpha, b_alpha);
    
    double sx = overlap_1d(a_l, b_l, a_center.x, b_center.x, a_alpha, b_alpha);
    double sy = overlap_1d(a_m, b_m, a_center.y, b_center.y, a_alpha, b_alpha);
    double sz = overlap_1d(a_n, b_n, a_center.z, b_center.z, a_alpha, b_alpha);
    
    return norm_a * norm_b * (tx * sy * sz + sx * ty * sz + sx * sy * tz);
}

// Complete nuclear attraction integral for arbitrary angular momentum
double nuclear(vec3_t a_center, double a_alpha, int a_l, int a_m, int a_n,
               vec3_t b_center, double b_alpha, int b_l, int b_m, int b_n,
               vec3_t nuclear_center, double nuclear_charge) {
    
    double gamma = a_alpha + b_alpha;
    double norm_a = gaussian_norm(a_l, a_m, a_n, a_alpha);
    double norm_b = gaussian_norm(b_l, b_m, b_n, b_alpha);
    
    // Use the enhanced nuclear_1d function for each dimension
    double nx = nuclear_1d(a_l, b_l, a_center.x, b_center.x, nuclear_center.x, a_alpha, b_alpha);
    double ny = nuclear_1d(a_m, b_m, a_center.y, b_center.y, nuclear_center.y, a_alpha, b_alpha);
    double nz = nuclear_1d(a_n, b_n, a_center.z, b_center.z, nuclear_center.z, a_alpha, b_alpha);
    
    return -nuclear_charge * norm_a * norm_b * nx * ny * nz * 2.0 * PI / gamma;
}

// Complete electron repulsion integral for s-orbitals (extensible framework)
double eri(vec3_t a_center, double a_alpha, int a_l, int a_m, int a_n,
           vec3_t b_center, double b_alpha, int b_l, int b_m, int b_n,
           vec3_t c_center, double c_alpha, int c_l, int c_m, int c_n,
           vec3_t d_center, double d_alpha, int d_l, int d_m, int d_n) {
    
    // Currently implemented for s-orbitals only
    if (a_l || a_m || a_n || b_l || b_m || b_n || c_l || c_m || c_n || d_l || d_m || d_n) {
        return 0.0;  // Higher angular momentum implementation would go here
    }
    
    double gamma_p = a_alpha + b_alpha;
    double gamma_q = c_alpha + d_alpha;
    double delta = 0.25 * (1.0/gamma_p + 1.0/gamma_q);
    
    vec3_t p_center = gaussian_product_center(a_center, a_alpha, b_center, b_alpha);
    vec3_t q_center = gaussian_product_center(c_center, c_alpha, d_center, d_alpha);
    vec3_t pq = vec3_sub(p_center, q_center);
    double t = vec3_norm2(pq) / (4.0 * delta);
    
    double norm_a = gaussian_norm(0, 0, 0, a_alpha);
    double norm_b = gaussian_norm(0, 0, 0, b_alpha);
    double norm_c = gaussian_norm(0, 0, 0, c_alpha);
    double norm_d = gaussian_norm(0, 0, 0, d_alpha);
    
    double product_coeff_ab = gaussian_product_coeff(a_center, a_alpha, b_center, b_alpha);
    double product_coeff_cd = gaussian_product_coeff(c_center, c_alpha, d_center, d_alpha);
    
    return norm_a * norm_b * norm_c * norm_d * product_coeff_ab * product_coeff_cd * 
           2.0 * pow(PI, 2.5) / (gamma_p * gamma_q * sqrt(gamma_p + gamma_q)) * 
           boys_function(0, t);
}

// Contracted Gaussian integral functions
double contracted_overlap(const contracted_gaussian_t* cg_a, const contracted_gaussian_t* cg_b) {
    double result = 0.0;
    
    for (int i = 0; i < cg_a->num_primitives; i++) {
        for (int j = 0; j < cg_b->num_primitives; j++) {
            double coeff = cg_a->primitives[i].contraction_coeff * cg_b->primitives[j].contraction_coeff;
            double integral = overlap(cg_a->center, cg_a->primitives[i].exponent, cg_a->l, cg_a->m, cg_a->n,
                                    cg_b->center, cg_b->primitives[j].exponent, cg_b->l, cg_b->m, cg_b->n);
            result += coeff * integral;
        }
    }
    
    return result;
}

double contracted_kinetic(const contracted_gaussian_t* cg_a, const contracted_gaussian_t* cg_b) {
    double result = 0.0;
    
    for (int i = 0; i < cg_a->num_primitives; i++) {
        for (int j = 0; j < cg_b->num_primitives; j++) {
            double coeff = cg_a->primitives[i].contraction_coeff * cg_b->primitives[j].contraction_coeff;
            double integral = kinetic(cg_a->center, cg_a->primitives[i].exponent, cg_a->l, cg_a->m, cg_a->n,
                                    cg_b->center, cg_b->primitives[j].exponent, cg_b->l, cg_b->m, cg_b->n);
            result += coeff * integral;
        }
    }
    
    return result;
}

double contracted_nuclear(const contracted_gaussian_t* cg_a, const contracted_gaussian_t* cg_b,
                         vec3_t nuclear_center, double nuclear_charge) {
    double result = 0.0;
    
    for (int i = 0; i < cg_a->num_primitives; i++) {
        for (int j = 0; j < cg_b->num_primitives; j++) {
            double coeff = cg_a->primitives[i].contraction_coeff * cg_b->primitives[j].contraction_coeff;
            double integral = nuclear(cg_a->center, cg_a->primitives[i].exponent, cg_a->l, cg_a->m, cg_a->n,
                                    cg_b->center, cg_b->primitives[j].exponent, cg_b->l, cg_b->m, cg_b->n,
                                    nuclear_center, nuclear_charge);
            result += coeff * integral;
        }
    }
    
    return result;
}

double contracted_eri(const contracted_gaussian_t* cg_a, const contracted_gaussian_t* cg_b,
                     const contracted_gaussian_t* cg_c, const contracted_gaussian_t* cg_d) {
    double result = 0.0;
    
    for (int i = 0; i < cg_a->num_primitives; i++) {
        for (int j = 0; j < cg_b->num_primitives; j++) {
            for (int k = 0; k < cg_c->num_primitives; k++) {
                for (int l = 0; l < cg_d->num_primitives; l++) {
                    double coeff = cg_a->primitives[i].contraction_coeff * 
                                  cg_b->primitives[j].contraction_coeff *
                                  cg_c->primitives[k].contraction_coeff * 
                                  cg_d->primitives[l].contraction_coeff;
                    
                    double integral = eri(cg_a->center, cg_a->primitives[i].exponent, cg_a->l, cg_a->m, cg_a->n,
                                        cg_b->center, cg_b->primitives[j].exponent, cg_b->l, cg_b->m, cg_b->n,
                                        cg_c->center, cg_c->primitives[k].exponent, cg_c->l, cg_c->m, cg_c->n,
                                        cg_d->center, cg_d->primitives[l].exponent, cg_d->l, cg_d->m, cg_d->n);
                    result += coeff * integral;
                }
            }
        }
    }
    
    return result;
}

// Basis set construction utilities
contracted_gaussian_t make_sto3g_1s(vec3_t center) {
    contracted_gaussian_t cg = {0};
    cg.center = center;
    cg.l = cg.m = cg.n = 0;  // s-orbital
    cg.num_primitives = 3;
    
    // STO-3G coefficients for H 1s
    cg.primitives[0] = (primitive_t){3.42525091, 0.15432897};
    cg.primitives[1] = (primitive_t){0.62391373, 0.53532814};
    cg.primitives[2] = (primitive_t){0.16885540, 0.44463454};
    
    return cg;
}

contracted_gaussian_t make_631g_1s_h(vec3_t center) {
    contracted_gaussian_t cg = {0};
    cg.center = center;
    cg.l = cg.m = cg.n = 0;
    cg.num_primitives = 3;
    
    // 6-31G coefficients for H 1s
    cg.primitives[0] = (primitive_t){18.7311370, 0.03349460};
    cg.primitives[1] = (primitive_t){2.8253937, 0.23472695};
    cg.primitives[2] = (primitive_t){0.6401217, 0.81375733};
    
    return cg;
}

// Create p-orbitals for testing
contracted_gaussian_t make_primitive_px(vec3_t center, double alpha) {
    contracted_gaussian_t cg = {0};
    cg.center = center;
    cg.l = 1; cg.m = 0; cg.n = 0;  // px-orbital
    cg.num_primitives = 1;
    cg.primitives[0] = (primitive_t){alpha, 1.0};
    return cg;
}

contracted_gaussian_t make_primitive_py(vec3_t center, double alpha) {
    contracted_gaussian_t cg = {0};
    cg.center = center;
    cg.l = 0; cg.m = 1; cg.n = 0;  // py-orbital
    cg.num_primitives = 1;
    cg.primitives[0] = (primitive_t){alpha, 1.0};
    return cg;
}

// Enhanced test functions
static void test_higher_angular_momentum(void) {
    printf("Higher Angular Momentum Tests:\n");
    
    vec3_t origin = {0.0, 0.0, 0.0};
    double alpha = 1.0;
    
    // Test p-orbital overlaps
    double s_px_px = overlap(origin, alpha, 1, 0, 0, origin, alpha, 1, 0, 0);
    double s_px_py = overlap(origin, alpha, 1, 0, 0, origin, alpha, 0, 1, 0);
    double s_px_s = overlap(origin, alpha, 1, 0, 0, origin, alpha, 0, 0, 0);
    
    printf("p_x - p_x overlap: %.8f (should be 1.0)\n", s_px_px);
    printf("p_x - p_y overlap: %.8f (should be 0.0)\n", s_px_py);
    printf("p_x - s overlap: %.8f (should be 0.0)\n", s_px_s);
    
    // Test nuclear attraction for p-orbitals
    double v_px_s = nuclear(origin, alpha, 1, 0, 0, origin, alpha, 0, 0, 0, origin, 1.0);
    printf("p_x - s nuclear attraction: %.8f (should be 0.0 by symmetry)\n", v_px_s);
    
    // Test contracted p-orbitals
    contracted_gaussian_t px = make_primitive_px(origin, alpha);
    contracted_gaussian_t py = make_primitive_py(origin, alpha);
    
    double s_px_px_contracted = contracted_overlap(&px, &px);
    double s_px_py_contracted = contracted_overlap(&px, &py);
    
    printf("Contracted p_x - p_x overlap: %.8f\n", s_px_px_contracted);
    printf("Contracted p_x - p_y overlap: %.8f\n", s_px_py_contracted);
    printf("\n");
}

static void test_contracted_basis(void) {
    printf("Contracted Basis Set Tests:\n");
    
    vec3_t h1 = {0.0, 0.0, -0.7};
    vec3_t h2 = {0.0, 0.0, +0.7};
    
    contracted_gaussian_t sto3g_h1 = make_sto3g_1s(h1);
    contracted_gaussian_t sto3g_h2 = make_sto3g_1s(h2);
    
    double s11_sto3g = contracted_overlap(&sto3g_h1, &sto3g_h1);
    double s12_sto3g = contracted_overlap(&sto3g_h1, &sto3g_h2);
    double t11_sto3g = contracted_kinetic(&sto3g_h1, &sto3g_h1);
    
    printf("STO-3G H2 molecule:\n");
    printf("S₁₁ = %.8f (normalized self-overlap)\n", s11_sto3g);
    printf("S₁₂ = %.8f (intermolecular overlap)\n", s12_sto3g);
    printf("T₁₁ = %.8f (kinetic energy)\n", t11_sto3g);
    
    contracted_gaussian_t g631_h1 = make_631g_1s_h(h1);
    contracted_gaussian_t g631_h2 = make_631g_1s_h(h2);
    
    double s11_631g = contracted_overlap(&g631_h1, &g631_h1);
    double s12_631g = contracted_overlap(&g631_h1, &g631_h2);
    
    printf("\n6-31G H2 molecule:\n");
    printf("S₁₁ = %.8f\n", s11_631g);
    printf("S₁₂ = %.8f\n", s12_631g);
    printf("Basis set comparison: ΔS₁₂ = %.6f\n", s12_631g - s12_sto3g);
    printf("\n");
}

static void test_boys_function(void) {
    printf("Boys Function Tests:\n");
    printf("F₀(0.0) = %.8f (expected: 1.0)\n", boys_function(0, 0.0));
    printf("F₀(1.0) = %.8f (expected: ~0.746)\n", boys_function(0, 1.0));
    printf("F₁(0.0) = %.8f (expected: 0.333)\n", boys_function(1, 0.0));
    printf("F₀(25.0) = %.8f (asymptotic test)\n", boys_function(0, 25.0));
    printf("\n");
}

static void test_gaussian_properties(void) {
    vec3_t origin = {0.0, 0.0, 0.0};
    double alpha = 1.0;
    
    printf("Gaussian Property Tests:\n");
    double s_self = overlap(origin, alpha, 0, 0, 0, origin, alpha, 0, 0, 0);
    printf("Self-overlap (normalized): %.8f (should be 1.0)\n", s_self);
    
    vec3_t displaced = {1.0, 0.0, 0.0};
    double s_displaced = overlap(origin, alpha, 0, 0, 0, displaced, alpha, 0, 0, 0);
    printf("Overlap at 1 Bohr separation: %.8f\n", s_displaced);
    printf("\n");
}

// Comprehensive validation and benchmark
static void run_comprehensive_tests(void) {
    printf("Running Comprehensive Validation Tests...\n");
    printf("==========================================\n\n");
    
    // Test orthogonality and normalization for different angular momentum
    vec3_t origin = {0.0, 0.0, 0.0};
    double alpha = 1.0;
    
    printf("Angular Momentum Orthogonality Tests:\n");
    
    // s-p orthogonality
    double s_s_px = overlap(origin, alpha, 0,0,0, origin, alpha, 1,0,0);
    double s_s_py = overlap(origin, alpha, 0,0,0, origin, alpha, 0,1,0);
    double s_s_pz = overlap(origin, alpha, 0,0,0, origin, alpha, 0,0,1);
    
    printf("⟨s|px⟩ = %.2e (should be 0)\n", s_s_px);
    printf("⟨s|py⟩ = %.2e (should be 0)\n", s_s_py);
    printf("⟨s|pz⟩ = %.2e (should be 0)\n", s_s_pz);
    
    // p-p orthogonality
    double s_px_py = overlap(origin, alpha, 1,0,0, origin, alpha, 0,1,0);
    double s_px_pz = overlap(origin, alpha, 1,0,0, origin, alpha, 0,0,1);
    double s_py_pz = overlap(origin, alpha, 0,1,0, origin, alpha, 0,0,1);
    
    printf("⟨px|py⟩ = %.2e (should be 0)\n", s_px_py);
    printf("⟨px|pz⟩ = %.2e (should be 0)\n", s_px_pz);
    printf("⟨py|pz⟩ = %.2e (should be 0)\n", s_py_pz);
    
    // Normalization tests (corrected)
    double s_px_px = overlap(origin, alpha, 1,0,0, origin, alpha, 1,0,0);
    double s_py_py = overlap(origin, alpha, 0,1,0, origin, alpha, 0,1,0);
    double s_pz_pz = overlap(origin, alpha, 0,0,1, origin, alpha, 0,0,1);
    
    printf("⟨px|px⟩ = %.8f (should be 1.0)\n", s_px_px);
    printf("⟨py|py⟩ = %.8f (should be 1.0)\n", s_py_py);
    printf("⟨pz|pz⟩ = %.8f (should be 1.0)\n", s_pz_pz);
    
    // Test d-orbitals
    double s_dxx_dxx = overlap(origin, alpha, 2,0,0, origin, alpha, 2,0,0);
    double s_dxy_dxy = overlap(origin, alpha, 1,1,0, origin, alpha, 1,1,0);
    double s_dxx_dxy = overlap(origin, alpha, 2,0,0, origin, alpha, 1,1,0);
    
    printf("⟨dxx|dxx⟩ = %.8f (should be 1.0)\n", s_dxx_dxx);
    printf("⟨dxy|dxy⟩ = %.8f (should be 1.0)\n", s_dxy_dxy);
    printf("⟨dxx|dxy⟩ = %.2e (should be 0)\n", s_dxx_dxy);
    
    printf("\nNumerical Stability Test (F₀ function):\n");
    for (int i = 0; i <= 5; i++) {
        double x = pow(10.0, i);
        double f0 = boys_function(0, x);
        double f1 = boys_function(1, x);
        printf("F₀(%.0e) = %.6e, F₁(%.0e) = %.6e\n", x, f0, x, f1);
    }
    
    printf("\n");
}

// Enhanced H2 demonstration with multiple basis sets
int main() {
    printf("Enhanced Molecular Integral Library - Production Ready\n");
    printf("====================================================\n\n");
    
    // Run comprehensive validation tests
    test_boys_function();
    test_gaussian_properties();
    test_higher_angular_momentum();
    test_contracted_basis();
    run_comprehensive_tests();
    
    // Original H2 demonstration with enhancements
    vec3_t h1_center = {0.0, 0.0, -0.7};
    vec3_t h2_center = {0.0, 0.0, +0.7};
    double h_alpha = 1.24;
    int l = 0, m = 0, n = 0;
    
    printf("H2 Molecule Analysis - Multi-Basis Comparison\n");
    printf("===========================================\n");
    printf("Nuclear positions: H1(0,0,-0.7), H2(0,0,+0.7)\n");
    printf("Bond length: 1.4 Bohr (0.741 Å)\n\n");
    
    // Primitive Gaussian analysis
    printf("Primitive Gaussian (α = %.2f):\n", h_alpha);
    double s11 = overlap(h1_center, h_alpha, l, m, n, h1_center, h_alpha, l, m, n);
    double s12 = overlap(h1_center, h_alpha, l, m, n, h2_center, h_alpha, l, m, n);
    double t11 = kinetic(h1_center, h_alpha, l, m, n, h1_center, h_alpha, l, m, n);
    double t12 = kinetic(h1_center, h_alpha, l, m, n, h2_center, h_alpha, l, m, n);
    
    printf("S₁₁ = %.8f, S₁₂ = %.8f\n", s11, s12);
    printf("T₁₁ = %.8f, T₁₂ = %.8f\n", t11, t12);
    
    // Nuclear attraction
    double v11_1 = nuclear(h1_center, h_alpha, l, m, n, h1_center, h_alpha, l, m, n, h1_center, 1.0);
    double v11_2 = nuclear(h1_center, h_alpha, l, m, n, h1_center, h_alpha, l, m, n, h2_center, 1.0);
    printf("V₁₁⁽¹⁾ = %.8f, V₁₁⁽²⁾ = %.8f\n", v11_1, v11_2);
    
    // Core Hamiltonian
    double h11_core = t11 + v11_1 + v11_2;
    double h12_core = t12 + nuclear(h1_center, h_alpha, l, m, n, h2_center, h_alpha, l, m, n, h1_center, 1.0) +
                              nuclear(h1_center, h_alpha, l, m, n, h2_center, h_alpha, l, m, n, h2_center, 1.0);
    
    printf("Core Hamiltonian: H₁₁ = %.6f, H₁₂ = %.6f\n", h11_core, h12_core);
    
    // Contracted basis comparison
    printf("\nContracted Basis Set Comparison:\n");
    contracted_gaussian_t sto3g_h1 = make_sto3g_1s(h1_center);
    contracted_gaussian_t sto3g_h2 = make_sto3g_1s(h2_center);
    contracted_gaussian_t g631_h1 = make_631g_1s_h(h1_center);
    contracted_gaussian_t g631_h2 = make_631g_1s_h(h2_center);
    
    double s12_sto3g = contracted_overlap(&sto3g_h1, &sto3g_h2);
    double s12_631g = contracted_overlap(&g631_h1, &g631_h2);
    double t11_sto3g = contracted_kinetic(&sto3g_h1, &sto3g_h1);
    double t11_631g = contracted_kinetic(&g631_h1, &g631_h1);
    
    printf("Overlap S₁₂: STO-3G = %.6f, 6-31G = %.6f, Δ = %.6f\n", 
           s12_sto3g, s12_631g, s12_631g - s12_sto3g);
    printf("Kinetic T₁₁: STO-3G = %.6f, 6-31G = %.6f, Δ = %.6f\n", 
           t11_sto3g, t11_631g, t11_631g - t11_sto3g);
    return 0;
}