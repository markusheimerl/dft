#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// Minimal structure for 3D vectors only
typedef struct { double x, y, z; } vec3_t;

// Mathematical constants
#define PI 3.14159265358979323846
#define SQRT_PI 1.77245385090551602729

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

// Pure functional factorial implementation
static double factorial_impl(int n, double acc) {
    return (n <= 1) ? acc : factorial_impl(n - 1, acc * n);
}

static double factorial(int n) {
    return (n < 0) ? 0.0 : factorial_impl(n, 1.0);
}

// Double factorial: n!! = n * (n-2) * (n-4) * ...
static double double_factorial_impl(int n, double acc) {
    return (n <= 1) ? acc : double_factorial_impl(n - 2, acc * n);
}

static double double_factorial(int n) {
    return (n < 0) ? 1.0 : double_factorial_impl(n, 1.0);
}

// Boys function F_n(x) using series expansion for small x, asymptotic for large x
static double boys_series(int n, double x, int term, double acc) {
    if (term > 50) return acc;  // Convergence check
    double term_val = pow(-x, term) / (factorial(term) * (2.0 * n + 2.0 * term + 1.0));
    if (fabs(term_val) < 1e-15) return acc;
    return boys_series(n, x, term + 1, acc + term_val);
}

static double boys_asymptotic(int n, double x) {
    return double_factorial(2 * n - 1) * sqrt(PI) / (2.0 * pow(2.0 * x, n + 0.5));
}

static double boys_function(int n, double x) {
    if (x < 1e-10) return 1.0 / (2.0 * n + 1.0);  // Limit as x -> 0
    if (x > 20.0) return boys_asymptotic(n, x);
    return boys_series(n, x, 0, 0.0);
}

// Gaussian normalization constant
static double gaussian_norm(int l, int m, int n, double alpha) {
    double norm = pow(2.0 * alpha / PI, 0.75) * 
                  pow(8.0 * alpha, (l + m + n) / 2.0) / 
                  sqrt(double_factorial(2 * l - 1) * 
                       double_factorial(2 * m - 1) * 
                       double_factorial(2 * n - 1));
    return norm;
}

// Gaussian product rule: product of two Gaussians is a Gaussian
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

// 1D overlap integral - analytical formula for Gaussian primitives
static double overlap_1d(int i, int j, double xa, double xb, double alpha, double beta) {
    if (i < 0 || j < 0) return 0.0;
    
    double gamma = alpha + beta;
    double xp = (alpha * xa + beta * xb) / gamma;
    double prefactor = exp(-alpha * beta * (xa - xb) * (xa - xb) / gamma) * sqrt(PI / gamma);
    
    // For s-orbitals (i=j=0)
    if (i == 0 && j == 0) {
        return prefactor;
    }
    
    // Recursive Obara-Saika for higher angular momentum
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

// 1D kinetic energy integral
static double kinetic_1d(int i, int j, double xa, double xb, double alpha, double beta) {
    // T_ij = β * [(2j+1)S_ij - 2β*S_i,j+2 - (1/2)*j*(j-1)*S_i,j-2]
    double s_ij = overlap_1d(i, j, xa, xb, alpha, beta);
    double s_ij_plus2 = overlap_1d(i, j + 2, xa, xb, alpha, beta);
    double s_ij_minus2 = (j >= 2) ? overlap_1d(i, j - 2, xa, xb, alpha, beta) : 0.0;
    
    double term1 = beta * (2 * j + 1) * s_ij;
    double term2 = -2.0 * beta * beta * s_ij_plus2;
    double term3 = (j >= 2) ? -0.5 * beta * j * (j - 1) * s_ij_minus2 : 0.0;
    
    return term1 + term2 + term3;
}

// Main overlap integral
double overlap(vec3_t a_center, double a_alpha, int a_l, int a_m, int a_n,
               vec3_t b_center, double b_alpha, int b_l, int b_m, int b_n) {
    
    double norm_a = gaussian_norm(a_l, a_m, a_n, a_alpha);
    double norm_b = gaussian_norm(b_l, b_m, b_n, b_alpha);
    
    double sx = overlap_1d(a_l, b_l, a_center.x, b_center.x, a_alpha, b_alpha);
    double sy = overlap_1d(a_m, b_m, a_center.y, b_center.y, a_alpha, b_alpha);
    double sz = overlap_1d(a_n, b_n, a_center.z, b_center.z, a_alpha, b_alpha);
    
    return norm_a * norm_b * sx * sy * sz;
}

// Kinetic energy integral
double kinetic(vec3_t a_center, double a_alpha, int a_l, int a_m, int a_n,
               vec3_t b_center, double b_alpha, int b_l, int b_m, int b_n) {
    
    double norm_a = gaussian_norm(a_l, a_m, a_n, a_alpha);
    double norm_b = gaussian_norm(b_l, b_m, b_n, b_alpha);
    
    // T = Tx * Sy * Sz + Sx * Ty * Sz + Sx * Sy * Tz
    double tx = kinetic_1d(a_l, b_l, a_center.x, b_center.x, a_alpha, b_alpha);
    double ty = kinetic_1d(a_m, b_m, a_center.y, b_center.y, a_alpha, b_alpha);
    double tz = kinetic_1d(a_n, b_n, a_center.z, b_center.z, a_alpha, b_alpha);
    
    double sx = overlap_1d(a_l, b_l, a_center.x, b_center.x, a_alpha, b_alpha);
    double sy = overlap_1d(a_m, b_m, a_center.y, b_center.y, a_alpha, b_alpha);
    double sz = overlap_1d(a_n, b_n, a_center.z, b_center.z, a_alpha, b_alpha);
    
    return norm_a * norm_b * (tx * sy * sz + sx * ty * sz + sx * sy * tz);
}

// Nuclear attraction integral
double nuclear(vec3_t a_center, double a_alpha, int a_l, int a_m, int a_n,
               vec3_t b_center, double b_alpha, int b_l, int b_m, int b_n,
               vec3_t nuclear_center, double nuclear_charge) {
    
    double gamma = a_alpha + b_alpha;
    vec3_t p_center = gaussian_product_center(a_center, a_alpha, b_center, b_alpha);
    vec3_t pc = vec3_sub(p_center, nuclear_center);
    double t = gamma * vec3_norm2(pc);
    
    double norm_a = gaussian_norm(a_l, a_m, a_n, a_alpha);
    double norm_b = gaussian_norm(b_l, b_m, b_n, b_alpha);
    double product_coeff = gaussian_product_coeff(a_center, a_alpha, b_center, b_alpha);
    
    // For s-orbitals only
    if (a_l == 0 && a_m == 0 && a_n == 0 && b_l == 0 && b_m == 0 && b_n == 0) {
        return -nuclear_charge * norm_a * norm_b * product_coeff * 
               2.0 * PI / gamma * boys_function(0, t);
    }
    
    return 0.0;  // Higher angular momentum not implemented
}

// Electron repulsion integral (s-orbitals only)
double eri(vec3_t a_center, double a_alpha, int a_l, int a_m, int a_n,
           vec3_t b_center, double b_alpha, int b_l, int b_m, int b_n,
           vec3_t c_center, double c_alpha, int c_l, int c_m, int c_n,
           vec3_t d_center, double d_alpha, int d_l, int d_m, int d_n) {
    
    // Only implement for s-orbitals
    if (a_l || a_m || a_n || b_l || b_m || b_n || c_l || c_m || c_n || d_l || d_m || d_n) {
        return 0.0;
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

// Test functions for validation
static void test_boys_function(void) {
    printf("Boys Function Tests:\n");
    printf("F₀(0.0) = %.8f (expected: 1.0)\n", boys_function(0, 0.0));
    printf("F₀(1.0) = %.8f (expected: ~0.746)\n", boys_function(0, 1.0));
    printf("F₁(0.0) = %.8f (expected: 0.333)\n", boys_function(1, 0.0));
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

// Demonstration: H2 molecule calculation with enhanced validation
int main() {
    printf("Molecular Integral Library - Functional Programming Style\n");
    printf("=========================================================\n\n");
    
    // Run validation tests
    test_boys_function();
    test_gaussian_properties();
    
    // H2 molecule: two hydrogen atoms
    vec3_t h1_center = {0.0, 0.0, -0.7};  // Bohrs
    vec3_t h2_center = {0.0, 0.0, +0.7};
    double h_alpha = 1.24;  // Exponent for H 1s orbital
    
    // All orbitals are 1s (l=m=n=0)
    int l = 0, m = 0, n = 0;
    
    printf("H2 Molecule Integrals (Bond length = 1.4 Bohr)\n");
    printf("Nuclear positions: H1(0,0,-0.7), H2(0,0,+0.7)\n");
    printf("Orbital exponent: α = %.2f\n\n", h_alpha);
    
    // Overlap integrals
    double s11 = overlap(h1_center, h_alpha, l, m, n, h1_center, h_alpha, l, m, n);
    double s12 = overlap(h1_center, h_alpha, l, m, n, h2_center, h_alpha, l, m, n);
    double s22 = overlap(h2_center, h_alpha, l, m, n, h2_center, h_alpha, l, m, n);
    
    printf("Overlap Integrals:\n");
    printf("S₁₁ = %12.8f  (should be ≈ 1.0)\n", s11);
    printf("S₁₂ = %12.8f  (intermolecular overlap)\n", s12);
    printf("S₂₂ = %12.8f  (should be ≈ 1.0)\n", s22);
    printf("Normalization check: |S₁₁ - 1.0| = %.2e\n", fabs(s11 - 1.0));
    printf("Symmetry check: |S₁₁ - S₂₂| = %.2e\n\n", fabs(s11 - s22));
    
    // Kinetic energy integrals
    double t11 = kinetic(h1_center, h_alpha, l, m, n, h1_center, h_alpha, l, m, n);
    double t12 = kinetic(h1_center, h_alpha, l, m, n, h2_center, h_alpha, l, m, n);
    double t22 = kinetic(h2_center, h_alpha, l, m, n, h2_center, h_alpha, l, m, n);
    
    printf("Kinetic Energy Integrals:\n");
    printf("T₁₁ = %12.8f  (α²/2 = %.8f)\n", t11, h_alpha * h_alpha / 2.0);
    printf("T₁₂ = %12.8f\n", t12);
    printf("T₂₂ = %12.8f\n", t22);
    printf("Symmetry check: |T₁₁ - T₂₂| = %.2e\n\n", fabs(t11 - t22));
    
    // Nuclear attraction integrals
    double v11_1 = nuclear(h1_center, h_alpha, l, m, n, h1_center, h_alpha, l, m, n, h1_center, 1.0);
    double v11_2 = nuclear(h1_center, h_alpha, l, m, n, h1_center, h_alpha, l, m, n, h2_center, 1.0);
    double v12_1 = nuclear(h1_center, h_alpha, l, m, n, h2_center, h_alpha, l, m, n, h1_center, 1.0);
    double v12_2 = nuclear(h1_center, h_alpha, l, m, n, h2_center, h_alpha, l, m, n, h2_center, 1.0);
    double v22_1 = nuclear(h2_center, h_alpha, l, m, n, h2_center, h_alpha, l, m, n, h1_center, 1.0);
    double v22_2 = nuclear(h2_center, h_alpha, l, m, n, h2_center, h_alpha, l, m, n, h2_center, 1.0);
    
    printf("Nuclear Attraction Integrals:\n");
    printf("V₁₁⁽¹⁾ = %12.8f  (orbital 1 with nucleus 1)\n", v11_1);
    printf("V₁₁⁽²⁾ = %12.8f  (orbital 1 with nucleus 2)\n", v11_2);
    printf("V₁₂⁽¹⁾ = %12.8f  (cross term with nucleus 1)\n", v12_1);
    printf("V₁₂⁽²⁾ = %12.8f  (cross term with nucleus 2)\n", v12_2);
    printf("Symmetry checks:\n");
    printf("  |V₁₁⁽²⁾ - V₂₂⁽¹⁾| = %.2e\n", fabs(v11_2 - v22_1));
    printf("  |V₁₂⁽¹⁾ - V₁₂⁽²⁾| = %.2e\n\n", fabs(v12_1 - v12_2));
    
    // Electron repulsion integrals
    double eri_1111 = eri(h1_center, h_alpha, l, m, n, h1_center, h_alpha, l, m, n,
                         h1_center, h_alpha, l, m, n, h1_center, h_alpha, l, m, n);
    double eri_1112 = eri(h1_center, h_alpha, l, m, n, h1_center, h_alpha, l, m, n,
                         h1_center, h_alpha, l, m, n, h2_center, h_alpha, l, m, n);
    double eri_1122 = eri(h1_center, h_alpha, l, m, n, h1_center, h_alpha, l, m, n,
                         h2_center, h_alpha, l, m, n, h2_center, h_alpha, l, m, n);
    double eri_1212 = eri(h1_center, h_alpha, l, m, n, h2_center, h_alpha, l, m, n,
                         h1_center, h_alpha, l, m, n, h2_center, h_alpha, l, m, n);
    
    printf("Electron Repulsion Integrals:\n");
    printf("(11|11) = %12.8f  (Coulomb self-repulsion)\n", eri_1111);
    printf("(11|12) = %12.8f  (Cross Coulomb)\n", eri_1112);
    printf("(11|22) = %12.8f  (Inter-center Coulomb)\n", eri_1122);
    printf("(12|12) = %12.8f  (Exchange integral)\n\n", eri_1212);
    
    // Core Hamiltonian matrix elements
    double h11 = t11 + v11_1 + v11_2;
    double h12 = t12 + v12_1 + v12_2;
    double h22 = t22 + v22_1 + v22_2;
    
    printf("Core Hamiltonian Matrix:\n");
    printf("H₁₁ = %12.8f\n", h11);
    printf("H₁₂ = %12.8f\n", h12);
    printf("H₂₂ = %12.8f\n", h22);
    printf("Symmetry check: |H₁₁ - H₂₂| = %.2e\n\n", fabs(h11 - h22));
    
    // Simple energy estimates
    printf("Simple Energy Estimates:\n");
    printf("Binding energy (crude): %.6f Hartree\n", h12 / s12);
    printf("Orbital energy estimate: %.6f Hartree\n\n", h11);
    return 0;
}