#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAX_ITER 100
#define TOLERANCE 1e-8
#define PI 3.14159265359

// Structure to hold atomic orbital parameters
typedef struct {
    double alpha;    // exponent
    double x, y, z;  // position
} orbital_t;

// Function prototypes
double overlap_integral(orbital_t *a, orbital_t *b);
double kinetic_integral(orbital_t *a, orbital_t *b);
double nuclear_attraction(orbital_t *a, orbital_t *b, double Zx, double Zy, double Zz);
double electron_repulsion(orbital_t *a, orbital_t *b, orbital_t *c, orbital_t *d);
double boys_function(double t, int n);
void normalize_orbitals(orbital_t *orbs, int n);
int hartree_fock(orbital_t *orbitals, double *nuclear_pos, int num_orbs);

// Overlap integral between two 1s Gaussian orbitals
double overlap_integral(orbital_t *a, orbital_t *b) {
    double dx = a->x - b->x;
    double dy = a->y - b->y;
    double dz = a->z - b->z;
    double r_sq = dx*dx + dy*dy + dz*dz;
    
    double p = a->alpha + b->alpha;
    double prefactor = pow(PI/p, 1.5);
    
    return prefactor * exp(-a->alpha * b->alpha * r_sq / p);
}

// Kinetic energy integral
double kinetic_integral(orbital_t *a, orbital_t *b) {
    double dx = a->x - b->x;
    double dy = a->y - b->y;
    double dz = a->z - b->z;
    double r_sq = dx*dx + dy*dy + dz*dz;
    
    double p = a->alpha + b->alpha;
    double mu = a->alpha * b->alpha / p;
    double prefactor = pow(PI/p, 1.5);
    
    return mu * (3.0 - 2.0 * mu * r_sq) * prefactor * exp(-mu * r_sq);
}

// Nuclear attraction integral
double nuclear_attraction(orbital_t *a, orbital_t *b, double Zx, double Zy, double Zz) {
    double dx = a->x - b->x;
    double dy = a->y - b->y;
    double dz = a->z - b->z;
    double r_sq = dx*dx + dy*dy + dz*dz;
    
    double p = a->alpha + b->alpha;
    double Px = (a->alpha * a->x + b->alpha * b->x) / p;
    double Py = (a->alpha * a->y + b->alpha * b->y) / p;
    double Pz = (a->alpha * a->z + b->alpha * b->z) / p;
    
    double PCx = Px - Zx;
    double PCy = Py - Zy;
    double PCz = Pz - Zz;
    double PC_sq = PCx*PCx + PCy*PCy + PCz*PCz;
    
    double prefactor = -2.0 * PI / p;
    double overlap = pow(PI/p, 1.5) * exp(-a->alpha * b->alpha * r_sq / p);
    
    if (PC_sq < 1e-10) {
        return prefactor * overlap;
    }
    
    double t = p * PC_sq;
    double boys0 = sqrt(PI / (4.0 * t)) * erf(sqrt(t));
    
    return prefactor * overlap * boys0;
}

// Simplified Boys function F_0(t)
double boys_function(double t, int n) {
    if (n != 0) return 0.0;  // Only implement F_0 for simplicity
    
    if (t < 1e-10) return 1.0;
    
    return sqrt(PI / (4.0 * t)) * erf(sqrt(t));
}

// Electron repulsion integral (simplified for 1s orbitals)
double electron_repulsion(orbital_t *a, orbital_t *b, orbital_t *c, orbital_t *d) {
    double dx1 = a->x - b->x, dy1 = a->y - b->y, dz1 = a->z - b->z;
    double dx2 = c->x - d->x, dy2 = c->y - d->y, dz2 = c->z - d->z;
    double r1_sq = dx1*dx1 + dy1*dy1 + dz1*dz1;
    double r2_sq = dx2*dx2 + dy2*dy2 + dz2*dz2;
    
    double p1 = a->alpha + b->alpha;
    double p2 = c->alpha + d->alpha;
    double p_total = p1 + p2;
    
    double P1x = (a->alpha * a->x + b->alpha * b->x) / p1;
    double P1y = (a->alpha * a->y + b->alpha * b->y) / p1;
    double P1z = (a->alpha * a->z + b->alpha * b->z) / p1;
    
    double P2x = (c->alpha * c->x + d->alpha * d->x) / p2;
    double P2y = (c->alpha * c->y + d->alpha * d->y) / p2;
    double P2z = (c->alpha * c->z + d->alpha * d->z) / p2;
    
    double P12x = P1x - P2x, P12y = P1y - P2y, P12z = P1z - P2z;
    double P12_sq = P12x*P12x + P12y*P12y + P12z*P12z;
    
    double prefactor = 2.0 * pow(PI, 2.5) / (p1 * p2 * sqrt(p_total));
    double exp_factor = exp(-a->alpha * b->alpha * r1_sq / p1 - 
                           c->alpha * d->alpha * r2_sq / p2);
    
    double t = p1 * p2 * P12_sq / p_total;
    double boys0 = boys_function(t, 0);
    
    return prefactor * exp_factor * boys0;
}

// Solve 2x2 generalized eigenvalue problem FC = SCE
void solve_2x2_generalized_eigenproblem(double F[2][2], double S[2][2], 
                                       double eigenvals[2], double C[2][2]) {
    // Solve |F - λS| = 0
    double det_S = S[0][0] * S[1][1] - S[0][1] * S[1][0];
    
    if (fabs(det_S) < 1e-12) {
        printf("Error: Singular overlap matrix!\n");
        return;
    }
    
    // Compute S^(-1) * F
    double S_inv[2][2];
    S_inv[0][0] = S[1][1] / det_S;
    S_inv[0][1] = -S[0][1] / det_S;
    S_inv[1][0] = -S[1][0] / det_S;
    S_inv[1][1] = S[0][0] / det_S;
    
    double SF[2][2];
    SF[0][0] = S_inv[0][0] * F[0][0] + S_inv[0][1] * F[1][0];
    SF[0][1] = S_inv[0][0] * F[0][1] + S_inv[0][1] * F[1][1];
    SF[1][0] = S_inv[1][0] * F[0][0] + S_inv[1][1] * F[1][0];
    SF[1][1] = S_inv[1][0] * F[0][1] + S_inv[1][1] * F[1][1];
    
    // Solve characteristic equation |SF - λI| = 0
    double trace = SF[0][0] + SF[1][1];
    double det = SF[0][0] * SF[1][1] - SF[0][1] * SF[1][0];
    double discriminant = trace * trace - 4.0 * det;
    
    if (discriminant < 0) {
        printf("Error: Complex eigenvalues!\n");
        return;
    }
    
    eigenvals[0] = 0.5 * (trace - sqrt(discriminant));  // Lower eigenvalue
    eigenvals[1] = 0.5 * (trace + sqrt(discriminant));  // Higher eigenvalue
    
    // Find eigenvector for lowest eigenvalue
    double lambda = eigenvals[0];
    double a = SF[0][0] - lambda;
    double b = SF[0][1];
    
    if (fabs(b) > 1e-12) {
        C[0][0] = 1.0;
        C[1][0] = -a / b;
    } else {
        C[0][0] = 1.0;
        C[1][0] = 0.0;
    }
    
    // Normalize eigenvector with respect to S
    double norm_sq = 0.0;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            norm_sq += C[i][0] * S[i][j] * C[j][0];
        }
    }
    double norm = 1.0 / sqrt(norm_sq);
    C[0][0] *= norm;
    C[1][0] *= norm;
    
    // Second eigenvector (not needed for ground state, but for completeness)
    C[0][1] = C[1][0];
    C[1][1] = -C[0][0];
}

// Main Hartree-Fock routine
int hartree_fock(orbital_t *orbitals, double *nuclear_pos, int num_orbs) {
    printf("Starting Hartree-Fock calculation for H2\n");
    printf("Number of basis functions: %d\n\n", num_orbs);
    
    // Matrices
    double S[2][2], H[2][2], F[2][2], P[2][2], P_old[2][2];
    double C[2][2], eigenvals[2];
    
    // Initialize density matrix to zero
    memset(P, 0, sizeof(P));
    
    // Calculate overlap matrix S
    printf("Calculating overlap matrix:\n");
    for (int i = 0; i < num_orbs; i++) {
        for (int j = 0; j < num_orbs; j++) {
            S[i][j] = overlap_integral(&orbitals[i], &orbitals[j]);
            printf("S[%d][%d] = %12.8f\n", i, j, S[i][j]);
        }
    }
    printf("\n");
    
    // Calculate core Hamiltonian H = T + V_ne
    printf("Calculating core Hamiltonian:\n");
    for (int i = 0; i < num_orbs; i++) {
        for (int j = 0; j < num_orbs; j++) {
            double T = kinetic_integral(&orbitals[i], &orbitals[j]);
            double V1 = nuclear_attraction(&orbitals[i], &orbitals[j], 
                                         nuclear_pos[0], nuclear_pos[1], nuclear_pos[2]);
            double V2 = nuclear_attraction(&orbitals[i], &orbitals[j], 
                                         nuclear_pos[3], nuclear_pos[4], nuclear_pos[5]);
            H[i][j] = T + V1 + V2;
            printf("H[%d][%d] = %12.8f (T=%8.4f, V1=%8.4f, V2=%8.4f)\n", 
                   i, j, H[i][j], T, V1, V2);
        }
    }
    printf("\n");
    
    // SCF iterations
    double energy_old = 0.0;
    for (int iter = 0; iter < MAX_ITER; iter++) {
        // Save old density matrix
        memcpy(P_old, P, sizeof(P));
        
        // Build Fock matrix F = H + G
        for (int i = 0; i < num_orbs; i++) {
            for (int j = 0; j < num_orbs; j++) {
                F[i][j] = H[i][j];
                
                // Add electron repulsion terms
                for (int k = 0; k < num_orbs; k++) {
                    for (int l = 0; l < num_orbs; l++) {
                        double coulomb = electron_repulsion(&orbitals[i], &orbitals[j], 
                                                          &orbitals[k], &orbitals[l]);
                        double exchange = electron_repulsion(&orbitals[i], &orbitals[k], 
                                                           &orbitals[j], &orbitals[l]);
                        F[i][j] += P[k][l] * (2.0 * coulomb - exchange);
                    }
                }
            }
        }
        
        // Solve generalized eigenvalue problem FC = SCE
        solve_2x2_generalized_eigenproblem(F, S, eigenvals, C);
        
        // Build new density matrix (occupy lowest orbital with 2 electrons)
        for (int i = 0; i < num_orbs; i++) {
            for (int j = 0; j < num_orbs; j++) {
                P[i][j] = 2.0 * C[i][0] * C[j][0];  // Factor of 2 for two electrons
            }
        }
        
        // Calculate energy
        double energy = 0.0;
        for (int i = 0; i < num_orbs; i++) {
            for (int j = 0; j < num_orbs; j++) {
                energy += 0.5 * P[i][j] * (H[i][j] + F[i][j]);
            }
        }
        
        // Add nuclear repulsion
        double R = sqrt(pow(nuclear_pos[3] - nuclear_pos[0], 2) +
                       pow(nuclear_pos[4] - nuclear_pos[1], 2) +
                       pow(nuclear_pos[5] - nuclear_pos[2], 2));
        energy += 1.0 / R;  // Nuclear repulsion (Z1*Z2/R with Z1=Z2=1)
        
        printf("Iteration %3d: Energy = %15.10f Hartree (ε₀ = %12.8f)\n", 
               iter+1, energy, eigenvals[0]);
        
        // Check convergence
        if (fabs(energy - energy_old) < TOLERANCE) {
            printf("\nSCF Converged!\n");
            printf("Final Energy: %15.10f Hartree\n", energy);
            printf("Final Energy: %15.10f eV\n", energy * 27.2114);
            printf("Orbital energies: ε₀ = %12.8f, ε₁ = %12.8f Hartree\n", 
                   eigenvals[0], eigenvals[1]);
            printf("Molecular orbital coefficients:\n");
            printf("MO 1: %8.5f * φ₁ + %8.5f * φ₂\n", C[0][0], C[1][0]);
            printf("MO 2: %8.5f * φ₁ + %8.5f * φ₂\n", C[0][1], C[1][1]);
            return iter + 1;
        }
        
        energy_old = energy;
    }
    
    printf("SCF did not converge in %d iterations\n", MAX_ITER);
    return -1;
}

int main() {
    printf("Two-Hydrogen Hartree-Fock Simulation\n");
    printf("=====================================\n\n");
    
    // Set up two hydrogen atoms
    orbital_t orbitals[2];
    
    // H atom 1 at origin
    orbitals[0].alpha = 1.24;  // STO-1G exponent for hydrogen
    orbitals[0].x = 0.0;
    orbitals[0].y = 0.0;
    orbitals[0].z = 0.0;
    
    // H atom 2 at distance R along z-axis
    double bond_length = 1.4;  // Bohr radii
    orbitals[1].alpha = 1.24;
    orbitals[1].x = 0.0;
    orbitals[1].y = 0.0;
    orbitals[1].z = bond_length;
    
    // Nuclear positions (same as orbital positions for minimal basis)
    double nuclear_pos[6] = {
        orbitals[0].x, orbitals[0].y, orbitals[0].z,  // H1
        orbitals[1].x, orbitals[1].y, orbitals[1].z   // H2
    };
    
    printf("System setup:\n");
    printf("H atom 1 at (%.3f, %.3f, %.3f) Bohr\n", 
           orbitals[0].x, orbitals[0].y, orbitals[0].z);
    printf("H atom 2 at (%.3f, %.3f, %.3f) Bohr\n", 
           orbitals[1].x, orbitals[1].y, orbitals[1].z);
    printf("Bond length: %.3f Bohr (%.3f Angstrom)\n\n", 
           bond_length, bond_length * 0.529177);
    
    // Run Hartree-Fock calculation
    int iterations = hartree_fock(orbitals, nuclear_pos, 2);
    
    if (iterations > 0) {
        printf("\nCalculation completed successfully in %d iterations.\n", iterations);
    } else {
        printf("\nCalculation failed to converge.\n");
    }
    
    return 0;
}