#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DIM 2
#define NELEC 2
#define EPSILON 1e-12

double ENUC = 1.323;

long eint(int a, int b, int c, int d) {
    long ab, cd, abcd;
    
    if (a > b) {
        ab = (long)a * (a + 1) / 2 + b;
    } else {
        ab = (long)b * (b + 1) / 2 + a;
    }
    
    if (c > d) {
        cd = (long)c * (c + 1) / 2 + d;
    } else {
        cd = (long)d * (d + 1) / 2 + c;
    }
    
    if (ab > cd) {
        abcd = ab * (ab + 1) / 2 + cd;
    } else {
        abcd = cd * (cd + 1) / 2 + ab;
    }
    
    return abcd;
}

double tei(int a, int b, int c, int d) {
    long key = eint(a, b, c, d);
    
    switch(key) {
        case 5:  return 0.7283;
        case 12: return 0.3418;
        case 14: return 0.2192;
        case 17: return 0.585;
        case 19: return 0.4368;
        case 20: return 0.9927;
        default: return 0.0;
    }
}

void makefock(double Hcore[DIM][DIM], double P[DIM][DIM], double F[DIM][DIM]) {
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            F[i][j] = Hcore[i][j];
            for (int k = 0; k < DIM; k++) {
                for (int l = 0; l < DIM; l++) {
                    F[i][j] = F[i][j] + P[k][l] * (tei(i+1,j+1,k+1,l+1) - 0.5*tei(i+1,k+1,j+1,l+1));
                }
            }
        }
    }
}

double currentenergy(double D[DIM][DIM], double Hcore[DIM][DIM], double F[DIM][DIM]) {
    double EN = 0.0;
    for (int mu = 0; mu < DIM; mu++) {
        for (int nu = 0; nu < DIM; nu++) {
            EN += 0.5 * D[mu][nu] * (Hcore[mu][nu] + F[mu][nu]);
        }
    }
    return EN;
}

double deltap(double D[DIM][DIM], double Dold[DIM][DIM]) {
    double DELTA = 0.0;
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            DELTA = DELTA + ((D[i][j] - Dold[i][j]) * (D[i][j] - Dold[i][j]));
        }
    }
    return sqrt(DELTA);
}

// Analytical eigenvalue/eigenvector solution for 2x2 symmetric matrix
void diagonalize_2x2(double mat[DIM][DIM], double eigenvalues[DIM], double eigenvectors[DIM][DIM]) {
    double a = mat[0][0];
    double b = mat[0][1];  // = mat[1][0] for symmetric matrix
    double c = mat[1][1];
    
    // Calculate eigenvalues using quadratic formula
    double trace = a + c;
    double det = a * c - b * b;
    double discriminant = trace * trace - 4 * det;
    
    if (discriminant < 0) {
        printf("Error: negative discriminant in eigenvalue calculation\n");
        exit(1);
    }
    
    eigenvalues[0] = (trace - sqrt(discriminant)) / 2.0;  // smaller eigenvalue
    eigenvalues[1] = (trace + sqrt(discriminant)) / 2.0;  // larger eigenvalue
    
    // Calculate eigenvectors
    // For eigenvalue λ₀ (smaller)
    if (fabs(b) > EPSILON) {
        double norm0 = sqrt((eigenvalues[0] - c) * (eigenvalues[0] - c) + b * b);
        eigenvectors[0][0] = (eigenvalues[0] - c) / norm0;
        eigenvectors[1][0] = b / norm0;
    } else if (fabs(eigenvalues[0] - a) < EPSILON) {
        eigenvectors[0][0] = 1.0;
        eigenvectors[1][0] = 0.0;
    } else {
        eigenvectors[0][0] = 0.0;
        eigenvectors[1][0] = 1.0;
    }
    
    // For eigenvalue λ₁ (larger)
    if (fabs(b) > EPSILON) {
        double norm1 = sqrt((eigenvalues[1] - c) * (eigenvalues[1] - c) + b * b);
        eigenvectors[0][1] = (eigenvalues[1] - c) / norm1;
        eigenvectors[1][1] = b / norm1;
    } else if (fabs(eigenvalues[1] - a) < EPSILON) {
        eigenvectors[0][1] = 1.0;
        eigenvectors[1][1] = 0.0;
    } else {
        eigenvectors[0][1] = 0.0;
        eigenvectors[1][1] = 1.0;
    }
    
    // Ensure orthogonality and correct signs to match NumPy convention
    // NumPy tends to prefer positive first components
    if (eigenvectors[0][0] < 0) {
        eigenvectors[0][0] = -eigenvectors[0][0];
        eigenvectors[1][0] = -eigenvectors[1][0];
    }
    if (eigenvectors[0][1] < 0) {
        eigenvectors[0][1] = -eigenvectors[0][1];
        eigenvectors[1][1] = -eigenvectors[1][1];
    }
}

void print_matrix(const char* name, double mat[DIM][DIM]) {
    printf("%s:\n", name);
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            printf("%12.8f ", mat[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

int main() {
    // Initialize matrices
    double S[DIM][DIM] = {{1.0000, 0.5017}, {0.5017, 1.0000}};
    double T[DIM][DIM] = {{0.6249, 0.2395}, {0.2395, 1.1609}};
    double V[DIM][DIM] = {{-2.2855, -1.5555}, {-1.5555, -3.4639}};
    
    double Hcore[DIM][DIM];
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            Hcore[i][j] = T[i][j] + V[i][j];
        }
    }
    
    printf("=== DEBUG INFORMATION ===\n");
    print_matrix("S matrix", S);
    print_matrix("Hcore matrix", Hcore);
    
    // Diagonalize S
    double SVAL[DIM];
    double SVEC[DIM][DIM];
    diagonalize_2x2(S, SVAL, SVEC);
    
    printf("S eigenvalues: %.4f %.4f\n", SVAL[0], SVAL[1]);
    print_matrix("S eigenvectors", SVEC);
    
    // Calculate S^(-1/2)
    double SVAL_minhalf[DIM][DIM] = {{0}};
    for (int i = 0; i < DIM; i++) {
        SVAL_minhalf[i][i] = pow(SVAL[i], -0.5);
    }
    
    printf("SVAL^(-1/2) diagonal: %.8f, %.8f\n\n", SVAL_minhalf[0][0], SVAL_minhalf[1][1]);
    
    // Calculate S_minhalf = SVEC * SVAL_minhalf * SVEC^T
    double temp[DIM][DIM], S_minhalf[DIM][DIM];
    
    // temp = SVEC * SVAL_minhalf
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            temp[i][j] = 0.0;
            for (int k = 0; k < DIM; k++) {
                temp[i][j] += SVEC[i][k] * SVAL_minhalf[k][j];
            }
        }
    }
    
    // S_minhalf = temp * SVEC^T
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            S_minhalf[i][j] = 0.0;
            for (int k = 0; k < DIM; k++) {
                S_minhalf[i][j] += temp[i][k] * SVEC[j][k];
            }
        }
    }
    
    print_matrix("S^(-1/2) matrix", S_minhalf);
    
    double P[DIM][DIM] = {{0}};
    double OLDP[DIM][DIM];
    double DELTA = 1.0;
    int count = 0;
    
    printf("=== SCF ITERATIONS ===\n");
    
    while (DELTA > 0.00001) {
        count++;
        
        // Make Fock matrix
        double F[DIM][DIM];
        makefock(Hcore, P, F);
        
        if (count == 1) {
            printf("\nIteration %d:\n", count);
            print_matrix("P matrix (density)", P);
            print_matrix("F matrix (Fock)", F);
        }
        
        // Calculate F' = S_minhalf^T * F * S_minhalf
        double temp1[DIM][DIM], Fprime[DIM][DIM];
        
        // temp1 = S_minhalf^T * F
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                temp1[i][j] = 0.0;
                for (int k = 0; k < DIM; k++) {
                    temp1[i][j] += S_minhalf[k][i] * F[k][j];
                }
            }
        }
        
        // Fprime = temp1 * S_minhalf
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                Fprime[i][j] = 0.0;
                for (int k = 0; k < DIM; k++) {
                    Fprime[i][j] += temp1[i][k] * S_minhalf[k][j];
                }
            }
        }
        
        if (count == 1) {
            print_matrix("F' matrix (transformed Fock)", Fprime);
        }
        
        // Diagonalize F'
        double E[DIM];
        double Cprime[DIM][DIM];
        diagonalize_2x2(Fprime, E, Cprime);
        
        if (count == 1) {
            printf("F' eigenvalues: %.8f %.8f\n", E[0], E[1]);
            print_matrix("F' eigenvectors (C')", Cprime);
        }
        
        // Back transform: C = S_minhalf * Cprime
        double C[DIM][DIM];
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                C[i][j] = 0.0;
                for (int k = 0; k < DIM; k++) {
                    C[i][j] += S_minhalf[i][k] * Cprime[k][j];
                }
            }
        }
        
        if (count == 1) {
            print_matrix("C matrix (back-transformed coefficients)", C);
        }
        
        // Make density matrix
        for (int mu = 0; mu < DIM; mu++) {
            for (int nu = 0; nu < DIM; nu++) {
                OLDP[mu][nu] = P[mu][nu];
                P[mu][nu] = 0.0;
                for (int m = 0; m < NELEC/2; m++) {
                    P[mu][nu] = P[mu][nu] + 2.0 * C[mu][m] * C[nu][m];
                }
            }
        }
        
        if (count == 1) {
            print_matrix("New P matrix", P);
        }
        
        DELTA = deltap(P, OLDP);
        
        double electronic_energy = currentenergy(P, Hcore, F);
        double total_energy = electronic_energy + ENUC;
        
        if (count == 1) {
            printf("Electronic energy: %.15f\n", electronic_energy);
            printf("Total energy: %.15f\n", total_energy);
            printf("DELTA: %.15f\n", DELTA);
        }
        
        printf("E = %.6f, N(SCF) = %d\n", total_energy, count);
    }
    
    double F_final[DIM][DIM];
    makefock(Hcore, P, F_final);
    double final_energy = currentenergy(P, Hcore, F_final) + ENUC;
    
    printf("SCF procedure complete, TOTAL E(SCF) = %.15f hartrees\n", final_energy);
    
    return 0;
}