#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <cblas.h>
#include <lapacke.h>

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
    // Copy Hcore to F
    cblas_dcopy(DIM * DIM, (double*)Hcore, 1, (double*)F, 1);
    
    // Add two-electron contributions
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            for (int k = 0; k < DIM; k++) {
                for (int l = 0; l < DIM; l++) {
                    F[i][j] += P[k][l] * (tei(i+1,j+1,k+1,l+1) - 0.5*tei(i+1,k+1,j+1,l+1));
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
    double diff[DIM][DIM];
    
    // Calculate difference matrix: diff = D - Dold
    cblas_dcopy(DIM * DIM, (double*)D, 1, (double*)diff, 1);
    cblas_daxpy(DIM * DIM, -1.0, (double*)Dold, 1, (double*)diff, 1);
    
    // Calculate Frobenius norm: sqrt(sum of squares)
    double norm = cblas_dnrm2(DIM * DIM, (double*)diff, 1);
    
    return norm;
}

// Diagonalize symmetric matrix using LAPACK, but match original analytical results
void diagonalize_symmetric(double mat[DIM][DIM], double eigenvalues[DIM], double eigenvectors[DIM][DIM]) {
    // Work arrays for LAPACK
    double work_matrix[DIM * DIM];
    double work_eigenvalues[DIM];
    
    // Copy matrix to work array (column-major for LAPACK)
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            work_matrix[j * DIM + i] = mat[i][j];  // Transpose for column-major
        }
    }
    
    // Call LAPACK symmetric eigenvalue solver
    int info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', DIM, work_matrix, DIM, work_eigenvalues);
    
    if (info != 0) {
        printf("Error: LAPACK dsyev failed with info = %d\n", info);
        exit(1);
    }
    
    // Copy eigenvalues
    for (int i = 0; i < DIM; i++) {
        eigenvalues[i] = work_eigenvalues[i];
    }
    
    // Copy eigenvectors back to row-major format
    // LAPACK returns eigenvectors as columns in column-major format
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            eigenvectors[i][j] = work_matrix[j * DIM + i];
        }
    }
    
    // Apply sign convention to match original implementation
    // Ensure first component of each eigenvector is positive
    for (int j = 0; j < DIM; j++) {
        if (eigenvectors[0][j] < 0) {
            for (int i = 0; i < DIM; i++) {
                eigenvectors[i][j] = -eigenvectors[i][j];
            }
        }
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
    
    // Hcore = T + V using BLAS
    cblas_dcopy(DIM * DIM, (double*)T, 1, (double*)Hcore, 1);
    cblas_daxpy(DIM * DIM, 1.0, (double*)V, 1, (double*)Hcore, 1);
    
    printf("=== DEBUG INFORMATION ===\n");
    print_matrix("S matrix", S);
    print_matrix("Hcore matrix", Hcore);
    
    // Diagonalize S
    double SVAL[DIM];
    double SVEC[DIM][DIM];
    diagonalize_symmetric(S, SVAL, SVEC);
    
    printf("S eigenvalues: %.4f %.4f\n", SVAL[0], SVAL[1]);
    print_matrix("S eigenvectors", SVEC);
    
    // Calculate S^(-1/2)
    double SVAL_minhalf[DIM][DIM] = {{0}};
    for (int i = 0; i < DIM; i++) {
        SVAL_minhalf[i][i] = pow(SVAL[i], -0.5);
    }
    
    printf("SVAL^(-1/2) diagonal: %.8f, %.8f\n\n", SVAL_minhalf[0][0], SVAL_minhalf[1][1]);
    
    // Calculate S_minhalf = SVEC * SVAL_minhalf * SVEC^T using BLAS
    double temp[DIM][DIM], S_minhalf[DIM][DIM];
    
    // temp = SVEC * SVAL_minhalf
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                DIM, DIM, DIM, 1.0, 
                (double*)SVEC, DIM, 
                (double*)SVAL_minhalf, DIM, 
                0.0, (double*)temp, DIM);
    
    // S_minhalf = temp * SVEC^T
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 
                DIM, DIM, DIM, 1.0, 
                (double*)temp, DIM, 
                (double*)SVEC, DIM, 
                0.0, (double*)S_minhalf, DIM);
    
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
        
        // Calculate F' = S_minhalf^T * F * S_minhalf using BLAS
        double temp1[DIM][DIM], Fprime[DIM][DIM];
        
        // temp1 = S_minhalf^T * F
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 
                    DIM, DIM, DIM, 1.0, 
                    (double*)S_minhalf, DIM, 
                    (double*)F, DIM, 
                    0.0, (double*)temp1, DIM);
        
        // Fprime = temp1 * S_minhalf
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                    DIM, DIM, DIM, 1.0, 
                    (double*)temp1, DIM, 
                    (double*)S_minhalf, DIM, 
                    0.0, (double*)Fprime, DIM);
        
        if (count == 1) {
            print_matrix("F' matrix (transformed Fock)", Fprime);
        }
        
        // Diagonalize F'
        double E[DIM];
        double Cprime[DIM][DIM];
        diagonalize_symmetric(Fprime, E, Cprime);
        
        if (count == 1) {
            printf("F' eigenvalues: %.8f %.8f\n", E[0], E[1]);
            print_matrix("F' eigenvectors (C')", Cprime);
        }
        
        // Back transform: C = S_minhalf * Cprime using BLAS
        double C[DIM][DIM];
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                    DIM, DIM, DIM, 1.0, 
                    (double*)S_minhalf, DIM, 
                    (double*)Cprime, DIM, 
                    0.0, (double*)C, DIM);
        
        if (count == 1) {
            print_matrix("C matrix (back-transformed coefficients)", C);
        }
        
        // Save old density matrix
        cblas_dcopy(DIM * DIM, (double*)P, 1, (double*)OLDP, 1);
        
        // Make new density matrix P = 2 * C_occ * C_occ^T using BLAS
        // For closed shell with NELEC/2 occupied orbitals
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 
                    DIM, DIM, NELEC/2, 2.0, 
                    (double*)C, DIM, 
                    (double*)C, DIM, 
                    0.0, (double*)P, DIM);
        
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