#include <stdio.h>
#include <string.h>
#include <math.h>
#include <lapacke.h>

#define N 2
#define MAX_ITER 50
#define THRESH 1e-8

// Integrals from Szabo & Ostlund Table 3.5
#define S12     0.6593
#define H11    -1.1204
#define H12    -0.9584

// Two-electron integrals (chemist's notation)
#define ERI_1111 0.7746
#define ERI_1112 0.5697
#define ERI_1122 0.4441
#define ERI_1212 0.3855
#define ERI_1222 0.5697

#define ENN    0.7143    // 1/R, R=1.4

int main() {
    double S[N][N] = {{1.0, S12},{S12,1.0}};
    double H[N][N] = {{H11, H12},{H12, H11}};
    double P[N][N] = {{0,0},{0,0}};
    double F[N][N], C[N][N], E[N];
    double Cprime[N][N], X[N][N]={0};
    double eig[N], work[N*N];
    double Eold = 0, Etot = 0, dE = 0;
    
    // Orthogonalize S (Lowdin)
    memcpy(work, S, sizeof(S));
    LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', N, work, N, eig);
    for(int i=0;i<N;i++) for(int j=0;j<N;j++)
        for(int k=0;k<N;k++)
            X[i][j] += work[i*N+k]/sqrt(eig[k])*work[j*N+k];

    for(int it=1; it<=MAX_ITER; ++it) {
        // Build two-electron part of Fock
        // G matrix: G(mu,nu) = sum P(lambda,sigma)[(m n|l s) - 0.5 (m s|l n)]
        double G[N][N] = {{0,0},{0,0}};

        // All explicit for 2x2
        // (mu,nu) = (0,0)
        G[0][0] = P[0][0]*(ERI_1111 - 0.5*ERI_1111)
                 +P[0][1]*(ERI_1112 - 0.5*ERI_1112)
                 +P[1][0]*(ERI_1112 - 0.5*ERI_1212)
                 +P[1][1]*(ERI_1122 - 0.5*ERI_1222);
        // (0,1) and (1,0)
        G[0][1] = P[0][0]*(ERI_1112 - 0.5*ERI_1212)
                 +P[0][1]*(ERI_1122 - 0.5*ERI_1222)
                 +P[1][0]*(ERI_1212 - 0.5*ERI_1112)
                 +P[1][1]*(ERI_1222 - 0.5*ERI_1122);
        G[1][0] = G[0][1];
        // (1,1)
        G[1][1] = P[0][0]*(ERI_1122 - 0.5*ERI_1222)
                 +P[0][1]*(ERI_1222 - 0.5*ERI_1122)
                 +P[1][0]*(ERI_1222 - 0.5*ERI_1122)
                 +P[1][1]*(ERI_1111 - 0.5*ERI_1111);

        // Fock matrix = Hcore + G
        for(int i=0;i<N;i++) for(int j=0;j<N;j++) F[i][j]=H[i][j]+G[i][j];
        
        // F' = X^T F X
        double FX[N][N] = {{0,0},{0,0}}, Fprime[N][N]={{0,0},{0,0}};
        for(int i=0;i<N;i++) for(int j=0;j<N;j++)
            for(int k=0;k<N;k++) FX[i][j] += F[i][k]*X[k][j];
        for(int i=0;i<N;i++) for(int j=0;j<N;j++)
            for(int k=0;k<N;k++) Fprime[i][j] += X[k][i]*FX[k][j];

        // Diagonalize Fprime
        memcpy(Cprime, Fprime, sizeof(Fprime));
        LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', N, Cprime[0], N, E);

        // C = X * Cprime
        for(int i=0;i<N;i++) for(int j=0;j<N;j++) {
            C[i][j]=0;
            for(int k=0;k<N;k++) C[i][j] += X[i][k]*Cprime[k][j];
        }

        // New density matrix: closed shell, only occ = 0
        double Pn[N][N] = {{0,0},{0,0}};
        for(int mu=0;mu<N;mu++) for(int nu=0;nu<N;nu++)
            Pn[mu][nu] = 2.0 * C[mu][0]*C[nu][0];

        // Electronic energy
        double Eelec=0;
        for(int mu=0;mu<N;mu++) for(int nu=0;nu<N;nu++)
            Eelec += 0.5*Pn[mu][nu]*(H[mu][nu]+F[mu][nu]);
        Etot = Eelec + ENN;

        dE = fabs(Etot-Eold);
        printf("Iter %2d: E = % .7f  dE = %.1e\n", it, Etot, dE);
        if(dE<THRESH) break;

        memcpy(P, Pn, sizeof(P));
        Eold=Etot;
    }

    printf("\nConverged SCF Energy for H2 (R=1.4 bohr) = % .7f hartree\n", Etot);
    return 0;
}