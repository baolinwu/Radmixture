#include <R.h>
#include <string.h>
#include <math.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>   // definitions of linear solvers
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include "admix.h"

/* 
  Admixture model
    for sample i SNP j: G_{ij} ~ Binom(2,p_{ij}) with p_{ij}=\sum_k Q_{ik}P_{kj}
    assume working linkage equilibrium assumption; some pruning is recommended (see Admixture manual for more details)
*/

/* compute binomial log likelihood */
// (1) ignore Binom coef const (2, g); (2) imputed SNP dosage OK.
double Lbinom(double p, double g){
  // gate-keeping: it happens for (p=0, g=1,2) or (p=1, g=0,1) ....
  if( (p<=0)&(g>0) ) return -50;
  if( (p>=1)&(g<2) ) return -50;
  
  double llk;
  if(g<=0){
    llk = 2*log(1-p); // (1-p)*(1-p);
  } 
  if(g>=2){
    llk = 2*log(p); // p*p;
  }
  if( (g>0)&&(g<2) ){
    llk = g*log(p)+(2-g)*log(1-p); // pow(p,g)*pow(1-p,2-g);
  }
  return llk;
}

/* Marginal log likelihood for one given SNP (across N samples) */
double Pllk(double *P_j, double *Q_ki, double *G_j, int *N, int *K){
  int i, k;
  double g, llk=0, p, *qi;
  
  for(i=0; i<*N; i++){
    g = G_j[i]; 
    if(g>=0){
      p = 0;  qi = Q_ki+(*K)*i;
      for(k=0; k<*K; k++){
        p += qi[k]*P_j[k];
      }
      p = fmax(fmin(p,1-1e-16),1e-16);  // safer for convergence; avoid being trapped at boundary
      llk += Lbinom(p,g);
    }
  }
  return llk;
}


/* Marginal log likelihood for all SNPs (across N samples) */
double LogLikelihood(double *P_kj, double *Q_ki, double *G, int *N, int *M, int *K){
  int j; double llk=0;
  for(j=0; j<*M; j++){
    llk += Pllk(P_kj+(*K)*j, Q_ki, G+(*N)*j, N, K);
  }
  return llk;
}


/* 
  linear equation solving with lapack : Q^{-1}V, Q (general; symmetric not needed). 
   Refer to R_ext/Lapack.h 
*/
// compute V <= Q^{-1}V; Q: KxK square matrix; ipiv: K-dim vector
void Lsolve(double *V, double *Q, int K, int *ipiv){
  char *trans = "N";
  int dim=1, info=-2, K1=K;
  F77_CALL(dgetrf)(&K1,&K1, Q, &K1, ipiv, &info);
  F77_CALL(dgetrs)(trans, &K1, &dim, Q, &K1, ipiv, V, &K1, &info  FCONE);
}




// update all Q's in one-step
void Qa_qn(double *Q1, double *Q2, double *Q3, double *Q4, double *Q5, double *Qa, // double *P_kj, double *G,
	   int *K, int *N, int *M, double *Q6){
  int i, k, ii, ik;
  // QN3 acc
  double u1, u2, u3, v1, v2, v3, UV[9], UX[3];  int ipiv[3];
  memset(UV, 0, 9*sizeof(double)); memset(UX, 0, 3*sizeof(double));
  for(i=0; i<*N; i++){
    ii = (*K)*i;
    for(k=0; k<*K; k++){
      ik = ii+k;
      u1 = Q2[ik]-Q1[ik];  u2 = v1 = Q3[ik]-Q2[ik]; u3 = v2 = Q4[ik]-Q3[ik]; v3 = Q5[ik]-Q4[ik];
      UX[0] += u1*v3/(*N);  UX[1] += u2*v3/(*N);  UX[2] += u3*v3/(*N);
      UV[0] += u1*u1/(*N) - u1*v1/(*N); UV[1] += u2*u1/(*N) - u2*v1/(*N); UV[2] += u3*u1/(*N) - u3*v1/(*N);
      UV[3] += u1*u2/(*N) - u1*v2/(*N); UV[4] += u2*u2/(*N) - u2*v2/(*N); UV[5] += u3*u2/(*N) - u3*v2/(*N);
      UV[6] += u1*u3/(*N) - u1*v3/(*N); UV[7] += u2*u3/(*N) - u2*v3/(*N); UV[8] += u3*u3/(*N) - u3*v3/(*N);
    }
  }
  Lsolve(UX, UV, 3, ipiv);
  for(i=0; i<*N; i++){
    ii = (*K)*i;
    for(k=0; k<*K; k++){
      ik = ii+k;
      v1 = Q3[ik]-Q2[ik]; v2 = Q4[ik]-Q3[ik]; v3 = Q5[ik]-Q4[ik];
      Q6[k] = Q5[ik] + v1*UX[0] + v2*UX[1] + v3*UX[2];
    }
    proj_q(Qa+ii, Q6, K);
  }
}


void Pa_qn(double *P1, double *P2, double *P3, double *P4, double *P5, double *Pa, // double *Q_ki, double *G,
	   int *K, int *N, int *M){
  int j, k, jj, jk; 
  // QN3 acc
  double u1, u2, u3, v1, v2, v3, UV[9], UX[3];  int ipiv[3];
  memset(UV, 0, 9*sizeof(double)); memset(UX, 0, 3*sizeof(double));
  for(j=0; j<*M; j++){
    jj = (*K)*j;
    for(k=0; k<*K; k++){
      jk = jj+k;
      u1 = P2[jk]-P1[jk];  u2 = v1 = P3[jk]-P2[jk]; u3 = v2 = P4[jk]-P3[jk]; v3 = P5[jk]-P4[jk];
      UX[0] += u1*v3/(*M);  UX[1] += u2*v3/(*M);  UX[2] += u3*v3/(*M);
      UV[0] += u1*u1/(*M) - u1*v1/(*M); UV[1] += u2*u1/(*M) - u2*v1/(*M); UV[2] += u3*u1/(*M) - u3*v1/(*M);
      UV[3] += u1*u2/(*M) - u1*v2/(*M); UV[4] += u2*u2/(*M) - u2*v2/(*M); UV[5] += u3*u2/(*M) - u3*v2/(*M);
      UV[6] += u1*u3/(*M) - u1*v3/(*M); UV[7] += u2*u3/(*M) - u2*v3/(*M); UV[8] += u3*u3/(*M) - u3*v3/(*M);
    }
  }
  Lsolve(UX, UV, 3, ipiv);
  double a1 = 0;
  for(j=0; j<*M; j++){
    jj = (*K)*j;
    for(k=0; k<*K; k++){
      jk = jj+k;
      v1 = P3[jk]-P2[jk]; v2 = P4[jk]-P3[jk]; v3 = P5[jk]-P4[jk];
      a1 = v1*UX[0] + v2*UX[1] + v3*UX[2];
      Pa[jk] = fmin(fmax(P5[jk] + a1, 1e-7), 1-1e-7);
    }
  }
}





// Update Q/P all in once
void QP_qn(double *Q1, double *Q2, double *Q3, double *Q4, double *Q5,
	   double *P1, double *P2, double *P3, double *P4, double *P5,
	   double *Qa, double *Pa,  int *K, int *N, int *M, double *Q6){
  int i, k, ii, ik, j, jj, jk;
  double u1, u2, u3, v1, v2, v3, UV[9], UX[3], NM=*N+*M;
  int ipiv[3];
  memset(UV, 0, 9*sizeof(double)); memset(UX, 0, 3*sizeof(double));
  for(i=0; i<*N; i++){
    ii = (*K)*i;
    for(k=0; k<*K; k++){
      ik = ii+k;
      u1 = Q2[ik]-Q1[ik];  u2 = v1 = Q3[ik]-Q2[ik]; u3 = v2 = Q4[ik]-Q3[ik]; v3 = Q5[ik]-Q4[ik];
      UX[0] += u1*v3/NM;  UX[1] += u2*v3/NM;  UX[2] += u3*v3/NM;
      UV[0] += u1*u1/NM - u1*v1/NM; UV[1] += u2*u1/NM - u2*v1/NM; UV[2] += u3*u1/NM - u3*v1/NM;
      UV[3] += u1*u2/NM - u1*v2/NM; UV[4] += u2*u2/NM - u2*v2/NM; UV[5] += u3*u2/NM - u3*v2/NM;
      UV[6] += u1*u3/NM - u1*v3/NM; UV[7] += u2*u3/NM - u2*v3/NM; UV[8] += u3*u3/NM - u3*v3/NM;
    }
  }
  for(j=0; j<*M; j++){
    jj = (*K)*j;
    for(k=0; k<*K; k++){
      jk = jj+k;
      u1 = P2[jk]-P1[jk];  u2 = v1 = P3[jk]-P2[jk]; u3 = v2 = P4[jk]-P3[jk]; v3 = P5[jk]-P4[jk];
      UX[0] += u1*v3/NM;  UX[1] += u2*v3/NM;  UX[2] += u3*v3/NM;
      UV[0] += u1*u1/NM - u1*v1/NM; UV[1] += u2*u1/NM - u2*v1/NM; UV[2] += u3*u1/NM - u3*v1/NM;
      UV[3] += u1*u2/NM - u1*v2/NM; UV[4] += u2*u2/NM - u2*v2/NM; UV[5] += u3*u2/NM - u3*v2/NM;
      UV[6] += u1*u3/NM - u1*v3/NM; UV[7] += u2*u3/NM - u2*v3/NM; UV[8] += u3*u3/NM - u3*v3/NM;
    }
  }
  Lsolve(UX, UV, 3, ipiv);
  for(i=0; i<*N; i++){
    ii = (*K)*i;
    for(k=0; k<*K; k++){
      ik = ii+k;
      v1 = Q3[ik]-Q2[ik]; v2 = Q4[ik]-Q3[ik]; v3 = Q5[ik]-Q4[ik];
      Q6[k] = Q5[ik] + v1*UX[0] + v2*UX[1] + v3*UX[2];
    }
    proj_q(Qa+ii, Q6, K);
  }
  for(j=0; j<*M; j++){
    jj = (*K)*j;
    for(k=0; k<*K; k++){
      jk = jj+k;
      v1 = P3[jk]-P2[jk]; v2 = P4[jk]-P3[jk]; v3 = P5[jk]-P4[jk];
      u1 = v1*UX[0] + v2*UX[1] + v3*UX[2];
      Pa[jk] = fmin(fmax(P5[jk] + u1, 1e-7), 1-1e-7);
    }
  }
}




