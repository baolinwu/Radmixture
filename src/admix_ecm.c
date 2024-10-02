#include <R.h>
#include <string.h>
#include <math.h>
#include <Rmath.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
/* 
  Admixture model
    Binomial model for genotype of sample i SNP j: G_{ij} ~ Binom(2,p_{ij}) with p_{ij}=\sum_k Q_{ik}P_{kj}
    assume working linkage equilibrium assumption; some LD pruning is recommended (see Admixture manual for more details)
*/


/************************************/
/*   EM nested within ECME          */
/* individual sample mixing proportion estimation: Q */
void Qem(double *Q_i, double *P_kj, double *G_ij, int *K, int *N, int *M, int ii, int itMax, double *qi){
  int j, k, it=0, j1;
  double gj, a1, *pj, delta, q0;
  
  while(it<itMax){
    it++;
    memset(qi, 0, (*K)*sizeof(double));
    for(j=0; j<*M; j++){
      j1 = (*N)*j+ii;  gj = G_ij[j1];
      if(gj>=0){
        a1 = 0;   pj = P_kj+(*K)*j;
	for(k=0; k<*K; k++){
	  a1 += Q_i[k]*pj[k];
	}
	for(k=0; k<*K; k++){
	  if(a1>0) qi[k] = qi[k]+gj*Q_i[k]*pj[k]/a1;
	  if(a1<1) qi[k] = qi[k]+(2-gj)*Q_i[k]*(1-pj[k])/(1-a1);
	}
      }
    }
    // better to re-normalize \sum_k Qk=1; otherwise sum >1 might happen.
    a1 = 0;
    for(k=0; k<*K; k++){
      a1 += qi[k];
    }
    gj=0; delta = 0;
    for(k=0; k<*K; k++){
      q0 = Q_i[k];
      Q_i[k] = fmax(qi[k]/a1, 1e-7); // Q_i[k] = qi[k]/a1; ... thresholding; improve global convergence
      gj += Q_i[k];
      delta += fabs(Q_i[k]-q0);
    }
    for(k=0; k<*K; k++){
      Q_i[k] /= gj;
    }
    if(delta<1e-5){
      it = itMax;
    }
  }
}

/* allele freq estimation */
void Pem(double *P_j, double *Q_ki, double *G_j, int *N, int *K, int itMax, double *pks){
  int i, k, it=0; 
  double gi, a1, *qi, *pk1, *pk2;
  pk1 = pks; pk2 = pks+(*K);
  
  while(it<itMax){
    it++;
    memset(pks, 0, 2*(*K)*sizeof(double));
    for(i=0; i<*N; i++){
      gi = G_j[i];
      if(gi>=0){
        a1 = 0;  qi = Q_ki+(*K)*i;
	for(k=0; k<*K; k++){
	  a1 += qi[k]*P_j[k];
	}
	for(k=0; k<*K; k++){
          if(a1>0) pk1[k] = pk1[k]+gi*qi[k]/a1;  
          if(a1<1) pk2[k] = pk2[k]+(2-gi)*qi[k]/(1-a1);
	}
      }
    }
    a1 = 0;
    for(k=0; k<*K; k++){
      gi = P_j[k];
      P_j[k] = P_j[k]*pk1[k]/(P_j[k]*pk1[k]+(1-P_j[k])*pk2[k]);
      P_j[k] = fmin(fmax(P_j[k],1e-7), 1-1e-7);  // .... thresholding; improve global convergence
      a1 += fabs(P_j[k]-gi);
    }
    if(a1<1e-5){
      it = itMax;
    }
  }
}

