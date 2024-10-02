#include <R.h>
#include <string.h>
#include <math.h>
#include <Rmath.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include "admix.h"


void Qadmixture_qp(double *Q_ki, double *P_kj, double *LLK, int *N, int *M, double *G_ij, int *K, int *Kem,
		   int *itMax, double *eps, int *verb){
  int i, k, it, KN=(*K)*(*N);
  double *Qs, *pks, delta;
  Qs = (double *) R_alloc(KN*6, sizeof(double));
  pks = (double *) R_alloc((*K)*2, sizeof(double));
  
  // ECM
  for(it=0; it<6; it++){
    for(i=0; i<*N; i++){
      Qem(Q_ki+(*K)*i, P_kj, G_ij, K,N,M,i, *Kem, pks);
    }
    if(it<5){
      memcpy(Qs+KN*it, Q_ki, KN*sizeof(double));
    }
  }
  /* inital llk */
  double llk0, llkq, llk1;
  llk1 = LogLikelihood(P_kj, Q_ki, G_ij, N, M, K);
  LLK[0] = llk1;
  if(*verb>0) Rprintf("Initial LLK = %10.9f \n", llk1);
  // QN
  it = 0; delta=1;
  
  double *Dm,*Aq,*bq,*d1q,*work;
  int *ivec;
  
  Dm = (double *) R_alloc((*K)*(*K), sizeof(double));
  Aq = (double *) R_alloc((*K)*(*K+1), sizeof(double));
  bq = (double *) R_alloc((*K)*2, sizeof(double));
  d1q = (double *) R_alloc(*K, sizeof(double));
  work = (double *) R_alloc(8*(*K)+2+(*K)*(*K+5)/2, sizeof(double));
  ivec = (int *) R_alloc(2*(*K)+4, sizeof(int));
  
  memset(Aq, 0, (*K)*(*K+1)*sizeof(double));
  for(k=0; k<*K; k++){
    Aq[(*K)*(k+1)+k] = 1;
    Aq[k] = 1;
  }

  // QP/ECM/QN
  while( (it<*itMax) && (delta>*eps) ){
    it++;
    llk0 = llk1;
    // QP
    memcpy(Qs+5*KN, Q_ki, KN*sizeof(double));
    for(i=0; i<*N; i++){
      Q_qp(pks, G_ij,Q_ki,P_kj, N,M,K, &i, Dm,Aq,bq, d1q, work, ivec);
    }
    llkq = LogLikelihood(P_kj, Q_ki, G_ij, N, M, K);
    if(llkq<llk0){
      // ECM
      memcpy(Q_ki, Qs+5*KN, KN*sizeof(double));
      for(i=0; i<*N; i++){
	Qem(Q_ki+(*K)*i, P_kj, G_ij, K, N, M, i, *Kem, pks);
      }
      llkq = LogLikelihood(P_kj, Q_ki, G_ij, N, M, K);
      // Rprintf("\t %d: ECM > QP, \n", it); 
    }
    
    // QN
    Qa_qn(Qs+KN*2, Qs+KN*3, Qs+KN*4, Qs+KN*5, Q_ki, Qs,  K, N, M, pks);
    // Pa_qn(Ps+KM*2, Ps+KM*3, Ps+KM*4, Ps+KM*5, P_kj, Ps,  K, N, M);
    // QN: 1-step; too greedy?
    // QP_qn(Qs+KN*2, Qs+KN*3, Qs+KN*4, Qs+KN*5, Q_ki,  Ps+KM*2, Ps+KM*3, Ps+KM*4, Ps+KM*5, P_kj,  Qs, Ps, K, N, M, pks);
    llk1 = LogLikelihood(P_kj, Qs, G_ij, N, M, K);
    if(llk1<llkq){
      llk1 = llkq;
    } else{
      memcpy(Q_ki, Qs, KN*sizeof(double));
      // Rprintf("\t %d: QN!; QP=%lf, QN=%lf \n", it, llkq, llk1); 
    }
    for(k=0; k<5; k++){
      memcpy(Qs+KN*k, Qs+KN*(k+1), KN*sizeof(double));
    }
    delta = fabs(llk1-llk0)/(1+fabs(llk0));
    if(*verb>0){
      Rprintf("%d: LLK= %10.9f, Rel Diff= %10.9f \n", it, llk1, delta);
    }
    LLK[it] = llk1;
  }
  itMax[0] = it+1;
  Rprintf(" maxIter= %d; maxLLK= %10.9f; convDiff= %10.9f \n", it, llk1, delta);
}






