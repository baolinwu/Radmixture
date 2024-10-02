#include <R.h>
#include <string.h>
#include <math.h>
#include <Rmath.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include "admix.h"



// project a K-dim vector Qa to the Simplex
//  s.t.: q_k\geq 1e-7; \sum_k q_k = 1
void proj_q(double *Q, double *Qa, int *K){
  double minp = 1e-7;
  double tau = 1-(*K)*minp;
  int k;
  for(k=0; k<*K; k++){
    Qa[k] = Qa[k]-minp;
  }
  memcpy(Q, Qa, (*K)*sizeof(double));  
  R_rsort(Qa, *K);
  int bget=0;
  double tsum=0, tmax;
  
  k = 0;
  while( (k<(*K-1))&&(bget<=0) ){
    tsum += Qa[(*K)-1-k];
    tmax = (tsum-tau)/(k+1);
    if(tmax>=Qa[(*K)-2-k]){
      bget = 1;
    }  else{
      k++;
    }
  }
  if( (bget<=0)&&(k>=(*K-1)) ){
    tmax = (tsum+Qa[0]-tau)/(*K);
  }
  for(k=0; k<*K; k++){
    Q[k] = fmax(Q[k]-tmax, 0);
  }
  for(k=0; k<*K; k++){
    Q[k] += minp;
  }
}




/*
  \min_{\beta} \beta^TD\beta/2 - \beta^Td
  s.t.  A^T\beta \geq b
  first m constraints are equality constraints.
*/
void solve_QP(double *sol, double *Dmat, double *dvec, double *Amat, double *bvec, int *meq, int*K, int *q,
	      double *work, int *ivec){
  int *iact, *nact, *iter, *ierr, q1=*q+1;
  ivec[0] = 0;
  ierr = ivec; iact = ivec+1;
  nact = ivec+q1; iter = ivec+q1+1;
  double *lagr, *crval;
  crval = work; lagr = work+1;
    
  F77_SUB(qpgen2)(Dmat, dvec,  K, K, sol,  lagr,crval,
		  Amat, bvec, K, q, meq,  iact, nact, iter, work+q1, ierr);
}

void Q_qp(double *sol, double *G, double *Q, double *P, int *n, int*m, int *K, int *i1,
	  double *Dm, double *Aq, double *bq, double *d1q, double *work, int*ivec){
  int j, k, k1;
  double pj, gj, pg;
  double *Pj, *Qi;
  Qi = Q+(*K)*(*i1);
  
  memset(d1q, 0, (*K)*sizeof(double));   memset(Dm, 0, (*K)*(*K)*sizeof(double));
  for(j=0; j<*m; j++){
    Pj = P + (*K)*j;
    gj = G[(*n)*j+*i1];
    pj = 0;
    for(k=0; k<*K; k++){
      pj += Qi[k]*Pj[k];
    }
    pg = (gj-2*pj)/pj/(1-pj);
    for(k=0; k<*K; k++){
      d1q[k] += pg*Pj[k];  //
    }
    pg = gj/pj/pj + (2-gj)/(1-pj)/(1-pj);
    for(k=0; k<*K; k++){
      Dm[(*K)*k+k] += pg*Pj[k]*Pj[k];
      if(k<(*K-1)){
	for(k1=(k+1); k1<*K; k1++){
	  pj = pg*Pj[k]*Pj[k1];
	  Dm[(*K)*k+k1] += pj;
	  Dm[(*K)*k1+k] += pj;
	}
      }
    }
  }
  bq[0] = 0;
  for(k=0; k<*K; k++){
    bq[k+1] = - Qi[k];
  }
  
  int meq = 1, q = *K+1;
  solve_QP(sol, Dm, d1q, Aq, bq, &meq, K,&q, work, ivec);
  // update
  if(ivec[0]<=0){
    for(k=0; k<*K; k++){
      bq[k] = Qi[k]+sol[k];
    }
    proj_q(Qi, bq, K);
  }
}


void P_qp(double *sol, double *G, double *Q, double *P, int *n, int*m, int *K, int *j1,
	  double *Dm, double *Ap, double *bp, double *d1p, double *work, int*ivec){
  int i, k, k1;
  double pi, gi, pg;
  double *Pj, *Qi;
  Pj = P+(*K)*(*j1);
  
  memset(d1p, 0, (*K)*sizeof(double));
  memset(Dm, 0, (*K)*(*K)*sizeof(double));
  for(i=0; i<*n; i++){
    Qi = Q + (*K)*i;
    gi = G[(*n)*(*j1)+i];
    pi = 0;
    for(k=0; k<*K; k++){
      pi += Qi[k]*Pj[k];
    }
    pg = (gi-2*pi)/pi/(1-pi);
    for(k=0; k<*K; k++){
      d1p[k] += pg*Qi[k];  //
    }
    pg = gi/pi/pi + (2-gi)/(1-pi)/(1-pi);
    // pg = 2/pi/(1-pi);
    for(k=0; k<*K; k++){
      Dm[(*K)*k+k] += pg*Qi[k]*Qi[k];
      if(k<(*K-1)){
	for(k1=(k+1); k1<*K; k1++){
	  pi = pg*Qi[k]*Qi[k1];
	  Dm[(*K)*k+k1] += pi;
	  Dm[(*K)*k1+k] += pi;
	}
      }
    }
  }
  for(k=0; k<*K; k++){
    bp[k] = - Pj[k];
    bp[k+(*K)] = Pj[k] - 1;
  }
  
  int meq = 0, q = 2*(*K);
  solve_QP(sol, Dm, d1p, Ap, bp, &meq, K,&q, work, ivec);
  // update
  if(ivec[0]<=0){
    for(k=0; k<*K; k++){
      Pj[k] = fmin(fmax(Pj[k]+sol[k], 1e-7), 1-1e-7);
    }
  }
}




void radmixture_qp(double *Q_ki, double *P_kj, double *LLK, int *N, int *M, double *G_ij, int *K, int *Kem,
		   int *itMax, double *eps, int *verb){
  int i, j, k, it, KN=(*K)*(*N), KM=(*K)*(*M);
  double *Ps, *Qs, *pks, delta;
  Ps = (double *) R_alloc(KM*6, sizeof(double));
  Qs = (double *) R_alloc(KN*6, sizeof(double));
  pks = (double *) R_alloc((*K)*2, sizeof(double));
  
  // ECM
  for(it=0; it<6; it++){
    for(i=0; i<*N; i++){
      Qem(Q_ki+(*K)*i, P_kj, G_ij, K,N,M,i, *Kem, pks);
    }
    for(j=0; j<*M; j++){
      Pem(P_kj+(*K)*j, Q_ki, G_ij+(*N)*j, N,K, *Kem, pks);
    }
    if(it<5){
      memcpy(Qs+KN*it, Q_ki, KN*sizeof(double));  memcpy(Ps+KM*it, P_kj, KM*sizeof(double));
    }
  }
  /* inital llk */
  double llk0, llkq, llk1;
  llk1 = LogLikelihood(P_kj, Q_ki, G_ij, N, M, K);
  LLK[0] = llk1;
  if(*verb>0) Rprintf("Initial LLK = %10.9f \n", llk1);
  // QN
  it = 0; delta=1;
  
  double *Dm,*Aq,*bq,*d1q,*work, *Ap;
  int *ivec;
  
  Dm = (double *) R_alloc((*K)*(*K), sizeof(double));
  Aq = (double *) R_alloc((*K)*(*K+1), sizeof(double));
  Ap = (double *) R_alloc((*K)*(*K)*2, sizeof(double));
  bq = (double *) R_alloc((*K)*2, sizeof(double));
  d1q = (double *) R_alloc(*K, sizeof(double));
  work = (double *) R_alloc(8*(*K)+2+(*K)*(*K+5)/2, sizeof(double));
  ivec = (int *) R_alloc(2*(*K)+4, sizeof(int));
  
  memset(Ap, 0, (*K)*(*K)*2*sizeof(double));   memset(Aq, 0, (*K)*(*K+1)*sizeof(double));
  for(k=0; k<*K; k++){
    Ap[(*K)*k+k] = 1;
    Ap[(*K)*(*K+k)+k] = -1;
    Aq[(*K)*(k+1)+k] = 1;
    Aq[k] = 1;
  }
  
  // QP/ECM/QN
  while( (it<*itMax) && (delta>*eps) ){
    it++;
    llk0 = llk1;
    // QP
    memcpy(Qs+5*KN, Q_ki, KN*sizeof(double));  memcpy(Ps+5*KM, P_kj, KM*sizeof(double));
    for(i=0; i<*N; i++){
      Q_qp(pks, G_ij,Q_ki,P_kj, N,M,K, &i, Dm,Aq,bq, d1q, work, ivec);
    }
    for(j=0; j<*M; j++){
      P_qp(pks, G_ij,Q_ki,P_kj, N,M,K, &j, Dm,Ap,bq, d1q, work, ivec);
    }
    llkq = LogLikelihood(P_kj, Q_ki, G_ij, N, M, K);
    if(llkq<llk0){
      // ECM
      memcpy(Q_ki, Qs+5*KN, KN*sizeof(double)); memcpy(P_kj, Ps+5*KM, KM*sizeof(double));
      for(i=0; i<*N; i++){
	Qem(Q_ki+(*K)*i, P_kj, G_ij, K, N, M, i, *Kem, pks);
      }
      for(j=0; j<*M; j++){
	Pem(P_kj+(*K)*j, Q_ki, G_ij+(*N)*j, N,K, *Kem, pks);
      }
      llkq = LogLikelihood(P_kj, Q_ki, G_ij, N, M, K);
    }
    
    // QN: 1-step
    QP_qn(Qs+KN*2, Qs+KN*3, Qs+KN*4, Qs+KN*5, Q_ki,  Ps+KM*2, Ps+KM*3, Ps+KM*4, Ps+KM*5, P_kj,  Qs, Ps, K, N, M, pks);
    llk1 = LogLikelihood(Ps, Qs, G_ij, N, M, K);
    if(llk1<llkq){
      llk1 = llkq;
    } else{
      memcpy(P_kj, Ps, KM*sizeof(double));  memcpy(Q_ki, Qs, KN*sizeof(double));
    }
    for(k=0; k<5; k++){
      memcpy(Qs+KN*k, Qs+KN*(k+1), KN*sizeof(double));  memcpy(Ps+KM*k, Ps+KM*(k+1), KM*sizeof(double));
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


