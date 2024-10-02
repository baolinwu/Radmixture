#include <R.h>
#include <Rinternals.h>
#include <string.h>
#include <math.h>
#include <Rmath.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include "admix.h"


SEXP Radmixture(SEXP Q, SEXP P, SEXP G, SEXP LLK, SEXP N, SEXP M, SEXP K, SEXP Kem,
		SEXP maxit, SEXP rtol, SEXP trace){
  double *Q_ki, *P_kj, *X, *pLLK, *eps;
  int *n, *m, *pK, *pKem, *itMax, *verb;
  SEXP ret;
  
  PROTECT(Q); PROTECT(P); PROTECT(LLK); PROTECT(maxit);
  
  //
  Q_ki = REAL(Q);
  P_kj = REAL(P);
  X = REAL(G);
  pLLK = REAL(LLK);
  n = INTEGER(N);
  m = INTEGER(M);
  pK = INTEGER(K);
  pKem = INTEGER(Kem);
  itMax = INTEGER(maxit);
  eps = REAL(rtol);
  verb = INTEGER(trace);
  
  radmixture_qp(Q_ki,P_kj,pLLK, n,m,X, pK,pKem, itMax, eps, verb);
  
  ret = PROTECT(allocVector(VECSXP,4));
  SET_VECTOR_ELT(ret,0,Q);
  SET_VECTOR_ELT(ret,1,P);
  SET_VECTOR_ELT(ret,2,LLK);
  SET_VECTOR_ELT(ret,3,maxit);
  UNPROTECT(5);
  return ret;
}

