#include <R.h>
#include <string.h>
#include <math.h>
#include <Rmath.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>


/* marginal log likelihood for allele j, sample i */
double Lbinom(double p, double g);
double Pllk(double *P_j, double *Q_ki, double *G_j, int *n, int *K);
double LogLikelihood(double *P_kj, double *Q_ki, double *G, int *N, int *M, int *K);
void Lsolve(double *V, double *Q, int K, int *ipiv);

/* individual mixing proportion estimation */
void Qem(double *Q_i, double *P_kj, double *G_ij, int *K, int *N, int *M, int ii, int itMax, double *qi);
/* allele freq estimation */
void Pem(double *P_j, double *Q_ki, double *G_j, int *N, int *K, int itMax, double *pks);






void F77_SUB(qpgen2)
  (double *dmat, double *dvec, int *fddmat, int *n,
   double *sol, double *lagr, double *crval,
   double *amat, double *bvec, int *fdamat, int *q,
   int *meq, int *iact, int *nact, int *iter,
   double *work, int *ierr);


// QN
void Pa_qn(double *P1, double *P2, double *P3, double *P4, double *P5, double *Pa,
	   int *K, int *N, int *M);
void Qa_qn(double *Q1, double *Q2, double *Q3, double *Q4, double *Q5, double *Qa, 
	   int *K, int *N, int *M, double *Q6);
void QP_qn(double *Q1, double *Q2, double *Q3, double *Q4, double *Q5,
	   double *P1, double *P2, double *P3, double *P4, double *P5,
	   double *Qa, double *Pa,  int *K, int *N, int *M, double *Q6);


// QP
void proj_q(double *Q, double *Qa, int *K);
void solve_QP(double *sol, double *Dmat, double *dvec, double *Amat, double *bvec, int *meq, int*K, int *q,
	      double *work, int *ivec);
void Q_qp(double *sol, double *G, double *Q, double *P, int *n, int*m, int *K, int *i1,
	  double *Dm, double *Aq, double *bq, double *d1q, double *work, int*ivec);
void P_qp(double *sol, double *G, double *Q, double *P, int *n, int*m, int *K, int *j1,
	  double *Dm, double *Ap, double *bp, double *d1p, double *work, int*ivec);

void radmixture_qp(double *Q_ki, double *P_kj, double *LLK, int *N, int *M, double *G_ij, int *K, int *Kem,
		   int *itMax, double *eps, int *verb);

// PQP
void Q_pqp(double *sol, double *G, double *Q, double *P, int *n, int*m, int *K, int *i1,
	   double *Dm, double *Aq, double *bq, double *d1q, double *work, int*ivec);
void Q_pem(double *Q_i, double *P_kj, double *G_ij, int *K, int *N, int *M, int ii, int itMax, double *qi);
double PenLLK(double *P_kj, double *Q_ki, double *G, int *N, int *M, int *K);
void radmixture_pqp(double *Q_ki, double *P_kj, double *LLK, int *N, int *M, double *G_ij, int *K, int *Kem,
		    int *itMax, double *eps, int *verb);

// ss



