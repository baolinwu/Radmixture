#' Estimate the ancestry proportions in an admixture model with known ancestry allele freqs
#'
#' Given n unrelated individuals with m SNP markers, coming from K ancestry populations,
#' the admixture model has parameters
#' (1) Q, an nxK matrix of ancestry proportions for n individuals; (2) P, a Kxm matrix of ancestry population allele freq for m SNPs.
#' The individual sample allele freq is modeled as \eqn{\theta=PQ} and the genotypes are modeled with a Binomial distribution with probability \eqn{\theta}.
#' Here we assume P is given and try to estimate Q for each individual. An EM algorithm with quasi-Newton acceleration is implemented.
#' 
#' @param G genotype matrix of n by m
#' @param P Kxm matrix of allele freqs for K ancestry populations
#' @param pinit list of starting values of Q (nxK ancestry proportions).  If not provided, a k-means clustering will be used to provide an intial estimate as starting values. 
#' @param Kem  number of EM iterations within each ECME step
#' @param maxit maximum number of iterations
#' @param rtol convergence threshold for relative change of log likelihood
#' @param trace whether monitoring convergence
#' @param ... additional inputs passed to kmeans() to compute the initial parameter values
#' @return
#' \describe{
#'   \item{Q}{  nxK matrix of ancestry proportions }
#'   \item{Cs}{ estimated ancestry indicators }
#'   \item{llk}{ sequence of computed log likelihoods over all iterations }
#' }
#'
#' @references
#' X. Wu and B. Wu (2024) Radmixture: an R package for fast and robust likelihood-based estimation of admixture model for global ancestry inference. tech rep.
#'
#' D.H. Alexander, J. Novembre, and K. Lange (2009) Fast model-based estimation of ancestry in unrelated individuals. Genome Research, 19:1655-1664.
#'
#' H. Zhou, D.H. Alexander, and K. Lange (2011). A quasi-Newton method for accelerating the convergence of iterative optimization algorithms. Statistics and Computing, 21:261-273.
#' @export
#' @examples
#' n=100; m = 1000; K=3; P = t(sapply(1:K, function(j) runif(m)))
#' Q = matrix(runif(n*K), n,K); Q = Q/rowSums(Q)
#' tha = Q%*%P; G = matrix(rbinom(n*m, 2, tha), n,m)
#' a1 = Qadmixture(G,P, trace=TRUE)
#' apply(Q, 2, range); apply(a1$Q, 2, range); cor(Q, a1$Q)
#' table(a1$Cs, apply(Q, 1, which.max))
#' par(mfrow=c(1,2))
#' barplot(t(a1$Q), col=rainbow(K), xlab="Individual #", ylab="Ancestry", border=NA)
#' matplot(a1$Q, Q); abline(0,1, col=K+1, lwd=2, lty=2)
#' Y = apply(Q, 1, which.max)
#' ml = nnet::multinom(Y~a1$Q[,-1])
#' Yc = apply(ml$fitted, 1, which.max)
#' table(Y, Yc); mean(Yc!=Y)
Qadmixture  <-  function(G, P, pinit=NULL, Kem=5, maxit=1e2,rtol=1e-7,trace=FALSE, ...){
    n = dim(G)[1]; m = dim(G)[2]; K = dim(P)[1]
    if(is.null(pinit)){
        cl = kmeans(scale(G), centers=P, ...)$cluster
        Q.ki = matrix(table(cl)/length(cl), K,n)
    } else{
        Q.ki = t(pinit$Q)
    }
    
    obj = .C("Qadmixture_qp",
             Q.ki = as.double(Q.ki), as.double(P), llk = double(maxit),
             as.integer(n), as.integer(m), as.double(G), as.integer(K), as.integer(Kem),
             itMax=as.integer(maxit), as.double(rtol), as.integer(trace), PACKAGE='Radmixture')
    Q.ki = matrix(obj$Q.ki, K,n);  Cs = apply(Q.ki, 2, which.max);  llk = obj$llk[1:min(obj$itMax,maxit)]
    return( list(Q=t(Q.ki), Cs=Cs, llk=llk) )
}

