#' An efficient and robust admixture model fitting to large-scale SNP data
#'
#' Given n unrelated individuals with m SNP markers, coming from K ancestry populations,
#' we implement a quasi-Newton (with 3 secant conditions) accelerated quadratic program/ECME algorithms to compute the maximum likelihood estimates of an admixture model, with parameters
#' (1) Q: an nxK matrix of ancestry proportions for n individuals; (2) P: a Kxm matrix of ancestry population allele frequencies across m SNPs.
#' The individual allele frequencies are modeled as \eqn{\theta=QP} and the genotypes are modeled with Binomial distributions with probability \eqn{\theta}.
#'
#' @param G genotype matrix of n by m
#' @param K number of ancestry populations
#' @param pinit list of starting values of parameters: Q (nxK ancestry proportions); P (Kxm ancestry population allele freq). If not provided, a k-means clustering will be used to provide an intial estimate as starting values.
#' @param Kem  number of EM iterations within each ECME step
#' @param maxit maximum number of iterations
#' @param rtol convergence threshold for relative change of log likelihood
#' @param trace whether monitoring convergence
#' @param ... additional inputs passed to kmeans() to compute the initial parameter values
#' @return
#' \describe{
#'   \item{Q}{  nxK matrix of ancestry proportions }
#'   \item{P}{ Kxm matrix of ancestry population allele frequencies }
#'   \item{Cs}{ estimated ancestry indicators }
#'   \item{llk}{ sequence of computed log likelihoods over all iterations }
#' }
#'
#' @export
#'
#' @references
#' Wu X and Wu B (2024) Radmixture: an R package for fast and robust likelihood-based estimation of admixture model for global ancestry inference. tech rep.
#'
#' Ko S, Chu BB, Peterson D, Okenwa C, Papp JC, Alexander DH, et al. (2023) Unsupervised discovery of ancestry-informative markers and genetic admixture proportions in biobank-scale datasets. Am J Hum Genet. 2023 Feb 2;110(2):314–25.
#'
#' D.H. Alexander, J. Novembre, and K. Lange (2009) Fast model-based estimation of ancestry in unrelated individuals. Genome Research, 19:1655-1664.
#'
#' H. Zhou, D.H. Alexander, and K. Lange (2011). A quasi-Newton method for accelerating the convergence of iterative optimization algorithms. Statistics and Computing, 21:261-273.
#' @examples
#' n=1000; K=4; m1 = 3000; m0 = 5000; m=m0+m1
#' P1 = t(sapply(1:K, function(j) runif(m1)))
#' P0 = matrix(runif(m0), K,m0, byrow=TRUE)
#' P = cbind(P1,P0)
#' Q = matrix(runif(n*K), n,K); Q = Q/rowSums(Q)
#' tha = Q%*%P; G = matrix(rbinom(n*m, 2, tha), n,m)
#' jj = which(apply(G, 2, sd)>0)
#' G = G[,jj]; P = P[,jj]
#' a1 = Radmixture(G,K=K, Kem=10,maxit=1e3,rtol=1e-7,trace=TRUE)
#' apply(Q, 2, range); apply(a1$Q, 2, range); cor(Q, a1$Q)
#' table(a1$Cs, apply(Q, 1, which.max))
#' par(mfrow=c(1,2))
#' barplot(t(a1$Q), col=rainbow(K), xlab="Individual #", ylab="Ancestry", border=NA)
#' ik = apply(abs(cor(Q, a1$Q)), 1, which.max); ik
#' cor(t(a1$P[ik,]), t(P))
#' matplot(a1$Q[,ik], Q); abline(0,1, col=K+1, lwd=2, lty=2)
#' Y = apply(Q, 1, which.max)
#' ml = nnet::multinom(Y~a1$Q[,-1], maxit=1e3)
#' Yc = apply(ml$fitted, 1, which.max)
#' table(Y, Yc); mean(Yc!=Y)
#' @export
Radmixture <- function(G, K=2, pinit=NULL, Kem=5, maxit=1e2,rtol=1e-7,trace=FALSE, ...) {
    n = dim(G)[1]; m = dim(G)[2]
    if(is.null(pinit)) {
        cl = kmeans(scale(G), centers=K, ...)$cluster
        ## cl = sample(1:K, size=n, rep=TRUE)
        Q.ki = matrix(table(cl)/length(cl), K,n)
        P.kj = matrix(0, K,m)
        for(k in 1:K) {
            P.kj[k,] = colMeans(G[cl==k,,drop=FALSE], na.rm=TRUE)/2
        }
    } else {
        Q.ki = t(pinit$Q); P.kj = pinit$P
    }

    obj = .Call("Radmixture",
                as.double(Q.ki), as.double(P.kj), as.double(G), double(maxit+1),
                as.integer(n), as.integer(m), as.integer(K), as.integer(Kem),
                as.integer(maxit), rtol, as.integer(trace))
    Q.ki = matrix(obj[[1]], K,n); Cs = apply(Q.ki, 2, which.max)
    P.kj = matrix(obj[[2]], K,m);
    llk = obj[[3]]
    itMax = obj[[4]]
    return( list(Q=t(Q.ki), P=P.kj, Cs=Cs, llk=llk[1:itMax]) )
}






## .C interface ## a bit slower
Radmixture0 <- function(G, K=2, pinit=NULL, Kem=5, maxit=1e2,rtol=1e-7,trace=FALSE, ...) {
    n = dim(G)[1]; m = dim(G)[2]
    if(is.null(pinit)) {
        cl = kmeans(scale(G), centers=K, ...)$cluster
        Q.ki = matrix(table(cl)/length(cl), K,n)
        P.kj = matrix(0, K,m)
        for(k in 1:K) {
            P.kj[k,] = colMeans(G[cl==k,,drop=FALSE], na.rm=TRUE)/2
        }
    } else {
        Q.ki = t(pinit$Q); P.kj = pinit$P
    }

    obj = .C("radmixture_qp",
             Q.ki = as.double(Q.ki), P.kj = as.double(P.kj), llk = double(maxit),
             as.integer(n), as.integer(m), as.double(G), as.integer(K), as.integer(Kem),
             itMax=as.integer(maxit), as.double(rtol), as.integer(trace), PACKAGE='Radmixture')
    Q.ki = matrix(obj$Q.ki, K,n); Cs = apply(Q.ki, 2, which.max)
    P.kj = matrix(obj$P.kj, K,m); ## llk = obj$llk
    return( list(Q=t(Q.ki), P=P.kj, Cs=Cs, llk=obj$llk[1:min(obj$itMax,maxit)]) )
}

