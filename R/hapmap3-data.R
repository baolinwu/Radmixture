#' Hapmap3 data used for illustration of ancestry inference based on the Admixture model
#'
#' Data are constructed from the Phase 3 of the HapMap Project. Selected individuals include:
#' (1) MEX: Mexican ancestry sampled in Los Angeles (n=50),
#' (2) ASW: African ancestry  sampled in the American Southwest (n=49),
#' (3) CEU: Utah residents with ancestry from Northern and Western Europe (n=112); and
#' (4) YRI: Yoruba in Ibadan, Nigeria (n=113).
#' A subset of the 1,440,616 available markers were chosen according to two criteria: (1) to minimize background linkage disequilibrium, adjacent markers must be no closer than 200 kb apart,
#' and (2) no more than 5\% of the genotypes must be missing.
#' Based on the genotypes for these markers for the unrelated individuals from the CEU, YRI, MEX, and ASW samples, a data set of 13,298 markers typed on 324 individuals are selected.
#'
#' @format A data list hapmap3 with two components
#' \describe{
#'   \item{Y}{ population indicators: 1,2,3,4 for (CEU, ASW, MEX, YRI)  }
#'   \item{G}{ 324x13928 matrix of genotypes }
#' }
#'
#' @references
#' Wu X, and Wu B (2024) Radmixture: an R package for fast and robust likelihood-based estimation of admixture model for global ancestry inference. tech rep.
#'
#' Ko S, Chu BB, Peterson D, Okenwa C, Papp JC, Alexander DH, et al. (2023) Unsupervised discovery of ancestry-informative markers and genetic admixture proportions in biobank-scale datasets. Am J Hum Genet. 2023 Feb 2;110(2):314–25. 
#' 
#' D.H. Alexander, J. Novembre, and K. Lange (2009) Fast model-based estimation of ancestry in unrelated individuals. Genome Research, 19:1655-1664.
#'
#' H. Zhou, D.H. Alexander, and K. Lange (2011) A quasi-Newton method for accelerating the convergence of iterative optimization algorithms. Statistics and Computing, 21:261-273.
#'
#' @source \url{https://dalexander.github.io/admixture/download.html}
#' 
#' @examples
#' data(hapmap3); G = hapmap3$G; Y = hapmap3$Y
#' jj = which(apply(G, 2, sd)>0); G = G[,jj]
#' h1 = Radmixture(G,K=4, Kem=5,maxit=1e2,rtol=1e-7,trace=TRUE, nstart=5)
#' apply(h1$Q, 2, range)
#' table(h1$Cs, Y)
#' ## pdf(file='hapmap3-Q.pdf', width=8, height=6)
#' barplot(t(h1$Q), col=rainbow(4), xlab="Individuals", ylab="Ancestry proportions (Q)", border=NA, ylim=c(-0.08,1))
#' yn = c(0, cumsum(table(Y)[c(2,1,3,4)]))
#' for(i in 1:4) lines(c(yn[i]+1,yn[i+1])*1.2, c(-0.02,-0.02)+i*0.002, lwd=2, col=1)
#' text(1.2*(yn[-5]/2+yn[-1]/2),-0.04+1:4*0.002, c('ASW', 'CEU', 'MEX', 'YRI'), col=1)
#' pairs(h1$Q, col=Y, pch=Y)
#' ml = nnet::multinom(Y~h1$Q[,-1])
#' logLik(ml)
#' cl = apply(ml$fitted, 1, which.max)
#' table(cl, Y)
"hapmap3"


