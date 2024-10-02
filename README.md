# Radmixture

An R package for population structure inference based on the admixture model. We provide very efficient and robust computation
of the maximum likelihood estimates. 

-----
## Sample codes for running the analysis
```
n=1000; K=4; m1 = 3000; m0 = 5000; m=m0+m1
P1 = t(sapply(1:K, function(j) runif(m1)))
P0 = matrix(runif(m0), K,m0, byrow=TRUE)
P = cbind(P1,P0)
Q = matrix(runif(n*K), n,K); Q = Q/rowSums(Q)
tha = Q%*%P; G = matrix(rbinom(n*m, 2, tha), n,m)
jj = which(apply(G, 2, sd)>0)
G = G[,jj]; P = P[,jj]
a1 = Radmixture(G,K=K, Kem=10,maxit=1e3,rtol=1e-7,trace=TRUE)
apply(Q, 2, range); apply(a1$Q, 2, range); cor(Q, a1$Q)
table(a1$Cs, apply(Q, 1, which.max))
par(mfrow=c(1,2))
barplot(t(a1$Q), col=rainbow(K), xlab="Individual #", ylab="Ancestry", border=NA)
ik = apply(abs(cor(Q, a1$Q)), 1, which.max); ik
cor(t(a1$P[ik,]), t(P))
matplot(a1$Q[,ik], Q); abline(0,1, col=K+1, lwd=2, lty=2)
Y = apply(Q, 1, which.max)
ml = nnet::multinom(Y~a1$Q[,-1], maxit=1e3)
Yc = apply(ml$fitted, 1, which.max)
table(Y, Yc); mean(Yc!=Y)
```
