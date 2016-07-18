#checking rlImportSamp with dissertation results
vinflate<-seq(from=nForPrior, to=25, by=2)
nForPrior<-3
data(phones)
phones<-data.frame(cbind(phones$year, phones$calls))
names(phones)<-c("year", "calls")
summary(phones)
phones$year<-phones$year-mean(phones$year) #center
fit3<-rlm(calls~year,data=phones[-c(1:nForPrior),], maxit=100, psi=psi.bisquare)
#abline(fit3, col=1)
summary(fit3)
fit3$s^2

phonesFit<-phones[-c(1:nForPrior),]
n<-nrow(phonesFit)
X<-cbind(rep(1, n), phonesFit$year)
y<-phonesFit$calls
phonesPrior<-phones[c(1:nForPrior),]




#priors
priorFit<-lm(calls~year, data=phonesPrior)
solve(t(model.matrix(priorFit))%*%model.matrix(priorFit))*summary(priorFit)$sigma^2
summary(priorFit)

mu0<-coef(priorFit)
mu0
sigma2Hat<-summary(priorFit)$sigma^2
alpha<-2
beta<-(alpha-1) #sigma2Hat*(alpha-1)
beta
beta/(alpha-1)
#curve(dinvgamma(x, alpha,beta), 0,10)
# integrate(dinvgamma, 0,s.obs^2,shape=alpha,scale=beta)
# beta^2/((alpha-1)^2*(alpha-2))
library(mvtnorm)
inflate<-2
Sigma0<-vinflate[inflate]*vcov(priorFit)

x=X;y=y; psi=psi.bisquare; scale.est='MAD'; mu0=mu0; Sigma0=Sigma0; alpha=alpha; beta=beta; smooth=1;N=1000;Nins=100; maxit=1000; k2=1.345



tst<-rlImportSamp(x=X,y=y, psi=psi.bisquare, scale.est='Huber', mu0=mu0, Sigma0=Sigma0, alpha=alpha, beta=beta, smooth=1,N=10000,Nins=10000, maxit=1000)
tst
plot(tst$w, cex=.1, pch=19)
