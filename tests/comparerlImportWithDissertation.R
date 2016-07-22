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

fit<-rlm(calls~year,data=phonesFit, maxit=100, psi=psi.bisquare)
summary(fit)
fit2<-rlm(X,y, maxit=100, psi=psi.bisquare)
summary(fit2)

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
source('~/Dropbox/school/osu/dissertationResearch/snm/rPackage/brlm/R/rlImportSamp.R')
inflate<-10
Sigma0<-vinflate[inflate]*vcov(priorFit)

y=y; psi=psi.bisquare; scale.est='MAD'; mu0=mu0; Sigma0=Sigma0; alpha=alpha; beta=beta; smooth=1;N=100;Nins=100; maxit=1000; k2=1.345



tst<-rlImportSamp(X=X,y=y, psi=psi.bisquare, scale.est=scale.est, mu0=mu0, Sigma0=Sigma0, alpha=alpha, beta=beta, smooth=1,N=1e4,Nins=1e4, maxit=10000)

plot(tst$w, cex=.1, pch=19)
plot(density(tst$impSamps[,1],weights=tst$w))
abline(v=b.obs[1], col=2)
plot(density(tst$impSamps[,2],weights=tst$w))
abline(v=b.obs[2], col=2)
plot(density(tst$impSamps[,3],weights=tst$w))
abline(v=s.obs, col=2)
abline(h=11)
range(tst$impSamps[,3])
weighted.mean(tst$impSamps[,1],tst$w)
weighted.mean(tst$impSamps[,2],tst$w)
b.obs
weighted.mean(tst$impSamps[,3],tst$w)
s.obs^2
