set.seed(1)
y<-rnorm(100)



tst<-fitOrderStat(y,k=4, theta.lims=c(-4,4),length.theta=500, sigma2.lims=c(.0001,10),length.sigma2=500, mu=0, tau=100, alpha=2, beta=1)
tst$jointPost
tst$thetaPost
plot(tst$thetaPost[,1],tst$thetaPost[,2], type='l')
plot(tst$sigma2Post[,1],tst$sigma2Post[,2], type='l')
head(tst$thetaPost)
tail(tst$thetaPost)
head(tst$sigma2Post)
tail(tst$sigma2Post)
sum(tst$thetaPost[,2])*diff(tst$thetaPost[,1])[1]
sum(tst$sigma2Post[,2])*diff(tst$sigma2Post[,1])[1]

tst$postMeans

set.seed(1)
y<-rnorm(10)
tst<-fitOrderStat(y,k=4, theta.lims=c(-5,5),length.theta=500, sigma2.lims=c(.0001,10),length.sigma2=500, mu=-1, tau=.1, alpha=1, beta=2)


tst$jointPost
tst$thetaPost
plot(tst$thetaPost[,1],tst$thetaPost[,2], type='l')
abline(v=mean(y))
abline(v=-1)
plot(tst$sigma2Post[,1],tst$sigma2Post[,2], type='l')
head(tst$thetaPost)
tail(tst$thetaPost)
head(tst$sigma2Post)
tail(tst$sigma2Post)
sum(tst$thetaPost[,2])*diff(tst$thetaPost[,1])[1]
sum(tst$sigma2Post[,2])*diff(tst$sigma2Post[,1])[1]
tst$postMeans


#fit Mixture:
n<-100
c<-10
theta_true<-0
sigma2_true<-1
p<-.1
c<-100
x<-rbinom(n,1,p=p)

y<-rnorm(n, mean=theta_true, sqrt(sigma2_true*ifelse(x, c, 1)))
hist(y)

tst<-fitMixture(y, theta.lims=c(-2,2), length.theta=500, sigma2.lims=c(.01,5), length.sigma2=500, mu=0, tau=100, alpha=1, beta=2, p=p, c=c)


tst$jointPost
tst$thetaPost
plot(tst$thetaPost[,1],tst$thetaPost[,2], type='l')
abline(v=mean(y))
abline(v=theta_true)
plot(tst$sigma2Post[,1],tst$sigma2Post[,2], type='l')
abline(v=var(y))
abline(v=sigma2_true)
head(tst$thetaPost)
tail(tst$thetaPost)
head(tst$sigma2Post)
tail(tst$sigma2Post)
sum(tst$thetaPost[,2])*diff(tst$thetaPost[,1])[1]
sum(tst$sigma2Post[,2])*diff(tst$sigma2Post[,1])[1]
tst$postMeans


###


#fits used in my dissertation
load("~/Dropbox/school/osu/dissertationResearch/snm/locationAndScale/dataAnalysis/newcomb/directSampling/wsNewcombDirectSamplingHuber2.RData")
set.seed(1)
fit<-rlDirectEval(y=newcomb, psi=psi.huber, scale.est='Huber', eta=23.6, tau=2.04, alpha=5, beta=10, mu_lims=c(20,32), sigma2_lims=c(0.001,100), length_mu=500, length_sigma2=500,smooth=1,N=1e5)
fit$H
set.seed(1)
fit2<-rlDirectEval(y=newcomb, psi=psi.huber, scale.est='Huber',eta=23.6, tau=2.04, alpha=5, beta=10, mu_lims=c(20,32), sigma2_lims=c(0.001,100), length_mu=100, length_sigma2=100,smooth=1,N=1e3,maxit=1000)


plot(fit$muPost[,1],fit$muPost[,2], type='l', col=4)
lines(theta.grid, post.theta.pdf, col=2)
lines(fit2$muPost[,1],fit2$muPost[,2])
plot(fit$sigma2Post[,1],fit$sigma2Post[,2], type='l')
lines(sigma2.grid, post.sigma2.pdf, col=2)
lines(fit2$sigma2Post[,1],fit2$sigma2Post[,2])
