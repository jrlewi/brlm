#[beta| mu0, Sigma0, a]~N(mu0, a*Sigma0) #mu0, Sigma0- fixed
#[b_i|beta, Sigma0, b]~iid N(beta,b^*Sigma0)

#b=b^*swSq, where swSq is fixed (swSq=sum w_i^2, from previous data set)

#Variance is partitioned between beta and (b_i, i=1,..,K) to connect to the non-hierarchical version of this model
# this reasoning implies a+b^*=1

#[b^*|v1, v2]~beta(v1, v2)

#priors are placed on mu_b*=v1/(v1+v2) and psi_b*=v1+v2 in the following way
#mu_b*~beta(alpha_mu*, beta_mu*) with parameters fixed so the mean is small....see notes and setting of prior parameters
#psi_b*~gamma(a_psib, b_psi_b)
#that is, the mean of b^*|v1,v2 is mu_b* and the variance is mu_b*(1-mu_b*)/(psi_b*+1)


#[sigma^2_i |mu_rho, K_rho, rho]~InG(a0,b0) a0, b0 fixed from previous data set (use the same values as the non-hiearchical fit). But sigma^2_i's not independent

# connect the sigma^2_i's via the formulation
# Z~N_K(0, Sigma_rho) with sigma^2_i=H^-1(phi(z_i)) where H is the cdf of IG(a0, b0) and Sigma_rho=rho*J+(1-rho)I

#to allow for data dependence in the correlation of the sigma^2_i

#[rho| mu_rho, psi_rho]~beta(alpha_rho, beta_rho) where alpha_rho and beta_rho are functions of mu_rho, psi_rho, so that the mean of the beta for rho is mu_rho and the variance is mu_rho(1-mu_rho)/(psi_rho+1)

#[mu_rho] ~beta(w1,w2), [psi_rho]~Gamma(a_psir, b_psir), fixed parameters


#finally, for the K groups i=1,...,K
# y_i=b_i*x_ij+e_ij j=1,..., n_i (n_i is the number of obs in the ith group)

# y_ij|b_i,sigma_i^2~N(b_i*x_ij,sigma_i^2)

# e_ij~(iid) N(0, sigma_i^2)

# mu_0 and Sigma0 are fixed a priori, so are a0 and b0
# NOTE: y and X are lists of the individual responses and design matrices for each group
#Model: obtain samples from [mu_b*, psi_b*, b*,beta, beta_i, i=1,...,K, Z,  mu_rho, psi_rho,rho|Y]
#strategy: write sampling functions for each full conditional, then combine into one function for one repetition
#fixed HyperParameters are
#mu0, Sigma0,Sigma0Inv, a0,b0,alpha_mustr, beta_mustr,a_psib, b_psib
#w1, w2,a_psir, b_psir
# other fixed values are swSq (sum of squared weights), p, nGroups, J_nGroups (nGroups by nGroups matrix of ones), I_nGroups (nGroups by nGroups Identity)


#fix mu_bstr and psi_bstr

#' sampling mu in hierarchical model of paper.
#' [mu_bstr|-]
#' function to compute shape parameters of beta from the mean mu<-a/(a+b) and psi=a+b (precision)
#' @rdname hier_samp_funs
#' @param mu = a/(a+b)
#' @param psi = a + b
#' @details a, b are the shape parameters of a beta distribution
fn.compute.ab<-function(mu, psi){
  a<-mu*psi
  b<-psi-a
  c(a,b)
}

#' quad1 and qaud2
#' @name quads
#' @details Used for sampling Beta in hierarchical model of paper. quad1 is the current value of (beta-mu0)'Sigma0Inv(beta-mu0). quad2 is the current value of sum_i^nGroups (beta_i-beta) Sigma0Inv (beta_i-beta)
#' @param Beta Beta
#' @param mu0 mu0
#' @param betalMat matrix of betal's
fn.compute.quad1<-function(Beta,mu0, Sigma0Inv){
  t(Beta-mu0)%*%Sigma0Inv%*%(Beta-mu0)
}
#' @rdname quads
fn.compute.quad2<-function(Beta, betalMat,Sigma0Inv){
  diff<-betalMat-Beta
  sum(apply(diff,MARGIN=2,FUN=function(x) t(x)%*%Sigma0Inv%*%x))
}


#' functions for bstar step in hierarchical model
#' @name bstar_step
#' @details functions for bstar step
#' @param bstar bstar
#' @param v1 v1
#' @param v2 v2
#' @param quad1 quad1
#' @param quad2 quad2
#' @param K K
#' @param p p
#' @param logbstarCur description
#' @param step_logbstar description
#' @param bstarProp description
#' @param bstarCur description
fn.loglike.bstar<-function(bstar,v1,v2, quad1, quad2, K,p)
{
  if(bstar>0 && bstar<1){
    #ab<-fn.compute.ab(mu_bstr, psi_bstr)
    #v1<-ab[1]; v2<-ab[2]
    (v1-.5*K*p-1)*log(bstar)+(v2-.5*p-1)*log(1-bstar)-.5*quad1/(1-bstar)-.5*swSq*quad2/bstar #-log(beta(v1,v2))
  } else {
    -Inf}
}
#' @rdname bstar_step
fn.proposal.logbstar<-function(logbstarCur, step_logbstar){
  rnorm(1, logbstarCur,step_logbstar)# sd=step)
}
#' @rdname bstar_step
fn.logpdf.prop.logbstar<-function(bstarProp){
  -log(bstarProp)
} #dnorm(log(bstarProp),log(bstarCur),step_logbstar, log=TRUE)-log(bstarProp)
#' @rdname bstar_step
fn.sample.bstar<-function(bstarCur,v1,v2,quad1,quad2, K,p, step_logbstar) {
  logbstarCur<-log(bstarCur)
  logbstarProp<-fn.proposal.logbstar(logbstarCur, step_logbstar)
  bstarProp<-exp(logbstarProp)
  logRat1<-fn.loglike.bstar(bstarProp, v1,v2, quad1, quad2, K,p)-fn.loglike.bstar(bstarCur, v1,v2, quad1, quad2,K,p)
  logRat2<-fn.logpdf.prop.logbstar(bstarCur)-fn.logpdf.prop.logbstar(bstarProp)
  mhRat6<-min(exp(logRat1+logRat2),1)
  if(runif(1)<mhRat6){
    return(bstarProp)
  } else {
    return(bstarCur)
  }
}



#' sampling for beta in hieararchical model
#' [Beta|-]
#'@param betalMat, matrix with columns equal to current beta_l's from each group
#' @param bstar: current value of bstar
#' @details Sigma0 is fixed,  nGroups (number of groups) and swSq are fixed, and used within this function.

fn.sample.Beta<-function(betalMat,bstar){
  a<-1-bstar
  b<-bstar/swSq
  prec<-(1/a+nGroups/b)
  SigmaB<-(Sigma0)/prec
  muB<-(rowSums(betalMat)/b+mu0/a)/prec
  mvrnorm(1, muB, SigmaB)
}


#' Sampling for beta_l's in hierarchical model
#' 5. [beta_l|-] l=1,..,nGroups
#' @name beta_ls
#' @param yl,Xl these are the responses and design matrix for the lth group
#' @param XtX.l X'X for lth group
#' @param XtX list of X'X for the l groups
#' @param y,X  list of yl and Xl
#' @param sigma2l: the sigma^2_l for the l group
#' @param bstar b=b^*/swSq
#' @param Beta Beta
#' @details Sigma0Inv: fixed  loop through for each group and do sampling for each group beta
fn.sample.beta.l<-function(yl, Xl,XtX.l,Beta,sigma2l,bstar, Sigma0Inv){
  b<-bstar/swSq
  Sigma.betal<-solve(XtX.l/(sigma2l)+Sigma0Inv/b)
  mu.l<-Sigma.betal%*%(crossprod(Xl,yl)/(sigma2l)+Sigma0Inv%*%Beta/b)
  mvrnorm(1, mu.l, Sigma.betal)
}
#' @rdname beta_ls
fn.one.rep.beta.l<-function(y, X,XtX,Beta,sigma2,bstar,betalMat, Sigma0Inv){
  for(l in 1:ncol(betalMat)){
    yl<-y[[l]]
    Xl<-X[[l]]
    XtX.l<-XtX[[l]]
    sigma2l<-sigma2[l]
    beta.l<-fn.sample.beta.l(yl, Xl,XtX.l,Beta,sigma2l,bstar, Sigma0Inv)
    betalMat[,l]<-beta.l
  }
  betalMat
}

#' Sampling mu_rho
#' [mu_rho|-]
#' @name mu_rho
#' @param mu_rho, mu_rhoCur
#' @param psi_rho descrip
#' @param rho descrip
#' @param mu_rho_step descrip
#' @details sampling mu_rho in hierarchical model
fn.loglik.mu_rho<-function(mu_rho,psi_rho, rho){
  #  if(mu_rho>0.001 && mu_rho<.999){
  #  if((mu_rho>0.001 && mu_rho<.999) && psi_rho>0){
  if(mu_rho>0 && mu_rho<1){
    ab<-fn.compute.ab(mu_rho,psi_rho)
    alpha_rho<-ab[1]
    beta_rho<-ab[2]
    lik<-(w1-1)*log(mu_rho)+(w2-1)*log(1-mu_rho)+(alpha_rho-1)*log(rho)+(beta_rho-1)*log(1-rho)-log(beta(alpha_rho,beta_rho))
  } else {
    lik<--Inf
  }
  return(lik)
}
#' @rdname  mu_rho
fn.proposal.mu_rho<-function(mu_rhoCur, mu_rho_step){
  rnorm(1,mu_rhoCur, mu_rho_step)
}
#' @rdname  mu_rho
fn.sample.mu_rho<-function(mu_rhoCur,psi_rho, rho, mu_rho_step){
  mu_rhoProp<-fn.proposal.mu_rho(mu_rhoCur, mu_rho_step)
  #mu_rhoProp<-fn.proposal.mu_rho()
  logRat1<-fn.loglik.mu_rho(mu_rhoProp,psi_rho, rho)-fn.loglik.mu_rho(mu_rhoCur,psi_rho, rho)
  logRat2<-0 #fn.logpdf.prop.mu_rho(mu_rhoCur,mu_rhoProp,mu_rho_step)-fn.logpdf.prop.mu_rho(mu_rhoProp, mu_rhoCur,mu_rho_step)
  mhRat7<-min(exp(logRat1+logRat2),1)
  if(runif(1)<mhRat7){
    return(mu_rhoProp)
  } else {
    return(mu_rhoCur)
  }
}
#' Sampling psi_rho: [psi_rho|-]
#' @name psi_rho_samp
#' @param psi_rho,psi_rhoCur desc
#' @param mu_rho desc
#' @param rho desc
#' @param psi_rho_step desc
fn.loglik.psi_rho<-function(psi_rho,mu_rho, rho){
  #if((mu_rho>0 && mu_rho<1) && psi_rho>0){(mu_rho>0.001 && mu_rho<.999) &&
  if(psi_rho>0){
    ab<-fn.compute.ab(mu_rho,psi_rho)
    alpha_rho<-ab[1]
    beta_rho<-ab[2]
    lik<-(a_psir-1)*log(psi_rho)-b_psir*psi_rho+(alpha_rho-1)*log(rho)+(beta_rho-1)*log(1-rho)-log(beta(alpha_rho,beta_rho))
  } else {
    lik<--Inf
  }
  return(lik)
}
#'@rdname psi_rho_samp
fn.proposal.psi_rho<-function(psi_rhoCur, psi_rho_step){ #a_psir, b_psir
  rnorm(1,psi_rhoCur, psi_rho_step)
}
#'@rdname psi_rho_samp
fn.logpdf.prop.psirho<-function(psi_rho){
  (a_psir-1)*log(psi_rho)-b_psir*psi_rho
} #
#'@rdname psi_rho_samp
fn.sample.psi_rho<-function(psi_rhoCur,mu_rho, rho, psi_rho_step){
  #psi_rhoProp<-fn.proposal.psi_rho()
  psi_rhoProp<-fn.proposal.psi_rho(psi_rhoCur, psi_rho_step)
  #print(psi_rhoProp)
  #print(c(mu_rho,rho))
  logRat1<-fn.loglik.psi_rho(psi_rhoProp,mu_rho, rho)-fn.loglik.psi_rho(psi_rhoCur,mu_rho, rho)
  logRat2<-0 #fn.logpdf.prop.psirho(psi_rhoCur)-fn.logpdf.prop.psirho(psi_rhoProp)
  mhRat1<-min(exp(logRat1+logRat2),1)
  if(runif(1)<mhRat1){
    return(psi_rhoProp)
  } else {
    return(psi_rhoCur)
  }
}

#' sampling for [rho|-] in hierchical model
#' @name rho_samp
#' @param rho desc
#' @param K  desc
#' @param rhoCur desc
#' @param rho_step desc
#' @param rho desc
#' @param Z desc
#' @param alpha_rho desc
#' @param beta_rho desc
#' @param mu_rho desc
#' @param psi_rho desc
#' @param Zdesc desc
#' @details J_nGroups must be defined, I_nGroups must be defined
#' function to effeciently caclulate Sigma_rhoInv using its special form and a well known matrix inverse result (from the fix-rank kriging chapter 8 of the spatial handbook book); note; the K variable added for the sampling of the z_i's later.
#' fn.compute.logDetSigmaRho used to effeciently caclulate log(det(Sigma_rho)) using its special form and Sylvester's determinant theorem: det(I_n+AB)=det(I_m+BA)
fn.compute.SigmaRhoInv<-function(rho, K){
  #nGroups<-nrow(I_nGroups)
  diag(K)/(1-rho)-matrix(1, K,K)/((1-rho)^2*(K/(1-rho)+1/rho))
}
#' @rdname rho_samp
fn.compute.logDetSigmaRho<-function(rho, K){
  if(rho>0 & rho <1){
    K*log(1-rho)+log(1+K*rho/(1-rho))} else {
      if(rho==0) {log(1)}
      else {if(rho==1){log(0)
      }
      }
    }
}
#' @rdname rho_samp
fn.proposal.rho<-function(rhoCur, rho_step){
  rnorm(1, rhoCur, rho_step)
}
#' @rdname rho_samp
fn.loglike.rho<-function(rho, Z,alpha_rho,beta_rho){
  if(rho<1 && rho>0){
    #if(rho<0.999 && rho>0.001){
    #Sigma_rho<-rho*J_nGroups+(1-rho)*I_nGroups #do i need this?
    K<-length(Z)
    logdetSigma_rho<-fn.compute.logDetSigmaRho(rho,K)
    Sigma_rhoInv<-fn.compute.SigmaRhoInv(rho,K)
    lik<-(alpha_rho-1)*log(rho)+(beta_rho-1)*log(1-rho)-.5*logdetSigma_rho-.5*t(Z)%*%Sigma_rhoInv%*%Z
  } else {
    lik<--Inf
  }
  lik
}
#' @rdname rho_samp
fn.sample.rho<-function(rhoCur, mu_rho, psi_rho, Z, rho_step){
  ab<-fn.compute.ab(mu_rho,psi_rho)
  alpha_rho<-ab[1]
  beta_rho<-ab[2]
  #rhoProp<-fn.proposal.rho()
  #rhoProp<-fn.proposal.rho(alpha_rho, beta_rho)
  rhoProp<-fn.proposal.rho(rhoCur, rho_step)
  # print(c(rhoCur, rhoProp))
  logRat1<-fn.loglike.rho(rhoProp, Z,alpha_rho,beta_rho)-fn.loglike.rho(rhoCur, Z,alpha_rho,beta_rho)
  logRat2<-0 #fn.logpdf.prop.rho(rhoCur,rhoProp,rho_step)-fn.logpdf.prop.rho(rhoProp,rhoCur,rho_step)
  mhRat2<-min(exp(logRat1+logRat2),1)
  if(runif(1)<mhRat2){
    return(rhoProp)
  } else {
    return(rhoCur)
  }
}


# Sampling Zl (i.e. sigma2 with dependence) in hierarchical model
#' [Zl|-]: sample one at a time
#' @name z_l
#' @param Z,sigma2 desc
#' @param a0,b0 desc
#' @param Xl,bl desc
#' @param zl,zlCur,zlProp desc
#' @param step desc
#' @param zminus, desc
#' @param yl, desc
#' @param fitsl, desc
#' @param condsd, desc
#' @param condMeanMult desc
#' @details function to compute simga2_is from a z vector (or a single z) from nGroups dimensional Normal distribution with marginal means 0 and vars 1
# a0,b0 fixed
fn.compute.sigma2<-function(Z,a0,b0){
  X<-qgamma(pnorm(Z),shape=a0, scale=1/b0) #rate=b0
  1/X
}
#' @rdname z_l
fn.compute.Z2<-function(sigma2,a0,b0){
  X<- 1-pgamma(1/sigma2, a0, 1/b0)
  qnorm(X)
}
#' @rdname z_l
fn.compute.Z<-function(sigma2,a0,b0){
  X<-integrate(function(x) {dinvgamma(x, shape=a0, scale = b0)},0, sigma2)$value #rate=b0
  qnorm(X)

}
#' @rdname z_l
fn.hat<-function(Xl, betal) {Xl%*%betal}
#' @rdname z_l
fn.logLike.zl<-function(zl,zminus, yl, fitsl, condsd,condMeanMult){
  #zl is the component I am updating, zminus is all other components; condsd is the condition sd for the conditional normal (zl|zminus): this depends only on current
  #condMeanMult%*%zminus is the conditional mean for zl|zminus
  condMean<-condMeanMult%*%zminus
  sigma2l<-fn.compute.sigma2(zl,a0,b0)
  dnorm(zl,condMean,condsd,log=TRUE)+sum(dnorm(yl, fitsl, sqrt(sigma2l),log=TRUE))
}
#' @rdname z_l
fn.proposal.zl<-function(zlCur, step=1){
  rnorm(1, zlCur,step)
}
#' @rdname z_l
fn.logpdf.proposal.Z<-function(zlProp,zlCur,step){
  dnorm(zlProp, zlCur,step,log=TRUE)
}
#' @rdname z_l
fn.sample.Zl<-function(zlCur,zminus,yl, fitsl, condsd,condMeanMult, step=1){
  zlProp<-fn.proposal.zl(zlCur, step)
  logRat1<-fn.logLike.zl(zlProp,zminus, yl, fitsl, condsd,condMeanMult)-fn.logLike.zl(zlCur,zminus, yl, fitsl, condsd,condMeanMult)
  logRat2<-0 #random walk: symmetric
  mhRat3<-min(exp(logRat1+logRat2),1)
  if(runif(1)<mhRat3){
    return(zlProp)
  } else {
    return(zlCur)
  }
}
#' @rdname z_l
fn.one.rep.Z<-function(Zcur, y,X, betalMat, rho,step=rep(1,length(Zcur))){
  Sigma_rhoInvMinus<-fn.compute.SigmaRhoInv(rho, K=length(Zcur)-1)    #inverse (1-pho)I+rhoJ for K-1 by K-1 dim matrices
  rho_vec<-rep(rho, length(Zcur)-1)
  condsd<-sqrt(1-t(rho_vec)%*%Sigma_rhoInvMinus%*%rho_vec)
  condMeanMult<-t(rho_vec)%*%Sigma_rhoInvMinus
  for(j in 1:length(Zcur)){
    zlCur<-Zcur[j]
    zminus<-Zcur[-j]
    yl<-y[[j]]
    Xl<-X[[j]]
    betal<-betalMat[,j]
    fitsl<-fn.hat(Xl, betal)
    zlCur<-fn.sample.Zl(zlCur,zminus,yl, fitsl, condsd,condMeanMult,step[j])
    Zcur[j]<-zlCur
  }
  Zcur
}


#----- Fitting functions for Normal Theory Model -------

#' MCMC functions for the hierarchical versions of Normal Theory Model, corresponding restricted version, and heavy-tailed version used in the paper.
#' @name hier_models
#' @param y,X lists of group level responses and design matrices
#' @param XtX list of X'X - for all the groups (input for efficiency)
#' @param v1,v2 desc
#' @param bstar desc
#' @param Beta desc
#' @param betalMat desc
#' @param Z desc
#' @param mu_rho desc
#' @param psi_rho desc
#' @param rho desc
#' @details fn.hier.one.rep is one rep for the full normal thoery model. heirNormTheoryLm uses fn.hier.one.rep for the complete MCMC of full normal thoert model. hierNormTheoryRestLm does the same for the restricted versions. The corresponding t-model versions are handedl by fn.one.rep.tHierModel and hier_TLm.
#'
#' Details of the model are describted in the paper.
#' @param sigma2Int is the vector of sigma2i initial values
#' @param mu0: prior mean of each beta
#' @param Sigma0 is the 'variance' matrix of beta b_i~N(mu0,b^*Sigma0)
#' @param a0,b0 prior parameters for sigma2
#' @param v1,v2 parameters for beta on b_start. mu_bstr is the mean and psi_bstr (originally had a prior) is the precision of the beta prior for b^*
#in detail b*~beta(v1,v2), mu_bstr=v1/(v1+v2), psi_bstr=v1+v2
#' @param mu_bstr,psi_bstr
#' @param swSq default to 1.
#' @param w1,w2,a_psir,b_psir parameters definining prior for rho. In detail: rho~beta(mean=mu_rho, precision=psi_rho), mu_rho~beta(w1,w2) and psi_rho~gamma(a_psir, b_psir)
#' @param nkeep,nburn number of MCMC iterations to keep, number for burn in.
#' @param step_logbstar,mu_rho_step,psi_rho_step,rho_step,step_Z tunning parameters for MH Steps

fn.hier.one.rep<-function(y,
                          X,
                          XtX,
                          v1,#mu_bstr,
                          v2,#psi_bstr,
                          bstar,
                          Beta,
                          betalMat, #first to update
                          Z,
                          mu_rho,
                          psi_rho,
                          rho,
                          step_logbstar,
                          mu_rho_step,
                          psi_rho_step,
                          rho_step,
                          step_Z,
                          Sigma0Inv
){ #y, X are list of group level responses
  #fn.one.rep.beta.l loops through for each beta.l
  #[beta_i|-]
  sigma2<-fn.compute.sigma2(Z,a0,b0)
  betalMat<-fn.one.rep.beta.l(y, X,XtX,Beta,sigma2,bstar,betalMat, Sigma0Inv)
  #[Z|-] i.e. the sigma2_is
  #fn.one.rep.Z loops through for all the z_i
  Z<-fn.one.rep.Z(Z, y,X, betalMat, rho,step=step_Z)
  sigma2<-fn.compute.sigma2(Z,a0,b0)
  #[Beta|-]
  Beta<-fn.sample.Beta(betalMat,bstar)

  # [mu_bstr|-]
  # mu_bstr<-fn.sample.mubstr(mu_bstr,psi_bstr, bstar,step_mubstr)

  #[psi_bstr|-]
  #  psi_bstr<-fn.sample.psibstr(psi_bstr,mu_bstr, bstar)

  #[bstar|-]
  quad1<-fn.compute.quad1(Beta, mu0, Sigma0Inv)
  quad2<-fn.compute.quad2(Beta, betalMat, Sigma0Inv)
  bstar<-fn.sample.bstar(bstar,v1,v2,quad1,quad2,K=nGroups,p, step_logbstar)
  #[mu_rho|-]
  mu_rho<-fn.sample.mu_rho(mu_rho,psi_rho, rho, mu_rho_step)

  #[psi_rho|-]
  psi_rho<-fn.sample.psi_rho(psi_rho,mu_rho, rho, psi_rho_step)

  #[rho|-]
  rho<-fn.sample.rho(rho, mu_rho, psi_rho, Z, rho_step)

  out<-list()
  out$Beta<-Beta
  out$betalMat<-betalMat
  out$Z<-Z
  out$sigma2<-sigma2
  #out$mu_bstr<-mu_bstr
  #out$psi_bstr<-psi_bstr
  out$bstar<-bstar
  out$mu_rho<-mu_rho
  out$psi_rho<-psi_rho
  out$rho<-rho
  out
}
#' @rdname hier_models
#' @export
hierNormTheoryLm<-function(y,
                           X,
                           nkeep=1e4, nburn=1e3,
                           mu0,
                           Sigma0,
                           a0,
                           b0,
                           mu_bstr,
                           psi_bstr,
                           swSq=1,
                           w1,
                           w2,
                           a_psir,
                           b_psir,
                           step_logbstar,
                           mu_rho_step,
                           psi_rho_step,
                           rho_step,
                           step_Z
){
  mu0<-mu0
  Sigma0<-Sigma0
  a0<-a0
  b0<-b0
  #alpha_mustr<-alpha_mustr
  #beta_mustr<-beta_mustr
  #a_psib<-a_psib
  #b_psib<-b_psib
  v12<-fn.compute.ab(mu_bstr,psi_bstr)
  v1<-v12[1]
  v2<-v12[2]
  swSq<-swSq
  w1<-w1; w2<-w2
  a_psir<-a_psir
  b_psir<-b_psir
  ############
  XtX<-lapply(X, FUN=function(X) t(X)%*%X)
  p<-length(mu0) #the number of reg. coefficients per group
  pTot<-length(X)*p
  #n<-length(unlist(y))
  nGroups<-pTot/p
  Sigma0Inv<-solve(Sigma0)
  total<-nkeep+nburn
  #initialize outputs
  #list of betaSamples: each element of list is betaSamples from corresponding group
  betaGroupSamples<-array(NA, c(p,nGroups, total))
  BetaSamples<-matrix(NA,total,p)
  sigma2Samples<-matrix(NA,total,nGroups)
  Zsamples<-matrix(NA,total,nGroups)
  #mu_bstrSamples<-numeric(total)
  #psi_bstrSamples<-numeric(total)
  bstarSamples<-numeric(total)
  mu_rhoSamples<-numeric(total)
  psi_rhoSamples<-numeric(total)
  rhoSamples<-numeric(total)
  tstbst<-matrix(NA, total, 4)


  #initial values
  #starting values of parameters
  #mu_bstr<-rbeta(1, alpha_mustr,beta_mustr)
  #psi_bstr<-rgamma(1, a_psib,b_psib)
  #v<-fn.compute.ab(mu_bstr,psi_bstr)
  bstar<-mu_bstr #rbeta(1, v[1],v[2]) qbeta(.5,v[1],v[2]) #
  Beta<-mvrnorm(1,mu=mu0, Sigma=(1-bstar)*Sigma0)
  betalMat<-matrix(Beta, p,nGroups)
  mu_rho<-rbeta(1,w1,w2)
  psi_rho<-rgamma(1, a_psir,b_psir)
  ab<-fn.compute.ab(mu_rho,psi_rho)
  #   rho<-rbeta(1,ab[1],ab[2])
  rho<-runif(1, .1,.9)
  Sigma_rho<-(1-rho)*diag(nGroups)+matrix(1, nGroups,nGroups)*rho

  Z<-mvrnorm(1, rep(0,nGroups), Sigma_rho)



  for(iter in 1:total){
    samp <- fn.hier.one.rep(y,
                          X,
                          XtX,
                          v1,#mu_bstr,
                          v2,#psi_bstr,
                          bstar,
                          Beta,
                          betalMat, #first to update
                          Z,
                          mu_rho,
                          psi_rho,
                          rho,
                          step_logbstar,
                          mu_rho_step,
                          psi_rho_step,
                          rho_step,
                          step_Z,
                          Sigma0Inv)
    #update temp values
    Beta<-samp$Beta
    betalMat<-samp$betalMat
    Z<-samp$Z
    sigma2<-samp$sigma2
    #mu_bstr<-samp$mu_bstr
    #psi_bstr<-samp$psi_bstr
    bstar<-samp$bstar
    mu_rho<-samp$mu_rho
    psi_rho<-samp$psi_rho
    rho<-samp$rho
    #update outputs
    BetaSamples[iter,]<-Beta
    betaGroupSamples[,,iter]<-betalMat
    sigma2Samples[iter,]<-sigma2
    Zsamples[iter,]<-Z
    #mu_bstrSamples[iter]<-mu_bstr
    #psi_bstrSamples[iter]<-psi_bstr
    bstarSamples[iter]<-bstar
    mu_rhoSamples[iter]<-mu_rho
    psi_rhoSamples[iter]<-psi_rho
    rhoSamples[iter]<-rho
  }
  out<-list()
  out$Beta<-BetaSamples[-c(1:nburn),]
  out$betal<-betaGroupSamples[,,-c(1:nburn)]
  out$sigma2s<-sigma2Samples[-c(1:nburn),]
  out$Z<-Zsamples[-c(1:nburn),]
  #out$mu_bstr<-mu_bstrSamples[-c(1:nburn)]
  #out$psi_bstr<-psi_bstrSamples[-c(1:nburn)]
  out$bstar<-bstarSamples[-c(1:nburn)]
  out$mu_rho<-mu_rhoSamples[-c(1:nburn)]
  out$psi_rho<-psi_rhoSamples[-c(1:nburn)]
  out$rho<-rhoSamples[-c(1:nburn)]
  hypers<-c(a0,b0,swSq,w1,w2,a_psir,b_psir, mu_bstr,psi_bstr)
  names(hypers)<-c("a0","b0","swSq",'w1',"w2","a_psir","b_psir", "mu_bstr","psi_bstr")
  out$hypers<-hypers
  out
}



#' Function to fit restricted likelihood version of hierarchical model in paper
#' @rdname hier_models
#' @inheritParams hierNormTheoryLm
#' @param regEst Regression estimator on which to condition . Either Huber or Tukey.
#' @param scaleEst Scale estimator on which to condition('Huber' is only option here)
#' @export
hierNormTheoryRestLm <- function(y,
                                X,
                                regEst='Huber',
                                scaleEst='Huber',
                                nkeep=1e4, nburn=1e3,
                                mu0,
                                Sigma0,
                                a0,
                                b0,
                                mu_bstr,
                                psi_bstr,
                                swSq=1,
                                w1,
                                w2,
                                a_psir,
                                b_psir,
                                maxit=400,
                                step_logbstar,
                                mu_rho_step,
                                psi_rho_step,
                                rho_step,
                                step_Z)
{
  #y is a list of the responses for each group
  #X is the design Matrix-a list of the design matrices Xi for each group
  #X[[i]] is the design matrix for each group

  mu0<-mu0
  Sigma0<-Sigma0
  a0<-a0
  b0<-b0
  #alpha_mustr<-alpha_mustr
  #beta_mustr<-beta_mustr
  #a_psib<-a_psib
  #b_psib<-b_psib
  swSq<-swSq
  w1<-w1; w2<-w2
  a_psir<-a_psir
  b_psir<-b_psir
  v12<-fn.compute.ab(mu_bstr,psi_bstr)
  v1<-v12[1]
  v2<-v12[2]
  ############
  projList=NULL #leave this as null; projection onto deviation space for each group
  XtX<-lapply(X, FUN=function(X) t(X)%*%X)
  p<-length(mu0) #the number of reg. coefficients per group
  pTot<-length(X)*p
  ni<-sapply(y, length)
  nGroups<-pTot/p
  Sigma0Inv<-solve(Sigma0)
  total<-nkeep+nburn
  #initialize outputs
  #list of betaSamples: each element of list is betaSamples from corresponding group
  betaGroupSamples<-array(NA, c(p,nGroups, total))
  BetaSamples<-matrix(NA,total,p)
  sigma2Samples<-matrix(NA,total,nGroups)
  #mu_bstrSamples<-numeric(total)
  #psi_bstrSamples<-numeric(total)
  bstarSamples<-numeric(total)
  mu_rhoSamples<-numeric(total)
  psi_rhoSamples<-numeric(total)
  rhoSamples<-numeric(total)
  yAccept<-matrix(NA,total,nGroups)
  logprop.curr<-matrix(NA,total,nGroups)
  yMhRatios<-matrix(NA,total,nGroups)

  #initial values
  #starting values of parameters
  #mu_bstr<-rbeta(1, alpha_mustr,beta_mustr)
  #psi_bstr<-rgamma(1, a_psib,b_psib)
  #v<-fn.compute.ab(mu_bstr,psi_bstr)
  bstar<-mu_bstr#rbeta(1, v[1],v[2])
  Beta<-mvrnorm(1,mu=mu0, Sigma=(1-bstar)*Sigma0)
  betalMat<-matrix(Beta, p,nGroups)
  mu_rho<-rbeta(1,w1,w2)
  psi_rho<-rgamma(1, a_psir,b_psir)
  ab<-fn.compute.ab(mu_rho,psi_rho)
  #rho<-rbeta(1,ab[1],ab[2])
  rho<-runif(1, .1, .9)
  Sigma_rho<-(1-rho)*diag(nGroups)+matrix(1, nGroups,nGroups)*rho
  Z<-mvrnorm(1, rep(0,nGroups), Sigma_rho)




  ############################
  #define the psi function
  ############################
  if(regEst=='Huber') {
    psi<-get('psi.huber') #internal
    fn.psi<-get('fn.psi.huber')

  } else {
    if(regEst=='Tukey'){
      psi<-get('psi.bisquare') #internal
      fn.psi<-get('fn.psi.bisquare')
    } else {stop("only set up for Huber or Tukey regression estimates")}}

  if(scaleEst!='Huber'){
    stop('scale estimate only set up for Hubers Prop2')
  }
  fn.chi<-fn.chi.prop2

  robustList<-list();length(robustList)<-nGroups
  bHatObsList<-list();length(bHatObsList)<-nGroups
  sigHatObsList<-list();length(sigHatObsList)<-nGroups

  #fit the robust linear model to each group separately
  #modified for small groups?
  for(groups in 1:nGroups){
    robust<-rlm(X[[groups]],y[[groups]],psi=psi, scale.est=scaleEst, maxit=maxit)
    robustList[[groups]]<-robust
    bHatObsList[[groups]]<-robust$coef
    sigHatObsList[[groups]]<-robust$s
  }

  #####################
  #choose starting value for yi
  #####################
  log.prop.den.curr<-numeric(nGroups)
  yCurr<-list(); length(yCurr)<-nGroups
  for(i in 1:nGroups){
    yCurr[[i]]<-numeric(ni[i])
  }
  if(is.null(projList)) {projList<-list(); length(projList)<-nGroups
  QtList<-list(); length(QtList)<-nGroups
  for(i in 1:nGroups){
    Q<-qr.Q(qr(X[[i]]))
    projList[[i]]<-diag(ni[i])-tcrossprod(Q,Q)
    QtList[[i]]<-t(Q)
  }
  }
  #finding starting yi's for each group
  for(i in 1:nGroups){
    y.prop <- rnorm(ni[i])
    Xi<-X[[i]]
    bHatObsi<-bHatObsList[[i]] #observed stats in each group
    sigHatObsi<-sigHatObsList[[i]]
    proji<-projList[[i]]
    Qti<-QtList[[i]]
    yCurri <- fn.comp.ystst(y.prop,Xi,l1obs=bHatObsi,s1obs=sigHatObsi,psi,scale.est=scaleEst,maxit)
    yCurri<-yCurri
    yCurr[[i]]<-yCurri
    log.prop.den.curr[i]<-log.prop.den2(yCurri,Xi, proji,Qti,bHatObsi, sigHatObsi,fn.psi, fn.chi,ni[i],p)
  }


  ##################
  #MCMC
  #################
  for(iter in 1:total){

    #[theta|everything] step
    samp<-fn.hier.one.rep(yCurr,
                          X,
                          XtX,
                          v1,#mu_bstr,
                          v2,#psi_bstr,
                          bstar,
                          Beta,
                          betalMat, # first to update
                          Z,
                          mu_rho,
                          psi_rho,
                          rho,
                          step_logbstar,
                          mu_rho_step,
                          psi_rho_step,
                          rho_step,
                          step_Z,
                          Sigma0Inv)
    #update temp values
    Beta<-samp$Beta
    betalMat<-samp$betalMat
    Z<-samp$Z
    sigma2<-samp$sigma2
    #mu_bstr<-samp$mu_bstr
    #psi_bstr<-samp$psi_bstr
    bstar<-samp$bstar
    mu_rho<-samp$mu_rho
    psi_rho<-samp$psi_rho
    rho<-samp$rho
    #update outputs
    BetaSamples[iter,]<-Beta
    betaGroupSamples[,,iter]<-betalMat
    sigma2Samples[iter,]<-sigma2
    # mu_bstrSamples[iter]<-mu_bstr
    #psi_bstrSamples[iter]<-psi_bstr
    bstarSamples[iter]<-bstar
    mu_rhoSamples[iter]<-mu_rho
    psi_rhoSamples[iter]<-psi_rho
    rhoSamples[iter]<-rho

    #[y|everything + robust statistics] step; each group updated separately
    for(gp in 1:nGroups){
      yicurr<-yCurr[[gp]]
      Xi<-X[[gp]]
      betaCuri<-betalMat[,gp]
      sigma2Curi<-sigma2[gp]
      bHatObsi<-bHatObsList[[gp]]
      sigHatObsi<-sigHatObsList[[gp]]
      log.prop.den.curri<-log.prop.den.curr[gp]
      proj<-projList[[gp]]
      Qt<-QtList[[gp]]
      ySample<-fn.one.rep.y2(yicurr,betaCuri,sqrt(sigma2Curi),bHatObsi, sigHatObsi,Xi, log.prop.den.curri, proj,Qt,fn.psi,fn.chi, psi,scaleEst,maxit)
      yCurr[[gp]]<-ySample[[1]]
      yAccept[iter,gp]<-ySample[[2]]
      log.prop.den.curr[gp]<-ySample[[3]]
      logprop.curr[iter,gp]<-ySample[[3]]
      yMhRatios[iter,gp]<-ySample[[4]]
    }
  }

  out<-list()
  out$Beta<-BetaSamples[-c(1:nburn),]
  out$betal<-betaGroupSamples[,,-c(1:nburn)]
  out$sigma2s<-sigma2Samples[-c(1:nburn),]
  # out$mu_bstr<-mu_bstrSamples[-c(1:nburn)]
  #  out$psi_bstr<-psi_bstrSamples[-c(1:nburn)]
  out$bstar<-bstarSamples[-c(1:nburn)]
  out$mu_rho<-mu_rhoSamples[-c(1:nburn)]
  out$psi_rho<-psi_rhoSamples[-c(1:nburn)]
  out$rho<-rhoSamples[-c(1:nburn)]
  #   out$yAccept<-yAccept #colMeans(yAccept)
  out$yAccept<-colMeans(yAccept)
  out$rlmFits<-robustList
  hypers<-c(a0,b0,swSq,w1,w2,a_psir,b_psir, mu_bstr,psi_bstr)
  names(hypers)<-c("a0","b0","swSq",'w1',"w2","a_psir","b_psir", 'mu_bstr', "psi_bstr")
  out$hypers<-hypers
  out
}



#' @rdname hier_models
#' @param betalsamples the array of betals: in the specific format: the 3rd dimension is the groups. columns represent samples, row represent slopes
#' @param sigma2lsamples: columns represnt groups, rows represent samples
#' @param yhold list of holdout samples
#' @param Xhold list of design matrices
#' @export
fn.compute.marginals.hierModelNormal<-function(betalsamples, sigma2lsamples, yhold,Xhold){


  fits<-function(betahats, X){X%*%betahats}

  betalSampsList <- lapply(1:dim(betalsamples)[3],function(x) betalsamples[,,x])
  names(betalSampsList) <- names(yhold)
  #holdout means

  nMuList<-mapply(fits, betalSampsList, Xhold)

  #sd's across samples for each holdout set
  nSigmaList<-split(sqrt(sigma2lsamples), rep(1:ncol(sigma2lsamples), each = nrow(sigma2lsamples)))
  names(nSigmaList)<-names(yhold)

  mapply(FUN=function(yh,mean,sd){
    rowMeans(dnorm(yh,mean,matrix(rep(sd, length(yh)), length(yh), length(sd),byrow=TRUE)))
  }
  ,yhold,nMuList,nSigmaList)

}



#' @rdname hier_models
#' @param y data
#' @param mean center of t distribution
#' @param sigma scale of t distribution
#' @param nu  degrees of freedom
#' @export
tdensity<-function(y, mean, sigma, nu){
  (gamma(.5*(nu+1))/(gamma(.5*nu)*sigma*sqrt(nu*pi)))*(1+((y-mean)/sigma)^2/nu)^(-.5*(nu+1))
}


#' @rdname hier_models
#' @param betalsamples the array of betals: in the specific format: the 3rd dimension is the groups. columns represent samples, row represent slopes
#' @param sigma2lsamples: columns represnt groups, rows represent samples
#' @param yhold list of holdout samples
#' @param Xhold list of design matrices
#' @param nu degrees of freedom in t model
#' @export
fn.compute.marginals.hierModelTmodel <- function(betalsamples, sigma2lsamples, yhold,Xhold){

  fits<-function(betahats, X){X%*%betahats}

  betalSampsList<-lapply(1:dim(betalsamples)[3],function(x) betalsamples[,,x])
  names(betalSampsList)<-names(yhold)
  #holdout means
  nMuList<-mapply(fits, betalSampsList, Xhold)

  #sd's across samples for each houldout set
  nSigmaList<-split(sqrt(sigma2lsamples), rep(1:ncol(sigma2lsamples), each = nrow(sigma2lsamples)))
  names(nSigmaList)<-names(yhold)

  mapply(FUN=function(yh,mean,sd){
    rowMeans(tdensity(yh,mean,matrix(rep(sd, length(yh)), length(yh), length(sd),byrow=TRUE), nu = nu))
  }
  ,yhold,nMuList,nSigmaList)

}






