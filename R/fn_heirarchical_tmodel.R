#--------------------
#Sampling and fitting functions for T-model
#-------------------

#reformulation of the data model: See BDA book
#y_ij~N(t(x_ij)beta_j, V_ij)
#V_ij~Inv-Chisq(nu, sigma_j^2)


#steps changing in the Gibbs sampler are [beta_j|-] and [Z|-]


#' Sampling for beta_l's and v_ls in hierarchical T-model
#' @rdname beta_ls_t
#' @param yl,Xl desc
#' @param betal desc
#' @param SigmaVlInv desc
#' @param Beta desc
#' @param bstar desc
#' @param betal
#' @param sigma2l
#' @details Sampling beta.l and v.l for the T-model. fn.one.rep.betasAndVs puts the two together to sample both for all groups.
fn.sample.beta.lTmodel<-function(yl, Xl,SigmaVlInv,Beta,bstar){
  b<-bstar/swSq
  Sigma.betal<-solve(t(Xl)%*%SigmaVlInv%*%Xl+Sigma0Inv/b)
  mu.l<-Sigma.betal%*%(t(Xl)%*%SigmaVlInv%*%yl+Sigma0Inv%*%Beta/b)
  mvrnorm(1, mu.l, Sigma.betal)
}

#' @rdname beta_ls_t
fn.sample.v.lTmodel<-function(yl, Xl,betal,sigma2l){
  nl<-length(yl)
  as.vector(((yl-Xl%*%betal)^2+nu*sigma2l)/rchisq(nl,nu+1))
}
#' @rdname beta_ls_t
fn.one.rep.betasAndVs<-function(y, X,Beta,sigma2,betalMat, vlList,bstar){
  for(l in 1:ncol(betalMat)){
    yl<-y[[l]]
    Xl<-X[[l]]
    vis<-vlList[[l]]
    SigmaVlInv<-diag(1/vis)
    sigma2l<-sigma2[l]
    beta.l<-fn.sample.beta.lTmodel(yl, Xl,SigmaVlInv,Beta,bstar)
    betalMat[,l]<-beta.l
    vsl<-fn.sample.v.lTmodel(yl, Xl,beta.l,sigma2l)
    vlList[[l]]<-vsl
  }
  out<-list()
  out$betalMat<-betalMat
  out$vlList<-vlList
  out
}

#--------------------
#[Z|-]
#--------------------


#function for the likelihood of one z.

#' Sampling for the z's and in hierarchical T-model
#' @rdname z_t desc
#' @param zl,zlCur desc
#' @param zminus desc
#' @param yl desc
#' @param vsl desc
#' @param condsd desc
#' @param betal desc
#' @param condMeanMult desc
#' @param step desc.
#' @details Sampling z for the hierarchical T-model.
fn.logLike.zlTmodel<-function(zl,zminus, yl,vsl, condsd,condMeanMult){
  #vsl is this vector of vij's
  #zl is the compenent I am updating, zminus is all other components; condvar is the condition sd for the conditional normal (zl|zminus): this depends only on current
  #condMeanMult%*%zminus is the conditional mean for zl|zminus
  condMean<-condMeanMult%*%zminus
  sigma2l<-fn.compute.sigma2(zl,a0,b0)
  nl<-length(yl)
  dnorm(zl,condMean,condsd,log=TRUE)+.5*nl*nu*log(sigma2l)-.5*nu*sigma2l*sum(1/vsl)
}
#'@rdname z_t
fn.proposal.zlTmodel<-function(zlCur, step=1){
  rnorm(1, zlCur,step)
}


#'@rdname z_t
fn.sample.ZlTmodel<-function(zlCur,zminus,yl,vsl, condsd,condMeanMult, step=1){
  zlProp<-fn.proposal.zlTmodel(zlCur, step)
  logRat1<-fn.logLike.zlTmodel(zlProp,zminus, yl,vsl, condsd,condMeanMult)-fn.logLike.zlTmodel(zlCur,zminus, yl,vsl,condsd,condMeanMult)
  logRat2<-0 #random walk: symmetric
  mhRat3<-min(exp(logRat1+logRat2),1)
  if(runif(1)<mhRat3){
    return(zlProp)
  } else {
    return(zlCur)
  }
}

#'@rdname z_t
fn.one.rep.ZTmodel<-function(Zcur, y,X, betalMat,vlList, rho,step=rep(1,length(Zcur))){
  Sigma_rhoInvMinus<-fn.compute.SigmaRhoInv(rho, K=length(Zcur)-1)    #inverse (1-pho)I+rhoJ for K-1 by K-1 dim matrices
  rho_vec<-rep(rho, length(Zcur)-1)
  condsd<-sqrt(1-t(rho_vec)%*%Sigma_rhoInvMinus%*%rho_vec)
  condMeanMult<-t(rho_vec)%*%Sigma_rhoInvMinus
  for(j in 1:length(Zcur)){
    zlCur<-Zcur[j]
    zminus<-Zcur[-j]
    yl<-y[[j]]
    Xl<-X[[j]]
    vsl<-vlList[[j]]
    betal<-betalMat[,j]
    zlCur<-fn.sample.ZlTmodel(zlCur,zminus,yl,vsl, condsd,condMeanMult,step[j])
    Zcur[j]<-zlCur
  }
  Zcur
}


#' @rdname hier_models
#' @inheritParams hierNormTheoryLm
#' @param vlList
fn.one.rep.tHierModel<-function(y,
                                X,
                                #XtX,
                                v1,#mu_bstr,
                                v2,#psi_bstr,
                                bstar,
                                Beta,
                                betalMat, #first to update
                                vlList,
                                Z,
                                mu_rho,
                                psi_rho,
                                rho,
                                step_logbstar,
                                mu_rho_step,
                                psi_rho_step,
                                rho_step,
                                step_Z
){ #y, X are list of group level responses
  #fn.one.rep.beta.l loops through for each beta.l
  #[beta_i|-]
  sigma2<-fn.compute.sigma2(Z,a0,b0)

  betalVSamp<-fn.one.rep.betasAndVs(y, X,Beta,sigma2,betalMat, vlList,bstar)
  betalMat<-betalVSamp$betalMat
  vlList<-betalVSamp$vlList

  #[Z|-] i.e. the sigma2_is
  #fn.one.rep.Z loops through for all the z_i
  Z<-fn.one.rep.ZTmodel(Z, y,X, betalMat,vlList, rho,step=step_Z)
  sigma2<-fn.compute.sigma2(Z,a0,b0)
  #[Beta|-]
  Beta<-fn.sample.Beta(betalMat,bstar)

  #[mu_bstr|-]
  # mu_bstr<-fn.sample.mubstr(mu_bstr,psi_bstr, bstar,step_mubstr)

  #[psi_bstr|-]
  #  psi_bstr<-fn.sample.psibstr(psi_bstr,mu_bstr, bstar)

  #[bstar|-]
  quad1<-fn.compute.quad1(Beta, mu0)
  quad2<-fn.compute.quad2(Beta, betalMat)
  bstar<-fn.sample.bstar(bstar,v1,v2,quad1,quad2,K=nGroups,p, step_logbstar)
  #[mu_rho|-]
  mu_rho<-fn.sample.mu_rho(mu_rho,psi_rho, rho, mu_rho_step)

  #[psi_rho|-]
  psi_rho<-fn.sample.psi_rho(psi_rho,mu_rho, rho,  psi_rho_step)

  #[rho|-]
  rho<-fn.sample.rho(rho, mu_rho, psi_rho, Z,rho_step)

  out<-list()
  out$Beta<-Beta
  out$betalMat<-betalMat
  out$Z<-Z
  out$sigma2<-sigma2
  out$vlList<-vlList
  #out$mu_bstr<-mu_bstr
  #out$psi_bstr<-psi_bstr
  out$bstar<-bstar
  out$mu_rho<-mu_rho
  out$psi_rho<-psi_rho
  out$rho<-rho
  out
}




#' @rdname hier_models
#' @inheritParams hierNormTheoryLm
#' @param nu fixed degrees of freedom for assumed t-distribution
#' @export
hier_TLm<-function(y,
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
                        nu, #degrees of freedom for t-dist.,
                        step_logbstar,
                        mu_rho_step,
                        psi_rho_step,
                        rho_step,
                        step_Z
){
  mu0<<-mu0
  Sigma0<<-Sigma0
  a0<<-a0
  b0<<-b0
  #alpha_mustr<<-alpha_mustr
  #beta_mustr<<-beta_mustr
  #a_psib<<-a_psib
  #b_psib<<-b_psib
  v12<<-fn.compute.ab(mu_bstr,psi_bstr)
  v1<<-v12[1]
  v2<<-v12[2]
  swSq<<-swSq
  w1<<-w1; w2<<-w2
  a_psir<<-a_psir
  b_psir<<-b_psir
  ############
  #XtX<-lapply(X, FUN=function(X) t(X)%*%X)
  p<<-length(mu0) #the number of reg. coefficients per group
  pTot<<-length(X)*p
  ni<<-sapply(y, length, simplify=TRUE)
  nGroups<<-pTot/p
  Sigma0Inv<<-solve(Sigma0)
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
  Vs<-matrix(NA, total, length(unlist(y)))

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

  vlList<-list(); length(vlList)<-nGroups
  for(i in 1:nGroups){
    vlList[[i]]<-nu*fn.compute.sigma2(Z[i],a0,b0)/rchisq(ni[i],nu)
  }

  for(iter in 1:total){
    samp<-fn.one.rep.tHierModel(y,
                                X,
                                #XtX,
                                v1,#mu_bstr,
                                v2,#psi_bstr,
                                bstar,
                                Beta,
                                betalMat, #first to update
                                vlList,
                                Z,
                                mu_rho,
                                psi_rho,
                                rho)
    #update temp values
    Beta<-samp$Beta
    betalMat<-samp$betalMat
    Z<-samp$Z
    sigma2<-samp$sigma2
    vlList<-samp$vlList
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
    Vs[iter,]<-unlist(vlList)
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
  out$Vs<-Vs[-c(1:nburn),]
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




