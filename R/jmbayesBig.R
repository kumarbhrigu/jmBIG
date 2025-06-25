




#############################################################
## function for joint model in BIG DATA using Jmbayes2 ######


#' @title Joint model for BIG data using JMbayes2
#' @description function for joint model in BIG DATA using \code{JMbayes2}
#' @param dtlong longitudinal dataset, which contains id,visit time,longitudinal measurements along with various covariates
#' @param dtsurv survival dataset corresponding to the longitudinal dataset, with survival status and survival time
#' @param longm fixed effect model for longitudinal response
#' @param survm survival model
#' @param samplesize sample size to divide the Big data
#' @param rd random effect model part
#' @param timeVar time variable in longitudinal model, included in the longitudinal data
#' @param nchain number of chain for MCMC
#' @param id name of id column in longitudinal dataset
#' @param niter number of iteration for MCMC chain
#' @param nburnin number of burnin sample for MCMC chain
#' @return returns a list containing various output which are useful for prediction.
#' @importFrom  JMbayes2 jm
#' @importFrom nlme lme
#' @import survival
#' @importFrom stats pnorm
#' @export
#' @references Rizopoulos, D., G. Papageorgiou, and P. Miranda Afonso. "JMbayes2: extended joint models for longitudinal and time-to-event data." R package version 0.2-4 (2022).
#' @examples
#'
#'  \donttest{
#' ##
#' library(survival)
#' library(nlme)
#' library(dplyr)
#' fit5<-jmbayesBig(dtlong=long2,dtsurv = surv2,longm=y~ x7+visit,survm=Surv(time,status)~x1+visit,
#' rd= ~ visit|id,timeVar='visit',nchain=1,samplesize=200,id='id')
#' ydt<-long2%>%filter(id%in%c(900))
#' cdt<-surv2[,'id']%>%filter(id%in%c(900))
#' newdata<-full_join(ydt,cdt,by='id')
#' P2<-predJMbayes(model<-fit5,ids<-c(900),newdata=newdata,process = 'event')
#' plot(P2$p1[[1]])
#' ##
#' }
#' @author Atanu Bhattacharjee, Bhrigu Kumar Rajbongshi and Gajendra Kumar Vishwakarma
#' @seealso   \link{jmcsBig},\link{jmstanBig},\link{joinRMLBig}

jmbayesBig<-function(dtlong,dtsurv,longm,survm,samplesize=50,rd,timeVar,nchain=1,id,niter=2000,nburnin=1000){

  cl<-match.call()

  if(!id%in%names(dtlong) ){
    stop("\n Longitudinal data must have column 'id' ")
  }
  if(!id%in%names(dtsurv) ){
    stop("\n Survival data must have column 'id' ")
  }
  if(!names(dtlong)[names(dtlong)==id]==names(dtsurv)[names(dtsurv)==id]){
    stop("\n'dtlong' and 'dtsurv' must have same id.")
  }
  if(!timeVar%in%names(dtlong)){
    stop("\n 'timeVar' should be in longitudinal dataset")
  }

  # Preparing the data

  dtlong<-dtlong
  dtsurv<-dtsurv
  if(names(dtlong)[names(dtlong)==id]=='id'){dtlong<-dtlong}else{
    dtlong<-dtlong; names(dtlong)[names(dtlong)==id]<-'id'}
  if(names(dtsurv)[names(dtsurv)==id]=='id'){dtsurv<-dtsurv}else{
    dtsurv<-dtsurv; names(dtsurv)[names(dtsurv)==id]<-'id'}

  nr<-nrow(dtsurv)
  samplesize<-samplesize# sample of size 50 from the survival data
  dtsurv1<-split(dtsurv,rep(1:ceiling(nr/samplesize),each=samplesize,length.out=nr))

  lid<-nr-(floor(nr/samplesize)*samplesize)
  if(lid!=0){
    maxdt<-length(dtsurv1)
    leftID<-dtsurv1[[maxdt]]$id
    dtsurv1[[maxdt]]<-rbind(dtsurv1[[maxdt-1]][(length(leftID)+1):samplesize,],dtsurv1[[maxdt]])
    #dtsurv1[[maxdt]]<-NULL
  }else{
    dtsurv1<-dtsurv1
  }
  f1<-function(x){
    return(max(unique(x$id)))
  }
  nk<-as.vector(unlist(lapply(dtsurv1,f1)))
  nk1<-c()
  for(j in 1:length(nk)){
    nk1[j]<-nrow(dtlong[dtlong$id<=nk[j],])
  }
  nk2<-c(0,nk1[-length(nk1)])
  nk3<-nk1-nk2
  dtlong1<-split(dtlong,rep(1:ceiling(nr/samplesize),times=as.vector(nk3)))
  dtlong1[[length(dtlong1)]]<-dtlong[dtlong$id%in%dtsurv1[[length(dtsurv1)]]$id,]

  rlist<-list();betalist<-list()
  gammalist<-list();alphalist<-list();vcovlist<-list()
  betaselist<-list();gammaselist<-list();alphaselist<-list()
  updatebeta<-0;updategamma<-0;updatealpha<-0;updatebetaSE<-0;
  updategammaSE<-0;updatealpahSE<-0

  for(i in 1:ceiling(nr/samplesize)){#floor(nr/samplesize)

    fm1<-lme(fixed=longm,random=rd,data=dtlong1[[i]])
    fm2<-coxph(survm,data=dtsurv1[[i]])

    mod1 <- jm(fm2,fm1,time_var=timeVar,n_chains=nchain,n_iter=niter,n_burnin=nburnin
    )
    rlist[[i]]<-mod1
    #uprlist[[i]]<-summary(mod1, probs = c(.025,.975))
    # here we have considered only the estimate part fro the mod1 object
    betalist[[i]]<-mod1$statistics$Mean$betas1
    gammalist[[i]]<-mod1$statistics$Mean$gammas
    alphalist[[i]]<-mod1$statistics$Mean$alphas
    betaselist[[i]]<-mod1$statistics$SE$betas1
    gammaselist[[i]]<-mod1$statistics$SE$gammas
    alphaselist[[i]]<-mod1$statistics$SE$alphas
    #vcovlist[[i]]<-mod1$statistics$post_vars

    if(i==1){
      updatebeta<-betalist[[1]]
      updategamma<-gammalist[[1]]
      updatealpha<-alphalist[[1]]
      updatebetaSE<-betaselist[[1]]
      updategammaSE<-gammaselist[[1]]
      updatealpahSE<-alphaselist[[1]]
    }
    updatebeta<-suppressWarnings((updatebeta+betalist[[i]])/2)
    #updategamma<-0;updatenu<-0
    updategamma<-suppressWarnings((updategamma+gammalist[[i]])/2)
    updatealpha<-(updatealpha+alphalist[[i]])/2
    #updatevcov<-(updatevcov+vcovlist[[i]])/2
    updatebetaSE<-(updatebetaSE+betaselist[[i]])/2
    updategammaSE<-(updategammaSE+gammaselist[[i]])/2
    updatealpahSE<-(updatealpahSE+alphaselist[[i]])/2
  }

  # update in mod1 object

  uprlist1<-rlist
  for(i in 1:length(rlist)){
    uprlist1[[i]]$statistics$Mean$betas1<-updatebeta
    uprlist1[[i]]$statistics$Mean$gammas<-updategamma
    uprlist1[[i]]$statistics$Mean$alphas<-updatealpha
    uprlist1[[i]]$statistics$SE$betas1<-updatebetaSE
    uprlist1[[i]]$statistics$SE$gammas<-updategammaSE
    uprlist1[[i]]$statistics$SE$alphas<-updatealpahSE

    # uprlist1[[i]]$vcov<-updatevcov
    #uprlist1[[i]]$ses$Event<-updatese[(length(mod1$ses[[1]])+1):length(updatese)]
    #uprlist1[[i]]$stan_summary<-updatest
  }

  mod11<-NULL
  mod11<-mod1
  #mlist<-updatestimate
  mod11$statistics$Mean$betas1<-updatebeta
  mod11$statistics$Mean$gammas<-updategamma
  mod11$statistics$Mean$alphas<-updatealpha
  mod11$statistics$SE$betas1<-updatebetaSE
  mod11$statistics$SE$gammas<-updategammaSE
  mod11$statistics$SE$alphas<-updatealpahSE
  #Results
  result<-list()
  result$call<-cl
  result$allmodel<-rlist
  result$uprlist<-uprlist1
  result$pseudoMod<-mod11
  result$nr<-nr
  result$samplesize<-samplesize
  class(result)<-'jmbayesBig'
  #message('Results for the joint model')
  result

}



