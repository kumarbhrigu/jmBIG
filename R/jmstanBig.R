#' @title  Joint model for BIG data using rstanarm
#' @description function for joint model in BIG DATA using \code{rstanarm} package
#' @param dtlong longitudinal dataset, which contains id,visit time,longitudinal measurements along with various covariates
#' @param dtsurv survival dataset corresponding to the longitudinal dataset, with survival status and survival time
#' @param longm model for longitudinal response
#' @param survm survival model
#' @param samplesize sample size to divide the Big data
#' @param time_var time variable in longitudinal model, included in the longitudinal data
#' @param id name of id column in longitudinal dataset
#' @param nchain number of chain for MCMC
#' @param refresh refresh rate for MCMC chain
#' @return returns a list containing various output which are useful for prediction.
#' @importFrom rstanarm stan_jm
#' @importFrom stats pnorm
#' @export
#' @references Goodrich, B., et al. "rstanarm: Bayesian applied regression modeling via Stan. R package version 2.17. 4." Online< http://mc-stan. org (2018).
#' @examples
#'
#'  \donttest{
#' ##
#' library(survival)
#' library(dplyr)
#' fit3<-jmstanBig(dtlong=long2,dtsurv = surv2,longm=y~ x7+visit+(1|id),
#' survm=Surv(time,status)~x1+visit,samplesize=200,time_var='visit',id='id')
#' P2<-postTraj(model<-fit3,m<-1,ids<-c(1,2,100))
#' pp1<-plot(P2$p1[[1]],plot_observed = TRUE)
#' pp2<-plot(P2$p1[[2]],plot_observed = TRUE)
#' pp3<-plot(P2$p1[[3]],plot_observed = TRUE)
#' ##
#' }
#'
#' @author Atanu Bhattacharjee, Bhrigu Kumar Rajbongshi and Gajendra Kumar Vishwakarma
#' @seealso   \link{jmbayesBig},\link{jmcsBig},\link{joinRMLBig}

jmstanBig<-function(dtlong,dtsurv,longm,survm,samplesize=50,time_var,id,nchain=1,refresh=2000){

  cl<-match.call()
  dtlong<-dtlong
  dtsurv<-dtsurv
  if(!id%in%names(dtlong) ){
    stop("\n Longitudinal data must have column 'id' ")
  }
  if(!id%in%names(dtsurv) ){
    stop("\n Survival data must have column 'id' ")
  }
  if(!names(dtlong)[names(dtlong)==id]==names(dtsurv)[names(dtsurv)==id]){
    stop("\n'dtlong' and 'dtsurv' must have same id.")
  }
  if(!time_var%in%names(dtlong)){
    stop("\n 'timeVar' should be in longitudinal dataset")
  }

  # Preparing the data


  if(names(dtlong)[names(dtlong)==id]=='id'){dtlong<-dtlong}else{
    dtlong<-dtlong; names(dtlong)[names(dtlong)==id]<-'id'}
  if(names(dtsurv)[names(dtsurv)==id]=='id'){dtsurv<-dtsurv}else{
    dtsurv<-dtsurv; names(dtsurv)[names(dtsurv)==id]<-'id'}

  nr<-nrow(dtsurv)
  samplesize=samplesize# sample of size 50 from the survival data
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

  rlist<-list();uprlist<-list()
  upselist<-list();upstlist<-list()
  updatestimate<-0;updatese<-0;updatest<-0

  for(i in 1:ceiling(nr/samplesize)){#floor(nr/samplesize)

    # Y ~ X1+X2+time+(1|id)
    # survival::Surv(Time, status) ~ X1+X2
    mod1 <- stan_jm(formulaLong =longm,
                    dataLong = dtlong1[[i]],
                    formulaEvent = survm,
                    dataEvent = dtsurv1[[i]],
                    time_var = time_var,
                    chains = nchain, refresh = refresh, seed = 12345)
    rlist[[i]]<-mod1
    #uprlist[[i]]<-summary(mod1, probs = c(.025,.975))
    # here we have considered only the estimate part fro the mod1 object
    uprlist[[i]]<-c(mod1$coefficients[[1]],mod1$coefficients$Event)
    upselist[[i]]<-c(mod1$ses[[1]],mod1$ses$Event)
    upstlist[[i]]<-mod1$stan_summary

    if(i==1){updatestimate<-uprlist[[1]]; updatese<-upselist[[1]]; updatest<-upstlist[[1]] }
    updatestimate<-suppressWarnings((updatestimate+uprlist[[i]])/2)
    updatese<-suppressWarnings((updatese+upselist[[i]])/2)
    updatest<-(updatest+upstlist[[i]])/2
  }

  # update in mod1 object

  uprlist1<-rlist
  for(i in 1:length(rlist)){
    uprlist1[[i]]$coefficients[[1]]<-updatestimate[1:length(mod1$coefficients[[1]])]
    uprlist1[[i]]$coefficients$Event<-updatestimate[(length(mod1$coefficients[[1]])+1):length(updatestimate)]
    uprlist1[[i]]$ses[[1]]<-updatese[1:length(mod1$ses[[1]])]
    uprlist1[[i]]$ses$Event<-updatese[(length(mod1$ses[[1]])+1):length(updatese)]
    uprlist1[[i]]$stan_summary<-updatest
  }

  mod11<-NULL
  mod11<-mod1
  mlist<-updatestimate
  mod11$coefficients[[1]]<-updatestimate[1:length(mod1$coefficients[[1]])]
  mod11$coefficients$Event<-updatestimate[(length(mod1$coefficients[[1]])+1):length(updatestimate)]
  mod11$ses[[1]]<-updatese[1:length(mod1$ses[[1]])]
  mod11$ses$Event<-updatese[(length(mod1$ses[[1]])+1):length(updatese)]
  mod11$stan_summary<-updatest
  #Results
  result<-list()
  result$call<-cl
  result$allmodel<-rlist
  result$uprlist<-uprlist1
  result$pseudoMod<-mod11
  result$nr<-nr
  result$samplesize<-samplesize
  class(result)<-'jmstanBig'
  #message('Results for the joint model')
  result


}





