#' @title  Joint model for BIG data using FastJM
#' @description function for joint model in BIG DATA using \code{FastJM}
#' @param dtlong longitudinal dataset, which contains id,visit time,longitudinal measurements along with various covariates
#' @param dtsurv survival dataset corresponding to the longitudinal dataset, with survival status and survival time
#' @param longm model for longitudinal response
#' @param survm survival model
#' @param samplesize sample size to divide the Big data
#' @param rd random effect part
#' @param id name of id column in longitudinal dataset
#' @return returns a list containing various output which are useful for prediction.
#' @importFrom FastJM jmcs
#' @importFrom stats pnorm
#' @export
#' @references Li, Shanpeng, et al. "Efficient Algorithms and Implementation of a Semiparametric Joint Model for Longitudinal and Competing Risk Data: With Applications to Massive Biobank Data." Computational and Mathematical Methods in Medicine 2022 (2022).
#' @examples
#'   \donttest{
#' ##
#' library(survival)
#' library(dplyr)
#' fit2<-jmcsBig(dtlong=data.frame(long2),dtsurv = data.frame(surv2),
#' longm=y~ x7+visit,survm=Surv(time,status)~x1+visit,rd= ~ visit|id,samplesize=200,id='id')
#' print(fit2)
#' ##
#'   }
#'
#'
#' @author Atanu Bhattacharjee, Bhrigu Kumar Rajbongshi and Gajendra Kumar Vishwakarma
#' @seealso   \link{jmbayesBig},\link{jmstanBig},\link{joinRMLBig}

jmcsBig<-function(dtlong,dtsurv,longm,survm,samplesize=50,rd,id){

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
  #if(!timeVar%in%names(dtlong)){
  #stop("\n 'timeVar' should be in longitudinal dataset")
  #}

  # Preparing the data

  dtlong<-dtlong
  dtsurv<-dtsurv
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

  rlist<-list();betalist<-list()
  gammalist<-list();nulist<-list();vcovlist<-list()
  updatebeta<-0;updategamma<-0;updatenu<-0;updatevcov<-0

  for(i in 1:ceiling(nr/samplesize)){#floor(nr/samplesize)

    # Y ~ X1+X2+time+(1|id)
    # survival::Surv(Time, status) ~ X1+X2
    # ydata=dtlong[dtlong$id<230,],cdata = dtsurv[dtsurv$id<230,],longm=Y ~ X1+X2+time,survm=Surv(Time, status) ~ X1+X2,rd=~1|id
    mod1 <- jmcs(long.formula =longm,
                 ydata = data.frame(dtlong1[[i]]),
                 surv.formula = survm,
                 cdata = data.frame(dtsurv1[[i]]),
                 random=rd
    )
    rlist[[i]]<-mod1
    #uprlist[[i]]<-summary(mod1, probs = c(.025,.975))
    # here we have considered only the estimate part fro the mod1 object
    betalist[[i]]<-mod1$beta
    gammalist[[i]]<-mod1$gamma1
    nulist[[i]]<-mod1$nu1
    vcovlist[[i]]<-mod1$vcov

    if(i==1){updatebeta<-betalist[[1]];updategamma<-gammalist[[1]]; updatenu<-nulist[[1]] ; updatevcov<-vcovlist[[1]]}
    updatebeta<-suppressWarnings((updatebeta+betalist[[i]])/2)
    #updategamma<-0;updatenu<-0
    updategamma<-suppressWarnings((updategamma+gammalist[[i]])/2)
    updatenu<-(updatenu+nulist[[i]])/2
    updatevcov<-(updatevcov+vcovlist[[i]])/2

  }

  # update in mod1 object

  uprlist1<-rlist
  for(i in 1:length(rlist)){
    uprlist1[[i]]$beta<-updatebeta
    uprlist1[[i]]$gamma1<-updategamma
    uprlist1[[i]]$nu1<-updatenu
    uprlist1[[i]]$vcov<-updatevcov
    #uprlist1[[i]]$ses$Event<-updatese[(length(mod1$ses[[1]])+1):length(updatese)]
    #uprlist1[[i]]$stan_summary<-updatest
  }

  mod11<-NULL
  mod11<-mod1
  #mlist<-updatestimate
  mod11$beta<-updatebeta
  mod11$gamma1<-updategamma
  mod11$nu1<-updatenu
  mod11$vcov<-updatevcov
  #Results
  result<-list()
  result$call<-cl
  result$allmodel<-rlist
  result$uprlist<-uprlist1
  result$pseudoMod<-mod11
  result$nr<-nr
  result$samplesize<-samplesize
  result$others<-list(dtlong=dtlong,dtsurv=dtsurv,longm=longm,survm=survm,rd=rd,id=id,samplesize=samplesize)
  class(result)<-'jmcsBig'
  result

}



