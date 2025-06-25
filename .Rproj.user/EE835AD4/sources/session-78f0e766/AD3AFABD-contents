#' @title Joint model for BIG data using joineRML
#' @description function for joint model in BIG DATA using \code{joineRML}
#' @param dtlong longitudinal dataset, which contains id,visit time,longitudinal measurements along with various covariates
#' @param dtsurv survival dataset corresponding to the longitudinal dataset, with survival status and survival time
#' @param longm model for longitudinal response
#' @param survm survival model
#' @param samplesize random effect part
#' @param rd random effect part
#' @param timeVar time variable in longitudinal model, included in the longitudinal data
#' @param id name of id column in longitudinal dataset
#' @return returns a list containing various output which are useful for prediction.
#' @importFrom joineRML mjoint
#' @importFrom stats pnorm
#' @export
#' @references Hickey, Graeme L., et al. "joineRML: a joint model and software package for time-to-event and multivariate longitudinal outcomes." BMC medical research methodology 18 (2018): 1-14.
#' @examples
#'   \donttest{
#' ##
#' library(survival)
#' library(dplyr)
#' fit4<-joinRMLBig(dtlong=long2,dtsurv = surv2,longm=y~ x7+visit,survm=Surv(time,status)~x1+visit,
#' rd=~ visit|id,timeVar='visit',samplesize=200,id='id')
#' P2<-predJRML(model<-fit4,ids<-c(10),dtlong=long2,dtsurv=surv2)
#' pp1<-plot(P2$plong[[1]])
#' pp1<-plot(P2$psurv[[1]])
#' ##
#'    }
#'
#'
#' @author Atanu Bhattacharjee, Bhrigu Kumar Rajbongshi and Gajendra Kumar Vishwakarma
#' @seealso   \link{jmbayesBig},\link{jmstanBig},\link{jmcsBig}

joinRMLBig<-function(dtlong,dtsurv,longm,survm,samplesize=50,rd,timeVar,id){

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
  updatebeta<-0;updategamma<-0;updateD<-0;updateH<-0
  Dlist<-list();Hlist<-list()
  for(i in 1:ceiling(nr/samplesize)){#floor(nr/samplesize)

    mod1 <- mjoint(formLongFixed=longm,formLongRandom=rd,formSurv=survm,data=dtlong1[[i]]
                   ,survData=dtsurv1[[i]],timeVar= timeVar)
    rlist[[i]]<-mod1
    #uprlist[[i]]<-summary(mod1, probs = c(.025,.975))
    # here we have considered only the estimate part fro the mod1 object
    betalist[[i]]<-mod1$coefficients$beta
    gammalist[[i]]<-mod1$coefficients$gamma
    Dlist[[i]]<-mod1$coefficients$D
    Hlist[[i]]<-mod1$Hessian

    if(i==1){updatebeta<-betalist[[1]]; updategamma<-gammalist[[1]]; updateD<-Dlist[[1]];updateH<-Hlist[[1]]}
    updatebeta<-suppressWarnings((updatebeta+betalist[[i]])/2)
    #updategamma<-0;updatenu<-0
    updategamma<-suppressWarnings((updategamma+gammalist[[i]])/2)
    updateD<-(updateD+Dlist[[i]])/2
    updateH<-(updateH+Hlist[[i]])/2

  }

  # update in mod1 object

  uprlist1<-rlist
  for(i in 1:length(rlist)){
    uprlist1[[i]]$coefficients$beta<-updatebeta
    uprlist1[[i]]$coefficients$gamma<-updategamma
    uprlist1[[i]]$coefficients$D<-updateD
    uprlist1[[i]]$Hessian<-updateH

  }

  mod11<-NULL
  mod11<-mod1
  #mlist<-updatestimate
  mod11$coefficients$beta<-updatebeta
  mod11$coefficients$gamma<-updategamma
  mod11$coefficients$D<-updateD
  mod11$Hessian<-updateH

  #Results
  result<-list()
  result$call<-cl
  result$allmodel<-rlist
  result$uprlist<-uprlist1
  result$pseudoMod<-mod11
  result$nr<-nr
  result$samplesize<-samplesize
  class(result)<-'joinRMLBig'
  #message('Results for the joint model')
  result
}



