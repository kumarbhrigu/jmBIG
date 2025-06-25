

#' @title  print.jmstanBig
#' @description
#' print method for class 'jmstanBig'
#'
#' @param x fitted object
#' @param ... others
#'
#' @return prints table containing various parameter estimates,
#'         SE, P- value for both survival and longitudinal submodel,
#'         if the model is bayesian it includes their credible interval too.
#' @importFrom rstanarm VarCorr
#' @export
#'
#' @examples
#'
#'  \donttest{
#' ##
#' library(survival)
#' library(dplyr)
#' mod1<-jmstanBig(dtlong=long2,
#'          dtsurv = surv2,
#'          longm=y~ x7+visit+(1|id),
#'          survm=Surv(time,status)~x1+visit,
#'          samplesize=200,
#'          time_var='visit',id='id')
#' print(mod1)
#' }
#' @method print jmstanBig
print.jmstanBig<-function(x,...){
  #x<-object
  digits<-3
  if(!inherits(x,'jmstanBig'))
    stop("\n Not a 'jmstanBig' object.\n")

  cat("\n Joint model for Big data using rstanarm")
  cat("\n Call: \n",paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("\n Total observation in longitudinal data:" ,x$nr ,'\n')
  cat("\n Chunk size:", x$samplesize, "\n")

  cat('\n Longitudinal process: \n')
  vname<-rownames(attr(x$pseudoMod$terms$Long1,'factors'))
  vname<-vname[vname!='id']
  vname[1]<-'Intercept'
  lstan<-x$pseudoMod$stan_summary[starts_with('Long1',vars =rownames(x$pseudoMod$stan_summary)),c(1,3,4,10)]
  lstan<-cbind(lstan,Zvalue=lstan[,1]/lstan[,2])
  f1<-function(x){
    return(if((x[1]>0))2*(1-pnorm((x[1]/x[2]),mean=0,sd=1))else 2*(pnorm((x[1]/x[2]),mean=0,sd=1)))
  }
  attr(lstan,'dimnames')[[2]]<-c('Mean','StDev','2.5%','97.5%','Zvalue')
  Pvalue<-apply(lstan,1,f1)
  lstan<-cbind(lstan,Pvalue)
  row.names(lstan)<-gsub("Long1\\|","",row.names(lstan))
  #dat<-data.frame(dat,row.names = c(names(x$pseudoMod$statistics$Mean$betas1),'sigma'))
  lstan<-round(lstan,digits=digits)
  print(lstan)
  cat('\n Survival process: \n ')

  sdat<-list()
  sstan<-x$pseudoMod$stan_summary[c(starts_with('Event',vars =rownames(x$pseudoMod$stan_summary)),
                                    ends_with('etavalue',vars =rownames(x$pseudoMod$stan_summary)))
                                  ,c(1,3,4,10)]
  sstan<-cbind(sstan,Zvalue=sstan[,1]/sstan[,2])
  attr(sstan,'dimnames')[[2]]<-c('Mean','StDev','2.5%','97.5%','Zvalue')
  sstan<-cbind(sstan,Pvalue=apply(sstan,1,f1))
  row.names(sstan)<-gsub('Event\\|','',row.names(sstan))
  sstan<-round(sstan,digits=digits)
  print(sstan)

  cat('Random effects covariance matrix:\n')
  D<-VarCorr(x$pseudoMod)
  #D<-round(D,digits=digits)
  print(D)

  invisible(x)
}



#' @title  print.jmcsBig
#' @description
#' print method for class 'jmcsBig'
#'
#' @param x fitted object
#' @param ... others
#'
#' @return prints table containing various parameter estimates,
#'         SE, P- value for both survival and longitudinal submodel,
#'         if the model is bayesian it includes their credible interval too.
#' @export
#'
#' @examples
#'  \donttest{
#' ##
#' library(survival)
#' library(dplyr)

#' ################################
#' mod2<-jmcsBig(dtlong=data.frame(long2),
#' dtsurv = data.frame(surv2),
#' longm=y~ x7+visit,
#' survm=Surv(time,status)~x1+visit,
#' rd= ~ visit|id,
#' samplesize=200,id='id')
#' print(mod2)
#'    }
#' @method print jmcsBig
print.jmcsBig<-function(x,...){
  #x<-object
  digits<-3
  if(!inherits(x,'jmcsBig'))
    stop("\n Not a 'jmcsBig' object.\n")
  cat("\n Joint model for Big data using FastJM")
  cat("\n Call: \n",paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("\n Total observation in longitudinal data:" ,x$nr ,'\n')
  cat("\n Chunk size:", x$samplesize, "\n")
  #cat('Random effects covariance matrix:\n')
  #D<-x$pseudoMod$D
  cat('\n Longitudinal process: \n')
  ldat<-list()
  ldat$Estimate<-c(x$pseudoMod$beta,x$pseudoMod$sigma)
  ldat$SE<-c(x$pseudoMod$sebeta,x$pseudoMod$sesigma)
  ldat$Zvalue<-ldat$Estimate/ldat$SE
  ldat<-data.frame(ldat)

  f1<-function(x){
    return(if((x[3]>0))2*(1-pnorm((x[3]),mean=0,sd=1))else 2*(pnorm((x[3]),mean=0,sd=1)))
  }

  Pvalue<-apply(ldat,1,f1)
  ldat<-cbind(ldat,Pvalue)
  row.names(ldat)[nrow(ldat)]<-'sigma^2'
  #dat<-data.frame(dat,row.names = c(names(x$pseudoMod$statistics$Mean$betas1),'sigma'))
  ldat<-round(ldat,digits=digits)
  print(ldat)
  cat('\n Survival process: \n ')

  sdat<-list()
  sdat$Estimate<-x$pseudoMod$gamma1
  sdat$SE<-x$pseudoMod$segamma1
  sdat$ZValue<-sdat$Estimate/sdat$SE
  sdat<-data.frame(sdat)
  sPvalue<-apply(sdat,1,f1)
  sdat<-cbind(sdat,Pvalue=sPvalue)
  sdat<-round(sdat,digits=digits)
  print(sdat)

  cat('\n Association parameters :\n')
  adat<-list()
  adat$Estimate<-x$pseudoMod$nu1
  adat$SE<-x$pseudoMod$senu1
  adat$Zvalue<-adat$Estimate/adat$SE
  adat<-data.frame(adat)
  aPvalue<-apply(adat,1,f1)
  adat<-cbind(adat,Pvalue=aPvalue)
  adat<-round(adat,digits = digits)
  random<-all.vars(x$pseudoMod$random)
  #if (length(x$pseudoMod$nu1) == 1) rownames(adat) <- c("(Intercept)_1", "(Intercept)_2")
  if (length(x$pseudoMod$nu1) == 2) rownames(adat) <- c("(Intercept)_1", paste0(random[1], "_1"))
  #if (length(x$pseudoMod$nu1) == 3) rownames(adat) <- c("(Intercept)_1", paste0(random[1], "_1"), paste0(random[2], "_1"),
  # "(Intercept)_2", paste0(random[1], "_2"), paste0(random[2], "_2"))
  adat<-round(adat,digits=digits)
  print(adat)

  cat('\n Variance Covariance matrix of Random effects:\n')
  rdat<-x$pseudoMod$Sig
  if(dim(rdat)[[2]]==2){
    rvar<-all.vars(x$pseudoMod$random)
    rownames(rdat)<-c('Intercept',paste0(rvar[[1]]))
    colnames(rdat)<-c('Intercept',paste0(rvar[[1]]))
  }
  rdat<-round(rdat,digits=digits)
  print(rdat)
  invisible(x)
}


#' @title  print.jmbayesBig
#' @description
#' print method for class 'jmbayesBig'
#'
#' @param x fitted object
#' @param ... others
#'
#' @return prints table containing various parameter estimates,
#'         SE, P- value for both survival and longitudinal submodel,
#'         if the model is bayesian it includes their credible interval too.
#' @export
#'
#' @examples
#'  \donttest{
#' ##
#' library(survival)
#' library(dplyr)
#'
#' #################################
#' mod3<-jmbayesBig(dtlong=long2,
#' dtsurv = surv2 ,
#' longm=y~ x7+visit,
#' survm=Surv(time,status)~x1+visit,
#' rd= ~ visit|id,
#' timeVar='visit',
#' nchain=1,
#' samplesize=200,
#' id='id')
#' print(mod3)
#'    }
#' @method print jmbayesBig
print.jmbayesBig<-function(x,...){
  #x<-object
  digits<-3
  if(!inherits(x,'jmbayesBig'))
    stop("\n Not a 'jmbayesBig' object.\n")

  cat("\n Joint model for Big data using jmbayes2")
  cat("\n Call: \n",paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("\n Total observation in longitudinal data:" ,x$nr ,'\n')
  cat("\n Chunk size:", x$samplesize, "\n")

  cat('\n Longitudinal process: \n')
  dat<-list()
  dat$Mean<-c(x$pseudoMod$statistics$Mean$betas1, x$pseudoMod$statistics$Mean$sigmas)
  dat$StDev<-c(x$pseudoMod$statistics$SD$betas1,x$pseudoMod$statistics$SD$sigmas)
  dat$'2.5%'<-c(x$pseudoMod$statistics$CI_low$betas1,x$pseudoMod$statistics$CI_low$sigmas)
  dat$'97.5%'<-c(x$pseudoMod$statistics$CI_upp$betas1,x$pseudoMod$statistics$CI_upp$sigmas)
  f1<-function(x){
    return(if((x[1]>0))2*(1-pnorm((x[1]/x[2]),mean=0,sd=1))else 2*(pnorm((x[1]/x[2]),mean=0,sd=1)))
  }

  dat$Pvalue<-apply(data.frame(dat$Mean,dat$StDev),1,f1)
  names(dat)[c(3,4)]<-c(as.character("2.5%"),as.character("97.5%"))
  dat<-data.frame(dat,row.names = c(names(x$pseudoMod$statistics$Mean$betas1),'sigma'),check.names = F)
  dat<-round(dat,digits=digits)
  print(dat)
  cat('\n Survival process: \n ')

  sdat<-list()
  sdat$Mean<-c(x$pseudoMod$statistics$Mean$gammas,x$pseudoMod$statistics$Mean$alphas)
  sdat$StDev<-c(x$pseudoMod$statistics$SD$gammas,x$pseudoMod$statistics$SD$alphas)
  sdat$'2.5%'<-c(x$pseudoMod$statistics$CI_low$gammas,x$pseudoMod$statistics$CI_low$alphas)
  sdat$'97.5%'<-c(x$pseudoMod$statistics$CI_upp$gammas,x$pseudoMod$statistics$CI_upp$alphas)
  sdat$Pvalue<-apply(data.frame(sdat$Mean,sdat$StDev),1,f1)
  names(sdat)[c(3,4)]<-c(as.character("2.5%"),as.character("97.5%"))
  sdat<-data.frame(sdat,row.names=c(names(x$pseudoMod$statistics$Mean$gammas),names(x$pseudoMod$statistics$Mean$alphas)),check.names = F)
  sdat<-round(sdat,digits=digits)
  print(sdat)
  cat('\n Random effects covariance matrix:\n')
  D<-x$pseudoMod$statistics$Mean$D
  #D<-round(D,digits=digits)
  print(D)

  invisible(x)
}

#' @title  print.joinRMLBig
#' @description
#' print method for class 'joinRMLBig'
#'
#' @param x fitted object
#' @param ... others
#'
#' @return prints table containing various parameter estimates,
#'         SE, P- value for both survival and longitudinal submodel,
#'         if the model is bayesian it includes their credible interval too.
#' @importFrom stats vcov
#' @export
#'
#' @examples
#'  \donttest{
#' ##
#' library(survival)
#' library(dplyr)
#' mod4<-joinRMLBig(dtlong=long2,
#' dtsurv = surv2,
#' longm=y~ x7+visit,
#' survm=Surv(time,status)~x1+visit,
#' rd=~ visit|id,
#' timeVar='visit',
#' samplesize=200,
#' id='id')
#' print(mod4)
#'    }
#' @method print joinRMLBig
print.joinRMLBig<-function(x,...){
  #x<-object
  digits<-3
  if(!inherits(x,'joinRMLBig'))
    stop("\n Not a 'joinRMLBig' object.\n")
  cat("\n Joint model for Big data using joineRML")
  cat("\n Call: \n",paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("\n Total observation in longitudinal data:" ,x$nr ,'\n')
  cat("\n Chunk size:", x$samplesize, "\n")
  #cat('Random effects covariance matrix:\n')
  #D<-x$pseudoMod$D
  vc<-vcov(x$pseudoMod)
  nd<-sum(1:dim(x$pseudoMod$coefficients$D)[[2]])
  nb<-(nd+1):(nd+length(x$pseudoMod$coefficients$beta)+1)
  ng<-(nd+length(x$pseudoMod$coefficients$beta)+1+1):dim(vc)[[2]]
  ese<-diag(vc)
  cat('\n Longitudinal process: \n')
  ldat<-list()
  ldat$Estimate<-c(x$pseudoMod$coefficients$beta,x$pseudoMod$coefficients$sigma2)
  ldat$SE<-c(ese[nb])
  ldat$Zvalue<-ldat$Estimate/ldat$SE
  ldat<-data.frame(ldat)

  f1<-function(x){
    return(if((x[3]>0))2*(1-pnorm((x[3]),mean=0,sd=1))else 2*(pnorm((x[3]),mean=0,sd=1)))
  }

  Pvalue<-apply(ldat,1,f1)
  ldat<-cbind(ldat,Pvalue)
  #row.names(ldat)<-
  #dat<-data.frame(dat,row.names = c(names(x$pseudoMod$statistics$Mean$betas1),'sigma'))
  ldat<-round(ldat,digits=digits)
  print(ldat)
  cat('\n Survival process: \n ')

  sdat<-list()
  sdat$Estimate<-x$pseudoMod$coefficients$gamma
  sdat$SE<-ese[ng]
  sdat$ZValue<-sdat$Estimate/sdat$SE
  sdat<-data.frame(sdat)
  sPvalue<-apply(sdat,1,f1)
  sdat<-cbind(sdat,Pvalue=sPvalue)
  sdat<-round(sdat,digits=digits)
  print(sdat)


  cat('\n Variance Covariance matrix of Random effects:\n')

  rdat<-x$pseudoMod$coefficients$D
  if(dim(rdat)[[2]]==2){
    rvar<-all.vars(x$pseudoMod$formLongRandom[[1]])
    rownames(rdat)<-c('Intercept',paste0(rvar[[1]]))
    colnames(rdat)<-c('Intercept',paste0(rvar[[1]]))
  }
  rdat<-round(rdat,digits=digits)
  print(rdat)
  invisible(x)
}


#' @export
#' @method print survfitJMCS
print.survfitJMCS<-function(x,...){
  object<-x
  if(!inherits(object,"survfitJMCS"))
    stop("\n Not a 'survfitJMCS' object.\n")

  print(object$P1)
  invisible(object)
}

#' @export
#' @method print cisurvfitJMCS
print.cisurvfitJMCS<-function(x,...){
  object<-x
  if(!inherits(object,"cisurvfitJMCS"))
    stop("\n Not a 'cisurvfitJMCS' object.\n")
  cat("\n Predicted survival proability data\n")
  print(object$bootCI)
  invisible(object)
}
