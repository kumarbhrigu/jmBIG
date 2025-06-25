



#' @title Prediction using \code{rstanarm}
#' @description prediction of the posterior trajectory for longitudinal marker while using \code{rstanarm} for Big data
#' @param model fitted model object
#'
#' @param m m for \code{posterior_traj} function
#' @param ids value of id
#' @param ... other parameter option, see \code{posterior_traj}
#' @return list of predicted values for the given id
#' @importFrom rstanarm posterior_traj
#' @importFrom dplyr arrange desc starts_with ends_with
#' @export postTraj
#' @examples
#'  \donttest{
#' ##
#' library(survival)
#' library(dplyr)
#' fit6<-jmstanBig(dtlong=long2,dtsurv = surv2,longm=y~ x7+visit+(1|id),
#'          survm=Surv(time,status)~x1+visit,samplesize=200,time_var='visit',id='id')
#' P2<-postTraj(model<-fit6,m<-1,ids<-c(1,2,100))
#' pp1<-plot(P2$p1[[1]],plot_observed = TRUE)
#' pp2<-plot(P2$p1[[2]],plot_observed = TRUE)
#' pp3<-plot(P2$p1[[3]],plot_observed = TRUE)
#' ##
#'    }
postTraj<-function(model,m,ids,...){
  model<-model$uprlist;m<-m;ids<-ids
  nlength<-ids
  mlist<-list()
  mlonglist<-list()
  for(j in 1:length(model)){
    mlonglist[[j]]<-model[[j]]$dataEvent
  }
  f1<-function(x,id){as.numeric(id%in%x$id)}

  for(i in 1:length(ids)){
    mlist[[i]]<-sapply(mlonglist,f1,id=ids[i])
  }
  mc<-matrix(unlist(mlist),nrow=length(ids),ncol=length(model),byrow=T)
  idm<-apply(mc,1,function(x){which(x!=0)})

  p1<-list()
  for(k in 1:length(ids)){
    Mod11<-model[[idm[k]]]
    p1[[k]]<-posterior_traj(Mod11, m = 1, ids = ids[[k]],...)
  }

  P1<-Reduce('rbind',p1)
  P1<-arrange(P1,desc(P1$id))
  result<-list()
  result$P1<-P1
  result$p1<-p1
  result
}

#' @title Prediction using \code{rstanarm}
#' @description posterior survival probability estimates from rstanarm for BIG data
#' @param model fitted model
#'
#' @param ids value of id
#' @param ... other parameter option, see \code{posterior_survfit}
#' @return list of predicted value for the given id
#' @importFrom dplyr arrange desc starts_with ends_with
#' @importFrom  rstanarm posterior_survfit
#' @export postSurvfit
#' @examples
#'  \donttest{
#' ##
#' library(survival)
#' library(dplyr)
#' jmstan<-jmstanBig(dtlong=long2,
#'          dtsurv = surv2,
#'          longm=y~ x7+visit+(1|id),
#'          survm=Surv(time,status)~x1+visit,
#'          samplesize=200,
#'          time_var='visit',id='id')
#' mod1<-jmstan
#' P2<-postSurvfit(model<-mod1,ids<-c(1,2,210))
#' pp1<-plot(P2$p1[[1]])
#' pp1
#' pp2<-plot(P2$p1[[2]])
#' pp2
#' pp3<-plot(P2$p1[[3]])
#' pp3
#'  ##
#'    }
postSurvfit<-function(model,ids,...){
  model<-model$uprlist;ids<-ids
  nlength<-ids
  mlist<-list()
  mlonglist<-list()
  for(j in 1:length(model)){
    mlonglist[[j]]<-model[[j]]$dataEvent
  }
  f1<-function(x,id){as.numeric(id%in%x$id)}

  for(i in 1:length(ids)){
    mlist[[i]]<-sapply(mlonglist,f1,id=ids[i])
  }
  mc<-matrix(unlist(mlist),nrow=length(ids),ncol=length(model),byrow=T)
  idm<-apply(mc,1,function(x){which(x!=0)})

  p1<-list()
  for(k in 1:length(ids)){
    Mod11<-model[[idm[k]]]
    p1[[k]]<-posterior_survfit(Mod11, ids = ids[[k]],...)
  }

  P1<-Reduce('rbind',p1)
  P1<-arrange(P1,desc(P1$id))
  result<-list()
  result$P1<-P1
  result$p1<-p1
  result
}

#' @title  Prediction using \code{FastJM}
#' @description prediction of survival probability using \code{FastJM} for BIG data
#' @param model fitted model object
#' @param ids value of id
#' @param u see \code{survfitjmcs}
#' @param method options are 'Laplace','GH'
#' @param obs.time vector which represents time variable in the longitudinal data
#' @return list of predicted value for the given id along with other information relevant for survival probability confidence plot
#' @importFrom  FastJM survfitjmcs
#' @importFrom dplyr arrange desc starts_with ends_with
#' @export survfitJMCS
#' @examples
#'   \donttest{
#' ##
#' library(survival)
#' library(dplyr)
#' jmcs1<-jmcsBig(dtlong=data.frame(long2),
#' dtsurv = data.frame(surv2),
#' longm=y~ x7+visit,
#' survm=Surv(time,status)~x1+visit,
#' rd= ~ visit|id,
#' samplesize=200,id='id')
#' mod2<-jmcs1
#' P2<-survfitJMCS(model<-mod2,ids<-c(5),u<-seq(surv2[surv2$id==5,]$time,
#' surv2[surv2$id==5,]$time+10,0.2),obs.time='time')
#' print(P2)
#' ##
#'   }
survfitJMCS<-function(model,ids,u,method='GH',obs.time){
  model1<-model
  model<-model$uprlist;ids<-ids
  nlength<-ids
  mlist<-list()
  mlonglist<-list()
  for(j in 1:length(model)){
    mlonglist[[j]]<-model[[j]]$cdata
  }
  f1<-function(x,id){as.numeric(id%in%x$id)}

  for(i in 1:length(ids)){
    mlist[[i]]<-sapply(mlonglist,f1,id=ids[i])
  }
  mc<-matrix(unlist(mlist),nrow=length(ids),ncol=length(model),byrow=T)
  idm<-apply(mc,1,function(x){which(x!=0)})

  p1<-list()
  for(k in 1:length(ids)){
    Mod11<-model[[idm[[k]][1]]]
    ynewdata<-Mod11$ydata[Mod11$ydata$id==ids[k],]
    cnewdata<-Mod11$cdata[Mod11$cdata$id==ids[k],]
    p1[[k]]<-survfitjmcs(Mod11,
                         ynewdata=ynewdata,
                         cnewdata=cnewdata,
                         u = u,seed=100,method=method,obs.time = obs.time)
  }

  P1<-Reduce('rbind',p1)
  #P1<-arrange(P1,desc(P1$id))
  result<-list()
  result$P1<-P1
  result$p1<-p1
  result$others<-list(jmcs_others=model1,obs.time=obs.time,ids=ids)
  class(result)<-"survfitJMCS"
  result
}


#' @title Prediction using \code{JMbayes2}
#' @description prediction of survival probability and longitudinal marker using \code{jmBayes2} for BIG data
#' @param model fitted model object
#'
#' @param ids value of id
#' @param process see \code{jm}
#' @param ... other parameter options, see \code{predict.jm}
#' @param newdata dataset having covariate information for the ids mentioned above.
#' @return list of predicted value for the given id
#' @importFrom stats predict
#' @importFrom dplyr arrange desc starts_with ends_with
#' @export
#' @examples
#'
#'  \donttest{
#' ##
#' library(survival)
#' library(nlme)
#' library(dplyr)
#' jmcs1<-jmbayesBig(dtlong=long2,
#' dtsurv = surv2 ,
#' longm=y~ x7+visit,
#' survm=Surv(time,status)~x1+visit,
#' rd= ~ visit|id,
#' timeVar='visit',
#' nchain=1,
#' samplesize=200,
#' id='id')
#' mod3<-jmcs1
#' ydt<-long2%>%filter(id%in%c(900))
#' names(ydt)
#' cdt<-surv2[,'id']%>%filter(id%in%c(900))
#' names(cdt)
#' newdata<-full_join(ydt,cdt,by='id')
#' P2<-predJMbayes(model<-mod3,ids<-c(900),newdata=newdata,process = 'event')
#' plot(P2$p1[[1]])
#'
#' ##
#'
#'    }
predJMbayes<-function(model,ids,process='longitudinal',newdata,...){
  model<-model$pseudoMod;ids<-ids
  nlength<-ids
  p1<-list()
  for(k in 1:length(ids)){
    Mod11<-model
    #nd<-newdata[newdata$id==ids[k],]
    p1[[k]]<-predict(Mod11,
                     newdata=newdata[newdata$id==ids[k],],
                     return_newdata=TRUE,process=process,...)
  }

  P1<-Reduce('rbind',p1)
  #P1<-arrange(P1,desc(P1$id))
  result<-list()
  result$P1<-P1
  result$p1<-p1
  result
}

#' @title Prediction using \code{joineRML}
#' @description prediction of survival probability and longitudinal marker using \code{joineRML} for BIG data
#' @param model fitted model object
#'
#' @param ids value of id
#' @param dtlong longitudinal data
#' @param dtsurv survival data
#' @param ... other parameter options, see \code{dynSurv}
#' @return list of predicted values for the given id
#' @importFrom joineRML dynLong dynSurv
#' @importFrom dplyr arrange desc starts_with ends_with
#' @export predJRML
#' @examples
#'    \donttest{
#' ##
#' library(survival)
#' library(dplyr)
#' jmcs1<-joinRMLBig(dtlong=long2,
#' dtsurv = surv2,
#' longm=y~ x7+visit,
#' survm=Surv(time,status)~x1+visit,
#' rd=~ visit|id,
#' timeVar='visit',
#' samplesize=200,
#' id='id')
#' mod4<-jmcs1
#' P2<-predJRML(model<-mod4,ids<-c(10),dtlong=long2,dtsurv=surv2)
#' plot(P2$plong[[1]])
#' plot(P2$psurv[[1]])
#' ##
#'    }
predJRML<-function(model,ids,dtlong,dtsurv,...){
  model<-model;ids<-ids
  nlength<-ids
  mlist<-list()

  p1<-list();p2<-list()
  for(k in 1:length(ids)){
    Mod11<-model$pseudoMod
    ydata<-dtlong[dtlong$id==ids[k],]
    cdata<-dtsurv[dtsurv$id==ids[k],]
    #cdata<-cdata[names(cdata)%in%c('id','Time','status')]
    #ND<-full_join(ydata,cdata,by='id')
    p1[[k]]<-dynLong(Mod11,
                     newdata=ydata)
    p2[[k]]<-dynSurv(Mod11,newdata=ydata,newSurvData=cdata,type='simulated',...)
  }

  P1<-Reduce('rbind',p1)
  P2<-Reduce('rbind',p2)
  #P1<-arrange(P1,desc(P1$id))
  result<-list()
  result$P1<-P1
  result$plong<-p1
  result$psurv<-p2
  result
}



#' @title Bootstrapped CI using \code{FastJM}
#' @description Bootstrapped CI for predicted survival probability
#' @param object a \code{survfitJMCS} object
#'
#' @return Bootstrap CI for the survival probability and other relevant information for predicted survival plot
#' @importFrom stats quantile
#' @export
#'
#' @examples
#'   \donttest{
#' ##
#' library(survival)
#' library(dplyr)
#' jmcs1<-jmcsBig(dtlong=data.frame(long2),
#' dtsurv = data.frame(surv2),
#' longm=y~ x7+visit,
#' survm=Surv(time,status)~x1+visit,
#' rd= ~ visit|id,
#' samplesize=200,id='id')
#' mod2<-jmcs1
#' P2<-survfitJMCS(model<-mod2,ids<-c(5),u<-seq(surv2[surv2$id==5,]$time,
#' surv2[surv2$id==5,]$time+10,0.2),obs.time='time')
#' bootci<-cisurvfitJMCS(P2)
#' print(bootci)
#' ##
#'   }
cisurvfitJMCS<-function(object){
  if(!inherits(object,"survfitJMCS"))
    stop("\n Not a 'survfitJMCS' object.\n")
  longdata<-object$others$jmcs_others$others$dtlong
  survdata<-object$others$jmcs_others$others$dtsurv
  id<-object$others$jmcs_others$others$id
  idNumber<-object$others$ids
  jmcsModel<-object$others$jmcs_others
  obs.time<-object$others$obs.time
  longm<-object$others$jmcs_others$others$longm
  survm<-object$others$jmcs_others$others$survm
  rd<-object$others$jmcs_others$others$rd
  samplesize<-object$others$jmcs_others$samplesize
  lastid<-max(survdata[[id]])



  bootstrap_longitudinal_survival <- function(longitudinal_data, survival_data, n_bootstrap = 10,id,idNumber){
    bootstrap_samples <- vector("list", length = n_bootstrap)
    unique_ids <- as.numeric(unique(longitudinal_data[[id]]))
    for (i in 1:n_bootstrap) {
      # Sample IDs with replacement
      bootstrap_ids <- sample(unique_ids[unique_ids!=idNumber], replace = TRUE)
      n_bootstrap_ids<-seq(1,length(unique(survival_data[[id]])))
      bootstrap_longitudinal_sample <- data.frame()
      bootstrap_survival_sample <- data.frame()
      count<-0

      for (j in 1:(length(bootstrap_ids))) {
        id_longitudinal_data <- longitudinal_data[longitudinal_data[[id]] == bootstrap_ids[j], ]
        id_longitudinal_data[[id]]<-j
        id_survival_data <- survival_data[survival_data[[id]] == bootstrap_ids[[j]], ]
        id_survival_data[[id]]<-j
        bootstrap_longitudinal_sample <- rbind(bootstrap_longitudinal_sample, id_longitudinal_data)
        bootstrap_survival_sample <- rbind(bootstrap_survival_sample, id_survival_data)
      }
      newdatalong<-longitudinal_data[longitudinal_data[id]==idNumber,]
      newdatalong[id]<-length(unique(longitudinal_data[[id]]))
      bootstrap_longitudinal_sample<-rbind(bootstrap_longitudinal_sample,newdatalong)
      newdatasurv<-survival_data[survival_data[id]==idNumber,]
      newdatasurv[id]<-length(unique(longitudinal_data[[id]]))
      bootstrap_survival_sample<-rbind(bootstrap_survival_sample,newdatasurv)

      bootstrap_samples[[i]] <- list(longitudinal = bootstrap_longitudinal_sample, survival = bootstrap_survival_sample)
    }
    return(bootstrap_samples)
  }
  bootstrapped_data <- bootstrap_longitudinal_survival(longdata,
                                                       survdata,
                                                       n_bootstrap = 10,id=id,
                                                       idNumber = idNumber)

  fit2<-list();P2<-list();P3<-list()
  for(i in 1:10){
    P2<-try(survfitJMCS(model=jmcsModel,ids=idNumber,
                     u<-seq(survdata[survdata[id]==idNumber,][[obs.time]],survdata[survdata[id]==idNumber,][[obs.time]]+10,0.2),
                     obs.time=obs.time),silent=TRUE)
    fit2[[i]]<-try(jmcsBig(dtlong=data.frame(bootstrapped_data[[i]]$longitudinal),
                        dtsurv = data.frame(bootstrapped_data[[i]]$survival),
                        longm=longm,survm=survm,
                        rd= rd,samplesize=200,id=id),silent=TRUE)

    P2[[i]]<-try(survfitJMCS(model<-fit2[[i]],ids<-c(lastid),
                          u<-seq(bootstrapped_data[[i]]$survival[bootstrapped_data[[i]]$survival[[id]]==lastid,][[obs.time]],
                                 bootstrapped_data[[i]]$survival[bootstrapped_data[[i]]$survival[[id]]==lastid,][[obs.time]]+10,0.2),
                          obs.time=obs.time),silent=TRUE)
    P3[[i]]<-try(P2[[i]]$P1$Pred[[as.character(lastid)]][,2],silent=TRUE)

  }

  CIdata<-Reduce('cbind',P3)
  qCIdata<-data.frame(Times=P2[[i]]$P1$Pred[[as.character(lastid)]]$times,LL=apply(CIdata,1,function(x){quantile(x,0.025)}),Med=apply(CIdata,1,function(x){quantile(x,0.5)}),UL=apply(CIdata,1,function(x){quantile(x,0.975)}))
  result<-list()
  result$lastid<-lastid
  result$bootCI<-qCIdata
  result$P2<-P2
  result$P3<-P3
  class(result)<-"cisurvfitJMCS"
  result

}


#' @title  Plot for \code{cisurvfitJMCS} object
#' @description prediction of survival probability and longitudinal marker using \code{FastJM} for BIG data
#' @param object fitted \code{survfitJMCS} object
#' @return Plot for predicted survival probability
#' @import ggplot2
#' @export
#' @examples
#'   \donttest{
#' ##
#' library(survival)
#' library(dplyr)
#' jmcs1<-jmcsBig(dtlong=data.frame(long2),
#' dtsurv = data.frame(surv2),
#' longm=y~ x7+visit,
#' survm=Surv(time,status)~x1+visit,
#' rd= ~ visit|id,
#' samplesize=200,id='id')
#' mod2<-jmcs1
#' P2<-survfitJMCS(model<-mod2,ids<-c(5),u<-seq(surv2[surv2$id==5,]$time,
#' surv2[surv2$id==5,]$time+10,0.2),obs.time='time')
#' P3<-cisurvfitJMCS(P2)
#' plot_cisurvfitJMCS(P3)
#' ##
#'   }
plot_cisurvfitJMCS<-function(object){
  object<-object
  if(!inherits(object,"cisurvfitJMCS"))
    stop("\n Not a 'cisurvfitJMCS' object.\n")
  P3<-object$P3
  P2<-object$P2
  lastid<-object$lastid
  qCIdata<-object$bootCI
  #CIdata<-Reduce('cbind',P3)
  #qCIdata<-data.frame(Times=P2[[1]]$P1$Pred[[as.character(lastid)]]$times,LL=apply(CIdata,1,function(x){quantile(x,0.025)}),Med=apply(CIdata,1,function(x){quantile(x,0.5)}),UL=apply(CIdata,1,function(x){quantile(x,0.975)}))
  survPlot<-ggplot(qCIdata, aes(x = Times, y = Med)) +
    geom_line() +
    geom_ribbon(aes(ymin = LL, ymax = UL), alpha = 0.3,fill='blue') +
    labs(title = "Survival Probabilities with Prediction Intervals",
         x = "Time",
         y = "Survival Probability")

  survPlot
  #invisible(object)
}


utils::globalVariables(c('ggplot',
                         'aes','Times',
                         'geom_line','Med','geom_line',
                         'geom_ribbon','LL','UL','labs',''))
