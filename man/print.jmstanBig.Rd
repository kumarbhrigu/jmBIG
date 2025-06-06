% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.R
\name{print.jmstanBig}
\alias{print.jmstanBig}
\title{print.jmstanBig}
\usage{
\method{print}{jmstanBig}(object, digits = 3, ...)
}
\arguments{
\item{object}{object}

\item{digits}{used for round the numeric values after decimal}

\item{...}{others}
}
\value{
prints table containing various parameter estimates,
        SE, P- value for both survival and longitudinal submodel,
        if the model is bayesian it includes their credible interval too.
}
\description{
print.jmstanBig
}
\examples{

 \donttest{
##
library(survival)
library(dplyr)
mod1<-jmstanBig(dtlong=long2,
         dtsurv = surv2,
         longm=y~ x7+visit+(1|id),
         survm=Surv(time,status)~x1+visit,
         samplesize=200,
         time_var='visit',id='id')
print(mod1)
}
}
