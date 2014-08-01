\name{env1.params}
\alias{env1.params}
\docType{data}
\title{Parameters to generate environmental exposure data}
\description{A table that contains parameters to set the environmental exposure}
\usage{data(env1.params)}
\format{
  A data frame with 4 observations on the following 12 variables.
  \describe{
    \item{\code{scenario.id}}{scenario ID}
    \item{\code{env1.model}}{model of the environmental exposure; 0 for binary and 1 for continuous}
    \item{\code{env1.prevalence}}{prevelance of the 'at risk' environmental exposure}
    \item{\code{env1.OR}}{odds ratio related to the environmental exposure (when exposure is binary)}
    \item{\code{env1.efkt}}{effect size related to the environmental exposure (when exposure is quantitative)}
    \item{\code{env1.mean}}{statistical mean of the normally distributed environment exposure}
    \item{\code{env1.sd}}{standard deviation of the normally distributed environment exposure}
    \item{\code{env1.low.lim}}{Lower limit of the uniformly distributed environment exposure}
    \item{\code{env1.up.lim}}{Upper limit of the uniformly distributed environment exposure}
    \item{\code{env1.sensitivity}}{sensitivity of the assessment of the binary exposure}
    \item{\code{env1.specificity}}{specificity of the assessment of the binary exposure}
    \item{\code{env1.reliability}}{reliability of the assessment of the quantative exposure}
  }
}

\examples{data(env1.params)}
