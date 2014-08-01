\name{env2.params}
\alias{env2.params}
\docType{data}
\title{Parameters to generate environmental exposure data}
\description{A table that contains parameters to set the environmental exposure}
\usage{data(env2.params)}
\format{
  A data frame with 4 observations on the following 12 variables.
  \describe{
    \item{\code{scenario.id}}{scenario ID}
    \item{\code{env2.model}}{model of the environmental exposure; 0 for binary and 1 for continuous}
    \item{\code{env2.prevalence}}{prevelance of the 'at risk' environmental exposure}
    \item{\code{env2.OR}}{odds ratio related to the environmental exposure (when exposure is binary)}
    \item{\code{env2.efkt}}{effect size related to the environmental exposure (when exposure is quantitative)}
    \item{\code{env2.mean}}{statistical mean of the normally distributed environment exposure}
    \item{\code{env2.sd}}{standard deviation of the normally distributed environment exposure}
    \item{\code{env2.low.lim}}{Lower limit of the uniformly distributed environment exposure}
    \item{\code{env2.up.lim}}{Upper limit of the uniformly distributed environment exposure}
    \item{\code{env2.sensitivity}}{sensitivity of the assessment of the binary exposure}
    \item{\code{env2.specificity}}{specificity of the assessment of the binary exposure}
    \item{\code{env2.reliability}}{reliability of the assessment of the quantative exposure}
  }
}

\examples{data(env2.params)}
