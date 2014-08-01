\name{get.observed.data.GxE}
\alias{get.observed.data.GxE}
\title{Generates exposure data with some error}
\usage{
get.observed.data.GxE(data = NULL, e1.error = NULL, e1.model = NULL,
  e1.prev = NULL, e1.sd = NULL, e1.reliability = NULL, e2.error = NULL,
  e2.model = NULL, e2.prev = NULL, e2.sd = NULL, e2.reliability = NULL)
}
\arguments{
  \item{data}{Input table of simulated data considered as
  true data}

  \item{e1.error}{Misclassification rates in the 1st
  environmental exposure assessment: 1-sensitivity and
  1-specificity}

  \item{e1.model}{Model of the 1st exposure: binary=0,
  quantitative-normal=1 or quantitative-uniform=2}

  \item{e1.prev}{Prevalence of the 1st environmental
  exposure}

  \item{e1.sd}{Standard Deviation of the 1st environmental
  exposure}

  \item{e1.reliability}{Reliability of the assessment of
  1st environmental exposure}

  \item{e2.error}{Misclassification rates in the 2nd
  environmental exposure assessment: 1-sensitivity and
  1-specificity}

  \item{e2.model}{Model of the 2nd exposure: binary=0,
  quantitative-normal=1 or quantitative-uniform=2}

  \item{e2.prev}{Prevalence of the 2nd environmental
  exposure}

  \item{e2.sd}{Standard Deviation of the 2nd environmental
  exposure}

  \item{e2.reliability}{Reliability of the assessment of
  2nd environmental exposure}
}
\value{
A matrix
}
\description{
Uses functions make.obs.geno and make.obs.env to generate
effect data with a set level of error
}
\details{

}
\author{
Gaye A
}
\keyword{internal}

