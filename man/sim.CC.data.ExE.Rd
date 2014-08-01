\name{sim.CC.data.ExE}
\alias{sim.CC.data.ExE}
\title{Simulates case and controls}
\usage{
sim.CC.data.ExE(n = NULL, ncases = NULL, ncontrols = NULL,
  max.sample.size = NULL, pheno.prev = NULL, e1.model = NULL,
  e1.prev = NULL, e1.mean = NULL, e1.sd = NULL, e1.low.lim = NULL,
  e1.up.lim = NULL, e1.OR = NULL, e2.model = NULL, e2.prev = NULL,
  e2.mean = NULL, e2.sd = NULL, e2.low.lim = NULL, e2.up.lim = NULL,
  e2.OR = NULL, i.OR = NULL, b.OR = NULL, ph.error = NULL)
}
\arguments{
  \item{n}{Number of observations to generate per
  iteration}

  \item{ncases}{Number of cases to simulate}

  \item{ncontrols}{Number of controls to simulate}

  \item{max.sample.size}{Maximum number of observations
  allowed}

  \item{pheno.prev}{Prevalence of the binary outcome}

  \item{e1.model}{Model of the 1st environmental exposure}

  \item{e1.prev}{Prevelance of the 1st environmental
  exposure}

  \item{e1.mean}{Mean of the 1st environmental exposure}

  \item{e1.sd}{Standard deviation of the 1st environmental
  exposure}

  \item{e1.low.lim}{Lower limit of the 1st environmental
  exposure}

  \item{e1.up.lim}{Upper limit of the 1st environmental
  exposure}

  \item{e1.OR}{Odds ratios of the 1st environmental
  exposure}

  \item{e2.model}{Model of the 2nd environmental exposure}

  \item{e2.prev}{Prevelance of the 2nd environmental
  exposure}

  \item{e2.mean}{Mean of the 2nd environmental exposure}

  \item{e2.sd}{Standard deviation of the 2nd environmental
  exposure}

  \item{e2.low.lim}{Lower limit of the 2nd environmental
  exposure}

  \item{e2.up.lim}{Upper limit of the 2nd environmental
  exposure}

  \item{e2.OR}{Odds ratios of the 2nd environmental
  exposure}

  \item{i.OR}{Odds ration of the interaction}

  \item{b.OR}{Baseline odds ratio for subject on 95 percent
  population centile versus 5 percentile. This parameter
  reflects the heterogeneity in disease risk arising from
  determinates that have not been measured or have not been
  included in the model.}

  \item{ph.error}{misclassification rates: 1-sensitivity
  and 1-specificity}
}
\value{
A matrix
}
\description{
Generates affected and non-affected subjects until the set
sample size is achieved.
}
\details{

}
\author{
Gaye A.
}
\keyword{internal}

