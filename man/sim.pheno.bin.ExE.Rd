\name{sim.pheno.bin.ExE}
\alias{sim.pheno.bin.ExE}
\title{Generates phenotype status}
\usage{
sim.pheno.bin.ExE(num.obs = NULL, disease.prev = NULL,
  environment1 = NULL, environment2 = NULL, interaction = NULL,
  subject.effect.data = NULL, env1.OR = NULL, env2.OR = NULL,
  int.OR = NULL)
}
\arguments{
  \item{num.obs}{Number of observations to generate per
  iteration}

  \item{disease.prev}{Prevalence of the binary outcome}

  \item{environment1}{Exposure data for the 1st
  determinant}

  \item{environment2}{Exposure data for the 2nd
  determinant}

  \item{interaction}{data}

  \item{subject.effect.data}{Subject effect data, reflects
  the heterogenity in baseline disease risk}

  \item{env1.OR}{Odds ratios of the 1st determinant}

  \item{env2.OR}{Odds ratios of the 2nd determinant}

  \item{int.OR}{Odds ration of the interaction}
}
\value{
A dataframe of phenotype
}
\description{
Generates affected and non-affected subjects
}
\details{

}
\author{
Gaye A.
}
\keyword{internal}

