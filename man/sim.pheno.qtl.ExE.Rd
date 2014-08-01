\name{sim.pheno.qtl.ExE}
\alias{sim.pheno.qtl.ExE}
\title{Simulates continuous outcome data}
\usage{
sim.pheno.qtl.ExE(numsubjects = NULL, pheno.mean = NULL, pheno.sd = NULL,
  environment1 = NULL, env1.efkt = NULL, environment2 = NULL,
  env2.efkt = NULL, interaction = NULL, int.efkt = NULL)
}
\arguments{
  \item{num.subjects}{Number of subjects to simulate}

  \item{pheno.mean}{statistical mean}

  \item{pheno.sd}{standard deviation}

  \item{environment1}{Exposure data for the first
  determinant}

  \item{env1.efkt}{Effects of the first determinant}

  \item{environment2}{Exposure data for the second
  determinant}

  \item{env2.efkt}{Effects of the second determinant}

  \item{interaction}{data}

  \item{int.efkt}{Interaction effect}
}
\value{
A dataframe of phenotype
}
\description{
Uses the effects data of the determinants to construct a
linear predictor(LP). The outcome is normally distributed
variable generated with a mean equal to LP and a standard
deviation of 1. Some error is then added to the simulated
outcome to obtained the observed outcome.
}
\details{

}
\author{
Gaye A.
}
\keyword{internal}

