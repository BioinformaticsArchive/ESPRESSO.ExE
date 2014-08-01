\name{get.critical.results.GxE}
\alias{get.critical.results.GxE}
\title{Summarizes the main results}
\usage{
get.critical.results.GxE(scenario = NULL, pheno.model = NULL,
  env1.model = NULL, env2.model = NULL, sample.sizes.required = NULL,
  empirical.power = NULL, modelled.power = NULL, mean.beta = NULL)
}
\arguments{
  \item{scenario}{Scenario number}

  \item{pheno.model}{Type of the outcome; 0 for binary and
  2 for continuous}

  \item{env1.model}{Model of the 1st enviromental
  explosure}

  \item{env2.model}{Model of the 2nd enviromental
  explosure}

  \item{sample.sizes.required}{Number of cases and controls
  or number of subjects required to achieve the desired
  power}

  \item{empirical.power}{Estimated empirical power}

  \item{modelled.power}{Calculated theoretical power}

  \item{mean.beta}{Mean beta value of each of the
  determinants}
}
\value{
A table containing the following items:
\code{genetic.model} Model of the genetic determinant
\code{environment.model} Model of the environmental
determinant \code{number.of.cases.required} Number of cases
required to achieve the desired power under binary outcome
model \code{number.of.controls.required} Number of controls
required to achieve the desired power under binary outcome
model \code{number.of.subjects.required} Number of subjects
required to achieve the desired power under a quantatative
outcome model \code{empirical.power} Estimated empirical
power under each model \code{modelled.power} Power achieved
under each model with specified sample size
\code{estimated.OR} Esitmated odds-ratios due to shrinkage
toward the null resulting from misclassification
\code{estimated.effect} Estitmated effect size if the
outocme is continuous
}
\description{
Gets the number of cases and controls or subjects and the
empirical and theoretical power under each model and prints
a summary on the screen
}
\details{

}
\author{
Gaye A.
}
\keyword{internal}

