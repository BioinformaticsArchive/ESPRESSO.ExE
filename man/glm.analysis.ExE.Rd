\name{glm.analysis.ExE}
\alias{glm.analysis.ExE}
\title{Carries out regression analysis}
\usage{
glm.analysis.ExE(pheno.model = NULL, observed.data = NULL)
}
\arguments{
  \item{pheno.model}{Type of outcome; 0=binary and
  1=continuous}

  \item{observed.data}{A dataframe that contains covariates
  and outcome data}
}
\value{
A vector containing the beta, standard-error and
z-statistic of each of the covariates
}
\description{
Fits a conventional unconditional logistic regression model
with a binary or continuous phentype as outcome and the
genetic, environmental, interaction determinants as
covariates.
}
\details{

}
\author{
Gaye A.
}
\keyword{internal}

